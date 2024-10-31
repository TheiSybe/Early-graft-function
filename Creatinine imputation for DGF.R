#### Select data to impute

NOTR_creat_1 <-subset(NOTR_DGF.df,Typeofdonor=='Deceased' 
                      & Recipientage >17
                      & Donorage >17
                      & Number_of_transplantations==1
                      & (Combinedtransplants=='LKi'
                         | Combinedtransplants=='RKi') 
                      & PNF== 0
                      & is.na(Transatl)
                      & PT_Earlydeath ==0
                      & !is.na(dDGF)
                      & Year_Of_Tx > 2013
                      & Year_Of_Tx < 2023) %>% 
  select( c("ID", "dDGF", "fictiefrecipientnummer",  "Follow_up_time",
            'Recipientage','Recipientsex','CreatinineD1', 'CreatinineD2',
            'CreatinineD3','CreatinineD4','CreatinineD5','CreatinineD6','CreatinineD7',
            "Recipient_weight" ,  "HLA_mismatch", 
            "Recipient_height" , 
            "Recipient_Diabetes", 
            "Warm_ischaemic_period_2" ,"DonorDeathGroup.dbd",
            "Cold_ischaemic_period" ,  "status2", 
            "Recipientsex",
            "Donor_type_deceased",
            "Donorsex", "Donorage", 
            "Smoking",	"Donorheight", "Donorweight",
            "Lastcreatinineµmoll", "Year_Of_Tx"))

## Identify creatinine columns
creat.cols <- c("CreatinineD1", "CreatinineD2", "CreatinineD3", "CreatinineD4",
                "CreatinineD5", "CreatinineD6", "CreatinineD7")

## visualize the amount of missing data
vis_miss(NOTR_creat_1 %>% select(-creat.cols, -fictiefrecipientnummer,-Follow_up_time, -status2,  -ID))
gg_miss_fct(subset(NOTR_creat_1, Donor_type_deceased =="DCD") %>% select(-creat.cols, -fictiefrecipientnummer,-Follow_up_time, -status2,  -ID),Year_Of_Tx)



# Subset to only non-DGF cases because that is what we are actually interested in
NOTR_creat <- subset(NOTR_creat_1, dDGF ==0)

# We had some discrepancies which were related to the year of data registration
gg_miss_fct(NOTR_creat %>% select(creat.cols,Year_Of_Tx), Year_Of_Tx)

# Check the missings. 
length(NOTR_creat[(rowSums(is.na(NOTR_creat[creat.cols]))>2),]$ID)

# We are only going to impute creatinine in cases that have at least two recorded values for creatinine (otherwise the slope might be completely random)

NOTR_creat<- NOTR_creat[rowSums(is.na(NOTR_creat[creat.cols])) <= 5 ,]
sum(is.na(gather(NOTR_creat, day, creatinine, CreatinineD1:CreatinineD7)$creatinine))
gg_miss_fct(NOTR_creat %>% select(creat.cols,Year_Of_Tx), Year_Of_Tx)

## Change ID variables character -> numeric to avoid factoring. Factored ID cols will make mice think they need to be used in imputation
NOTR_creat$ID <- as.numeric(NOTR_creat$ID)
NOTR_creat$fictiefrecipientnummer <- as.numeric(NOTR_creat$fictiefrecipientnummer)


# Identify the other character variables in the data frame
character_variables <- sapply(NOTR_creat, is.character)

# Convert the character variables to factors
NOTR_creat[, character_variables] <- lapply(NOTR_creat[, character_variables], as.factor)


### To impute
init = mice(NOTR_creat, maxit=0) 
meth = init$method
pred = mice::make.predictorMatrix(NOTR_creat)

## Set predictionMatrix
pred[ ,c("fictiefrecipientnummer", "ID", creat.cols, "Follow_up_time","dDGF")] <- 0

# change variables that should not be imputed
meth[c(creat.cols, "dDGF",  "Follow_up_time", "status2", "Year_Of_Tx")] <- ""

# MICE
imp.qp <- mice(NOTR_creat, meth = meth, predictorMatrix = pred, seed =1, maxit= 20, m=20)

# nothing wrong?
imp.qp$loggedEvents
summary(imp.qp)

# now for the actual creatinine prediction: loop mixed effects model through the imputed sets
# Initialize a list to store the predicted creatinine columns
# Before using this methos
predicted_creatinine_list <- list()
imp_df <- list()
imp_df.long <-list()
creat.lme <-list()

## here we loop the mized model through the 20 imputation of the dataset to create 20 different predicted values for creatinine at all time points.
for (i in 1:imp.qp$m) {
  imp_df[[i]] <- mice::complete(imp.qp, action = i)
  
  imp_df.long[[i]] <- gather(imp_df[[i]], day, creatinine, CreatinineD1:CreatinineD7, factor_key = TRUE)
  imp_df.long[[i]]$day.num <- as.numeric(factor(imp_df.long[[i]]$day, 
                                                levels = c('CreatinineD1', 'CreatinineD2', 'CreatinineD3', 'CreatinineD4',
                                                           'CreatinineD5', 'CreatinineD6', 'CreatinineD7')))
  
  creat.lme[[i]] <- lmer(
    log(creatinine) ~ ns(day.num, 3) + Recipientsex + Recipientage + Recipient_weight + Recipient_height +
      Recipient_Diabetes + Warm_ischaemic_period_2 + Cold_ischaemic_period +
      Donorage + Smoking + Donor_type_deceased +  status2 + DonorDeathGroup.dbd +
      Donorheight + Donorsex + Donorweight + Lastcreatinineµmoll +
      (ns(day.num, 3) | ID), REML = F, data = imp_df.long[[i]], na.action = na.omit
  )
  
  # Append the predicted values to the list
  predicted_creatinine_list[[i]] <- predict(creat.lme[[i]], allow.new.levels = TRUE, newdata = imp_df.long[[i]])
  
}

#### Create long df to add columns 
NOTR_creat.long <- gather(NOTR_creat, day, creatinine, CreatinineD1:CreatinineD7, factor_key=TRUE)
NOTR_creat.long$day.num <- as.numeric(factor(NOTR_creat.long$day, 
                                             levels = c('CreatinineD1','CreatinineD2','CreatinineD3', 'CreatinineD4',
                                                        'CreatinineD5', 'CreatinineD6',
                                                        'CreatinineD7')))

# Add the predicted values as new columns to creat_df.long 
for (i in 1:imp.qp$m) {
  col_name <- paste("creatinine_predicted",i)
  NOTR_creat.long[[col_name]] <- predicted_creatinine_list[[i]]
}


# Calculate the mean of the predictions from multiple imputations

NOTR_creat.long$creatinine_combined <-
  exp(rowMeans(NOTR_creat.long[, grep("creatinine_predicted", colnames(NOTR_creat.long))], na.rm = TRUE))

table(is.na(NOTR_creat.long$creatinine_combined), is.na(NOTR_creat.long$creatinine))



### Check the fit of the model
NOTR_creat.long <- NOTR_creat.long %>% arrange(ID) 

# For the actually imputed cases so if more than 3 values of creatinine are missing how well does the predicted values fit the slope
prbn <- NOTR_creat.long %>%
  group_by(ID) %>%
  filter(sum(is.na(creatinine))>3)

prbn2 <- prbn[1:154,]

# see how it fits the observed data
prbn3 <- NOTR_creat.long[1000:1400,]


entryplot <- prbn3

ggplot(
  entryplot,
  aes(
    x = day.num,
    y = creatinine_combined,
    color = factor(ID),
    group = factor(ID)
  )
) +
  geom_point() +
  geom_line() +
  geom_point(
    data = entryplot,
    aes(x = day.num, y = creatinine),
    shape = 2
  ) +
  facet_wrap( ~ factor(ID), nrow = 4)


### create the imputed value (observed if observed, predicted if missing)
NOTR_creat.long$creatinine_imputed <-
  ifelse(
    is.na(NOTR_creat.long$creatinine),
    NOTR_creat.long$creatinine_combined,
    NOTR_creat.long$creatinine
  ) #final
summary(NOTR_creat.long$creatinine_imputed)

# Sort the data frame by ID and day
NOTR_creat.long <- NOTR_creat.long %>%
  arrange(ID, day.num)


### Clean the data remove all unnecessary creatinine columns (the observed only with missing, the predicted only, and the exp(predicted))
drop <- c( "creatinine_slope", "creatinine", "creatinine_combined",
          names(NOTR_creat.long[grep("creatinine_predicted", colnames(NOTR_creat.long))]),"day.num", "eGFR")

# drop the unnecessary columns 
NOTR_creat.long.drop = NOTR_creat.long[,!(names(NOTR_creat.long) %in% drop)]

NOTR_creatdata.short <- spread(NOTR_creat.long.drop, "day", "creatinine_imputed")

## create df for the imputed values (so these are the observed values plus the missings filled in by the predicted values) 
# we will merge this one with the observed only (so with missings) df
# to make sure that we can still identify the imputed values we should add "imputed" to our column names
NOTR_creatdata.short.imp <- NOTR_creatdata.short %>% rename_with(~paste0(., "_imputed"), CreatinineD1:CreatinineD7)

### Clean the data remove all unnecessary creatinine columns we are only interested in the observed creatinine values
keep <- c("ID", "creatinine", "day")
NOTR_creat.long.drop2 = NOTR_creat.long[,(names(NOTR_creat.long) %in% keep)]

## 
NOTR_creatdata.short.observed <- spread(NOTR_creat.long.drop2, "day", "creatinine")

NOTR_creatdata.short.observed.imp  <- merge(NOTR_creatdata.short.observed, NOTR_creatdata.short, by.x = "ID", by.y = "ID")

NOTR_creatdata.short.observed.imp <- NOTR_creatdata.short.observed.imp %>% group_by (ID)


NOTR_creatdata.short.dgf <-subset(NOTR_creat_1, dDGF == 1)

#### Merge  imputed creat with the NOTR_Data set

## set to characters to match merging set
NOTR_creatdata.short.observed.imp$ID <- as.character(NOTR_creatdata.short.observed.imp$ID)
NOTR_creatdata.short.observed.imp$fictiefrecipientnummer <- as.character(NOTR_creatdata.short.observed.imp$fictiefrecipientnummer)

## bind the imputed set with the DGF cases (we excluded them before imputation)
df <- bind_rows(NOTR_creatdata.short.observed.imp, NOTR_creatdata.short.dgf)

## merge the selected set with the whole dataset (so we effectively add cases that had more than 5 missings in creatinine and no DGF)
df <- merge(select(df,ID, dDGF,  Cold_ischaemic_period , CreatinineD1,CreatinineD2,CreatinineD3, CreatinineD4,
                   CreatinineD5, CreatinineD6,
                   CreatinineD7, CreatinineD1_imputed,CreatinineD2_imputed,CreatinineD3_imputed, CreatinineD4_imputed,
                   CreatinineD5_imputed, CreatinineD6_imputed,
                   CreatinineD7_imputed), NOTR_DGF.df[,!names(NOTR_DGF.df) %in% c("dDGF", "Cold_ischaemic_period",'CreatinineD1','CreatinineD2','CreatinineD3', 'CreatinineD4',
                                                                'CreatinineD5', 'CreatinineD6',
                                                                'CreatinineD7')], by.x = "ID", by.y = "ID")
## check to see if the numbers are still ok
table(df$dDGF)


## create the early graft function variable (DGF, SGF, IGF)
df <- df %>%
  mutate(impaired_start = as.factor(ifelse(
    CreatinineD3_imputed > 0.9* CreatinineD2_imputed
    & CreatinineD2_imputed > 0.9 * CreatinineD1_imputed 
    & dDGF == 0, 1, 0)))
  mutate( early_graft_function = as.factor(case_when(impaired_start == 0 & dDGF == 0 ~  "IGF", 
                                                     impaired_start == 1  ~ "SGF",
                                                     dDGF == 1  ~ "DGF"
  )))%>%
  mutate( early_graft_function = relevel(as.factor(early_graft_function), ref = "IGF")) 


##### Rename to use later
NOTR_DGF_creatinine.df<-df

