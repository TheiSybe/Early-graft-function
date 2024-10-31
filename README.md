# Early-graft-function

## README R script for the analyses published in Differential impact of early graft function in DBD and DCD ##
These scripts belong to doi: (https://doi.org/10.1016/j.ajt.2024.09.030)

*There are three files available in this project

# Datacleaning file
In this file we sort out some difficulties when working with a registry data set.
The rationale for the values that we have created (not just renamed) can be found in the supplementary datafile that was published with the manuscript.

EDIT 08-10-24: unfortunately, the variable names in this script are no longer in use by the NOTR. 

# Creatinine imputation file
In some cases, we couldn't determine if the recipient had full function or slow graft function because one or two creatinine values were missing in the first three days after transplantatation.
We decided that we would need at least two measurements of plasma creatinine within the first 7 days after transplantatation to determine a slope over those 7 days.
If we know the slope, we can predict the values we need.

Because we also have some missing values in the covariates that we want to use to give some more information to the model, we have to solve that first.
First we used the MICE package to impute the covariates we would like to use in the model. This generates 20 imputed data sets. 
We then looped a mixed effects model, which we adapted to the non-linear creatinine data with natural cubic splines (splines package), through all of those sets.

From the examples  and the fit of  the model that we provided in the supplements we can conclude that the predictions seem to be very acurate. 

This file shows the model we have used to approximate the course of plasma creatinine for every individual recipient in the first week after transplantation.


# Analysis and plotting file
Cumulative incidence competing risk analyses and cox regression analysis. Mixed model to analyze and plot the eGFR data.



For questions regarding this script \feel free to contact me:
t.steenvoorden@amsterdamumc.nl 
