## ----load packages, message=FALSE---------------------------------------------
## load required libraries 
#install.packages("corrplot")
library(corrplot)
#install.packages("ggplot2")
library(ggplot2)
#install.packages("tidyverse")
library(tidyverse)
#install.packages("randomForest")
library(randomForest) 

## ----load data----------------------------------------------------------------
## read in data and only consider complete data 
## this drops 327 individuals
## some tree methods handle missing data but we will not deal with that here
nhanes <- na.omit(read.csv(paste0(here::here(),"/Data/studypop.csv")))

## ----format outcome-----------------------------------------------------------
## our y variable - ln transformed and scaled mean telomere length
## some tree methods make assumptions about the error variance and some do not
## generally safer to transform the dependent variable to reduce heteroskedasticity
lnLTL_z <- scale(log(nhanes$TELOMEAN)) 

## ----format exposures---------------------------------------------------------
## exposures matrix
mixture <- with(nhanes, cbind(LBX074LA, LBX099LA, LBX118LA, LBX138LA, LBX153LA, LBX170LA, LBX180LA, LBX187LA, LBX194LA, LBXHXCLA, LBXPCBLA,
                              LBXD03LA, LBXD05LA, LBXD07LA,
                              LBXF03LA, LBXF04LA, LBXF05LA, LBXF08LA)) 
# Most tree models are invariant the transformations of the predictors (they don't make a difference)
# here we transform it to make the plots consistent with other methods so they can be compared
mixture   <- apply(mixture, 2, log)
mixture <- scale(mixture)
colnames(mixture) <- c(paste0("PCB",c(74, 99, 118, 138, 153, 170, 180, 187, 194, 169, 126)), 
                           paste0("Dioxin",1:3), paste0("Furan",1:4)) 
exposure_names <- colnames(mixture)

## ----format covariates--------------------------------------------------------
## our X matrix
covariates <- with(nhanes, cbind(age_cent, male, bmi_cat3, edu_cat, race_cat,
                                 LBXWBCSI, LBXLYPCT, LBXMOPCT, 
                                 LBXNEPCT, LBXEOPCT, LBXBAPCT, ln_lbxcot)) 

## ----regress out covariates---------------------------------------------------
# regress covariates out
lnLTL_z_residuals <- lm(lnLTL_z~covariates)$residuals

## ----fit random forest--------------------------------------------------------
# fit the random forest model
set.seed(1000)
fit_rf <- randomForest(y=lnLTL_z_residuals,
                       x=mixture,
                       ntree=1000,
                       mtry=6,   # number of variables used in each tree
                       importance = TRUE)  # assess mixture component importance

## ----predict with random forests----------------------------------------------
pred_rf <- predict(fit_rf)

# view predicted vs observed
plot(pred_rf~lnLTL_z_residuals, 
     main="Predicted vs observed with random forest")

## ----variable importance with random forests quick plot-----------------------
varImpPlot(fit_rf, main = "Random forest important variables", type=1)

## ----variable importance with random forests----------------------------------
# extract variable importance from the fit random forest model
rf_var_imortance <- importance(fit_rf)
rf_var_imortance

## ----random forest partial dependence plot for Furan 1------------------------
partialPlot(x = fit_rf, 
            pred.data = mixture, 
            x.var = "Furan1",
            main = "Funan 1 partial dependence with random forest", 
            xlab = "Exposure (z-score of log-transformed exposure)",
            ylab = "Estimate (mean response)") 

