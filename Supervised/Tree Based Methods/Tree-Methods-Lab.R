## ----load packages, message=FALSE---------------------------------------------
## load required libraries 
#install.packages("corrplot")
library(corrplot)
#install.packages("ggplot2")
library(ggplot2)
#install.packages("tidyverse")
library(tidyverse)
#install.packages("BART")
library(BART)
#install.packages("randomForest")
library(randomForest) 

## ----load bartmix, message=FALSE----------------------------------------------
# install.packages("remotes")
# remotes::install_github("AnderWilson/bartmix")
library(bartmix)

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

## ----fit BART model-----------------------------------------------------------
# fit the BART model
set.seed(1000)
fit_bart <- gbart(x.train=cbind(mixture,covariates), 
                  y.train=lnLTL_z,
                  nskip=2000,    # MCMC iterations that are discarded as burning
                  ndpost=2000)   # MCMC iterations after burning that are retained for 

load(paste0(here::here(), "/Supervised/Tree Based Methods/bart_fits_for_lab.rda"))

## ----fitted values------------------------------------------------------------
plot(fit_bart$yhat.train.mean~lnLTL_z, 
     main="fitted vs observed with BART")

## ----BART variable importance-------------------------------------------------
# variable importance for the mixture components with BART
varcount <- data.frame(exposure=exposure_names,
                       mean=fit_bart$varcount.mean[exposure_names],
                       lower=apply(fit_bart$varcount[,exposure_names],2,quantile,0.025),
                       upper=apply(fit_bart$varcount[,exposure_names],2,quantile,0.975))

# visualize the variable importance with BART
ggplot(varcount, aes(x=exposure, y=mean, ymin=lower, ymax=upper)) +
  geom_errorbar(width=0) +
  geom_point() +
  theme_minimal() + 
  ggtitle("Variable Inclusion Count with BART") + 
  xlab("Exposure") + 
  ylab("Count") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

## ----BART partial dependence for Furan 1, eval=FALSE--------------------------
#  # partial dependence of BART
#  funan1_pd <- partialdependence1(fit_bart,
#                                  data=cbind(mixture,covariates),
#                                  exposures = "Furan1",
#                                  L=50)

## ----BART partial dependence for Furan 1 visualization------------------------
ggplot(funan1_pd, aes(x=x, y=mean, ymin=lower, ymax=upper)) + 
  geom_ribbon(fill="grey70") + 
  geom_line() + 
  theme_minimal() + 
  ggtitle("Funan 1 partial dependence with BART") + 
  xlab("Exposure (z-score of log-transformed exposure)") + 
  ylab("Estimate (mean response)") 

## ----BART partial dependence for all components, echo=FALSE, eval=FALSE-------
#  all_pd <- partialdependence1(fit_bart,
#                              data=cbind(mixture,covariates),
#                              exposures = exposure_names,
#                              L=50)

## ----BART partial dependence for all components visualization-----------------
plt <- all_pd %>%
mutate(exposure = fct_recode(exposure, "PCB 74" = "PCB74",
                                "PCB 99" = "PCB99",
                                "PCB 118" = "PCB118",
                                "PCB 138" = "PCB138",
                                "PCB 153" = "PCB153",
                                "PCB 170" = "PCB170",
                                "PCB 180" = "PCB180",
                                "PCB 187" = "PCB187",
                                "PCB 194" = "PCB194",
                                "1,2,3,6,7,8-hxcdd" = "Dioxin1",
                                "1,2,3,4,6,7,8-hpcdd" = "Dioxin2",
                               "1,2,3,4,6,7,8,9-ocdd" =  "Dioxin3",
                               "2,3,4,7,8-pncdf" =  "Furan1",
                               "1,2,3,4,7,8-hxcdf" =  "Furan2",
                               "1,2,3,6,7,8-hxcdf" =  "Furan3",
                               "1,2,3,4,6,7,8-hxcdf" =  "Furan4",
                               "PCB 169" =  "PCB169",
                                "PCB 126" = "PCB126")) %>% 
  ggplot(aes(x=x, y=mean, ymin=lower, ymax=upper)) +
  facet_wrap(~exposure) +
  geom_ribbon(fill="grey70") +
  geom_line() +
  theme_minimal() +
  ggtitle("Partial dependence with BART") +
  xlab("Exposure (z-score of log-transformed exposure)") +
  ylab("Estimate (mean response)")
plt

## ----two way partial dependence with BART, eval=FALSE-------------------------
#  pd_2way <- partialdependence2(fit_bart,
#                                  data=cbind(mixture,covariates),
#                                  var = "PCB169",
#                                var2 = "Furan1",
#                                qtls = c(0.1,0.25,0.5,0.75,0.9),
#                                  L=20)

## ----two way partial dependence with BART visualize---------------------------
pd_2way$qtl <- as.factor(pd_2way$qtl)
ggplot(pd_2way, aes(x=x, y=mean, 
                    color=qtl,
                    linetype=qtl)) + 
  geom_line() + 
  theme_minimal() + 
  ggtitle("PCB169 partial dependence by quantile of Furan 1 with BART") + 
  xlab("Exposure (z-score of log-transformed exposure)") + 
  ylab("Estimate (mean response)") 

## ----BART total mixture effect, eval=FALSE------------------------------------
#  totalmix <- totalmixtureeffect(fit_bart,
#                                 data=cbind(mixture,covariates),
#                                 exposures = exposure_names,
#                                 qtls = seq(0.2,0.8,0.05))

## ----BART total mixture effect visualization----------------------------------
ggplot(totalmix, aes(x=quantile, y=mean, ymin=lower, ymax=upper)) + 
  geom_errorbar(width=0) + geom_point() + theme_minimal() +
  ggtitle("Total mixture effect with BART") + 
  xlab("Quantile") + 
  ylab("Estimate (mean response)")

## ----subset the data to male and female---------------------------------------
combind_data <- cbind(mixture,covariates)
male_data_subset <- combind_data[which(combind_data[,"male"]==1),]
female_data_subset <- combind_data[which(combind_data[,"male"]==0),]

# dimension and distribution of male in combined data
dim(combind_data)
table(combind_data[,"male"])
# dimension and distribution of male in male only data
dim(male_data_subset)
table(male_data_subset[,"male"])
# dimension and distribution of male in female only data
dim(female_data_subset)
table(female_data_subset[,"male"])

## ----estimate partial dependence for each subset, eval=FALSE------------------
#  # estimate among males with male subset of data
#  funan1_pd_male_subset <- partialdependence1(fit_bart,
#                                   data=male_data_subset,  # subset of data provided
#                                   exposures = "Furan1",
#                                   L=50)
#  
#  # estimate among females with female subset of data
#  funan1_pd_female_subset <- partialdependence1(fit_bart,
#                                       data=female_data_subset,   # subset of data provided
#                                       exposures = "Furan1",
#                                       L=50)

## ----graph subset specific partial dependence functions-----------------------
funan1_pd_male_subset$subgroup <- "Male"
funan1_pd_female_subset$subgroup <- "Female"
funan1_pd_mf_subset <- rbind(funan1_pd_male_subset,funan1_pd_female_subset)
ggplot(funan1_pd_mf_subset, aes(x=x, y=mean, ymin=lower, ymax=upper)) + 
  geom_ribbon(fill="grey70") + 
  geom_line() + 
  theme_minimal() + 
  ggtitle("Funan 1 partial dependence with BART") + 
  xlab("Exposure (z-score of log-transformed exposure)") + 
  ylab("Estimate (mean response)")  +
  facet_grid(.~subgroup)

## ----create datasets where all individuals are assigned one level of the modifier----
combind_data_allmale <- combind_data
combind_data_allmale[,"male"] <- 1

combind_data_allfemale <- combind_data
combind_data_allfemale[,"male"] <- 0



# dimension and distribution of male in combined data
dim(combind_data)
table(combind_data[,"male"])
# dimension and distribution of male in male only data
dim(combind_data_allmale)
table(combind_data_allmale[,"male"])
# dimension and distribution of male in female only data
dim(combind_data_allfemale)
table(combind_data_allfemale[,"male"])

## ----estimate partial dependence functions on each dataset, eval=FALSE--------
#  funan1_pd_male_all <- partialdependence1(fit_bart,
#                                       data=combind_data_allmale,
#                                       exposures = "Furan1",
#                                       L=50)
#  
#  funan1_pd_female_all <- partialdependence1(fit_bart,
#                                         data=combind_data_allfemale,
#                                         exposures = "Furan1",
#                                         L=50)

## ----graph subset specific partial dependence functions using alternative approach----
funan1_pd_male_all$subgroup <- "Male"
funan1_pd_female_all$subgroup <- "Female"
funan1_pd_mf_all <- rbind(funan1_pd_male_all,funan1_pd_female_all)
ggplot(funan1_pd_mf_all, aes(x=x, y=mean, ymin=lower, ymax=upper)) + 
  geom_ribbon(fill="grey70") + 
  geom_line() + 
  theme_minimal() + 
  ggtitle("Funan 1 partial dependence with BART") + 
  xlab("Exposure (z-score of log-transformed exposure)") + 
  ylab("Estimate (mean response)")  +
  facet_grid(.~subgroup)

