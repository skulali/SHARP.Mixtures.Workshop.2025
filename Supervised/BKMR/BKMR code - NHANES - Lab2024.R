################################################
###  BMKR code; Columbia mixtures workshop   ###
###  developed by Katrina Devick             ###
###  8/10/18                                 ###
###  edited for speed for lab session        ###
###  by Brent Coull                          ###
###  8/22/18                                 ###
###  Last updated by BC: 7/27/2022           ###
################################################

## load required libraries 
#install.packages("bkmr")
library(bkmr)
#install.packages("corrplot")
library(corrplot)
#install.packages("zoo")
library(zoo)
#install.packages("ggplot2")
library(ggplot2)
#install.packages("tidyverse")
library(tidyverse)
#install.packages("Matrix")
library(Matrix)
#install.packages("fields") # need to select "no" when asked question during install
library(fields)
#install.packages("bindrcpp")
library(bindrcpp)


################################################
###         Data Manipulation                ###
################################################

## read in data and only consider complete data 
## this drops 327 individuals, but BKMR does not handle missing data
nhanes <- na.omit(read.csv("Data/studypop.csv"))

## center/scale continous covariates and create indicators for categorical covariates
nhanes$age_z         <- scale(nhanes$age_cent)         ## center and scale age
nhanes$agez_sq       <- nhanes$age_z^2                 ## square this age variable
nhanes$bmicat2       <- as.numeric(nhanes$bmi_cat3==2) ## 25 <= BMI < 30
nhanes$bmicat3       <- as.numeric(nhanes$bmi_cat3==3) ## BMI >= 30 (BMI < 25 is the reference)
nhanes$educat1       <- as.numeric(nhanes$edu_cat==1)  ## no high school diploma
nhanes$educat3       <- as.numeric(nhanes$edu_cat==3)  ## some college or AA degree
nhanes$educat4       <- as.numeric(nhanes$edu_cat==4)  ## college grad or above (reference is high schol grad/GED or equivalent)
nhanes$otherhispanic <- as.numeric(nhanes$race_cat==1) ## other Hispanic or other race - including multi-racial
nhanes$mexamerican   <- as.numeric(nhanes$race_cat==2) ## Mexican American 
nhanes$black         <- as.numeric(nhanes$race_cat==3) ## non-Hispanic Black (non-Hispanic White as reference group)
nhanes$wbcc_z        <- scale(nhanes$LBXWBCSI)
nhanes$lymphocytes_z <- scale(nhanes$LBXLYPCT)
nhanes$monocytes_z   <- scale(nhanes$LBXMOPCT)
nhanes$neutrophils_z <- scale(nhanes$LBXNEPCT)
nhanes$eosinophils_z <- scale(nhanes$LBXEOPCT)
nhanes$basophils_z   <- scale(nhanes$LBXBAPCT)
nhanes$lncotinine_z  <- scale(nhanes$ln_lbxcot)         ## to access smoking status, scaled ln cotinine levels


## our y variable - ln transformed and scaled mean telomere length
lnLTL_z <- scale(log(nhanes$TELOMEAN)) 

## our Z matrix
mixture <- with(nhanes, cbind(LBX074LA, LBX099LA, LBX118LA, LBX138LA, LBX153LA, LBX170LA, LBX180LA, LBX187LA, LBX194LA, LBXHXCLA, LBXPCBLA,
                              LBXD03LA, LBXD05LA, LBXD07LA,
                              LBXF03LA, LBXF04LA, LBXF05LA, LBXF08LA)) 
lnmixture   <- apply(mixture, 2, log)
lnmixture_z <- scale(lnmixture)
colnames(lnmixture_z) <- c(paste0("PCB",c(74, 99, 118, 138, 153, 170, 180, 187, 194, 169, 126)), 
                           paste0("Dioxin",1:3), paste0("Furan",1:4)) 


## our X matrix
covariates <- with(nhanes, cbind(age_z, agez_sq, male, bmicat2, bmicat3, educat1, educat3, educat4, 
                                 otherhispanic, mexamerican, black, wbcc_z, lymphocytes_z, monocytes_z, 
                                 neutrophils_z, eosinophils_z, basophils_z, lncotinine_z)) 


### create knots matrix for Gaussian predictive process (to speed up BKMR with large datasets)
set.seed(10) #if you want exact reproducibility, you need to also set a seed for selecting the knots because this also
##happens at random. 

#generates 50 vectors describing this 18 dimensional space. 
knots50     <- fields::cover.design(lnmixture_z, nd = 50)$design  
#save(knots50, file="Supervised/BKMR/saved_model/NHANES_knots50.RData") 


################################################
###         Fit Models                       ###
################################################

load("Supervised/BKMR/saved_model/NHANES_knots50.RData")

##### fit BKMR models WITH Gaussian predictive process using 50 knots

### Group VS fit with all exposures using GPP and 50 knots 
set.seed(1000)

## Note: the commented out statement is what generated the saved model fits. 
#fit_gvs_knots50 <-  kmbayes(y=lnLTL_z, Z=lnmixture_z, X=covariates, iter=100000, verbose=TRUE, varsel=TRUE, 
#groups=c(rep(1,times=2), 2, rep(1,times=6), rep(3,times=2), rep(2,times=7)), knots=knots50)

## For this lab we are only going to generate 100 MCMC samples to get a sense
##   for what the program outputs.  

  
temp <-  kmbayes(y=lnLTL_z, Z=lnmixture_z, X=covariates, iter=100, verbose=TRUE, varsel=TRUE, 
                             groups=c(rep(1,times=2), 2, rep(1,times=6), rep(3,times=2), rep(2,times=7)), knots=knots50)

## The following statement saved the model fits using 100,000 MCMC samples. We won't save
##    our fit based on 100 samples. Rather we will load in the results from 100,000. 
#save(Supervised/BKMR/saved_model/fit_gvs_knots50,file="bkmr_NHANES_gvs_knots50.RData")

load("Supervised/BKMR/saved_model/bkmr_NHANES_gvs_knots50.RData")

summary(fit_gvs_knots50)

## obtain posterior inclusion probabilities (PIPs)
## PIP -> refers to the probability of inclusion for the given variable
ExtractPIPs(fit_gvs_knots50)

##############################################
###        PLOTS                           ###
##############################################

### correlation matrix
cor.Z <- cor(lnmixture_z, use="complete.obs")

#pdf(file="Supervised/BKMR/figures_pdf/cor_nhanes.pdf", width=12, height=12)
corrplot.mixed(cor.Z, upper = "ellipse", lower.col="black")
#dev.off()


###############################################


### change this for each model you fit and then rerun the code from here to the bottom
modeltoplot      <- fit_gvs_knots50   ## name of model object
modeltoplot.name <- "fit_gvs_knots50" ## name of model for saving purposes
plot.name        <- "gvs_knots50"     ## part that changed in plot name 
Z                <- lnmixture_z        ## Z matrix to match what was used in model

### values to keep after burnin/thin
### I throw out the beginning samples (burn-in) and keep every 10th or 50th 
### sample, as follows. This is because we start the algorithm at initial starting 
###values that might not be close to the final posterior distribution of the parameters of interest. 
###It takes a little while for the stochastic process to converge to the true posterior distribution. 

sel<-seq(50001,100000,by=50)

### assess convergence - the iterations with MCMC have to be in the 
###tens of thousands to see convergence of data with actual distribution - 
###with traceplots 
TracePlot(fit = modeltoplot, par = "beta", sel=sel)
TracePlot(fit = modeltoplot, par = "sigsq.eps", sel=sel) #sigma squared is the residual

par(mfrow=c(2,2))
TracePlot(fit = modeltoplot, par = "r", comp = 1, sel=sel) # these won't look quite as 
#random around the center point because bounce back and forth from 0 to estimate
TracePlot(fit = modeltoplot, par = "r", comp = 2, sel=sel)
TracePlot(fit = modeltoplot, par = "r", comp = 3, sel=sel)
TracePlot(fit = modeltoplot, par = "r", comp = 4, sel=sel)
par(mfrow=c(1,1))

# Remove fit because it takes A LOT of RAM in Rcloud
rm(fit_gvs_knots50)
rm(modeltoplot)

#### create dataframes for ggplot (this takes a little while to run)

## We will not run these 6 commands that generate all the data for the summaries of
## the exposure-response function. Instead we will read in the previously generated results 
## just to save a little time. 

#pred.resp.univar <- PredictorResponseUnivar(fit = modeltoplot, sel=sel, method="approx")
#pred.resp.bivar  <- PredictorResponseBivar(fit = modeltoplot,  min.plot.dist = 1, sel=sel, method="approx")
#pred.resp.bivar.levels <- PredictorResponseBivarLevels(pred.resp.df = pred.resp.bivar, Z = Z,
#                                                          both_pairs = TRUE, qs = c(0.25, 0.5, 0.75))
#risks.overall <- OverallRiskSummaries(fit = modeltoplot, qs = seq(0.25, 0.75, by = 0.05), q.fixed = 0.5, 
#method = "approx",sel=sel)
#risks.singvar <- SingVarRiskSummaries(fit = modeltoplot, qs.diff = c(0.25, 0.75),
#                                        q.fixed = c(0.25, 0.50, 0.75), method = "approx")
#risks.int <- SingVarIntSummaries(fit = modeltoplot, qs.diff = c(0.25, 0.75), qs.fixed = c(0.25, 0.75))

## Save the data need to plot the summaries of interest generated 
#save(pred.resp.univar, pred.resp.bivar, pred.resp.bivar.levels, risks.overall, risks.singvar, risks.int, 
#file=paste0("Supervised/BKMR/saved_model/", modeltoplot.name,"_plots.RData"))

# Load in the results, which were computed previously and saved using the command just above this
load(paste0("Supervised/BKMR/saved_model/", modeltoplot.name,"_plots.RData"))

### run and save ggplots for each bkmr model
#pdf(file=paste0("Supervised/BKMR/figures_pdf/univar_",plot.name,".pdf"), width=15, height=15)
pred.resp.univar %>% 
  mutate(variable = fct_recode(variable, "PCB 74" = "PCB74",
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
  ggplot(aes(z, est, ymin = est - 1.96*se, ymax = est + 1.96*se)) + 
  geom_smooth(stat = "identity") + labs(y = "Estimate", x = "Exposure") + 
  facet_wrap(~ variable) + theme_bw() +
  theme(strip.background = element_rect(fill = "white"))
#dev.off()
# didn't specify linearity (didn't force it to be linear), 
#but the model estimated them to be linear, which provides some insight.


#### Bivariable Exposure-Response Functions
#pdf(file=paste0("Supervised/BKMR/figures_pdf/bivar_levels_",plot.name,".pdf"), width=30, height=30)
pred.resp.bivar.levels %>% 
  mutate(variable1 = fct_recode(variable1, "PCB 74" = "PCB74",
                                "PCB 99" = "PCB99",
                                "PCB 118" = "PCB118",
                                "PCB 138" = "PCB138",
                                "PCB 153" = "PCB153",
                                "PCB 170" = "PCB170",
                                "PCB 180" = "PCB180",
                                "PCB 187" = "PCB187",
                                "PCB 194" = "PCB194",
                                "1,2,3,6,7,8- hxcdd" = "Dioxin1",
                                "1,2,3,4,6,7,8- hpcdd" = "Dioxin2",
                                "1,2,3,4,6,7,8,9- ocdd" =  "Dioxin3",
                                "2,3,4,7,8- pncdf" =  "Furan1",
                                "1,2,3,4,7,8- hxcdf" =  "Furan2",
                                "1,2,3,6,7,8- hxcdf" =  "Furan3",
                                "1,2,3,4,6,7,8- hxcdf" =  "Furan4",
                                "PCB 169" =  "PCB169",
                                "PCB 126" = "PCB126"),
         variable2 = fct_recode(variable2, "PCB 74" = "PCB74",
                                "PCB 99" = "PCB99",
                                "PCB 118" = "PCB118",
                                "PCB 138" = "PCB138",
                                "PCB 153" = "PCB153",
                                "PCB 170" = "PCB170",
                                "PCB 180" = "PCB180",
                                "PCB 187" = "PCB187",
                                "PCB 194" = "PCB194",
                                "1,2,3,6,7,8- hxcdd" = "Dioxin1",
                                "1,2,3,4,6,7,8- hpcdd" = "Dioxin2",
                                "1,2,3,4,6,7,8,9- ocdd" =  "Dioxin3",
                                "2,3,4,7,8- pncdf" =  "Furan1",
                                "1,2,3,4,7,8- hxcdf" =  "Furan2",
                                "1,2,3,6,7,8- hxcdf" =  "Furan3",
                                "1,2,3,4,6,7,8- hxcdf" =  "Furan4",
                                "PCB 169" =  "PCB169",
                                "PCB 126" = "PCB126")) %>% 
  ggplot(aes(z1, est)) + 
  geom_smooth(aes(col = quantile), stat = "identity") + 
  facet_grid(variable2 ~ variable1, scales = "free", space = "free",
             labeller = labeller(variable1 = label_wrap_gen(5),
                                 variable2 = label_wrap_gen(5))) +
  ggtitle("h(Exposure 1 | Quantiles of Exposure 2)") +
  xlab("Exposure 1") + theme_bw() + labs(col = "Quantile", y = "Estimate") +
  theme(strip.background = element_rect(fill = "white"),
        panel.spacing = unit(0.05, "lines"),
        legend.position = "bottom")
#dev.off()

#overall mixture effect - estimated difference in the mean response (telomere length) when 
#comparing all exposures being at lowest percentile and highest percentile
#pdf(file=paste0("Supervised/BKMR/figures_pdf/overallrisks_",plot.name,".pdf"), width=10, height=10)
ggplot(risks.overall, aes(quantile, est, ymin = est - 1.96*sd, ymax = est + 1.96*sd)) + 
  geom_hline(yintercept = 00, linetype = "dashed", color = "gray") + 
  geom_pointrange() + theme_bw() +
  labs(x = "Quantile", y = "Estimate")
#dev.off()

#pdf(file=paste0("Supervised/BKMR/figures_pdf/singvar_",plot.name,".pdf"), width=5, height=10)
risks.singvar %>% mutate(variable = fct_recode(variable, "PCB 74" = "PCB74",
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
  ggplot(aes(variable, est, ymin = est - 1.96*sd,  ymax = est + 1.96*sd, col = q.fixed)) + 
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "gray") +  
  geom_pointrange(position = position_dodge(width = 0.75)) +  coord_flip() + 
  theme_bw() +
  labs(x = "", y = "Estimate", col = "Fixed Quantile")
#dev.off()


#### Single Variable Interaction Terms
#pdf(file=paste0("Supervised/BKMR/figures_pdf/interactplot_",plot.name,".pdf"), width=5, height=10)
risks.int %>% mutate(variable = fct_recode(variable, "PCB 74" = "PCB74",
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
  ggplot(aes(variable, est, ymin = est - 1.96*sd, ymax = est + 1.96*sd)) + 
  geom_pointrange(position = position_dodge(width = 0.75)) + 
  geom_hline(yintercept = 0, lty = 2, col = "gray") + coord_flip() + theme_bw() +
  labs(x = "", y = "Estimate")
#dev.off()

## fit simple linear regression models just to double check we are not finding artifact in the data
## univariate Furan1 model 
summary(lm(lnLTL_z ~lnmixture_z[,"Furan1"] + covariates))

cor(lnmixture_z[,c("PCB169","PCB126","Furan1")], use="complete.obs")
## Furan1, PCB126 and PCB169 main effects
summary(lm(lnLTL_z ~lnmixture_z[,c("PCB169","PCB126","Furan1")] + covariates))

## Furan1:PCB126 and Furan1:PCB169 interaction 
summary(lm(lnLTL_z ~lnmixture_z[,"Furan1"] + lnmixture_z[,"PCB169"] + lnmixture_z[,"PCB126"] + 
             lnmixture_z[,"Furan1"]:lnmixture_z[,"PCB169"]+ + lnmixture_z[,"Furan1"]:lnmixture_z[,"PCB126"]+covariates))

