# Workshop on mixtures, Columbia July 2022

#rm(list = ls())
# install the packages needed for the workshop
#install.packages("gWQS")
# import the packages just installed
library(gWQS)

# define the path  
directory_path_out = paste0(here::here(),"/Supervised/WQS/Output/")
# import the dataset
dataset = read.csv(paste0(here::here(),"/Data/studypop.csv"))
dim(dataset)
dataset = dataset[complete.cases(dataset),]
dim(dataset)

# define the chemicals to include in the mixture
mixture=names(dataset)[grep("LA",names(dataset))]
#mixture = c("LBX074LA", "LBX099LA", "LBX118LA", "LBX138LA", "LBX153LA", "LBX170LA", "LBX180LA", "LBX187LA",
#             "LBX194LA", "LBXD03LA", "LBXD05LA", "LBXD07LA", "LBXF03LA", "LBXF04LA", "LBXF05LA", "LBXF08LA",
#             "LBXHXCLA", "LBXPCBLA")
mixture

################################################################
# correlation among the chemicals
#install.packages("corrplot")
library("corrplot")

dataset_ = dataset[complete.cases(dataset),]
forcorr = data.frame(dataset_[,mixture])
corrplot(cor(forcorr), method="circle", sig.level = c(.001, .01, .05), type="upper")
###############################################################

# log-transform the outcome
dataset$log_TELOMEAN = log(dataset$TELOMEAN)
summary(dataset)
hist(dataset$TELOMEAN, xlab="TELOMEAN", col="blue")
hist(dataset$log_TELOMEAN, xlab="log(TELOMEAN)", col='green')

###############################################################
############################# single chemical analyses
library(broom)
library(tidyr)
library(dplyr)
### log2 transformed
datalog = dataset
datalog[, mixture] <- lapply(dataset[, mixture], log, 2)

##  Quantile-transformation
dataquant = dataset
dataquant[, mixture] <- lapply(dataset[, mixture], ntile, 10)

#######################################################
# Linear Model w/ one biomarker 
# create a list of markers of interest
biomarkers <- mixture
# linear model
lm_fit <-
  lapply(biomarkers, function(x){
    temp = datalog[,c(mixture, "log_TELOMEAN", 
                      "LBXWBCSI", "LBXLYPCT", "LBXMOPCT","LBXEOPCT", "LBXBAPCT","LBXNEPCT",
                      "age_cent", "age_sq", "race_cat", "bmi_cat3", "ln_lbxcot", "edu_cat", "male")] 
    temp= temp[complete.cases(temp),]
    glm(
      formula = substitute(log_TELOMEAN ~ i + 
                             LBXWBCSI + LBXLYPCT + LBXMOPCT + LBXEOPCT + LBXBAPCT + LBXNEPCT + age_cent + age_sq + 
                             as.factor(race_cat) + as.factor(bmi_cat3) + ln_lbxcot + as.factor(edu_cat) + male,
                           list(i = as.name(x))),
      data = temp,
      family = gaussian())
  })
warnings()
# creates a table with multiple lists
# number of lists created = number of biomarkers chosen
# extract coefficients with confidence interval
lm_coefs <- lapply(lm_fit, tidy, conf.int = T)
names(lm_coefs) <- biomarkers

# combine all lists into one big table
# gives estimates and CI for biomarkers and variables
lm_coefs <- bind_rows(lm_coefs, .id = "biomarkers")

# subset estimates and CI for _only the biomarkers_
lm_coefs <- filter(lm_coefs, term %in% biomarkers)
warnings()
lm_coefs

write.csv(lm_coefs, file=paste0(directory_path_out, "single_log2.csv"))
#write.csv(lm_coefs, file=paste0(directory_path_out, "single_decile.csv"))
####################################################################################
## multiple linear regression
####################################################################################
mixture
multipleGLM = glm(log_TELOMEAN ~ 
                    LBX074LA + LBX099LA + LBX118LA + LBX138LA + LBX153LA + LBX170LA + LBX180LA + LBX187LA + LBX194LA +
                    LBXD03LA + LBXD05LA + LBXD07LA + LBXF03LA + LBXF04LA + LBXF05LA + LBXF08LA + LBXPCBLA + LBXHXCLA +
                    LBXWBCSI + LBXLYPCT + LBXMOPCT + LBXEOPCT + LBXBAPCT + LBXNEPCT + age_cent + age_sq + 
                    as.factor(race_cat) + as.factor(bmi_cat3) + ln_lbxcot + as.factor(edu_cat) + male,  data = datalog)
summary(multipleGLM)
matrix_coef <- summary(multipleGLM)$coefficients  # Extract coefficients in matrix
matrix_coef   

write.csv(matrix_coef, file=paste0(directory_path_out, "multiplereg.csv"))


####################################################################################

# covariate only model
cov_only = glm(log_TELOMEAN ~ LBXWBCSI + LBXLYPCT + LBXMOPCT + LBXEOPCT + LBXBAPCT + LBXNEPCT + age_cent + age_sq + 
                as.factor(race_cat) + as.factor(bmi_cat3) + ln_lbxcot + as.factor(edu_cat) + male,  data = dataset)
summary(cov_only)

#####################################################################################
#  WQS regression
#####################################################################################

# fit an unadjusted model to look at the association between the mixture and the outcome
# TELOMEAN = Mean Telomere Length
## there's a signal only in the negative direction
#summary(dataset[,mixture])
results1 = gwqs(log_TELOMEAN ~ wqs, mix_name = mixture, data = dataset, q = 10, validation = 0.6, 
                b = 100, b1_pos = FALSE, b1_constr = FALSE, family = "gaussian", seed = 123)
gwqs_barplot(results1)
gwqs_scatterplot(results1)
gwqs_fitted_vs_resid(results1)
gwqs_summary_tab(results1)
sum(results1$bres$b1 < 0) 
## all 100 bootstraps had negative beta estmates in the GLM

summary(results1$fit)
results1$final_weights

# adjusting for covariates:
# blood data: LBXWBCSI LBXLYPCT LBXMOPCT LBXEOPCT LBXBAPCT LBXNEPCT
# demographics: age_cent age_sq race_cat bmi_cat3 ln_lbxcot edu_cat male
dataset$race_cat=as.factor(dataset$race_cat)
dataset$bmi_cat3=as.factor(dataset$bmi_cat3)
dataset$edu_cat=as.factor(dataset$edu_cat)
is.factor(dataset$race_cat)
is.factor(dataset$bmi_cat3)
is.factor(dataset$edu_cat)
is.factor(dataset$male)

# positive direction -  when adjusted the association is in the positive direction: F08 F03 PCB
result2 = gwqs(log_TELOMEAN ~ wqs + 
                 LBXWBCSI + LBXLYPCT + LBXMOPCT + LBXEOPCT + LBXBAPCT + LBXNEPCT + age_cent + age_sq + 
                 race_cat + bmi_cat3 + ln_lbxcot + edu_cat + male, 
               mix_name = mixture, data = dataset, q = 10, 
               validation = 0.6, b = 100, b1_pos = TRUE, b1_constr = FALSE, family = "gaussian", 
               seed = 123)  
## signal="abst","t2"(default),"expt"
## can require stricter convergence:  control = list(trace = FALSE, maxit = 2000, reltol = 1e-10) 
## different optimization algorithms, optim.method="BFGS"
## plan_strategy = "multisession"
gwqs_barplot(result2)
gwqs_summary_tab(result2)
gwqs_scatterplot(result2)
sum(result2$bres$b1 > 0)   ## 95 of 100 are positive beta estimates
sum(result2$bres$b1 < 0)

# negative direction - although constrained to be negative, only 23 of 100 bootstraps resulted in a negative coeff; 
#                      final estimate is positive and sig
result3 = gwqs(log_TELOMEAN ~ wqs+LBXWBCSI + LBXLYPCT + LBXMOPCT + LBXEOPCT + LBXBAPCT + LBXNEPCT + age_cent + age_sq + 
                 race_cat + bmi_cat3 + ln_lbxcot + edu_cat + male, 
               mix_name = mixture, 
               data = dataset, q = 10, 
               validation = 0.6,  b = 100, b1_pos = FALSE, b1_constr = TRUE, family = "gaussian", 
               seed = 123)
gwqs_summary_tab(result3)
gwqs_barplot(result3)
gwqs_scatterplot(result3)
sum(result3$bres$b1 < 0)  #23 of 100 are negative
sum(result3$bres$b1 > 0)

#compare to random subset WQS - F08 F03 PCB D05 D07 pos and sig
result4 = gwqs(log_TELOMEAN ~ wqs +
                 LBXWBCSI + LBXLYPCT + LBXMOPCT + LBXEOPCT + LBXBAPCT + LBXNEPCT + age_cent + age_sq + 
                 race_cat + bmi_cat3 + ln_lbxcot + edu_cat + male, 
               mix_name = mixture, data = dataset, q = 10, 
               validation = 0.6,  b = 500, b1_pos = TRUE, b1_constr = TRUE, family = "gaussian", 
               seed = 123,rs=TRUE, n_var=5) # 5 for set of 18

gwqs_summary_tab(result4) 
gwqs_barplot(result4)
gwqs_scatterplot(result4)
gwqs_weights_tab(result4)
sum(result4$bres$b1 > 0)

# WQS with interaction - b12 NS (p=0.378)
result2int = gwqs(log_TELOMEAN ~ wqs*male +
                 LBXWBCSI + LBXLYPCT + LBXMOPCT + LBXEOPCT + LBXBAPCT + LBXNEPCT + age_cent + age_sq + 
                 race_cat + bmi_cat3 + ln_lbxcot + edu_cat + male, 
               mix_name = mixture, data = dataset, q = 10, 
               validation = 0.6,  b = 100, b1_pos = TRUE, b1_constr = TRUE, family = "gaussian", 
               seed = 123)

gwqs_summary_tab(result2int) 
gwqs_barplot(result2int)
gwqs_scatterplot(result2int)
gwqs_weights_tab(result2int)


# run the wqs model using the stratified variables in the mixtures  
dataset$sex = as.factor(dataset$male) 
is.factor(dataset$male)
is.factor(dataset$sex) 
## females(0): F08 F03 PCB 194 180 138 D05
## males(1): F03 F08 D07 PCB (~D05)
result5 = gwqs(log_TELOMEAN ~ wqs +
                 LBXWBCSI + LBXLYPCT + LBXMOPCT + LBXEOPCT + LBXBAPCT + LBXNEPCT + age_cent + age_sq + 
                 race_cat + bmi_cat3 + ln_lbxcot + edu_cat + male, 
               stratified="sex", mix_name = mixture, data = dataset, q = 10, 
               validation = 0.6,  b = 100, b1_pos = TRUE, b1_constr = TRUE, family = "gaussian", 
               seed = 123)

gwqs_summary_tab(result5)  
gwqs_barplot(result5) 
gwqs_scatterplot(result5) 
gwqs_fitted_vs_resid(result5) 
gwqs_weights_tab(result5) 


# run the wqs model using the stratified variables in the mixtures with interaction
## the interaction is NS (p=0.851)
result6 = gwqs(log_TELOMEAN ~ wqs*male +
                 LBXWBCSI + LBXLYPCT + LBXMOPCT + LBXEOPCT + LBXBAPCT + LBXNEPCT + age_cent + age_sq + 
                 race_cat + bmi_cat3 + ln_lbxcot + edu_cat + male, 
               stratified="sex", mix_name = mixture, data = dataset, q = 10, 
               validation = 0.6,  b = 100, b1_pos = TRUE, b1_constr = TRUE, family = "gaussian", 
               seed = 123)

gwqs_summary_tab(result6) 
gwqs_scatterplot(result6) 
gwqs_fitted_vs_resid(result6) 
gwqs_barplot(result6) 


#####################################################################
## divergent plots
## load libraries
library(gapminder)
library(tidyverse)
library(stringr)

forplot = result6$final_weights  
#forplot = result5rh$final_weights  ## specify the data to be used in the plot; i.e., WQS str or WQS str int model

forplot$strata = str_sub(forplot$mix_name,-1)
forplot$chem = str_sub(forplot$mix_name,1,nchar(as.character(forplot$mix_name))-2)
forplot$strataname = ifelse(forplot$strata == "1", "Male", "Female")

resultstr <- forplot %>%
  mutate(mean_weight = ifelse(strata == "1",
                              mean_weight,
                              -1*mean_weight)) 
resultstr
dim(resultstr)
## calculate breaks values
breaks_values <- pretty(resultstr$mean_weight)

ggplot(data=resultstr,
       aes(reorder(x = chem, mean_weight), y = mean_weight, fill = strataname))+
  geom_bar(stat = "identity")+
  coord_flip()+
  scale_y_continuous(breaks = breaks_values,
                     labels = abs(breaks_values))+
  theme_minimal()+
  geom_hline(yintercept = 1/nrow(resultstr), linetype='dotted', col='black')+  
  geom_hline(yintercept = -1/nrow(resultstr), linetype='dotted', col='black')+
  geom_hline(yintercept = 0, col='black')+  
  labs(x="Component", y="Mean Weight", fill="Strata")

wt_perc_males = round(100*sum(resultstr[resultstr$strata == 1,]$mean_weight),1)
wt_perc_males

################################ takes FOREVER ##########################
#repeated holdout WQS
# positive direction
result5rh = gwqsrh(log_TELOMEAN ~ wqs +
                 LBXWBCSI + LBXLYPCT + LBXMOPCT + LBXEOPCT + LBXBAPCT + LBXNEPCT + age_cent + age_sq + 
                 race_cat + bmi_cat3 + ln_lbxcot + edu_cat + sex, 
                 stratified="sex", mix_name = mixture, data = dataset, q = 10, 
               validation = 0.6, b = 100, b1_pos = TRUE, b1_constr = TRUE, family = "gaussian", 
               seed = 123, rh=30,
               plan_strategy = "multisession")
gwqs_summary_tab(result5rh)
gwqsrh_boxplot(result5rh)
gwqs_weights_tab(result5rh)
summary(result5rh$fit)
result5rh$final_weights
result5rh$final_weights$mean_weight = result5rh$final_weights$Estimate


write.csv(result5rh$wmat, paste0(directory_path_out, "labstWQSwtsrh30.csv"))


result6rh = gwqsrh(log_TELOMEAN ~ wqs*sex +
                     LBXWBCSI + LBXLYPCT + LBXMOPCT + LBXEOPCT + LBXBAPCT + LBXNEPCT + age_cent + age_sq + 
                     race_cat + bmi_cat3 + ln_lbxcot + edu_cat + sex, 
                   stratified="sex", mix_name = mixture, data = dataset, q = 10, 
                   validation = 0.6, b = 100, b1_pos = TRUE, b1_constr = TRUE, family = "gaussian", 
                   seed = 123, rh=30,
                   plan_strategy = "multisession")
gwqs_summary_tab(result6rh)
gwqsrh_boxplot(result6rh)
gwqs_weights_tab(result6rh)
summary(result6rh$fit)
result6rh$final_weights
result6rh$final_weights$mean_weight = result6rh$final_weights$Estimate


write.csv(result6rh$wmat, paste0(directory_path_out, "labstintWQSwtsrh30.csv"))


#### boxplots ################################################################
## for WQS st rh model 
boxplot(result5rh$coefmat[,"wqs"],
       xlab="WQS",  ylab = "Beta Estimates, rh=30", 
       col = c("orange"),  border = "brown", horizontal = FALSE,  notch = TRUE)


## for WQS st int rh model
WQS = result6rh$coefmat[,"wqs"]
WQS_SEX = result6rh$coefmat[,"wqs:sex1"]

boxplot(WQS, WQS_SEX,
        at = c(1,2),
        names = c("WQS", "WQS_SEX"),
        ylab = "Beta Estimates, rh=30",
        col = c("orange","red"),
        border = "brown",
        horizontal = FALSE,
        notch = TRUE
)

##############################################################
# plots for repeated holdout weights with multiple chemical classes
# Author: Marlene Stratmann
# Date: 19.04.2022
# Title: GGPlot for RHWQS 

library(ggplot2)
library(tidyverse)
library(dplyr)
library(magrittr)
library(reshape2)

weightswide = result6rh$wmat
weights = melt(weightswide)
names(weights) = c("Rep", "name", "mean_weight")
head(weights)

weights <- weights  %>% mutate( group = str_sub(name, start=10), mix_name=str_sub(name, start=4, end=6))
head(weights)

# add the categories to the chemicals
weights <- weights %>%
  mutate(
    type_chemical = ifelse(
        mix_name == "F03" | 
        mix_name == "F04" | 
        mix_name == "F05" |
        mix_name == "F08",
      "Furans",
      ifelse(
          mix_name == "D03" |
          mix_name == "D05" | 
          mix_name == "D07" ,
        "Dioxins",
        ifelse(
            mix_name == "PCB" |
            mix_name == "HXC" | 
            mix_name == "138" | 
            mix_name == "153" | 
            mix_name == "118" | 
            mix_name == "194" | 
            mix_name == "180" | 
            mix_name == "170" | 
            mix_name == "187" | 
            mix_name == "074" | 
            mix_name == "099" ,
          "PCBs", 0
        ))))

# Factor the chemicals in order to make the plot
weights<- weights %>%
  mutate(type_chemical = as.factor(type_chemical))
head(weights)

# Create colour palettes inspired by Tanner et al paper
awesome_colours <- c(
  `Dioxins`     = "#991F85",
  `Furans`      = "#0000FF",
  `PCBs`        = "#138C89")

################################################
#reorder the levels so it looks better in the plot. 
#GGPlot does not accept if you change the order beforehand
plot <- weights %>% filter(group == "0" ) %>% 
  mutate(mix_name = fct_relevel(mix_name,
                                "D07", "D05", "D03")) %>%
  mutate(mix_name = fct_relevel(mix_name,
                                "F08", "F05", "F04", "F03"    )) %>%
  mutate(mix_name = fct_relevel(mix_name,
                                "PCB", "HXC", "138", "153", "118", "194", "180", "170", "187", "074", "099")) %>%
  mutate(type_chemical = fct_relevel(type_chemical,
                                     "PFAS", "Pesticides", "PCBs")) %>%
  ggplot(aes(x = mix_name, #Now I can start building the plot
             y = mean_weight*100, 
             fill = type_chemical)) + 
  geom_boxplot(alpha = 0.5, 
               outlier.shape = NA, 
               show.legend = FALSE) + #makes the boxes more transparent and removes the points for the outliers
  scale_x_discrete(guide = guide_axis(angle = -45), limits = rev) + # I rotate the legend of the x-axis to 45 degrees
  theme_classic() + # Just a simple and minimal theme for the legends and the background of the plot
  geom_hline(yintercept = 3, 
             color = "black", 
             size = 0.45) + # add the threshold line here
  labs(y = "Weight (%)", 
       x = "Chemical",
       #      title = "Identification of chemicals of concern with repeated holdout validation",
       caption = "The black line indicates the threshold guideline for chemicals of concern",
       fill = "Type of chemical") +
  theme(legend.position = "bottom", 
        legend.title = element_text(size = 10),
        legend.key.width = unit(6.5, "pt")) + # This puts the legends on the bottom and closer to each other in the legend points
  geom_point(size = 1, 
             alpha = 0.4, 
             show.legend = TRUE, 
             aes(color = type_chemical)) +
  scale_fill_manual(name = "Type of chemical", 
                    labels = c("Dioxins", "Furans",  "PCBs"), 
                    values = awesome_colours) + # apply the custom colour palette to the boxplot
  scale_color_manual(name = "Type of chemical", 
                     labels = c("Dioxins", "Furans", "PCBs"), 
                     values = awesome_colours) + #make the points small and a little transparent and hide the legend
  guides(color = guide_legend(override.aes = list(size = 5, 
                                                  shape = 15, 
                                                  alpha = 0.8)))

plot

# Extend the plot with the diamonds 
meanwts <- weights %>% filter(group == "0") %>%
  group_by(mix_name) %>%
  mutate(meanwt = mean(mean_weight))
meanwts

awesome_plot <- plot +
  geom_point(data = meanwts, 
             size = 2, 
             alpha = 0.8, 
             show.legend = FALSE, 
             shape = 5, #Diamond shape
             aes(x = mix_name, 
                 y = meanwt * 100, 
                 stroke = 1))  # stroke makes the borders thicker

awesome_plot

#################################################################################
## violin plots for betas
head(result6rh$coefmat)
betas = result6rh$coefmat[, c(2, 21)]  # wqs and interaction columns

b1 = betas[,1]
b12 = betas[,2]
b2 = betas[,1] + betas[,2]

round(quantile(b1, probs=c(0.025, 0.05, 0.5, 0.95, 0.975)), digits=3)
round(quantile(b12, probs=c(0.025, 0.05, 0.5, 0.95, 0.975)), digits=3)
round(quantile(b2, probs=c(0.025, 0.05, 0.1, 0.5,0.9, 0.95, 0.975)), digits=3)

#install.packages("vioplot") 
library(vioplot)
vioplot(b1, b2, b12, names=c("WQS_girls", "WQS_boys", "WQS*sex"), col=c("green", "green", "blue"),
        ylab="Beta Estimates")





