#--- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
#--- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# Comparison of 4 precision medicine strategies - main analysis  ---------------
#--- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
#--- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---


#--- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# Paths ------------------------------------------------------------------------
#--- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---


setwd(" ")

path_data              <- " "
plot_path              <- " "
additional_code        <- " "
study_code             <- " "
function_path          <- " " 

path_summary           <- c("path_summary","path_data", "plot_path", "additional_code", "study_code", "function_path")



#--- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# Packages ---------------------------------------------------------------------
#--- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

# install.packages("data.table")                                                # for rbindlist and melt function function
library(data.table)

# install.packages("haven")                                                     # for opeing STATA data files
library(haven)

# install.packages("rms")                                                       # for lm with rcs 
library(rms)

# install.packages("dplyr")                                                     # group analysis and reshaping of data
library(dplyr)

# install.packages("tidyr")                                                     # group analysis and reshaping of data
library(tidyr)

# install.packages("lme4")                                                      # mixed effect models (lmer function)
library(lme4)

# install.packages("Matrix")                                                    # apparently needed to run lmer command
library(Matrix) 

# install.packages("performance")                                               # r2_nakagawa function 
library(performance)

# install.packages("ggplot2")                                                   # plots 
library(ggplot2)

# install.packages("rstatix")                                                   # for 3 group t-test comparisons
library(rstatix)

# install.packages("emmeans")                                                   # emmeans function
library(emmeans)

# install.packages("cowplot")                                                   # to plot in a grid (https://cran.r-project.org/web/packages/cowplot/vignettes/introduction.html)
library("cowplot")

# install.packages("Cairo")                                                     # for saveing nice plots
library(Cairo)

# install.packages("lmerTest")                                                  # to get p values from lmer models
library(lmerTest)

# install.packages("car")                                                       # for contrasting regression coefficients with linearHypothesis()
library(car)
    
# install.packages("effectsize")                                                # for function standardize()
library(effectsize)

# install.packages("performance")                                               # R^2 for mixed model
library(performance)


#--- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# Load in data -----------------------------------------------------------------
#--- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

#
## Main dataset of TriMaster cohort --------------------------------------------
#

load(paste0(path_data, "/data_wcluster.Rdata"))                                 

#
## 5 Drug Model, Dennis et al. (2025) ------------------------------------------
#

load(" ")




#
## Load and merge data for smoking ---------------------------------------------
#

smoke_data_merge_orig <- data.frame(read_dta(paste0(" ")))
smoke_data_merge      <- data.frame(study_id = as.vector(smoke_data_merge_orig$study_id), smoke = smoke_data_merge_orig$v1_Smoker)
# table(is.na(smoke_data_merge$smoke))

data <- merge(data_wcluster, smoke_data_merge, by = "study_id", all.x = TRUE)

rm(data_wcluster)
rm(smoke_data_merge_orig)
rm(smoke_data_merge)


#
## Load in genetic cluster data and merge --------------------------------------
#

load(paste0(path_data, "/genetic_cluser_dataRdata"))

table(genetic_cluser_data$study_id%in%data$study_id)

data_merged <- merge(data, genetic_cluser_data, by = "study_id", all.x = T)
data        <- data_merged

rm(data_merged)
rm(genetic_cluser_data)



#--- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# Data preparation -------------------------------------------------------------
#--- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

# Exclude participants from the 1_SAID cluster in the analysis
data_clusterexcluded <- data[!data$clusterID_name == "1_SAID", ]                # 10 GAD positive individuals excluded

# generate a few variables needed later
data_clusterexcluded$drugline <- data_clusterexcluded$vs_SU + data_clusterexcluded$vs_MFN + 1
# data_clusterexcluded$vs_DIET_only[is.na(data_clusterexcluded$drugline)]

data_clusterexcluded$ncurrtx  <- data_clusterexcluded$drugline

variables_to_include <- c("study_id", "clusterID_name", "validhba1cA", "validhba1cB", "validhba1cC", 
                          "TreatmentOrder", "vs_HbA1c_result", "vs_Age_at_diagnosis", "vs_Gender", 
                          "vs_eGFR", "vs_BMI", "vs_ALT_result", "CHOLmmolLV1", "HDLmmolLV1", "CREATumolLV1",
                          "v1_Blood_Pressure_1sys", "v1_Blood_Pressure_1dia", 
                          "v2v1datediff", "v3v2datediff", "v4v3datediff", "smoke",
                          "vs_Ethnic_Group", "ageatscreening", "drugline", "ncurrtx", "HOMA_B", "HOMA_IR", 
                          "Drug1rank", "Drug1rank2", "Drug2rank", "Drug2rank2", "Drug3rank", "Drug3rank2",
                          "incp1", "incp2", "incp3", 
                          "t2d_betacell_pineg", "t2d_betacell_pipos", "t2d_bodyfat", "t2d_lipodystrophy", 
                          "t2d_liverlipid", "t2d_metabollicsyn", "t2d_obesity", "t2d_residualglycaemic")

# reshape dataset 
data_long_temp            <- melt(setDT(data_clusterexcluded[, variables_to_include]), 
                                  id.vars = c("study_id", "clusterID_name", "TreatmentOrder", "vs_HbA1c_result",
                                              "vs_Age_at_diagnosis", "vs_Gender", 
                                              "vs_eGFR", "vs_BMI", "vs_ALT_result", "CHOLmmolLV1", "HDLmmolLV1", "CREATumolLV1",
                                              "v1_Blood_Pressure_1sys", "v1_Blood_Pressure_1dia", 
                                              "v2v1datediff", "v3v2datediff", "v4v3datediff",  "smoke",
                                              "vs_Ethnic_Group", "ageatscreening", "drugline", "ncurrtx", 
                                              "HOMA_B", "HOMA_IR", 
                                              "Drug1rank", "Drug1rank2", "Drug2rank", "Drug2rank2", "Drug3rank", "Drug3rank2",
                                              "incp1", "incp2", "incp3", 
                                              "t2d_betacell_pineg", "t2d_betacell_pipos", "t2d_bodyfat", 
                                              "t2d_lipodystrophy", "t2d_liverlipid", "t2d_metabollicsyn", "t2d_obesity", "t2d_residualglycaemic"),
                                  variable.name = c("hba1c"))

# removed all participants that have not valid hba1c recorded (will not have a study_id in long data format)
data_long            <- data_long_temp[!is.na(data_long_temp$study_id), ]                  

colnames(data_long)  <- c("study_id", "clusterID_name", "TreatmentOrder", "prehba1c", "Age_at_diagnosis", "sex",
                          "preegfr", "prebmi", "prealt", "pretotalcholesterol", "prehdl", "precreat", "presys", "predia", 
                          "v2v1datediff", "v3v2datediff", "v4v3datediff", "smoke", 
                          "ethnicity", "ageatscreening", "drugline", "ncurrtx", "HOMA_B", "HOMA_IR",
                          "Drug1rank", "Drug1rank2", "Drug2rank", "Drug2rank2", "Drug3rank", "Drug3rank2", 
                          "tol1", "tol2", "tol3", 
                          "t2d_betacell_pineg", "t2d_betacell_pipos", "t2d_bodyfat", 
                          "t2d_lipodystrophy", "t2d_liverlipid", "t2d_metabollicsyn", "t2d_obesity", "t2d_residualglycaemic", 
                          "hba1c_tx", "posthba1cfinal")

# treatment variable
data_long$tx_taken                                 <- rep(NA, times = dim(data_long)[1])
data_long$tx_taken[data_long$hba1c_tx == "validhba1cA"] <- "A"
data_long$tx_taken[data_long$hba1c_tx == "validhba1cB"] <- "B"
data_long$tx_taken[data_long$hba1c_tx == "validhba1cC"] <- "C"

data_long$drugclass                                      <- rep(NA, times = dim(data_long)[1])
data_long$drugclass[data_long$hba1c_tx == "validhba1cA"] <- "TZD"
data_long$drugclass[data_long$hba1c_tx == "validhba1cB"] <- "DPP4"
data_long$drugclass[data_long$hba1c_tx == "validhba1cC"] <- "SGLT2"

# period variable 
data_long$period                            <- rep(NA, times = dim(data_long)[1])
data_long$period[data_long$tx_taken == "A"] <- unlist(gregexpr("A", data_long$TreatmentOrder[data_long$tx_taken == "A"]))
data_long$period[data_long$tx_taken == "B"] <- unlist(gregexpr("B", data_long$TreatmentOrder[data_long$tx_taken == "B"]))
data_long$period[data_long$tx_taken == "C"] <- unlist(gregexpr("C", data_long$TreatmentOrder[data_long$tx_taken == "C"]))

# create drugclass as factor with all 5 treatments
data_long$drugclass                         <- factor(data_long$drugclass, levels = c("DPP4", "GLP1", "SGLT2", "SU", "TZD"))

# create drugline as factor
data_long$drugline                          <- factor(data_long$drugline, levels = c("2", "3", "4", "5+"))

# create ncurrtx as factor
data_long$ncurrtx                           <- data_long$ncurrtx - 1            # as per 5 drug model definition (Dennis et al. 2025, Lancet)
data_long$ncurrtx                           <- factor(data_long$ncurrtx, levels = c("1", "2", "3", "4+"))


# change coding of ethnicity variable to fit with model
table(data_long$ethnicity)
data_long$ethnicity[data_long$ethnicity == "A"] <- "White"
data_long$ethnicity[data_long$ethnicity == "B"] <- "White"
data_long$ethnicity[data_long$ethnicity == "C"] <- "White"
data_long$ethnicity[data_long$ethnicity == "E"] <- "Mixed"
data_long$ethnicity[data_long$ethnicity == "H"] <- "South Asian"
data_long$ethnicity[data_long$ethnicity == "J"] <- "South Asian"
data_long$ethnicity[data_long$ethnicity == "K"] <- "South Asian"
data_long$ethnicity[data_long$ethnicity == "L"] <- "South Asian"
data_long$ethnicity[data_long$ethnicity == "N"] <- "Black"
data_long$ethnicity[data_long$ethnicity == "S"] <- "Other"
data_long$ethnicity[data_long$ethnicity == "Z"] <- "Missing"

data_long$ethnicity                             <- factor(data_long$ethnicity, levels = c("White", "South Asian", "Black", "Mixed", "Other", "Missing"))

# deprivation index variable
data_long$imd5                            <- rep("3", times = dim(data_long)[1])
data_long$imd5                            <- factor(data_long$imd5, levels = c("1 (least)", "2", "3", "4", "5 (most)"))

# age at treatment variable
data_long$agetx                           <- data_long$ageatscreening  

# sex variable needs to be Female/Male, female is baseline
data_long$sex[data_long$sex == "M"]        <- "Male"
data_long$sex[data_long$sex == "F"]        <- "Female"

str(data_long$sex)
data_long$sex                              <- as.factor(data_long$sex)

# T2D duration variable
data_long$t2dmduration                     <- data_long$agetx - data_long$Age_at_diagnosis

summary(data_long$t2dmduration)

# recode smoke variable
data_long$smoke[data_long$smoke == 0]   <- "Non-smoker"
data_long$smoke[data_long$smoke == 1]   <- "Active smoker"
data_long$smoke[is.na(data_long$smoke)] <- "Not recorded"
table(data_long$smoke)
data_long$smoke                         <- factor(data_long$smoke, levels = c("Non-smoker", "Active smoker", "Ex-smoker", "Not recorded"))
str(data_long$smoke)
table(data_long$smoke)

# hba1cmonth variable
data_long$hba1cmonth                    <- rep(NA, times = dim(data_long)[1])

study_id_unique                         <- unique(data_long$study_id)
for(i in 1: length(study_id_unique)){
  
  subset_participant <- subset(data_long, data_long$study_id == study_id_unique[i])
  
  data_long$hba1cmonth[data_long$study_id%in%study_id_unique[i]] <- as.vector(unlist(subset_participant[1, c("v2v1datediff", "v3v2datediff", "v4v3datediff")]))[subset_participant$period]
  
  print(i)
}

data_long$hba1cmonth                     <- data_long$hba1cmonth /30.4


# just checking:
#subset(data_long, data_long$study_id == "TM010002")
#subset(data_long, data_long$study_id == "TM030300")


data_long_cc      <- data_long[!is.na(data_long$posthba1cfinal), ]
table(table(data_long_cc$study_id))
# 50 participants have only one valid hba1c recorded
# 96 participants have two valid hba1c recorded
# 309 participants have 3 3 valid hba1c recorded


# Check for missing values 
contains_any_na = sapply(data_long_cc, function(x) any(is.na(x)))
names(data_long_cc)[contains_any_na]



#
## Set up study cohort ---------------------------------------------------------
#

subset_3treatments <- data_long_cc[data_long_cc$study_id%in%names(which(table(data_long_cc$study_id) == 3)), ]


# Some data preparation for further analysis
subset_3treatments$hba1c_diff <- subset_3treatments$posthba1cfinal - subset_3treatments$prehba1c

subset_3treatments$treatment_class <- rep(NA , times = dim(subset_3treatments)[1])
subset_3treatments$treatment_class[subset_3treatments$drugclass == "SGLT2"] <- rep("SGLT2", times = as.numeric(table(subset_3treatments$drugclass == "SGLT2")["TRUE"]))
subset_3treatments$treatment_class[subset_3treatments$drugclass == "DPP4"]  <- rep("DPP4", times = as.numeric(table(subset_3treatments$drugclass == "DPP4")["TRUE"]))
subset_3treatments$treatment_class[subset_3treatments$drugclass == "TZD"]   <- rep("TZD", times = as.numeric(table(subset_3treatments$drugclass == "TZD")["TRUE"]))
subset_3treatments$treatment_class                                          <- factor(subset_3treatments$treatment_class, levels = c("SGLT2", "DPP4", "TZD"))



#--- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# Recalibration of the 5 drug model --------------------------------------------
#--- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

 
rm(list=setdiff(ls(), c("data_long_cc", "orig_model", "subset_3treatments", path_summary)))


# Load closet testing procedure (for continuous outcome) 
source(paste0(function_path, "/function_closettesting_procedure_continuous.R"))


# Application of testing procedure

test_results_SGLT2_cohort <- closedtest_continuous_function(cohort         = "SGLT2 cohort", 
                                                            dataset        = subset(subset_3treatments, subset_3treatments$drugclass == "SGLT2"), 
                                                            original_model = orig_model, 
                                                            outcome_name   = "posthba1cfinal", 
                                                            p_value        = 0.05)


test_results_DPP4_cohort  <- closedtest_continuous_function(cohort         = "DPP4 cohort", 
                                                            dataset        = subset(subset_3treatments, subset_3treatments$drugclass == "DPP4"), 
                                                            original_model = orig_model, 
                                                            outcome_name   = "posthba1cfinal", 
                                                            p_value        = 0.05)


test_results_TZD_cohort   <- closedtest_continuous_function(cohort         = "TZD cohort", 
                                                            dataset        = subset(subset_3treatments, subset_3treatments$drugclass == "TZD"), 
                                                            original_model = orig_model, 
                                                            outcome_name   = "posthba1cfinal", 
                                                            p_value        = 0.05)



#--- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# Analysis of overall association to treatment response ------------------------
#--- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

# Analysis: Are methods associated with differential treatment response?

#
## Clinical cluster ------------------------------------------------------------
#

# full model:
model_cluster <- lmer(hba1c_diff ~ treatment_class * clusterID_name + period + (1 | study_id), data = subset_3treatments)
summary(model_cluster)

# reduced model:
model_cluster_reduced <- lmer(hba1c_diff ~ treatment_class + period + (1 | study_id), data = subset_3treatments)

# computing the overall p-value
anova(model_cluster_reduced, model_cluster)


#
## 5 drug model ----------------------------------------------------------------
#

#
### Calculation of the model predicted HbA1c outcome ---------------------------
#

# for prediction period = 1
prediction_data        <- subset_3treatments
prediction_data$period <- rep(0, times = dim(subset_3treatments)[1])                  


# prediction for calibrated model 
# models decisions for the different treatment groups 
test_results_SGLT2_cohort$name_chosen_model
test_results_DPP4_cohort$name_chosen_model
test_results_TZD_cohort$name_chosen_model

subset_3treatments$hba1cprediction_model                                          <- rep(NA, times = dim(subset_3treatments)[1])
subset_3treatments$hba1cprediction_model[subset_3treatments$drugclass == "SGLT2"] <- predict_with_modelchoice_function(test_results_SGLT2_cohort, prediction_data[subset_3treatments$drugclass == "SGLT2", ])
subset_3treatments$hba1cprediction_model[subset_3treatments$drugclass == "DPP4"]  <- predict_with_modelchoice_function(test_results_DPP4_cohort,  prediction_data[subset_3treatments$drugclass == "DPP4", ])
subset_3treatments$hba1cprediction_model[subset_3treatments$drugclass == "TZD"]   <- predict_with_modelchoice_function(test_results_TZD_cohort,   prediction_data[subset_3treatments$drugclass == "TZD", ])


# predictions for uncalibrated model 
subset_3treatments$hba1cprediction_modeluncal        <- rep(NA, times = dim(subset_3treatments)[1])
subset_3treatments$hba1cprediction_modeluncal        <- predict(orig_model, newdata = prediction_data)


#
### Calculation of the model associations --------------------------------------
#

subset_3treatments$hba1c_diff_modelpredicted  <- subset_3treatments$hba1cprediction_model - subset_3treatments$prehba1c

# full model 
model_model <- lmer(hba1c_diff ~ treatment_class * hba1c_diff_modelpredicted + period + (1 | study_id), data = subset_3treatments)
summary(model_model)


# computing the overall p-value
model_model_reduced <- lmer(hba1c_diff ~ treatment_class + period + (1 | study_id), data = subset_3treatments)
anova(model_model_reduced, model_model)


#
## Genetic clusters ------------------------------------------------------------
#

gcluster_names       <- c("t2d_betacell_pineg", "t2d_betacell_pipos", "t2d_bodyfat", "t2d_lipodystrophy", "t2d_liverlipid",
                          "t2d_metabollicsyn", "t2d_obesity", "t2d_residualglycaemic")

gcluster_zscore_names <- paste0(gcluster_names, "_zscore")

interaction           <- paste0("treatment_class*", gcluster_zscore_names)

# Standardization of the genetic risk scores
subset_3treatments[, gcluster_zscore_names] <- NA
subset_3treatments[, gcluster_zscore_names] <- data.frame(apply(subset_3treatments[,  ..gcluster_names], 2, function(x) (x-mean(x, na.rm = T))/sd(x, na.rm = T)))



gcluster_model <- lmer(hba1c_diff ~ treatment_class * t2d_betacell_pineg_zscore + 
                         treatment_class * t2d_betacell_pipos_zscore +
                         treatment_class * t2d_bodyfat_zscore + 
                         treatment_class * t2d_lipodystrophy_zscore +
                         treatment_class * t2d_liverlipid_zscore +
                         treatment_class * t2d_metabollicsyn_zscore +
                         treatment_class * t2d_obesity_zscore +
                         treatment_class * t2d_residualglycaemic_zscore + 
                         + period + (1 | study_id), data = subset_3treatments[complete.cases(subset_3treatments[ , t2d_betacell_pipos_zscore]), ])

summary(gcluster_model)

# reduced model 
gcluster_model_reduced <- lmer(hba1c_diff ~ treatment_class + period + (1 | study_id), data = subset_3treatments[complete.cases(subset_3treatments[ , t2d_betacell_pipos_zscore]), ])


# computing the overall p-value
anova(gcluster_model_reduced, gcluster_model)

#
### Calculation of genetic cluster calculations for plot 2 ---------------------
#







#
## DDR tree --------------------------------------------------------------------
#


load(" ")

names(tri_ddr_dim)[names(tri_ddr_dim) == "dat$study_id"] <- "study_id"

data_DDRmerge            <- merge(subset_3treatments, tri_ddr_dim, by = "study_id")

data_DDRmerge$dim1_z     <- standardize(data_DDRmerge$new_test_dim1)
data_DDRmerge$dim2_z     <- standardize(data_DDRmerge$new_test_dim2)
data_DDRmerge$hba1c_diff <- data_DDRmerge$posthba1cfinal - data_DDRmerge$prehba1c

# full model 
model_DDR <- lmer(hba1c_diff ~ treatment_class*dim1_z +
                    treatment_class*dim2_z + 
                    period + (1 | study_id), data = data_DDRmerge)
summary(model_DDR)
confint(model_DDR)

# reduced model 
model_DDR_reduced <- lmer(hba1c_diff ~ period + (1 | study_id), data = data_DDRmerge)

# computing of the overall p-values
anova(model_DDR_reduced, model_DDR)




#--- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# Concordant/ discordant treatment analysis ------------------------------------
#--- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---


#
## Predictions for cluster strategy --------------------------------------------
#

predictionmodel_cluster              <- lm(posthba1cfinal ~ drugclass + period +
                                                clusterID_name + 
                                                drugclass*clusterID_name, 
                                              data = subset_3treatments)

subset_3treatments$hba1cprediction_cluster <- as.vector(predict(predictionmodel_cluster, newdata = prediction_data))

# (find predictions for model strategy above)


#
## Find optimal treatment (min of predicted HbA1c) and worst treatment ---------
#

# Model calibrated
subset_3treatments <- subset_3treatments %>% group_by(study_id) %>% mutate(min_predHbA1c_model = min(hba1cprediction_model))
subset_3treatments <- subset_3treatments %>% group_by(study_id) %>% mutate(max_predHbA1c_model = max(hba1cprediction_model))

# Model uncalibrated
subset_3treatments <- subset_3treatments %>% group_by(study_id) %>% mutate(min_predHbA1c_modeluncal = min(hba1cprediction_modeluncal))
subset_3treatments <- subset_3treatments %>% group_by(study_id) %>% mutate(max_predHbA1c_modeluncal = max(hba1cprediction_modeluncal))

# Clinical cluster
subset_3treatments <- subset_3treatments %>% group_by(study_id) %>% mutate(min_predHbA1c_cluster = min(hba1cprediction_cluster))
subset_3treatments <- subset_3treatments %>% group_by(study_id) %>% mutate(max_predHbA1c_cluster = max(hba1cprediction_cluster))

# just checking: 
# View(subset(data_long_cc, data_long_cc$study_id == "TM010002"))


#
## Find concordant/ discordant treatment ---------------------------------------
#

# Strategy
# concordant: individuals who received their best treatment according to the cluster or model
# discordant: individuals who received treatment of predicted least good hba1c 

# Model strategy, calibrated
subset_3treatments$group_model <- rep(NA, times = dim(subset_3treatments)[1])
subset_3treatments$group_model[subset_3treatments$min_predHbA1c_model == subset_3treatments$hba1cprediction_model] <- "C"
subset_3treatments$group_model[subset_3treatments$max_predHbA1c_model == subset_3treatments$hba1cprediction_model] <- "D" 

# Model strategy, undcalibrated
subset_3treatments$group_modeluncal <- rep(NA, times = dim(subset_3treatments)[1])
subset_3treatments$group_modeluncal[subset_3treatments$min_predHbA1c_modeluncal == subset_3treatments$hba1cprediction_modeluncal] <- "C"
subset_3treatments$group_modeluncal[subset_3treatments$max_predHbA1c_modeluncal == subset_3treatments$hba1cprediction_modeluncal] <- "D" 

# Cluster strategy
subset_3treatments$group_cluster <- rep(NA, times = dim(subset_3treatments)[1])
subset_3treatments$group_cluster[subset_3treatments$min_predHbA1c_cluster == subset_3treatments$hba1cprediction_cluster] <- "C"
subset_3treatments$group_cluster[subset_3treatments$max_predHbA1c_cluster == subset_3treatments$hba1cprediction_cluster] <- "D"


#
## Choose the analysis cohort from this dataset --------------------------------
#

# Strategy:
# All participants who have taken all three treatments

at_least_three_treatments <- names(which(table(subset_3treatments$study_id) > 2))   # study IDs of participants who have at least 3 treatments
data_analysis_temp        <- subset_3treatments[subset_3treatments$study_id%in%at_least_three_treatments, ]

data_analysis_model      <- data_analysis_temp[complete.cases(data_analysis_temp[ , "group_model"]), ]
data_analysis_modeluncal <- data_analysis_temp[complete.cases(data_analysis_temp[ , "group_modeluncal"]), ]
data_analysis_cluster    <- data_analysis_temp[complete.cases(data_analysis_temp[ , "group_cluster"]), ]


# tidy up the code 
data_analysis_model$group_cluster       <- NULL
data_analysis_model$group_modeluncal    <- NULL

data_analysis_modeluncal$group_cluster  <- NULL
data_analysis_modeluncal$group_model    <- NULL

data_analysis_cluster$group_model       <- NULL
data_analysis_cluster$group_modeluncal  <- NULL




#--- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
#--- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# Main visual Items for paper --------------------------------------------------
#--- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
#--- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

#
## Data preparation ------------------------------------------------------------
# 

subset_3treatments$drugclass  <- factor(subset_3treatments$drugclass, levels = c("SGLT2", "DPP4", "TZD"))

# merge information about model concordant treatment in study data (for calibrated model plots)
model_c_tx                   <- subset(data_analysis_model, data_analysis_model$group_model == "C")[ , c("study_id", "drugclass")]
colnames(model_c_tx)         <- c("study_id", "concordant_tx")
subset_3treatments_wmodelC   <- merge(subset_3treatments, model_c_tx, by = "study_id", all.x = T )

# merge information about model concordant treatment in study data (for uncalibrated model plots)
modeluncal_c_tx                  <- subset(data_analysis_modeluncal, data_analysis_modeluncal$group_modeluncal == "C")[ , c("study_id", "drugclass")]
colnames(modeluncal_c_tx)        <- c("study_id", "concordant_tx")
subset_3treatments_wuncalmodelC  <- merge(subset_3treatments, modeluncal_c_tx, by = "study_id", all.x = T )


#--- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# Table 1: Baseline characteristics of the selected TriMaster cohort -----------
#--- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---


source(paste0(additional_code, "/table1_study_cohorts.R"))


#--- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# Figure 1, Part 1: Main results of genetic clusters ---------------------------
#--- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---


source(paste0(additional_code, "/Figure1Part1.R"))
# plot_results_geneticclusterss.R


#--- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# Figure 1, Part 2: Main results of DDRTree ------------------------------------
#--- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---


source(paste0(additional_code, "/Figure1Part2.R"))
# separate code from John: "Trimaster_DDRTree.R" 


#--- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# Figure 2, Part 1: Main results of clinical clusters  -------------------------
#--- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---


source(paste0(additional_code, "/Figure2Part1.R"))


#--- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# Figure 2, Part 2: Main results of prediction model  --------------------------
#--- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---


source(paste0(additional_code, "/Figure2Part2.R"))


#--- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# Figure 3: Main results comparing clinical cluster and prediction model  ------
#--- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---


source(paste0(additional_code, "/Figure3.R"))


#--- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# Figure S3: Results of clinical clusters  -------------------------------------
#--- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---


source(paste0(additional_code, "/FigureS3.R"))

Figure3_model <- lmer(hba1c_diff ~ group_model + period + (1 | study_id), data = data_analysis_model)
confint(Figure3_model)

Figure3_cluster <- lmer(hba1c_diff ~ group_cluster + period + (1 | study_id), data = data_analysis_cluster)
confint(Figure3_cluster)

                                                            
#--- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# Figure S4: Results of prediction model  --------------------------------------
#--- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---


source(paste0(additional_code, "/FigureS4.R"))


#--- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# Figure S5: Results of prediction model, uncalibrated -------------------------
#--- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---


source(paste0(additional_code, "/FigureS5.R"))


# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --
# Table S4 --------------------------------------------------------------------- 
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --


source(paste0(additional_code, "/TableS4.R"))


