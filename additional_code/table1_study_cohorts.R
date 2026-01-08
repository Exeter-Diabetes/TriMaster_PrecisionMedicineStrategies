#
# Table 1 for study cohort -----------------------------------------------------
#

#
# Packages ---------------------------------------------------------------------
#


# install.packages("table1")                                                    # table 1 cohort description
library(table1)


#
# Setting ----------------------------------------------------------------------
#


table1_cohorts_variables <- c("drugclass","prehba1c", "prebmi", "preegfr", "prehdl", "prealt", "precreat","presys", "predia",
                              "pretotalcholesterol", "Age_at_diagnosis", "sex", "ncurrtx", "HOMA_B", "HOMA_IR", "smoke","t2dmduration", "ethnicity",
                              "imd5", "hba1cmonth")


#
# Table 1: subset_3treatments -------------------------------------------------- 
#


data_table1                     <- subset_3treatments[!duplicated(subset_3treatments$study_id), ]

label(data$prehba1c)            <- "HbA1c (mmol/mol)"
label(data$prebmi)              <- "BMI (kg/m2)"
label(data$preegfr)             <- "eGFR (mL/min/1.73m2)"
label(data$prehdl)              <- "HDL (mg/dL) (1.6 mmol/L)"
label(data$prealt)              <- "ALT (U/L)"
label(data$pretotalcholesterol) <- "total cholesterol"
label(data$precreat)            <- "Creatinine"
label(data$presys)              <- "Systolic blood pressure"
label(data$predia)              <- "Diastolic blood pressure"


label(data$Age_at_diagnosis)    <- "Age (years)"
label(data$sex)                 <- "Gender"
label(data$ncurrtx)             <- "n concurrent tx"
label(data$HOMA_B)              <- "HOMA2-B (%)"
label(data$HOMA_IR)             <- "HOMA2-IR"
label(data$smoke)               <- "Smoke status"
label(data$ethnicity)           <- "Etnicity"
label(data$imd5)                <- "Deprivation index"
label(data$t2dmduration)        <- "T2D duration"
label(data$hba1cmonth)          <- "HbA1c month"


formula_data_description  <- as.formula(paste("~ ", paste0(c(table1_cohorts_variables), collapse =  " + ")))

table_data_description    <- table1(formula_data_description, droplevels = TRUE, render.continuous = "Mean (SD)", data = data_table1)
print(table_data_description)



