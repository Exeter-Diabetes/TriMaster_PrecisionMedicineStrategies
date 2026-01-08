
##### TRIMASTER PREP #####

#library(haven)
library(tidyverse)
library(rms)
library(viridis)
#library(caret)
library(data.table)
library(lme4)
library(splines)
library(ggthemes)
library(MuMIn)
library(tableone)
library(stargazer)
library(cowplot)
library(haven)
library(patchwork)

setwd("")
path_data  <- ""

load(paste0(path_data, "/data_wcluster.Rdata"))                                 # patient data
smoke_data_merge_orig <- data.frame(read_dta(paste0("")))
smoke_data_merge      <- data.frame(study_id = as.vector(smoke_data_merge_orig$study_id), smoke = smoke_data_merge_orig$v1_Smoker)
dat <- merge(data_wcluster, smoke_data_merge, by = "study_id", all.x = TRUE)
rm(data_wcluster)
rm(smoke_data_merge_orig)
rm(smoke_data_merge)


head(dat)

dat <- dat %>% filter(randomised==1)

#Trimaster groups
table(dat$obese)
table(dat$egfr90)

#Included patients
table(dat$incp1)
table(dat$incp2)
table(dat$incp3)



#Work out HbA1c on each drug
dat <- dat %>% mutate(drug1=Stage1,drug2=Stage2,drug3=Stage3)

dat <- dat %>% mutate(validhba1cA=if_else(drug1=="A" & incp1==1,hba1c1,
                                          if_else(drug2=="A" & incp2==1,hba1c2,
                                                  if_else(drug3=="A" & incp3==1,hba1c3,NA_real_))))
dat <- dat %>% mutate(validhba1cB=if_else(drug1=="B" & incp1==1,hba1c1,
                                          if_else(drug2=="B" & incp2==1,hba1c2,
                                                  if_else(drug3=="B" & incp3==1,hba1c3,NA_real_))))
dat <- dat %>% mutate(validhba1cC=if_else(drug1=="C" & incp1==1,hba1c1,
                                          if_else(drug2=="C" & incp2==1,hba1c2,
                                                  if_else(drug3=="C" & incp3==1,hba1c3,NA_real_))))

#how many have valid on both drugs (for each pairwise comparison)?

dat <- dat %>% mutate(combvalid =
                        ifelse(is.na(validhba1cA) & is.na(validhba1cB) & is.na(validhba1cC),"none",
                               ifelse(!is.na(validhba1cA) & is.na(validhba1cB) & is.na(validhba1cC),"A only",
                                      ifelse(is.na(validhba1cA) & !is.na(validhba1cB) & is.na(validhba1cC),"B only",
                                             ifelse(is.na(validhba1cA) & is.na(validhba1cB) & !is.na(validhba1cC),"C only",
                                                    ifelse(!is.na(validhba1cA) & !is.na(validhba1cB) & is.na(validhba1cC),"A&B only",
                                                           ifelse(is.na(validhba1cA) & !is.na(validhba1cB) & !is.na(validhba1cC),"B&C only",
                                                                  ifelse(!is.na(validhba1cA) & is.na(validhba1cB) & !is.na(validhba1cC),"A&C only",
                                                                         ifelse(!is.na(validhba1cA) & !is.na(validhba1cB) & !is.na(validhba1cC),"All three",
                                                                                NA)))))))))
describe(dat$combvalid)                      
table(dat$combvalid,dat$baseline)                      

dat <- dat %>% dplyr::select(study_id, drug1,drug2,drug3,validhba1cA,validhba1cB,validhba1cC, 
                      egfr90,obese,
                      Age_at_diagnosis=vs_Age_at_diagnosis,
                      ageatscreening,
                      sex=vs_Gender,
                      prebmi=bmi,
                      preegfr=egfr,
                      prealt=vs_ALT_result,
                      pretotalcholesterol=CHOLmmolLV1,
                      prehdl=HDLmmolLV1,
                      prehba1c=vs_HbA1c_result,
                      smoke,
                      v2v1datediff, v3v2datediff, v4v3datediff, 
                      ethnicity=vs_Ethnic_Group, 
                      HOMA_B, 
                      HOMA_IR,
                      pretrig = TRIGmmolLV1,
                      precreat = CREATumolLV1,
                      presys = v1_Blood_Pressure_1sys,
                      predia = v1_Blood_Pressure_1dia,
                      vs_MFN, 
                      vs_SU,
                      combvalid,
                      baseline,
                      incp1,incp2,incp3,
                      clusterID_name, 
                      TreatmentOrder)
head(dat)
describe(dat)

# change coding of ethnicity variable to fit with model
table(dat$ethnicity)
dat$ethnicity[dat$ethnicity == "A"] <- "White"
dat$ethnicity[dat$ethnicity == "B"] <- "White"
dat$ethnicity[dat$ethnicity == "C"] <- "White"
dat$ethnicity[dat$ethnicity == "E"] <- "Mixed"
dat$ethnicity[dat$ethnicity == "F"] <- "Mixed"
dat$ethnicity[dat$ethnicity == "H"] <- "South Asian"
dat$ethnicity[dat$ethnicity == "J"] <- "South Asian"
dat$ethnicity[dat$ethnicity == "K"] <- "South Asian"
dat$ethnicity[dat$ethnicity == "L"] <- "South Asian"
dat$ethnicity[dat$ethnicity == "N"] <- "Black"
dat$ethnicity[dat$ethnicity == "R"] <- "Other"
dat$ethnicity[dat$ethnicity == "S"] <- "Other"
dat$ethnicity[dat$ethnicity == "Z"] <- "Missing"

dat$ethnicity                             <- factor(dat$ethnicity, levels = c("White", "South Asian", "Black", "Mixed", "Other", "Missing"))

# deprivation index variable
dat$imd5                            <- rep("3", times = dim(dat)[1])
dat$imd5                            <- factor(dat$imd5, levels = c("1 (least)", "2", "3", "4", "5 (most)"))

# age at treatment variable
dat$agetx                           <- dat$ageatscreening  

# sex variable needs to be Female/Male, female is baseline
dat$sex[dat$sex == "M"]        <- "Male"
dat$sex[dat$sex == "F"]        <- "Female"

str(dat$sex)
dat$sex                              <- as.factor(dat$sex)

# T2D duration variable
dat$t2dmduration                     <- dat$agetx - dat$Age_at_diagnosis
summary(dat$t2dmduration)

# recode smoke variable
dat$smoke[dat$smoke == 0]   <- "Non-smoker"
dat$smoke[dat$smoke == 1]   <- "Active smoker"
dat$smoke[is.na(dat$smoke)] <- "Not recorded"
table(dat$smoke)
dat$smoke                         <- factor(dat$smoke, levels = c("Non-smoker", "Active smoker", "Ex-smoker", "Not recorded"))
str(dat$smoke)
table(dat$smoke)



head(dat)
describe(dat)

#Define period for each drug
dat <- dat %>% mutate(period.sglt2=ifelse(drug1=="C","1",
                                          ifelse(drug2=="C","2",
                                                 ifelse(drug3=="C","3",
                                                        NA))))
describe(dat$period.sglt2)

dat <- dat %>% mutate(period.dpp4=ifelse(drug1=="B","1",
                                         ifelse(drug2=="B","2",
                                                ifelse(drug3=="B","3",
                                                       NA))))
describe(dat$period.dpp4)

dat <- dat %>% mutate(period.pio=ifelse(drug1=="A","1",
                                        ifelse(drug2=="A","2",
                                               ifelse(drug3=="A","3",
                                                      NA))))
describe(dat$period.pio)

#Drugline (number of current treatments + treatment initiated)
dat$drugline <- dat$vs_MFN + dat$vs_SU + 1
# create drugline as factor
dat$drugline                          <- factor(dat$drugline, levels = c("2", "3", "4", "5+"))
table(dat$drugline)

#ncurrtx (number of concomitant treatments)
dat$ncurrtx <- dat$vs_MFN + dat$vs_SU
head(dat)
# create ncurrtx as factor
dat$ncurrtx                           <- factor(dat$ncurrtx, levels = c("1", "2", "3", "4+"))
table(dat$ncurrtx)

#Drugclass dummy
dat$drugclass <- "DPP4"

#Subset n=309
load("C:/Users/lg704/OneDrive - University of Exeter/Documents - Trimaster Analysis/LauraAnalysis/analysis_cluster/study_data_output/studyid_3tx.Rdata")
dat <- merge(dat,studyid_3tx,by="study_id")

#################
#DDR TREE MAPPING TO SCOTTISH TREE
#################

load("")
load("")
load("")

dat$prehba1c.pc <- (dat$prehba1c/10.929) + 2.15

data.ddr <- dat %>% relocate(sex, agetx, prehdl, pretotalcholesterol, pretrig, prehba1c.pc, prebmi, presys, predia, prealt, precreat,
                             validhba1cA, validhba1cB, validhba1cC)
data.ddr <- data.ddr[1:14]
data.ddr <- data.ddr %>% dplyr::filter(complete.cases(data.ddr[,c(1:11)]))

ddrtree_map<-function(data){
  
  library(mgcv)
  library(ggplot2)
  
  #data <- pca_ext_formapping
  data <- data[1:11]

  # prepare the data
  names(data)[1]<-'sex.x';names(data)[2]<-'clu_age_at_diag';names(data)[3]<-'hdl';names(data)[4]<-'chol';
  names(data)[5]<-'trigs';names(data)[6]<-'hba1c_dcct';names(data)[7]<-'bmi_final';names(data)[8]<-'sbp';
  names(data)[9]<-'dbp';names(data)[10]<-'alt';names(data)[11]<-'creat'
  
  print(names(data))
  
  # Data class
  data[,c(2:9)]<-apply(data[,c(2:9)],2,function(x) as.numeric(as.character(x)))
  data[,c(10,11)]<-apply(data[,c(10,11)],2,function(x) as.integer(as.character(x)))
  #data$sex.x<-factor(data$sex.x,levels=c('M','F'))
  data$sex.x<-factor(ifelse(data$sex.x=="Male","M","F"))
  head(data)
  
  # get the position predict
  dim1<-predict.gam(fit_dim1,newdata = data)
  dim2<-predict.gam(fit_dim2,newdata = data)
  #return(test)
  predict_data<-as.data.frame(cbind(dim1,dim2))
  #return(predict_data)
  
  ## Re position the predicted points based on Euclidean distance 
  for(i in 1:dim(predict_data)[1]){
    distance<-((nl1$data_dim_1-predict_data$dim1[i])^2+(nl1$data_dim_2-predict_data$dim2[i])^2)^0.5
    j<-which.min(distance)
    #print(j)
    predict_data$new_test_dim1[i]<-nl1$data_dim_1[j]
    predict_data$new_test_dim2[i]<-nl1$data_dim_2[j]
  }
  
  ### output
  k<-ggplot(nl1,aes(x=data_dim_1,y=data_dim_2))+geom_point(color='grey')+
    geom_point(data=predict_data,aes(x=new_test_dim1,y=new_test_dim2),color='red')+
    xlab('Dimension 1')+ylab('Dimension 2')+theme_minimal()+
    labs(title = "Positions of new individuals")+
    theme(plot.title = element_text(size = 15, face = "bold"))
  print(k)
  return(predict_data)
  return(k)
}

TRI_mapped <- ddrtree_map(data.ddr)
k<-ggplot(nl1,aes(x=data_dim_1,y=data_dim_2))+geom_point(color='grey')+
  geom_point(data=ddrtree_map(data.ddr),aes(x=new_test_dim1,y=new_test_dim2),color='red')+
  xlab('Dimension 1')+ylab('Dimension 2')+theme_minimal()+
  labs(title = "Positions of new individuals")+
  theme(plot.title = element_text(size = 15, face = "bold"))

TRI_mapped <- cbind(TRI_mapped,data.ddr)
TRI_mapped$prehba1c <- 10.929*(TRI_mapped$prehba1c.pc-2.15)
TRI_mapped$TZD <- TRI_mapped$validhba1cA - TRI_mapped$prehba1c
TRI_mapped$SGLT2 <- TRI_mapped$validhba1cC - TRI_mapped$prehba1c
TRI_mapped$DPP4 <- TRI_mapped$validhba1cB - TRI_mapped$prehba1c
hist(TRI_mapped$TZD)
hist(TRI_mapped$SGLT2)
hist(TRI_mapped$DPP4)

tzd.data <- TRI_mapped %>% dplyr::filter(complete.cases(TRI_mapped$TZD))
sglt2.data <- TRI_mapped %>% dplyr::filter(complete.cases(TRI_mapped$SGLT2))
dpp4.data <- TRI_mapped %>% dplyr::filter(complete.cases(TRI_mapped$DPP4))

#residuals
m.tzd <- lm(TZD~rcs(prehba1c,3),data=tzd.data)
m.sglt2 <- lm(SGLT2~rcs(prehba1c,3),data=sglt2.data)
m.dpp4 <- lm(DPP4~rcs(prehba1c,3),data=dpp4.data)
tzd.data$resid <- resid(m.tzd)
sglt2.data$resid <- resid(m.sglt2)
dpp4.data$resid <- resid(m.dpp4)

TRI_mapped$SGLT2DPP4 <- TRI_mapped$SGLT2 - TRI_mapped$DPP4
TRI_mapped$SGLT2TZD <- TRI_mapped$SGLT2 - TRI_mapped$TZD
TRI_mapped$DPP4TZD <- TRI_mapped$DPP4 - TRI_mapped$TZD
head(TRI_mapped)
hist(TRI_mapped$SGLT2DPP4)
hist(TRI_mapped$SGLT2TZD)
hist(TRI_mapped$DPP4TZD)

sglt2.dpp4.data <- TRI_mapped %>% dplyr::filter(complete.cases(TRI_mapped$SGLT2DPP4))
sglt2.tzd.data <- TRI_mapped %>% dplyr::filter(complete.cases(TRI_mapped$SGLT2TZD))
dpp4.tzd.data <- TRI_mapped %>% dplyr::filter(complete.cases(TRI_mapped$DPP4TZD))

#Differential response: SGLT2 v DPP4
  psglt2dpp4 <-  ggplot(nl1,aes(x=data_dim_1,y=data_dim_2))+geom_point(color='grey')+
    geom_point(data=sglt2.dpp4.data,aes(x=new_test_dim1,y=new_test_dim2,color=SGLT2DPP4))+
    scale_color_gradient(low="green",high="darkmagenta", limits = c(-50, 25), name = " ")+
    xlab('Dimension 1')+ylab('Dimension 2')+theme_minimal()+
    labs(title = "Positions of new individuals")+
    theme(plot.title = element_text(size = 15, face = "bold")) + ggtitle("SGLT2i - DPP4")
  
  m1 <- lm(SGLT2DPP4~ns(new_test_dim1,3)+ns(new_test_dim2,3),data=sglt2.dpp4.data)
  summary(m1)
  anova(m1)
 
#Differential response: TZD v DPP4
  pdpp4tzd <-  ggplot(nl1,aes(x=data_dim_1,y=data_dim_2))+geom_point(color='grey')+
    geom_point(data=dpp4.tzd.data,aes(x=new_test_dim1,y=new_test_dim2,color=DPP4TZD))+
    scale_color_gradient(low="green",high="darkmagenta", limits = c(-50, 25), name = " ")+
    xlab('Dimension 1')+ylab('Dimension 2')+theme_minimal()+
    labs(title = "Positions of new individuals")+
    theme(plot.title = element_text(size = 15, face = "bold")) + ggtitle("DPP4i - TZD")
  
  m1 <- lm(DPP4TZD~ns(new_test_dim1,3)+ns(new_test_dim2,3),data=dpp4.tzd.data)
  summary(m1)
  anova(m1)
 
#Differential response: SGLT2 v TZD
  psglt2tzd <-  ggplot(nl1,aes(x=data_dim_1,y=data_dim_2))+geom_point(color='grey')+
    geom_point(data=sglt2.tzd.data,aes(x=new_test_dim1,y=new_test_dim2,color=SGLT2TZD))+
    scale_color_gradient(low="green",high="darkmagenta", limits = c(-50, 25), name = "treatment response")+ #  c(-50, 25)
    xlab('Dimension 1')+ylab('Dimension 2')+theme_minimal()+
    labs(title = "Positions of new individuals")+
    theme(plot.title = element_text(size = 15, face = "bold")) + ggtitle("SGLT2i - TZD")

  m1 <- lm(SGLT2TZD~ns(new_test_dim1,3)+ns(new_test_dim2,3),data=sglt2.tzd.data)
  summary(m1)
  anova(m1)

#save
  pdf.options(reset=TRUE, onefile=FALSE)
  pdf("/ddr_phenotype_TRI_diffdrugresponse_n309.pdf",width=15,height=5)
  psglt2dpp4 + pdpp4tzd + psglt2tzd + plot_layout(nrow=1)
  dev.off()
  


#Save dim data
tri_ddr_dim <- TRI_mapped %>% select(new_test_dim1, new_test_dim2) 
tri_ddr_dim <- cbind(dat$study_id,tri_ddr_dim)
save(tri_ddr_dim, file="/tri_ddr_dim.Rdata")

head(sglt2.dpp4.data)

#
# Calculation of the p-values --------------------------------------------------
#

variables_to_include <- c("study_id", "validhba1cA", "validhba1cB", "validhba1cC", "prehba1c")


# reshape dataset 
data_long           <- melt(setDT(dat[, variables_to_include]), 
                            id.vars = c("study_id", "prehba1c"), 
                            variable.name = c("hba1c"))

# treatment variable
data_long$tx_taken                                 <- rep(NA, times = dim(data_long)[1])
data_long$tx_taken[data_long$hba1c == "validhba1cA"] <- "A"
data_long$tx_taken[data_long$hba1c == "validhba1cB"] <- "B"
data_long$tx_taken[data_long$hba1c == "validhba1cC"] <- "C"

data_long$drugclass                                      <- rep(NA, times = dim(data_long)[1])
data_long$drugclass[data_long$hba1c == "validhba1cA"] <- "TZD"
data_long$drugclass[data_long$hba1c == "validhba1cB"] <- "DPP4"
data_long$drugclass[data_long$hba1c == "validhba1cC"] <- "SGLT2"

names(data_long)[names(data_long) == "value"] <- "hba1coutcome"


# merge with DDRtree results 

names(tri_ddr_dim)[names(tri_ddr_dim) == "dat$study_id"] <- "study_id"

merged_data <- merge(data_long, tri_ddr_dim, by=c("study_id"))


# subsets data in treatment comparisons 

subset_SGLT2DPP4                <- subset(merged_data, merged_data$drugclass == "SGLT2" | merged_data$drugclass == "DPP4")
subset_SGLT2DPP4$achieved_hba1c <- subset_SGLT2DPP4$hba1coutcome - subset_SGLT2DPP4$prehba1c

subset_DPP4TZD                <- subset(merged_data, merged_data$drugclass == "DPP4" | merged_data$drugclass == "TZD")
subset_DPP4TZD$achieved_hba1c <- subset_DPP4TZD$hba1coutcome - subset_DPP4TZD$prehba1c

subset_SGLT2TZD                <- subset(merged_data, merged_data$drugclass == "SGLT2" | merged_data$drugclass == "TZD")
subset_SGLT2TZD$achieved_hba1c <- subset_SGLT2TZD$hba1coutcome - subset_SGLT2TZD$prehba1c


# calculation of hba1cdifference

diff_in_diff_SGLT2DPP4_data <- subset_SGLT2DPP4 %>%
  group_by(study_id) %>%
  mutate(diff_in_diff = lead(achieved_hba1c) - achieved_hba1c) %>%
  fill(diff_in_diff)


diff_in_diff_DPP4TZD_data <- subset_DPP4TZD %>%
  group_by(study_id) %>%
  mutate(diff_in_diff = lead(achieved_hba1c) - achieved_hba1c) %>%
  fill(diff_in_diff)


diff_in_diff_SGLT2TZD_data <- subset_SGLT2TZD %>%
  group_by(study_id) %>%
  mutate(diff_in_diff = lead(achieved_hba1c) - achieved_hba1c) %>%
  fill(diff_in_diff)


# regression to derive the p-values:

library(splines)

model_SGLT2DPP4 <- lm(diff_in_diff ~ ns(new_test_dim1,1) + ns(new_test_dim2,1), data = unique(diff_in_diff_SGLT2DPP4_data[, c("diff_in_diff", "new_test_dim1", "new_test_dim2")]))
summary(model_SGLT2DPP4)

model_DPP4TZD <- lm(diff_in_diff ~ ns(new_test_dim1,1) + ns(new_test_dim2,1), data = unique(diff_in_diff_DPP4TZD_data[, c("diff_in_diff", "new_test_dim1", "new_test_dim2")]))
summary(model_DPP4TZD)

model_SGLT2TZD <- lm(diff_in_diff ~ ns(new_test_dim1,1) + ns(new_test_dim2,1), data = unique(diff_in_diff_SGLT2TZD_data[, c("diff_in_diff", "new_test_dim1", "new_test_dim2")]))
summary(model_SGLT2TZD)






