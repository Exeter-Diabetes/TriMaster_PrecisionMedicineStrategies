# ------------------------------------------------------------------------------
# Description of generated cluster in TRIMASTER data ---------------------------
# ------------------------------------------------------------------------------


#
# Paths ------------------------------------------------------------------------
#

setwd("C:/Users/lg704/OneDrive - University of Exeter/Documents - Trimaster Analysis/LauraAnalysis/analysis_cluster/study_data_output")

path_cluster_results   <- "C:/Users/lg704/OneDrive - University of Exeter/Documents - Trimaster Analysis/LauraAnalysis/analysis_cluster/02_TRIMASTER_cluster_descriptive_results/"

save_plots             <- "C:/Users/lg704/OneDrive - University of Exeter/Documents - Trimaster Analysis/LauraAnalysis/analysis_cluster/result_plots/"



#
# Packages ---------------------------------------------------------------------
#

# install.packages("ggplot2")                                                   # for plots 
library(ggplot2)

# install.packages("cowplot")                                                   # plot all figures in one panel 
library(cowplot)

# install.packages("Cairo")                                                     # to save the plots in a nice format
library(Cairo)

# install.packages("table1")                                                    # table 1 cohort description
library(table1)

#
# Load data --------------------------------------------------------------------
#


load("data_wcluster.Rdata")                                                     # patient data
# clusterID_name indicates the cluster membership with name of cluster
# cluster_IDall indicates cluster (numberic)


load("studyid_3tx.Rdata")

#
# Data preparation -------------------------------------------------------------
#

# include only N = 309 participants of study cohort 


data             <- data_wcluster[data_wcluster$study_id%in%studyid_3tx$study_id, ] 

rm(data_wcluster)
rm(studyid_3tx)


#
# Cluster characteristics visualization ----------------------------------------
# 

# plots are made on the same scale as Ahlqvist  


cluster_colors <- c("#9bd6e8", "#91bea7", "#fab79c", "#c6a1c2")                 # color for cluster 1: "#9eaed6"
x_axis_lables  <- c("SIDD", "SIRD", "MOD", "MARD")

# HbA1c plot
hba1c <- ggplot(data = data, aes(x = clusterID_name, y = vs_HbA1c_result)) + 
                geom_boxplot(fill = cluster_colors) + 
                stat_summary(fun = mean, colour = "black", geom = "point", shape = 18, size = 5) +
                ylab("HbA1c (mmol/mol)") + xlab("") + theme_bw() +
                scale_y_continuous(breaks = c(50, 100, 150), limits = c(35, 175)) +
                theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), 
                axis.text.x = element_text(angle = 45, vjust = 1.0, hjust = 1, size = 25, face = "bold"),
                axis.text.y = element_text(size = 25, face = "bold"), axis.title.y = element_text(size = 25)) +
                scale_x_discrete(labels = x_axis_lables)


plot(hba1c)

# BMI plot
bmi <- ggplot(data = data, aes(x = clusterID_name, y = vs_BMI)) + 
              geom_boxplot(fill = cluster_colors) + 
              stat_summary(fun = mean, colour = "black", geom = "point", shape = 18, size = 5) +
              ylab("BMI (kg/m2)") + xlab("") + theme_bw() +
              scale_y_continuous(breaks = c(10,20,30,40,50,60,70), limits = c(10,70)) +
              theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), 
              axis.text.x = element_text(angle = 45, vjust = 1.0, hjust = 1, size = 25, face = "bold"),
              axis.text.y = element_text(size = 25, face = "bold"), axis.title.y = element_text(size = 25)) +
              scale_x_discrete(labels = x_axis_lables) 

plot(bmi)


# Age plot 
age <- ggplot(data = data, aes(x = clusterID_name, y = vs_Age_at_diagnosis)) + 
       geom_boxplot(fill = cluster_colors) + 
       stat_summary(fun = mean, colour = "black", geom = "point", shape = 18, size = 5) +
       ylab("Age (years)") + xlab("") + theme_bw() +
       scale_y_continuous(breaks = c(25,50,75), limits = c(12.5,83)) +
        theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), 
        axis.text.x = element_text(angle = 45, vjust = 1.0, hjust = 1, size = 25, face = "bold"),
        axis.text.y = element_text(size = 25, face = "bold"), axis.title.y = element_text(size = 25)) +
        scale_x_discrete(labels = x_axis_lables) 

plot(age)


# HOMA2-B plot
HOMAB <- ggplot(data = data, aes(x = clusterID_name, y = HOMA_B)) + 
         geom_boxplot(fill = cluster_colors) + 
         stat_summary(fun = mean, colour = "black", geom = "point", shape = 18, size = 5) +
         ylab("HOMA2-B (%)") + xlab("") + theme_bw() +
         scale_y_continuous(breaks = c(0,100,200,300), limits = c(0,300)) +
         theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
         axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), 
         axis.text.x = element_text(angle = 45, vjust = 1.0, hjust = 1, size = 25, face = "bold"),
         axis.text.y = element_text(size = 25, face = "bold"), axis.title.y = element_text(size = 25)) +
         scale_x_discrete(labels = x_axis_lables) 

plot(HOMAB)


# HOMA-IR plot
HOMAIR <- ggplot(data = data, aes(x = clusterID_name, y = HOMA_IR)) + 
          geom_boxplot(fill = cluster_colors) + 
          stat_summary(fun = mean, colour = "black", geom = "point", shape = 18, size = 5) +
          ylab("HOMA2-IR") + xlab("") + theme_bw() +
          scale_y_continuous(breaks = c(0,10,20), limits = c(0,26)) +
          theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), 
          axis.text.x = element_text(angle = 45, vjust = 1.0, hjust = 1, size = 25, face = "bold"),
          axis.text.y = element_text(size = 25, face = "bold"), axis.title.y = element_text(size = 25)) +
          scale_x_discrete(labels = x_axis_lables) 

plot(HOMAIR)


#
# Plot all figures together ----------------------------------------------------
#

plot_grid(hba1c, bmi, age, HOMAB, HOMAIR, ncol = 5)




Cairo(file = paste0(save_plots,"cluster_description.png"), 
      type = "png",
      units = "in", 
      width = 20, #10
      height = 8, # 16 
      pointsize = 12, 
      dpi = 300)


plot_grid(hba1c, bmi, age, HOMAB, HOMAIR, ncol = 5)


dev.off()




#
# Cluster description table ----------------------------------------------------
#


variables_of_interest <- c("vs_HbA1c_result", "vs_BMI", "vs_Age_at_diagnosis", "HOMA_B", "HOMA_IR",
                           "vs_Gender")

label(data$vs_HbA1c_result)     <- "HbA1c (mmol/mol)"
label(data$vs_BMI)              <- "BMI (kg/m2)"
label(data$vs_Age_at_diagnosis) <- "Age (years)"
label(data$HOMA_B)              <- "HOMA2-B (%)"
label(data$HOMA_IR)             <- "HOMA2-IR"
label(data$vs_Gender)           <- "Gender"

formula_data_description  <- as.formula(paste("~ ", paste0(c(variables_of_interest), collapse =  " + "), "| clusterID_name "))
table_data_description    <- table1(formula_data_description, droplevels = TRUE, render.continuous = "Mean (SD)", data = data)
print(table_data_description)


#
# Save all figures -------------------------------------------------------------
#

Cairo(file = paste0(path_cluster_results, "cluster_description_TRIMASTER_hba1c.png"), 
      type = "png",
      units = "in", 
      width = 8, #10
      height = 10, # 16 
      pointsize = 12, 
      dpi = 300)

plot(hba1c)

dev.off()


Cairo(file = paste0(path_cluster_results, "cluster_description_TRIMASTER_bmi.png"), 
      type = "png",
      units = "in", 
      width = 8, #10
      height = 10, # 16 
      pointsize = 12, 
      dpi = 300)

plot(bmi)

dev.off()


Cairo(file = paste0(path_cluster_results, "cluster_description_TRIMASTER_age.png"), 
      type = "png",
      units = "in", 
      width = 8, #10
      height = 10, # 16 
      pointsize = 12, 
      dpi = 300)

plot(age)

dev.off()


Cairo(file = paste0(path_cluster_results, "cluster_description_TRIMASTER_HOMAB.png"), 
      type = "png",
      units = "in", 
      width = 8, #10
      height = 10, # 16 
      pointsize = 12, 
      dpi = 300)

plot(HOMAB)

dev.off()



Cairo(file = paste0(path_cluster_results, "cluster_description_TRIMASTER_HOMAIR.png"), 
      type = "png",
      units = "in", 
      width = 8, #10
      height = 10, # 16 
      pointsize = 12, 
      dpi = 300)

plot(HOMAIR)

dev.off()



Cairo(file = paste0(path_cluster_results, "cluster_description_TRIMASTER.jpeg"), 
      type = "jpeg",
      units = "in", 
      width = 20, #20
      height = 8, # 8 
      pointsize = 12, 
      dpi = 300)
plot_grid(hba1c, bmi, age, HOMAB, HOMAIR, ncol = 5)

dev.off()





