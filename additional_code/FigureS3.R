#
## Part A ----------------------------------------------------------------------
#

model_FigureA <- lmer(hba1c_diff ~ treatment_class * clusterID_name + prehba1c + period + (1 | study_id), data = subset_3treatments)

emm_FigureA      <- emmeans(model_FigureA, ~ treatment_class | clusterID_name, cov.reduce = mean)
data_FigureA     <- as.data.frame(emm_FigureA[1:12,1])


FigureA_cluster2 <- data_FigureA[which(data_FigureA$clusterID_name == "2_SIDD"), ] %>% 
  ggplot(aes(x = clusterID_name, y = emmean, fill = treatment_class))+
  geom_col( position = "dodge", width = 0.5, alpha = 0.5, color = "black", linewidth = 0.1)+
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                position =  position_dodge(width = 0.5), width = 0.2, linewidth = 1)+
  scale_fill_manual(values = c("#DA9700", "#005686", "#C33D13"), name = NULL) +
  scale_y_continuous(limits = c(-17.5, 0)) + 
  xlab("") +
  ylab ("") + # Achieved HbA1c change (mmol/mol)
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        legend.position = "none", axis.text.x=element_blank(), axis.title=element_text(size = 15),
        axis.text = element_text(size = 15))

FigureA_cluster3 <- data_FigureA[which(data_FigureA$clusterID_name == "3_SIRD"), ] %>% 
  ggplot(aes(x = clusterID_name, y = emmean, fill = treatment_class))+
  geom_col( position = "dodge", width = 0.5, alpha = 0.5, color = "black", linewidth = 0.1)+
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                position =  position_dodge(width = 0.5), width = 0.2, linewidth = 1)+
  scale_fill_manual(values = c("#DA9700", "#005686", "#C33D13"), name = NULL) +
  scale_y_continuous(limits = c(-17.5, 0)) + 
  xlab("") +
  ylab (" ") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        legend.position = "none", axis.text.x=element_blank(), axis.title=element_text(size = 15),
        axis.text = element_text(size = 15))

FigureA_cluster4 <- data_FigureA[which(data_FigureA$clusterID_name == "4_MOD"), ] %>% 
  ggplot(aes(x = clusterID_name, y = emmean, fill = treatment_class))+
  geom_col( position = "dodge", width = 0.5, alpha = 0.5, color = "black", linewidth = 0.1)+
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                position =  position_dodge(width = 0.5), width = 0.2, linewidth = 1)+
  scale_fill_manual(values = c("#DA9700", "#005686", "#C33D13"), name = NULL) +
  scale_y_continuous(limits = c(-17.5, 0)) + 
  xlab("") +
  ylab (" ") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        legend.position = "none", axis.text.x=element_blank(), axis.title=element_text(size = 15),
        axis.text = element_text(size = 15))

FigureA_cluster5 <- data_FigureA[which(data_FigureA$clusterID_name == "5_MARD"), ] %>% 
  ggplot(aes(x = clusterID_name, y = emmean, fill = treatment_class))+
  geom_col( position = "dodge", width = 0.5, alpha = 0.5, color = "black", linewidth = 0.1)+
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                position =  position_dodge(width = 0.5), width = 0.2, linewidth = 1)+
  scale_fill_manual(values = c("#DA9700", "#005686", "#C33D13"), name = NULL) +
  scale_y_continuous(limits = c(-17.5, 0)) + 
  xlab("") +
  ylab (" ") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        legend.position = "none", axis.text.x=element_blank(), axis.title=element_text(size = 15),
        axis.text = element_text(size = 15))



plot_grid(FigureA_cluster2, FigureA_cluster3, FigureA_cluster4, FigureA_cluster5, ncol = 4)



#
## pairwise t-test -------------------------------------------------------------
#

pairwise_comparison_Figure <- contrast(emm_FigureA, method = "pairwise", by = "clusterID_name", adjust = "bonferroni")

#
## Part B ----------------------------------------------------------------------
#

FigureB_data           <- data.frame(pairwise_comparison_Figure)[c(1, 2, 4, 5, 7, 8, 10, 11), ]
FigureB_data$drugclass <- c("DPP4", "TZD", "DPP4", "TZD", "DPP4", "TZD", "DPP4", "TZD") 
FigureB_data$CI_low    <- FigureB_data$estimate - 1.96*FigureB_data$SE
FigureB_data$CI_high   <- FigureB_data$estimate + 1.96*FigureB_data$SE

FigureB_plotdata          <- FigureB_data
FigureB_plotdata$estimate <- -1*FigureB_data$estimate
FigureB_plotdata$SE       <- -1*FigureB_data$SE
FigureB_plotdata$CI_low   <- -1*FigureB_data$CI_low
FigureB_plotdata$CI_high  <- -1*FigureB_data$CI_high


FigureB_cluster2 <- FigureB_plotdata[which(FigureB_plotdata$clusterID_name == "2_SIDD"), ]%>% 
  ggplot(aes(x = clusterID_name, y = estimate, fill = drugclass))+
  geom_col( position = "dodge", width = 0.5, alpha = 0.5, color = "black", linewidth = 0.1)+
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high),
                position =  position_dodge(width = 0.5), width = 0.2, linewidth = 1)+
  scale_fill_manual(values = c("#005686", "#C33D13"), name = NULL) +
  geom_hline(yintercept = 0, color = "#DA9700", linewidth = 0.8) + 
  scale_y_continuous(limits = c(-6, 6)) + 
  xlab("") +
  ylab (" ") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        legend.position = "none", axis.text.x=element_blank(), axis.title=element_text(size = 15),
        axis.text = element_text(size = 15))

FigureB_cluster3 <- FigureB_plotdata[which(FigureB_plotdata$clusterID_name == "3_SIRD"), ]%>% 
  ggplot(aes(x = clusterID_name, y = estimate, fill = drugclass))+
  geom_col( position = "dodge", width = 0.5, alpha = 0.5, color = "black", linewidth = 0.1)+
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high),
                position =  position_dodge(width = 0.5), width = 0.2, linewidth = 1)+
  scale_fill_manual(values = c("#005686", "#C33D13"), name = NULL) +
  geom_hline(yintercept = 0, color = "#DA9700", linewidth = 0.8) + 
  scale_y_continuous(limits = c(-6, 6)) + 
  xlab("") +
  ylab (" ") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        legend.position = "none", axis.text.x=element_blank(), axis.title=element_text(size = 15),
        axis.text = element_text(size = 15))

FigureB_cluster4 <- FigureB_plotdata[which(FigureB_plotdata$clusterID_name == "4_MOD"), ]%>% 
  ggplot(aes(x = clusterID_name, y = estimate, fill = drugclass))+
  geom_col( position = "dodge", width = 0.5, alpha = 0.5, color = "black", linewidth = 0.1)+
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high),
                position =  position_dodge(width = 0.5), width = 0.2, linewidth = 1)+
  scale_fill_manual(values = c("#005686", "#C33D13"), name = NULL) +
  geom_hline(yintercept = 0, color = "#DA9700", linewidth = 0.8) + 
  scale_y_continuous(limits = c(-6, 6)) + 
  xlab("") +
  ylab (" ") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        legend.position = "none", axis.text.x=element_blank(), axis.title=element_text(size = 15),
        axis.text = element_text(size = 15))

FigureB_cluster5 <- FigureB_plotdata[which(FigureB_plotdata$clusterID_name == "5_MARD"), ]%>% 
  ggplot(aes(x = clusterID_name, y = estimate, fill = drugclass))+
  geom_col( position = "dodge", width = 0.5, alpha = 0.5, color = "black", linewidth = 0.1)+
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high),
                position =  position_dodge(width = 0.5), width = 0.2, linewidth = 1)+
  scale_fill_manual(values = c("#005686", "#C33D13"), name = NULL) +
  geom_hline(yintercept = 0, color = "#DA9700", linewidth = 0.8) + 
  scale_y_continuous(limits = c(-6, 6)) + 
  xlab("") +
  ylab (" ") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        legend.position = "none", axis.text.x=element_blank(), axis.title=element_text(size = 15),
        axis.text = element_text(size = 15))


plot_grid(FigureB_cluster2, FigureB_cluster3, FigureB_cluster4, FigureB_cluster5, ncol = 4)


#
## Combination of all plots ----------------------------------------------------
#


Cairo(file = paste0(plot_path,"/main_results_adjusted_cluster.png" ), 
      type = "png",
      units = "in", 
      width = 18,
      height = 8, 
      pointsize = 12, 
      dpi = 72)

plot_grid(FigureA_cluster2, FigureA_cluster3, FigureA_cluster4, FigureA_cluster5, 
          FigureB_cluster2, FigureB_cluster3, FigureB_cluster4, FigureB_cluster5, 
          ncol = 4)

dev.off()