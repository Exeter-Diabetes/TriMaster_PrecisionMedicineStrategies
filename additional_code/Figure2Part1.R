
#
## Part A ----------------------------------------------------------------------
#

Figure_data_full       <- subset_3treatments[, c("clusterID_name", "drugclass", "hba1c_diff")]

FigureA_data           <- Figure_data_full %>% 
  group_by(clusterID_name, drugclass) %>% 
  summarise(mean= mean(hba1c_diff),
            se = 1.96* sd(hba1c_diff)/sqrt(n()))


FigureA_cluster2 <- FigureA_data[which(FigureA_data$clusterID_name == "2_SIDD"), ] %>%  
  ggplot(aes(x = clusterID_name, y = mean, fill = drugclass))+
  geom_col( position = "dodge", width = 0.5, alpha = 0.5, color = "black", linewidth = 0.1)+
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position =  position_dodge(width = 0.5), width = 0.2, linewidth = 1)+
  scale_fill_manual(values = c("#DA9700","#005686", "#C33D13"), name = NULL) +
  scale_y_continuous(limits = c(-22.5, 0)) + 
  xlab("") +
  ylab (" ") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        legend.position = "none", axis.text.x=element_blank(), axis.title=element_text(size = 15),
        axis.text = element_text(size = 15))


FigureA_cluster3 <- FigureA_data[which(FigureA_data$clusterID_name == "3_SIRD"), ] %>%  
  ggplot(aes(x = clusterID_name, y = mean, fill = drugclass))+
  geom_col( position = "dodge", width = 0.5, alpha = 0.5, color = "black", linewidth = 0.1)+
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position =  position_dodge(width = 0.5), width = 0.2, linewidth = 1)+
  scale_y_continuous(limits = c(-22.5, 0)) + 
  xlab("") +
  ylab (" ") + 
  scale_fill_manual(values = c("#DA9700","#005686", "#C33D13"), name = NULL) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        legend.position = "none", axis.text.x=element_blank(), axis.title=element_text(size = 15),
        axis.text = element_text(size = 15))


FigureA_cluster4 <- FigureA_data[which(FigureA_data$clusterID_name == "4_MOD"), ] %>%  
  ggplot(aes(x = clusterID_name, y = mean, fill = drugclass))+
  geom_col( position = "dodge", width = 0.5, alpha = 0.5, color = "black", linewidth = 0.1)+
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position =  position_dodge(width = 0.5), width = 0.2, linewidth = 1)+
  scale_y_continuous(limits = c(-22.5, 0)) + 
  xlab("") +
  ylab (" ") + 
  scale_fill_manual(values = c("#DA9700","#005686", "#C33D13"), name = NULL) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        legend.position = "none", axis.text.x=element_blank(), axis.title=element_text(size = 15),
        axis.text = element_text(size = 15))


FigureA_cluster5 <- FigureA_data[which(FigureA_data$clusterID_name == "5_MARD"), ] %>%  
  ggplot(aes(x = clusterID_name, y = mean, fill = drugclass))+
  geom_col( position = "dodge", width = 0.5, alpha = 0.5, color = "black", linewidth = 0.1)+
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position =  position_dodge(width = 0.5), width = 0.2, linewidth = 1)+
  scale_y_continuous(limits = c(-22.5, 0)) + 
  xlab("") +
  ylab (" ") + 
  scale_fill_manual(values = c("#DA9700","#005686", "#C33D13"), name = NULL) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        legend.position = "none", axis.text.x=element_blank(), axis.title=element_text(size = 15),
        axis.text = element_text(size = 15))


plot_grid(FigureA_cluster2, FigureA_cluster3, FigureA_cluster4, FigureA_cluster5, ncol = 4)



#
## pairwise t-test -------------------------------------------------------------
#


Figure_SIDD2_diff <- subset(subset_3treatments, subset_3treatments$clusterID_name == "2_SIDD")[ , c("treatment_class", "hba1c_diff")] %>% 
  pairwise_t_test(hba1c_diff ~ treatment_class, p.adjust.method = "bonferroni", paired = TRUE, pool.sd = FALSE, detailed = TRUE)

Figure_SIRD3_diff <- subset(subset_3treatments, subset_3treatments$clusterID_name == "3_SIRD")[ , c("treatment_class", "hba1c_diff")] %>% 
  pairwise_t_test(hba1c_diff ~ treatment_class, p.adjust.method = "bonferroni", paired = TRUE, pool.sd = FALSE, detailed = TRUE)

Figure_MOD4_diff <- subset(subset_3treatments, subset_3treatments$clusterID_name == "4_MOD")[ , c("treatment_class", "hba1c_diff")] %>% 
  pairwise_t_test(hba1c_diff ~ treatment_class, p.adjust.method = "bonferroni", paired = TRUE, pool.sd = FALSE, detailed = TRUE)

Figure_MARD5_diff <- subset(subset_3treatments, subset_3treatments$clusterID_name == "5_MARD")[ , c("treatment_class", "hba1c_diff")] %>% 
  pairwise_t_test(hba1c_diff ~ treatment_class, p.adjust.method = "bonferroni", paired = TRUE, pool.sd = FALSE, detailed = TRUE)


#
## Part B ----------------------------------------------------------------------
# 


FigureB_data          <- data.frame(clusterID_name = c("2_SIDD", "2_SIDD", "3_SIRD", "3_SIRD", "4_MOD", "4_MOD", "5_MARD", "5_MARD"), 
                                      drugclass      = c("DPP4", "TZD", "DPP4", "TZD", "DPP4", "TZD", "DPP4", "TZD"), 
                                      mean           = as.numeric(c(Figure_SIDD2_diff$estimate[1:2],  Figure_SIRD3_diff$estimate[1:2],  Figure_MOD4_diff$estimate[1:2],  Figure_MARD5_diff$estimate[1:2])), 
                                      CI_low         = as.numeric(c(Figure_SIDD2_diff$conf.low[1:2],  Figure_SIRD3_diff$conf.low[1:2],  Figure_MOD4_diff$conf.low[1:2],  Figure_MARD5_diff$conf.low[1:2])), 
                                      CI_high        = as.numeric(c(Figure_SIDD2_diff$conf.high[1:2], Figure_SIRD3_diff$conf.high[1:2], Figure_MOD4_diff$conf.high[1:2], Figure_MARD5_diff$conf.high[1:2])))


FigureB_plot_data           <- FigureB_data
FigureB_plot_data$mean      <- -1*FigureB_data$mean
FigureB_plot_data$CI_high   <- -1*FigureB_data$CI_high
FigureB_plot_data$CI_low    <- -1*FigureB_data$CI_low


FigureB_cluster2 <- FigureB_plot_data[which(FigureB_plot_data$clusterID_name == "2_SIDD"), ]%>% 
  ggplot(aes(x = clusterID_name, y = mean, fill = drugclass))+
  geom_col( position = "dodge", width = 0.5, alpha = 0.5, color = "black", linewidth = 0.1)+
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high),
                position =  position_dodge(width = 0.5), width = 0.2, linewidth = 1)+
  scale_fill_manual(values = c("#005686", "#C33D13"), name = NULL) +
  geom_hline(yintercept = 0, color = "#DA9700", linewidth = 0.8) + 
  scale_y_continuous(limits = c(-5.5, 5.5)) + 
  xlab("") +
  ylab (" ") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        legend.position = "none", axis.text.x=element_blank(), axis.title=element_text(size = 15),
        axis.text = element_text(size = 15))

FigureB_cluster3 <- FigureB_plot_data[which(FigureB_plot_data$clusterID_name == "3_SIRD"), ]%>% 
  ggplot(aes(x = clusterID_name, y = mean, fill = drugclass))+
  geom_col( position = "dodge", width = 0.5, alpha = 0.5, color = "black", linewidth = 0.1)+
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high),
                position =  position_dodge(width = 0.5), width = 0.2, linewidth = 1)+
  scale_fill_manual(values = c("#005686", "#C33D13"), name = NULL) +
  geom_hline(yintercept = 0, color = "#DA9700", linewidth = 0.8) + 
  scale_y_continuous(limits = c(-5.5, 5.5)) + 
  xlab("") +
  ylab (" ") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        legend.position = "none", axis.text.x=element_blank(), axis.title=element_text(size = 15),
        axis.text = element_text(size = 15))

FigureB_cluster4 <- FigureB_plot_data[which(FigureB_plot_data$clusterID_name == "4_MOD"), ]%>% 
  ggplot(aes(x = clusterID_name, y = mean, fill = drugclass))+
  geom_col( position = "dodge", width = 0.5, alpha = 0.5, color = "black", linewidth = 0.1)+
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high),
                position =  position_dodge(width = 0.5), width = 0.2, linewidth = 1)+
  scale_fill_manual(values = c("#005686", "#C33D13"), name = NULL) +
  geom_hline(yintercept = 0, color = "#DA9700", linewidth = 0.8) + 
  scale_y_continuous(limits = c(-5.5, 5.5)) + 
  xlab("") +
  ylab (" ") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        legend.position = "none", axis.text.x=element_blank(), axis.title=element_text(size = 15),
        axis.text = element_text(size = 15))

FigureB_cluster5 <- FigureB_plot_data[which(FigureB_plot_data$clusterID_name == "5_MARD"), ]%>% 
  ggplot(aes(x = clusterID_name, y = mean, fill = drugclass))+
  geom_col( position = "dodge", width = 0.5, alpha = 0.5, color = "black", linewidth = 0.1)+
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high),
                position =  position_dodge(width = 0.5), width = 0.2, linewidth = 1)+
  scale_fill_manual(values = c("#005686", "#C33D13"), name = NULL) +
  geom_hline(yintercept = 0, color = "#DA9700", linewidth = 0.8) + 
  scale_y_continuous(limits = c(-5.5, 5.5)) + 
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


Cairo(file = paste0(plot_path,"/clinical_cluster_SAresults_unadjusted.png"), 
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
