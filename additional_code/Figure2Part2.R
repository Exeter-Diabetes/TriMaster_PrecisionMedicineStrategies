#
## Part A ----------------------------------------------------------------------
#


Figure_data_full <- subset_3treatments_wmodelC[, c("concordant_tx", "treatment_class", "hba1c_diff")]

FigureA_data     <- Figure_data_full %>% 
  group_by(concordant_tx, treatment_class) %>% 
  summarise(mean= mean(hba1c_diff),
            se = 1.96* sd(hba1c_diff)/sqrt(n()))


FigureA_DPP4 <- FigureA_data[which(FigureA_data$concordant_tx  == "DPP4"), ] %>%  
  ggplot(aes(x = concordant_tx, y = mean, fill = treatment_class))+
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


FigureA_SGLT2 <- FigureA_data[which(FigureA_data$concordant_tx  == "SGLT2"), ] %>%  
  ggplot(aes(x = concordant_tx, y = mean, fill = treatment_class))+
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

FigureA_TZD <- FigureA_data[which(FigureA_data$concordant_tx  == "TZD"), ] %>%  
  ggplot(aes(x = concordant_tx, y = mean, fill = treatment_class))+
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



#
## pairwise t-test -------------------------------------------------------------
#


Figure_DPP4_diff <- subset(subset_3treatments_wmodelC, subset_3treatments_wmodelC$concordant_tx == "DPP4")[ , c("treatment_class", "hba1c_diff")] %>% 
  pairwise_t_test(hba1c_diff ~ treatment_class, p.adjust.method = "bonferroni", paired = TRUE, pool.sd = FALSE, detailed = TRUE)

Figure_SGLT2_diff <- subset(subset_3treatments_wmodelC, subset_3treatments_wmodelC$concordant_tx == "SGLT2")[ , c("treatment_class", "hba1c_diff")] %>% 
  pairwise_t_test(hba1c_diff ~ treatment_class, p.adjust.method = "bonferroni", paired = TRUE, pool.sd = FALSE, detailed = TRUE)

Figure_TZD_diff <- subset(subset_3treatments_wmodelC, subset_3treatments_wmodelC$concordant_tx == "TZD")[ , c("treatment_class", "hba1c_diff")] %>% 
  pairwise_t_test(hba1c_diff ~ treatment_class, p.adjust.method = "bonferroni", paired = TRUE, pool.sd = FALSE, detailed = TRUE)



#
## Part B ----------------------------------------------------------------------
#


FigureB_data <- data.frame(clusterID_name = c("DPP4", "DPP4", "SGLT2", "SGLT2", "TZD", "TZD"), 
                             drugclass      = c("DPP4", "TZD", "DPP4", "TZD", "DPP4", "TZD"), 
                             mean           = as.numeric(c(Figure_DPP4_diff$estimate[1:2],  Figure_SGLT2_diff$estimate[1:2],  Figure_TZD_diff$estimate[1:2])), 
                             CI_low         = as.numeric(c(Figure_DPP4_diff$conf.low[1:2],  Figure_SGLT2_diff$conf.low[1:2],  Figure_TZD_diff$conf.low[1:2])), 
                             CI_high        = as.numeric(c(Figure_DPP4_diff$conf.high[1:2], Figure_SGLT2_diff$conf.high[1:2], Figure_TZD_diff$conf.high[1:2])))

FigureB_plotdata         <- FigureB_data
FigureB_plotdata$mean    <- -1*FigureB_data$mean
FigureB_plotdata$CI_low  <- -1*FigureB_data$CI_low
FigureB_plotdata$CI_high <- -1*FigureB_data$CI_high



FigureB_DPP4 <- FigureB_plotdata[which(FigureB_plotdata$clusterID_name == "DPP4"), ]%>% 
  ggplot(aes(x = clusterID_name, y = mean, fill = drugclass))+
  geom_col( position = "dodge", width = 0.5, alpha = 0.5, color = "black", linewidth = 0.1)+
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high),
                position =  position_dodge(width = 0.5), width = 0.2, linewidth = 1)+
  scale_fill_manual(values = c("#005686", "#C33D13"), name = NULL) +
  geom_hline(yintercept = 0, color = "#DA9700", linewidth = 0.8) + 
  scale_y_continuous(limits = c(-6.5, 5.6)) + 
  xlab("") +
  ylab (" ") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        legend.position = "none", axis.text.x=element_blank(), axis.title=element_text(size = 15),
        axis.text = element_text(size = 15))


FigureB_SGLT2 <- FigureB_plotdata[which(FigureB_plotdata$clusterID_name == "SGLT2"), ]%>% 
  ggplot(aes(x = clusterID_name, y = mean, fill = drugclass))+
  geom_col( position = "dodge", width = 0.5, alpha = 0.5, color = "black", linewidth = 0.1)+
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high),
                position =  position_dodge(width = 0.5), width = 0.2, linewidth = 1)+
  scale_fill_manual(values = c("#005686", "#C33D13"), name = NULL) +
  geom_hline(yintercept = 0, color = "#DA9700", linewidth = 0.8) + 
  scale_y_continuous(limits = c(-6.5, 5.6)) + 
  xlab("") +
  ylab (" ") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        legend.position = "none", axis.text.x=element_blank(), axis.title=element_text(size = 15),
        axis.text = element_text(size = 15))


FigureB_TZD <- FigureB_plotdata[which(FigureB_plotdata$clusterID_name == "TZD"), ]%>% 
  ggplot(aes(x = clusterID_name, y = mean, fill = drugclass))+
  geom_col( position = "dodge", width = 0.5, alpha = 0.5, color = "black", linewidth = 0.1)+
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high),
                position =  position_dodge(width = 0.5), width = 0.2, linewidth = 1)+
  scale_fill_manual(values = c("#005686", "#C33D13"), name = NULL) +
  geom_hline(yintercept = 0, color = "#DA9700", linewidth = 0.8) + 
  scale_y_continuous(limits = c(-6.5, 5.6)) + 
  xlab("") +
  ylab (" ") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        legend.position = "none", axis.text.x=element_blank(), axis.title=element_text(size = 15),
        axis.text = element_text(size = 15))



#
## Combination of all plots ----------------------------------------------------
#


Cairo(file = paste0(plot_path,"/model_SAresults_unadjusted_calibrated.png" ), 
      type = "png",
      units = "in", 
      width = 18,
      height = 8, 
      pointsize = 12, 
      dpi = 72)

plot_grid(FigureA_DPP4, FigureA_SGLT2, FigureA_TZD,
          FigureB_DPP4, FigureB_SGLT2, FigureB_TZD,
          ncol = 3)

dev.off()

