#
## Part A ----------------------------------------------------------------------
#

model_FigureA <- lmer(hba1c_diff ~ treatment_class * concordant_tx + prehba1c + period + (1 | study_id), data = subset_3treatments_wmodelC)

emm_FigureA   <- emmeans(model_FigureA, ~ treatment_class | concordant_tx, cov.reduce = mean)
data_FigureA  <- as.data.frame(emm_FigureA[1:9,1])


FigureA_DPP4 <- data_FigureA[which(data_FigureA$concordant_tx == "DPP4"), ] %>% 
  ggplot(aes(x = concordant_tx, y = emmean, fill = treatment_class))+
  geom_col( position = "dodge", width = 0.5, alpha = 0.5, color = "black", linewidth = 0.1)+
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                position =  position_dodge(width = 0.5), width = 0.2, linewidth = 1)+
  scale_fill_manual(values = c("#DA9700", "#005686", "#C33D13"), name = NULL) +
  scale_y_continuous(limits = c(-17, 0)) + 
  xlab("") +
  ylab ("") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        legend.position = "none", axis.text.x=element_blank(), axis.title=element_text(size = 15),
        axis.text = element_text(size = 15))


FigureA_SGLT2 <- data_FigureA[which(data_FigureA$concordant_tx == "SGLT2"), ] %>% 
  ggplot(aes(x = concordant_tx, y = emmean, fill = treatment_class))+
  geom_col( position = "dodge", width = 0.5, alpha = 0.5, color = "black", linewidth = 0.1)+
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                position =  position_dodge(width = 0.5), width = 0.2, linewidth = 1)+
  scale_fill_manual(values = c("#DA9700", "#005686", "#C33D13"), name = NULL) +
  scale_y_continuous(limits = c(-17, 0)) + 
  xlab("") +
  ylab ("") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        legend.position = "none", axis.text.x=element_blank(), axis.title=element_text(size = 15),
        axis.text = element_text(size = 15))


FigureA_TZD <- data_FigureA[which(data_FigureA$concordant_tx == "TZD"), ] %>% 
  ggplot(aes(x = concordant_tx, y = emmean, fill = treatment_class))+
  geom_col( position = "dodge", width = 0.5, alpha = 0.5, color = "black", linewidth = 0.1)+
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                position =  position_dodge(width = 0.5), width = 0.2, linewidth = 1)+
  scale_fill_manual(values = c("#DA9700", "#005686", "#C33D13"), name = NULL) +
  scale_y_continuous(limits = c(-17, 0)) + 
  xlab("") +
  ylab ("") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        legend.position = "none", axis.text.x=element_blank(), axis.title=element_text(size = 15),
        axis.text = element_text(size = 15))


plot_grid(FigureA_DPP4, FigureA_SGLT2, FigureA_TZD, ncol = 3)


#
## pairwise t-test -------------------------------------------------------------
#

pairwise_comparison_Figure <- contrast(emm_FigureA, method = "pairwise", by = "concordant_tx", adjust = "bonferroni")

#
# Part B -----------------------------------------------------------------------
#

FigureB_data           <- data.frame(pairwise_comparison_Figure)[c(1, 2, 4, 5, 7, 8), ]
FigureB_data$drugclass <- c("DPP4", "TZD", "DPP4", "TZD", "DPP4", "TZD") 
FigureB_data$CI_low    <- FigureB_data$estimate - 1.96*FigureB_data$SE
FigureB_data$CI_high   <- FigureB_data$estimate + 1.96*FigureB_data$SE

FigureB_plotdata            <- FigureB_data
FigureB_plotdata$estimate   <- -1*FigureB_data$estimate
FigureB_plotdata$CI_low     <- -1*FigureB_data$CI_low
FigureB_plotdata$CI_high    <- -1*FigureB_data$CI_high 


FigureB_DPP4 <- FigureB_plotdata[which(FigureB_plotdata$concordant_tx == "DPP4"), ]%>% 
  ggplot(aes(x = concordant_tx, y = estimate, fill = drugclass))+
  geom_col( position = "dodge", width = 0.5, alpha = 0.5, color = "black", linewidth = 0.1)+
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high),
                position =  position_dodge(width = 0.5), width = 0.2, linewidth = 1)+
  scale_fill_manual(values = c("#005686", "#C33D13"), name = NULL) +
  geom_hline(yintercept = 0, color = "#DA9700", linewidth = 0.8) + 
  scale_y_continuous(limits = c(-5.7, 6.5)) + 
  xlab("") +
  ylab (" ") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        legend.position = "none", axis.text.x=element_blank(), axis.title=element_text(size = 15),
        axis.text = element_text(size = 15))


FigureB_SGLT2 <- FigureB_plotdata[which(FigureB_plotdata$concordant_tx == "SGLT2"), ]%>% 
  ggplot(aes(x = concordant_tx, y = estimate, fill = drugclass))+
  geom_col( position = "dodge", width = 0.5, alpha = 0.5, color = "black", linewidth = 0.1)+
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high),
                position =  position_dodge(width = 0.5), width = 0.2, linewidth = 1)+
  scale_fill_manual(values = c("#005686", "#C33D13"), name = NULL) +
  geom_hline(yintercept = 0, color = "#DA9700", linewidth = 0.8) + 
  scale_y_continuous(limits = c(-5.7, 6.5)) + 
  xlab("") +
  ylab (" ") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        legend.position = "none", axis.text.x=element_blank(), axis.title=element_text(size = 15),
        axis.text = element_text(size = 15))


FigureB_TZD <- FigureB_plotdata[which(FigureB_plotdata$concordant_tx == "TZD"), ]%>% 
  ggplot(aes(x = concordant_tx, y = estimate, fill = drugclass))+
  geom_col( position = "dodge", width = 0.5, alpha = 0.5, color = "black", linewidth = 0.1)+
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high),
                position =  position_dodge(width = 0.5), width = 0.2, linewidth = 1)+
  scale_fill_manual(values = c("#005686", "#C33D13"), name = NULL) +
  geom_hline(yintercept = 0, color = "#DA9700", linewidth = 0.8) + 
  scale_y_continuous(limits = c(-5.7, 6.5)) + 
  xlab("") +
  ylab (" ") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        legend.position = "none", axis.text.x=element_blank(), axis.title=element_text(size = 15),
        axis.text = element_text(size = 15))


plot_grid(FigureB_DPP4, FigureB_SGLT2, FigureB_TZD, ncol = 3)


#
## Combination of all plots ----------------------------------------------------
#

Cairo(file = paste0(plot_path,"/main_results_adjusted_model.png" ), 
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

