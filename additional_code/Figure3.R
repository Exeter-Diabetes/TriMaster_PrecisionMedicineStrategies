
Figure_cluster_data <- data_analysis_cluster %>% 
  group_by(group_cluster) %>% 
  summarise(mean= mean(hba1c_diff),
            se = 1.96* sd(hba1c_diff)/sqrt(n()))

Figure_cluster_data           <- data.frame(Figure_cluster_data)
colnames(Figure_cluster_data) <- c("groups", "mean", "se")


Figure_model_data <- data_analysis_model %>% 
  group_by(group_model) %>% 
  summarise(mean= mean(hba1c_diff),
            se = 1.96* sd(hba1c_diff)/sqrt(n()))

Figure_model_data           <- data.frame(Figure_model_data)
colnames(Figure_model_data) <- c("groups", "mean", "se")


Figure_diff_results         <- rbind(Figure_cluster_data, Figure_model_data)
Figure_diff_results$methods <- c("cluster", "cluster", "model", "model")


Figure <- Figure_diff_results %>% 
  ggplot(aes(x = methods, y = mean, fill = groups))+
  geom_col( position = "dodge", width = 0.5, alpha = 0.5, color = "black", linewidth = 0.1)+
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position =  position_dodge(width = 0.5), width = 0.2, linewidth = 1)+
  scale_fill_manual(values = c("#059212", "brown"), name = NULL) +
  scale_y_continuous(limits = c(-16, 0)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))



Cairo(file = paste0(plot_path,"/summary_results_cluster_model.png" ), 
      type = "png",
      units = "in", 
      width = 14,
      height = 8, 
      pointsize = 12, 
      dpi = 72)

plot(Figure)

dev.off()


#
## pairwise t-test -------------------------------------------------------------
#

# adjusted for multiple comparisons
pairwise.t.test(data_analysis_cluster$hba1c_diff, data_analysis_cluster$group_cluster, paired = TRUE)$p.value*2
pairwise.t.test(data_analysis_model$hba1c_diff, data_analysis_model$group_model, paired = TRUE)$p.value*2

