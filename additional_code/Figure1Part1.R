
drug_order <- list(c("SGLT2", "DPP4", "TZD"), c("DPP4", "SGLT2", "TZD"))

gcluster_associations                   <- data.frame(matrix(NA, nrow = 8*3, ncol = 6))
colnames(gcluster_associations)         <- c("genetic cluster","tx_comparison", "estimate", "CI_lower", "CI_higher", "pvalues")
gcluster_associations$tx_comparison     <- c(rep("SGLT2 vs DPP4", times = 8), rep("DPP4 vs TZD", times = 8), rep("SGLT2 vs TZD", times = 8))
gcluster_associations$`genetic cluster` <- c(rep(gcluster_names, times = 3))


## first two comparisons -------------------------------------------------------

subset_3treatments$treatment_class <- factor(subset_3treatments$treatment_class, levels = drug_order[[1]])      
  
gcluster_model1 <- lmer(hba1c_diff ~ treatment_class  * t2d_betacell_pineg_zscore    + period + (1 | study_id), data = subset_3treatments)
gcluster_model2 <- lmer(hba1c_diff ~ treatment_class  * t2d_betacell_pipos_zscore    + period + (1 | study_id), data = subset_3treatments)
gcluster_model3 <- lmer(hba1c_diff ~ treatment_class  * t2d_bodyfat_zscore           + period + (1 | study_id), data = subset_3treatments)
gcluster_model4 <- lmer(hba1c_diff ~ treatment_class  * t2d_lipodystrophy_zscore     + period + (1 | study_id), data = subset_3treatments)
gcluster_model5 <- lmer(hba1c_diff ~ treatment_class  * t2d_liverlipid_zscore        + period + (1 | study_id), data = subset_3treatments)
gcluster_model6 <- lmer(hba1c_diff ~ treatment_class  * t2d_metabollicsyn_zscore     + period + (1 | study_id), data = subset_3treatments)
gcluster_model7 <- lmer(hba1c_diff ~ treatment_class  * t2d_obesity_zscore           + period + (1 | study_id), data = subset_3treatments)
gcluster_model8 <- lmer(hba1c_diff ~ treatment_class  * t2d_residualglycaemic_zscore + period + (1 | study_id), data = subset_3treatments)

# Point estimates
gcluster_associations[gcluster_associations$`genetic cluster` == "t2d_betacell_pineg", ][c(1,3), "estimate"] <- 
summary(gcluster_model1)$coef[c("treatment_classDPP4:t2d_betacell_pineg_zscore", "treatment_classTZD:t2d_betacell_pineg_zscore"), "Estimate"]

gcluster_associations[gcluster_associations$`genetic cluster` == "t2d_betacell_pipos", ][c(1,3), "estimate"] <- 
summary(gcluster_model2)$coef[c("treatment_classDPP4:t2d_betacell_pipos_zscore", "treatment_classTZD:t2d_betacell_pipos_zscore"), "Estimate"]

gcluster_associations[gcluster_associations$`genetic cluster` == "t2d_bodyfat", ][c(1,3), "estimate"] <- 
summary(gcluster_model3)$coef[c("treatment_classDPP4:t2d_bodyfat_zscore", "treatment_classTZD:t2d_bodyfat_zscore"), "Estimate"]

gcluster_associations[gcluster_associations$`genetic cluster` == "t2d_lipodystrophy", ][c(1,3), "estimate"] <- 
summary(gcluster_model4)$coef[c("treatment_classDPP4:t2d_lipodystrophy_zscore", "treatment_classTZD:t2d_lipodystrophy_zscore"), "Estimate"]

gcluster_associations[gcluster_associations$`genetic cluster` == "t2d_liverlipid", ][c(1,3), "estimate"] <- 
summary(gcluster_model5)$coef[c("treatment_classDPP4:t2d_liverlipid_zscore", "treatment_classTZD:t2d_liverlipid_zscore"), "Estimate"]

gcluster_associations[gcluster_associations$`genetic cluster` == "t2d_metabollicsyn", ][c(1,3), "estimate"] <- 
summary(gcluster_model6)$coef[c("treatment_classDPP4:t2d_metabollicsyn_zscore", "treatment_classTZD:t2d_metabollicsyn_zscore"), "Estimate"]

gcluster_associations[gcluster_associations$`genetic cluster` == "t2d_obesity", ][c(1,3), "estimate"] <- 
summary(gcluster_model7)$coef[c("treatment_classDPP4:t2d_obesity_zscore", "treatment_classTZD:t2d_obesity_zscore"), "Estimate"]

gcluster_associations[gcluster_associations$`genetic cluster` == "t2d_residualglycaemic", ][c(1,3), "estimate"] <- 
summary(gcluster_model8)$coef[c("treatment_classDPP4:t2d_residualglycaemic_zscore", "treatment_classTZD:t2d_residualglycaemic_zscore"), "Estimate"]


# p-values
gcluster_associations[gcluster_associations$`genetic cluster` == "t2d_betacell_pineg", ][c(1,3), "pvalues"] <- 
  summary(gcluster_model1)$coef[c("treatment_classDPP4:t2d_betacell_pineg_zscore", "treatment_classTZD:t2d_betacell_pineg_zscore"), "Pr(>|t|)"]

gcluster_associations[gcluster_associations$`genetic cluster` == "t2d_betacell_pipos", ][c(1,3), "pvalues"] <- 
  summary(gcluster_model2)$coef[c("treatment_classDPP4:t2d_betacell_pipos_zscore", "treatment_classTZD:t2d_betacell_pipos_zscore"), "Pr(>|t|)"]

gcluster_associations[gcluster_associations$`genetic cluster` == "t2d_bodyfat", ][c(1,3), "pvalues"] <- 
  summary(gcluster_model3)$coef[c("treatment_classDPP4:t2d_bodyfat_zscore", "treatment_classTZD:t2d_bodyfat_zscore"), "Pr(>|t|)"]

gcluster_associations[gcluster_associations$`genetic cluster` == "t2d_lipodystrophy", ][c(1,3), "pvalues"] <- 
  summary(gcluster_model4)$coef[c("treatment_classDPP4:t2d_lipodystrophy_zscore", "treatment_classTZD:t2d_lipodystrophy_zscore"), "Pr(>|t|)"]

gcluster_associations[gcluster_associations$`genetic cluster` == "t2d_liverlipid", ][c(1,3), "pvalues"] <- 
  summary(gcluster_model5)$coef[c("treatment_classDPP4:t2d_liverlipid_zscore", "treatment_classTZD:t2d_liverlipid_zscore"), "Pr(>|t|)"]

gcluster_associations[gcluster_associations$`genetic cluster` == "t2d_metabollicsyn", ][c(1,3), "pvalues"] <- 
  summary(gcluster_model6)$coef[c("treatment_classDPP4:t2d_metabollicsyn_zscore", "treatment_classTZD:t2d_metabollicsyn_zscore"), "Pr(>|t|)"]

gcluster_associations[gcluster_associations$`genetic cluster` == "t2d_obesity", ][c(1,3), "pvalues"] <- 
  summary(gcluster_model7)$coef[c("treatment_classDPP4:t2d_obesity_zscore", "treatment_classTZD:t2d_obesity_zscore"), "Pr(>|t|)"]

gcluster_associations[gcluster_associations$`genetic cluster` == "t2d_residualglycaemic", ][c(1,3), "pvalues"] <- 
  summary(gcluster_model8)$coef[c("treatment_classDPP4:t2d_residualglycaemic_zscore", "treatment_classTZD:t2d_residualglycaemic_zscore"), "Pr(>|t|)"]


# Confidence intervals 
gcluster_associations[gcluster_associations$`genetic cluster` == "t2d_betacell_pineg", ][c(1,3), "CI_lower"] <- 
confint(gcluster_model1)[c("treatment_classDPP4:t2d_betacell_pineg_zscore", "treatment_classTZD:t2d_betacell_pineg_zscore"), ][ , "2.5 %"]
gcluster_associations[gcluster_associations$`genetic cluster` == "t2d_betacell_pineg", ][c(1,3), "CI_higher"] <- 
confint(gcluster_model1)[c("treatment_classDPP4:t2d_betacell_pineg_zscore", "treatment_classTZD:t2d_betacell_pineg_zscore"), ][ , "97.5 %"]

gcluster_associations[gcluster_associations$`genetic cluster` == "t2d_betacell_pipos", ][c(1,3), "CI_lower"] <-
confint(gcluster_model2)[c("treatment_classDPP4:t2d_betacell_pipos_zscore", "treatment_classTZD:t2d_betacell_pipos_zscore"), ][ , "2.5 %"]
gcluster_associations[gcluster_associations$`genetic cluster` == "t2d_betacell_pipos", ][c(1,3), "CI_higher"] <-
confint(gcluster_model2)[c("treatment_classDPP4:t2d_betacell_pipos_zscore", "treatment_classTZD:t2d_betacell_pipos_zscore"), ][ , "97.5 %"]

gcluster_associations[gcluster_associations$`genetic cluster` == "t2d_bodyfat", ][c(1,3), "CI_lower"] <-
confint(gcluster_model3)[c("treatment_classDPP4:t2d_bodyfat_zscore", "treatment_classTZD:t2d_bodyfat_zscore"), ][ , "2.5 %"]
gcluster_associations[gcluster_associations$`genetic cluster` == "t2d_bodyfat", ][c(1,3), "CI_higher"] <-
confint(gcluster_model3)[c("treatment_classDPP4:t2d_bodyfat_zscore", "treatment_classTZD:t2d_bodyfat_zscore"), ][ , "97.5 %"]

gcluster_associations[gcluster_associations$`genetic cluster` == "t2d_lipodystrophy", ][c(1,3), "CI_lower"] <-
confint(gcluster_model4)[c("treatment_classDPP4:t2d_lipodystrophy_zscore", "treatment_classTZD:t2d_lipodystrophy_zscore"), ][ , "2.5 %"]
gcluster_associations[gcluster_associations$`genetic cluster` == "t2d_lipodystrophy", ][c(1,3), "CI_higher"] <-
confint(gcluster_model4)[c("treatment_classDPP4:t2d_lipodystrophy_zscore", "treatment_classTZD:t2d_lipodystrophy_zscore"), ][ , "97.5 %"]

gcluster_associations[gcluster_associations$`genetic cluster` == "t2d_liverlipid", ][c(1,3), "CI_lower"] <-
confint(gcluster_model5)[c("treatment_classDPP4:t2d_liverlipid_zscore", "treatment_classTZD:t2d_liverlipid_zscore"), ][ , "2.5 %"]
gcluster_associations[gcluster_associations$`genetic cluster` == "t2d_liverlipid", ][c(1,3), "CI_higher"] <-
confint(gcluster_model5)[c("treatment_classDPP4:t2d_liverlipid_zscore", "treatment_classTZD:t2d_liverlipid_zscore"), ][ , "97.5 %"]

gcluster_associations[gcluster_associations$`genetic cluster` == "t2d_metabollicsyn", ][c(1,3), "CI_lower"] <-
confint(gcluster_model6)[c("treatment_classDPP4:t2d_metabollicsyn_zscore", "treatment_classTZD:t2d_metabollicsyn_zscore"), ][ , "2.5 %"]
gcluster_associations[gcluster_associations$`genetic cluster` == "t2d_metabollicsyn", ][c(1,3), "CI_higher"] <-
confint(gcluster_model6)[c("treatment_classDPP4:t2d_metabollicsyn_zscore", "treatment_classTZD:t2d_metabollicsyn_zscore"), ][ , "97.5 %"]

gcluster_associations[gcluster_associations$`genetic cluster` == "t2d_obesity", ][c(1,3), "CI_lower"] <-
confint(gcluster_model7)[c("treatment_classDPP4:t2d_obesity_zscore", "treatment_classTZD:t2d_obesity_zscore"), ][ , "2.5 %"]
gcluster_associations[gcluster_associations$`genetic cluster` == "t2d_obesity", ][c(1,3), "CI_higher"] <-
confint(gcluster_model7)[c("treatment_classDPP4:t2d_obesity_zscore", "treatment_classTZD:t2d_obesity_zscore"), ][ , "97.5 %"]

gcluster_associations[gcluster_associations$`genetic cluster` == "t2d_residualglycaemic", ][c(1,3), "CI_lower"] <-
confint(gcluster_model8)[c("treatment_classDPP4:t2d_residualglycaemic_zscore", "treatment_classTZD:t2d_residualglycaemic_zscore"), ][ , "2.5 %"]
gcluster_associations[gcluster_associations$`genetic cluster` == "t2d_residualglycaemic", ][c(1,3), "CI_higher"] <-
confint(gcluster_model8)[c("treatment_classDPP4:t2d_residualglycaemic_zscore", "treatment_classTZD:t2d_residualglycaemic_zscore"), ][ , "97.5 %"]



## third comparison ------------------------------------------------------------

subset_3treatments$treatment_class <- factor(subset_3treatments$treatment_class, levels = drug_order[[2]])      

gcluster_model1 <- lmer(hba1c_diff ~ treatment_class  * t2d_betacell_pineg_zscore    + period + (1 | study_id), data = subset_3treatments)
gcluster_model2 <- lmer(hba1c_diff ~ treatment_class  * t2d_betacell_pipos_zscore    + period + (1 | study_id), data = subset_3treatments)
gcluster_model3 <- lmer(hba1c_diff ~ treatment_class  * t2d_bodyfat_zscore           + period + (1 | study_id), data = subset_3treatments)
gcluster_model4 <- lmer(hba1c_diff ~ treatment_class  * t2d_lipodystrophy_zscore     + period + (1 | study_id), data = subset_3treatments)
gcluster_model5 <- lmer(hba1c_diff ~ treatment_class  * t2d_liverlipid_zscore        + period + (1 | study_id), data = subset_3treatments)
gcluster_model6 <- lmer(hba1c_diff ~ treatment_class  * t2d_metabollicsyn_zscore     + period + (1 | study_id), data = subset_3treatments)
gcluster_model7 <- lmer(hba1c_diff ~ treatment_class  * t2d_obesity_zscore           + period + (1 | study_id), data = subset_3treatments)
gcluster_model8 <- lmer(hba1c_diff ~ treatment_class  * t2d_residualglycaemic_zscore + period + (1 | study_id), data = subset_3treatments)


# Point estimates
gcluster_associations[gcluster_associations$`genetic cluster` == "t2d_betacell_pineg", ][c(2), "estimate"] <- 
  summary(gcluster_model1)$coef["treatment_classTZD:t2d_betacell_pineg_zscore", "Estimate"]

gcluster_associations[gcluster_associations$`genetic cluster` == "t2d_betacell_pipos", ][c(2), "estimate"] <- 
  summary(gcluster_model2)$coef["treatment_classTZD:t2d_betacell_pipos_zscore", "Estimate"]

gcluster_associations[gcluster_associations$`genetic cluster` == "t2d_bodyfat", ][c(2), "estimate"] <- 
  summary(gcluster_model3)$coef["treatment_classTZD:t2d_bodyfat_zscore", "Estimate"]

gcluster_associations[gcluster_associations$`genetic cluster` == "t2d_lipodystrophy", ][c(2), "estimate"] <- 
  summary(gcluster_model4)$coef["treatment_classTZD:t2d_lipodystrophy_zscore", "Estimate"]

gcluster_associations[gcluster_associations$`genetic cluster` == "t2d_liverlipid", ][c(2), "estimate"] <- 
  summary(gcluster_model5)$coef["treatment_classTZD:t2d_liverlipid_zscore", "Estimate"]

gcluster_associations[gcluster_associations$`genetic cluster` == "t2d_metabollicsyn", ][c(2), "estimate"] <- 
  summary(gcluster_model6)$coef["treatment_classTZD:t2d_metabollicsyn_zscore", "Estimate"]

gcluster_associations[gcluster_associations$`genetic cluster` == "t2d_obesity", ][c(2), "estimate"] <- 
  summary(gcluster_model7)$coef["treatment_classTZD:t2d_obesity_zscore", "Estimate"]

gcluster_associations[gcluster_associations$`genetic cluster` == "t2d_residualglycaemic", ][c(2), "estimate"] <- 
  summary(gcluster_model8)$coef["treatment_classTZD:t2d_residualglycaemic_zscore", "Estimate"]

# p-values 
gcluster_associations[gcluster_associations$`genetic cluster` == "t2d_betacell_pineg", ][c(2), "pvalues"] <- 
  summary(gcluster_model1)$coef["treatment_classTZD:t2d_betacell_pineg_zscore", "Pr(>|t|)"]

gcluster_associations[gcluster_associations$`genetic cluster` == "t2d_betacell_pipos", ][c(2), "pvalues"] <- 
  summary(gcluster_model2)$coef["treatment_classTZD:t2d_betacell_pipos_zscore", "Pr(>|t|)"]

gcluster_associations[gcluster_associations$`genetic cluster` == "t2d_bodyfat", ][c(2), "pvalues"] <- 
  summary(gcluster_model3)$coef["treatment_classTZD:t2d_bodyfat_zscore", "Pr(>|t|)"]

gcluster_associations[gcluster_associations$`genetic cluster` == "t2d_lipodystrophy", ][c(2), "pvalues"] <- 
  summary(gcluster_model4)$coef["treatment_classTZD:t2d_lipodystrophy_zscore", "Pr(>|t|)"]

gcluster_associations[gcluster_associations$`genetic cluster` == "t2d_liverlipid", ][c(2), "pvalues"] <- 
  summary(gcluster_model5)$coef["treatment_classTZD:t2d_liverlipid_zscore", "Pr(>|t|)"]

gcluster_associations[gcluster_associations$`genetic cluster` == "t2d_metabollicsyn", ][c(2), "pvalues"] <- 
  summary(gcluster_model6)$coef["treatment_classTZD:t2d_metabollicsyn_zscore", "Pr(>|t|)"]

gcluster_associations[gcluster_associations$`genetic cluster` == "t2d_obesity", ][c(2), "pvalues"] <- 
  summary(gcluster_model7)$coef["treatment_classTZD:t2d_obesity_zscore", "Pr(>|t|)"]

gcluster_associations[gcluster_associations$`genetic cluster` == "t2d_residualglycaemic", ][c(2), "pvalues"] <- 
  summary(gcluster_model8)$coef["treatment_classTZD:t2d_residualglycaemic_zscore", "Pr(>|t|)"]


# confidence intervals 

gcluster_associations[gcluster_associations$`genetic cluster` == "t2d_betacell_pineg", ][c(2), "CI_lower"] <- 
  confint(gcluster_model1)["treatment_classTZD:t2d_betacell_pineg_zscore", ]["2.5 %"]
gcluster_associations[gcluster_associations$`genetic cluster` == "t2d_betacell_pineg", ][c(2), "CI_higher"] <- 
  confint(gcluster_model1)["treatment_classTZD:t2d_betacell_pineg_zscore", ]["97.5 %"]

gcluster_associations[gcluster_associations$`genetic cluster` == "t2d_betacell_pipos", ][c(2), "CI_lower"] <- 
  confint(gcluster_model2)["treatment_classTZD:t2d_betacell_pipos_zscore", ]["2.5 %"]
gcluster_associations[gcluster_associations$`genetic cluster` == "t2d_betacell_pipos", ][c(2), "CI_higher"] <- 
  confint(gcluster_model2)["treatment_classTZD:t2d_betacell_pipos_zscore", ]["97.5 %"]

gcluster_associations[gcluster_associations$`genetic cluster` == "t2d_bodyfat", ][c(2), "CI_lower"] <- 
  confint(gcluster_model3)["treatment_classTZD:t2d_bodyfat_zscore", ]["2.5 %"]
gcluster_associations[gcluster_associations$`genetic cluster` == "t2d_bodyfat", ][c(2), "CI_higher"] <- 
  confint(gcluster_model3)["treatment_classTZD:t2d_bodyfat_zscore", ]["97.5 %"]

gcluster_associations[gcluster_associations$`genetic cluster` == "t2d_lipodystrophy", ][c(2), "CI_lower"] <- 
  confint(gcluster_model4)["treatment_classTZD:t2d_lipodystrophy_zscore", ]["2.5 %"]
gcluster_associations[gcluster_associations$`genetic cluster` == "t2d_lipodystrophy", ][c(2), "CI_higher"] <- 
  confint(gcluster_model4)["treatment_classTZD:t2d_lipodystrophy_zscore", ]["97.5 %"]

gcluster_associations[gcluster_associations$`genetic cluster` == "t2d_liverlipid", ][c(2), "CI_lower"] <- 
  confint(gcluster_model5)["treatment_classTZD:t2d_liverlipid_zscore", ]["2.5 %"]
gcluster_associations[gcluster_associations$`genetic cluster` == "t2d_liverlipid", ][c(2), "CI_higher"] <- 
  confint(gcluster_model5)["treatment_classTZD:t2d_liverlipid_zscore", ]["97.5 %"]

gcluster_associations[gcluster_associations$`genetic cluster` == "t2d_metabollicsyn", ][c(2), "CI_lower"] <- 
  confint(gcluster_model6)["treatment_classTZD:t2d_metabollicsyn_zscore", ]["2.5 %"]
gcluster_associations[gcluster_associations$`genetic cluster` == "t2d_metabollicsyn", ][c(2), "CI_higher"] <- 
  confint(gcluster_model6)["treatment_classTZD:t2d_metabollicsyn_zscore", ]["97.5 %"]

gcluster_associations[gcluster_associations$`genetic cluster` == "t2d_obesity", ][c(2), "CI_lower"] <- 
  confint(gcluster_model7)["treatment_classTZD:t2d_obesity_zscore", ]["2.5 %"]
gcluster_associations[gcluster_associations$`genetic cluster` == "t2d_obesity", ][c(2), "CI_higher"] <- 
  confint(gcluster_model7)["treatment_classTZD:t2d_obesity_zscore", ]["97.5 %"]

gcluster_associations[gcluster_associations$`genetic cluster` == "t2d_residualglycaemic", ][c(2), "CI_lower"] <- 
  confint(gcluster_model8)["treatment_classTZD:t2d_residualglycaemic_zscore", ]["2.5 %"]
gcluster_associations[gcluster_associations$`genetic cluster` == "t2d_residualglycaemic", ][c(2), "CI_higher"] <- 
  confint(gcluster_model8)["treatment_classTZD:t2d_residualglycaemic_zscore", ]["97.5 %"]


## revert drug order back to normal --- --- --- --- --- --- --- --- --- --- --- 


subset_3treatments$treatment_class <- factor(subset_3treatments$treatment_class, levels = drug_order[[1]])      


# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --
# Make plot for genetic cluster ------------------------------------------------
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --


genetic_cluster_colors <- c("orchid3", "chartreuse4", "goldenrod")
gcluster_associations$`genetic cluster` <- factor(gcluster_associations$`genetic cluster`, 
                                                  levels = rev(gcluster_names))


gc_plot_SGLT2vsDPP4 <- 
ggplot(gcluster_associations[gcluster_associations$tx_comparison == "SGLT2 vs DPP4", ], aes(`genetic cluster`, estimate)) +
  geom_col(aes(fill = tx_comparison), position = position_dodge(0.6), 
           color = "black", alpha = 0.7, width = 0.6) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_higher), 
                position = position_dodge(0.6), width = 0.15, size = 0.8) +
  scale_fill_manual(values = genetic_cluster_colors[1]) +
  theme_minimal(base_size = 20) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),  
    panel.grid.minor = element_blank(),     
    axis.text.y = element_blank(), legend.position = "none" # element_text(size = 22),
  ) + 
  labs(x = NULL, y = " ") +
  #coord_cartesian(ylim = c(min(gcluster_associations[, c(3,4,5)]), max(gcluster_associations[, c(3,4,5)]))) +
  geom_vline(xintercept = (0:8)+0.5) +
  geom_hline(yintercept = 0, color = "black", size = 1.2) +
  ylim(-2.1, 2.1) + 
  geom_vline(xintercept = (0:8)+0.5, color = "gray") +
  coord_flip()


gc_plot_DPP4vsTZD <- 
  ggplot(gcluster_associations[gcluster_associations$tx_comparison == "DPP4 vs TZD", ], aes(`genetic cluster`, estimate)) +
  geom_col(aes(fill = tx_comparison), position = position_dodge(0.6), 
           color = "black", alpha = 0.7, width = 0.6) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_higher), 
                position = position_dodge(0.6), width = 0.15, size = 0.8) +
  scale_fill_manual(values = genetic_cluster_colors[2]) +
  theme_minimal(base_size = 20) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),  
    panel.grid.minor = element_blank(),     
    axis.text.y = element_blank(), legend.position = "none"
  ) + 
  labs(x = NULL, y = " ") +
  #coord_cartesian(ylim = c(min(gcluster_associations[, c(3,4,5)]), max(gcluster_associations[, c(3,4,5)]))) +
  geom_vline(xintercept = (0:8)+0.5) +
  geom_hline(yintercept = 0, color = "black", size = 1.2) +
  ylim(-2.1, 2.1) + 
  geom_vline(xintercept = (0:8)+0.5, color = "gray") +
  coord_flip()



gc_plot_SGLT2vsTZD <- 
  ggplot(gcluster_associations[gcluster_associations$tx_comparison == "SGLT2 vs TZD", ], aes(`genetic cluster`, estimate)) +
  geom_col(aes(fill = tx_comparison), position = position_dodge(0.6), 
           color = "black", alpha = 0.7, width = 0.6) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_higher), 
                position = position_dodge(0.6), width = 0.15, size = 0.8) +
  scale_fill_manual(values = genetic_cluster_colors[3]) +
  theme_minimal(base_size = 20) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),  
    panel.grid.minor = element_blank(),     
    axis.text.y = element_blank(), legend.position = "none"
  ) + 
  labs(x = NULL, y = " ") +
  #coord_cartesian(ylim = c(min(gcluster_associations[, c(3,4,5)]), max(gcluster_associations[, c(3,4,5)]))) +
  geom_vline(xintercept = (0:8)+0.5) +
  geom_hline(yintercept = 0, color = "black", size = 1.2) +
  ylim(-2.1, 2.1) + 
  geom_vline(xintercept = (0:8)+0.5, color = "gray") +
  coord_flip()



Cairo(file = paste0(plot_path,"/genetic_cluster_association_results.jpeg" ), 
      type = "jpeg",
      units = "in", 
      width = 18,
      height = 8, 
      pointsize = 12, 
      dpi = 300)

plot_grid(gc_plot_SGLT2vsDPP4, gc_plot_DPP4vsTZD, gc_plot_SGLT2vsTZD, nrow = 1)

dev.off()




