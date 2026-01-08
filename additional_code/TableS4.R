#
## Cluster Strategy ------------------------------------------------------------
#


cluster_c_tx           <- subset(data_analysis_cluster, data_analysis_cluster$group_cluster == "C")[, c("study_id", "drugclass")]
colnames(cluster_c_tx) <- c("study_id", "concordant_tx")

subset_3treatments_merged_cluster             <- merge(subset_3treatments, cluster_c_tx, by = "study_id", all.x = T )
subset_3treatments_merged_cluster$group_cluster[is.na(subset_3treatments_merged_cluster$group_cluster)] <- "D"

table(subset_3treatments_merged_cluster$concordant_tx)/3


# Concordant: DPP4

data_cluster_concordant_DPP4 <- subset(subset_3treatments_merged_cluster, subset_3treatments_merged_cluster$concordant_tx =="DPP4")

# DPP4i (best) vs. SGLT2i
coef(summary(lmer(posthba1cfinal ~ treatment_class + period + (1|study_id), data = data_cluster_concordant_DPP4[-which(data_cluster_concordant_DPP4$treatment_class == "TZD"), ])))
confint(lmer(posthba1cfinal ~ treatment_class + period + (1|study_id), data = data_cluster_concordant_DPP4[-which(data_cluster_concordant_DPP4$treatment_class == "TZD"), ]))

# DPP4i (best) vs TZD
coef(summary(lmer(posthba1cfinal ~ treatment_class + period + (1|study_id), data = data_cluster_concordant_DPP4[-which(data_cluster_concordant_DPP4$treatment_class == "SGLT2"), ])))
confint(lmer(posthba1cfinal ~ treatment_class + period + (1|study_id), data = data_cluster_concordant_DPP4[-which(data_cluster_concordant_DPP4$treatment_class == "SGLT2"), ]))

# DPP4i (best) vs other pooled
coef(summary(lmer(posthba1cfinal ~ group_cluster + period + (1|study_id), data = data_cluster_concordant_DPP4)))
confint(lmer(posthba1cfinal ~ group_cluster + period + (1|study_id), data = data_cluster_concordant_DPP4))

# Concordant: SGLT2

data_cluster_concordant_SGLT2                 <- subset(subset_3treatments_merged_cluster, subset_3treatments_merged_cluster$concordant_tx =="SGLT2")
data_cluster_concordant_SGLT2$treatment_class <- factor(data_cluster_concordant_SGLT2$treatment_class, levels = c("SGLT2", "DPP4", "TZD"))

# SGLT2i (best) vs. DPP4i
coef(summary(lmer(posthba1cfinal ~ treatment_class + period + (1|study_id), data = data_cluster_concordant_SGLT2[-which(data_cluster_concordant_SGLT2$treatment_class == "TZD"), ])))
confint(lmer(posthba1cfinal ~ treatment_class + period + (1|study_id), data = data_cluster_concordant_SGLT2[-which(data_cluster_concordant_SGLT2$treatment_class == "TZD"), ]))

# SGLT2i (best) vs. TZD
coef(summary(lmer(posthba1cfinal ~ treatment_class + period + (1|study_id), data = data_cluster_concordant_SGLT2[-which(data_cluster_concordant_SGLT2$treatment_class == "DPP4"), ])))
confint(lmer(posthba1cfinal ~ treatment_class + period + (1|study_id), data = data_cluster_concordant_SGLT2[-which(data_cluster_concordant_SGLT2$treatment_class == "DPP4"), ]))

# SGLT2i (best) vs other pooled
coef(summary(lmer(posthba1cfinal ~ group_cluster + period + (1|study_id), data = data_cluster_concordant_SGLT2)))
confint(lmer(posthba1cfinal ~ group_cluster + period + (1|study_id), data = data_cluster_concordant_SGLT2))


# Concordant: TZD

data_cluster_concordant_TZD                 <- subset(subset_3treatments_merged_cluster, subset_3treatments_merged_cluster$concordant_tx =="TZD")
data_cluster_concordant_TZD$treatment_class <- factor(data_cluster_concordant_TZD$treatment_class, levels = c("TZD", "DPP4", "SGLT2"))

# TZD (best) vs. DPP4i
coef(summary(lmer(posthba1cfinal ~ treatment_class + period + (1|study_id), data = data_cluster_concordant_TZD[-which(data_cluster_concordant_TZD$treatment_class == "SGLT2"), ])))
confint(lmer(posthba1cfinal ~ treatment_class + period + (1|study_id), data = data_cluster_concordant_TZD[-which(data_cluster_concordant_TZD$treatment_class == "SGLT2"), ]))

# TZD (best) vs. SGLT2i
coef(summary(lmer(posthba1cfinal ~ treatment_class + period + (1|study_id), data = data_cluster_concordant_TZD[-which(data_cluster_concordant_TZD$treatment_class == "DPP4"), ])))
confint(lmer(posthba1cfinal ~ treatment_class + period + (1|study_id), data = data_cluster_concordant_TZD[-which(data_cluster_concordant_TZD$treatment_class == "DPP4"), ]))

# TZD (best) vs. other pooled
coef(summary(lmer(posthba1cfinal ~ group_cluster + period + (1|study_id), data = data_cluster_concordant_TZD)))
confint(lmer(posthba1cfinal ~ group_cluster + period + (1|study_id), data = data_cluster_concordant_TZD))

# concordant: overall 
coef(summary(lmer(posthba1cfinal ~ group_cluster + period + (1|study_id), data = data_analysis_cluster)))
confint(lmer(posthba1cfinal ~ group_cluster + period + (1|study_id), data = data_analysis_cluster))


#
## Model Strategy (calibrated) -------------------------------------------------
#


model_c_tx           <- subset(data_analysis_model, data_analysis_model$group_model == "C")[, c("study_id", "drugclass")]
colnames(model_c_tx) <- c("study_id", "concordant_tx")

subset_3treatments_merged_model             <- merge(subset_3treatments, model_c_tx, by = "study_id", all.x = T )
subset_3treatments_merged_model$group_model[is.na(subset_3treatments_merged_model$group_model)] <- "D"

table(subset_3treatments_merged_model$concordant_tx)/3


# Concordant: DPP4
data_model_concordant_DPP4 <- subset(subset_3treatments_merged_model, subset_3treatments_merged_model$concordant_tx =="DPP4")

# DPP4i (best) vs. SGLT2i
coef(summary(lmer(posthba1cfinal ~ treatment_class + period + (1|study_id), data = data_model_concordant_DPP4[-which(data_model_concordant_DPP4$treatment_class == "TZD"), ])))
confint(lmer(posthba1cfinal ~ treatment_class + period + (1|study_id), data = data_model_concordant_DPP4[-which(data_model_concordant_DPP4$treatment_class == "TZD"), ]))

# DPP4i (best) vs. TZD
coef(summary(lmer(posthba1cfinal ~ treatment_class + period + (1|study_id), data = data_model_concordant_DPP4[-which(data_model_concordant_DPP4$treatment_class == "SGLT2"), ])))
confint(lmer(posthba1cfinal ~ treatment_class + period + (1|study_id), data = data_model_concordant_DPP4[-which(data_model_concordant_DPP4$treatment_class == "SGLT2"), ]))

# DPP4i (best) vs. pooled other
coef(summary(lmer(posthba1cfinal ~ group_model + period + (1|study_id), data = data_model_concordant_DPP4)))
confint(lmer(posthba1cfinal ~ group_model + period + (1|study_id), data = data_model_concordant_DPP4))


# Concordant: SGLT2
data_model_concordant_SGLT2                 <- subset(subset_3treatments_merged_model, subset_3treatments_merged_model$concordant_tx =="SGLT2")
data_model_concordant_SGLT2$treatment_class <- factor(data_model_concordant_SGLT2$treatment_class, levels = c("SGLT2", "DPP4", "TZD"))

# SGLT2i (best) vs. DPP4i
coef(summary(lmer(posthba1cfinal ~ treatment_class + period + (1|study_id), data = data_model_concordant_SGLT2[-which(data_model_concordant_SGLT2$treatment_class == "TZD"), ])))
confint(lmer(posthba1cfinal ~ treatment_class + period + (1|study_id), data = data_model_concordant_SGLT2[-which(data_model_concordant_SGLT2$treatment_class == "TZD"), ]))

# SGLT2i (best) vs. TZD
coef(summary(lmer(posthba1cfinal ~ treatment_class + period + (1|study_id), data = data_model_concordant_SGLT2[-which(data_model_concordant_SGLT2$treatment_class == "DPP4"), ])))
confint(lmer(posthba1cfinal ~ treatment_class + period + (1|study_id), data = data_model_concordant_SGLT2[-which(data_model_concordant_SGLT2$treatment_class == "DPP4"), ]))

# SGLT2i (best) vs. pooled other
coef(summary(lmer(posthba1cfinal ~ group_model + period + (1|study_id), data = data_model_concordant_SGLT2)))
confint(lmer(posthba1cfinal ~ group_model + period + (1|study_id), data = data_model_concordant_SGLT2))


# Concordant: TZD
data_model_concordant_TZD                 <- subset(subset_3treatments_merged_model, subset_3treatments_merged_model$concordant_tx =="TZD")
data_model_concordant_TZD$treatment_class <- factor(data_model_concordant_TZD$treatment_class, levels = c("TZD", "DPP4", "SGLT2"))

# TZD (best) vs. DPP4i
coef(summary(lmer(posthba1cfinal ~ treatment_class + period + (1|study_id), data = data_model_concordant_TZD[-which(data_model_concordant_TZD$treatment_class == "SGLT2"), ])))
confint(lmer(posthba1cfinal ~ treatment_class + period + (1|study_id), data = data_model_concordant_TZD[-which(data_model_concordant_TZD$treatment_class == "SGLT2"), ]))

# TZD (best) vs. SGLT2i
coef(summary(lmer(posthba1cfinal ~ treatment_class + period + (1|study_id), data = data_model_concordant_TZD[-which(data_model_concordant_TZD$treatment_class == "DPP4"), ])))
confint(lmer(posthba1cfinal ~ treatment_class + period + (1|study_id), data = data_model_concordant_TZD[-which(data_model_concordant_TZD$treatment_class == "DPP4"), ]))

# TZD (best) vs. pooled other
coef(summary(lmer(posthba1cfinal ~ group_model + period + (1|study_id), data = data_model_concordant_TZD)))
confint(lmer(posthba1cfinal ~ group_model + period + (1|study_id), data = data_model_concordant_TZD))

# concordant: overall 
coef(summary(lmer(posthba1cfinal ~ group_model + period + (1|study_id), data = data_analysis_model)))
confint(lmer(posthba1cfinal ~ group_model + period + (1|study_id), data = data_analysis_model))


#
## Model Strategy (uncalibrated) ---------------------------------------------------------------
#

modeluncal_c_tx           <- subset(data_analysis_modeluncal, data_analysis_modeluncal$group_modeluncal == "C")[, c("study_id", "drugclass")]
colnames(modeluncal_c_tx) <- c("study_id", "concordant_tx")

subset_3treatments_merged_modeluncal             <- merge(subset_3treatments, modeluncal_c_tx, by = "study_id", all.x = T )
subset_3treatments_merged_modeluncal$group_modeluncal_cd[is.na(subset_3treatments_merged_modeluncal$group_modeluncal)] <- "D"

table(subset_3treatments_merged_modeluncal$concordant_tx)/3


# Concordant: DPP4
data_modeluncal_concordant_DPP4 <- subset(subset_3treatments_merged_modeluncal, subset_3treatments_merged_modeluncal$concordant_tx =="DPP4")

# DPP4i (best) vs. SGLT2i
coef(summary(lmer(posthba1cfinal ~ treatment_class + period + (1|study_id), data = data_modeluncal_concordant_DPP4[-which(data_modeluncal_concordant_DPP4$treatment_class == "TZD"), ])))
confint(lmer(posthba1cfinal ~ treatment_class + period + (1|study_id), data = data_modeluncal_concordant_DPP4[-which(data_modeluncal_concordant_DPP4$treatment_class == "TZD"), ]))

# DPP4i (best) vs. TZD
coef(summary(lmer(posthba1cfinal ~ treatment_class + period + (1|study_id), data = data_modeluncal_concordant_DPP4[-which(data_modeluncal_concordant_DPP4$treatment_class == "SGLT2"), ])))
confint(lmer(posthba1cfinal ~ treatment_class + period + (1|study_id), data = data_modeluncal_concordant_DPP4[-which(data_modeluncal_concordant_DPP4$treatment_class == "SGLT2"), ]))

# DPP4i (best) vs. pooled other
coef(summary(lmer(posthba1cfinal ~ group_modeluncal + period + (1|study_id), data = data_modeluncal_concordant_DPP4)))
confint(lmer(posthba1cfinal ~ group_modeluncal + period + (1|study_id), data = data_modeluncal_concordant_DPP4))


# Concordant: SGLT2
data_modeluncal_concordant_SGLT2                 <- subset(subset_3treatments_merged_modeluncal, subset_3treatments_merged_modeluncal$concordant_tx =="SGLT2")
data_modeluncal_concordant_SGLT2$treatment_class <- factor(data_modeluncal_concordant_SGLT2$treatment_class, levels = c("SGLT2", "DPP4", "TZD"))

# SGLT2i (best) vs. DPP4i
coef(summary(lmer(posthba1cfinal ~ treatment_class + period + (1|study_id), data = data_modeluncal_concordant_SGLT2[-which(data_modeluncal_concordant_SGLT2$treatment_class == "TZD"), ])))
confint(lmer(posthba1cfinal ~ treatment_class + period + (1|study_id), data = data_modeluncal_concordant_SGLT2[-which(data_modeluncal_concordant_SGLT2$treatment_class == "TZD"), ]))

# SGLT2i (best) vs. TZD
coef(summary(lmer(posthba1cfinal ~ treatment_class + period + (1|study_id), data = data_modeluncal_concordant_SGLT2[-which(data_modeluncal_concordant_SGLT2$treatment_class == "DPP4"), ])))
confint(lmer(posthba1cfinal ~ treatment_class + period + (1|study_id), data = data_modeluncal_concordant_SGLT2[-which(data_modeluncal_concordant_SGLT2$treatment_class == "DPP4"), ]))

# SGLT2i (best) vs. pooled other
coef(summary(lmer(posthba1cfinal ~ group_modeluncal + period + (1|study_id), data = data_modeluncal_concordant_SGLT2)))
confint(lmer(posthba1cfinal ~ group_modeluncal + period + (1|study_id), data = data_modeluncal_concordant_SGLT2))


# Concordant: TZD
data_modeluncal_concordant_TZD                 <- subset(subset_3treatments_merged_modeluncal, subset_3treatments_merged_modeluncal$concordant_tx =="TZD")
data_modeluncal_concordant_TZD$treatment_class <- factor(data_modeluncal_concordant_TZD$treatment_class, levels = c("TZD", "DPP4", "SGLT2"))

# TZD (best) vs. SGLT2i
coef(summary(lmer(posthba1cfinal ~ treatment_class + period + (1|study_id), data = data_modeluncal_concordant_TZD[-which(data_modeluncal_concordant_TZD$treatment_class == "SGLT2"), ])))
confint(lmer(posthba1cfinal ~ treatment_class + period + (1|study_id), data = data_modeluncal_concordant_TZD[-which(data_modeluncal_concordant_TZD$treatment_class == "SGLT2"), ]))

# TZD (best) vs. DPP4i
coef(summary(lmer(posthba1cfinal ~ treatment_class + period + (1|study_id), data = data_modeluncal_concordant_TZD[-which(data_modeluncal_concordant_TZD$treatment_class == "DPP4"), ])))
confint(lmer(posthba1cfinal ~ treatment_class + period + (1|study_id), data = data_modeluncal_concordant_TZD[-which(data_modeluncal_concordant_TZD$treatment_class == "DPP4"), ]))

# TZD (best) vs. pooled other
coef(summary(lmer(posthba1cfinal ~ group_modeluncal + period + (1|study_id), data = data_modeluncal_concordant_TZD)))
confint(lmer(posthba1cfinal ~ group_modeluncal + period + (1|study_id), data = data_modeluncal_concordant_TZD))

# concordant: overall 
coef(summary(lmer(posthba1cfinal ~ group_modeluncal + period + (1|study_id), data = data_analysis_modeluncal)))
confint(lmer(posthba1cfinal ~ group_modeluncal + period + (1|study_id), data = data_analysis_modeluncal))
