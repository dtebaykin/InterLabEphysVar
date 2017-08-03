# Project: Interlab ephys variance modeling
# Author: Dmitry Tebaykin
# Contact: dmitry.tebaykin@gmail.com
# Version: v1 (2 August 2017)

### Compare NeuroElectro data to never-before-seen AIBS data before and after model adjustment.
### Run requirements: DataProcessing.R, FeatureSelection.R


# Comparison: The frequency of features being picked per EP
# Try correcting NeuroElectro data for Allen Institute for Brain Science brain data
aibs_aggregated_ephys_v3 <- read.csv("~/Dropbox/dmitry/ephys_datasets/aibs_aggregated_ephys_v3.csv")
aibs_neuroelectro_cell_types_mapping <- read.csv("~/Dropbox/dmitry/ephys_datasets/aibs_neuroelectro_cell_types_mapping.csv")

aibs_aggregated_ephys_v3$NeuronName = apply(aibs_aggregated_ephys_v3, 1, function(x) {
  aibs_neuroelectro_cell_types_mapping[which(x[2] == aibs_neuroelectro_cell_types_mapping$transgene), "NeuronName"]
})

test_data = subset(aibs_aggregated_ephys_v3, NeuronName != "")
names(test_data)[names(test_data) == 'input_resistance_mohm'] <- 'rin'
names(test_data)[names(test_data) == 'vrest'] <- 'rmp'
test_data$apamp = test_data$peak_v_long_square - test_data$threshold_v_long_square
names(test_data)[names(test_data) == 'half_width'] <- 'aphw'
test_data$ahpamp = test_data$threshold_v_long_square - test_data$fast_trough_v_long_square
names(test_data)[names(test_data) == 'threshold_i_long_square'] <- 'rheo'
names(test_data)[names(test_data) == 'threshold_v_long_square'] <- 'apthr'
test_data$cap = (test_data$tau / test_data$rin) * 1000
test_data$maxfreq = 1 / (test_data$avg_isi / 1000)
test_data$adratio = 1 / ((test_data$adaptation * 10) +1)

# Only plot RMP and APhw
ephys_aibs <- c('rmp','aphw')

test_data = aggregate(test_data[,ephys_aibs], by = list(test_data$NeuronName), median, na.rm = TRUE)
colnames(test_data) = c("NeuronName", setdiff(colnames(test_data), "Group.1"))

test = matrix(0, nrow(test_data), length(features) - 1)
test = as.data.frame(test)
colnames(test) = setdiff(features, c("NeuronName"))
test_data = cbind(test_data, test)

rm(test)

test_data$AnimalAge = 57.5
test_data$Species = as.factor("Mice")
test_data$Strain = as.factor("C57BL")
test_data$RecTemp = 34
test_data$PubYear = 2016
test_data$external_0_Mg = 1
test_data$external_0_Na = 153.25
test_data$external_0_K = 2.5
test_data$external_0_Cl = 135
test_data$external_0_Ca = 2
test_data$external_0_glucose = 12.5
#test_data$external_0_kynur = 1
#test_data$external_0_picro = 0.1
test_data$internal_0_K = 130
test_data$internal_0_Na = 20.3
test_data$internal_0_HEPES = 10
test_data$internal_0_Cl = 4
test_data$internal_0_ATP = 4
test_data$internal_0_GTP = 0.3
test_data$internal_0_EGTA = 0.3
test_data$internal_0_Mg = 4

test_data$fold = 2


# Predict NE data using other folds, switch base 
# Shifted EP to AIBS = obs_ne - pred_ne + shift_aibs
# shift_aibs is a model trained on all NE data, predict EP using AIBS metadata per NT
result = data.frame("ep" = character(0), "NeuronName" = character(0), "Ne_obs" = numeric(0), "NE_corr" = numeric(0))
for (ep in ephys_aibs) {
  rf_data = fit_data[!is.na(fit_data[,ep]),]
  rf_data$fold = 1
  sel_feat = colnames(selected_features[ep, which(selected_features[ep,] >= 8)])
  rf_data = rbind(rf_data[, c(ep, sel_feat, "fold")], test_data[!is.na(test_data[,ep]), c(ep, sel_feat, "fold")])
  frml = as.formula(paste(ep, " ~ ", paste(sel_feat, collapse = "+")))
  rf_model = cforest(frml, rf_data[rf_data$fold == 1,], controls = cforest_control(teststat = "max",
                                                                                   testtype = "Teststatistic",
                                                                                   mincriterion = x,
                                                                                   savesplitstats = FALSE,
                                                                                   ntree = y, mtry = z, replace = T,
                                                                                   trace = T))
  obs_ne = rf_data[rf_data$fold == 1 & rf_data$NeuronName %in% as.character(rf_data[rf_data$fold == 2, "NeuronName"]),]
  pred_ne = predict(rf_model, newdata = obs_ne)
  
  shift_aibs_base = rf_data[rf_data$fold == 2,]
  shift_aibs_base$base = predict(rf_model, newdata = shift_aibs_base)
  
  shift_aibs = obs_ne[, c(ep, "NeuronName")]
  shift_aibs$base = apply(shift_aibs, 1, function(x) {
    subset(shift_aibs_base, NeuronName == x[2])$base
  })
  
  corrected_y = obs_ne[,ep] - pred_ne + shift_aibs$base
  
  result = rbind(result, data.frame("ep" = ep, "NeuronName" = shift_aibs$NeuronName, "NE_obs" = obs_ne[,ep], "NE_corr" = as.vector(corrected_y)))
}

# After metadata adjustment compare MSEs
plot_data = result

plot_data[,log10_ephys] <- 10 ** plot_data[,log10_ephys]

plot_data = melt(plot_data, id = c("ep", "NeuronName"))
plot_data$NeuronName = droplevels(plot_data$NeuronName)
toBind = subset(test_data, fold == 2)
toBind = toBind[, c(ephys_aibs, "NeuronName")]
toBind = melt(toBind, id = "NeuronName")
colnames(toBind) = c("NeuronName", "ep", "value")
toBind$variable = "AIBS"
plot_data = rbind(plot_data, toBind)
rm(toBind)
colnames(plot_data) = c("ep", "NeuronName", "Dataset", "value")
plot_data$ep = droplevels(plot_data$ep)
plot_data = subset(plot_data, NeuronName %in% c("Neocortex basket cell", "Neocortex pyramidal cell layer 2-3"))
plot_data$NeuronName = droplevels(plot_data$NeuronName)
levels(plot_data$NeuronName) = c("Basket", "L2/3 pyramidal")
levels(plot_data$ep) = c("Vrest", "APhw")
plot_data$shape = "mean"

toBind = data.frame()
toBind = aibs_aggregated_ephys_v3[,c("NeuronName", "vrest", "half_width")]
colnames(toBind) = c("NeuronName", "Vrest", "APhw")
toBind = melt(toBind)
toBind = subset(toBind, NeuronName %in% c("Neocortex basket cell", "Neocortex pyramidal cell layer 2-3"))
toBind$NeuronName = droplevels(toBind$NeuronName)
levels(toBind$NeuronName) = c("Basket", "L2/3 pyramidal")
colnames(toBind) = c("NeuronName", "ep", "value")
toBind$Dataset = "AIBS"
toBind$shape = "point"
plot_data = rbind(plot_data, toBind)
rm(toBind)

p = ggplot(plot_data, aes(y = value, x = NeuronName, colour = Dataset)) + 
  geom_boxplot(aes(fill = Dataset), position = position_dodge(width = 0.8), alpha = 0.3, outlier.alpha = 0)  + 
  geom_point(aes(shape = shape), na.rm=TRUE, position = position_jitterdodge(jitter.width = 0.3, jitter.height = 0, dodge.width = 0.7)) +
  scale_shape_manual(values = c(1, 16)) +
  facet_wrap( ~ ep, scales = "free_y") + 
  labs(x = "Neuron Type", y = "Electrophysiological measurements") +
  scale_color_manual(values=c("#000000", "#e29b20", "#349438")) +
  scale_fill_manual(values=c("#000000", "#e29b20", "#349438")) +
  theme_few(15) + 
  theme(axis.text.x = element_text(angle = 90,hjust =1, vjust =0.5))
print(p)
ggsave(plot=p,height=6,width=10,dpi=1000, filename= paste0(getwd(), "/Plots/AIBS_validation.pdf"), useDingbats=FALSE)

test = aggregate(value ~ ep + NeuronName + Dataset, plot_data, mean)
test_AIBS = subset(test, Dataset == "AIBS")
test_NE = setdiff(test, test_AIBS)
list_ep = c("rmp", "aphw")
for (i in unique(test$ep)) {
  ep = list_ep[which(unique(test$ep) == i)]
  for (j in unique(test$NeuronName)) {
    test_NE[which(test_NE$ep == i & test_NE$NeuronName == j), "value"] = abs(test_NE[which(test_NE$ep == i & test_NE$NeuronName == j), "value"] - test_AIBS[which(test_AIBS$ep == i & test_AIBS$NeuronName == j), "value"])
    # median(data_models[!is.na(data_models[,ep]),ep])
  }
}

ggplot(test_NE, aes(x = NeuronName, y = value, color = Dataset)) + 
  geom_point(size = 4) +
  facet_wrap( ~ ep, scales = "free_y") + 
  labs(x = "Neuron Type", y = "Electrophysiological measurements") +
  scale_color_manual(values=c("#246ccf", "#349438")) +
  theme_bw(15) + 
  theme(axis.text.x = element_text(angle = 90,hjust =1, vjust =0.5))