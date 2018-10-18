# Project: Interlab ephys variance modeling
# Author: Dmitry Tebaykin
# Contact: dmitry.tebaykin@gmail.com
# Version: v1 (2 August 2017)

### Feature selection and random forest script. Ephys models are created here.
### Run requirements: DataProcessing.R and CVfunc.R

### Find optimal parameters for cforest. Tried the following value ranges:
# x in 0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8
# y in 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 500
# z in 5:20

# optimal parameters for cforest: (replace = T) <- R^2 goes down by 0.1 if we do not sample with replacement
x = 0.05
y = 110
z = 10

# Variable importance
n_varimp_results = list()
n_aicc_results = list()
n_choice_results = list()
for (nrun in 1:10) {
  run_data = fit_data[sample(nrow(fit_data)),]
  run_data$fold = 0
  pmids <- unique(run_data$Pmid)
  for (k in 1:10) {
    split <- pmids[seq(k, length(pmids), 10)]
    run_data[run_data$Pmid %in% split,]$fold = k
  }
  
  varimp_results = list()
  aicc_results = list()
  choice_results = list()
  for (fold in 1:10) {
    run_result = list()
    for (ep in ephys_interest) {
      print(paste0("Fold: ", fold, ", ep: ", ep))
      
      rf_data = run_data[!is.na(run_data[,ep]),]
      frml =  as.formula(paste(ep, " ~ ", paste(features, collapse = "+")))
      run_result[ep] = cforest(frml, rf_data[rf_data$fold != fold, ], controls = cforest_control(teststat = "max",
                                                                                                 testtype = "Teststatistic",
                                                                                                 mincriterion = x,
                                                                                                 savesplitstats = FALSE,
                                                                                                 ntree = y, mtry = z, replace = T,
                                                                                                 trace = T))
    }
    print("Step 1 done")
    
    varimp_result = data.frame(matrix(0, 0, length(features)))
    for (ep in ephys_interest) {
      varimp_result = rbind(varimp_result, varimp(run_result[[ep]]))
    }
    rownames(varimp_result) = ephys_interest
    colnames(varimp_result) = features
    print("Step 2 done")
    
    varimp_results[[fold]] = varimp_result
    
    # Compare models using AIC: 2k - 2ln(L), substituting RMSE for L, assumption: real SDs of response variables are the same (within each model).
    # For limited num samples: AICc = AIC + 2 * k * (k + 1) / (n - k - 1)
    aicc_result = data.frame(matrix(0, 0, length(features)))
    for (i in 1:length(ephys_interest)) {
      ep = ephys_interest[i]
      rf_data = run_data[!is.na(run_data[,ep]),]
      sigma = sd(rf_data[rf_data$fold != fold, ep])
      n = nrow(rf_data)
      
      temp_aicc = c()  
      for (k in 1:length(features)) {
        print(paste0("EP: ", ep, ", Formula #", k))
        frml = as.formula(paste(ep, " ~ ", paste(colnames(sort(varimp_result[i,], decreasing = T)[1:k]), collapse = " + ")) )
        
        rf = cforest(frml, rf_data[rf_data$fold != fold, ], controls = cforest_control(teststat = "max",
                                                                                       testtype = "Teststatistic",
                                                                                       mincriterion = x,
                                                                                       savesplitstats = FALSE,
                                                                                       ntree = y, mtry = z, replace = T,
                                                                                       trace = T))
        rf_pred = predict(rf, newdata = rf_data[rf_data$fold != fold, ], type = "response")
        logL = (-n/2)*log(2 * pi * sigma^2) + ( -1 / (2 * sigma^2) * mse(rf_data[rf_data$fold != fold, ep], rf_pred) * n) 
        aicc = 2 * k - 2 * logL + 2 * k * (k + 1) / (n - k - 1)
        temp_aicc = c(temp_aicc, aicc)
      }
      aicc_result = rbind(aicc_result, temp_aicc)
    }
    rownames(aicc_result) = ephys_interest
    colnames(aicc_result) = 1:length(features)
    print("Step 3 done")
    
    #train and test new models
    aicc_choice = apply(aicc_result, 1, which.min)
    choice_result = data.frame(ep = ephys_interest, custom_rf = numeric(11), all = numeric(11))
    for (i in 1:length(ephys_interest)) {
      ep = ephys_interest[i]
      rf_data = run_data[!is.na(run_data[,ep]),]
      
      frml = as.formula(paste(ep, " ~ ", paste(colnames(sort(varimp_result[i,], decreasing = T)[1:aicc_choice[[ep]]]), collapse = " + ")) )
      print(frml)
      
      frml_all = paste(ep, " ~ ", paste(features, collapse = " + "))
      frml_all = gsub("NeuronName", "1", frml)
      
      for(level in unique(rf_data$NeuronName)){
        frml_all = paste(frml_all, paste("NeuronName", make.names(level), sep = "_"), sep = "+")
      }
      frml_all = as.formula(gsub("1\\+", "", frml_all))
      
      runList = convertNN(rf_data, frml)
      rf_data = runList[[1]]
      frml = as.formula(runList[[2]])
      
      rf = randomForest(frml, rf_data[rf_data$fold != fold, ], na.action = na.omit)
      rf_pred = predict(rf, newdata = rf_data[rf_data$fold == fold, ], type = "response")
      choice_result[i, 2] = get_R2(rf_data[rf_data$fold == fold, ep], rf_pred)
      
      rf = randomForest(frml_all, rf_data[rf_data$fold != fold, ], na.action = na.omit)
      rf_pred = predict(rf, newdata = rf_data[rf_data$fold == fold, ], type = "response")
      choice_result[i, 3] = get_R2(rf_data[rf_data$fold == fold, ep], rf_pred)
    }
    print(paste0("Run: ", nrun, "Step 4 done"))
    aicc_results[[fold]] = aicc_result
    choice_results[[fold]] = choice_result
  }
  n_aicc_results[[nrun]] = aicc_results
  n_choice_results[[nrun]] = choice_results
  n_varimp_results[[nrun]] = varimp_results
}

### Display results of the feature selected custom models

# Number of times each feature has been chosen per EP
plot_data = data.frame(matrix(0,length(ephys_interest), length(features)))
rownames(plot_data) = ephys_interest
colnames(plot_data) = features
for (k in 1:10) {
  for (i in 1:length(ephys_interest)) {
    aicc_index = apply(aicc_results[[k]][i,], 1, which.min)
    cols = colnames(sort(varimp_results[[k]][i,], decreasing = T)[1: aicc_index])
    plot_data[i, cols] = plot_data[i, cols] + 1
  }
}
colnames(plot_data) = new_features
rownames(plot_data) = response_vars_plot_names
pheatmap(plot_data, fontsize = 14, cluster_rows = F, cluster_cols = F)
colnames(plot_data) = features

# Save features that were selected for each EP
selected_features = plot_data

# Relative AICc per EP. 
plot_data = as.data.frame(t(scale(t(aicc_results[[1]]))))
plot_data$fold = 1
for (k in 2:10) {
  toBind = as.data.frame(t(scale(t(aicc_results[[k]]))))
  toBind$fold = k
  plot_data = rbind(plot_data, toBind)
}
plot_data$ep = ephys_interest
plot_data = melt(plot_data, id = c("ep", "fold"))

ggplot(plot_data, aes(x = variable, y = value, color = fold, group = fold)) + 
  geom_point() + 
  geom_line() + 
  facet_wrap(~ep, ncol = 3) + 
  labs(x = "Features, in order of importance", y = "AICc score") +
  theme_bw(10)

new_feat_levels = colnames(sort(varimp_results[[k]]["rin",], decreasing = T))
for (i in 1:length(solns)) {
  new_feat_levels[which(new_feat_levels == solns[i])] = new_solns[i]
}
levels(plot_data$variable) = new_feat_levels
ggplot(subset(plot_data, ep == "rin"), aes(x = variable, y = value, color = fold, group = fold)) + 
  geom_point() + 
  geom_line() + 
  labs(x = "Features, in order of importance", y = "AICc score") +
  theme_bw(15) +
  theme(axis.text.x = element_text(angle = 90,hjust =1, vjust =0.5))