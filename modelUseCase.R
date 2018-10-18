### Example use case script: Modeling RMP with respect to solution components.
### Apply NeuroElectro ephys prediction models to remove as much explained variance as possible
### from RMP and plot univariate relationships of the resulting RMP values and solution components
# Author: Dmitry Tebaykin
# Contact: dmitry.tebaykin@gmail.com

# Load required and potentially useful packages
libs <- c('ggplot2','RColorBrewer','ggbiplot','reshape2', 'devtools',
          'formula.tools', 'randomForest', 'cowplot', 'cluster','caret',
          'pheatmap', 'ggthemes', 'dplyr','mlbench','party','Metrics',
          'splines2','data.table')
biocLite <- source("http://bioconductor.org/biocLite.R")

for (L in libs){
  if(!require(L, character.only = T)) {
    install.packages(L, dep = T, quiet = T)
    
    if (!require(L, character.only = T)) {
      suppressWarnings(biocLite(L))
    }
  }
}

sourceFiles <- c('CVfunc.R')
for (s in sourceFiles) {
  source(s)
}
rm(L, libs, s, sourceFiles)
theme_set(theme_cowplot())

load("Data/ephys_models.RData")

# Generate univariate plots after shifting the RMP values with the proper model
ep = "rmp"

# Create data frame with the ephys property of interest (RMP in this case)
rf_data = fit_data[!is.na(fit_data[,ep]),]

# Reshuffle the data rows for fold creation
rf_data = rf_data[sample(nrow(rf_data)),]

# Randomly assign folds to data rows
rf_data$fold = 0
# Enforce unique Pmid's for each fold
pmids <- unique(rf_data$Pmid)
for (f in 1:10) {
  split <- pmids[seq(f, length(pmids), 10)]
  rf_data[rf_data$Pmid %in% split,]$fold = f
}
rm(f, pmids, split)

for (i in 1:length(relevant_solns)) {
  # Create data copy for the current run
  run_data = rf_data
  run_data$pred = 0
  
  # Use features that were chosen at least 8 times during the FeatureSelection.R 10-fold CV process
  # The value "8" above is arbitrary and can be adjusted on an individual basis
  ep_plotname = response_vars_plot_names[grep(ep, ephys_interest)]
  feat_plotname = new_features[grep(relevant_solns[i], features)]
  sel_feat = colnames(selected_features[ep_plotname, which(selected_features[ep_plotname,] >= 8)])
  
  print(paste0("New model run: ", ep_plotname, " ~ ", feat_plotname))
  
  # Drop the feature we are currently plotting if it is part of the model - we don't want to regress out its effect
  sel_feat = setdiff(sel_feat, relevant_solns[i])
  
  # Build the formula
  frml = as.formula(paste(ep, " ~ ", paste(sel_feat, collapse = "+")))
  
  for (k in 1:10) {
    print(paste0("Modeling: ", ep_plotname, " ~ ", feat_plotname, ". Fold: ", k, " out of 10"))
    rf_model = cforest(frml, run_data[run_data$fold != k,], controls = cforest_control(teststat = "max",
                                                                 testtype = "Teststatistic",
                                                                 mincriterion = x,
                                                                 savesplitstats = FALSE,
                                                                 ntree = y, mtry = z, replace = T,
                                                                 trace = T))
    # Predict the ephys prop using the optimal model
    run_data[run_data$fold == k,]$pred = predict(rf_model, newdata = run_data[run_data$fold == k,])
  }
  
  # Remove model-explained variance from the EP, using median EP value as the base value to adjust around
  run_data$ep_shifted = run_data[,ep] - run_data$pred + rep(median(run_data[, ep]), nrow(run_data))
  
  # Plot the univariate relationships
  plot_data = run_data
  if (ep %in% log10_ephys) {
    plot_data$ep_shifted = 10**plot_data$ep_shifted
  }
  
  plot_data$group = "All"
  plot_data[plot_data$NeuronName == "Hippocampus CA1 pyramidal cell",]$group = "CA1"
  plot_data[plot_data$NeuronName == "Neocortex basket cell",]$group = "Basket"
  
  temp = subset(plot_data, group != "All")
  temp$group = "All"
  plot_data = rbind(plot_data, temp)
  
  # Remove textmining outliers for plotting
  minThresh = quantile(plot_data[,relevant_solns[i]], 0.02)
  maxThresh = quantile(plot_data[,relevant_solns[i]], 0.98)
  test = plot_data[which(plot_data[,relevant_solns[i]] >= minThresh &
                        plot_data[,relevant_solns[i]] <= maxThresh),]

  print(paste0("Saving plot to: ", paste0(getwd(), "/Plots/univariatePlots/", ep, "Vs", relevant_solns[i], ".pdf")))
  p = ggplot(plot_data, aes_string(x = relevant_solns[i], y = "ep_shifted", color = "group")) +
    geom_point(shape = 19, alpha = 0.2, size = 2) +
    stat_smooth(method = "lm", size = 1.5, level = 0.95, se = F) +
    
    scale_color_manual(values = c("#000000", "#8A458A", "darkblue")) +
    ggtitle(paste0("All NTs, adj.: ", ep_plotname, " ~ ", feat_plotname)) +
    labs(x = feat_plotname, y = ep_plotname) +
    theme_few(15) +
    facet_wrap(~group, scales = "free_x")
  ggsave(plot=p,height=4, width=8,dpi=200, filename = paste0(getwd(), "/Plots/univariatePlots/", ep, "Vs", relevant_solns[i], ".pdf"), useDingbats=FALSE)
  
  # Stats to go along with the figures
  frml_uni = as.formula(paste("ep_shifted ~ ", relevant_solns[i]))
  print(paste("All p-val: ", lmp(lm(frml_uni, data = subset(plot_data, group == "All")))))
  print(paste("All cor: ", cor(subset(plot_data, group == "All")$ep_shifted, subset(plot_data, group == "All")[,relevant_solns[i]])))
  
  print(paste("CA1 p-val: " ,lmp(lm(frml_uni, data = subset(plot_data, group == "CA1")))))
  print(paste("CA1 cor: ", cor(subset(plot_data, group == "CA1")$ep_shifted, subset(plot_data, group == "CA1")[,relevant_solns[i]])))
  
  print(paste("Basket p-val: ", lmp(lm(frml_uni, data = subset(plot_data, group == "Basket")))))
  print(paste("Basket cor: ", cor(subset(plot_data, group == "Basket")$ep_shifted, subset(plot_data, group == "Basket")[,relevant_solns[i]])))  
  print("Model complete.")
}
rm(k, i, ep_plotname, sel_feat, rf_model, temp, p)
