# Project: Interlab ephys variance modeling
# Author: Dmitry Tebaykin
# Contact: dmitry.tebaykin@gmail.com
# Version: v1 (2 August 2017)

### Univariate analysis of ephys properties and major ion reversal potentials
### Run DataProcessing.R first

# Overhead for Reversal Potentials
fit_data <- subset(fit_data, is.finite(Veq_Ca))
fit_data <- subset(fit_data, is.finite(Veq_Mg))
fit_data <- subset(fit_data, is.finite(Veq_Na))
fit_data <- subset(fit_data, is.finite(Veq_K))
fit_data <- subset(fit_data, is.finite(Veq_Cl))

### Check univariate relationships between all features and EPs (after regressing out the effect of other features)
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

rf_result <- data.frame(matrix(NA, ncol = length(ions), nrow = length(ephys_interest)))
colnames(rf_result) <- paste0("Veq_", ions)
rownames(rf_result) <- ephys_interest

### Generate univariate plots: each EP for each E_ion
for (ep in ephys_interest) {
  for (ion in ions) {
    if ( !((ep == "rmp" & ion == "Na") | (ep == "aphw" & ion == 'Cl'))) {
      next
    }
    ft = paste0("Veq_", ion)
    frml = paste(colnames(selected_features[ep, which(selected_features[ep,] >= 9)]), collapse = "+")
    cat(paste0("EP: ", ep, ", Feature: ", ft, "\n"))
    
    run_data = fit_data[!is.na(fit_data[,ep]),]
    run_data = run_data[complete.cases(run_data[,c(ep, features)]),]
    run_data = run_data[sample(nrow(run_data)),]
    
    rf_data = run_data
    
    rf_data$group = "All"
    rf_data[rf_data$NeuronName == "Hippocampus CA1 pyramidal cell",]$group = "CA1"
    rf_data[rf_data$NeuronName == "Neocortex basket cell",]$group = "Basket"
    
    frml_corr = paste0(ep, "~", gsub("^\\+|\\+$", "", gsub( paste0(
      paste0("internal_0_", ion), "|",
      paste0("external_0_", ion), "|",
      "RecTemp"),
      "", frml)))
    frml_corr = gsub("\\+\\+", "\\+", frml_corr)
    frml_uni = as.formula(paste("ep_shifted ~ ", ft))
    
    if (grepl("NeuronName", frml_corr)) {
      rf_formula = gsub("NeuronName", "1", frml_corr)
      rf_data$NeuronName = droplevels(rf_data$NeuronName)
      
      for(level in unique(rf_data$NeuronName)){
        newname = gsub("\\.\\.", "\\.", level)
        rf_data[paste("NeuronName", make.names(newname), sep = "_")] <- ifelse(rf_data$NeuronName == level, 1, 0)
        rf_formula = paste(rf_formula, paste("NeuronName", make.names(newname), sep = "_"), sep = "+")
      }
      rf_data$NeuronName = NULL
      frml_corr = gsub("1\\+", "", rf_formula)
    }
    frml_corr = as.formula(frml_corr)
    
    rf = randomForest(frml_corr, rf_data, importance = F, na.action = na.omit)
    rf_pred = predict(rf, newdata = rf_data, type = "response")
    rf_data$ep_shifted = rf_data[,ep] - rf_pred + rep(median(rf_data[, ep]), nrow(rf_data))
    
    plot_data = rf_data
    
    # Filter out Eion < 1st quantile and > 3rd quantile. They likely do not represent real biological reversal potentials
    if (ion == "Na") {
      plot_data = subset(plot_data, Veq_Na > 55 & Veq_Na < 110)
    }
    if (ion == "Cl") {
      plot_data = subset(plot_data, Veq_Cl > -82 & Veq_Cl < -49)
    }
    if (ion == "K") {
      plot_data = subset(plot_data, Veq_K > -106 & Veq_K < -97)
    }
    
    if (ep %in% log10_ephys) {
      plot_data$ep_shifted = 10**plot_data$ep_shifted
    }
    
    rf_result[ep, ft] = -log10(lmp(lm(frml_uni, data = plot_data)))
    
    temp = subset(plot_data, group != "All")
    temp$group = "All"
    plot_data = rbind(plot_data, temp)
    
    p = ggplot(plot_data, aes_string(x = ft, y = "ep_shifted", color = "group")) +
      geom_point(shape = 19, alpha = 0.35, size = 2) +
      stat_smooth(method = "lm", size = 1.5, level = 0.95, se = F) +
      
      scale_color_manual(values = c("#000000", "#8A458A", "darkblue")) +
      ggtitle(paste0("All NTs, adj.: ", ep, " ~ ", ft)) +
      labs(x = ft, y = ep) +
      theme_few(15) +
      facet_wrap(~group, scales = "free_x")
    ggsave(plot=p,height=4, width=8,dpi=200, filename = paste0(getwd(), "/Plots/univariatePlots/", ep, "Vs", ft, ".pdf"), useDingbats=FALSE)
    
    print(paste("All p-val: ",lmp(lm(frml_uni, data = subset(plot_data, group == "All")))))
    print(paste("All cor: ", cor(subset(plot_data, group == "All")$ep_shifted, subset(plot_data, group == "All")[,ft])))
    
    print(paste("CA1 p-val: ",lmp(lm(frml_uni, data = subset(plot_data, group == "CA1")))))
    print(paste("CA1 cor: ", cor(subset(plot_data, group == "CA1")$ep_shifted, subset(plot_data, group == "CA1")[,ft])))
    
    print(paste("Basket p-val: ",lmp(lm(frml_uni, data = subset(plot_data, group == "Basket")))))
    print(paste("Basket cor: ", cor(subset(plot_data, group == "Basket")$ep_shifted, subset(plot_data, group == "Basket")[,ft])))
  }
} # end of function

pheatmap(rf_result, show_colnames = T, fontsize = 14, cluster_rows = F, cluster_cols = F, scale = "none",
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100))