# Project: Interlab ephys variance modeling
# Author: Dmitry Tebaykin
# Contact: dmitry.tebaykin@gmail.com
# Version: v1 (2 August 2017)

### Compare different models and generate figures
### Feature selection, 100 runs with 90% train : 10% test data split. Next time paralellize it: mclapply or foreach
# Plot rf real vs predicted values for each ephys property
plot_rf_vals <- function(vals, rf_res) {
  plot_data <- as.data.frame(vals[, 1:(length(ephys_interest)*2)])
  colnames(plot_data) <- gsub("_\\d*$", "", colnames(plot_data))
  for (i in 2:nrun - 1) {
    startI = length(ephys_interest)*2*i + 1
    endI = length(ephys_interest)*2*(i + 1)
    toBind <- vals[, startI : endI]
    colnames(toBind) <- colnames(plot_data)
    plot_data <- rbind(plot_data, toBind)
  }
  
  #   toAttach <- data.frame(matrix(0, ncol = 1, nrow = nrun*2000))
  #   for (i in 0:9) {
  #     toAttach[(i*2000 + 1): ((i+1)*2000),] <- splits[i + 1]
  #   }
  #colnames(toAttach) = "splits"
  #plot_data <- cbind(plot_data, toAttach)
  #plot_data$splits <- as.factor(plot_data$splits)
  #levels(plot_data$splits) <- c(0.1, 0.14, 0.16, 0.16, 0.18, 0.18)
  #droplevels(plot_data$splits)
  
  # log10(Input Resistance), MΩ. R^2 = 0.28
  g1 = ggplot(plot_data, aes(y = log10(rin_pred), x = log10(rin_real))) + geom_point() + geom_abline(aes(slope = 1, intercept = 0)) + labs(x = "Measured values: log10(Rin), MΩ", y = "Predicted values: log10(Rin), MΩ") + xlim(1.5, 3) + ylim(1.5, 3) + ggtitle("") + theme_few(15)
  #plot(g1)
  g2 = ggplot(plot_data, aes(y = rmp_pred, x = rmp_real)) + geom_abline(aes(slope = 1, intercept = 0), color = "grey") + geom_point(alpha = 0.25) + geom_smooth(method = "lm", se = FALSE, color = "darkred") + xlim(-90, -40) + ylim(-90, -40) + labs(x = "Measured Vrest, mV", y = "Predicted Vrest, mV") + ggtitle(paste0("R^2: ", format(mean(as.data.frame(t(rf_res))$rmp), digits = 2, nsmall = 2))) + theme_few(15)
  plot(g2)
  ggsave(plot=g2,height=6,width=5,dpi=200, filename= paste0(getwd(), "/Plots/Real vs predicted RMP.pdf"), useDingbats=FALSE)
  
  g3 = ggplot(plot_data, aes(y = apthr_pred, x = apthr_real)) + geom_point() + geom_abline(aes(slope = 1, intercept = 0)) + xlim(-60, -20) + ylim(-60, -20) + ggtitle(paste0("R^2: ", format(mean(as.data.frame(t(rf_res))$apthr), digits = 2, nsmall = 2))) + theme_bw(15)
  g4 = ggplot(plot_data, aes(y = apamp_pred, x = apamp_real)) + geom_point() + geom_abline(aes(slope = 1, intercept = 0)) + xlim(50, 100) + ylim(50, 100) + ggtitle(paste0("R^2: ", format(mean(as.data.frame(t(rf_res))$apamp), digits = 2, nsmall = 2))) + theme_bw(15)
  g5 = ggplot(plot_data, aes(y = aphw_pred, x = aphw_real)) + geom_point() + geom_abline(aes(slope = 1, intercept = 0)) + xlim(-0.5, 0.5) + ylim(-0.5, 0.5) + ggtitle(paste0("R^2: ", format(mean(as.data.frame(t(rf_res))$aphw), digits = 2, nsmall = 2))) + theme_few(15)
  g6 = ggplot(plot_data, aes(y = tau_pred, x = tau_real)) + geom_point() + geom_abline(aes(slope = 1, intercept = 0)) + xlim(0, 2) + ylim(0, 2) + ggtitle(paste0("R^2: ", format(mean(as.data.frame(t(rf_res))$tau), digits = 2, nsmall = 2))) + theme_bw(15)
  g7 = ggplot(plot_data, aes(y = ahpamp_pred, x = ahpamp_real)) + geom_point() + geom_abline(aes(slope = 1, intercept = 0)) + xlim(0, 30) + ylim(0, 30) + ggtitle(paste0("R^2: ", format(mean(as.data.frame(t(rf_res))$ahpamp), digits = 2, nsmall = 2))) + theme_bw(15)
  
  #print(multiplot(g1, g2, g3, g4, g5, g6, g7, cols = 3))
}

## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

compare_models <- function() {
  for (ep in ephys_interest) {
    frmls = list(as.formula(paste(ep, " ~ NeuronName")))
    
    # Custom models for each ep, pick features that are chosen at least 90% of the time (from basic meta), NeuronName workaround for capacitance
    frml = paste(ep, " ~ ", paste(features[colnames(selected_features[ep, which(selected_features[ep,] >= 9)]) %in% basic_metadata], collapse = "+"))
    if (!grepl("NeuronName", frml)) {
      frml = paste(frml, "NeuronName")
    }
    frmls[[2]] = as.formula(frml)   
    
    # Custom models for each ep, pick features that are chosen at least 90% of the time, NeuronName workaround for capacitance
    frml = paste(ep, " ~ ", paste(colnames(selected_features[ep, which(selected_features[ep,] >= 9)]), collapse = "+"))
    if (!grepl("NeuronName", frml)) {
      frml = paste(frml, "NeuronName")
    }
    frmls[[3]] = as.formula(frml)
    
    run_data = fit_data[!is.na(fit_data[,ep]),]
    run_data = run_data[sample(nrow(run_data)),]
    
    # Force unique pmid's into splits
    run_data$fold = 0
    pmids <- unique(run_data$Pmid)
    for (i in 1:10) {
      split <- pmids[seq(i, length(pmids), 10)]
      run_data[run_data$Pmid %in% split,]$fold = i
    }
    
    # Force unique NTs and "Other" to not appear in test fold
    run_data[which(run_data$NeuronName == "Other"),]$fold <- 0
    run_data[which(run_data$NeuronName %in% names(table(run_data$NeuronName)[table(run_data$NeuronName) == 1])),]$fold <- 0
    
    # 10x cross-validation length(frmls)
    for (i in 1:length(frmls)) {
      for (k in 1:10) {
        print(paste0("Run: ", runnum, ", EP: ", ep, ", Model #", i, ", Fold #", k))
        
        if (grepl("NeuronName", frmls[[i]]) & (i < 7 | i > 12 )) {
          rf_formula = gsub("NeuronName", "1", frmls[[i]])
          rf_data = run_data
          rf_data$NeuronName = droplevels(rf_data$NeuronName)
          
          for(level in unique(rf_data$NeuronName)){
            rf_data[paste("NeuronName", make.names(level), sep = "_")] <- ifelse(rf_data$NeuronName == level, 1, 0)
            rf_formula = paste(rf_formula, paste("NeuronName", make.names(level), sep = "_"), sep = "+")
          }
          rf_data$NeuronName = NULL
          rf_formula = gsub("1\\+", "", rf_formula)
        } else {
          rf_formula = frmls[[i]]
          rf_data = run_data
        }
        
        rf_fit = randomForest(as.formula(rf_formula), data = rf_data[rf_data$fold != k,], importance = F)
        rf_pred = predict(rf_fit, newdata = rf_data[rf_data$fold == k,], type = "response")
        
        rf_result[[i]][paste0(ep, "_r2"), (runnum-1)*10 + k] <<- get_R2(run_data[run_data$fold == k, ep], rf_pred)
        rf_result[[i]][paste0(ep, "_mse"), (runnum-1)*10 + k] <<- mse(run_data[run_data$fold == k, ep], rf_pred)
        
        rf_vals[[i]][paste0(ep, "_real_", k)] <<- c(run_data[run_data$fold == k, ep], rep(NaN, 2000 - length(rf_pred)))
        rf_vals[[i]][paste0(ep, "_pred_", k)] <<- c(rf_pred, rep(NaN, 2000 - length(rf_pred)))
      } # end of 10x cross validation
    } 
  } 
  
  runnum <<- runnum + 1
}

### Run start for compare_models()
nrun = 10
rf_result = list()
rf_vals = list()
rf_vals_names = c()
for (k in 1:nrun) {
  rf_vals_names = c(rf_vals_names, c(paste0(ephys_interest, "_real_", k), paste0(ephys_interest, "_pred_", k)))
}

for (i in 1:3) {
  rf_temp = data.frame(matrix(0, ncol = nrun*10, nrow = length(ephys_interest) * 2))
  colnames(rf_temp) <- 1:nrun*10
  rownames(rf_temp) <- c(paste0(ephys_interest, "_r2"), paste0(ephys_interest, "_mse"))
  rf_result[[i]] <- rf_temp
  
  rf_temp = data.frame(matrix(NaN, ncol = length(ephys_interest)*2*nrun, nrow = 2000))
  colnames(rf_temp) = rf_vals_names
  rf_vals[[i]] <- rf_temp
}
rm(rf_temp, rf_vals_names)

runnum = 1
for (i in 1: nrun) {
  compare_models()
}

# Plot real vs predicted vals
plot_rf_vals(rf_vals[[3]], rf_result[[3]][1:length(ephys_interest),])


# Plot means of each model
models = c("Neuron Type", "NT and sel.basic.meta", "Selected all feat.")
plot_data <- as.data.frame(t(rf_result[[1]]))
plot_data$model = models[1]
for (i in 2:length(models)) {
  toBind = as.data.frame(t(rf_result[[i]]))
  toBind$model =  models[i]
  plot_data <- rbind(plot_data, toBind)
}

plot_data$model = factor(plot_data$model, levels = models)
plot_data[,grepl("mse", colnames(plot_data))] = apply(plot_data[,grepl("mse", colnames(plot_data))], 2, function(x) {
  scale(x)
})

# Convert to mean R^2 and MSE per run
rownames(plot_data) = 1:nrow(plot_data)
temp = plot_data[0,]
for (i in 1:(nrow(plot_data) / 10)) {
  temp[i, 1:22] = colMeans(plot_data[((i-1)*10 + 1) : (i*10), 1:22])
  temp[i, 23] = plot_data[i*10, 23]
}
plot_data = temp
rm(temp)
plot_data <- melt(plot_data, id = "model")

cbbPalette <- c("#F0E442", "#56B4E9", "#E69F00") # Solutions vs basic metadata

plot_data <- plot_data[grepl("r2", plot_data$variable),]
plot_data$variable = droplevels(plot_data$variable)
levels(plot_data$variable) = response_vars_plot_names

p = ggplot(plot_data, aes(x = variable, y = value)) + 
  geom_boxplot(aes(fill = model), position = position_dodge(width = 0.9)) +
  scale_fill_manual(values=cbbPalette) +
  ylim(-0.5,0.5) +
  labs(x = "Electrophysiological properties", y = "Out-of-sample R2") +
  theme_few(15)
print(p)
ggsave(plot=p,height=5,width=12,dpi=200, filename= paste0(getwd(), "/Plots/R2 model comparison.pdf"), useDingbats=FALSE)

