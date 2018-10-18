# Project: Interlab ephys variance modeling
# Author: Dmitry Tebaykin
# Contact: dmitry.tebaykin@gmail.com
# Version: v2 (17 October 2018)

### This script loads and processes the data. It must be run before any further analysis can take place.
### It requires ephys_data_curated.csv dataset file from the /Data folder

### Load required libraries
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

### Load the data
ions <- c("Na", "Ca", "Mg", "Cl", "K")
rev_ions <- c("Veq_Na", "Veq_Ca", "Veq_Mg", "Veq_Cl", "Veq_K")
solns <- c("external_0_Na", "external_0_K", "external_0_Cl", "external_0_Ca", "external_0_Mg", "external_0_glucose", "external_0_HEPES", "external_0_EGTA", "external_0_EDTA", "external_0_BAPTA", "external_0_ATP", "external_0_GTP",  
           "internal_0_Na", "internal_0_K", "internal_0_Cl", "internal_0_Ca", "internal_0_Mg", "internal_0_glucose", "internal_0_HEPES", "internal_0_EGTA", "internal_0_EDTA", "internal_0_BAPTA", "internal_0_ATP", "internal_0_GTP", "internal_0_Cs", "external_0_Cs")
relevant_solns <- c('external_0_Mg','external_0_Ca','external_0_Na','external_0_Cl','external_0_K',
                    'internal_0_Mg','internal_0_Ca','internal_0_Na','internal_0_Cl','internal_0_K',
                    "external_0_glucose", "internal_0_HEPES",  
                    "internal_0_EGTA", "internal_0_ATP", "internal_0_GTP")

log10_ephys <- c("cap", "rin", "tau", "rheo", "aphw", "maxfreq")

features = c('NeuronName','Species','Strain',
             'RecTemp','AnimalAge','PubYear', 
             'external_0_Mg','external_0_Ca','external_0_Na','external_0_Cl','external_0_K',
             'internal_0_Mg','internal_0_Ca','internal_0_Na','internal_0_Cl','internal_0_K',
             "external_0_glucose", "internal_0_HEPES",  
             "internal_0_EGTA", "internal_0_ATP", "internal_0_GTP")

new_features <- c("Neuron Type", "Species", "Strain", "Recording Temp", "Animal Age", "Publication Year",
                  "External [Mg]",
                  "External [Ca]",
                  "External [Na]",
                  "External [Cl]",
                  "External [K]",
                  "Internal [Mg]",
                  "Internal [Ca]",
                  "Internal [Na]",
                  "Internal [Cl]",
                  "Internal [K]",
                  "External [glucose]",
                  "Internal [HEPES]",
                  "Internal [EGTA]",
                  "Internal [ATP]",
                  "Internal [GTP]")

ephys_interest <- c("rin", "rmp","apthr","apamp", "aphw", "tau", "ahpamp", "rheo","maxfreq", "cap", "adratio")
response_vars_plot_names = c('Rin', 'Vrest', 'APthr', 'APamp', 'APhw', 'Tau', 'AHPamp', 'Rheo', 'FRmax', 'Cm', 'SFA')

ne_all_curated <- read.delim("Data/ephys_data_curated.csv", stringsAsFactors=FALSE)

### Filter and process the data
ne_all_filter_solns <- ne_all_curated[apply(ne_all_curated[, solns], 1, function(y) !all(is.na(y))),]

# Assign Na's: a small number to avoid 'divide by zero' issue when computing Veq
rev_data <- ne_all_filter_solns
rev_data[is.na(rev_data$AnimalAge), "AnimalAge"] <- median(rev_data$AnimalAge, na.rm = TRUE)
rev_data[is.na(rev_data$RecTemp), "RecTemp"] <- median(rev_data$RecTemp, na.rm = TRUE)
rev_data[,solns][is.na(rev_data[, solns])] <- 0.000001

# Constants (at rest), P's are 'typical' values from neuroscience course / textbook / the internet
Far <- 9.6485 * 10000 / 1000
K <- 273.15
R <- 8.314
P_K <- 1
P_Na <- 0.05
P_Cl <- 0.45

### Calculate reversal potentials
for (ion in ions) {
  new_col <- data.frame(a = R * (K + rev_data[,"RecTemp"]) / Far / 2 * log(rev_data[, paste("external_0_", ion, sep = "")] / rev_data[, paste("internal_0_", ion, sep = "")]))
  colnames(new_col) <- paste("Veq_", ion, sep = "")
  rev_data <- cbind(rev_data, new_col)
}
rm(new_col, ion)

# Filter for in vitro, patch-clamp studies in rats, mice or guinea pigs
rev_data_filtered <- subset(rev_data, AnimalAge > 0)
rev_data_filtered[,log10_ephys] <- log10(rev_data_filtered[,log10_ephys])
rev_data_filtered <- subset(rev_data_filtered,  PrepType == 'in vitro'
                            & (Species == 'Rats' | Species == 'Mice' | Species == 'Guinea Pigs'))
rev_data_filtered <- subset(rev_data_filtered,  internal_0_Cs < 0.001 & external_0_Cs < 0.001 ) # current clamp recordings
rev_data_filtered <- subset(rev_data_filtered,  ElectrodeType == "Patch-clamp")
rev_data_filtered <- subset(rev_data_filtered, external_0_HEPES < 0.001) # remove cell cultures that were incorrectly curated as in vitro
rev_data_filtered <- droplevels(rev_data_filtered)
rev_data_filtered[,log10_ephys] <- 10 ** rev_data_filtered[,log10_ephys]

# Create data for models to be used in feature_selection.R
data_models <- rev_data_filtered
data_models$Species <- factor(data_models$Species)

# Fix RMP based on JxnPotential and JxnOffset
data_models$rmp = as.numeric(apply(data_models, 1 , function(x) {
  if (x['JxnPotential'] == "Corrected") {
    if (is.na(x['JxnOffset'])) {
      as.numeric(x['rmp']) + median(abs(data_models[,'JxnOffset']), na.rm = TRUE)
    } else {
      as.numeric(x['rmp']) + abs(as.numeric(x['JxnOffset']))
    }  
  } else {
    x['rmp']
  }
}))

# Fix AP threshold based on JxnPotential and JxnOffset
data_models$apthr = as.numeric(apply(data_models, 1 , function(x) {
  if (x['JxnPotential'] == "Corrected") {
    if (is.na(x['JxnOffset'])) {
      as.numeric(x['apthr']) + median(abs(data_models[,'JxnOffset']), na.rm = TRUE)
    } else {
      as.numeric(x['apthr']) + abs(as.numeric(x['JxnOffset']))
    }  
  } else {
    x['apthr']
  }
}))

data_models_x <- as.data.frame(unclass(data_models[,c(features, "Pmid")]))
fit_data <- cbind(data_models_x, data_models[,ephys_interest])

# Clean-up
rm(data_models_x, rev_data, ne_all_filter_solns, rev_data_filtered)
