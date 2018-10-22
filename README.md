# InterLabEphysVar
Modeling sources of inter-laboratory variability in electrophysiological properties of mammalian neurons

## Abstract

Patch-clamp electrophysiology is widely used to characterize neuronal electrical phenotypes. However, there are no standard experimental conditions for in vitro whole-cell patch-clamp electrophysiology, complicating direct comparisons between datasets. Here, we sought to understand how basic experimental conditions differ among labs and how these differences might impact measurements of electrophysiological parameters. We curated the compositions of external bath solutions (ACSF), internal pipette solutions, and other methodological details from 509 published neurophysiology articles studying rodent neurons. We found that very few articles used the exact same experimental solutions as any other and some solution differences stem from recipe inheritance from adviser to advisee as well as changing trends over the years. Next, we used statistical models to understand how the use of different experimental conditions impacts downstream electrophysiological measurements such as resting potential and action potential width. While these experimental condition features could explain up to 43% of the study-to-study variance in electrophysiological parameters, the majority of the variability was left unexplained. Our results suggest that there are likely additional experimental factors that contribute to cross-laboratory electrophysiological variability, and identifying and addressing these will be important to future efforts to assemble consensus descriptions of neurophysiological phenotypes for mammalian cell types.

## Setup instructions
1) Clone (or download) the repository files
2) Set InterLabEphysVar as working directory

## Using the existing models
Use modelUseCase.R, the file is fully commented, installs and loads all required libraries and it loads a pre-calculated .RData file. The current use case creates scatterplots and calcualtes relevant stats of Vrest (after removing explained variance from other metadata) versus individual solution components. The plots can be located in Plots/univariatePlots.

## Full feature selection run
Start with DataProcessing.R, adjust features as necessary:
features - list of features to use for feature selection during model creation
new_features - the same features as above, these names are used for data plotting

Next, check ephys_interest - does it include all the ephys properties that should be modeled?
Finally, check line 73: ne_all_curated <- read.delim("Data/ephys_data_curated.csv", stringsAsFactors=FALSE)
Either replace the existing spreadsheet or point the script to a new one. The filtering steps should be robust to any dataset and can be left as-is, unless the user desires to change them explicitly. 

Next step: run FeatureSelection.R - be prepared for a ~24 hour runtime, it does ten 10-fold feature selection CVs for each ephys property. In the paper I ran 100 10-fold CVs, that took several days to run on a computational server. When the run is done - feature_selection table will hold the desired information. The script will make several figures regarding model performance.

Other R files correspond to specific areas of the paper and mostly include figure generation and data wrangling code. AIBS.R contains collaboration code with Allen Institute for Brain Science.