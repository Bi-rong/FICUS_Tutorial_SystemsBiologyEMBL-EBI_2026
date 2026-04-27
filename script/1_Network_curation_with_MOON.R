# installing and loading FICUS ------------------------------------------------
path <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

#if (!requireNamespace("devtools", quietly = TRUE)){
#  install.packages("devtools")}
#devtools::install_github("https://github.com/saezlab/ficus")
library(ficus)

# retrieve the small patient cohorts ------------------------------------------
load('../data/cohortA.RData')
load('../data/cohortB.RData')

# setting some variables for MOON ---------------------------------------------
PKN_path <- '../data/clean_omnipath_PKN.RData' # data for initial PKN
min_size_PKN = 20 # minimal size for PKN 
significant_input_threshold <- 2 # threshold to filter TF activities 
n_steps <- 6 # number of steps during network pruning
use_subset = F # if desired, can select subset of proteins for network

# thresholds for soft network reduction
primary_threshold1 = 0.5
secondary_threshold1 = 0.25

# CollecTRI net used in MOON 
load('../data/Collectri_PROGENy_networks.RData') 
net$confidence <- NA
net <- net[,c('source', 'confidence', 'target', 'mor')]

# running MOON ----------------------------------------------------------------
# this section takes around 5 minutes to run  
output_folder <- '../output/MOON_cohortA/' 
start <- 1
end <- length(cosmos_inputs_A) 
run_MOON(cosmos_inputs=cosmos_inputs_A, 
         significant_input_threshold=significant_input_threshold,
         PKN_path=PKN_path, 
         n_steps=n_steps, 
         net=net,
         use_subset=use_subset, 
         min_size_PKN=min_size_PKN, 
         primary_threshold=primary_threshold1, 
         secondary_threshold=secondary_threshold1, 
         output_folder=output_folder, 
         start=start, 
         end=end)




