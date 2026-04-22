# FICUS Practical SystemsBiology2026
This file includes both explanation and the code for the practical. Separate R scripts are also available for each section containing only the code. This practical includes a demonstration of the main FICUS functions, allowing you to go from omics data to logic models. Two small patient cohorts are provided, one for demonstration purposes and one for you to use during the assignments. The assignments here are more guidelines on what can be tuned in the current framework, and what you might want to consider when using the tools. They are nothing more than a place where you can try out the FICUS tools, change variables, investigate the outputs and compare across data sets. The aim is to familiarize you with the methods used in FICUS, so if useful, you could consider applying it to your own datasets. 

```ruby
# installing and loading FICUS
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("https://github.com/saezlab/ficus")
library(ficus)
```

## (0) Patient cohorts 
Each of the cohorts are based on a published data set: cohort A consists of 20 random patients from the SU2C-MARK lung cancer cohort [(Ravi et al., 2023) ](https://www.nature.com/articles/s41588-023-01355-5) while cohort B consists of 20 random patients from the [TCGA-KIRC cohort](https://gdc.cancer.gov/about-data/publications/kirc_2013).  

```ruby
# retrieve TF activities and annotations of the small patient cohorts 
library(RCurl)
cohortA <- read.csv(text = getURL("https://raw.github.com/aronlindberg/dummy.csv"))
cohortB <- read.csv(text = getURL("https://raw.github.com/aronlindberg/dummy.csv"))
```

## (1) Network curation with MOON
To derive patient-specific protein networks and corresponding functional scores for protiens, we'll provide the following inputs to MOON: (1) transcription factor (TF) activities, derived from bulk transcriptomics using [CollecTRI](https://github.com/saezlab/CollecTRI), and (2) a general Omnipath PKN. TF activities have already been loaded for ```cohortA``` and ```cohortB```, which we can use with the Omnipath PKN to run MOON. While running MOON, coherence checks are also performed and a soft network reduction to create smaller networks more computationally feasible for follow-up steps but still of significant size which we leverage when creating patient subgroups. 

```ruby
# load Omnipath PKN

# setting some variables for MOON
PKN_path <- '../data/clean_omnipath_PKN.RData' # data for initial PKN
min_size_PKN = 20 # minimal size for PKN 
significant_input_threshold <- 2 # threshold to filter TF activities 
n_steps <- 6 # number of steps during network pruning
use_subset = F # if desired, can select subset of proteins for network

# thresholds for soft network reduction
primary_threshold1 = 0.5
secondary_threshold1 = 0.25

# running MOON
cosmos_inputs <- cohortA[['CosmosInputs']]
start <- 1
end <- length(cosmos_inputs)
run_MOON(cosmos_inputs, significant_input_threshold, PKN_path, n_steps, 
         use_subset, min_size_PKN, primary_threshold1, secondary_threshold1, 
         output_folder, start, end)
```

The ```run_MOON``` function saves the several files, either required for the conversion to a logic model, or for analysis purposes (**(*)** stands for the cell line or patient id):
- **allSIF.RData** : collection of all patient-specific protein networks, can be used for easy comparison and analysis of the curated networks
- **(*)_BeforeReduction_SIF_decouplerino_full.csv** : patient-specific protein network before soft network reduction 
- **(*)_SIF_decouplerino_full.csv** : final patient-specific protein network 
- **(*)_ATT_decouplerino_full.csv** : functional scoring of each protein in the patient-specific protein network 
- **(*).RData** : R object with patient-specific data 
- **00_PKNsizes.RData** : overview of all protein network sizes (number of edges)

## (2) Network clustering 


## (3) Converting MOON outputs to CellNOpt inputs 

## (4) Model optimization 

## (5) In silico knockout screenings 

## (6) Post-session remarks 



