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
# retrieve the small patient cohorts 
load('../data/cohortA.RData')
load('../data/cohortB.RData')
```

## (1) Network curation with MOON
To derive patient-specific protein networks and corresponding functional scores for protiens, we'll provide the following inputs to MOON: (1) transcription factor (TF) activities, derived from bulk transcriptomics using [CollecTRI](https://github.com/saezlab/CollecTRI), and (2) a general Omnipath PKN. TF activities have already been loaded for ```cohortA``` and ```cohortB```, which we can use with the Omnipath PKN to run MOON. While running MOON, coherence checks are also performed and a soft network reduction to create smaller networks more computationally feasible for follow-up steps but still of significant size which we leverage when creating patient subgroups. Please note that you'll need to create and define the folder for output files. 

```ruby
# setting some variables for MOON
PKN_path <- '../data/clean_omnipath_PKN.RData' # data for initial PKN
min_size_PKN = 20 # minimal size for PKN 
significant_input_threshold <- 2 # threshold to filter TF activities 
n_steps <- 6 # number of steps during network pruning
use_subset = F # if desired, can select subset of proteins for network
output_folder <- ## ADD EXISTING OUTPUT FOLDER

# thresholds for soft network reduction
primary_threshold1 = 0.5
secondary_threshold1 = 0.25

# CollecTRI net used in MOON 
load('../data/Collectri_PROGENy_networks.RData') 
net$confidence <- NA
net <- net[,c('source', 'confidence', 'target', 'mor')]

# running MOON
# this section takes around 5 minutes to run 
start <- 1
end <- length(cosmos_inputs_A) # 11:49
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
```

The ```run_MOON``` function saves the several files, either required for the conversion to a logic model, or for analysis purposes (**(*)** stands for the cell line or patient id):
- **allSIF.RData** : collection of all patient-specific protein networks, can be used for easy comparison and analysis of the curated networks
- **(*)_BeforeReduction_SIF_decouplerino_full.csv** : patient-specific protein network before soft network reduction 
- **(*)_SIF_decouplerino_full.csv** : final patient-specific protein network 
- **(*)_ATT_decouplerino_full.csv** : functional scoring of each protein in the patient-specific protein network
- **(*)_MOONscores.png** : visualization of all patient-specific MOON scores, including the interval that was selected for further processing
- **(*).RData** : R object with patient-specific data 
- **00_PKNsizes.RData** : overview of all protein network sizes (number of edges)
- **PKNsizes.png** : visualization of PKN sizes after reduction

### Assignment 
After you've run the demo code above for cohortA, you can inspect the outputs, including functional scoring (ATT) and protein networks (SIF). The SIF files can also be loaded into a software such as Cytoscape for visualization. Note that only a soft network reduction is used, you can try changing the ```primary_threshold1``` and ```secondary_threshold1``` to get a feel for how these thresholds affect the network sizes. 

Important is that if you use higher thresholds, the effect on the network becomes relatively patient-specific. You can also test how the variation in network sizes changes based on these threshold values. If you use another data set, note that for a soft network reduction, a generalized value can be used. If you wish to do already stronger network filtering before the clustering, it might be necessary to set these thresholds more manually to retain consistency across patients. 

## (2) Network clustering 
If you're satisfied with the patient-specific networks, you can use these networks to create patient subgroups. We introduce here a simple hierarchical clustering approach based on similarity between network edges, also referred to as the 'Jaccard' clustering metric. As we saved quite a lot of intermediate results in the previous section, we can load them directly. Note that the clustering scales with the number of patients, as it requires calculation of all pairwise Jaccard indices. For 20 patients, it'll only take a few minutes at most, but this can scale up quite rapidly if you have hundreds of patients.   

```ruby
# load data from previous run 
load('../output/MOON/allSIF.RData')
nPatients <- length(unique(allSIF$patient))
all_edges <- unique(allSIF[,c('source', 'target')])

# perform clustering 
# this section takes a few minutes 
cluster_metric = 'Jaccard'
output_folder = ## ADD EXISTING OUTPUT FOLDER
nClusters = 3

res <- cluster_networks(nPatients, all_edges, allSIF,
                        output_folder, nClusters, 
                        metric = cluster_metric) 
save(res, file=paste(output_folder, 'Heatmap_', cluster_metric, '_data.RData', sep=''))

```

The ```cluster_networks``` will save three outputs (**(*)** stands for the number of clusters you've chosen for the hierarchical clustering):
- **Automatic_Clusters_Jaccard_nClusters=(*)** : overview of all patients and corresponding cluster index that was assigned to them
- **Heatmap_Jaccard.png**: visualization of pairwise Jaccard indices and corresponding clusters
- **Heatmap_Jaccard_data.RData**: includes variable ```res``` which can be used to reproduce the heatmap and in case you want to set a different number of clusters, prevents you from having to calculate all pairwise Jaccard indices if you only want to create a different number of clusters 

### Assignment 
You can run the clustering for both cohorts A and B and compare how the Jaccard indices differ, and how the clusters are formed. 

In general, any approach that clusters the patient-specific network would fit in this step. You could take a moment to reflect (if you work with networks in your research) whether you have specific topology-related or biologically-related features that could be relevent for patient stratification. Think of number of interactions per protein, the presence of a protein subset, connectivity, or if you have access to patient annotations, whether e.g. cancer subtype or tumor grade could be used for patient subgroup formation. 

## (3) Converting MOON outputs to CellNOpt inputs 
Until now, we've created patient-specific networks and subgroups. With our clusters formed, we can start preparing the group-specific logic-ODE model inputs: (1) a prior knowledge network (PKN) combining the patient-specific protein networks, and (2) training data based on the patient-specific functional scoring of corresponding proteins. 

Let's start with the PKN, which will be saved as a SIF file for the logic-ODE model in CellNOpt. A straightforward way of creating a PKN representing the subgroup is by simply aggregating all patient-specific networks (intersect between all networks). One important issue is however that this significantly increases the number of edges in the network. To optimize an ODE model with thousands of edges (the number of parameters scales directly with the network size) becomes rapidly infeasible. We therefore include a strong network reduction step in the preparation of the PKN. 



part 1 : preparing combined PKN 
part 2 : preparing MIDAS (training data)

Assignment: inspect SIF and MIDAS, try out with other data set, try out different strengths in network reduction 

## (4) Model optimization 

## (5) In silico knockout screenings 

## (6) Post-session remarks 



