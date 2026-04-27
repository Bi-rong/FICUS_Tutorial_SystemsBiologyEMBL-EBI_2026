# FICUS Practical SystemsBiology2026
Welcome to this introduction to FICUS! Next to this markdown with additional explanations for each step, there are two other R scripts important for the session:
1. **00_MAIN.R** : Combines all the code from this markdown into one script, you can run this script in parallel to the explanations below. The suggestions is to run it line-by-line in Rstudio (CTRL + ENTER) so you can follow step by step what is happening, instead of running the entire script all at once.
2. **00_ASSIGNMENTS.R** : While the previous script contains code for one example cohort, this script replaced specific lines with comments to indicate which variables and paths you can set yourself to get more hands-on experience in running FICUS. The main idea is that you can use another example cohort from this repository to run the code, or if you're interested, you can also run data from your own project. Additionally, this script also encourages you to further analyze the results, so you can take a deeper dive into what you're actually getting from FICUS. 

This practical focuses on the main FICUS functions, allowing you to go from omics data to logic models. Two small patient cohorts are provided, one for demonstration purposes and one for you to use during the assignments. Each of the cohorts are based on a published data set: cohort A consists of 20 random patients from the SU2C-MARK lung cancer cohort [(Ravi et al., 2023) ](https://www.nature.com/articles/s41588-023-01355-5) while cohort B consists of 20 random patients from the [TCGA-KIRC cohort](https://gdc.cancer.gov/about-data/publications/kirc_2013).  

As briefly mentioned above, the assignments here are more guidelines on what can be tuned in the current framework, and what you might want to consider when using the tools. They are nothing more than a place where you can try out the FICUS tools, change variables, investigate the outputs and compare across data sets. The aim is to familiarize you with the methods used in FICUS, so if useful, you could consider applying it to your own datasets. 

If you've finished going through the installation guide in the README, didn't bump into weird errors and don't have any more questions, you're ready to get started!

## (1) Network curation with MOON
To derive patient-specific protein networks and corresponding functional scores for protiens, we'll provide the following inputs to MOON: (1) transcription factor (TF) activities, derived from bulk transcriptomics using [CollecTRI](https://github.com/saezlab/CollecTRI), and (2) a general Omnipath PKN. TF activities have already been loaded for ```cohortA``` and ```cohortB```, which we can use with the Omnipath PKN to run MOON. While running MOON, coherence checks are also performed together with a soft network reduction. The goal of this reduction is to create smaller protein networks that are computationally feasible for follow-up steps but still of significant size which we leverage when creating patient subgroups.

```ruby
# Set directory and load libraries and data ----------------------------------
path <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

library(ficus)
library(reshape2)
library(CellNOptR)
library(CNORode)

# retrieve the small patient cohorts 
load('../data/cohortA.RData')
```

```ruby
## STEP 1 : NETWORK CURATION WITH MOON ----------------------------------------
# setting some variables for MOON 
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

# running MOON 
# this section takes around 5 minutes to run  
output_folder <- '../output/MOON/' 
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
# load data from previous step
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

Let's start with the PKN, which will be saved as a SIF file for the logic-ODE model in CellNOpt. A straightforward way of creating a PKN representing the subgroup is by simply aggregating all patient-specific networks (union between all networks). One important issue is however that this significantly increases the number of edges in the network. Optimizing an ODE model with thousands of edges (the number of parameters scales directly with the network size) becomes rapidly infeasible. We therefore include a strong network reduction step in the preparation of the PKN, with thresholds ```primary_threshold2``` and ```secondary_threshold2```.

In addition to the strict two-threshold network reduction, we also do some filtering based on occurrences of nodes and edges across patients. Firstly, we filter away proteins for which MOON activity scores are unknown for a large majority of the patients (```NA_threshold```). Secondly, we remove edges if they are not present across many patients (```NA_edge_threshold```). The edge filtering is important as an absence of such edge in a patient-specific MOON network might indicate that the edge was removed before during the consistency check in MOON. Such edges would be difficult to model, especially in combination with other patients.

```ruby
# load data from previous step 
load(## PATH TO "allSIF.RData") # load RData with all SIF files 

# clusters : define cluster_metric and nClusters corresponding with the chosen values during clustering
cluster_metric = 'Jaccard' 
nClusters = 3 
clusters = read.csv(
  paste(## PATH TO "paste('Automatic_Clusters_', cluster_metric, '_nClusters=', 
        nClusters, '.csv', sep='')")

# part 1 : preparing combined PKN =========================================================================
# (1) Define some variables 
output_folder = ## PATH TO MOON OUTPUTS 

# threshold for protein selection : fraction of patients that don't have NA for this protein
NA_threshold = 0.8

# threshold for edge selection : fraction of patient that need to have this edge before selection 
NA_edge_threshold = 0.70

# set thresholds for network reduction
primary_threshold2 = 3
secondary_threshold2 = 2

# select cluster to prepare for logic-ODE model 
# you can use heatmap visualization from previous step to determine which 
# cluster could be suitable 
nIndex = ## CHOOSE CLUSTER INDEX TO PREPARE FOR LOGIC-MODEL

# get dataframe with nIndex and nPatients
nPat_nInd <- data.frame()
for (i in unique(clusters$cluster)){
  nPat_nInd <- rbind(nPat_nInd, data.frame(
    'nPatients'=nrow(clusters[clusters$cluster ==i,]),
    'nIndex'=i
  ))
}
write.csv(nPat_nInd, paste(output_folder, 'cluster_annots.csv', sep=''),
          row.names = F)

# (2) get PKN for clusters 
# get union of all protein networks 
# this includes filtering SIF based on how often an edge occurs across patients 
SIF <- compile_PKN(clusters, allSIF, NA_edge_threshold)

# get average MOON scores 
MOON_scores <- aggregate_MOON_scores(clusters, SIF, nIndex, output_folder,
                                     NA_threshold)

# simplify PKN with strict two-threshold reduction
filename = paste(output_folder, 'nClust=',
                 nClusters, '_nIndex=', nIndex, '.RData',
                 sep='')
simplify_PKN(MOON_scores, SIF[[nIndex]], nIndex, primary_threshold2, secondary_threshold2,
             filename, nClusters)
```
You could say that with our PKN, we now have the **model structure** of the logic-ODE model (i.e. model parameters), but to tune the model, we still need to prepare the **training data**. An advantage from the MOON outputs and the criteria we used to filter the network is that we have functional scores for a large set of proteins in the network, which we aim to fully leverage in the conversion to dynamic models. 

For the training data, important steps include (1) defining patients as conditions, and (2) converting functional scores to CellNOpt compatible ranges. Moreover, we'd like to use a crossvalidation approach for model optimization, meaning we also need to prepare the different folds, each with different patients. The framework is specified to prepare MIDAS CSV files. If you're interested in the MIDAS format, you can refer to the [CNOdocs of CellNOpt section 2.2.1](https://saezlab.github.io/CellNOptR/6_CNODocs/). 

```ruby
# part 2 : preparing training data ==============================================
# (1) define some variables
# input preparation variables 
scale_threshold = 2 # threshold for capping values 
retain_perc <- 1 # fraction of proteins to optimize
n_proteins <- 1000 # choose large number if you want to include all proteins
k_folds <- 5 # number of folds in crossvalidation 

# needed to ensure getting the same stimuli 
set.seed(42)

# load files 
filename = ## PATH TO WHERE COMBINED PKN WAS SAVED #paste('../output/', output_folder,
            #'/nClust=', nClusters, '_nIndex=', nIndex, '.RData', sep='')
load(filename)
SIF <- solution_network$SIF

# get all MOON activity scores for each patient in the cluster 
temp_output <- ## SET OUTPUT FOLDER 
output <- get_ATT(sample_names, SIF, output_folder, nClusters, nIndex, 
                  scale_threshold, temp_output, normalize=T)
allATT <- output[['allATT']]
allInputs <- output[['allInputs']]

# (2) prepare inputs for CellNOpt - CNORode 
temp_output <- paste(output_folder, 'model_inputs', sep='')
prep_MIDAS(sample_names, SIF, nIndex, allInputs, allATT, temp_output,
           active_threshold = NA, inactive_threshold = NA,
           formalism = 'ODE', reduce_inputs = F)

# split into train-validation with crossvalidation 
MIDAS_filename <- paste('ODE_nPatients=', length(sample_names), 
                        '_nIndex=', nIndex,  sep='')
prep_crossvalidation(sample_names, nIndex, temp_output, MIDAS_filename, k_folds, 
                     retain_perc = retain_perc, 
                     n_proteins = n_proteins, formalism = 'ODE')

print(paste('::: Final PKN size:', nrow(SIF)))
print(paste('::: CNORode inputs are ready!', sep=''))
```
Similar to the previous steps, we got several output files from this step. The most important ones for the model optimization (**(XX)** is specific to the cluster):
- **MIDAS_ODE_nPatients=(XX)_nIndex=(XX).csv** : MIDAS file with all patients, this file is used to create the different folds. The MIDAS files corresponding with the different folds (train, validation) have a similar naming convention, just with an additional indication whether it's train or validation, and which fold they correspond to.
- **SIF_ODE_nPatients=(XX)_nIndex=(XX).SIF** : SIF file containing the combined PKN of all patients in the subgroup, can be used directly for CellNOpt

### Assignment
You can use this opportunity to inspect the SIF and MIDAS files, or prepare model inputs with the other data set. You can also play around with the variables that affect how strongly you're filtering the network: for instance, how do these influence the number of fully connected components within the network? 

## (4) Model optimization 
With all inputs ready, it's time to optimize the model and get a dynamic logic-ODE model. How to run CellNOpt with the crossvalidation folds is provided below. Running the model optimization can take however quite some time, which will not be feasible for this session. An example subgroup was therefore used to retrieve optimized parameters, which are provided in this repository (data) for both cohorts. These outputs can be used further for follow-up analyses and simulations, as illustrated in the next section. 

```ruby
# loading packages -------------------------------------------------------------
path <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

library(CellNOptR)
library(CNORode)
library(ficus)

# optimization parameters -----------------------------------------------------
folder <- 'MOON/model_inputs'
cluster_metric <- 'Jaccard'

# load overview CSV with all hyperparameters used to generate training data 
overview <- read.csv(paste('../output/MOON/cluster_annots.csv', sep=''), check.names = F)[1,]

maxTime <- 500
maxNumSteps <- 1e5
startRuns <- 1  
nRuns <- 1
k_folds <- 5
lam = 0

# starting with which cluster and hyperparameters
startRow = 1
endRow = nrow(overview)

# for plotting of optimized data
plotParams=list(margin=0.1, width=20, height=15, cmap_scale=1, cex=1, ymin=0,
                maxrow=50)

# set optimization parameters 
paramsSSm = init_pars_ssm()

# optimization ----------------------------------------------------------------
for (row in startRow:endRow){
  all_res <- list()
  
  current_row <- overview[row,]
  nIndex <- current_row$nIndex
  nPat <- current_row$nPatients
  retain_percentage <- current_row$retain_percentage
  nClust <- current_row$nClusters
  print('-----------')
  
  general_idx <- 1
  
  for (k in 1:k_folds){
    # load model inputs 
    pknmodel = readSIF(paste('../output/', folder, '/SIF_ODE_nPatients=', 
                             nPat, '_nIndex=', nIndex, '.SIF', sep=''))
    
    MIDAS_train_filename <- paste('../output/', folder,  
                                  '/MIDAS_train_fold=', k, 
                                  '_ODE_nPatients=', nPat, '_nIndex=', nIndex, '.csv', sep='')
    
    cnolist = CNOlist(MIDAS_train_filename, verbose = F)
    model <- preprocessing(cnolist, pknmodel, compression = F, expansion = F)
    
    # set initial parameters 
    initial_parameters=createLBodeContPars(model, LB_n = 1, LB_k = 0,
                                           LB_tau = 0, UB_n = 4, UB_k = 1, UB_tau = 1, default_n = 3,
                                           default_k = 0.5, default_tau = 0.01, opt_n = FALSE, opt_k = TRUE,
                                           opt_tau = TRUE, random = TRUE)
    
    for (iter in startRuns:nRuns){
      print(paste('Metric=', cluster_metric, ', nClusters=', nClust, ', nIndex=', nIndex, 
                  ', nPatients=', nPat, ', retain=', retain_percentage, ', Iteration=', iter, ', Fold=', k, sep=''))
      print(paste('Optimizing', nrow(cnolist@cues), 'patients and',ncol(cnolist@signals$`0`), 'proteins'))
      
      # optimization 
      # plot the simulated data using the initial parameters guess
      simulated_data_initial_parameters=plotLBodeFitness(cnolist = cnolist, 
                                                         model = model, 
                                                         transfer_function=paramsSSm$transfer_function,
                                                         ode_parameters=initial_parameters)
      
      optimized_parameters=parEstimationLBode(cnolist, model, method="essm", 
                                              ode_parameters=initial_parameters, paramsSSm=paramsSSm)
      all_res[[general_idx]] <- optimized_parameters
      general_idx <- general_idx + 1 
      
      # evaluation (training)
      simulated_data_optimized=plotLBodeFitness(cnolist, model, 
                                                transfer_function=paramsSSm$transfer_function, 
                                                ode_parameters=optimized_parameters, 
                                                plotParams = plotParams)  
      
    }
  }
  
  # save results 
  save(all_res, paramsSSm, plotParams, nClust, nPat, nIndex, retain_percentage,
       nRuns, k_folds, 
       file = paste('../output/Results_continuous_model.RData', sep=''))
}
```

(SOMETHING ABOUT OUTPUT)

## (5) In silico knockout screenings 
Now that we have the models, we can think about how to leverage the dynamic logic-ODE models for analyses and simulations. An advantage of the logic-ODE is that each protein and interaction is linked to a parameter, setting the corresponding parameter to 0 can be interpreted as removing or knocking out the protein or interaction. 


## (6) Post-session remarks 



