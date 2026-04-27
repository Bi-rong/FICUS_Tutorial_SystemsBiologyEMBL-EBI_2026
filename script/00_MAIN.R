# installing and loading FICUS ------------------------------------------------
path <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

#if (!requireNamespace("devtools", quietly = TRUE)){
#  install.packages("devtools")}
#devtools::install_github("https://github.com/saezlab/ficus")
library(ficus)
library(reshape2)
library(CellNOptR)
library(CNORode)

## STEP 1 : NETWORK CURATION WITH MOON ----------------------------------------
# retrieve the small patient cohorts 
load('../data/cohortA.RData')
load('../data/cohortB.RData')

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


## STEP 2 : NETWORK CLUSTERING ------------------------------------------------
# load data from previous step 
load('../output/MOON/allSIF.RData')
nPatients <- length(unique(allSIF$patient))
all_edges <- unique(allSIF[,c('source', 'target')])

# perform clustering 
# this section takes a few minutes
cluster_metric = 'Jaccard'
output_folder = '../output/clustering/'
nClusters = 3

res <- cluster_networks(nPatients, all_edges, allSIF,
                        output_folder, nClusters, 
                        metric = cluster_metric) 
save(res, file=paste(output_folder, 'Heatmap_', cluster_metric, '_data.RData', sep=''))





## STEP 3 : CONVERTING MOON OUTPUTS INTO CELLNOPT INPUTS ----------------------
# load data from previous step 
load('../output/MOON/allSIF.RData') # load RData with all SIF files 

# clusters
cluster_metric = 'Jaccard' 
nClusters = 3 
clusters = read.csv(
  paste('../output/clustering/Automatic_Clusters_', cluster_metric, '_nClusters=', 
        nClusters, '.csv', sep=''))


# part 1 : preparing combined PKN 
# (1) Define some variables 
output_folder = '../output/MOON/' # folder where you also saved MOON outputs

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
nIndex = 1

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




# part 2 : preparing training data 
# (1) define some variables
# input preparation variables 
scale_threshold = 2 # threshold for capping values 
retain_perc <- 1 # fraction of proteins to optimize
n_proteins <- 1000 # choose large number if you want to include all proteins
k_folds <- 5 # number of folds in crossvalidation 

# needed to ensure getting the same stimuli 
set.seed(42)

# load files 
filename = paste(output_folder, '/nClust=',
                 nClusters, '_nIndex=', nIndex, '.RData',
                 sep='')
load(filename)
SIF <- solution_network$SIF

# get all MOON activity scores for each patient in the cluster 
temp_output <- paste(output_folder, '/', sep='')
output <- get_ATT(sample_names, SIF, output_folder, nClusters, nIndex, 
                  scale_threshold, temp_output, normalize=T)
allATT <- output[['allATT']]
allInputs <- output[['allInputs']]

# (2) prepare inputs for CellNOpt - CNORode 
temp_output <- output_folder
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



## STEP 5 : IN SILICO KNOCKOUT SCREENINGS -------------------------------------
folder <- 'MOON'
clusters <- read.csv(paste('../output/clustering/Automatic_Clusters_Jaccard_nClusters=3.csv', sep=''),
                     check.names = F)

metric <- 'Jaccard'
nPatients <- 12 
one_ind <- 1

# load model inputs
pknmodel = readSIF(paste('../output/', folder, '/SIF_ODE_nPatients=', 
                         nPatients, '_nIndex=', one_ind, '.SIF', sep=''))
MIDAS_filename <- paste('../output/', folder, 
                        '/MIDAS_ODE_nPatients=', nPatients, '_nIndex=', one_ind, '.csv', sep='')
cnolist = CNOlist(MIDAS_filename, verbose = F)
MIDAS <- read.csv(MIDAS_filename, check.names=F)
patient_names <- colnames(MIDAS)[grepl('CellLine', colnames(MIDAS))]
patient_names <- gsub('TR:', '', patient_names)
patient_names <- gsub(':CellLine', '', patient_names)
patient_names <- data.frame('Condition'=1:length(patient_names),
                            'Patient'=patient_names)

# get optimized model(s)
load('../data/trained_models/cohortA_cluster1.RData')

all_simulations <- data.frame()
for (current_idx in 1:length(all_res)){
  print(current_idx)
  
  res <- all_res[[current_idx]]
  model <- preprocessing(cnolist, pknmodel, compression = FALSE, expansion = FALSE,
                         verbose=F)
  
  simulated <- plotLBodeFitness(cnolist = cnolist, 
                                model = model, 
                                transfer_function=paramsSSm$transfer_function,
                                ode_parameters=res)
  
  temp <- simulated[[2]]
  colnames(temp) <- colnames(data.frame(cnolist@signals$`10`))
  simulation_df <- melt(temp)
  colnames(simulation_df) <- c('Condition', 'Protein', 'Activity')
  simulation_df <- merge(simulation_df, patient_names, by='Condition')
  simulation_df$run <- current_idx  
  
  all_simulations <- rbind(all_simulations, simulation_df)
}

all_changed = list()
proteins <- unique(all_simulations$Protein)[1:5] # run for 5 proteins as example
knockout_simulations <- data.frame()
changed_idx = 1 
for (protein in proteins){ 
  for (current_idx in 1:length(all_res)){
    optimized_parameters <- all_res[[current_idx]]
    
    print(paste(protein, current_idx))
    
    temp <- simulate_knockout(protein, optimized_parameters, 
                              cnolist, model, paramsSSm, plotParams, 
                              patient_names)
    one_knockout <- temp[['simulation']]
    all_changed[[changed_idx]] <- temp[['changed']]
    changed_idx = changed_idx + 1 
    
    one_knockout$run <- current_idx
    one_knockout$knockout <- protein
    
    knockout_simulations <- rbind(knockout_simulations, one_knockout)
  }
}

# replace the "parent_of" nodes so it's better for visualization 
old_names <- c()
new_names <- c()
idx = 1 
for (one_protein in unique(all_simulations$Protein)){
  if (grepl('parent_of', one_protein)){
    old_names <- c(old_names, one_protein) 
    new_names <- c(new_names, paste('parent', idx, sep=''))
    
    idx = idx + 1 
  }
}
knockout_simulations <- replace_values(old_names, new_names, 'knockout', 
                                       knockout_simulations)

# save results results of complete in silico knockout screening 
save(knockout_simulations, all_simulations, old_names, new_names, 
     file = paste('../output/', folder, '/knockout_results.RData', sep=''))
