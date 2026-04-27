# loading FICUS ---------------------------------------------------------------
path <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

library(ficus)

# load data from previous step -------------------------------------------------
load('../output/MOON_cohortA/allSIF.RData') # load RData with all SIF files 

# clusters
cluster_metric = 'Jaccard' 
nClusters = 3 
clusters = read.csv(
  paste('../output/clustering_cohortA/Automatic_Clusters_', cluster_metric, '_nClusters=', 
        nClusters, '.csv', sep=''))


# part 1 : preparing combined PKN ---------------------------------------------
# (1) Define some variables 
output_folder = '../output/MOON_cohortA/' # folder where you also saved MOON outputs

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




# part 2 : preparing training data -------------------------------------------
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
