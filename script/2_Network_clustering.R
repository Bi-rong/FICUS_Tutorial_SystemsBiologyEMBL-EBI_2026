# loading FICUS ---------------------------------------------------------------
path <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

library(ficus)

# load data from previous step -------------------------------------------------
load('../output/MOON_cohortA/allSIF.RData')
nPatients <- length(unique(allSIF$patient))
all_edges <- unique(allSIF[,c('source', 'target')])

# perform clustering ----------------------------------------------------------
# this section takes around (13:33)
cluster_metric = 'Jaccard'
output_folder = '../output/clustering_cohortA/'
nClusters = 3

res <- cluster_networks(nPatients, all_edges, allSIF,
                        output_folder, nClusters, 
                        metric = cluster_metric) 
save(res, file=paste(output_folder, 'Heatmap_', cluster_metric, '_data.RData', sep=''))
