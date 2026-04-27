# loading packages -------------------------------------------------------------
path <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

library(ficus)
library(reshape2)

folder <- 'MOON_cohortA'
clusters <- read.csv(paste('../output/clustering_cohortA/Automatic_Clusters_Jaccard_nClusters=3.csv', sep=''),
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
load('../output/Results_continuous_model.RData')

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
proteins <- unique(all_simulations$Protein)[:5] # run for 5 proteins as example
knockout_simulations <- data.frame()
for (protein in proteins){ 
  for (current_idx in 1:length(all_res)){
    optimized_parameters <- all_res[[current_idx]]
    
    print(paste(protein, current_idx))
    
    one_knockout <- simulate_knockout(protein, optimized_parameters, 
                                      cnolist, model, paramsSSm, plotParams, 
                                      patient_names)
    
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
