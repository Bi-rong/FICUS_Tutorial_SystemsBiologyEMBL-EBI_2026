# loading packages -------------------------------------------------------------
path <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

library(CellNOptR)
library(CNORode)
library(ficus)


# optimization parameters -----------------------------------------------------
folder <- 'MOON_cohortA/'
cluster_metric <- 'Jaccard'

# load overview CSV with all hyperparameters used to generate training data 
overview <- read.csv(paste('../output/MOON_cohortA/cluster_annots.csv', sep=''), check.names = F)[1,]

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


