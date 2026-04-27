path <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

# install CellNOptR 
if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")}
BiocManager::install("CellNOptR")

# install CNORode
BiocManager::install("CNORode")
BiocManager::install("MEIGOR")

# install cosmosR
devtools::install_github("saezlab/cosmosR")

# install ficus
if (!requireNamespace("devtools", quietly = TRUE)){
  install.packages("devtools")}
devtools::install('../../ficus-main')
#devtools::install_deps('../../ficus-main')
#devtools::install_github("https://github.com/saezlab/ficus") # install from repository

# test if libraries were installed correctly
library(CellNOptR)
library(CNORode)
library(ficus)