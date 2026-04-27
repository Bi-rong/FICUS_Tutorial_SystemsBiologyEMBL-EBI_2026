# FICUS_Tutorial_SystemsBiologyEMBL-EBI_2026
Repository with course material for the session "From omics data to networks to mechanistic models" of the Systems Biology course EMBL-EBI 2026. 

This repository contains the assignments (incl. data sets) and slides used for the session. The assignments provide examples of how FICUS and the integrated methods can be used to go from omics data to dynamic logic models. This tutorial will be more general. The main github of FICUS (ADD LINK) includes extended vignettes, which include more details also on the different variables, so please refer to them if you're interested. 

(INSTRUCTIONS ON WHAT TO DOWNlOAD ETC.)


https://drive.google.com/file/d/1B12r4hdZBSy_6WuvZO-_M79YkqJMgCya/view?usp=sharing

```ruby
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

# test if libraries were installed correctly
library(CellNOptR)
library(CNORode)
library(ficus)

```
