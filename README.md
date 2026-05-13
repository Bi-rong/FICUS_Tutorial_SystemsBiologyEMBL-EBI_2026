# FICUS_Tutorial_SystemsBiologyEMBL-EBI_2026
_Wang et al. 2026 Manuscript in preparation_ <br>
<br>
Repository with course material for the session "From omics data to networks to mechanistic models" of the Systems Biology course EMBL-EBI 2026. 

This repository contains the assignments (incl. data sets) and slides used for the session. The assignments provide examples of how FICUS and the integrated methods can be used to go from omics data to dynamic logic models. The main github of FICUS (link will be added if repository goes public) includes extended vignettes, with more details the different variables that you can tune in the framework, so please refer to them if you're interested. 

## Installation 
Below is a step-by-step description of how to get the tutorial up and running. The first part is specific to the VM environment that you can use during the course. If you're planning to run this tutorial on your own device, you most likely can skip the first step. 

### Step 1 : Activating environment 
In the VM, the first step is to activate the correct environment so you can start RStudio. In the terminal, you can type the following to get started in RStudio:
```ruby
activate conda environment
rstudio
```

### Step 2 : Downloading repositories 
Download this repository as a ZIP file. There are two ways of downloading FICUS: either via Github directly or through a local ZIP file. The ZIP file can be retrieved from: https://drive.google.com/file/d/1B12r4hdZBSy_6WuvZO-_M79YkqJMgCya/view?usp=sharing.

Extract both ZIP files into the same folder. 

### Step 3 : Installing FICUS 
Open the [Installation.R](https://github.com/Bi-rong/FICUS_Tutorial_SystemsBiologyEMBL-EBI_2026/blob/main/script/00_Installation.R) file in RStudio, and run the script to install all packages. The VM already comes with some packages installed, so these will be skipped when running this script. If you're running this on your own device, this script helps in downloading the necessary libraries to run FICUS. 

To give a brief summary, Installation.R installs three packages from Bioconductor (CellNOptR, CNORode and MEIGOR), and one from Github (cosmosR). 
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
```

All other dependencies should be installed automatically when installing FICUS. As mentioned before, you can install FICUS through the downloaded, extracted FICUS repository:
```ruby
devtools::install('../../ficus-main')
```
Or download directly from the Github repository:
```ruby
devtools::install_github("https://github.com/saezlab/ficus")
```

You can load libraries to test if everything has been installed correctly:
```ruby
library(CellNOptR)
library(CNORode)
library(ficus)
library(reshape2) # check if dependencies are installed correctly
```

After this step, you should be ready to go and get started with the practical! Please refer to the [main markdown file](https://github.com/Bi-rong/FICUS_Tutorial_SystemsBiologyEMBL-EBI_2026/blob/main/script/FICUS%20Practical%20SystemsBiology2026.md) to start going through the example code and assignments.  
