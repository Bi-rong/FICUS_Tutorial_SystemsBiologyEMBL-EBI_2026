# FICUS Practical SystemsBiology2026
This file includes both explanation and the code for the practical. Separate R scripts are also available for each section containing only the code. 

```ruby
# installing and loading FICUS
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("https://github.com/saezlab/ficus")
library(ficus)
```

## (0) Patient cohorts 
This practical will use two small patient cohorts, one for demonstration purposes and one for you to use during the assignments. Each of the cohorts are based on a published data set: cohort A consists of 20 random patients from the SU2C-MARK lung cancer cohort [(Ravi et al., 2023) ](https://www.nature.com/articles/s41588-023-01355-5) while cohort B consists of 20 random patients from the [TCGA-KIRC cohort](https://gdc.cancer.gov/about-data/publications/kirc_2013).  

```ruby
# retrieve TF activities and annotations of the small patient cohorts 
library(RCurl)
cohortA <- read.csv(text = getURL("https://raw.github.com/aronlindberg/dummy.csv"))
cohortB <- read.csv(text = getURL("https://raw.github.com/aronlindberg/dummy.csv"))
```

## (1) Network curation with MOON
From a bulk transcriptomics set, it is possible to derive transcription factor activities using methods such as [CollecTRI](https://github.com/saezlab/CollecTRI). 


## (2) Network clustering 


## (3) Converting MOON outputs to CellNOpt inputs 

## (4) Model optimization 

## (5) In silico knockout screenings 

## (6) Post-session remarks 



