# installing and loading FICUS ------------------------------------------------
if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")
devtools::install_github("https://github.com/saezlab/ficus")
library(ficus)


# retrieve TF activities and annotations of the small patient cohorts ----------
library(RCurl)
#cohortA <- read.csv(text = getURL("https://raw.github.com/aronlindberg/dummy.csv"))
cohortB <- read.csv(text = getURL("https://raw.github.com/aronlindberg/dummy.csv"))
