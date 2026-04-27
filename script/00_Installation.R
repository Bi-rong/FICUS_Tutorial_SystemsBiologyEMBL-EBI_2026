path <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)
install_local('../../Downloads/ficus.zip', subdir = NULL, quiet = FALSE)