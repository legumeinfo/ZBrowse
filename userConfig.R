# User-configurable global variables
userConfig <- list()

# Default values for communication with GCV
# GCV base URL (may override in organism files)
userConfig$default_gcv_url <- "https://gcv-microservices.lis.ncgr.org/lis"
# for micro-synteny
userConfig$bcName <- "GCV"
userConfig$neighbors <- 20
userConfig$matched <- 4
userConfig$intermediate <- 5
# for macro-synteny
userConfig$macroMatched <- 20
userConfig$macroIntermediate <- 10
userConfig$macroMask <- 10
userConfig$macroDistance <- "levenshtein" # for Jaccard distance use something like "jaccard:1:false"

# TODO: move default colors here?
