# User-configurable global variables
userConfig <- list()

# GCV microservices base URL (may override in organism files)
userConfig$default_gcv_microservices_url <- "https://gcv-microservices.lis.ncgr.org/lis"
# GCV client URL (for Broadcast Channel communication)
# if NULL, it will default to <url_protocol>//<url_hostname>:<url_port>/gcv2 from session$clientData fields,
# or you may explicitly set it here
userConfig$gcv_client_url <- NULL
userConfig$bcName <- "GCV"

# for micro-synteny
userConfig$neighbors <- 20
userConfig$matched <- 4
userConfig$intermediate <- 5
# for macro-synteny
userConfig$macroMatched <- 20
userConfig$macroIntermediate <- 10
userConfig$macroMask <- 10
userConfig$macroDistance <- "levenshtein" # for Jaccard distance use something like "jaccard:1:false"

# TODO: move default colors here?
