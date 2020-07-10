# --------------------------------------------------------------
# Build a GWAS data frame from files accessible by HTTP
# --------------------------------------------------------------

library(jsonlite)
library(stringi)
library(RCurl)

# --------------------------------------------------------------

gwas.sources <- data.frame(
  name = c("CyVerse", "DSCensor"),
  # use more basic URLs
  url = c("http://de.cyverse.org/dl/d/F61A306C-92D2-4595-8226-A195D46EBB50/FT10.gwas", "http://dev.lis.ncgr.org:50020"),
  status = FALSE,
  stringsAsFactors = FALSE
)
row.names(gwas.sources) <- gwas.sources$name
gwas.sources$status <- sapply(gwas.sources$url, url.exists)

# The following structure makes it easy to add any organism whose
# GWAS files are in the same format as those for Medicago truncatula.
# TODO: generalize for other formats, if necessary.

gwas.filenames <- list()
gwas.filenames[["Arabidopsis thaliana"]] <- c(
  "http://de.cyverse.org/dl/d/F61A306C-92D2-4595-8226-A195D46EBB50/FT10.gwas",
  "http://de.cyverse.org/dl/d/64C5BD48-CF10-4833-96B4-3C74CEC47257/FT16.gwas",
  "http://de.cyverse.org/dl/d/BBB7DCAB-87FE-4C9E-8F81-DDF8FFEFF806/FT22.gwas",
  "http://de.cyverse.org/dl/d/57E47B19-FF8C-47F6-AB32-9A111AB1F2A7/Trichome avg C.gwas",
  "http://de.cyverse.org/dl/d/9B68EAA3-D105-49B7-B2C1-E1690E8BAF23/Trichome avg JA.gwas"
  # ...
)

# Get the GWAS filenames from DSCensor
gwasBaseUrl <- "http://dev.lis.ncgr.org:50020/api/v1/nodes/labels/"
organisms.gwas <- c("Medicago truncatula") #, "Soybean", "Pigeonpea", "Cowpea")
for (o.gwas in organisms.gwas) {
  ss <- strsplit(org.Genus_species[[o.gwas]], split = " ")[[1]]
  genus <- ss[1]
  species <- ss[2]

  # If the connection is valid, return the query as a list containing a data frame,
  # otherwise return an empty list
  query <- list()
  if (gwas.sources["DSCensor", "status"]) tryCatch(
    query <- fromJSON(paste0(gwasBaseUrl, tolower(genus), ":", species, ":gwas"))
  )
  if (length(query) == 0 || length(query$data$url) == 0) {
    gwas.filenames[[paste(ss[1], ss[2])]] <- c()
  } else {
    oo <- order(sapply(query$data$url, basename))
    gwas.filenames[[paste(ss[1], ss[2])]] <- query$data$url[oo]
  }
}

gwas.traits <- list()
gwas.traits[["Arabidopsis thaliana"]] <- stri_match(basename(gwas.filenames[["Arabidopsis thaliana"]]), regex = ".*(?=.gwas)")[, 1]
for (o.gwas in organisms.gwas) {
  if (length(gwas.filenames[[o.gwas]]) == 0) {
    gwas.traits[[o.gwas]] <- c()
  } else {
    gwas.traits[[o.gwas]] <- stri_match(basename(gwas.filenames[[o.gwas]]), regex = ".*(?=_results.gwas)")[, 1]
  }
}

# TODO: standardize column names
gwas.cols <- list()
gwas.cols[["Arabidopsis thaliana"]] <- c("Chromosome", "Position", "Trait", "P.Value", "negLogP", "MAF")
for (o.gwas in organisms.gwas) {
  # TODO: this works for Medicago truncatula;
  # either make sure it works for other species, or generalize it
  gwas.cols[[o.gwas]] <- c("Chromosome", "pos", "P.value")
}

# Start with an empty data frame
init.gwas <- function(o.gwas) {
  organism <- stri_match(o.gwas, regex = ".*(?= GWAS)")[, 1]
  if (organism == "Medicago truncatula") {
    # TODO: this works for Medicago truncatula;
    # either make sure it works for other species, or generalize it
    df.gwas <- data.frame(Chromosome = "1", pos = 1L, Trait = "-", P.value = 0.1, stringsAsFactors = FALSE)
  } else if (organism == "Arabidopsis thaliana") {
    df.gwas <- data.frame(Chromosome = "1", Position = 1L, Trait = "-", P.Value = 0.1, negLogP = 1.0, MAF = 0.01, stringsAsFactors = FALSE)
  }
  df.gwas <- df.gwas[-1, ]
  df.gwas
}

load.gwas.remote <- function(organism, filename, trait) {
  filename <- gsub(" ", "%20", filename)
  if (!gwas.sources["CyVerse", "status"]) {
    print(paste("No connection to", filename))
    return()
  }

  t0 <- proc.time()[3]
  cat(paste("Loading", trait, "data ... "))
  nid <- paste0("load.remote.", stri_replace_all(trait, "_", regex = " "))
  showNotification(paste0("Loading remote trait ", trait, ". Please wait."),
    duration = NULL, id = nid, type = "message")

  filenames <- gwas.filenames[[organism]]
  traits <- gwas.traits[[organism]]
  cols <- gwas.cols[[organism]]

  if (organism == "Medicago truncatula") {
    # TODO: this works for Medicago truncatula;
    # either make sure it works for other species, or generalize it
    df.gwas <- read.table(file = url(filename, method = "libcurl"), header = TRUE, sep = "\t", quote = "\"", stringsAsFactors = FALSE)[, cols]
    df.gwas$Trait <- trait

    # Clean up: change Chromosome format from "chr1" to "1", and move Trait to the third column
    df.gwas$Chromosome <- sapply(df.gwas$Chromosome, FUN = function(chr) stri_sub(chr, 4))
    df.gwas <- df.gwas[, c("Chromosome", "pos", "Trait", "P.value")]

  } else if (organism == "Arabidopsis thaliana") {
    df.gwas <- read.table(file = url(filename, method = "libcurl"), header = TRUE, sep = ",", quote = "\"", stringsAsFactors = FALSE)[, cols]
    df.gwas$Chromosome <- as.character(df.gwas$Chromosome)
    df.gwas$Trait <- trait
  }

  cat(sprintf("Done. (%2.1f seconds)\n", proc.time()[3] - t0))
  removeNotification(nid)
  df.gwas
}

# load.gwas.local <- function(organism, filename, trait) {
#   # TODO ...
# }

# Tests
# df.gwas <- load.gwas.local("Medicago truncatula", ..., "floweringdate")
# write.csv(df.gwas, file = "Medicago_truncatula_GWAS.csv", quote = FALSE, row.names = FALSE)

# --------------------------------------------------------------
