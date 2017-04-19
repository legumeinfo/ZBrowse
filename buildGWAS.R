# --------------------------------------------------------------
# Build a GWAS data frame from files accessible by HTTP
# --------------------------------------------------------------

library(stringi)

# --------------------------------------------------------------

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
gwas.filenames[["Medicago truncatula"]] <- c(
  "http://de.cyverse.org/dl/d/8F20C8BF-BEEC-4801-BCFF-41534832958B/floweringdate_results.gwas",
  "http://de.cyverse.org/dl/d/DAAF68FD-3A80-4E5D-83D5-2245C7E0D747/height_results.gwas",
  "http://de.cyverse.org/dl/d/03FCF51A-8501-487F-A12E-DF4F8CDB0446/noda_results.gwas",
  "http://de.cyverse.org/dl/d/D7082A0C-5945-47DB-B5B0-3301C78BB12D/nodb_results.gwas",
  "http://de.cyverse.org/dl/d/2C0AC931-4C07-47F9-B236-DDB581497979/occupancyA_results.gwas",
  "http://de.cyverse.org/dl/d/EE69C622-4432-479C-82CF-66BF33C1C38C/occupancyB_results.gwas",
  "http://de.cyverse.org/dl/d/F2485FCD-D5D3-4F55-B693-F00551B4FF47/totalnod_results.gwas",
  "http://de.cyverse.org/dl/d/6030E8B2-4BD9-4393-87D6-6E67D60DA7E8/trichomes_results.gwas"
)

gwas.traits <- list()
gwas.traits[["Arabidopsis thaliana"]] <- stri_match(basename(gwas.filenames[["Arabidopsis thaliana"]]), regex = ".*(?=.gwas)")[, 1]
gwas.traits[["Medicago truncatula"]] <- stri_match(basename(gwas.filenames[["Medicago truncatula"]]), regex = ".*(?=_results.gwas)")[, 1]

# TODO: standardize column names
gwas.cols <- list()
gwas.cols[["Arabidopsis thaliana"]] <- c("Chromosome", "Position", "Trait", "P.Value", "negLogP", "MAF")
gwas.cols[["Medicago truncatula"]] <- c("Chromosome", "pos", "P.value")

# Start with an empty data frame
init.gwas <- function(organism.gwas) {
  organism <- stri_match(organism.gwas, regex = ".*(?= GWAS)")[, 1]
  if (organism == "Medicago truncatula") {
    df.gwas <- data.frame(Chromosome = "1", pos = 1L, Trait = "-", P.value = 0.1, stringsAsFactors = FALSE)
  } else if (organism == "Arabidopsis thaliana") {
    df.gwas <- data.frame(Chromosome = "1", Position = 1L, Trait = "-", P.Value = 0.1, negLogP = 1.0, MAF = 0.01, stringsAsFactors = FALSE)
  }
  df.gwas <- df.gwas[-1, ]
  df.gwas
}

load.gwas.remote <- function(organism, filename, trait) {
  t0 <- proc.time()[3]
  cat(paste("Loading", trait, "data ... "))

  filenames <- gwas.filenames[[organism]]
  traits <- gwas.traits[[organism]]
  cols <- gwas.cols[[organism]]

  if (organism == "Medicago truncatula") {
    df.gwas <- read.table(file = url(filename, method = "libcurl"), header = TRUE, sep = "\t", quote = "\"", stringsAsFactors = FALSE)[, cols]
    df.gwas$Trait <- trait

    # Clean up: change Chromosome format from "chr1" to "1", and move Trait to the third column
    df.gwas$Chromosome <- sapply(df.gwas$Chromosome, FUN = function(chr) stri_sub(chr, 4))
    df.gwas <- df.gwas[, c("Chromosome", "pos", "Trait", "P.value")]

  } else if (organism == "Arabidopsis thaliana") {
    df.gwas <- read.table(file = url(filename, method = "libcurl"), header = TRUE, sep = ",", quote = "\"", stringsAsFactors = FALSE)[, cols]
    df.gwas$Trait <- trait
  }

  cat(sprintf("Done. (%2.1f seconds)\n", proc.time()[3] - t0))
  df.gwas
}

# load.gwas.local <- function(organism, filename, trait) {
#   # TODO ...
# }

# Tests
# df.gwas <- load.gwas.local("Medicago truncatula", ..., "floweringdate")
# write.csv(df.gwas, file = "Medicago_truncatula_GWAS.csv", quote = FALSE, row.names = FALSE)

# --------------------------------------------------------------
