# --------------------------------------------------------------
# Build a GWAS data frame from files accessible by HTTP
# --------------------------------------------------------------

library(stringi)

# --------------------------------------------------------------

# The following structure makes it easy to add any organism whose
# GWAS files are in the same format as those for Medicago truncatula.
# TODO: generalize for other formats, if necessary.

gwas.filenames <- list()
gwas.filenames[["Medicago truncatula GWAS"]] <- c(
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
gwas.traits[["Medicago truncatula GWAS"]] <- c("floweringdate", "height", "noda", "nodb", "occupancyA", "occupancyB", "totalnod", "trichomes")

gwas.cols <- list()
gwas.cols[["Medicago truncatula GWAS"]] <- c("Chromosome", "pos", "P.value")

build.gwas <- function(organism.gwas) {
  t0 <- proc.time()[3]
  cat(paste("Constructing", organism.gwas, "results ... "))

  filenames <- gwas.filenames[[organism.gwas]]
  traits <- gwas.traits[[organism.gwas]]
  cols <- gwas.cols[[organism.gwas]]

  df.gwas <- read.table(file = url(filenames[1], method = "libcurl"), header = TRUE, sep = "\t", quote = "\"", stringsAsFactors = FALSE)[, cols]
  df.gwas$Trait <- traits[1]
  for (i in 2:length(filenames)) {
    df.i <- read.table(file = url(filenames[i], method = "libcurl"), header = TRUE, sep = "\t", quote = "\"", stringsAsFactors = FALSE)[, cols]
    df.i$Trait <- traits[i]
    df.gwas <- rbind(df.gwas, df.i)
  }

  # Clean up: convert format like "chr1	940235	5.73345464453568e-06" to "1 940235  5.241584"
  df.gwas$Chromosome <- sapply(df.gwas$Chromosome, FUN = function(chr) stri_sub(chr, 4))
  df.gwas$negLogP <- -log10(df.gwas$P.value)
  names(df.gwas)[2] <- "bp"
  df.gwas <- df.gwas[, -3] # remove P.value

  cat(sprintf("Done. (%2.1f seconds)\n", proc.time()[3] - t0))
  df.gwas
}

# Tests
# df.gwas <- build.gwas("Medicago truncatula GWAS")
# write.csv(df.gwas, file = "Medicago_truncatula_GWAS.csv", quote = FALSE, row.names = FALSE)

# --------------------------------------------------------------
