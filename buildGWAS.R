# --------------------------------------------------------------
# Build a GWAS data frame from files accessible by HTTP
# --------------------------------------------------------------

library(stringi)

# --------------------------------------------------------------

# The following structure makes it easy to add any organism whose
# GWAS files are in the same format as those for Medicago truncatula.
# TODO: generalize for other formats, if necessary.

gwas.filenames <- list()
gwas.filenames[["Arabidopsis thaliana GWAS"]] <- c(
  "http://de.cyverse.org/dl/d/F61A306C-92D2-4595-8226-A195D46EBB50/FT10.gwas",
  "http://de.cyverse.org/dl/d/64C5BD48-CF10-4833-96B4-3C74CEC47257/FT16.gwas",
  "http://de.cyverse.org/dl/d/BBB7DCAB-87FE-4C9E-8F81-DDF8FFEFF806/FT22.gwas",
  "http://de.cyverse.org/dl/d/57E47B19-FF8C-47F6-AB32-9A111AB1F2A7/Trichome avg C.gwas",
  "http://de.cyverse.org/dl/d/9B68EAA3-D105-49B7-B2C1-E1690E8BAF23/Trichome avg JA.gwas"
  # ...
)
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
gwas.traits[["Arabidopsis thaliana GWAS"]] <- stri_match(basename(gwas.filenames[["Arabidopsis thaliana GWAS"]]), regex = ".*(?=.gwas)")[, 1]
gwas.traits[["Medicago truncatula GWAS"]] <- stri_match(basename(gwas.filenames[["Medicago truncatula GWAS"]]), regex = ".*(?=_results.gwas)")[, 1]

# TODO: standardize column names
gwas.cols <- list()
gwas.cols[["Arabidopsis thaliana GWAS"]] <- c("Chromosome", "Position", "Trait", "P.Value", "negLogP", "MAF")
gwas.cols[["Medicago truncatula GWAS"]] <- c("Chromosome", "pos", "P.value")

build.gwas <- function(organism.gwas) {
  t0 <- proc.time()[3]
  cat(paste("Constructing", organism.gwas, "results ... "))

  filenames <- gwas.filenames[[organism.gwas]]
  traits <- gwas.traits[[organism.gwas]]
  cols <- gwas.cols[[organism.gwas]]

  if (organism.gwas == "Medicago truncatula GWAS") {
    df.gwas <- read.table(file = url(filenames[1], method = "libcurl"), header = TRUE, sep = "\t", quote = "\"", stringsAsFactors = FALSE)[, cols]
    df.gwas$Trait <- traits[1]
    for (i in 2:length(filenames)) {
      df.i <- read.table(file = url(filenames[i], method = "libcurl"), header = TRUE, sep = "\t", quote = "\"", stringsAsFactors = FALSE)[, cols]
      df.i$Trait <- traits[i]
      df.gwas <- rbind(df.gwas, df.i)
    }

    # Clean up: change Chromosome format from "chr1" to "1", and move Trait to the third column
    df.gwas$Chromosome <- sapply(df.gwas$Chromosome, FUN = function(chr) stri_sub(chr, 4))
    df.gwas <- df.gwas[, c("Chromosome", "pos", "Trait", "P.value")]

    # Start with an empty data frame
    # df.gwas <- data.frame(Chromosome = "1", pos = 1L, Trait = "-", P.Value = 0.1, stringsAsFactors = FALSE)
    # df.gwas <- df.gwas[-1, ]

  } else if (organism.gwas == "Arabidopsis thaliana GWAS") {
    # Start with an empty data frame
    df.gwas <- data.frame(Chromosome = "1", Position = 1L, Trait = "-", P.Value = 0.1, negLogP = 1.0, MAF = 0.01, stringsAsFactors = FALSE)
    df.gwas <- df.gwas[-1, ]
  }

  cat(sprintf("Done. (%2.1f seconds)\n", proc.time()[3] - t0))
  df.gwas
}

# Tests
# df.gwas <- build.gwas("Medicago truncatula GWAS")
# write.csv(df.gwas, file = "Medicago_truncatula_GWAS.csv", quote = FALSE, row.names = FALSE)

# --------------------------------------------------------------
