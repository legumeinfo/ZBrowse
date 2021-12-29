# --------------------------------------------------------------

# set path to the data directory
# setwd(".../ZZBrowse/www/config/data")

combineGWASandQTL <- function(species) {
  fin.gwas <- sprintf("%s GWAS", species)
  fin.qtl <- sprintf("%s QTL", species)
  fin.both <- sprintf("%s GWAS-QTL", species)

  df.gwas <- read.csv(fin.gwas)
  df.gwas$identifier <- NA
  df.gwas$start_pos <- NA
  df.gwas$end_pos <- NA
  df.qtl <- read.csv(fin.qtl)
  df.qtl <- df.qtl[, c("markers", "chromosome", "center_pos", "trait_id", "trait", "val", "publication", "identifier", "start_pos", "end_pos")]
  names(df.qtl) <- names(df.gwas)
  df.both <- rbind(df.gwas, df.qtl)
  write.csv(df.both, fin.both, row.names = FALSE, na = "")
}

# uncomment to run
# combineGWASandQTL("Common Bean")
# combineGWASandQTL("Cowpea")
# combineGWASandQTL("Peanut")

# --------------------------------------------------------------
