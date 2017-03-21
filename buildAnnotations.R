# --------------------------------------------------------------
# Build an annotations data frame from a tab-indexed, gzipped GFF3 file accessible by HTTP
# --------------------------------------------------------------

library(Rsamtools)
library(stringi)

# --------------------------------------------------------------

extract.gff.attribute <- function(text, s) {
  i <- regexpr(s, text)
  if (length(i) == 0) return("")

  isc <- regexpr(";", text)
  isc <- isc[isc > i]
  if (length(isc) == 0) return(stri_sub(text, i + nchar(s) + 1))

  stri_sub(text, i + nchar(s) + 1, isc[1] - 1)
}

# --------------------------------------------------------------

pp.special <- c(
  "%20", "%21", "%22", "%23", "%24", "%25", "%26", "%27",
  "%28", "%29", "%2A", "%2B", "%2C", "%2D", "%2E", "%2F",
  "%3A", "%3B", "%3C", "%3D", "%3E", "%3F"
)
rr.special <- c(
  " ", "!", "&quot;", "#", "$", "%", "&amp;", "'",
    # use "'" because &apos; does not work in some browsers
  "(", ")", "*", "+", ",", "-", ".", "/",
  ":", ";", "&lt;", "=", "&gt;", "?"
)
replace.special.characters <- function(text) {
  stri_replace_all_fixed(text, pp.special, rr.special, vectorize_all = FALSE)
}

# --------------------------------------------------------------

# The following assumes that column 3 (type) = "gene".
# If not, include and filter on that column.
gff.cols <- c(1, 4, 5, 7, 9) # see names(df.annot) below

# TODO: generalize for other possible formats of column 9 (attributes)
build.annotations <- function(filename, chrLengths) {
  t0 <- proc.time()[3]
  cat(paste("Constructing", key, "annotations ... "))

  chrs <- paste0("chr", 1:length(chrLengths))
  pp <- GRanges(seqnames = chrs, ranges = IRanges(start = 1, end = chrLengths))

  df.annot <- unlist(scanTabix(filename, param = pp), use.names = FALSE)
  df.annot <- stri_split_fixed(df.annot, "\t")
  df.annot <- lapply(df.annot, FUN = function(x) x[gff.cols])
  df.annot <- as.data.frame(do.call(rbind, df.annot), stringsAsFactors = FALSE)
#  df.annot <- df.annot[df.annot[, 2] == "gene", -2] # uncomment if type column is not already filtered
  names(df.annot) <- c("chromosome", "transcript_start", "transcript_end", "strand", "attributes")
  df.annot$chromosome <- sapply(df.annot$chromosome, FUN = function(chr) stri_sub(chr, 4))
  df.annot$transcript_start <- as.integer(df.annot$transcript_start)
  df.annot$transcript_end <- as.integer(df.annot$transcript_end)
  df.annot$id <- sapply(df.annot$attributes, FUN = function(s) extract.gff.attribute(s, "ID"))
#  df.annot$name <- sapply(df.annot$attributes, FUN = function(s) extract.gff.attribute(s, "Name"))
  df.annot$name <- paste0("medtr.", df.annot$id)
  df.annot$description <- sapply(df.annot$attributes, FUN = function(s) replace.special.characters(extract.gff.attribute(s, "Note")))
  df.annot$attributes <- NULL # remove attributes column

  cat(sprintf("Done. (%2.1f seconds)\n", proc.time()[3] - t0))
  df.annot
}

# Tests
# mt.chromosome.lengths <- c(52991155, 45729672, 55515152, 56582383, 43630510, 35275713, 49172423, 45569985)
# df.annot <- build.annotations(
#   "http://data.comparative-legumes.org/genomes/medtr/Mt4.0v1_genes_20130731_1800.just_genes.gff3.gz",
#   mt.chromosome.lengths
# )
# write.csv(df.annot, file = "Medicago_truncatula_annotations.csv", quote = 7, row.names = FALSE)

# --------------------------------------------------------------
