# --------------------------------------------------------------
# Build an annotations data frame from a tab-indexed, gzipped GFF3 file accessible by HTTP
# --------------------------------------------------------------

library(Rsamtools)
library(stringi)

# --------------------------------------------------------------

extract.gff.attribute <- function(text, s) {
  stri_match_first(text, regex = sprintf("%s=(.*?)(;|$)", s))[, 2]
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

# GFF column indices for ("chromosome", "type", "transcript_start", "transcript_end", "strand", "attributes")
gff.cols <- c(1, 3, 4, 5, 7, 9)

build.annotations <- function(key, filename, chrLengths, chrPrefix) {
  t0 <- proc.time()[3]
  cat(paste("Constructing", key, "annotations ... "))

  # We expect that the number part of the chromosome name will have enough leading zeros to allow alphabetic sorting.
  # This is only necessary for this parsing step, afterward we will refer to chromosomes by number
  # unless otherwise specified in the organism file (as for cowpea and soybean).
  num.chromosomes <- length(chrLengths)
  num.digits <- 1 + floor(log10(num.chromosomes))
  chrs <- sprintf(sprintf("%s%%0%dd", chrPrefix, num.digits), 1:num.chromosomes)
  pp <- GRanges(seqnames = chrs, ranges = IRanges(start = 1, end = chrLengths))

  if (startsWith(filename, "https:")) {
    # since R cannot handle https directly
    annotations.dir <- "./annotations/"
    gff <- paste0(annotations.dir, chrPrefix, ".gz")
    #allow this to be cached
    if (! file.exists(gff)) { 
        download.file(filename, gff, method = "wget", quiet = TRUE)
        download.file(paste0(filename, ".tbi"), paste0(gff, ".tbi"), method = "wget", quiet = TRUE)
    }
    df.annot <- unlist(scanTabix(gff, param = pp), use.names = FALSE)
  } else {
    df.annot <- unlist(scanTabix(filename, param = pp), use.names = FALSE)
  }
  df.annot <- stri_split_fixed(df.annot, "\t")
  df.annot <- lapply(df.annot, FUN = function(x) x[gff.cols])
  df.annot <- as.data.frame(do.call(rbind, df.annot), stringsAsFactors = FALSE)
  df.annot <- df.annot[df.annot[, 2] == "gene", -2] # filter and then remove the type column
  names(df.annot) <- c("chromosome", "transcript_start", "transcript_end", "strand", "attributes")
  if (key %in% c("Cowpea", "Soybean")) {
    # "Vu01" or "Gm01" instead of "1", so start two characters back
    df.annot$chromosome <- sapply(df.annot$chromosome, FUN = function(chr) stri_sub(chr, nchar(chrPrefix) - 1))
  } else {
    df.annot$chromosome <- sapply(df.annot$chromosome, FUN = function(chr) stri_sub(chr, nchar(chrPrefix) + 1))
  }
  df.annot$transcript_start <- as.integer(df.annot$transcript_start)
  df.annot$transcript_end <- as.integer(df.annot$transcript_end)
  df.annot$id <- sapply(df.annot$attributes, FUN = function(s) extract.gff.attribute(s, "ID"))
  if (key == "Medicago truncatula") {
    # it has no Name field, so construct one
    #df.annot$name <- paste0("medtr.", df.annot$id)
    # the Names in the new file have additional prefixing that miust be stripped off before the linkout service will 
    # function properly
    df.annot$name <- gsub("^medtr\\.A17_HM341\\.", "medtr.",
            sapply(df.annot$attributes, FUN = function(s) extract.gff.attribute(s, "Name")))
  } else {
    df.annot$name <- sapply(df.annot$attributes, FUN = function(s) extract.gff.attribute(s, "Name"))
  }
  df.annot$description <- sapply(df.annot$attributes, FUN = function(s) replace.special.characters(extract.gff.attribute(s, "Note")))
  df.annot$attributes <- NULL # remove attributes column

  cat(sprintf("Done. (%2.1f seconds)\n", proc.time()[3] - t0))
  df.annot
}

# --------------------------------------------------------------
