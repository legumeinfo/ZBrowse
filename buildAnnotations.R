# --------------------------------------------------------------
# Build an annotations data frame from a tab-indexed, gzipped GFF3 file accessible by HTTP
# --------------------------------------------------------------
library(Rsamtools)
library(stringi)

source("common.R")
# --------------------------------------------------------------

extract.gff.attribute <- function(text, s) {
  stri_match_first(text, regex = sprintf("%s=(.*?)(;|$)", s))[, 2]
}

# GFF column indices for ("chromosome", "type", "transcript_start", "transcript_end", "strand", "attributes")
gff.cols <- c(1, 3, 4, 5, 7, 9)

# --------------------------------------------------------------

build.annotations <- function(key, filename, chrLengths, annotChrFormat) {
  t0 <- proc.time()[3]
  cat(paste("Constructing", key, "annotations ... "))

  chrs <- sprintf(annotChrFormat, 1:length(chrLengths))
  pp <- GRanges(seqnames = chrs, ranges = IRanges(start = 1, end = chrLengths))

  if (startsWith(filename, "https:")) {
    # since R cannot handle https directly
    gff <- stri_match(annotChrFormat, regex = "(.+)\\.")[, 2]
    gff <- paste0("./annotations/", gff, ".gff.gz")
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
  df.annot$chromosome <- trailingChromosomeName(df.annot$chromosome, organism = key)
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
  df.annot$description <- sapply(df.annot$attributes, FUN = function(s) URLdecode(extract.gff.attribute(s, "Note")))
  df.annot$attributes <- NULL # remove attributes column

  cat(sprintf("Done. (%2.1f seconds)\n", proc.time()[3] - t0))
  df.annot
}

# --------------------------------------------------------------
