# --------------------------------------------------------------
# Build an annotations data frame from a tab-indexed, gzipped GFF3 file accessible by HTTP
# --------------------------------------------------------------
source("common.R")
# --------------------------------------------------------------

build.annotations <- function(key, filename, chrLengths, annotChrFormat) {
  t0 <- proc.time()[3]
  cat(paste("Constructing", key, "annotations ... "))

  if (startsWith(filename, "https:")) {
    # since R cannot handle https directly
    gff <- stri_match(annotChrFormat, regex = "(.+)\\.")[, 2]
    gff <- paste0("./annotations/", gff, ".gff.gz")
    #allow this to be cached
    if (! file.exists(gff)) { 
        download.file(filename, gff, method = "wget", quiet = TRUE)
        download.file(paste0(filename, ".tbi"), paste0(gff, ".tbi"), method = "wget", quiet = TRUE)
    }
    df.annot <- readGFF(gff, columns = c("seqid", "start", "end", "strand"), tags = c("ID", "Name", "Note"), filter = list(type = "gene"))
  } else {
    df.annot <- readGFF(filename, columns = c("seqid", "start", "end", "strand"), tags = c("ID", "Name", "Note"), filter = list(type = "gene"))
  }
  names(df.annot) <- c("chromosome", "transcript_start", "transcript_end", "strand", "id", "name", "description")
  df.annot$chromosome <- trailingChromosomeName(df.annot$chromosome, organism = key)
  # Add prefixes ("gensp.") to the Name field if necessary
  # (only pigeonpea and soybean already have them)
  gensp.dot <- paste0(tolower(org.Gensp[[key]]), ".")
  if (!startsWith(df.annot$name[1], gensp.dot)) {
    df.annot$name <- paste0(gensp.dot, df.annot$name)
  }
  # Convert the Note/description column from a list to a URL-decoded character vector,
  # replacing missing values with a blank string
  df.annot$description <- sapply(df.annot$description, function(desc) {
    ifelse(identical(desc, character(0)), "", desc)
  })
  df.annot$description <- sapply(df.annot$description, URLdecode, USE.NAMES = FALSE)
  # Finally, convert back from S4Vectors DataFrame to a basic data frame
  # to prevent crashing on genomic linkages
  df.annot <- as.data.frame(df.annot)

  cat(sprintf("Done. (%2.1f seconds)\n", proc.time()[3] - t0))
  df.annot
}

# --------------------------------------------------------------
