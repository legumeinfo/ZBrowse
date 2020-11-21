# --------------------------------------------------------------
# Build a QTL data frame from files accessible by HTTP
# --------------------------------------------------------------

# LIS Data Store information
lis.datastore.info[["Cowpea QTL"]] <- lis.datastore.info[["Cowpea GWAS"]]

read.qtl.lis.datastore <- function(fin.expt, fin.marker) {
  # first read QTL metadata
  zz <- gzcon(url(fin.expt, "r"))
  ll <- readLines(zz)
  close(zz)
  # read metadata before line beginning "#identifier"
  src.name <- src.url <- ""
  i <- 1
  while (!startsWith(ll[i], "#")) {
    ss <- strsplit(ll[i], split = "\t")[[1]]
    if (ss[1] == "Name") src.name <- ss[2]
    else if (ss[1] == "DOI") src.url <- paste0("https://doi.org/", ss[2])
    else if (ss[1] == "PMID") src.url <- paste0("https://pubmed.ncbi.nlm.nih.gov/", ss[2], "/")
    i <- i + 1
  }
  # read the rest
  df.expt <- read.csv(textConnection(ll[i:length(ll)]), header = TRUE, sep = '\t', stringsAsFactors = FALSE)
  names(df.expt) <- c("identifier", "trait")
  df.expt$publication <- paste0("<a href='", src.url, "' target=_blank>", src.name, "</a>")

  # then read QTL markers (first two columns)
  zz <- gzcon(url(fin.marker, "r"))
  ll <- readLines(zz)
  close(zz)
  df.marker <- read.csv(textConnection(ll), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  df.marker <- df.marker[, 1:2]
  names(df.marker) <- c("identifier", "marker")

  # finally, merge them
  df.qtl <- merge(df.expt, df.marker, by = "identifier")
  df.qtl
}

merge.qtl <- function(df.qtl, df.gff) {
  # df.qtl is the QTL data frame merged from expt and marker files:
  #   identifier,trait,publication,marker
  # df.gff is the processed GFF file:
  #   marker,chromosome,position

  # Merge on marker to convert to marker,identifier,trait,publication,chromosome,position
  df.1 <- merge(df.qtl, df.gff, by = "marker", sort = FALSE)
  # Determine start and end position for each QTL identifier
  df.2 <- df.1[, c("identifier", "position")]
  f3 <- function(r) c(min(r), max(r))
  df.3 <- do.call(function(...) data.frame(..., stringsAsFactors = FALSE),
    aggregate(position ~ identifier, data = df.2, FUN = f3))
  names(df.3) <- c("identifier", "start_pos", "end_pos")
  # Merge relevant columns
  df.4 <- unique(df.1[, c("identifier", "trait", "publication", "chromosome")])
  df.merged <- merge(df.4, df.3, by = "identifier", sort = FALSE)

  # Clean up column order (identifier, trait, chromosome, start_pos, end_pos, publication)
  df.merged <- df.merged[, c(1, 2, 4, 5, 6, 3)]
  df.merged
}

build.qtl.from.lis.datastore <- function(key) {
  nid <- paste0("load.", gsub(" ", ".", key))
  showNotification(paste("Loading", key, "data. Please wait."), duration = NULL, id = nid, type = "message")

  # TODO: Discover GFF and QTL files from DSCensor (instead of as below)
  df.qtl <- init.qtl(key)
  gffBaseUrl <- paste(url_dscensor, "api/v1/nodes/labels/", sep = "/")
  qtlBaseUrl <- paste(url_lis, "data/public/Vigna_unguiculata/mixed.qtl.KF1G/", sep = "/")
  qtlPrefix <- "vigun.mixed.qtl.KF1G."
  qtlFileNumbers <- c("22691139", "25620880", "26450274", "27658053", "29356213", "29674702", "30143525")
  exptFiles <- paste0(qtlBaseUrl, qtlPrefix, qtlFileNumbers, ".expt.tsv.gz")
  markerFiles <- paste0(qtlBaseUrl, qtlPrefix, qtlFileNumbers, ".marker.tsv.gz")

  # GFF (marker) files
  query.mrk <- fromJSON(paste0(gffBaseUrl, "mrk:", lis.datastore.info[[key]]$mrkFilter))
  if (length(query.mrk$data) > 0) {
    df.mrk <- query.mrk$data[endsWith(query.mrk$data$url, ".gff3.gz"), ]
    if (nrow(df.mrk) > 0) {
      for (i in 1:nrow(df.mrk)) {
        df.gff <- read.gff3.lis.datastore(df.mrk[i, "url"])
        df.gff <- scrub.gff(df.gff, lis.datastore.info[[key]])

        # Associated QTL files
        # query.qtl <- fromJSON(paste0(qtlBaseUrl, "qtl:", lis.datastore.info[[key]]$mrkFilter))
        # if (length(query.qtl$data) == 0) next

        # Read the QTL files, merge with the GFF data, and append the results
        for (j in 1:length(exptFiles)) {
          df.f <- read.qtl.lis.datastore(exptFiles[j], markerFiles[j])
          df.f2 <- merge.qtl(df.f, df.gff)
          if (nrow(df.qtl) == 0) {
            df.qtl <- df.f2
          } else {
            df.qtl <- rbind(df.qtl, df.f2)
          }
        }
      }
    }
  }
  # Add spurious column for QTL y axis position
  df.qtl$val <- 1.0
  chr2n <- list()
  for (i in 1:nrow(df.qtl)) {
    chr <- df.qtl$chromosome[i]
    if (chr %in% names(chr2n)) {
      df.qtl$val[i] <- df.qtl$val[i] + chr2n[[chr]]
      chr2n[[chr]] <- chr2n[[chr]] + 0.1
    } else {
      chr2n[[chr]] <- 0.1
    }
  }
  removeNotification(nid)
  # deduplicate the results (TODO: do we need to?)
  unique(df.qtl)
}

# Start with an empty data frame
init.qtl <- function(o.qtl) {
  # organism <- stri_match(o.qtl, regex = ".*(?= QTL)")[, 1]
  df.qtl <- data.frame(identifier = "QTL", chromosome = "1", start_pos = 1L, end_pos = 2L, trait = "-", stringsAsFactors = FALSE)
  df.qtl <- df.qtl[-1, ]
  df.qtl
}

# --------------------------------------------------------------
