# --------------------------------------------------------------
# Build a QTL data frame from files accessible by HTTP
# --------------------------------------------------------------

# LIS Data Store information
lis.datastore.chrRegex[["Common Bean QTL"]] <- lis.datastore.chrRegex[["Common Bean GWAS"]]
lis.datastore.chrRegex[["Cowpea QTL"]] <- lis.datastore.chrRegex[["Cowpea GWAS"]]
lis.datastore.chrRegex[["Peanut QTL"]] <- lis.datastore.chrRegex[["Peanut GWAS"]]
lis.datastore.qtlUrls <- readLines(paste0(lis.datastore.localDir, "datasets-qtl.txt"))

read.qtl.lis.datastore <- function(fin.qtl) {
  tmp <- tempfile()
  download.file(fin.qtl, tmp, method = "wget", quiet = TRUE)
  df.qtl <- read.csv(gzfile(tmp), header = TRUE, sep = '\t', stringsAsFactors = FALSE)
  df.qtl <- df.qtl[, 1:2]
  names(df.qtl) <- c("identifier", "trait_id")

  # QTL -> marker
  fin.qtlmrk <- gsub("\\.qtl\\.tsv\\.gz", "\\.qtlmrk\\.tsv\\.gz", fin.qtl)
  download.file(fin.qtlmrk, tmp, method = "wget", quiet = TRUE)
  df.qtlmrk <- read.csv(gzfile(tmp), header = TRUE, sep = '\t', stringsAsFactors = FALSE)
  df.qtlmrk <- df.qtlmrk[, 1:2]
  names(df.qtlmrk) <- c("identifier", "marker")
  df.qtl <- merge(df.qtlmrk, df.qtl, by = "identifier", sort = FALSE)

  # trait id -> trait ontology -> trait name
  fin.obo <- gsub("\\.qtl\\.tsv\\.gz", "\\.obo\\.tsv\\.gz", fin.qtl)
  download.file(fin.obo, tmp, method = "wget", quiet = TRUE)
  df.obo <- read.csv(gzfile(tmp), header = TRUE, sep = '\t', stringsAsFactors = FALSE)
  unlink(tmp)

  names(df.obo) <- c("trait_id", "ontology_code")
  df.to <- df.obo[startsWith(df.obo$ontology_code, "TO:"), ]
  ontologyId <- sapply(df.qtl$trait_id, function(id) {
    # prioritize TO
    oids <- df.to$ontology_code[df.to$trait_id == id]
    if (length(oids) == 0) oids <- df.obo$ontology_code[df.obo$trait_id == id]
    ifelse(length(oids) == 0, "", oids[1]) # choose the first one
  }, USE.NAMES = FALSE)
  trait <- sapply(ontologyId, function(oid) {
    traitName <- allOntologies$name[allOntologies$id == oid] # there should be at most one
    ifelse(length(traitName) == 0, "", stri_trans_totitle(traitName, type = "sentence"))
  }, USE.NAMES = FALSE)
  # for missing ones, use the trait id
  bb <- which(nchar(trait) == 0)
  trait[bb] <- df.qtl$trait_id[bb]
  df.qtl$trait <- trait

  df.qtl$publication <- read.metadata(fin.qtl)$publication
  df.qtl
}

merge.qtl <- function(df.qtl, df.gff) {
  # df.qtl is the QTL data frame merged from qtlmrk, qtl, obo files:
  #   identifier,marker,trait_id,trait,publication
  # df.gff is the processed GFF file:
  #   marker,chromosome,position

  # Merge on marker to convert to marker,identifier,trait_id,trait,publication,chromosome,position
  df.1 <- merge(df.qtl, df.gff, by = "marker", sort = FALSE)
  # Determine start and end position for each QTL identifier
  f2 <- function(r) c(min(r), max(r), (min(r) + max(r)) %/% 2)
  df.2 <- do.call(function(...) data.frame(..., stringsAsFactors = FALSE),
    aggregate(position ~ identifier, data = df.1, FUN = f2))
  names(df.2) <- c("identifier", "start_pos", "end_pos", "center_pos")
  df.2$markers <- sapply(df.2$identifier, function(id) {
    markers <- sort(df.1$marker[df.1$identifier == id])
    numMarkers <- length(markers)
    ifelse(numMarkers <= 3, paste(markers, collapse = ", "), sprintf("%d markers", numMarkers))
  }, USE.NAMES = FALSE)
  # Merge relevant columns
  df.3 <- unique(df.1[, c("identifier", "trait_id", "trait", "publication", "chromosome")])
  df.merged <- merge(df.3, df.2, by = "identifier", sort = FALSE)

  # Clean up column order
  df.merged <- df.merged[, c("identifier", "trait_id", "trait", "markers", "chromosome", "start_pos", "end_pos", "center_pos", "publication")]
  df.merged
}

build.qtl.from.lis.datastore <- function(key) {
  nid <- paste0("load.", gsub(" ", ".", key))
  showNotification(paste("Loading", key, "data. Please wait."), duration = NULL, id = nid, type = "message")

  # Read GFF and QTL files from file (instead of DSCensor)
  df.qtl <- init.qtl(key)

  # GFF (marker) files
  organism <- stri_match(key, regex = ".*(?= QTL)")[, 1]
  org.filter <- paste0("/", gsub(" ", "/", org.Genus_species[[organism]]), "/")
  markerUrls <- readLines(paste0(lis.datastore.localDir, "markers.txt"))
  markerUrls <- markerUrls[!startsWith(markerUrls, "#")]
  markerUrls <- markerUrls[grepl(org.filter, markerUrls)]
  if (length(markerUrls) > 0) {
    ll.mrk <- lapply(markerUrls, function(u) {
      df.mrku <- read.gff3.lis.datastore(u)
      df.mrku <- scrub.gff(df.mrku, lis.datastore.chrRegex[[key]])
      df.mrku
    })
    df.gff <- do.call(rbind, ll.mrk)
  }

  # QTL files
  qtlUrls <- lis.datastore.qtlUrls
  qtlUrls <- qtlUrls[!startsWith(qtlUrls, "#")]
  qtlUrls <- qtlUrls[grepl(org.filter, qtlUrls)]
  if (length(qtlUrls) > 0) {
    ll.qtl <- lapply(qtlUrls, function(u) {
      df.q <- NULL
      tryCatch({
        df.q <- read.qtl.lis.datastore(u)
      }, error = function(e) {
        print(e)
      })
      df.q
    })
    nn <- sapply(ll.qtl, is.null)
    if (any(nn)) ll.qtl <- ll.qtl[-which(nn)]
    if (length(ll.qtl) > 0) {
      df.f <- do.call(rbind, ll.qtl)
      # Merge QTL and marker data frames
      df.qtl <- merge.qtl(df.f, df.gff)
    }
  }
  # Add spurious column for QTL y axis position, as we will assign it dynamically in each chart
  df.qtl$val <- 0

  removeNotification(nid)
  # deduplicate the results (TODO: do we need to?)
  unique(df.qtl)
}

# Start with an empty data frame
init.qtl <- function(o.qtl) {
  # organism <- stri_match(o.qtl, regex = ".*(?= QTL)")[, 1]
  df.qtl <- data.frame(identifier = "QTL", chromosome = "1", start_pos = 1L, end_pos = 2L, center_pos = 1L, trait = "-", stringsAsFactors = FALSE)
  df.qtl <- df.qtl[-1, ]
  df.qtl
}

# --------------------------------------------------------------
