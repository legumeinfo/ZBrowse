# --------------------------------------------------------------
# Build a GWAS data frame from files accessible by HTTP
# --------------------------------------------------------------
library(stringi)
library(yaml)
# --------------------------------------------------------------
# LIS Data Store information
lis.datastore.chrRegex <- list()
lis.datastore.chrRegex[["Common Bean GWAS"]] <- "phavu.G19833.gnm\\d{1}.(Chr\\d+)" # was gnm1
lis.datastore.chrRegex[["Cowpea GWAS"]] <- "vigun.IT97K-499-35.gnm1.(Vu\\d+)"
lis.datastore.chrRegex[["Mung bean GWAS"]] <- "vigra.VC1973A.gnm6.(Vr\\d+)"
lis.datastore.chrRegex[["Peanut GWAS"]] <- "arahy.Tifrunner.gnm1.(Arahy.\\d+)"
lis.datastore.chrRegex[["Soybean GWAS"]] <- "glyma.Wm82.gnm\\d{1}.(Gm\\d+)"
lis.datastore.mrkRegex <- "Name=(\\w[^;]+);?"
lis.datastore.localDir <- "www/config/lis-datastore/"
lis.datastore.gwasUrls <- readLines(paste0(lis.datastore.localDir, "datasets-gwas.txt"))

# Ontology terms
traitOntology <- read.csv(paste0(lis.datastore.localDir, "traitOntology.tsv"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
chebiOntology <- read.csv(paste0(lis.datastore.localDir, "chebiOntology.tsv"), header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE)
plantOntology <- read.csv(paste0(lis.datastore.localDir, "plantOntology.tsv"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
soyOntology <- read.csv(paste0(lis.datastore.localDir, "soyOntology.tsv"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
cropOntologyPhavu <- read.csv(paste0(lis.datastore.localDir, "cropOntology335_CommonBean.tsv"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
cropOntologyGlyma <- read.csv(paste0(lis.datastore.localDir, "cropOntology336_Soybean.tsv"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
cropOntologyVigun <- read.csv(paste0(lis.datastore.localDir, "cropOntology340_Cowpea.tsv"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
cropOntologyCajca <- read.csv(paste0(lis.datastore.localDir, "cropOntology341_Pigeonpea.tsv"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
cropOntologyVigra <- read.csv(paste0(lis.datastore.localDir, "cropOntology346_Mungbean.tsv"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
allOntologies <- rbind(traitOntology, chebiOntology, plantOntology, soyOntology,
  cropOntologyPhavu, cropOntologyGlyma, cropOntologyVigun, cropOntologyCajca, cropOntologyVigra)

# LIS Data Store methods
read.gff3.lis.datastore <- function(fin) {
  # download a temporary file as R cannot handle https
  tmp <- tempfile()
  download.file(fin, tmp, method = "wget", quiet = TRUE)
  df.gff <- readGFF(tmp, columns = c("seqid", "start"), tags = "Name", filter = list(type = "genetic_marker"))
  unlink(tmp)
  df.gff
}

scrub.gff <- function(df.gff.in, chrRegex) {
  # Extract chromosome
  chromosome <- sapply(df.gff.in$seqid, function(s) stri_match_first(s, regex = chrRegex)[, 2])
  # Restore internal commas (if any) to marker name
  marker <- sapply(df.gff.in$Name, function(name) paste(name, collapse = ","))
  # output
  df.gff <- data.frame(marker, chromosome, position = df.gff.in$start, row.names = NULL, stringsAsFactors = FALSE)
  df.gff[!is.na(df.gff$chromosome), ]
}

# Read metadata from README file in same directory as input file
read.metadata <- function(fin) {
  directory <- dirname(fin)
  collection <- basename(directory)
  fReadme <- sprintf("%s/README.%s.yml", directory, collection)
  tmp <- tempfile()
  download.file(fReadme, tmp, method = "wget", quiet = TRUE)
  yy <- read_yaml(tmp)
  unlink(tmp)

  metadata <- list()
  metadata$publication <- sprintf("<a href='https://doi.org/%s' title='%s' target=_blank>%s</a>",
    yy$publication_doi, yy$synopsis, yy$identifier)
  metadata
}

read.gwas.lis.datastore <- function(fin) {
  tmp <- tempfile()
  download.file(fin, tmp, method = "wget", quiet = TRUE)
  df.gwas <- read.csv(gzfile(tmp), header = FALSE, skip = 1, sep = '\t', stringsAsFactors = FALSE)

  # trait id -> trait ontology -> trait name
  fin.obo <- gsub("\\.([^\\.]+)\\.tsv\\.gz", "\\.obo\\.tsv\\.gz", fin)
  download.file(fin.obo, tmp, method = "wget", quiet = TRUE)
  df.obo <- read.csv(gzfile(tmp), header = TRUE, sep = '\t', stringsAsFactors = FALSE)
  unlink(tmp)

  names(df.obo) <- c("trait_id", "ontology_code")
  df.to <- df.obo[startsWith(df.obo$ontology_code, "TO:"), ]
  ontologyId <- sapply(df.gwas[, 1], function(id) {
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
  trait[bb] <- df.gwas[bb, 1]
  df.gwas <- data.frame(df.gwas[, 1], trait, df.gwas[, 2], df.gwas[, 3], stringsAsFactors = FALSE)
  names(df.gwas) <- c("trait_id", "trait", "marker", "p_value") # clean up column names

  df.gwas$publication <- read.metadata(fin)$publication
  df.gwas
}

merge.gwas <- function(df.gwas, df.gff) {
  # df.gwas is the GWAS data frame downloaded from the data store:
  #   trait_id,trait,marker,p_value,publication
  # df.gff is the processed GFF file:
  #   marker,chromosome,position

  # Merge on marker to convert to marker,trait_id,trait,p_value,publication,chromosome,position
  df.merged <- merge(df.gwas, df.gff, by = "marker", sort = FALSE)

  # Clean up column order
  df.merged <- df.merged[, c("marker", "chromosome", "position", "trait_id", "trait", "p_value", "publication")]
  df.merged
}

build.gwas.from.lis.datastore <- function(key) {
  nid <- paste0("load.", gsub(" ", ".", key))
  showNotification(paste("Loading", key, "data. Please wait."), duration = NULL, id = nid, type = "message")

  # Read GFF and GWAS files from file (instead of DSCensor)
  df.gwas <- init.gwas(key)

  # GFF (marker) files
  organism <- stri_match(key, regex = ".*(?= GWAS)")[, 1]
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

  # GWAS files
  gwasUrls <- lis.datastore.gwasUrls
  gwasUrls <- gwasUrls[!startsWith(gwasUrls, "#")]
  gwasUrls <- gwasUrls[grepl(org.filter, gwasUrls)]
  if (length(gwasUrls) > 0) {
    ll.gwas <- lapply(gwasUrls, function(u) {
      df.g <- NULL
      tryCatch({
        df.g <- read.gwas.lis.datastore(u)
      }, error = function(e) {
        print(e)
      })
      df.g
    })
    nn <- sapply(ll.gwas, is.null)
    if (any(nn)) ll.gwas <- ll.gwas[-which(nn)]
    if (length(ll.gwas) > 0) {
      df.f <- do.call(rbind, ll.gwas)
      # Merge GWAS and marker data frames
      df.gwas <- merge.gwas(df.f, df.gff)
    }
  }

  removeNotification(nid)
  # deduplicate the results
  unique(df.gwas)
}

# --------------------------------------------------------------

# Remote GWAS datasets from CyVerse (for A. thaliana and M. truncatula)
gwas.filenames <- list()
gwas.filenames[["Arabidopsis thaliana"]] <- paste(url_cyverse,
  c("dl/d/F61A306C-92D2-4595-8226-A195D46EBB50/FT10.gwas",
    "dl/d/64C5BD48-CF10-4833-96B4-3C74CEC47257/FT16.gwas",
    "dl/d/BBB7DCAB-87FE-4C9E-8F81-DDF8FFEFF806/FT22.gwas",
    "dl/d/57E47B19-FF8C-47F6-AB32-9A111AB1F2A7/Trichome avg C.gwas",
    "dl/d/9B68EAA3-D105-49B7-B2C1-E1690E8BAF23/Trichome avg JA.gwas"
    # ...
  ), sep = "/")
gwas.filenames[["Medicago truncatula"]] <- paste(url_cyverse,
  c("dl/d/D7082A0C-5945-47DB-B5B0-3301C78BB12D/nodb_results.gwas",
    "dl/d/DAAF68FD-3A80-4E5D-83D5-2245C7E0D747/height_results.gwas",
    "dl/d/8F20C8BF-BEEC-4801-BCFF-41534832958B/floweringdate_results.gwas",
    "dl/d/2C0AC931-4C07-47F9-B236-DDB581497979/occupancyA_results.gwas",
    "dl/d/EE69C622-4432-479C-82CF-66BF33C1C38C/occupancyB_results.gwas",
    "dl/d/03FCF51A-8501-487F-A12E-DF4F8CDB0446/noda_results.gwas",
    "dl/d/6030E8B2-4BD9-4393-87D6-6E67D60DA7E8/trichomes_results.gwas",
    "dl/d/F2485FCD-D5D3-4F55-B693-F00551B4FF47/totalnod_results.gwas"
    # ...
  ), sep = "/")

gwas.traits <- list()
gwas.traits[["Arabidopsis thaliana"]] <- stri_match(basename(gwas.filenames[["Arabidopsis thaliana"]]), regex = ".*(?=.gwas)")[, 1]
gwas.traits[["Medicago truncatula"]] <- stri_match(basename(gwas.filenames[["Medicago truncatula"]]), regex = ".*(?=_results.gwas)")[, 1]

# TODO: standardize column names
gwas.cols <- list()
gwas.cols[["Arabidopsis thaliana"]] <- c("Chromosome", "Position", "Trait", "P.Value", "negLogP", "MAF")
gwas.cols[["Medicago truncatula"]] <- c("Chromosome", "pos", "P.value")

gwas.sources <- data.frame(
  name = c("CyVerse", "LIS Data Store"),
  url = c(gwas.filenames[[1]][1], lis.datastore.gwasUrls[1]),
  status = FALSE,
  stringsAsFactors = FALSE
)
row.names(gwas.sources) <- gwas.sources$name

# Start with an empty data frame
init.gwas <- function(o.gwas) {
  organism <- stri_match(o.gwas, regex = ".*(?= GWAS)")[, 1]
  if (organism == "Medicago truncatula") {
    df.gwas <- data.frame(Chromosome = "1", pos = 1L, Trait = "-", P.value = 0.1, stringsAsFactors = FALSE)
  } else if (organism == "Arabidopsis thaliana") {
    df.gwas <- data.frame(Chromosome = "1", Position = 1L, Trait = "-", P.Value = 0.1, negLogP = 1.0, MAF = 0.01, stringsAsFactors = FALSE)
  } else {
    df.gwas <- data.frame(chromosome = "1", position = 1L, trait = "-", p_value = 0.1, stringsAsFactors = FALSE)
  }
  df.gwas <- df.gwas[-1, ]
  df.gwas
}

load.gwas.remote <- function(organism, filename, trait) {
  filename <- URLencode(filename)
  if (!gwas.sources["CyVerse", "status"]) {
    print(paste("No connection to", filename))
    return()
  }

  t0 <- proc.time()[3]
  cat(paste("Loading", trait, "data ... "))
  nid <- paste0("load.remote.", stri_replace_all(trait, "_", regex = " "))
  showNotification(paste0("Loading remote trait ", trait, ". Please wait."),
    duration = NULL, id = nid, type = "message")

  filenames <- gwas.filenames[[organism]]
  traits <- gwas.traits[[organism]]
  cols <- gwas.cols[[organism]]

  if (organism == "Medicago truncatula") {
    df.gwas <- read.table(file = url(filename, method = "libcurl"), header = TRUE, sep = "\t", quote = "\"", stringsAsFactors = FALSE)[, cols]
    df.gwas$Trait <- trait
    # Clean up: move Trait to the third column
    df.gwas <- df.gwas[, c("Chromosome", "pos", "Trait", "P.value")]

  } else if (organism == "Arabidopsis thaliana") {
    df.gwas <- read.table(file = url(filename, method = "libcurl"), header = TRUE, sep = ",", quote = "\"", stringsAsFactors = FALSE)[, cols]
    df.gwas$Chromosome <- as.character(df.gwas$Chromosome)
    df.gwas$Trait <- trait
  }

  cat(sprintf("Done. (%2.1f seconds)\n", proc.time()[3] - t0))
  removeNotification(nid)
  df.gwas
}

# --------------------------------------------------------------
