# --------------------------------------------------------------
# Scan for GWAS and QTL files in the LIS data store, as in the specifications:
# https://github.com/legumeinfo/datastore-specifications/tree/main/Genus/species/[gwas|qtl|markers]

# TODO: investigate why this works with curl >= 7.59.0 but not <= 7.54.0
# (https://curl.se/changes.html - not sure at which version it starts working)
# --------------------------------------------------------------
library(RCurl)
library(readr)
library(stringi)

# set path to the data directory
# setwd(".../ZZBrowse/www/config/lis-datastore")
# --------------------------------------------------------------

url.data.store <- "https://data.legumeinfo.org/"
genus.species <- c("Arachis/hypogaea/", "Cajanus/cajan/",
  "Glycine/max/", "Medicago/truncatula/", "Phaseolus/vulgaris/",
  "Vigna/radiata/", "Vigna/unguiculata/")

buildFileList <- function(gsDirs, folders, types, mtypes, fout) {
  rgx.xx <- paste0(gsub("\\/$", "", url.data.store), "(.+)")
  if (file.exists(fout)) file.remove(fout)
  ftmp <- tempfile()
  for (gs in gsDirs) {
    for (ff in folders) {
      url.i <- paste0(url.data.store, gs, ff, "/")
      if (!url.exists(url.i)) next
      system(sprintf('curl -s %s -o %s', url.i, ftmp))
      txt <- read_file(ftmp)
      xx <- stri_match_first(url.i, regex = rgx.xx)[, 2]
      rgx <- paste0('<a href="', xx, '([^<>]+)">([^<>]+)</a>')
      dd <- URLdecode(stri_match_all(txt, regex = rgx)[[1]][, 2])
      urls.j <- paste0(url.i, dd)
      for (url.j in urls.j) {
        if (!url.exists(url.j)) next
        system(sprintf('curl -s %s -o %s', url.j, ftmp))
        txt <- read_file(ftmp)
        xx <- stri_match_first(url.i, regex = rgx.xx)[, 2]
        rgx <- paste0('<a href="', xx, '([^<>]+)">([^<>]+)</a>')
        dd2 <- URLdecode(stri_match_all(txt, regex = rgx)[[1]][, 2])
        bb <- rep(FALSE, length(dd2))
        for (tt in types) {
          bb <- bb | endsWith(dd2, tt)
        }
        dd2 <- dd2[bb]
        # check that each data file in dd2 has a counterpart of each type in mtypes
        bb <- rep(TRUE, length(dd2))
        for (k in 1:length(dd2)) {
          for (mtype in mtypes) {
            dd2.m <- paste0(url.i, gsub(types[1], mtype, dd2[k]))
            bb[k] <- bb[k] && url.exists(dd2.m)
          }
        }
        dd2 <- dd2[bb]
        if (length(dd2) > 0) write(paste0(url.i, dd2), fout, append = TRUE)
      }
    }
  }
  unlink(ftmp)
}

# GWAS datasets
buildFileList(gsDirs = genus.species,
  folders = c("gwas"),
  types = c(".result.tsv.gz"),
  mtypes = c(".obo.tsv.gz"),
  fout = "datasets-gwas.txt"
)

# QTL datasets
buildFileList(gsDirs = genus.species,
  folders = c("qtl"),
  types = c(".qtlmrk.tsv.gz"),
  mtypes = c(".obo.tsv.gz"),
  fout = "datasets-qtl.txt"
)

# Marker files
buildFileList(gsDirs = genus.species,
  folders = c("markers"),
  types = c(".gff3.gz"),
  mtypes = character(0),
  fout = "markers.txt"
)

# --------------------------------------------------------------

# Download data files from a list of datasets created by buildFileList()
downloadFileList <- function(fin, dirOut) {
  ff <- readLines(fin)
  ff <- ff[!startsWith(ff, "#")]
  if (!dir.exists(dirOut)) dir.create(dirOut)
  setwd(dirOut)
  for (fi in ff) {
    system(sprintf("curl -s %s -o %s", URLencode(fi), basename(fi)))
  }
}

# Download GWAS files
downloadFileList(fin = "datasets-gwas.txt", dirOut = "gwas")

# Download QTL files
downloadFileList(fin = "datasets-qtl.txt", dirOut = "qtl")

# --------------------------------------------------------------

# Trait ontology
buildTraitOntology <- function() {
  ll <- readLines("https://raw.githubusercontent.com/Planteome/plant-trait-ontology/master/to.obo")
  ii <- which(grepl("^id: TO:", ll))
  to.id <- stri_match_first(ll[ii], regex = "^id: (.+)$")[, 2]
  to.name <- stri_match_first(ll[ii + 1], regex = "^name: (.+)$")[, 2]
  df.to <- data.frame(id = to.id, name = to.name, stringsAsFactors = FALSE)
  write.table(df.to, "traitOntology.tsv", sep = "\t", row.names = FALSE)
}
buildTraitOntology()

# CHEBI ontology
buildChebiOntology <- function() {
  tmp <- tempfile()
  download.file("https://ftp.ebi.ac.uk/pub/databases/chebi/ontology/chebi.obo.gz", tmp, method = "wget", quiet = TRUE)
  ll <- readLines(gzfile(tmp))
  unlink(tmp)
  ii <- which(grepl("^id: CHEBI:", ll))
  chebi.id <- stri_match_first(ll[ii], regex = "^id: (.+)$")[, 2]
  chebi.name <- stri_match_first(ll[ii + 1], regex = "^name: (.+)$")[, 2]
  df.chebi <- data.frame(id = chebi.id, name = chebi.name, stringsAsFactors = FALSE)
  # no quotes, as it contains double quotes
  write.table(df.chebi, "chebiOntology.tsv", quote = FALSE, sep = "\t", row.names = FALSE)
}
buildChebiOntology()

# Plant ontology
buildPlantOntology <- function() {
  ll <- readLines("https://data.bioontology.org/ontologies/PO/submissions/25/download?apikey=8b5b7825-538d-40e0-9e9e-5ab9274a9aeb")
  ii <- which(grepl("^id: PO:", ll))
  po.id <- stri_match_first(ll[ii], regex = "^id: (.+)$")[, 2]
  po.name <- stri_match_first(ll[ii + 1], regex = "^name: (.+)$")[, 2]
  df.po <- data.frame(id = po.id, name = po.name, stringsAsFactors = FALSE)
  write.table(df.po, "plantOntology.tsv", sep = "\t", row.names = FALSE)
}
buildPlantOntology()

# Soy ontology
buildSoyOntology <- function() {
  ll <- readLines("http://data.agroportal.lirmm.fr/ontologies/SOY/submissions/7/download?apikey=1de0a270-29c5-4dda-b043-7c3580628cd5")
  ii <- which(grepl("^id: SOY:", ll))
  soy.id <- stri_match_first(ll[ii], regex = "^id: (.+)$")[, 2]
  soy.name <- stri_match_first(ll[ii + 1], regex = "^name: (.+)$")[, 2]
  df.soy <- data.frame(id = soy.id, name = soy.name, stringsAsFactors = FALSE)
  write.table(df.soy, "soyOntology.tsv", sep = "\t", row.names = FALSE)
}
buildSoyOntology()

# Crop ontologies (one per species, from Excel spreadsheets at https://cropontology.org)
scrubCropOntology <- function(fin) {
  df.co <- read.csv(fin, header = TRUE, sep = "\t")
  df.co <- unique(df.co)
  names(df.co) <- c("id", "name")
  write.table(df.co, fin, sep = "\t", row.names = FALSE)
}
scrubCropOntology("cropOntology335_CommonBean.tsv")
scrubCropOntology("cropOntology336_Soybean.tsv")
scrubCropOntology("cropOntology340_Cowpea.tsv")
scrubCropOntology("cropOntology341_Pigeonpea.tsv")
scrubCropOntology("cropOntology346_Mungbean.tsv")

# --------------------------------------------------------------
