# --------------------------------------------------------------
# Build a GWAS data frame from files accessible by HTTP
# --------------------------------------------------------------

library(jsonlite)
library(stringi)
library(RCurl)

# --------------------------------------------------------------
# LIS Data Store methods

# TODO: add column names
read.gz.url <- function(fin, skip = 0, header = FALSE) {
  zz <- gzcon(url(fin, "r"))
  xx <- read.csv(textConnection(readLines(zz)), header = header, sep = '\t', as.is = TRUE, skip = skip)
  close(zz)
  xx
}

# returns the number of lines before the first line starting with ch
determine.skip <- function(fin, ch = "#", nmax = 16) {
  zz <- gzcon(url(fin, "r"))
  xx <- readLines(zz, n = nmax)
  close(zz)
  bb <- startsWith(xx, ch)
  if (!any(bb)) return(NA)
  which(bb)[1] - 1
}

scrub.gff <- function(gff) {
  # gff is the GFF data frame downloaded from the data store: V1-V9
  v1 <- gff[, 1]
  position <- gff[, 4] # start position (note v5 is end position which may be different)
  v9 <- gff[, 9]
  # TODO: keep other columns?

  # Extract chromosome and marker
  # Note: for other species the format (and thus the parsing) may be different.
  chromosome <- sapply(v1, function(s) {
    # stri_match_first(s, regex = "\\d+$")[1]
    # as.integer(stri_sub(s, nchar(s) - 1))
    stri_match_first(s, regex = "glyma.Wm82.gnm2.(Gm\\d+)")[2]
  })
  marker <- sapply(v9, function(s) {
    stri_match_first(s, regex = "ID=glyma.Wm82.gnm2.(\\S[^;]+);?")[2]
  })

  # output
  df.out <- data.frame(marker, chromosome, position, row.names = NULL)
  df.out[!is.na(df.out$chromosome), ]
}

merge.gwas <- function(df.gwas, df.gff) {
  # df.gwas is the GWAS data frame downloaded from the data store:
  #   X.identifier,phenotype,marker,pvalue
  # TODO: keep other columns?
  # df.gff is the processed GFF file:
  #   marker,chromosome,position

  # Merge on marker to convert to chromosome,phenotype(=trait),position,p_value
  df.merged <- merge(df.gwas, df.gff, by = "marker", sort = FALSE)

  # Clean up column order and names
  df.merged <- df.merged[, c(1, 5, 6, 3, 4)]
  names(df.merged) <- c("marker", "chromosome", "position", "phenotype", "p_value")
  df.merged
}

build.gwas.from.lis.datastore <- function(gwasDir, gffFile) {
  nid <- "load.gwas"
  showNotification("Loading Soybean GWAS data. Please wait.", duration = NULL, id = nid, type = "message")

  df.gff <- read.gz.url(gffFile, skip = 2)
  df.gff <- scrub.gff(df.gff)

  # TODO: equivalent of list.files, or discover from DSCensor
  #gwasFiles <- list.files(gwasDir, pattern = "gwas.tsv.gz$")
  gwasFiles <- paste("https://legumeinfo.org/data/public/Glycine_max/mixed.gwas.1W14", c(
    "glyma.mixed.gwas.1W14.KGK20170707-1.gwas.tsv.gz",
    "glyma.mixed.gwas.1W14.KGK20170711-1.gwas.tsv.gz",
    "glyma.mixed.gwas.1W14.KGK20170714-1.gwas.tsv.gz",
    "glyma.mixed.gwas.1W14.KGK20170803-1.gwas.tsv.gz",
    "glyma.mixed.gwas.1W14.KGK20170808-1.gwas.tsv.gz",
    "glyma.mixed.gwas.1W14.KGK20170814-1.gwas.tsv.gz",
    "glyma.mixed.gwas.1W14.KGK20170908-1.gwas.tsv.gz",
    "glyma.mixed.gwas.1W14.KGK20170915-1.gwas.tsv.gz",
    "glyma.mixed.gwas.1W14.KGK20171006-1.gwas.tsv.gz",
    "glyma.mixed.gwas.1W14.KGK20171018-1.gwas.tsv.gz",
    "glyma.mixed.gwas.1W14.KGK20171027-1.gwas.tsv.gz",
    "glyma.mixed.gwas.1W14.LBC20180516-2.gwas.tsv.gz",
    "glyma.mixed.gwas.1W14.LBC20180516-3.gwas.tsv.gz",
    "glyma.mixed.gwas.1W14.LBC20180521-1.gwas.tsv.gz",
    "glyma.mixed.gwas.1W14.LBC20180521-5.gwas.tsv.gz",
    "glyma.mixed.gwas.1W14.LBC20180601-3.gwas.tsv.gz",
    "glyma.mixed.gwas.1W14.LBC20180625-2.gwas.tsv.gz",
    "glyma.mixed.gwas.1W14.LBC20180625-3.gwas.tsv.gz",
    "glyma.mixed.gwas.1W14.SST20180209-1.gwas.tsv.gz"
  ), sep = "/")

  # Read the GWAS files, merge with the GFF data, and append the results
  df.gwas <- NULL
  for (f in gwasFiles) {
    nskip <- determine.skip(f)
    df.f <- read.gz.url(f, skip = nskip, header = TRUE)
    df.f2 <- merge.gwas(df.f, df.gff)
    if (is.null(df.gwas)) {
      df.gwas <- df.f2
    } else {
      df.gwas <- rbind(df.gwas, df.f2)
    }
  }
  removeNotification(nid)

  # sort and deduplicate the output
  df.gwas <- df.gwas[order(df.gwas$chromosome, df.gwas$position), ]
  unique(df.gwas)
}

# --------------------------------------------------------------

gwas.sources <- data.frame(
  name = c("CyVerse", "DSCensor", "LIS Data Store"),
  # use more basic URLs
  url = c("http://de.cyverse.org/dl/d/F61A306C-92D2-4595-8226-A195D46EBB50/FT10.gwas",
          "http://dev.lis.ncgr.org:50020",
          "https://legumeinfo.org/data/public/Glycine_max/mixed.gwas.1W14/glyma.mixed.gwas.1W14.KGK20170707-1.gwas.tsv.gz"),
  status = FALSE,
  stringsAsFactors = FALSE
)
row.names(gwas.sources) <- gwas.sources$name
gwas.sources$status <- sapply(gwas.sources$url, url.exists)

# The following structure makes it easy to add any organism whose
# GWAS files are in the same format as those for Medicago truncatula.
# TODO: generalize for other formats, if necessary.

gwas.filenames <- list()
gwas.filenames[["Arabidopsis thaliana"]] <- c(
  "http://de.cyverse.org/dl/d/F61A306C-92D2-4595-8226-A195D46EBB50/FT10.gwas",
  "http://de.cyverse.org/dl/d/64C5BD48-CF10-4833-96B4-3C74CEC47257/FT16.gwas",
  "http://de.cyverse.org/dl/d/BBB7DCAB-87FE-4C9E-8F81-DDF8FFEFF806/FT22.gwas",
  "http://de.cyverse.org/dl/d/57E47B19-FF8C-47F6-AB32-9A111AB1F2A7/Trichome avg C.gwas",
  "http://de.cyverse.org/dl/d/9B68EAA3-D105-49B7-B2C1-E1690E8BAF23/Trichome avg JA.gwas"
  # ...
)

# Get the GWAS filenames from DSCensor
gwasBaseUrl <- "http://dev.lis.ncgr.org:50020/api/v1/nodes/labels/"
organisms.gwas <- c("Medicago truncatula") #, "Soybean", "Pigeonpea", "Cowpea")
for (o.gwas in organisms.gwas) {
  ss <- strsplit(org.Genus_species[[o.gwas]], split = " ")[[1]]
  genus <- ss[1]
  species <- ss[2]

  # If the connection is valid, return the query as a list containing a data frame,
  # otherwise return an empty list
  query <- list()
  if (gwas.sources["DSCensor", "status"]) tryCatch(
    query <- fromJSON(paste0(gwasBaseUrl, tolower(genus), ":", species, ":gwas"))
  )
  if (length(query) == 0 || length(query$data$url) == 0) {
    gwas.filenames[[paste(ss[1], ss[2])]] <- c()
  } else {
    oo <- order(sapply(query$data$url, basename))
    gwas.filenames[[paste(ss[1], ss[2])]] <- query$data$url[oo]
  }
}

gwas.traits <- list()
gwas.traits[["Arabidopsis thaliana"]] <- stri_match(basename(gwas.filenames[["Arabidopsis thaliana"]]), regex = ".*(?=.gwas)")[, 1]
for (o.gwas in organisms.gwas) {
  if (length(gwas.filenames[[o.gwas]]) == 0) {
    gwas.traits[[o.gwas]] <- c()
  } else {
    gwas.traits[[o.gwas]] <- stri_match(basename(gwas.filenames[[o.gwas]]), regex = ".*(?=_results.gwas)")[, 1]
  }
}

# TODO: standardize column names
gwas.cols <- list()
gwas.cols[["Arabidopsis thaliana"]] <- c("Chromosome", "Position", "Trait", "P.Value", "negLogP", "MAF")
for (o.gwas in organisms.gwas) {
  # TODO: this works for Medicago truncatula;
  # either make sure it works for other species, or generalize it
  gwas.cols[[o.gwas]] <- c("Chromosome", "pos", "P.value")
}

# Start with an empty data frame
init.gwas <- function(o.gwas) {
  organism <- stri_match(o.gwas, regex = ".*(?= GWAS)")[, 1]
  if (organism == "Medicago truncatula") {
    # TODO: this works for Medicago truncatula;
    # either make sure it works for other species, or generalize it
    df.gwas <- data.frame(Chromosome = "1", pos = 1L, Trait = "-", P.value = 0.1, stringsAsFactors = FALSE)
  } else if (organism == "Arabidopsis thaliana") {
    df.gwas <- data.frame(Chromosome = "1", Position = 1L, Trait = "-", P.Value = 0.1, negLogP = 1.0, MAF = 0.01, stringsAsFactors = FALSE)
  }
  df.gwas <- df.gwas[-1, ]
  df.gwas
}

load.gwas.remote <- function(organism, filename, trait) {
  filename <- gsub(" ", "%20", filename)
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
    # TODO: this works for Medicago truncatula;
    # either make sure it works for other species, or generalize it
    df.gwas <- read.table(file = url(filename, method = "libcurl"), header = TRUE, sep = "\t", quote = "\"", stringsAsFactors = FALSE)[, cols]
    df.gwas$Trait <- trait

    # Clean up: change Chromosome format from "chr1" to "1", and move Trait to the third column
    df.gwas$Chromosome <- sapply(df.gwas$Chromosome, FUN = function(chr) stri_sub(chr, 4))
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

# load.gwas.local <- function(organism, filename, trait) {
#   # TODO ...
# }

# Tests
# df.gwas <- load.gwas.local("Medicago truncatula", ..., "floweringdate")
# write.csv(df.gwas, file = "Medicago_truncatula_GWAS.csv", quote = FALSE, row.names = FALSE)

# --------------------------------------------------------------
