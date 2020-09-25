# --------------------------------------------------------------
# Build a GWAS data frame from files accessible by HTTP
# --------------------------------------------------------------

library(jsonlite)
library(stringi)
library(RCurl)

# --------------------------------------------------------------
# LIS Data Store information
lis.datastore.info <- list()
lis.datastore.info[["Cowpea GWAS"]] <- list(
  mrkFilter = "vigna:unguiculata",
  chrRegex = "vigun.IT97K-499-35.gnm1.(Vu\\d+)",
  mrkRegex = "ID=vigun.IT97K-499-35.gnm1.(\\S[^;]+);?"
)
lis.datastore.info[["Soybean GWAS"]] <- list(
  mrkFilter = "glycine:max",
  chrRegex = "glyma.Wm82.gnm2.(Gm\\d+)",
  mrkRegex = "ID=glyma.Wm82.gnm2.(\\S[^;]+);?"
)

# LIS Data Store methods
read.gff3.lis.datastore <- function(fin) {
  zz <- gzcon(url(fin, "r"))
  ll <- readLines(zz)
  close(zz)
  # remove header lines with leading "##"
  ll <- ll[!startsWith(ll, "#")]
  df.gff <- read.csv(textConnection(ll), header = FALSE, sep = '\t', as.is = TRUE)
  df.gff
}

scrub.gff <- function(df.gff.in, what) {
  # df.gff.in is the GFF data frame downloaded from the data store: V1-V9
  v1 <- df.gff.in[, 1]
  position <- df.gff.in[, 4] # start position (note v5 is end position which may be different)
  v9 <- df.gff.in[, 9]
  # TODO: keep other columns?

  # Extract chromosome and marker
  chromosome <- sapply(v1, function(s) {
    stri_match_first(s, regex = what$chrRegex)[2]
  })
  marker <- sapply(v9, function(s) {
    stri_match_first(s, regex = what$mrkRegex)[2]
  })

  # output
  df.gff <- data.frame(marker, chromosome, position, row.names = NULL)
  df.gff[!is.na(df.gff$chromosome), ]
}

read.gwas.lis.datastore <- function(fin) {
  zz <- gzcon(url(fin, "r"))
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
  df.gwas <- read.csv(textConnection(ll[i:length(ll)]), header = TRUE, sep = '\t', as.is = TRUE)
  df.gwas$publication <- paste0("<a href='", src.url, "' target=_blank>", src.name, "</a>")
  df.gwas
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
  # Ordering by column index should be correct, sometimes the names differ (like "phenotype" v. "trait")
  df.merged <- df.merged[, c(1, 6, 7, 3, 4, 5)]
  names(df.merged) <- c("marker", "chromosome", "position", "phenotype", "p_value", "publication")
  df.merged
}

build.gwas.from.lis.datastore <- function(key) {
  nid <- paste0("load.", gsub(" ", ".", key))
  showNotification(paste("Loading", key, "data. Please wait."), duration = NULL, id = nid, type = "message")

  # Discover GFF and GWAS files from DSCensor
  df.gwas <- init.gwas(key)
  gwasBaseUrl <- "http://dev.lis.ncgr.org:50020/api/v1/nodes/labels/"

  # GFF (marker) files
  query.mrk <- fromJSON(paste0(gwasBaseUrl, "mrk:", lis.datastore.info[[key]]$mrkFilter))
  if (length(query.mrk$data) > 0) {
    df.mrk <- query.mrk$data[endsWith(query.mrk$data$url, ".gff3.gz"), ]
    if (nrow(df.mrk) > 0) {
      for (i in 1:nrow(df.mrk)) {
        df.gff <- read.gff3.lis.datastore(df.mrk[i, "url"])
        df.gff <- scrub.gff(df.gff, lis.datastore.info[[key]])

        # Associated GWAS files
        query.gwas <- fromJSON(paste0(gwasBaseUrl, "gwas:", lis.datastore.info[[key]]$mrkFilter))
        if (length(query.gwas$data) == 0) next
        gwasFiles <- query.gwas$data$url[endsWith(query.gwas$data$url, ".tsv.gz")]

        # Read the GWAS files, merge with the GFF data, and append the results
        for (f in gwasFiles) {
          df.f <- read.gwas.lis.datastore(f)
          df.f2 <- merge.gwas(df.f, df.gff)
          if (nrow(df.gwas) == 0) {
            df.gwas <- df.f2
          } else {
            df.gwas <- rbind(df.gwas, df.f2)
          }
        }
      }
    }
  }
  removeNotification(nid)
  # deduplicate the results
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
  } else {
    df.gwas <- data.frame(chromosome = "1", position = 1L, phenotype = "-", p_value = 0.1, stringsAsFactors = FALSE)
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

# --------------------------------------------------------------
