# GCV microservices used in ZZBrowse
# https://github.com/legumeinfo/gcv-microservices

library(RCurl) # for basicTextGatherer()

genesMicroservice <- function(url, genes) {
  pf <- list(genes = genes)
  btg <- basicTextGatherer()
  curlPerform(url = sprintf("%s/genes", url), postfields = toJSON(pf), writefunction = btg$update)
  btg.text <- btg$value()
  if (!startsWith(btg.text, "{\"genes\":")) {
    return(list(error = sprintf("%s in genes microservice", btg.text)))
  }
  results <- fromJSON(btg.text, simplifyVector = FALSE)
  results
}

chromosomeMicroservice <- function(url, chromosome) {
  btg <- basicTextGatherer()
  curlPerform(url = sprintf("%s/chromosome?chromosome=%s", url, chromosome), writefunction = btg$update)
  btg.text <- btg$value()
  if (!startsWith(btg.text, "{\"chromosome\":")) {
    return(list(error = sprintf("%s in chromosome microservice (%s)", btg.text, chromosome)))
  }
  results <- fromJSON(btg.text, simplifyVector = FALSE)
  results
}

microSyntenySearchMicroservice <- function(url, families, matched, intermediate) {
  pf <- list(query = families, matched = unbox(matched), intermediate = unbox(intermediate))
  btg <- basicTextGatherer()
  curlPerform(url = sprintf("%s/micro-synteny-search", url), postfields = toJSON(pf), writefunction = btg$update)
  btg.text <- btg$value()
  if (!startsWith(btg.text, "{\"tracks\":")) {
    return(list(error = sprintf("%s in micro-synteny-search microservice", btg.text)))
  }
  results <- fromJSON(btg.text, simplifyVector = FALSE)
  results
}

macroSyntenyBlocksMicroservice <- function(url, families, matched, intermediate, mask, targets = NULL, distanceMetric = "levenshtein", chr1name = "", org1 = "", org2 = "") {
  pf <- list(chromosome = families, matched = unbox(matched), intermediate = unbox(intermediate), mask = unbox(mask), optionalMetrics = distanceMetric)
  if (!is.null(targets)) pf$targets <- targets
  btg <- basicTextGatherer()
  curlPerform(url = sprintf("%s/macro-synteny-blocks", url), postfields = toJSON(pf), writefunction = btg$update)
  btg.text <- btg$value()
  if (!startsWith(btg.text, "{\"blocks\":")) {
    return(list(error = sprintf("%s in macro-synteny-blocks microservice (%s %s - %s)", btg.text, org1, chr1name, org2)))
  }
  results <- fromJSON(btg.text, simplifyVector = FALSE)

  # Process the results
  chr1 <- trailingInteger(chr1name)
  # Construct blocks for species 2
  bb2 <- do.call(rbind, lapply(results$blocks, function(b) {
    do.call(rbind, lapply(b$blocks, function(bi) {
      # Note that bi$i, bi$j are the gene indices (0-based, so add 1 for R) of the species 1 block,
      # and bi$fmin, bi$fmax are the positions (on b$chromosome) of the species 2 block.
      data.frame(chromosome = trailingInteger(b$chromosome), i = bi$i + 1, j = bi$j + 1,
        fmin = bi$fmin, fmax = bi$fmax, orientation = bi$orientation, distance = bi$optionalMetrics[[1]],
        stringsAsFactors = FALSE)
    }))
  }))
  if (is.null(bb2)) return(list(NULL, NULL))
  cumbp2 <- c(0, cumsum(as.numeric(chrSize[[org2]])))
  bb2$cumfmin <- bb2$fmin + cumbp2[bb2$chromosome]
  bb2$cumfmax <- bb2$fmax + cumbp2[bb2$chromosome]
  bb2$chr1 <- chr1
  bb2 <- bb2[order(bb2$cumfmin, bb2$cumfmax), ]
  # Determine number of genes in species 2 blocks, for normalization
  df.annot2 <- org.annotGeneLoc[[org2]]
  gg.chr <- trailingInteger(df.annot2$chromosome)
  gg.fmin <- df.annot2$transcript_start
  gg.fmax <- df.annot2$transcript_end
  bb2$n2 <- apply(bb2[, c("chromosome", "fmin", "fmax")], 1, function(b) {
    sum(gg.chr == b["chromosome"] & gg.fmin >= b["fmin"] & gg.fmax <= b["fmax"])
  })

  # Construct blocks for species 1
  cumbp1 <- c(0, cumsum(as.numeric(chrSize[[org1]])))
  bb1 <- bb2[, c("i", "j")]
  bb1$chromosome <- chr1
  df.annot1 <- subset(org.annotGeneLoc[[org1]], trailingInteger(chromosome) == chr1)
  bb1$fmin <- df.annot1[bb1$i, "transcript_start"]
  bb1$fmax <- df.annot1[bb1$j, "transcript_end"]
  bb1$cumfmin <- bb1$fmin + cumbp1[chr1]
  bb1$cumfmax <- bb1$fmax + cumbp1[chr1]
  # Determine number of genes in species 1 blocks, for normalization
  bb2$n1 <- apply(bb1[, c("fmin", "fmax")], 1, function(b) {
    sum(df.annot1$transcript_start >= b["fmin"] & df.annot1$transcript_end <= b["fmax"])
  })
  # sort
  bb1 <- bb1[order(bb1$cumfmin, bb1$cumfmax), ]
  
  # Match and assign colors
  bb1$color <- bb2$color <- macrosyntenyColors[chr1]
  list(bb1, bb2)
}
