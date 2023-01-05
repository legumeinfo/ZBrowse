readOntologyTerms <- function() {
  ontologyDir <- "www/config/lis-datastore/"
  traitOntology <- read.csv(paste0(ontologyDir, "traitOntology.tsv"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  chebiOntology <- read.csv(paste0(ontologyDir, "chebiOntology.tsv"), header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE)
  plantOntology <- read.csv(paste0(ontologyDir, "plantOntology.tsv"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  soyOntology <- read.csv(paste0(ontologyDir, "soyOntology.tsv"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  cropOntologyPhavu <- read.csv(paste0(ontologyDir, "cropOntology335_CommonBean.tsv"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  cropOntologyGlyma <- read.csv(paste0(ontologyDir, "cropOntology336_Soybean.tsv"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  cropOntologyVigun <- read.csv(paste0(ontologyDir, "cropOntology340_Cowpea.tsv"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  cropOntologyCajca <- read.csv(paste0(ontologyDir, "cropOntology341_Pigeonpea.tsv"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  cropOntologyVigra <- read.csv(paste0(ontologyDir, "cropOntology346_Mungbean.tsv"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  rbind(traitOntology, chebiOntology, plantOntology, soyOntology,
    cropOntologyPhavu, cropOntologyGlyma, cropOntologyVigun, cropOntologyCajca, cropOntologyVigra)
}

allOntologies <- readOntologyTerms()
