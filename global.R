#Load GenomicRange packages locally, order here matters
pkgs <- c("BiocGenerics","S4Vectors","IRanges","XVector","Rsamtools","GenomeInfoDb","GenomicRanges")
if (!require("BiocManager")) install.packages("BiocManager", dependencies = TRUE)

for(p in pkgs){
  if(system.file(package=p) == ""){
    BiocManager::install(p, update = FALSE)
  }  
  suppressPackageStartupMessages(library(p, quietly=TRUE, character.only=TRUE))
}
  

library(shiny)
library(plyr)
library(rCharts)
library(xtable)
library(tools)
#install shinyIncubator like this (gives progress bar):
#install.packages("devtools")
#library(devtools)
#devtools::install_github("shiny-incubator", "rstudio")
#library(shinyIncubator)
#options(shiny.maxRequestSize=-1)
library(jsonlite)
library(shinyjs)
library(rintrojs)

addResourcePath('datatables','www/DataTables/')
addResourcePath('tabletools','www/TableTools/')
addResourcePath('highcharts','www/highcharts/')

#find the non-numeric values in a vector
which.nonnum <- function(x) {
  badNum <- is.na(suppressWarnings(as.numeric(as.character(x))))
  which(badNum & !is.na(x))
}

# For constructing the annotations data frame on the fly
# (must go before constructing organism-specific properties)
source("./buildAnnotations.R")

# Add each organism's properties from its configuration file
# instead of using hardcoded arrays
chrSize<-list()
chrName<-list()
# Species full taxonomic name (Genus species)
org.Genus_species <- list()
# Species abbreviation (G.species)
org.G.species <- list()
# 5-letter species abbreviation (Gensp)
org.Gensp <- list()
# Annotations file
org.annotGeneLoc<-list()
# Chromosome prefix for annotations
org.annotChrPrefix <- list()
# Base URL for Services API genomic linkage queries
org.gcvUrlBase <- list()
# Chromosome name format, for sending chromosome-based queries to the Genome Context Viewer
org.gcvChrFormat <- list()

# For constructing gene mouseover text (annotTable, annotTableReverse) in zChart
org.tag_strand <- list()
org.strand_fwd <- list()
org.strand_rev <- list()
org.tag_start <- list()
org.tag_end <- list()
org.urlFormat <- list()
org.tag_url <- list()
org.tag_name <- list()
org.tag_chr <- list()
org.tag_desc <- list()

files<-list.files(path="./organisms/")
for(i in 1:length(files)){
  if(tools::file_ext(files[i]) == "txt"){
    filename=""
    filename=paste("./organisms/",files[i],sep="")
    conn=file(filename,open="r")
    data<-readLines(conn)
    
    key<-data[1]
    #added ability to specify chrom names in organisms file by separating with a :, otherwise it just assumes they are alphanumerically sorted    
    if(length(grep(":",data[2]))){
      features <- strsplit(data[2], ",")[[1]]
      names <- read.table(text=features,sep=":",stringsAsFactors = FALSE) #this will die if any names are missing
      #order so that numeric come first in order, then named chrs
      if(length(which.nonnum(names$V1))>0){
        numNames <- names[-which.nonnum(names$V1),]
      }else{
        numNames <- names
      }
      numNames <- numNames[order(as.numeric(numNames$V1)),]
      nonNumNames <- names[which.nonnum(names$V1),]
      nonNumNames <- nonNumNames[order(nonNumNames$V1),]
      namesOrdered <- rbind(numNames,nonNumNames)
      value <- namesOrdered$V2
      name <- namesOrdered$V1
    }else{
      value<-unlist(c(lapply(strsplit(data[2], ","), as.numeric)))
      name <- as.character(1:length(value))
    }    
    chrSize[key]<-list(value)
    chrName[key]<-list(name)

    ss.org.names <- strsplit(data[3], split = ",")[[1]]
    org.Genus_species[key] <- ss.org.names[1]
    org.G.species[key] <- ss.org.names[2]
    org.Gensp[key] <- ss.org.names[3]
    annotFilename <- data[4]
    org.annotChrPrefix[key] <- data[5]
    if (stri_endswith_fixed(annotFilename, "gff3.gz")) {
      chromosome.lengths <- chrSize[key][[1]]
      locValue <- build.annotations(key, annotFilename, chromosome.lengths, org.annotChrPrefix[key])
    } else {
      locValue<-read.table(annotFilename,sep=",",head=TRUE,stringsAsFactors = FALSE,quote = c("\""))
    }
    chr2i <- suppressWarnings(as.integer(locValue$chromosome))
    if (!any(is.na(chr2i))) locValue$chromosome <- as.character(chr2i)
    org.annotGeneLoc[key]<-list(locValue)
    org.gcvUrlBase[key] <- data[6]
    org.gcvChrFormat[key] <- data[7]

    ss.org.annot <- strsplit(data[8], split = ",")[[1]]
    org.tag_strand[key] <- ss.org.annot[1]
    org.strand_fwd[key] <- ss.org.annot[2]
    org.strand_rev[key] <- ss.org.annot[3]
    org.tag_start[key] <- ss.org.annot[4]
    org.tag_end[key] <- ss.org.annot[5]
    org.urlFormat[key] <- ss.org.annot[6]
    org.tag_url[key] <- ss.org.annot[7]
    org.tag_name[key] <- ss.org.annot[8]
    org.tag_chr[key] <- ss.org.annot[9]
    org.tag_desc[key] <- ss.org.annot[10]

    close(conn)
  }
}

# For constructing the GWAS data frame on the fly
# (must go after creating organism-specific properties, as it uses org.Genus_species)
source("./buildGWAS.R")

helpPopup <- function(title, content, placement=c('right', 'top', 'left', 'bottom'),
                      trigger=c('click', 'hover', 'focus', 'manual')) {  
  tagList(
    singleton(tags$head(tags$script("$(function() { $(\"[data-toggle='popover']\").popover(); })"))),
    tags$a(href = "#", `data-toggle` = "popover", title = title, `data-content` = content,
           `data-placement` = match.arg(placement, several.ok=TRUE)[1],
           `data-trigger` = match.arg(trigger, several.ok=TRUE)[1], tags$i(class="icon-question-sign"))
  )
}

helpModal <- function(title, link, content) {
  html <- sprintf("<div id='%s' class='modal hide fade in' style='display: none; '>
<div class='modal-header'><a class='close' data-dismiss='modal' href='#'>&times;</a>
<h3>%s</h3>
</div>
<div class='modal-body'>%s</div>
</div>
<a data-toggle='modal' href='#%s' class='icon-question-sign'></a>", link, title, content, link)
  Encoding(html) <- 'UTF-8'
  HTML(html)
}

# binding for a text input that only updates when the return key is pressed
returnTextInput <- function(inputId, label, value = "") {
  tagList(
    singleton(tags$head(tags$script(src = "js/returnTextInputBinding.js"))),
    tags$label(label, `for` = inputId),
    tags$input(id = inputId, type = "text", value = value, class = "returnTextInput")
  )
}

#functions from shiny leaflet example, sets up div columns in page layout
row <- function(...) {
  tags$div(class="row", ...)
}

col <- function(width, ...) {
  tags$div(class=paste0("span", width), ...)
}

#initiate colors for plot data (these 10 are the highchart defaults)
colors <- c('#2f7ed8','#0d233a', '#8bbc21','#910000','#1aadce','#492970','#f28f43','#77a1e5','#c42525','#a6c96a')

#a list of 20 colors from http://stackoverflow.com/questions/470690/how-to-automatically-generate-n-distinct-colors
moreColors <- c(
  '#FFB300', #Vivid Yellow
  '#803E75', #Strong Purple
  '#FF6800', #Vivid Orange
  '#A6BDD7', #Very Light Blue
  '#C10020', #Vivid Red
  '#CEA262', #Grayish Yellow
  '#817066', #Medium Gray
  '#007D34', #Vivid Green
  '#F6768E', #Strong Purplish Pink
  '#00538A', #Strong Blue
  '#FF7A5C', #Strong Yellowish Pink
  '#53377A', #Strong Violet
  '#FF8E00', #Vivid Orange Yellow
  '#B32851', #Strong Purplish Red
  '#F4C800', #Vivid Greenish Yellow
  '#7F180D', #Strong Reddish Brown
  '#93AA00', #Vivid Yellowish Green
  '#593315', #Deep Yellowish Brown
  '#F13A13', #Vivid Reddish Orange
  '#232C16'  #Dark Olive Green
)

allColors <- c(moreColors,colors)

# background colors for comparing two organisms
bgColors <- c("lightblue", "lightsalmon")

# locally configurable global variables
source("private.R")
