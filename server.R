shinyServer(function(input, output, session) {
  # Construct the full name of a named list item for species j (j = 1 or 2).
  # for example, input$datasets == input[[jth_ref("datasets", 1)]]
  # and input$datasets2 == input[[jth_ref("datasets", 2)]]
  jth_ref <- function(name, j) {
    ifelse(j == 1, name, paste0(name, j))
  }

  #Load any saved datasets
  values <- reactiveValues()
  dataPath <- "./www/config/data/"
  dataFiles <- list.files(dataPath,recursive=T)
  # Append names of any data we will create on the fly
  legumeInfo.gwas <- c("Arabidopsis thaliana GWAS", "Medicago truncatula GWAS")
  legumeInfo.organisms <- c("Arabidopsis thaliana", "Medicago truncatula")
  dataFiles <- c(dataFiles, legumeInfo.gwas)
  for(i in dataFiles){
    if (i %in% legumeInfo.gwas) {
      values[[i]] <- init.gwas(i)
    } else {
      values[[i]] <- read.table(paste0(dataPath,i),sep=",",stringsAsFactors=FALSE,head=TRUE)
    }
  }  
  values$datasetlist <- dataFiles
#  values[["ionomics"]] <- aggTable
#  values$datasetlist <- dataFiles
  values$datasetToOrganism <- NULL # map each dataset to an organism

  # what to do when the user changes the jth dataset selection
  datasetChanged <- function(j) {
    if (is.null(input[[jth_ref("datasets", j)]])) return()
    
    # Make sure the organism corresponds to the selected dataset
    # (invoke values$organism first to force a single reaction)
    if (is.null(values[[jth_ref("organism", j)]])) {
      values[[jth_ref("organism", j)]] = "Corn"
    } else {
      isolate({
        values[[jth_ref("organism", j)]] <- values$datasetToOrganism[[input[[jth_ref("datasets", j)]]]]
      })
    }
    
    # Uncheck the Append to Current Dataset checkbox and clear any previously selected GWAS files
    updateSelectizeInput(session, jth_ref("gwasTraits", j), choices = gwas.traits[[values[[jth_ref("organism", j)]]]], selected = NULL)
    # As there is no updateFileInput() method,
    # send a custom message (to clear the progress bar) and clear values$needsToUploadFiles
    session$sendCustomMessage(type = "resetFileInputHandler", jth_ref("uploadfile", j))
    values[[jth_ref("needsToUploadFiles", j)]] <- FALSE
    updateCheckboxInput(session, jth_ref("appendSNPs", j), value = FALSE)
    # Clear all genomic linkages if either organism changes
    values$glSelectedGene <- NULL
    values$glGenes <- values$glGenes2 <- values$glColors <- NULL
  }
  # This should be the first code block to detect a change in input$datasets
  observe(datasetChanged(1))
  # This should be the first code block to detect a change in input$datasets2
  observe(datasetChanged(2))

  # The user selected one or more local GWAS files (not yet loaded)
  gwasFilesSelected <- function(j) {
    input[[jth_ref("uploadfile", j)]]
    values[[jth_ref("needsToUploadFiles", j)]] <- TRUE
  }
  observe(gwasFilesSelected(1))
  observe(gwasFilesSelected(2))

  #handles what displays in the sidebar based on what tab is selected
  manageTabSelected <- "input.datatabs == 'Manage'"
  dataTableTabSelected <- "input.datatabs == 'Table'"
  wholeGenomeTabSelected <- "input.datatabs == 'WhGen'"
  #the second part of the statement [chromosomeTabSelected] is what is allowing the detection of changing panels due to a click event, i [original ZBrowse author] couldn't figure out how to update input$datatabs with javascript
  chromosomeTabSelected <- "input.datatabs == 'Chrom' || $('li.active a').first().html()==='Chromosome View'"
  annotationsTabSelected <- "input.datatabs == 'Annot'"
  wholeGenomeOrChromosomeTabSelected <- paste(wholeGenomeTabSelected, "||", chromosomeTabSelected)
  annotationsOrChromosomeTabSelected <- paste(annotationsTabSelected, "||", chromosomeTabSelected)

  createManageSidebar <- function(j) {
    list(
      wellPanel(
        style = paste0("background-color: ", bgColors[j], ";"),
        uiOutput(jth_ref("datasets", j))
      ),
      wellPanel(
        style = paste0("background-color: ", bgColors[j], ";"),
        radioButtons(inputId = jth_ref("dataType", j), label = "Load data (Max. 5MB):", c(".csv" = "csv", ".rda" = "rda", "examples" = "examples"), selected = "csv"),
        conditionalPanel(condition = paste0("input.", jth_ref("dataType", j), " != 'examples'"),
          conditionalPanel(condition = paste0("input.", jth_ref("dataType", j), " == 'csv'"),
            checkboxInput(jth_ref('header', j), 'Header', TRUE),
            radioButtons(jth_ref('sep', j), '', c(Comma=',', Semicolon=';', Tab='\t'), ',')
          ),
          selectizeInput(jth_ref("gwasTraits", j), "Remote Trait Files:", choices = NULL, multiple = TRUE),
          fileInput(jth_ref("uploadfile", j), "Local Trait Files:", multiple = TRUE),
          checkboxInput(jth_ref("appendSNPs", j), "Append to current dataset", FALSE),
          actionButton(jth_ref("loadTraits", j), "Load Data")
        ),
        conditionalPanel(condition = paste0("input.", jth_ref("dataType", j), " == 'examples'"),
          actionButton(jth_ref('loadExampleData', j), 'Load examples')
        )
      ),
      wellPanel(
        style = paste0("background-color: ", bgColors[j], ";"),
        h6("Once your file is finished uploading, press the Save Dataset button below and reload ZBrowse."),
        actionButton(jth_ref('saveDatasetButton', j), 'Save Current Dataset'),
        conditionalPanel(condition = paste0("input.", jth_ref("saveDatasetButton", j), " > 0"),
          h5("Dataset successfully saved!")
        )
      )
    )
  }
  
  createDataTableSidebar <- function(j) {
    wellPanel(
      style = paste0("background-color: ", bgColors[j], ";"),
      uiOutput(jth_ref("columns", j)),
      tags$br(),tags$br(),tags$br(),tags$br(),tags$br(),tags$br(), #add some space between selection columns and subset search
      # uiOutput("view_order"), checkboxInput("view_order_desc", "DESC", value = FALSE),
      returnTextInput(jth_ref("dv_select", j), "Subset (e.g., RMIP > 20 & Location == 'FL06')", '')
    )
  }
  
  createTraitsSidebar <- function(j) {
    wellPanel(
      style = paste0("background-color: ", bgColors[j], ";"),
      uiOutput(jth_ref("traitColBoxes", j)),
      uiOutput(jth_ref("legend", j)),
      uiOutput(jth_ref("overlaps", j)),
      conditionalPanel(condition = paste0("input.", jth_ref("overlaps", j), "==true"),
        uiOutput(jth_ref("overlapSize", j)),
        uiOutput(jth_ref("numOverlapTraits", j))
      )#end conditional for plotting only overlaps
      #submitButton("Update View"),
    )
  }
  createChromosomeSidebar <- function(j) {
    wellPanel(
      style = paste0("background-color: ", bgColors[j], ";"),
      uiOutput(jth_ref("selectChr", j))
    )
  }
  createAnnotationWindowSidebar <- function(j) {
    wellPanel(
      style = paste0("background-color: ", bgColors[j], ";"),
      h5("Annotation window options:"),
      h6("Click a point or type a basepair value:"),
      uiOutput(jth_ref("selectedOut", j)),
      uiOutput(jth_ref("windowOut", j))
    )
  }
  createDownloadAnnotationsSidebar <- function(j) {
    wellPanel(
      style = paste0("background-color: ", bgColors[j], ";"),
      helpText(h5(p(paste("Download a CSV of the annotations in the selected window.")))),
      downloadButton(jth_ref('downloadAnnot', j), 'Download')
    )
  }
  createGenomicLinkageSidebar <- function() {
    wellPanel(
      h5("Genomic Linkage options:"),
      checkboxInput('boolGenomicLinkage', 'ON', FALSE),
      conditionalPanel("input.boolGenomicLinkage == true",
        wellPanel(
          uiOutput("selectedGene"),
          numericInput("neighbors", "Neighbors:", min = 1, max = 20, value = 10),
          numericInput("matched", "Matched:", min = 1, max = 20, value = 6),
          numericInput("intermediate", "Intermediate:", min = 1, max = 10, value = 3),
          style = paste0("background-color: ", bgColors[1], ";")
        ),
        wellPanel(
          uiOutput("relatedRegions"),
          style = paste0("background-color: ", bgColors[2], ";")
        )
      )
    )
  }

  output$ui_All <- renderUI({
    list(
      conditionalPanel(condition = manageTabSelected,
        createManageSidebar(1),
        createManageSidebar(2),
        tags$script("Shiny.addCustomMessageHandler('resetFileInputHandler', function(x) {
          var id = '#' + x;
          var idProgress = id + '_progress';
          var idBar = id + ' .bar';
          $(idProgress).css('visibility', 'hidden');
          $(idBar).css('width', '0%');
          // TODO: clear the loaded file(s) (none of the following methods work, for Javascript security reasons)
          /*
          $(id).val(''); // method 1
          $(id).replaceWith($(id).val('').clone(true)); // method 2
          $(id).replaceWith($(id) = $(id).clone(true)); // method 3
          // method 4
          $(id).wrap('<form></form>').closest('form').reset();
          $(id).unwrap();
          $(id).stopPropagation();
          $(id).preventDefault();
          */
        });"),
        helpModal('Manage','manage',includeMarkdown("tools/manage.md")),HTML('<p style="font-size:10px;">Powered by <a href="http://www.rstudio.com/shiny/">Shiny</a>, <a href="http://rcharts.io/">rCharts</a> and <a href="http://www.highcharts.com">Highcharts</a></p>')             
      ),#end conditional Manage

      conditionalPanel(condition = dataTableTabSelected,
        createDataTableSidebar(1),
        createDataTableSidebar(2),
        helpModal('Data Table View','view',includeMarkdown("tools/manage.md"))      
      ),#end conditional Table

      # Sidebar panels for the remaining tabs, in conditional order
      conditionalPanel(condition = wholeGenomeOrChromosomeTabSelected,
        helpText(h5(p("Interactive Graphs for GWAS Data"))),
        createTraitsSidebar(1)
      ),
      conditionalPanel(condition = chromosomeTabSelected,
        createChromosomeSidebar(1)
      ),
      conditionalPanel(condition = annotationsOrChromosomeTabSelected,
        createAnnotationWindowSidebar(1)
      ),
      conditionalPanel(condition = annotationsTabSelected,
        createDownloadAnnotationsSidebar(1)
      ),
      conditionalPanel(condition = chromosomeTabSelected,
        createGenomicLinkageSidebar() # goes between the sidebars for organism 1 (above) and organism 2 (below)
      ),
      conditionalPanel(condition = wholeGenomeOrChromosomeTabSelected,
        createTraitsSidebar(2)
      ),
      conditionalPanel(condition = chromosomeTabSelected,
        createChromosomeSidebar(2)
      ),
      conditionalPanel(condition = annotationsOrChromosomeTabSelected,
        createAnnotationWindowSidebar(2)
      ),
      conditionalPanel(condition = annotationsTabSelected,
        createDownloadAnnotationsSidebar(2)
      ),
      # end sidebar panels for remaining tabs

      conditionalPanel(condition = wholeGenomeOrChromosomeTabSelected,
        helpModal('Browser Help','browser',includeMarkdown("tools/manage.md"))        
      )#add help button for browser tabs
    )#end list
  }) #end ui_All
  outputOptions(output, "ui_All", suspendWhenHidden=FALSE)

  # find the appropriate UI
  output$ui_finder <- renderUI({
#    if(is.null(input$datatabs)){      
#      get("ui_All")()
#    }else{
#      get(paste0('ui_',input$datatabs))()
#    }
#     if(input$tool == "data") {
#       if(!is.null(input$datatabs)) get(paste0('ui_',input$datatabs))()
#     } else {
#       if(!is.null(input$tool)) get(paste0('ui_',input$tool))()
#     }
  })  
  outputOptions(output, "ui_finder", suspendWhenHidden=FALSE)

  # create the menu for selecting the jth dataset
  renderDatasets <- function(j, defaultDataset) {
    if (is.null(input[[jth_ref("loadTraits", j)]])) return()
    
    # Load any requested GWAS files
    isolate({
      values[[jth_ref("uploadFiles", j)]] <- input[[jth_ref("uploadfile", j)]]
      values[[jth_ref("gwasTraits", j)]] <- input[[jth_ref("gwasTraits", j)]]
      
      # Local GWAS files
      inFile <- NULL
      if (values[[jth_ref("needsToUploadFiles", j)]]) inFile <- values[[jth_ref("uploadFiles", j)]]
      if (!is.null(inFile)) {
        # iterating through the files to upload
        for (i in 1:(dim(inFile)[1])) {
          loadUserData(inFile[i, 'name'], inFile[i, 'datapath'], j)
        }
        values[[jth_ref("needsToUploadFiles", j)]] <- FALSE # since we just loaded them
      }
      # Remote GWAS files
      inTraits <- values[[jth_ref("gwasTraits", j)]]
      if (!is.null(inTraits)) {
        for (trait in inTraits) {
          trait.url <- gwas.filenames[[values[[jth_ref("organism", j)]]]][which(gwas.traits[[values[[jth_ref("organism", j)]]]] == trait)]
          loadRemoteData(trait, trait.url, j)
        }
      }
    })
    
    dat <- isolate(input[[jth_ref("datasets", j)]])
    if (!is.null(dat)) {
      appendSNPs <- isolate(input[[jth_ref("appendSNPs", j)]])
      if (appendSNPs || (is.null(inFile) && is.null(inTraits))) {
        val <- dat
      } else {
        val <- values$datasetlist[1]
      }
    }else{
      val <- defaultDataset
    }
    
    # Drop-down selection of data set
    selectInput(inputId = jth_ref("datasets", j), label = paste0("Dataset ", j, ":"), choices = values$datasetlist, selected = values$datasetlist[values$datasetlist == val], multiple = FALSE, selectize = FALSE)
  }
  output$datasets <- renderUI(renderDatasets(1, "Medicago truncatula GWAS"))
  output$datasets2 <- renderUI(renderDatasets(2, "Arabidopsis thaliana GWAS"))
  
  createAnnotTable <- function(j) {
    if (is.null(input[[jth_ref("datasets", j)]])) return()
    centerBP <- as.numeric(input[[jth_ref("selected", j)]][[1]])
    winHigh <- centerBP + input[[jth_ref("window", j)]][1]
    winLow <- centerBP - input[[jth_ref("window", j)]][1]
    if (winLow < 0) { winLow <- 0 }
    thisChrAnnot <- subset(annotGeneLoc[values[[jth_ref("organism", j)]]][[1]], chromosome == input[[jth_ref("chr", j)]])
    thisAnnot <- thisChrAnnot[thisChrAnnot$transcript_start >= winLow & thisChrAnnot$transcript_end <= winHigh, ]
    thisAnnot
  }
  
  #Returns the nicely formatted preview table
  createHtmlDataExample <- function(j) {
    if(is.null(input[[jth_ref("datasets", j)]])) return()

    dat <- getdata(j)

    # necessary when deleting a dataset
    if(is.null(dat) || nrow(dat) == 0) return()

    # Show only the first 10 rows
    nr <- min(10,nrow(dat))
    dat <- data.frame(dat[1:nr,, drop = FALSE])

    #dat <- date2character_dat(dat) #may be needed to print table if there is a data column

    html <- print(xtable::xtable(dat), type='html', print.results = FALSE)
    html <- sub("<TABLE border=1>","<table class='table table-condensed table-hover'>", html)
    Encoding(html) <- 'UTF-8'
    html
  }
  output$htmlDataExample <- reactive(createHtmlDataExample(1))
  output$htmlDataExample2 <- reactive(createHtmlDataExample(2))
  
  createNrowDataset <- function(j) {
    if(is.null(input[[jth_ref("datasets", j)]])) return()
    dat <- getdata(j)
    if(is.null(dat) || nrow(dat) == 0) return()
    nr <- nrow(dat)
    
    if(nr>2500){
      paste0('<p>First 10 rows shown of ',nr,' rows. See Data Table tab for details.<br>More than 2500 rows found, only the top 2500 will be plotted.</p>')
    }else{
      paste0('<p>First 10 rows shown of ',nr,' rows. See Data Table tab for details.</p>')
    }
  }
  output$nrowDataset <- reactive(createNrowDataset(1))
  output$nrowDataset2 <- reactive(createNrowDataset(2))

  createColumnSettingsPanel <- function(j) {
    tags$div(
      class = "container",
      style = paste0("background-color: ", bgColors[j], ";"),

      htmlOutput(jth_ref("htmlDataExample", j)),
      htmlOutput(jth_ref("nrowDataset", j)),

      #tags$p(tags$br()),
      row(
        col(3, tags$br()),
        col(7, h4('Select appropriate columns to be used for plotting.'))
        #HTML('<h4>Select appropriate columns to be used for plotting.</h4>'),
      ),
      tags$hr(),
      row(
        #col(2, tags$br()),
        col(2,uiOutput(jth_ref("chrColumn", j)),uiOutput(jth_ref("bpColumn", j))),
        col(2,uiOutput(jth_ref("plotAll", j)),uiOutput(jth_ref("traitColumns", j))),
        col(2,uiOutput(jth_ref("yAxisColumn", j)),uiOutput(jth_ref("logP", j))),
        col(2,uiOutput(jth_ref("axisLimBool", j)),uiOutput(jth_ref("axisLim", j))),
        col(2,actionButton(jth_ref("SubmitColsButton", j),"Submit"))
      ),
      tags$hr(),
      row(
        col(7, uiOutput(jth_ref("supportInterval", j)))#
      ),
      row(
        col(2, uiOutput(jth_ref("SIbpStart", j))),
        col(2, uiOutput(jth_ref("SIyAxisColumn", j))),
        col(2, uiOutput(jth_ref("SIaxisLimBool", j)),uiOutput(jth_ref("SIaxisLim", j)))
      )
    )
  }
  output$ui_data_tabs <- renderUI({
    tabsetPanel(id = "datatabs",      
      tabPanel(title="Manage",value="Manage",
        createColumnSettingsPanel(1),
        createColumnSettingsPanel(2)
      ),
      tabPanel(title="Data Table",value="Table",
        wellPanel(dataTableOutput("dataviewer"), style = paste0("background-color: ", bgColors[1], ";")),
        wellPanel(dataTableOutput("dataviewer2"), style = paste0("background-color: ", bgColors[2], ";"))
      ),
      tabPanel(title="Whole Genome View",value="WhGen",
        wellPanel(showOutput("gChart", "highcharts"), style = paste0("background-color: ", bgColors[1], ";")),
        wellPanel(showOutput("gChart2", "highcharts"), style = paste0("background-color: ", bgColors[2], ";"))
      ),
      tabPanel(title="Chromosome View",value="Chrom",
        wellPanel(showOutput("pChart", "highcharts"), showOutput("zChart", "highcharts"),
          tags$script('Shiny.addCustomMessageHandler("customMsg", function(bandOpts){
            chartXAxis = $("#pChart").highcharts().xAxis[0]
            chartXAxis.removePlotBand()
            chartXAxis.addPlotBand(bandOpts)
          })'),
          style = paste0("background-color: ", bgColors[1], ";")
        ),
        wellPanel(showOutput("zChart2", "highcharts"), showOutput("pChart2", "highcharts"),
          tags$script('Shiny.addCustomMessageHandler("customMsg2", function(bandOpts){
            chartXAxis = $("#pChart2").highcharts().xAxis[0]
            chartXAxis.removePlotBand()
            chartXAxis.addPlotBand(bandOpts)
          })'),
          style = paste0("background-color: ", bgColors[2], ";")
        )
      ),
      tabPanel(title="Annotations Table",value="Annot",
        wellPanel(dataTableOutput("annotViewer"), style = paste0("background-color: ", bgColors[1], ";")),
        wellPanel(dataTableOutput("annotViewer2"), style = paste0("background-color: ", bgColors[2], ";"))
      )
    )#end tabsetPanel
  })#end data tabs
  outputOptions(output, "ui_data_tabs", suspendWhenHidden=FALSE)

  annotViewer.options <- list(orderClasses = TRUE, bCaseInsensitive = TRUE,
    lengthMenu = c(15, 50, 100, 200, 500), pageLength = 15,
    "dom" = 'T<"clear">lfrtip',
    "oTableTools" = list(
      "sSwfPath" = "/tabletools/swf/copy_csv_xls_pdf.swf",
      "aButtons" = list(
        "copy",
        "print",
        list("sExtends" = "collection",
          "sButtonText" = "Save",
          "aButtons" = c("csv","xls","pdf")
        )
      )
    )                    
  )
  output$annotViewer <- renderDataTable({
    createAnnotTable(1)
  }, options = annotViewer.options)
  output$annotViewer2 <- renderDataTable({
    createAnnotTable(2)
  }, options = annotViewer.options)

  dataViewer.options <- annotViewer.options # they happen to be identical, but we could make them different
  createDataViewer <- function(j) {
    if(is.null(input[[jth_ref("datasets", j)]]) || is.null(input[[jth_ref("columns", j)]])) return()
    
    dat <- getdata(j)
    #dat <- date2character()
    
    if(!all(input[[jth_ref("columns", j)]] %in% colnames(dat))) return()
    
    if(input[[jth_ref("dv_select", j)]] != '') {
      selcom <- input[[jth_ref("dv_select", j)]]
      selcom <- gsub(" ", "", selcom)
      
      seldat <- try(do.call(subset, list(dat,parse(text = selcom))), silent = TRUE)
      
      if(!is(seldat, 'try-error')) {
        if(is.data.frame(seldat)) {
          dat <- seldat
          seldat <- NULL
        }
      }
    }
    
    dat <- data.frame(dat[, input[[jth_ref("columns", j)]], drop = FALSE])
    dat
    
    # html <- print(xtable::xtable(dat), type='html', print.results = FALSE)
    # html <- sub("<TABLE border=1>","<table class='table table-condensed table-hover'>", html)
    # html
  }
  output$dataviewer <-renderDataTable(createDataViewer(1), options = dataViewer.options)
  output$dataviewer2 <-renderDataTable(createDataViewer(2), options = dataViewer.options)
  
  createDownloadAnnot <- function(j) {
    downloadHandler(
      filename = function() {paste0("AnnotationsAround.chr",input[[jth_ref("chr", j)]],".",input[[jth_ref("selected", j)]][[1]],"bp.",values[[jth_ref("organism", j)]],".csv")},
      content = function(file) {write.csv(createAnnotTable(j),file,row.names=F)}
    )
  }
  output$downloadAnnot <- createDownloadAnnot(1)
  output$downloadAnnot2 <- createDownloadAnnot(2)

  createColumns <- function(j) {
    cols <- varnames(j)    
    selectInput(jth_ref("columns", j), "Select columns to show:", choices = as.list(cols), selected = cols, multiple = TRUE)
  }
  output$columns <- renderUI(createColumns(1))
  output$columns2 <- renderUI(createColumns(2))
  
  createAxisLimBool <- function(j) {
    if(is.null(input[[jth_ref("datasets", j)]])){return()}
    if(input[[jth_ref("datasets", j)]] %in% datasetProp()$dataset){
      val = datasetProp()$axisLim[datasetProp()$dataset == input[[jth_ref("datasets", j)]]]
    } else if (values[[jth_ref("organism", j)]] %in% legumeInfo.organisms) {
      val <- FALSE
    }else{
      val = TRUE}
    checkboxInput(jth_ref('axisLimBool', j), 'Set Y-axis Limits?', val)
  }
  output$axisLimBool <- renderUI(createAxisLimBool(1))
  output$axisLimBool2 <- renderUI(createAxisLimBool(2))

  createLogP <- function(j) {
    if(is.null(input[[jth_ref("datasets", j)]])){return()}
    if(input[[jth_ref("datasets", j)]] %in% datasetProp()$dataset){
      val = datasetProp()$logP[datasetProp()$dataset == input[[jth_ref("datasets", j)]]]
    } else if (values[[jth_ref("organism", j)]] %in% legumeInfo.organisms) {
      val <- TRUE
    }else{
      val = FALSE}
    checkboxInput(jth_ref('logP', j), 'Take -log10 of column?', val)
  }
  output$logP <- renderUI(createLogP(1))
  output$logP2 <- renderUI(createLogP(2))

  createChrColumn <- function(j) {
    if(is.null(input[[jth_ref("datasets", j)]])){return()}
    cols <- varnames(j)    
    if(input[[jth_ref("datasets", j)]] %in% datasetProp()$dataset){
      selected = datasetProp()$chrColumn[datasetProp()$dataset == input[[jth_ref("datasets", j)]]]
    }else{
      selected = names(cols[1])
    }    
    selected <- stri_match(selected, regex = ".*(?=\\ \\{)")[, 1]
    selectizeInput(jth_ref("chrColumn", j), "Chromosome Column:", choices = as.list(cols), selected = selected, multiple = FALSE, options = list(dropdownParent="body"))
  }
  output$chrColumn <- renderUI(createChrColumn(1))
  output$chrColumn2 <- renderUI(createChrColumn(2))

  createBpColumn <- function(j) {
    if(is.null(input[[jth_ref("datasets", j)]])){return()}
    cols <- varnames(j)
    if(input[[jth_ref("datasets", j)]] %in% datasetProp()$dataset){
      selected = datasetProp()$bpColumn[datasetProp()$dataset == input[[jth_ref("datasets", j)]]]
    }else{
      selected = names(cols[2])
    }
    selected <- stri_match(selected, regex = ".*(?=\\ \\{)")[, 1]
    selectizeInput(jth_ref("bpColumn", j), "Base Pair Column:", choices = as.list(cols), selected = selected, multiple = FALSE, options = list(dropdownParent="body"))
  }
  output$bpColumn <- renderUI(createBpColumn(1))
  output$bpColumn2 <- renderUI(createBpColumn(2))

  createTraitColumns <- function(j) {
    if(is.null(input[[jth_ref("datasets", j)]])){return()}
    cols <- varnames(j)
    if(input[[jth_ref("datasets", j)]] %in% datasetProp()$dataset){
      selected = unlist(strsplit(datasetProp()$traitCol[datasetProp()$dataset == input[[jth_ref("datasets", j)]]],";"))
    } else if (values[[jth_ref("organism", j)]] %in% legumeInfo.organisms) {
      selected <- names(cols[3])
    }else{
      selected = names(cols[3:4])
    }
    selected <- sapply(selected, FUN = function(s) stri_match(selected, regex = ".*(?=\\ \\{)")[, 1], USE.NAMES = FALSE)
    conditionalPanel(condition = paste0("input.", jth_ref("plotAll", j), "==false"),
      selectizeInput(jth_ref("traitColumns", j), "Group by these trait column(s):", choices = as.list(cols), selected = selected, multiple = TRUE, options = list(dropdownParent="body"))
    )        
  }
  output$traitColumns <- renderUI(createTraitColumns(1))
  output$traitColumns2 <- renderUI(createTraitColumns(2))

  createPlotAll <- function(j) {
    if(is.null(input[[jth_ref("datasets", j)]])){return()}
    if(input[[jth_ref("datasets", j)]] %in% datasetProp()$dataset){      
      val = datasetProp()$plotAll[datasetProp()$dataset == input[[jth_ref("datasets", j)]]]
    }else{
      val = FALSE
    }    
    checkboxInput(jth_ref('plotAll', j), 'All data is the same trait', val)
  }
  output$plotAll <- renderUI(createPlotAll(1))
  output$plotAll2 <- renderUI(createPlotAll(2))

  createYAxisColumn <- function(j) {
    if(is.null(input[[jth_ref("datasets", j)]])){return()}
    cols <- varnames(j)
    if(input[[jth_ref("datasets", j)]] %in% datasetProp()$dataset){      
      selected = datasetProp()$yAxisColumn[datasetProp()$dataset == input[[jth_ref("datasets", j)]]]
    } else if (values[[jth_ref("organism", j)]] %in% legumeInfo.organisms) {
      selected <- names(cols[4])
    }else{
      #selected = names(cols[10])
      selected = as.character(cols[10])
    }    
    selected <- stri_match(selected, regex = ".*(?=\\ \\{)")[, 1]
    selectizeInput(jth_ref("yAxisColumn", j), "Y-axis column:", choices = as.list(cols), selected = selected, multiple = FALSE, options = list(dropdownParent="body"))
  }
  output$yAxisColumn <- renderUI(createYAxisColumn(1))
  outputOptions(output, "yAxisColumn", suspendWhenHidden=FALSE)
  output$yAxisColumn2 <- renderUI(createYAxisColumn(2))
  outputOptions(output, "yAxisColumn2", suspendWhenHidden=FALSE)
  
  createAxisLim <- function(j) {
    if(is.null(input[[jth_ref("datasets", j)]])){return()}
    if(input[[jth_ref("datasets", j)]] %in% datasetProp()$dataset){
      min = datasetProp()$axisMin[datasetProp()$dataset == input[[jth_ref("datasets", j)]]]
      max = datasetProp()$axisMax[datasetProp()$dataset == input[[jth_ref("datasets", j)]]]
    }else{
      min = 0
      max = 1
    }
    conditionalPanel(condition = paste0("input.", jth_ref("axisLimBool", j), "==true"),
      numericInput(jth_ref("axisMin", j), "Min:", value=min),
      numericInput(jth_ref("axisMax", j), "Max:", value=max)
    )
  }
  output$axisLim <- renderUI(createAxisLim(1))
  output$axisLim2 <- renderUI(createAxisLim(2))

  createSupportInterval <- function(j) {
    if(is.null(input[[jth_ref("datasets", j)]])){return()}
    if(input[[jth_ref("datasets", j)]] %in% datasetProp()$dataset){
      val = datasetProp()$supportInterval[datasetProp()$dataset == input[[jth_ref("datasets", j)]]]
    }else{
      val = FALSE}
    checkboxInput(jth_ref('supportInterval', j), 'Plot base pair intervals (e.g., Joint linkage support intervals)?', val)
  }
  output$supportInterval <- renderUI(createSupportInterval(1))
  output$supportInterval2 <- renderUI(createSupportInterval(2))

  createSIbpStart <- function(j) {
    if(is.null(input[[jth_ref("datasets", j)]])){return()}
    cols <- varnames(j)
    if(input[[jth_ref("datasets", j)]] %in% datasetProp()$dataset){
      selected = datasetProp()$SIbpStart[datasetProp()$dataset == input[[jth_ref("datasets", j)]]]
      selectedEnd = datasetProp()$SIbpEnd[datasetProp()$dataset == input[[jth_ref("datasets", j)]]]
    }else{
      selected = names(cols[2])
      selectedEnd = names(cols[2])
    }
    conditionalPanel(condition = paste0("input.", jth_ref("supportInterval", j), "==true"),
      selectizeInput(jth_ref("SIbpStart", j), "Interval Base Pair Start:", choices = as.list(cols), selected = selected, multiple = FALSE, options = list(dropdownParent="body")),
      selectizeInput(jth_ref("SIbpEnd", j), "Interval Base Pair End:", choices = as.list(cols), selected = selectedEnd, multiple = FALSE, options = list(dropdownParent="body"))
    )
  }
  output$SIbpStart <- renderUI(createSIbpStart(1))
  output$SIbpStart2 <- renderUI(createSIbpStart(2))

  createSIyAxisColumn <- function(j) {
    if(is.null(input[[jth_ref("datasets", j)]])){return()}
    cols <- varnames(j)
    if(input[[jth_ref("datasets", j)]] %in% datasetProp()$dataset){  
      selected = datasetProp()$SIyAxisColumn[datasetProp()$dataset == input[[jth_ref("datasets", j)]]]
    }else{
      #selected = names(cols[10])
      selected = as.character(cols[10])
    }        
    conditionalPanel(condition = paste0("input.", jth_ref("supportInterval", j), "==true"),
      selectizeInput(jth_ref("SIyAxisColumn", j), "Support Interval Y-axis column:", choices = as.list(cols), selected = selected, multiple = FALSE, options = list(dropdownParent="body"))
    )
  }
  output$SIyAxisColumn <- renderUI(createSIyAxisColumn(1))
  outputOptions(output, "SIyAxisColumn", suspendWhenHidden=FALSE)  
  output$SIyAxisColumn2 <- renderUI(createSIyAxisColumn(2))
  outputOptions(output, "SIyAxisColumn2", suspendWhenHidden=FALSE)
  
  createSIaxisLimBool <- function(j) {
    if(is.null(input[[jth_ref("datasets", j)]])){return()}
    if(input[[jth_ref("datasets", j)]] %in% datasetProp()$dataset){
      val = datasetProp()$SIaxisLimBool[datasetProp()$dataset == input[[jth_ref("datasets", j)]]]
    } else if (values[[jth_ref("organism", j)]] %in% legumeInfo.organisms) {
      val <- FALSE
    }else{
      val = TRUE
    }
    conditionalPanel(condition = paste0("input.", jth_ref("supportInterval", j), "==true"),
      checkboxInput(jth_ref('SIaxisLimBool', j), 'Set Support Interval Y-axis Limits?', val)
    )
  }
  output$SIaxisLimBool <- renderUI(createSIaxisLimBool(1))
  output$SIaxisLimBool2 <- renderUI(createSIaxisLimBool(2))

  createSIaxisLim <- function(j) {
    if(is.null(input[[jth_ref("datasets", j)]])){return()}
    if(input[[jth_ref("datasets", j)]] %in% datasetProp()$dataset){
      min = datasetProp()$axisMin[datasetProp()$dataset == input[[jth_ref("datasets", j)]]]
      max = datasetProp()$axisMax[datasetProp()$dataset == input[[jth_ref("datasets", j)]]]
    }else{
      min = 0
      max = 1
    }
    conditionalPanel(condition = paste0("input.", jth_ref("supportInterval", j), "==true && input.", jth_ref("SIaxisLimBool", j), "==true"),
      numericInput(jth_ref("SIaxisMin", j), "Min:", value=min),
      numericInput(jth_ref("SIaxisMax", j), "Max:", value=max)
    )
  }
  output$SIaxisLim <- renderUI(createSIaxisLim(1))
  output$SIaxisLim2 <- renderUI(createSIaxisLim(2))

  #builds list of multiple selection boxes for traits that have multiple columns in dataset
  createTraitColBoxes <- function(j) {
    if(is.null(input[[jth_ref("datasets", j)]])){return()}
    if(input[[jth_ref("plotAll", j)]] == TRUE){return()}
    lapply(input[[jth_ref("traitColumns", j)]], function(i) {
      traits <- c("Select All ",sort(unique(values[[input[[jth_ref("datasets", j)]]]][,i])))
      selectizeInput(inputId=jth_ref(i, j), label=paste0("Select ",i),traits,
        selected=traits[2],
        multiple=TRUE, options = list(dropdownParent="body",plugins=list("remove_button")))
    })
  }
  output$traitColBoxes <- renderUI(createTraitColBoxes(1))
  output$traitColBoxes2 <- renderUI(createTraitColBoxes(2))

  updateTraitsMenu <- function(j) {
    lapply(input[[jth_ref("traitColumns", j)]], function(i){
      if("Select All " %in% input[[jth_ref(i, j)]]){
        selected_choices <- sort(unique(values[[input[[jth_ref("datasets", j)]]]][,i]))
        updateSelectizeInput(session, jth_ref(i, j), selected = selected_choices)
      }
    })
  }
  observe(updateTraitsMenu(1))
  observe(updateTraitsMenu(2))

  #checkbox to suppress plot legend
  createLegend <- function(j) {
    if(is.null(input[[jth_ref("datasets", j)]])){return()}
    checkboxInput(jth_ref('legend', j), 'Suppress Legend', FALSE)
  }
  output$legend <- renderUI(createLegend(1))
  outputOptions(output, "legend", suspendWhenHidden=FALSE)
  output$legend2 <- renderUI(createLegend(2))
  outputOptions(output, "legend2", suspendWhenHidden=FALSE)
  
  #checkbox for whether to filter for only overlapping SNPs
  createOverlaps <- function(j) {
    if(is.null(input[[jth_ref("datasets", j)]])){return()}
#    if(input[[jth_ref("plotAll", j)]] == TRUE){return()}
    checkboxInput(jth_ref('overlaps', j), 'Show only overlapping SNPs', FALSE)
  }
  output$overlaps <- renderUI(createOverlaps(1))
  outputOptions(output, "overlaps", suspendWhenHidden=FALSE)
  output$overlaps2 <- renderUI(createOverlaps(2))
  outputOptions(output, "overlaps2", suspendWhenHidden=FALSE)
  
  #how many traits must overlap to be included in output, 1 means traits that overlap with themselves will be included
  createNumOverlapTraits <- function(j) {
    numericInput(jth_ref("numOverlaps", j), "Minimum number of overlaps?", value=2,min=1)
  }
  output$numOverlapTraits <- renderUI(createNumOverlapTraits(1))
  output$numOverlapTraits2 <- renderUI(createNumOverlapTraits(2))

  #how big is the window when calculating whether two snps overlap
  createOverlapSize <- function(j) {
    numericInput(inputId=jth_ref("overlapSize", j), label="Overlap size around each point:",min=1,max=.5e6,value=10000)
  }
  output$overlapSize <- renderUI(createOverlapSize(1))
  output$overlapSize2 <- renderUI(createOverlapSize(2))

  createSelectChr <- function(j) {
    if (is.null(values[[jth_ref("organism", j)]])) { return() }
    selectInput(jth_ref("chr", j), "Chromosome:", chrName[values[[jth_ref("organism", j)]]][[1]], selectize = FALSE)
  }
  output$selectChr <- renderUI(createSelectChr(1))
  outputOptions(output, "selectChr", suspendWhenHidden=FALSE)
  output$selectChr2 <- renderUI(createSelectChr(2))
  outputOptions(output, "selectChr2", suspendWhenHidden=FALSE)
  
  createSelectedOut <- function(j) {
    numericInput(jth_ref("selected", j), "", value=100000)
  }
  output$selectedOut <- renderUI(createSelectedOut(1))
  outputOptions(output, "selectedOut", suspendWhenHidden=FALSE)
  output$selectedOut2 <- renderUI(createSelectedOut(2))
  outputOptions(output, "selectedOut2", suspendWhenHidden=FALSE)

  createWindowOut <- function(j) {
    sliderInput(inputId=jth_ref("window", j), label="Window size around selected point:",min=1000,max=.5e6,value=2.5e5)
  }
  output$windowOut <- renderUI(createWindowOut(1))
  outputOptions(output, "windowOut", suspendWhenHidden=FALSE)
  output$windowOut2 <- renderUI(createWindowOut(2))
  outputOptions(output, "windowOut2", suspendWhenHidden=FALSE)
  
  #returns datasets from uploaded file
  getdata <- function(j) {
    if (is.null(input[[jth_ref("datasets", j)]])) { return() }
    values[[input[[jth_ref("datasets", j)]]]]
  }

  #builds list of column names and type in dataset
  varnames <- function(j) {
    if (is.null(input[[jth_ref("datasets", j)]])) return()
    dat <- getdata_class(j)
    vars <- names(dat)
    names(vars) <- paste(vars, " {", dat, "}", sep = "")
    vars
  }

  getdata_class <- function(j) {
    if (is.null(input[[jth_ref("datasets", j)]])) return()
    cls <- sapply(getdata(j), function(x) class(x)[1])
    gsub("ordered","factor", cls)
  }

  # Function to handle loading of data from a file or rObject, for the jth dataset
  loadUserData <- function(filename, uFile, j) {
    ext <- file_ext(filename)
    objname <- sub(paste(".", ext, sep = ""), "", basename(filename))
    ext <- tolower(ext)

    if (ext == 'rda' || ext == 'rdata') {
      # objname will hold the name of the object(s) inside the R datafile
      robjname <- load(uFile)
      if (length(robjname) > 1) {
        loaded.values <- data.frame(get(robjname[1]))
      } else {
        loaded.values <- data.frame(get(robjname))  # only work with data.frames
      }
    }

    appendSNPs <- isolate(input[[jth_ref("appendSNPs", j)]])
    if (!appendSNPs) {
      # add new datasets to the datasetToOrganism map
      values$datasetToOrganism[[objname]] <- values[[jth_ref("organism", j)]]
    }
    if (length(values[['datasetlist']]) == 0 || values[['datasetlist']][1] == '') {
      values[['datasetlist']] <- c(objname)
    } else if (!appendSNPs) {
      values[['datasetlist']] <- unique(c(objname, values[['datasetlist']]))
    }

    if (ext == 'sav') {
      loaded.values <- as.data.frame(as.data.set(spss.system.file(uFile)))
    } else if (ext == 'dta') {
      loaded.values <- read.dta(uFile)
    } else {
      loaded.values <- read.csv(uFile, header = input[[jth_ref("header", j)]], sep = input[[jth_ref("sep", j)]], stringsAsFactors = FALSE)
    }

    if (appendSNPs) {
      values[[input[[jth_ref("datasets", j)]]]]$totalBP <- NULL
      names(loaded.values) <- names(values[[input[[jth_ref("datasets", j)]]]])
      values[[input[[jth_ref("datasets", j)]]]] <- rbind(values[[input[[jth_ref("datasets", j)]]]], loaded.values)
    } else {
      values[[objname]] <- loaded.values
    }
  }
  
  # Load data from .csv files at a remote URL, for the jth dataset
  loadRemoteData <- function(trait, traitUrl, j) {
    ext <- file_ext(traitUrl)
    objname <- sub(paste(".", ext, sep = ""), "", basename(traitUrl))

    appendSNPs <- isolate(input[[jth_ref("appendSNPs", j)]])
    if (!appendSNPs) {
      # add new datasets to the datasetToOrganism map
      values$datasetToOrganism[[objname]] <- values[[jth_ref("organism", j)]]
    }
    if (length(values[['datasetlist']]) == 0 || values[['datasetlist']][1] == '') {
      values[['datasetlist']] <- c(objname)
    } else if (!appendSNPs) {
      values[['datasetlist']] <- unique(c(objname, values[['datasetlist']]))
    }

    loaded.values <- load.gwas.remote(values[[jth_ref("organism", j)]], traitUrl, trait)

    if (appendSNPs) {
      values[[input[[jth_ref("datasets", j)]]]]$totalBP <- NULL
      names(loaded.values) <- names(values[[input[[jth_ref("datasets", j)]]]])
      values[[input[[jth_ref("datasets", j)]]]] <- rbind(values[[input[[jth_ref("datasets", j)]]]], loaded.values)
    } else {
      values[[objname]] <- loaded.values
    }
  }

  #add a totalBP column to an input dataset if not already present
  calculateTotalBP <- function(j) {
    if("totalBP" %in% colnames(values[[input[[jth_ref("datasets", j)]]]])){
      
    }else{
#      progress <- Progress$new(session, min=1, max=1)
#      on.exit(progress$close())
      
#      progress$set(message = 'Preparing data for plotting.',
#                   detail = 'This may take a few seconds...')
#       withProgress(session, min=1, max=400, expr={
#         for(i in 1:400) {
#           setProgress(message = 'Preparing Data for Plotting',
#                       detail = 'This may take a few seconds...',
#                       value=i)
#           #print(i)
#           Sys.sleep(0.1)
#         }
#       })     
      cumBP<-c(0,cumsum(as.numeric(chrSize[values[[jth_ref("organism", j)]]][[1]])))
      #to order by desired chromosome add factor levels in the desired order to the chrColumn, any chr names that differ in gwas file compared
      #to organism file will turn into NA
      values[[input[[jth_ref("datasets", j)]]]][,input[[jth_ref("chrColumn", j)]]] <- factor(values[[input[[jth_ref("datasets", j)]]]][,input[[jth_ref("chrColumn", j)]]],levels=chrName[values[[jth_ref("organism", j)]]][[1]])
      values[[input[[jth_ref("datasets", j)]]]] <- values[[input[[jth_ref("datasets", j)]]]][order(values[[input[[jth_ref("datasets", j)]]]][,input[[jth_ref("chrColumn", j)]]],values[[input[[jth_ref("datasets", j)]]]][,input[[jth_ref("bpColumn", j)]]]),]
      numeachchr<-aggregate(values[[input[[jth_ref("datasets", j)]]]][,input[[jth_ref("bpColumn", j)]]],list(values[[input[[jth_ref("datasets", j)]]]][,input[[jth_ref("chrColumn", j)]]]),length)
#      adjust<-rep(cumBP[1],numeachchr$x[numeachchr$Group.1==1])            
      adjust <- numeric()
      for (i in 1:(length(cumBP)-1)){#max(unique(values[[input[[jth_ref("datasets", j)]]]][,input[[jth_ref("chrColumn", j)]]]))){
        if(length(numeachchr$x[numeachchr$Group.1==chrName[values[[jth_ref("organism", j)]]][[1]][i]])==0){next;}
        adjust<-c(adjust,rep(cumBP[i],numeachchr$x[numeachchr$Group.1==chrName[values[[jth_ref("organism", j)]]][[1]][i]]))
      }
      #newval <- values[[input[[jth_ref("datasets", j)]]]][600,input[[jth_ref("bpColumn", j)]]]+adjust[600]      
      values[[input[[jth_ref("datasets", j)]]]]$totalBP <- values[[input[[jth_ref("datasets", j)]]]][,input[[jth_ref("bpColumn", j)]]]+adjust
      
      #values[[input[[jth_ref("datasets", j)]]]] <- adply(values[[input[[jth_ref("datasets", j)]]]],1,function(x){data.frame(totalBP=sum(x[[input[[jth_ref("bpColumn", j)]]]],chrSize$bp[chrSize$chr %in% if(x[[input[[jth_ref("chrColumn", j)]]]]==1) 0 else c(1:(x[[input[[jth_ref("chrColumn", j)]]]]-1))]))})
    }
   if(input[[jth_ref("supportInterval", j)]] == TRUE){
      if("SIbpStartTotal" %in% colnames(values[[input[[jth_ref("datasets", j)]]]])){
        
      }else{
         
        cumBP<-c(0,cumsum(as.numeric(chrSize[values[[jth_ref("organism", j)]]][[1]])))
        values[[input[[jth_ref("datasets", j)]]]][,input[[jth_ref("chrColumn", j)]]] <- factor(values[[input[[jth_ref("datasets", j)]]]][,input[[jth_ref("chrColumn", j)]]],levels=chrName[values[[jth_ref("organism", j)]]][[1]])
        values[[input[[jth_ref("datasets", j)]]]] <- values[[input[[jth_ref("datasets", j)]]]][order(values[[input[[jth_ref("datasets", j)]]]][,input[[jth_ref("chrColumn", j)]]],values[[input[[jth_ref("datasets", j)]]]][,input[[jth_ref("SIbpStart", j)]]]),]
        numeachchr<-aggregate(values[[input[[jth_ref("datasets", j)]]]][,input[[jth_ref("SIbpStart", j)]]],list(values[[input[[jth_ref("datasets", j)]]]][,input[[jth_ref("chrColumn", j)]]]),length)
        adjust <- numeric()
        for (i in 1:(length(cumBP)-1)){#max(unique(values[[input[[jth_ref("datasets", j)]]]][,input[[jth_ref("chrColumn", j)]]]))){
          if(length(numeachchr$x[numeachchr$Group.1==chrName[values[[jth_ref("organism", j)]]][[1]][i]])==0){next;}
          adjust<-c(adjust,rep(cumBP[i],numeachchr$x[numeachchr$Group.1==chrName[values[[jth_ref("organism", j)]]][[1]][i]]))
        }
        values[[input[[jth_ref("datasets", j)]]]]$SIbpStartTotal <- values[[input[[jth_ref("datasets", j)]]]][,input[[jth_ref("SIbpStart", j)]]]+adjust    
      }
    
      if("SIbpEndTotal" %in% colnames(values[[input[[jth_ref("datasets", j)]]]])){
        
      }else{
        
        cumBP<-c(0,cumsum(as.numeric(chrSize[values[[jth_ref("organism", j)]]][[1]])))
        values[[input[[jth_ref("datasets", j)]]]][,input[[jth_ref("chrColumn", j)]]] <- factor(values[[input[[jth_ref("datasets", j)]]]][,input[[jth_ref("chrColumn", j)]]],levels=chrName[values[[jth_ref("organism", j)]]][[1]])
        values[[input[[jth_ref("datasets", j)]]]] <- values[[input[[jth_ref("datasets", j)]]]][order(values[[input[[jth_ref("datasets", j)]]]][,input[[jth_ref("chrColumn", j)]]],values[[input[[jth_ref("datasets", j)]]]][,input[[jth_ref("SIbpEnd", j)]]]),]
        numeachchr<-aggregate(values[[input[[jth_ref("datasets", j)]]]][,input[[jth_ref("SIbpEnd", j)]]],list(values[[input[[jth_ref("datasets", j)]]]][,input[[jth_ref("chrColumn", j)]]]),length)
        adjust <- numeric()
        for (i in 1:(length(cumBP)-1)){#max(unique(values[[input[[jth_ref("datasets", j)]]]][,input[[jth_ref("chrColumn", j)]]]))){
          if(length(numeachchr$x[numeachchr$Group.1==chrName[values[[jth_ref("organism", j)]]][[1]][i]])==0){next;}
          adjust<-c(adjust,rep(cumBP[i],numeachchr$x[numeachchr$Group.1==chrName[values[[jth_ref("organism", j)]]][[1]][i]]))
        }
        values[[input[[jth_ref("datasets", j)]]]]$SIbpEndTotal <- values[[input[[jth_ref("datasets", j)]]]][,input[[jth_ref("SIbpEnd", j)]]]+adjust    
      }
    } #end SI total bp calculation
  }#end calculateTotalBP

  create_pChart <- function(j) {
    #this function makes the chromsomeview chart  
    #subset whole chart based on selection
    chromChart <- values[[input[[jth_ref("datasets", j)]]]]
    chromChart <- chromChart[chromChart[,input[[jth_ref("chrColumn", j)]]]==input[[jth_ref("chr", j)]],]
    
    if(input[[jth_ref("plotAll", j)]]==FALSE){
      for(i in input[[jth_ref("traitColumns", j)]]){
        chromChart <- chromChart[chromChart[,i] %in% input[[jth_ref(i, j)]],]
      }    
      if(length(input[[jth_ref("traitColumns", j)]]) > 1){
        chromChart$trait <- do.call(paste,c(chromChart[,input[[jth_ref("traitColumns", j)]]],sep="_"))
      }else{
        chromChart$trait <- chromChart[,input[[jth_ref("traitColumns", j)]]]
      }
    }else{
      chromChart$trait <- input[[jth_ref("datasets", j)]]
    }
    
    #Separate Support Interval data from GWAS data, if support, GWAS data is assumed to be anything that has an NA in the SIbpStart column
    if(input[[jth_ref("supportInterval", j)]] == TRUE){
      SIchart <- chromChart[!(is.na(chromChart[,input[[jth_ref("SIbpStart", j)]]])),]
      chromChart <- chromChart[is.na(chromChart[,input[[jth_ref("SIbpStart", j)]]]),]
    }        
    #check if there is any data for the selected traits
    chromChart <- chromChart[!(is.na(chromChart[,input[[jth_ref("bpColumn", j)]]])),]
    chromChart <- chromChart[!(is.na(chromChart[,input[[jth_ref("yAxisColumn", j)]]])),]
    
    #if checked, filter for only overlapping SNPs
    if(!is.null(input[[jth_ref("overlaps", j)]]) & input[[jth_ref("overlaps", j)]] == TRUE){
      chromChart <- findGWASOverlaps(chromChart, j)
    }
    
    if(nrow(chromChart)==0){ #nothing is in the window, but lets still make a data.frame
      chromChart <- values[[input[[jth_ref("datasets", j)]]]][1,]
      chromChart[,input[[jth_ref("yAxisColumn", j)]]] <- -1    
      if(length(input[[jth_ref("traitColumns", j)]]) > 1){
        chromChart$trait <- do.call(paste,c(chromChart[,input[[jth_ref("traitColumns", j)]]],sep="_"))
      }else{
        chromChart$trait <- chromChart[,input[[jth_ref("traitColumns", j)]]]
      }             
    }    
    colorTable <- getColorTable(j)
    
    #take -log10 of y-axis column if requested
    if(input[[jth_ref("logP", j)]] == TRUE && chromChart[1,input[[jth_ref("yAxisColumn", j)]]] != -1){
      chromChart[,input[[jth_ref("yAxisColumn", j)]]] <- -log(chromChart[,input[[jth_ref("yAxisColumn", j)]]],10)
    }
    
    #check if there is too much data (>2500 data points), trim to 2500
    if(nrow(chromChart)>2500){
      cutVal <- sort(chromChart[,input[[jth_ref("yAxisColumn", j)]]],decreasing = T)[2500]
      chromChart <- chromChart[chromChart[,input[[jth_ref("yAxisColumn", j)]]] >= cutVal,]
    }
    
    #calculate window for plotband
    pbWin <- isolate({
      center <- as.numeric(input[[jth_ref("selected", j)]][[1]])
      winHigh <- center + input[[jth_ref("window", j)]][1]
      winLow <- center - input[[jth_ref("window", j)]][1]
      list(winLow=winLow,winHigh=winHigh)
    })
    
    pkTable <- data.frame(x=chromChart[,input[[jth_ref("bpColumn", j)]]],y=chromChart[,input[[jth_ref("yAxisColumn", j)]]],trait=chromChart$trait,
                          name=sprintf("Base Pair: %1$s<br/>Chromosome: %2$s<br/>",
                                       prettyNum(chromChart[,input[[jth_ref("bpColumn", j)]]], big.mark = ","),
                                       chromChart[,input[[jth_ref("chrColumn", j)]]]
                          ),
                          url="http://danforthcenter.org",
                          chr=chromChart[,input[[jth_ref("chrColumn", j)]]],
                          bp=chromChart[,input[[jth_ref("bpColumn", j)]]],stringsAsFactors=FALSE)
    pkSeries <- lapply(split(pkTable, pkTable$trait), function(x) {
      res <- lapply(split(x, rownames(x)), as.list)
      names(res) <- NULL
      res <- res[order(sapply(res, function(x) x$x))]
      return(res)
    })
    
    #build JL series
    if(input[[jth_ref("supportInterval", j)]]==TRUE){
      if(nrow(SIchart)==0){ #nothing is in the window, but lets still make a data.frame
        SIchart <- values[[input[[jth_ref("datasets", j)]]]][1,]
        SIchart[,input[[jth_ref("SIyAxisColumn", j)]]] <- -1    
        if(length(input[[jth_ref("traitColumns", j)]]) > 1){
          SIchart$trait <- do.call(paste,c(SIchart[,input[[jth_ref("traitColumns", j)]]],sep="_"))
        }else{
          SIchart$trait <- SIchart[,input[[jth_ref("traitColumns", j)]]]
        }
      }
      SIchart$loc_el <- SIchart$trait
      SIchart$trait <- paste(SIchart$trait,"Int",sep="_")
      SIchart <- SIchart[order(SIchart[[input[[jth_ref("SIbpStart", j)]]]]),]
      jlTable <- adply(SIchart,1,function(x) {data.frame(x=c(x[[input[[jth_ref("SIbpStart", j)]]]],x[[input[[jth_ref("SIbpEnd", j)]]]],x[[input[[jth_ref("SIbpEnd", j)]]]]),y=c(x[[input[[jth_ref("SIyAxisColumn", j)]]]],x[[input[[jth_ref("SIyAxisColumn", j)]]]],NA),trait=x$trait,
                                                         name=sprintf("<table cellpadding='4' style='line-height:1.5'><tr><td align='left'>Interval: %1$s-%2$s<br>Chromosome: %3$s</td></tr></table>",
                                                                      prettyNum(x[[input[[jth_ref("SIbpStart", j)]]]], big.mark = ","),
                                                                      prettyNum(x[[input[[jth_ref("SIbpEnd", j)]]]], big.mark = ","),
                                                                      x[[input[[jth_ref("chrColumn", j)]]]]
                                                         ),loc_el=x$loc_el,bp=x[[input[[jth_ref("bpColumn", j)]]]],chr=x[[input[[jth_ref("chrColumn", j)]]]],stringsAsFactors=FALSE
      )}#end jlTable and function
      )#end adply
      jlTable <- jlTable[,c("x","y","trait","name","loc_el","bp","chr")]
    }#end build jlTable if support intervals        
    
    a <- rCharts::Highcharts$new()
    a$LIB$url <- 'highcharts/' #use the local copy of highcharts, not the one installed by rCharts
    a$xAxis(title = list(text = "Base Pairs"),startOnTick=TRUE,min=1,max=chrSize[values[[jth_ref("organism", j)]]][[1]][as.numeric(input[[jth_ref("chr", j)]])],endOnTick=FALSE,
            plotBands = list(list(from=pbWin$winLow,to=pbWin$winHigh,color='rgba(68, 170, 213, 0.4)')))
    
    if(input[[jth_ref("axisLimBool", j)]] == TRUE){
      a$yAxis(title=list(text=input[[jth_ref("yAxisColumn", j)]]),min=input[[jth_ref("axisMin", j)]],max=input[[jth_ref("axisMax", j)]],startOnTick=FALSE)
    }else{
      a$yAxis(title=list(text=input[[jth_ref("yAxisColumn", j)]]),startOnTick=FALSE)      
    }    
    
    if(input[[jth_ref("supportInterval", j)]]==TRUE){
      if(input[[jth_ref("SIaxisLimBool", j)]] == TRUE){
        a$yAxis(title=list(text=input[[jth_ref("SIyAxisColumn", j)]]),min=input[[jth_ref("SIaxisMin", j)]],max=input[[jth_ref("SIaxisMax", j)]],gridLineWidth=0,minorGridLineWidth=0,startOnTick=FALSE,opposite=TRUE,replace=FALSE)
      }else{
        a$yAxis(title=list(text=input[[jth_ref("SIyAxisColumn", j)]]),gridLineWidth=0,minorGridLineWidth=0,startOnTick=FALSE,opposite=TRUE,replace=FALSE)   
      }
      
      if(SIchart[1,input[[jth_ref("SIyAxisColumn", j)]]] != -1){
        d_ply(jlTable,.(trait),function(x){
          a$series(
            data = toJSONArray2(x,json=F,names=T),
            type = "line",
            name = unique(x$trait),
            yAxis=1,
            color = colorTable$color[colorTable$trait == as.character(unique(x$loc_el))])})            
      }
    }
    
    if(chromChart[1,input[[jth_ref("yAxisColumn", j)]]] != -1){
      invisible(sapply(pkSeries, function(x) {if(length(x)==0){return()};a$series(data = x, type = "scatter", turboThreshold=5000, name = paste0(x[[1]]$trait), color = colorTable$color[colorTable$trait == as.character(x[[1]]$trait)])}))
    }
    a$chart(zoomType="x", alignTicks=FALSE,events=list(click = "#!function(event) {this.tooltip.hide();}!#"))
    a$title(text=paste(input[[jth_ref("datasets", j)]],"Results for Chromosome",input[[jth_ref("chr", j)]],sep=" "))
    a$subtitle(text="Rollover for more info. Drag chart area to zoom. Click point for zoomed annotated plot.")
    
    js <- jth_ref("", j)
    doClickOnPoint <- sprintf("#! function(){$('input#selected%s').val(this.options.bp); $('input#selected%s').trigger('change');} !#", js, js)
    doClickOnLine <- doClickOnPoint # they happen to be identical, but we could make them different
    a$plotOptions(
      scatter = list(
        cursor = "pointer",
        point = list(
          events = list(
            click = doClickOnPoint
          )
        ),
        marker = list(
          symbol = "circle",
          radius = 5
        ),
        tooltip = list(
          headerFormat = "<b>{series.name}</b><br/>{point.key}<br/>Y-value: {point.y}<br/>",
          pointFormat = "",
          followPointer = TRUE
        )
      ),
      line = list(
        lineWidth = 10,
        dashStyle = 'Solid',
        cursor = "pointer",
        point = list(
          events = list(
            click = doClickOnLine
          )
        ),
        marker = list(
          enabled = FALSE,
          states = list(hover = list(enabled=FALSE))
        )
      ),
      spline = list(
        lineWidth = 3,
        cursor = "pointer"
      )
    )
    a$exporting(enabled=TRUE,filename='chromChart',sourceWidth=2000)
    a$set(dom = jth_ref('pChart', j))
    return(a)
  }
  output$pChart <- renderChart(create_pChart(1))
  output$pChart2 <- renderChart(create_pChart(2))
  
  #Genome wide chart
  create_gChart <- function(j) {
    calculateTotalBP(j)
    
    #subset whole chart based on selection
    genomeChart <- values[[input[[jth_ref("datasets", j)]]]]
    if(input[[jth_ref("plotAll", j)]] == FALSE){
      for(i in input[[jth_ref("traitColumns", j)]]){
        genomeChart <- genomeChart[genomeChart[,i] %in% input[[jth_ref(i, j)]],]
      }    
      if(length(input[[jth_ref("traitColumns", j)]]) > 1){
        genomeChart$trait <- do.call(paste,c(genomeChart[,input[[jth_ref("traitColumns", j)]]],sep="_"))    
      }else{
        genomeChart$trait <- genomeChart[,input[[jth_ref("traitColumns", j)]]]
      }
    }else{
      genomeChart$trait <- input[[jth_ref("datasets", j)]]
    }

    #Separate Support Interval data from GWAS data, if support, GWAS data is assumed to be anything that has an NA in the SIbpStart column
    if(input[[jth_ref("supportInterval", j)]] == TRUE){
      SIchart <- genomeChart[!(is.na(genomeChart[,input[[jth_ref("SIbpStart", j)]]])),]
      genomeChart <- genomeChart[is.na(genomeChart[,input[[jth_ref("SIbpStart", j)]]]),]
    }
    
    #filter genomeChart for only rows that have a base pair and yaxis value
    genomeChart <- genomeChart[!(is.na(genomeChart[,input[[jth_ref("bpColumn", j)]]])),]
    genomeChart <- genomeChart[!(is.na(genomeChart[,input[[jth_ref("yAxisColumn", j)]]])),]
    
    #if checked, filter for only overlapping SNPs
    if(!is.null(input[[jth_ref("overlaps", j)]]) & input[[jth_ref("overlaps", j)]] == TRUE){
      genomeChart <- findGWASOverlaps(genomeChart, j)
    }
    
    #check if there is any data for the selected traits
    if(nrow(genomeChart)==0){ #nothing is in the window, but lets still make a data.frame
      genomeChart <- values[[input[[jth_ref("datasets", j)]]]][1,]
      genomeChart[,input[[jth_ref("yAxisColumn", j)]]] <- -1    
      if(length(input[[jth_ref("traitColumns", j)]]) > 1){
        genomeChart$trait <- do.call(paste,c(genomeChart[,input[[jth_ref("traitColumns", j)]]],sep="_"))    
      }else{
        genomeChart$trait <- genomeChart[,input[[jth_ref("traitColumns", j)]]]
      }             
    }
    
    #take -log10 of y-axis column if requested
    if(input[[jth_ref("logP", j)]] == TRUE && genomeChart[1,input[[jth_ref("yAxisColumn", j)]]] != -1){
      genomeChart[,input[[jth_ref("yAxisColumn", j)]]] <- -log(genomeChart[,input[[jth_ref("yAxisColumn", j)]]],10)
    }
    
    #check if there is too much data (>2500 data points), trim to 2500
    if(nrow(genomeChart)>2500){
      cutVal <- sort(genomeChart[,input[[jth_ref("yAxisColumn", j)]]],decreasing = T)[2500]
      genomeChart <- genomeChart[genomeChart[,input[[jth_ref("yAxisColumn", j)]]] >= cutVal,]
    }
    
    colorTable <- getColorTable(j)
    genomeTable <- data.frame(x=genomeChart$totalBP,y=genomeChart[,input[[jth_ref("yAxisColumn", j)]]],trait=genomeChart$trait,
                              #                               name=sprintf("<table cellpadding='4' style='line-height:1.5'><tr><th>%1$s</th></tr><tr><td align='left'>RMIP: %2$s<br>Location: %3$s<br>Base Pairs: %4$s<br>SNP: %5$s<br>Chromosome: %6$s</td></tr></table>",
                              #                               name=sprintf("<table cellpadding='4' style='line-height:1.5'><tr><th>%1$s</th></tr><tr><td align='left'>Y-value: %2$s<br>Base Pairs: %3$s<br>Chromosome: %4$s</td></tr></table>",
                              name=sprintf("Base Pair: %1$s<br/>Chromosome: %2$s<br/>",
                                           #                                            genomeChart$trait,
                                           #                                            genomeChart[,input[[jth_ref("yAxisColumn", j)]]],
                                           #                                            genomeChart$loc,
                                           prettyNum(genomeChart[,input[[jth_ref("bpColumn", j)]]], big.mark = ","),
                                           #                                            genomeChart$SNP,
                                           genomeChart[,input[[jth_ref("chrColumn", j)]]]
                              ),
                              url="http://danforthcenter.org",
                              chr=genomeChart[,input[[jth_ref("chrColumn", j)]]],
                              bp=genomeChart[,input[[jth_ref("bpColumn", j)]]],stringsAsFactors=FALSE)
    genomeSeries <- lapply(split(genomeTable, genomeTable$trait), function(x) {
      res <- lapply(split(x, rownames(x)), as.list)
      names(res) <- NULL
      res <- res[order(sapply(res, function(x) x$x))]
      return(res)
    })
    #     
    #build JL series
    if(input[[jth_ref("supportInterval", j)]]==TRUE){
      if(nrow(SIchart)==0){ #nothing is in the window, but lets still make a data.frame
        SIchart <- values[[input[[jth_ref("datasets", j)]]]][1,]
        SIchart[,input[[jth_ref("SIyAxisColumn", j)]]] <- -1    
        if(length(input[[jth_ref("traitColumns", j)]]) > 1){
          SIchart$trait <- do.call(paste,c(SIchart[,input[[jth_ref("traitColumns", j)]]],sep="_"))    
        }else{
          SIchart$trait <- SIchart[,input[[jth_ref("traitColumns", j)]]]
        }             
      }
      SIchart$loc_el <- SIchart$trait
      SIchart$trait <- paste(SIchart$trait,"Int",sep="_")
      SIchart <- SIchart[order(SIchart$SIbpStartTotal),]
      jlTable <- adply(SIchart,1,function(x) {data.frame(x=c(x$SIbpStartTotal,x$SIbpEndTotal,x$SIbpEndTotal),y=c(x[[input[[jth_ref("SIyAxisColumn", j)]]]],x[[input[[jth_ref("SIyAxisColumn", j)]]]],NA),trait=x$trait,
                                                         #name=sprintf("<table cellpadding='4' style='line-height:1.5'><tr><th>%1$s</th></tr><tr><td align='left'>Y-value: %2$.2f <br>Interval: %3$s-%4$s<br>Chromosome: %5$s</td></tr></table>",
                                                         name=sprintf("<table cellpadding='4' style='line-height:1.5'><tr><td align='left'>Interval: %1$s-%2$s<br>Chromosome: %3$s</td></tr></table>",
                                                                      #                                                                 x$trait,
                                                                      #                                                                 x[[input[[jth_ref("SIyAxisColumn", j)]]]],
                                                                      prettyNum(x[[input[[jth_ref("SIbpStart", j)]]]], big.mark = ","),
                                                                      prettyNum(x[[input[[jth_ref("SIbpEnd", j)]]]], big.mark = ","),
                                                                      x[[input[[jth_ref("chrColumn", j)]]]]
                                                         ),loc_el=x$loc_el,bp=x[[input[[jth_ref("bpColumn", j)]]]],chr=x[[input[[jth_ref("chrColumn", j)]]]],stringsAsFactors=FALSE
                                                         #                                                   
                                                         #                                                   totalBP=x$totalBP,
                                                         #                                                   chr=x$Chromosome,stringsAsFactors=FALSE
      )}#end jlTable and function
      )#end adply
      jlTable <- jlTable[,c("x","y","trait","name","loc_el","bp","chr")]
      #jlTable <- jlTable[order(jlTable$x),]
    }#end build jlTable if support intervals
    
    #build list for where to put plotbands for this organism
    bigList <- list()
    cumBP<-c(0,cumsum(as.numeric(chrSize[values[[jth_ref("organism", j)]]][[1]])))
    for(i in 1:(length(cumBP)-1)){
      if(i %% 2 == 0 ){ #even
        bigList[[length(bigList)+1]] <- list(from=cumBP[i]+1,to=cumBP[i+1],label=list(text=chrName[values[[jth_ref("organism", j)]]][[1]][i],style=list(color="#6D869F"),verticalAlign="bottom"))
      }else{ #odd
        bigList[[length(bigList)+1]] <- list(from=cumBP[i]+1,to=cumBP[i+1],color='rgba(68, 170, 213, 0.1)',label=list(text=chrName[values[[jth_ref("organism", j)]]][[1]][i],style=list(color="#6D869F"),verticalAlign="bottom"))
      }
    }    
    
    c <- rCharts::Highcharts$new()
    c$LIB$url <- 'highcharts/'
    c$xAxis(title = list(text = "Chromosome",margin=15),startOnTick=TRUE,min=0,max=sum(as.numeric(chrSize[values[[jth_ref("organism", j)]]][[1]])),endOnTick=FALSE,labels=list(enabled=FALSE),tickWidth=0,
            plotBands = bigList)   
    
    if(input[[jth_ref("axisLimBool", j)]] == TRUE){       
      c$yAxis(title=list(text=input[[jth_ref("yAxisColumn", j)]]),min=input[[jth_ref("axisMin", j)]],max=input[[jth_ref("axisMax", j)]],startOnTick=FALSE)
    }else{
      c$yAxis(title=list(text=input[[jth_ref("yAxisColumn", j)]]),startOnTick=FALSE)      
    }
    
    if(input[[jth_ref("supportInterval", j)]]==TRUE){
      if(input[[jth_ref("SIaxisLimBool", j)]] == TRUE){
        c$yAxis(title=list(text=input[[jth_ref("SIyAxisColumn", j)]]),min=input[[jth_ref("SIaxisMin", j)]],max=input[[jth_ref("SIaxisMax", j)]],gridLineWidth=0,minorGridLineWidth=0,startOnTick=FALSE,opposite=TRUE,replace=FALSE)   
      }else{
        c$yAxis(title=list(text=input[[jth_ref("SIyAxisColumn", j)]]),gridLineWidth=0,minorGridLineWidth=0,startOnTick=FALSE,opposite=TRUE,replace=FALSE)
      }
      
      if(SIchart[1,input[[jth_ref("SIyAxisColumn", j)]]] != -1){
        d_ply(jlTable,.(trait),function(x){
          c$series(
            data = toJSONArray2(x,json=F,names=T),
            type = "line",
            name = unique(x$trait),
            dashStyle = 'Solid',
            marker = list(enabled=F),
            yAxis=1,           
            color = colorTable$color[colorTable$trait == as.character(unique(x$loc_el))])})            
      }
    }
    if(genomeChart[1,input[[jth_ref("yAxisColumn", j)]]] != -1){
      invisible(sapply(genomeSeries, function(x) {if(length(x)==0){return()};c$series(data = x, turboThreshold=5000,type = "scatter", color = colorTable$color[colorTable$trait == as.character(x[[1]]$trait)], name = paste0(x[[1]]$trait))}))
    }
    
    c$chart(zoomType="x",alignTicks=FALSE,events=list(click = "#!function(event) {this.tooltip.hide();}!#"))
    c$title(text=paste(input[[jth_ref("datasets", j)]]," Results",sep=" "))
    c$subtitle(text="Rollover for more info. Drag chart area to zoom. Click point to switch to chromosome and annotation view.")
    
    js <- jth_ref("", j)
    doClickOnPoint <- sprintf(
      paste(
        "#! function(){$('select#chr%s').val(this.options.chr); $('select#chr%s').trigger('change'); $('input#selected%s').val(this.options.bp);",
        "$('input#selected%s').trigger('change'); $('ul#datatabs li').eq(0).removeClass('active');",
        "$('ul#datatabs li').eq(1).removeClass('active'); $('ul#datatabs li').eq(2).removeClass('active');",
        "$('ul#datatabs li').eq(4).removeClass('active');",
        "$('ul#datatabs li').eq(3).addClass('active');",
        "$('#pChart%s').trigger('change');$('#pChart%s').trigger('shown');",
        "$('.tab-content div').toggleClass(function(){if(this.getAttribute('data-value')=='Chrom' || this.getAttribute('data-value')=='WhGen'){return 'active';}else{return '';}});",
        "$('.tab-content div').trigger('change');$('ul#datatabs li').trigger('change');} !#"
      ), js, js, js, js, js, js)
    doClickOnLine <- doClickOnPoint # they happen to be identical, but we could make them different
    c$plotOptions(
      scatter = list(
        cursor = "pointer",
        point = list(
          events = list(
            #click = "#! function() { window.open(this.options.url); } !#")), #open webpage
            #click = "#! function(event) {alert(this.name);} !#")), #display popup
            #click = "#! function(event) {console.log(this);} !#")), #write object to log
            #click = "#! function(){$('input#selected').val(134); $('input#selected').trigger('change');} !#")),
            click = doClickOnPoint
          )
        ),
        marker = list(
          symbol = "circle",
          radius = 5
        ),
        tooltip = list(
          headerFormat = "<b>{series.name}</b><br/>{point.key}<br/>Y-value: {point.y}<br/>",
          pointFormat = "",
          followPointer = TRUE
        )
          ),
      line = list(
        lineWidth = 10,
        cursor = "pointer",
        point = list(
          events = list(
            #click = "#! function() { window.open(this.url); } !#")), #open webpage
            #click = "#! function(event) {alert(this.url);} !#")), #display popup
            #click = "#! function(event) {console.log(this);} !#")), #write object to log
            #click = "#! function(){$('input#selected').val(134); $('input#selected').trigger('change');} !#")),
            #click = "#! function(){$('input#selected').val(this.options.bp); $('input#selected').trigger('change');} !#")),
            #click = "#! function(){$('select#chr').val(this.options.chr); $('select#chr').trigger('change'); $('input#selected').val(this.options.bp); $('input#selected').trigger('change');  $('ul#datatabs li').eq(2).removeClass('active'); $('ul#datatabs li').eq(3).addClass('active'); $('.tab-content div').toggleClass('active'); $('#pChart').trigger('shown')} !#")),
            click = doClickOnLine
          )
        ),
        marker = list(
          enabled = FALSE,
          states = list(hover = list(enabled=FALSE)),
          symbol = "square"
        )
          )            
        )#end plotOptions        
    #c$tooltip(useHTML = T, formatter = "#! function() { return this.point.name; } !#")
    #c$tooltip(formatter = "#! function() { return this.point.name; } !#")
    c$exporting(enabled=TRUE,filename='genomeChart',sourceWidth=2000)
    if(!is.null(input[[jth_ref("legend", j)]]) & input[[jth_ref("legend", j)]] == TRUE){
      c$legend(enabled=FALSE)
    }
    
    c$credits(enabled=TRUE)
    c$set(dom = jth_ref('gChart', j))
    return(c)
  }
  output$gChart <- renderChart(create_gChart(1))
  output$gChart2 <- renderChart(create_gChart(2))

  create_zChart <- function(j) {
    if (is.null(input[[jth_ref("selected", j)]])) return()
    
    centerBP <- as.numeric(input[[jth_ref("selected", j)]][[1]])
    winHigh <- centerBP + input[[jth_ref("window", j)]][1]
    winLow <- centerBP - input[[jth_ref("window", j)]][1]
    if (winLow < 0) {winLow <- 0}
    
    zoomChart <- values[[input[[jth_ref("datasets", j)]]]]
    zoomChart <- zoomChart[zoomChart[,input[[jth_ref("chrColumn", j)]]]==input[[jth_ref("chr", j)]],]    
    
    if(input[[jth_ref("plotAll", j)]] == FALSE){
      for(i in input[[jth_ref("traitColumns", j)]]){
        zoomChart <- zoomChart[zoomChart[,i] %in% input[[jth_ref(i, j)]],]
      }
      
      if(length(input[[jth_ref("traitColumns", j)]]) > 1){
        zoomChart$trait <- do.call(paste,c(zoomChart[,input[[jth_ref("traitColumns", j)]]],sep="_"))    
      }else{
        zoomChart$trait <- zoomChart[,input[[jth_ref("traitColumns", j)]]]
      }
    }else{
      zoomChart$trait <- input[[jth_ref("datasets", j)]]
    }
    
    #Separate Support Interval data from GWAS data, if support, GWAS data is assumed to be anything that has an NA in the SIbpStart column
    if(input[[jth_ref("supportInterval", j)]] == TRUE){
      SIchart <- zoomChart[!(is.na(zoomChart[,input[[jth_ref("SIbpStart", j)]]])),]
      zoomChart <- zoomChart[is.na(zoomChart[,input[[jth_ref("SIbpStart", j)]]]),]
      #not sure the below logic works for subsetting SIchart, probably not necessary anyways, since there are usually very few SI rows for one chromosome anyways (e.g. small overhead)
      #SIchart <- SIchart[((SIchart[,input[[jth_ref("SIbpStart", j)]]] <= winHigh & SIchart[,input[[jth_ref("SIbpStart", j)]]] >= winLow) | (SIchart[,input[[jth_ref("SIbpEnd", j)]]] <= winHigh & SIchart[,input[[jth_ref("SIbpEnd", j)]]] >= winLow)),]
    }    
    
    zoomChart <- zoomChart[(zoomChart[,input[[jth_ref("bpColumn", j)]]] <= winHigh) & (zoomChart[,input[[jth_ref("bpColumn", j)]]] >= winLow),]    
    
    #filter for only rows that have a base pair value
    zoomChart <- zoomChart[!(is.na(zoomChart[,input[[jth_ref("bpColumn", j)]]])),]
    zoomChart <- zoomChart[!(is.na(zoomChart[,input[[jth_ref("yAxisColumn", j)]]])),]
    
    #if checked, filter for only overlapping SNPs
    if(!is.null(input[[jth_ref("overlaps", j)]]) & input[[jth_ref("overlaps", j)]] == TRUE){
      zoomChart <- findGWASOverlaps(zoomChart, j)
    }                    
    
    if(nrow(zoomChart)==0){ #nothing is in the window, but lets still make a data.frame
      zoomChart <- values[[input[[jth_ref("datasets", j)]]]][1,]
      zoomChart[,input[[jth_ref("yAxisColumn", j)]]] <- -1    
      if(length(input[[jth_ref("traitColumns", j)]]) > 1){
        zoomChart$trait <- do.call(paste,c(zoomChart[,input[[jth_ref("traitColumns", j)]]],sep="_"))
      }else{
        zoomChart$trait <- zoomChart[,input[[jth_ref("traitColumns", j)]]]
      }                   
    }
    colorTable <- getColorTable(j)
    
    #take -log10 of y-axis column if requested
    if(input[[jth_ref("logP", j)]] == TRUE && zoomChart[1,input[[jth_ref("yAxisColumn", j)]]] != -1){
      zoomChart[,input[[jth_ref("yAxisColumn", j)]]] <- -log(zoomChart[,input[[jth_ref("yAxisColumn", j)]]],10)
    }                
    
    #check if there is too much data (>2500 data points), trim to 2500
    if(nrow(zoomChart)>2500){
      cutVal <- sort(zoomChart[,input[[jth_ref("yAxisColumn", j)]]],decreasing = T)[2500]
      zoomChart <- zoomChart[zoomChart[,input[[jth_ref("yAxisColumn", j)]]] >= cutVal,]
    }                
    
    zoomTable <- data.frame(x=zoomChart[,input[[jth_ref("bpColumn", j)]]],y=zoomChart[,input[[jth_ref("yAxisColumn", j)]]],trait=zoomChart$trait,
#                                         name=sprintf("<table cellpadding='4' style='line-height:1.5'><tr><th>%1$s</th></tr><tr><td align='left'>RMIP: %2$s<br>Location: %3$s<br>Base Pairs: %4$s<br>SNP: %5$s<br>Chromosome: %6$s</td></tr></table>",
                            name=sprintf("<table cellpadding='4' style='line-height:1.5'><tr><th>%1$s</th></tr><tr><td align='left'>Y-axis value: %2$s<br>Base Pairs: %3$s<br>Chromosome: %4$s</td></tr></table>",                                         
                                         zoomChart$trait,
                                         zoomChart[,input[[jth_ref("yAxisColumn", j)]]],
                                         #zoomChart$loc,
                                         prettyNum(zoomChart[,input[[jth_ref("bpColumn", j)]]], big.mark = ","),
                                         #zoomChart$SNP,
                                         zoomChart[,input[[jth_ref("chrColumn", j)]]]
                            ),
                            url="http://danforthcenter.org",
                            chr=zoomChart[,input[[jth_ref("chrColumn", j)]]],
                            bp=zoomChart[,input[[jth_ref("bpColumn", j)]]])
    zoomSeries <- lapply(split(zoomTable, zoomTable$trait), function(x) {
      res <- lapply(split(x, rownames(x)), as.list)
      names(res) <- NULL
      res <- res[order(sapply(res, function(x) x$x))]
      return(res)
    })
    
#     #build the sliding window GWAS summary line
#     if(input$SlidingWinCheckbox==TRUE){
#       winTable <- fullWinTable[fullWinTable$chr==input$chr & fullWinTable$el %in% els & fullWinTable$loc %in% input$locs]
#     }else{
#       winTable <- data.frame()
#     }      
#     
#     if(nrow(winTable)==0){ #make a dummy point if there are none in this plot
#       winTable <- fullWinTable[1,]
#       winTable$chr <- input$chr
#       winTable$env <- input$locs[1]
#       winTable$el <- els[1]
#       winTable$sumRMIP <- -5
#     }
#     
#     winTable$loc_el <- paste(winTable$loc,winTable$el,"Window",sep="_")
#     winTableSeries <- adply(winTable,1,function(x) {data.table(x=x$bp,y=x$sumRMIP,element=x$loc_el,
#                                                                name=sprintf("<table cellpadding='4' style='line-height:1.5'><tr><th>%1$s</th></tr><tr><td align='left'>sumRMIP: %2$s<br>Num Points: %3$s<br>Loc: %4$s<br>Base Pair: %5$s<br>cM: %6$s<br>Chromosome: %7$s</td></tr></table>",
#                                                                             x$el,
#                                                                             x$sumRMIP,
#                                                                             x$numPoints,
#                                                                             x$loc,
#                                                                             prettyNum(x$bp, big.mark = ","),
#                                                                             round(x$cM,digits=3),
#                                                                             x$chr
#                                                                ),
#                                                                bp=x$bp,
#                                                                chr=x$chr                                                    
#     )})
#     #this removes columns not needed (using data.table which is why its complicated)
#     winTableSeries[,colnames(winTableSeries)[which(!colnames(winTableSeries) %in% c("x","y","name","element","bp","chr"))] := NULL,with=FALSE]
#     #winTableSeries <- winTableSeries[order(winTableSeries$bp,winTableSeries$x),] #should be ordered from the file
#     #build JL series
    if(input[[jth_ref("supportInterval", j)]]==TRUE){
      if(nrow(SIchart) == 0){ #make a dummy table, but we won't plot the series anyways
        SIchart <- values[[input[[jth_ref("datasets", j)]]]][1,]
        SIchart[,input[[jth_ref("SIyAxisColumn", j)]]] <- -1    
        if(length(input[[jth_ref("traitColumns", j)]]) > 1){
          SIchart$trait <- do.call(paste,c(SIchart[,input[[jth_ref("traitColumns", j)]]],sep="_"))
        }else{
          SIchart$trait <- SIchart[,input[[jth_ref("traitColumns", j)]]]
        }                      
      }     
      SIchart$loc_el <- SIchart$trait
      SIchart$trait <- paste(SIchart$trait,"Int",sep="_")
      SIchart <- SIchart[order(SIchart[[input[[jth_ref("SIbpStart", j)]]]]),]
      jlTable <- adply(SIchart,1,function(x) {data.frame(x=c(x[[input[[jth_ref("SIbpStart", j)]]]],x[[input[[jth_ref("SIbpEnd", j)]]]],x[[input[[jth_ref("SIbpEnd", j)]]]]),y=c(x[[input[[jth_ref("SIyAxisColumn", j)]]]],x[[input[[jth_ref("SIyAxisColumn", j)]]]],NA),trait=x$trait,
                                                         name=sprintf("<table cellpadding='4' style='line-height:1.5'><tr><td align='left'>Interval: %1$s-%2$s<br>Chromosome: %3$s</td></tr></table>",
#                                                                      x$trait,
#                                                                      x[[input[[jth_ref("SIyAxisColumn", j)]]]],
                                                                      prettyNum(x[[input[[jth_ref("SIbpStart", j)]]]], big.mark = ","),
                                                                      prettyNum(x[[input[[jth_ref("SIbpEnd", j)]]]], big.mark = ","),
                                                                      x[[input[[jth_ref("chrColumn", j)]]]]
                                                         ),loc_el=x$loc_el,bp=x[[input[[jth_ref("bpColumn", j)]]]],chr=x[[input[[jth_ref("chrColumn", j)]]]],stringsAsFactors=FALSE
#                                                   
#                                                   totalBP=x$totalBP,
#                                                   chr=x$Chromosome,stringsAsFactors=FALSE
      )}#end jlTable and function
      )#end adply
      jlTable <- jlTable[,c("x","y","trait","name","loc_el","bp","chr")]
      #jlTable <- jlTable[order(jlTable$x),]
    }#end if support interval
#     
#     
#     jl$loc_el <- paste(jl$env,jl$el,"JL",sep="_")      
#     jlTable <- adply(jl,1,function(x) {data.frame(x=c(x$lowerCIbp,x$upperCIbp,x$upperCIbp),y=c(x$F,x$F,NA),element=x$loc_el,url="http://danforthcenter.org",
#                                                   name=sprintf("<table cellpadding='4' style='line-height:1.5'><tr><th>%1$s</th></tr><tr><td align='left'>F: %2$.2f  -logP: %3$.2f<br>Location: %4$s<br>Base Pair: %5$s<br>SI: %6$s-%7$s<br>SNP: %8$s<br>Chromosome: %9$s</td></tr></table>",       
#                                                                x$el,
#                                                                x$F,
#                                                                x$negLogP,
#                                                                x$env,
#                                                                prettyNum(x$bp, big.mark = ","),
#                                                                prettyNum(x$lowerCIbp, big.mark = ","),
#                                                                prettyNum(x$upperCIbp, big.mark = ","),
#                                                                x$Name,
#                                                                x$Chromosome
#                                                   ),
#                                                   bp=x$bp,
#                                                   chr=x$Chromosome,stringsAsFactors=FALSE
#     )}
#     )
#     jlTable <- jlTable[,c("x","y","name","element","bp","chr","url")]
#     
#     jlTable <- jlTable[order(jlTable$bp,jlTable$x),]
#     #     jlSeries <- lapply(split(jlTable, jlTable$element), function(x) {
#     #       res <- lapply(split(x, rownames(x)), as.list)
#     #       names(res) <- NULL
#     #       #res <- res[order(sapply(res, function(x) x$x))] #
#     #       return(res)
#     #     })
    
    #build annotation series
    #thisChrAnnot <- subset(annotGeneLoc,chromosome==input[[jth_ref("chr", j)]])
    thisChrAnnot <- subset(annotGeneLoc[values[[jth_ref("organism", j)]]][[1]],chromosome==input[[jth_ref("chr", j)]])
    thisAnnot <- thisChrAnnot[thisChrAnnot$transcript_start >= winLow & thisChrAnnot$transcript_end <= winHigh,]
    if(nrow(thisAnnot)==0){ #nothing is in the window, but lets still make a data.frame (actually make it big just to hopefully pick up one row from each strand...)
      thisAnnot <- thisChrAnnot[1:100,]
    }
    thisAnnot <- thisAnnot[order(thisAnnot$transcript_start),]
    
#    urlBase <- 'http://www.maizesequence.org/Zea_mays/Transcript/ProteinSummary?db=core;t='
    urlBase <- 'http://maizegdb.org/cgi-bin/displaygenemodelrecord.cgi?id='
    soyurlBase <- 'http://www.soybase.org/sbt/search/search_results.php?category=FeatureName&search_term='
    araburlBase <- 'http://arabidopsis.org/servlets/TairObject?type=locus&name='
    sorgurlBase <- 'http://phytozome.jgi.doe.gov/pz/portal.html#!gene?search=1&detail=1&searchText=transcriptid:'
    legumeInfo_urlBase <- 'https://legumeinfo.org/gene_links/'
    
    annotYvalReverse <- 0.02
    #if(input[[jth_ref("axisLimBool", j)]] == TRUE){annotYvalReverse <- input[[jth_ref("axisMin", j)]] + 0.01}
    annotYvalForward <- annotYvalReverse + 0.04
    if (values[[jth_ref("organism", j)]] == "Corn") {
      annotTable <- adply(thisAnnot[thisAnnot$transcript_strand==1,],1,function(x) {
        data.frame(x=c(x$transcript_start,x$transcript_end,x$transcript_end),y=c(annotYvalForward,annotYvalForward,NA),url=paste0(urlBase,x$transcript_id),
          name=sprintf("<table cellpadding='4' style='line-height:1.5'><tr><th>%1$s</th></tr><tr><td align='left'>Location: %2$s-%3$s<br>Chromosome: %4$s, Strand: %5$s<br>%6$s</td></tr></table>",
            x$translation_id,
            prettyNum(x$transcript_start, big.mark = ","),
            prettyNum(x$transcript_end, big.mark = ","),
            x$chromosome,
            x$transcript_strand,
            x$V2
          ),
          marker=c(NA,"Arrow",NA),
          stringsAsFactors=FALSE
        )
      })
      annotTableReverse <- adply(thisAnnot[thisAnnot$transcript_strand==-1,],1,function(x) {
        data.frame(x=c(x$transcript_start,x$transcript_end,x$transcript_end),y=c(annotYvalReverse,annotYvalReverse,NA),url=paste0(urlBase,x$transcript_id),
          name=sprintf("<table cellpadding='4' style='line-height:1.5'><tr><th>%1$s</th></tr><tr><td align='left'>Location: %2$s-%3$s<br>Chromosome: %4$s, Strand: %5$s<br>%6$s</td></tr></table>",
            x$translation_id,
            prettyNum(x$transcript_start, big.mark = ","),
            prettyNum(x$transcript_end, big.mark = ","),
            x$chromosome,
            x$transcript_strand,
            x$V2
          ),
          marker=c("Arrow",NA,NA),
          stringsAsFactors=FALSE
        )
      })

    } else if(values[[jth_ref("organism", j)]] == "Soybean") { #strand is '+' or '-'
      annotTable <- adply(thisAnnot[thisAnnot$strand=="+",],1,function(x) {
        data.frame(x=c(x$transcript_start,x$transcript_end,x$transcript_end),y=c(annotYvalForward,annotYvalForward,NA),url=paste0(soyurlBase,x$transcript_id),
          name=sprintf("<table cellpadding='4' style='line-height:1.5'><tr><th>%1$s</th></tr><tr><td align='left'>Location: %2$s-%3$s, Protein Length: %4$s<br>Chromosome: %5$s, Strand: %6$s<br>Top TAIR Hit Desc.: %7$s<br>Top Uniref Hit Desc.: %8$s</td></tr></table>",
            x$transcript_id,
            prettyNum(x$transcript_start, big.mark = ","),
            prettyNum(x$transcript_end, big.mark = ","),
            x$Protein.Length,
            x$chromosome,                                                                           
            x$strand,
            x$TopTAIRHitDescription,
            x$TopUniref100DescriptionExtraSmall
          ),
          stringsAsFactors=FALSE
        )
      })
      annotTableReverse <- adply(thisAnnot[thisAnnot$strand=="-",],1,function(x) {
        data.frame(x=c(x$transcript_start,x$transcript_end,x$transcript_end),y=c(annotYvalReverse,annotYvalReverse,NA),url=paste0(soyurlBase,x$transcript_id),
          name=sprintf("<table cellpadding='4' style='line-height:1.5'><tr><th>%1$s</th></tr><tr><td align='left'>Location: %2$s-%3$s, Protein Length: %4$s<br>Chromosome: %5$s, Strand: %6$s<br>Top TAIR Hit Desc.: %7$s<br>Top Uniref Hit Desc.: %8$s</td></tr></table>",
            x$transcript_id,
            prettyNum(x$transcript_start, big.mark = ","),
            prettyNum(x$transcript_end, big.mark = ","),
            x$Protein.Length,
            x$chromosome,                                                                           
            x$strand,
            x$TopTAIRHitDescription,
            x$TopUniref100DescriptionExtraSmall
          ),
          stringsAsFactors=FALSE
        )
      })

    } else if(values[[jth_ref("organism", j)]] %in% c("Arabidopsis", "Arabidopsis thaliana")) { #strand is '+' or '-'
      annotTable <- adply(thisAnnot[thisAnnot$strand=="+",],1,function(x) {
        data.frame(x=c(x$transcript_start,x$transcript_end,x$transcript_end),y=c(annotYvalForward,annotYvalForward,NA),url=paste0(araburlBase,x$Locus),
          name=sprintf("<table cellpadding='4' style='line-height:1.5'><tr><th>%1$s</th></tr><tr><td align='left'>Location: %2$s-%3$s<br>Chromosome: %4$s, Strand: %5$s<br>Short Desc.: %6$s</td></tr></table>",
            x$name,
            prettyNum(x$transcript_start, big.mark = ","),
            prettyNum(x$transcript_end, big.mark = ","),
            x$chromosome,                                                                           
            x$strand,
            x$short_description
            # x$Curator_summary
          ),
          stringsAsFactors=FALSE
        )
      })
      annotTableReverse <- adply(thisAnnot[thisAnnot$strand=="-",],1,function(x) {
        data.frame(x=c(x$transcript_start,x$transcript_end,x$transcript_end),y=c(annotYvalReverse,annotYvalReverse,NA),url=paste0(araburlBase,x$Locus),
          name=sprintf("<table cellpadding='4' style='line-height:1.5'><tr><th>%1$s</th></tr><tr><td align='left'>Location: %2$s-%3$s<br>Chromosome: %4$s, Strand: %5$s<br>Short Desc.: %6$s</td></tr></table>",
            x$name,
            prettyNum(x$transcript_start, big.mark = ","),
            prettyNum(x$transcript_end, big.mark = ","),
            x$chromosome,                                                                           
            x$strand,
            x$short_description
            # x$Curator_summary
          ),
          stringsAsFactors=FALSE
        )
      })

    } else if (values[[jth_ref("organism", j)]] == "Medicago truncatula") { # strand is '+' or '-'
      annotTable <- adply(thisAnnot[thisAnnot$strand=="+",],1,function(x) {
        data.frame(x=c(x$transcript_start,x$transcript_end,x$transcript_end),y=c(annotYvalForward,annotYvalForward,NA),url=paste0(legumeInfo_urlBase, x$name, "/json"),
          name=sprintf("<table cellpadding='4' style='line-height:1.5'><tr><th>%1$s</th></tr><tr><td align='left'>Location: %2$s-%3$s<br>Chromosome: %4$s, Strand: %5$s<br>Desc: %6$s</td></tr></table>",
            x$name,
            prettyNum(x$transcript_start, big.mark = ","),
            prettyNum(x$transcript_end, big.mark = ","),
            x$chromosome,
            x$strand,
            x$description
          ),
          stringsAsFactors=FALSE
        )
      })
      annotTableReverse <- adply(thisAnnot[thisAnnot$strand=="-",],1,function(x) {
        data.frame(x=c(x$transcript_start,x$transcript_end,x$transcript_end),y=c(annotYvalReverse,annotYvalReverse,NA),url=paste0(legumeInfo_urlBase, x$name, "/json"),
          name=sprintf("<table cellpadding='4' style='line-height:1.5'><tr><th>%1$s</th></tr><tr><td align='left'>Location: %2$s-%3$s<br>Chromosome: %4$s, Strand: %5$s<br>Desc: %6$s</td></tr></table>",
            x$name,
            prettyNum(x$transcript_start, big.mark = ","),
            prettyNum(x$transcript_end, big.mark = ","),
            x$chromosome,
            x$strand,
            x$description
          ),
          stringsAsFactors=FALSE
        )
      })

    } else { #} if(values[[jth_ref("organism", j)]] == "Sorghum"){#strand is '+' or '-'
      annotTable <- adply(thisAnnot[thisAnnot$strand=="+",],1,function(x) {
        data.frame(x=c(x$transcript_start,x$transcript_end,x$transcript_end),y=c(annotYvalForward,annotYvalForward,NA),url=paste0(sorgurlBase,x$ID),
          name=sprintf("<table cellpadding='4' style='line-height:1.5'><tr><th>%1$s</th></tr><tr><td align='left'>Location: %2$s-%3$s<br>Chromosome: %4$s, Strand: %5$s<br>Desc: %6$s<br>Top TAIR Hit: %7$s<br>Top Rice Hit: %8$s</td></tr></table>",
            x$name,
            #1,#x$V2.1,
            prettyNum(x$transcript_start, big.mark = ","),
            prettyNum(x$transcript_end, big.mark = ","),
            x$chromosome,                                                                           
            x$strand,
            x$defLine,
            x$bestArabHitDefline,
            x$bestRiceHitDefline
            # x$Curator_summary
          ),
          stringsAsFactors=FALSE
        )
      })
      annotTableReverse <- adply(thisAnnot[thisAnnot$strand=="-",],1,function(x) {
        data.frame(x=c(x$transcript_start,x$transcript_end,x$transcript_end),y=c(annotYvalReverse,annotYvalReverse,NA),url=paste0(sorgurlBase,x$ID),
          name=sprintf("<table cellpadding='4' style='line-height:1.5'><tr><th>%1$s</th></tr><tr><td align='left'>Location: %2$s-%3$s<br>Chromosome: %4$s, Strand: %5$s<br>Desc: %6$s<br>Top TAIR Hit: %7$s<br>Top Rice Hit: %8$s</td></tr></table>",
            x$name,
            #1,#x$V2.1,
            prettyNum(x$transcript_start, big.mark = ","),
            prettyNum(x$transcript_end, big.mark = ","),
            x$chromosome,                                                                           
            x$strand,
            x$defLine,
            x$bestArabHitDefline,
            x$bestRiceHitDefline
            # x$Curator_summary
          ),
          stringsAsFactors=FALSE
        )
      })
    }
    #annotTable <- annotTable[,c("x","y","name","url","marker")]
    
    annotTable <- annotTable[,c("x","y","name","url")]
    #annotTable <- annotTable[order(annotTable$x),]
    
    #annotTableReverse <- annotTableReverse[,c("x","y","name","url","marker")]
    if(nrow(annotTableReverse)==0){
      annotTableReverse <- data.frame(x=character(0),y=character(0),name=character(0),url=character(0),stringsAsFactors = FALSE)
    }
    annotTableReverse <- annotTableReverse[,c("x","y","name","url")]
    #annotTableReverse <- annotTableReverse[order(annotTableReverse$x),]
    
    annotArray <- toJSONArray2(annotTable, json = F, names = T)
    #     for(i in 1:length(annotArray)){ #use this to add a symbol before or after the gene track
    #       if(is.na(annotArray[[i]]$marker)){
    #         annotArray[[i]]$marker <- NULL
    #       }else{
    #         annotArray[[i]]$marker <- NULL
    #         annotArray[[i]]$marker$symbol <- "url(./forwardArrow.svg)"
    #       }
    #     }

    annotArrayReverse <- toJSONArray2(annotTableReverse, json = F, names = T)
    #     for(i in 1:length(annotArrayReverse)){
    #       if(is.na(annotArrayReverse[[i]]$marker)){
    #         annotArrayReverse[[i]]$marker <- NULL
    #       }else{
    #         annotArrayReverse[[i]]$marker <- NULL
    #         annotArrayReverse[[i]]$marker$symbol <- "url(./reverseArrow.svg)"
    #       }
    #     }
    
    
    b <- rCharts::Highcharts$new()
    b$LIB$url <- 'highcharts/'
    #b$xAxis(title=list(text=zoomTitle), min= -0, max=1,startOnTick = FALSE,reversed=FALSE)
    #b$yAxis(title = list(text = "Base Pairs"),min=winLow,max=winHigh,startOnTick=TRUE,endOnTick=TRUE)
    #b$xAxis(list(title=list(text="RMIP"),min=-0,max=1,startOnTick=FALSE,reversed=FALSE),list(title=list(text="F-score"),min=0,max=1,startOnTick=FALSE,opposite=TRUE,reversed=FALSE))      
    b$chart(zoomType="x",alignTicks=FALSE,events=list(click = "#!function(event) {this.tooltip.hide();}!#"))
    b$xAxis(title = list(text = "Base Pairs"),startOnTick=FALSE,min=winLow,max=winHigh,endOnTick=FALSE)      
    #     if(winTable$sumRMIP[1] != -5){
    #       b$yAxis(list(title=list(text="RMIP"),min=-0,max=1,startOnTick=FALSE),list(title=list(text="F-score"),min=0,startOnTick=FALSE,opposite=TRUE,gridLineWidth=0),
    #               list(title=list(text="sumRMIP"),min=-0,startOnTick=FALSE,opposite=TRUE,gridLineWidth=0))          
    #     }else{
    #      b$yAxis(list(title=list(text="RMIP"),min=-0,max=1,startOnTick=FALSE),list(title=list(text="F-score"),min=0,startOnTick=FALSE,opposite=TRUE,gridLineWidth=0))
    #    }
    if(input[[jth_ref("axisLimBool", j)]] == TRUE){
      b$yAxis(title=list(text=input[[jth_ref("yAxisColumn", j)]]),min=input[[jth_ref("axisMin", j)]],max=input[[jth_ref("axisMax", j)]],startOnTick=FALSE)
      #create a hidden axis to put the gene track on, all the options are setting to hide everything from the axis 
      b$yAxis(labels=list(enabled=FALSE),title=list(text=NULL),min=0,max=1,lineWidth=0,gridLineWidth=0,minorGridLineWidth=0,lineColor="transparent",minorTickLength=0,tickLength=0,startOnTick=FALSE,opposite=TRUE,replace=FALSE)
    }else{      
      b$yAxis(title=list(text=input[[jth_ref("yAxisColumn", j)]]),startOnTick=FALSE) 
      #create a hidden axis to put the gene track on, all the options are setting to hide everything from the axis
      b$yAxis(labels=list(enabled=FALSE),title=list(text=NULL),min=0,max=1,lineWidth=0,gridLineWidth=0,minorGridLineWidth=0,lineColor="transparent",minorTickLength=0,tickLength=0,startOnTick=FALSE,opposite=TRUE,replace=FALSE)
    }
    
    if(input[[jth_ref("supportInterval", j)]]==TRUE){
      if(input[[jth_ref("SIaxisLimBool", j)]] == TRUE){
        b$yAxis(title=list(text=input[[jth_ref("SIyAxisColumn", j)]]),min=input[[jth_ref("SIaxisMin", j)]],max=input[[jth_ref("SIaxisMax", j)]],gridLineWidth=0,minorGridLineWidth=0,startOnTick=FALSE,opposite=TRUE,replace=FALSE)
      }else{
        b$yAxis(title=list(text=input[[jth_ref("SIyAxisColumn", j)]]),gridLineWidth=0,minorGridLineWidth=0,startOnTick=FALSE,opposite=TRUE,replace=FALSE)   
      }
      
      if(SIchart[1,input[[jth_ref("SIyAxisColumn", j)]]] != -1){
        d_ply(jlTable,.(trait),function(x){
          b$series(
            data = toJSONArray2(x,json=F,names=T),
            type = "line",
            dashStyle = 'Solid',
            lineWidth = 10,
            name = unique(x$trait),
            yAxis=2,           
            color = colorTable$color[colorTable$trait == as.character(unique(x$loc_el))])})            
      }
    }
    
    #b$series(
    #  data = toJSONArray2(annotTable, json = F, names = F),
    #  type = "columnrange",
    #  pointWidth=15,
    #  pointPadding=0,
    #  pointPlacement=0,
    #  name="Genes"    
    #)    
    
    #invisible(sapply(jlSeries, function(x) {b$series(data = x, type = "line", yAxis=1,name = x[[1]]$element)}))
    #     if(jl$F[1] != -5){
    #       d_ply(jlTable,.(element),function(x){
    #         b$series(
    #           data = toJSONArray2(x,json=F,names=T),
    #           type = "line",
    #           name = unique(x$element),
    #           yAxis=1,
    #           visible = if(input$JLCheckbox==FALSE) FALSE else TRUE,
    #           color = colorTable$color[colorTable$loc_el == as.character(unique(x$element))])})
    #     }
    #     
    #     if(winTable$sumRMIP[1] != -5){      
    #       d_ply(winTableSeries,.(element),function(x){
    #         b$series(
    #           data = toJSONArray2(x,json=F,names=T),
    #           type = "line",
    #           name = unique(x$element),
    #           yAxis = 2,
    #           marker = list(enabled = FALSE),
    #           lineWidth = 3,
    #           turboThreshold=5500,
    #           color = colorTable$color[colorTable$loc_el == as.character(sub("_Window","",unique(x$element)))]
    #         )
    #       })
    #     }
    #     
    if(zoomChart[1,input[[jth_ref("yAxisColumn", j)]]] != -1){
      invisible(sapply(zoomSeries, function(x) {if(length(x)==0){return()};b$series(data = x, type = "scatter", showInLegend = FALSE, color = colorTable$color[colorTable$trait == as.character(x[[1]]$trait)], name = paste0(x[[1]]$trait))}))
    }
    
    b$series(
      data = annotArray,
      type = "line",
      showInLegend = FALSE,
      name = "Forward Genes",
      id = "forward-genes",
      zIndex = 1,
      color = "#53377A",
      marker = list(symbol = "circle", enabled = FALSE),
      yAxis = 1
    )    
    
    b$series(
      data = annotArrayReverse,
      type = "line",
      showInLegend = FALSE,
      name = "Reverse Genes",
      id = "reverse-genes",
      zIndex = 1,
      color = "#53377A",
      marker = list(symbol = "circle", enabled = FALSE),
      yAxis = 1
    )      
    
    if (!is.null(values[[jth_ref("glGenes", j)]])) {
      apply(values[[jth_ref("glGenes", j)]], 1, FUN = function(g) {
        g.strand <- as.integer(g$strand)
        yh <- -1
        if (g.strand == 1) {
          yh <- annotYvalForward
          sid <- "forward-genes"
        } else if (g.strand == -1) {
          yh <- annotYvalReverse
          sid <- "reverse-genes"
        }
        if (yh > 0 && g$chr == input[[jth_ref("chr", j)]]) {
          g.data <- vector("list", 2)
          g.data[[1]]$x <- as.integer(g$fmin)
          g.data[[2]]$x <- as.integer(g$fmax)
          g.data[[1]]$y <- g.data[[2]]$y <- yh
          b$series(
            type = "line",
            data = g.data,
            color = g$color,
            linkedTo = sid,
            lineWidth = 12,
            yAxis = 1,
            zIndex = 0
          )
        }
      })
    }
    
    #b$series(data = annotSeries[[1]], type = "columnrange", pointWidth=15,pointPadding=0,pointPlacement=0,name = "Genes")    
    b$chart(zoomType="x",alignTicks=FALSE,events=list(click = "#!function(event) {this.tooltip.hide();}!#"))
    #b$title(text=paste("NAM GWAS Results",sep=" "))
    #b$subtitle(text="Rollover for more info. Drag chart area to zoom. Click point for zoomed annotated plot.")

    # User clicked on a point -> display trait in popup
    doClickOnPoint <- "#! function(event) { alert(this.trait); } !#"
    # User clicked on a line -> various possible responses:
    microSyntenySearch <- paste(
      "$.ajax({",
        "url: 'https://' + url2 + '/lis_context_server/services/v1/micro-synteny-search/',",
        "dataType: 'json',",
        "data: JSON.stringify({",
          "query: families1,",
          "matched: $('input#matched').val(),",
          "intermediate: $('input#intermediate').val()",
        "}),",
        "type: 'POST',",
        "success: function(response2) {",
          "obj2 = JSON.parse(response2);",
          # Send information about neighboring and related genes back to the chart
          "Shiny.onInputChange('genomicLinkages', {",
            "results1: obj1,",
            "results2: obj2",
          "});",
        "},",
        "error: function(errmsg2) { alert('FAIL2: ' + errmsg2.responseText); }",
      "});"
    )
    geneToQueryTrack <- paste(
      # Query the Genome Context Viewer for genomic linkages
      "var url1 = '';",
      "var url2 = '';",
      "var speciesName2 = '';",
      "var geneString = '';",
      "mt0 = this.url.search('medtr');",
      "mt1 = this.url.search('/json');",
      "if (mt0 >= 0) {",
        "geneString = this.url.substring(mt0, mt1);",
        "url1 = 'legumeinfo.org';",
        "url2 = 'legumefederation.org';",
        "speciesName2 = 'A.thaliana';",
      "} else if (this.url.search('arabidopsis.org') >= 0) {",
        "at0 = this.url.search('name=');",
        "geneString = 'arath.Col.' + this.url.substring(at0 + 5);",
        "url1 = 'legumefederation.org';",
        "url2 = 'legumeinfo.org';",
        "speciesName2 = 'M.truncatula';",
      "} else {",
        "return;",
      "}",
      "$.ajax({",
        "url: 'https://' + url1 + '/lis_context_server/services/v1/gene-to-query-track/',",
        "dataType: 'json',",
        "data: JSON.stringify({",
          "gene: geneString,",
          "neighbors: $('input#neighbors').val()",
        "}),",
        "type: 'POST',",
        "success: function(response) {",
          "obj1 = JSON.parse(response);",
          "families1 = Array.from(obj1.genes, x => x.family);",
          microSyntenySearch,
        "},",
        "error: function(errmsg) { alert('FAIL: ' + errmsg.responseText); }",
      "});"
    )
    provideMultipleURLs <- paste(
      # From the JSON at this.url, extract the URLs related to this gene.
      # Note that this.url = legumeInfo_urlBase + geneString + '/json'
      #  legumeInfo_urlBase currently has 34 characters (see above)
      #  and geneString = <5-character species abbreviation>.geneName
      # And for now, add the gene family phylogram URL by hand.
      "$.getJSON(this.url, function(data) {",
        "var geneString = this.url.substring(34, this.url.indexOf('/json'));",
        "var geneName = geneString.substring(6);",
        "var content = '';",
        "if (data.length == 0) {",
          "content = '<p>No ' + geneName + ' links found.</p>';",
        "} else {",
          "$.each(data, function(i, obj) {",
            "content = content + '<p><a href=' + obj.href + ' target=_blank>' + obj.text + '</a></p>';",
            "if (i == 0) {",
              "var urlPhylogram = 'http://legumeinfo.org/chado_gene_phylotree_v2?gene_name=' + geneString;",
              "var textPhylogram = 'View LIS gene family phylogram page for : ' + geneName;",
              "content = content + '<p><a href=' + urlPhylogram + ' target=_blank>' + textPhylogram + '</a></p>';",
            "}",
          "});",
        "}",
        "var $div = $('<div></div>');",
        "$div.html(content);",
        "$div.dialog({",
          "title: geneName + ' Links',",
          "width: 512,",
          "height: 'auto',",
          "modal: true",
        "});",
      "});"
    )
    bGenomicLinkage <- ifelse(j == 1, 1, 0)
    doClickOnLine <- sprintf(paste(
      "#! function() {",
        "if (%d && $('input#boolGenomicLinkage').prop('checked')) {",
          geneToQueryTrack,
        "} else if (this.url.includes('legumeinfo.org')) {",
          provideMultipleURLs,
        "} else {",
          # for all other species
          "window.open(this.url);", #open webpage
        "}",
      "} !#"
    ), bGenomicLinkage)

    # Create the gene family legend. For organism 2, show only the gene families on the currently selected chromosome.
    glFamilies <- names(values$glColors)
    if (j == 2 && !is.null(values$glGenes2) && !(is.null(input$relatedRegions) || length(input$relatedRegions) == 0)) {
      glGenes2 <- values$glGenes2
      # parse from the format "chr[Chr] [minBP]-[maxBP] Mbp"
      ss <- strsplit(input$relatedRegions, split = " ")[[1]]
      chr <- as.integer(stri_sub(ss[1], 4))
      glFamilies <- intersect(glFamilies, glGenes2$family[glGenes2$chr == chr])
    }
    doClickOnColumn <- paste(
      "#! function() {",
        "window.open('https://legumeinfo.org/chado_phylotree/' + this.name);", # go to the gene family's web page
        "return false;", # and disable toggling the legend item
      "} !#"
    )
    sapply(glFamilies, FUN = function(f) {
      b$series(
        data = list(),
        type = "column",
        name = f,
        color = values$glColors[[f]],
        events = list(legendItemClick = doClickOnColumn)
      )
    })

    b$plotOptions(
      scatter = list(
        cursor = "pointer",
        point = list(
          events = list(
            #click = "#! function() { window.open(this.options.url); } !#")), #open webpage
            click = doClickOnPoint
          )
        ),
        #click = "#! function(event) {console.log(this);} !#")), #write object to log
        #click = "#! function(){$('input#selected').val(134); $('input#selected').trigger('change');} !#")),
        #click = "#! function(){$('input#selected').val(this.options.bp); $('input#selected').trigger('change');} !#")),
        marker = list(
          symbol = "circle",
          radius = 5
        )
      ),
      line = list(
        lineWidth = 6,
        cursor = "pointer",
        #stickyTracking=FALSE,
        point = list(
          events = list(
            click = doClickOnLine
          )
        )
      ),
      #click = "#! function(event) {alert(this.url);} !#")), #display popup
      #click = "#! function(event) {console.log(this);} !#")), #write object to log
      #click = "#! function(){$('input#selected').val(134); $('input#selected').trigger('change');} !#")),
      #click = "#! function(){$('input#selected').val(this.options.bp); $('input#selected').trigger('change');} !#")),
      marker = list(
        enabled = FALSE,
        radius = 1,
        #symbol = "url(./sun.png)",
        states = list(hover = list(enabled=FALSE))
      )
    )

    #it seems almost impossible to get the tooltip to hover along the chart with this version of highcharts (4.0.1), perhaps a question to stackoverflow could solve it.
    #see an example of the problem here: http://jsfiddle.net/N5ymb/
    #one hack/fix would be to add dummy points to the middle of the line that show up when moused over
    b$tooltip(snap=5, useHTML = T, formatter = "#! function() { return this.point.name; } !#") #followTouchMove = T, shared=T, followPointer = T
    b$exporting(enabled=TRUE,filename='zoomChart',sourceWidth=2000)
    b$credits(enabled=TRUE)
    b$set(dom = jth_ref('zChart', j))
    return(b)
  }
  output$zChart <- renderChart(create_zChart(1))
  output$zChart2 <- renderChart(create_zChart(2))

  #highcharts test chart
  output$testChart <- renderChart({
    h1 <- hPlot(x = "Wr.Hnd", y = "NW.Hnd", data = MASS::survey, type = c("line", 
                                                                          "bubble", "scatter"), group = "Clap", size = "Age")
    h1$set(dom = 'testChart')
    return(h1)     
  })
  
  #return color table (which sets colors for series in charts) up to date with all combinations of traits
  #this uses allColors set up in global, colors will repeat after 30series
  getColorTable <- function(j) {
    traitVals <- list()
    if(input[[jth_ref("plotAll", j)]] == FALSE){
      for(i in input[[jth_ref("traitColumns", j)]]){
        traitVals[[i]] <- input[[jth_ref(i, j)]]
      }
      
      traits <- do.call(paste,c(expand.grid(traitVals),sep="_"))     
      if(length(traits)==0){return(NULL)}
      
      colorTable <- data.frame(trait=traits,color=rep(allColors,ceiling(length(traits)/30))[1:length(traits)])
    }else{
      colorTable <- data.frame(trait=input[[jth_ref("datasets", j)]],color=allColors[1])
    }
    colorTable
  }

  findGWASOverlaps <- function(genomeChart, j) {
    if(is.null(input[[jth_ref("overlapSize", j)]])){return(genomeChart[1,])}
    tableIn <- genomeChart
    tableIn$winStart <- tableIn[, input[[jth_ref("bpColumn", j)]]] - input[[jth_ref("overlapSize", j)]]
    tableIn$winStop <- tableIn[, input[[jth_ref("bpColumn", j)]]] + input[[jth_ref("overlapSize", j)]]
    
    allGr <- GRanges(tableIn[,input[[jth_ref("chrColumn", j)]]], IRanges(start=tableIn$winStart,end=tableIn$winStop))

    tableIn$group <- subjectHits(findOverlaps(allGr, reduce(allGr)))
    
    #just groups that have more than one unique SNP
    gwasDataOverlap <- tableIn[tableIn$group %in% as.data.frame(table(tableIn$group))[as.data.frame(table(tableIn$group))$Freq>1,"Var1"],]
    
    #just groups that have more than one unique phenotype
    gwasDataOverlapDiffPheno <- ddply(gwasDataOverlap,.(group),function(x){if(nrow(unique(as.data.frame(x[,"trait"])))>=input[[jth_ref("numOverlaps", j)]]){x}else{x[0,]}})

    return(gwasDataOverlapDiffPheno)    
  }

  handleSubmitColsButton <- function(j) {
    if(is.null(input[[jth_ref("SubmitColsButton", j)]]) || input[[jth_ref("SubmitColsButton", j)]] == 0){return()}
    isolate({
      currDatasetProp <- datasetProp()
      #print("before")
      #print(currDatasetProp)
      if(as.character(input[[jth_ref("datasets", j)]]) %in% currDatasetProp$dataset){
        currDatasetProp <- currDatasetProp[currDatasetProp$dataset != as.character(input[[jth_ref("datasets", j)]]),]
      }
      #print("after")
      #print(currDatasetProp)
      cols <- varnames(j)
      #print("data.frame")
      #print(data.frame(dataset=input$datasets,chrColumn=names(cols[cols==input$chrColumn]),bpColumn=names(cols[cols==input$bpColumn]),
      #                 traitCol=paste(names(cols[cols %in% input$traitColumns]),collapse=";"),yAxisColumn=names(cols[cols==input$yAxisColumn]),axisLim=input$axisLimBool,axisMin=input$axisMin,axisMax=input$axisMax,stringsAsFactors=FALSE))
      currDatasetProp <-  rbind(currDatasetProp,data.frame(dataset=input[[jth_ref("datasets", j)]],chrColumn=names(cols[cols==input[[jth_ref("chrColumn", j)]]]),bpColumn=names(cols[cols==input[[jth_ref("bpColumn", j)]]]),
        traitCol=paste(names(cols[cols %in% input[[jth_ref("traitColumns", j)]]]),collapse=";"),yAxisColumn=names(cols[cols==input[[jth_ref("yAxisColumn", j)]]]),
        logP=input[[jth_ref("logP", j)]],axisLim=input[[jth_ref("axisLimBool", j)]],axisMin=input[[jth_ref("axisMin", j)]],axisMax=input[[jth_ref("axisMax", j)]],organism=values[[jth_ref("organism", j)]],plotAll=input[[jth_ref("plotAll", j)]],
        supportInterval=input[[jth_ref("supportInterval", j)]],SIyAxisColumn=input[[jth_ref("SIyAxisColumn", j)]],SIbpStart=input[[jth_ref("SIbpStart", j)]],SIbpEnd=input[[jth_ref("SIbpEnd", j)]],
        SIaxisLimBool=input[[jth_ref("SIaxisLimBool", j)]],SIaxisMin=input[[jth_ref("SIaxisMin", j)]],SIaxisMax=input[[jth_ref("SIaxisMax", j)]],stringsAsFactors=FALSE))      
      #print("rbind")
      #print(currDatasetProp)
      write.table(file="./www/config/datasetProperties.csv",x=currDatasetProp,col.names=TRUE,row.names=FALSE,sep=",")
      updateTabsetPanel(session, "datatabs", selected = "WhGen")
    })
#    if(input$selected != 1e5){
#      updateTabsetPanel(session, "datatabs", selected = "WhGen")  
#    }
  }
  observe(handleSubmitColsButton(1))
  observe(handleSubmitColsButton(2))

  saveDataset <- function(j) {
    if (is.null(input[[jth_ref("saveDatasetButton", j)]]) || input[[jth_ref("saveDatasetButton", j)]] == 0) { return() }
    isolate({
      if(!file.exists(paste0("./www/config/data/",input[[jth_ref("datasets", j)]]))){
        write.table(getdata(j),paste0("./www/config/data/",input[[jth_ref("datasets", j)]]),sep=",",col.names=TRUE,row.names=FALSE)
      }
    })
  }
  observe(saveDataset(1))
  observe(saveDataset(2))
  
   datasetProp <- function(){
     dp <- read.table("./www/config/datasetProperties.csv",sep=",",head=TRUE,stringsAsFactors=FALSE)
     if (is.null(values$datasetToOrganism)) {
       # populate the datasetToOrganism map from stored values
       values$datasetToOrganism <- list()
       for (i in 1:nrow(dp)) {
         values$datasetToOrganism[[dp$dataset[i]]] <- dp$organism[i]
       }
     }
     return(dp)
   }
   #this provides the functionality to update the plotband window after a user clicks without rerendering the whole plot
   #From:
   #http://stackoverflow.com/questions/20247759/add-highcharts-plotband-after-render-in-r-shiny-rcharts/20249933?noredirect=1#20249933
   updatePlotbandWindow <- function(j) {
     center <- as.numeric(input[[jth_ref("selected", j)]][[1]])
     winHigh <- center + input[[jth_ref("window", j)]][1]
     winLow <- center - input[[jth_ref("window", j)]][1]
     #eventually I would use winLow/winHigh to change the plotband range
     band = list(from = winLow, to = winHigh, color = "rgba(68, 170, 213, 0.4)")
     #print(band)
     session$sendCustomMessage(type = jth_ref("customMsg", j), band)
   }
   observe(updatePlotbandWindow(1))
   observe(updatePlotbandWindow(2))

  output$selectedGene <- renderUI({
    h5(values$glSelectedGene)
  })
  output$relatedRegions <- renderUI({
    selectInput("relatedRegions", "Related Regions:", choices = NULL, selectize = FALSE)
  })
  observe({
    if (is.null(input$relatedRegions) || length(input$relatedRegions) == 0) {
      updateNumericInput(session, "selected2", value = input$selected2)
    } else {
      # parse from the format "chr[Chr] [minBP]-[maxBP] Mbp"
      ss <- strsplit(input$relatedRegions, split = " ")[[1]]
      chr <- as.integer(stri_sub(ss[1], 4))
      ss2 <- strsplit(ss[2], split = "-")[[1]]
      centerBP <- as.integer(1.0e6*mean(as.numeric(ss2)))
      updateSelectInput(session, "chr2", selected = chr)
      updateNumericInput(session, "selected2", value = centerBP)
    }
  })

  clearGenomicLinkages <- function() {
    values$glGenes <- values$glGenes2 <- values$glColors <- NULL
    updateSelectInput(session, "relatedRegions", choices = character(0))
  }

  # If the user unchecks Genomic Linkage options, clear any existing genomic linkage query results
  observe({
    glOn <- input$boolGenomicLinkage
    if (!(is.null(glOn) || glOn)) {
      values$glSelectedGene <- NULL
      clearGenomicLinkages()
    }
  })

  observe({
    # Handle and display genomic linkage query results
    if (is.null(input$genomicLinkages)) return()

    # Parse neighboring genes from species 1
    results1 <- input$genomicLinkages$results1
    values$glSelectedGene <- results1$genes[[(length(results1$genes) + 1) %/% 2]]$name
    glGenes <- data.frame(matrix(unlist(results1$genes), nrow = length(results1$genes), byrow = TRUE),
      stringsAsFactors = FALSE)[, 3:6]
    glGenes$chr <- as.integer(stri_match(results1$chromosome_name, regex = "(?i)(?<=\\.chr)\\d+$")[, 1])
    names(glGenes) <- c("family", "fmin", "fmax", "strand", "chr")
    glGenes <- glGenes[nchar(glGenes$family) > 0, ]
    if (nrow(glGenes) == 0) {
      # could reach here if none of the (2*neighbors + 1) genes has a family id
      clearGenomicLinkages()
      return()
    }

    # Convert (for example) "Medicago truncatula" to "M.truncatula"
    ss.org2 <- strsplit(values$organism2, split = " ")[[1]]
    abbrSpeciesName2 <- paste(stri_sub(ss.org2[1], 1, 1), ss.org2[2], sep = ".")
    # Parse related genes from species 2
    results2 <- input$genomicLinkages$results2
    if (length(results2$groups) == 0) {
      clearGenomicLinkages()
      return()
    }
    for (i in 1:length(results2$groups)) {
      results2$groups[[i]]$id <- i
    }
    glGenes2 <- do.call(rbind, lapply(results2$groups, FUN = function(gr) {
      gr.chr <- as.integer(stri_match(gr$chromosome_name, regex = "(?i)(?<=\\.chr)\\d+$")[, 1])
      if (gr$species_name == abbrSpeciesName2 && !is.na(gr.chr)) {
        gr.genes <- data.frame(matrix(unlist(gr$genes), nrow = length(gr$genes), byrow = TRUE),
          stringsAsFactors = FALSE)[, 3:6]
        gr.genes$chr <- gr.chr
        gr.genes$id <- gr$id
        names(gr.genes) <- c("family", "fmin", "fmax", "strand", "chr", "id")
        gr.genes <- gr.genes[nchar(gr.genes$family) > 0, ]
        gr.genes
      }
    }))

    # Highlight families common to both genomes
    families <- intersect(glGenes$family, glGenes2$family)
    nf <- length(families)
    if (nf == 0) {
      clearGenomicLinkages()
      return()
    }
    # Create nf colors
    fc <- rainbow(nf, end = 5/6) # TODO: a more clearly distinguishable set of colors
    familyColors <- list()
    lapply(1:nf, FUN = function(i) familyColors[[families[i]]] <<- stri_sub(fc[i], 1, 7))

    # Construct the chart data
    glGenes <- glGenes[glGenes$family %in% families, ]
    glGenes$color <- familyColors[glGenes$family]
    glGenes2 <- glGenes2[glGenes2$family %in% families, ]
    glGenes2$color <- familyColors[glGenes2$family]

    # Construct the related regions (each corresponds to a group from results2$groups)
    glGroupIds <- unique(glGenes2$id)
    glRelatedRegions <- do.call(rbind.data.frame, compact(lapply(results2$groups, FUN = function(gr) {
      if (gr$id %in% glGroupIds) {
        gr.chr <- as.integer(stri_match(gr$chromosome_name, regex = "(?i)(?<=\\.chr)\\d+$")[, 1])
        gr.minBP <- gr$genes[[1]]$fmin
        gr.maxBP <- gr$genes[[length(gr$genes)]]$fmax
        list(region = sprintf("chr%d %3.2f-%3.2f Mbp", gr.chr, gr.minBP*1.0e-6, gr.maxBP*1.0e-6),
          chr = gr.chr, minBP = gr.minBP, maxBP = gr.maxBP)
      }
    })))
    # Sort the related regions
    glRelatedRegions <- glRelatedRegions[with(glRelatedRegions, order(chr, minBP)), ]

    # Recenter the window around the selected gene
    centerBP1 <- (as.integer(glGenes$fmin[1]) + as.integer(glGenes$fmax[nrow(glGenes)])) %/% 2
    updateNumericInput(session, "selected", value = centerBP1)
    # Update the charts and other output
    values$glGenes <- glGenes
    values$glGenes2 <- glGenes2
    values$glColors <- familyColors
    updateSelectInput(session, "relatedRegions", choices = glRelatedRegions$region)
  })

})#end server
