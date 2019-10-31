source("gChart.R")
source("pChart.R")
source("zChart.R")

shinyServer(function(input, output, session) {
  #Load any saved datasets
  values <- reactiveValues()
  dataPath <- "./www/config/data/"
  dataFiles <- list.files(dataPath,recursive=T)
  # Append names of any data we will create on the fly:
  # Add your organism to legumeInfo.gwas if its GWAS files live on a server instead of locally.
  legumeInfo.gwas <- c("Arabidopsis thaliana GWAS", "Medicago truncatula GWAS")
  # TODO: Do we really need legumeInfo.organisms?
  legumeInfo.organisms <- c("Arabidopsis thaliana", "Medicago truncatula", "Soybean", "Cowpea", "Pigeonpea")
  dataFiles <- c(dataFiles, legumeInfo.gwas)
  for(i in dataFiles){
    if (i %in% legumeInfo.gwas) {
      values[[i]] <- init.gwas(i)
    } else {
      values[[i]] <- read.table(paste0(dataPath,i),sep=",",stringsAsFactors=FALSE,head=TRUE)
    }
  }  
  values$datasetlist <- dataFiles
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
    trait.choices <- gwas.traits[[values[[jth_ref("organism", j)]]]]
    if (is.null(trait.choices)) trait.choices <- character(0)
    updateSelectizeInput(session, jth_ref("gwasTraits", j), choices = trait.choices, selected = NULL)
    # As there is no updateFileInput() method,
    # send a custom message (to clear the progress bar) and clear values$needsToUploadFiles
    session$sendCustomMessage(type = "resetFileInputHandler", jth_ref("uploadfile", j))
    values[[jth_ref("needsToUploadFiles", j)]] <- FALSE
    updateCheckboxInput(session, jth_ref("appendSNPs", j), value = TRUE)
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
      tags$div(id = jth_ref("tour-dataset", j), wellPanel(
        style = paste0("background-color: ", bgColors[j], ";"),
        uiOutput(jth_ref("datasets", j))
      )),
      tags$div(id = jth_ref("tour-loaddata", j), wellPanel(
        style = paste0("background-color: ", bgColors[j], ";"),
        radioButtons(inputId = jth_ref("dataType", j), label = "Load data (Max. 5MB):", c(".csv" = "csv", ".rda" = "rda", "examples" = "examples"), selected = "csv"),
        conditionalPanel(condition = paste0("input.", jth_ref("dataType", j), " != 'examples'"),
          conditionalPanel(condition = paste0("input.", jth_ref("dataType", j), " == 'csv'"),
            checkboxInput(jth_ref('header', j), 'Header', TRUE),
            radioButtons(jth_ref('sep', j), '', c(Comma=',', Semicolon=';', Tab='\t'), ',')
          ),
          selectizeInput(jth_ref("gwasTraits", j), "Remote Trait Files:", choices = NULL, multiple = TRUE),
          fileInput(jth_ref("uploadfile", j), "Local Trait Files:", multiple = TRUE),
          checkboxInput(jth_ref("appendSNPs", j), "Append to current dataset (combine multiple traits)", TRUE),
          actionButton(jth_ref("loadTraits", j), "Load Data")
        ),
        conditionalPanel(condition = paste0("input.", jth_ref("dataType", j), " == 'examples'"),
          actionButton(jth_ref('loadExampleData', j), 'Load examples')
        )
      )) #,
      # wellPanel(
      #   style = paste0("background-color: ", bgColors[j], ";"),
      #   h6("Once your file is finished uploading, press the Save Dataset button below and reload ZBrowse."),
      #   actionButton(jth_ref('saveDatasetButton', j), 'Save Current Dataset'),
      #   conditionalPanel(condition = paste0("input.", jth_ref("saveDatasetButton", j), " > 0"),
      #     h5("Dataset successfully saved!")
      #   )
      # )
    )
  }
  observeEvent(input$appendSNPs, shinyjs::disable("appendSNPs"))
  observeEvent(input$appendSNPs, shinyjs::disable("appendSNPs2"))
  
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
    tags$div(id = "tour-genLink", wellPanel(
      h5("Genomic Linkage options:"),
      checkboxInput('boolGenomicLinkage', 'ON', FALSE),
      conditionalPanel("input.boolGenomicLinkage == true",
        h5("Broadcast Channel options:"),
        checkboxInput('boolBroadcastToBC', 'Broadcast', TRUE),
        checkboxInput('boolListenToBC', 'Listen', TRUE),
        tags$div(id = "tour-genLink-1", wellPanel(
          uiOutput("selectedGene"),
          numericInput("neighbors", "Neighbors:", min = 1, max = 20, value = 20),
          numericInput("matched", "Matched:", min = 1, max = 20, value = 4),
          numericInput("intermediate", "Intermediate:", min = 1, max = 10, value = 5),
          conditionalPanel("input.boolBroadcastToBC == true",
            actionLink("viewInGCV", "View in GCV")
          ),
          style = paste0("background-color: ", bgColors[1], ";")
        )),
        tags$div(id = "tour-genLink-2", wellPanel(
          uiOutput("relatedRegions"),
          style = paste0("background-color: ", bgColors[2], ";")
        ))
      )
    ))
  }

  displayConnectionStatus <- function() {
    apply(gwas.sources, FUN = function(ci) {
      # id for each connection's label
      nn <- paste("gwasSource", ci["name"], sep = "_")
      # server part comes first
      output[[nn]] <- renderText({
        ci["status"] <- url.exists(ci["url"])
        if (ci["status"]) {
          cc <- "green"
          ss <- ""
        } else {
          cc <- "red"
          ss <- " unavailable"
        }
        tt <- 60000 # msec between connection tests
        invalidateLater(tt, session)
        paste0("<span style=\"color:", cc, "\">", ci["name"], ss, "</span>")
      })
      # then the ui part
      htmlOutput(nn)
    }, MARGIN = 1)
  }
  output$ui_All <- renderUI({
    list(
      conditionalPanel(condition = manageTabSelected,
        actionLink("zbrowseTour", "Start Tour"),
        wellPanel(h5("Connection Status"), displayConnectionStatus()),
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
        helpModal('Manage','manage', includeMarkdown("tools/manage.md")),
        HTML(paste('<p style="font-size:10px;">Powered by',
          '<a href="http://www.rstudio.com/shiny/">Shiny</a>,',
          '<a href="http://rcharts.io/">rCharts</a>,',
          '<a href="http://www.highcharts.com">Highcharts</a>,',
          'and <a href="https://github.com/carlganz/rintrojs">rintrojs</a>',
          '</p>'))
      ),#end conditional Manage

      conditionalPanel(condition = dataTableTabSelected,
        tags$div(id = "tour-datatableSidebar", createDataTableSidebar(1)),
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
    thisChrAnnot <- subset(org.annotGeneLoc[values[[jth_ref("organism", j)]]][[1]], chromosome == input[[jth_ref("chr", j)]])
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
    tags$div(id = jth_ref("tour-columnSettings", j), wellPanel(
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
    ))
  }
  output$ui_data_tabs <- renderUI({
    tabsetPanel(id = "datatabs",      
      tabPanel(title="Manage",value="Manage",
        createColumnSettingsPanel(1),
        createColumnSettingsPanel(2)
      ),
      tabPanel(title="Data Table",value="Table",
        tags$div(id = "tour-datatable",
          wellPanel(dataTableOutput("dataviewer"), style = paste0("background-color: ", bgColors[1], ";"))
        ),
        wellPanel(dataTableOutput("dataviewer2"), style = paste0("background-color: ", bgColors[2], ";"))
      ),
      tabPanel(title="Whole Genome View",value="WhGen",
        tags$div(id = "tour-wholegenome",
          wellPanel(showOutput("gChart", "highcharts"), style = paste0("background-color: ", bgColors[1], ";"))
        ),
        wellPanel(showOutput("gChart2", "highcharts"), style = paste0("background-color: ", bgColors[2], ";"))
      ),
      tabPanel(title="Chromosome View",value="Chrom",
        wellPanel(
          tags$div(id = "tour-pChart", showOutput("pChart", "highcharts")),
          tags$div(id = "tour-zChart", showOutput("zChart", "highcharts")),
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
        tags$div(id = "tour-annotations",
          wellPanel(dataTableOutput("annotViewer"), style = paste0("background-color: ", bgColors[1], ";"))
        ),
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
  }, options = annotViewer.options, escape = FALSE)
  output$annotViewer2 <- renderDataTable({
    createAnnotTable(2)
  }, options = annotViewer.options, escape = FALSE)

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
  output$dataviewer <-renderDataTable(createDataViewer(1), options = dataViewer.options, escape = FALSE)
  output$dataviewer2 <-renderDataTable(createDataViewer(2), options = dataViewer.options, escape = FALSE)
  
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
      traits <- c("Select All", "Deselect All", sort(unique(values[[input[[jth_ref("datasets", j)]]]][,i])))
      selectizeInput(inputId=jth_ref(i, j), label=paste0("Select ",i),traits,
        selected = traits[1], # 1 for Select All, 2 for Deselect All, 3 for the first trait, etc
        multiple=TRUE, options = list(dropdownParent="body",plugins=list("remove_button")))
    })
  }
  output$traitColBoxes <- renderUI(createTraitColBoxes(1))
  output$traitColBoxes2 <- renderUI(createTraitColBoxes(2))

  updateTraitsMenu <- function(j) {
    lapply(input[[jth_ref("traitColumns", j)]], function(i){
      if ("Select All" %in% input[[jth_ref(i, j)]]) {
        selected_choices <- sort(unique(values[[input[[jth_ref("datasets", j)]]]][,i]))
        updateSelectizeInput(session, jth_ref(i, j), selected = selected_choices)
      } else if ("Deselect All" %in% input[[jth_ref(i, j)]]) {
        updateSelectizeInput(session, jth_ref(i, j), selected = character(0))
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

  #this function makes the chromsomeview chart  
  #subset whole chart based on selection
  output$pChart <- renderChart(create_pChart(1, input, values))
  output$pChart2 <- renderChart(create_pChart(2, input, values))
  
  #Genome wide chart
  output$gChart <- renderChart(create_gChart(1, input, values))
  output$gChart2 <- renderChart(create_gChart(2, input, values))

  output$zChart <- renderChart(create_zChart(1, input, values))
  output$zChart2 <- renderChart(create_zChart(2, input, values))

  #highcharts test chart
  output$testChart <- renderChart({
    h1 <- hPlot(x = "Wr.Hnd", y = "NW.Hnd", data = MASS::survey, type = c("line", 
                                                                          "bubble", "scatter"), group = "Clap", size = "Age")
    h1$set(dom = 'testChart')
    return(h1)     
  })
  
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
      # write.table(file="./www/config/datasetProperties.csv",x=currDatasetProp,col.names=TRUE,row.names=FALSE,sep=",")
      updateTabsetPanel(session, "datatabs", selected = "WhGen")
    })
#    if(input$selected != 1e5){
#      updateTabsetPanel(session, "datatabs", selected = "WhGen")  
#    }
  }
  observe(handleSubmitColsButton(1))
  observe(handleSubmitColsButton(2))

  # saveDataset <- function(j) {
  #   if (is.null(input[[jth_ref("saveDatasetButton", j)]]) || input[[jth_ref("saveDatasetButton", j)]] == 0) { return() }
  #   isolate({
  #     if(!file.exists(paste0("./www/config/data/",input[[jth_ref("datasets", j)]]))){
  #       write.table(getdata(j),paste0("./www/config/data/",input[[jth_ref("datasets", j)]]),sep=",",col.names=TRUE,row.names=FALSE)
  #     }
  #   })
  # }
  # observe(saveDataset(1))
  # observe(saveDataset(2))
  
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
    values$glGenesGlobalPlot <- NULL
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

isolate({
    # Parse neighboring genes from species 1
    results1 <- input$genomicLinkages$results1
    values$glSelectedGene <- results1$genes[[(length(results1$genes) + 1) %/% 2]]$name
    glGenes <- data.frame(matrix(unlist(results1$genes), nrow = length(results1$genes), byrow = TRUE),
      stringsAsFactors = FALSE)
    names(glGenes) <- names(results1$genes[[1]])
    glGenes <- glGenes[, c("family", "fmin", "fmax", "strand")]
    glGenes$chr <- trailingInteger(results1$chromosome_name)
    glGenes <- glGenes[nchar(glGenes$family) > 0, ]
    if (nrow(glGenes) == 0) {
      # could reach here if none of the (2*neighbors + 1) genes has a family id
      clearGenomicLinkages()
      return()
    }

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
      gr.chr <- trailingInteger(gr$chromosome_name)
      if (paste(substr(gr$genus,1,1),gr$species,sep=".") == org.G.species[values$organism2] && !is.na(gr.chr)) {
        gr.genes <- data.frame(matrix(unlist(gr$genes), nrow = length(gr$genes), byrow = TRUE),
          stringsAsFactors = FALSE)
        names(gr.genes) <- names(gr$genes[[1]])
        gr.genes <- gr.genes[, c("family", "fmin", "fmax", "strand")]
        gr.genes$chr <- gr.chr
        gr.genes$id <- gr$id
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
    fc <- getRainbowColors(nf)
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
        gr.chr <- trailingInteger(gr$chromosome_name)
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
  })

  observe({
    if (is.null(input$genesGlobalPlot)) return()
    if (!(input$boolGenomicLinkage && input$boolListenToBC)) return()
isolate({
    clearGenomicLinkages()

    # Parse the genes
    glGenes <- input$genesGlobalPlot$genes
    if (length(glGenes) == 0) return()
    nn <- names(glGenes[[1]]) # the column names
    glGenes <- data.frame(matrix(unlist(glGenes), nrow = length(glGenes), byrow = TRUE),
      stringsAsFactors = FALSE)
    names(glGenes) <- nn
    glGenes <- glGenes[, c("family", "fmin", "fmax", "strand")]
    glGenes$chr <- trailingInteger(input$genesGlobalPlot$chromosome)
    centerBP <- as.numeric(input$selected[[1]])
    winHigh <- centerBP + input$window[1]
    winLow <- centerBP - input$window[1]
    glGenes <- glGenes[as.integer(glGenes$fmax) >= winLow & as.integer(glGenes$fmin) <= winHigh, ]
    glGenes <- glGenes[order(as.integer(glGenes$fmin)), ] # sort as they may not be in order
    # assign colors
    ff <- unique(glGenes$family)
    nf <- length(ff)
    fc <- getRainbowColors(nf)
    familyColors <- list()
    lapply(1:nf, FUN = function(i) {
      if (i <= nf) familyColors[[ff[i]]] <<- stri_sub(fc[i], 1, 7)
    })
    glGenes$color <- familyColors[glGenes$family]

    if (nrow(glGenes) > 0) {
      values$glGenes <- glGenes
      values$glColors <- familyColors
    }
})
  })

  # Handle Broadcast Channel messages from the Genome Context Viewer
  # (move to a separate file?...)
  observeEvent(input$bc_gcv, {
    if (!(input$boolGenomicLinkage && input$boolListenToBC)) return()
    # TODO: move to servicesAPI.R ?
    if (input$bc_gcv$type == 'select') {
      # Parse the message
      flags = 0
# or handle each one separately...
# TODO: add a SNPs field
      if (!is.null(input$bc_gcv$targets$organism)) flags <- flags + 1
      if (!is.null(input$bc_gcv$targets$chromosome)) flags <- flags + 2
      if (!is.null(input$bc_gcv$targets$genes)) flags <- flags + 4
      if (!is.null(input$bc_gcv$targets$family)) flags <- flags + 8
      if (!is.null(input$bc_gcv$targets$extent)) flags <- flags + 16
      if (!is.null(input$bc_gcv$targets$block)) flags <- flags + 32
      isMacroSyntenyRow <- (flags == 3)
      isMacroSyntenyBlock <- (flags == 33)
      isMacroSyntenyOrganism <- (flags == 1)
      isMicroSyntenyRow <- (flags == 23)
      isMicroSyntenyGene <- isMicroSyntenyFamily <- FALSE
      if (flags == 12) {
        isMicroSyntenyGene <- (length(input$bc_gcv$targets$genes) == 1)
        isMicroSyntenyFamily <- !isMicroSyntenyGene
      }
      # TODO:
      # isDotPlot <- isMicroSyntenyGene # they have the same fields
      # Multi-view fields:
      # isMacroSyntenyCircosInnerBlock <- isMacroSyntenyRow # they have the same fields
      # isMacroSyntenyCircosOuterBlock <- (flags == 32)

      # Do something
      if (isMacroSyntenyRow) {
        j <- 0
        org <- input$bc_gcv$targets$organism
        # Note - the following tests only work for "Genus species" organism names
        if (org == values$organism) {
          j <- 1
        } else if (org == values$organism2) {
          j <- 2
        }
        if (j > 0) {
          # Extract the chromosome number
          chr <- trailingInteger(input$bc_gcv$targets$chromosome)
          # Adjust the Chromosome window to match the selection
          updateTabsetPanel(session, "datatabs", selected = "Chrom")
          updateSelectInput(session, jth_ref("chr", j), selected = chr)
        }

      } else if (isMacroSyntenyBlock) {
        org <- input$bc_gcv$targets$organism
        blk <- input$bc_gcv$targets$block
        # TODO: ...

      } else if (isMacroSyntenyOrganism) {
        org <- input$bc_gcv$targets$organism
        # TODO: ...

      } else if (isMicroSyntenyRow) {
        if (input$bc_gcv$targets$organism == values$organism) {
          # Extract the chromosome number
          chr <- trailingInteger(input$bc_gcv$targets$chromosome)
          # Range of base pairs
          bpMin <- input$bc_gcv$targets$extent[[1]]
          bpMax <- input$bc_gcv$targets$extent[[2]]
          centerBP <- (bpMax + bpMin) %/% 2
          widthBP <- (bpMax - bpMin) %/% 2 + 1000 # give it 1000 BPs on either side to ensure visibility
          # Adjust the Chromosome window to match the selection
          updateTabsetPanel(session, "datatabs", selected = "Chrom")
          updateSelectInput(session, "chr", selected = chr)
          updateNumericInput(session, "selected", value = centerBP)
          updateSliderInput(session, "window", value = widthBP)

          # TODO: Do something with the genes?
          # genes <- input$bc_gcv$targets$genes
        }

      } else if (isMicroSyntenyGene) {
        gene <- input$bc_gcv$targets$genes[[1]]
        fam <- input$bc_gcv$targets$family
# print(gene)
# print(fam)
# print("-----------------")
        # TODO: (possibilities)
        # Center the gene in its organism's window
        # Highlight the gene
        # Highlight other genes in the same family (in a different color?)

      } else if (isMicroSyntenyFamily) {
        # Check for singleton and orphan genes
        genes <- input$bc_gcv$targets$genes
        fam <- input$bc_gcv$targets$family
        isSingletons <- startsWith(fam, "singleton")
        isOrphans <- (fam == "")

        if (isSingletons) {
          # Highlight all singleton genes, for both organisms
          chr <- organismToChromosomeName(values$organism, as.integer(input$chr))
          # genes <- genes[-1] # ?
          # Parse "singleton,phytozome_10_2.xxxxxxxx,phytozome_10_2.yyyyyyyy,..."
          fam <- strsplit(fam, split = ",")[[1]][-1]
          fam <- paste0("'", fam, "'", collapse = ",")
          fam <- sprintf("[%s]", fam)
          query <- getGlobalPlot(chr, fam)
          runjs(query)

        } else if (isOrphans) {
          # Highlight all genes with no family, for both organisms
          chr <- organismToChromosomeName(values$organism, as.integer(input$chr))
          fam <- sprintf("['%s']", fam)
          # TODO: (can the GCV return genes with no family?)
          query <- getGlobalPlot(chr, fam)
          runjs(query)

        } else {
          # Highlight all genes in the selected family, for both organisms
          chr <- organismToChromosomeName(values$organism, as.integer(input$chr))
          fam <- sprintf("['%s']", fam)
          query <- getGlobalPlot(chr, fam)
          runjs(query)
          # TODO: organism 2
        }
      }
    }
  })
  # The messages GCV currently implements conform to the following schema:
    # {
    #   type: "select" | "deselect",
    #   targets: {
    #     organism?: String,  // an organism identifier of the form "<genus> <species>"
    #     chromosome?: String,  // a chromosome identifier
    #     genes?: String[],  // an array of gene identifiers
    #     family?: String,  // a gene family identifier
    #     extent?: Integer[2],  // a length 2 array representing a genomic interval
    #     block?: {  // an object representing a pairwise synteny block
    #       source: {
    #         chromosome: String,
    #         locus: Integer[2]
    #       },
    #       reference: {
    #         chromosome: String,
    #         locus: Integer[2]
    #       },
    #       orientation: "+" | "-"
    #     }
    #   }
    # }

  # Turn off (deselect) any old highlighted targets
  clearGCVHighlights <- function() {
    if (!is.null(values$gcvTargets)) {
      robj <- list(type = "deselect", targets = values$gcvTargets)
      msg <- sprintf("bc.postMessage(JSON.parse('%s'));", toJSON(robj, auto_unbox = TRUE))
      runjs(msg)
    }
  }

  # View the current genomic linkage query in the Genome Context Viewer
  # (this does not involve Broadcast Channel)
  observeEvent(input$viewInGCV, isolate({
    if (!is.null(values$glSelectedGene)) {
      gcvQuery <- sprintf("window.open('/gcv/legfed_v1_0/search/lis/%s?neighbors=%d&matched=%d&intermediate=%d&regexp=%s', 'gcv');",
        values$glSelectedGene, input$neighbors, input$matched, input$intermediate, tolower(org.Gensp[values$organism2]))
      runjs(gcvQuery)
    }
  }))

  # Post a Broadcast Channel message to the Genome Context Viewer
  # (which handles the same kind of messages it sends)
  observe({
    if (is.null(input$gcvGeneFamily)) return()
    isolate(clearGCVHighlights())

    # User selected a family (in the Chromosome View) to highlight in the Genome Context Viewer
    gcvTargets = list(family = input$gcvGeneFamily)
    robj <- list(type = "select", targets = gcvTargets)
    msg <- sprintf("bc.postMessage(JSON.parse('%s'));", toJSON(robj, auto_unbox = TRUE))
    values$gcvTargets <- gcvTargets
    runjs(msg)

    # Then reset it to enable selecting the same family in the future
    runjs("Shiny.onInputChange('gcvGeneFamily', null);")
  })

  # Tour based on rintrojs
  df.tour <- read.csv(file = "tour.csv", header = TRUE)
  observeEvent(input$zbrowseTour, {
    introjs(session,
      options = list(steps = df.tour, showBullets = FALSE, showStepNumbers = FALSE, skipLabel = "End Tour"),
      # allow switching tabs
      events = list(onbeforechange = readCallback("switchTabs"))
    )
  })

})#end server
