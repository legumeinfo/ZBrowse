shinyServer(function(input, output, session) {

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
  
  # This should be the first code block to detect a change in input$datasets
  observe({
    if (is.null(input$datasets)) return()

    # Make sure the organism corresponds to the selected dataset
    # (invoke values$organism first to force a single reaction)
    if (is.null(values$organism)) {
      values$organism = "Corn"
    } else {
      isolate({
        values$organism <- values$datasetToOrganism[[input$datasets]]
      })
    }

    # Uncheck the Append to Current Dataset checkbox and clear any previously selected GWAS files
    updateSelectizeInput(session, "gwasTraits", choices = gwas.traits[[values$organism]], selected = NULL)
    # As there is no updateFileInput() method,
    # send a custom message (to clear the progress bar) and clear values$needsToUploadFiles
    session$sendCustomMessage(type = "resetFileInputHandler", "uploadfile")
    values$needsToUploadFiles <- FALSE
    updateCheckboxInput(session, "appendSNPs", value = FALSE)
    # Clear all genomic linkages if either organism changes
    values$glSelectedGene <- NULL
    values$glGenes1 <- values$glGenes2 <- NULL
  })
  # This should be the first code block to detect a change in input$datasets2
  observe({
    if (is.null(input$datasets2)) return()
    
    # Make sure the organism corresponds to the selected dataset
    # (invoke values$organism2 first to force a single reaction)
    if (is.null(values$organism2)) {
      values$organism2 = "Corn"
    } else {
      isolate({
        values$organism2 <- values$datasetToOrganism[[input$datasets2]]
      })
    }
    
    # Uncheck the Append to Current Dataset checkbox and clear any previously selected GWAS files
    updateSelectizeInput(session, "gwasTraits2", choices = gwas.traits[[values$organism2]], selected = NULL)
    # As there is no updateFileInput() method,
    # send a custom message (to clear the progress bar) and clear values$needsToUploadFiles2
    session$sendCustomMessage(type = "resetFileInputHandler", "uploadfile2")
    values$needsToUploadFiles2 <- FALSE
    updateCheckboxInput(session, "appendSNPs2", value = FALSE)
    # Clear all genomic linkages if either organism changes
    values$glSelectedGene <- NULL
    values$glGenes1 <- values$glGenes2 <- NULL
  })

  # The user selected one or more local GWAS files (not yet loaded)
  observe({
    input$uploadfile
    values$needsToUploadFiles <- TRUE
  })
  observe({
    input$uploadfile2
    values$needsToUploadFiles2 <- TRUE
  })

  #handles what displays in the sidebar based on what tab is selected
  output$ui_All <- renderUI({
    list(
      conditionalPanel(condition = "input.datatabs == 'Manage'",
         wellPanel(
           style = paste0("background-color: ", bgColors[1], ";"),
           uiOutput("datasets")
         ),                       
         wellPanel(
           style = paste0("background-color: ", bgColors[1], ";"),
           radioButtons(inputId = "dataType", label = "Load data (Max. 5MB):", c(".csv" = "csv", ".rda" = "rda", "examples" = "examples"), selected = "csv"),
           conditionalPanel(condition = "input.dataType != 'examples'",
                            conditionalPanel(condition = "input.dataType == 'csv'",
                                             checkboxInput('header', 'Header', TRUE),
                                             radioButtons('sep', '', c(Comma=',', Semicolon=';', Tab='\t'), ',')                                             
                            ),
                            selectizeInput("gwasTraits", "Remote Trait Files:", choices = NULL, multiple = TRUE),
                            fileInput("uploadfile", "Local Trait Files:", multiple = TRUE),
                            checkboxInput("appendSNPs", "Append to current dataset", FALSE),
                            actionButton("loadTraits", "Load Data"),
tags$script("Shiny.addCustomMessageHandler('resetFileInputHandler', function(x) {      
  var id = '#' + x;
  var idProgress = id + '_progress';
  var idBar = id + ' .bar';
  $(idProgress).css('visibility', 'hidden');
  $(idBar).css('width', '0%');
  // TODO: clear the loaded file(s) (none of the following work, for Javascript security reasons)
  //$(id).val('');
  //$(id).replaceWith($(id).val('').clone(true));
  //$(id).replaceWith($(id) = $(id).clone(true));
  //$(id).wrap('<form></form>').closest('form').reset();
    //$(id).unwrap();
    //$(id).stopPropagation();
    //$(id).preventDefault();
});")
           ),      
           conditionalPanel(condition = "input.dataType == 'examples'",
                            actionButton('loadExampleData', 'Load examples')
           )
         ),
         wellPanel(
           style = paste0("background-color: ", bgColors[1], ";"),
            h6("Once your file is finished uploading, press the Save Dataset button below and reload Zbrowse."),
           actionButton('saveDatasetButton', 'Save Current Dataset'),
           conditionalPanel(condition = "input.saveDatasetButton > 0",
                            h5("Dataset successfully saved!")
           )
         ),
        wellPanel(
          style = paste0("background-color: ", bgColors[2], ";"),
          uiOutput("datasets2")
        ),
        wellPanel(
          style = paste0("background-color: ", bgColors[2], ";"),
          radioButtons(inputId = "dataType2", label = "Load data (Max. 5MB):", c(".csv" = "csv", ".rda" = "rda", "examples" = "examples"), selected = "csv"),
          conditionalPanel(condition = "input.dataType2 != 'examples'",
            conditionalPanel(condition = "input.dataType2 == 'csv'",
              checkboxInput('header2', 'Header', TRUE),
              radioButtons('sep2', '', c(Comma=',', Semicolon=';', Tab='\t'), ',')
            ),
            selectizeInput("gwasTraits2", "Remote Trait Files:", choices = NULL, multiple = TRUE),
            fileInput("uploadfile2", "Local Trait Files:", multiple = TRUE),
            checkboxInput("appendSNPs2", "Append to current dataset", FALSE),
            actionButton("loadTraits2", "Load Data")
          ),
          conditionalPanel(condition = "input.dataType2 == 'examples'",
            actionButton('loadExampleData2', 'Load examples')
          )
        ),
        wellPanel(
          style = paste0("background-color: ", bgColors[2], ";"),
          h6("Once your file is finished uploading, press the Save Dataset button below and reload Zbrowse."),
          actionButton('saveDatasetButton2', 'Save Current Dataset'),
          conditionalPanel(condition = "input.saveDatasetButton2 > 0",
            h5("Dataset successfully saved!")
          )
        ),
        helpModal('Manage','manage',includeMarkdown("tools/manage.md")),HTML('<p style="font-size:10px;">Powered by <a href="http://www.rstudio.com/shiny/">Shiny</a>, <a href="http://rcharts.io/">rCharts</a> and <a href="http://www.highcharts.com">Highcharts</a></p>')             
      ),#end conditional Manage
      conditionalPanel(condition = "input.datatabs == 'Table'",
         wellPanel(
           style = paste0("background-color: ", bgColors[1], ";"),
           uiOutput("columns"),
           tags$br(),tags$br(),tags$br(),tags$br(),tags$br(),tags$br(), #add some space between selection columns and subset search
           # uiOutput("view_order"), checkboxInput("view_order_desc", "DESC", value = FALSE),
           returnTextInput("dv_select", "Subset (e.g., RMIP > 20 & Location == 'FL06')", '')
         ),
        wellPanel(
          style = paste0("background-color: ", bgColors[2], ";"),
          uiOutput("columns2"),
          tags$br(),tags$br(),tags$br(),tags$br(),tags$br(),tags$br(), #add some space between selection columns and subset search
          # uiOutput("view_order"), checkboxInput("view_order_desc", "DESC", value = FALSE),
          returnTextInput("dv_select2", "Subset (e.g., RMIP > 20 & Location == 'FL06')", '')
        ),
        helpModal('Data Table View','view',includeMarkdown("tools/manage.md"))      
      ),#end conditional Table
      #the second part of the statement is what is allowing the detection of changing panels due to a click event, i couldn't figure out how to update input$datatabs with javascript
      conditionalPanel(condition = "input.datatabs == 'panel1' || input.datatabs == 'panel2' || $('li.active a').first().html()==='Chromosome View'",
         helpText(h5(p("Interactive Graphs for GWAS Data"))),
         wellPanel(
            style = paste0("background-color: ", bgColors[1], ";"),
            uiOutput("traitColBoxes"),
            uiOutput("legend"),
            uiOutput("overlaps"),
            conditionalPanel(condition = "input.overlaps==true",
                             uiOutput("overlapSize"),
                             uiOutput("numOverlapTraits")
                             )#end conditional for plotting only overlaps
         )
         #submitButton("Update View"),

      ),#end conditional for traitCols (genome view and chromosomeview)
      conditionalPanel(condition = "input.datatabs == 'panel2' || $('li.active a').first().html()==='Chromosome View'",
         wellPanel(
            style = paste0("background-color: ", bgColors[1], ";"),
            uiOutput("selectChr")
         )
      ),#end conitional for chromsomeview (panel2)
      conditionalPanel(condition = "input.datatabs == 'Annot' || input.datatabs == 'panel2' || $('li.active a').first().html()==='Chromosome View'",
         wellPanel(
            style = paste0("background-color: ", bgColors[1], ";"),
            h5("Annotation window options:"),
            h6("Click a point or type a basepair value:"),
            uiOutput("selectedOut"),
            uiOutput("windowOut")                        
         )
      ),#end conditional panel for Annotation plot and Table
      conditionalPanel(condition = "input.datatabs == 'Annot'",
         style = paste0("background-color: ", bgColors[1], ";"),
         helpText(h5(p(paste("Download a CSV of the annotations in the selected window.")))),
         wellPanel(                        
           style = paste0("background-color: ", bgColors[1], ";"),
           downloadButton('downloadAnnot','Download')
         )
      ),#end conditional panel for Annotation Table
      # Genomic Linkages panel
      conditionalPanel(condition = "input.datatabs == 'panel2' || $('li.active a').first().html()==='Chromosome View'",
        wellPanel(
          h5("Genomic Linkage options:"),
          checkboxInput('boolGenomicLinkage', 'ON', FALSE),
          conditionalPanel("input.boolGenomicLinkage == true",
            wellPanel(
              uiOutput("selectedGene"),
              numericInput("neighbors", "Neighbors:", min = 1, max = 20, value=10),
              numericInput("matched", "Matched:", min = 1, max = 20, value=6),
              numericInput("intermediate", "Intermediate:", min = 1, max = 10, value=3),
              style = paste0("background-color: ", bgColors[1], ";")
            ),
            wellPanel(
              uiOutput("relatedRegions"),
              style = paste0("background-color: ", bgColors[2], ";")
            )
          )
        )
      ),
conditionalPanel(condition = "input.datatabs == 'panel1' || input.datatabs == 'panel2' || $('li.active a').first().html()==='Chromosome View'",
  wellPanel(
    style = paste0("background-color: ", bgColors[2], ";"),
    uiOutput("traitColBoxes2"),
    uiOutput("legend2"),
    uiOutput("overlaps2"),
    conditionalPanel(condition = "input.overlaps2==true",
      uiOutput("overlapSize2"),
      uiOutput("numOverlapTraits2")
    )#end conditional for plotting only overlaps
  )
),#end conditional for traitCols (genome view and chromosomeview)
conditionalPanel(condition = "input.datatabs == 'panel2' || $('li.active a').first().html()==='Chromosome View'",
  wellPanel(
    style = paste0("background-color: ", bgColors[2], ";"),
    uiOutput("selectChr2")
  )
),#end conitional for chromsomeview (panel2)
conditionalPanel(condition = "input.datatabs == 'Annot' || input.datatabs == 'panel2' || $('li.active a').first().html()==='Chromosome View'",
  wellPanel(
    style = paste0("background-color: ", bgColors[2], ";"),
    h5("Annotation window options:"),
    h6("Click a point or type a basepair value:"),
    uiOutput("selectedOut2"),
    uiOutput("windowOut2")                        
  )
),#end conditional panel for Annotation plot and Table
conditionalPanel(condition = "input.datatabs == 'Annot'",
  style = paste0("background-color: ", bgColors[2], ";"),
  helpText(h5(p(paste("Download a CSV of the annotations in the selected window.")))),
  wellPanel(                        
    style = paste0("background-color: ", bgColors[2], ";"),
    downloadButton('downloadAnnot2','Download')
  )
),#end conditional panel for Annotation Table 2
      conditionalPanel(condition = "input.datatabs == 'panel1' || input.datatabs == 'panel2' || $('li.active a').first().html()==='Chromosome View'",      
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

  output$datasets <- renderUI({   
    if (is.null(input$loadTraits)) return()

    # Load any requested GWAS files
    isolate({
      values$uploadFiles <- input$uploadfile
      values$gwasTraits <- input$gwasTraits

      # Local GWAS files
      inFile <- NULL
      if (values$needsToUploadFiles) inFile <- values$uploadFiles
      if(!is.null(inFile)) {
        # iterating through the files to upload
        for(i in 1:(dim(inFile)[1])) {
          loadUserData(inFile[i,'name'], inFile[i,'datapath'])
          # unlink(inFile[i,'datapath'], recursive = FALSE, force = TRUE)
        }
        values$needsToUploadFiles <- FALSE # since we just loaded them
      }
      # Remote GWAS files
      inTraits <- values$gwasTraits
      if (!is.null(inTraits)) {
        for (trait in inTraits) {
          trait.url <- gwas.filenames[[values$organism]][which(gwas.traits[[values$organism]] == trait)]
          loadRemoteData(trait, trait.url)
        }
      }
    })

    dat <- isolate(input$datasets)
    if (!is.null(dat)) {
      appendSNPs <- isolate(input$appendSNPs)
      if (appendSNPs || (is.null(inFile) && is.null(inTraits))) {
        val <- dat
      } else {
        val <- values$datasetlist[1]
      }
    }else{
      val <- "Medicago truncatula GWAS"
    }
        
    # Drop-down selection of data set
    # selectInput(inputId = "datasets", label = "Datasets:", choices = datasets, selected = datasets[1], multiple = FALSE)
    selectInput(inputId = "datasets", label = "Dataset 1:", choices = values$datasetlist, selected = values$datasetlist[values$datasetlist==val], multiple = FALSE, selectize=FALSE)
  })
  output$datasets2 <- renderUI({   
    if (is.null(input$loadTraits2)) return()

    # Load any requested GWAS files
    isolate({
      values$uploadFiles2 <- input$uploadfile2
      values$gwasTraits2 <- input$gwasTraits2

      # Local GWAS files
      inFile <- NULL
      if (values$needsToUploadFiles2) inFile <- values$uploadFiles2
      if(!is.null(inFile)) {
        # iterating through the files to upload
        for(i in 1:(dim(inFile)[1])) {
          loadUserData2(inFile[i,'name'], inFile[i,'datapath'])
        }
        values$needsToUploadFiles2 <- FALSE # since we just loaded them
      }
      # Remote GWAS files
      inTraits <- values$gwasTraits2
      if (!is.null(inTraits)) {
        for (trait in inTraits) {
          trait.url <- gwas.filenames[[values$organism2]][which(gwas.traits[[values$organism2]] == trait)]
          loadRemoteData2(trait, trait.url)
        }
      }
    })

    dat <- isolate(input$datasets2)
    if (!is.null(dat)) {
      appendSNPs <- isolate(input$appendSNPs2)
      if (appendSNPs || (is.null(inFile) && is.null(inTraits))) {
        val <- dat
      } else {
        val <- values$datasetlist[1]
      }
    }else{
      val <- "Arabidopsis thaliana GWAS"
    }

    # Drop-down selection of data set
    selectInput(inputId = "datasets2", label = "Dataset 2:", choices = values$datasetlist, selected = values$datasetlist[values$datasetlist==val], multiple = FALSE, selectize=FALSE)
  })
  
  reactiveAnnotTable <- reactive({
    if(is.null(input$datasets)) return()
    centerBP <- as.numeric(input$selected[[1]])
    winHigh <- centerBP+input$window[1]
    winLow <- centerBP-input$window[1]
    if(winLow < 0){winLow <- 0}    
    thisChrAnnot <- subset(annotGeneLoc[values$organism][[1]],chromosome==input$chr)
    thisAnnot <- thisChrAnnot[thisChrAnnot$transcript_start >= winLow & thisChrAnnot$transcript_end <= winHigh,]    
    thisAnnot  
  })
  reactiveAnnotTable2 <- reactive({
    if(is.null(input$datasets2)) return()
    centerBP <- as.numeric(input$selected2[[1]])
    winHigh <- centerBP+input$window2[1]
    winLow <- centerBP-input$window2[1]
    if(winLow < 0){winLow <- 0}    
    thisChrAnnot <- subset(annotGeneLoc[values$organism2][[1]],chromosome==input$chr2)
    thisAnnot <- thisChrAnnot[thisChrAnnot$transcript_start >= winLow & thisChrAnnot$transcript_end <= winHigh,]    
    thisAnnot  
  })
  
  #Returns the nicely formatted preview table
  output$htmlDataExample <- reactive({
    if(is.null(input$datasets)) return()
    
    dat <- getdata()

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
    
  })
  output$htmlDataExample2 <- reactive({
    if(is.null(input$datasets2)) return()
    
    dat <- getdata2()
    
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
  })

  output$nrowDataset <- reactive({
    if(is.null(input$datasets)) return()
    dat <- getdata()
    if(is.null(dat) || nrow(dat) == 0) return()
    nr <- nrow(dat)
    
    if(nr>2500){
      paste0('<p>First 10 rows shown of ',nr,' rows. See Data Table tab for details.<br>More than 2500 rows found, only the top 2500 will be plotted.</p>')
    }else{
      paste0('<p>First 10 rows shown of ',nr,' rows. See Data Table tab for details.</p>')
    }
  })
  output$nrowDataset2 <- reactive({
    if(is.null(input$datasets2)) return()
    dat <- getdata2()
    if(is.null(dat) || nrow(dat) == 0) return()
    nr <- nrow(dat)
    
    if(nr>2500){
      paste0('<p>First 10 rows shown of ',nr,' rows. See Data Table tab for details.<br>More than 2500 rows found, only the top 2500 will be plotted.</p>')
    }else{
      paste0('<p>First 10 rows shown of ',nr,' rows. See Data Table tab for details.</p>')
    }
  })

  output$ui_data_tabs <- renderUI({
    tabsetPanel(id = "datatabs",      
      tabPanel(title="Manage",value="Manage",
               tags$div(
                class = "container",
                style = paste0("background-color: ", bgColors[1], ";"),
                
                htmlOutput("htmlDataExample"),
                htmlOutput("nrowDataset"),
                
                #tags$p(tags$br()),
                row(
                  col(3, tags$br()),
                  col(7, h4('Select appropriate columns to be used for plotting.'))
                #HTML('<h4>Select appropriate columns to be used for plotting.</h4>'),
                ),
                tags$hr(),
                row(
                   #col(2, tags$br()),
                   col(2,uiOutput("chrColumn"),uiOutput("bpColumn")),
                   col(2,uiOutput("plotAll"),uiOutput("traitColumns")),
                   col(2,uiOutput("yAxisColumn"),uiOutput("logP")),
                   col(2,uiOutput("axisLimBool"),uiOutput("axisLim")),
                   col(2,actionButton("SubmitColsButton","Submit"))
                ),
                tags$hr(),
                row(
                  col(7, uiOutput("supportInterval"))#
                ),
                row(
                  col(2, uiOutput("SIbpStart")),
                  col(2, uiOutput("SIyAxisColumn")),
                  col(2, uiOutput("SIaxisLimBool"),uiOutput("SIaxisLim"))
                )
              ),
#                   col(
#                     4,
#                     uiOutput("traitColumns"),
#                     uiOutput("yAxisColumn")
#                   )
#                 )                                                
# #               #HTML(dataDescriptionOutput())
#                )
tags$div(
  class = "container",
  style = paste0("background-color: ", bgColors[2], ";"),

  htmlOutput("htmlDataExample2"),
  htmlOutput("nrowDataset2"),
  
  #tags$p(tags$br()),
  row(
    col(3, tags$br()),
    col(7, h4('Select appropriate columns to be used for plotting.'))
    #HTML('<h4>Select appropriate columns to be used for plotting.</h4>'),
  ),
  tags$hr(),
  row(
    #col(2, tags$br()),
    col(2,uiOutput("chrColumn2"),uiOutput("bpColumn2")),
    col(2,uiOutput("plotAll2"),uiOutput("traitColumns2")),
    col(2,uiOutput("yAxisColumn2"),uiOutput("logP2")),
    col(2,uiOutput("axisLimBool2"),uiOutput("axisLim2")),
    col(2,actionButton("SubmitColsButton2","Submit"))
  ),
  tags$hr(),
  row(
    col(7, uiOutput("supportInterval2"))#
  ),
  row(
    col(2, uiOutput("SIbpStart2")),
    col(2, uiOutput("SIyAxisColumn2")),
    col(2, uiOutput("SIaxisLimBool2"),uiOutput("SIaxisLim2"))
  )
)
      ),
      tabPanel(title="Data Table",value="Table",
        wellPanel(dataTableOutput("dataviewer"), style = paste0("background-color: ", bgColors[1], ";")),
        wellPanel(dataTableOutput("dataviewer2"), style = paste0("background-color: ", bgColors[2], ";"))
      ),
      tabPanel(title="Whole Genome View",value="panel1",
        wellPanel(showOutput("gChart", "highcharts"), style = paste0("background-color: ", bgColors[1], ";")),
        wellPanel(showOutput("gChart2", "highcharts"), style = paste0("background-color: ", bgColors[2], ";"))
      ),
      tabPanel(title="Chromosome View",value="panel2",
        wellPanel(showOutput("pChart", "highcharts"), showOutput("zChart", "highcharts"),
          tags$script('Shiny.addCustomMessageHandler("customMsg", function(bandOpts){
            chartXAxis = $("#pChart").highcharts().xAxis[0]
            chartXAxis.removePlotBand()
            chartXAxis.addPlotBand(bandOpts)
          })'),
          style = paste0("background-color: ", bgColors[1], ";")
        ),
        wellPanel(showOutput("pChart2", "highcharts"), showOutput("zChart2", "highcharts"),
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
  output$annotViewer <- renderDataTable({
#    if(is.null(input$datasets)) return()
#    centerBP <- as.numeric(input$selected[[1]])
#    winHigh <- centerBP+input$window[1]
#    winLow <- centerBP-input$window[1]
#    if(winLow < 0){winLow <- 0}    
#    thisChrAnnot <- subset(annotGeneLoc[input$organism][[1]],chromosome==input$chr)
#    thisAnnot <- thisChrAnnot[thisChrAnnot$transcript_start >= winLow & thisChrAnnot$transcript_end <= winHigh,]        
#    thisAnnot
    reactiveAnnotTable()
  }, options = list(orderClasses = TRUE, bCaseInsensitive = TRUE,
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
  )#end annotation table
output$annotViewer2 <- renderDataTable({
  reactiveAnnotTable2()
}, options = list(orderClasses = TRUE, bCaseInsensitive = TRUE,
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
)#end annotation table 2

  output$dataviewer <-renderDataTable({    
    if(is.null(input$datasets) || is.null(input$columns)) return()
    
    dat <- getdata()
    #dat <- date2character()
    
    if(!all(input$columns %in% colnames(dat))) return()
    
    if(input$dv_select != '') {
      selcom <- input$dv_select
      selcom <- gsub(" ", "", selcom)
      
      seldat <- try(do.call(subset, list(dat,parse(text = selcom))), silent = TRUE)
      
      if(!is(seldat, 'try-error')) {
        if(is.data.frame(seldat)) {
          dat <- seldat
          seldat <- NULL
        }
      }
    }
    
    dat <- data.frame(dat[, input$columns, drop = FALSE])
    dat
    
    # html <- print(xtable::xtable(dat), type='html', print.results = FALSE)
    # html <- sub("<TABLE border=1>","<table class='table table-condensed table-hover'>", html)
    # html
    
  }, options = list(orderClasses = TRUE, bCaseInsensitive = TRUE,
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
  )#end dataviewer
output$dataviewer2 <-renderDataTable({    
  if(is.null(input$datasets2) || is.null(input$columns2)) return()
  
  dat <- getdata2()
  #dat <- date2character()
  
  if(!all(input$columns2 %in% colnames(dat))) return()
  
  if(input$dv_select2 != '') {
    selcom <- input$dv_select2
    selcom <- gsub(" ", "", selcom)
    
    seldat <- try(do.call(subset, list(dat,parse(text = selcom))), silent = TRUE)
    
    if(!is(seldat, 'try-error')) {
      if(is.data.frame(seldat)) {
        dat <- seldat
        seldat <- NULL
      }
    }
  }
  
  dat <- data.frame(dat[, input$columns2, drop = FALSE])
  dat
  
  # html <- print(xtable::xtable(dat), type='html', print.results = FALSE)
  # html <- sub("<TABLE border=1>","<table class='table table-condensed table-hover'>", html)
  # html
  
}, options = list(orderClasses = TRUE, bCaseInsensitive = TRUE,
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
)#end dataviewer2
  
  output$downloadAnnot <- downloadHandler(
    filename = function() {paste0("AnnotationsAround.chr",input$chr,".",input$selected[[1]],"bp.",values$organism,".csv")},
    content = function(file) {write.csv(reactiveAnnotTable(),file,row.names=F)}
  )
  output$downloadAnnot2 <- downloadHandler(
    filename = function() {paste0("AnnotationsAround.chr",input$chr2,".",input$selected2[[1]],"bp.",values$organism2,".csv")},
    content = function(file) {write.csv(reactiveAnnotTable2(),file,row.names=F)}
  )
  
  output$columns <- renderUI({
    cols <- varnames()    
    selectInput("columns", "Select columns to show:", choices = as.list(cols), selected = cols, multiple = TRUE)
  })
  output$columns2 <- renderUI({
    cols <- varnames2()    
    selectInput("columns2", "Select columns to show:", choices = as.list(cols), selected = cols, multiple = TRUE)
  })
  
  output$axisLimBool <- renderUI({
    if(is.null(input$datasets)){return()}
    if(input$datasets %in% datasetProp()$dataset){
      val = datasetProp()$axisLim[datasetProp()$dataset == input$datasets]
    } else if (values$organism %in% legumeInfo.organisms) {
      val <- FALSE
    }else{
      val = TRUE}
    checkboxInput('axisLimBool', 'Set Y-axis Limits?', val)
  })
  output$axisLimBool2 <- renderUI({
    if(is.null(input$datasets2)){return()}
    if(input$datasets2 %in% datasetProp()$dataset){
      val = datasetProp()$axisLim[datasetProp()$dataset == input$datasets2]
    } else if (values$organism2 %in% legumeInfo.organisms) {
      val <- FALSE
    }else{
      val = TRUE}
    checkboxInput('axisLimBool2', 'Set Y-axis Limits?', val)
  })
  
  output$logP <- renderUI({
    if(is.null(input$datasets)){return()}
    if(input$datasets %in% datasetProp()$dataset){
      val = datasetProp()$logP[datasetProp()$dataset == input$datasets]
    } else if (values$organism %in% legumeInfo.organisms) {
      val <- TRUE
    }else{
      val = FALSE}
    checkboxInput('logP', 'Take -log10 of column?', val)
  })
  output$logP2 <- renderUI({
    if(is.null(input$datasets2)){return()}
    if(input$datasets2 %in% datasetProp()$dataset){
      val = datasetProp()$logP[datasetProp()$dataset == input$datasets2]
    } else if (values$organism2 %in% legumeInfo.organisms) {
      val <- TRUE
    }else{
      val = FALSE}
    checkboxInput('logP2', 'Take -log10 of column?', val)
  })

  output$chrColumn <- renderUI({
    if(is.null(input$datasets)){return()}
    cols <- varnames()    
    if(input$datasets %in% datasetProp()$dataset){
      selected = datasetProp()$chrColumn[datasetProp()$dataset == input$datasets]
    }else{
      selected = names(cols[1])
    }    
    selected <- stri_match(selected, regex = ".*(?=\\ \\{)")[, 1]
    selectizeInput("chrColumn", "Chromosome Column:", choices = as.list(cols), selected = selected, multiple = FALSE, options = list(dropdownParent="body"))
  })
  output$chrColumn2 <- renderUI({
    if(is.null(input$datasets2)){return()}
    cols <- varnames2()    
    if(input$datasets2 %in% datasetProp()$dataset){
      selected = datasetProp()$chrColumn[datasetProp()$dataset == input$datasets2]
    }else{
      selected = names(cols[1])
    }    
    selected <- stri_match(selected, regex = ".*(?=\\ \\{)")[, 1]
    selectizeInput("chrColumn2", "Chromosome Column:", choices = as.list(cols), selected = selected, multiple = FALSE, options = list(dropdownParent="body"))
  })

  output$bpColumn <- renderUI({
    if(is.null(input$datasets)){return()}
    cols <- varnames()
    if(input$datasets %in% datasetProp()$dataset){
      selected = datasetProp()$bpColumn[datasetProp()$dataset == input$datasets]
    }else{
      selected = names(cols[2])
    }
    selected <- stri_match(selected, regex = ".*(?=\\ \\{)")[, 1]
    selectizeInput("bpColumn", "Base Pair Column:", choices = as.list(cols), selected = selected, multiple = FALSE, options = list(dropdownParent="body"))
  })
  output$bpColumn2 <- renderUI({
    if(is.null(input$datasets2)){return()}
    cols <- varnames2()
    if(input$datasets2 %in% datasetProp()$dataset){
      selected = datasetProp()$bpColumn[datasetProp()$dataset == input$datasets2]
    }else{
      selected = names(cols[2])
    }
    selected <- stri_match(selected, regex = ".*(?=\\ \\{)")[, 1]
    selectizeInput("bpColumn2", "Base Pair Column:", choices = as.list(cols), selected = selected, multiple = FALSE, options = list(dropdownParent="body"))
  })
  
  output$traitColumns <- renderUI({
    if(is.null(input$datasets)){return()}
    cols <- varnames()
    if(input$datasets %in% datasetProp()$dataset){
      selected = unlist(strsplit(datasetProp()$traitCol[datasetProp()$dataset == input$datasets],";"))
    } else if (values$organism %in% legumeInfo.organisms) {
      selected <- names(cols[3])
    }else{
      selected = names(cols[3:4])
    }
    selected <- sapply(selected, FUN = function(s) stri_match(selected, regex = ".*(?=\\ \\{)")[, 1], USE.NAMES = FALSE)
    conditionalPanel(condition = "input.plotAll==false",
                     selectizeInput("traitColumns", "Group by these trait column(s):", choices = as.list(cols), selected = selected, multiple = TRUE, options = list(dropdownParent="body"))
    )        
  })  
  output$traitColumns2 <- renderUI({
    if(is.null(input$datasets2)){return()}
    cols <- varnames2()
    if(input$datasets2 %in% datasetProp()$dataset){
      selected = unlist(strsplit(datasetProp()$traitCol[datasetProp()$dataset == input$datasets2],";"))
    } else if (values$organism2 %in% legumeInfo.organisms) {
      selected <- names(cols[3])
    }else{
      selected = names(cols[3:4])
    }
    selected <- sapply(selected, FUN = function(s) stri_match(selected, regex = ".*(?=\\ \\{)")[, 1], USE.NAMES = FALSE)
    conditionalPanel(condition = "input.plotAll2==false",
      selectizeInput("traitColumns2", "Group by these trait column(s):", choices = as.list(cols), selected = selected, multiple = TRUE, options = list(dropdownParent="body"))
    )        
  })
  
  output$plotAll <- renderUI({
    if(is.null(input$datasets)){return()}
    if(input$datasets %in% datasetProp()$dataset){      
      val = datasetProp()$plotAll[datasetProp()$dataset == input$datasets]
    }else{
      val = FALSE
    }    
    checkboxInput('plotAll', 'All data is the same trait', val)    
  })
  output$plotAll2 <- renderUI({
    if(is.null(input$datasets2)){return()}
    if(input$datasets2 %in% datasetProp()$dataset){      
      val = datasetProp()$plotAll[datasetProp()$dataset == input$datasets2]
    }else{
      val = FALSE
    }    
    checkboxInput('plotAll2', 'All data is the same trait', val)    
  })
  
  output$yAxisColumn <- renderUI({
    if(is.null(input$datasets)){return()}
    cols <- varnames()       
    if(input$datasets %in% datasetProp()$dataset){      
      selected = datasetProp()$yAxisColumn[datasetProp()$dataset == input$datasets]
    } else if (values$organism %in% legumeInfo.organisms) {
      selected <- names(cols[4])
    }else{
      #selected = names(cols[10])
      selected = as.character(cols[10])
    }    
    selected <- stri_match(selected, regex = ".*(?=\\ \\{)")[, 1]
    selectizeInput("yAxisColumn", "Y-axis column:", choices = as.list(cols), selected = selected, multiple = FALSE, options = list(dropdownParent="body"))
  })
  outputOptions(output, "yAxisColumn", suspendWhenHidden=FALSE)
  output$yAxisColumn2 <- renderUI({
    if(is.null(input$datasets2)){return()}
    cols <- varnames2()       
    if(input$datasets2 %in% datasetProp()$dataset){      
      selected = datasetProp()$yAxisColumn[datasetProp()$dataset == input$datasets2]
    } else if (values$organism2 %in% legumeInfo.organisms) {
      selected <- names(cols[4])
    }else{
      #selected = names(cols[10])
      selected = as.character(cols[10])
    }    
    selected <- stri_match(selected, regex = ".*(?=\\ \\{)")[, 1]
    selectizeInput("yAxisColumn2", "Y-axis column:", choices = as.list(cols), selected = selected, multiple = FALSE, options = list(dropdownParent="body"))
  })
  outputOptions(output, "yAxisColumn2", suspendWhenHidden=FALSE)
  
  output$axisLim <- renderUI({    
    if(is.null(input$datasets)){return()}
    if(input$datasets %in% datasetProp()$dataset){
      min = datasetProp()$axisMin[datasetProp()$dataset == input$datasets]
      max = datasetProp()$axisMax[datasetProp()$dataset == input$datasets]
    }else{
      min = 0
      max = 1
    }
    conditionalPanel(condition = "input.axisLimBool==true",
                     numericInput("axisMin","Min:",value=min),
                     numericInput("axisMax","Max:",value=max)
    )    
  })  
  output$axisLim2 <- renderUI({    
    if(is.null(input$datasets2)){return()}
    if(input$datasets2 %in% datasetProp()$dataset){
      min = datasetProp()$axisMin[datasetProp()$dataset == input$datasets2]
      max = datasetProp()$axisMax[datasetProp()$dataset == input$datasets2]
    }else{
      min = 0
      max = 1
    }
    conditionalPanel(condition = "input.axisLimBool2==true",
      numericInput("axisMin2","Min:",value=min),
      numericInput("axisMax2","Max:",value=max)
    )    
  })
  
  output$supportInterval <- renderUI({
    if(is.null(input$datasets)){return()}
    if(input$datasets %in% datasetProp()$dataset){
      val = datasetProp()$supportInterval[datasetProp()$dataset == input$datasets]
    }else{
      val = FALSE}
    checkboxInput('supportInterval', 'Plot base pair intervals (e.g., Joint linkage support intervals)?', val)
  })
  output$supportInterval2 <- renderUI({
    if(is.null(input$datasets2)){return()}
    if(input$datasets2 %in% datasetProp()$dataset){
      val = datasetProp()$supportInterval[datasetProp()$dataset == input$datasets2]
    }else{
      val = FALSE}
    checkboxInput('supportInterval2', 'Plot base pair intervals (e.g., Joint linkage support intervals)?', val)
  })
  
  output$SIbpStart <- renderUI({
    if(is.null(input$datasets)){return()}
    cols <- varnames()
    if(input$datasets %in% datasetProp()$dataset){
      selected = datasetProp()$SIbpStart[datasetProp()$dataset == input$datasets]
      selectedEnd = datasetProp()$SIbpEnd[datasetProp()$dataset == input$datasets]
    }else{
      selected = names(cols[2])
      selectedEnd = names(cols[2])
    }
    conditionalPanel(condition = "input.supportInterval==true",
        selectizeInput("SIbpStart", "Interval Base Pair Start:", choices = as.list(cols), selected = selected, multiple = FALSE, options = list(dropdownParent="body")),
        selectizeInput("SIbpEnd", "Interval Base Pair End:", choices = as.list(cols), selected = selectedEnd, multiple = FALSE, options = list(dropdownParent="body"))
    )
  })   
  output$SIbpStart2 <- renderUI({
    if(is.null(input$datasets2)){return()}
    cols <- varnames2()
    if(input$datasets2 %in% datasetProp()$dataset){
      selected = datasetProp()$SIbpStart[datasetProp()$dataset == input$datasets2]
      selectedEnd = datasetProp()$SIbpEnd[datasetProp()$dataset == input$datasets2]
    }else{
      selected = names(cols[2])
      selectedEnd = names(cols[2])
    }
    conditionalPanel(condition = "input.supportInterval2==true",
      selectizeInput("SIbpStart2", "Interval Base Pair Start:", choices = as.list(cols), selected = selected, multiple = FALSE, options = list(dropdownParent="body")),
      selectizeInput("SIbpEnd2", "Interval Base Pair End:", choices = as.list(cols), selected = selectedEnd, multiple = FALSE, options = list(dropdownParent="body"))
    )
  })
  
  output$SIyAxisColumn <- renderUI({
    if(is.null(input$datasets)){return()}
    cols <- varnames()       
    if(input$datasets %in% datasetProp()$dataset){  
      selected = datasetProp()$SIyAxisColumn[datasetProp()$dataset == input$datasets]
    }else{
      #selected = names(cols[10])
      selected = as.character(cols[10])
    }        
    conditionalPanel(condition = "input.supportInterval==true",                     
                     selectizeInput("SIyAxisColumn", "Support Interval Y-axis column:", choices = as.list(cols), selected = selected, multiple = FALSE, options = list(dropdownParent="body"))
    )
  })
  outputOptions(output, "SIyAxisColumn", suspendWhenHidden=FALSE)  
  output$SIyAxisColumn2 <- renderUI({
    if(is.null(input$datasets2)){return()}
    cols <- varnames2()       
    if(input$datasets2 %in% datasetProp()$dataset){  
      selected = datasetProp()$SIyAxisColumn[datasetProp()$dataset == input$datasets2]
    }else{
      #selected = names(cols[10])
      selected = as.character(cols[10])
    }        
    conditionalPanel(condition = "input.supportInterval2==true",                     
      selectizeInput("SIyAxisColumn2", "Support Interval Y-axis column:", choices = as.list(cols), selected = selected, multiple = FALSE, options = list(dropdownParent="body"))
    )
  })
  outputOptions(output, "SIyAxisColumn2", suspendWhenHidden=FALSE)
  
  output$SIaxisLimBool <- renderUI({
    if(is.null(input$datasets)){return()}
    if(input$datasets %in% datasetProp()$dataset){
      val = datasetProp()$SIaxisLimBool[datasetProp()$dataset == input$datasets]
    } else if (values$organism %in% legumeInfo.organisms) {
      val <- FALSE
    }else{
      val = TRUE
    }
    conditionalPanel(condition = "input.supportInterval==true", 
      checkboxInput('SIaxisLimBool', 'Set Support Interval Y-axis Limits?', val)
    )
  })  
  output$SIaxisLimBool2 <- renderUI({
    if(is.null(input$datasets2)){return()}
    if(input$datasets2 %in% datasetProp()$dataset){
      val = datasetProp()$SIaxisLimBool[datasetProp()$dataset == input$datasets2]
    } else if (values$organism2 %in% legumeInfo.organisms) {
      val <- FALSE
    }else{
      val = TRUE
    }
    conditionalPanel(condition = "input.supportInterval2==true", 
      checkboxInput('SIaxisLimBool2', 'Set Support Interval Y-axis Limits?', val)
    )
  })
  
  output$SIaxisLim <- renderUI({    
    if(is.null(input$datasets)){return()}
    if(input$datasets %in% datasetProp()$dataset){
      min = datasetProp()$axisMin[datasetProp()$dataset == input$datasets]
      max = datasetProp()$axisMax[datasetProp()$dataset == input$datasets]
    }else{
      min = 0
      max = 1
    }
    conditionalPanel(condition = "input.supportInterval==true && input.SIaxisLimBool==true",
                     numericInput("SIaxisMin","Min:",value=min),
                     numericInput("SIaxisMax","Max:",value=max)
    )    
  })  
  output$SIaxisLim2 <- renderUI({    
    if(is.null(input$datasets2)){return()}
    if(input$datasets2 %in% datasetProp()$dataset){
      min = datasetProp()$axisMin[datasetProp()$dataset == input$datasets2]
      max = datasetProp()$axisMax[datasetProp()$dataset == input$datasets2]
    }else{
      min = 0
      max = 1
    }
    conditionalPanel(condition = "input.supportInterval2==true && input.SIaxisLimBool2==true",
      numericInput("SIaxisMin2","Min:",value=min),
      numericInput("SIaxisMax2","Max:",value=max)
    )    
  })
  
  #builds list of multiple selection boxes for traits that have multiple columns in dataset
  output$traitColBoxes <- renderUI({
    if(is.null(input$datasets)){return()}
    if(input$plotAll == TRUE){return()}
    lapply(input$traitColumns, function(i) {      
      traits <- c("Select All",sort(unique(values[[input$datasets]][,i])))
      selectizeInput(inputId=i, label=paste0("Select ",i),traits,
                  selected=traits[2],
                  multiple=TRUE, options = list(dropdownParent="body",plugins=list("remove_button")))
    })
  })
  output$traitColBoxes2 <- renderUI({
    if(is.null(input$datasets2)){return()}
    if(input$plotAll2 == TRUE){return()}
    lapply(input$traitColumns2, function(i) {      
      traits <- c("Select All ",sort(unique(values[[input$datasets2]][,i])))
      selectizeInput(inputId=paste0(i, "2"), label=paste0("Select ",i),traits,
        selected=traits[2],
        multiple=TRUE, options = list(dropdownParent="body",plugins=list("remove_button")))
    })
  })
  
  observe({
    lapply(input$traitColumns, function(i){
      if("Select All" %in% input[[i]]){
        selected_choices <- sort(unique(values[[input$datasets]][,i]))
        updateSelectizeInput(session, i, selected = selected_choices)
      }
    })
  })
  observe({
    lapply(input$traitColumns2, function(i){
      i2 <- paste0(i, "2")
      if("Select All " %in% input[[i2]]){
        selected_choices <- sort(unique(values[[input$datasets2]][,i]))
        updateSelectizeInput(session, i2, selected = selected_choices)
      }
    })
  })
  
  #checkbox to suppress plot legend
  output$legend <- renderUI({
    if(is.null(input$datasets)){return()}
    checkboxInput('legend', 'Suppress Legend', FALSE)
  })
  outputOptions(output, "legend", suspendWhenHidden=FALSE)
  output$legend2 <- renderUI({
    if(is.null(input$datasets2)){return()}
    checkboxInput('legend2', 'Suppress Legend', FALSE)
  })
  outputOptions(output, "legend2", suspendWhenHidden=FALSE)
  
  #checkbox for whether to filter for only overlapping SNPs
  output$overlaps <- renderUI({
    if(is.null(input$datasets)){return()}
#    if(input$plotAll == TRUE){return()}
    checkboxInput('overlaps', 'Show only overlapping SNPs', FALSE)
  })  
  outputOptions(output, "overlaps", suspendWhenHidden=FALSE)
  output$overlaps2 <- renderUI({
    if(is.null(input$datasets2)){return()}
    #    if(input$plotAll2 == TRUE){return()}
    checkboxInput('overlaps2', 'Show only overlapping SNPs', FALSE)
  })  
  outputOptions(output, "overlaps2", suspendWhenHidden=FALSE)
  
  #how many traits must overlap to be included in output, 1 means traits that overlap with themselves will be included
  output$numOverlapTraits <- renderUI({
    numericInput("numOverlaps", "Minimum number of overlaps?", value=2,min=1)    
  })
  output$numOverlapTraits2 <- renderUI({
    numericInput("numOverlaps2", "Minimum number of overlaps?", value=2,min=1)    
  })
  
  #how big is the window when calculating whether two snps overlap
  output$overlapSize <- renderUI({
    numericInput(inputId="overlapSize",label="Overlap size around each point:",min=1,max=.5e6,value=10000)
  })
  output$overlapSize2 <- renderUI({
    numericInput(inputId="overlapSize2",label="Overlap size around each point:",min=1,max=.5e6,value=10000)
  })
  
  output$selectChr <- renderUI({
    if(is.null(values$organism)){return()}
    selectInput("chr", "Chromosome:",chrName[values$organism][[1]],selectize = FALSE)
    #selectInput("chr", "Chromosome:",1:length(chrSize[input$organism][[1]]),selectize = FALSE)
  })
  outputOptions(output, "selectChr", suspendWhenHidden=FALSE)
  output$selectChr2 <- renderUI({
    if(is.null(values$organism2)){return()}
    selectInput("chr2", "Chromosome:",chrName[values$organism2][[1]],selectize = FALSE)
  })
  outputOptions(output, "selectChr2", suspendWhenHidden=FALSE)
  
  output$selectedOut <- renderUI({
    numericInput("selected", "", value=100000)
  })
  outputOptions(output, "selectedOut", suspendWhenHidden=FALSE)
  output$selectedOut2 <- renderUI({
    numericInput("selected2", "", value=100000)
  })
  outputOptions(output, "selectedOut2", suspendWhenHidden=FALSE)

  output$windowOut <- renderUI({
    #sliderInput(inputId="window",label="Window size around selected point:",min=-1e6,max=1e6,value=c(-7.5e5,7.5e5))
    sliderInput(inputId="window",label="Window size around selected point:",min=1000,max=.5e6,value=2.5e5)
  })
  outputOptions(output, "windowOut", suspendWhenHidden=FALSE)
  output$windowOut2 <- renderUI({
    sliderInput(inputId="window2",label="Window size around selected point:",min=1000,max=.5e6,value=2.5e5)
  })
  outputOptions(output, "windowOut2", suspendWhenHidden=FALSE)
  
  #returns datasets from uploaded file
  getdata <- reactive({
    if(is.null(input$datasets)){return()}
    values[[input$datasets]]
  })
  getdata2 <- reactive({
    if(is.null(input$datasets2)){return()}
    values[[input$datasets2]]
  })
  
  #builds list of column names and type in dataset
  varnames <- reactive({
    if(is.null(input$datasets)) return()
    dat <- getdata_class()
    vars <- names(dat)
    names(vars) <- paste(vars, " {", dat, "}", sep = "")
    vars
  })
  varnames2 <- reactive({
    if(is.null(input$datasets2)) return()
    dat <- getdata_class2()
    vars <- names(dat)
    names(vars) <- paste(vars, " {", dat, "}", sep = "")
    vars
  })
  
  getdata_class <- reactive({
    if(is.null(input$datasets)) return()
    cls <- sapply(getdata(), function(x) class(x)[1])
    gsub("ordered","factor", cls)
  })
  getdata_class2 <- reactive({
    if(is.null(input$datasets2)) return()
    cls <- sapply(getdata2(), function(x) class(x)[1])
    gsub("ordered","factor", cls)
  })
  
  #Function to handle loading of data from a file or rObject
  loadUserData <- function(filename, uFile) {  
    ext <- file_ext(filename)
    # objname <- robjname <- sub(paste(".",ext,sep = ""),"",basename(filename))
    objname <- sub(paste(".",ext,sep = ""),"",basename(filename))
    ext <- tolower(ext)
    
    if(ext == 'rda' || ext == 'rdata') {
      # objname will hold the name of the object(s) inside the R datafile
      robjname <- load(uFile)
      
      if(length(robjname) > 1) {
        
        #loaded.values <- data.frame(get(robjname[-which(robjname == "description")]))
        loaded.values <- data.frame(get(robjname[1]))
        #values[[paste0(objname,"_descr")]] <- get("description")
        
      } else {
        
        loaded.values <- data.frame(get(robjname))  # only work with data.frames
      }
    }
    
    appendSNPs <- isolate(input$appendSNPs)
    if (!appendSNPs) {
      # add new datasets to the datasetToOrganism map
      values$datasetToOrganism[[objname]] <- values$organism
    }
    if(length(values[['datasetlist']]) == 0 || values[['datasetlist']][1] == '') {
      values[['datasetlist']] <- c(objname)
    } else if (!appendSNPs) {
      values[['datasetlist']] <- unique(c(objname,values[['datasetlist']]))
    }
    
    if(ext == 'sav') {
      loaded.values <- as.data.frame(as.data.set(spss.system.file(uFile)))
    } else if(ext == 'dta') {
      loaded.values <- read.dta(uFile)
    } else{
      loaded.values <- read.csv(uFile, header=input$header, sep=input$sep,stringsAsFactors=FALSE)
    }

    if (appendSNPs) {
      values[[input$datasets]]$totalBP <- NULL
      names(loaded.values) <- names(values[[input$datasets]])
      values[[input$datasets]] <- rbind(values[[input$datasets]], loaded.values)
    } else {
      values[[objname]] <- loaded.values
    }
  }
  loadUserData2 <- function(filename, uFile) {  
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

    appendSNPs <- isolate(input$appendSNPs2)
    if (!appendSNPs) {
      # add new datasets to the datasetToOrganism map
      values$datasetToOrganism[[objname]] <- values$organism2
    }
    if (length(values[['datasetlist']]) == 0 || values[['datasetlist']][1] == '') {
      values[['datasetlist']] <- c(objname)
    } else if (!appendSNPs) {
      values[['datasetlist']] <- unique(c(objname, values[['datasetlist']]))
    }

    if (ext == 'sav') {
      loaded.values <- as.data.frame(as.data.set(spss.system.file(uFile)))
    } else if(ext == 'dta') {
      loaded.values <- read.dta(uFile)
    } else{
      loaded.values <- read.csv(uFile, header = input$header2, sep = input$sep2, stringsAsFactors=FALSE)
    }

    if (appendSNPs) {
      values[[input$datasets2]]$totalBP <- NULL
      names(loaded.values) <- names(values[[input$datasets2]])
      values[[input$datasets2]] <- rbind(values[[input$datasets2]], loaded.values)
    } else {
      values[[objname]] <- loaded.values
    }
  }
  
  # Load data from .csv files at a remote URL
  loadRemoteData <- function(trait, traitUrl) {
    ext <- file_ext(traitUrl)
    objname <- sub(paste(".", ext, sep = ""), "", basename(traitUrl))

    appendSNPs <- isolate(input$appendSNPs)
    if (!appendSNPs) {
      # add new datasets to the datasetToOrganism map
      values$datasetToOrganism[[objname]] <- values$organism
    }
    if(length(values[['datasetlist']]) == 0 || values[['datasetlist']][1] == '') {
      values[['datasetlist']] <- c(objname)
    } else if (!appendSNPs) {
      values[['datasetlist']] <- unique(c(objname, values[['datasetlist']]))
    }

    loaded.values <- load.gwas.remote(values$organism, traitUrl, trait)

    if (appendSNPs) {
      values[[input$datasets]]$totalBP <- NULL
      names(loaded.values) <- names(values[[input$datasets]])
      values[[input$datasets]] <- rbind(values[[input$datasets]], loaded.values)
    } else {
      values[[objname]] <- loaded.values
    }
  }
  loadRemoteData2 <- function(trait, traitUrl) {
    ext <- file_ext(traitUrl)
    objname <- sub(paste(".", ext, sep = ""), "", basename(traitUrl))

    appendSNPs <- isolate(input$appendSNPs2)
    if (!appendSNPs) {
      # add new datasets to the datasetToOrganism map
      values$datasetToOrganism[[objname]] <- values$organism2
    }
    if (length(values[['datasetlist']]) == 0 || values[['datasetlist']][1] == '') {
      values[['datasetlist']] <- c(objname)
    } else if (!appendSNPs) {
      values[['datasetlist']] <- unique(c(objname, values[['datasetlist']]))
    }

    loaded.values <- load.gwas.remote(values$organism2, traitUrl, trait)

    if (appendSNPs) {
      values[[input$datasets2]]$totalBP <- NULL
      names(loaded.values) <- names(values[[input$datasets2]])
      values[[input$datasets2]] <- rbind(values[[input$datasets2]], loaded.values)
    } else {
      values[[objname]] <- loaded.values
    }
  }
  
  #add a totalBP column to an input dataset if not already present
  calculateTotalBP <- reactive({ 
    if("totalBP" %in% colnames(values[[input$datasets]])){
      
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
      cumBP<-c(0,cumsum(as.numeric(chrSize[values$organism][[1]])))
      #to order by desired chromosome add factor levels in the desired order to the chrColumn, any chr names that differ in gwas file compared
      #to organism file will turn into NA
      values[[input$datasets]][,input$chrColumn] <- factor(values[[input$datasets]][,input$chrColumn],levels=chrName[values$organism][[1]])
      values[[input$datasets]] <- values[[input$datasets]][order(values[[input$datasets]][,input$chrColumn],values[[input$datasets]][,input$bpColumn]),]
      numeachchr<-aggregate(values[[input$datasets]][,input$bpColumn],list(values[[input$datasets]][,input$chrColumn]),length)
#      adjust<-rep(cumBP[1],numeachchr$x[numeachchr$Group.1==1])            
      adjust <- numeric()
      for (i in 1:(length(cumBP)-1)){#max(unique(values[[input$datasets]][,input$chrColumn]))){
        if(length(numeachchr$x[numeachchr$Group.1==chrName[values$organism][[1]][i]])==0){next;}
        adjust<-c(adjust,rep(cumBP[i],numeachchr$x[numeachchr$Group.1==chrName[values$organism][[1]][i]]))
      }
      #newval <- values[[input$datasets]][600,input$bpColumn]+adjust[600]      
      values[[input$datasets]]$totalBP <- values[[input$datasets]][,input$bpColumn]+adjust
      
      #values[[input$datasets]] <- adply(values[[input$datasets]],1,function(x){data.frame(totalBP=sum(x[[input$bpColumn]],chrSize$bp[chrSize$chr %in% if(x[[input$chrColumn]]==1) 0 else c(1:(x[[input$chrColumn]]-1))]))})
    }
   if(input$supportInterval == TRUE){
      if("SIbpStartTotal" %in% colnames(values[[input$datasets]])){
        
      }else{
         
        cumBP<-c(0,cumsum(as.numeric(chrSize[values$organism][[1]])))
        values[[input$datasets]][,input$chrColumn] <- factor(values[[input$datasets]][,input$chrColumn],levels=chrName[values$organism][[1]])
        values[[input$datasets]] <- values[[input$datasets]][order(values[[input$datasets]][,input$chrColumn],values[[input$datasets]][,input$SIbpStart]),]
        numeachchr<-aggregate(values[[input$datasets]][,input$SIbpStart],list(values[[input$datasets]][,input$chrColumn]),length)
        adjust <- numeric()
        for (i in 1:(length(cumBP)-1)){#max(unique(values[[input$datasets]][,input$chrColumn]))){
          if(length(numeachchr$x[numeachchr$Group.1==chrName[values$organism][[1]][i]])==0){next;}
          adjust<-c(adjust,rep(cumBP[i],numeachchr$x[numeachchr$Group.1==chrName[values$organism][[1]][i]]))
        }
        values[[input$datasets]]$SIbpStartTotal <- values[[input$datasets]][,input$SIbpStart]+adjust    
      }
    
      if("SIbpEndTotal" %in% colnames(values[[input$datasets]])){
        
      }else{
        
        cumBP<-c(0,cumsum(as.numeric(chrSize[values$organism][[1]])))
        values[[input$datasets]][,input$chrColumn] <- factor(values[[input$datasets]][,input$chrColumn],levels=chrName[values$organism][[1]])
        values[[input$datasets]] <- values[[input$datasets]][order(values[[input$datasets]][,input$chrColumn],values[[input$datasets]][,input$SIbpEnd]),]
        numeachchr<-aggregate(values[[input$datasets]][,input$SIbpEnd],list(values[[input$datasets]][,input$chrColumn]),length)
        adjust <- numeric()
        for (i in 1:(length(cumBP)-1)){#max(unique(values[[input$datasets]][,input$chrColumn]))){
          if(length(numeachchr$x[numeachchr$Group.1==chrName[values$organism][[1]][i]])==0){next;}
          adjust<-c(adjust,rep(cumBP[i],numeachchr$x[numeachchr$Group.1==chrName[values$organism][[1]][i]]))
        }
        values[[input$datasets]]$SIbpEndTotal <- values[[input$datasets]][,input$SIbpEnd]+adjust    
      }
    } #end SI total bp calculation
  })#end calculateTotalBP
  calculateTotalBP2 <- reactive({ 
    if("totalBP" %in% colnames(values[[input$datasets2]])){
      
    }else{
      cumBP<-c(0,cumsum(as.numeric(chrSize[values$organism2][[1]])))
      #to order by desired chromosome add factor levels in the desired order to the chrColumn, any chr names that differ in gwas file compared
      #to organism file will turn into NA
      values[[input$datasets2]][,input$chrColumn2] <- factor(values[[input$datasets2]][,input$chrColumn2],levels=chrName[values$organism2][[1]])
      values[[input$datasets2]] <- values[[input$datasets2]][order(values[[input$datasets2]][,input$chrColumn2],values[[input$datasets2]][,input$bpColumn2]),]
      numeachchr<-aggregate(values[[input$datasets2]][,input$bpColumn2],list(values[[input$datasets2]][,input$chrColumn2]),length)
      adjust <- numeric()
      for (i in 1:(length(cumBP)-1)){#max(unique(values[[input$datasets2]][,input$chrColumn2]))){
        if(length(numeachchr$x[numeachchr$Group.1==chrName[values$organism2][[1]][i]])==0){next;}
        adjust<-c(adjust,rep(cumBP[i],numeachchr$x[numeachchr$Group.1==chrName[values$organism2][[1]][i]]))
      }
      #newval <- values[[input$datasets2]][600,input$bpColumn2]+adjust[600]      
      values[[input$datasets2]]$totalBP <- values[[input$datasets2]][,input$bpColumn2]+adjust
      
      #values[[input$datasets2]] <- adply(values[[input$datasets2]],1,function(x){data.frame(totalBP=sum(x[[input$bpColumn2]],chrSize$bp[chrSize$chr %in% if(x[[input$chrColumn2]]==1) 0 else c(1:(x[[input$chrColumn2]]-1))]))})
    }
    if(input$supportInterval2 == TRUE){
      if("SIbpStartTotal" %in% colnames(values[[input$datasets2]])){
        
      }else{
        
        cumBP<-c(0,cumsum(as.numeric(chrSize[values$organism2][[1]])))
        values[[input$datasets2]][,input$chrColumn2] <- factor(values[[input$datasets2]][,input$chrColumn2],levels=chrName[values$organism2][[1]])
        values[[input$datasets2]] <- values[[input$datasets2]][order(values[[input$datasets2]][,input$chrColumn2],values[[input$datasets2]][,input$SIbpStart2]),]
        numeachchr<-aggregate(values[[input$datasets2]][,input$SIbpStart2],list(values[[input$datasets2]][,input$chrColumn2]),length)
        adjust <- numeric()
        for (i in 1:(length(cumBP)-1)){#max(unique(values[[input$datasets2]][,input$chrColumn2]))){
          if(length(numeachchr$x[numeachchr$Group.1==chrName[values$organism2][[1]][i]])==0){next;}
          adjust<-c(adjust,rep(cumBP[i],numeachchr$x[numeachchr$Group.1==chrName[values$organism2][[1]][i]]))
        }
        values[[input$datasets2]]$SIbpStartTotal <- values[[input$datasets2]][,input$SIbpStart2]+adjust    
      }
      
      if("SIbpEndTotal" %in% colnames(values[[input$datasets2]])){
        
      }else{
        
        cumBP<-c(0,cumsum(as.numeric(chrSize[values$organism2][[1]])))
        values[[input$datasets2]][,input$chrColumn2] <- factor(values[[input$datasets2]][,input$chrColumn2],levels=chrName[values$organism2][[1]])
        values[[input$datasets2]] <- values[[input$datasets2]][order(values[[input$datasets2]][,input$chrColumn2],values[[input$datasets2]][,input$SIbpEnd2]),]
        numeachchr<-aggregate(values[[input$datasets2]][,input$SIbpEnd2],list(values[[input$datasets2]][,input$chrColumn2]),length)
        adjust <- numeric()
        for (i in 1:(length(cumBP)-1)){#max(unique(values[[input$datasets2]][,input$chrColumn2]))){
          if(length(numeachchr$x[numeachchr$Group.1==chrName[values$organism2][[1]][i]])==0){next;}
          adjust<-c(adjust,rep(cumBP[i],numeachchr$x[numeachchr$Group.1==chrName[values$organism2][[1]][i]]))
        }
        values[[input$datasets2]]$SIbpEndTotal <- values[[input$datasets2]][,input$SIbpEnd2]+adjust    
      }
    } #end SI total bp calculation
  })#end calculateTotalBP2
  
  output$pChart <- renderChart({
    #this function makes the chromsomeview chart  
    #subset whole chart based on selection
    chromChart <- values[[input$datasets]]
    chromChart <- chromChart[chromChart[,input$chrColumn]==input$chr,]
    
    if(input$plotAll==FALSE){
      for(i in input$traitColumns){
        chromChart <- chromChart[chromChart[,i] %in% input[[i]],]
      }    
      if(length(input$traitColumns) > 1){
        chromChart$trait <- do.call(paste,c(chromChart[,input$traitColumns],sep="_"))    
      }else{
        chromChart$trait <- chromChart[,input$traitColumns]
      }
    }else{
      chromChart$trait <- input$datasets
    }
            
    #Separate Support Interval data from GWAS data, if support, GWAS data is assumed to be anything that has an NA in the SIbpStart column
    if(input$supportInterval == TRUE){
      SIchart <- chromChart[!(is.na(chromChart[,input$SIbpStart])),]
      chromChart <- chromChart[is.na(chromChart[,input$SIbpStart]),]
    }        
    #check if there is any data for the selected traits
    chromChart <- chromChart[!(is.na(chromChart[,input$bpColumn])),]
    chromChart <- chromChart[!(is.na(chromChart[,input$yAxisColumn])),]
    
    #if checked, filter for only overlapping SNPs
    if(!is.null(input$overlaps) & input$overlaps == TRUE){
      chromChart <- findGWASOverlaps(chromChart)
    }            
    
    if(nrow(chromChart)==0){ #nothing is in the window, but lets still make a data.frame
      chromChart <- values[[input$datasets]][1,]
      chromChart[,input$yAxisColumn] <- -1    
      if(length(input$traitColumns) > 1){
        chromChart$trait <- do.call(paste,c(chromChart[,input$traitColumns],sep="_"))    
      }else{
        chromChart$trait <- chromChart[,input$traitColumns]
      }             
    }    
    colorTable <- colorTable() 
    
    #take -log10 of y-axis column if requested
    if(input$logP == TRUE && chromChart[1,input$yAxisColumn] != -1){
      chromChart[,input$yAxisColumn] <- -log(chromChart[,input$yAxisColumn],10)
    }        
    
    #check if there is too much data (>2500 data points), trim to 2500
    if(nrow(chromChart)>2500){
      cutVal <- sort(chromChart[,input$yAxisColumn],decreasing = T)[2500]
      chromChart <- chromChart[chromChart[,input$yAxisColumn] >= cutVal,]
    }        
    
    #calculate window for plotband
    pbWin <- isolate({
      center <- as.numeric(input$selected[[1]])
      winHigh <- center + input$window[1]
      winLow <- center - input$window[1]
      list(winLow=winLow,winHigh=winHigh)
    })    
    
    pkTable <- data.frame(x=chromChart[,input$bpColumn],y=chromChart[,input$yAxisColumn],trait=chromChart$trait,
                          #                          name=sprintf("<table cellpadding='4' style='line-height:1.5'><tr><th>%1$s</th></tr><tr><td align='left'>Y-value: %2$s<br>Location: %3$s<br>Base Pair: %4$s<br>SNP: %5$s<br>Chromosome: %6$s</td></tr></table>",
                          #name=sprintf("<table cellpadding='4' style='line-height:1.5'><tr><th>%1$s</th></tr><tr><td align='left'>Y-value: %2$s<br>Base Pairs: %3$s<br>Chromosome: %4$s</td></tr></table>",
                          name=sprintf("Base Pair: %1$s<br/>Chromosome: %2$s<br/>",
#                                       chromChart$trait,
#                                       chromChart[,input$yAxisColumn],
                                       #                                       pk$loc,
                                       prettyNum(chromChart[,input$bpColumn], big.mark = ","),
                                       #                                       pk$SNP,
                                       chromChart[,input$chrColumn]
                          ),
                          url="http://danforthcenter.org",
                          chr=chromChart[,input$chrColumn],
                          bp=chromChart[,input$bpColumn],stringsAsFactors=FALSE)
    pkSeries <- lapply(split(pkTable, pkTable$trait), function(x) {
      res <- lapply(split(x, rownames(x)), as.list)
      names(res) <- NULL
      res <- res[order(sapply(res, function(x) x$x))]
      return(res)
    })

    #build JL series
    if(input$supportInterval==TRUE){
      if(nrow(SIchart)==0){ #nothing is in the window, but lets still make a data.frame
        SIchart <- values[[input$datasets]][1,]
        SIchart[,input$SIyAxisColumn] <- -1    
        if(length(input$traitColumns) > 1){
          SIchart$trait <- do.call(paste,c(SIchart[,input$traitColumns],sep="_"))    
        }else{
          SIchart$trait <- SIchart[,input$traitColumns]
        }             
      }
      SIchart$loc_el <- SIchart$trait
      SIchart$trait <- paste(SIchart$trait,"Int",sep="_")
      SIchart <- SIchart[order(SIchart[[input$SIbpStart]]),]
      jlTable <- adply(SIchart,1,function(x) {data.frame(x=c(x[[input$SIbpStart]],x[[input$SIbpEnd]],x[[input$SIbpEnd]]),y=c(x[[input$SIyAxisColumn]],x[[input$SIyAxisColumn]],NA),trait=x$trait,
                                                         name=sprintf("<table cellpadding='4' style='line-height:1.5'><tr><td align='left'>Interval: %1$s-%2$s<br>Chromosome: %3$s</td></tr></table>",
#                                                                      x$trait,
#                                                                      x[[input$SIyAxisColumn]],
                                                                      prettyNum(x[[input$SIbpStart]], big.mark = ","),
                                                                      prettyNum(x[[input$SIbpEnd]], big.mark = ","),
                                                                      x[[input$chrColumn]]
                                                         ),loc_el=x$loc_el,bp=x[[input$bpColumn]],chr=x[[input$chrColumn]],stringsAsFactors=FALSE
                                                         #                                                   
                                                         #                                                   totalBP=x$totalBP,
                                                         #                                                   chr=x$Chromosome,stringsAsFactors=FALSE
      )}#end jlTable and function
      )#end adply
      jlTable <- jlTable[,c("x","y","trait","name","loc_el","bp","chr")]
      #jlTable <- jlTable[order(jlTable$x),]
    }#end build jlTable if support intervals        
    
    a <- rCharts::Highcharts$new()
    a$LIB$url <- 'highcharts/' #use the local copy of highcharts, not the one installed by rCharts
    a$xAxis(title = list(text = "Base Pairs"),startOnTick=TRUE,min=1,max=chrSize[values$organism][[1]][as.numeric(input$chr)],endOnTick=FALSE,
            plotBands = list(list(from=pbWin$winLow,to=pbWin$winHigh,color='rgba(68, 170, 213, 0.4)')))
#    a$xAxis(title = list(text = "Base Pairs"),startOnTick=TRUE,min=1,max=30000000,endOnTick=TRUE)
    
    if(input$axisLimBool == TRUE){
      a$yAxis(title=list(text=input$yAxisColumn),min=input$axisMin,max=input$axisMax,startOnTick=FALSE)
    }else{
      a$yAxis(title=list(text=input$yAxisColumn),startOnTick=FALSE)      
    }    

    if(input$supportInterval==TRUE){
      if(input$SIaxisLimBool == TRUE){
        a$yAxis(title=list(text=input$SIyAxisColumn),min=input$SIaxisMin,max=input$SIaxisMax,gridLineWidth=0,minorGridLineWidth=0,startOnTick=FALSE,opposite=TRUE,replace=FALSE)
      }else{
        a$yAxis(title=list(text=input$SIyAxisColumn),gridLineWidth=0,minorGridLineWidth=0,startOnTick=FALSE,opposite=TRUE,replace=FALSE)   
      }
      
      if(SIchart[1,input$SIyAxisColumn] != -1){
        d_ply(jlTable,.(trait),function(x){
          a$series(
            data = toJSONArray2(x,json=F,names=T),
            type = "line",
            name = unique(x$trait),
            yAxis=1,           
            color = colorTable$color[colorTable$trait == as.character(unique(x$loc_el))])})            
      }
    }

   if(chromChart[1,input$yAxisColumn] != -1){    
      invisible(sapply(pkSeries, function(x) {if(length(x)==0){return()};a$series(data = x, type = "scatter", turboThreshold=5000, name = paste0(x[[1]]$trait), color = colorTable$color[colorTable$trait == as.character(x[[1]]$trait)])}))
    }
    a$chart(zoomType="x", alignTicks=FALSE,events=list(click = "#!function(event) {this.tooltip.hide();}!#"))
    a$title(text=paste(input$datasets,"Results for Chromosome",input$chr,sep=" "))
    a$subtitle(text="Rollover for more info. Drag chart area to zoom. Click point for zoomed annotated plot.")
    
    a$plotOptions(
      scatter = list(
        cursor = "pointer",
        #           tooltip = list(
        #             pointFormat = "#! function() { return this.point.name; } !#"
        #           ),          
        point = list(
          events = list(
            #click = "#! function() { window.open(this.options.url); } !#")), #open webpage
            #click = "#! function(event) {alert(this.url);} !#")), #display popup
            #click = "#! function(event) {console.log(this);} !#")), #write object to log
            click = "#! function(){$('input#selected').val(this.options.bp); $('input#selected').trigger('change');} !#")),
        marker = list(
          symbol = "circle",
          radius = 5
#          states = list(hover = list(enabled = TRUE))
        ),
        tooltip = list(
          headerFormat = "<b>{series.name}</b><br/>{point.key}<br/>Y-value: {point.y}<br/>",
          pointFormat = "",
          followPointer = TRUE
        )
 #       states = list(hover = list(marker = list(enabled = FALSE)))
#        tooltip = list(
#          followPointer = TRUE
#        )
      ),
      line = list(
        lineWidth = 10,
        dashStyle = 'Solid',
        cursor = "pointer",
#        stickyTracking=FALSE,
        point = list(
          events = list(
            #click = "#! function() { window.open(this.url); } !#")), #open webpage
            #click = "#! function(event) {alert(this.url);} !#")), #display popup
            #click = "#! function(event) {console.log(this);} !#")), #write object to log
            #click = "#! function(){$('select#chr').val(this.options.chr); $('select#chr').trigger('change'); $('input#selected').val(this.options.bp); $('input#selected').trigger('change');  $('ul#methodtabs li').eq(0).removeClass('active'); $('ul#methodtabs li').eq(1).addClass('active'); $('.tab-content div').toggleClass('active'); $('#pChart').trigger('shown')} !#")),
            click = "#! function(){$('input#selected').val(this.options.bp); $('input#selected').trigger('change');} !#")),
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
    #a$tooltip(useHTML = T, formatter = "#! function() { return this.point.name; } !#")
    #a$tooltip(snap = 5, useHTML = T, formatter = "#! function() { return this.point.name; } !#")
    a$exporting(enabled=TRUE,filename='chromChart',sourceWidth=2000)
    a$set(dom = 'pChart')
    return(a)
    
  })#end pchart
  output$pChart2 <- renderChart({
    #this function makes the chromsomeview chart  
    #subset whole chart based on selection
    chromChart <- values[[input$datasets2]]
    chromChart <- chromChart[chromChart[,input$chrColumn2]==input$chr2,]
    
    if(input$plotAll2==FALSE){
      for(i in input$traitColumns2){
        i2 <- paste0(i, "2")
        chromChart <- chromChart[chromChart[,i] %in% input[[i2]],]
      }    
      if(length(input$traitColumns2) > 1){
        chromChart$trait <- do.call(paste,c(chromChart[,input$traitColumns2],sep="_"))
      }else{
        chromChart$trait <- chromChart[,input$traitColumns2]
      }
    }else{
      chromChart$trait <- input$datasets2
    }
    
    #Separate Support Interval data from GWAS data, if support, GWAS data is assumed to be anything that has an NA in the SIbpStart column
    if(input$supportInterval2 == TRUE){
      SIchart <- chromChart[!(is.na(chromChart[,input$SIbpStart2])),]
      chromChart <- chromChart[is.na(chromChart[,input$SIbpStart2]),]
    }        
    #check if there is any data for the selected traits
    chromChart <- chromChart[!(is.na(chromChart[,input$bpColumn2])),]
    chromChart <- chromChart[!(is.na(chromChart[,input$yAxisColumn2])),]
    
    #if checked, filter for only overlapping SNPs
    if(!is.null(input$overlaps2) & input$overlaps2 == TRUE){
      chromChart <- findGWASOverlaps2(chromChart)
    }
    
    if(nrow(chromChart)==0){ #nothing is in the window, but lets still make a data.frame
      chromChart <- values[[input$datasets2]][1,]
      chromChart[,input$yAxisColumn2] <- -1    
      if(length(input$traitColumns2) > 1){
        chromChart$trait <- do.call(paste,c(chromChart[,input$traitColumns2],sep="_"))
      }else{
        chromChart$trait <- chromChart[,input$traitColumns2]
      }             
    }    
    colorTable <- colorTable2()
    
    #take -log10 of y-axis column if requested
    if(input$logP2 == TRUE && chromChart[1,input$yAxisColumn2] != -1){
      chromChart[,input$yAxisColumn2] <- -log(chromChart[,input$yAxisColumn2],10)
    }
    
    #check if there is too much data (>2500 data points), trim to 2500
    if(nrow(chromChart)>2500){
      cutVal <- sort(chromChart[,input$yAxisColumn2],decreasing = T)[2500]
      chromChart <- chromChart[chromChart[,input$yAxisColumn2] >= cutVal,]
    }
    
    #calculate window for plotband
    pbWin <- isolate({
      center <- as.numeric(input$selected2[[1]])
      winHigh <- center + input$window2[1]
      winLow <- center - input$window2[1]
      list(winLow=winLow,winHigh=winHigh)
    })
    
    pkTable <- data.frame(x=chromChart[,input$bpColumn2],y=chromChart[,input$yAxisColumn2],trait=chromChart$trait,
      name=sprintf("Base Pair: %1$s<br/>Chromosome: %2$s<br/>",
        prettyNum(chromChart[,input$bpColumn2], big.mark = ","),
        chromChart[,input$chrColumn2]
      ),
      url="http://danforthcenter.org",
      chr=chromChart[,input$chrColumn2],
      bp=chromChart[,input$bpColumn2],stringsAsFactors=FALSE)
    pkSeries <- lapply(split(pkTable, pkTable$trait), function(x) {
      res <- lapply(split(x, rownames(x)), as.list)
      names(res) <- NULL
      res <- res[order(sapply(res, function(x) x$x))]
      return(res)
    })
    
    #build JL series
    if(input$supportInterval2==TRUE){
      if(nrow(SIchart)==0){ #nothing is in the window, but lets still make a data.frame
        SIchart <- values[[input$datasets2]][1,]
        SIchart[,input$SIyAxisColumn2] <- -1    
        if(length(input$traitColumns2) > 1){
          SIchart$trait <- do.call(paste,c(SIchart[,input$traitColumns2],sep="_"))
        }else{
          SIchart$trait <- SIchart[,input$traitColumns2]
        }
      }
      SIchart$loc_el <- SIchart$trait
      SIchart$trait <- paste(SIchart$trait,"Int",sep="_")
      SIchart <- SIchart[order(SIchart[[input$SIbpStart2]]),]
      jlTable <- adply(SIchart,1,function(x) {data.frame(x=c(x[[input$SIbpStart2]],x[[input$SIbpEnd2]],x[[input$SIbpEnd2]]),y=c(x[[input$SIyAxisColumn2]],x[[input$SIyAxisColumn2]],NA),trait=x$trait,
        name=sprintf("<table cellpadding='4' style='line-height:1.5'><tr><td align='left'>Interval: %1$s-%2$s<br>Chromosome: %3$s</td></tr></table>",
          prettyNum(x[[input$SIbpStart2]], big.mark = ","),
          prettyNum(x[[input$SIbpEnd2]], big.mark = ","),
          x[[input$chrColumn2]]
        ),loc_el=x$loc_el,bp=x[[input$bpColumn2]],chr=x[[input$chrColumn2]],stringsAsFactors=FALSE
      )}#end jlTable and function
      )#end adply
      jlTable <- jlTable[,c("x","y","trait","name","loc_el","bp","chr")]
    }#end build jlTable if support intervals        
    
    a <- rCharts::Highcharts$new()
    a$LIB$url <- 'highcharts/' #use the local copy of highcharts, not the one installed by rCharts
    a$xAxis(title = list(text = "Base Pairs"),startOnTick=TRUE,min=1,max=chrSize[values$organism2][[1]][as.numeric(input$chr2)],endOnTick=FALSE,
      plotBands = list(list(from=pbWin$winLow,to=pbWin$winHigh,color='rgba(68, 170, 213, 0.4)')))
    
    if(input$axisLimBool2 == TRUE){
      a$yAxis(title=list(text=input$yAxisColumn2),min=input$axisMin2,max=input$axisMax2,startOnTick=FALSE)
    }else{
      a$yAxis(title=list(text=input$yAxisColumn2),startOnTick=FALSE)      
    }    
    
    if(input$supportInterval2==TRUE){
      if(input$SIaxisLimBool2 == TRUE){
        a$yAxis(title=list(text=input$SIyAxisColumn2),min=input$SIaxisMin2,max=input$SIaxisMax2,gridLineWidth=0,minorGridLineWidth=0,startOnTick=FALSE,opposite=TRUE,replace=FALSE)
      }else{
        a$yAxis(title=list(text=input$SIyAxisColumn2),gridLineWidth=0,minorGridLineWidth=0,startOnTick=FALSE,opposite=TRUE,replace=FALSE)   
      }
      
      if(SIchart[1,input$SIyAxisColumn2] != -1){
        d_ply(jlTable,.(trait),function(x){
          a$series(
            data = toJSONArray2(x,json=F,names=T),
            type = "line",
            name = unique(x$trait),
            yAxis=1,
            color = colorTable$color[colorTable$trait == as.character(unique(x$loc_el))])})            
      }
    }
    
    if(chromChart[1,input$yAxisColumn2] != -1){
      invisible(sapply(pkSeries, function(x) {if(length(x)==0){return()};a$series(data = x, type = "scatter", turboThreshold=5000, name = paste0(x[[1]]$trait), color = colorTable$color[colorTable$trait == as.character(x[[1]]$trait)])}))
    }
    a$chart(zoomType="x", alignTicks=FALSE,events=list(click = "#!function(event) {this.tooltip.hide();}!#"))
    a$title(text=paste(input$datasets2,"Results for Chromosome",input$chr2,sep=" "))
    a$subtitle(text="Rollover for more info. Drag chart area to zoom. Click point for zoomed annotated plot.")
    
    a$plotOptions(
      scatter = list(
        cursor = "pointer",
        point = list(
          events = list(
            click = "#! function(){$('input#selected2').val(this.options.bp); $('input#selected2').trigger('change');} !#")),
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
            click = "#! function(){$('input#selected2').val(this.options.bp); $('input#selected2').trigger('change');} !#")),
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
    a$set(dom = 'pChart2')
    return(a)
  })#end pchart2
  
  #Genome wide chart
  output$gChart <- renderChart({
       
    calculateTotalBP()

    #subset whole chart based on selection
    genomeChart <- values[[input$datasets]]
    if(input$plotAll == FALSE){
      for(i in input$traitColumns){
        genomeChart <- genomeChart[genomeChart[,i] %in% input[[i]],]
      }    
      if(length(input$traitColumns) > 1){
        genomeChart$trait <- do.call(paste,c(genomeChart[,input$traitColumns],sep="_"))    
      }else{
        genomeChart$trait <- genomeChart[,input$traitColumns]
      }
    }else{
        genomeChart$trait <- input$datasets
    }
    
    #Separate Support Interval data from GWAS data, if support, GWAS data is assumed to be anything that has an NA in the SIbpStart column
    if(input$supportInterval == TRUE){
      SIchart <- genomeChart[!(is.na(genomeChart[,input$SIbpStart])),]
      genomeChart <- genomeChart[is.na(genomeChart[,input$SIbpStart]),]
    }
    
    #filter genomeChart for only rows that have a base pair and yaxis value
    genomeChart <- genomeChart[!(is.na(genomeChart[,input$bpColumn])),]
    genomeChart <- genomeChart[!(is.na(genomeChart[,input$yAxisColumn])),]
    
    #if checked, filter for only overlapping SNPs
    if(!is.null(input$overlaps) & input$overlaps == TRUE){
      genomeChart <- findGWASOverlaps(genomeChart)
    }
    
    #check if there is any data for the selected traits
    if(nrow(genomeChart)==0){ #nothing is in the window, but lets still make a data.frame
       genomeChart <- values[[input$datasets]][1,]
       genomeChart[,input$yAxisColumn] <- -1    
       if(length(input$traitColumns) > 1){
         genomeChart$trait <- do.call(paste,c(genomeChart[,input$traitColumns],sep="_"))    
       }else{
         genomeChart$trait <- genomeChart[,input$traitColumns]
       }             
    }
    
    #take -log10 of y-axis column if requested
    if(input$logP == TRUE && genomeChart[1,input$yAxisColumn] != -1){
      genomeChart[,input$yAxisColumn] <- -log(genomeChart[,input$yAxisColumn],10)
    }
    
    #check if there is too much data (>2500 data points), trim to 2500
    if(nrow(genomeChart)>2500){
      cutVal <- sort(genomeChart[,input$yAxisColumn],decreasing = T)[2500]
      genomeChart <- genomeChart[genomeChart[,input$yAxisColumn] >= cutVal,]
    }
    
    colorTable <- colorTable()
     genomeTable <- data.frame(x=genomeChart$totalBP,y=genomeChart[,input$yAxisColumn],trait=genomeChart$trait,
#                               name=sprintf("<table cellpadding='4' style='line-height:1.5'><tr><th>%1$s</th></tr><tr><td align='left'>RMIP: %2$s<br>Location: %3$s<br>Base Pairs: %4$s<br>SNP: %5$s<br>Chromosome: %6$s</td></tr></table>",
#                               name=sprintf("<table cellpadding='4' style='line-height:1.5'><tr><th>%1$s</th></tr><tr><td align='left'>Y-value: %2$s<br>Base Pairs: %3$s<br>Chromosome: %4$s</td></tr></table>",
                                name=sprintf("Base Pair: %1$s<br/>Chromosome: %2$s<br/>",
#                                            genomeChart$trait,
#                                            genomeChart[,input$yAxisColumn],
#                                            genomeChart$loc,
                                            prettyNum(genomeChart[,input$bpColumn], big.mark = ","),
#                                            genomeChart$SNP,
                                            genomeChart[,input$chrColumn]
                               ),
                               url="http://danforthcenter.org",
                               chr=genomeChart[,input$chrColumn],
                               bp=genomeChart[,input$bpColumn],stringsAsFactors=FALSE)
     genomeSeries <- lapply(split(genomeTable, genomeTable$trait), function(x) {
       res <- lapply(split(x, rownames(x)), as.list)
       names(res) <- NULL
       res <- res[order(sapply(res, function(x) x$x))]
       return(res)
     })
#     
      #build JL series
    if(input$supportInterval==TRUE){
      if(nrow(SIchart)==0){ #nothing is in the window, but lets still make a data.frame
        SIchart <- values[[input$datasets]][1,]
        SIchart[,input$SIyAxisColumn] <- -1    
        if(length(input$traitColumns) > 1){
          SIchart$trait <- do.call(paste,c(SIchart[,input$traitColumns],sep="_"))    
        }else{
          SIchart$trait <- SIchart[,input$traitColumns]
        }             
      }
      SIchart$loc_el <- SIchart$trait
      SIchart$trait <- paste(SIchart$trait,"Int",sep="_")
      SIchart <- SIchart[order(SIchart$SIbpStartTotal),]
       jlTable <- adply(SIchart,1,function(x) {data.frame(x=c(x$SIbpStartTotal,x$SIbpEndTotal,x$SIbpEndTotal),y=c(x[[input$SIyAxisColumn]],x[[input$SIyAxisColumn]],NA),trait=x$trait,
                                                    #name=sprintf("<table cellpadding='4' style='line-height:1.5'><tr><th>%1$s</th></tr><tr><td align='left'>Y-value: %2$.2f <br>Interval: %3$s-%4$s<br>Chromosome: %5$s</td></tr></table>",
                                                    name=sprintf("<table cellpadding='4' style='line-height:1.5'><tr><td align='left'>Interval: %1$s-%2$s<br>Chromosome: %3$s</td></tr></table>",
#                                                                 x$trait,
#                                                                 x[[input$SIyAxisColumn]],
                                                                 prettyNum(x[[input$SIbpStart]], big.mark = ","),
                                                                 prettyNum(x[[input$SIbpEnd]], big.mark = ","),
                                                                 x[[input$chrColumn]]
                                                    ),loc_el=x$loc_el,bp=x[[input$bpColumn]],chr=x[[input$chrColumn]],stringsAsFactors=FALSE
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
    cumBP<-c(0,cumsum(as.numeric(chrSize[values$organism][[1]])))
    for(i in 1:(length(cumBP)-1)){
      if(i %% 2 == 0 ){ #even
        bigList[[length(bigList)+1]] <- list(from=cumBP[i]+1,to=cumBP[i+1],label=list(text=chrName[values$organism][[1]][i],style=list(color="#6D869F"),verticalAlign="bottom"))
      }else{ #odd
        bigList[[length(bigList)+1]] <- list(from=cumBP[i]+1,to=cumBP[i+1],color='rgba(68, 170, 213, 0.1)',label=list(text=chrName[values$organism][[1]][i],style=list(color="#6D869F"),verticalAlign="bottom"))
      }
    }    
    
     c <- rCharts::Highcharts$new()
    c$LIB$url <- 'highcharts/'
    c$xAxis(title = list(text = "Chromosome",margin=15),startOnTick=TRUE,min=0,max=sum(as.numeric(chrSize[values$organism][[1]])),endOnTick=FALSE,labels=list(enabled=FALSE),tickWidth=0,
            plotBands = bigList)   
    
     if(input$axisLimBool == TRUE){       
       c$yAxis(title=list(text=input$yAxisColumn),min=input$axisMin,max=input$axisMax,startOnTick=FALSE)
     }else{
       c$yAxis(title=list(text=input$yAxisColumn),startOnTick=FALSE)      
     }
     
     if(input$supportInterval==TRUE){
       if(input$SIaxisLimBool == TRUE){
         c$yAxis(title=list(text=input$SIyAxisColumn),min=input$SIaxisMin,max=input$SIaxisMax,gridLineWidth=0,minorGridLineWidth=0,startOnTick=FALSE,opposite=TRUE,replace=FALSE)   
       }else{
         c$yAxis(title=list(text=input$SIyAxisColumn),gridLineWidth=0,minorGridLineWidth=0,startOnTick=FALSE,opposite=TRUE,replace=FALSE)
       }
       
       if(SIchart[1,input$SIyAxisColumn] != -1){
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
     if(genomeChart[1,input$yAxisColumn] != -1){
       invisible(sapply(genomeSeries, function(x) {if(length(x)==0){return()};c$series(data = x, turboThreshold=5000,type = "scatter", color = colorTable$color[colorTable$trait == as.character(x[[1]]$trait)], name = paste0(x[[1]]$trait))}))
     }
     
     c$chart(zoomType="x",alignTicks=FALSE,events=list(click = "#!function(event) {this.tooltip.hide();}!#"))
     c$title(text=paste(input$datasets," Results",sep=" "))
     c$subtitle(text="Rollover for more info. Drag chart area to zoom. Click point to switch to chromosome and annotation view.")
     
     c$plotOptions(
        scatter = list(
          cursor = "pointer",
          point = list(
            events = list(
              #click = "#! function() { window.open(this.options.url); } !#")), #open webpage
              #click = "#! function(event) {alert(this.name);} !#")), #display popup
              #click = "#! function(event) {console.log(this);} !#")), #write object to log
              #click = "#! function(){$('input#selected').val(134); $('input#selected').trigger('change');} !#")),
              click = "#! function(){$('select#chr').val(this.options.chr); $('select#chr').trigger('change'); $('input#selected').val(this.options.bp); 
                                     $('input#selected').trigger('change'); $('ul#datatabs li').eq(0).removeClass('active'); 
                                     $('ul#datatabs li').eq(1).removeClass('active'); $('ul#datatabs li').eq(2).removeClass('active');
                                     $('ul#datatabs li').eq(4).removeClass('active');
                                     $('ul#datatabs li').eq(3).addClass('active'); 
                                     $('#pChart').trigger('change');$('#pChart').trigger('shown');
                                     $('.tab-content div').toggleClass(function(){if(this.getAttribute('data-value')=='panel2' || this.getAttribute('data-value')=='panel1'){return 'active';}else{return '';}});
                                     $('.tab-content div').trigger('change');$('ul#datatabs li').trigger('change');} !#")), 
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
             click = "#! function(){$('select#chr').val(this.options.chr); $('select#chr').trigger('change'); $('input#selected').val(this.options.bp); 
                                    $('input#selected').trigger('change'); $('ul#datatabs li').eq(0).removeClass('active'); 
                                    $('ul#datatabs li').eq(1).removeClass('active'); $('ul#datatabs li').eq(2).removeClass('active');
                                    $('ul#datatabs li').eq(4).removeClass('active');
                                    $('ul#datatabs li').eq(3).addClass('active'); 
                                    $('#pChart').trigger('change');$('#pChart').trigger('shown');
                                    $('.tab-content div').toggleClass(function(){if(this.getAttribute('data-value')=='panel2' || this.getAttribute('data-value')=='panel1'){return 'active';}else{return '';}});
                                    $('.tab-content div').trigger('change');$('ul#datatabs li').trigger('change');} !#")),             
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
     if(!is.null(input$legend) & input$legend == TRUE){
       c$legend(enabled=FALSE)
     }
     
     c$credits(enabled=TRUE)
     c$set(dom = 'gChart')     
     return(c)
#    h1 <- hPlot(x = "Wr.Hnd", y = "NW.Hnd", data = MASS::survey, type = c("line", 
#                                                                          "bubble", "scatter"), group = "Clap", size = "Age")
#    h1$set(dom = 'gChart')
#    return(h1)    
   })#end gchart
  output$gChart2 <- renderChart({
    calculateTotalBP2()
    
    #subset whole chart based on selection
    genomeChart <- values[[input$datasets2]]
    if(input$plotAll2 == FALSE){
      for(i in input$traitColumns2){
        i2 <- paste0(i, "2")
        genomeChart <- genomeChart[genomeChart[,i] %in% input[[i2]],]
      }    
      if(length(input$traitColumns2) > 1){
        genomeChart$trait <- do.call(paste,c(genomeChart[,input$traitColumns2],sep="_"))
      }else{
        genomeChart$trait <- genomeChart[,input$traitColumns2]
      }
    }else{
      genomeChart$trait <- input$datasets2
    }
    
    #Separate Support Interval data from GWAS data, if support, GWAS data is assumed to be anything that has an NA in the SIbpStart column
    if(input$supportInterval2 == TRUE){
      SIchart <- genomeChart[!(is.na(genomeChart[,input$SIbpStart2])),]
      genomeChart <- genomeChart[is.na(genomeChart[,input$SIbpStart2]),]
    }
    
    #filter genomeChart for only rows that have a base pair and yaxis value
    genomeChart <- genomeChart[!(is.na(genomeChart[,input$bpColumn2])),]
    genomeChart <- genomeChart[!(is.na(genomeChart[,input$yAxisColumn2])),]
    
    #if checked, filter for only overlapping SNPs
    if(!is.null(input$overlaps2) & input$overlaps2 == TRUE){
      genomeChart <- findGWASOverlaps2(genomeChart)
    }
    
    #check if there is any data for the selected traits
    if(nrow(genomeChart)==0){ #nothing is in the window, but lets still make a data.frame
      genomeChart <- values[[input$datasets2]][1,]
      genomeChart[,input$yAxisColumn2] <- -1    
      if(length(input$traitColumns2) > 1){
        genomeChart$trait <- do.call(paste,c(genomeChart[,input$traitColumns2],sep="_"))    
      }else{
        genomeChart$trait <- genomeChart[,input$traitColumns2]
      }             
    }
    
    #take -log10 of y-axis column if requested
    if(input$logP2 == TRUE && genomeChart[1,input$yAxisColumn2] != -1){
      genomeChart[,input$yAxisColumn2] <- -log(genomeChart[,input$yAxisColumn2],10)
    }
    
    #check if there is too much data (>2500 data points), trim to 2500
    if(nrow(genomeChart)>2500){
      cutVal <- sort(genomeChart[,input$yAxisColumn2],decreasing = T)[2500]
      genomeChart <- genomeChart[genomeChart[,input$yAxisColumn2] >= cutVal,]
    }
    
    colorTable <- colorTable2()
    genomeTable <- data.frame(x=genomeChart$totalBP,y=genomeChart[,input$yAxisColumn2],trait=genomeChart$trait,
      name=sprintf("Base Pair: %1$s<br/>Chromosome: %2$s<br/>",
        prettyNum(genomeChart[,input$bpColumn2], big.mark = ","),
        genomeChart[,input$chrColumn2]
      ),
      url="http://danforthcenter.org",
      chr=genomeChart[,input$chrColumn2],
      bp=genomeChart[,input$bpColumn2],stringsAsFactors=FALSE)
    genomeSeries <- lapply(split(genomeTable, genomeTable$trait), function(x) {
      res <- lapply(split(x, rownames(x)), as.list)
      names(res) <- NULL
      res <- res[order(sapply(res, function(x) x$x))]
      return(res)
    })
    #     
    #build JL series
    if(input$supportInterval2==TRUE){
      if(nrow(SIchart)==0){ #nothing is in the window, but lets still make a data.frame
        SIchart <- values[[input$datasets2]][1,]
        SIchart[,input$SIyAxisColumn2] <- -1    
        if(length(input$traitColumns2) > 1){
          SIchart$trait <- do.call(paste,c(SIchart[,input$traitColumns2],sep="_"))    
        }else{
          SIchart$trait <- SIchart[,input$traitColumns2]
        }             
      }
      SIchart$loc_el <- SIchart$trait
      SIchart$trait <- paste(SIchart$trait,"Int",sep="_")
      SIchart <- SIchart[order(SIchart$SIbpStartTotal),]
      jlTable <- adply(SIchart,1,function(x) {data.frame(x=c(x$SIbpStartTotal,x$SIbpEndTotal,x$SIbpEndTotal),y=c(x[[input$SIyAxisColumn2]],x[[input$SIyAxisColumn2]],NA),trait=x$trait,
        name=sprintf("<table cellpadding='4' style='line-height:1.5'><tr><td align='left'>Interval: %1$s-%2$s<br>Chromosome: %3$s</td></tr></table>",
          prettyNum(x[[input$SIbpStart2]], big.mark = ","),
          prettyNum(x[[input$SIbpEnd2]], big.mark = ","),
          x[[input$chrColumn2]]
        ),loc_el=x$loc_el,bp=x[[input$bpColumn2]],chr=x[[input$chrColumn2]],stringsAsFactors=FALSE
      )}#end jlTable and function
      )#end adply
      jlTable <- jlTable[,c("x","y","trait","name","loc_el","bp","chr")]
      #jlTable <- jlTable[order(jlTable$x),]
    }#end build jlTable if support intervals
    
    #build list for where to put plotbands for this organism
    bigList <- list()
    cumBP<-c(0,cumsum(as.numeric(chrSize[values$organism2][[1]])))
    for(i in 1:(length(cumBP)-1)){
      if(i %% 2 == 0 ){ #even
        bigList[[length(bigList)+1]] <- list(from=cumBP[i]+1,to=cumBP[i+1],label=list(text=chrName[values$organism2][[1]][i],style=list(color="#6D869F"),verticalAlign="bottom"))
      }else{ #odd
        bigList[[length(bigList)+1]] <- list(from=cumBP[i]+1,to=cumBP[i+1],color='rgba(68, 170, 213, 0.1)',label=list(text=chrName[values$organism2][[1]][i],style=list(color="#6D869F"),verticalAlign="bottom"))
      }
    }    
    
    c <- rCharts::Highcharts$new()
    c$LIB$url <- 'highcharts/'
    c$xAxis(title = list(text = "Chromosome",margin=15),startOnTick=TRUE,min=0,max=sum(as.numeric(chrSize[values$organism2][[1]])),endOnTick=FALSE,labels=list(enabled=FALSE),tickWidth=0,
      plotBands = bigList)   
    
    if(input$axisLimBool2 == TRUE){       
      c$yAxis(title=list(text=input$yAxisColumn2),min=input$axisMin2,max=input$axisMax2,startOnTick=FALSE)
    }else{
      c$yAxis(title=list(text=input$yAxisColumn2),startOnTick=FALSE)      
    }
    
    if(input$supportInterval2==TRUE){
      if(input$SIaxisLimBool2 == TRUE){
        c$yAxis(title=list(text=input$SIyAxisColumn2),min=input$SIaxisMin2,max=input$SIaxisMax2,gridLineWidth=0,minorGridLineWidth=0,startOnTick=FALSE,opposite=TRUE,replace=FALSE)   
      }else{
        c$yAxis(title=list(text=input$SIyAxisColumn2),gridLineWidth=0,minorGridLineWidth=0,startOnTick=FALSE,opposite=TRUE,replace=FALSE)
      }
      
      if(SIchart[1,input$SIyAxisColumn2] != -1){
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
    if(genomeChart[1,input$yAxisColumn2] != -1){
      invisible(sapply(genomeSeries, function(x) {if(length(x)==0){return()};c$series(data = x, turboThreshold=5000,type = "scatter", color = colorTable$color[colorTable$trait == as.character(x[[1]]$trait)], name = paste0(x[[1]]$trait))}))
    }
    
    c$chart(zoomType="x",alignTicks=FALSE,events=list(click = "#!function(event) {this.tooltip.hide();}!#"))
    c$title(text=paste(input$datasets2," Results",sep=" "))
    c$subtitle(text="Rollover for more info. Drag chart area to zoom. Click point to switch to chromosome and annotation view.")
    
    c$plotOptions(
      scatter = list(
        cursor = "pointer",
        point = list(
          events = list(
            click = "#! function(){$('select#chr2').val(this.options.chr); $('select#chr2').trigger('change'); $('input#selected2').val(this.options.bp); 
            $('input#selected2').trigger('change'); $('ul#datatabs li').eq(0).removeClass('active'); 
            $('ul#datatabs li').eq(1).removeClass('active'); $('ul#datatabs li').eq(2).removeClass('active');
            $('ul#datatabs li').eq(4).removeClass('active');
            $('ul#datatabs li').eq(3).addClass('active'); 
            $('#pChart2').trigger('change');$('#pChart2').trigger('shown');
            $('.tab-content div').toggleClass(function(){if(this.getAttribute('data-value')=='panel2' || this.getAttribute('data-value')=='panel1'){return 'active';}else{return '';}});
            $('.tab-content div').trigger('change');$('ul#datatabs li').trigger('change');} !#")), 
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
            click = "#! function(){$('select#chr2').val(this.options.chr); $('select#chr2').trigger('change'); $('input#selected2').val(this.options.bp); 
            $('input#selected2').trigger('change'); $('ul#datatabs li').eq(0).removeClass('active'); 
            $('ul#datatabs li').eq(1).removeClass('active'); $('ul#datatabs li').eq(2).removeClass('active');
            $('ul#datatabs li').eq(4).removeClass('active');
            $('ul#datatabs li').eq(3).addClass('active'); 
            $('#pChart2').trigger('change');$('#pChart2').trigger('shown');
            $('.tab-content div').toggleClass(function(){if(this.getAttribute('data-value')=='panel2' || this.getAttribute('data-value')=='panel1'){return 'active';}else{return '';}});
            $('.tab-content div').trigger('change');$('ul#datatabs li').trigger('change');} !#")),             
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
    if(!is.null(input$legend2) & input$legend2 == TRUE){
      c$legend(enabled=FALSE)
    }
    
    c$credits(enabled=TRUE)
    c$set(dom = 'gChart2')     
    return(c)
  })#end gchart2

  output$zChart <- renderChart({
    if(is.null(input$selected)) return()

    centerBP <- as.numeric(input$selected[[1]])
    winHigh <- centerBP+input$window[1]
    winLow <- centerBP-input$window[1]
    if(winLow < 0){winLow <- 0}
    
    zoomChart <- values[[input$datasets]]
    zoomChart <- zoomChart[zoomChart[,input$chrColumn]==input$chr,]    
    
    if(input$plotAll == FALSE){
      for(i in input$traitColumns){
        zoomChart <- zoomChart[zoomChart[,i] %in% input[[i]],]
      }
      
      if(length(input$traitColumns) > 1){
        zoomChart$trait <- do.call(paste,c(zoomChart[,input$traitColumns],sep="_"))    
      }else{
        zoomChart$trait <- zoomChart[,input$traitColumns]
      }
    }else{
      zoomChart$trait <- input$datasets
    }
    
    #Separate Support Interval data from GWAS data, if support, GWAS data is assumed to be anything that has an NA in the SIbpStart column
    if(input$supportInterval == TRUE){
      SIchart <- zoomChart[!(is.na(zoomChart[,input$SIbpStart])),]
      zoomChart <- zoomChart[is.na(zoomChart[,input$SIbpStart]),]
      #not sure the below logic works for subsetting SIchart, probably not necessary anyways, since there are usually very few SI rows for one chromosome anyways (e.g. small overhead)
      #SIchart <- SIchart[((SIchart[,input$SIbpStart] <= winHigh & SIchart[,input$SIbpStart] >= winLow) | (SIchart[,input$SIbpEnd] <= winHigh & SIchart[,input$SIbpEnd] >= winLow)),]
    }    
    
    zoomChart <- zoomChart[(zoomChart[,input$bpColumn] <= winHigh) & (zoomChart[,input$bpColumn] >= winLow),]    
    
    #filter for only rows that have a base pair value
    zoomChart <- zoomChart[!(is.na(zoomChart[,input$bpColumn])),]
    zoomChart <- zoomChart[!(is.na(zoomChart[,input$yAxisColumn])),]    
    
    #if checked, filter for only overlapping SNPs
    if(!is.null(input$overlaps) & input$overlaps == TRUE){
      zoomChart <- findGWASOverlaps(zoomChart)
    }                    
    
    if(nrow(zoomChart)==0){ #nothing is in the window, but lets still make a data.frame
      zoomChart <- values[[input$datasets]][1,]
      zoomChart[,input$yAxisColumn] <- -1    
      if(length(input$traitColumns) > 1){
        zoomChart$trait <- do.call(paste,c(zoomChart[,input$traitColumns],sep="_"))    
      }else{
        zoomChart$trait <- zoomChart[,input$traitColumns]
      }                   
    }
    colorTable <- colorTable() 

    #take -log10 of y-axis column if requested
    if(input$logP == TRUE && zoomChart[1,input$yAxisColumn] != -1){
      zoomChart[,input$yAxisColumn] <- -log(zoomChart[,input$yAxisColumn],10)
    }                
    
    #check if there is too much data (>2500 data points), trim to 2500
    if(nrow(zoomChart)>2500){
      cutVal <- sort(zoomChart[,input$yAxisColumn],decreasing = T)[2500]
      zoomChart <- zoomChart[zoomChart[,input$yAxisColumn] >= cutVal,]
    }                
    
    zoomTable <- data.frame(x=zoomChart[,input$bpColumn],y=zoomChart[,input$yAxisColumn],trait=zoomChart$trait,
#                                         name=sprintf("<table cellpadding='4' style='line-height:1.5'><tr><th>%1$s</th></tr><tr><td align='left'>RMIP: %2$s<br>Location: %3$s<br>Base Pairs: %4$s<br>SNP: %5$s<br>Chromosome: %6$s</td></tr></table>",
                                         name=sprintf("<table cellpadding='4' style='line-height:1.5'><tr><th>%1$s</th></tr><tr><td align='left'>Y-axis value: %2$s<br>Base Pairs: %3$s<br>Chromosome: %4$s</td></tr></table>",                                         
                                         zoomChart$trait,
                                         zoomChart[,input$yAxisColumn],
                                         #zoomChart$loc,
                                         prettyNum(zoomChart[,input$bpColumn], big.mark = ","),
                                         #zoomChart$SNP,
                                         zoomChart[,input$chrColumn]
                            ),
                            url="http://danforthcenter.org",
                            chr=zoomChart[,input$chrColumn],
                            bp=zoomChart[,input$bpColumn])
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
     if(input$supportInterval==TRUE){
       if(nrow(SIchart) == 0){ #make a dummy table, but we won't plot the series anyways
         SIchart <- values[[input$datasets]][1,]
         SIchart[,input$SIyAxisColumn] <- -1    
         if(length(input$traitColumns) > 1){
           SIchart$trait <- do.call(paste,c(SIchart[,input$traitColumns],sep="_"))    
         }else{
           SIchart$trait <- SIchart[,input$traitColumns]
         }                      
       }     
       SIchart$loc_el <- SIchart$trait
       SIchart$trait <- paste(SIchart$trait,"Int",sep="_")
       SIchart <- SIchart[order(SIchart[[input$SIbpStart]]),]
       jlTable <- adply(SIchart,1,function(x) {data.frame(x=c(x[[input$SIbpStart]],x[[input$SIbpEnd]],x[[input$SIbpEnd]]),y=c(x[[input$SIyAxisColumn]],x[[input$SIyAxisColumn]],NA),trait=x$trait,
                                                          name=sprintf("<table cellpadding='4' style='line-height:1.5'><tr><td align='left'>Interval: %1$s-%2$s<br>Chromosome: %3$s</td></tr></table>",
                                                                       #                                                                      x$trait,
                                                                       #                                                                      x[[input$SIyAxisColumn]],
                                                                       prettyNum(x[[input$SIbpStart]], big.mark = ","),
                                                                       prettyNum(x[[input$SIbpEnd]], big.mark = ","),
                                                                       x[[input$chrColumn]]
                                                          ),loc_el=x$loc_el,bp=x[[input$bpColumn]],chr=x[[input$chrColumn]],stringsAsFactors=FALSE
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
    #thisChrAnnot <- subset(annotGeneLoc,chromosome==input$chr)
    thisChrAnnot <- subset(annotGeneLoc[values$organism][[1]],chromosome==input$chr)
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
    #if(input$axisLimBool == TRUE){annotYvalReverse <- input$axisMin+0.01}
    annotYvalForward <- annotYvalReverse + 0.04
    if(values$organism == "Corn"){
      annotTable <- adply(thisAnnot[thisAnnot$transcript_strand==1,],1,function(x) {data.frame(x=c(x$transcript_start,x$transcript_end,x$transcript_end),y=c(annotYvalForward,annotYvalForward,NA),url=paste0(urlBase,x$transcript_id),
                                                              name=sprintf("<table cellpadding='4' style='line-height:1.5'><tr><th>%1$s</th></tr><tr><td align='left'>Location: %2$s-%3$s<br>Chromosome: %4$s, Strand: %5$s<br>%6$s</td></tr></table>",
                                                                           x$translation_id,
                                                                           prettyNum(x$transcript_start, big.mark = ","),
                                                                           prettyNum(x$transcript_end, big.mark = ","),
                                                                           x$chromosome,
                                                                           x$transcript_strand,
                                                                           x$V2
                                                              ),
                                                              marker=c(NA,"Arrow",NA),
                                                              stringsAsFactors=FALSE)})
      
      annotTableReverse <- adply(thisAnnot[thisAnnot$transcript_strand==-1,],1,function(x) {data.frame(x=c(x$transcript_start,x$transcript_end,x$transcript_end),y=c(annotYvalReverse,annotYvalReverse,NA),url=paste0(urlBase,x$transcript_id),
                                                              name=sprintf("<table cellpadding='4' style='line-height:1.5'><tr><th>%1$s</th></tr><tr><td align='left'>Location: %2$s-%3$s<br>Chromosome: %4$s, Strand: %5$s<br>%6$s</td></tr></table>",
                                                                           x$translation_id,
                                                                           prettyNum(x$transcript_start, big.mark = ","),
                                                                           prettyNum(x$transcript_end, big.mark = ","),
                                                                           x$chromosome,
                                                                           x$transcript_strand,
                                                                           x$V2
                                                              ),
                                                              marker=c("Arrow",NA,NA),
                                                              stringsAsFactors=FALSE)})
            
    }else if(values$organism == "Soybean"){#strand is '+' or '-'
      annotTable <- adply(thisAnnot[thisAnnot$strand=="+",],1,function(x) {data.frame(x=c(x$transcript_start,x$transcript_end,x$transcript_end),y=c(annotYvalForward,annotYvalForward,NA),url=paste0(soyurlBase,x$transcript_id),
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
                                                              stringsAsFactors=FALSE)})
      
      annotTableReverse <- adply(thisAnnot[thisAnnot$strand=="-",],1,function(x) {data.frame(x=c(x$transcript_start,x$transcript_end,x$transcript_end),y=c(annotYvalReverse,annotYvalReverse,NA),url=paste0(soyurlBase,x$transcript_id),
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
                                                                                     stringsAsFactors=FALSE)})
    }else if(values$organism %in% c("Arabidopsis", "Arabidopsis thaliana")){#strand is '+' or '-'
      annotTable <- adply(thisAnnot[thisAnnot$strand=="+",],1,function(x) {data.frame(x=c(x$transcript_start,x$transcript_end,x$transcript_end),y=c(annotYvalForward,annotYvalForward,NA),url=paste0(araburlBase,x$Locus),
                                                              name=sprintf("<table cellpadding='4' style='line-height:1.5'><tr><th>%1$s</th></tr><tr><td align='left'>Location: %2$s-%3$s<br>Chromosome: %4$s, Strand: %5$s<br>Short Desc.: %6$s</td></tr></table>",
                                                                           x$name,
                                                                           prettyNum(x$transcript_start, big.mark = ","),
                                                                           prettyNum(x$transcript_end, big.mark = ","),
                                                                           x$chromosome,                                                                           
                                                                           x$strand,
                                                                           x$short_description
                                                                          # x$Curator_summary
                                                              ),
                                                              stringsAsFactors=FALSE)})   
      
      annotTableReverse <- adply(thisAnnot[thisAnnot$strand=="-",],1,function(x) {data.frame(x=c(x$transcript_start,x$transcript_end,x$transcript_end),y=c(annotYvalReverse,annotYvalReverse,NA),url=paste0(araburlBase,x$Locus),
                                                                                     name=sprintf("<table cellpadding='4' style='line-height:1.5'><tr><th>%1$s</th></tr><tr><td align='left'>Location: %2$s-%3$s<br>Chromosome: %4$s, Strand: %5$s<br>Short Desc.: %6$s</td></tr></table>",
                                                                                                  x$name,
                                                                                                  prettyNum(x$transcript_start, big.mark = ","),
                                                                                                  prettyNum(x$transcript_end, big.mark = ","),
                                                                                                  x$chromosome,                                                                           
                                                                                                  x$strand,
                                                                                                  x$short_description
                                                                                                  # x$Curator_summary
                                                                                     ),
                                                                                     stringsAsFactors=FALSE)})
    } else if (values$organism == "Medicago truncatula") { # strand is '+' or '-'
      annotTable <- adply(thisAnnot[thisAnnot$strand=="+",],1,function(x) {data.frame(x=c(x$transcript_start,x$transcript_end,x$transcript_end),y=c(annotYvalForward,annotYvalForward,NA),url=paste0(legumeInfo_urlBase, x$name, "/json"),
        name=sprintf("<table cellpadding='4' style='line-height:1.5'><tr><th>%1$s</th></tr><tr><td align='left'>Location: %2$s-%3$s<br>Chromosome: %4$s, Strand: %5$s<br>Desc: %6$s</td></tr></table>",
          x$name,
          prettyNum(x$transcript_start, big.mark = ","),
          prettyNum(x$transcript_end, big.mark = ","),
          x$chromosome,
          x$strand,
          x$description
        ),
        stringsAsFactors=FALSE)})

      annotTableReverse <- adply(thisAnnot[thisAnnot$strand=="-",],1,function(x) {data.frame(x=c(x$transcript_start,x$transcript_end,x$transcript_end),y=c(annotYvalReverse,annotYvalReverse,NA),url=paste0(legumeInfo_urlBase, x$name, "/json"),
        name=sprintf("<table cellpadding='4' style='line-height:1.5'><tr><th>%1$s</th></tr><tr><td align='left'>Location: %2$s-%3$s<br>Chromosome: %4$s, Strand: %5$s<br>Desc: %6$s</td></tr></table>",
          x$name,
          prettyNum(x$transcript_start, big.mark = ","),
          prettyNum(x$transcript_end, big.mark = ","),
          x$chromosome,
          x$strand,
          x$description
        ),
        stringsAsFactors=FALSE)})
    }else{#} if(input$organism == "Sorghum"){#strand is '+' or '-'
      annotTable <- adply(thisAnnot[thisAnnot$strand=="+",],1,function(x) {data.frame(x=c(x$transcript_start,x$transcript_end,x$transcript_end),y=c(annotYvalForward,annotYvalForward,NA),url=paste0(sorgurlBase,x$ID),
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
                                                                                      stringsAsFactors=FALSE)})   
      
      annotTableReverse <- adply(thisAnnot[thisAnnot$strand=="-",],1,function(x) {data.frame(x=c(x$transcript_start,x$transcript_end,x$transcript_end),y=c(annotYvalReverse,annotYvalReverse,NA),url=paste0(sorgurlBase,x$ID),
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
                                                                                             stringsAsFactors=FALSE)})
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
    if(input$axisLimBool == TRUE){
      b$yAxis(title=list(text=input$yAxisColumn),min=input$axisMin,max=input$axisMax,startOnTick=FALSE)
      #create a hidden axis to put the gene track on, all the options are setting to hide everything from the axis 
      b$yAxis(labels=list(enabled=FALSE),title=list(text=NULL),min=0,max=1,lineWidth=0,gridLineWidth=0,minorGridLineWidth=0,lineColor="transparent",minorTickLength=0,tickLength=0,startOnTick=FALSE,opposite=TRUE,replace=FALSE)
    }else{      
      b$yAxis(title=list(text=input$yAxisColumn),startOnTick=FALSE) 
      #create a hidden axis to put the gene track on, all the options are setting to hide everything from the axis
      b$yAxis(labels=list(enabled=FALSE),title=list(text=NULL),min=0,max=1,lineWidth=0,gridLineWidth=0,minorGridLineWidth=0,lineColor="transparent",minorTickLength=0,tickLength=0,startOnTick=FALSE,opposite=TRUE,replace=FALSE)
    }

    if(input$supportInterval==TRUE){
      if(input$SIaxisLimBool == TRUE){
        b$yAxis(title=list(text=input$SIyAxisColumn),min=input$SIaxisMin,max=input$SIaxisMax,gridLineWidth=0,minorGridLineWidth=0,startOnTick=FALSE,opposite=TRUE,replace=FALSE)
      }else{
        b$yAxis(title=list(text=input$SIyAxisColumn),gridLineWidth=0,minorGridLineWidth=0,startOnTick=FALSE,opposite=TRUE,replace=FALSE)   
      }
      
      if(SIchart[1,input$SIyAxisColumn] != -1){
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
    if(zoomChart[1,input$yAxisColumn] != -1){
      invisible(sapply(zoomSeries, function(x) {if(length(x)==0){return()};b$series(data = x, type = "scatter", color = colorTable$color[colorTable$trait == as.character(x[[1]]$trait)], name = paste0(x[[1]]$trait))}))
    }
    
    b$series(
      data = annotArray,
      type = "line",
      name = "Forward Genes",
      id = "forward-genes",
      zIndex = 1,
      color = "#53377A",
      yAxis = 1
    )    

    b$series(
      data = annotArrayReverse,
      type = "line",
      name = "Reverse Genes",
      id = "reverse-genes",
      zIndex = 1,
      color = "#53377A",
      yAxis = 1
    )      

    if (!is.null(values$glGenes1)) {
      apply(values$glGenes1, 1, FUN = function(g) {
        g.strand <- as.integer(g$strand)
        yh <- -1
        if (g.strand == 1) {
          yh <- annotYvalForward
          sid <- "forward-genes"
        } else if (g.strand == -1) {
          yh <- annotYvalReverse
          sid <- "reverse-genes"
        }
        if (yh > 0 && g$chr == input$chr) {
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
    b$plotOptions(
      scatter = list(
        cursor = "pointer",
        point = list(
          events = list(
            #click = "#! function() { window.open(this.options.url); } !#")), #open webpage
            click = "#! function(event) {alert(this.trait);} !#")), #display popup
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
            click = paste0(
              "#! function() {",
                "if ($('input#boolGenomicLinkage').prop('checked')) {",
                  "var url1 = '';
                  var url2 = '';
                  var speciesName2 = '';
                  var geneString = '';
                  mt0 = this.url.search('medtr');
                  mt1 = this.url.search('/json');
                  if (mt0 >= 0) {
                    geneString = this.url.substring(mt0, mt1);
                    url1 = 'legumeinfo.org';
                    url2 = 'legumefederation.org';
                    speciesName2 = 'A.thaliana';
                  } else if (this.url.search('arabidopsis.org') >= 0) {
                    at0 = this.url.search('name=');
                    geneString = 'arath.Col.' + this.url.substring(at0 + 5);
                    url1 = 'legumefederation.org';
                    url2 = 'legumeinfo.org';
                    speciesName2 = 'M.truncatula';
                  } else {
                    return;
                  }",
                  "$.ajax({
                    url: 'https://' + url1 + '/lis_context_server/services/v1/gene-to-query-track/',
                    dataType: 'json',
                    data: JSON.stringify({
                      gene: geneString,
                      neighbors: $('input#neighbors').val()
                    }),
                    type: 'POST',
                    success: function(response) {
                      obj1 = JSON.parse(response);
                      families1 = Array.from(obj1.genes, x => x.family);
                      $.ajax({
                        url: 'https://' + url2 + '/lis_context_server/services/v1/micro-synteny-search/',
                        dataType: 'json',
                        data: JSON.stringify({
                          query: families1,
                          matched: $('input#matched').val(),
                          intermediate: $('input#intermediate').val()
                        }),
                        type: 'POST',
                        success: function(response2) {
                          obj2 = JSON.parse(response2);
                          // Send information about neighboring and related genes back to the chart
                          Shiny.onInputChange('genomicLinkages', {
                            results1: obj1,
                            results2: obj2
                          });
                        },
                        error: function(errmsg2) { alert('FAIL2: ' + errmsg2.responseText); }
                      });
                    },
                    error: function(errmsg) { alert('FAIL: ' + errmsg.responseText); }
                  });",

                "} else {",
                  "if (this.url.includes('legumeinfo.org')) {",
                    # From the JSON at this.url, extract the URLs related to this gene.
                    # Note that this.url = legumeInfo_urlBase + geneString + '/json'
                    #  legumeInfo_urlBase currently has 34 characters (see above)
                    #  and geneString = <5-character species abbreviation>.geneName
                    # And for now, add the gene family phylogram URL by hand.
                    "$.getJSON(this.url, function(data) {
                      var geneString = this.url.substring(34, this.url.indexOf('/json'));
                      var geneName = geneString.substring(6);
                      var content = '';
                      if (data.length == 0) {
                        content = '<p>No ' + geneName + ' links found.</p>';
                      } else {
                        $.each(data, function(i, obj) {
                          content = content + '<p><a href=' + obj.href + ' target=_blank>' + obj.text + '</a></p>';
                          if (i == 0) {
                            var urlPhylogram = 'http://legumeinfo.org/chado_gene_phylotree_v2?gene_name=' + geneString;
                            var textPhylogram = 'View LIS gene family phylogram page for : ' + geneName;
                            content = content + '<p><a href=' + urlPhylogram + ' target=_blank>' + textPhylogram + '</a></p>';
                          }
                        });
                      }

                      var $div = $('<div></div>');
                      $div.html(content);
                      $div.dialog({
                        title: geneName + ' Links',
                        width: 512,
                        height: 'auto',
                        modal: true
                      });
                    });",

                  "} else {",
                    # for all other species
                    "window.open(this.url);", #open webpage
                  "}",
                "}",
              "} !#"
        ))),
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
    )        
    #it seems almost impossible to get the tooltip to hover along the chart with this version of highcharts (4.0.1), perhaps a question to stackoverflow could solve it.
    #see an example of the problem here: http://jsfiddle.net/N5ymb/
    #one hack/fix would be to add dummy points to the middle of the line that show up when moused over
    b$tooltip(snap=5, useHTML = T, formatter = "#! function() { return this.point.name; } !#") #followTouchMove = T, shared=T, followPointer = T
    b$exporting(enabled=TRUE,filename='zoomChart',sourceWidth=2000)
    b$credits(enabled=TRUE)
    b$set(dom = 'zChart')
    return(b)
  })#end zchart
  output$zChart2 <- renderChart({
    if(is.null(input$selected2)) return()
    
    centerBP <- as.numeric(input$selected2[[1]])
    winHigh <- centerBP+input$window2[1]
    winLow <- centerBP-input$window2[1]
    if(winLow < 0){winLow <- 0}
    
    zoomChart <- values[[input$datasets2]]
    zoomChart <- zoomChart[zoomChart[,input$chrColumn2]==input$chr2,]
    
    if(input$plotAll2 == FALSE){
      for(i in input$traitColumns2){
        i2 <- paste0(i, "2")
        zoomChart <- zoomChart[zoomChart[,i] %in% input[[i2]],]
      }
      
      if(length(input$traitColumns2) > 1){
        zoomChart$trait <- do.call(paste,c(zoomChart[,input$traitColumns2],sep="_"))
      }else{
        zoomChart$trait <- zoomChart[,input$traitColumns2]
      }
    }else{
      zoomChart$trait <- input$datasets2
    }
    
    #Separate Support Interval data from GWAS data, if support, GWAS data is assumed to be anything that has an NA in the SIbpStart column
    if(input$supportInterval2 == TRUE){
      SIchart <- zoomChart[!(is.na(zoomChart[,input$SIbpStart2])),]
      zoomChart <- zoomChart[is.na(zoomChart[,input$SIbpStart2]),]
      #not sure the below logic works for subsetting SIchart, probably not necessary anyways, since there are usually very few SI rows for one chromosome anyways (e.g. small overhead)
      #SIchart <- SIchart[((SIchart[,input$SIbpStart] <= winHigh & SIchart[,input$SIbpStart] >= winLow) | (SIchart[,input$SIbpEnd] <= winHigh & SIchart[,input$SIbpEnd] >= winLow)),]
    }    
    
    zoomChart <- zoomChart[(zoomChart[,input$bpColumn2] <= winHigh) & (zoomChart[,input$bpColumn2] >= winLow),]
    
    #filter for only rows that have a base pair value
    zoomChart <- zoomChart[!(is.na(zoomChart[,input$bpColumn2])),]
    zoomChart <- zoomChart[!(is.na(zoomChart[,input$yAxisColumn2])),]    
    
    #if checked, filter for only overlapping SNPs
    if(!is.null(input$overlaps2) & input$overlaps2 == TRUE){
      zoomChart <- findGWASOverlaps2(zoomChart)
    }
    
    if(nrow(zoomChart)==0){ #nothing is in the window, but lets still make a data.frame
      zoomChart <- values[[input$datasets2]][1,]
      zoomChart[,input$yAxisColumn2] <- -1    
      if(length(input$traitColumns2) > 1){
        zoomChart$trait <- do.call(paste,c(zoomChart[,input$traitColumns2],sep="_"))
      }else{
        zoomChart$trait <- zoomChart[,input$traitColumns2]
      }                   
    }
    colorTable <- colorTable2()
    
    #take -log10 of y-axis column if requested
    if(input$logP2 == TRUE && zoomChart[1,input$yAxisColumn2] != -1){
      zoomChart[,input$yAxisColumn2] <- -log(zoomChart[,input$yAxisColumn2],10)
    }
    
    #check if there is too much data (>2500 data points), trim to 2500
    if(nrow(zoomChart)>2500){
      cutVal <- sort(zoomChart[,input$yAxisColumn2],decreasing = T)[2500]
      zoomChart <- zoomChart[zoomChart[,input$yAxisColumn2] >= cutVal,]
    }                
    
    zoomTable <- data.frame(x=zoomChart[,input$bpColumn2],y=zoomChart[,input$yAxisColumn2],trait=zoomChart$trait,
      name=sprintf("<table cellpadding='4' style='line-height:1.5'><tr><th>%1$s</th></tr><tr><td align='left'>Y-axis value: %2$s<br>Base Pairs: %3$s<br>Chromosome: %4$s</td></tr></table>",
        zoomChart$trait,
        zoomChart[,input$yAxisColumn2],
        prettyNum(zoomChart[,input$bpColumn2], big.mark = ","),
        zoomChart[,input$chrColumn2]
      ),
      url="http://danforthcenter.org",
      chr=zoomChart[,input$chrColumn2],
      bp=zoomChart[,input$bpColumn2])
    zoomSeries <- lapply(split(zoomTable, zoomTable$trait), function(x) {
      res <- lapply(split(x, rownames(x)), as.list)
      names(res) <- NULL
      res <- res[order(sapply(res, function(x) x$x))]
      return(res)
    })
    
    # build JL series
    if(input$supportInterval2==TRUE){
      if(nrow(SIchart) == 0){ #make a dummy table, but we won't plot the series anyways
        SIchart <- values[[input$datasets2]][1,]
        SIchart[,input$SIyAxisColumn2] <- -1    
        if(length(input$traitColumns2) > 1){
          SIchart$trait <- do.call(paste,c(SIchart[,input$traitColumns2],sep="_"))
        }else{
          SIchart$trait <- SIchart[,input$traitColumns2]
        }                      
      }     
      SIchart$loc_el <- SIchart$trait
      SIchart$trait <- paste(SIchart$trait,"Int",sep="_")
      SIchart <- SIchart[order(SIchart[[input$SIbpStart2]]),]
      jlTable <- adply(SIchart,1,function(x) {data.frame(x=c(x[[input$SIbpStart2]],x[[input$SIbpEnd2]],x[[input$SIbpEnd2]]),y=c(x[[input$SIyAxisColumn2]],x[[input$SIyAxisColumn2]],NA),trait=x$trait,
        name=sprintf("<table cellpadding='4' style='line-height:1.5'><tr><td align='left'>Interval: %1$s-%2$s<br>Chromosome: %3$s</td></tr></table>",
          prettyNum(x[[input$SIbpStart2]], big.mark = ","),
          prettyNum(x[[input$SIbpEnd2]], big.mark = ","),
          x[[input$chrColumn2]]
        ),loc_el=x$loc_el,bp=x[[input$bpColumn2]],chr=x[[input$chrColumn2]],stringsAsFactors=FALSE
      )}#end jlTable and function
      )#end adply
      jlTable <- jlTable[,c("x","y","trait","name","loc_el","bp","chr")]
      #jlTable <- jlTable[order(jlTable$x),]
    }#end if support interval
    
    #build annotation series
    thisChrAnnot <- subset(annotGeneLoc[values$organism2][[1]],chromosome==input$chr2)
    thisAnnot <- thisChrAnnot[thisChrAnnot$transcript_start >= winLow & thisChrAnnot$transcript_end <= winHigh,]
    if(nrow(thisAnnot)==0){ #nothing is in the window, but lets still make a data.frame (actually make it big just to hopefully pick up one row from each strand...)
      thisAnnot <- thisChrAnnot[1:100,]
    }
    thisAnnot <- thisAnnot[order(thisAnnot$transcript_start),]
    
    urlBase <- 'http://maizegdb.org/cgi-bin/displaygenemodelrecord.cgi?id='
    soyurlBase <- 'http://www.soybase.org/sbt/search/search_results.php?category=FeatureName&search_term='
    araburlBase <- 'http://arabidopsis.org/servlets/TairObject?type=locus&name='
    sorgurlBase <- 'http://phytozome.jgi.doe.gov/pz/portal.html#!gene?search=1&detail=1&searchText=transcriptid:'
    legumeInfo_urlBase <- 'https://legumeinfo.org/gene_links/'
    
    annotYvalReverse <- 0.02
    annotYvalForward <- annotYvalReverse + 0.04
    if(values$organism2 == "Corn"){
      annotTable <- adply(thisAnnot[thisAnnot$transcript_strand==1,],1,function(x) {data.frame(x=c(x$transcript_start,x$transcript_end,x$transcript_end),y=c(annotYvalForward,annotYvalForward,NA),url=paste0(urlBase,x$transcript_id),
        name=sprintf("<table cellpadding='4' style='line-height:1.5'><tr><th>%1$s</th></tr><tr><td align='left'>Location: %2$s-%3$s<br>Chromosome: %4$s, Strand: %5$s<br>%6$s</td></tr></table>",
          x$translation_id,
          prettyNum(x$transcript_start, big.mark = ","),
          prettyNum(x$transcript_end, big.mark = ","),
          x$chromosome,
          x$transcript_strand,
          x$V2
        ),
        marker=c(NA,"Arrow",NA),
        stringsAsFactors=FALSE)})
      
      annotTableReverse <- adply(thisAnnot[thisAnnot$transcript_strand==-1,],1,function(x) {data.frame(x=c(x$transcript_start,x$transcript_end,x$transcript_end),y=c(annotYvalReverse,annotYvalReverse,NA),url=paste0(urlBase,x$transcript_id),
        name=sprintf("<table cellpadding='4' style='line-height:1.5'><tr><th>%1$s</th></tr><tr><td align='left'>Location: %2$s-%3$s<br>Chromosome: %4$s, Strand: %5$s<br>%6$s</td></tr></table>",
          x$translation_id,
          prettyNum(x$transcript_start, big.mark = ","),
          prettyNum(x$transcript_end, big.mark = ","),
          x$chromosome,
          x$transcript_strand,
          x$V2
        ),
        marker=c("Arrow",NA,NA),
        stringsAsFactors=FALSE)})
      
    }else if(values$organism2 == "Soybean"){#strand is '+' or '-'
      annotTable <- adply(thisAnnot[thisAnnot$strand=="+",],1,function(x) {data.frame(x=c(x$transcript_start,x$transcript_end,x$transcript_end),y=c(annotYvalForward,annotYvalForward,NA),url=paste0(soyurlBase,x$transcript_id),
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
        stringsAsFactors=FALSE)})
      
      annotTableReverse <- adply(thisAnnot[thisAnnot$strand=="-",],1,function(x) {data.frame(x=c(x$transcript_start,x$transcript_end,x$transcript_end),y=c(annotYvalReverse,annotYvalReverse,NA),url=paste0(soyurlBase,x$transcript_id),
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
        stringsAsFactors=FALSE)})
    }else if(values$organism2 %in% c("Arabidopsis", "Arabidopsis thaliana")){#strand is '+' or '-'
      annotTable <- adply(thisAnnot[thisAnnot$strand=="+",],1,function(x) {data.frame(x=c(x$transcript_start,x$transcript_end,x$transcript_end),y=c(annotYvalForward,annotYvalForward,NA),url=paste0(araburlBase,x$Locus),
        name=sprintf("<table cellpadding='4' style='line-height:1.5'><tr><th>%1$s</th></tr><tr><td align='left'>Location: %2$s-%3$s<br>Chromosome: %4$s, Strand: %5$s<br>Short Desc.: %6$s</td></tr></table>",
          x$name,
          prettyNum(x$transcript_start, big.mark = ","),
          prettyNum(x$transcript_end, big.mark = ","),
          x$chromosome,                                                                           
          x$strand,
          x$short_description
          # x$Curator_summary
        ),
        stringsAsFactors=FALSE)})   
      
      annotTableReverse <- adply(thisAnnot[thisAnnot$strand=="-",],1,function(x) {data.frame(x=c(x$transcript_start,x$transcript_end,x$transcript_end),y=c(annotYvalReverse,annotYvalReverse,NA),url=paste0(araburlBase,x$Locus),
        name=sprintf("<table cellpadding='4' style='line-height:1.5'><tr><th>%1$s</th></tr><tr><td align='left'>Location: %2$s-%3$s<br>Chromosome: %4$s, Strand: %5$s<br>Short Desc.: %6$s</td></tr></table>",
          x$name,
          prettyNum(x$transcript_start, big.mark = ","),
          prettyNum(x$transcript_end, big.mark = ","),
          x$chromosome,                                                                           
          x$strand,
          x$short_description
          # x$Curator_summary
        ),
        stringsAsFactors=FALSE)})
    } else if (values$organism2 == "Medicago truncatula") { # strand is '+' or '-'
      annotTable <- adply(thisAnnot[thisAnnot$strand=="+",],1,function(x) {data.frame(x=c(x$transcript_start,x$transcript_end,x$transcript_end),y=c(annotYvalForward,annotYvalForward,NA),url=paste0(legumeInfo_urlBase, x$name, "/json"),
        name=sprintf("<table cellpadding='4' style='line-height:1.5'><tr><th>%1$s</th></tr><tr><td align='left'>Location: %2$s-%3$s<br>Chromosome: %4$s, Strand: %5$s<br>Desc: %6$s</td></tr></table>",
          x$name,
          prettyNum(x$transcript_start, big.mark = ","),
          prettyNum(x$transcript_end, big.mark = ","),
          x$chromosome,
          x$strand,
          x$description
        ),
        stringsAsFactors=FALSE)})
      
      annotTableReverse <- adply(thisAnnot[thisAnnot$strand=="-",],1,function(x) {data.frame(x=c(x$transcript_start,x$transcript_end,x$transcript_end),y=c(annotYvalReverse,annotYvalReverse,NA),url=paste0(legumeInfo_urlBase, x$name, "/json"),
        name=sprintf("<table cellpadding='4' style='line-height:1.5'><tr><th>%1$s</th></tr><tr><td align='left'>Location: %2$s-%3$s<br>Chromosome: %4$s, Strand: %5$s<br>Desc: %6$s</td></tr></table>",
          x$name,
          prettyNum(x$transcript_start, big.mark = ","),
          prettyNum(x$transcript_end, big.mark = ","),
          x$chromosome,
          x$strand,
          x$description
        ),
        stringsAsFactors=FALSE)})
    }else{#} if(input$organism == "Sorghum"){#strand is '+' or '-'
      annotTable <- adply(thisAnnot[thisAnnot$strand=="+",],1,function(x) {data.frame(x=c(x$transcript_start,x$transcript_end,x$transcript_end),y=c(annotYvalForward,annotYvalForward,NA),url=paste0(sorgurlBase,x$ID),
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
        stringsAsFactors=FALSE)})   
      
      annotTableReverse <- adply(thisAnnot[thisAnnot$strand=="-",],1,function(x) {data.frame(x=c(x$transcript_start,x$transcript_end,x$transcript_end),y=c(annotYvalReverse,annotYvalReverse,NA),url=paste0(sorgurlBase,x$ID),
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
        stringsAsFactors=FALSE)})
    }
    
    annotTable <- annotTable[,c("x","y","name","url")]
    
    if(nrow(annotTableReverse)==0){
      annotTableReverse <- data.frame(x=character(0),y=character(0),name=character(0),url=character(0),stringsAsFactors = FALSE)
    }
    annotTableReverse <- annotTableReverse[,c("x","y","name","url")]
    
    annotArray <- toJSONArray2(annotTable, json = F, names = T)
    annotArrayReverse <- toJSONArray2(annotTableReverse, json = F, names = T)
    
    b <- rCharts::Highcharts$new()
    b$LIB$url <- 'highcharts/'
    b$chart(zoomType="x",alignTicks=FALSE,events=list(click = "#!function(event) {this.tooltip.hide();}!#"))
    b$xAxis(title = list(text = "Base Pairs"),startOnTick=FALSE,min=winLow,max=winHigh,endOnTick=FALSE)      
    if(input$axisLimBool2 == TRUE){
      b$yAxis(title=list(text=input$yAxisColumn2),min=input$axisMin2,max=input$axisMax2,startOnTick=FALSE)
      #create a hidden axis to put the gene track on, all the options are setting to hide everything from the axis 
      b$yAxis(labels=list(enabled=FALSE),title=list(text=NULL),min=0,max=1,lineWidth=0,gridLineWidth=0,minorGridLineWidth=0,lineColor="transparent",minorTickLength=0,tickLength=0,startOnTick=FALSE,opposite=TRUE,replace=FALSE)
    }else{      
      b$yAxis(title=list(text=input$yAxisColumn2),startOnTick=FALSE) 
      #create a hidden axis to put the gene track on, all the options are setting to hide everything from the axis
      b$yAxis(labels=list(enabled=FALSE),title=list(text=NULL),min=0,max=1,lineWidth=0,gridLineWidth=0,minorGridLineWidth=0,lineColor="transparent",minorTickLength=0,tickLength=0,startOnTick=FALSE,opposite=TRUE,replace=FALSE)
    }
    
    if(input$supportInterval2==TRUE){
      if(input$SIaxisLimBool2 == TRUE){
        b$yAxis(title=list(text=input$SIyAxisColumn2),min=input$SIaxisMin2,max=input$SIaxisMax2,gridLineWidth=0,minorGridLineWidth=0,startOnTick=FALSE,opposite=TRUE,replace=FALSE)
      }else{
        b$yAxis(title=list(text=input$SIyAxisColumn2),gridLineWidth=0,minorGridLineWidth=0,startOnTick=FALSE,opposite=TRUE,replace=FALSE)   
      }
      
      if(SIchart[1,input$SIyAxisColumn2] != -1){
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
    
    if(zoomChart[1,input$yAxisColumn2] != -1){
      invisible(sapply(zoomSeries, function(x) {if(length(x)==0){return()};b$series(data = x, type = "scatter", color = colorTable$color[colorTable$trait == as.character(x[[1]]$trait)], name = paste0(x[[1]]$trait))}))
    }
    
    b$series(
      data = annotArray,
      type = "line",
      name = "Forward Genes",
      id = "forward-genes",
      zIndex = 1,
      color = "#53377A",
      yAxis = 1
    )
    
    b$series(
      data = annotArrayReverse,
      type = "line",
      name = "Reverse Genes",
      id = "reverse-genes",
      zIndex = 1,
      color = "#53377A",
      yAxis = 1
    )
    
    if (!is.null(values$glGenes2)) {
      apply(values$glGenes2, 1, FUN = function(g) {
        g.strand <- as.integer(g$strand)
        yh <- -1
        if (g.strand == 1) {
          yh <- annotYvalForward
          sid <- "forward-genes"
        } else if (g.strand == -1) {
          yh <- annotYvalReverse
          sid <- "reverse-genes"
        }
        if (yh > 0 && g$chr == input$chr2) {
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

    b$chart(zoomType="x",alignTicks=FALSE,events=list(click = "#!function(event) {this.tooltip.hide();}!#"))
    b$plotOptions(
      scatter = list(
        cursor = "pointer",
        point = list(
          events = list(
            click = "#! function(event) {alert(this.trait);} !#")), #display popup
        marker = list(
          symbol = "circle",
          radius = 5
        )
      ),
      line = list(
        lineWidth = 6,
        cursor = "pointer",
        point = list(
          events = list(
            click = paste0(
              "#! function() {",
                "if (this.url.includes('legumeinfo.org')) {",
                  # From the JSON at this.url, extract the URLs related to this gene.
                  # Note that this.url = legumeInfo_urlBase + geneString + '/json'
                  #  legumeInfo_urlBase currently has 34 characters (see above)
                  #  and geneString = <5-character species abbreviation>.geneName
                  # And for now, add the gene family phylogram URL by hand.
                  "$.getJSON(this.url, function(data) {
                    var geneString = this.url.substring(34, this.url.indexOf('/json'));
                    var geneName = geneString.substring(6);
                    var content = '';
                    if (data.length == 0) {
                      content = '<p>No ' + geneName + ' links found.</p>';
                    } else {
                      $.each(data, function(i, obj) {
                        content = content + '<p><a href=' + obj.href + ' target=_blank>' + obj.text + '</a></p>';
                        if (i == 0) {
                          var urlPhylogram = 'http://legumeinfo.org/chado_gene_phylotree_v2?gene_name=' + geneString;
                          var textPhylogram = 'View LIS gene family phylogram page for : ' + geneName;
                          content = content + '<p><a href=' + urlPhylogram + ' target=_blank>' + textPhylogram + '</a></p>';
                        }
                      });
                    }

                    var $div = $('<div></div>');
                    $div.html(content);
                    $div.dialog({
                      title: geneName + ' Links',
                      width: 512,
                      height: 'auto',
                      modal: true
                    });
                  });",

                "} else {",
                  # for all other species
                  "window.open(this.url);", #open webpage
                "}",
              "} !#"
          ))),
        marker = list(
          enabled = FALSE,
          radius = 1,
          states = list(hover = list(enabled=FALSE))
        )
      )            
    )
    #it seems almost impossible to get the tooltip to hover along the chart with this version of highcharts (4.0.1), perhaps a question to stackoverflow could solve it.
    #see an example of the problem here: http://jsfiddle.net/N5ymb/
    #one hack/fix would be to add dummy points to the middle of the line that show up when moused over
    b$tooltip(snap=5, useHTML = T, formatter = "#! function() { return this.point.name; } !#")
    b$exporting(enabled=TRUE,filename='zoomChart',sourceWidth=2000)
    b$credits(enabled=TRUE)
    b$set(dom = 'zChart2')
    return(b)
  })#end zchart2
  
  #highcharts test chart
  output$testChart <- renderChart({
    h1 <- hPlot(x = "Wr.Hnd", y = "NW.Hnd", data = MASS::survey, type = c("line", 
                                                                          "bubble", "scatter"), group = "Clap", size = "Age")
    h1$set(dom = 'testChart')
    return(h1)     
  })
  
  #return color table (which sets colors for series in charts) up to date with all combinations of traits
  #this uses allColors set up in global, colors will repeat after 30series
   colorTable <- reactive({     
     traitVals <- list()
     if(input$plotAll == FALSE){
       for(i in input$traitColumns){
         traitVals[[i]] <- input[[i]]
       }     
       
       traits <- do.call(paste,c(expand.grid(traitVals),sep="_"))     
       if(length(traits)==0){return(NULL)}
       
       colorTable <- data.frame(trait=traits,color=rep(allColors,ceiling(length(traits)/30))[1:length(traits)])
     }else{
       colorTable <- data.frame(trait=input$datasets,color=allColors[1])
     }
     colorTable
   })
   colorTable2 <- reactive({     
     traitVals <- list()
     if(input$plotAll2 == FALSE){
       for(i in input$traitColumns2){
         i2 <- paste0(i, "2")
         traitVals[[i]] <- input[[i2]]
       }
       
       traits <- do.call(paste,c(expand.grid(traitVals),sep="_"))     
       if(length(traits)==0){return(NULL)}
       
       colorTable <- data.frame(trait=traits,color=rep(allColors,ceiling(length(traits)/30))[1:length(traits)])
     }else{
       colorTable <- data.frame(trait=input$datasets2,color=allColors[1])
     }
     colorTable
   })
   
  findGWASOverlaps <- function(genomeChart){  
    if(is.null(input$overlapSize)){return(genomeChart[1,])}
    tableIn <- genomeChart
    tableIn$winStart <- tableIn[,input$bpColumn]-input$overlapSize
    tableIn$winStop <- tableIn[,input$bpColumn]+input$overlapSize
    
    allGr <- GRanges(tableIn[,input$chrColumn], IRanges(start=tableIn$winStart,end=tableIn$winStop))

    tableIn$group <- subjectHits(findOverlaps(allGr, reduce(allGr)))
    
    #just groups that have more than one unique SNP
    gwasDataOverlap <- tableIn[tableIn$group %in% as.data.frame(table(tableIn$group))[as.data.frame(table(tableIn$group))$Freq>1,"Var1"],]
    
    #just groups that have more than one unique phenotype
    gwasDataOverlapDiffPheno <- ddply(gwasDataOverlap,.(group),function(x){if(nrow(unique(as.data.frame(x[,"trait"])))>=input$numOverlaps){x}else{x[0,]}})

    return(gwasDataOverlapDiffPheno)    
  }
  findGWASOverlaps2 <- function(genomeChart){  
    if(is.null(input$overlapSize2)){return(genomeChart[1,])}
    tableIn <- genomeChart
    tableIn$winStart <- tableIn[,input$bpColumn2]-input$overlapSize2
    tableIn$winStop <- tableIn[,input$bpColumn2]+input$overlapSize2
    
    allGr <- GRanges(tableIn[,input$chrColumn2], IRanges(start=tableIn$winStart,end=tableIn$winStop))
    
    tableIn$group <- subjectHits(findOverlaps(allGr, reduce(allGr)))
    
    #just groups that have more than one unique SNP
    gwasDataOverlap <- tableIn[tableIn$group %in% as.data.frame(table(tableIn$group))[as.data.frame(table(tableIn$group))$Freq>1,"Var1"],]
    
    #just groups that have more than one unique phenotype
    gwasDataOverlapDiffPheno <- ddply(gwasDataOverlap,.(group),function(x){if(nrow(unique(as.data.frame(x[,"trait"])))>=input$numOverlaps2){x}else{x[0,]}})
    
    return(gwasDataOverlapDiffPheno)    
  }
  
   observe({     
    if(is.null(input$SubmitColsButton) || input$SubmitColsButton == 0){return()}
    isolate({
      currDatasetProp <- datasetProp()
      #print("before")
      #print(currDatasetProp)
      if(as.character(input$datasets) %in% currDatasetProp$dataset){
        currDatasetProp <- currDatasetProp[currDatasetProp$dataset != as.character(input$datasets),]
      }
      #print("after")
      #print(currDatasetProp)
      cols <- varnames()
      #print("data.frame")
      #print(data.frame(dataset=input$datasets,chrColumn=names(cols[cols==input$chrColumn]),bpColumn=names(cols[cols==input$bpColumn]),
      #                 traitCol=paste(names(cols[cols %in% input$traitColumns]),collapse=";"),yAxisColumn=names(cols[cols==input$yAxisColumn]),axisLim=input$axisLimBool,axisMin=input$axisMin,axisMax=input$axisMax,stringsAsFactors=FALSE))
      currDatasetProp <-  rbind(currDatasetProp,data.frame(dataset=input$datasets,chrColumn=names(cols[cols==input$chrColumn]),bpColumn=names(cols[cols==input$bpColumn]),
                                                  traitCol=paste(names(cols[cols %in% input$traitColumns]),collapse=";"),yAxisColumn=names(cols[cols==input$yAxisColumn]),
                                                  logP=input$logP,axisLim=input$axisLimBool,axisMin=input$axisMin,axisMax=input$axisMax,organism=values$organism,plotAll=input$plotAll,
                                                  supportInterval=input$supportInterval,SIyAxisColumn=input$SIyAxisColumn,SIbpStart=input$SIbpStart,SIbpEnd=input$SIbpEnd,
                                                  SIaxisLimBool=input$SIaxisLimBool,SIaxisMin=input$SIaxisMin,SIaxisMax=input$SIaxisMax,stringsAsFactors=FALSE))      
      #print("rbind")
      #print(currDatasetProp)
      write.table(file="./www/config/datasetProperties.csv",x=currDatasetProp,col.names=TRUE,row.names=FALSE,sep=",")
      updateTabsetPanel(session, "datatabs", selected = "panel1")
     })
#    if(input$selected != 1e5){
#      updateTabsetPanel(session, "datatabs", selected = "panel1")  
#    }
   })
   observe({     
     if(is.null(input$SubmitColsButton2) || input$SubmitColsButton2 == 0){return()}
     isolate({
       currDatasetProp <- datasetProp()
       if(as.character(input$datasets2) %in% currDatasetProp$dataset){
         currDatasetProp <- currDatasetProp[currDatasetProp$dataset != as.character(input$datasets2),]
       }
       cols <- varnames2()
       currDatasetProp <-  rbind(currDatasetProp,data.frame(dataset=input$datasets2,chrColumn=names(cols[cols==input$chrColumn2]),bpColumn=names(cols[cols==input$bpColumn2]),
         traitCol=paste(names(cols[cols %in% input$traitColumns2]),collapse=";"),yAxisColumn=names(cols[cols==input$yAxisColumn2]),
         logP=input$logP2,axisLim=input$axisLimBool2,axisMin=input$axisMin2,axisMax=input$axisMax2,organism=values$organism2,plotAll=input$plotAll2,
         supportInterval=input$supportInterval2,SIyAxisColumn=input$SIyAxisColumn2,SIbpStart=input$SIbpStart2,SIbpEnd=input$SIbpEnd2,
         SIaxisLimBool=input$SIaxisLimBool2,SIaxisMin=input$SIaxisMin2,SIaxisMax=input$SIaxisMax2,stringsAsFactors=FALSE))      
       write.table(file="./www/config/datasetProperties.csv",x=currDatasetProp,col.names=TRUE,row.names=FALSE,sep=",")
       updateTabsetPanel(session, "datatabs", selected = "panel1")
     })
   })
   
   observe({
     if(is.null(input$saveDatasetButton) || input$saveDatasetButton == 0){return()}
       isolate({
          if(!file.exists(paste0("./www/config/data/",input$datasets))){
            write.table(getdata(),paste0("./www/config/data/",input$datasets),sep=",",col.names=TRUE,row.names=FALSE) 
          }
       })
   })
   observe({
     if(is.null(input$saveDatasetButton2) || input$saveDatasetButton2 == 0){return()}
     isolate({
       if(!file.exists(paste0("./www/config/data/",input$datasets2))){
         write.table(getdata2(),paste0("./www/config/data/",input$datasets2),sep=",",col.names=TRUE,row.names=FALSE) 
       }
     })
   })
   
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
   observe({
    center <- as.numeric(input$selected[[1]])
    winHigh <- center + input$window[1]
    winLow <- center - input$window[1]
    #eventually I would use winLow/winHigh to change the plotband range
    band = list(from = winLow, to = winHigh, color = "rgba(68, 170, 213, 0.4)")
    #print(band)
    session$sendCustomMessage(type = "customMsg", band)
   })  
   observe({
     center <- as.numeric(input$selected2[[1]])
     winHigh <- center + input$window2[1]
     winLow <- center - input$window2[1]
     #eventually I would use winLow/winHigh to change the plotband range
     band = list(from = winLow, to = winHigh, color = "rgba(68, 170, 213, 0.4)")
     session$sendCustomMessage(type = "customMsg2", band)
   })
   
  output$selectedGene <- renderUI({
    h5(values$glSelectedGene)
  })
  output$relatedRegions <- renderUI({
    selectInput("relatedRegions", "Related Regions:", choices = NULL, selectize = FALSE)
  })
  observe({
    if (is.null(input$relatedRegions) || length(input$relatedRegions) == 0) return()

    ss <- strsplit(input$relatedRegions, split = " ")[[1]]
    chr <- as.integer(stri_sub(ss[1], 4))
    ss2 <- strsplit(ss[2], split = "-")[[1]]
    centerBP <- as.integer(1.0e6*mean(as.numeric(ss2)))
    updateSelectInput(session, "chr2", selected = chr) #input$relatedRegions$chr)
    updateNumericInput(session, "selected2", value = centerBP) #(input$relatedRegions$minBP + input$relatedRegions$maxBP) %/% 2)
  })

  observe({
    # Handle and display genomic linkage query results
    if (is.null(input$genomicLinkages)) return()

    # Parse neighboring genes from species 1
    results1 <- input$genomicLinkages$results1
    values$glSelectedGene <- results1$genes[[(length(results1$genes) + 1) %/% 2]]$name
    glChr1 <- as.integer(stri_match(results1$chromosome_name, regex = "(?i)(?<=\\.chr)\\d+$")[, 1])
    glGenes1 <- data.frame(matrix(unlist(results1$genes), nrow = length(results1$genes), byrow = TRUE),
      stringsAsFactors = FALSE)[, 3:6]
    glGenes1$chr <- glChr1
    names(glGenes1) <- c("family", "fmin", "fmax", "strand", "chr")
    glGenes1 <- glGenes1[nchar(glGenes1$family) > 0, ]

    # Convert (for example) "Medicago truncatula" to "M.truncatula"
    ss.org2 <- strsplit(values$organism2, split = " ")[[1]]
    abbrSpeciesName2 <- paste(stri_sub(ss.org2[1], 1, 1), ss.org2[2], sep = ".")
    # Parse related genes from species 2
    results2 <- input$genomicLinkages$results2
    if (length(results2$groups) == 0) {
      values$glGenes1 <- values$glGenes2 <- NULL
      updateSelectInput(session, "relatedRegions", choices = NULL)

    } else {
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
      families <- intersect(glGenes1$family, glGenes2$family)
      nf <- length(families)
      if (nf == 0) return()
      # Create nf colors
      fc <- rainbow(nf, end = 5/6) # TODO: a more clearly distinguishable set of colors
      familyColors <- vector("list", nf)
      for (i in 1:nf) familyColors[[families[i]]] <- stri_sub(fc[i], 1, 7)

      # Construct the chart data
      glGenes1 <- glGenes1[glGenes1$family %in% families, ]
      glGenes1$color <- familyColors[glGenes1$family]
      glGenes2 <- glGenes2[glGenes2$family %in% families, ]
      glGenes2$color <- familyColors[glGenes2$family]
      values$glGenes1 <- glGenes1
      values$glGenes2 <- glGenes2

      # Construct the related regions (each corresponds to a group from results2$groups)
      glGroupIds <- unique(glGenes2$id)
      if (length(glGroupIds) == 0) {
        glRelatedRegions <- NULL
      } else {
        glRelatedRegions <- unlist(compact(lapply(results2$groups, FUN = function(gr) {
          if (gr$id %in% glGroupIds) {
            gr.chr <- as.integer(stri_match(gr$chromosome_name, regex = "(?i)(?<=\\.chr)\\d+$")[, 1])
            gr.minBP <- gr$genes[[1]]$fmin
            gr.maxBP <- gr$genes[[length(gr$genes)]]$fmax
            sprintf("chr%d %3.2f-%3.2f Mbp", gr.chr, gr.minBP*1.0e-6, gr.maxBP*1.0e-6)
          }
        })))
      }
      updateSelectInput(session, "relatedRegions", choices = glRelatedRegions)
    }
  })

#  observe({
#     print(input$datatabs)     
#  })
})#end server
