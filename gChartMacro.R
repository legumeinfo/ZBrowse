source("common.R")

#Genome wide chart
create_gChartMacro <- function(j, input, values) {
  nid <- jth_ref("notify.create_gChartMacro", j)
  showNotification(paste0("Creating Whole Genome macro-synteny chart for ", values[[jth_ref("organism", j)]], ". Please wait."),
    duration = NULL, id = nid, type = "message")

  calculateTotalBP(j, input, values)

  #build list for where to put plotbands for this organism
  bigList <- list()
  title.margin <- 15
  cumBP<-c(0,cumsum(as.numeric(chrSize[[values[[jth_ref("organism", j)]]]])))
  for(i in 1:(length(cumBP)-1)){
    chrLabel <- chrName[[values[[jth_ref("organism", j)]]]][i]
    # split character names like "Arahy.01"
    k.dot <- unlist(gregexpr("\\.", chrLabel))
    if (k.dot > 1) {
      chrLabel <- paste(substring(chrLabel, 1, k.dot - 1), substring(chrLabel, k.dot), sep = "<br>")
      title.margin <- 30
    }
    bigList[[i]] <- list(
      from = cumBP[i] + 1, to = cumBP[i + 1],
      label = list(
        text = chrLabel,
        style = list(color = "#6D869F", fontSize = "10px"),
        verticalAlign = "bottom"
      )
    )
    # different color for odd plot bands
    if (i %% 2 != 0) bigList[[i]]$color <- "rgba(68, 170, 213, 0.1)"
  }    
  
  c <- rCharts::Highcharts$new()
  c$LIB$url <- 'highcharts/'
  c$xAxis(title = list(text = "Chromosome", margin = title.margin), startOnTick = TRUE,
    min = 0, max = sum(as.numeric(chrSize[[values[[jth_ref("organism", j)]]]])),
    endOnTick = FALSE, labels = list(enabled = FALSE), tickWidth = 0, plotBands = bigList)

  # Display macro-synteny blocks
  blocks <- values$pairwiseBlocks[[j]]
  if (is.null(blocks)) {
    # make an empty data frame in order to display the x axis
    columns <- c("i", "j", "chromosome")
    blocks <- data.frame(matrix(nrow = 0, ncol = length(columns)), stringsAsFactors = FALSE)
    names(blocks) <- columns
  }
  macroDistanceMetric <- isolate(input$macroDistance)
  chr1 <- trailingInteger(input$macroChromosome) # to filter by macro-synteny species 1 chromosome
  if (j == 1) {
    if (!is.na(chr1)) blocks <- blocks[blocks$chromosome == chr1, ]
    c$yAxis(labels=list(enabled=FALSE),title=list(text=NULL),min=0,max=1,lineWidth=0,gridLineWidth=0,minorGridLineWidth=0,lineColor="transparent",minorTickLength=0,tickLength=0,endOnTick=FALSE)
    blocks2 <- values$pairwiseBlocks[[2]]
  } else {
    if (!is.na(chr1)) blocks <- blocks[blocks$chr1 == chr1, ]
    ylab2 <- paste(macroDistanceMetric, "distance")
    if (macroDistanceMetric == "Levenshtein") ylab2 <- paste("Normalized", ylab2)
    c$yAxis(title=list(text=ylab2),min=0,max=1,reversed=TRUE,lineWidth=0,gridLineWidth=0,minorGridLineWidth=0,lineColor="transparent",minorTickLength=0,tickLength=0,endOnTick=FALSE)
    blocks1 <- values$pairwiseBlocks[[1]]
  }
  apply(blocks, 1, function(r) {
    r <- data.frame(as.list(r), stringsAsFactors = FALSE) # to avoid "$ operator is invalid for atomic vectors" warning
    if (j == 1 && nrow(blocks) > 0) {
      yh <- 0.5
    } else {
      yh <- as.numeric(r$distance)
      if (macroDistanceMetric == "Levenshtein") yh <- yh*2/(as.numeric(r$n1) + as.numeric(r$n2))
    }
    r.data <- vector("list", 2)
    r.data[[1]]$x <- r.data[[1]]$min <- r.data[[2]]$min <- as.numeric(r$cumfmin)
    r.data[[2]]$x <- r.data[[1]]$max <- r.data[[2]]$max <- as.numeric(r$cumfmax)
    r.data[[1]]$y <- r.data[[2]]$y <- yh
    r.distance <- sprintf("%4.3f", as.numeric(r$distance))
    if (macroDistanceMetric == "Levenshtein") r.distance <- sprintf("%d", as.integer(r$distance))
    if (j == 1) {
      # there may be multiple block2s
      block2 <- blocks2[as.integer(blocks2$chr1) == as.integer(r$chromosome) & as.integer(blocks2$i) == as.integer(r$i) & as.integer(blocks2$j) == as.integer(r$j), ]
      r.data[[1]]$minSrc <- r.data[[2]]$minSrc <- as.numeric(block2$cumfmin)
      r.data[[1]]$maxSrc <- r.data[[2]]$maxSrc <- as.numeric(block2$cumfmax)
    } else if (j == 2) {
      # there should be only one block1
      block1 <- blocks1[as.integer(blocks1$chromosome) == as.integer(r$chr1) & as.integer(blocks1$i) == as.integer(r$i) & as.integer(blocks1$j) == as.integer(r$j), ]
      r.data[[1]]$minRef <- r.data[[2]]$minRef <- as.numeric(block1$cumfmin)
      r.data[[1]]$maxRef <- r.data[[2]]$maxRef <- as.numeric(block1$cumfmax)
    }
    c$series(
      type = "line",
      data = r.data,
      color = r$color,
      lineWidth = 6,
      showInLegend = FALSE,
      tooltip = list(
        headerFormat = ifelse(j == 1,
          sprintf("<b>Macro-synteny</b><br>%s chromosome %d<br>Blocks %d-%d<br>Location %s-%s",
            values$organism, as.integer(r$chromosome), as.integer(r$i), as.integer(r$j),
            prettyNum(as.integer(r$fmin), big.mark = ","), prettyNum(as.integer(r$fmax), big.mark = ",")),
          sprintf("<b>Macro-synteny</b><br>%s chromosome %d (%d genes)<br>Location %s-%s Orientation: %s<br>with %s chromosome %d (%d genes)<br>Blocks %d-%d<br>%s distance: %s",
            values$organism2, as.integer(r$chromosome), as.integer(r$n2), prettyNum(as.integer(r$fmin), big.mark = ","), prettyNum(as.integer(r$fmax), big.mark = ","), r$orientation,
            values$organism, as.integer(r$chr1), as.integer(r$n1), as.integer(r$i), as.integer(r$j), macroDistanceMetric, r.distance)
        ),
        pointFormat = '',
        followPointer = TRUE
      ),
      # put shortest blocks on top (note zIndex < 0 for macrosynteny blocks)
      zIndex = as.integer(r$fmin) - as.integer(r$fmax)
    )
  })

  chartHeight <- ifelse(j == 1, 150, 300)
  c$chart(height=chartHeight,zoomType="x",alignTicks=FALSE,events=list(click = "#!function(event) {this.tooltip.hide();}!#"))
  c$title(text = paste(values[[jth_ref("organism", j)]], "Macro-Synteny"))
  c$subtitle(text = "Rollover for more info. Drag chart area to zoom.")
  doClickOnLine <- paste(
    "#! function() {",
      sprintf("this.options['j'] = %d;", j),
      "Shiny.onInputChange('set_gChartMacro', this.options);",
    "} !#"
  )
  c$plotOptions(
    line = list(
      cursor = "pointer",
      point = list(
        events = list(
          click = doClickOnLine
        )
      ),
      marker = list(
        enabled = FALSE,
        states = list(hover = list(enabled=FALSE)),
        symbol = "square"
      )
    )
  )
  c$exporting(enabled=TRUE,filename='genomeChartMacro',sourceWidth=2000)

  removeNotification(nid)

  c$credits(enabled=TRUE)
  c$set(dom = jth_ref('gChartMacro', j))
  return(c)
}
