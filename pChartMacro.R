source("common.R")

create_pChartMacro <- function(j, input, values) {
  nid <- jth_ref("notify.create_pChartMacro", j)
  showNotification(paste0("Creating Chromosome macro-synteny chart for ", values[[jth_ref("organism", j)]], ". Please wait."),
    duration = NULL, id = nid, type = "message")

  #calculate window for plotband
  pbWin <- isolate({
    center <- as.numeric(input[[jth_ref("selected", j)]])
    winHigh <- center + input[[jth_ref("window", j)]]
    winLow <- center - input[[jth_ref("window", j)]]
    list(winLow=winLow,winHigh=winHigh)
  })

  a <- rCharts::Highcharts$new()
  a$LIB$url <- 'highcharts/' #use the local copy of highcharts, not the one installed by rCharts
  chrNumber <- trailingInteger(input[[jth_ref("chr", j)]])
  a$xAxis(title = list(text = "Base Pairs"),startOnTick=TRUE,min=1,max=chrSize[[values[[jth_ref("organism", j)]]]][chrNumber],endOnTick=FALSE,
          plotBands = list(list(from=pbWin$winLow,to=pbWin$winHigh,color='rgba(68, 170, 213, 0.4)')))

  # Display macro-synteny blocks
  blocks <- values$pairwiseBlocks[[j]]
  if (is.null(blocks)) {
    # make an empty data frame in order to display the x axis
    columns <- c("i", "j", "chromosome")
    blocks <- data.frame(matrix(nrow = 0, ncol = length(columns)), stringsAsFactors = FALSE)
    names(blocks) <- columns
  }
  blocks <- blocks[blocks$chromosome == trailingInteger(input[[jth_ref("chr", j)]]), ]
  macroDistanceMetric <- isolate(input$macroDistance)
  if (j == 1) {
    a$yAxis(labels=list(enabled=FALSE),title=list(text=NULL),min=0,max=1,lineWidth=0,gridLineWidth=0,minorGridLineWidth=0,lineColor="transparent",minorTickLength=0,tickLength=0,endOnTick=FALSE)
    blocks2 <- values$pairwiseBlocks[[2]]
    blocks2 <- blocks2[blocks2$chromosome == trailingInteger(input$chr2), ]
  } else {
    ylab2 <- paste(macroDistanceMetric, "distance")
    if (macroDistanceMetric == "Levenshtein") ylab2 <- paste("Normalized", ylab2)
    a$yAxis(title=list(text=ylab2),min=0,max=1,reversed=TRUE,lineWidth=0,gridLineWidth=0,minorGridLineWidth=0,lineColor="transparent",minorTickLength=0,tickLength=0,endOnTick=FALSE)
    blocks1 <- values$pairwiseBlocks[[1]]
    blocks1 <- blocks1[blocks1$chromosome == trailingInteger(input$chr), ]
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
    r.data[[1]]$x <- r.data[[1]]$min <- r.data[[2]]$min <- as.numeric(r$fmin)
    r.data[[2]]$x <- r.data[[1]]$max <- r.data[[2]]$max <- as.numeric(r$fmax)
    r.data[[1]]$y <- r.data[[2]]$y <- yh
    r.distance <- sprintf("%4.3f", as.numeric(r$distance))
    if (macroDistanceMetric == "Levenshtein") r.distance <- sprintf("%d", as.integer(r$distance))
    if (j == 1) {
      # there may be multiple block2s
      block2 <- blocks2[as.integer(blocks2$chr1) == as.integer(r$chromosome) & as.integer(blocks2$i) == as.integer(r$i) & as.integer(blocks2$j) == as.integer(r$j), ]
      if (length(block2$fmin) > 0) {
        r.data[[1]]$minSrc <- r.data[[2]]$minSrc <- as.numeric(block2$fmin)
        r.data[[1]]$maxSrc <- r.data[[2]]$maxSrc <- as.numeric(block2$fmax)
      }
    } else if (j == 2) {
      # there should be only one block1
      block1 <- blocks1[as.integer(blocks1$chromosome) == as.integer(r$chr1) & as.integer(blocks1$i) == as.integer(r$i) & as.integer(blocks1$j) == as.integer(r$j), ]
      if (length(block1$fmin) > 0) {
        r.data[[1]]$minRef <- r.data[[2]]$minRef <- as.numeric(block1$fmin)
        r.data[[1]]$maxRef <- r.data[[2]]$maxRef <- as.numeric(block1$fmax)
      }
    }
    a$series(
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
  a$chart(height=chartHeight,zoomType="x", alignTicks=FALSE,events=list(click = "#!function(event) {this.tooltip.hide();}!#"))
  a$title(text=paste(values[[jth_ref("organism", j)]],"Macro-Synteny for Chromosome",input[[jth_ref("chr", j)]],sep=" "))
  a$subtitle(text="Rollover for more info. Drag chart area to zoom.")
  doClickOnLine <- paste(
    "#! function() {",
      sprintf("this.options['j'] = %d;", j),
      "Shiny.onInputChange('set_pChartMacro', this.options);",
    "} !#"
  )
  a$plotOptions(
    line = list(
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

  removeNotification(nid)

  a$exporting(enabled=TRUE,filename='chromChartMacro',sourceWidth=2000)
  a$credits(enabled=TRUE)
  a$set(dom = jth_ref('pChartMacro', j))
  return(a)
}
