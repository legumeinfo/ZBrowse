source("common.R")

create_pChartMacro <- function(j, input, values) {
  nid <- jth_ref("notify.create_pChartMacro", j)
  showNotification(paste0("Creating Chromosome macro-synteny chart for ", values[[jth_ref("organism", j)]], ". Please wait."),
    duration = NULL, id = nid, type = "message")

  a <- rCharts::Highcharts$new()
  a$LIB$url <- 'highcharts/' #use the local copy of highcharts, not the one installed by rCharts
  chrNumber <- trailingInteger(input[[jth_ref("chr", j)]])
  a$xAxis(title = list(text = "Base Pairs"),startOnTick=TRUE,min=1,max=chrSize[[values[[jth_ref("organism", j)]]]][chrNumber],endOnTick=FALSE)

  # Display macro-synteny blocks
  blocks <- values$pairwiseBlocks[[j]]
  if (!is.null(blocks)) {
    blocks <- blocks[blocks$chromosome == trailingInteger(input[[jth_ref("chr", j)]]), ]
    if (j == 1) {
      a$yAxis(labels=list(enabled=FALSE),title=list(text=NULL),min=0,max=1,lineWidth=0,gridLineWidth=0,minorGridLineWidth=0,lineColor="transparent",minorTickLength=0,tickLength=0,endOnTick=FALSE)
    } else {
      ylab2 <- paste(input$macroDistance, "distance")
      if (input$macroDistance == "Levenshtein") ylab2 <- paste("Normalized", ylab2)
      a$yAxis(title=list(text=ylab2),min=0,max=1,reversed=TRUE,lineWidth=0,gridLineWidth=0,minorGridLineWidth=0,lineColor="transparent",minorTickLength=0,tickLength=0,endOnTick=FALSE)
    }
    apply(blocks, 1, function(r) {
      r <- data.frame(as.list(r), stringsAsFactors = FALSE) # to avoid "$ operator is invalid for atomic vectors" warning
      if (j == 1) {
        yh <- 0.5
      } else {
        yh <- as.numeric(r$distance)
        if (input$macroDistance == "Levenshtein") yh <- yh*2/(as.numeric(r$n1) + as.numeric(r$n2))
      }
      r.data <- vector("list", 2)
      r.data[[1]]$x <- as.numeric(r$fmin)
      r.data[[2]]$x <- as.numeric(r$fmax)
      r.data[[1]]$y <- r.data[[2]]$y <- yh
      r.distance <- sprintf("%4.3f", as.numeric(r$distance))
      if (input$macroDistance == "Levenshtein") r.distance <- sprintf("%d", as.integer(r$distance))
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
              values$organism, as.integer(r$chr1), as.integer(r$n1), as.integer(r$i), as.integer(r$j), input$macroDistance, r.distance)
          ),
          pointFormat = '',
          followPointer = TRUE
        ),
        # put shortest blocks on top (note zIndex < 0 for macrosynteny blocks)
        zIndex = as.integer(r$fmin) - as.integer(r$fmax)
      )
    })
  }

  chartHeight <- ifelse(j == 1, 150, 300)
  a$chart(height=chartHeight,zoomType="x", alignTicks=FALSE,events=list(click = "#!function(event) {this.tooltip.hide();}!#"))
  a$title(text=paste(values[[jth_ref("organism", j)]],"Macro-Synteny for Chromosome",input[[jth_ref("chr", j)]],sep=" "))
  a$subtitle(text="Rollover for more info. Drag chart area to zoom.")
  a$plotOptions(
    line = list(
      dashStyle = 'Solid',
      cursor = "pointer",
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
