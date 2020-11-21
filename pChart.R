source("common.R")

create_pChart <- function(j, input, values) {
  nid <- jth_ref("notify.create_pChart", j)
  showNotification(paste0("Creating Chromosome chart for ", values[[jth_ref("organism", j)]], ". Please wait."),
    duration = NULL, id = nid, type = "message")

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
    chromChart <- findGWASOverlaps(chromChart, j, input)
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
  colorTable <- getColorTable(j, input)
  
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
    if (dynamic.interval.height) {
      SIchart$h <- rank(SIchart$trait, ties = "first")
    } else {
      SIchart$h <- SIchart[[input[[jth_ref("SIyAxisColumn", j)]]]]
    }
    SIchart <- SIchart[order(SIchart[[input[[jth_ref("SIbpStart", j)]]]]),]
    jlTable <- adply(SIchart,1,function(x) {data.frame(x=c(x[[input[[jth_ref("SIbpStart", j)]]]],x[[input[[jth_ref("SIbpEnd", j)]]]],x[[input[[jth_ref("SIbpEnd", j)]]]]),y=c(x$h,x$h,NA),trait=x$trait,
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
  chrNumber <- trailingInteger(input[[jth_ref("chr", j)]])
  a$xAxis(title = list(text = "Base Pairs"),startOnTick=TRUE,min=1,max=chrSize[values[[jth_ref("organism", j)]]][[1]][chrNumber],endOnTick=FALSE,
          plotBands = list(list(from=pbWin$winLow,to=pbWin$winHigh,color='rgba(68, 170, 213, 0.4)')))
  
  if(input[[jth_ref("axisLimBool", j)]] == TRUE){
    a$yAxis(title=list(text=input[[jth_ref("yAxisColumn", j)]]),min=input[[jth_ref("axisMin", j)]],max=input[[jth_ref("axisMax", j)]],startOnTick=FALSE)
  }else{
    a$yAxis(title=list(text=input[[jth_ref("yAxisColumn", j)]]),min=0,startOnTick=FALSE)
  }    
  
  if(input[[jth_ref("supportInterval", j)]]==TRUE){
    if(input[[jth_ref("SIaxisLimBool", j)]] == TRUE){
      a$yAxis(visible=FALSE,title=list(text=input[[jth_ref("SIyAxisColumn", j)]]),min=input[[jth_ref("SIaxisMin", j)]],max=input[[jth_ref("SIaxisMax", j)]],gridLineWidth=0,minorGridLineWidth=0,startOnTick=FALSE,opposite=TRUE,replace=FALSE)
    }else{
      a$yAxis(visible=FALSE,title=list(text=input[[jth_ref("SIyAxisColumn", j)]]),min=0,gridLineWidth=0,minorGridLineWidth=0,startOnTick=FALSE,opposite=TRUE,replace=FALSE)
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

  removeNotification(nid)

  a$exporting(enabled=TRUE,filename='chromChart',sourceWidth=2000)
  a$set(dom = jth_ref('pChart', j))
  return(a)
}
