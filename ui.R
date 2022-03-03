shinyUI(pageWithSidebar(  
  
  headerPanel(
  introjsUI(), # to enable rintrojs
  singleton(tags$head(tags$title("ZZBrowse")))),
  sidebarPanel(
    includeCSS('www/style.css'),
    h4(HTML("<a href='https://legumeinfo.org' target='_blank'><img src='lis-6044923.png' width='40px' height='40px'></a> ZZBrowse")),
#    wellPanel(
#      uiOutput("datasets")
#    ),
    #uiOutput("ui_Manage")
    uiOutput("ui_All"),
    width=2
  ),
  mainPanel(
    tagList(
      useShinyjs(),
      singleton(tags$head(tags$script(src='DataTables/js/jquery.dataTables.js',type='text/javascript'))),
      singleton(tags$head(tags$script(src='TableTools/js/TableTools.js',type='text/javascript'))),
      singleton(tags$head(tags$script(src='TableTools/js/ZeroClipboard.js',type='text/javascript'))),
      singleton(tags$head(tags$link(href='TableTools/css/TableTools.css',rel='stylesheet',type='text/css'))),
      singleton(tags$head(tags$script(src='highcharts/js/highcharts.js',type='text/javascript'))),
      # singleton(tags$head(tags$script(src='https://code.highcharts.com/7.0.0/highcharts.js',type='text/javascript'))),
      # For jQuery dialogs
      singleton(tags$head(tags$script(src='jquery-ui/jquery-ui.js',type='text/javascript'))),
      singleton(tags$head(tags$link(href='jquery-ui/jquery-ui.css',rel='stylesheet',type='text/css')))
    ),    
    #progressInit(),    
    uiOutput("ui_data_tabs"),
    width=10
    #tableOutput('contents')
  )
))
