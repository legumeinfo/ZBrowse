shinyUI(pageWithSidebar(  
  
  headerPanel(singleton(tags$head(tags$title("Zbrowse")))),
  sidebarPanel(
    includeCSS('www/style.css'),
#    wellPanel(
#      uiOutput("datasets")
#    ),
    #uiOutput("ui_Manage")
    uiOutput("ui_All"),
    width=2
  ),
  mainPanel(
    tagList( # The four core files: 3 JS files and 1 CSS file --
      useShinyjs(),
#      singleton(tags$head(tags$script(src='js/highcharts.js',type='text/javascript'))),
      singleton(tags$head(tags$script(src='DataTables/js/jquery.dataTables.js',type='text/javascript'))),
      singleton(tags$head(tags$script(src='TableTools/js/TableTools.js',type='text/javascript'))),
      singleton(tags$head(tags$script(src='TableTools/js/ZeroClipboard.js',type='text/javascript'))),
      singleton(tags$head(tags$link(href='TableTools/css/TableTools.css',rel='stylesheet',type='text/css'))),
      #singleton(tags$head(tags$script(src='http://code.highcharts.com/highcharts.js',type='text/javascript')))
      # For jQuery dialogs
      singleton(tags$head(tags$script(src='jquery-ui/jquery-ui.js',type='text/javascript'))),
      singleton(tags$head(tags$link(href='jquery-ui/jquery-ui.css',rel='stylesheet',type='text/css'))),
      # For Broadcast Channel
      singleton(tags$head(tags$script(src='js/BroadcastChannel-GCV.js',type='text/javascript')))
    ),    
    #progressInit(),    
    uiOutput("ui_data_tabs"),
    width=10
    #tableOutput('contents')
  )
))
