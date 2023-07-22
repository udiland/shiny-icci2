library(ggplot2)
library(shiny)
library(shinyFiles)
library(htmlwidgets)
library(shinycssloaders)
library(DT)
library(data.table)

ui <- fluidPage(

  tags$head(
    tags$link(
      rel="stylesheet",
      type="text/css",
      href="stile.css"),

    tags$div(class="tau_logo",
             tags$img(src="logo.png")
    ),

    tags$div(tags$span(id="faculty_name", class="another_gray","The George S. Wise"),

             tags$span(class="black", "Faculty of Life Sciences")),

    tags$span("Intitute for Cereal Crops Research"),

    tags$p(class="universisy-name", "Tel Aviv University")


  ),

  tags$hr(),
  
  tags$body(
    
    h1("Plot common snps"),
    
    fluidRow(
      column(3,
             # create a selection pan for input files
             shinyDirButton("folder", label = "Choose folder", "Choose folder"),
             
             
             # get a title for the plot
             textInput("title","Enter plot title", ""),
             
             # sellect number of runs in which snp appear
             numericInput("n_runs", "Show snps common in more than:", value = 1),
            
              # create a run button fo the manhattan plot
             actionButton("run", "Run"),
             
             # create a button to download the plot
             downloadButton("downloadPlot", "Download plot")
             
      ),
      
      # create a display for the graphs output with spinner
      column(9,
             
             plotOutput('plot'),
             # for debugging
             verbatimTextOutput('x4')
             
             # conditionalPanel(
             #   condition = "input.run > 0",
             #   withSpinner(plotOutput('plot'))
             # )
      )
  
      ), # close row
      
    # show table
    DTOutput("snp_table"),
    
    # create a button to download the table
    downloadButton("downloadTable", "Download table")
    
    ))

    
server <- function(input, output) {
  
  #set a container for reactive values (to create a reactive value use: values$name, than it will be available through all applications)
  values <- reactiveValues()
  
  #values$files <- list.files("/Users/udila/Documents/ICCI/shiny/common_snps/tables")

  #values$InputDir <- "/Users/udila/Documents/ICCI/shiny/common_snps/tables"

  volumes = getVolumes()

   observe({
     #shinyDirChoose(input, 'folder', roots=c(home = '~'), filetypes=c("","csv"))
     shinyDirChoose(input, 'folder', roots=volumes(), filetypes=c("","csv"))
     if(!is.null(input$folder)){
   
       values$InputDir <- parseDirPath(volumes, input$folder)
       #output$x4 <- renderText(InputDir)
       values$files <- list.files(values$InputDir, pattern="*.csv")
       }
     })
  

  #
  # # # create pie plot
  analysis <- observeEvent(input$run, {


    # create an empty data frame with column names
    tbl <- data.frame(matrix(ncol = 4, nrow = 0))
    colnames(tbl) <- c("SNP", "model", "N_models", "Run")

  # go over each file, read it and output one table
  for (f in values$files){
    #res <- read.csv(paste0("/Users/udila/Documents/ICCI/shiny/common_snps/tables/", f))
    
    res <- read.csv(paste0(values$InputDir, "/", f))
    res <- res[c("SNP", "model", "N_models")]

    # split the file name and take all but the extension
    f_split <- strsplit(f, "[.]")[[1]]
    res$Run <- paste(f_split[1:length(f_split) -1], collapse = ".")

    # use setDT (conversion to data.frame) to coerce same snps to one column and create a column with models name
    res <- setDT(res)[,.(models = paste(model, collapse = ",")), by = c("SNP", "N_models", "Run")]
    tbl <- rbind(tbl,res)
    }

  # summarise by Runs
  all_snps <- setDT(tbl)[,.(Runs = paste(Run, collapse = ",")), by = c("SNP")]

  # add a column that count the number of runs a snps is in
  all_snps$n_runs <- sapply(all_snps$Runs, function(x) length(strsplit(x, ",")[[1]]))
  
  values$tbl <- all_snps

  # if the number of runs selected is bigger than the maximum give a message, else plot the data
  if (input$n_runs > max(all_snps$n_runs) - 1){
    
    plt <- NULL
    
    output$plot <- renderPlot(plt)
    
    output$x4 <- renderPrint("No snps to display")
    
  } else if (input$n_runs <= max(all_snps$n_runs)) {
    
    output$x4 <- NULL
    
  plt <- ggplot(all_snps[all_snps$n_runs > input$n_runs], aes(x="", y=Runs, fill=SNP)) +
    geom_bar(stat = "identity") +
    coord_polar("y", start = 0) +
    theme_void() +
    ggtitle(input$title) +
    theme(plot.title = element_text(hjust = 0.5, size = 22))

  output$plot <- renderPlot(plt)
  
  # save the plot
  output$downloadPlot <- downloadHandler(
    filename = paste(input$title,".pdf", sep = ""),
    content = function(file) {
      ggsave(file, plot=plt)
    })
  } # close else if
  
  # display the table
  output$snp_table <- renderDT(all_snps)
  
  # download the table
  output$downloadTable <- downloadHandler(
    filename = paste(input$title, ".csv", sep = ""),
    content = function(file) {
      write.csv(all_snps, file, row.names = F)
    })

  }) # close reactive
  
}

shinyApp(ui = ui, server = server)

