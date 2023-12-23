library(shiny)
library(htmlwidgets)
library(shinyFiles)
library(shinycssloaders)
library(ggplot2)
library(dplyr)
library(ggnewscale)
options(shiny.maxRequestSize=1000*1024^2)

# get the function to create the plot
source("plot_snp_kmr_function.R")

ui <- fluidPage(
  # set up the header with TAU logo
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
    
    h1("Plot SNPs and Kmers"),
    
    fluidRow(
      column(3,
             # create a selection pan for input files
             fileInput(inputId = "gwas", label = "Upload GAPIT results", accept = ".csv"),
             # upload kmers data
             fileInput(inputId = "kmers", label = "Upload Kmers data", accept = ".csv"),
             h3("Use the file created in the 'direct' pipline, e.g. 'data_for_plot_pass5.csv'"),
             # sellect p-value cutoff
             selectInput("pv", "P-value cutoff", c("0.05", "0.01", "0.001"), selected = "0.001"),
             # create a run button fo the plot
             actionButton("runPlot", "Run"),
             
             # create a button to download the manhattan plot
             downloadButton("downloadPlot", "Download plot")
             
      ),
      column(9,
             
             conditionalPanel(
               condition = "input.runPlot > 0",
               shinycssloaders::withSpinner(plotOutput('plot'))
             )
      ))))

server <- function(input, output){
  observeEvent(input$runPlot,{
    
    #get the gwas results file:
    file <- input$gwas
    ext <- tools::file_ext(file$datapath)
    gwasResults <- read.csv(file$datapath)
    
    #get the kmers results file:
    file <- input$kmers
    ext <- tools::file_ext(file$datapath)
    kmersResults <- read.csv(file$datapath)
    
    # create the plot
    plt <- plot_snp_kmers(kmersResults, gwasResults, input$pv)
  
  output$plot <- renderPlot(plot_snp_kmers(kmersResults, gwasResults, input$pv))
  })
  
  # download plot
  output$downloadPlot <- downloadHandler(
    filename = "plot.png",
    content = function(file) {
      ggsave(file, plot)}
    )
    
}

shinyApp(ui = ui, server = server)

  