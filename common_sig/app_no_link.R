library(shiny)
library(shinyFiles)
library(shinycssloaders)
library(DT)
library(data.table)
library(dplyr)

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
    
    h1("Search for mutual SNPs"),
    h2("selecet a folder with GAPIT output files from different runs"),
    
    fluidRow(
      column(3,
             # create a selection pan for input files
             shinyDirButton("folder", label = "Choose folder", "Choose folder"),
	     #fileInput(inputId = "folder", label = "Choose folder"),             
             
             # select model
             textInput("model","Enter model", "Type: Blink, MLM etc."),
             
             # sellect number of runs in which snp appear
             numericInput("logpv", "Filter SNPs with LOD higher than:", value = 6),
             
             # create a run button fo the manhattan plot
             actionButton("run", "Run")
             
      ),
      
      # create a display for the graphs output with spinner
      column(9,
        # show table
        conditionalPanel(
          condition = "input.run > 0",
          withSpinner(DTOutput("snp_table"))
        )
      ) # close column
      
    ), # close row
    
    # create a button to download the table
    downloadButton("downloadTable", "Download table")
    
  ))


server <- function(input, output) {
  
  #set a container for reactive values (to create a reactive value use: values$name, than it will be available through all applications)
  values <- reactiveValues()
  
  volumes = getVolumes()
  
  observe({
    # choose from home in the server
    #shinyDirChoose(input, 'folder', roots=c(home='~'), filetypes=c("csv"))
    # choose from root
    shinyDirChoose(input, 'folder', roots=volumes(), filetypes=c("","csv"))
    if(!is.null(input$folder)){
      
      values$InputDir <- parseDirPath(volumes, input$folder)
      #output$x4 <- renderText(InputDir)
      values$files <- list.files(values$InputDir, pattern="*.csv", full.names = T)
    }
  })
  
  
  #
  analysis <- observeEvent(input$run, {
    
    # create an empty data frame with column names
    tbl <- data.table(matrix(ncol = 7, nrow = 0))
    colnames(tbl) <- c("SNP", "Chromosome", "Position", "P.value", "maf", "FDR_Adjusted_P.values", "Trait")
    # create list of files with the user selected model
    files = c()
    for (i in seq_along(values$files)){
      if (grepl(input$model, values$files[i])){
        files[i] <- values$files[i]}
    }
    # go over each file, read it and output one table
    for (f in files){

      res <- read.table(f, header = T, sep = ",")
      res <- res[c("SNP", "Chromosome", "Position", "P.value", "maf", "FDR_Adjusted_P.values")]
      # filter for LOD
      res <- res[-log10(res$P.value) >= input$logpv ,]
      
      # add column with trait and combine with the rest of the data (only if there are results in the table)
      if (dim(res)[1] != 0){
      trait <- strsplit(f, "[.]")[[1]][3]
      res$Trait <- trait

      tbl <- rbind(tbl,res)}
    }
    
    # summarise by snp, chrom, position,maf
    snps <- tbl %>% group_by(SNP, Chromosome, Position, maf) %>% summarize(Trait = paste0(Trait, collapse = ','),
                                                n_trait = length(strsplit(Trait, ",")[[1]]))

    # if the number of runs selected is bigger than the maximum give a message, else plot the data
    if (dim(snps)[1] == 0){
      
      output$x4 <- renderPrint("No snps to display")
      
    } else (output$x4 <- NULL)
    
    # display the table
    output$snp_table <- renderDT(snps)
    
    # download the table
    output$downloadTable <- downloadHandler(
      filename = paste(input$title, ".csv", sep = ""),
      content = function(file) {
        write.csv(snps, file, row.names = F)
      })
    
  }) 
}

shinyApp(ui = ui, server = server)


# paste("http://icci-1.tau.ac.il:8000/cgi-bin/hgTracks?db=hub_112_Longissima&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=",
#       "chr", don$Chromosome, "%3A" ,don$Position, "%2D", don$Position,
#       "&hgsid=329_aFWaA0H7kp7ROLdvGJNjECGBoAin", sep="")

