library(shiny)
library(shinyFiles)
library(dplyr)
library(htmlwidgets)
library(plotly)
library(ggplot2)
library(shinycssloaders)
options(shiny.maxRequestSize=500*1024^2) 



ui <- fluidPage(
  
  
  # create a selection pan for input files
  fileInput(inputId = "glm", label = "Upload GLM results", accept = ".csv"),
  
  fileInput(inputId = "mlm", label = "Upload MLM results", accept = ".csv"),
  
  fileInput(inputId = "mlmm", label = "Upload MLMM results", accept = ".csv"),
  
  fileInput(inputId = "farm", label = "Upload FarmCPU results", accept = ".csv"),
  
  fileInput(inputId = "blink", label = "Upload BLINK results", accept = ".csv"),
  
#  selectInput("n_models", "SNP common to number of models", c("all", "2", "3", "4", "5"), selected = "all"),
  
  # create a button to download
  downloadButton("downloadData", "Download"),
  
  # create a run button
  actionButton("run", "Run"),
  
  
  mainPanel(
    conditionalPanel(
      condition = "input.run > 0",
      withSpinner(plotOutput("plot")))
  )
)


server <- function(input, output) {
  
  gwas_res <- observeEvent(input$run, {
  
  file <- input$glm
  ext <- tools::file_ext(file$datapath)
  glm <- read.csv(file$datapath)
  
  
  glm <- glm[glm['P.value'] < 0.01 ,]
  
  glm$Chromosome[glm$Chromosome == "UN"] <- "Un"
  
  GLM <- glm %>% 
    
    # Compute chromosome size
    group_by(Chromosome) %>% 
    summarise(chr_len=max(Position)) %>% 
    
    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
    select(-chr_len) %>%
    
    # Add this info to the initial dataset
    left_join(glm, ., by=c("Chromosome"="Chromosome")) %>%
    
    # Add a cumulative position of each SNP
    arrange(Chromosome, Position) %>%
    mutate( BPcum=Position+tot)
  
 axisdf = GLM %>% group_by(Chromosome) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
  
  
  #use this calculation of chromosome start position for the other positions
  chr_start_pos <- GLM %>% group_by(Chromosome) %>% summarise(tot = min(tot))
  
  #leave only FDR sig snp
  GLM <- GLM[GLM$FDR_Adjusted_P.values < 0.05 ,]
  
  #add model column
  GLM$model <- "GLM"
  
  GLM$text <- paste("snp: ", GLM$SNP, "\nPosition: ", GLM$Position, 
                    "\nChromosome: ", GLM$Chromosome, "\nFDR:", -log10(GLM$FDR_Adjusted_P.values) %>% round(2),
                    "\nMAF:", GLM$maf %>% round(4), "\nModel:", GLM$model, sep = "")
  
  GLM$link <- paste("http://icci.tau.ac.il:8000/cgi-bin/hgTracks?db=hub_112_Longissima&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=", 
                    "chr", GLM$Chromosome, "%3A" ,GLM$Position, "%2D", GLM$Position,
                    "&hgsid=329_aFWaA0H7kp7ROLdvGJNjECGBoAin", sep="")
  
  
  #####################################################################################################################
  # read data of next model
  
  file <- input$mlm
  ext <- tools::file_ext(file$datapath)
  MLM <- read.csv(file$datapath)

  MLM$Chromosome[MLM$Chromosome == "UN"] <- "Un"
  
  # leave only FDR sig snps
  MLM <- MLM[MLM$FDR_Adjusted_P.values < 0.05 ,]
  
  # make cumulative position the same as the first model
  MLM <- merge(MLM, chr_start_pos, by="Chromosome")
  MLM <- MLM %>% group_by(Chromosome) %>% mutate(BPcum = Position+tot)
  
  
  #add model column
  MLM$model <- "MLM"
  
  MLM$text <- paste("snp: ", MLM$SNP, "\nPosition: ", MLM$Position, 
                    "\nChromosome: ", MLM$Chromosome, "\nFDR:", -log10(MLM$FDR_Adjusted_P.values) %>% round(2),
                    "\nMAF:", MLM$maf %>% round(4), "\nModel:", MLM$model, sep = "")
  
  MLM$link <- paste("http://icci.tau.ac.il:8000/cgi-bin/hgTracks?db=hub_112_Longissima&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=", 
                    "chr", MLM$Chromosome, "%3A" ,MLM$Position, "%2D", MLM$Position,
                    "&hgsid=329_aFWaA0H7kp7ROLdvGJNjECGBoAin", sep="")
  
  
  # merge the this model with previous results
  all <- rbind(GLM,MLM)
  
  
  #################################################################################################################
  #MLMM
  
  # read data of next model
  file <- input$mlmm
  ext <- tools::file_ext(file$datapath)
  MLMM <- read.csv(file$datapath)
  
  MLMM$Chromosome[MLMM$Chromosome == "UN"] <- "Un"
  
  # leave only FDR sig snps
  MLMM <- MLMM[MLMM$FDR_Adjusted_P.values < 0.05 ,]
  
  # make cumulative position the same as the first model
  MLMM <- merge(MLMM, chr_start_pos, by="Chromosome")
  MLMM <- MLMM %>% group_by(Chromosome) %>% mutate(BPcum = Position+tot)
  
  
  #add model column
  MLMM$model <- "MLMM"
  
  MLMM$text <- paste("snp: ", MLMM$SNP, "\nPosition: ", MLMM$Position, 
                     "\nChromosome: ", MLMM$Chromosome, "\nFDR:", -log10(MLMM$FDR_Adjusted_P.values) %>% round(2),
                     "\nMAF:", MLMM$maf %>% round(4), "\nModel:", MLMM$model, sep = "")
  
  MLMM$link <- paste("http://icci.tau.ac.il:8000/cgi-bin/hgTracks?db=hub_112_Longissima&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=", 
                     "chr", MLMM$Chromosome, "%3A" ,MLMM$Position, "%2D", MLMM$Position,
                     "&hgsid=329_aFWaA0H7kp7ROLdvGJNjECGBoAin", sep="")
  
  
  # merge the this model with previous results
  all <- rbind(all,MLMM)
  
  
  #################################################################################################################
  #FarmCPU
  
  # read data of next model
  file <- input$farm 
  ext <- tools::file_ext(file$datapath)
  farm <- read.csv(file$datapath)

  farm$Chromosome[farm$Chromosome == "UN"] <- "Un"
  
  # leave only FDR sig snps
  farm <- farm[farm$FDR_Adjusted_P.values < 0.05 ,]
  
  # make cumulative position the same as the first model
  farm <- merge(farm, chr_start_pos, by="Chromosome")
  farm <- farm %>% group_by(Chromosome) %>% mutate(BPcum = Position+tot)
  
  
  #add model column
  farm$model <- "FramCPU"
  
  farm$text <- paste("snp: ", farm$SNP, "\nPosition: ", farm$Position, 
                     "\nChromosome: ", farm$Chromosome, "\nFDR:", -log10(farm$FDR_Adjusted_P.values) %>% round(2),
                     "\nMAF:", farm$maf %>% round(4), "\nModel:", farm$model, sep = "")
  
  farm$link <- paste("http://icci.tau.ac.il:8000/cgi-bin/hgTracks?db=hub_112_Longissima&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=", 
                     "chr", farm$Chromosome, "%3A" ,farm$Position, "%2D", farm$Position,
                     "&hgsid=329_aFWaA0H7kp7ROLdvGJNjECGBoAin", sep="")
  
  
  # merge the this model with previous results
  all <- rbind(all,farm)
  
  
  #################################################################################################################
  
  #BLINK
  
  # read data of next model
  file <- input$blink 
  ext <- tools::file_ext(file$datapath)
  blink <- read.csv(file$datapath)
  
  blink$Chromosome[blink$Chromosome == "UN"] <- "Un"
  
  # leave only FDR sig snps
  blink <- blink[blink$FDR_Adjusted_P.values < 0.05 ,]
  
  # make cumulative position the same as the first model
  blink <- merge(blink, chr_start_pos, by="Chromosome")
  blink <- blink %>% group_by(Chromosome) %>% mutate(BPcum = Position+tot)
  
  
  #add model column
  blink$model <- "BLINK"
  
  blink$text <- paste("snp: ", blink$SNP, "\nPosition: ", blink$Position, 
                      "\nChromosome: ", blink$Chromosome, "\nFDR:", -log10(blink$FDR_Adjusted_P.values) %>% round(2),
                      "\nMAF:", blink$maf %>% round(4), "\nModel:", blink$model, sep = "")
  
  blink$link <- paste("http://icci.tau.ac.il:8000/cgi-bin/hgTracks?db=hub_112_Longissima&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=", 
                      "chr", blink$Chromosome, "%3A" ,blink$Position, "%2D", blink$Position,
                      "&hgsid=329_aFWaA0H7kp7ROLdvGJNjECGBoAin", sep="")
  
  
  # merge the this model with previous results
  all <- rbind(all,blink)
  
  
  #################################################################################################################
  
  #mark duplication
  n_dup  <- all %>% group_by(Position) %>% filter(n()>1) %>% summarise(N  = n())
  all <- merge(all, n_dup, by = "Position", all.x= T)
  all$N[is.na(all$N)] <- 1
  
  
  # create the graph 
  graph <- observeEvent(input$run, {
  
  # make the data to plot according to user input
#  n_models = input$n_models
#  if (n_models == "all"){
#    dat = all
#  } else {dat = all[all$N == as.integer(n_models) ,]}
  
  
  # plot only snp that apear in ? models
  p <- ggplot(all, aes(x=BPcum, y=-log10(FDR_Adjusted_P.values), text=text, customdata=link, 
                                   shape=(model))) +
    
    # Show all points  (duplicate kmers are circles and uniq are triengels)
    geom_jitter(aes(color=as.factor(Chromosome)), size=3, height = 0.9, width = 0.9) +
    
    
    scale_color_manual(values = rep(c("blue", "orange"), 4 )) +
    

    # custom X axis:
    scale_x_continuous( label = axisdf$Chromosome, breaks= axisdf$center ) +
    
    # remove space between plot area and x axis and set the y axis limits 
    scale_y_continuous(expand = c(0, 0), limits = c(min(-log10(all$FDR_Adjusted_P.values)) - 1, max(-log10(all$FDR_Adjusted_P.values)) + 0.5)) +   
    
    # Add horizontal line
    #geom_hline(yintercept = 5, linetype="dashed", color = "red", size=2) +
    
    # Custom the theme:
    theme_bw() +
    theme( 
      legend.position="none",
      legend.title = element_blank(),
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.title.x = element_blank()
    ) + ylab("FDR [-log10]")
  
  
  p2 <- ggplotly(p, tooltip = c('text'))
  
  p2$x$data[[1]]$customdata <- all$link
  
  js <- "
    function(el) {
    el.on('plotly_click', function(d) {
    var url = d.points[0].customdata;
    //url
    window.open(url);
    });
    }"
  
  
  p3 <- onRender(p2, js)
  
  output$plot <- renderPlotly(plotly::ggplotly(p3))
  
  
  })
})
  
  # download the graph
  output$downloadData <- downloadHandler(
    filename = function() {
      paste(input$dataset, ".html", sep = "")
    },
    content = function(file) {
      htmlwidgets::saveWidget(renderPlotly(plotly::ggplotly(p3)), file)
    })
  
}

shinyApp(ui = ui, server = server)
  
  