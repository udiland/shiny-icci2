library(shiny)
library(shinyFiles)
library(dplyr)
library(htmlwidgets)
library(plotly)
library(ggplot2)
library(shinycssloaders)
library("DT")
options(shiny.maxRequestSize=2000*1024^2)



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
    
    tags$p(class="universisy-name", "Tel Aviv Universisy")
    
    
  ),
  
  tags$hr(),
    
    h1("Find common SNPs"),
  
  # create a selection pan for input files
  fileInput(inputId = "glm", label = "Upload GLM results", accept = ".csv"),
  
  fileInput(inputId = "mlm", label = "Upload MLM results", accept = ".csv"),
  
  fileInput(inputId = "mlmm", label = "Upload MLMM results", accept = ".csv"),
  
  fileInput(inputId = "farm", label = "Upload FarmCPU results", accept = ".csv"),
  
  fileInput(inputId = "blink", label = "Upload BLINK results", accept = ".csv"),
  
  selectInput("n_models", "SNP common to number of models", c("any", ">1", ">2", ">3", ">4"), selected = "any"),
  
  # get a title for the manhattan plot
  textInput("title","Enter plot title", ""),
  
  
  # create a button to download
  downloadButton("downloadPlot", "Download Plot"),
  
  
  # create a run button
  actionButton("run", "Run"),
  
  
  mainPanel(
    conditionalPanel(
      condition = "input.run > 0",
      withSpinner(plotlyOutput("plot")))
  ),
  
  # show results table
  DTOutput("res"),
  
  # create a button to download the table
  downloadButton("downloadData", "Download Table"),
  
  verbatimTextOutput('message')
  
)

server <- function(input, output) {
  
  values <- reactiveValues()
  
  gwas_res <- observeEvent(input$run, {
    all <- data.frame()
    
    
    file <- input$glm
    ext <- tools::file_ext(file$datapath)
    try(glm <- read.csv(file$datapath), silent = T)
    
    #glm <- read.csv("/Users/udila/Documents/ICCI/gwas_snp/presentation/SR/all_chr/oneFile/GAPIT.GLM.phenotype.GWAS.Results.csv")
    if (exists("GLM")){
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
    
    if (dim(GLM)[1] > 0){
    #add model column
    GLM$model <- "GLM"
    
    GLM$text <- paste("snp: ", GLM$SNP, "\nPosition: ", GLM$Position,
                      "\nChromosome: ", GLM$Chromosome, "\nFDR:", -log10(GLM$FDR_Adjusted_P.values) %>% round(2),
                      "\nMAF:", GLM$maf %>% round(4), "\nModel:", GLM$model, sep = "")
    
    GLM$link <- paste("http://icci-1.tau.ac.il:8000/cgi-bin/hgTracks?db=hub_112_Longissima&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=",
                      "chr", GLM$Chromosome, "%3A" ,GLM$Position, "%2D", GLM$Position,
                      "&hgsid=329_aFWaA0H7kp7ROLdvGJNjECGBoAin", sep="")
    }
    
      all <- rbind(all,GLM)
    }
    #####################################################################################################################
    # read data of next model
    
    file <- input$mlm
    ext <- tools::file_ext(file$datapath)
    try(MLM <- read.csv(file$datapath), silent = T)
    
    #MLM <- read.csv("/Users/udila/Documents/ICCI/gwas_snp/presentation/SR/all_chr/oneFile/GAPIT.MLM.phenotype.GWAS.Results.csv")
    if (exists("GLM")){
    
    MLM$Chromosome[MLM$Chromosome == "UN"] <- "Un"
    
    # leave only FDR sig snps
    MLM <- MLM[MLM$FDR_Adjusted_P.values <= 0.05 ,]
    
    
    if (dim(MLM)[1] > 0){
    # make cumulative position the same as the first model
    
    MLM <- merge(MLM, chr_start_pos, by="Chromosome")
    MLM <- MLM %>% group_by(Chromosome) %>% mutate(BPcum = Position+tot)
    
    #add model column
    MLM$model <- "MLM"
    
    MLM$text <- paste("snp: ", MLM$SNP, "\nPosition: ", MLM$Position,
                      "\nChromosome: ", MLM$Chromosome, "\nFDR:", -log10(MLM$FDR_Adjusted_P.values) %>% round(2),
                      "\nMAF:", MLM$maf %>% round(4), "\nModel:", MLM$model, sep = "")
    
    MLM$link <- paste("http://icci-1.tau.ac.il:8000/cgi-bin/hgTracks?db=hub_112_Longissima&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=",
                      "chr", MLM$Chromosome, "%3A" ,MLM$Position, "%2D", MLM$Position,
                      "&hgsid=329_aFWaA0H7kp7ROLdvGJNjECGBoAin", sep="")
    }
    
    # merge the this model with previous results

      all <- rbind(all,MLM)
    }else{
      mlm <- read.csv(file$datapath)
      mlm <- mlm[mlm['P.value'] < 0.01 ,]
      
      mlm$Chromosome[mlm$Chromosome == "UN"] <- "Un"
      
      MLM <- mlm %>%
        
        # Compute chromosome size
        group_by(Chromosome) %>%
        summarise(chr_len=max(Position)) %>%
        
        # Calculate cumulative position of each chromosome
        mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
        select(-chr_len) %>%
        
        # Add this info to the initial dataset
        left_join(mlm, ., by=c("Chromosome"="Chromosome")) %>%
        
        # Add a cumulative position of each SNP
        arrange(Chromosome, Position) %>%
        mutate( BPcum=Position+tot)
      
      axisdf = MLM %>% group_by(Chromosome) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
      
      
      #use this calculation of chromosome start position for the other positions
      chr_start_pos <- MLM %>% group_by(Chromosome) %>% summarise(tot = min(tot))
      
      #leave only FDR sig snp
      MLM <- MLM[MLM$FDR_Adjusted_P.values < 0.05 ,]
      
      if (dim(MLM)[1] > 0){
        #add model column
        MLM$model <- "MLM"
        
        MLM$text <- paste("snp: ", MLM$SNP, "\nPosition: ", MLM$Position,
                          "\nChromosome: ", MLM$Chromosome, "\nFDR:", -log10(MLM$FDR_Adjusted_P.values) %>% round(2),
                          "\nMAF:", MLM$maf %>% round(4), "\nModel:", MLM$model, sep = "")
        
        MLM$link <- paste("http://icci-1.tau.ac.il:8000/cgi-bin/hgTracks?db=hub_112_Longissima&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=",
                          "chr", MLM$Chromosome, "%3A" ,MLM$Position, "%2D", MLM$Position,
                          "&hgsid=329_aFWaA0H7kp7ROLdvGJNjECGBoAin", sep="")
      }
      
      all <- rbind(all,MLM)
    }
    
    #################################################################################################################
    #MLMM
    
    # read data of next model
    file <- input$mlmm
    ext <- tools::file_ext(file$datapath)
    try(MLMM <- read.csv(file$datapath))
    
    #MLMM <- read.csv("/Users/udila/Documents/ICCI/gwas_snp/presentation/SR/all_chr/oneFile/GAPIT.MLMM.phenotype.GWAS.Results.csv")
    if (exists("MLMM")){
    
    MLMM$Chromosome[MLMM$Chromosome == "UN"] <- "Un"
    
    # leave only FDR sig snps
    MLMM <- MLMM[MLMM$FDR_Adjusted_P.values <= 0.05 ,]
    
    # make cumulative position the same as the first model
    MLMM <- merge(MLMM, chr_start_pos, by="Chromosome")
    MLMM <- MLMM %>% group_by(Chromosome) %>% mutate(BPcum = Position+tot)
    
    
    #add model column
    MLMM$model <- "MLMM"
    
    MLMM$text <- paste("snp: ", MLMM$SNP, "\nPosition: ", MLMM$Position,
                       "\nChromosome: ", MLMM$Chromosome, "\nFDR:", -log10(MLMM$FDR_Adjusted_P.values) %>% round(2),
                       "\nMAF:", MLMM$maf %>% round(4), "\nModel:", MLMM$model, sep = "")
    
    MLMM$link <- paste("http://icci-1.tau.ac.il:8000/cgi-bin/hgTracks?db=hub_112_Longissima&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=",
                       "chr", MLMM$Chromosome, "%3A" ,MLMM$Position, "%2D", MLMM$Position,
                       "&hgsid=329_aFWaA0H7kp7ROLdvGJNjECGBoAin", sep="")
    
    
    # merge the this model with previous results
      all <- rbind(all,MLMM)
    }
    #################################################################################################################
    #FarmCPU
    
    # read data of next model
    file <- input$farm
    ext <- tools::file_ext(file$datapath)
    try(farm <- read.csv(file$datapath))
    
    if (exists("farm")){
    
    farm$Chromosome[farm$Chromosome == "UN"] <- "Un"
    
    #farm <- read.csv("/Users/udila/Documents/ICCI/gwas_snp/presentation/SR/all_chr/oneFile/GAPIT.FarmCPU.phenotype.GWAS.Results.csv")
    
    
    # leave only FDR sig snps
    farm <- farm[farm$FDR_Adjusted_P.values <= 0.05 ,]
    
    # make cumulative position the same as the first model
    farm <- merge(farm, chr_start_pos, by="Chromosome")
    farm <- farm %>% group_by(Chromosome) %>% mutate(BPcum = Position+tot)
    
    
    #add model column
    farm$model <- "FramCPU"
    
    farm$text <- paste("snp: ", farm$SNP, "\nPosition: ", farm$Position,
                       "\nChromosome: ", farm$Chromosome, "\nFDR:", -log10(farm$FDR_Adjusted_P.values) %>% round(2),
                       "\nMAF:", farm$maf %>% round(4), "\nModel:", farm$model, sep = "")
    
    farm$link <- paste("http://icci-1.tau.ac.il:8000/cgi-bin/hgTracks?db=hub_112_Longissima&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=",
                       "chr", farm$Chromosome, "%3A" ,farm$Position, "%2D", farm$Position,
                       "&hgsid=329_aFWaA0H7kp7ROLdvGJNjECGBoAin", sep="")
    
    
    # merge the this model with previous results
      all <- rbind(all,farm)
    }
    #################################################################################################################
    
    #BLINK
    
    # read data of next model
    file <- input$blink
    ext <- tools::file_ext(file$datapath)
    try(blink <- read.csv(file$datapath))
    
    #blink <- read.csv("/Users/udila/Documents/ICCI/gwas_snp/presentation/SR/all_chr/oneFile/GAPIT.Blink.phenotype.GWAS.Results.csv")
    if (exists("blink")){
    
    blink$Chromosome[blink$Chromosome == "UN"] <- "Un"
    
    # leave only FDR sig snps
    blink <- blink[blink$FDR_Adjusted_P.values <= 0.05 ,]
    
    # make cumulative position the same as the first model
    blink <- merge(blink, chr_start_pos, by="Chromosome")
    blink <- blink %>% group_by(Chromosome) %>% mutate(BPcum = Position+tot)
    
    
    #add model column
    blink$model <- "BLINK"
    
    blink$text <- paste("snp: ", blink$SNP, "\nPosition: ", blink$Position,
                        "\nChromosome: ", blink$Chromosome, "\nFDR:", -log10(blink$FDR_Adjusted_P.values) %>% round(2),
                        "\nMAF:", blink$maf %>% round(4), "\nModel:", blink$model, sep = "")
    
    blink$link <- paste("http://icci-1.tau.ac.il:8000/cgi-bin/hgTracks?db=hub_112_Longissima&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=",
                        "chr", blink$Chromosome, "%3A" ,blink$Position, "%2D", blink$Position,
                        "&hgsid=329_aFWaA0H7kp7ROLdvGJNjECGBoAin", sep="")
    
    
    # merge the this model with previous results
      all <- rbind(all,blink)
    }
    #################################################################################################################
    
    #mark duplication
    n_dup  <- all %>% group_by(Position) %>% filter(n()>1) %>% summarise(N_models  = n())
    all <- merge(all, n_dup, by = "Position", all.x= T)
    all$N_models[is.na(all$N_models)] <- 1
    
    # select data
    if (input$n_models == "any"){
      dat <- all
    }
    if (input$n_models == ">1"){
      dat <- all[all$N_models > 1 ,]
    }
    
    if (input$n_models == ">2"){
      dat <- all[all$N_models > 2 ,]
    }
    
    if (input$n_models == ">3"){
      dat <- all[all$N_models > 3 ,]
    }
    
    if (input$n_models == ">4"){
      dat <- all[all$N_models > 4 ,]
    }
    
    # if the table has results plot them otherwise print message
    if (dim(dat)[1] != 0) {
    # create the graph
    p <- ggplot(dat, aes(x=BPcum, y=-log10(FDR_Adjusted_P.values), text=text, 
                         customdata=link, shape=model, color=as.factor(Chromosome)))+
      
      geom_jitter(size=3, height = 0.5, width = 0.9) +
      
      scale_color_manual(values = rep(c("blue", "orange"), 4 )) +
      
      # custom X axis:
      scale_x_continuous( label = axisdf$Chromosome, breaks= axisdf$center ) +
      
      # remove space between plot area and x axis and set the y axis limits
      scale_y_continuous(expand = c(0, 0), limits = c(min(-log10(dat$FDR_Adjusted_P.values)) - 0.5, max(-log10(dat$FDR_Adjusted_P.values)) + 0.5)) +
      
      theme_bw() +
      theme(
        #legend.position="none",
        legend.title = element_blank(),
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title.x = element_blank()
      ) +
      ylab("FDR [-log10]") +
      ggtitle(input$title) +
      theme(plot.title = element_text(hjust = 0.5, size = 22)) +
      guides(color='none')
     
    
    
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
    
    
    output$plot <- renderPlotly(ggplotly(p3))
    
    } else{
            output$plot <- renderPlot(ggplot())
            output$message = renderPrint({
            "Nothing to show" })
          }
    
    
    output$downloadPlot <- downloadHandler(
      filename = paste(input$title,".html", sep = ""),
      content = function(file) {
        saveWidget(ggplotly(p3), file)
      })
    
    # prepare the table for display
    col_rm <- c("tot", "BPcum", "text", "link")
    
    dat <- dat[, !(names(dat) %in% col_rm)]
    
    output$res <-  renderDT(dat)
    
    # add the table to reactive values 
    values$tbl <- dat
  
  
  
 output$downloadData <- downloadHandler(
   filename = paste(input$title, ".csv", sep = ""),
   content = function(file) {
     write.csv(values$tbl, file, row.names = F)
     })
  })
 }

shinyApp(ui = ui, server = server)

