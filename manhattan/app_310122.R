library(shiny)
library(shinyFiles)
library(dplyr)
library(htmlwidgets)
library(plotly)
library(ggplot2)
library(shinycssloaders)
library(data.table)
options(shiny.maxRequestSize=1000*1024^2)
library("DT")



ui <- fluidPage(
  
  # add university logo and try to style with css file (it didnt work but didnt fail...)
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

tags$body(
  
        h1("Plot Manhhatan"),
  
        fluidRow(
        column(3,
        # create a selection pan for input files
          fileInput(inputId = "gwas", label = "Upload GAPIT results"),
  
  
        # get a title for the manhattan plot
          textInput("titleManhattan","Enter title for Manhattan plot", ""),
  
        # sellect p-value cutoff
          selectInput("pv", "P-value cutoff", c("0.05", "0.01", "0.001"), selected = "0.001"),
  
        # create a run button fo the manhattan plot
          actionButton("runManhattan", "Run"),
  
        # create a button to download the manhattan plot
          downloadButton("downloadData", "Download Manhattan plot")
        ),
  
      # create a display for the graphs output with spinner
      column(9,   
             
            conditionalPanel(
            condition = "input.runManhattan > 0",
            withSpinner(plotlyOutput('plot'))
              )
            )
      ), # close row
  
      # add a line separa tor section 2
      tags$hr(),
  
      # Title section 2
      h1("Plot Phenotype"),
      fluidRow(
      column(3,
  
        # upload phenotype table
        fileInput(inputId = "phenotype", label = "Upload phenotype table"),
  
  
        #select chromosome
        selectInput(inputId = "chr", label = "Choose chromosome", choices = c(1,2,3,4,5,6,7,'Un'), multiple = F, selected = 'Un'),
  
        #select start positions
        numericInput(inputId = "start", label = "Choose start position", 0),
  
        #select end positions
        numericInput(inputId = "end", label = "Choose end position", 1e9),
  

        #choose plot type
        radioButtons("plotype", "Choose geno-pheno plot type", c("Violin", "Boxplot"), "violin"),
  
  
        # get a title for the geno-pheon plot
        textInput("titleGenePheno","Enter title for geno-pheno plot", ""),
  
        # create a run button fo the phenotype table and plot
        actionButton("runPhenotype", "Run"),
  
        # create a button to download table
        downloadButton("downloadTable", "Download geno-pheon table"),
  
        # create a button to download figure
        downloadButton("downloadPlot", "Download geno-pheon plot"),
  
        ),
  
    column(9,
    # show genotype-phenotype plot
    plotOutput("genopheno"),  
        )
    
  ), # close row
  
  # Show table in next line
  fluidRow(
    conditionalPanel(
      condition = "input.runPhenotype > 0",
      withSpinner(DTOutput('tbl'))
    ),
  
  ), #close row
  
  fluidRow(
  # show pheno geno table
  DTOutput("tbl2"),
  
  ),
  
  
  # for debugging
  verbatimTextOutput('x4')
  

) # close body tag

) # close fluid page


server <- function(input, output) {
  
  #set a container for reactive values (to create a reactive value use: values$name, than it will be available through all applications)
  values <- reactiveValues()
  
  # create manhattan plot
  analysis <- observeEvent(input$runManhattan, {
    
    #get the gwas results file:
    file <- input$gwas
    ext <- tools::file_ext(file$datapath)
    gwasResults <- read.csv(file$datapath)
    
    
    # Select  sig snps define by user
    gwasResults <- gwasResults[gwasResults['P.value'] < as.numeric(input$pv) ,]
    
    
    
    # change "UN" to "Un"
    gwasResults$Chromosome[gwasResults$Chromosome == "UN"] <- "Un"
    
    
    don <- gwasResults %>%
      
      # Compute chromosome size
      group_by(Chromosome) %>%
      summarise(chr_len=max(Position)) %>%
      
      # Calculate cumulative position of each chromosome
      mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
      select(-chr_len) %>%
      
      # Add this info to the initial dataset
      left_join(gwasResults, ., by=c("Chromosome"="Chromosome")) %>%
      
      # Add a cumulative position of each SNP
      arrange(Chromosome, Position) %>%
      mutate( BPcum=Position+tot)
    
    
    #Then we need to prepare the X axis. Indeed we do not want to display the cumulative position of SNP in bp, but just show the chromosome name instead.
    
    axisdf = don %>% group_by(Chromosome) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
    
    
    # use plotly to add data to points
    
    
    don$text <- paste("snp: ", don$rs, "\nPosition: ", don$Position,
                      "\nChromosome: ", don$Chromosome, "\nFDR:",don$FDR_Adjusted_P.values %>% round(6),
                      "\nMAF:", don$maf %>% round(4), sep = "")
    
    don$link <- paste("http://icci.tau.ac.il:8000/cgi-bin/hgTracks?db=hub_112_Longissima&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=",
                      "chr", don$Chromosome, "%3A" ,don$Position, "%2D", don$Position,
                      "&hgsid=329_aFWaA0H7kp7ROLdvGJNjECGBoAin", sep="")
    
    
    p <- ggplot(don, aes(x=BPcum, y=-log10(P.value), text=text, customdata=link, shape=FDR_Adjusted_P.values < 0.05)) +
      
      # Show all points 
      geom_jitter(aes(color=as.factor(Chromosome)), size=2) +
      
      scale_color_manual(values = rep(c("blue", "orange"), 4 )) +
      
      # set the shapes
      scale_shape_manual(values = c(16,8)) +
      
      # custom X axis:
      scale_x_continuous( label = axisdf$Chromosome, breaks= axisdf$center ) +
      
      # remove space between plot area and x axis and set the y axis limits
      scale_y_continuous(expand = c(0, 0), limits = c(min(-log10(don$P.value)), max(-log10(don$P.value)) + 0.5)) +
      
        
      # Custom the theme:
      theme_bw() +
      theme(
        legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title.x=element_blank()
      ) +
      # add a line that separate FDR significant snps (take the p-value coresponding to the mninmal FDR > 0.05)
      geom_hline(yintercept = -log10(min(don$P.value[don$FDR_Adjusted_P.values > 0.05])), linetype="dashed", color = "red", size=1) +
      # add title from user
      ggtitle(input$titleManhattan) +
      # position title
      theme(plot.title = element_text(hjust = 0.5)) 
   
            
    
    p2 <- ggplotly(p, tooltip = c('text'))
    
    p2$x$data[[1]]$customdata <- don$link
    
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
    
    
    
    # download the graph
    output$downloadData <- downloadHandler(
      filename = "Manhattan.html",
      content = function(file) {
        htmlwidgets::saveWidget(plotly::ggplotly(p3), file)
      })
    
  })
  

pheno <- req(reactive({
  #get the phenotpye file:
  file <- input$phenotype
  ext <- tools::file_ext(file$datapath)
  pheno_tbl <- read.table(file$datapath, sep='\t', head=T)
  
  #give new column name
  colnames(pheno_tbl)[1] <- "accession"
  colnames(pheno_tbl)[2] <- "phenotype"
  pheno_tbl
}))



hapmap <- req(reactive({
    # get genotype data
  req(if(input$chr == "Un"){
    hapmap <- read.table("/storage/udiland/gwas/snp/data/geno/all_samples/LD_pruned_100/long_chr8.hmp.txt", sep = '\t', comment.char = "", head=T)
    
   # read.table("/Users/udila/Downloads/long_chr8.hmp.txt", sep = '\t', comment.char = "", head=T)
  
  }
  else(
    read.table(paste("/storage/udiland/gwas/snp/data/geno/all_samples/LD_pruned_100/long_chr",input$chr,".hmp.txt", sep = ""), sep = '\t', comment.char = "", head=T)
  ))
  
}))

gwas <- reactive({
  #get the gwas results file:
  file <- input$gwas
  ext <- tools::file_ext(file$datapath)
  gwasResults <- read.csv(file$datapath)
  
  # change "UN" to "Un"
  gwasResults$Chromosome[gwasResults$Chromosome == "UN"] <- "Un"
  gwasResults
})
  


observeEvent(input$runPhenotype, {

  
  #show the table of selected range range
  selected_pos <- hapmap()[hapmap()$pos >= input$start & hapmap()$pos <= input$end ,]
  
  # get only first 5 columns and merge with gwas results
  selected_pos <- selected_pos[,1:5]
  
  
  output_tbl <- merge(selected_pos, gwas(), by.x="pos", by.y = "Position", all.x=T)
  
  output_tbl <- output_tbl[, -c(2,7)]
  
  # display the table
  output$tbl <- renderDT(output_tbl, selection = 'single')
  
  # make the table as reactive value
  values$tbl <- output_tbl
  
  })
  

  output$genopheno <- renderPlot({
    

    # take genotype data
    position <- req(values$tbl[input$tbl_rows_selected ,]$pos)
    
    
    geno <- hapmap()[hapmap()$pos == position ,][c(12 : length(hapmap()))]
    
    geno[2,] <- names(geno)
    
    genotype <- transpose(geno)
    
    #genotype <- as.data.frame(genotype)
    
    colnames(genotype) <- c("genotype", "accession")
    
    
    res <- merge(pheno(), genotype, by="accession")
    
    # add allele column by converting the hapmap genotype coding
    for(i in 1:dim(res)[1]){
      
      if (res$genotype[i] == "A"){
        res$allele[i] <- "A/A"
        
      }else if (res$genotype[i] == "G"){
        res$allele[i] <- "G/G"
        
      }else if (res$genotype[i] == "C"){
        res$allele[i] <- "C/C"
        
      }else if (res$genotype[i] == "T"){
        res$allele[i] <- "T/T"
        
      }else if (res$genotype[i] == "R"){
        res$allele[i] <- "A/G"
        
      }else if (res$genotype[i] == "Y"){
        res$allele[i] <- "C/T"
        
      }else if (res$genotype[i] == "S"){
        res$allele[i] <- "G/C"
        
      }else if (res$genotype[i] == "W"){
        res$allele[i] <- "A/T"
        
      }else if (res$genotype[i] == "K"){
        res$allele[i] <- "G/T" 
        
      }else if (res$genotype[i] == "M"){
        res$allele[i] <- "A/C"
        
      }else if (res$genotype[i] == "B"){
        res$allele[i] <- "C/G/T"
        
      }else if (res$genotype[i] == "D"){
        res$allele[i] <- "A/G/T"
        
      }else if (res$genotype[i] == "H"){
        res$allele[i] <- "A/C/T"
        
      }else if (res$genotype[i] == "V"){
        res$allele[i] <- "A/C/G"
        
      }else if (res$genotype[i] == "N"){
        res$allele[i] <- "N"
        
      } else if (res$genotype[i] == "-" || res$genotype[i] == ".") {
        res$allele[i] <- "gap"}
      
    }
    
    
    # download the geno pheno table
    output$downloadTable <- downloadHandler(
      filename = "table.csv",
      content = function(file) {
        write.csv(res, file, row.names = F)
      })
    
    #colnames(res)[2] <- "genotype"
    
#for debugging
#    output$x4 = renderPrint({
#      s = input$tbl_rows_selected
#      if (length(s)) {
#        cat('These rows were selected:\n\n')
#        cat(s, sep = ', ')
#      }
#    })
    
    
# display the pheno geno data
output$tbl2 <- renderDT(res)
  

  # remove the 'N' allele for ploting
  res <- res[!res$allele == 'N' ,]
  
  
  geom <- switch(input$plotype, Violin = geom_violin,  Boxplot = geom_boxplot)
    
  # plot the genotype phenotype data
  values$plot <- ggplot(res, aes(x=allele, y=phenotype, fill=allele)) + geom() + theme_classic()+ ggtitle(input$titleGenePheno)
  
  values$plot
  
  
    })
  
  # download the graph
  output$downloadPlot <- downloadHandler(
    filename = "plot.png",
    content = function(file) {
      ggsave(file, plot=values$plot)
      
    })
  
  }



shinyApp(ui = ui, server = server)

