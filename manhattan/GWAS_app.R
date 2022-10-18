#devtools::install_github('ropensci/plotly')
#install.packages("shinycssloaders")
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
  fileInput(inputId = "gwas", label = "Upload GAPIT results"),
  
  # sellect p-value cutoff
  selectInput("pv", "P-value cutoff", c("0.05", "0.01", "0.001"), selected = "0.001"),
  
  #select chromosome
  #selectInput(inputId = "chr", label = "Choose chromosome", choices = c(1,2,3,4,5,6,7,'Un', 'All'), multiple = F, selected='All'),
  

  # create a button to download
  downloadButton("downloadData", "Download"),
  
  # create a run button
  actionButton("run", "Run"),
  
  # create a display for the graphs output with spinner
  mainPanel(
	conditionalPanel(
	condition = "input.run > 0",
	shinycssloaders::withSpinner(plotlyOutput('plot')))
  )

)



server <- function(input, output) {
    
library(plotly) 

analysis <- observeEvent(input$run, {
    
    #get the phenotpye file:
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
                      "\nChromosome: ", don$Chromosome, "\nPV score:", -log10(don$P.value) %>% round(2),
                      "\nMAF:", don$maf %>% round(4), sep = "")
    
    don$link <- paste("http://icci.tau.ac.il:8000/cgi-bin/hgTracks?db=hub_112_Longissima&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=", 
                      "chr", don$Chromosome, "%3A" ,don$Position, "%2D", don$Position,
                      "&hgsid=329_aFWaA0H7kp7ROLdvGJNjECGBoAin", sep="")
    
    
    p <- ggplot(don, aes(x=BPcum, y=-log10(P.value), text=text, customdata=link, shape=FDR_Adjusted_P.values < 0.05)) +

      # Show all points  (duplicate kmers are circles and uniq are triengels)
      geom_jitter(aes(color=as.factor(Chromosome)), size=2) +

      scale_color_manual(values = rep(c("blue", "orange"), 4 )) +

      # set the shapes
      scale_shape_manual(values = c(16,8)) +

      # custom X axis:
      scale_x_continuous( label = axisdf$Chromosome, breaks= axisdf$center ) +

      # remove space between plot area and x axis and set the y axis limits
      scale_y_continuous(expand = c(0, 0), limits = c(min(-log10(don$P.value)), max(-log10(don$P.value)) + 0.5)) +

      # Add horizontal line
      #geom_hline(yintercept = 5, linetype="dashed", color = "red", size=2) +

      # Custom the theme:
      theme_bw() +
      theme(
        legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title.x=element_blank()
        )
    
    
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
    
    
    man_plot <- reactive(onRender(p2, js))
    
  

    output$plot <- renderPlotly(plotly::ggplotly(man_plot()))

  })


  # download the graph
  output$downloadData <- downloadHandler(
    filename = function() {
      paste(input$gwas, ".html", sep = "")
    },
    content = function(file) {
      htmlwidgets::saveWidget(plotly::ggplotly(man_plot()), file)
    })
  

}


shinyApp(ui = ui, server = server)
