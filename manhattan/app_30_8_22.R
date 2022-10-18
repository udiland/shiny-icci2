library(shiny)
library(shinyFiles)
library(dplyr)
library(htmlwidgets)
library(plotly)
library(ggplot2)
library(shinycssloaders)
library(data.table)
options(shiny.maxRequestSize=1000*1024^2)
library(DT)
library(leaflet)
library(mapview)
library(ggthemes)
source("./all2geno.R")
library(gridExtra)
library(data.table)
pdf(NULL)

# set the map parameters
Lat <- 32.33349
Lon <- 34.92786
Zoom <- 7

# read geographical data
geo_data <- read.csv("./geo_data2.csv")

# add noise to show overlapping regions
geo_data$latitude <- jitter(geo_data$latitude, factor = 2)
geo_data$longitude <- jitter(geo_data$longitude, factor = 2)

# make sample icon for the map  # no need when using circles
#iconURL <- "/Users/udila/Documents/ICCI/shiny/test_map/pick.png"
#leafletIcon <- makeIcon(
#  iconUrl = iconURL,
#  iconWidth = 16, iconHeight = 16
#)

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
               shinycssloaders::withSpinner(plotlyOutput('plot'))
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

             # choose if the genotype data is of filtered or pruned
             radioButtons("chr_data", "Genomic data:",
               c("LD pruned" = "pruned",
                 "All SNPs" = "filtered"),
                selected="pruned"),
             
             
             #select chromosome
             selectInput(inputId = "chr", label = "Choose chromosome", choices = c('All',1,2,3,4,5,6,7,'Un'), multiple = F, selected = 'All'),
             
             #select start positions
             numericInput(inputId = "start", label = "Choose start position", 0),
             
             #select end positions
             numericInput(inputId = "end", label = "Choose end position", 1e9),
             
             
             #choose plot type
             radioButtons("plotype", "Choose geno-pheno plot type", c("Violin", "Boxplot"), "Violin"),
             
             # get a y-axis label for the geno-pheon plot
             textInput("YLabelGenePheno","Enter Y axis label (phenotypic measurement)", ""),
             
             # get a title for the geno-pheon plot
             textInput("titleGenePheno","Enter title for geno-pheno plot", ""),
             
             # create a run button fo the phenotype table and plot
             actionButton("runPhenotype", "Run"),
             
             # create a button to download violin/boxplot
             downloadButton("downloadPlot", "Download  plot")
             
             
      ),
      
      column(9,
             # show genotype-phenotype plot
             plotOutput("genopheno"),
      )
      
    ), # close row
    
    # Show snps table
    conditionalPanel(
      condition = "input.runPhenotype > 0",
      withSpinner(DTOutput('tbl'))
    ),
    
    # draw line
    tags$hr(),
    
    fluidRow(
      column(5,
             # show pheno geno table
             DTOutput("tbl2"),
             
             # create a button to download table
             downloadButton("downloadTable", "Download geno-pheon table")
             
      ),
      
      # for debugging
      verbatimTextOutput('x4'),
      
      #DTOutput("tbl3"),
      column(6,
             # show map
             leafletOutput(outputId = "leafletMap"),
             
             # create a button to download the manhattan plot
             downloadButton("downloadMap", "Download Map")
      )
    ) # close row
    
  ) # close fluid page
  
)

server <- function(input, output) {
  
  #set a container for reactive values (to create a reactive value use: values$name, than it will be available through all applications)
  values <- reactiveValues()
  
  # create manhattan plot
  analysis <- observeEvent(input$runManhattan, {
    
    #get the gwas results file:
    file <- input$gwas
    ext <- tools::file_ext(file$datapath)
    gwasResults <- fread(file$datapath)
    
    
    # Select  sig snps define by user
    gwasResults <- gwasResults[P.value < as.numeric(input$pv) ]
    
    
    
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
                      "\nChromosome: ", don$Chromosome, "\nFDR:",don$`FDR_Adjusted_P-values` %>% round(6),
                      "\nMAF:", don$maf %>% round(4), sep = "")
    
    don$link <- paste("http://icci-1.tau.ac.il:8000/cgi-bin/hgTracks?db=hub_112_Longissima&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=",
                      "chr", don$Chromosome, "%3A" ,don$Position, "%2D", don$Position,
                      "&hgsid=329_aFWaA0H7kp7ROLdvGJNjECGBoAin", sep="")
    
    
    p <- ggplot(don, aes(x=BPcum, y=-log10(P.value), text=text, customdata=link, shape=`FDR_Adjusted_P-values` < 0.05)) +
      
      # Show all points  (duplicate kmers are circles and uniq are triengels)
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
      
      # add a line that separate FDR significant snps
      geom_hline(yintercept = -log10(max(don$P.value[don$`FDR_Adjusted_P-values` < 0.05])), linetype="dashed",
                 color = "red", size=1) +
      
      # add title from the user
      ggtitle(input$titleManhattan) +
      theme(plot.title = element_text(hjust = 0.5, size = 22))
    
    
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
      filename = paste(input$titleManhattan,".html", sep = ""),
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
    if(input$chr_data == "pruned"){
	if(input$chr == "Un"){
		read.table("/storage/udiland/gwas/snp/data/geno/all_samples/LD_pruned_100/long_chr8.hmp.txt", sep = "\t",comment.char = "",  head=T)
      
		 #hapmap <- read.table("/Users/udila/Downloads/long_chr8.hmp.txt", sep = '\t', comment.char = "", head=T)
      
		} else if(input$chr == "All"){
		read.table("/storage/udiland/gwas/snp/data/geno/all_samples/LD_pruned_100/long_all_pruned.hmp.txt", sep= "\t", comment.char = "", head=T)
      
		} else(
		read.table(paste("/storage/udiland/gwas/snp/data/geno/all_samples/LD_pruned_100/long_chr",input$chr,".hmp.txt", sep=""), sep = "\t", comment.char = "", head=T)
		)
    }
   if(input$chr_data == "filtered"){
        if(input$chr == "Un"){
                hapmap <- read.table("/storage/udiland/gwas/snp/data/geno/all_samples/long_chr8.hmp.txt", sep = "\t", head=T, comment.char = "")

                } else if(input$chr == "All"){
                        read.table("/storage/udiland/gwas/snp/data/geno/all_samples/long_all_pruned.hmp.txt", sep= "\t", comment.char = "", head=T)

                } else(
                        read.table(paste("/storage/udiland/gwas/snp/data/geno/all_samples/long_chr",input$chr,".hmp.txt", sep=""), sep = "\t", comment.char = "", head=T)
                )
   }
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
    values$selected_pos <- selected_pos[,1:5]

#    req(values$selected_pos)

#    output_tbl <- merge(values$selected_pos, gwas(), by.x="pos", by.y = "Position", all.x=T)
    # bypass the selected positions
    output_tbl <- merge(hapmap()[,1:5], gwas(), by.x="pos", by.y = "Position", all.x=T)    

    output_tbl <- output_tbl[, -c(2,7)]
    
    # make the table as reactive value for later use
    values$tbl <- output_tbl
    
    # display the table
    output$tbl <- renderDT(output_tbl)
    
  })
  
  # save selected rows from snp table
  snp_rows <- reactive(input$tbl_rows_selected)
  
  
  # if only one snp is selected plot the geno_pheno
  observe({
    if (length(snp_rows()) == 1 ){
      
      ## take genotype data
      # get snp position
      position <- req(values$tbl[snp_rows() ,]$pos)
      
      # get snp name
      snp_name <- values$tbl[snp_rows() ,]$SNP
      
      # get genotype data from snp table using the snp position (subset from column number 12..)
      geno <- hapmap()[hapmap()$pos == position ,][c(12 : length(hapmap()))]
      
      geno[2,] <- names(geno)
      
      genotype <- transpose(geno)
      
      
      colnames(genotype) <- c("allele", "accession")
      
      # merge phenotype and genotype data
      res <- merge(pheno(), req(genotype), by="accession")
      
      # add snp name
      res$SNP <- snp_name
      
      # add genotype column by converting the hapmap allele coding
      for(i in 1:dim(res)[1]){
        
        if (res$allele[i] == "A"){
          res$genotype[i] <- "A/A"
          
        }else if (res$allele[i] == "G"){
          res$genotype[i] <- "G/G"
          
        }else if (res$allele[i] == "C"){
          res$genotype[i] <- "C/C"
          
        }else if (res$allele[i] == "T"){
          res$genotype[i] <- "T/T"
          
        }else if (res$allele[i] == "R"){
          res$genotype[i] <- "A/G"
          
        }else if (res$allele[i] == "Y"){
          res$genotype[i] <- "C/T"
          
        }else if (res$allele[i] == "S"){
          res$genotype[i] <- "G/C"
          
        }else if (res$allele[i] == "W"){
          res$genotype[i] <- "A/T"
          
        }else if (res$allele[i] == "K"){
          res$genotype[i] <- "G/T"
          
        }else if (res$allele[i] == "M"){
          res$genotype[i] <- "A/C"
          
        }else if (res$allele[i] == "B"){
          res$genotype[i] <- "C/G/T"
          
        }else if (res$allele[i] == "D"){
          res$genotype[i] <- "A/G/T"
          
        }else if (res$allele[i] == "H"){
          res$genotype[i] <- "A/C/T"
          
        }else if (res$allele[i] == "V"){
          res$genotype[i] <- "A/C/G"
          
        }else if (res$allele[i] == "N"){
          res$genotype[i] <- "N"
          
        } else if (res$allele[i] == "-" || res$allele[i] == ".") {
          res$genotype[i] <- "gap"}
        
      }
      
      # remove the allel column
      res <- res[c("SNP", "accession", "genotype","phenotype")]
      
      # download the geno pheno table
      output$downloadTable <- downloadHandler(
        filename = "table.csv",
        content = function(file) {
          write.csv(res, file, row.names = F)
        })
      
      #for debugging
      #    output$x4 = renderPrint({
      #      s = input$tbl_rows_selected
      #      if (length(s)) {
      #        cat('These rows were selected:\n\n')
      #        cat(s, sep = ', ')
      #      }
      #    })
      
      
      # remove the 'N' genotype for ploting
      res <- res[!res$genotype == 'N' ,]
      
      # add geographical data
      res <- merge(res, geo_data, by = "accession", all.x=TRUE)
      
      # columns to display:
      res <- res[c("SNP", "accession", "genotype","phenotype", "site_general", "het")]
      
      # display the pheno geno data
      rownames(res) <- 1:dim(res)[1]
      output$tbl2 <- renderDT(res)
      
      #save as reactive value
      values$res <- res[c("SNP", "accession", "genotype","phenotype")]
      
      # select type of plot from the user
      geom <- switch(input$plotype, Violin = geom_violin,  Boxplot = geom_boxplot)
      
      # calculate means
      means = aggregate(phenotype ~ genotype, res, mean)
      means$phenotype <- round(means$phenotype ,2)
      
      # plot the genotype phenotype data
      plot <- ggplot(res, aes(x=genotype, y=phenotype, fill=genotype)) + geom() +
        theme_classic()+ ggtitle(input$titleGenePheno) +
        ylab(input$YLabelGenePheno) +
        stat_summary(fun=mean, geom="point", shape=20, size=10, color="darkred") +
        geom_text(data = means, aes(label = phenotype, y = phenotype + 0.5), size = 6, show.legend = FALSE) +
        theme(axis.title.x=element_blank()) + theme(legend.position="none") +
        theme(plot.title = element_text(hjust = 0.5, size = 22))
      
      #plot geno_pheno
      output$genopheno <- renderPlot(plot)
      
      
       output$downloadPlot <- downloadHandler(
          filename = "plot.png",
          content = function(file) {
            ggsave(file, plot)
         })
      
      
      
    } # close if for one snp selection
    
    # if more than one row is selected:
    if (length(snp_rows()) > 1){
      
      ## take genotype data
      # get snp position
      positions <- req(values$tbl[snp_rows() ,]$pos)
      
      # get snp name
      snp_names <- values$tbl[snp_rows() ,]$SNP
      
      # get genotype data from snp table using the snp position (subset from column number 12..)
      geno <- hapmap()[hapmap()$pos %in% positions ,][c(4, 12 : length(hapmap()))]
      
      # transpose and change to data frame
      genotype <- t(geno)
      
      genotype <- as.data.frame(genotype)
      
      #take column names from 1st row
      colnames(genotype) <- genotype[1,]
      
      # remove 1st row
      genotype <- genotype[-1 ,]
      
      # set column of accessions from row names
      genotype$accession <- rownames(genotype)
      
      # change rownames to numbers
      rownames(genotype) <- 1:dim(genotype)[1]
      
      
      # add snp name
      all <- merge(genotype, values$res, by ="accession")
      
      # order by phenotype
      all <- all[order(all$phenotype) ,]
      
      all$accession <- factor(all$accession, levels = unique(all$accession))
      
      # take the phenotype to make a different plot
      pheno_col <- subset(all, select = c(accession, phenotype))
      
      pheno_col <- as.data.table(pheno_col)
      
      pheno_col <- melt(pheno_col)
      
      #take genotype columns
      all <- subset(all, select= -c(genotype, SNP, phenotype))
      
      #snp columns only
      cols <- colnames(all)[!colnames(all) %in% "accession"]
      
      # go over the columns and
      all[cols] <- apply(all[cols], c(1,2), allel2geno)
      
      # change to data.table
      all <- as.data.table(all)
      
      dat_to_plot <- melt(all, id.vars = "accession")
      
      # show table
      output$tbl2 <- renderDT(all)
      
      # plot the accessions genotype
      snp <- ggplot(dat_to_plot, aes(variable, accession)) + geom_raster(aes(fill=value)) + theme_bw() +
        theme(axis.text.y = element_blank(),
              axis.ticks.y=element_blank(),
              axis.title.y = element_blank(),
              axis.title.x=element_blank()) +
        scale_fill_calc() +
        scale_x_discrete(expand=c(0,0))
      
      # plot sorted phenotype values as color
      phenot <- ggplot(pheno_col, aes(variable, accession)) + geom_tile(aes(fill = value)) + theme_bw()+
        theme(axis.title.x=element_blank(),
              
              axis.ticks.x=element_blank(),
              axis.title.y=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank()) +
        scale_x_discrete(expand=c(0,0))
      
      # layout matrix
      lay = rbind(c(1,1,1,1,2),
                  c(1,1,1,1,2),
                  c(1,1,1,1,2))
      
      plt  <- arrangeGrob(snp, phenot, layout_matrix = lay)
      
      output$genopheno <- renderPlot(grid.arrange(snp, phenot, layout_matrix = lay))
      
      
      
      #for debugging
      #output$x4 = renderPrint({length(snp_rows())})
      
      output$downloadPlot <- downloadHandler(
        filename = "plot.png",
        content = function(file) {
          ggsave(file, plt)
                  })
      
      
    } # close if statment
    
  }) # close observe
  
  # download the geno pheno graph
  #  output$downloadPlot <- downloadHandler(
  #    filename = "plot.png",
  #    content = function(file) {
  #      ggsave(file, plot=values$snp_plot)
  
  #    })
  
  
  # create a table of selected rows from geno_pheno and the geographical data
  geno_pheno_geo <- reactive(req(merge(req(values$res[input$tbl2_rows_selected ,]), geo_data, by="accession", all.x=T)))
  
  #show the table to user
  output$tbl3 <- renderDT(geno_pheno_geo())
  
  # creat color pallet for m phenotype
  col <- reactive(colorFactor(palette = 'RdYlGn', unique(geno_pheno_geo()$phenotype)))
  
  #create a leaflet map reactive object
  map <- reactive({
    leaflet(data = geno_pheno_geo()) %>%
      setView(lat = Lat, lng = Lon, zoom = Zoom) %>%
      addTiles() %>%
      addProviderTiles(providers$Esri.WorldStreetMap)
  })
  
  # show the map and add data to each point
  output$leafletMap <- renderLeaflet(req(map()) %>% addCircleMarkers(data = geno_pheno_geo(),
                                                                     lng = geno_pheno_geo()$longitude,
                                                                     lat = geno_pheno_geo()$latitude,
                                                                     radius = 5,
                                                                     popup = paste("Accession:", geno_pheno_geo()$accession, "<br>",
                                                                                   "Phenotype:", geno_pheno_geo()$phenotype, "<br>",
                                                                                   "Site:", geno_pheno_geo()$site_general, "<be>"),
                                                                     label = geno_pheno_geo()$site_general,
                                                                     color = col()(geno_pheno_geo()$phenotype),
                                                                     fillOpacity = 1) %>%
                                       addLegend("bottomright", pal = col(),
                                                 values = unique(geno_pheno_geo()$phenotype),
                                                 opacity = 1))  #%>% leafletOutput(height=500) )#title = input$titleGenePheno)#,

  
  
  # download the map
  map_to_save <- reactive(map() %>% addCircleMarkers(data = geno_pheno_geo(),
                                                     lng = geno_pheno_geo()$longitude,
                                                     lat = geno_pheno_geo()$latitude,
                                                     radius=5,
                                                     color = col()(geno_pheno_geo()$phenotype),
                                                     fillOpacity = 1) %>%
                            addLegend("bottomright", pal = col(),
                                      values = unique(geno_pheno_geo()$phenotype),
                                      opacity = 1)
  )
  #title = input$titleGenePheno)
  
  
  output$downloadMap <- downloadHandler(
    filename = "Map.pdf"
    
    , content = function(file) {
      mapshot(map_to_save()
              , file = file
              #, cliprect = "viewport" # the clipping rectangle matches the height & width from the viewing port
              #, selfcontained = FALSE # when this was not specified, the function for produced a PDF of two pages: one of the leaflet map, the other a blank page.
      )} # end of content() function
  ) # end of downloadHandler() funct
  
  
}



shinyApp(ui = ui, server = server)


