library(shiny)
library(shinyFiles)
library(htmlwidgets)
library(ggplot2)
library(GenomicFeatures)
library(genetics)
library(ggbio)
library(DT)
library(data.table)
library(dplyr)
options(shiny.maxRequestSize=1000*1024^2)


# read genes data (table name is: longissima_genes)
load("./genes_longissima.Rdata")

# get a GRangesList object of all genes (used to plot genes)
load("./long_genes_grlist.R")

# import functions: getResults, read.genes, makeLDplot, get_loci_tbl
source("./utils.R")

source("./get_loci_tbls.R")


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
  tags$h1("Plot genes associated with significant SNPs"),
  tags$hr(),

  tags$body(

  fluidRow(
    column(3,
           
  # create a selection pan for input files
  shinyDirButton(id = "folder", label = "Choose folder", title = "Choose GAPIT results folder"),

  #select gwas threshold
  numericInput(inputId = "pv", label = "GWAS threshold", 1e-3),

  #select lod score threshold
  numericInput(inputId = "lod", label = "Lod threshold", 5.3),

  #select MAF threshold
  numericInput(inputId = "maf", label = "MAF", 5e-2),

  #select locus width
  numericInput(inputId = "locus", label = "Locus width", 1e4),
  

  # create a run button for creating the loci table
  actionButton("get_table", "Create loci table"),
  
  verbatimTextOutput(outputId = "x4")
  
    ),

column(9,
       
       conditionalPanel(
         condition = "input.get_table > 0",
         # show table
         shinycssloaders::withSpinner(DTOutput("table"))

                        )
     )
  ),


tags$hr(),
  fluidRow(column(3,
                  # choose if the genotype data is of filtered or pruned
                  radioButtons("chr_data", "Genomic data:",
                               choices=c("LD pruned" = "pruned",
                                 "All SNPs" = "filtered"),
                               selected="All SNPs"),


                  #select chromosome
                  selectInput(inputId = "chr", label = "Choose chromosome", choices = c(1,2,3,4,5,6,7,'Un'), multiple = F, selected = 5),
                  # choose lod filter
                  numericInput(inputId = "lodfilter", label = "Lod filter", 0),
                  # select lod max value
                  numericInput(inputId = "lodmax", label = "Lod max", 17),
                  # base pairs up and down stream of the genes of interest
                  numericInput(inputId = "winup", label = "WinUp", 500),
                  numericInput(inputId = "windown", label = "WinDown", 500),
                  # write plot title
                  textInput(inputId = "plot_title", label = "Enter plot title", value = ""),
                  # create a run button for creating the loci plot
                  actionButton("get_plot", "Create plot"),
                  # create a button to download the plot
                  downloadButton("downloadPlot", label = "Download plot")
                  

    ),
          column(6, 
                 conditionalPanel(
                   condition = "input.get_plot > 0",
                   # show table
                   shinycssloaders::withSpinner(plotOutput("loci"))
                                  )
                  )
                 
    ) # close row
  ) # close tag body
)

  
server <- function(input, output) {
  
  values <- reactiveValues()
  
  volumes = getVolumes()
  
  observe({
    shinyDirChoose(input, 'folder', roots=volumes(), filetypes=c("","csv"))
    if(!is.null(input$folder)){
      
      values$InputDir <- parseDirPath(volumes, input$folder)
      values$files <- list.files(values$InputDir, pattern="*.csv")
    }
  })
  
  observeEvent(input$get_table, {
  
  # read all gwas data in selected folder
  gwas <- getResults(values$InputDir, p_value_thresh = input$pv)
  
  # create a table of loci 
  genes_tbl <- get_loci_tbl(gwas, longissima_genes, lod.thresh = input$lod, locus.range = input$locus, maf.thr = input$maf)[[2]]
  
  if (identical(genes_tbl, data.table())){
      output$x4 <- renderText('No genes in range, try to increase range or lower LOD')
      output$table <- renderDT(genes_tbl)
      }else{
  
  # browser()
  
  # stack the columns of Run id
  genes_tbl <- data.frame(genes_tbl[,c(1:15)], stack(genes_tbl[16:ncol(genes_tbl)]))

  genes_tbl <- genes_tbl[!(genes_tbl$values == '.') ,]

  # group bi loci id (this table will be displaied)
  loci_tbl <- genes_tbl %>% group_by(locusID) %>% summarise(chrom=chrom, start=min(pos_i), end=max(pos_f),
                                                            best_SNP_lod = max(best_SNP_lod),
                                                            n_genes_in_loci=n(),
                                                            genes_names=paste0(gene, collapse=", "),
                                                            Model= values) %>% distinct()


  # save to reactive value

  values$loci_tbl <- loci_tbl

  values$genes_tbl <- genes_tbl

  output$table <- renderDT(loci_tbl, selection = 'single')
      } # close the else
  })
  
  
  
  snps <- req(reactive({
    # get genotype data
    if(input$chr_data == "pruned"){
      if(input$chr == "Un"){
        fread("/home/udiland/gwas/snp/data/geno/all_samples/LD_pruned_100/long_chr8.hmp.txt", sep = '\t', head=T)
        #fread("/Users/udila/Downloads/long_chr8.hmp.txt", sep = '\t', head=T)
        
      } else(
        fread(paste0("/home/udiland/gwas/snp/data/geno/all_samples/LD_pruned_100/long_chr",input$chr,".hmp.txt", sep = ""), sep = '\t', head=T)
      )
    }
    if(input$chr_data == "filtered"){
      if(input$chr == "Un"){
        fread("/home/udiland/gwas/snp/data/geno/all_samples/long_chr8.hmp.txt", sep = '\t', head=T)
        
      } else(
        fread(paste0("/home/udiland/gwas/snp/data/geno/all_samples/long_chr",input$chr,".hmp.txt", sep = ""), sep = '\t', head=T)
        # load snp data
        #fread("/Users/udila/Downloads/long_chr5.hmp.txt", sep = '\t', header = T)
      )
    }
  }))
  
  # store the number of row that the user selected
  loci_row <- reactive(input$table_rows_selected)
  
  # plot LD based on loci selected from table
  observeEvent(input$get_plot, {
  # get the loci number selected
  loci_number <- req(values$loci_tbl[loci_row() ,]$locusID)
  
  # read genes and models for that loci
  gen_traits <- read.genes(values$genes_tbl, loci_number)
  print(gen_traits$genes)
  
  
  run <- req(values$loci_tbl[loci_row() ,]$Model)
  
  genes <- gen_traits$genes
  
  genes2plot <- genes_grrangelst[genes]
  
  #find all files, which are Results from GWAS
  traitcsv <- dir(values$InputDir, pattern = "*.Results.csv")
  
    pos <- grep(run,traitcsv) # get the index of the run in the runs vector taken from the results folder
    
    # read gwas results
    gwas <- fread(file.path(values$InputDir, traitcsv[pos])) 
    
    print(paste0("Ploting: ", run, collapse = " "))
    
    if (input$plot_title == ""){
      plotTitle<-paste0(run,"_",names(gen_traits$genes))
    } else (plotTitle <- input$plot_title)
    
    plot <- makeLDplot(gwas, longissima_genes, snps(), genes2plot, plotTitle, input$lodfilter, input$lodmax, input$winup, input$windown)
    
    if (identical(class(plot), "character")){
      output$x4 <- renderText(plot)}
      
   # plot output
   output$loci <- renderPlot(plot)

  # download the graph
   output$downloadPlot <- downloadHandler(
    filename = paste(plotTitle,".png", sep = ""),
    content = function(file) {
      ggbio::ggsave(file, plot, width = 10, height = 8, units = "cm")})
  
  }) # close obsereEvent
}

shinyApp(ui = ui, server = server)


