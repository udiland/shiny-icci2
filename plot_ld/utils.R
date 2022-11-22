"
This function is used to create a list of gene names and name of gapit run for specific locus from 
a table that has all genes and euns fro all loci
"

read.genes <-function(gene_table, locus){
  ## function that takes the gene info table created with gwas output (in 'build_loci_gene_tbls.R')
  # and locusID number
  # and generate a list with 2 elements, $genes - name of genes in loci and $Runs - name of Run (GAPIT algo) 
  gen.name <-c()
  genID <-c()
  Runs <-c()
  # make a vector of genes in locus
  genID <- append(genID, as.vector(gene_table[which(gene_table[,1]==locus),15])) 
  gen.name <- append(gen.name, as.vector(gene_table[which(gene_table[,1]==locus),15])) 
  # gen.name <-append(gen.name,(LG[which(LG[,1]==locus),16])) 
  names(gen.name)<- genID
  
  # collect the different type of runs (usually different algorithms of GAPIT)
  for (i in 16:ncol(gene_table)){ 
    Runs <-append(Runs,as.vector(gene_table[which(gene_table[,1]==locus),i]))
  }
  Runs <-unique(Runs)
  Runs <-Runs[Runs!="."]
  gen_traits <-list("genes"=gen.name,"Runs"=Runs)
  
  return(gen_traits)
}

################################################################################

"
This script will read all gapit gwas results in a specified folder, 
and output an R list with each item as a result table of a different run found in the input folder
can be runs of different algorithms of te same data or different gwas runs
(the results are named after the the algorithm and name of treatment and are filtered by p.value)
"
getResults <- function(path, p_value_thresh=0.001) {
  # create a vector of gwas results file names
  files <- dir(path, pattern="GWAS.Results.csv")
  results  <- list()
  for(f in files) {
    message(f)
    # make a regex object and extract from file name the name of algorithm used and treatment,
    # e.g blink.ME1_Control or MLM.ME1_Control etc
    z <- regexpr("GAPIT\\.\\.?(.*)\\.GWAS.Results", f, perl=TRUE)
    stopifnot(z != -1)
    run <- substring(f, attr(z, 'capture.start')[1],
                     attr(z, 'capture.start')[1] + attr(z, 'capture.length')[1] - 1)
    # read gwas results
    res <- fread(file.path(path, f))
    # add lod score column 
    res$lod <- -log10(res$P.value)
    res$SNP <- rownames(res)
    res$Run <- run
    results[[run]] <- res[res$P.value < p_value_thresh, ]
  }
  return(results)
}

########################################################################################

"
This function read genotypes in hapmap format (one letter) and output the full genotype two letters
for example: 'A' -> 'A/A' or 'W' -> 'A/T'
"
allel2geno <- function(x){
  if (x == "A"){
    return("A/A")
    
  }else if (x == "G"){
    return("G/G")
    
  }else if (x == "C"){
    return("C/C")
    
  }else if (x == "T"){
    return("T/T")
    
  }else if (x == "R"){
    return("A/G")
    
  }else if (x == "Y"){
    return("C/T")
    
  }else if (x == "S"){
    return("G/C")
    
  }else if (x == "W"){
    return("A/T")
    
  }else if (x == "K"){
    return("G/T")
    
  }else if (x == "M"){
    return("A/C")
    
  }else {
    return(NA)
    
  }
  
}

##################################################################################

get_ld <- function(snps_gwas, snps_hmp){
  # oredr gwas by lod score (high to low)
  snps_gwas <- snps_gwas[order(-snps_gwas$lod) ,]
  # take snps names
  snps_in_loci <- snps_gwas$SNP
  
  # select only snps in loci
  snps <- snps_hmp[snps_hmp$`rs#` %in% snps_in_loci ,]
  
  # make an empty matrix that could store the new genotype data
  res <- data.frame(matrix(nrow = nrow(snps), ncol = (ncol(snps)- 11)))
  
  # go over each row in the snps data, transform each cell from hapmap format
  for (i in 1:nrow(snps)){
    res[i,] <- sapply(snps[i, 12:ncol(snps)], allel2geno)}
  
  # transpose the table so each sample (accession) is a column
  res_t <- as.data.frame(t(res))
  # give sample name
  colnames(res_t) <- snps$`rs#`     
  
  # order the table by the value of lod score of each accession (column)
  # get the index of samples in res_t by the location in the snps_gwas table that is sorted by lod score 
  col_numbers <-match( snps_gwas$SNP ,colnames(res_t))
  res_t <- res_t[, c(col_numbers)]
  
  data <- makeGenotypes(res_t)
  
  print("Calculating LD")
  
  ld_object <- LD(data)
  
  print("LD calculated")
  
#  r_square <- ld_object$`R^2`
   r_square <- ld_object$`D'`
  

    # return the value of LD compare with best snp (first row of matrix)
  return(r_square[1,])
}
  
########################################################################################

# make plot
makeLDplot <- function(gwas, longissima_genes, snp, genes, plotTitle,
                       lodFilter, lodMax, winUp, winDown)
{    
  
  gwas$lod <- -log10(gwas$P.value)
  gwas <- gwas[gwas$lod > lodFilter, ]
  gwas <- gwas[order(gwas$lod),] # put less significant SNPs first
  
  # store vector of lod scores
  lod <- gwas$lod
  
  # get genes coordination
  genes_info <- longissima_genes[ match(names(genes),longissima_genes$gene), c("chrom","pos_i","pos_f", "strand", "gene")]
  
  # #change gwas chrom name so it match with all gene data chrom names
  # if (gwas$Chromosome[0] == "UN") {
  #    gwas$Chromosome <- "chrUn"} else {  
  #    gwas$Chromosome <- paste("chr", gwas$Chromosome, sep = "")}
  
  # take all snps rows in the range of [minimum gene start - window down, maximum gene end + window up]
  rows_snps_in_loci <- which(gwas[,Chromosome] == genes_info[1,1] & gwas[,Position] >= min(genes_info$pos_i) - winDown &
                               gwas[,Position] <= max(genes_info$pos_f) + winUp)
  
  # if there are less than 2 snps there is nothing to plot
  if (length(rows_snps_in_loci) <= 1){
    return("Less than 2 snps in loci ")
  }
  
  # make a table of snp positions (in Mb) and lod score
  #data_for_plot <- as.data.frame(cbind(gwas$Pos[rows_snps_in_loci], lod[rows_snps_in_loci], gwas$FDR_Adjusted_P.values[rows_snps_in_loci]))    
  
  #colnames(data_for_plot) <- c("Position", "lod", "FDR")
  data_for_plot <- gwas[rows_snps_in_loci, c("Position", "lod", "FDR_Adjusted_P-values")]
  
  # maximum lod score
  best_snp_lod_index <- which.max(lod[rows_snps_in_loci])
  #stopifnot(best_snp_lod_index == nrow(data_for_plot))
  
  # take relevant snps from gwas results
  snps_gwas <- gwas[rows_snps_in_loci, ]
  
  # get the LD R^2 calculations for the best snp
  ld <- get_ld(snps_gwas, snp)
  # give the ref snp value of 1
  ld[1] <- 1
  # if some snp hase no value st it to 0
  ld[is.na(ld)] <- 0
  

  
  # take the data of best snps
  best_snp <- data_for_plot[best_snp_lod_index ,]
  
  # make a table with gene sart gene end for plot shading (geom_rect)
  rects <- data.frame(xstart = genes_info[,2], xend = genes_info[,3])
  
  # calculate ld value of fdr
  lod_sig_thresh <- min(data_for_plot[data_for_plot$`FDR_Adjusted_P-values` < 0.05 , ]$lod)
  
  # use ggbio package to plot genes
  plt <- autoplot(genes, fill="steelBlue") + 
    geom_arrow(arrow.rate=0.1, type='closed', angle = 30, length = unit(0.25, "cm")) +
    theme_pack_panels()
  
  p <- ggplot() +
       geom_rect(data = rects, aes(ymin = 0, ymax = as.numeric(data_for_plot[best_snp_lod_index,2]+3), xmin = xstart, xmax=xend), alpha=0.3) + 
       scale_color_gradient2(low="pink",high="red", mid="orange", midpoint=0.5) +
       geom_point(data = data_for_plot[-best_snp_lod_index,], aes(x=Position, y=lod, col = ld[2:length(ld)]), size=3) + 
       geom_point(data = best_snp, x=best_snp$Position, y=best_snp$lod, color="red", size=3, shape=8) +
       theme_classic() +
       labs(color= expression("D`")) +
       theme(axis.text.x=element_blank()) + 
       ylab("LOD Score")
  
  if (lod_sig_thresh != Inf){p <- p + geom_hline(yintercept = lod_sig_thresh, linetype='dotted', col = 'black', size=1)}
  
  
  xlim <- c(min(data_for_plot$Position) - 100, max(data_for_plot$Position) + 100)
  
  tracks_all <- tracks(list(p, plt), xlim = xlim, heights = c(5,3), main = plotTitle)
  
  return(tracks_all)
  
}



