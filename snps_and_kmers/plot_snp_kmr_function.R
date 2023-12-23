

plot_snp_kmers <- function(kmer_tbl, gwas_tbl, pv_limit=0.001){
  
  # filter GWAS results
  gwasResults <- gwas_tbl[gwas_tbl['P.value'] < pv_limit ,]	
  
  # take only kmers tha were aligned
  kmers <- kmer_tbl[kmer_tbl$state != "not_aligned" ,]
  
  # take only kmers that intersect with scaffold
  kmers <- kmers %>% group_by(name_kmer) %>% filter(startsWith(as.character(name_contig), "NODE"))
  
  # mark partial match of alignment
  kmers$match <- ifelse(kmers$match == "full", "full", "partial")
  
  # mark if the kmer is match uniqly
  kmers$state <- ifelse(kmers$state == "duplicate", "duplicate", "uniqe")
  
  # take only uniqly aligned kmers
  kmers <- kmers[kmers$state == "uniqe" ,]
  
  # change the name of chromosomes
  kmers$chr_kmer[kmers$chr_kmer == 'chr1S'] <- '1S'
  
  kmers$chr_kmer[kmers$chr_kmer == 'chr2S'] <- '2S'
  
  kmers$chr_kmer[kmers$chr_kmer == 'chr3S'] <- '3S'
  
  kmers$chr_kmer[kmers$chr_kmer == 'chr4S'] <- '4S'
  
  kmers$chr_kmer[kmers$chr_kmer == 'chr5S'] <- '5S'
  
  kmers$chr_kmer[kmers$chr_kmer == 'chr6S'] <- '6S'
  
  kmers$chr_kmer[kmers$chr_kmer == 'chr7S'] <- '7S'
  
  kmers$chr_kmer[kmers$chr_kmer == 'chrUn'] <- 'UN'
  
  # take relevant columns
  kmers <- kmers[, c(1,3,4,17)]
  
  kmers <- as.data.frame(kmers)
  
  head(gwasResults)
  
  colnames(kmers)[2] <- "Chromosome"
  
  colnames(kmers)[3] <- "Position"
  
  res <- merge(gwasResults, kmers, by = c("Chromosome", "Position"), all = T)
  
  
  don <- res %>%
    
    # Compute chromosome size
    group_by(Chromosome) %>%
    summarise(chr_len=max(Position)) %>%
    
    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
    select(-chr_len) %>%
    
    # Add this info to the initial dataset
    left_join(res, ., by=c("Chromosome"="Chromosome")) %>%
    
    # Add a cumulative position of each SNP
    arrange(Chromosome, Position) %>%
    mutate(BPcum=Position+tot)
  
  
  axisdf = don %>% group_by(Chromosome) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
  
  # add a column name 'size' that hold the number of same value 'logPV' in a window of 10 positions
  don <- don %>% group_by(Chromosome, grp = rep(row_number(), length.out = n(), each = 10), logPV) %>% mutate(size=n())
  
  # create the plot
  p <- ggplot(data = don, aes(x=BPcum, y=-log10(P.value)), size=size) +
    
    # plot snps
    geom_point(aes(col=as.factor(don$Chromosome))) +
    
    scale_color_manual(values = rep(c("darkblue", "blue"), 4)) +
    
    new_scale_colour() +
    
    # plot the kmers
    geom_point(aes(x=BPcum, y=logPV, size=size, col=as.factor(don$Chromosome))) +
    
    scale_color_manual(values = rep(c("darkgreen", "green"), 4)) +
    
    # custom X axis:
    scale_x_continuous(label = axisdf$Chromosome, breaks= axisdf$center) +
    
    # add a line that separate FDR significant snps
    geom_hline(yintercept = -log10(max(na.exclude(don$P.value[don$FDR_Adjusted_P.values < 0.05]))), linetype="dashed",
               color = "red", size=1, alpha=0.5) +
    
    # remove space between plot area and x axis and set the y axis limits
    scale_y_continuous(expand = c(0, 0), limits = c(min(-log10(don$P.value)), max(na.omit(don$logPV) + 0.5))) +
    
    
    # Custom the theme:
    theme_bw() +
    theme(
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.title.x=element_blank()
    ) 
  
  return(p)
}

# gwas <- read.csv("GAPIT.MLM.FE1_DTA.GWAS.Results.csv")
# kmers_tbl <- read.csv("data_for_plot_pass5.csv")
# 
# plot_snp_kmers(kmers_tbl, gwas)
