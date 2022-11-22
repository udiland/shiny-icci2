


get_loci_tbl <- function(gwas, longissima_genes, lod.thresh=5.3, locus.range=1e4, maf.thr=0.05){
  
  # combine all tables to one 
  z <- do.call('rbind', gwas)
  
  # sort by chrom and position
  z <- z[order(Chromosome, Position)]
  
  # filter snps by p-value (lod > threshold)
  z <- z[lod > lod.thresh]
  
  #additional filtering of maf
  z <- z[maf > maf.thr ]
  
  
  # Define loci, group SNPs to loci if they are in the same area as a predefined window (loci_range)
  # if a anp is in the same window of 'loci_range' (default 1e4) assign them to the same group
  
  
  # first define each snp to different loci (loci=number of snps)
  loci <- rep(1, nrow(z))
  # than give them same loci designation if they are in same loci_range
  for(i in 2:nrow(z)) {
    if(z[i, 'Chromosome'] == z[i-1, 'Chromosome'] &
       abs(z[i, 'Position'] - z[i-1, 'Position'] ) < locus.range)
    {
      loci[i] <- loci[i-1]
    } else {
      loci[i] <- loci[i-1] + 1
    }
  }
  
  # reset rownames
  rownames(z) <- NULL
  z <- data.frame(z)
  
  # add the group number to the table
  z$locusID <- loci
  
  # fix chron names
  z$Chromosome <- paste0("chr", z$Chromosome)
  
  #fix chrUN
  z$Chromosome <- ifelse(z$Chromosome == "chrUN", "chrUn", z$Chromosome)
  
  # for each locus get the gene names (give a list that each item is a locus and the content is a vector of gene names)
  genes <- lapply(unique(z$locusID), function(loci) {
    w <- z[z$locusID == loci, ]
    chr <- w$Chromosome[1]
    x  <- min(w$Position)
    y  <- max(w$Position)
    # collect the genes that are in a window of loci_range up or down the locus 
    j <- which(chr == longissima_genes$chrom & x - longissima_genes$pos_f < locus.range & longissima_genes$pos_i - y < locus.range)
    longissima_genes$gene[j]
  })
  
  
  # check if no gene has found
  if (identical(unique(genes)[[1]], character(0))){
    return(list(0, data.table()))
  }
  

  new_loci_lst <- list(1)# list with one elemnt equal to 1
  k <- 1
  for(i in 2:length(genes)) {
    if(any(genes[[i]] %in% unlist(genes[ new_loci_lst[[k]] ]))) {
      new_loci_lst[[k]] <- c( new_loci_lst[[k]] , i )
    } else {
      k <- k+1
      new_loci_lst[[k]] <- i
    }
  }

  # create new genes list for the new loci ID
  new_genes <- list()
  lociID <- rep(0, nrow(z))
  for(i in 1:length(new_loci_lst)) {
    new_genes[[i]] <- unique(unlist(genes[new_loci_lst[[i]]]))
    lociID[ z$locusID %in% new_loci_lst[[i]] ] <- i
  }
  
  # store new locus id in the snps data table
  genes_in_loci_lst <- new_genes
  
  # check that the there are no multiple genes
  stopifnot(length(unlist(genes_in_loci_lst)) == length(unique(unlist(genes_in_loci_lst))))
  z$locusID <- lociID
  rm(lociID, new_genes, k, i)

  # sort snps table by locusID and lod score and rename it
  z <- z[order(z$locusID, -z$lod), ]
  loci_nfo <- z
  # how many genes in each locus?
  ngenes <- sapply(genes_in_loci_lst, length)
  
  # take the genes in loci create a table of: chrom, start, end with gene names as rownames 
  genes_in_loci <- unlist(genes_in_loci_lst)
  genes_coord_tbl <- longissima_genes[match(genes_in_loci, longissima_genes$gene), c('chrom', 'pos_i', 'pos_f')]
  rownames(genes_coord_tbl) <- genes_in_loci
  
  
  
  # create an NA matrix with number of rows correspond to number of genes found in all 
  # loci and number of columns correspond to number of snp`s
  snp_gene_Distance_mtx <- matrix(NA, nrow(genes_coord_tbl), nrow(z))
  rownames(snp_gene_Distance_mtx) <- genes_in_loci
  
  # for each gene loop over all snp`s `
  for(i in 1:nrow(genes_coord_tbl)) {
    for(j in 1:nrow(z)) {
      if(z$Chromosome[j] != genes_coord_tbl$chrom[i])
        next
      if(z$Position[j] >= genes_coord_tbl$pos_f[i]) # if the snp is after the gene store the distance from the end of the gene
        snp_gene_Distance_mtx[i,j] <- genes_coord_tbl$pos_f[i] - z$Position[j]
      else if(z$Position[j] <= genes_coord_tbl$pos_i[i]) # if the is snp before the gene store the distance from the start
        snp_gene_Distance_mtx[i,j] <- genes_coord_tbl$pos_i[i] - z$Position[j]
      else # else, the snps is in the gene so the distance is 0
        snp_gene_Distance_mtx[i,j] <- 0
    }
  }
  
  # go over each locus and store the best snp (highest lod) and nearest snp
  # associated to genes in that locus and store to lists, 'best' and 'nearest'
  nearest <- list()
  best    <- list()
  for(i in 1:length(genes_in_loci_lst)) {
    j <- which(z$locusID == i) # rows numbers in snp table for snp that have the same locusID
    best[[i]] <- rep(j[1], length(genes_in_loci_lst[[i]])) 
    d <- snp_gene_Distance_mtx[genes_in_loci_lst[[i]], j, drop=F]
    nearest[[i]] <- j[apply(d, 1, function(w) which.min(abs(w)))]
  }
  
  # change to vector
  nearest <- unlist(nearest)
  best    <- unlist(best)
  
  # create a table of genes with data of best and nearest snp
  gene_nfo <- data.frame(
    locusID = rep(1:length(genes_in_loci_lst), times=ngenes),
    no_SNP  = rep(table( z$locusID ), times=ngenes),
    best_SNP      = z$SNP[best],
    best_SNP_lod  = z$lod[best],
    best_SNP_dist = mapply(function(i,j) snp_gene_Distance_mtx[i,j], 1:length(best), best),
    nearest_SNP   = z$SNP[nearest],
    nearest_SNP_lod = z$lod[nearest],
    nearest_SNP_dist = mapply(function(i,j) snp_gene_Distance_mtx[i,j], 1:length(nearest), nearest),
    best_eq_nearest = best == nearest
  )
  
  # add data on each gene (from the gff)
  tr = longissima_genes[match(genes_in_loci, longissima_genes$gene), ]
  gene_nfo <- cbind(gene_nfo, tr)
  
  
  # this section add the 'Run' names to the gene_info data. (the name of the GAPIT run)
  runs_by_loci = split(loci_nfo$Run, loci_nfo$locusID) # get a list for each locus the names of the 'RUN`s'
  runs_by_loci = sapply(runs_by_loci, unique, simplify=FALSE) # remove duplicates
  
  m = max(sapply(runs_by_loci, length)) # mximum nuber of runs
  
  # create a matrix with each row as loci and column the different runs (add '.' to fill empty cells)
  if(m == 1){
    runs_by_loci <- cbind(Peak.ID=unlist(runs_by_loci))
  } else {
    runs_by_loci = t(sapply(runs_by_loci, function(w) if(length(w) == m) w else c(w, rep(".", m-length(w)))))
    colnames(runs_by_loci) <- sprintf("Run.ID_%s", 1:m)
  }
  
  
  # combine with genes info table
  runs_by_loci = runs_by_loci[gene_nfo$locusID,, drop=FALSE]
  gene_nfo <- cbind(gene_nfo, runs_by_loci)
  
  return(list(loci_nfo, gene_nfo))

}
  
#a <- get_loci_tbl(gwas, longissima_genes, locus.range = 1e6)
  