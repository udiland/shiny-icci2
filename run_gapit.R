# Get GAPIT
source("/home/username/GAPIT.library.R")
source("/home/username/gapit_functions.txt")

# upload phenotype file
myY <- read.table("/home/username/../phenotype.tsv", sep='\t', head=T)

myY[,2] <- as.numeric(myY[,2])


# Read genomic data - this will load the LD fiiltered file with snps in all chromosomes
myG <- read.table("/home/udiland/gwas/snp/data/geno/all_samples/LD_pruned_100/long_all_pruned.hmp.txt", sep='\t', comment.char = "", head=F)


# Run GAPIT3
myGAPIT <- GAPIT(
  Y=myY,
  PCA.total=2,
  G=myG,
# the models to run
  model=c("GLM","MLM", "MLMM", "FarmCPU", "Blink")
)
