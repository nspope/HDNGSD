snp2vcf <- function()
{

  library(haplodiplo)

  # call SNPs here
  pos    <- read.table("clades.snp.pos.gz")
  samps  <- read.csv("../pepo_samples.csv")
  wonky  <- samps$Sample %in% c("PP0809", "PP1197", "PP0807") #contaminated, duplicate, duplicate
  samps  <- samps[!wonky,]

  gno    <- array(readBin("clades.snp.fglf", size=4, what="numeric", n=3*nrow(pos)*nrow(samps)), c(3,nrow(samps),nrow(pos)))
  cnt    <- array(readBin("clades.snp.scnt", size=2, what="integer", n=4*nrow(pos)*nrow(samps)), c(4,nrow(samps),nrow(pos)))

  Hd <- Haplodiplo$new(...)
  #need to filter here just to dbl check

  #make VCF output here
  chrlen <- read.table("Epru.chrlen")

  #thinit
  OK <- samps$State %in% c("MA", "OH", "PA", "NY", "ON", "QC", "MS", "DE") & samps$Ploidy == 1
  gno    <- gno[,OK,]
  samps  <- samps[OK,]

  header <- 
'#fileformat=VCFv4.2
##source="HDNGSD"
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'

  chrlen <- apply(chrlen, 1, function(x) paste0('##contig=<ID="', x[1], '",length=', x[2]))

  colnm <-
  '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT'

  colnm <- paste(c(colnm, as.character(samps$Sample)), collapse="\t")

  gno <- apply(gno, 1, function(x) paste(paste0(x, "/", x), collapse="\t"))

  body <- cbind(position[,1:2], ".", position[,3:4], ".", ".", ".", gno)
  body <- apply(body, 1, function(x) paste(x, collapse="\t"))

  sink("clades.snp.vcf")
  cat(header, "\n")
  cat(chrlen, "\n")
  cat(colnm, "\n")
  cat(gno, "\n")
  cat(body, "\n")
  sink()
  
}
