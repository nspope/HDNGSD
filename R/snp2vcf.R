snp2vcf <- function()
{

  library(haplodiplo)

  # call SNPs here
  pos    <- read.table("clades.snp.pos.gz")
  samps  <- read.csv("../pepo_samples.csv")
  glf    <- array(readBin("clades.snp.fglf", size=4, what="numeric", n=3*nrow(pos)*nrow(samps)), c(3,nrow(samps),nrow(pos)))
  cnt    <- array(readBin("clades.snp.scnt", size=2, what="integer", n=4*nrow(pos)*nrow(samps)), c(4,nrow(samps),nrow(pos)))

  # remove the wonky bees and diploids
  wonky  <- samps$Sample %in% c("PP0809", "PP1197", "PP0807") | samps$Ploidy==2 #contaminated, duplicate, duplicate; diploids
  samps  <- samps[!wonky,]
  glf    <- glf[,!wonky,]
  cnt    <- cnt[,!wonky,]
  Hd <- Haplodiplo$new(glf, cnt, pos, samps$Ploidy)

  # priors and stuff
  samps$State <- factor(as.character(samps$State))
  strat <- as.numeric(samps$State)-1
  mind  <- rep(1, length(unique(strat)))
  load("clades.sfs1d.fold.RData")
  Sfs1d <- Sfs1d[levels(samps$State)] 

  # filter individual genotypes
  Hd$hetbias(1:Hd$samples()-1, c(1e-5, 1e-2), rep(0.0010,Hd$samples()))

  # do I want to be more stringent about missing data?

  # calculate genotype probabilities ... should also remove fixed sites?
  Hd$callsnps(1:Hd$samples()-1, strat, Sfs1d, -1, 0.05)

  #make VCF output here
  chrlen <- read.table("Epru.chrlen")

  #thinit
  OK <- samps$State %in% c("MA", "OH", "PA", "NY", "ON", "QC", "MS", "DE") & samps$Ploidy == 1
  gno    <- gno[,OK,]
  samps  <- samps[OK,]

  header <- 
'#fileformat=VCFv4.2
##source="HDNGSDv0.1"
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'

  chrlen <- apply(chrlen, 1, function(x) paste0('##contig=<ID="', x[1], '",length=', x[2]))

  colnm <-
  '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT'

  colnm <- paste(c(colnm, as.character(samps$Sample)), collapse="\t")

  gno <- apply(gno, 1, function(x) paste(paste0(x, "/", x), collapse="\t"))

  body <- cbind(position[,1:2], ".", position[,3:4], ".", ".", ".", gno)
  body <- apply(body, 1, function(x) paste(x, collapse="\t"))

  sink("clades.east.snp.vcf")
  cat(header, "\n")
  cat(chrlen, "\n")
  cat(colnm, "\n")
  cat(gno, "\n")
  cat(body, "\n")
  sink()
  
}
