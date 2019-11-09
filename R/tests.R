glf2array <- function(name, nsite, nsamp, isfloat=FALSE)
{
  array(readBin(name, "double", size=if(isfloat) 4 else 8, n=3*nsamp*nsite), dim=c(3,nsamp,nsite))
}

count2array <- function(name)
{
  in1 <- read.table(name, header=T)
  nsamp <- ncol(in1)/4
  nsite <- nrow(in1)
  out <- array(c(t(as.matrix(in1))), c(4,nsamp,nsite))
  attr(out, "nsamp") <- nsamp
  attr(out, "nsite") <- nsite
  out
}

pos2array <- function(name)
{
  pos <- read.table(name)
  contigs <- levels(pos[,1])
  pos[,1] <- as.numeric(pos[,1])-1
  pos[,2] <- pos[,2]-1
  pos[,3] <- match(pos[,3], c("A","C","G","T"))-1
  pos[,4] <- match(pos[,4], c("A","C","G","T"))-1
  pos     <- as.matrix(pos)
  attr(pos, "contigs") <- contigs
  pos
}

readglf <- function(prefix, isfloat=FALSE)
{
  counts <- count2array(paste0(prefix,".counts.gz"))
  glf <- glf2array(paste0(prefix,".glf"), attr(counts, "nsite"), attr(counts, "nsamp"), isfloat)
  pos <- pos2array(paste0(prefix,".glf.pos.gz"))
  list(counts=counts,glf=glf,pos=pos)
}

outline <- function()
{
  library(haplodiplo)
  library(Matrix)

  contig <- "ScQjuvi_472"
  outdir <- "output"

  dir.create(outdir)

  data <- readglf(contig)
  samps <- read.csv("../pepo_samples.csv")
  samps <- rbind(samps, data.frame(BAM=NA, Sample="XJ002", Species=NA, Sex="M", Ploidy=1, Lat=NA, Lon=NA, State="MEX", County=NA))
  samps <- rbind(samps, data.frame(BAM=NA, Sample="XJ005", Species=NA, Sex="M", Ploidy=1, Lat=NA, Lon=NA, State="PA", County=NA))

  Hd <- Haplodiplo$new(data$glf, data$counts, data$pos, samps$Ploidy)

  Hd$hetbias(0:95, c(1e-6, 1e-3), rep(0.0010,96))

  strat <- as.numeric(samps$State)-1
  mind  <- rep(2, length(unique(strat)))
  mind[which(levels(samps$State)=="MEX")] <- 1
  Hd$minchrom(0:95, strat, mind)

  saf <- list()
  for (state in levels(samps$State))
    saf[[state]] <- Hd$saf(0:(Hd$sites()-1), which(samps$State==state)-1)
  pos <- data.frame(Hd$position())
  pos[,1] <- contig
  pos[,3] <- c("A","C","G","T")[pos[,3]+1]
  pos[,4] <- c("A","C","G","T")[pos[,4]+1]

  save(saf, pos, file=paste0(outdir, "/", contig, ".saf.RData"))

  Hd$minmaf(0:95, 0, 1)

  glf <- Hd$likelihoods()
  counts <- Hd$counts()
  pos <- data.frame(Hd$position())

  pos[,1] <- contig
  pos[,3] <- c("A","C","G","T")[pos[,3]+1]
  pos[,4] <- c("A","C","G","T")[pos[,4]+1]

  writeBin(c(glf), paste0(outdir, "/", contig, ".variable.fglf"), size=4)
  writeBin(as.integer(c(counts)), paste0(outdir, "/", contig, ".variable.scnt"), size=2)
  write.table(pos, paste0(outdir, "/", contig, ".variable.pos"), row.names=FALSE, col.names=FALSE, quote=F, sep="\t")

  #test that write was OK
  array(readBin(paste0(outdir, "/", contig, ".variable.fglf"), "double", size=4, n=96*3*nrow(pos)), c(3,96,nrow(pos)))[,,2]
  array(readBin(paste0(outdir, "/", contig, ".variable.scnt"), "integer", size=2, n=96*4*nrow(pos)), c(4,96,nrow(pos)))[,,2]
  read.table(paste0(outdir, "/", contig, ".variable.pos"))[2,]

}

tests <- function()
{
  library(haplodiplo)
  library(Matrix)

  load("test.RData")

  pos[,1] <- as.numeric(pos[,1])
  pos[,3] <- match(pos[,3], c("A","C","G","T"))-1
  pos[,4] <- match(pos[,4], c("A","C","G","T"))-1
  pos <- as.matrix(pos[,1:4])

  # introduce a couple missing sites
  glf[,1,3] <- 0
  glf[,4,3] <- 0
  glf[,3,8] <- 0

  Hd <- Haplodiplo$new(glf[,,1:10], counts[,,1:10], pos[1:10,], ploidy)

  # are GLF correctly read in?
  all(round(Hd$likelihoods()[,,2],5) == round(apply(glf[,,2], 2, function(x) exp(x)/sum(exp(x))),5))

  # are counts correctly read in?
  all(Hd$majorminor()[,,2] == counts[4:3,,2])

  # are missing data correctly detected?
  all(which(Hd$missing()==1, arr.ind=T) == rbind(c(1,3),c(4,3),c(3,8)))

  # are missing counts correctly eliminated?
  all(Hd$counts()[,3,8] == 0)

  # are missing tallies calculated correctly?
  all(Hd$missing_sites()[c(1,3,4)] == c(1,1,1))
  all(Hd$missing_samples()[c(3,8)] == c(2,1))

  # filter by strata
  Hd$minsample(0:12, c(0,1,1,0,0,1,1,1,1,1,1,1,1), c(2,1))
  Hd$sites()
  str(Hd$likelihoods())
  str(Hd$counts())

  # filter by MAF
  Hd$minmaf(0:12, 0.10, 1.0)
  Hd$sites()
  str(Hd$likelihoods())
  str(Hd$counts())

  # hetbias filter for 
  Hd <- Haplodiplo$new(glf, counts, pos, ploidy)
  ooo <- glf[,Hd$ploidy() == 1,]
  goo <- Hd$counts()[,Hd$ploidy() == 1,]
  het <- apply(ooo,c(2,3),function(x) if(abs(x[1]-x[2])<0.0001 & abs(x[2]-x[3])<0.0001) 9 else which.max(x)-1) == 1
  ooh <- goo[1,,][het] + goo[2,,][het]
  ohh <- ifelse(goo[1,,][het] < goo[2,,][het], goo[1,,][het], goo[2,,][het])
  sum(1-pbinom(ohh-1,ooh,0.0004) < 1e-5) == Hd$hetbias(0:12,c(1e-5,0),rep(0.0004,13))

  ooo <- glf[,Hd$ploidy() == 2,]
  goo <- Hd$counts()[,Hd$ploidy() == 2,]
  het <- apply(ooo,c(2),function(x) if(abs(x[1]-x[2])<0.0001 & abs(x[2]-x[3])<0.0001) 9 else which.max(x)-1) == 1
  ooh <- goo[1,][het] + goo[2,][het]
  ohh <- ifelse(goo[1,][het] > goo[2,][het], goo[1,][het], goo[2,][het])
  sum(1-pbinom(ohh-1, ooh, 0.5) < 1e-3) == Hd$hetbias(0:12,c(0,1e-3),rep(0.0004,13))

  # does polarization work?
  Hd <- Haplodiplo$new(glf[,,1:10], counts[,,1:10], pos[1:10,], ploidy)

  Hd$polarize(c(pos[1:5,3],pos[6:10,4]))

  # general filtering proceedure:
  # Hetbias ==> missing data
  # (do SAF stuff) ==> frequency/SNP
  Hd <- Haplodiplo$new(glf, counts, pos, ploidy)
  Hd$hetbias(0:12,c(1e-1,1e-3),rep(0.0004,13))
  Hd$minsample(0:12, rep(0,13), c(13))

  # does other stuff work
  Hd$frequencies(0:(Hd$sites()-1), 0:12)
  Hd$genotypes(0:(Hd$sites()-1), 0:12)
  Hd$paralogs(0:(Hd$sites()-1), 0:12)
  #Hd$maxhet(0:12, 0.75) #will remove all hets in this example!

  # SFS estimation
  saf0 <- Hd$saf(0:(Hd$sites()-1), 0:12)

  n_blocks <- 20
  block <- rep(0:19, each=ceiling(Hd$sites()/n_blocks))[1:Hd$sites()]
  samps <- sfs2d(saf0, saf0, block, 5, TRUE)

  sfs1d(saf0, block, 0, FALSE)
  sfs1d(saf0, block, 0, TRUE)

  sfs2d(saf0, saf0, block, 0, FALSE)
  sfs2d(saf0, saf0, block, 0, TRUE)

  sfs3d(saf0, saf0, saf0, block, 0, FALSE)
  sfs3d(saf0, saf0, saf0, block, 0, TRUE)

  # Fst
  saf0 <- Hd$saf(0:(Hd$sites()-1), 0:6)
  saf1 <- Hd$saf(0:(Hd$sites()-1), 6:12)
  myunfold <- sfs2d(saf0, saf1, block, 0, FALSE)
  myfold <- sfs2d(saf0, saf1, block, 0, TRUE)
  head(FST(saf0, saf1, myunfold))
  head(FST(saf0, saf1, myfold))


  # SNP calling
  glf[,1,1] <- 0
  Hd <- Haplodiplo$new(glf, counts, pos, ploidy)

  n_blocks <- 20
  block <- rep(0:19, each=ceiling(Hd$sites()/n_blocks))[1:Hd$sites()]
  saf0 <- Hd$saf(0:(Hd$sites()-1), 0:6)
  sfs_fold <- sfs1d(saf0, block, 0, TRUE)
  sfs_unfold <- sfs1d(saf0, block, 0, FALSE)

  stats<-Hd$postsnp(1:Hd$sites()-1, 0:6, sfs_unfold)
  stats<-Hd$postsnp(1:Hd$sites()-1, 0:6, sfs_fold)
  stats$mean[1]

  Hd$callsnps(0:6, rep(0,7), list(sfs_fold), 0, 1)
  Hd$retain(0:2)
  Hd$genotypes(TRUE)
  Hd$genotypes(FALSE)
}
