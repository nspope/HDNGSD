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
  n_blocks <- 20
  block <- rep(0:19, each=ceiling(Hd$sites()/n_blocks))[1:Hd$sites()]
  saf0 <- Hd$saf(0:(Hd$sites()-1), 0:6)
  saf1 <- Hd$saf(0:(Hd$sites()-1), 6:12)
  myunfold <- sfs2d(saf0, saf1, block, 0, FALSE)
  myfold <- sfs2d(saf0, saf1, block, 0, TRUE)
  whoa = FST(saf0, saf1, myunfold[,,1])
  head(FST(saf0, saf1, myfold[,,1]))
  sum(whoa[,1])/sum(whoa[,2])


  # SNP calling
  glf[,1,1] <- 0
  Hd <- Haplodiplo$new(glf, counts, pos, ploidy)

  n_blocks <- 20
  block <- rep(0:19, each=ceiling(Hd$sites()/n_blocks))[1:Hd$sites()]
  saf0 <- Hd$saf(1:Hd$sites()-1, 0:6)
  sfs_fold <- sfs1d(saf0, block, 0, TRUE)
  sfs_unfold <- sfs1d(saf0, block, 0, FALSE)

  stats<-Hd$postsnp(1:Hd$sites()-1, 0:6, sfs_unfold)
  stats<-Hd$postsnp(1:Hd$sites()-1, 0:6, sfs_fold)
  stats$mean[1]

  Hd$callsnps(0:6, rep(0,7), list(sfs_fold), 0, 1)
  Hd$retain(0:2)
  Hd$genotypes(TRUE)
  Hd$genotypes(FALSE)

  saf0 <- Hd$saf(1:Hd$sites()-1, 1:Hd$samples()-1)
  sfs_fold <- sfs1d(saf0, block, 0, TRUE)
  Hd$callsnps(1:Hd$samples()-1, rep(0,Hd$samples()), list(sfs_fold), 0, 0.01)

}

.test_shf <- function()
{
  library(haplodiplo)
  library(Matrix)

  # haplotype likelihood (haploid)
  glf <- array(0, c(3,4,3))
  glf[1,1,1] <- 1; glf[3,2,1] <- 1; glf[1,3,1] <- 1; glf[3,4,1] <- 1
  glf[3,1,2] <- 1; glf[3,2,2] <- 1; glf[3,3,2] <- 1; glf[3,4,2] <- 1
  glf[3,1,3] <- 1; glf[3,2,3] <- 1; glf[1,3,3] <- 1; glf[3,4,3] <- 1
  counts <- array(0, c(4,4,3))
  counts[1,1,1] <- 1; counts[2,2,1] <- 1; counts[1,3,1] <- 1; counts[2,4,1] <- 1
  counts[2,1,2] <- 1; counts[2,2,2] <- 1; counts[2,3,2] <- 1; counts[2,4,2] <- 1
  counts[2,1,3] <- 1; counts[2,2,3] <- 1; counts[1,3,3] <- 1; counts[2,4,3] <- 1
  pos <- cbind(0, c(1,2,3), rep(0,3), rep(1,3))
  Hd <- Haplodiplo$new(log(glf), counts, pos, rep(1,4))
  Hd$shf(rbind(c(0,1),c(0,2),c(1,2)), 0:3)

  # diploid likelihood
  glf <- array(0, c(3,4,3))
  glf[1,1,1] <- 1; glf[3,2,1] <- 1; glf[1,3,1] <- 1; glf[2,4,1] <- 1 # aa AA aa aA
  glf[3,1,2] <- 1; glf[3,2,2] <- 1; glf[3,3,2] <- 1; glf[2,4,2] <- 1 # BB BB BB bB
  glf[3,1,3] <- 1; glf[2,2,3] <- 1; glf[1,3,3] <- 1; glf[3,4,3] <- 1 # CC cC cc CC
  counts <- array(0, c(4,4,3))
  counts[1,1,1] <- 2; counts[2,2,1] <- 2;   counts[1,3,1] <- 2; counts[1:2,4,1] <- 1
  counts[2,1,2] <- 2; counts[2,2,2] <- 2;   counts[2,3,2] <- 2; counts[1:2,4,2] <- 1
  counts[2,1,3] <- 2; counts[1:2,2,3] <- 1; counts[1,3,3] <- 2; counts[2,4,3] <- 2
  pos <- cbind(0, c(1,2,3), rep(0,3), rep(1,3))
  Hd <- Haplodiplo$new(log(glf), counts, pos, rep(2,4))
  Hd$shf(rbind(c(0,1),c(0,2),c(1,2)), 0:3)

  # haplodiploid likelihood
  glf <- array(0, c(3,4,3))
  glf[1,1,1] <- 1; glf[3,2,1] <- 1; glf[1,3,1] <- 1; glf[2,4,1] <- 1 # a AA aa aA
  glf[3,1,2] <- 1; glf[3,2,2] <- 1; glf[3,3,2] <- 1; glf[2,4,2] <- 1 # B BB BB bB
  glf[3,1,3] <- 1; glf[2,2,3] <- 1; glf[1,3,3] <- 1; glf[3,4,3] <- 1 # C cC cc CC
  counts <- array(0, c(4,4,3))
  counts[1,1,1] <- 1; counts[2,2,1] <- 2;   counts[1,3,1] <- 2; counts[1:2,4,1] <- 1
  counts[2,1,2] <- 1; counts[2,2,2] <- 2;   counts[2,3,2] <- 2; counts[1:2,4,2] <- 1
  counts[2,1,3] <- 1; counts[1:2,2,3] <- 1; counts[1,3,3] <- 2; counts[2,4,3] <- 2
  pos <- cbind(0, c(1,2,3), rep(0,3), rep(1,3))
  Hd <- Haplodiplo$new(log(glf), counts, pos, c(1,rep(2,3)))
  Hd$shf(rbind(c(0,1),c(0,2),c(1,2)), 0:3)

  # haplotype likelihood (haploid, missing data)
  glf <- array(0, c(3,4,3))
  glf[1:3,1,1] <- 1/3; glf[3,2,1] <- 1; glf[1,3,1] <- 1; glf[3,4,1] <- 1 # N A a A
  glf[1:3,1,2] <- 1; glf[3,2,2] <- 1; glf[3,3,2] <- 1; glf[3,4,2] <- 1   # N B B B
  glf[3,1,3] <- 1; glf[3,2,3] <- 1; glf[1,3,3] <- 1; glf[3,4,3] <- 1     # C C c C
  counts <- array(0, c(4,4,3))
  counts[1,1,1] <- 0; counts[2,2,1] <- 1; counts[1,3,1] <- 1; counts[2,4,1] <- 1
  counts[2,1,2] <- 0; counts[2,2,2] <- 1; counts[2,3,2] <- 1; counts[2,4,2] <- 1
  counts[2,1,3] <- 0; counts[2,2,3] <- 1; counts[1,3,3] <- 1; counts[2,4,3] <- 1
  pos <- cbind(0, c(1,2,3), rep(0,3), rep(1,3))
  Hd <- Haplodiplo$new(log(glf), counts, pos, rep(1,4))
  Hd$shf(rbind(c(0,1),c(0,2),c(1,2)), 0:3)
}

.test_hfs1d <- function()
{
  library(haplodiplo)
  library(Matrix)
  # haplotype likelihood (haploid, no missing data)
  set.seed(1)
  nsite <- 200
  glf <- array(0, c(3,9,nsite))
  counts <- array(0, c(4,9,nsite))
  for (i in 1:nsite)
  {
#    for (j in 1:3)
#    {
#      k <- sample(1:3, 1)
#      glf[k,j,i] <- 1
#      if (k == 1)
#        counts[1,j,i] <- 2
#      if (k == 2)
#        counts[1:2,j,i] <- 1
#      if (k == 3)
#        counts[2,j,i] <- 2
#    }
    for (j in 1:9)
    {
      k <- sample(c(1,3), 1)
      glf[k,j,i] <- 1
      if (k == 1)
        counts[1,j,i] <- 1
      if (k == 3)
        counts[2,j,i] <- 1
    }
  }
  pos <- cbind(0, 1:nsite, rep(0,nsite), rep(1,nsite))
  Hd <- Haplodiplo$new(log(glf), counts, pos, rep(1,9))
  pairs <- t(combn(0:99,2))
  shf1 <- Hd$shf(pairs, 0:8)
  block <- rep(0,nrow(pairs))
  hfs1 <- hfs1d(shf1$SHF, shf1$config, block, 0, max(shf1$config))
  hfs2 <- hfs1d(shf1$SHF, shf1$config, block, 0, 4)
  proj <- matrix(0, length(hfs2[[1]]), length(hfs1[[1]]))
  for (i in 1:length(hfs2[[1]])) for(j in 1:length(hfs1[[1]]))
    proj[i,j] <- extraDistr::dmvhyper(
      c(4-sum(hfs2[[2]][,i]),hfs2[[2]][,i]),
      c(max(shf1$config)-sum(hfs1[[2]][,j]),hfs1[[2]][,j]),
      4)

  # because each haplotype is selected uniformly at random
  # the number of counts per bin should follow multinomial sampling probabilities
  mcoef1 <- rep(NA, length(hfs1[[1]]))
  for (i in 1:length(hfs1[[1]]))
    mcoef1[i] <- multicool::multinom(c(max(shf1$config)-sum(hfs1[[2]][,i]),hfs1[[2]][,i]), counts=TRUE)
  plot(mcoef1/sum(mcoef1), hfs1[[1]]); abline(0,1)#this works

  mcoef2 <- rep(NA, length(hfs2[[1]]))
  for (i in 1:length(hfs2[[1]]))
    mcoef2[i] <- multicool::multinom(c(4-sum(hfs2[[2]][,i]),hfs2[[2]][,i]), counts=TRUE)
  plot(mcoef2/sum(mcoef2), hfs2[[1]]); abline(0,1)#this doesn't
  plot(mcoef2/sum(mcoef2), proj %*% hfs1[[1]]); abline(0,1)#but this does

  #how about...
  yow <- as(proj %*% shf1[[1]], "dgCMatrix")
  hfs2.2 <- hfs1d(yow, hfs2$config, block, 0, 4)
  plot(mcoef2/sum(mcoef2), hfs2.2[[1]])
  #nope

  EMstep <- function(par) { rowMeans(apply(yow, 2, function(x) x*par/sum(x*par))) }
  plot(mcoef2/sum(mcoef2), EMstep(EMstep(EMstep(EMstep(EMstep(rep(1,length(hfs2.2[[1]]))))))))
  #nope

  #what do we have
  #the SHF is p(X|G) p(G|a)
  #to project... we have to ask, what is p(a|a_proj)
  #... arg. we know that:
  #p(a_proj|a) = hypergeometric probability
  #p(a|a_proj) = p(a_proj|a) p(a)/sum p(a_proj|a) p(a)
  EMstep2 <- function(par) { 
    #rowMeans(apply(apply(proj, 2, function(x) x*par/sum(x*par)) %*% shf1[[1]], 2, function(z) z/sum(z)))
    rowMeans(apply(t(apply(proj, 1, function(x) x/sum(x))) %*% shf1[[1]], 2, function(z) (z*par)/sum(z*par)))
  }
  plot(mcoef2/sum(mcoef2), EMstep2(EMstep2(EMstep2(EMstep2(EMstep2(rep(1,length(hfs2.2[[1]]))))))))

  yow2 <- as(apply(proj, 2, function(x) (x/mcoef2)/sum(x/mcoef2)) %*% shf1[[1]], "dgCMatrix")
  hfs2.3 <- hfs1d(yow2, hfs2$config, block, 0, 4)

  # I think we just have to project after the fact, which totally blows.


}
