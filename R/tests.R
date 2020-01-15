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

.simulate_glf <- function(seed=0, nsnp=1000, ploidy=c(1,2,1,2), pseg=0.1, perr=0.01,depth=10)
{
  set.seed(seed)

  pder <- rbeta(nsnp, 1, 5)

  glik <- function(cnt, err)
    c(dbinom(cnt[1], sum(cnt), 1-err), dbinom(cnt[1], sum(cnt), 0.5), dbinom(cnt[1], sum(cnt), err))

  tglf <- array(0, dim=c(3,length(ploidy),nsnp))
  glf <- array(0, dim=c(3,length(ploidy),nsnp))
  counts <- array(0, dim=c(4,length(ploidy),nsnp))
  for(i in 1:nsnp)
  {
    if (runif(1) > pseg) #nonsegregating
    {
      for(j in 1:length(ploidy))
      {
        tglf[1,j,i] <- 1
        counts[1:2,j,i] <- rpois(2,depth*c(1-perr,perr))
        glf[,j,i] <- glik(counts[,j,i], perr)
      }
    }
    else #nonsegregating
    {
      for(j in 1:length(ploidy))
      {
        if (ploidy[j]==2) #diploid
        {
          tglf[sample(1:3,1,prob=c((1-pder[i])^2,2*pder[i]*(1-pder[i]),pder[i]^2)),j,i] <- 1
          if (tglf[1,j,i]) 
            counts[1:2,j,i] <- rpois(2,depth*c(1-perr,perr))
          else if (tglf[2,j,i]) #het
            counts[1:2,j,i] <- rpois(2,depth*0.5)
          else
            counts[1:2,j,i] <- rpois(2,depth*c(perr,1-perr))
          glf[,j,i] <- glik(counts[,j,i], perr)
        } else { #haploid
          tglf[sample(c(1,3),1,prob=c(1-pder[i],pder[i])),j,i] <- 1
          if (tglf[1,j,i]) 
            counts[1:2,j,i] <- rpois(2,depth*c(1-perr,perr))
          else
            counts[1:2,j,i] <- rpois(2,depth*c(perr,1-perr))
          glf[,j,i] <- glik(counts[,j,i], perr)
        }
      }
    }
  }

  return(list(tglf=tglf, glf=glf, counts=counts, ploidy=ploidy))
}

.test_shfwindow <- function(seed=0, perr=0., depth=100, nsnp=10)
{
  library(haplodiplo)
  library(Matrix)

  # high depth, haploid, no error
  ploidy <- c(1,1,1,1,1,1)
  w1 <- .simulate_glf (seed=seed+1, ploidy=ploidy, perr=perr, pseg=0.25, nsnp=nsnp, depth=depth)
  w2 <- .simulate_glf (seed=seed+2, ploidy=ploidy, perr=perr, pseg=0.25, nsnp=nsnp, depth=depth)
  w3 <- .simulate_glf (seed=seed+3, ploidy=ploidy, perr=perr, pseg=0.25, nsnp=nsnp, depth=depth)
  p1 <- cbind(0, 1:nsnp, 0, 1)
  p2 <- cbind(0, 50001:(50000+nsnp), 0, 1)
  p3 <- cbind(1, 1:nsnp, 0, 1)
  glf <- array(c(w1$glf, w2$glf, w3$glf), c(3,6,nsnp*3))
  counts <- array(c(w1$counts, w2$counts, w3$counts), c(4,6,nsnp*3))
  pos <- rbind(p1,p2,p3)

  Hd <- Haplodiplo$new(log(glf), counts, pos, ploidy)
  rmap <- cbind(c(0,0,0,1,1,1), c(0,40000,80000,0,40000,80000), c(1.,1.,1.,1.,1.,1), c(0.0,0.04,0.08,0.0,0.04,0.08))
  bins <- cbind(c(0.0,0.01),c(0.01,Inf)) #what if bin is empty; how does Shf struct do

  shfwindow <- Hd$shfwindow(0:5, rmap, bins, 0.5)
  hfs_w1w2 <- hfs1d(shfwindow$SHF[[2]], shfwindow$weights[[2]])
  hfs_w1w1_w2w2_w3w3 <- hfs1d(shfwindow$SHF[[1]], shfwindow$weights[[1]])

  #to check...count haplotypes in simulated data
  thfs_w1w2 <- array(0, c(7,7,7))
  for(i in 1:nsnp)
    for(j in 1:nsnp)
    {
      bin <- c("Ab"=0, "aB"=0, "AB"=0)
      for(k in 1:6)
      {
        if(w1$tglf[3,k,i] & w2$tglf[1,k,j])
          bin["Ab"] = bin["Ab"] + 1
        else if(w1$tglf[1,k,i] & w2$tglf[3,k,j])
          bin["aB"] = bin["aB"] + 1
        else if(w1$tglf[3,k,i] & w2$tglf[3,k,j])
          bin["AB"] = bin["AB"] + 1
      }
      thfs_w1w2[bin[1]+1,bin[2]+1,bin[3]+1] <- thfs_w1w2[bin[1]+1,bin[2]+1,bin[3]+1] + 1
    }

  thfs_w1w1 <- array(0, c(7,7,7))
  for(i in 1:(nsnp-1))
    for(j in (i+1):nsnp)
    {
      bin <- c("Ab"=0, "aB"=0, "AB"=0)
      for(k in 1:6)
      {
        if(w1$tglf[3,k,i] & w1$tglf[1,k,j])
          bin["Ab"] = bin["Ab"] + 1
        else if(w1$tglf[1,k,i] & w1$tglf[3,k,j])
          bin["aB"] = bin["aB"] + 1
        else if(w1$tglf[3,k,i] & w1$tglf[3,k,j])
          bin["AB"] = bin["AB"] + 1
      }
      thfs_w1w1[bin[1]+1,bin[2]+1,bin[3]+1] <- thfs_w1w1[bin[1]+1,bin[2]+1,bin[3]+1] + 1
    }

  thfs_w2w2 <- array(0, c(7,7,7))
  for(i in 1:(nsnp-1))
    for(j in (i+1):nsnp)
    {
      bin <- c("Ab"=0, "aB"=0, "AB"=0)
      for(k in 1:6)
      {
        if(w2$tglf[3,k,i] & w2$tglf[1,k,j])
          bin["Ab"] = bin["Ab"] + 1
        else if(w2$tglf[1,k,i] & w2$tglf[3,k,j])
          bin["aB"] = bin["aB"] + 1
        else if(w2$tglf[3,k,i] & w2$tglf[3,k,j])
          bin["AB"] = bin["AB"] + 1
      }
      thfs_w2w2[bin[1]+1,bin[2]+1,bin[3]+1] <- thfs_w2w2[bin[1]+1,bin[2]+1,bin[3]+1] + 1
    }

  thfs_w3w3 <- array(0, c(7,7,7))
  for(i in 1:(nsnp-1))
    for(j in (i+1):nsnp)
    {
      bin <- c("Ab"=0, "aB"=0, "AB"=0)
      for(k in 1:6)
      {
        if(w3$tglf[3,k,i] & w3$tglf[1,k,j])
          bin["Ab"] = bin["Ab"] + 1
        else if(w3$tglf[1,k,i] & w3$tglf[3,k,j])
          bin["aB"] = bin["aB"] + 1
        else if(w3$tglf[3,k,i] & w3$tglf[3,k,j])
          bin["AB"] = bin["AB"] + 1
      }
      thfs_w3w3[bin[1]+1,bin[2]+1,bin[3]+1] <- thfs_w3w3[bin[1]+1,bin[2]+1,bin[3]+1] + 1
    }

  thfs_w1w1_w2w2_w3w3 <- thfs_w1w1 + thfs_w2w2 + thfs_w3w3

  thfs_w1w1_w2w2_w3w3_lin <- rep(0, ncol(shfwindow$config[[1]]))
  thfs_w1w2_lin <- rep(0, ncol(shfwindow$config[[1]]))
  for (i in 1:ncol(shfwindow$config[[1]]))
  {
    bin <- shfwindow$config[[1]][,i]
    thfs_w1w1_w2w2_w3w3_lin[i] <- thfs_w1w1_w2w2_w3w3[bin[1]+1,bin[2]+1,bin[3]+1]
    thfs_w1w2_lin[i] <- thfs_w1w2[bin[1]+1,bin[2]+1,bin[3]+1]
  }
  thfs_w1w1_w2w2_w3w3_lin <- thfs_w1w1_w2w2_w3w3_lin / sum(thfs_w1w1_w2w2_w3w3_lin)
  thfs_w1w2_lin <- thfs_w1w2_lin / sum(thfs_w1w2_lin)
  
  par(mfrow=c(2,1))
  plot(thfs_w1w1_w2w2_w3w3_lin, hfs_w1w1_w2w2_w3w3); abline(0,1)
  plot(thfs_w1w2_lin, hfs_w1w2); abline(0,1)

  #browser()

  print("stop");
}

.test_robinsonhill <- function(seed=0)
{
  set.seed(seed)

  # test single population
  nsim <- 1e6
  p0 <- runif(nsim, 0.2, 0.8)
  q0 <- runif(nsim, 0.2, 0.8)
  Fst <- 0.1
  p <- q <- D <- f <- HFS <- list()
  mn <- c(0.3, -0.2, 0.)
  sd <- c(0.1, 0.001, 0.3)
  chr <- c(4,5,6)
  for(i in 1:3)
  {
    p[[i]] <- rbeta(nsim, (1-Fst)/Fst * p0, (1-Fst)/Fst * (1-p0) )
    q[[i]] <- rbeta(nsim, (1-Fst)/Fst * q0, (1-Fst)/Fst * (1-q0) )
    D[[i]] <- truncnorm::rtruncnorm(nsim, pmax(-p[[i]]*q[[i]], -(1-p[[i]])*(1-q[[i]])), pmin(p[[i]]*(1-q[[i]]), (1-p[[i]])*q[[i]]), mn[i], sd[i])
    f[[i]] <- matrix(NA,4,nsim); rownames(f[[i]]) <- c("ab", "Ab", "aB", "AB")
    f[[i]]["ab",] <- (1-p[[i]])*(1-q[[i]]) + D[[i]]
    f[[i]]["Ab",] <- p[[i]]*(1-q[[i]]) - D[[i]]
    f[[i]]["aB",] <- (1-p[[i]])*q[[i]] - D[[i]]
    f[[i]]["AB",] <- p[[i]]*q[[i]] + D[[i]]
    HFS[[i]] <- apply(f[[i]], 2, function(x) paste(rmultinom(1, chr[i], x),collapse=""))
  }
  HFS3d <- table(HFS[[1]], HFS[[2]], HFS[[3]])
  HFS3d <- HFS3d/sum(HFS3d)

  configs <- list()
  for(i in 1:3)
    configs[[i]] <- sapply(dimnames(HFS3d)[[i]], function(x) as.numeric(strsplit(x, "")[[1]][2:4]))

  est <- rbind(basis(HFS3d, configs[[1]], configs[[2]], configs[[3]]),
               basis(cycle_cube(HFS3d), configs[[2]], configs[[3]], configs[[1]]),
               basis(cycle_cube(cycle_cube(HFS3d)), configs[[3]], configs[[1]], configs[[2]]))

  tru <- rbind(
         c("Di" = mean(D[[1]]), 
           "DiDi" = mean(D[[1]]^2), 
           "DiDj" = mean(D[[1]]*D[[2]]), 
           "DiDk" = mean(D[[1]]*D[[3]]), 
           "DiZii" = mean(D[[1]]*(1-2*p[[1]])*(1-2*q[[1]])),
           "DiZij" = mean(D[[1]]*(1-2*p[[1]])*(1-2*q[[2]])),
           "DiZji" = mean(D[[1]]*(1-2*p[[2]])*(1-2*q[[1]])),
           "DiZik" = mean(D[[1]]*(1-2*p[[1]])*(1-2*q[[3]])),
           "DiZki" = mean(D[[1]]*(1-2*p[[3]])*(1-2*q[[1]])),
           "DiZjk" = mean(D[[1]]*(1-2*p[[2]])*(1-2*q[[3]])),
           "DiZkj" = mean(D[[1]]*(1-2*p[[3]])*(1-2*q[[2]]))
           ),
         c("Di" = mean(D[[2]]), 
           "DiDi" = mean(D[[2]]^2), 
           "DiDj" = mean(D[[2]]*D[[3]]), 
           "DiDk" = mean(D[[2]]*D[[1]]), 
           "DiZii" = mean(D[[2]]*(1-2*p[[2]])*(1-2*q[[2]])),
           "DiZij" = mean(D[[2]]*(1-2*p[[2]])*(1-2*q[[3]])),
           "DiZji" = mean(D[[2]]*(1-2*p[[3]])*(1-2*q[[2]])),
           "DiZik" = mean(D[[2]]*(1-2*p[[2]])*(1-2*q[[1]])),
           "DiZki" = mean(D[[2]]*(1-2*p[[1]])*(1-2*q[[2]])),
           "DiZjk" = mean(D[[2]]*(1-2*p[[3]])*(1-2*q[[1]])),
           "DiZkj" = mean(D[[2]]*(1-2*p[[1]])*(1-2*q[[3]]))
           ),
         c("Di" = mean(D[[3]]), 
           "DiDi" = mean(D[[3]]^2), 
           "DiDj" = mean(D[[3]]*D[[1]]), 
           "DiDk" = mean(D[[3]]*D[[2]]), 
           "DiZii" = mean(D[[3]]*(1-2*p[[3]])*(1-2*q[[3]])),
           "DiZij" = mean(D[[3]]*(1-2*p[[3]])*(1-2*q[[1]])),
           "DiZji" = mean(D[[3]]*(1-2*p[[1]])*(1-2*q[[3]])),
           "DiZik" = mean(D[[3]]*(1-2*p[[3]])*(1-2*q[[2]])),
           "DiZki" = mean(D[[3]]*(1-2*p[[2]])*(1-2*q[[3]])),
           "DiZjk" = mean(D[[3]]*(1-2*p[[1]])*(1-2*q[[2]])),
           "DiZkj" = mean(D[[3]]*(1-2*p[[2]])*(1-2*q[[1]]))
           ))
  
  plot(tru, est); abline(0,1)
  
  cat("stop\n")
}

