#
NULL

Rcpp::loadModule("Haplodiplo", TRUE)

#simulate data
sim <- function(seed, sites, ploidy, K, errmax, pmiss)
{
  set.seed(seed)
  err <- matrix(runif(length(ploidy)*sites, 0, errmax), length(ploidy), sites)
  Q <- matrix(runif(length(ploidy)*K), K, length(ploidy))
  Q <- matrix(apply(Q, 2, function(x) x/sum(x)), K, length(ploidy))
  F <- matrix(runif(sites*K), K, sites)
  h <- t(Q) %*% F
  gno <- h
  for (i in 1:nrow(h))
    gno[i,] <- rbinom(ncol(h), size=ploidy[i], h[i,])
  ogno <- gno
  for (i in 1:nrow(ogno))
    ogno[i,] <- abs(gno[i,]-rbinom(ncol(gno), ploidy[i], err[i,]))
  # likelihood
  lik0 <- lik1 <- lik2 <- matrix(0, length(ploidy), sites)
  for (i in 1:nrow(gno))
    if (ploidy[i] == 1)
    {
      miss <- which(rbinom(ncol(lik0), size=1, pmiss)==1)
      lik0[i,] <- dbinom(abs(0-ogno[i,]), 1, err[i,])
      lik2[i,] <- dbinom(abs(1-ogno[i,]), 1, err[i,])
      lik0[i,miss] <- lik2[i,miss] <- 1/2
    } else if (ploidy[i] == 2) {
      miss <- which(rbinom(ncol(lik0), size=1, pmiss)==1)
      lik0[i,] <- dbinom(abs(0-ogno[i,]), 2, err[i,])
      lik1[i,] <- dbinom(abs(1-ogno[i,]), 2, err[i,])
      lik2[i,] <- dbinom(abs(2-ogno[i,]), 2, err[i,])
      lik0[i,miss] <- lik1[i,miss] <- lik2[i,miss] <- 1/3
    }
  lik <- array(0, dim=c(3, length(ploidy), sites))
  lik[1,,] <- lik0
  lik[2,,] <- lik1
  lik[3,,] <- lik2

  #remove fixed
  bad <- (colSums(lik0)==0 & colSums(lik1)==0) | (colSums(lik0)==0 & colSums(lik2)==0) | (colSums(lik1)==0 & colSums(lik2)==0)
  gno <- gno[,!bad]
  ogno <- ogno[,!bad]
  F <- F[,!bad]
  err <- err[,!bad]
  h <- h[,!bad]
  lik0 <- lik0[,!bad]
  lik1 <- lik1[,!bad]
  lik2 <- lik2[,!bad]
  lik <- lik[,,!bad]

  #this would be the output from ANGSD
  likf <- c()
  for(i in 1:length(ploidy))
    likf <- cbind(likf, lik0[i,], lik1[i,], lik2[i,])

  return(list(likf=likf, lik=lik, gno=gno, obsgno=ogno, ploidy=ploidy, h=h, Q=Q, F=F, err=err))
}

#ddd <- c(rep(1,10),rep(2,10))
#sdf <- sim(1, ddd, sites=5000, K=2, errmax=0.0, pmiss=0.3)
#set.seed(1)
#subsit<-sample(0:(dim(sdf$lik)[3]-1), 3000)
#subsam<-sample(0:(length(ddd)-1), 10)
#writeBin(c(t(sdf$likf)), "test/test_data.gl", size=8)
#head(t(matrix(readBin("testDat.gl", "double", n=5000*3*7),3*7,5000)))
#library(Rcpp)
#set.seed(1)
#hi <- Haplodiplo$new(sdf$lik, ddd)
##test with all
#set.seed(1)
#moo1 <- hi$admixture(0:(dim(sdf$lik)[3]-1), 0:(length(ddd)-1), 1, 4)
#set.seed(1)
#moo2  <- hi$admixture(0:(dim(sdf$lik)[3]-1), 0:(length(ddd)-1), 2, 4)
#set.seed(1)
#moo3 <- hi$admixture(0:(dim(sdf$lik)[3]-1), 0:(length(ddd)-1), 3, 4)
#set.seed(1)
#moo4 <- hi$admixture(0:(dim(sdf$lik)[3]-1), 0:(length(ddd)-1), 4, 4)
#par(mfrow=c(1,2))
#plot(sdf$F[2:1,], moo$F)
#plot(sdf$Q, moo$Q)

##test with subset of sites, samples
#set.seed(1)
#fr <- hi$frequencies(0:(dim(sdf$lik)[3]-1), 0:(length(ddd)-1))
#moo2 <- hi$admixture(subsit, subsam, matrix(length(subsam)*2,2,length(subsam)))
#par(mfrow=c(1,2))
#plot(sdf$F[,subsit+1], moo2$F)
#plot(sdf$Q[,subsam+1], moo2$Q)
#
#set.seed(1)
#hi2 <- haplodiplo$new(sdf$lik[,subsam+1,subsit+1], ddd[subsam+1])
#moo3 <- hi2$admixture(0:(length(subsit)-1), 0:(length(subsam)-1), 2, 0)
#par(mfrow=c(1,2))
#plot(moo2$F, moo3$F)
#plot(moo2$Q, moo3$Q)
#OK

#frequencies
#ddd2 <- c(rep(1,10),rep(2,10))
#ahh <- sim(1, ddd2, sites=100, K=1, errmax=0.1, pmiss=0.1)
#ho <- haplodiplo$new(ahh$lik, ddd2)
#boo <- ho$frequencies(0:(dim(ahh$lik)[3]-1), 0:(length(ddd2)-1))

#saf
#ddd2 <- c(rep(1,2),rep(2,2))
#ahh <- sim(1, ddd2, sites=4, K=1, errmax=0.1, pmiss=0.2)
#ho <- haplodiplo$new(ahh$lik, ddd2)
#boo <- ho$saf(0:(dim(ahh$lik)[3]-1), 0:(length(ddd2)-1))

# table(replicate(1000, sum(apply(ahh$lik[,,3],2,function(x) which(rmultinom(1,size=1, prob=x)>0)-1))))/1000

# genotype likelihoods for angsd equiv test
# wha <- read.table("../angsd_equiv/example.beagle.gz", header=T)
# lik <- array(NA, c(3,2,4))
# for(i in 1:nrow(wha))
#   lik[,,i] <- unlist(wha[i,4:9])
# ho <- haplodiplo$new(lik, c(2,2))
# ho$saf(0:3, 0:1)
# #more sites
# lik2 <- array(NA, c(3,2,4*100))
# for(i in 1:nrow(wha))
#   for(r in 1:100)
#     lik2[,,(r-1)*4 + i] <- unlist(wha[i,4:9])
# ho <- haplodiplo$new(lik2, c(2,2))
# ho$saf(0:399, 0:1)

# TODO should come up with test for haplodiploid. 

# #covariance
# hi <- haplodiplo$new(sdf$lik, ddd) 
# fr <- hi$frequencies(0:4999, 0:19)
# ok <- fr$pval < 0.9
# foo0 <- eigen(hi$covar(which(ok)-1, 0:19, sdf$h[,ok], fr$freq[ok])[[1]])$vectors #this will EXPLODE if super low frequency alleles are included
# fr2 <- t(sapply(1:20, function(x) fr$freq))
# foo1 <- eigen(hi$covar(0:4999, 0:19, fr2, fr$freq)[[1]])$vectors #this seems far less likely to explode
# plot(sdf$Q[1,], foo0[,1])
# plot(sdf$Q[1,], foo1[,1])
# plot(foo0[,1], foo0[,2])
# plot(foo1[,1], foo1[,2])


####### real data
# load("debug_pepo.RData"); hi <- Master$new(lik, rep(1,96))
# hi$paralogs(0:9, 0:95)
