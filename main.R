rm(list=ls(all=TRUE))
setwd('C:/TempGit/seqStationarityTestECF')
library(Rcpp)

source('models.R')
sourceCpp('tstat.cpp')

T <- 100
L <- 1
m <- 1

model <- 1

a <- c(.5,1,1.5)
N <- T*(1+L)

# Bootstrap block and subsampling size
KT <- max(5, sqrt(log10(T)))
lambda <- function(t) {
  abst <- abs(t)
  (abst <= .5) + 2*(1 - abst)*((abst > .5) & (abst <= 1))
}

mc <- 1000
alpha <- .05

na <- length(a)

Tn <- matrix(0, nrow=mc, ncol=na)
Tn.star <- matrix(0, nrow=mc, ncol=na)

rej0 <- rep(0,na)
names(rej0) <- a

set.seed(20104480)

startTime <- Sys.time()
for (i in 1:mc) {
  x <- genData(N, T, model)
  # plot(x,type='l')
  
  R <- as.vector(acf(x[1:T], lag.max=2*KT, plot=FALSE)[[1]])
  tmp <- which(abs(R[1:KT]) < 2*sqrt(log10(T)/T))
  if (length(tmp) > 0) {
    M <- 2*(tmp[1] - 1)
  } else {
    M <- 2*(KT - 1)
  }
  ghat <- sum(lambda((-M:M)/M)*R[abs(-M:M)+1])
  Ghat <- sum(lambda((-M:M)/M)*R[abs(-M:M)+1]*abs(-M:M))
  D.CB <- 4/3*ghat^2
  blocksize <- max(1, round((2*Ghat^2/D.CB*T)^(1/3)))
  
  x.star <- c()
  while (length(x.star) < N) {
    idx.start <- sample.int(T, 1, replace=TRUE)
    idx.len <- rgeom(1, 1/blocksize) + 1
    idx <- idx.start:(idx.start + idx.len)
    idx <- (idx - 1) %% T + 1
    
    x.star <- c(x.star, x[idx])
  }
  
  rho.hat <- acf(x,lag.max=1,plot=FALSE)$acf[[2]]
  pT <- min(floor((3*T/2)^(1/3)*(2*rho.hat/(1-rho.hat^2))^(2/3)),floor(8*(T/100)^(1/3)))
  
  Tn[i,] <- apply(tstats(x,T,m,a),1,max)
  Tn.star[i,] <- apply(tstats(x.star,T,m,a),1,max)
  
  if (i %% 20 == 0) {
    critvals <- apply(Tn.star[1:i,], 2, function(x) { quantile(x, 1-alpha) })
    
    for (k in 1:na) {
      rej0[k] <- sum(critvals[k] < Tn[1:i,k])
    }
    
    cat('\n')
    rej <- round(100*rej0/i,1)
    print(rej)
  }
  
  elapsed = difftime(Sys.time(),startTime,units='mins')
  remaining = round((mc-i)*elapsed/i)
  cat('\r', floor(i/mc*100),'% complete (',remaining,' min remaining)                        ',sep='')
  flush.console()
}

cat('\n')


