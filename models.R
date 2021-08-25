genData <- function(N,T,model,parms=0) {
  burnin <- ceiling(T/2)
  obs <- burnin + N
  L <- N/T-1
  
  eps <- rnorm(obs)
  x <- rep(0, obs)
  
  if (model == 1) {
    x <- eps
  } else if (model == 2) {
    x[1] <- eps[1]
    for (i in 2:obs) {
      x[i] <- .5*x[i-1] + eps[i]
    }
  } else if (model == 3) {
    h <- rep(0, obs)
    h[1] <- .2
    x[1] <- sqrt(h[1])*eps[1]
    for (i in 2:obs) {
      h[i] <- .2 + .3*x[i-1]^2
      x[i] <- sqrt(h[i])*eps[i]
    }
  } else if (model == 4) {
    h <- rep(0, obs)
    h[1] <- .1
    x[1] <- sqrt(h[1])*eps[1]
    for (i in 2:obs) {
      h[i] <- .1 + .3*x[i-1]^2 + .3*h[i-1]
      x[i] <- sqrt(h[i])*eps[i]
    }
  } else if (model == 5) {
    h <- rep(0, obs)
    h[1] <- .1
    x[1] <- sqrt(h[1])*eps[1]
    for (i in 2:obs) {
      h[i] <- .1 + .7*eps[i-1]^2 + .3*h[i-1]
      x[i] <- sqrt(h[i])*eps[i]
    }
  } else if (model == 6) {
    eta <- rnorm(obs, 0, .1)
    beta <- rep(0, obs)
    x[1] <- eps[1]
    for (i in 2:obs) {
      beta[i] <- 0.5*beta[i-1] + eta[i]
      x[i] <- beta[i]*x[i-1] + eps[i]
    }
  } else if (model == 7) {
    x <- rt(obs, 1)
  } else if (model == 101) {
    x <- eps
    change.pt <- floor(burnin + T*(1 + L*runif(1, 0, 4/5)))
    x[change.pt:obs] <- x[change.pt:obs] + 1
  } else if (model == 102) {
    x <- eps
    change.pt <- floor(burnin + T*(1 + L*runif(1, 0, 4/5)))
    x[change.pt:obs] <- 2*x[change.pt:obs]
  } else if (model == 103) {
    change.pt <- floor(burnin + T*(1 + L*runif(1, 0, 4/5)))
    x <- c(1+sqrt(2)*eps[1:(change.pt-1)], eps[change.pt:obs]^2)
  } else if (model == 104) {
    tmp.x <- eps[(T+1):N]*exp(.5 - abs(.5 - (1:(L*T))/(L*T)))
    x <- c(rep(0, burnin), eps[1:T], tmp.x)
  } else if (model == 105) {
    library(stabledist)
    change.pt <- floor(burnin + T*(1 + L*runif(1, 0, 4/5)))
    xi <- rstable(obs - change.pt + 1, 1.5, 0.25)
    x <- c(eps[1:(change.pt-1)], xi)
  } else if (model == 106) {
    x <- rep(0, obs)
    change.pt <- burnin + floor(T*(1 + L/2))
    rho0 <- 0
    rhoA <- .5
    for (i in 2:(change.pt-1)) {
      x[i] <- rho0*x[i-1] + sqrt(1 - rho0^2)*eps[i]
    }
    for (i in change.pt:obs) {
      x[i] <- rhoA*x[i-1] + sqrt(1 - rhoA^2)*eps[i]
    }
  } else if (model == 107) {
    w <- rep(0, obs)
    change.pt <- burnin + floor(T*(1 + L/2))
    a <- 10
    rho0 <- 0.8
    x <- z <- rnorm(obs)
    for (i in seq(2, obs, by=2)) {
      x[i] <- (z[i] + a*z[i-1])/sqrt(1 + a^2)
    }
    x <- c(z[1:(change.pt-1)], x[change.pt:obs])
  } else if (model == 108) {
    w <- rep(0, obs)
    change.pt <- burnin + floor(T*(1 + L/2))
    a <- 10
    rho0 <- 0.8
    for (i in 2:obs) {
      w[i] <- rho0*w[i-1] + sqrt(1 - rho0^2)*eps[i]
    }
    x <- z <- rnorm(obs)
    for (i in seq(2, obs, by=2)) {
      x[i] <- (z[i] + a*w[i])/sqrt(1 + a^2)
    }
    x <- c(z[1:(change.pt-1)], x[change.pt:obs])
  }
  x <- x[-(1:burnin)]
}
