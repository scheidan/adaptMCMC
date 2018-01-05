## =======================================================
## Project: adaptMCMC
##
## Description: test cases for adaptMCMC package
##
## andreas.scheidegger@eawag.ch
## =======================================================


library(devtools)

package.path = "adaptMCMC/"     # path must point to the folder containing the WaMaSim files

## simulate a new package installation
load_all(package.path)


## run R CMD check
check(package.path)

## build(package.path)


## -------------------------------------------------------
## multi variate normal

d <- 5                                  # number of dimensions
(means <- 2^(1:d -1 ))


## covariance matrix
S <- matrix(0, ncol=d, nrow=d)
diag(S) <- 1:d
S[lower.tri(S)] <- 1:sum(lower.tri(S))

Sigma <- S%*%t(S)

library(mvtnorm)

p.log <- function(x) {
  dmvnorm(x, means, Sigma, log=TRUE)
}

p.log.list <- function(x) {
  if(x[1]<0) {
    return (list(log.density=dmvnorm(x, means, Sigma, log=TRUE), x=x))
  } else {
    return (list(log.density=dmvnorm(x, means, Sigma, log=TRUE), x=x, extra="positive!"))
  } 
}

p.log.error <- function(x) {
  c(dmvnorm(x, means, Sigma, log=TRUE), 1)
}


## ----------------------
## sampling one chain

n <- 2500
burn.in <- n/2

samp <- MCMC(p.log, n, init=rep(0,d), acc.rate=0.234, adapt=TRUE, showProgressBar=T)
samp <- MCMC.add.samples(samp, 500)
str(samp)

samp <- MCMC(p.log.list, n, init=rep(0,d), acc.rate=0.234, adapt=TRUE, showProgressBar=T)
samp <- MCMC.add.samples(samp, 500)
str(samp)

samp <- MCMC(p.log.error, n, init=rep(0,d), acc.rate=0.234, adapt=TRUE, showProgressBar=T)


means
colMeans(samp$samples[-(1:burn.in),])

Sigma
round(var(samp$samples[-(1:burn.in),]),1)

samp$acceptance.rate

plot(convert.to.coda(samp))


## ----------------------
## sampling parallel


n <- 2500
burn.in <- n/2

samp <- MCMC.parallel(p.log, n, n.chain=3, n.cpu=3, init=rep(0,d),
                      acc.rate=0.234, adapt=TRUE, packages='mvtnorm')

str(samp)


samp <- MCMC.parallel(p.log.list, n, n.chain=3, n.cpu=3, init=rep(0,d),
                      acc.rate=0.234, adapt=TRUE, packages='mvtnorm')

str(samp)


