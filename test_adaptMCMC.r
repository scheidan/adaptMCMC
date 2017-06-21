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
## check(package.path)


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

## ----------------------
## sampling one chain

n <- 2500
burn.in <- n/2

samp <- MCMC(p.log, n, init=rep(0,d), acc.rate=0.234, adapt=TRUE, showProgressBar=F)

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

means

Sigma

samp[[1]]$acceptance.rate

plot(convert.to.coda(samp))


