## =======================================================
## Project: adaptMCMC
##
## Description: test cases for adaptMCMC package
##
## andreas.scheidegger@eawag.ch
## =======================================================


library(adaptMCMC)


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
## sampling oen chain

n <- 25000
burn.in <- n/2

samp <- MCMC(p.log, n, init=rep(0,d), acc.rate=0.234, adapt=TRUE)

## means
colMeans(samp$samples[-(1:burn.in),])

## Sigma
round(var(samp$samples[-(1:burn.in),]),1)


samp$acceptance.rate

plot(convert.to.coda(samp))


## ----------------------
## sampling parallel


samp <- MCMC.parallel(p.log, n, n.chain=5, n.cpu=5, init=rep(0,d),
                      acc.rate=0.234, adapt=TRUE, packages='mvtnorm')

means

Sigma

samp[[1]]$acceptance.rate

plot(convert.to.coda(samp))
