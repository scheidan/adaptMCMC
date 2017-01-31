adaptMCMC
=========

R package that provides an implementation of the generic adaptive Monte Carlo Markov chain sampler proposed by Vihola (2011).

# Geting started

```R
library(adaptMCMC)

## ----------------------
## Define (non-normalized) log density

## log-pdf to sample from
p.log <- function(x) {
  B <- 0.03                              # controls 'bananacity'
  -x[1]^2/200 - 1/2*(x[2]+B*x[1]^2-100*B)^2
}


## ----------------------
## generate samples

## 1) non-adaptive sampling
samp.1 <- MCMC(p.log, n=200, init=c(0, 1), scale=c(1, 0.1),
adapt=FALSE)

## 2) adaptive sampling
samp.2 <- MCMC(p.log, n=200, init=c(0, 1), scale=c(1, 0.1),
adapt=TRUE, acc.rate=0.234)


## ----------------------
## summarize results

str(samp.2)
summary(samp.2$samples)

## covariance of last jump distribution
samp.2$cov.jump
```