\name{MCMC}
\alias{MCMC}

\title{
(Adaptive) Metropolis Sampler
}

\description{
Implementation of the robust adaptive Metropolis sampler of Vihola (2012).
}

\usage{
MCMC(p, n, init, scale = rep(1, length(init)),
    adapt = !is.null(acc.rate), acc.rate = NULL, gamma = 0.5,
    list = TRUE, showProgressBar=interactive(), n.start = 0, ...)
}

\arguments{
  \item{p}{
    function that returns a value proportional to the log probability
    density to sample from. Alternatively it can be a function that returns a list
    with at least one element named \code{log.density}. See details below.
  }
  \item{n}{
    number of samples.
  }
  \item{init}{
    vector with initial values.
  }
  \item{scale}{
    vector with the variances \emph{or} covariance matrix of the jump distribution.
  }
  \item{adapt}{
    if \code{TRUE}, adaptive sampling is used, if \code{FALSE} classic metropolis
    sampling, if a positive integer the adaption stops after \code{adapt} iterations.
  }
  \item{acc.rate}{
    desired acceptance rate (ignored if \code{adapt=FALSE})
  }
  \item{gamma}{
    controls the speed of adaption. Should be between 0.5 and 1. A lower
    gamma leads to faster adaption.
  }
  \item{list}{
    logical. If \code{TRUE} a list is returned otherwise only a matrix with the samples.
  }
  \item{showProgressBar}{
    logical. If \code{TRUE} a progress bar is shown.
  }
  \item{n.start}{
    iteration where the adaption starts. Only internally used.
  }
  \item{\dots}{
    further arguments passed to \code{p}.
  }
}
\details{

  The algorithm tunes the covariance matrix of the (normal) jump
  distribution to achieve the desired acceptance rate. Classic
  (non-adaptive) Metropolis sampling can be obtained by setting \code{adapt=FALSE}.

  Note, due to the calculation for the adaption steps the sampler is
  rather slow. However, with a suitable jump distribution good mixing can
  be observed with less samples. This is crucial if
  the computation of \code{p} is slow.

  In some cases the function \code{p} may not only calculate the log
  density but return a list containing also other values. For example
  if \code{p} is a log posterior one may be also interested to store
  the corresponding prior and likelihood values. The function must
  either return always a scalar or always a list, however, the length of
  the list may vary.


}
\value{
  If \code{list=FALSE} a matrix is with the samples.

  If \code{list=TRUE} a list is returned with the following components:
  \item{samples}{matrix with samples}
  \item{log.p}{vector with the (unnormalized) log density for each sample}
  \item{n.sample}{number of generated samples}
  \item{acceptance.rate}{acceptance rate}
  \item{adaption}{either logical if adaption was used or not, or the
    number of adaption steps.}
  \item{sampling.parameters}{a list with further sampling
    parameters. Mainly used by \code{MCMC.add.samples()}}.
  \item{extra.values}{A list containing additional return values provided by
    \code{p}. Only if \code{p} provides a list.}
}

\references{
Vihola, M. (2012) Robust adaptive Metropolis algorithm with coerced acceptance rate.
Statistics and Computing, 22(5), 997-1008. doi:10.1007/s11222-011-9269-5.
}

\author{Andreas Scheidegger, \email{andreas.scheidegger@eawag.ch} or
  \email{scheidegger.a@gmail.com}.

  Thanks to David Pleydell and Venelin for spotting errors and providing
  improvements.
}

\note{Due to numerical errors it may happen that the computed
  covariance matrix is not positive definite. In such a case the nearest positive
  definite matrix is calculated with \code{nearPD()}
  of the package \pkg{Matrix}.
}

\seealso{
  \code{\link{MCMC.parallel}}, \code{\link{MCMC.add.samples}}

  The package \code{HI} provides an adaptive rejection Metropolis sampler
  with the function \code{\link[HI]{arms}}. See also
  \code{\link[MHadaptive]{Metro_Hastings}} of the \code{MHadaptive} package.
}
\examples{
## ----------------------
## Banana shaped distribution

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


## ----------------------
## plot density and samples

x1 <- seq(-15, 15, length=80)
x2 <- seq(-15, 15, length=80)
d.banana <- matrix(apply(expand.grid(x1, x2), 1,  p.log), nrow=80)

par(mfrow=c(1,2))
image(x1, x2, exp(d.banana), col=cm.colors(60), asp=1, main="no adaption")
contour(x1, x2, exp(d.banana), add=TRUE, col=gray(0.6))
lines(samp.1$samples, type='b', pch=3)

image(x1, x2, exp(d.banana), col=cm.colors(60), asp=1, main="with adaption")
contour(x1, x2, exp(d.banana), add=TRUE, col=gray(0.6))
lines(samp.2$samples, type='b', pch=3)



## ----------------------
## function returning extra information in a list


p.log.list <- function(x) {
  B <- 0.03                              # controls 'bananacity'
  log.density <- -x[1]^2/200 - 1/2*(x[2]+B*x[1]^2-100*B)^2

  result <- list(log.density=log.density)

  ## under some conditions one may want to return other infos
  if(x[1]<0) {
    result$message <- "Attention x[1] is negative!"
    result$x <- x[1]
  }

  result
}

samp.list <- MCMC(p.log.list, n=200, init=c(0, 1), scale=c(1, 0.1),
                  adapt=TRUE, acc.rate=0.234)

## the additional values are stored under `extras.values`
head(samp.list$extras.values)

}