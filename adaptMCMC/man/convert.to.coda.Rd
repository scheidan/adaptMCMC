\name{convert.to.coda}
\alias{convert.to.coda}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Converts chain(s) into \pkg{coda} objects.
}
\description{
  Converts chain(s) produced by \code{MCMC} or \code{MCMC.parallel} into
  \pkg{coda} objects.
}
\usage{
convert.to.coda(sample)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sample}{
    output of \code{MCMC} or \code{MCMC.parallel}.
  }
}

\details{
  Converts chain(s) produced by \code{MCMC} or \code{MCMC.parallel} so
  that they can be used with functions of the \pkg{coda} package.
}

\value{
 An object of the class \code{mcmc} or \code{mcmc.list}.
}

\author{Andreas Scheidegger, \email{andreas.scheidegger@eawag.ch} or
  \email{scheidegger.a@gmail.com}}

\seealso{
  \code{\link{MCMC}}, \code{\link[coda]{mcmc}}, \code{\link[coda]{mcmc.list}}
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
## generate 200  samples

samp <- MCMC(p.log, n=200, init=c(0, 1), scale=c(1, 0.1),
               adapt=TRUE, acc.rate=0.234)


## ----------------------
## convert in object of class 'mcmc'
samp.coda <- convert.to.coda(samp)

class(samp.coda)

## ----------------------
## use functions of package 'coda'

require(coda)

plot(samp.coda)
cumuplot(samp.coda)
}
