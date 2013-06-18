\name{MCMC.add.samples}
\alias{MCMC.add.samples}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Add samples to an existing chain.
}
\description{
  Add samples to an existing chain produced by \code{MCMC} or \code{MCMC.parallel}.
}
\usage{
MCMC.add.samples(MCMC.object, n.update, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{MCMC.object}{
    a list produced by \code{MCMC}  or \code{MCMC.parallel} with option \code{list = TRUE}.
  }
  \item{n.update}{
    number of additional samples.
  }
  \item{\dots}{
    further arguments passed to \code{p}.
  }
}

\details{
  Only objects generated with the option \code{list = TRUE} can be
  updated.

 A list of chains produced by \code{MCMC.parallel} can be
 updated. However, the calculations are \emph{not} performed in parallel
 (i.e. only a single CPU is used).
}

\value{
  A updated version of \code{MCMC.object}.
}

\author{Andreas Scheidegger, \email{andreas.scheidegger@eawag.ch} or
  \email{scheidegger.a@gmail.com}}


\seealso{
  \code{\link{MCMC}}, \code{\link{MCMC.parallel}}
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
               adapt=TRUE, acc.rate=0.234, list=TRUE)


## ----------------------
## add 200 to the existing chain
samp <- MCMC.add.samples(samp, n.update=200)

str(samp)
}

