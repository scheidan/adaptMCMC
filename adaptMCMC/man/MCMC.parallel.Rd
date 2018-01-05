\name{MCMC.parallel}
\alias{MCMC.parallel}

\title{
Parallel computation of MCMC()
}
\description{
A wrapper function to generate several independent Markov chains by
stetting up cluster on a multi-core machine. The function is
based on the \pkg{parallel} package.
}
\usage{
MCMC.parallel(p, n, init, n.chain = 4, n.cpu, packages = NULL, dyn.libs=NULL,
    scale = rep(1, length(init)),  adapt = !is.null(acc.rate),
    acc.rate = NULL, gamma = 2/3, list = TRUE, ...)
}

\arguments{
  \item{p}{
    function that returns a value proportional to the log probability
    density to sample from. Alternatively the function can return a list
    with at least one element named \code{log.density}.
  }
  \item{n}{
    number of samples.
  }
  \item{init}{
    vector with initial values.
  }
  \item{n.chain}{
    number of independent chains.
  }
  \item{n.cpu}{
    number of CPUs that should be used in parallel.
  }
  \item{packages}{
    vector with name of packages to load into each instance. (Typically,
    all packages on which \code{p} depends.)
  }
  \item{dyn.libs}{
    vector with name of dynamic link libraries (shared objects) to load into each
    instance. The libraries must be located in the working directory.
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
    logical. If \code{TRUE} a list of lits is returned otherwise a list of matrices with the samples.
  }
  \item{\dots}{
    further arguments passed to \code{p}
  }

}

\details{
  This function is just a wrapper to use \code{MCMC} in parallel. It is
  based on \pkg{parallel}. Obviously, the application of this function
  makes only sense on a multi-core machine.
}
\value{
  A list with a list or matrix for each chain. See \code{\link{MCMC}} for details.
}

\author{Andreas Scheidegger, \email{andreas.scheidegger@eawag.ch} or
  \email{scheidegger.a@gmail.com}
}

\seealso{
  \code{\link{MCMC}}
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
## compute 4 independent chains on 2 CPU's (if available) in parallel

samp <- MCMC.parallel(p.log, n=200, init=c(x1=0, x2=1),
    n.chain=4, n.cpu=2, scale=c(1, 0.1),
    adapt=TRUE, acc.rate=0.234)

str(samp)

}
