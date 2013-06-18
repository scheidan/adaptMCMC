\name{adaptMCMC-package}
\alias{adaptMCMC-package}
\alias{adaptMCMC}
\docType{package}
\title{
Generic adaptive Monte Carlo Markov Chain sampler
}
\description{
This package provides an implementation of the generic adaptive Monte
Carlo Markov chain sampler proposed by Vihola (2011).
}
\details{
\tabular{ll}{
Package: \tab adaptMCMC\cr
Type: \tab Package\cr
Version: \tab 1.1\cr
Date: \tab 2012-06-22\cr
License: \tab GPL (>= 2)\cr
LazyLoad: \tab yes\cr
}

The workhorse function is \code{\link{MCMC}}. Chains can be updated with
 \code{\link{MCMC.add.samples}}. \code{\link{MCMC.parallel}} is a
 wrapper to generate independent chains on several CPU's in parallel
 using \pkg{parallel}. \pkg{coda}-functions can be used after conversion
 with \code{\link{convert.to.coda}}.
}

\author{Andreas Scheidegger, \email{andreas.scheidegger@eawag.ch} or
  \email{scheidegger.a@gmail.com}}

\references{
Vihola, M. (2011) Robust adaptive Metropolis algorithm with
coerced acceptance rate. Statistics and Computing.
[online] \url{http://www.springerlink.com/content/672270222w79h431/},
(Accessed December 8, 2011)
}

\seealso{
\code{\link{MCMC}}, \code{\link{MCMC.add.samples}},
\code{\link{MCMC.parallel}}, \code{\link{convert.to.coda}}

  The package \code{HI} provides an adaptive rejection Metropolis sampler
  with the function \code{\link[HI]{arms}}. See also
  \code{\link[MHadaptive]{Metro_Hastings}} of the \code{MHadaptive} package.
}
% \examples{
% ##
% }
