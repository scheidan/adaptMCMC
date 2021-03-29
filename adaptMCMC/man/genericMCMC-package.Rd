\name{adaptMCMC-package}
\alias{adaptMCMC-package}
\alias{adaptMCMC}
\docType{package}
\title{
Generic adaptive Monte Carlo Markov Chain sampler
}
\description{
Enables sampling from arbitrary distributions if the log density is known up to a constant; a common situation in the context of Bayesian inference. The implemented sampling algorithm was proposed by Vihola (2012) and achieves often a high efficiency by tuning the proposal distributions to a user defined acceptance rate.
}
\details{
\tabular{ll}{
Package: \tab adaptMCMC\cr
Type: \tab Package\cr
Version: \tab 1.4\cr
Date: \tab 2021-03-29\cr
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
Vihola, M. (2012) Robust adaptive Metropolis algorithm with coerced acceptance rate.
Statistics and Computing, 22(5), 997-1008. doi:10.1007/s11222-011-9269-5.
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
