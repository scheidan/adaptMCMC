## =======================================================
## (Adaptive) Metropolis Sampler
##
## Implementation of the RAM (robust adaptive Metropolis)
## sampler of
## Vihola, M. (2011) Robust adaptive Metropolis algorithm with
## coerced acceptance rate. Statistics and Computing.
## [online] http://www.springerlink.com/content/672270222w79h431/
## (Accessed December 8, 2011).

## Version 1.1.3

## June 21, 2017 -- Andreas Scheidegger
## =======================================================


MCMC <- function(p, n, init, scale=rep(1, length(init)),
                 adapt=!is.null(acc.rate), acc.rate=NULL, gamma=2/3, list=TRUE,
                 showProgressBar=interactive(), n.start=0, ...) {

  ## checks
  if(adapt & !is.numeric(acc.rate)) stop('Argument "acc.rate" is missing!')
  if(gamma<=0.5 | gamma>1) stop('Argument "gamma" must be in (0.5, 1]!')


  ## number of adaption steps
  if(is.numeric(adapt)) n.adapt <- adapt
  if(adapt==TRUE) n.adapt <- Inf
  if(adapt==FALSE) n.adapt <- 0

  ## number of parameter
  d <- length(init)

  ## matrix to store MC chain
  X <- matrix(NA, ncol=d, nrow=n)
  colnames(X) <- names(init)
  X[1,] <- init

  ## vector to store log densities p(x)
  p.val <- rep(NA, n)

  val <- p(X[1,], ...)
  if(is.list(val)) {
    returns.list <- TRUE
    extras <- list()                    # list to store additional return values of p
    if(!"log.density" %in% names(val)) {
      stop("The list returned by 'p' must contain an element named 'log.density!'")
    }
    if(length(val$log.density)>1) stop("The list element 'log.density' must be a scalar value!")
    
    p.val[1] <- val$log.density
    extras[[1]] <- val["log.density" != names(val)]
    
  } else {
    returns.list <- FALSE
    if(length(val)>1) stop("The function 'p' must return a scalar value or a named list! See ?MCMC.!")
    p.val[1] <- val
  }

  ## initial S
  if(d>1) {
    if(length(scale)==d) {
      M <- diag(scale)
    } else {
      M <- scale
    }
  } else {
    M <- matrix(scale)
  }

  ## check
  if(ncol(M) != length(init)) stop("Length or dimension of 'init' and 'scale' do not match!")

  S <-  t(chol(M))

  ## initialize progress bar
  cat('  generate', n, 'samples \n')
  if(showProgressBar){
    pb <- txtProgressBar(min=0, max=n, style=3)
  }
  update.step <- max(5, floor(n/100))

  k <- 0
  for(i in 2:n){

    if(showProgressBar && i %% update.step == 0) {
      setTxtProgressBar(pb, i)
    }

    ## proposal value
    U <- rt(d, df=d)
    X.prop <- c( X[i-1,] + S %*% U )
    names(X.prop) <- names(init)

    ## calculate density at X.prop
    val <- p(X.prop, ...)
    if(returns.list) {
      p.val.prop <- val$log.density
      extras.prop <- val["log.density" != names(val)]
    } else {
      p.val.prop <- val
    }
    
    ## acceptance probability
    alpha <- min(1, exp( p.val.prop - p.val[i-1] )) # for log density

    if(!is.finite(alpha)) alpha <- 0    # if zero divided by zero

    ## accept with P=alpha
    if(runif(1)<alpha) {
      X[i,] <- X.prop                   # accept
      p.val[i] <- p.val.prop
      if(returns.list) {
        extras[[i]] <- extras.prop
      }
      k <- k+1
    } else {
      X[i,] <- X[i-1,]                  # or not
      p.val[i] <- p.val[i-1]
      if(returns.list) {
        extras[[i]] <- extras[[i-1]]
      }
    }

    ## compute new S
    ii <- i+n.start
    if(ii < n.adapt) {
      
      ## ramcmc package performs rank 1 cholesky update/downdate as required
      S <- ramcmc::adapt_S(S, U, alpha, ii, acc.rate, gamma)

    }
  }

  if(showProgressBar){
    close(pb)                             # close progress bar
  }
  
  ## calculate accpetance rate
  acceptance.rate <- round(k/(n-1), 3)

  if(list) {
    res <- list(samples=X,
                log.p=p.val,
                cov.jump=S %*% t(S),
                n.sample=n,
                acceptance.rate=acceptance.rate,
                adaption=adapt,
                sampling.parameters=list(sample.density=p,
                                         acc.rate=acc.rate,
                                         gamma=gamma)
                )
    if(returns.list) {
      res$extra.values = extras
    }
    return(res)
  } else {
    cat("Acceptance rate:", acceptance.rate, "\n")
    return(X)
  }
}


## ----------------------
## Adds more samples to an existing chain
## does work with objet generated with MCMC.parallel()
## but updating is not computed  parallel.

MCMC.add.samples <- function(MCMC.object, n.update, ...) {

  ## if single chain
  if(!is.null(names(MCMC.object))) {
    if(is.matrix(MCMC.object)) stop("Only MCMC objects generated with option 'list=TRUE' can be updated!")

    ## get values from last sampling
    p <- MCMC.object$sampling.parameters$sample.density
    init <- MCMC.object$samples[nrow(MCMC.object$samples),]
    scale <- MCMC.object$cov.jump
    acc.rate <- MCMC.object$sampling.parameters$acc.rate
    gamma <- MCMC.object$sampling.parameters$gamma
    n.before <- MCMC.object$n.sample    # number of existing samples
    adapt <- MCMC.object$adaption

    ## generate new samples
    samp.update <- MCMC(p=p, n=n.update, init=init, scale=scale,  adapt=adapt, acc.rate=acc.rate,
                        gamma=gamma, list=TRUE, n.start=n.before, ...)

    ## update old sampling object
    MCMC.object$cov.jump <- samp.update$cov.jump
    m <- c(MCMC.object$n.sample, samp.update$n.sample)
    MCMC.object$acceptance.rate <-  1/sum(m)*(m[1]*MCMC.object$acceptance.rate + m[2]*samp.update$acceptance.rate)
    MCMC.object$n.sample <- MCMC.object$n.sample + n.update

    MCMC.object$samples <- rbind(MCMC.object$samples, samp.update$samples)
    MCMC.object$log.p <- c(MCMC.object$log.p, samp.update$log.p)
    if("extra.values" %in% names(MCMC.object)) {
      MCMC.object$extra.values <- c(MCMC.object$extra.values, samp.update$extra.values)
    }

    ## return the updated list
    return(MCMC.object)
  }

  ## if list of chains
  if(is.null(names(MCMC.object))) {
    ## recursive call of MCMC.add.samples() to update single chains
    MCMC.object <- lapply(MCMC.object, function(x) MCMC.add.samples(x, n.update=n.update, ...))
    return(MCMC.object)
  }
}

## ----------------------
## Wrapper for parallel calculation of independent chains
MCMC.parallel <- function(p, n, init, n.chain=4, n.cpu, packages=NULL, dyn.libs=NULL,
                          scale=rep(1, length(init)), adapt=!is.null(acc.rate),
                          acc.rate=NULL, gamma=2/3, list=TRUE, ...) {

  ## initialisation of (local) parallel computing
  cl <- makeCluster(min(n.cpu, detectCores()))

  ## stop parallel computing on exit
  on.exit({ stopCluster(cl); print("Cluster stopped.")})

  ## export complete work space of master
  varlist <- unique(c(ls(), ls(envir=.GlobalEnv), ls(envir=parent.env(environment()))))
  clusterExport(cl, varlist=varlist, envir=environment())

  ## init random generators
  clusterSetRNGStream(cl)

  ## export 'packages', 'dyn.libs', and current working directory
  wd <- getwd()
  clusterExport(cl, varlist=c("packages", "dyn.libs", "wd"), envir=environment())

  ## wrapper function to be called in parallel
  MCMC.wrap <- function(x, ...) {
    if(!is.null(packages)) sapply(packages, function(x) require(x, character.only=TRUE))
    ## load and unload dynamic libraries
    if (!is.null(dyn.libs)) {
      sapply(dyn.libs, function(x) dyn.load(paste(wd, x, sep = "/")))
      on.exit( sapply(dyn.libs, function(x) dyn.unload(paste(wd, x, sep = "/"))) )
    }
    MCMC(...)
  }


  ## sample n chains in parallel
  result <- clusterApply(cl, 1:n.chain, MCMC.wrap, p=p, n=n, init=init,
                         scale=scale,  adapt=adapt, acc.rate=acc.rate,
                         gamma=gamma, list=list, ...)

  return(result)

}


## ----------------------
## converts a sample into coda object

convert.to.coda <- function(sample) {
  ## if single chain
  if(!is.null(names(sample))) {
    if(is.matrix(sample)) {
      obj <- coda::mcmc(sample)
    }
    if(is.list(sample)) {
      obj <- coda::mcmc(sample$samples)
    }
    return(obj)
  } else {

    ## if sample is a list of chains
    if(is.matrix(sample[[1]])) {
      obj <- as.mcmc.list(lapply(sample, coda::mcmc))
    }
    if(is.list(sample[[1]])) {
      obj <- as.mcmc.list(lapply(sample, function(x) {coda::mcmc(x$samples)}))
    }
    return(obj)
  }
}
