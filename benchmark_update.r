# This file benchmarks the cholesky update code:
# naive pure R vs ramcmc::update_s
# For the same set of inputs both methods should have the same resulting matrix
# The methods are also timed for a variety of dimensions

# Function to be compared against ramcmc::adapt_S
# this code is the update from the version of adaptMCMC on CRAN as of 2021-03-26
r_adapt_S <- function(S, U, alpha, ii, acc.rate, gamma) {
  
  d <- length(U)
  
  adapt.rate <-  min(1, d*ii^(-gamma))
  M <- S %*% (diag(d) + adapt.rate*(alpha - acc.rate) * U%*%t(U)/sum(U^2)) %*% t(S)
  
  ## check if M is positive definite. If not, use nearPD().
  eig <- eigen(M, only.values = TRUE)$values
  tol <- ncol(M)*max(abs(eig))*.Machine$double.eps
  
  if( !isSymmetric(M) | is.complex(eig) | !all(Re(eig)>tol) ){
    ## nearPD() computes the 'nearest' positive definite matrix
    M <- as.matrix(Matrix::nearPD(M)$mat)
  }
  
  S <- t(chol(M))
  
  S
}

# Define the run
n <- 1000
acc.rate <- 0.234
gamma <- 0.67
dims <- c(2,3,5, 10, 20, 100, 200, 300)

set.seed(123)

# Arrays to store times
r.times <- ramcmc.times <- c()

for (dim in dims) {
  # Generate some fake proposals and alphas to pass to each method
  U <- matrix(rnorm(n*dim), nrow=n, ncol=dim)
  alpha <- runif(n)
  
  # Each method starts with identity scale matrix
  S.r <- S.ramcmc <- diag(dim)
  
  # Pure R method
  r.start <- Sys.time()
  for (i in 1:n) {
    S.r <- r_adapt_S(S.r, U[i,], alpha[i], i, acc.rate, gamma)
  }
  r.end <- Sys.time()
  
  # ramcmc method
  ramcmc.start <- Sys.time()
  for (i in 1:n) {
    S.ramcmc <- ramcmc::adapt_S(S.ramcmc, U[i,], alpha[i], i, acc.rate, gamma)
  }
  ramcmc.end <- Sys.time()
  
  # How long did each take?
  cat("dim =", dim, "R time per update:", (r.end - r.start) / n, "\n")
  cat("dim =", dim, "ramcmc time per update:", (ramcmc.end - ramcmc.start) / n, "\n")
  r.times <- c(r.times, (r.end - r.start) / n)
  ramcmc.times <- c(ramcmc.times, (ramcmc.end - ramcmc.start) / n)
  
  # Equal within numerical tolerance?
  cat("dim =", dim, "Equal result?", 
      all.equal(S.r[lower.tri(S.r, diag=T)], S.ramcmc[lower.tri(S.ramcmc, diag=T)]), 
      "\n")
  
  cat("\n")
}

# Plot the times for comparison
barplot(rbind(r.times, ramcmc.times), beside=T, names.arg=dims,
        legend.text=c("R", "ramcmc"), args.legend=list(x="topleft"),
        xlab="Target dimension", ylab="Time per update (sec)",
        main="Comparison of Update Times")

barplot(r.times / ramcmc.times, names.arg=dims, 
        xlab="Target dimension", ylab="Speedup (R / ramcmc)",
        main="Update Time Ratio (Speedup)")
