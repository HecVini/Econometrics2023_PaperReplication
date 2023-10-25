# This function evaluates the GLN4 cdf distribution. See appendix of Dennis, Quintero and Sieg (2019).

cdfGLN4 <- function(x, mediancdforig, mu, sigma, r, beta) {
  # We are changing the scale of theta and beta everywhere
  beta <- 1000 * beta
  
  # Initialize variables
  uu <- mu
  ss <- sigma
  rr <- r
  lengthx <- length(x)
  v <- 1 / r
  menor <- 0
  mayor <- 0
  
  # Create vectors to store the results
  B <- numeric(lengthx)
  T <- numeric(lengthx)
  F <- numeric(lengthx)
  M <- numeric(lengthx)
  
  # Loop through the input vector
  for (i in 1:lengthx) {
    # If x[i] + beta is less than mediancdforig, the cdf is less than 0.5
    if (x[i] + beta < mediancdforig) {
      menor <- 1
      B[i] <- ((mu - log(x[i] + beta) / sigma)^r) / r
      T[i] <- B[i]
      F[i] <- gamma(v) * pgamma(B[i], shape = v, lower.tail = FALSE) / (2 * gamma(v))
    }
    
    # If x[i] + beta is equal to mediancdforig, the cdf is equal to 0.5
    if (x[i] + beta == mediancdforig) {
      med <- 1
      F[i] <- 0.5
    }
    
    # If x[i] + beta is greater than mediancdforig, the cdf is greater than 0.5
    if (x[i] + beta > mediancdforig) {
      mayor <- 1
      M[i] <- ((log(x[i] + beta) - mu) / sigma)^r / r
      T[i] <- M[i]
      F[i] <- 0.5 + gamma(v) * pgamma(M[i], shape = v) / (2 * gamma(1/r))
    }
  }
  
  # Return the cdf vector
  return(F)
}


