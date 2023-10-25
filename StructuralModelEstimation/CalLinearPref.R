# This calculates the preference parameters implied by a single type model. 
# See appendix in Dennis, Quintero and Sieg (2019).

CalLinearPref <- function(Py, Pv) {
  # INCOME t=1
  mu <- Py[1]
  sigma <- Py[2]
  beta <- Py[4]
  
  # VALUE t=1
  omega <- Pv[1]
  tau <- Pv[2]
  theta <- Pv[4]
  
  # We are changing the scale of theta and beta everywhere
  theta <- 1000 * theta
  beta <- 1000 * beta
  
  # PREFERENCES
  b1 <- sigma / tau
  a1 <- mu - (b1) * omega
  A1 <- exp(a1)
  
  # Create U
  U <- list()
  U$alpha <- 1 / (b1 - 1)
  U$phi <- 1 / A1
  U$eta <- theta
  U$gamma <- 1 - b1
  
  return(U)
}
