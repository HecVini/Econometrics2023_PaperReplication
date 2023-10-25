# This function estimates the quantile of the GLN4 distribution. See appendix of Dennis, Quintero and Sieg (2019).

GLN4quantiles <- function(p, mediancdforig, mu, sigma, r, beta) {
  
  # Set the optimization options
  opts <- list(Algorithm = "interior-point", MaxFunEvals = 2000000000, MaxIter = 2000000000, TolX = 1e-30, TolFun = 1e-20)
  
  # Call the function
  source("cdfGLN4.R")
  
  # Define the objective function
  GLN4quant <- function(q, beta) {
    if (-1000 * beta < exp(q)) {
      F <- (cdfGLN4(exp(q), mediancdforig, mu, sigma, r, beta) - p)^2
    } else {
      F <- 1e15  # A large value to discourage q values that are too small
    }
    return(F)
  }
  
  qini <- max(-1000 * beta + 100, 5)
  # tive que usar outra fun??o de minimizacao que nao tem aqueles controles
  result <- optim(par = log(qini), fn = GLN4quant, beta = beta, method = "L-BFGS-B", control = opts)
  q_hat <- exp(result$par)
  fval <- result$value
  
  return(list(q_hat, fval))
}


