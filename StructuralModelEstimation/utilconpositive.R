#Calculates if utility of housing component is well defined for the normalization in equation (26)
utilconpositive <- function(up, hg) {
  # Order: alpha, phi, eta, gamma, kappa
  c <- 1 - up$phi * ((hg + up$eta) ^ up$gamma) - 0.001 # Compute nonlinear inequalities at x.
  # ceq <- []
  # Compute nonlinear equalities at x.
  return(c)
}