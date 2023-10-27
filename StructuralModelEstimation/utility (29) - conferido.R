# Evaluate the utility function in equation 26 for an income of y, a
# consumption of housing of h at prices v
utility <- function(y,h,phi,gamma,alpha,eta,kappa, v){
  A <- log(1-phi*((h+eta)^gamma))
  B <- (1/alpha)*log(y-v-kappa)
  U <- A + B
  return(U)
}


