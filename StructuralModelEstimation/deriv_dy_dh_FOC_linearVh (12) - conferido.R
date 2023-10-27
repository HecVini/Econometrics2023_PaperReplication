# Derivative dy_dh from the FOCs, when pricing function is linear
deriv_dy_dh_FOC_linearVh <- function(h,p,alpha,phi,gammap,eta){
  S <- -p/(alpha*phi*gammap)
  M <- (h+eta)^(-gammap)
  L <- -alpha*phi*gammap -phi + (1 - gammap) * M
  return(S*L)
  
}