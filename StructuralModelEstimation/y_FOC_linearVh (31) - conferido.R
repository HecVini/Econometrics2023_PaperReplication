#Income as a function of quality, obtained from the FOCs, when pricing function is linear
y_FOC_linearVh <- function(h,p,alpha,phi,gammap,eta,kappa){
  S <- p/(alpha*phi*gammap)
  A <- ((h+eta)^(-gammap))-phi
  M <- S*(h+eta)*A
  y <- kappa + p*h - M
  return(y)
}