#Income as a function of quality, obtained from the FOCs
y_FOC_Vh <- function(h,alpha,phi,gammap,eta,kappa,v,vp){
  heta <- (h+eta)
  S <- 1 - phi*(heta^gammap)
  L <- -phi*gammap*alpha*(heta^(gammap-1))
  M <- (vp*(S))/(L)
  y <- kappa + v + M
  return(y)
}
