deriv_dy_dh_FOC <- function(h,alpha,phi,gammap,eta,vp,vpp){
  #     vp = Vprime_GLN4(h,mu,sigma,tau,omega,phi,gammap,alpha,eta);
  #     vpp = Vdoubleprime_GLN4(h,mu,sigma,tau,omega,phi,gammap,alpha,eta);
  heta <- h+eta
  TT <- vp*(1-phi*((heta)^gammap))
  TTp <- vpp*(1-phi*(heta^gammap)) + vp*(-phi*gammap*(heta^(gammap-1)))
  B <- -phi*gammap*alpha*(heta^(gammap-1))
  Bp <- -phi*gammap*alpha*(gammap-1)*(heta^(gammap-2))
  Dydh <- ((TTp*B - Bp*TT)/(B^2)) + vp
 
  return(Dydh) 
}
