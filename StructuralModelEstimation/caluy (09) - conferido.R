# Calculate the quality h consumed at income y at prices v in quality grid hg
caluy <- function(y,uph,v,hg){

# What is the quality that a household with income y would consume?
  yh <- Ycutoff(hg,uph,v)
  
  diff_squared <- (yh - y)^2
  b <- which.min(diff_squared)
  a <- min(diff_squared)
  h <- hg(b)
  vh <- v(b)
  
  #What is the utility of this quality at this rent price and this income
  phi <- uph.phi
  gamma <- uph.gamma
  alpha <- uph.alpha
  eta <- uph.eta
  kappa <- uph.kappa
  
  u <- utility(y,h,phi,gamma,alpha,eta,kappa,vh);
  return(c(u,h,vh)) 
}
  