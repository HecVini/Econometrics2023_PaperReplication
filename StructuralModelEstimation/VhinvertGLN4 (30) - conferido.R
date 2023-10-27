VhinvertGLN4 <- function(v,muy,sigma,tau,omega,phi,gammap,alpha,eta,theta){
  #We are changing the scale of theta everywhere
  theta <- 1000*theta
  a <- muy-(sigma/tau)*omega
  A <- exp(a)
  b <- (sigma/tau)
  vtheta <- v+theta
  F <- 1-((vtheta^(1-b))/A)
  S <- F^(1/(alpha*(b-1)))
  Z <- ((1-S)/phi)^(1/gammap)
  h <- Z-eta;
  
  return(h)
}
