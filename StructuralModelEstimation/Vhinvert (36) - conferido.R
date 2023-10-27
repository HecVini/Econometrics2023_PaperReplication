Vhinvert <- function(v,con,muy,sigy,tau,omg,phi,gam,alpha,eta){
  #keyboard
  #v,con,
  #v(1,1),muy,sigy,tau,omg,phi,gam,alpha
  tht <- 0
  a <- muy-(sigy/tau)*omg
  A <- exp(a)
  b <- (sigy/tau)
  F <- 1-((v^(1-b))/A)
  Ffirst <- F(1,1)
  C <- exp((b-1)*con)
  Z <- (F/C)^(1/((b-1)*alpha))
  h <- ((1-Z)/phi)^(1/gam)-eta
  return(h)
}
