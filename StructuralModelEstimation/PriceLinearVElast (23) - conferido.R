#Calculate price quality consumed elasticity for given type when a linear pricing function is
#assumed
PriceLinearVElast <- function(h,U){
  
  #Recover Parameters
  phi <- U$phi
  gammap <- U$gamma
  alpha <- U$alph
  eta <- U$eta
  kappa  <-  U$kappa
#Elasticities
#Corresponding income
  v  <-  h
  y  <-  Ycutoff(h,U,v)
# Period 1
  Dp1h  <-  dp_dh_FOC_linearVh(h[1:(length(vetor)-0)],eta,alpha,phi,gammap,y,kappa)
  Dhp1  <-  1/Dp1h               
  p1  <-  p_FOC_linearVh(h[1:(length(vetor)-0)],y,alpha,phi,gammap,eta,kappa)
  ehp1  <-  Dhp1*(p1/h[1:(length(vetor)-0)])
  return(ehp1)
}
