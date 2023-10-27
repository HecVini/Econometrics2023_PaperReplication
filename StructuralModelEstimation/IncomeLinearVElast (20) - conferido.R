#library pracma e signal
IncomeLinearVElast <- function(h,U){
  
  install.packages("pracma", "signal")
  library("pracma")
  
  phi <- U$phi
  gammap <- U$gamma
  alpha <- U$alpha
  eta <- U$eta
  kappa <- U$kappa
  
  d  <-  2
  p  <-  polyfit(h, v, d) #obtém os coeficientes do polinomio - é polyfit no R tbm
  
  v  <-  polyval(p,h) #calcula o polinomio para todos os valores de h - tbm é polyval
  
  #v'(h)
  pd  <-  polyder(p) #obtém os coeficientes da derivada do polinomio - tbm é polyder
  vp  <-  polyval(pd,h)  #calcula a derivada
  #v''(h) - mesma coisa para a segunda derivada.
  pdd  <-  polyder(pd) 
  vpp  <-  polyval(pdd,h)
  
  
    
  Dy1h  <-  deriv_dy_dh_FOC_linearVh(h,p,alpha,phi,gammap,eta)
  Dhy1  <-  1/Dy1h
  y1  <-  y_FOC_linearVh(h,p,alpha,phi,gammap,eta,kappa)
  ehy  <-  Dhy1*(y1/h)
  return(ehy)
}
