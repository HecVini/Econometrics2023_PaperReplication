#library pracma e signal
IncomeElast <- function(h, v, U){
  phi <- U$phi
  gammap <- U$gamma
  alpha <- U$alpha
  eta <- U$eta
  kappa  <-  U$kappa

install.packages("pracma", "signal")
library("pracma")

polytool(h,v,d)
d  <-  2
p  <-  polyfit(h, v, d) #obtém os coeficientes do polinomio - é polyfit no R tbm

v  <-  polyval(p,h); #calcula o polinomio para todos os valores de h - tbm é polyval

#v'(h)
pd  <-  polyder(p) #obtém os coeficientes da derivada do polinomio - tbm é polyder
vp  <-  polyval(pd,h)  #calcula a derivada
#v''(h) - mesma coisa para a segunda derivada.
pdd  <-  polyder(pd) 
vpp  <-  polyval(pdd,h)



Dy1h  <-  deriv_dy_dh_FOC(h,alpha,phi,gammap,eta,vp,vpp)
Dhy1  <-  1/Dy1h
y1  <-  y_FOC_Vh(h,alpha,phi,gammap,eta,kappa,v,vp)
ehy  <-  Dhy1*(y1/h)
return(ehy)
}
