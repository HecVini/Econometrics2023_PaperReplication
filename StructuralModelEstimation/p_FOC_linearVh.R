#Função que determina a FOC quando se assume que a equação de precificação é linear
#Variáveis: vetor de qualidades (h), vetor de rendas (ym), parâmetros da função de utilidade (alpha,phi,gammap,eta,kappa)

p_FOC_linearVh<-function(h,ym,alpha,phi,gammap,eta,kappa){
  heta = (h + eta)
  S = -alpha*phi*gammap*h + (heta)**(1-gammap) - phi*heta
  M = -alpha*phi*gammap*(ym-kappa)
  return (M/S)
}

#Teste
h<-c(60,110,150)
ym<-c(1000,2500,3000)
up<-c('eta'=1,'alpha'=0.7,'phi'=0.5,'gamma'=-0.2,'kappa'=0.8)
print(p_FOC_linearVh(h,ym, up['alpha'],up['phi'],up['gamma'],up['eta'],up['kappa']))
      
