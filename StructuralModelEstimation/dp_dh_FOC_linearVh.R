#Função que determina a FOC quando se assume que a quação de precificação é linear
#Variáveis: vetor de qualidades (h),vetor de rendas (y), parâmetros da função de utilidade (alpha,phi,gammap,eta,kappa)

dp_dh_FOC_linearVh<-function(h,eta,alpha,phi,gammap,y,kappa){
  heta = (h + eta)
  S = -alpha*phi*gammap*h + (heta)**(1-gammap) - phi*heta
  Sp = -alpha*phi*gammap + (1-gammap)*(heta)**(-gammap) - phi
  M = -alpha*phi*gammap*(y-kappa)
  return (-M*Sp/(S**2))
}

#Teste
h<-c(60,110,150)
y<-c(1000,2500,3000)
up<-c('eta'=1,'alpha'=0.7,'phi'=0.5,'gamma'=-0.2,'kappa'=0.8)
print(dp_dh_FOC_linearVh(h,up['eta'], up['alpha'],up['phi'],up['gamma'],y,up['kappa']))
