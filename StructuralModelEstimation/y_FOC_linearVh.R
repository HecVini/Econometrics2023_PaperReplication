#Função que determina renda em função da qualidade com base nas FOCs quando a função de precificação é linear
#Variáveis: vetor de qualidades (h), vetor de preços (p), parâmetros da função de utilidade (alpha,phi,gammap,eta,kappa)

y_FOC_linearVh<-function(h,p,alpha,phi,gammap,eta,kappa){
  S = p/(alpha*phi*gammap)
  A = ((h+eta)**(-gammap))-phi
  M = S*(h+eta)*A;
  return (kappa + p*h - M)
}

#Teste
h<-c(60,110,150)
p<-c(100,250,300)
up<-c('eta'=1,'alpha'=0.7,'phi'=0.5,'gamma'=-0.2,'kappa'=0.8)
print(p_FOC_linearVh(h,ym, up['alpha'],up['phi'],up['gamma'],up['eta'],up['kappa']))