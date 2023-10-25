#Função que calcula os cutoffs na renda que tornam o agente indiferente entre qualidades j e j+1
#Variáveis: h (vetor de qualidades), v(vetor de aluguéis) e up (lista nomeada dos valores dos parâmetros da função de utilidade)

Ycutoff<-function(h,up,v){
  s=log(1-up['phi']*((h+up['eta'])**up['gamma']))
  E<-c()
  L<-c()
  M<-c()
  ycf<-c()
  for (i in 1:(length(h)-1)){
    E<-c(E, exp((s[i+1]-s[i])*up['alpha']))
    L<-c(L, 1-E[i])
    M<-c(M, v[i]-E[i]*v[i+1])
    ycf<-c(ycf, (M[i]/L[i])+up['kappa'])
  }
  return(unname(ycf))
}

#Teste
h<-c(60,110,150)
v<-c(50,100,130)
up<-c('eta'=1,'alpha'=0.7,'phi'=0.5,'gamma'=-0.2,'kappa'=0.8)
print(Ycutoff(h,up,v))
