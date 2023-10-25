# Função que determina os quantis da distribuição de cutoffs de renda (?)
# Variáveis: y (vetor de rendas), v(vetor de aluguéis) e up (lista nomeada dos valores dos parâmetros da função de utilidade)

Ycutoffinv<-function(y, up, v){
  hconsumed<-function(hvar){
    if((Ycutoff(exp(hvar),up,v)>0) & is.numeric(exp(hvar))){
      return (sum((Ycutoff(exp(hvar),up,v)-y[1:length(y)-1])**2))
    }
    else{
      return (100000000)
    }
  }
  
  hini <-mean(v[1:length(v)/2])
  hinip1<-mean(v[length(v)/2:length(v)])
  quant<-optim(log(c(hini, hinip1)),hconsumed)$par
  return(c('h_hat'=exp(quant[1]),'fval'=quant[2]))
  }

#Teste
y<-c(1000,2500,3000)
v<-c(50,100,130)
up<-c('eta'=1,'alpha'=0.7,'phi'=0.5,'gamma'=-0.2,'kappa'=0.8)
Ycutoffinv(y,up,v)

