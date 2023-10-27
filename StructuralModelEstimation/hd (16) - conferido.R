# Calculate the 
hd <-function(pxy,pxxy,pxxxy,si,sii,siii,pi,ycf){
#Calculate the demand for each quality j = 1:J implied by the income
#cutoffs ycf. ycf gives the cutoff points that make households 
#indifferent between quality j and j-1 given qualities J and prices v

  mux <- pxy(1)
  sigmax <- pxy(2) 
  rx <- pxy(3) 
  betax <- pxy(4)
  muxx <- pxxy(1) 
  sigmaxx <- pxxy(2) 
  rxx <- pxxy(3) 
  betaxx <- pxxy(4)
  muxxx <- pxxxy(1) 
  sigmaxxx <- pxxxy(2) 
  rxxx <- pxxxy(3) 
  betaxxx <- pxxxy(4)

  #Get the income distribution for each observed type evaluated at the income
  #cutoffs
  GLNyxj <- cdfGLN4(ycf[2:(length(ycf)-0)],exp(mux),mux,sigmax,rx,betax)
  GLNyxxj <- cdfGLN4(ycf[2:(length(ycf)-0)],exp(muxx),muxx,sigmaxx,rxx,betaxx)
  GLNyxxxj <- cdfGLN4(ycf[2:(length(ycf)-0)],exp(muxxx),muxxx,sigmaxxx,rxxx,betaxxx)
  
  
  GLNyxj_1 <- cdfGLN4(ycf[1:length(ycf)-1],exp(mux),mux,sigmax,rx,betax)
  GLNyxxj_1 <- cdfGLN4(ycf[1:length(ycf)-1],exp(muxx),muxx,sigmaxx,rxx,betaxx)
  GLNyxxxj_1 <- cdfGLN4(ycf[1:length(ycf)-1],exp(muxxx),muxxx,sigmaxxx,rxxx,betaxxx)


  #Get the income distribution for each unobserved type evaluated at the income
  #cutoffs
  # contemporaneous
  s <- c(si, sii, siii)
  FF <- c(t(GLNyxj), t(GLNyxxj), t(GLNyxxxj))
  D <- s*t(pi)
  N <- diag(t(t(s)*pi))*FF 
  Fi <-  N/D
                
  # lagged
  F_1 <- c(t(GLNyxj_1), t(GLNyxxj_1), t(GLNyxxxj_1))
  N_1 <- diag(t(t(s)*pi))*F_1
  Fi_1 <-  N_1/D
      
  # demands are the difference between the two cumulative income
  # distributions
  hd(1) <-  Fi(1)
  hd(2:length(ycf)) <-  Fi - Fi_1
  return(hd)      
}

