# Calculate the compensating variation required to keep utility that guarantees the same utility for a household 
# that moves from the reference m metro area to the base metro area.
# yphat is the resulting total income and cv is the compensating variation (the additional amount)
library(neldermead)
CalcCV <- function(y,uph,v,vm,hg){  

	opts = optimset(MaxFunEvals = 200, MaxIter = 200, TolFun = 1e-5, TolX = 1e-20)

	c(utest, htest, vtest)   <- caluy(y,uph,v,hg)

	c(um,hm,vhm)   <- caluy(y,uph,vm,hg)

	CV <- function(yp){

	  c(uh,hh,vh)   <- caluy(exp(yp),uph,v,hg) 

	  if (imag(um)==0 & imag(uh)==0){
	    F <- (um-uh).^2
	  } else {
	    F <- 1000000000000000
	  }
	  return(FF)
	}	

	ypini <- y
	[yp_hat ,fval] <- fminsearch(CV,log(ypini),opts)
	yphat <- exp(yp_hat)
	cv <- yphat - y

	phi <- uph$phi
	gamma <- uph$gamma
	alpha <- uph$alpha
	eta <- uph$eta
	kappa <- uph$kappa

[uhat,hhat,vhat] <- caluy(yphat,uph,v,hg)
utility(yphat,hh,phi,gamma,alpha,eta,kappa,vh)    
return(c(yphat, cv,fval))
}
