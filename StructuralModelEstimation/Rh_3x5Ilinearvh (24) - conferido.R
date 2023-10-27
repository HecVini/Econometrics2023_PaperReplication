# Calculate demand and supplu for each observed and unobserved type in t=1 for the base metro area, imposing a
# linear pricing function. 
Rh_3x5Ilinearvh <- function(h,pxy,pxxy,pxxxy,upi,upii,upiii,upiv,upv,si,sii,siii,p1,p2,p3){
  
  #This function calculates the quality distribution when we assume that pricing function is linear in the first period. 
  
  #Linear pricing function
  vp <- h
  
  #Income cutoffs for each type
  ycfi <- Ycutoff(h,upi,vp)
  ycfii <- Ycutoff(h,upii,vp)        
  ycfiii <- Ycutoff(h,upiii,vp)        
  ycfiv <- Ycutoff(h,upiv,vp)
  ycfv <- Ycutoff(h,upv,vp)        
  
  #Conditions
  monoyi <- monotonic(Re(ycfi), direction = 'inc') 
  monoyii <- monotonic(Re(ycfii), direction = 'inc') 
  monoyiii <- monotonic(Re(ycfiii), direction = 'inc') 
  monoyiv <- monotonic(Re(ycfiv), direction = 'inc') 
  monoyv <- monotonic(Re(ycfv), direction = 'inc') 
  
  realyi <- !is.complex(ycfi)
  realyii <- !is.complex(ycfii)
  realyiii <- !is.complex(ycfiii)
  realyiv <- !is.complex(ycfiv)
  realyv <- !is.complex(ycfv)
  
  budgeti <- ycfi-vp[1:(length(vetor)-1)]-upi.kappa
  budgetii <- ycfii-vp[1:(length(vetor)-1)]-upii.kappa
  budgetiii <- ycfiii-vp[1:(length(vetor)-1)]-upiii.kappa
  budgetiv <- ycfiv-vp[1:(length(vetor)-1)]-upiv.kappa
  budgetv <- ycfv-vp[1:(length(vetor)-1)]-upv.kappa
  
  
  if (realyi & realyii & realyiii & realyiv & realyv & 
  monoyi & monoyii & monoyiii & monoyiv & monoyv & 
  ycfi>0 & ycfii>0 & ycfiii>0 & ycfiv>0 & ycfv>0 & 
  budgeti>0 & budgetii>0 & budgetiii>0 & budgetiv>0 & budgetv>0){
  #I don't need the income distr to get the income cutoffs
      #but I do need it to get F. The demand is the difference in the Fs
      pi <- c(p1[1], p2[1], p3[1])
      pii <- c(p1[2], p2[2], p3[2])
      piii <- c(p1[3], p2[3], p3[3])
      piv <- c(p1[4], p2[4], p3[4])
      pv <- c(p1[5], p2[5], p3[5])
      
      #demand
      hdi <- hd(pxy,pxxy,pxxxy,si,sii,siii,pi,ycfi)
      hdii <- hd(pxy,pxxy,pxxxy,si,sii,siii,pii,ycfii)
      hdiii <- hd(pxy,pxxy,pxxxy,si,sii,siii,piii,ycfiii)
      hdiv <- hd(pxy,pxxy,pxxxy,si,sii,siii,piv,ycfiv)
      hdv <- hd(pxy,pxxy,pxxxy,si,sii,siii,pv,ycfv)
      
      #demand CDFs
      Hdcdi[1] <- hdi[1] 
      Hdcdii[1] <- hdii[1] 
      Hdcdiii[1] <- hdiii[1] 
      Hdcdiv[1] <- hdiv[1] 
      Hdcdv[1] <- hdv[1] 
      
      for (j in 2:length(h)-1){
          Hdcdi[j] <- Hdcdi[j-1]+hdi[j]
          Hdcdii[j] <- Hdcdii[j-1]+hdii[j]
          Hdcdiii[j] <- Hdcdiii[j-1]+hdiii[j]
          Hdcdiv[j] <- Hdcdiv[j-1]+hdiv[j]
          Hdcdv[j] <- Hdcdv[j-1]+hdv[j]
      }
      
      
      #Aggregate the CDFs (get aggregate demand distribution)
      H <- c(Hdcdi, Hdcdii, Hdcdiii, Hdcdiv, Hdcdv)
      s <- c(si, sii, siii)
      
      ax1 <- (tr(s)*p1)
      bx1 <- ax1[1,]
      
      ax2 <- (tr(s)*p2)
      bx2 <- ax2[2,]
      
      ax3 <- (tr(s)*p3)
      bx3 <- ax3[3,]
      
      Hx <- bx1*H
      Hxx <- bx2*H
      Hxxx <- bx3*H
      
      AgHdcdf <- Hx + Hxx + Hxxx
  
  #also return the observed type demands to construct moments in the main
  #objective function
  
      cond=1
    }else{
      AgHdcdf <- vector()
      Hdcdi <- vector()
      Hdcdii <- vector()
      Hdcdiii <- vector()
      Hx  <- vector() 
      Hxx <- vector() 
      Hxxx <- vector()
      cond=0
    }
  Rh <- AgHdcdf
  return(c(Rh,Hdcdi,Hdcdii,Hdcdiii,Hdcdiv,Hdcdv, Hx, Hxx, Hxxx,AgHdcdf,cond))  
    
}

