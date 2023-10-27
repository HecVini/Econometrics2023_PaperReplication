#Calculate Supplies and Demand for each quality level for the 2 metro area
# for the second or reference metro area
Hinv3x5INorm1m <- function(h,vp2m, pxym, pxxym, pxxxym, pxytm, pxxytm, pxxxytm, upi,upii,upiii,upiv,upv,sxm,sxxm,sxxxm,p1,p2,p3,zeta,VIm,vp1m,VI,Rh){
  
  # Calculate supply and demand 
  
  ## Check monotonicity
  monov1 <- monotonic(vp1m, direction = 'inc')
  
  if  (monov1){
    ## Demand for period 1
    c(Hdcdim,Hdcdiim,Hdcdiiim,Hdcdivm,Hdcdvm,Hxm, Hxxm, Hxxxm,AgHdcdfm,cond1m) <- Rh_3x5Ivh(h,pxym,pxxym,pxxxym,upi,upii,upiii,upiv,upv,sxm,sxxm,sxxxm, p1, p2, p3, vp1m)
  
    if (cond1m ==1){
  
      ## Check monotonicity
      monov2 <- monotonic(vp2m, direction = 'inc')
    
      if  (monov2){
        ## Demand for period 2
    
        c(Hdtcdim,Hdtcdiim,Hdtcdiiim,Hdtcdivm,Hdtcdvm,Htxm, Htxxm, Htxxxm,AgHdtcdfm,cond2m) <- Rh_3x5Ivh(h,pxytm,pxxytm,pxxxytm,upi,upii,upiii,upiv,upv,sxm,sxxm,sxxxm,p1, p2, p3, vp2m)
    
        if (cond2m == 1){
      
            ## Supply for period 1
            ## of periods between observations
            lengthper <- 4
            #Base metro
            Vbase1$k1ypolyf <- VI$uf
            Vbase1$V <- VI$V 
      
            #Ref Metro
            Vref1$k2ypolyf <- VIm$uf
            Vref1$Vt <- VIm$V
            
            # linear function for price of first period in the
            # refence metro
            vbase <- t(h)
            vref <- vp2m                    
            basesupply <- Rh
            
            c(Rh1m, N1m, gratehm) <- RhSl2PeriodsIm(h,vbase,vref,Vref1,zeta,basesupply,lengthper,Vbase1)    
            
            ## Supply for period 2
             #Base metro
            Vbase2$k1ypolyf <- VIm$uf
            Vbase2$V <- VIm$V 
                        
            #Ref Metro
            Vref2$k2ypolyf <- VIm$uft 
            Vref2$Vt <- VIm$Vt
            
            vbase <- vp1m
            vref <- vp2m                    
            basesupply <- Rh1m
            
            c(Rh2m, N2m, gratehtm) <- RhSl2PeriodsIm(h,vbase,vref,Vref2,zeta,basesupply,lengthper,Vbase2)
          
        } else {
            Rh1m <- vector() 
            N1m <- vector() 
            gratehm <- vector() 
            Rh2m <- vector() 
            N2m <- vector() 
            gratehtm <- vector()  
        }
      } else {
    
          Hdtcdim <- vector() 
          Hdtcdiim <- vector() 
          Hdtcdiiim <- vector() 
          Hdtcdivm <- vector() 
          Hdtcdvm <- vector()
          Htxm  <- vector() 
          Htxxm <- vector() 
          Htxxxm <- vector()
          AgHdtcdfm <- vector() 
          
          Rh1m <- vector() 
          N1m <- vector() 
          gratehm <- vector()  
          Rh2m <- vector() 
          N2m <- vector() 
          gratehtm <- vector()  
          
          cond2m <- 0
      }
       
    } else {
  
      Hdtcdim <- vector()
      Hdtcdiim <- vector() 
      Hdtcdiiim <- vector() 
      Hdtcdivm <- vector() 
      Hdtcdvm <- vector()
      Htxm  <- vector() 
      Htxxm <- vector() 
      Htxxxm <- vector()
      AgHdtcdfm <- vector() 
      Rh1m <- vector() 
      N1m <- vector() 
      gratehm <- vector()  
      Rh2m <- vector() 
      N2m <- vector() 
      gratehtm <- vector()  
      cond2m <- 0
    }
          
  } else { 
          
  Hdcdim <- vector() 
  Hdcdiim <- vector() 
  Hdcdiiim <- vector() 
  Hdcdivm <- vector() 
  Hdcdvm <- vector()
  Hdtcdim <- vector()
  Hdtcdiim <- vector() 
  Hdtcdiiim <- vector() 
  Hdtcdivm <- vector() 
  Hdtcdvm <- vector()
  
  Hxm  <- vector() 
  Hxxm <- vector() 
  Hxxxm <- vector()
  Htxm  <- vector() 
  Htxxm <- vector() 
  Htxxxm <- vector()
          
  AgHdcdfm <- vector() 
  AgHdtcdfm <- vector() 
  
  Rh1m <- vector() 
  N1m <- vector() 
  gratehm <- vector()  
  Rh2m <- vector() 
  N2m <- vector() 
  gratehtm <- vector()  
  cond2m <- 0
  cond1m <- 0
              
  }
  return(c(Rh1m,cond1m,cond2m,Hdcdim,Hdcdiim,Hdcdiiim,Hdcdivm,Hdcdvm,Hxm,Hxxm,Hxxxm,AgHdcdfm,Hdtcdim,Hdtcdiim,Hdtcdiiim,Hdtcdivm,Hdtcdvm,Htxm,Htxxm,Htxxxm,AgHdtcdfm,Rh2m))
}
