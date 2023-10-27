#Calculate Supplies and Demand for each quality level for the 1 metro area
#model or for the base metro area in the multiple metro areas model
Hinv3x5INorm1 <- function(h,vp2,pxy,pxxy,pxxxy, pxyt, pxxyt, pxxxyt, upi,upii,upiii,upiv,upv,sx,sxx,sxxx,p1,p2,p3,zeta,VI){

  # Rh from the normalization that sets v(h) = v
  ## Calculate supply and demand in t=1
  c(Rh,Hdcdi,Hdcdii,Hdcdiii,Hdcdiv,Hdcdv,Hx,Hxx,Hxxx,AgHdcdf,cond1) <-
  Rh_3x5Ilinearvh(h,pxy,pxxy,pxxxy,upi,upii,upiii,upiv,upv,sx,sxx,sxxx, p1, p2, p3)                     
  
  if (cond1 ==1){
    ## Check monotonicity
    monov <- monotonic(vp2,direction = 'inc')
  
    if  (monov) {
      ## Demand for period 2
      c(Hdtcdi,Hdtcdii,Hdtcdiii,Hdtcdiv,Hdtcdv,Htx, Htxx, Htxxx,AgHdtcdf,cond2) <- Rh_3x5Ivh(h,pxyt,pxxyt,pxxxyt,upi,upii,upiii,upiv,upv,sx,sxx,sxxx, p1, p2, p3, vp2)
  
      if (cond2 ==1){
        ## Supply for period 2
        
        ## of periods between observations
        lengthper <- 4
  
        VI$k1ypolyf <- VI$uf
        VI$k2ypolyf <- VI$uft 
  
        # linear function for price of first period
        vp1 <- t(h)
        c(Rh2, N2, grateh) <-RhSl2PeriodsI(h,vp1,vp2,VI,zeta,AgHdcdf,lengthper)
      } else {
        Rh2= vector()
        N2= vector() 
        grateh = vector()  
      }      
    } else { 
  
          Hdtcdi = vector()
          Hdtcdii = vector() 
          Hdtcdiii = vector() 
          Hdtcdiv = vector() 
          Hdtcdv = vector()
          Htx  = vector() 
          Htxx = vector() 
          Htxxx = vector()
          AgHdtcdf = vector() 
          Rh2= vector() 
          N2 = vector() 
          grateh = vector()  
          cond2=0
    }
    
  } else {
      
      Hdtcdi = vector()
      Hdtcdii = vector() 
      Hdtcdiii = vector() 
      Hdtcdiv = vector() 
      Hdtcdv = vector()
      Htx  = vector() 
      Htxx = vector() 
      Htxxx = vector()
      AgHdtcdf = vector() 
      Rh2= vector() 
      N2 = vector() 
      grateh = vector()  
      cond2=0
  }

  return(c(Rh,cond1,cond2,Hdcdi,Hdcdii,Hdcdiii,Hdcdiv,Hdcdv,Hx, Hxx, Hxxx,AgHdcdf,Hdtcdi,Hdtcdii,Hdtcdiii,Hdtcdiv,Hdtcdv,Htx, Htxx, Htxxx,AgHdtcdf, Rh2))
}
