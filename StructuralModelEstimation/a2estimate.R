#Função que estima os parâmetros do modelo usando ambas as regiões metropolitanas e suas correlações


library(dplyr)

a2Estimate<-function(YI,VI,PI,Pop, Corr){
  sx<-Pop$I$k1
  sxt<-Pop$I$k1
  sxx<-Pop$I$k2
  sxxt<-Pop$I$k2
  sxxx<-Pop$I$k3
  sxxxt<-Pop$I$k3
  
  typex<-'Type 1'
  typexx<-'Type 2'
  typexxx<-'Type 3'
  
  inicut<-2
  begy<-2
  beg<-3
  hg<-VI$V(beg:length(VI$V)-inicut)
  Ycut<-YI$Y(begy:length(VI$Y)-inicut)
  Ytcut<-YI$Yt(begy:length(VI$Yt)-inicut)
  
  px<-list()
  pxx<-list()
  pxxx<-list()
  px$y<-PI$GLN4ui$y
  pxx$y<-PI$GLN4uii$y
  pxxx$y<-PI$GLN4uiii$y
  px$yt<-PI$GLN4ui$yt
  pxx$yt<-PI$GLN4uii$yt
  pxxx$yt<-PI$GLN4uiii$yt
  
  
  mu1x<-px$y[1]
  sigma1x<-px$y[2]
  r1x<-px$y[3]
  beta1x<-px$y[4]
  F1xob <- cdfGLN4(YI$Y[begy:length(YI$Y)-inicut],exp(mu1x),mu1x,sigma1x,r1x,beta1x)
  
  mu2x<-px$yt[1] 
  sigma2x<-px$yt[2] 
  r2x<-px$yt[3]  
  beta2x<-px$yt[4]
  F2xob <- cdfGLN4(YI$Yt[begy:length(YI$Yt)-inicut],exp(mu2x),mu2x,sigma2x,r2x,beta2x)
  
  mu1xx<-pxx$y[1]  
  sigma1xx<-pxx$y[2]   
  r1xx<-pxx$y[3] 
  beta1xx<-pxx$y[4]
  F1xxob <- cdfGLN4(YI$Y[begy:length(YI$Y)-inicut],exp(mu1xx),mu1xx,sigma1xx,r1xx,beta1xx)
  
  mu2xx<-pxx$yt[1]
  sigma2xx<-pxx$yt[2] 
  r2xx<-pxx$yt[3]  
  beta2xx<-pxx$yt[4]
  F2xxob <- cdfGLN4(YI$Yt[begy:length(YI$Yt)-inicut],exp(mu2xx),mu2xx,sigma2xx,r2xx,beta2xx)
  
  mu1xxx<-pxxx$y[1] 
  sigma1xxx<-pxxx$y[2]  
  r1xxx<-pxxx$y[3]  
  beta1xxx<-pxxx$y[4]
  F1xxxob <- cdfGLN4(YI$Y[begy:length(YI$Y)-inicut],exp(mu1xxx),mu1xxx,sigma1xxx,r1xxx,beta1xxx)
  
  mu2xxx<-pxxx$yt[1] 
  sigma2xxx<-pxxx$yt[2] 
  r2xxx<-pxxx$yt[3]
  beta2xxx<-pxxx$yt[4]
  F2xxxob <- cdfGLN4(YI$Yt[begy:length(YI$Yt)-inicut],exp(mu2xxx),mu2xxx,sigma2xxx,r2xxx,beta2xxx)
  
  
  px$v<-PI$GLN4ui$v
  pxx$v<-PI$GLN4uii$v
  pxxx$v<-PI$GLN4uiii$v
  px$vt<-PI$GLN4ui$vt
  pxx$vt<-PI$GLN4uii$vt
  pxxx$vt<-PI$GLN4uiii$vt
  
  omega1x<-px$v[1]
  tau1x<-px$v[2]
  m1x<-px$v[3]
  theta1x<-px$v[4]
  G1xob <- cdfGLN4(VI$V[begy:length(VI$V)-inicut],exp(omega1x),omega1x,tau1x,m1x,theta1x)
  
  omega2x<-px$vt[1] 
  tau2x<-px$vt[2] 
  m2x<-px$vt[3]  
  theta2x<-px$vt[4]
  G2xob <- cdfGLN4(VI$Vt[begy:length(VI$Vt)-inicut],exp(omega2x),omega2x,tau2x,m2x,theta2x)
  
  omega1xx<-pxx$v[1]  
  tau1xx<-pxx$v[2]   
  m1xx<-pxx$v[3] 
  theta1xx<-pxx$v[4]
  G1xxob <- cdfGLN4(VI$V[begy:length(VI$V)-inicut],exp(omega1xx),omega1xx,tau1xx,m1xx,theta1xx)
  
  omega2xx<-pxx$vt[1]
  tau2xx<-pxx$vt[2] 
  m2xx<-pxx$vt[3]  
  theta2xx<-pxx$vt[4]
  G2xxob <- cdfGLN4(VI$Vt[begy:length(VI$Vt)-inicut],exp(omega2xx),omega2xx,tau2xx,m2xx,theta2xx)
  
  omega1xxx<-pxxx$v[1] 
  tau1xxx<-pxxx$v[2]  
  m1xxx<-pxxx$v[3]  
  theta1xxx<-pxxx$v[4]
  G1xxxob <- cdfGLN4(VI$V[begy:length(VI$V)-inicut],exp(omega1xxx),omega1xxx,tau1xxx,m1xxx,theta1xxx)
  
  omega2xxx<-pxxx$vt[1] 
  tau2xxx<-pxxx$vt[2] 
  m2xxx<-pxxx$vt[3]
  theta2xxx<-pxxx$vt[4]
  G2xxxob <- cdfGLN4(VI$Vt[begy:length(VI$Vt)-inicut],exp(omega2xxx),omega2xxx,tau2xxx,m2xxx,theta2xxx)
  
  pr<-matrix(nrow=3,ncol=5)
  pr[1,1:4] <- c(0.4, 0.05, 0.03, 0.3)
  pr[1,5] <- 1-sum(pr[1,1:4])
  pr[2,1:4] <- c(0.1, 0.4, 0.01, 0.05)
  pr[2,5] <- 1-sum(pr[2,1:4])
  pr[3,1:4] <- c(0.05, 0.2, 0.5, 0.1)
  pr[3,5] <- 1-sum(pr[3,1:4])
  
  test_prob_sumx1 <- sum(pr[1,])
  test_prob_sumx2 <- sum(pr[2,])
  test_prob_sumx3 <- sum(pr[3,])
  
  pr0<-log(pr)
  
  p0i<- PI$GLN4ui
  u1<-CalLinearPref(p0i$y,p0i$v)
  u1$kappa<-  0
  u2<- u1
  u3<-u1
  u4<- u1
  u5 <- u1
  
  u2$alpha <- u2$alpha+runif()/10
  u3$alpha <- u3$alpha+runif()/10
  u4$alpha <- u4$alpha+runif()/10
  u5$alpha <- u5$alpha+runif()/10
  
  u2$gamma <- u2$gamma-runif()/10
  u3$gamma <- u3$gamma-runif()/10
  u4$gamma <- u4$gamma-runif()/10
  u5$gamma <- u5$gamma-runif()/10
  
  
  v01 <- log(hg)
  v02 <- log(hg)
  

  testyi <- Ycutoff(hg,u1,exp(v01))
  testyii <- Ycutoff(hg,u2,exp(v01))
  testyiii <- Ycutoff(hg,u3,exp(v01))
  testyiv <- Ycutoff(hg,u4,exp(v01))
  testyiv <- Ycutoff(hg,u5,exp(v01))
  
  testyit <- Ycutoff(hg,u1,exp(v02))
  testyiit <- Ycutoff(hg,u2,exp(v02))
  testyiiit <- Ycutoff(hg,u3,exp(v02))
  testyivt <- Ycutoff(hg,u4,exp(v02))
  testyvt <- Ycutoff(hg,u5,exp(v02))
  
  utilconpositive(u1,hg)
  utilconpositive(u2,hg)
  utilconpositive(u3,hg)
  utilconpositive(u4,hg)
  utilconpositive(u5,hg)
  

  u01 <- c(u1$alpha, (u1$phi)/100, (u1$eta)/1000, u1$gamma*10, (u1$kappa)/1000)
  u02 <- c(u2$alpha, (u2$phi)/100, (u2$eta)/1000, u2$gamma*10, (u2$kappa)/1000)
  u03 <- c(u3$alpha, (u3$phi)/100, (u3$eta)/1000, u3$gamma*10, (u3$kappa)/1000)
  u04 <- c(u4$alpha, (u4$phi)/100, (u4$eta)/1000, u4$gamma*10, (u4$kappa)/1000)
  u05 <- c(u5$alpha, (u5$phi)/100, (u5$eta)/1000, u5$gamma*10, (u5$kappa)/1000)
  

  zeta0 <- log(0.065)
  
  Q<-function(p_hat){
    upi$alpha <- p_hat[1] 
    upi$phi <- p_hat[2]*100 
    upi$eta <- p_hat[3]*1000
    upi$gamma <- p_hat[4]/10 
    upi$kappa <- p_hat[5]*1000    
    
    upii$alpha <- p_hat[6]
    upii$phi <- p_hat[7]*100
    upii$eta <- p_hat[8]*1000
    upii$gamma <- p_hat[9]/10
    upii$kappa <- p_hat[10]*1000    
    
    upiii$alpha <- p_hat[11]
    upiii$phi <- p_hat[12]*100 
    upiii$eta <- p_hat[13]*1000
    upiii$gamma <- p_hat[14]/10
    upiii$kappa <- p_hat[15]*1000    
    
    upiv$alpha <- p_hat[16]
    upiv$phi <- p_hat[17]*100 
    upiv$eta <- p_hat[18]*1000
    upiv$gamma <- p_hat[19]/10
    upiv$kappa <- p_hat[20]*1000    
    
    upv$alpha <- p_hat[21] 
    upv$phi <- p_hat[22]*100 
    upv$eta <- p_hat[23]*1000
    upv$gamma <- p_hat[24]/10 
    upv$kappa <- p_hat[25]*1000    
    
    vp1 <- hg
    vp2 <- exp(p_hat[26:39])
    
    p1 <- exp(p_hat[40:43])
    p1 <- append(p1,1-sum(p1))
    
    p2 <- exp(p_hat[44:47])
    p2 <- append(p2,1-sum(p2))
    
    p3 <- exp(p_hat[48:51])
    p3 <- append(p3,1-sum(p3))
    
    zeta <- exp(p_hat[52])
    
    
    cons1 <- utilconpositive(upi,hg)    
    cons2 <- utilconpositive(upii,hg)
    cons3 <- utilconpositive(upiii,hg)
    cons4 <- utilconpositive(upiv,hg)
    cons5 <- utilconpositive(upv,hg)
    
    
    if (cons1>0 & cons2>0 & cons3>0 & cons4>0 & cons5>0 & upi$alpha>0 & upii$alpha>0 & 
        upiii$alpha>0  & upiv$alpha>0 & upv$alpha>0 & upi$gamma<0 & upii$gamma<0 & upiii$gamma<0 & 
        upiv$gamma<0 & upv$gamma<0 & upi$phi>0 & upii$phi>0 & upiii$phi>0 & upiv$phi>0 & upv$phi>0 & 
        upi$eta>0 & upii$eta>0 & upiii$eta>0 & upiii$eta>0 & upiii$eta>0 & sum(p1)>0 & sum(p2)>0 & 
        sum(p3)>0 & all(p1 > 0) & all(p2 > 0) & all(p3 > 0)){
       Hinv<-Hinv3x5INorm1(hg,vp2,px$y,pxx$y,pxxx$y,px$yt,pxx$yt,pxxx$yt,upi,upii,upiii,upiv,upv,sx,sxx,sxxx,p1,p2,p3,zeta,VI)
       for (i in names(Hinv)){
         assign(i,unlist(Hinv[i]))
       }
       if(cond1 & cond2){
         eq <- 10*(AgHdt-Rh2)
         eqs <- sum(eq**2)
         eqslow <- sum(eq[1:7]**2)
         eqsup <- sum(eq[8:length(eq)]**2)
         
         dif1x <- 10*(cdfGLN4(vp1[1:length(vp1)-1],exp(omega1x),omega1x,tau1x,m1x,theta1x)-t(Hx))
         verr1x <- sum(dif1x**2)
         verr1xlow <- sum(dif1x[1:7]**2)
         verr1xup <- sum(dif1x[8:length(dif1x)]**2)
         
         dif1xx <- 10*(cdfGLN4(vp1[1:length(vp1)-1],exp(omega1xx),omega1xx,tau1xx,m1xx,theta1xx)-t(Hxx))
         verr1xx <- sum(dif1xx**2)                
         verr1xxlow <- sum(dif1xx[1:7]**2)
         verr1xxup <- sum(dif1xx[8:length(dif1xx)]**2)
         
         dif1xxx <- 10*(cdfGLN4(vp1[1:length(vp1)-1],exp(omega1xxx),omega1xxx,tau1xxx,m1xxx,theta1xxx)-t(Hxxx))
         verr1xxx <- sum(dif1xxx**2)                
         verr1xxxlow <- sum(dif1xxx[1:7]**2)
         verr1xxxup <- sum(dif1xxx[8:length(dif1xxx)]**2)
         
         dif2x <- 10*(cdfGLN4(vp2[1:length(vp2)-1],exp(omega2x),omega2x,tau2x,m2x,theta2x)-t(Htx))
         verr2x <- sum(dif2x**2)
         verr2xlow <- sum(dif2x[1:7]**2)
         verr2xup <- sum(dif2x[8:length(dif2x)]**2)
         
         dif2xx <- 10*(cdfGLN4(vp2[1:length(vp2)-1],exp(omega2xx),omega2xx,tau2xx,m2xx,theta2xx)-t(Htxx))
         verr2xx <- sum(dif2xx**2)                
         verr2xxlow <- sum(dif2xx[1:7]**2)
         verr2xxup <- sum(dif2xx[8:length(dif2xx)]**2)
         
         dif2xxx <- 10*(cdfGLN4(vp2[1:length(vp2)-1],exp(omega2xxx),omega2xxx,tau2xxx,m2xxx,theta2xxx)-t(Htxxx))
         verr2xxx <- sum(dif2xxx**2)                
         verr2xxxlow <- sum(dif2xxx[1:7]**2)
         verr2xxxup <- sum(dif2xxx[8:length(dif2xxx)]**2)
         
         CorrMoment <- Yvcorrelations(sx, sxx, sxxx, p1, p2, p3, hg, upi,upii,upiii,upiv,upv, vp1, vp2, Corr)
         F <- 10*((eqslow + eqsup) + (verr1xlow + verr1xup +verr1xxlow + verr1xxup+verr1xxxlow + verr1xxxup +
                                        verr2xlow + verr2xup +verr2xxlow + verr2xxup +verr2xxxlow + verr2xxxup+ CorrMoment)) 
         }else{
      F <- 1000000000000000000000000
      }}else{
        F <- 1000000000000000000000000
      }
  }       
  
  opts<-list(abstol=1e-4,reltol=1e-4,maxit=5000)
  
  iter <-50
  for (i in 1:iter){
    hat <- optim(c(u01, u02, u03, u04, u05, v02, pr0[1,1:length(pr0)-1], pr0[2,1:length(pr0)-1], pr0[3,1:length(pr0)-1], zeta0), Q,control=opts)
    phat<-hat$par
    Ehat<-hat$value
    u01 <- c(phat[1], phat[2], phat[3], phat[4], phat[5])
    u02 <- c(phat[6], phat[7], phat[8], phat[9], phat[10])
    u03 <- c(phat[11], phat[12], phat[13], phat[14], phat[15])
    u04 <- c(phat[16], phat[17], phat[18], phat[19], phat[20])
    u05 <- c(phat[21], phat[22], phat[23], phat[24], phat[25])
    v02 <- phat[26:39]
    
    pr0[1,1:4] <- phat[40:43]
    pr0_nolog[1,1:4] <- exp(phat[40:43])
    pr0_nolog[1,5] <- 1 - sum(pr0_nolog[1,1:4])
    pr0[1,5] <- log(pr0_nolog[1,5])
    
    pr0[2,1:4] <- phat[44:47]
    pr0_nolog[2,1:4] <- exp(phat[44:47])
    pr0_nolog[2,5] <- 1 - sum(pr0_nolog[2,1:4])
    pr0[2,5] <- log(pr0_nolog[2,5])
    
    pr0[3,1:4] <- phat[48:51]
    pr0_nolog[3,1:4] <- exp(phat[48:51])
    pr0_nolog[3,5] <- 1 - sum(pr0_nolog[3,1:4])
    pr0[3,5] <- log(pr0_nolog[3,5])
    
    zeta0 <- phat[52]
    
    
    pr0test <- exp(pr0)
    sum(exp(t(pr0)))
    
    iter <-1+iter
  }
  
  hat <- optim(c(u01, u02, u03, u04, u05, v02, pr0[1,1:length(pr0)-1], pr0[2,1:length(pr0)-1], pr0[3,1:length(pr0)-1], zeta0), Q,control=opts)
  phat<-hat$par
  Ehat<-hat$value
  
  upih$alpha <- phat[1]
  upih$phi <- phat[2]*100  
  upih$eta <- phat[3]*1000 
  upih$gamma <- phat[4]/10  
  upih$kappa <- phat[5]*1000  
  upiih$alpha <- phat[6] 
  upiih$phi <- phat[7]*100
  upiih$eta <- phat[8]*1000 
  upiih$gamma <- phat[9]/10  
  upiih$kappa <- phat[10]*1000  
  upiiih$alpha <- phat[11]
  upiiih$phi <- phat[12]*100 
  upiiih$eta <- phat[13]*1000 
  upiiih$gamma <- phat[14]/10  
  upiiih$kappa <- phat[15]*1000  
  upivh$alpha <- phat[16] 
  upivh$phi <- phat[17]*100 
  upivh$eta <- phat[18]*1000 
  upivh$gamma <- phat[19]/10  
  upivh$kappa <- phat[20]*1000  
  upvh$alpha <- phat[21] 
  upvh$phi <- phat[22]*100 
  upvh$eta <- phat[23]*1000 
  upvh$gamma <- phat[24]/10  
  upvh$kappa <- phat[25]*1000  
  
  v1h <-hg
  v2h <- exp(phat[26:39])
  
  prh[1,1:4] <- exp(phat[40:43])
  prh[1,5] <- 1-sum(prh[1,1:4])
  
  prh[2,1:4] <- exp(phat[44:47])
  prh[2,5] <- 1-sum(prh[2,1:4])
  
  prh[3,1:4] <- exp(phat[48:51])
  prh[3,5] <- 1-sum(prh[3,1:4])
  
  zetah <- exp(phat[52])
  
  cons1h <- utilconpositive(upih,hg)    
  cons2h <- utilconpositive(upiih,hg)
  cons3h <- utilconpositive(upiiih,hg)
  cons4h <- utilconpositive(upiiih,hg)
  cons5h <- utilconpositive(upiiih,hg)
  
  hmini <- ((1/upih$phi)**(1/upih$gamma))-upih$eta
  hminii <- ((1/upiih$phi)**(1/upiih$gamma))-upiih$eta
  hminiii <- ((1/upiiih$phi)**(1/upiiih$gamma))-upiiih$eta
  hminiv <- ((1/upivh$phi)**(1/upivh$gamma))-upivh$eta
  hminv <- ((1/upvh$phi)**(1/upvh$gamma))-upvh$eta
  
  hminieta <- ((1/upih$phi)**(1/upih$gamma))-upih$eta
  hminiieta <- ((1/upiih$phi)**(1/upiih$gamma))-upiih$eta
  hminiiieta <- ((1/upiiih$phi)**(1/upiiih$gamma))-upiiih$eta
  hminiveta <- ((1/upivh$phi)**(1/upivh$gamma))-upivh$eta
  hminveta <- ((1/upvh$phi)**(1/upvh$gamma))-upvh$eta
  
  umini_m <- log(1-upih$phi*((hmini+upih$eta)**upih$gamma))
  uminii_m <- log(1-upiih$phi*((hminii+upiih$eta)**upiih$gamma))
  uminiii_m <- log(1-upiiih$phi*((hminiii+upiiih$eta)**upiiih$gamma))
  uminiv_m <- log(1-upivh$phi*((hminiv+upivh$eta)**upivh$gamma))
  uminv_m <- log(1-upvh$phi*((hminv+upvh$eta)**upvh$gamma))
  
  umini_0 <- log(1-upih$phi*((0+upih$eta)**upih$gamma))
  uminii_0 <- log(1-upiih$phi*((0+upiih$eta)**upiih$gamma))
  uminiii_0 <- log(1-upiiih$phi*((0+upiiih$eta)**upiiih$gamma))
  uminiv_0 <- log(1-upivh$phi*((0+upivh$eta)**upivh$gamma))
  uminv_0 <- log(1-upvh$phi*((0+upvh$eta)**upvh$gamma))
  
  umini_mo <- log(1-upih$phi*((min(hg)+upih$eta)**upih$gamma))
  uminii_mo <- log(1-upiih$phi*((min(hg)+upiih$eta)**upiih$gamma))
  uminiii_mo <- log(1-upiiih$phi*((min(hg)+upiiih$eta)**upiiih$gamma))
  uminiv_mo <- log(1-upivh$phi*((min(hg)+upivh$eta)**upivh$gamma))
  uminv_mo <- log(1-upvh$phi*((min(hg)+upvh$eta)**upvh$gamma))
  
  yih <- Ycutoff(hg,upih,v1h)
  yimidbin[1] <- (0+yih[1])/2
  yimidbin[2:length(yih)] <- (yih[1:length(yih)-1]+yih[2:length(yih)])/2
  
  yiih <- Ycutoff(hg,upiih,v1h)
  yiimidbin[1] <- (0+yiih[1])/2
  yiimidbin[2:length(yiih)] <- (yiih[1:length(yiih)-1]+yiih[2:length(yiih)])/2
  
  yiiih <- Ycutoff(hg,upiiih,v1h)
  yiiimidbin[1] <- (0+yiiih[1])/2
  yiiimidbin[2:length(yiiih)] <- (yiiih[1:length(yiiih)-1]+yiiih[2:length(yiiih)])/2
  
  yivh <- Ycutoff(hg,upivh,v1h)
  yivmidbin[1] <- (0+yivh[1])/2
  yivmidbin[2:length(yivh)] <- (yivh[1:length(yivh)-1]+yivh[2:length(yivh)])/2
  
  yvh <- Ycutoff(hg,upvh,v1h)
  yvmidbin[1] <- (0+yvh[1])/2
  yvmidbin[2:length(yvh)] <- (yvh[1:length(yvh)-1]+yvh[2:length(yvh)])/2
  
  
  ytih <- Ycutoff(hg,upih,v2h)
  ytimidbin[1] <- (0+ytih[1])/2
  ytimidbin[2:length(ytih)] <- (ytih[1:length(ytih)-1]+ytih[2:length(ytih)])/2
  
  ytiih <- Ycutoff(hg,upiih,v2h)
  ytiimidbin[1] <- (0+ytiih[1])/2
  ytiimidbin[2:length(ytiih)] <- (ytiih[1:length(ytiih)-1]+ytiih[2:length(ytiih)])/2
  
  ytiiih <- Ycutoff(hg,upiiih,v2h)
  ytiiimidbin[1] <- (0+ytiiih[1])/2
  ytiiimidbin[2:length(ytiiih)] <- (ytiiih[1:length(ytiiih)-1]+ytiiih[2:length(ytiiih)])/2
  
  ytivh <- Ycutoff(hg,upivh,v2h)
  ytivmidbin[1] <- (0+ytivh[1])/2
  ytivmidbin[2:length(ytivh)] <- (ytivh[1:length(ytivh)-1]+ytivh[2:length(ytivh)])/2
  
  ytvh <- Ycutoff(hg,upvh,v2h)
  ytvmidbin[1] <- (0+ytvh[1])/2
  ytvmidbin[2:length(ytvh)] <- (ytvh[1:length(ytvh)-1]+ytvh[2:length(ytvh)])/2
  
  
  shi <- v1h[1:length(v1h)-1]/yimidbin
  shii <- v1h[1:length(v1h)-1]/yiimidbin
  shiii <- v1h[1:length(v1h)-1]/yiiimidbin
  shiv <- v1h[1:length(v1h)-1]/yivmidbin
  shv <- v1h[1:length(v1h)-1]/yvmidbin
  
  shti <- v2h[1:length(v2h)-1]/ytimidbin
  shtii <- v2h[1:length(v2h)-1]/ytiimidbin
  shtiii <- v2h[1:length(v2h)-1]/ytiiimidbin
  shtiv <- v2h[1:length(v2h)-1]/ytivmidbin
  shtv <- v2h[1:length(v2h)-1]/ytvmidbin
  
  ehy1i <- IncomeLinearVElast(hg,upih)
  ehy1ii <- IncomeLinearVElast(hg,upiih)
  ehy1iii <- IncomeLinearVElast(hg,upiiih)
  ehy1iv <- IncomeLinearVElast(hg,upivh)
  ehy1v <- IncomeLinearVElast(hg,upvh)
  
  ehy2i <- IncomeElast(hg,v2h, upih)
  ehy2ii <- IncomeElast(hg,v2h,upiih)
  ehy2iii <- IncomeElast(hg,v2h,upiiih)
  ehy2iv <- IncomeElast(hg,v2h,upivh)
  ehy2v <- IncomeElast(hg,v2h,upvh)
  
  ehp1i <- PriceLinearVElast(hg,upih) 
  ehp1ii <- PriceLinearVElast(hg,upiih) 
  ehp1iii <- PriceLinearVElast(hg,upiiih) 
  ehp1iv <- PriceLinearVElast(hg,upivh) 
  ehp1v <- PriceLinearVElast(hg,upvh) 
  
  PIh<-list()
  
  PIh$upi <- upih
  PIh$upii <- upiih
  PIh$upiii <- upiiih
  PIh$upiv <- upivh
  PIh$upv <- upvh
  PIh$v2 <- v2h
  PIh$pr <- prh
  PIh$zeta <- zetah
  
  PIh$elast$ehy1i <- ehy1i
  PIh$elast$ehy1ii <- ehy1ii
  PIh$elast$ehy1iii <- ehy1iii
  PIh$elast$ehy1iv <- ehy1iv
  PIh$elast$ehy1v <- ehy1v
  
  PIh$elast$ehy2i <- ehy2i
  PIh$elast$ehy2ii <- ehy2ii
  PIh$elast$ehy2iii <- ehy2iii
  PIh$elast$ehy2iv <- ehy2iv
  PIh$elast$ehy2v <- ehy2v
  
  
  PIh$elast$ehp1i <- ehp1i
  PIh$elast$ehp1ii <- ehp1ii
  PIh$elast$ehp1iii <- ehp1iii
  PIh$elast$ehp1iv <- ehp1iv
  PIh$elast$ehp1v <- ehp1v
  
  return(list(PIh,Ehat))
}
