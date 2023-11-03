#install.packages('dplyr')
#install.packages('neldermead')
#install.packages("pracma")
#install.packages("signal")
#install.packages("MonoInc")
#install.packages('Rfast')
library(dplyr)
library(neldermead)
library(pracma)
library(signal)
library(MonoInc)
library(Rfast)



# This function loads the data to estimate the model in Dennis, Quintero and Sieg (2019) with 1 metro area.

a1loadinitialdata_1ma <- function(WPmetro) {
  YIk <- WPmetro$YIk
  VIk <- WPmetro$VIk
  PIk <- WPmetro$PIk
  Pop <- WPmetro$Pop
  Corr <- WPmetro$Corr

  return(list(YIk = YIk, VIk = VIk, PIk = PIk, Pop = Pop, Corr = Corr))
}

#Função que estima os parâmetros do modelo usando ambas as regiões metropolitanas e suas correlações

a2Estimate<-function(YI,VI,PI,Pop, Corr){
  start_time<-Sys.time()
  sx<-Pop$I$k1
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
  hg<-VI$V[beg:(length(VI$V)-inicut)]
  Ycut<-YI$Y[begy:(length(YI$Y)-inicut)]
  Ytcut<-YI$Yt[begy:(length(YI$Yt)-inicut)]

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
  F1xob <- cdfGLN4(YI$Y[begy:(length(YI$Y)-inicut)],exp(mu1x),mu1x,sigma1x,r1x,beta1x)

  mu2x<-px$yt[1]
  sigma2x<-px$yt[2]
  r2x<-px$yt[3]
  beta2x<-px$yt[4]
  F2xob <- cdfGLN4(YI$Yt[begy:(length(YI$Yt)-inicut)],exp(mu2x),mu2x,sigma2x,r2x,beta2x)

  mu1xx<-pxx$y[1]
  sigma1xx<-pxx$y[2]
  r1xx<-pxx$y[3]
  beta1xx<-pxx$y[4]
  F1xxob <- cdfGLN4(YI$Y[begy:(length(YI$Y)-inicut)],exp(mu1xx),mu1xx,sigma1xx,r1xx,beta1xx)

  mu2xx<-pxx$yt[1]
  sigma2xx<-pxx$yt[2]
  r2xx<-pxx$yt[3]
  beta2xx<-pxx$yt[4]
  F2xxob <- cdfGLN4(YI$Yt[begy:(length(YI$Yt)-inicut)],exp(mu2xx),mu2xx,sigma2xx,r2xx,beta2xx)

  mu1xxx<-pxxx$y[1]
  sigma1xxx<-pxxx$y[2]
  r1xxx<-pxxx$y[3]
  beta1xxx<-pxxx$y[4]
  F1xxxob <- cdfGLN4(YI$Y[begy:(length(YI$Y)-inicut)],exp(mu1xxx),mu1xxx,sigma1xxx,r1xxx,beta1xxx)

  mu2xxx<-pxxx$yt[1]
  sigma2xxx<-pxxx$yt[2]
  r2xxx<-pxxx$yt[3]
  beta2xxx<-pxxx$yt[4]
  F2xxxob <- cdfGLN4(YI$Yt[begy:(length(YI$Yt)-inicut)],exp(mu2xxx),mu2xxx,sigma2xxx,r2xxx,beta2xxx)


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
  G1xob <- cdfGLN4(VI$V[begy:(length(VI$V)-inicut)],exp(omega1x),omega1x,tau1x,m1x,theta1x)

  omega2x<-px$vt[1]
  tau2x<-px$vt[2]
  m2x<-px$vt[3]
  theta2x<-px$vt[4]
  G2xob <- cdfGLN4(VI$Vt[begy:(length(VI$Vt)-inicut)],exp(omega2x),omega2x,tau2x,m2x,theta2x)

  omega1xx<-pxx$v[1]
  tau1xx<-pxx$v[2]
  m1xx<-pxx$v[3]
  theta1xx<-pxx$v[4]
  G1xxob <- cdfGLN4(VI$V[begy:(length(VI$V)-inicut)],exp(omega1xx),omega1xx,tau1xx,m1xx,theta1xx)

  omega2xx<-pxx$vt[1]
  tau2xx<-pxx$vt[2]
  m2xx<-pxx$vt[3]
  theta2xx<-pxx$vt[4]
  G2xxob <- cdfGLN4(VI$Vt[begy:(length(VI$Vt)-inicut)],exp(omega2xx),omega2xx,tau2xx,m2xx,theta2xx)

  omega1xxx<-pxxx$v[1]
  tau1xxx<-pxxx$v[2]
  m1xxx<-pxxx$v[3]
  theta1xxx<-pxxx$v[4]
  G1xxxob <- cdfGLN4(VI$V[begy:(length(VI$V)-inicut)],exp(omega1xxx),omega1xxx,tau1xxx,m1xxx,theta1xxx)

  omega2xxx<-pxxx$vt[1]
  tau2xxx<-pxxx$vt[2]
  m2xxx<-pxxx$vt[3]
  theta2xxx<-pxxx$vt[4]
  G2xxxob <- cdfGLN4(VI$Vt[begy:(length(VI$Vt)-inicut)],exp(omega2xxx),omega2xxx,tau2xxx,m2xxx,theta2xxx)

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

  u2$alpha <- u2$alpha+runif(1)/10
  u3$alpha <- u3$alpha+runif(1)/10
  u4$alpha <- u4$alpha+runif(1)/10
  u5$alpha <- u5$alpha+runif(1)/10

  u2$gamma <- u2$gamma-runif(1)/10
  u3$gamma <- u3$gamma-runif(1)/10
  u4$gamma <- u4$gamma-runif(1)/10
  u5$gamma <- u5$gamma-runif(1)/10


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
  
#  u1<-list('alpha'=1.22,'phi'=9.29,'eta'=0.42,'gamma'=-1.77,'kappa'=0)
#  u2<-list('alpha'=1.61,'phi'=4.32,'eta'=1.19,'gamma'=-1.22,'kappa'=0)
#  u3<-list('alpha'=2.12,'phi'=1.15,'eta'=5.71,'gamma'=-1.618,'kappa'=0)
#  u4<-list('alpha'=1.11,'phi'=9.11,'eta'=0.99,'gamma'=-1.51,'kappa'=0)
#  u5<-list('alpha'=1.33,'phi'=7.77,'eta'=1.98,'gamma'=-1.99,'kappa'=0)
  
  u01 <- c(u1$alpha, (u1$phi)/100, (u1$eta)/1000, u1$gamma*10, (u1$kappa)/1000)
  u02 <- c(u2$alpha, (u2$phi)/100, (u2$eta)/1000, u2$gamma*10, (u2$kappa)/1000)
  u03 <- c(u3$alpha, (u3$phi)/100, (u3$eta)/1000, u3$gamma*10, (u3$kappa)/1000)
  u04 <- c(u4$alpha, (u4$phi)/100, (u4$eta)/1000, u4$gamma*10, (u4$kappa)/1000)
  u05 <- c(u5$alpha, (u5$phi)/100, (u5$eta)/1000, u5$gamma*10, (u5$kappa)/1000)


  zeta0 <- log(0.065)

  Q<-function(p_hat){
    upi<-list()
    upi$alpha <- p_hat[1]
    upi$phi <- p_hat[2]*100
    upi$eta <- p_hat[3]*1000
    upi$gamma <- p_hat[4]/10
    upi$kappa <- p_hat[5]*1000

    upii<-list()
    upii$alpha <- p_hat[6]
    upii$phi <- p_hat[7]*100
    upii$eta <- p_hat[8]*1000
    upii$gamma <- p_hat[9]/10
    upii$kappa <- p_hat[10]*1000

    upiii<-list()
    upiii$alpha <- p_hat[11]
    upiii$phi <- p_hat[12]*100
    upiii$eta <- p_hat[13]*1000
    upiii$gamma <- p_hat[14]/10
    upiii$kappa <- p_hat[15]*1000

    upiv<-list()
    upiv$alpha <- p_hat[16]
    upiv$phi <- p_hat[17]*100
    upiv$eta <- p_hat[18]*1000
    upiv$gamma <- p_hat[19]/10
    upiv$kappa <- p_hat[20]*1000

    upv<-list()
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


    if (all(cons1>0) & all(cons2>0) & all(cons3>0) & all(cons4>0) & all(cons5>0) & upi$alpha>0 & upii$alpha>0 &
        upiii$alpha>0  & upiv$alpha>0 & upv$alpha>0 & upi$gamma<0 & upii$gamma<0 & upiii$gamma<0 &
        upiv$gamma<0 & upv$gamma<0 & upi$phi>0 & upii$phi>0 & upiii$phi>0 & upiv$phi>0 & upv$phi>0 &
        upi$eta>0 & upii$eta>0 & upiii$eta>0 & upiii$eta>0 & upiii$eta>0 & sum(p1)>0 & sum(p2)>0 &
        sum(p3)>0 & all(p1 > 0) & all(p2 > 0) & all(p3 > 0)){
      Hinv<-Hinv3x5INorm1(hg,vp2,px$y,pxx$y,pxxx$y,px$yt,pxx$yt,pxxx$yt,upi,upii,upiii,upiv,upv,sx,sxx,sxxx,p1,p2,p3,zeta,VI)
      nomes<-c('Rh','cond1','cond2','Hdcdi','Hdcdii','Hdcdiii','Hdcdiv','Hdcdv','Hx','Hxx','Hxxx','AgHd','Hdtcdi','Hdtcdii','Hdtcdiii','Hdtcdiv','Hdtcdv','Htx','Htxx','Htxxx','AgHdt','Rh2')
      for (i in 1:length(names(Hinv))){
        assign(nomes[i],unlist(Hinv[[names(Hinv)[i]]]))
      }
      if(cond1 & cond2){
        eq <- 10*(AgHdt-Rh2)
        eqs <- sum(eq**2)
        eqslow <- sum(eq[1:7]**2)
        eqsup <- sum(eq[8:length(eq)]**2)

        dif1x <- 10*(cdfGLN4(vp1[1:(length(vp1)-1)],exp(omega1x),omega1x,tau1x,m1x,theta1x)-t(Hx))
        verr1x <- sum(dif1x**2)
        verr1xlow <- sum(dif1x[1:7]**2)
        verr1xup <- sum(dif1x[8:length(dif1x)]**2)

        dif1xx <- 10*(cdfGLN4(vp1[1:(length(vp1)-1)],exp(omega1xx),omega1xx,tau1xx,m1xx,theta1xx)-t(Hxx))
        verr1xx <- sum(dif1xx**2)
        verr1xxlow <- sum(dif1xx[1:7]**2)
        verr1xxup <- sum(dif1xx[8:length(dif1xx)]**2)

        dif1xxx <- 10*(cdfGLN4(vp1[1:(length(vp1)-1)],exp(omega1xxx),omega1xxx,tau1xxx,m1xxx,theta1xxx)-t(Hxxx))
        verr1xxx <- sum(dif1xxx**2)
        verr1xxxlow <- sum(dif1xxx[1:7]**2)
        verr1xxxup <- sum(dif1xxx[8:length(dif1xxx)]**2)

        dif2x <- 10*(cdfGLN4(vp2[1:(length(vp2)-1)],exp(omega2x),omega2x,tau2x,m2x,theta2x)-t(Htx))
        verr2x <- sum(dif2x**2)
        verr2xlow <- sum(dif2x[1:7]**2)
        verr2xup <- sum(dif2x[8:length(dif2x)]**2)

        dif2xx <- 10*(cdfGLN4(vp2[1:(length(vp2)-1)],exp(omega2xx),omega2xx,tau2xx,m2xx,theta2xx)-t(Htxx))
        verr2xx <- sum(dif2xx**2)
        verr2xxlow <- sum(dif2xx[1:7]**2)
        verr2xxup <- sum(dif2xx[8:length(dif2xx)]**2)

        dif2xxx <- 10*(cdfGLN4(vp2[1:(length(vp2)-1)],exp(omega2xxx),omega2xxx,tau2xxx,m2xxx,theta2xxx)-t(Htxxx))
        verr2xxx <- sum(dif2xxx**2)
        verr2xxxlow <- sum(dif2xxx[1:7]**2)
        verr2xxxup <- sum(dif2xxx[8:length(dif2xxx)]**2)

        CorrMoment <- Yvcorrelations(sx, sxx, sxxx, p1, p2, p3, hg, upi,upii,upiii,upiv,upv, vp1, vp2, Corr)
        F <- 10*((eqslow + eqsup) + (verr1xlow + verr1xup +verr1xxlow + verr1xxup+verr1xxxlow + verr1xxxup +
                                       verr2xlow + verr2xup +verr2xxlow + verr2xxup +verr2xxxlow + verr2xxxup+CorrMoment))
      }else{
        F <- 1000000000000000000000000
      }}else{
        F <- 1000000000000000000000000
      }
    return(F)
  }

  opts<-list(abstol=1e-4,reltol=1e-4,maxit=1)
  
  iter <-50
  for (i in 1:iter){
    hat <- optim(c(u01, u02, u03, u04, u05, v02, pr0[1,1:(length(pr0[1,])-1)], pr0[2,1:(length(pr0[2,])-1)], pr0[3,(1:length(pr0[3,])-1)], zeta0), Q,control=opts)
    phat<-hat$par
    Ehat<-hat$value
    u01 <- c(phat[1], phat[2], phat[3], phat[4], phat[5])
    u02 <- c(phat[6], phat[7], phat[8], phat[9], phat[10])
    u03 <- c(phat[11], phat[12], phat[13], phat[14], phat[15])
    u04 <- c(phat[16], phat[17], phat[18], phat[19], phat[20])
    u05 <- c(phat[21], phat[22], phat[23], phat[24], phat[25])
    v02 <- phat[26:39]

    pr0_nolog<-matrix(nrow = 3, ncol=5)
    
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
    print(iter)
  }

  hat <- optim(c(u01, u02, u03, u04, u05, v02, pr0[1,(1:length(pr0[1,])-1)], pr0[2,(1:length(pr0[2,])-1)], pr0[3,1:(length(pr0[3,])-1)], zeta0), Q,control=opts)
  phat<-hat$par
  Ehat<-hat$value
  
  upih<-list()
  upiih<-list()
  upiiih<-list()
  upivh<-list()
  upvh<-list()
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
  
  prh<-matrix(nrow = 3, ncol=5)
  
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
  yimidbin <- (0+yih[1])/2
  yimidbin[2:length(yih)] <- (yih[1:(length(yih)-1)]+yih[2:length(yih)])/2

  yiih <- Ycutoff(hg,upiih,v1h)
  yiimidbin <- (0+yiih[1])/2
  yiimidbin[2:length(yiih)] <- (yiih[1:(length(yiih)-1)]+yiih[2:length(yiih)])/2

  yiiih <- Ycutoff(hg,upiiih,v1h)
  yiiimidbin <- (0+yiiih[1])/2
  yiiimidbin[2:length(yiiih)] <- (yiiih[1:(length(yiiih)-1)]+yiiih[2:length(yiiih)])/2

  yivh <- Ycutoff(hg,upivh,v1h)
  yivmidbin <- (0+yivh[1])/2
  yivmidbin[2:length(yivh)] <- (yivh[1:(length(yivh)-1)]+yivh[2:length(yivh)])/2

  yvh <- Ycutoff(hg,upvh,v1h)
  yvmidbin <- (0+yvh[1])/2
  yvmidbin[2:length(yvh)] <- (yvh[1:(length(yvh)-1)]+yvh[2:length(yvh)])/2


  ytih <- Ycutoff(hg,upih,v2h)
  ytimidbin <- (0+ytih[1])/2
  ytimidbin[2:length(ytih)] <- (ytih[1:(length(ytih)-1)]+ytih[2:length(ytih)])/2

  ytiih <- Ycutoff(hg,upiih,v2h)
  ytiimidbin <- (0+ytiih[1])/2
  ytiimidbin[2:length(ytiih)] <- (ytiih[1:(length(ytiih)-1)]+ytiih[2:length(ytiih)])/2

  ytiiih <- Ycutoff(hg,upiiih,v2h)
  ytiiimidbin <- (0+ytiiih[1])/2
  ytiiimidbin[2:length(ytiiih)] <- (ytiiih[1:(length(ytiiih)-1)]+ytiiih[2:length(ytiiih)])/2

  ytivh <- Ycutoff(hg,upivh,v2h)
  ytivmidbin <- (0+ytivh[1])/2
  ytivmidbin[2:length(ytivh)] <- (ytivh[1:(length(ytivh)-1)]+ytivh[2:length(ytivh)])/2

  ytvh <- Ycutoff(hg,upvh,v2h)
  ytvmidbin <- (0+ytvh[1])/2
  ytvmidbin[2:length(ytvh)] <- (ytvh[1:(length(ytvh)-1)]+ytvh[2:length(ytvh)])/2


  shi <- v1h[1:(length(v1h)-1)]/yimidbin
  shii <- v1h[1:(length(v1h)-1)]/yiimidbin
  shiii <- v1h[1:(length(v1h)-1)]/yiiimidbin
  shiv <- v1h[1:(length(v1h)-1)]/yivmidbin
  shv <- v1h[1:(length(v1h)-1)]/yvmidbin

  shti <- v2h[1:(length(v2h)-1)]/ytimidbin
  shtii <- v2h[1:(length(v2h)-1)]/ytiimidbin
  shtiii <- v2h[1:(length(v2h)-1)]/ytiiimidbin
  shtiv <- v2h[1:(length(v2h)-1)]/ytivmidbin
  shtv <- v2h[1:(length(v2h)-1)]/ytvmidbin

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
  
  print(Sys.time()-start_time)
  return(list('PIh'=PIh,'Ehat'=Ehat))
}

#Função que estima os parâmetros do modelo usando ambas as regiões metropolitanas e suas correlações


library(dplyr)

a3Estimate<-function(YI,VI,PI,Pop,YIm,VIm,PIm,Popm, Corr){
  sx<-Pop$I$k1m
  sxm<-Popm$I$k1m
  sxx<-Pop$I$k2m
  sxxm<-Popm$I$k2m
  sxxx<-Pop$I$k3m
  sxxxm<-Popm$I$k3m

  typex<-'Type 1'
  typexx<-'Type 2'
  typexxx<-'Type 3'
  typexm<-'Type 1'
  typexxm<-'Type 2'
  typexxxm<-'Type 3'

  inicut<-2
  begy<-2
  beg<-3
  hg<-VI$V[beg:(length(VI$V)-inicut)]
  Ycut<-YI$Y[begy:(length(YI$Y)-inicut)]
  Ytcut<-YI$Yt[begy:(length(YI$Yt)-inicut)]

  Ycutm<-YIm$Y[begy:(length(YIm$Y)-inicut)]
  Ytcutm<-YIm$Yt[begy:(length(YIm$Yt)-inicut)]

  px<-list()
  pxx<-list()
  pxxx<-list()
  px$y<-PI$GLN4ui$y
  pxx$y<-PI$GLN4uii$y
  pxxx$y<-PI$GLN4uiii$y
  px$yt<-PI$GLN4ui$yt
  pxx$yt<-PI$GLN4uii$yt
  pxxx$yt<-PI$GLN4uiii$yt

  pxm<-list()
  pxxm<-list()
  pxxxm<-list()
  pxm$y<-PIm$GLN4ui$y
  pxxm$y<-PIm$GLN4uii$y
  pxxxm$y<-PIm$GLN4uiii$y
  pxm$yt<-PIm$GLN4ui$yt
  pxxm$yt<-PIm$GLN4uii$yt
  pxxxm$yt<-PIm$GLN4uiii$yt

  mu1x<-px$y[1]
  sigma1x<-px$y[2]
  r1x<-px$y[3]
  beta1x<-px$y[4]
  F1xob <- cdfGLN4(YI$Y[begy:(length(YI$Y)-inicut)],exp(mu1x),mu1x,sigma1x,r1x,beta1x)

  mu2x<-px$yt[1]
  sigma2x<-px$yt[2]
  r2x<-px$yt[3]
  beta2x<-px$yt[4]
  F2xob <- cdfGLN4(YI$Yt[begy:(length(YI$Yt)-inicut)],exp(mu2x),mu2x,sigma2x,r2x,beta2x)

  mu1xx<-pxx$y[1]
  sigma1xx<-pxx$y[2]
  r1xx<-pxx$y[3]
  beta1xx<-pxx$y[4]
  F1xxob <- cdfGLN4(YI$Y[begy:(length(YI$Y)-inicut)],exp(mu1xx),mu1xx,sigma1xx,r1xx,beta1xx)

  mu2xx<-pxx$yt[1]
  sigma2xx<-pxx$yt[2]
  r2xx<-pxx$yt[3]
  beta2xx<-pxx$yt[4]
  F2xxob <- cdfGLN4(YI$Yt[begy:(length(YI$Yt)-inicut)],exp(mu2xx),mu2xx,sigma2xx,r2xx,beta2xx)

  mu1xxx<-pxxx$y[1]
  sigma1xxx<-pxxx$y[2]
  r1xxx<-pxxx$y[3]
  beta1xxx<-pxxx$y[4]
  F1xxxob <- cdfGLN4(YI$Y[begy:(length(YI$Y)-inicut)],exp(mu1xxx),mu1xxx,sigma1xxx,r1xxx,beta1xxx)

  mu2xxx<-pxxx$yt[1]
  sigma2xxx<-pxxx$yt[2]
  r2xxx<-pxxx$yt[3]
  beta2xxx<-pxxx$yt[4]
  F2xxxob <- cdfGLN4(YI$Yt[begy:(length(YI$Yt)-inicut)],exp(mu2xxx),mu2xxx,sigma2xxx,r2xxx,beta2xxx)


  mu1xm<-pxm$y[1]
  sigma1xm<-pxm$y[2]
  r1xm<-pxm$y[3]
  beta1xm<-pxm$y[4]
  F1xobm <- cdfGLN4(YIm$Y[begy:(length(YIm$Y)-inicut)],exp(mu1xm),mu1xm,sigma1xm,r1xm,beta1xm)

  mu2xm<-pxm$yt[1]
  sigma2xm<-pxm$yt[2]
  r2xm<-pxm$yt[3]
  beta2xm<-pxm$yt[4]
  F2xobm <- cdfGLN4(YIm$Yt[begy:(length(YIm$Yt)-inicut)],exp(mu2xm),mu2xm,sigma2xm,r2xm,beta2xm)

  mu1xxm<-pxxm$y[1]
  sigma1xxm<-pxxm$y[2]
  r1xxm<-pxxm$y[3]
  beta1xxm<-pxxm$y[4]
  F1xxobm <- cdfGLN4(YIm$Y[begy:(length(YIm$Y)-inicut)],exp(mu1xxm),mu1xxm,sigma1xxm,r1xxm,beta1xxm)

  mu2xxm<-pxxm$yt[1]
  sigma2xxm<-pxxm$yt[2]
  r2xxm<-pxxm$yt[3]
  beta2xxm<-pxxm$yt[4]
  F2xxobm <- cdfGLN4(YIm$Yt[begy:(length(YIm$Yt)-inicut)],exp(mu2xxm),mu2xxm,sigma2xxm,r2xxm,beta2xxm)

  mu1xxxm<-pxxxm$y[1]
  sigma1xxxm<-pxxxm$y[2]
  r1xxxm<-pxxxm$y[3]
  beta1xxxm<-pxxxm$y[4]
  F1xxxobm <- cdfGLN4(YIm$Y[begy:(length(YIm$Y)-inicut)],exp(mu1xxxm),mu1xxxm,sigma1xxxm,r1xxxm,beta1xxxm)

  mu2xxxm<-pxxxm$yt[1]
  sigma2xxxm<-pxxxm$yt[2]
  r2xxxm<-pxxxm$yt[3]
  beta2xxxm<-pxxxm$yt[4]
  F2xxxobm <- cdfGLN4(YIm$Yt[begy:(length(YIm$Yt)-inicut)],exp(mu2xxxm),mu2xxxm,sigma2xxxm,r2xxxm,beta2xxxm)
  

  px$v<-PI$GLN4ui$v
  pxx$v<-PI$GLN4uii$v
  pxxx$v<-PI$GLN4uiii$v
  px$vt<-PI$GLN4ui$vt
  pxx$vt<-PI$GLN4uii$vt
  pxxx$vt<-PI$GLN4uiii$vt

  pxm$v<-PIm$GLN4ui$v
  pxxm$v<-PIm$GLN4uii$v
  pxxxm$v<-PIm$GLN4uiii$v
  pxm$vt<-PIm$GLN4ui$vt
  pxxm$vt<-PIm$GLN4uii$vt
  pxxxm$vt<-PIm$GLN4uiii$vt

  omega1x<-px$v[1]
  tau1x<-px$v[2]
  m1x<-px$v[3]
  theta1x<-px$v[4]
  G1xob <- cdfGLN4(VI$V[begy:(length(VI$V)-inicut)],exp(omega1x),omega1x,tau1x,m1x,theta1x)

  omega2x<-px$vt[1]
  tau2x<-px$vt[2]
  m2x<-px$vt[3]
  theta2x<-px$vt[4]
  G2xob <- cdfGLN4(VI$Vt[begy:(length(VI$Vt)-inicut)],exp(omega2x),omega2x,tau2x,m2x,theta2x)

  omega1xx<-pxx$v[1]
  tau1xx<-pxx$v[2]
  m1xx<-pxx$v[3]
  theta1xx<-pxx$v[4]
  G1xxob <- cdfGLN4(VI$V[begy:(length(VI$V)-inicut)],exp(omega1xx),omega1xx,tau1xx,m1xx,theta1xx)

  omega2xx<-pxx$vt[1]
  tau2xx<-pxx$vt[2]
  m2xx<-pxx$vt[3]
  theta2xx<-pxx$vt[4]
  G2xxob <- cdfGLN4(VI$Vt[begy:(length(VI$Vt)-inicut)],exp(omega2xx),omega2xx,tau2xx,m2xx,theta2xx)

  omega1xxx<-pxxx$v[1]
  tau1xxx<-pxxx$v[2]
  m1xxx<-pxxx$v[3]
  theta1xxx<-pxxx$v[4]
  G1xxxob <- cdfGLN4(VI$V[begy:(length(VI$V)-inicut)],exp(omega1xxx),omega1xxx,tau1xxx,m1xxx,theta1xxx)

  omega2xxx<-pxxx$vt[1]
  tau2xxx<-pxxx$vt[2]
  m2xxx<-pxxx$vt[3]
  theta2xxx<-pxxx$vt[4]
  G2xxxob <- cdfGLN4(VI$Vt[begy:(length(VI$Vt)-inicut)],exp(omega2xxx),omega2xxx,tau2xxx,m2xxx,theta2xxx)

  omega1xm<-pxm$v[1]
  tau1xm<-pxm$v[2]
  m1xm<-pxm$v[3]
  theta1xm<-pxm$v[4]
  G1xobm <- cdfGLN4(VIm$V[begy:(length(VIm$V)-inicut)],exp(omega1xm),omega1xm,tau1xm,m1xm,theta1xm)

  omega2xm<-pxm$vt[1]
  tau2xm<-pxm$vt[2]
  m2xm<-pxm$vt[3]
  theta2xm<-pxm$vt[4]
  G2xobm <- cdfGLN4(VIm$Vt[begy:(length(VIm$Vt)-inicut)],exp(omega2xm),omega2xm,tau2xm,m2xm,theta2xm)

  omega1xxm<-pxxm$v[1]
  tau1xxm<-pxxm$v[2]
  m1xxm<-pxxm$v[3]
  theta1xxm<-pxxm$v[4]
  G1xxobm <- cdfGLN4(VIm$V[begy:(length(VIm$V)-inicut)],exp(omega1xxm),omega1xxm,tau1xxm,m1xxm,theta1xxm)

  omega2xxm<-pxxm$vt[1]
  tau2xxm<-pxxm$vt[2]
  m2xxm<-pxxm$vt[3]
  theta2xxm<-pxxm$vt[4]
  G2xxobm <- cdfGLN4(VIm$Vt[begy:(length(VIm$Vt)-inicut)],exp(omega2xxm),omega2xxm,tau2xxm,m2xxm,theta2xxm)

  omega1xxx<-pxxxm$v[1]
  tau1xxxm<-pxxxm$v[2]
  m1xxxm<-pxxxm$v[3]
  theta1xxxm<-pxxxm$v[4]
  G1xxxobm <- cdfGLN4(VI$V[begy:(length(VIm$V)-inicut)],exp(omega1xxx),omega1xxx,tau1xxx,m1xxx,theta1xxx)

  omega1xxxm<-pxxxm$vt[1]
  tau2xxxm<-pxxxm$vt[2]
  m1xxxm<-pxxxm$vt[3]
  theta1xxxm<-pxxxm$vt[4]
  G1xxxobm <- cdfGLN4(VIm$V[begy:(length(VIm$Vt)-inicut)],exp(omega1xxxm),omega1xxxm,tau1xxxm,m1xxxm,theta1xxxm)
  
  omega2xxx<-pxxxm$v[1]
  tau2xxxm<-pxxxm$v[2]
  m2xxxm<-pxxxm$v[3]
  theta2xxxm<-pxxxm$v[4]
  G2xxxobm <- cdfGLN4(VI$Vt[begy:(length(VIm$V)-inicut)],exp(omega2xxx),omega2xxx,tau2xxx,m2xxx,theta2xxx)
  
  omega2xxxm<-pxxxm$vt[1]
  tau2xxxm<-pxxxm$vt[2]
  m2xxxm<-pxxxm$vt[3]
  theta2xxxm<-pxxxm$vt[4]
  G2xxxobm <- cdfGLN4(VIm$Vt[begy:(length(VIm$Vt)-inicut)],exp(omega2xxxm),omega2xxxm,tau2xxxm,m2xxxm,theta2xxxm)
  
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

  u2$alpha <- u2$alpha+runif(1)/10
  u3$alpha <- u3$alpha+runif(1)/10
  u4$alpha <- u4$alpha+runif(1)/10
  u5$alpha <- u5$alpha+runif(1)/10

  u2$gamma <- u2$gamma-runif(1)/10
  u3$gamma <- u3$gamma-runif(1)/10
  u4$gamma <- u4$gamma-runif(1)/10
  u5$gamma <- u5$gamma-runif(1)/10


  v01 <- log(hg)
  v02 <- log(hg)
  v01m <- log(hg)
  v02m <- log(hg)


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

  #  u1<-list('alpha'=1.22,'phi'=9.29,'eta'=0.42,'gamma'=-1.77,'kappa'=0)
  #  u2<-list('alpha'=1.61,'phi'=4.32,'eta'=1.19,'gamma'=-1.22,'kappa'=0)
  #  u3<-list('alpha'=2.12,'phi'=1.15,'eta'=5.71,'gamma'=-1.618,'kappa'=0)
  #  u4<-list('alpha'=1.11,'phi'=9.11,'eta'=0.99,'gamma'=-1.51,'kappa'=0)
  #  u5<-list('alpha'=1.33,'phi'=7.77,'eta'=1.98,'gamma'=-1.99,'kappa'=0)
  
  u01 <- c(u1$alpha, (u1$phi)/100, (u1$eta)/1000, u1$gamma*10, (u1$kappa)/1000)
  u02 <- c(u2$alpha, (u2$phi)/100, (u2$eta)/1000, u2$gamma*10, (u2$kappa)/1000)
  u03 <- c(u3$alpha, (u3$phi)/100, (u3$eta)/1000, u3$gamma*10, (u3$kappa)/1000)
  u04 <- c(u4$alpha, (u4$phi)/100, (u4$eta)/1000, u4$gamma*10, (u4$kappa)/1000)
  u05 <- c(u5$alpha, (u5$phi)/100, (u5$eta)/1000, u5$gamma*10, (u5$kappa)/1000)


  zeta0 <- log(0.065)

  Q<-function(p_hat){
    upi<-list()
    upii<-list()
    upiii<-list()
    upiv<-list()
    upv<-list()
    
    
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
    vp1m <- exp(p_hat[40:53])
    vp2m <- exp(p_hat[54:67])

    p1 <- exp(p_hat[68:71])
    p1 <- append(p1,1-sum(p1))

    p2 <- exp(p_hat[72:75])
    p2 <- append(p2,1-sum(p2))

    p3 <- exp(p_hat[76:79])
    p3 <- append(p3,1-sum(p3))

    zeta <- exp(p_hat[80])


    cons1 <- utilconpositive(upi,hg)
    cons2 <- utilconpositive(upii,hg)
    cons3 <- utilconpositive(upiii,hg)
    cons4 <- utilconpositive(upiv,hg)
    cons5 <- utilconpositive(upv,hg)


    if (all(cons1>0) & all(cons2>0) & all(cons3>0) & all(cons4>0) & all(cons5>0) & upi$alpha>0 & upii$alpha>0 &
        upiii$alpha>0  & upiv$alpha>0 & upv$alpha>0 & upi$gamma<0 & upii$gamma<0 & upiii$gamma<0 &
        upiv$gamma<0 & upv$gamma<0 & upi$phi>0 & upii$phi>0 & upiii$phi>0 & upiv$phi>0 & upv$phi>0 &
        upi$eta>0 & upii$eta>0 & upiii$eta>0 & upiii$eta>0 & upiii$eta>0 & sum(p1)>0 & sum(p2)>0 &
        sum(p3)>0 & all(p1 > 0) & all(p2 > 0) & all(p3 > 0)){
      Hinv<-Hinv3x5INorm1(hg,vp2,px$y,pxx$y,pxxx$y,px$yt,pxx$yt,pxxx$yt,upi,upii,upiii,upiv,upv,sx,sxx,sxxx,p1,p2,p3,zeta,VI)
      nomes<-c('Rh','cond1','cond2','Hdcdi','Hdcdii','Hdcdiii','Hdcdiv','Hdcdv','Hx','Hxx','Hxxx','AgHd','Hdtcdi','Hdtcdii','Hdtcdiii','Hdtcdiv','Hdtcdv','Htx','Htxx','Htxxx','AgHdt','Rh2')
      for (i in 1:length(names(Hinv))){
        assign(nomes[i],unlist(Hinv[[names(Hinv)[i]]]))
      }
      Hinvm<-Hinv3x5INorm1m(hg,vp2m,pxm$y,pxxm$y,pxxxm$y,pxm$yt,pxxm$yt,pxxxm$yt,upi,upii,upiii,upiv,upv,sxm,sxxm,sxxxm,p1,p2,p3,zeta,VIm,vp1m,VI,Rh)
      nomes<-c('Rhm','cond1m','cond2m','Hdcdim','Hdcdiim','Hdcdiiim','Hdcdivm','Hdcdvm','Hxm','Hxxm','Hxxxm','AgHdm','Hdtcdim','Hdtcdiim','Hdtcdiiim','Hdtcdivm','Hdtcdvm','Htxm','Htxxm','Htxxxm','AgHdtm','Rh2m')
      for (i in 1:length(names(Hinvm))){
        assign(nomes[i],unlist(Hinvm[[names(Hinvm)[i]]]))
      }
      if(cond1 & cond2 & cond1m & cond2m){

        eq1m <- 10*(AgHdtm-Rhm)
        eq1ms <- sum(eq1m**2)
        eq1mslow <- sum(eq1m[1:7]**2)
        eq1msup <- sum(eq1m[8:length(eq1m)]**2)

        eq2m <- 10*(AgHdtm-Rh2m)
        eq2ms <- sum(eq2m**2)
        eq2mslow <- sum(eq2m[1:7]**2)
        eq2msup <- sum(eq2m[8:length(eq2m)]**2)

        eq2 <- 10*(AgHdt-Rh2)
        eq2s <- sum(eq2**2)
        eq2slow <- sum(eq2[1:7]**2)
        eq2sup <- sum(eq2[8:length(eq2)]**2)

        dif1x <- 10*(cdfGLN4(vp1[1:(length(vp1)-1)],exp(omega1x),omega1x,tau1x,m1x,theta1x)-t(Hx))
        verr1x <- sum(dif1x**2)
        verr1xlow <- sum(dif1x[1:7]**2)
        verr1xup <- sum(dif1x[8:length(dif1x)]**2)

        dif1xx <- 10*(cdfGLN4(vp1[1:(length(vp1)-1)],exp(omega1xx),omega1xx,tau1xx,m1xx,theta1xx)-t(Hxx))
        verr1xx <- sum(dif1xx**2)
        verr1xxlow <- sum(dif1xx[1:7]**2)
        verr1xxup <- sum(dif1xx[8:length(dif1xx)]**2)

        dif1xxx <- 10*(cdfGLN4(vp1[1:(length(vp1)-1)],exp(omega1xxx),omega1xxx,tau1xxx,m1xxx,theta1xxx)-t(Hxxx))
        verr1xxx <- sum(dif1xxx**2)
        verr1xxxlow <- sum(dif1xxx[1:7]**2)
        verr1xxxup <- sum(dif1xxx[8:length(dif1xxx)]**2)

        dif2x <- 10*(cdfGLN4(vp2[1:(length(vp2)-1)],exp(omega2x),omega2x,tau2x,m2x,theta2x)-t(Htx))
        verr2x <- sum(dif2x**2)
        verr2xlow <- sum(dif2x[1:7]**2)
        verr2xup <- sum(dif2x[8:length(dif2x)]**2)

        dif2xx <- 10*(cdfGLN4(vp2[1:(length(vp2)-1)],exp(omega2xx),omega2xx,tau2xx,m2xx,theta2xx)-t(Htxx))
        verr2xx <- sum(dif2xx**2)
        verr2xxlow <- sum(dif2xx[1:7]**2)
        verr2xxup <- sum(dif2xx[8:length(dif2xx)]**2)

        dif2xxx <- 10*(cdfGLN4(vp2[1:length(vp2)-1],exp(omega2xxx),omega2xxx,tau2xxx,m2xxx,theta2xxx)-t(Htxxx))
        verr2xxx <- sum(dif2xxx**2)
        verr2xxxlow <- sum(dif2xxx[1:7]**2)
        verr2xxxup <- sum(dif2xxx[8:length(dif2xxx)]**2)

        dif1xm <- 10*(cdfGLN4(vp1m[1:(length(vp1m)-1)],exp(omega1xm),omega1xm,tau1xm,m1xm,theta1xm)-t(Hxm))
        verr1xm <- sum(dif1xm**2)
        verr1xmlow <- sum(dif1xm[1:7]**2)
        verr1xmup <- sum(dif1xm[8:length(dif1xm)]**2)

        dif1xxm <- 10*(cdfGLN4(vp1m[1:(length(vp1m)-1)],exp(omega1xxm),omega1xxm,tau1xxm,m1xxm,theta1xxm)-t(Hxxm))
        verr1xxm <- sum(dif1xxm**2)
        verr1xxmlow <- sum(dif1xxm[1:7]**2)
        verr1xxmup <- sum(dif1xxm[8:length(dif1xxm)]**2)

        dif1xxxm <- 10*(cdfGLN4(vp1m[1:(length(vp1m)-1)],exp(omega1xxxm),omega1xxxm,tau1xxxm,m1xxxm,theta1xxxm)-t(Hxxxm))
        verr1xxxm <- sum(dif1xxxm**2)
        verr1xxxmlow <- sum(dif1xxxm[1:7]**2)
        verr1xxxmup <- sum(dif1xxxm[8:length(dif1xxxm)]**2)

        dif2xm <- 10*(cdfGLN4(vp2m[1:(length(vp2m)-1)],exp(omega2xm),omega2xm,tau2xm,m2xm,theta2xm)-t(Htxm))
        verr2xm <- sum(dif2xm**2)
        verr2xmlow <- sum(dif2xm[1:7]**2)
        verr2xmup <- sum(dif2xm[8:length(dif2xm)]**2)

        dif2xxm <- 10*(cdfGLN4(vp2m[1:(length(vp2m)-1)],exp(omega2xxm),omega2xxm,tau2xxm,m2xxm,theta2xxm)-t(Htxxm))
        verr2xxm <- sum(dif2xxm**2)
        verr2xxmlow <- sum(dif2xxm[1:7]**2)
        verr2xxmup <- sum(dif2xxm[8:length(dif2xxm)]**2)

        dif2xxxm <- 10*(cdfGLN4(vp2m[1:(length(vp2m)-1)],exp(omega2xxxm),omega2xxxm,tau2xxxm,m2xxxm,theta2xxxm)-t(Htxxxm))
        verr2xxxm <- sum(dif2xxxm**2)
        verr2xxxmlow <- sum(dif2xxxm[1:7]**2)
        verr2xxxmup <- sum(dif2xxxm[8:length(dif2xxxm)]**2)

        CorrMoment <- Yvcorrelations(sx, sxx, sxxx, p1, p2, p3, hg, upi,upii,upiii,upiv,upv, vp1, vp2, Corr)
        F <- 10*((eq1mslow + eq1msup)+(eq2mslow + eq2msup) +(eq2slow + eq2sup) +
                   (verr1xlow + verr1xup+verr1xmlow + verr1xmup +verr1xxlow + verr1xxup +verr1xxmlow + verr1xxmup
                    +verr1xxxlow + verr1xxxup +verr1xxxmlow + verr1xxxmup +verr2xlow + verr2xup + verr2xmlow + verr2xmup
                    + verr2xxlow + verr2xxup +verr2xxmlow + verr2xxmup +verr2xxxlow + verr2xxxup +verr2xxxmlow + verr2xxxmup +CorrMoment))
      }else{
        F <- 1000000000000000000000000
      }}else{
        F <- 1000000000000000000000000
      }
  }

  opts<-list(abstol=1e-4,reltol=1e-4,maxit=1)
  
  iter <-50
  for (i in 1:iter){
    hat <- optim(c(u01, u02, u03, u04, u05, v02,v01m,v02m, pr0[1,1:(length(pr0[1,])-1)], pr0[2,1:(length(pr0[2,])-1)], pr0[3,1:(length(pr0[3,])-1)], zeta0), Q,control=opts)
    phat<-hat$par
    Ehat<-hat$value
    u01 <- c(phat[1], phat[2], phat[3], phat[4], 0)
    u02 <- c(phat[6], phat[7], phat[8], phat[9], 0)
    u03 <- c(phat[11], phat[12], phat[13], phat[14], 0)
    u04 <- c(phat[16], phat[17], phat[18], phat[19], 0)
    u05 <- c(phat[21], phat[22], phat[23], phat[24], 0)
    v02 <- phat[26:39]
    v01m <- phat[40:53]
    v02m <- phat[54:67]
    
    pr0_nolog<-matrix(nrow=3,ncol=5)
    pr0[1,1:4] <- phat[68:71]
    pr0_nolog[1,1:4] <- exp(phat[68:71])
    pr0_nolog[1,5] <- 1 - sum(pr0_nolog[1,1:4])
    pr0[1,5] <- log(pr0_nolog[1,5])

    pr0[2,1:4] <- phat[72:75]
    pr0_nolog[2,1:4] <- exp(phat[72:75])
    pr0_nolog[2,5] <- 1 - sum(pr0_nolog[2,1:4])
    pr0[2,5] <- log(pr0_nolog[2,5])

    pr0[3,1:4] <- phat[76:79]
    pr0_nolog[3,1:4] <- exp(phat[76:79])
    pr0_nolog[3,5] <- 1 - sum(pr0_nolog[3,1:4])
    pr0[3,5] <- log(pr0_nolog[3,5])

    zeta0 <- phat[80]


    pr0test <- exp(pr0)
    sum(exp(t(pr0)))
    iter<-iter-1
    print(iter)
  }
  
  hat <- optim(c(u01, u02, u03, u04, u05, v02,v01m,v02m, pr0[1,1:(length(pr0[1,])-1)], pr0[2,(1:length(pr0[2,])-1)], pr0[3,(1:length(pr0[3,])-1)], zeta0), Q,control=opts)
  phat<-hat$par
  Ehat<-hat$value
  
  upih<-list()
  upiih<-list()
  upiiih<-list()
  upivh<-list()
  upvh<-list()
  
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
  v1hm <- exp(phat[40:53])
  v2hm <- exp(phat[54:67])
  
  prh<-matrix(nrow=3,ncol=5)
  prh[1,1:4] <- exp(phat[68:71])
  prh[1,5] <- 1-sum(prh[1,1:4])

  prh[2,1:4] <- exp(phat[72:75])
  prh[2,5] <- 1-sum(prh[2,1:4])

  prh[3,1:4] <- exp(phat[76:79])
  prh[3,5] <- 1-sum(prh[3,1:4])

  zetah <- exp(phat[80])

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
  yimidbin <- (0+yih[1])/2
  yimidbin[2:length(yih)] <- (yih[1:(length(yih)-1)]+yih[2:length(yih)])/2

  yiih <- Ycutoff(hg,upiih,v1h)
  yiimidbin <- (0+yiih[1])/2
  yiimidbin[2:length(yiih)] <- (yiih[1:(length(yiih)-1)]+yiih[2:length(yiih)])/2

  yiiih <- Ycutoff(hg,upiiih,v1h)
  yiiimidbin <- (0+yiiih[1])/2
  yiiimidbin[2:length(yiiih)] <- (yiiih[1:(length(yiiih)-1)]+yiiih[2:length(yiiih)])/2

  yivh <- Ycutoff(hg,upivh,v1h)
  yivmidbin <- (0+yivh[1])/2
  yivmidbin[2:length(yivh)] <- (yivh[1:(length(yivh)-1)]+yivh[2:length(yivh)])/2

  yvh <- Ycutoff(hg,upvh,v1h)
  yvmidbin <- (0+yvh[1])/2
  yvmidbin[2:length(yvh)] <- (yvh[1:(length(yvh)-1)]+yvh[2:length(yvh)])/2


  ytih <- Ycutoff(hg,upih,v2h)
  ytimidbin <- (0+ytih[1])/2
  ytimidbin[2:length(ytih)] <- (ytih[1:(length(ytih)-1)]+ytih[2:length(ytih)])/2

  ytiih <- Ycutoff(hg,upiih,v2h)
  ytiimidbin <- (0+ytiih[1])/2
  ytiimidbin[2:length(ytiih)] <- (ytiih[1:(length(ytiih)-1)]+ytiih[2:length(ytiih)])/2

  ytiiih <- Ycutoff(hg,upiiih,v2h)
  ytiiimidbin <- (0+ytiiih[1])/2
  ytiiimidbin[2:length(ytiiih)] <- (ytiiih[1:(length(ytiiih)-1)]+ytiiih[2:length(ytiiih)])/2

  ytivh <- Ycutoff(hg,upivh,v2h)
  ytivmidbin <- (0+ytivh[1])/2
  ytivmidbin[2:length(ytivh)] <- (ytivh[1:(length(ytivh)-1)]+ytivh[2:length(ytivh)])/2

  ytvh <- Ycutoff(hg,upvh,v2h)
  ytvmidbin <- (0+ytvh[1])/2
  ytvmidbin[2:length(ytvh)] <- (ytvh[1:(length(ytvh)-1)]+ytvh[2:length(ytvh)])/2


  shi <- v1h[1:(length(v1h)-1)]/yimidbin
  shii <- v1h[1:(length(v1h)-1)]/yiimidbin
  shiii <- v1h[1:(length(v1h)-1)]/yiiimidbin
  shiv <- v1h[1:(length(v1h)-1)]/yivmidbin
  shv <- v1h[1:(length(v1h)-1)]/yvmidbin

  shti <- v2h[1:(length(v2h)-1)]/ytimidbin
  shtii <- v2h[1:(length(v2h)-1)]/ytiimidbin
  shtiii <- v2h[1:(length(v2h)-1)]/ytiiimidbin
  shtiv <- v2h[1:(length(v2h)-1)]/ytivmidbin
  shtv <- v2h[1:(length(v2h)-1)]/ytvmidbin

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
  
  yimidbinm<-c()
  yiimidbinm<-c()
  yiiimidbinm<-c()
  yivmidbinm<-c()
  yvmidbinm<-c()
  ytimidbinm<-c()
  ytiimidbinm<-c()
  ytiiimidbinm<-c()
  ytivmidbinm<-c()
  ytvmidbinm<-c()
  
  
  yihm <- Ycutoff(hg,upih,v1hm)
  yimidbinm[1] <- (0+yihm[1])/2
  yimidbinm[2:length(yihm)] <- (yihm[1:(length(yihm)-1)]+yihm[2:length(yihm)])/2

  yiihm <- Ycutoff(hg,upiih,v1hm)
  yiimidbinm[1] <- (0+yiihm[1])/2
  yiimidbinm[2:length(yiihm)] <- (yiihm[1:(length(yiihm)-1)]+yiihm[2:length(yiihm)])/2

  yiiihm <- Ycutoff(hg,upiiih,v1hm)
  yiiimidbinm[1] <- (0+yiiihm[1])/2
  yiiimidbinm[2:length(yiiihm)] <- (yiiihm[1:(length(yiiihm)-1)]+yiiihm[2:length(yiiihm)])/2

  yivhm <- Ycutoff(hg,upivh,v1hm)
  yivmidbinm[1] <- (0+yivhm[1])/2
  yivmidbinm[2:length(yivhm)] <- (yivhm[1:(length(yivhm)-1)]+yivhm[2:length(yivhm)])/2

  yvhm <- Ycutoff(hg,upvh,v1hm)
  yvmidbinm[1] <- (0+yvhm[1])/2
  yvmidbinm[2:length(yvhm)] <- (yvhm[1:(length(yvhm)-1)]+yvhm[2:length(yvhm)])/2


  ytihm <- Ycutoff(hg,upih,v2hm)
  ytimidbinm[1] <- (0+ytihm[1])/2
  ytimidbinm[2:length(ytihm)] <- (ytihm[1:(length(ytihm)-1)]+ytihm[2:length(ytihm)])/2

  ytiihm <- Ycutoff(hg,upiih,v2hm)
  ytiimidbinm[1] <- (0+ytiihm[1])/2
  ytiimidbinm[2:length(ytiihm)] <- (ytiihm[1:(length(ytiihm)-1)]+ytiihm[2:length(ytiihm)])/2

  ytiiihm <- Ycutoff(hg,upiiih,v2hm)
  ytiiimidbinm[1] <- (0+ytiiihm[1])/2
  ytiiimidbinm[2:length(ytiiihm)] <- (ytiiihm[1:(length(ytiiihm)-1)]+ytiiihm[2:length(ytiiihm)])/2

  ytivhm <- Ycutoff(hg,upivh,v2hm)
  ytivmidbinm[1] <- (0+ytivhm[1])/2
  ytivmidbinm[2:length(ytivhm)] <- (ytivhm[1:(length(ytivhm)-1)]+ytivhm[2:length(ytivhm)])/2

  ytvhm <- Ycutoff(hg,upvh,v2hm)
  ytvmidbinm[1] <- (0+ytvhm[1])/2
  ytvmidbinm[2:length(ytvhm)] <- (ytvhm[1:(length(ytvhm)-1)]+ytvhm[2:length(ytvhm)])/2


  shim <- v1hm[1:(length(v1hm)-1)]/yimidbinm
  shiim <- v1hm[1:(length(v1hm)-1)]/yiimidbinm
  shiiim <- v1hm[1:(length(v1hm)-1)]/yiiimidbinm
  shivm <- v1hm[1:(length(v1hm)-1)]/yivmidbinm
  shvm <- v1hm[1:(length(v1hm)-1)]/yvmidbinm

  shtim <- v2hm[1:(length(v2hm)-1)]/ytimidbinm
  shtiim <- v2hm[1:(length(v2hm)-1)]/ytiimidbinm
  shtiiim <- v2hm[1:(length(v2hm)-1)]/ytiiimidbinm
  shtivm <- v2hm[1:(length(v2hm)-1)]/ytivmidbinm
  shtvm <- v2hm[1:(length(v2hm)-1)]/ytvmidbinm

  ehy1im <- IncomeElast(hg,v1hm,upih)
  ehy1iim <- IncomeElast(hg,v1hm,upiih)
  ehy1iiim <- IncomeElast(hg,v1hm,upiiih)
  ehy1ivm <- IncomeElast(hg,v1hm,upivh)
  ehy1vm <- IncomeElast(hg,v1hm,upvh)

  ehy2im <- IncomeElast(hg,v2hm, upih)
  ehy2iim <- IncomeElast(hg,v2hm,upiih)
  ehy2iiim <- IncomeElast(hg,v2hm,upiiih)
  ehy2ivm <- IncomeElast(hg,v2hm,upivh)
  ehy2vm <- IncomeElast(hg,v2hm,upvh)

  ycv<-c()
  CV<-c()
  fvalCV<-c()
  for (i in 1:length(Ycut)){
    CVi<-CalcCV(Ycut[i],upih,v2h,v2hm, hg)
    ycv<-append(ycv,CVi$yphat)
    CV<-append(CV,CVi$cv)
    fvalCV<-append(fvalCV,CVi$fval)
  }

  PIh<-list()

  PIh$upi <- upih
  PIh$upii <- upiih
  PIh$upiii <- upiiih
  PIh$upiv <- upivh
  PIh$upv <- upvh

  PIh$v2 <- v2h
  PIh$v1m <- v1hm
  PIh$v2m <- v2hm

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

  PIh$elast$ehy1im <- ehy1im
  PIh$elast$ehy1iim <- ehy1iim
  PIh$elast$ehy1iiim <- ehy1iiim
  PIh$elast$ehy1ivm <- ehy1ivm
  PIh$elast$ehy1vm <- ehy1vm

  PIh$elast$ehy2im <- ehy2im
  PIh$elast$ehy2iim <- ehy2iim
  PIh$elast$ehy2iiim <- ehy2iiim
  PIh$elast$ehy2ivm <- ehy2ivm
  PIh$elast$ehy2vm <- ehy2vm

  PIh$elast$ehp1i <- ehp1i
  PIh$elast$ehp1ii <- ehp1ii
  PIh$elast$ehp1iii <- ehp1iii
  PIh$elast$ehp1iv <- ehp1iv
  PIh$elast$ehp1v <- ehp1v

  PIh$CV$ycv<-ycv[1:length(ycv)-1]
  PIh$CV$CV<-CV[1:length(CV)-1]

  return(list('PIh'=PIh,'Ehat'=Ehat))
}

# This function loads the data to estimate the model in Dennis, Quintero and Sieg (2019) with 1 metro area.

a3loadinitialdata_2ma <- function(WPmetro, WPmetrom, WPAgr) {
  YIk <- WPmetro$YIk
  VIk <- WPmetro$VIk
  PIk <- WPmetro$PIk
  Pop <- WPmetro$Pop

  YIkm <- WPmetrom$YIk
  VIkm <- WPmetrom$Vik
  PIkm <- WPmetrom$PIk
  Popm <- WPmetrom$Pop

  Corra <- WPAgr$Corr

  return(list(YIk = YIk, VIk = VIk, PIk = PIk, Pop = Pop,
              YIkm = YIkm, VIkm = VIkm, PIkm = PIkm, Popm = Popm, Corra = Corra))
}

# Calculate the compensating variation required to keep utility that guarantees the same utility for a household
# that moves from the reference m metro area to the base metro area.
# yphat is the resulting total income and cv is the compensating variation (the additional amount)
library(neldermead)
CalcCV <- function(y,uph,v,vm,hg){

  opts<-list(abstol=1e-4,reltol=1e-4,maxit=5000)

  test<- caluy(y,uph,v,hg)
  nomes<-c('utest','htest','vtest')
  for (i in 1:length(names(test))){
    assign(nomes[i],unlist(test[[names(test)[i]]]))
  }

  testm<- caluy(y,uph,vm,hg)
  nomes<-c('um','hm','vhm')
  for (i in 1:length(names(testm))){
    assign(nomes[i],unlist(testm[[names(testm)[i]]]))
  }

  CV <- function(yp){

    testh<- caluy(exp(yp),uph,v,hg)
    nomes<-c('uh','hh','vh')
    for (i in 1:length(names(testh))){
      assign(nomes[i],unlist(testh[[names(testh)[i]]]))
    }

    if (is.numeric(um)==0 & is.numeric(uh)==0){
      F <- (um-uh)**2
    } else {
      F <- 1000000000000000
    }
    return(F)
  }

  ypini <- y
  yp <- optim(log(ypini),CV,method='SANN',control=opts)
  yp_hat<-yp$par
  fval<-yp$value
  yphat <- exp(yp_hat)
  cv <- yphat - y

  phi <- uph$phi
  gamma <- uph$gamma
  alpha <- uph$alpha
  eta <- uph$eta
  kappa <- uph$kappa

  caly<- caluy(yphat,uph,v,hg)
  nomes<-c('uhat','hhat','vhat')
  for (i in 1:length(names(caly))){
    assign(nomes[i],unlist(caly[[names(caly)[i]]]))
  }
  return(list('yphat'=yphat, 'cv'=cv,'fval'=fval))
}

# This calculates the preference parameters implied by a single type model.
# See appendix in Dennis, Quintero and Sieg (2019).

CalLinearPref <- function(Py, Pv) {
  # INCOME t=1
  mu <- Py[1]
  sigma <- Py[2]
  beta <- Py[4]

  # VALUE t=1
  omega <- Pv[1]
  tau <- Pv[2]
  theta <- Pv[4]

  # We are changing the scale of theta and beta everywhere
  theta <- 1000 * theta
  beta <- 1000 * beta

  # PREFERENCES
  b1 <- sigma / tau
  a1 <- mu - (b1) * omega
  A1 <- exp(a1)

  # Create U
  U <- list()
  U$alpha <- 1 / (b1 - 1)
  U$phi <- 1 / A1
  U$eta <- theta
  U$gamma <- 1 - b1

  return(U)
}

# Calculate the quality h consumed at income y at prices v in quality grid hg
caluy <- function(y,uph,v,hg){

  # What is the quality that a household with income y would consume?
  yh <- Ycutoff(hg,uph,v)

  diff_squared <- (yh - y)**2
  b <- which.min(diff_squared)
  a <- min(diff_squared)
  h <- hg[b]
  vh <- v[b]

  #What is the utility of this quality at this rent price and this income
  uph<-list()
  phi <- uph$phi
  gamma <- uph$gamma
  alpha <- uph$alpha
  eta <- uph$eta
  kappa <- uph$kappa

  u <- utility(y,h,phi,gamma,alpha,eta,kappa,vh);

  return(list('u'=u,'h'=h,'vh'=vh))
}

# This function evaluates the GLN4 cdf distribution. See appendix of Dennis, Quintero and Sieg (2019).

cdfGLN4 <- function(x, mediancdforig, mu, sigma, r, beta) {
  # We are changing the scale of theta and beta everywhere
  beta <- 1000 * beta

  # Initialize variables
  uu <- mu
  ss <- sigma
  rr <- r
  lengthx <- length(x)
  v <- 1 / r
  menor <- 0
  mayor <- 0

  # Create vectors to store the results
  B <- numeric(lengthx)
  T <- numeric(lengthx)
  F <- numeric(lengthx)
  M <- numeric(lengthx)
  
  
  # Loop through the input vector
  for (i in 1:lengthx) {
    # If x[i] + beta is less than mediancdforig, the cdf is less than 0.5
    if (x[i] + beta < mediancdforig) {
      menor <- 1
      B[i] <- (((mu - log(x[i] + beta)) / sigma)**r) / r
      T[i] <- B[i]
      F[i] <- gamma(v) * pgamma(B[i], shape = v, lower.tail = FALSE) / (2 * gamma(v))

    }

    # If x[i] + beta is equal to mediancdforig, the cdf is equal to 0.5
    if (x[i] + beta == mediancdforig) {
      med <- 1
      F[i] <- 0.5
    }

    # If x[i] + beta is greater than mediancdforig, the cdf is greater than 0.5
    if (x[i] + beta > mediancdforig) {
      mayor <- 1
      M[i] <- ((log(x[i] + beta) - mu) / sigma)**r / r
      T[i] <- M[i]
      F[i] <- 0.5 + gamma(v) * pgamma(M[i], shape = v) / (2 * gamma(1/r))
    }

  }

  # Return the cdf vector

  return(F)
}

deriv_dy_dh_FOC <- function(h,alpha,phi,gammap,eta,vp,vpp){
  #     vp = Vprime_GLN4(h,mu,sigma,tau,omega,phi,gammap,alpha,eta);
  #     vpp = Vdoubleprime_GLN4(h,mu,sigma,tau,omega,phi,gammap,alpha,eta);
  heta <- h+eta
  T <- vp*(1-phi*((heta)**gammap))
  Tp <- vpp*(1-phi*(heta**gammap)) + vp*(-phi*gammap*(heta**(gammap-1)))
  B <- -phi*gammap*alpha*(heta**(gammap-1))
  Bp <- -phi*gammap*alpha*(gammap-1)*(heta**(gammap-2))
  Dydh <- ((Tp*B - Bp*T)/(B**2)) + vp

  return(Dydh)
}

#Função que determina a FOC quando se assume que a quação de precificação é linear
#Variáveis: vetor de qualidades (h),vetor de rendas (y), parâmetros da função de utilidade (alpha,phi,gammap,eta,kappa)

dp_dh_FOC_linearVh<-function(h,eta,alpha,phi,gammap,y,kappa){
  heta = (h + eta)
  S = -alpha*phi*gammap*h + (heta)**(1-gammap) - phi*heta
  Sp = -alpha*phi*gammap + (1-gammap)*(heta)**(-gammap) - phi
  M = -alpha*phi*gammap*(y-kappa)
  return (-M*Sp/(S**2))
}

# Derivative dy_dh from the FOCs, when pricing function is linear
deriv_dy_dh_FOC_linearVh <- function(h,p,alpha,phi,gammap,eta){
  S <- -p/(alpha*phi*gammap)
  M <- (h+eta)**(-gammap)
  L <- -alpha*phi*gammap -phi + (1 - gammap) * M
  return(S*L)

}

# This function estimates the quantile of the GLN4 distribution. See appendix of Dennis, Quintero and Sieg (2019).

GLN4quantiles <- function(p, mediancdforig, mu, sigma, r, beta) {

  # Set the optimization options
  opts<-list(abstol=1e-4,reltol=1e-4,maxit=5000)
  
  # Define the objective function
  GLN4quant <- function(q, beta) {
    if (-1000 * beta < exp(q)) {
      F <- (cdfGLN4(exp(q), mediancdforig, mu, sigma, r, beta) - p)**2
    } else {
      F <- 1e15  # A large value to discourage q values that are too small
    }
    return(F)
  }

  qini <- max(-1000 * beta + 100, 5)
  # tive que usar outra fun??o de minimizacao que nao tem aqueles controles
  result <- optim(par = log(qini), fn = GLN4quant, beta = beta, method = "L-BFGS-B", control = opts)
  q_hat <- exp(result$par)
  fval <- result$value

  return(list('q_hat'=q_hat,'fval'=fval))
}

#Função do modelo que determina que aluguéis com base na distribuição GLN4

GLN4u<-function(Y,V,P0,Solver,metro){
  final<-list(P=list(),E=list())
  opts<-list(abstol=1e-10,reltol=1e-10,maxit=1e9)

  #Renda t=1

  GLNobjy<-function(paramy){
    return(sum(cdfGLN4(Y$Y,exp(paramy[1]),paramy[1],exp(paramy[2]),exp(paramy[3]),paramy[4])-t(Y$Ycd))**2)
  }

  y<-optim(c(P0$y[1],log(P0$y[2]),log(PO$y[3]),P0$y[4]),GLNobjy, control=opts)
  final$P$y<-y$par
  final$E$y<-y$value

  mu1=final$P$y[1]
  sigma1=exp(final$P$y[2])
  final$P$y[2]=exp(final$P$y[2])
  r1=exp(final$P$y[3])
  final$P$y[3]=exp(final$P$y[3])
  beta1=final$P$y[4]

  GLNy<-cdfGLN4(Y$Y,exp(mu1),mu1,sigma1,r1,beta1)
  # Create a figure
  h <- fig_window()

  # Plot the data
  plot(log(Y$Y), GLNy)

  # Add a scatter plot
  points(log(Y$Y), Y$Ycd, col = "black", pch = 20, cex = 0.8)

  # Add axis labels and a title
  xlabel("log(Y_1)")
  title(sprintf("Income, Period 1. GLN4 Unconstrained fit."))

  # Set the font size of all text elements in the figure
  text_obj <- find_objects(h, type = "text")
  for (i in 1:length(text_obj)) {
    set_text_properties(text_obj[[i]], cex = 20)
  }

  GLNobjv<-function(paramv){
    return(sum((cdfGLN4(V$V,exp(paramv[1]),paramv[1],exp(paramv[2]),exp(paramv[3]),paramv[4])-t(V$Vcd))**2))
  }

  #Funções para encontrar mínimo em cada coluna de matriz e seus índices
  

  GLNobjvPh<-function(paramv){
    if(-1000*paramv[4]<min_value(cbind(t(V$V),t(V$Vt)))){
      return (sum((cdfGLN4(V$V,exp(paramv[1]),paramv[1],exp(paramv[2]),exp(paramv[3]),paramv[4]-t(V$Vcd))**2)))
    }
    else{
      return (100000000)
    }
  }

  v<-optim(c(P0$v[1],log(P0$v[2]),log(PO$v[3]),P0$v[4]),GLNobjv, control=opts)
  final$P$v<-V$par
  final$E$v<-V$value

  omega1<-final$P$v[1]
  tau1<-exp(final$P$v[2])
  final$P$v[2]<-exp(final$P$v[2])
  m1<-exp(final$P$v[3])
  final$P$v[3]<-exp(final$P$v[3])
  theta1<-final$P$v(4)

  GLN4v = cdfGLN4(V$V,exp(omega1),omega1,tau1,m1,theta1)
  # Create a figure
  h <- fig_window()

  # Plot the data
  plot(log(V$V), GLN4v)

  # Add a scatter plot
  points(log(V$V), V$Vcd, col = "black", pch = 20, cex = 0.8)

  # Add axis labels and a title
  xlabel("log(V_1)")
  title(sprintf("Value, Period 1. GLN4 Unconstrained fit."))

  # Set the font size of all text elements in the figure
  text_obj <- find_objects(h, type = "text")
  for (i in 1:length(text_obj)) {
    set_text_properties(text_obj[[i]], cex = 20)
  }
  #Renda t=2

  GLNobjyt<-function(paramyt){
    return(sum(cdfGLN4(Y$Yt,exp(paramyt[1]),paramyt[1],exp(paramyt[2]),exp(paramyt[3]),paramyt[4])-t(Y$Ytcd))**2)
  }

  yt<-optim(c(P0$yt[1],log(P0$yt[2]),log(PO$yt[3]),P0$yt[4]),GLNobjyt, control=opts)
  final$P$yt<-yt$par
  final$E$yt<-yt$value

  mu2=final$P$yt[1]
  sigma2=exp(final$P$yt[2])
  final$P$yt[2]=exp(final$P$yt[2])
  r2=exp(final$P$yt[3])
  final$P$yt[3]=exp(final$P$yt[3])
  beta2=final$P$yt[4]

  GLNyt<-cdfGLN4(Y$Yt,exp(mu2),mu2,sigma2,r2,beta2)
  # Create a figure
  h <- fig_window()

  # Plot the data
  plot(log(Y$Yt), GLNyt)

  # Add a scatter plot
  points(log(Y$Yt), Y$Ytcd, col = "black", pch = 20, cex = 0.8)

  # Add axis labels and a title
  xlabel("log(Y_2)")
  title(sprintf("Income, Period 2. GLN4 Unconstrained fit."))

  # Set the font size of all text elements in the figure
  text_obj <- find_objects(h, type = "text")
  for (i in 1:length(text_obj)) {
    set_text_properties(text_obj[[i]], cex = 20)
  }

  GLNobjvt<-function(paramvt){
    return(sum((cdfGLN4(V$Vt,exp(paramvt[1]),paramvt[1],exp(paramvt[2]),exp(paramvt[3]),paramvt[4])-t(V$Vtcd))**2))
  }

  if (metro=='NewYork'){
    vt<-optim(c(P0$vt[1],log(P0$vt[2]),log(PO$vt[3]),P0$vt[4]),GLNobjvtPh, control=opts)
  }
  else {
    vt<-optim(c(P0$vt[1],log(P0$vt[2]),log(PO$vt[3]),P0$vt[4]),GLNobjvt, control=opts)
  }
  final$P$vt<-vt$par
  final$E$vt<-vt$value

  omega2<-final$P$vt[1]
  tau2<-exp(final$P$vt[2])
  final$P$vt[2]<-exp(final$P$vt[2])
  m2<-exp(final$P$vt[3])
  final$P$vt[3]<-exp(final$P$vt[3])
  theta2<-final$P$vt(4)

  GLN4vt = cdfGLN4(V$Vt,exp(omega2),omega2,tau2,m2,theta2)
  GLN4v = cdfGLN4(V$V,exp(omega1),omega1,tau1,m1,theta1)
  # Create a figure
  h <- fig_window()

  # Plot the data
  plot(log(V$Vt), GLNvt)

  # Add a scatter plot
  points(log(V$Vt), V$Vtcd, col = "black", pch = 20, cex = 0.8)

  # Add axis labels and a title
  xlabel("log(V_2)")
  title(sprintf("Value, Period 2. GLN4 Unconstrained fit."))

  # Set the font size of all text elements in the figure
  text_obj <- find_objects(h, type = "text")
  for (i in 1:length(text_obj)) {
    set_text_properties(text_obj[[i]], cex = 20)
  }

  return (final)
}

# Calculate the
hd <-function(pxy,pxxy,pxxxy,si,sii,siii,pi,ycf){
  #Calculate the demand for each quality j = 1:J implied by the income
  #cutoffs ycf. ycf gives the cutoff points that make households
  #indifferent between quality j and j-1 given qualities J and prices v
  hd <- c()
  mux <- pxy[1]
  sigmax <- pxy[2]
  rx <- pxy[3]
  betax <- pxy[4]
  muxx <- pxxy[1]
  sigmaxx <- pxxy[2]
  rxx <- pxxy[3]
  betaxx <- pxxy[4]
  muxxx <- pxxxy[1]
  sigmaxxx <- pxxxy[2]
  rxxx <- pxxxy[3]
  betaxxx <- pxxxy[4]

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
  F <- matrix(c(GLNyxj, GLNyxxj, GLNyxxxj), ncol = length(ycf) - 1, byrow = TRUE)
  D <- sum(s*pi) #produto escalar, correto
  N <- (s*pi)%*%F #o termo entre parenteses é a diagonal da matriz, um vetor 1x3, que deve ser multiplicado por um 3x12
  Fi <-  N/D
  Fi <- c(Fi)

  # lagged
  F_1 <- matrix(c(GLNyxj_1, GLNyxxj_1, GLNyxxxj_1), ncol = length(ycf) - 1, byrow = TRUE)
  N_1 <- (s*pi)%*%F_1
  Fi_1 <-  N_1/D
  # demands are the difference between the two cumulative income
  # distributions
  hd[1] <-  Fi[1]
  hd[2:length(ycf)] <-  Fi - Fi_1

  return(hd)
}




#Calculate Supplies and Demand for each quality level for the 1 metro area
#model or for the base metro area in the multiple metro areas model
#Calculate Supplies and Demand for each quality level for the 1 metro area
#model or for the base metro area in the multiple metro areas model
Hinv3x5INorm1 <- function(h,vp2,pxy,pxxy,pxxxy, pxyt, pxxyt, pxxxyt, upi,upii,upiii,upiv,upv,sx,sxx,sxxx,p1,p2,p3,zeta,VI){
  
  Hdtcdi = c()
  Hdtcdii = c() 
  Hdtcdiii = c() 
  Hdtcdiv = c() 
  Hdtcdv = c()
  Htx  = c() 
  Htxx = c() 
  Htxxx = c()
  AgHdtcdf = c() 
  Rh2= c() 
  #N2 = c() 
  #grateh = c()  
  
  # Rh from the normalization that sets v(h) = v
  ## Calculate supply and demand in t=1
  Rhlinear <-Rh_3x5Ilinearvh(h,pxy,pxxy,pxxxy,upi,upii,upiii,upiv,upv,sx,sxx,sxxx, p1, p2, p3)
  nomes<-c('Rh','Hdcdi','Hdcdii','Hdcdiii','Hdcdiv','Hdcdv','Hx','Hxx','Hxxx','AgHdcdf','cond1')
  for (i in 1:length(names(Rhlinear))){
    assign(nomes[i],unlist(Rhlinear[[names(Rhlinear)[i]]]))
  }
  if (cond1 ==1){
    ## Check monotonicity
    monov <- monotonic(vp2,direction = 'inc')
    
    if  (monov) {
      ## Demand for period 2
      Rh_3 <-Rh_3x5Ivh(h,pxyt,pxxyt,pxxxyt,upi,upii,upiii,upiv,upv,sx,sxx,sxxx, p1, p2, p3, vp2)
      nomes<-c('Hdtcdi','Hdtcdii','Hdtcdiii','Hdtcdiv','Hdtcdv','Htx', 'Htxx', 'Htxxx','AgHdtcdf','cond2')
      for (i in 1:length(names(Rh_3))){
        assign(nomes[i],unlist(Rh_3[[names(Rh_3)[i]]]))
      }

      if (cond2 ==1){
        ## Supply for period 2
        
        ## of periods between observations
        lengthper <- 4
        
        VI$k1ypolyf <- VI$uf
        VI$k2ypolyf <- VI$uft 
        
        # linear function for price of first period
        #Joao
        #aqui entrava t(h) e o resultado era esquisito (um vetor saía transposto) em RhSl2PeriodsI
        vp1 <- t(h)
        Rhsl2P <-RhSl2PeriodsI(h,vp1,vp2,VI,zeta,AgHdcdf,lengthper)
        nomes<-c('Rh2', 'N2', 'grateh')
        for (i in 1:length(names(Rhsl2P))){
          assign(nomes[i],unlist(Rhsl2P[[names(Rhsl2P)[i]]]))
        }   
      } else {
        Rh2<- c()
        N2<- c() 
        grateh<- c()  
      }      
    } else { 
      
      Hdtcdi = c()
      Hdtcdii = c() 
      Hdtcdiii = c() 
      Hdtcdiv = c() 
      Hdtcdv = c()
      Htx  = c() 
      Htxx = c() 
      Htxxx = c()
      AgHdtcdf = c() 
      Rh2= c() 
      N2 = c() 
      grateh = c()  
      cond2=0
    }
    
  } else {
    
    Hdtcdi = c()
    Hdtcdii = c() 
    Hdtcdiii = c() 
    Hdtcdiv = c() 
    Hdtcdv = c()
    Htx  = c() 
    Htxx = c() 
    Htxxx = c()
    AgHdtcdf = c() 
    Rh2= c() 
    N2 = c() 
    grateh = c()  
    cond2=0
  }
  
  return(list('Rh'=Rh,'cond1'=cond1,'cond2'=cond2,'Hdcdi'=Hdcdi,'Hdcdii'=Hdcdii,'Hdcdiii'=Hdcdiii,'Hdcdiv'=Hdcdiv,'Hdcdv'=Hdcdv,'Hx'=Hx,'Hxx'=Hxx,'Hxxx'=Hxxx,'AgHdcdf'=AgHdcdf,
              'Hdtci'=Hdtcdi,'Hdtcii'=Hdtcdii,'Hdtciii'=Hdtcdiii,'Hdtciv'=Hdtcdiv,'Hdtcv'=Hdtcdv,'Htx'=Htx, 'Htxx'=Htxx,'Htxxx'=Htxxx,'AgHdtcdf'=AgHdtcdf, 'Rh2'=Rh2))
}


#Calculate Supplies and Demand for each quality level for the 2 metro area
# for the second or reference metro area
Hinv3x5INorm1m <- function(h,vp2m, pxym, pxxym, pxxxym, pxytm, pxxytm, pxxxytm, upi,upii,upiii,upiv,upv,sxm,sxxm,sxxxm,p1,p2,p3,zeta,VIm,vp1m,VI,Rh){

  # Calculate supply and demand

  ## Check monotonicity
  monov1 <- monotonic(vp1m, direction = 'inc')

  if  (monov1){
    ## Demand for period 1
    Rh_3 <- Rh_3x5Ivh(h,pxym,pxxym,pxxxym,upi,upii,upiii,upiv,upv,sxm,sxxm,sxxxm, p1, p2, p3, vp1m)
    nomes<-c('Hdcdim','Hdcdiim','Hdcdiiim','Hdcdivm','Hdcdvm','Hxm', 'Hxxm', 'Hxxxm','AgHdcdfm','cond1m')
    for (i in 1:length(names(Rh_3))){
      assign(nomes[i],unlist(Rh_3[[names(Rh_3)[i]]]))
    }
    
    if (cond1m ==1){

      ## Check monotonicity
      monov2 <- monotonic(vp2m, direction = 'inc')

      if  (monov2){
        ## Demand for period 2

        Rh_3x <- Rh_3x5Ivh(h,pxytm,pxxytm,pxxxytm,upi,upii,upiii,upiv,upv,sxm,sxxm,sxxxm,p1, p2, p3, vp2m)
        nomes<-c('Hdtcdim','Hdtcdiim','Hdtcdiiim','Hdtcdivm','Hdtcdvm','Htxm', 'Htxxm', 'Htxxxm','AgHdtcdfm','cond2m')
        for (i in 1:length(names(Rh_3x))){
          assign(nomes[i],unlist(Rh_3x[[names(Rh_3x)[i]]]))
        }
        
        if (cond2m == 1){

          ## Supply for period 1
          ## of periods between observations
          lengthper <- 4
          #Base metro
          Vbase1<-list()
          Vbase1$k1ypolyf <- VI$uf
          Vbase1$V <- VI$V

          #Ref Metro
          Vref1<-list()
          Vref1$k2ypolyf <- VIm$uf
          Vref1$Vt <- VIm$V

          # linear function for price of first period in the
          # refence metro
          vbase <- t(h)
          vref <- vp2m
          basesupply <- Rh

          Rhsl1 <- RhSl2PeriodsIm(h,vbase,vref,Vref1,zeta,basesupply,lengthper,Vbase1)
          nomes<-c('Rh1m', 'N1m', 'gratehm')
          for (i in 1:length(names(Rhsl1))){
            assign(nomes[i],unlist(Rhsl1[[names(Rhsl1)[i]]]))
          }

          ## Supply for period 2
          #Base metro
          Vbase2<-list()
          Vbase2$k1ypolyf <- VIm$uf
          Vbase2$V <- VIm$V

          #Ref Metro
          Vref2<-list()
          Vref2$k2ypolyf <- VIm$uft
          Vref2$Vt <- VIm$Vt

          vbase <- vp1m
          vref <- vp2m
          basesupply <- Rh1m

          Rhsl2 <- RhSl2PeriodsIm(h,vbase,vref,Vref2,zeta,basesupply,lengthper,Vbase2)
          nomes<-c('Rh2m', 'N2m', 'gratehtm')
          for (i in 1:length(names(Rhsl2))){
            assign(nomes[i],unlist(Rhsl2[[names(Rhsl2)[i]]]))
          }
          

        } else {
          Rh1m <- c()
          N1m <- c()
          gratehm <- c()
          Rh2m <- c()
          N2m <- c()
          gratehtm <- c()
        }
      } else {

        Hdtcdim <- c()
        Hdtcdiim <- c()
        Hdtcdiiim <- c()
        Hdtcdivm <- c()
        Hdtcdvm <- c()
        Htxm  <- c()
        Htxxm <- c()
        Htxxxm <- c()
        AgHdtcdfm <- c()

        Rh1m <- c()
        N1m <- c()
        gratehm <- c()
        Rh2m <- c()
        N2m <- c()
        gratehtm <- c()

        cond2m <- 0
      }

    } else {

      Hdtcdim <- c()
      Hdtcdiim <- c()
      Hdtcdiiim <- c()
      Hdtcdivm <- c()
      Hdtcdvm <- c()
      Htxm  <- c()
      Htxxm <- c()
      Htxxxm <- c()
      AgHdtcdfm <- c()
      Rh1m <- c()
      N1m <- c()
      gratehm <- c()
      Rh2m <- c()
      N2m <- c()
      gratehtm <- c()
      cond2m <- 0
    }

  } else {

    Hdcdim <- c()
    Hdcdiim <- c()
    Hdcdiiim <- c()
    Hdcdivm <- c()
    Hdcdvm <- c()
    Hdtcdim <- c()
    Hdtcdiim <- c()
    Hdtcdiiim <- c()
    Hdtcdivm <- c()
    Hdtcdvm <- c()

    Hxm  <- c()
    Hxxm <- c()
    Hxxxm <- c()
    Htxm  <- c()
    Htxxm <- c()
    Htxxxm <- c()

    AgHdcdfm <- c()
    AgHdtcdfm <- c()

    Rh1m <- c()
    N1m <- c()
    gratehm <- c()
    Rh2m <- c()
    N2m <- c()
    gratehtm <- c()
    cond2m <- 0
    cond1m <- 0

  }
  return(list('Rh1m'=Rh1m,'cond1m'=cond1m,'cond2m'=cond2m,'Hdcdim'=Hdcdim,'Hdcdiim'=Hdcdiim,'Hdcdiiim'=Hdcdiiim,'Hdcdivm'=Hdcdivm,'Hdcdvm'=Hdcdvm,
              'Hxm'=Hxm,'Hxxm'=Hxxm,'Hxxxm'=Hxxxm,'AgHdcdfm'=AgHdcdfm,'Hdtcdim'=Hdtcdim,'Hdtcdiim'=Hdtcdiim,'Hdtcdiiim'=Hdtcdiiim,'Hdtcdivm'=Hdtcdivm,
              'Hdtcdvm'=Hdtcdvm,'Htxm'=Htxm,'Htxxm'=Htxxm,'Htxxxm'=Htxxxm,'AgHdtcdfm'=AgHdtcdfm,'Rh2m'=Rh2m))
}

#library pracma e signal
IncomeElast <- function(h, v, U){
  phi <- U$phi
  gammap <- U$gamma
  alpha <- U$alpha
  eta <- U$eta
  kappa  <-  U$kappa

  #polytool(h,v,d)
  d  <-  2
  p  <-  polyfit(h, v, d) #obtém os coeficientes do polinomio - é polyfit no R tbm

  v  <-  polyval(p,h); #calcula o polinomio para todos os valores de h - tbm é polyval

  #v'(h)
  pd  <-  polyder(p) #obtém os coeficientes da derivada do polinomio - tbm é polyder
  vp  <-  polyval(pd,h)  #calcula a derivada
  #v''(h) - mesma coisa para a segunda derivada.
  pdd  <-  polyder(pd)
  vpp  <-  polyval(pdd,h)



  Dy1h  <-  deriv_dy_dh_FOC(h,alpha,phi,gammap,eta,vp,vpp)
  Dhy1  <-  1/Dy1h
  y1  <-  y_FOC_Vh(h,alpha,phi,gammap,eta,kappa,v,vp)
  ehy  <-  Dhy1*(y1/h)
  return(ehy)
}

#library pracma e signal
IncomeLinearVElast <- function(h,U){
  phi <- U$phi
  gammap <- U$gamma
  alpha <- U$alpha
  eta <- U$eta
  kappa <- U$kappa

  p<-1

  Dy1h  <-  deriv_dy_dh_FOC_linearVh(h,p,alpha,phi,gammap,eta)
  Dhy1  <-  1/Dy1h
  y1  <-  y_FOC_linearVh(h,p,alpha,phi,gammap,eta,kappa)
  ehy  <-  Dhy1*(y1/h)
  return(ehy)
}

#Função que determina a FOC quando se assume que a equação de precificação é linear
#Variáveis: vetor de qualidades (h), vetor de rendas (ym), parâmetros da função de utilidade (alpha,phi,gammap,eta,kappa)

p_FOC_linearVh<-function(h,ym,alpha,phi,gammap,eta,kappa){
  heta = (h + eta)
  S = -alpha*phi*gammap*h + (heta)**(1-gammap) - phi*heta
  M = -alpha*phi*gammap*(ym-kappa)
  return (M/S)
}

#Calculate price quality consumed elasticity for given type when a linear pricing function is
#assumed
PriceLinearVElast <- function(h,U){

  #Recover Parameters
  phi <- U$phi
  gammap <- U$gamma
  alpha <- U$alph
  eta <- U$eta
  kappa  <-  U$kappa
  #Elasticities
  #Corresponding income
  v  <-  h
  y  <-  Ycutoff(h,U,v)
  # Period 1
  Dp1h  <-  dp_dh_FOC_linearVh(h[1:(length(h)-1)],eta,alpha,phi,gammap,y,kappa)
  Dhp1  <-  1/Dp1h
  p1  <-  p_FOC_linearVh(h[1:(length(h)-1)],y,alpha,phi,gammap,eta,kappa)
  ehp1  <-  Dhp1*(p1/h[1:(length(h)-1)])
  return(ehp1)
}

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

  budgeti <- ycfi-vp[1:(length(h)-1)]-upi$kappa
  budgetii <- ycfii-vp[1:(length(h)-1)]-upii$kappa
  budgetiii <- ycfiii-vp[1:(length(h)-1)]-upiii$kappa
  budgetiv <- ycfiv-vp[1:(length(h)-1)]-upiv$kappa
  budgetv <- ycfv-vp[1:(length(h)-1)]-upv$kappa
  
  if (realyi && realyii && realyiii && realyiv && realyv && monoyi && monoyii && monoyiii && monoyiv && monoyv && all(ycfi>0) && all(ycfii>0) && all(ycfiii>0) && all(ycfiv>0) && all(ycfv>0) && all(budgeti>0) && all(budgetii>0) && all(budgetiii>0) && all(budgetiv>0) && all(budgetv>0)){
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
    Hdcdi = c() 
    Hdcdii = c()
    Hdcdiii = c()
    Hdcdiv = c()
    Hdcdv = c() 
    
    Hdcdi[1] <- hdi[1]
    Hdcdii[1] <- hdii[1]
    Hdcdiii[1] <- hdiii[1]
    Hdcdiv[1] <- hdiv[1]
    Hdcdv[1] <- hdv[1]
    
    for (j in 2:(length(h)-1)){
      Hdcdi[j] <- Hdcdi[j-1]+hdi[j]
      Hdcdii[j] <- Hdcdii[j-1]+hdii[j]
      Hdcdiii[j] <- Hdcdiii[j-1]+hdiii[j]
      Hdcdiv[j] <- Hdcdiv[j-1]+hdiv[j]
      Hdcdv[j] <- Hdcdv[j-1]+hdv[j]
    }

    #Aggregate the CDFs (get aggregate demand distribution)
    H <- as.matrix(rbind(Hdcdi, Hdcdii, Hdcdiii, Hdcdiv, Hdcdv))
    s <- c(si, sii, siii)
    
    ax1 <- (as.matrix(s)%*%p1)
    bx1 <- ax1[1,]

    ax2 <- (as.matrix(s)%*%p2)
    bx2 <- ax2[2,]

    ax3 <- (as.matrix(s)%*%p3)
    bx3 <- ax3[3,]
    
    Hx <- as.vector(bx1%*%H)
    Hxx <- as.vector(bx2%*%H)
    Hxxx <- as.vector(bx3%*%H)
    
    AgHdcdf <- Hx + Hxx + Hxxx
    #also return the observed type demands to construct moments in the main
    #objective function

    cond=1
  }else{
    AgHdcdf <- c()
    Hdcdi <- c()
    Hdcdii <- c()
    Hdcdiii <- c()
    Hx  <- c()
    Hxx <- c()
    Hxxx <- c()
    cond=0
  }
  Rh <- AgHdcdf
  return(list('Rh' = Rh,'Hdcdi' = Hdcdi,'Hdcdii' = Hdcdii, 'Hdcdiii' = Hdcdiii, 'Hdcdiv' = Hdcdiv,
  'Hdcdv' = Hdcdv, 'Hx' = Hx, 'Hxx' = Hxx, 'Hxxx' = Hxxx, 'AgHdcdf' = AgHdcdf, 'cond' = cond))
}

# Calculate aggregate demand and demand for each observed and unobserved type
Rh_3x5Ivh <- function(h,pxy,pxxy,pxxxy,upi,upii,upiii,upiv,upv,si,sii,siii,p1,p2,p3, vp){

  #This function calculates the quality distribution for an arbitrary set of prices.

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

  budgeti <- ycfi-vp[1:(length(vp)-1)]-upi$kappa
  budgetii <- ycfii-vp[1:(length(vp)-1)]-upii$kappa
  budgetiii <- ycfiii-vp[1:(length(vp)-1)]-upiii$kappa
  budgetiv <- ycfiv-vp[1:(length(vp)-1)]-upiv$kappa
  budgetv <- ycfv-vp[1:(length(vp)-1)]-upv$kappa


  if (realyi && realyii && realyiii && realyiv && realyv &&
      monoyi && monoyii && monoyiii && monoyiv && monoyv &&
      all(ycfi>0) && all(ycfii>0) && all(ycfiii>0) && all(ycfiv>0) && all(ycfv>0) &&
      all(budgeti>0) && all(budgetii>0) && all(budgetiii>0) && all(budgetiv>0) && all(budgetv>0)){
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

    Hdcdi <- c()
    Hdcdii <- c()
    Hdcdiii <- c()
    Hdcdiv <- c()
    Hdcdv <- c()
    
    
    #demand CDFs
    Hdcdi[1] <- hdi[1]
    Hdcdii[1] <- hdii[1]
    Hdcdiii[1] <- hdiii[1]
    Hdcdiv[1] <- hdiv[1]
    Hdcdv[1] <- hdv[1]

    for (j in (2:(length(h)-1))){
      Hdcdi[j] <- Hdcdi[j-1]+hdi[j]
      Hdcdii[j] <- Hdcdii[j-1]+hdii[j]
      Hdcdiii[j] <- Hdcdiii[j-1]+hdiii[j]
      Hdcdiv[j] <- Hdcdiv[j-1]+hdiv[j]
      Hdcdv[j] <- Hdcdv[j-1]+hdv[j]
    }
    

    #Aggregate the CDFs (get aggregate demand distribution)
    H <- as.matrix(rbind(Hdcdi, Hdcdii, Hdcdiii, Hdcdiv, Hdcdv))
    s <- c(si, sii, siii)

    ax1 <- (as.matrix(s)%*%p1)
    bx1 <- ax1[1,]

    ax2 <- (as.matrix(s)%*%p2)
    bx2 <- ax2[2,]

    ax3 <- (as.matrix(s)%*%p3)
    bx3 <- ax3[3,]

    Hx <- as.vector(bx1%*%H)
    Hxx <- as.vector(bx2%*%H)
    Hxxx <- as.vector(bx3%*%H)

    AgHdcdf <- Hx + Hxx + Hxxx

    cond<-1
  } else {
    AgHdcdf = c()
    Hdcdi =c()
    Hdcdii = c()
    Hdcdiii = c()
    Hdcdiv = c()
    Hdcdv = c()
    Hx  = c()
    Hxx = c()
    Hxxx = c()
    cond=0
  }
  return(list('Hdcdi' = Hdcdi,'Hdcdii' = Hdcdii,'Hdcdiii' = Hdcdiii, 'Hdcdiv' = Hdcdiv, 'Hdcdv' = Hdcdv,
              'Hx' = Hx, 'Hxx' = Hxx, 'Hxxx' = Hxxx, 'AgHdcdf' = AgHdcdf, 'cond' = cond))
}

# Calculate Supply determined by the ratio of values for each quality
# and the elasticity of supply zeta

RhSl2PeriodsI <- function(h, vp1, vp2, V, zeta, G1, Pt) {
  # This function calculates the supply when supply is determined by a supply
  # function that is a function of the changes in price in each quality. The
  # function is specified for the discrete model with heterogeneous types.

  # PRICING FUNCTION RATIO
  # Linear pricing function for the first period
  v1 <- vp1
  v2 <- vp2

  # In equilibrium for t=1
  R1 <- G1

  V1 <- c()
  V2 <- c()
  test <- c()
  test2 <- c()

  # Corresponding values
  for (j in 1:length(h)) {
    # t=1
    d1 <- (v1[j] - V$V[2:length(V$V)])**2
    b1 <- which.min(d1)

    # V$k1ypolyf gives the user cost estimated as in equation 32. This
    # was calculated before the estimation algorithm, as this only
    # needs to be done once
    if (b1 > length(V$k1ypolyf)) {
      k1 <- V$k1ypolyf[length(V$k1ypolyf)]
    } else {
      k1 <- V$k1ypolyf[b1]
    }

    V1[j] <- v1[j] / k1
    test[j] <- v1[j] - V$V[b1]

    # t=2
    d2 <- (v2[j] - V$Vt[2:length(V$Vt)])**2
    b2 <- which.min(d2)

    if (b2 > length(V$k2ypolyf)) {
      k2 <- V$k2ypolyf[length(V$k2ypolyf)]
    } else {
      k2 <- V$k2ypolyf[b2]
    }

    V2[j] <- v2[j] / k2
    test2[j] <- v2[j] - V$Vt[b2]
  }

  # t=1
  # Calculate approx pdf
  #G1 não está se comportando como uma cdf.
  #O que era uma matriz acaba por virar um vetor
  g1 <- c()
  g1[1] <- G1[1]
  for (i in 2:length(G1)) {
    g1[i] <- G1[i] - G1[i - 1]
  }
  g1e <- 1 - G1[length(G1)]
  # CDF
  # The distribution of houses with quality h is G1 in base period/place; this means that,
  # normalizing the population to be 1 in the 1st period, we have a g1(1) = q1(1)
  # houses of quality h(1) and so on.
  r1 <- g1
  r1e <- g1e
  R1e <- R1[length(R1)] + r1e
  N1 <- 1
  q1 <- r1 * N1
  q1e <- r1e * N1
  q1t <- c(q1, q1e)
  

  # t=2
  # Calculate approx pdf
  # We have q2 = q1(V2/V1)**(zeta*Pt) houses of quality h1 in the second period.
  vratioV <- (V2 / V1)
  vratioh <- (V2 / V1)**(Pt * zeta)
  q2 <- q1 * (vratioh[1:(length(vratioh) - 1)])
  q2e <- q1e * vratioh[length(vratioh)]
  q2t <- q1t * vratioh
  N2 <- sum(q2t)
  r2 <- q2 / N2
  r2e <- q2e / N2
  r2t <- q2t / N2
  
  # CDF
  R2 <- vector(length = 13)
  R2[1] <- r2[1]
  for (i in 2:(length(h) - 1)) {
    R2[i] <- R2[i - 1] + r2[i]
  }
  R2e <- R2[length(R2)] + r2e
  R2t <- c(R2, R2e)

  # Growth Rate
  vratioV <- (V2 / V1)
  vratioh <- (V2 / V1)**(Pt * zeta)

  grateV <- (vratioV)**(1 / Pt) - 1
  grateh <- (vratioh)**(1 / Pt) - 1

  meangrateV <- sum(grateV * r2t)
  meangrateh <- sum(grateh * r2t)
  
  return(list('R2' = R2, 'N2' = N2, 'grateh' = grateh))
}

# Calculate Supply for the 2 metro areas model for the second or reference metro area,
# determined by the ratio of values for each quality and
# the elasticity of supply zeta


RhSl2PeriodsIm <- function(h, vp1m, vp2m, Vref, zeta, G1, Pt, Vbase) {
  # This function calculates the supply when supply is determined by a supply
  # function that is a function of the changes in price in each quality. The
  # function is specified for the discrete model with heterogeneous types.

  # PRICING FUNCTION RATIO
  # Linear pricing function for the first period
  v1 <- vp1m
  v2 <- vp2m
  V1<-c()
  test<-c()
  V2<-c()
  test2<-c()

  # Corresponding values
  for (j in 1:length(h)) {
    # t=1
    d1 <- (v1[j] - Vbase$V[2:length(Vbase$V)])**2
    b1 <- which.min(d1)

    # k1ypolyf and k2ypolyf give the user cost estimated as in equation 32. This
    # was calculated before the estimation algorithm, as this only
    # needs to be done once
    if (b1 > length(Vbase$k1ypolyf)) {
      k1 <- Vbase$k1ypolyf[length(Vbase$k1ypolyf)]
    } else {
      k1 <- Vbase$k1ypolyf[b1]
    }

    V1[j] <- v1[j] / k1
    test[j] <- v1[j] - Vbase$V[b1]

    # t=2
    d2 <- (v2[j] - Vref$Vt[2:length(Vref$Vt)])**2
    b2 <- which.min(d2)

    if (b2 > length(Vref$k2ypolyf)) {
      k2 <- Vref$k2ypolyf[length(Vref$k2ypolyf)]
    } else {
      k2 <- Vref$k2ypolyf[b2]
    }

    V2[j] <- v2[j] / k2
    test2[j] <- v2[j] - Vref$Vt[b2]
  }

  # t=1
  # Calculate approx pdf
  g1 <- rep(0, length(G1))
  g1[1] <- G1[1]
  for (i in 2:length(G1)) {
    g1[i] <- G1[i] - G1[i - 1]
  }
  g1e <- 1 - G1[length(G1)]

  # CDF
  # the distribution of houses with quality h is G1 in the base period/place; this means that,
  # normalizing the population to be 1 in the 1st period, we have a g1(1) = q1(1)
  # houses of quality h(1) and so on.
  R1 <- G1
  r1 <- g1
  r1e <- g1e
  R1e <- R1[length(R1)] + r1e
  N1 <- 1
  q1 <- r1 * N1
  q1e <- r1e * N1
  q1t <- c(q1, q1e)

  # t=2
  # Calculate approx pdf
  # We have q2 = q1(V2/V1)**(zeta*Pt) houses of quality h1 in the second period.
  vratioV <- (V2 / V1)
  vratioh <- (V2 / V1)**(Pt * zeta)
  q2 <- q1 * (vratioh[1:(length(vratioh) - 1)])
  q2e <- q1e * vratioh[length(vratioh)]
  q2t <- q1t * vratioh
  N2m <- sum(q2t)
  r2 <- q2 / N2m
  r2e <- q2e / N2m
  r2t <- q2t / N2m

  # CDF
  R2m <- vector(length = 13)
  R2m[1] <- r2[1]
  for (i in 2:(length(h) - 1)) {
    R2m[i] <- R2m[i - 1] + r2[i]
  }
  R2e <- R2m[length(R2m)] + r2e
  R2t <- c(R2m, R2e)

  # Growth Rate
  vratioV <- (V2 / V1)
  vratioh <- (V2 / V1)**(Pt * zeta)

  grateV <- (vratioV)**(1 / Pt) - 1
  gratehm <- (vratioh)**(1 / Pt) - 1

  meangrateV <- sum(grateV * r2t)
  meangrateh <- sum(gratehm * r2t)

  return(list('R2m' = R2m, 'N2m' = N2m, 'gratehm' = gratehm))
}

#Calculates if utility of housing component is well defined for the normalization in equation (26)
utilconpositive <- function(up, hg) {
  # Order: alpha, phi, eta, gamma, kappa
  c <- 1 - up$phi * ((hg + up$eta)**up$gamma) - 0.001 # Compute nonlinear inequalities at x.
  # ceq <- []
  # Compute nonlinear equalities at x.
  return(c)
}

# Evaluate the utility function in equation 26 for an income of y, a
# consumption of housing of h at prices v
utility <- function(y,h,phi,gamma,alpha,eta,kappa, v){
  A <- log(1-phi*((h+eta**gamma)))
  B <- (1/alpha)*log(y-v-kappa)
  U <- A + B
  return(U)
}

Vhinvert <- function(v,con,muy,sigy,tau,omg,phi,gam,alpha,eta){
  #keyboard
  #v,con,
  #v(1,1),muy,sigy,tau,omg,phi,gam,alpha
  tht <- 0
  a <- muy-(sigy/tau)*omg
  A <- exp(a)
  b <- (sigy/tau)
  F <- 1-((v**(1-b))/A)
  Ffirst <- F(1,1)
  C <- exp((b-1)*con)
  Z <- (F/C)**(1/((b-1)*alpha))
  h <- ((1-Z)/phi)**(1/gam)-eta
  return(h)
}

VhinvertGLN4 <- function(v,muy,sigma,tau,omega,phi,gammap,alpha,eta,theta){
  #We are changing the scale of theta everywhere
  theta <- 1000*theta
  a <- muy-(sigma/tau)*omega
  A <- exp(a)
  b <- (sigma/tau)
  vtheta <- v+theta
  F <- 1-((vtheta**(1-b))/A)
  S <- F**(1/(alpha*(b-1)))
  Z <- ((1-S)/phi)**(1/gammap)
  h <- Z-eta;

  return(h)
}

#Income as a function of quality, obtained from the FOCs, when pricing function is linear
y_FOC_linearVh <- function(h,p,alpha,phi,gammap,eta,kappa){
  S <- p/(alpha*phi*gammap)
  A <- ((h+eta)**(-gammap))-phi
  M <- S*(h+eta)*A
  y <- kappa + p*h - M
  return(y)
}

#Income as a function of quality, obtained from the FOCs
y_FOC_Vh <- function(h,alpha,phi,gammap,eta,kappa,v,vp){
  heta <- (h+eta)
  S <- 1 - phi*(heta**gammap)
  L <- -phi*gammap*alpha*(heta**(gammap-1))
  M <- (vp*(S))/(L)
  y <- kappa + v + M
  return(y)
}

#Função que calcula os cutoffs na renda que tornam o agente indiferente entre qualidades j e j+1
#Variáveis: h (vetor de qualidades), v(vetor de aluguéis) e up (lista nomeada dos valores dos parâmetros da função de utilidade)

Ycutoff<-function(h,up,v){
  s=log(1-up$phi*((h+up$eta)**up$gamma))
  E<-c()
  L<-c()
  M<-c()
  ycf<-c()
  for (i in 1:(length(h)-1)){
    E<-c(E, exp((s[i+1]-s[i])*up$alpha))
    L<-c(L, 1-E[i])
    M<-c(M, v[i]-E[i]*v[i+1])
    ycf<-c(ycf, (M[i]/L[i])+up$kappa)
  }
  return(unname(ycf))
}

# Função que determina os quantis da distribuição de cutoffs de renda (?)
# Variáveis: y (vetor de rendas), v(vetor de aluguéis) e up (lista nomeada dos valores dos parâmetros da função de utilidade)

Ycutoffinv<-function(y, up, v){
  opts<-list(abstol=1e-10,reltol=1e-10,maxit=1e9)
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
  quant<-optim(log(c(hini, hinip1)),hconsumed,control=opts)
  return(c('h_hat'=exp(quantpar),'fval'=quant$value))
}

# This function calculates the correlation between income and rent for each
# unobserved type, and the corresponding implied correlations for the observed types.
# Then, it calculates moments based on the difference between the correlations implied by the model and the correlations implied by the data.

Yvcorrelations<- function(sx, sxx, sxxx, p1, p2, p3, h, upi,upii,upiii,upiv,upv, vp1, vp2, Corr){
  ycfx<-c()
  ytcfx<-c()
  vcfx<-c()
  vtcfx<-c()
  ycfxx<-c()
  ytcfxx<-c()
  vcfxx<-c()
  vtcfxx<-c()
  ycfxxx<-c()
  ytcfxxx<-c()
  vcfxxx<-c()
  vtcfxxx<-c()
  
  min_value<-function(A){
    if (is.vector(A)){
      return(min(A))
    }
    mv<-c()
    for (i in 1:ncol(A)){
      mv<-c(mv,min(A[,i]))
    }
    return(mv)
  }
  min_index<-function(A){
    mv<-min_value(A)
    if (is.vector(A)){
      return(match(mv,A))
    }
    mi<-c()
    for (i in 1:ncol(A)){
      mi<-c(mi,match(mv[i],A[,i]))
    }
    return(mi)
  }
  
  s<-c(sx,sxx,sxxx)
  sdiag<-diag(s)
  
  Pr<-rbind(p1,p2,p3)
  L<-t(Pr)
  
  pr_xi<-L%*%sdiag
  
  testprob<-sum(sum(pr_xi))
  
  pr_xi_vector<-c(pr_xi[,1], pr_xi[,2], pr_xi[,3])
  cum_prix_vector<-c()
  cum_prix_vector[1]<-pr_xi_vector[1]
  for (i in 2:length(pr_xi_vector)){
    cum_prix_vector[i]<-cum_prix_vector[i-1] + pr_xi_vector[i]
  }
  
  stypes<-length(s)
  cum_prix<-cbind(cum_prix_vector[1:length(p1)], cum_prix_vector[(length(p1)+1):(length(p1)*(stypes-1))], cum_prix_vector[(length(p1)*(stypes-1)+1):(length(p1)*stypes)])
  numdraw_x <- 0
  numdraw_xx <- 0
  numdraw_xxx <- 0

  r<-1000
  counter<-1
  while(min(c(numdraw_x, numdraw_xx, numdraw_xxx))<r){
    draw<-runif(1)
    
    locate<-(draw-cum_prix)**2
    a<-colMins(locate, value=TRUE)
    b<-colMins(locate)
    c<-min(a)
    xdraw<-min_index(a)
    idraw<-b[xdraw]
    counter<-counter+1
    
    if(xdraw==1){
      numdraw_x<-numdraw_x+1
    } else if(xdraw==2){
      numdraw_xx<-numdraw_xx+1
    } else if(xdraw==3){
      numdraw_xxx<-numdraw_xxx+1
    }
    
    if(idraw==1){
      ycf <- Ycutoff(h,upi,vp1)
      ytcf <- Ycutoff(h,upi,vp2)
    } else if (idraw==2){
      ycf <- Ycutoff(h,upii,vp1)
      ytcf <- Ycutoff(h,upii,vp2)
    } else if (idraw==3){
      ycf <- Ycutoff(h,upiii,vp1)
      ytcf <- Ycutoff(h,upii,vp2)
    } else if (idraw==4){
      ycf <- Ycutoff(h,upiv,vp1)
      ytcf <- Ycutoff(h,upii,vp2)
    } else if (idraw==5){
      ycf <- Ycutoff(h,upv,vp1)
      ytcf <- Ycutoff(h,upii,vp2)
    }
    
    htrim<-length(h)-1
    if(xdraw==1){
      ycfx<-c(ycfx,ycf)
      ytcfx<-c(ytcfx,ytcf)
      vcfx<-c(vcfx,vp1[1:htrim])
      vtcfx<-c(vtcfx,vp2[1:htrim])
    }
    else if(xdraw==2){
      ycfxx<-c(ycfxx,ycf)
      ytcfxx<-c(ytcfxx,ytcf)
      vcfxx<-c(vcfxx,vp1[1:htrim])
      vtcfxx<-c(vtcfxx,vp2[1:htrim])
    }
    else if(xdraw==3){
      ycfxxx<-c(ycfxxx,ycf)
      ytcfxxx<-c(ytcfxxx,ytcf)
      vcfxxx<-c(vcfxxx,vp1[1:htrim])
      vtcfxxx<-c(vtcfxxx,vp2[1:htrim])
    }
  }
  

  corryv_x <- cor(ycfx,vcfx)
  corrtyv_x <- cor(ytcfx,vtcfx)
  
  corryv_xx <- cor(ycfxx,vcfxx)
  corrtyv_xx <- cor(ytcfxx,vtcfxx)
  
  corryv_xxx <- cor(ycfxxx,vcfxxx)
  corrtyv_xxx <- cor(ytcfxxx,vtcfxxx)
  
  difx<- (corryv_x - Corr[['corr_x1']])**2
  diftx<- (corrtyv_x - Corr[['corrt_x1']])**2
  
  difxx<- (corryv_xx - Corr[['corr_x2']])**2
  diftxx<- (corrtyv_xx - Corr[['corrt_x2']])**2
  
  difxxx<- (corryv_xxx - Corr[['corr_x3']])**2
  diftxxx<- (corrtyv_xxx - Corr[['corrt_x3']])**2
  return(100* sum(c(difx, diftx, difxx, diftxx, difxxx, diftxxx)))
}
