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
  plot(log(Y.Y), GLNy)
  
  # Add a scatter plot
  points(log(Y.Y), Y.Ycd, col = "black", pch = 20, cex = 0.8)
  
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
  plot(log(V.V), GLN4v)
  
  # Add a scatter plot
  points(log(V.V), V.Vcd, col = "black", pch = 20, cex = 0.8)
  
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
  plot(log(Y.Yt), GLNyt)
  
  # Add a scatter plot
  points(log(Y.Yt), Y.Ytcd, col = "black", pch = 20, cex = 0.8)
  
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
  plot(log(V.Vt), GLNvt)
  
  # Add a scatter plot
  points(log(V.Vt), V.Vtcd, col = "black", pch = 20, cex = 0.8)
  
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
