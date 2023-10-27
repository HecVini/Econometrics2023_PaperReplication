# This function calculates the correlation between income and rent for each
# unobserved type, and the corresponding implied correlations for the observed types. 
# Then, it calculates moments based on the difference between the correlations implied by the model and the correlations implied by the data. 

Yvcorrelations <- function(sx, sxx, sxxx, p1, p2, p3, h, upi,upii,upiii,upiv,upv, vp1, vp2, Corr){

# - Given some probabilities mapping observed types to unobserved types, 
# I can calculate the implied joint probability of each one of the combinations for unobserved type i and observed type
# k.

  s <- c(sx, sxx, sxxx)
  sdiag <- diag(s) 
  
  Pr <- c(p1, p2, p3)
  L <- t(Pr)
  
  pr_xi <- L* sdiag
  
  
  # - These probabilities must add to 1. 
  testprob <- sum(sum(pr_xi))
  
  
  # Draw a renter using the probabilities of that matrix
  # Calculate the cumulative prob of this joint probability
  # place in a vector
  pr_xi_vector <- c(pr_xi[,1], pr_xi[,2], pr_xi[,3])
  
  cum_prix_vector(1) <- pr_xi_vector(1)
  
  for (i in 2:length(pr_xi_vector)){ #3*5 = 15
    cum_prix_vector(i) <- cum_prix_vector(i-1) + pr_xi_vector(i)
  }
  
  
  # Return to a matrix
  stypes <- length(s)
  cum_prix <- c(t(cum_prix_vector(1:length(p1))), t(cum_prix_vector(length(p1)+1:length(p1)*(stypes-1))), t(cum_prix_vector(length(p1)*(stypes-1)+1:length(p1)*stypes)))
  
  
  numdraw_x <- 0
  numdraw_xx <- 0
  numdraw_xxx <- 0
  
  #Draw numbers until there are r observations for each observed type
  r <-10000
  counter <-1
  while (min(c(numdraw_x, numdraw_xx, numdraw_xxx)<r)){
    # random number
    draw <-rand
  
    # Locate in the matrix cum_prix
  
    locate <-(draw - cum_prix)^2
    c(a, b) <-min(locate)
  
    c(c, xdraw) <-min(a)
  
    idraw <-b(xdraw)
  
    count[counter,] <-c(xdraw, idraw)
    counter <-counter+1
  
    # Count how many draws we have for observed type
    if (xdraw==1){
      numdraw_x <-numdraw_x+1
    }else if (xdraw==2){
      numdraw_xx <-numdraw_xx+1
    }else if (xdraw==3){
      numdraw_xxx <-numdraw_xxx+1
    }
  
  # For each quality in the grid and current corresponding price, find the
  # corresponding income 
  
    if (idraw == 1){
      ycf <-Ycutoff(h,upi,vp1)
      ytcf <-Ycutoff(h,upi,vp2)
      
    }else if (idraw == 2){
      ycf <-Ycutoff(h,upii,vp1)
      ytcf <-Ycutoff(h,upii,vp2)
    }else if (idraw == 3){
      ycf <-Ycutoff(h,upiii,vp1)
      ytcf <-Ycutoff(h,upii,vp2)
    }else if (idraw == 4){
      ycf <-Ycutoff(h,upiv,vp1)
      ytcf <-Ycutoff(h,upii,vp2)
    }else if (idraw == 5){
      ycf <-Ycutoff(h,upv,vp1)
      ytcf <-Ycutoff(h,upii,vp2)
    }
  
  #Assign it to the observed type according to the draw with did in the
  #cum_prix matrix
    htrim <-length(h)-1
    if (xdraw==1){
      ycfx(((numdraw_x-1)*htrim)+1:htrim*numdraw_x) <-ycf
      ytcfx(((numdraw_x-1)*htrim)+1:htrim*numdraw_x) <-ytcf
      vcfx(((numdraw_x-1)*htrim)+1:htrim*numdraw_x) <-vp1(1:htrim)
      vtcfx(((numdraw_x-1)*htrim)+1:htrim*numdraw_x) <-vp2(1:htrim)
    } else if (xdraw==2){
      ycfxx(((numdraw_xx-1)*htrim)+1:htrim*numdraw_xx) <-ycf
      ytcfxx(((numdraw_xx-1)*htrim)+1:htrim*numdraw_xx) <-ytcf
      vcfxx(((numdraw_xx-1)*htrim)+1:htrim*numdraw_xx) <-vp1(1:htrim)
      vtcfxx(((numdraw_xx-1)*htrim)+1:htrim*numdraw_xx) <-vp2(1:htrim)
    } else if (xdraw==3){
      ycfxxx(((numdraw_xxx-1)*htrim)+1:htrim*numdraw_xxx) <-ycf
      ytcfxxx(((numdraw_xxx-1)*htrim)+1:htrim*numdraw_xxx) <-ytcf
      vcfxxx(((numdraw_xxx-1)*htrim)+1:htrim*numdraw_xxx) <-vp1(1:htrim)
      vtcfxxx(((numdraw_xxx-1)*htrim)+1:htrim*numdraw_xxx) <-vp2(1:htrim)
    }
  }
  
  #Calculate correlations for each observed type
  
  corryv_x <-cor(ycfx,vcfx)
  corrtyv_x <-cor(ytcfx,vtcfx)
  
  corryv_xx <-cor(ycfxx,vcfxx)
  corrtyv_xx <-cor(ytcfxx,vtcfxx)
  
  corryv_xxx <-cor(ycfxxx,vcfxxx)
  corrtyv_xxx <-cor(ytcfxxx,vtcfxxx)
  
  
  difx <- (corryv_x(1,1) - Corr$corr_x1)^2
  diftx <- (corrtyv_x(1,1) - Corr$corrt_x1)^2
  
  difxx <- (corryv_xx(1,1) - Corr$corr_x2)^2
  diftxx <- (corrtyv_xx(1,1) - Corr$corrt_x2)^2
  
  difxxx <- (corryv_xxx(1,1) - Corr$corr_x3)^2
  diftxxx <- (corrtyv_xxx(1,1) - Corr$corrt_x3)^2
  
  M <-100* sum(c(difx, diftx, difxx, diftxx, difxxx, diftxxx))
  return(M)
}
