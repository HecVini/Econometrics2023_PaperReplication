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
  
  V1 <- vp1
  V2 <- vp2
  test <- vp1
  test2 <- vp2
  
  # Corresponding values
  for (j in 1:length(h)) {
    # t=1
    d1 <- (v1[j] - Vbase$V[2:length(Vbase$V)])^2
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
    d2 <- (v2[j] - Vref$Vt[2:length(Vref$Vt)])^2
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
  # We have q2 = q1(V2/V1)^(zeta*Pt) houses of quality h1 in the second period.
  vratioV <- (V2 / V1)
  vratioh <- (V2 / V1)^(Pt * zeta)
  q2 <- q1 * (vratioh[1:(length(vratioh) - 1)])
  q2e <- q1e * vratioh[length(vratioh)]
  q2t <- q1t * vratioh
  N2m <- sum(q2t)
  r2 <- q2 / N2m
  r2e <- q2e / N2m
  r2t <- q2t / N2m
  
  # CDF
  R2m <- rep(0, length(h))
  R2m[1] <- r2[1]
  for (i in 2:(length(h) - 1)) {
    R2m[i] <- R2m[i - 1] + r2[i]
  }
  R2e <- R2m[length(R2m)] + r2e
  R2t <- c(R2m, R2e)
  
  # Growth Rate
  vratioV <- (V2 / V1)
  vratioh <- (V2 / V1)^(Pt * zeta)
  
  grateV <- (vratioV)^(1 / Pt) - 1
  gratehm <- (vratioh)^(1 / Pt) - 1
  
  meangrateV <- sum(grateV * r2t)
  meangrateh <- sum(gratehm * r2t)
  
  return(list(R2m = R2m, N2m = N2m, gratehm = gratehm))
}
