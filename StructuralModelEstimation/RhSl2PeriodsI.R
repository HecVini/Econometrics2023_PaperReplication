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
  
  V1 <- vp1
  V2 <- vp2
  test <- vp1
  test2 <- vp2
  
  # Corresponding values
  for (j in 1:length(h)) {
    # t=1
    d1 <- (v1[j] - V$V[2:length(V$V)])^2
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
    d2 <- (v2[j] - V$Vt[2:length(V$Vt)])^2
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
  g1 <- rep(0, length(G1))
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
  # We have q2 = q1(V2/V1)^(zeta*Pt) houses of quality h1 in the second period.
  vratioV <- (V2 / V1)
  vratioh <- (V2 / V1)^(Pt * zeta)
  q2 <- q1 * (vratioh[1:(length(vratioh) - 1)])
  q2e <- q1e * vratioh[length(vratioh)]
  q2t <- q1t * vratioh
  N2 <- sum(q2t)
  r2 <- q2 / N2
  r2e <- q2e / N2
  r2t <- q2t / N2
  
  # CDF
  R2 <- rep(0, length(h))
  R2[1] <- r2[1]
  for (i in 2:(length(h) - 1)) {
    R2[i] <- R2[i - 1] + r2[i]
  }
  R2e <- R2[length(R2)] + r2e
  R2t <- c(R2, R2e)
  
  # Growth Rate
  vratioV <- (V2 / V1)
  vratioh <- (V2 / V1)^(Pt * zeta)
  
  grateV <- (vratioV)^(1 / Pt) - 1
  grateh <- (vratioh)^(1 / Pt) - 1
  
  meangrateV <- sum(grateV * r2t)
  meangrateh <- sum(grateh * r2t)
  
  return(list(R2 = R2, N2 = N2, grateh = grateh))
}
