np_fit <- function(L,R,order = 3,equal_space = T,nknot, myknots,conv_cri = 1e-9){
  if(equal_space == T){
    #the number of interior knots
    ik <- nknot-2
    #equal space knots
    #knots at the right end + 0.01 because all the basis equals to 0 at the end points if not
    mx <-  max(setdiff(c(R,L),Inf)) + 0.01
    knots <- seq(0,mx,length.out = nknot)
  }else{
    knots <- myknots
    ik <- length(knots) - 2
  }

  #Degree for Ispline
  dgr <- order
  #number of parameters
  K <- dgr + ik
  #number of observations
  n <- length(L)
  #delta
  d_1 <- rep(0,n)
  d_2 <- rep(0,n)
  d_3 <- rep(0,n)
  d_0 <- as.numeric(L == R)
  for(i in which(d_0 ==0)){
    if (L[i] == 0){
      d_1[i] <- 1
    }
    else if(R[i] == Inf){
      d_3[i] <- 1
    }
    else {
      d_2[i] <- 1
    }
  }

  Ml <- Mspline(L,dgr,knots)
  Mr <- Mspline(R,dgr,knots)
  Il <- Ispline(L,dgr,knots)
  Ir <- Ispline(R,dgr,knots)
  #rate for interval censored data
  I <- matrix(0,K,n)
  I[,d_2 == 1] <- Ir[,d_2 == 1] - Il[,d_2 == 1]

  #set an initial value
  gama <- rep(1,K)
  df <- rep(1,K)
  ite <- 0
  Ez <- rep(0,n)
  Ew <- rep(0,n)
  Epsi <- rep(0,n)
  Ephi <- rep(0,n)
  EU <- matrix(0,n,K)

  while( t(df)%*%df > conv_cri & ite < 20000){
    #exp(x_i*beta) is n*1 vector
    exb <- rep(1,n)
    #Lambda_0(R) and Lambda_0(L), n*1 vector
    Bsr <- t(Ir)%*%gama
    Bsl <- t(Il)%*%gama
    #Lambda_0(.)exp(x*beta) + 1, n*1 vector
    Bre <- Bsr*exb + 1
    Ble <- Bsl*exb + 1

    #Ez does not equal to zero only when d_1=1
    Ez[d_1 == 1] <- Bre[d_1 ==1]
    #Ew does not equal to zero only when d_2=1
    Ew[d_2 == 1] <- Bre[d_2 ==1]/Ble[d_2 ==1]
    #Ephi and Epsi are n*1 vectors
    Ephi[d_0 == 1 | d_3 ==1] <- 1/Ble[d_0 == 1 | d_3 ==1]
    Ephi[d_1 == 1] <- (Bre[d_1 == 1] + 1)/Bre[d_1 == 1]
    Ephi[d_2 == 1] <- (Bre + Ble)[d_2 == 1]/(Ble*Bre)[d_2 == 1]
    #Epsi does not equal to zero only when delta_0 = 1
    Epsi[d_0 == 1] <- 1/Ble[d_0 ==1]
    #EZ and EW are n*K matrix
    EZ <- t(Ir*gama)*as.vector(Ez/Bsr)
    EW <- t(I*gama)*as.vector(Ew/(Bsr - Bsl)^d_2)
    #EU is a n*K matrix, the ^del is used to get rid of the 0 denomator in right censored data
    EU[d_0 == 1,] <- t(Ml[,d_0 == 1]*gama)/as.vector(t(Ml[,d_0 == 1])%*%gama)

    #set the equations needed to be solved
    #first define the constant terms
    #A is a K*1 vector
    A <- apply( EU[d_0 == 1,],2,sum) + apply( EZ[d_1 == 1,],2,sum) + apply( EW[d_2 == 1,],2,sum)
    #B is a n*K matrix
    B <- t(Il)*Epsi + ( t(Il)*(d_0 + d_3) + t(Ir)*(d_1 + d_2) )*Ephi

    gamanew <-  A/( as.vector(t(B)%*%exb) )

    df <- gamanew - gama
    gama <- as.vector(gamanew)
    ite <- ite + 1

  }

  #define function log??likelihood AND calculate AIC and BIC
  llhd <- sum(log(((t(Ml)%*%gama)/( t(Il)%*%gama + 1 )^2 )^d_0 * (1 - 1/( t(Ir)%*%gama + 1 ))^(d_1) * ( 1/( t(Il)%*%gama + 1) - 1/( t(Ir)%*%gama + 1))^(d_2) * ( 1/( t(Il)%*%gama + 1) )^(d_3)))
  AIC <- 2*(P+K) - 2*llhd
  BIC <- (P+K)*log(n) - 2*llhd

  #########################################################
  # plots:odds, survival,hazard
  ########################################################
  #grids of time points
  tgrids <- seq(0,mx,0.1)
  #calculate baseline odds
  b <- Ispline(tgrids,dgr,knots)
  m <- Mspline(tgrids,dgr,knots)
  odds <- t(as.matrix(gama)) %*% b
  ff <- t(as.matrix(gama))%*%m

  #check baseline hazard rate
  hzd <- ff/(1+odds)
  bsl_hz = ggplot() + geom_line(data = as.data.frame( cbind(tgrids,t(hzd) ) ),aes(tgrids,t(hzd) )) + labs(x="t",y="h(t)")
  #check baseline odds
  bsl_odds = ggplot() + geom_line(data = as.data.frame( cbind(tgrids,t(odds) ) ),aes(tgrids,t(odds) )) + labs(x="t",y="Baseline odds")
  #check baseline survival
  sur <- 1/(1+odds)
  bsl_surv = ggplot() + geom_line(data = as.data.frame( cbind(tgrids,t(sur) ) ),aes(tgrids,t(sur) )) + labs(x="t",y="S(t)")

  #Return the result
  list(spline_coef = gamanew,knots= knots, AIC = AIC, BIC = BIC,Baseline_Surv = bsl_surv, Baseline_hazard = bsl_hz, Baseline_odds = bsl_odds)

}
