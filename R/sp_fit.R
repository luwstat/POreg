
sp_fit <- function(L,R,x,order = 3,equal_space = T,nknot, myknots, diagnosis = TRUE,conv_cri = 1e-9){
###################Functions will be used in the EM algorithm###############
  library(nleqslv)
  library(Matrix)
  library(ggplot2)
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
  #iv is the X's
  iv <- as.matrix(x)
  P <- ncol(iv)
  n <- nrow(iv)
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
  bt <- rep(1,P)
  gama <- rep(1,K)
  df <- rep(1,P)
  ite <- 0
  Ez <- rep(0,n)
  Ew <- rep(0,n)
  Epsi <- rep(0,n)
  Ephi <- rep(0,n)
  EU <- matrix(0,n,K)

  while( t(df)%*%df > conv_cri & ite < 20000){
    #exp(x_i*beta) is n*1 vector
    exb <- exp(iv%*%bt)
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
    A <- apply( matrix(EU[d_0 == 1,],ncol = K),2,sum) + apply( matrix(EZ[d_1 == 1,],ncol = K),2,sum) + apply( EW[d_2 == 1,],2,sum)
    #B is a n*K matrix
    B <- t(Il)*Epsi + ( t(Il)*(d_0 + d_3) + t(Ir)*(d_1 + d_2) )*Ephi
    one <- rep(1,K)
    btf <- function(x){
      y <- numeric(length(bt))
      for(h in 1:length(bt)){
        y[h] <- (d_0 + Ez + Ew)%*%iv[,h] - (one%*%((t(B)*A)/as.vector(t(B)%*%exp(iv%*%x)))) %*% (iv[,h]*exp(iv%*%x))
      }
      y
    }


    btstart <- rep(0,length(bt))
    #solve and get the updated bt
    sol <- nleqslv(btstart,btf,method="Newton")
    btnew <- sol$x

    gamanew <-  A/( as.vector(t(B)%*%exp(iv%*%btnew)) )

    df <- btnew - bt
    bt <- as.vector(btnew)
    gama <- as.vector(gamanew)
    ite <- ite + 1

  }

  ################calculate covariance matrix##################
  #First calculate all the expectationes and variances given the theta_hat
  ######Expectations######
  #x_i*beta is n*1 vector
  xb <- iv%*%bt
  #exp(x_i*beta) is n*1 vector
  exb <- exp(xb)
  #Lambda_0(R) and Lambda_0(L), n*1 vector
  Bsr <- t(Ir)%*%gama
  Bsl <- t(Il)%*%gama
  #Lambda_0(.)exp(x*beta) + 1, n*1 vector
  Bre <- Bsr*exb + 1
  Ble <- Bsl*exb + 1

  #Ez not zero when d_1=1
  Ez[d_1 == 1] <- Bre[d_1 ==1]
  #Ew not zero when d_2=1
  Ew[d_2 == 1] <- Bre[d_2 ==1]/Ble[d_2 ==1]
  #Ephi and epsi are n*1 vectors
  Ephi[d_0 == 1 | d_3 ==1] <- 1/Ble[d_0 == 1 | d_3 ==1]
  Ephi[d_1 == 1] <- (Bre[d_1 == 1] + 1)/Bre[d_1 == 1]
  Ephi[d_2 == 1] <- (Bre + Ble)[d_2 == 1]/(Ble*Bre)[d_2 == 1]
  #Epsi only when delta_0 = 1
  Epsi[d_0 == 1] <- 1/Ble[d_0 ==1]
  #EZ and EW are n*K matrix
  EZ <- t(Ir*gama)*as.vector(Ez/Bsr^d_1)
  EW <- t(I*gama)*as.vector(Ew/(Bsr - Bsl)^d_2)
  #EU is a n*K matrix, the ^del is used to get rid of the 0 denomator
  EU[d_0 == 1,] <- t(Ml[,d_0 == 1]*gama)/as.vector(t(Ml[,d_0 == 1])%*%gama)

  #########variances#############
  #variance phi is n*1 vector
  Vphi <- rep(0,n)
  Vphi[d_0 == 1 | d_3 ==1] <- Ephi[d_0 == 1 | d_3 ==1]^2
  Vphi[d_1 == 1] <- Bre[d_1 == 1]^(-2) + 1
  Vphi[d_2 == 1] <- Ble[d_2 == 1]^(-2) + Bre[d_2 == 1]^(-2)
  #variance psi is n*1 vector
  Vpsi <- rep(0,n)
  Vpsi[d_0 == 1]  <- Epsi[d_0 == 1]^2
  #Varicnce of z not zero when d_1=1 (vector of Var(z_i))
  Vz <- rep(0,n)
  Vz[d_1 == 1] <- Bre[d_1 == 1]*(Bre[d_1 == 1] - 1)
  #Varicnce of w not zero when d_2=1 (vector of Var(w_i))
  Vw <- rep(0,n)
  Vw[d_2 == 1] <- Ew[d_2 == 1]*(Ew[d_2 == 1] - 1  )
  #Variance of Z and W are n*K matrix
  VZ <- EZ*(1 - t(Ir*gama)/as.vector(Bsr)^d_1) + Vz*(t(Ir*gama)/as.vector(Bsr)^d_1)^2
  VW <- EW*(1 - t(I*gama)/as.vector(Bsr - Bsl)^d_2) + Vw*(t(I*gama)/as.vector(Bsr - Bsl)^d_2)^2
  #VU variance of U_il is a n*K matrix, the il_th entry is var(u_il)
  VU <- EU*(1- EU)

  ####part 1 Q########################################
  Q <- matrix(0,(K+P),(K+P))
  #same A and B as before (in the EM part). A is a K*1 vector, B is a n*K matrix
  A <- apply( matrix(EU[d_0 == 1,],ncol = K),2,sum) + apply( matrix(EZ[d_1 == 1,],ncol = K),2,sum) + apply( matrix(EW[d_2 == 1,],ncol = K),2,sum)
  B <- t(Il)*Epsi + ( t(Il)*(d_0 + d_3) + t(Ir)*(d_1 + d_2) )*Ephi
  #elements in xx^t as a vector, four entries as a group
  e <- as.vector(t(do.call(cbind,replicate(P,iv,simplify = F)))) * rep(as.vector(t(iv)),each = P)
  #coefficients for each matrix, repeat each elememt p times so it match the vector above
  co <-  rep(B%*%gama*exb, each = P^2)
  #fill in the corresponding place in Q
  Q[seq(1,P),seq(1,P)] <- - matrix(apply(matrix(e*co,P^2,n),1,sum),P,P)

  Q[seq(P+1,K+P),seq(1,P)] <- - t(B)%*%(iv*as.vector(exb))

  Q[seq(1,P),seq(P+1,P+K)] <- t( Q[seq(P+1,K+P),seq(1,P)])

  diag(Q[seq(1+P,P+K),seq(1+P,P+K)]) <- - gama^(-2)*A

  ####part 2 VC ###########################################
  vc <- matrix(0,(K+P),(K+P))
  D_i <- Bsl*(d_0 + d_3) + Bsr*(d_1 + d_2)
  D_il <- t(Il)*(d_0 + d_3) + t(Ir)*(d_1 + d_2)

  #cov(l_c/beta,l_c/beta )
  e_covc <- Vpsi*(Bsl*exb*d_0)^2 + Vphi*(D_i*exb)^2 + Vz*d_1 + Vw*d_2 - 2*(Bre - 1)^2*d_1 - 2*((Bre - Ble)/Ble^2)*(Bre - 1)*d_2
  covc <-  rep(e_covc, each = P^2)
  vc[seq(1,P),seq(1,P)] <- matrix(apply(matrix(e*covc,P^2,n),1,sum),P,P)

  #coefficients for cov(l_c/beta,l_c/gama )
  co1 <- t(Il)*as.vector(Bsl*exb^2)*Vpsi*d_0 + D_il*as.vector(D_i*exb^2)*Vphi
  co2 <- t(Ir)*as.vector(d_1*exb*Bre) + t(I)*as.vector( Ew*d_2*exb/Ble )
  co3 <- t(Ir)*as.vector(d_1*Bsr*exb^2) + t(Ir)*as.vector(d_2*(Bsr - Bsl)*exb^2*Ble^(-2)) + t(Ir)*as.vector(Bsr*exb^2*d_1) + t(I)*as.vector(d_2*Bsr*exb^2/Ble^2)
  co_t <- co1 + co2 - co3

  vc[seq(1,P),seq(P+1,P+K)] <- t(iv)%*%co_t
  vc[seq(P+1,K+P),seq(1,P)] <- t(vc[seq(1,P),seq(P+1,P+K)])
  #coefficients for cov(l_c/gama, l_c/gama) on diagnal
  dg1 <- t(VU*d_0 + VW*d_2 + VZ*d_1)/gama^2
  dg2 <- ( t(Il^2)*Vpsi*d_0 + Vphi*D_il^2 )*as.vector(exb^2)
  dg3 <- t(Ir)*t(I)*d_2*as.vector(exb^2/(Ble^2)) + t(Ir^2)*d_1*as.vector(exb^2)

  diag(vc[seq(P+1,K+P),seq(P+1,K+P)]) <- apply(t(dg1) + dg2 -2*dg3,2,sum)
  #coefficients for cov(l_c/gama_l, l_c/gama_k) off diagnal
  for(l in 1:K){
    for(m in 1:K){
      if(l != m){
        part_1 <- -(Ml[l,]*Ml[m,]*d_0/as.vector(t(Ml)%*%gama)^(2*d_0)) + d_1*Ir[l,]*Ir[m,]*(Vz-Ez)/Bsr^2 + d_2*I[l,]*I[m,]*(Vw - Ew)/(Bsr - Bsl)^(2*d_2)
        part_2 <- (Il[l,]*Il[m,]*d_0*Vpsi + D_il[,l]*D_il[,m]*Vphi)*exb^2
        part_3 <- d_1*Ir[m,]*Ir[l,]*exb^2 + d_2*Ir[m,]*I[l,]*exb^2/Ble^2
        part_4 <- d_1*Ir[l,]*Ir[m,]*exb^2 + d_2*Ir[l,]*I[m,]*exb^2/Ble^2
        vc[P+l,P+m] <- sum( part_1 + part_2 - part_3 - part_4 )
      }
    }
  }

  v <- -(Q + vc)

  ####part3
  tol <- 1e-7
  if(sum(gama < tol) == 0){
    vv <- v
  } else {
    rg=1:K
    index=rg[gama < tol]+P
    vv=v[-index,-index]
  }

  if( rcond(vv) > .Machine$double.eps ){
    se_theta = sqrt(diag(solve(vv))[1:P])
  } else {
    se_theta = sqrt(diag(solve(vv + diag(1e-4,nrow(vv),ncol(vv))))[1:P])
  }


  CI_lower <- bt - 1.96*se_theta
  CI_upper <- bt + 1.96*se_theta
  ci <- cbind(CI_lower,CI_upper)

  #define function log??likelihood AND calculate AIC and BIC
  llhd <- sum(log(((t(Ml)%*%gama) * exp(iv%*%bt)/( t(Il)%*%gama * exp(iv%*%bt) + 1 )^2 )^d_0 * (1 - 1/( t(Ir)%*%gama * exp(iv%*%bt) + 1 ))^(d_1) * ( 1/( t(Il)%*%gama * exp(iv%*%bt) + 1) - 1/( t(Ir)%*%gama * exp(iv%*%bt) + 1))^(d_2) * ( 1/( t(Il)%*%gama * exp(iv%*%bt) + 1) )^(d_3)))
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

  if(diagnosis == TRUE){
    ####################### model diagnosis with logistic distribution ############################
    epsl <- log(as.vector(t(as.matrix(gama))%*%Il)) + iv%*%bt
    epsr <- log(as.vector(t(as.matrix(gama))%*%Ir)) + iv%*%bt

    library(icenReg)
    estR <- ic_np(cbind(epsl,epsr))
    diplot <- plot(estR,ylab = "Survival function",xlab = "v")
    curve(1-plogis(x,0,1),add = TRUE,col = "red")
    legend(min(c(epsl,epsr)),0.2,legend = c('S(v)','S(v)_hat'), col = c('red','black'))
    A <- getSCurves(estR)

    jps <- c(A$Tbull_ints[,1],A$Tbull_ints[,2])
    st1 <- A$S_curves$baseline[findInterval(jps,A$Tbull_ints[,1])]
    st2 <- A$S_curves$baseline[findInterval(jps,A$Tbull_ints[,1])-1]
    xx1 = (1-plogis(jps)) - st1
    xx2 = (1-plogis(jps))- st2
    STATISTIC <-  max(abs(c(xx1,xx2)))
  }else{
    STATISTIC = NULL
    diplot = NULL
  }

  #Return the result
  list(coefficient_est = as.data.frame(cbind(bt,se_theta,ci),row.names = colnames(x)), spline_coef = gama, knots = knots, AIC = AIC, BIC = BIC, SK_Statistic = STATISTIC, Baseline_Surv = bsl_surv, Baseline_hazard = bsl_hz, Baseline_odds = bsl_odds )

}

