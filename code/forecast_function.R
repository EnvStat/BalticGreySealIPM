
# Note, in the code we use "bias" to denote hunting selectivity. The latter term 
# is used in the paper. The difference between the paper and code arises from 
# the reason that we used the less good "bias" in the original submission.


# inverse logit transform
logistic <- function(x) { return(1 / (1 + exp(-x))) }


# density dependent birth rate
b.N <-function (b_max, theta_0, theta_1, N) {
  b_N = exp(log(b_max) - theta_0 * (expm1(theta_1 * sum(N))))
  return(b_N)
}

# ODE for hunting 
dH_dtau <- function(tau, y,vars) {
  
  n0=vars[1]
  k=vars[2]
  E_1=vars[3]
  E_2=vars[4]
  mu=vars[5]
  dH = n0*exp(-(E_1+E_2)*(k*tau-tau^2/2) - mu*tau)*E_1*(k-tau)
  return(dH)
}

simulate<-function (N0, years.sim, Q_fi.sim, Q_sw.sim, h_S.sim, h_G.sim, bias_sw,bias_fi, g_sw, g_fi,fixed_K=0,q=1,
                         params, data) {
  tau_h = 8/12
  
  
  with(params, {
    t <- length(years.sim)
    b.old=b
    phi = S_diag
    tau_p = 0.5
    tau_h_0 = 1.5/12
    tau_h = 8/12
    a=6
    tau_h_m = 5/12
    N <- matrix(0, 12, t)
    H_sw <- rep(0,t)
    H_fi <- rep(0,t)
    U <- matrix(NA, 12 * 4, t)
    rho <- matrix(NA, 12, 4)
    A = matrix(0, 12, 12)
    A[2:6, 1:5] <- A[8:12, 7:11] <- diag(1, 5)
    A[6, 6] <- A[12, 12] <- 1
    if(fixed_K==1){
    theta1= log(1-log(bK/b0_av)/(theta0))/120000
    }
    b0 <- b0min + (b0max - b0min)/(1 + exp(-(alpha + beta * 
                                               (w * h_S.sim + (1 - w) * h_G.sim))))
    b <- rep(NA, length(years.sim))
    N[, 1] <- N0
    for (i in 1:t) {
      if (i == 1) {
        b[i] <-b.old[length(b.old)]
      }
      if (i > 1) {
        b[i] <- b.N(b0[i], theta0, theta1, N[, i - 1])
        
        N[, i] <- A %*% U[1:12, i - 1]
        newborns <- rbinom(1, N[6, i], b[i])
        newborn.females <- rbinom(1, newborns, 0.5)
        N[1, i] = N[1, i] + newborn.females
        N[7, i] = N[7, i] + newborns - newborn.females
      }
      
      if (sum(N[, i]) <= 0) {
        break
      }
      phi_tau_0 <- phi^tau_h_0
      
      
      H_tot_sw=q*Q_sw.sim[i]
      H_tot_fi=q*Q_fi.sim[i]
      
      N_tot[i]=sum(N[,i])
      

      if (sum(N[, i] %in% c(0)) > 0) {
        rho[which(N[, i] == 0), 3] = 0
        rho[which(N[, i] == 0), 4] = 0
        rho[, 1] <-  (1-rho[, 3]-rho[, 4])* phi
        rho[, 2] <- 1 - rowSums(rho[, c(1, 3:4)])
        if (min(rho) < 0) {
          rho[rho < 0] = 0
          
        }
      }
      else {
        
        rho[, 3] <- H_tot_sw*exp(bias_sw)/(sum(exp(bias_sw)*N[,i]))
        rho[, 4] <-   H_tot_fi*exp(bias_fi)/(sum(exp(bias_fi)*N[,i]))
        
        rho[,1] <- (1-rho[, 3]-rho[, 4])* phi
        rho[, 2] <- 1 - rowSums(rho[, c(1, 3:4)])
        
        
        if (sum(is.na(rho))>0){
          
          rho <- matrix(0, 12, 4)
        }
        else{
          
          if (min(rho) < 0) {
            rho <- matrix(0, 12, 4)
          }
        }
      }
      U_stacked <- matrix(0, nrow = 12, ncol = 4)
      for (j in 1:nrow(U_stacked)) {
        if (N[j, i] > 0) {
          if (sum(rho[j,])>0){
            U_stacked[j, ] <- t(rmultinom(1, N[j, i], rho[j, ]))
          }
        }
      }
      U[, i] <- c(U_stacked)
      H_sw[i]=sum(U_stacked[,3])
      H_fi[i]=sum(U_stacked[,4])
    }
    
    M = matrix(c(b, t(N), H_sw, H_fi), ncol = 15, nrow = length(years.sim))
    
    return(M)
  })
}




simulate_growth<-function (n_start, fixed_K=0,params, data) {
  
  with(params, {
    t <- 2
    phi = S_diag
    a=6
 
    A = matrix(0, 12, 12)
    A[2:6, 1:5] <- A[8:12, 7:11] <- diag(1, 5)
    A[6, 6] <- A[12, 12] <- 1
    
    T=A%*% diag(phi)
    
    
    if(fixed_K==1){
      theta1= log(1-log(bK/b0_av)/(theta0))/120000
    }
    
    b0 <- b0min + (b0max - b0min)/(1 + exp(-(alpha )))
    N<- rep(1/12,12)*n_start
    
    b <- b.N(b0, theta0, theta1, N)

    for ( j in 1:20){
      N <-T%*%N
      N[1]=b/2*N[6]
      N[7]=b/2*N[6]
      N=N*n_start/sum(N)
    }
    N <- T%*%N
    N[1]=b/2*N[6]
    N[7]=b/2*N[6]
    g=(sum(N)/n_start-1)*100
    return(g)
  })
}




