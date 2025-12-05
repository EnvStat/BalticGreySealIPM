#Daire Carroll University of Gothenburg 2022
#Grey seal PVA script, functions for predictive modelling

#t = number of cycles
#A
#K = carrying capacity 
#P = initial population as an age structured vector
#B = age specific birthrates as a vector
#S = age specific survival as a vector
#S_G = age specific survival as a vector (good year)
#S_B  = age specific survival as a vector (bad year)
#Yr_P = probability of a bad survival year per cycle
#R = random variation term 
#H = size of hunt, either proportion or set quota
#H_sk = hunt skew, a vector expressing how the hunt should be applied to the population
#proportion = TRUE/FALSE if hunt should be applied as a proportion of the population (TRUE) or as a set quota (FALSE) 
#stru = TRUE/FALSE if an age structured population should be returned for each cycle (TRUE) or a vector of the sums of populaiotn size for each cycle (FALSE) 

structure = function(B,S,N){
  B = B*S/2
  L = matrix(B,ncol = A)
  for(i in 1:(A-1)){
    x = numeric(A)
    x[i] = S[i]
    L = rbind(L, x)
  }
  ev = eigen(L)
  ev = as.numeric(ev$vectors[,1])
  P = c()
  for(i in 1:A){
    P[i] = ev[i]/sum(ev)*N
  } 
  return(P)
} #find dominant left eigen vector for given demographic scenario, with birthrates B (vector), survival S (vector) and popualiton size

eig_val = function(B,S){
  A = length(B)
  B = B*S/2
  L = matrix(B,ncol = A)
  for(i in 1:(A-1)){
    x = numeric(A)
    x[i] = S[i]
    L = rbind(L, x)
  }
  ev = eigen(L)
  return(as.numeric(ev$value[1]))
} #find dominant left eigen value for given demographic scenario

PopMod = function(t, A, K, P, B, S, R, H, H_sk,E,E_sk, proportion, stru){  ##Added by Vanko et al. added E. E_sk otherwise as not running (the function does not end up using it)
  
  ###housekeeping and troubleshooting begin###
  
  if(t==0){		
    
    if(stru == TRUE){
      
      return(P)
      
    }else{
      
      return(sum(P))
      
    }
    
  } 
  
  null_args = list(H_sk,E_sk)
  null_args_names = c("H_sk","E_sk")
  for(i in 1:length(null_args)){
    if(is.null(null_args[[i]])){
      assign(paste0(null_args_names[i]), 0)
    }
  }
  
  if(length(P)!=A || (length(S))!=A || length(B)!=A){
    stop('Lengths of P, S, & B must = A')
  }
  
  ###housekeeping and troubleshooting end###
  
  ###function for Leslie matrix assembly begin###
  
  const_leslie = function(B,S){
    
    B = B*S/2
    
    L = matrix(B,ncol = A) 
    
    for(i in 1:(A-1)){
      x = numeric(A)
      x[i] = S[i]
      L = rbind(L, x)
    }	
    
    return(L)
    
  }
  
  ###function for Leslie matrix assembly end###
  
  P_out = matrix(0, nrow = length(P),ncol = t+1)
  P_out[,1] = P
  
  ###carrying capacity begin### 
  
  Keffect = function(K,ev1,TP){ 
    
    if(is.null(K)){
      return(1)
    }else{
      return(((K+(ev1-1)*TP)/K)^-1)
    }  
  } #calculate the effect of carrying capactity at t based on population size at t-1, currently for best L1 matrix dominent eigen value
  
  ###carrying capacity end###
  
  ###hunting function start###
  
  even_spread = function(PI){
    remaining = sqrt(sum(PI[which(PI < 0)])^2)
    PI[which(PI < 0)] = 0               
    PI = PI - remaining*(PI/sum(PI))
    return(PI)
  }  
  
  if(proportion == TRUE){
    
    apply_H = function(PI,H,H_sk){
      
      basic_H = function(PI,H){
        
        PI = PI-PI*H
        return(PI)
        
      } #in cases where there is no skew hunting is applied evenly across age classes
      
      if(sum(H_sk) == 0){ 
        PI = basic_H(PI,H)
        return(PI) #with no hunting skew, basic hunting is applied across classes
      }else{
        
        PI  = PI-sum(PI)*H*H_sk
        PI = even_spread(PI)
        
        return(PI) #in cases where hunting is applied as a vector, when an age class is depleted, the quota is met with individuals from other classes
        
      }
      
    }
    
  
    }else if(proportion == FALSE){
    
    apply_H = function(PI,H,H_sk){
      
      if(sum(PI)<H){
        PI[1:length(PI)] = 0
        return(PI)
      }
      
      basic_H = function(PI,H){
        
        PI = PI-H*(PI/sum(PI))
        return(PI)
        
      } #in cases where there is no skew hunting is applied evenly across age classes
      
      PI = even_spread(PI)
      
      if(sum(H_sk) == 0){ 
        
        PI = basic_H(PI,H)
        
        return(PI) #with no hunting skew, basic hunting is applied across classes
        
      }else{
        
        PI = PI - H*H_sk
        PI = even_spread(PI)
        
        return(PI) #in cases where hunting is applied as a vector, when an age class is depleted, the quota is met with individuals from other classes
        
      }
      
    }
  }
  
  check = function(P){
    P[P < 0] = 0
    return(P)
  }
  
  ###hunting function end###
  
  ###stochastic###
  
  stocst = function(R,L){
    stoc = (1 - (rbeta(1,5,5)*2-1)*R)*L
    return(stoc)    
  } #random variation applied to L
  
  ###stochastic###
  
  ###first year start###
  
  TP = sum(P)
  L = const_leslie(B,S)
  ev = eigen(L)
  ev1 = as.numeric(ev$value[1]) #dominant eigan value (lambda), also 'long term population growth'
  
  Kef = Keffect(K,ev1,TP)
  LK = L*Kef #apply effect of carrying capacity to Leslie matrix
  LK = stocst(R,LK)
  PI = LK%*%P	 #population at t = 1 
  
  PI = apply_H(PI,H,H_sk)
  
  P_out[,2] = PI
  
  ###first year end###
  
  if(t==1){
    
    
    if(stru == TRUE){  
      
      return(P_out) 
      
    }else{
      
      yr_sum = c()
      
      for(i in 1:ncol(P_out)){
        
        N = sum(P_out[,i])
        yr_sum[i] = N
        
      }
      
      return(yr_sum)
      
    }
    
  }else{			
    
    for(i in 2:t){
      
      TP = sum(PI)
      L = const_leslie(B,S)
      ev = eigen(L)
      ev1 = as.numeric(ev$value[1]) #dominant eigan value (lambda), also 'long term population growth'
      
      Kef = Keffect(K,ev1,TP)
      LK = L*Kef #apply effect of carrying capacity to Leslie matrix

      LK = stocst(R,LK)
      
      PI = LK%*%PI	 #population at t = 1 
      
      PI = apply_H(PI,H,H_sk)
      
      P_out[,(i+1)] = PI 
      
    }
    
  }
  
  if(stru == TRUE){
    
    return(P_out) #to get all years plus age structure
    
  }else {
    
    yr_sum = c()
    
    for(i in 1:ncol(P_out)){
      
      N = sum(P_out[,i])
      yr_sum[i] = N
      
    }
    
    return(yr_sum)
    
  }
  
} #Population model, random variation dependent on stocasticity term R only

PopMod_gyby = function(t, A, K_G, K_B, K_Yr_P, P, B, S_G, S_B, Yr_P, R, H, H_sk, E, E_sk,proportion, stru){
  
  ###housekeeping and troubleshooting begin###
  
  null_args = list(H_sk,E_sk)
  null_args_names = c("H_sk","E_sk")
  for(i in 1:length(null_args)){
    if(is.null(null_args[[i]])){
      assign(paste0(null_args_names[i]), 0)
    }
  }
  
  if(length(P)!=A || (length(S_G))!=A || (length(S_B))!=A || length(B)!=A){
    stop('Lengths of P, S, & B must = A')
  }
  
  ###housekeeping and troubleshooting end###
  
  if(t==0){		
    
    if(stru == TRUE){
      
      return(P)
      
    }else{
      
      return(sum(P))
      
    }
    
  } 
  
  ###function for Leslie matrix assembly begin###
  
  const_leslie = function(B,S){
    
    B = B*S/2
    
    L = matrix(B,ncol = A) 
    
    for(i in 1:(A-1)){
      x = numeric(A)
      x[i] = S[i]
      L = rbind(L, x)
    }	
    
    return(L)
    
  }
  
  ###function for Leslie matrix assembly end###
  
  P_out = matrix(0, nrow = length(P),ncol = t+1)
  P_out[,1] = P
  
  ###carrying capacity begin### 
  
  Keffect = function(K,ev1,TP){ 
    
    if(is.null(K)){
      return(1)
    }else{
      return(((K+(ev1-1)*TP)/K)^-1)
    }  
  } #calculate the effect of carrying capactity at t based on population size at t-1, currently for best L1 matrix dominent eigen value
  
  ###carrying capacity end###
  
  ###hunting and entanglement function start###
  
  even_spread = function(PI){
    remaining = sqrt(sum(PI[which(PI < 0)])^2)
    PI[which(PI < 0)] = 0               
    PI = PI - remaining*(PI/sum(PI))
    return(PI)
  }  
  
  if(proportion == TRUE){
    
    apply_H = function(PI,H,H_sk){
      
      basic_H = function(PI,H){
        
        PI = PI-PI*H
        return(PI)
        
      } #in cases where there is no skew hunting is applied evenly across age classes
      
      if(sum(H_sk) == 0){ 
        PI = basic_H(PI,H)
        return(PI) #with no hunting skew, basic hunting is applied across classes
      }else{
        
        PI = PI = PI-sum(PI)*H*H_sk
        PI = even_spread(PI)
        
        return(PI) #in cases where hunting is applied as a vector, when an age class is depleted, the quota is met with individuals from other classes
        
      }
      
    }
    
  }else if(proportion == FALSE){
    
    apply_H = function(PI,H,H_sk){
      
      if(sum(PI)<H){
        PI[1:length(PI)] = 0
        return(PI)
      }
      
      basic_H = function(PI,H){
        
        PI = PI-H*(PI/sum(PI))
        return(PI)
        
      } #in cases where there is no skew hunting is applied evenly across age classes
      
      PI = even_spread(PI)
      
      if(sum(H_sk) == 0){ 
        
        PI = basic_H(PI,H)
        
        return(PI) #with no hunting skew, basic hunting is applied across classes
        
      }else{
        
        PI = PI - H*H_sk
        PI = even_spread(PI)
        
        return(PI) #in cases where hunting is applied as a vector, when an age class is depleted, the quota is met with individuals from other classes
        
      }
      
    }
  }
  
  check = function(P){
    P[P < 0] = 0
    return(P)
  }
  
  ###hunting and entanglement function end###
  
  ###good year bad year selection function begin###
  
  gyby = function(S_G,S_B,Yr_P){
    cho = runif(1)
    if(cho > Yr_P){
      return(S_G)
    }else{
      return(S_B)
    }
  } 
  
  ###good year bad year selection function end###
  
  stocst = function(R,L){
    stoc = (1 - (rbeta(1,5,5)*2-1)*R)*L
    return(stoc)    
  } #random variation applied to L
  
  ###first year start###
  
  S_yr = gyby(S_G,S_B,Yr_P) #gyby survuval
  K_yr = gyby(K_G,K_B,K_Yr_P) #gyby carrying capacity
  TP = sum(P)
  L = const_leslie(B,S_yr)
  ev = eigen(L)
  ev1 = as.numeric(ev$value[1]) #dominant eigan value (lambda), also 'long term population growth'
  
  Kef = Keffect(K_yr,ev1,TP)
  LK = L*Kef #apply effect of carrying capacity to Leslie matrix
  LK = stocst(R,LK)
  PI = LK%*%P	 #population at t = 1 
  
  PI = apply_H(PI,H,H_sk)
  
  P_out[,2] = PI
  
  ###first year end###
  
  if(t==1){
    
    if(stru == TRUE){  
      
      return(P_out) 
      
    }else{
      
      yr_sum = c()
      
      for(i in 1:ncol(P_out)){
        
        N = sum(P_out[,i])
        yr_sum[i] = N
        
      }
      
      return(yr_sum)
      
    }
    
  }else{			
    
    for(i in 2:t){
      
      S_yr = gyby(S_G,S_B,Yr_P)
      K_yr = gyby(K_G,K_B,K_Yr_P)
      
      TP = sum(PI)
      L = const_leslie(B,S_yr)
      ev = eigen(L)
      ev1 = as.numeric(ev$value[1]) #dominant eigan value (lambda), also 'long term population growth'
      
      Kef = Keffect(K_yr,ev1,TP)
      LK = L*Kef #apply effect of carrying capacity to Leslie matrix
   
      LK = stocst(R,LK)
      
      PI = LK%*%PI	 #population at t = 1 
      
      PI = apply_H(PI,H,H_sk)
      
      P_out[,(i+1)] = PI 
      
    }
    
  }
  
  if(stru == TRUE){
    
    return(P_out) #to get all years plus age structure
    
  }else {
    
    yr_sum = c()
    
    for(i in 1:ncol(P_out)){
      
      N = sum(P_out[,i])
      yr_sum[i] = N
      
    }
    
    return(yr_sum)
    
  }
  
} #Population model, stocastic good year bad year dynamics
