#Daire Carroll University of Gothenburg 2022
#Grey seal PVA script, function to estimate generation time based on Jonasson 2022 DOI: 10.1086/719667

generationtime = function(F,S){

  n = length(F)
  B = F*S/2
  L = matrix(B,ncol = n)
  for(i in 1:(n-1)){
    x = numeric(n)
    x[i] = S[i]
      L = rbind(L, x)
  }
  ev = eigen(L)
  lambda = as.numeric(ev$value[1])

  p = rep(NA,n)
  
  p[1]=F[1]
  p[2:n]=S[1:n-1]*F[2:n]
  R0=sum(p)
  j=c(1:n)
  v1=j*p
  v2=lambda^(-j)
  A=sum(v1*v2)/sum(v2*p)
  mu=sum(j*p)/R0
  T=log(R0)/log(lambda)
  return(data.frame(A,T,mu,lambda))
}



