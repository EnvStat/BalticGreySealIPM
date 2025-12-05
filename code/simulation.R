#### Read in ####

source(paste(code_folder,"readin.R",sep="")) # read in stan fit, samples, observed data

#### Set up ####

source(paste(code_folder,"forecast_function.R",sep="")) #contains simualting function
source(paste(code_folder,"simulation_setup.R",sep="")) # read in hunting quota, hunting bias and herring quality scenarios
##### Population size  and birth rate in modelled period 2003-2025 ####
N.past <- apply(samples$N[idx,,1:(t+1)], c(1,3), sum)
N.past.q<-t(apply(N.past,2, quantile, quantiles))



N.past.df <- data.frame(Scenario = rep(as.factor('Historical'), nrow(N.past.q)), 
                        Hunting_bias = rep(as.factor('Historical'), nrow(N.past.q)), 
                        Herring_quality = rep(as.factor('Historical'), nrow(N.past.q)), 
                        year = (2002):2025,
                        lb = c(N.past.q[,1]), 
                        md = c(N.past.q[,2]), 
                        ub = c(N.past.q[,3]))





b.past <- samples$b[idx,]
b.past.q<-t(apply(b.past,2, quantile, quantiles))

b.past.df <- data.frame(Scenario = rep(as.factor('Historical'), nrow(b.past.q)), 
                        year = (2002):2025,
                        lb = c(b.past.q[,1]), 
                        md = c(b.past.q[,2]), 
                        ub = c(b.past.q[,3]))
#### Simulation ####

N.sim <- array(NA, dim=c(nrow(Q_fi.scenarios), nrow(bias_fi.scenarios),nrow(h_S.scenarios), length(idx), length(years.sim))) 
b.sim <- array(NA, dim=c( nrow(h_S.scenarios), length(idx), length(years.sim.birth))) 
N.composition.sim <- array(NA, dim=c( length(years.sim.composition), 12, length(idx))) 


for(i in 1:length(idx)) {
  # extract posterior samples for parameters
  params.i <- lapply(samples, function(x) {
    if(ncol(as.matrix(x)) > 1) {x[idx[i],]} else {x[idx[i]]}
  })
  for(j in 1:nrow(Q_fi.scenarios)) {
    for(k in 1:nrow(bias_fi.scenarios)) {
      for(l in 1:nrow(h_S.scenarios)) {
        
        results <- simulate(N0 = round(samples$N[idx[i],,t+1]), years.sim = years.sim,         
                            bias_sw=bias_sw.scenarios[k,],
                            bias_fi=bias_fi.scenarios[k,],
                            Q_fi.sim = Q_fi.scenarios[j,],
                            Q_sw.sim = Q_sw.scenarios[j,],
                            h_S.sim = h_S.scenarios[l,] ,
                            h_G.sim = h_G.scenarios[l,] , 
                            params = params.i, data = data)
        N.sim[j,k,l,i,] <- rowSums(results[,2:13])
        if (j==1 & k==1){b.sim[l,i,]=results[1:length(years.sim.birth),1]}
        if(j==2 & k==1 & l==1){  N.composition.sim[,,i] <- results[1:length(years.sim.composition),2:13]}
        
      }
    }
  }
}


N=apply(N.sim,c(1,2,3,5), quantile, quantiles)  #get the quantiles from the whole simulated data
dimnames(N)[[2]] <- scenario.breaks[2:7]
dimnames(N)[[3]] <- bias.breaks
dimnames(N)[[4]] <- her.breaks
dimnames(N)[[5]]<- years.sim

N.df <- data.frame(melt(N[1,,,,])) 
colnames(N.df)<-c("Scenario", "Hunting_bias", "Herring_quality", "year", "lb")
N.df$md=melt(N[2,,,,])$value
N.df$ub=melt(N[3,,,,])$value
N.df[,1]=as.factor(N.df[,1])
N.df[,2]=as.factor(N.df[,2])

N.full.hunt=rbind(N.past.df, N.df)
N.full.hunt$group=paste(N.full.hunt$Scenario, N.full.hunt$Hunting_bias, N.full.hunt$Herring_quality)


write.csv(N.full.hunt,paste(simulations_folder, "Simulation.csv", sep="")) # Median and 95% CI of total pop size 2025-2080 for all scenarios


N.composition.sim.median=apply(N.composition.sim, c(1,2), quantile, 0.5)
write.csv(N.composition.sim.median,paste(simulations_folder, "Simulation_composition.csv",sep=""),row.names = FALSE) # Demographic group specific median and 95% CI of pop size 2025-2080 for 2024 harvest size and historcal bias and herring

b=apply(b.sim, c(1,3), quantile, quantiles)
b.df <- data.frame(Scenario = rep(as.factor(her.breaks), each=length(years.sim.birth)), 
                   year = rep(years.sim.birth, nrow(h_S.scenarios)),
                   lb = c(t(b[1,,])), 
                   md = c(t(b[2,,])), 
                   ub = c(t(b[3,,])))
write.csv(b.df,paste(simulations_folder,"Simulation_birth.csv",sep=""), row.names = FALSE) #Median and 95% CI of birth rate 2025-2080 under no hunting and varying herirng quality


N.sim.80.bias=array(NA, dim=c(length(B.list), length(idx)))
for ( j in 1:length(B.list)){
  N.sim.80.bias[j,]=N.sim[Q.list[j],B.list[j],1,,length(years.sim)]
}

N.sim.80.her=array(NA, dim=c(length(H.list), length(idx)))
for ( j in 1:length(H.list)){
  N.sim.80.her[j,]=N.sim[4,1,H.list[j],,length(years.sim)]
}

write.csv(as.data.frame(N.sim.80.bias), paste(simulations_folder,"Simulation_2080_hunting_bias.csv",sep=""), row.names = FALSE) #Posterior distribution in 2080 under historical herring and varying harvest size and bias
write.csv(as.data.frame(N.sim.80.her),paste(simulations_folder, "Simulation_2080_hunting_her.csv",sep=""), row.names = FALSE)  #Posterior distribution in 2080 under harvest size of 3600, historical bias and varying herring quality

N.sim.80.3600=apply(N.sim[4,,,,length(years.sim)]>100, c(1,2), sum)/length(idx)*100
N.sim.80.4800=apply(N.sim[5,,,,length(years.sim)]>100, c(1,2), sum)/length(idx)*100


write.csv(as.data.frame(N.sim.80.3600), paste(simulations_folder,"Simulation_2080_3600.csv",sep=""), row.names = FALSE) #Posterior distribution in 2080 under historical herring and varying harvest size and bias
write.csv(as.data.frame(N.sim.80.4800),paste(simulations_folder, "Simulation_2080_4800.csv",sep=""), row.names = FALSE)  #Posterior distribution in 2080 under harvest size of 3600, historical bias and varying herring quality
