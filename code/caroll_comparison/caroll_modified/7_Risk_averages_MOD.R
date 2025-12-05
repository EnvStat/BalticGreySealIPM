#Daire Carroll University of Gothenburg 2022
#Grey seal PVA script 6, estimating risks and averages for given scenarios (outcomes of 6_Scenarios.R)
library(reshape2)
library(rstan)


setwd(caroll_code_folder)

source("1_Ordered_scripts_MOD.R", echo = TRUE)

gt1 = generationtime(F1,S_1_Hunt)[1]
gt2 = generationtime(F2,S_2_Hunt)[1]
gt3 = generationtime(F3,S_3_Hunt)[1]  ##Added by Vanko et al.

setwd(caroll_data_folder)
my.data = read.csv("Areial_Counts.csv", sep = ",", header = TRUE, fileEncoding = 'UTF-8-BOM')
my.data2 = my.data[1:12,]
attach(my.data2)

setwd(caroll_code_folder)

histfun = function(S,K,P0,F){
  hist_pop = matrix(nrow = length(Estimate), ncol = rep)
  for(i in 1:rep){
    
    hist_pop[,i] = PopMod(t = (length(Estimate) - 1), 
                                     A = A, 
                                     K = K, 
                                     P = structure(F,S,P0), 
                                     B = F, 
                                     S = S, 
                                     R = R, 
                                     H = 0, 
                                     H_sk = H_sk,
                                     proportion = TRUE, 
                                     E = 0, 
                                     E_sk = E_sk, 
                                     stru = FALSE)
    
  } 
  return(hist_pop)
} #simulate hitorical popualiton growth

CI = function(P,rep){
  m = rowSums(P)/rep
  i96 = 1.96*( sqrt(rowSums((P-m)^2)/rep) )
  return(i96)
}

rep = 100 #1000
subs = histfun(S1,K1,(P2009_1),F1)
hist_K1_av = c(rep(NA,9), rowSums(subs)/ncol(subs), rep(NA, (t+1)))
hist_K1_ci = c(rep(NA,9), CI(subs,rep), rep(NA, (t+1)))

subs = histfun(S2,K2,(P2009_2),F2)
hist_K2_av = c(rep(NA,9), rowSums(subs)/ncol(subs), rep(NA, (t+1)))
hist_K2_ci = c(rep(NA,9), CI(subs,rep), rep(NA, (t+1)))

hists = data.frame(cbind(hist_K1_av,hist_K1_ci,hist_K2_av,hist_K2_ci))

detach(my.data2)
attach(my.data)
  
#####################################################

E_Range = c(0) 
H_sk_range = c(1) 
K_range = c(K1,K2, K3)  ##Added by Vanko et al.

gt_range = c(gt1,gt2, gt3)  ##Added by Vanko et al.
ice = c(1,2)

thresh_crit = (K1*0.3)
thresh_prec = (K1*0.7)

rsks_avs = function(dat,H_Range,gyby){
  
  store = data.frame(c(rev(Estimate)*0.59,rep(NA,t+1)),  c(rev(Year),c(2021:(2021+t))))
  names = c("Data","Year")
  sens = c()
  parms_store = c()
  
  q_prob_crit = c()#probability of being less than critical reference level at any point
  q_prob_prec = c()#probability of being less than precautionary reference level at any point
  iucn_10 = c() #90% decline in 3 generations
  iucn_70 = c() #70% decline in 3 generations
  iucn_50 = c() #50% decline in 3 generations
  
  if(gyby == TRUE){
    fer = c(0)
    K_range = c(K1,K2,K3, K.stan)
    for(i in c(1:length(H_Range))){
      H = H_Range[i]
      for(j in c(1:1)){
        fer_r = fer[j]
        for(k in c(1:length(H_sk_range))){
          H_skew = H_sk_range[k]
          for(l in c(1:length(K_range))){
            K = K_range[l]
            for(m in c(1:length(ice))){
              
              ic = ice[m]
              
              gt = gt_range[l]
              
              subs = dat[,which(dat[1,] == H & dat[2,] == fer_r & dat[3,] == H_skew & dat[4,] == ic & dat[5,] == K)]
              
              rep = ncol(subs)
              
              parms = c(H,fer_r,H_skew,ic,K)
              
              parms_store = rbind(parms_store,parms)
              
              av = rowSums(subs)/ncol(subs)
              
              er = CI(subs,ncol(subs)) 
              
              store = cbind(store, c(rep(NA,length(Year)),av[6:length(av)]), c(rep(NA,length(Year)), er[6:length(av)]))
              
              names = c(names, paste0("H",H,"fer_r",fer_r,"H_sk",H_skew,"Ice",ic,"K",K,"av",sep = "") , paste0("H",H,"fer_r",fer_r,"H_sk",H_skew,"Ice",ic,"K",K,"CI",sep = ""))
              
              names(store) = names
              
              sens = c(sens, paste0("H",H,"fer_r",fer_r,"H_sk",H_skew,"Ice",ic,"K",K,sep = ""))
              
              message(paste0("H",H,"fer_r",fer_r,"H_sk",H_skew,"Ice",ic,"K",K,sep = ""))
              
              q_crit = 0
              q_prec = 0
              iu_10 = 0
              iu_70 = 0
              iu_50 = 0
              
              for(m in c(1:ncol(subs))){
                
                cur = subs[,m][5:nrow(subs)]
                
                if(cur[length(na.omit(cur))] < thresh_crit){
                  r = 1
                }else{
                  r = 0
                }
                q_crit = q_crit + r ### 
                
                if(cur[length(na.omit(cur))] < thresh_prec){
                  r = 1
                }else{
                  r = 0
                }
                q_prec = q_prec + r
                
                if(min(na.omit(cur[1:(3*round(as.numeric(gt)))])) < P01*0.9){
                  r = 1
                }else{
                  r = 0
                }
                iu_10 = iu_10 + r
                
                if(min(na.omit(cur[1:(3*round(as.numeric(gt)))])) < P01*0.3){
                  r = 1
                }else{
                  r = 0
                }
                
                iu_70 = iu_70 + r
                
                if(min(na.omit(cur[1:(3*round(as.numeric(gt)))])) < P01*0.5){
                  r = 1
                }else{
                  r = 0
                }
                
                iu_50 = iu_50 + r
                
              }
              
              iu_10 = iu_10/rep
              iu_70 = iu_70/rep
              iu_50 = iu_50/rep
              q_crit = q_crit/rep 
              q_prec = q_prec/rep
              
              q_prob_crit = c(q_prob_crit,q_crit)
              q_prob_prec = c(q_prob_prec,q_prec)
              iucn_10 = c(iucn_10,iu_10) 
              iucn_70 = c(iucn_70,iu_70) 
              iucn_50 = c(iucn_50,iu_50) 
              
            }
          }
        }
      }
    } 
    store = cbind(store,hists)
    risks = data.frame(sens, q_prob_crit, q_prob_prec, iucn_10,iucn_70,iucn_50)
    risks = cbind(parms_store,risks)
    names(risks) = c("H","E","Skew","Ice","K","Scenario","crit_prob","prec_prob","iucn_10","iucn_70","iucn_50")
    
    newList = list("Risks" = risks, "Averages" = store)
    return(newList)
    
  }else if(gyby == FALSE){
    
    K_range = c(K1,K2,K3)
    H_lab=c("0000", H_Range[2:length(H_Range)])
    for(i in c(1:length(H_Range))){
      H = H_Range[i]
      for(j in c(1:length(E_Range))){
        E = E_Range[j]
        for(k in c(1:length(H_sk_range))){
          H_skew = H_sk_range[k]
          for(l in c(1:length(K_range))){
            fer = rbind(c(0), c(0), c(0),c(0))
            K = K_range[l]
            fer = fer[l,]
          
              
              f = 0
              print(f)
              print(H)
              print(K)
              print(H_sk)
              gt = gt_range[l]
              
              subs = dat[,which(dat[1,] == H & dat[2,] == E & dat[3,] == H_skew & dat[4,] == K & dat[5,] == f)]
              
              rep = ncol(subs)
              
              parms = c(H,E,H_skew,K,f)
              
              parms_store = rbind(parms_store,parms)
              
              av = rowSums(subs)/ncol(subs)
              
              er = CI(subs,ncol(subs))
              
              store = cbind(store, c(rep(NA,length(Year)),av[6:length(av)]), c(rep(NA,length(Year)), er[6:length(av)]))
              
              names = c(names, paste0("H",H_lab[i],"E",E,"H_sk",H_skew,"K",K,"fer",f,"av",sep = "") , paste0("H",H_lab[i],"E",E,"H_sk",H_skew,"K",K,"fer",f,"CI",sep = ""))
              
              names(store) = names
              
              sens = c(sens, paste0("H",H_lab[i],"E",E,"H_sk",H_skew,"K",K,"fer",f,sep = ""))
              
              message(paste0("H",H,"E",E,"H_sk",H_skew,"K",K,"fer",f,sep = ""))
              
              q_crit = 0
              q_prec = 0
              iu_10 = 0
              iu_70 = 0
              iu_50 = 0
              
              for(n in c(1:ncol(subs))){
               
                cur = subs[,n][7:nrow(subs)]
                
                if(cur[length(na.omit(cur))] < thresh_crit){
                  r = 1
                }else{
                  r = 0
                }
                q_crit = q_crit + r ### 
                
                if(cur[length(na.omit(cur))] < thresh_prec){
                  r = 1
                }else{
                  r = 0
                }
                q_prec = q_prec + r
                
                if(min(na.omit(cur[1:(3*round(as.numeric(gt)))])) < P01*0.9){
                  r = 1
                }else{
                  r = 0
                }
                iu_10 = iu_10 + r
                
                if(min(na.omit(cur[1:(3*round(as.numeric(gt)))])) < P01*0.3){
                  r = 1
                }else{
                  r = 0
                }
                
                iu_70 = iu_70 + r
                
                if(min(na.omit(cur[1:(3*round(as.numeric(gt)))])) < P01*0.5){
                  r = 1
                }else{
                  r = 0
                }
                
                iu_50 = iu_50 + r
                
              }
              
              iu_10 = iu_10/rep
              iu_70 = iu_70/rep
              iu_50 = iu_50/rep
              q_crit = q_crit/rep 
              q_prec = q_prec/rep
              
              q_prob_crit = c(q_prob_crit,q_crit)
              q_prob_prec = c(q_prob_prec,q_prec)
              iucn_10 = c(iucn_10,iu_10) 
              iucn_70 = c(iucn_70,iu_70) 
              iucn_50 = c(iucn_50,iu_50) 
              
            
          }
        }
      }
    } 
    
    store = cbind(store,hists)
    risks = data.frame(sens, q_prob_crit, q_prob_prec, iucn_10,iucn_70,iucn_50)
    risks = cbind(parms_store,risks)
    names(risks) = c("H","E","Skew","K","Fer","Scenario","crit_prob","prec_prob","iucn_10","iucn_70","iucn_50")
    
    newList = list("Risks" = risks, "Averages" = store)
    return(newList)
    
  }  
  
}  #estimate risks and averages for given data set of simulations (dat), under conditions of number of variable parameters (no_parms), hunting range (H_Range_P, H_Range_Q), if the simulation contains good year bad year dynamics (GYBY = TRUE/FALSE) 

setwd(caroll_code_folder)

##Added by Vanko et al.
dat = read.csv("Caroll_simulation.csv", sep = ",", header = TRUE, fileEncoding = 'UTF-8-BOM')
dat[5,which(dat[5,] == 0.333333333333333)] = 0.333
dat[5,which(dat[5,] == 0.33519 )] = 0.335
dat[4,] =round(dat[4,])

H_Range_Q = c(0,2400,3600) #seals hunted, quota


H_Range_Q = H_Range_Q*0.43 #Hunt_Ratio #0.5   females hunted, quota
Fer_Q.r = rsks_avs(dat,H_Range_Q,FALSE)


###### Here on simulation with Vanko et al. model + compariosn and plot of results with the two different models #####



setwd(code_folder)
source("forecast_function.R") #contains simualting function
source("readin.R") # read in stan fit, samples, observed data

years.sim= 2021:2100
quantiles <- c(0.025, 0.5, 0.975)


# Set up hunting bag for the 3 scenarios
Q_fi.scenarios.comp <- matrix(NA,3 , length(years.sim))
Q_fi.scenarios.comp[,1] <- 0 # known Finnish quota for 2023
Q_fi.scenarios.comp[1,] <- 0
Q_fi.scenarios.comp[2,] <- 1200
Q_fi.scenarios.comp[3,] <- 1800


Q_sw.scenarios.comp <- matrix(NA,3 , length(years.sim))
Q_sw.scenarios.comp[,1] <- 0 
Q_sw.scenarios.comp[1,] <- 0
Q_sw.scenarios.comp[2,] <- 1200
Q_sw.scenarios.comp[3,] <- 1800


bias_sw=log(colMeans(samples$rho[,1,,3]/apply(samples$rho[,1,,3],1, mean)))
bias_fi=log(colMeans(samples$rho[,1,,4]/apply(samples$rho[,1,,4],1, mean)))


n0=55997.36 #Start from same initial population in 2020 as Caroll et al in the K=120,000 case, calculated as  sum(P03)/0.59

N.sim <- array(NA, dim=c(nrow(Q_fi.scenarios.comp), length(idx), length(years.sim))) 

#run simulation
for(i in 1:length(idx)) {
  params.i <- lapply(samples, function(x) {
    if(ncol(as.matrix(x)) > 1) {x[idx[i],]} else {x[idx[i]]}
  })
  N0.stan=round(samples$N[idx[i],,22+1])
  
  for(j in 1:nrow(Q_fi.scenarios.comp)) {
    results <- simulate(N0 = n0/sum(N0.stan)*N0.stan, years.sim = years.sim,         
                          bias_sw=bias_sw,
                          bias_fi=bias_fi,
                          q=1,
                          Q_fi.sim = Q_fi.scenarios.comp[j,],
                          Q_sw.sim = Q_sw.scenarios.comp[j,],
                          h_S.sim=rep(0,length(years.sim)),
                          h_G.sim=rep(0,length(years.sim)), 
                          params = params.i, data = data)
      
      N.sim[j,i,] <- rowSums(results[,2:13])
    }
  }



# Put into dataframe and plot

N=apply(N.sim,c(1,3), quantile, quantiles)
dimnames(N)[[2]] <-c("0", "2400", "3600")
dimnames(N)[[3]]<- as.numeric(years.sim)

N.df <- data.frame(melt(N[1,,]))
colnames(N.df)<-c("Scenario",  "year", "lb")
N.df$md=melt(N[2,,])$value
N.df$ub=melt(N[3,,])$value
N.df[,1]=as.factor(N.df[,1])

N.df$K="This study"
N.df=N.df[,c("year", "Scenario", "ub", "md", "lb","K" )]

options("scipen"=100, "digits"=4)


# Caroll's results -> format to plot
Av=Fer_Q.r$Averages
Av[,3:ncol(Av)]=Av[,3:ncol(Av)]/0.59
Av.melt=melt(Av, id="Year")
Av.melt$Hunting=as.character(as.numeric(substr(Av.melt$variable,2,5))/0.43)
Av.melt$K=as.numeric(substr(Av.melt$variable,14,18))/0.59
Av.melt=Av.melt[substr(Av.melt$variable,23,24)=="av" & Av.melt$K < 125000,]
Av.melt=Av.melt[Av.melt$Year<=2100,]
Av.melt$K=paste("Caroll et al. (2024),",as.character(Av.melt$K))
Av.melt$Scenario=Av.melt$Hunting
Av.melt$year=as.numeric(Av.melt$Year)
Av.melt$lb=NA
Av.melt$md=Av.melt$value
Av.melt$ub=NA
Av.melt=Av.melt[,c("year", "Scenario", "ub", "md", "lb","K" )]


model.breaks= c("Caroll et al. (2024), 120000","Caroll et al. (2024), 100000", "Caroll et al. (2024), 88000", "This study")
model.linetype<- c("Caroll et al. (2024), 120000"=2,"Caroll et al. (2024), 100000"=3, "Caroll et al. (2024), 88000"=6, "This study"=1)
model.linewidth<- c("Caroll et al. (2024), 120000"=0.8,"Caroll et al. (2024), 100000"=0.8, "Caroll et al. (2024), 88000"=0.8, "This study"=0.8)
model.alpha<- c("Caroll et al. (2024), 120000"=0,"Caroll et al. (2024), 100000"=0, "Caroll et al. (2024), 88000"=0, "This study"=0.07)

scenario.colors <- c('Historical'='black', 
                     '0'='red', 
                     '2400'='blue', 
                     '3600'='orange', 
                     '4800'='green3',
                     "1271"="cyan3"
) #Use same colors as for analysis plots


p=ggplot(rbind(Av.melt, N.df), aes(x=year, y=md,  fill=Scenario,  linetype=K, linewidth = K))+
  geom_line(aes(col=Scenario))+
  geom_ribbon(aes(ymin=lb, ymax=ub, alpha=K),linetype=1) +
  labs(y='Population size', x='Year') +
  xlim(2015,2100)+
  theme_classic()+
  scale_linetype_manual("Model & carrying capacity", values=model.linetype, breaks=model.breaks) +
  scale_linewidth_manual("Model & carrying capacity", values=model.linewidth, breaks=model.breaks, guide="none") +
  scale_alpha_manual("", values=model.alpha, breaks=model.breaks, guide="none") +
  scale_color_manual("Hunted Seals", values=scenario.colors) +
  scale_fill_manual("Hunted Seals", values=scenario.colors, guide="none") +
  theme(legend.position = c(0.2,0.7), legend.direction = "vertical", legend.box = "vertical",text=element_text(size=20),
        legend.background = element_rect(fill="transparent" ,colour = "transparent"))
 
p


ggsave(filename= paste(figure_folder,"comp_Caroll_sim.png",sep=""),
       plot=p,   width=12, height=10)


df=rbind(Av.melt, N.df)
p=ggplot(df[df$Scenario==0,], aes(x=year, y=md,   linetype=K, linewidth = K), fill=scenario.colors[2])+
  geom_line(aes(col=Scenario))+
  geom_ribbon(aes(ymin=lb, ymax=ub, alpha=K),linetype=1, fill=scenario.colors[2]) +
  labs(y='Population size', x='Year') +
  xlim(2015,2100)+
  theme_classic()+
  scale_linetype_manual("Model & carrying capacity", values=model.linetype, breaks=model.breaks, guide="none") +
  scale_linewidth_manual("Model & carrying capacity", values=model.linewidth, breaks=model.breaks, guide="none") +
  scale_alpha_manual("", values=model.alpha, breaks=model.breaks, guide="none") +
  scale_color_manual("Hunted Seals", values=scenario.colors, guide="none") +
 # scale_fill_manual("Hunted Seals", values=scenario.colors, guide="none") +
  theme(legend.position = c(0.2,0.55), legend.direction = "vertical", legend.box = "vertical",text=element_text(size=20),
        legend.background = element_rect(fill="transparent" ,colour = "transparent"))

p

ggsave(filename= paste(figure_folder,"comp_Caroll_sim_0.png",sep=""),
       plot=p,   width=12, height=6)

p=ggplot(df[df$Scenario==2400,], aes(x=year, y=md,   linetype=K, linewidth = K), fill=scenario.colors[3])+
  geom_line(aes(col=Scenario))+
  geom_ribbon(aes(ymin=lb, ymax=ub, alpha=K),linetype=1, fill=scenario.colors[3]) +
  labs(y='Population size', x='Year') +
  xlim(2015,2100)+
  theme_classic()+
  scale_linetype_manual("Model & carrying capacity", values=model.linetype, breaks=model.breaks, guide="none") +
  scale_linewidth_manual("Model & carrying capacity", values=model.linewidth, breaks=model.breaks, guide="none") +
  scale_alpha_manual("", values=model.alpha, breaks=model.breaks, guide="none") +
  scale_color_manual("Hunted Seals", values=scenario.colors, guide="none") +
  # scale_fill_manual("Hunted Seals", values=scenario.colors, guide="none") +
  theme(legend.position = c(0.2,0.55), legend.direction = "vertical", legend.box = "vertical",text=element_text(size=20),
        legend.background = element_rect(fill="transparent" ,colour = "transparent"))

p

ggsave(filename= paste(figure_folder,"comp_Caroll_sim_2400.png",sep=""),
       plot=p,   width=12, height=6)

p=ggplot(df[df$Scenario==3600,], aes(x=year, y=md,   linetype=K, linewidth = K), fill=scenario.colors[4])+
  geom_line(aes(col=Scenario))+
  geom_ribbon(aes(ymin=lb, ymax=ub, alpha=K),linetype=1, fill=scenario.colors[4]) +
  labs(y='Population size', x='Year') +
  xlim(2015,2100)+
  theme_classic()+
  scale_linetype_manual("Model & carrying capacity", values=model.linetype, breaks=model.breaks)+
  scale_linewidth_manual("Model & carrying capacity", values=model.linewidth, breaks=model.breaks, guide="none") +
  scale_alpha_manual("", values=model.alpha, breaks=model.breaks, guide="none") +
  scale_color_manual("Hunted Seals", values=scenario.colors, guide="none") +
  # scale_fill_manual("Hunted Seals", values=scenario.colors, guide="none") +
  theme(legend.position = c(0.2,0.6), legend.direction = "vertical", legend.box = "vertical",text=element_text(size=20),
        legend.background = element_rect(fill="transparent" ,colour = "transparent"))

p

ggsave(filename= paste(figure_folder,"comp_Caroll_sim_3600.png",sep=""),
       plot=p,   width=12, height=6)
