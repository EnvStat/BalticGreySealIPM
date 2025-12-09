##### Read in ####
library(rstan)
source(paste(caroll_code_folder,"1_Ordered_scripts_MOD.R", sep=""),echo = TRUE)
source(paste(code_folder,"readin.R",sep="")) # read in stan fit, samples, observed data

fit=readRDS("~/milena/grey_seal_stan_fit/fit_long_num.rds")

samples <- rstan::extract(fit, pars=c('N','N_tot',"n0",
                                      "D",
                                      'b0max', "b0min_sc",'b0min','alpha',"w", 'beta',"alpha_sc",
                                      "b", "p", "b0_av",
                                      "theta0_sc",'theta0','theta1',
                                      "mu_m",
                                      "phi_a", "phi_pup", "phi_sc", "c","v0", 'S_diag',
                                      "mu", "r",   
                                      "rho",
                                      "g_sw", "g_fi","g_bc","bK",
                                      "H_sw", "H_fi",
                                      "H_sw_tot","H_fi_tot",
                                      "K", 
                                      "pi_s", "pi_c", "kappa"
))

P03=structure(F3,S3,P01)
P03.5=c(P03[1:5],sum(P03[6:A]))

P.f=rowMeans(apply(samples$N, c(2,3), mean))[1:6]
P.f=P.f*sum(P03.5)/sum(P.f)

P.stan=c(P.f[1:5], P.f[6]*P03[6:A]/sum(P03[6:A]))
E_Range = c(0) #increased entanglement ,0.01,0.02,0.03

S.f=apply(samples$S_diag, 2, mean)[1:6]
S.stan=c(S.f, rep(S.f[6], A-6))

F.stan=c(rep(0,5),rep(mean(samples$b), A-5))

K.stan=round(mean(samples$K)*0.59)


A=47
K=K.stan
F=F.stan
P0=P.stan
S=S.stan

K=K3
F=F3
S=S3
P0=P03

rep=50
len=150

ss=200

idx <- 1:1000
seq(from=ss/rep, to=ss, length=rep) # indices for posterior samples utilized in the simulation


years.sim=2024:(2024+len-1)
N_start=50000
#### Pop size ####
###### Caroll model with Bayesian post pred params #####

P0=structure(F,S,N_start*0.59)
rep=1000

N1_tot.sim=array(NA, dim=c(rep, len))
N1.sim=array(NA, dim=c(rep, A,len))
for (j in 1:rep){
  P = PopMod(t = 149, 
             A = A, 
             K = 120000*0.59, 
             P = P0, 
             B = F, 
             S = S, 
             R = R, 
             H = 0, 
             H_sk = H_sk,
             E = 0, 
             E_sk = E_sk, 
             proportion = FALSE, 
             stru = TRUE)
  N1_tot.sim[j, ]=colSums(P)
  N1.sim[j,,]=P
}

###### Bayesian model ####
source("~/milena/forecast_function.R")

N2_tot.sim=array(NA, dim=c(length(idx), len))
N2.sim=array(NA, dim=c(length(idx), len,12))

str.stan=apply(samples$N[,,1],2, mean)
str.stan=str.stan/sum(str.stan)

N02=N_start*str.stan

for(j in 1:length(idx)) {
  
  
  params.i <- lapply(samples, function(x) {
    if(ncol(as.matrix(x)) > 1) {x[idx[j],]} else {x[idx[j]]}
  })

 results <- simulate(N0 = N02, years.sim = years.sim, 
                     fixed_K=1,
                               bias_sw=bias_sw,bias_fi=bias_fi,q=1,
                               Q_fi.sim =  rep(0, length(years.sim)),
                               Q_sw.sim = rep(0, length(years.sim)),
                               h_S.sim=rep(0,length(years.sim)),
                               h_G.sim=rep(0,length(years.sim)), 
                               params = params.i, data = data)
  
  N2_tot.sim[j,]=rowSums(results[,2:13])
  N2.sim[j,,]=results[,2:13]

}


###### Quantiles and plot ####

qs=c(0.025,0.5,0.975)

N_tot1.q=as.data.frame(t(apply(N1_tot.sim/0.59, 2, quantile, qs)))
N_tot2.q=as.data.frame(t(apply(N2_tot.sim, 2, quantile, qs)))

t=seq(from=1, to=len)
colnames(N_tot1.q)<- c("lb", "m", "ub")
N_tot1.q$t=t
N_tot1.q$mod="Carrol et al., 2024"
N_tot1.q$lb=NA
N_tot1.q$ub=NA

colnames(N_tot2.q)<- c("lb", "m", "ub")
N_tot2.q$t=t
N_tot2.q$mod="This study"
font_size=20

model.color<- c("Carrol et al., 2024"="magenta3", "This study"="green4")
N_tot.q=rbind(N_tot1.q, N_tot2.q)
legend.title="Models"
p=ggplot(N_tot.q[N_tot.q$t<80,], aes(x=t,  fill=mod))+geom_line(aes(y=m, col=mod))+
  geom_ribbon(aes( ymin=lb, ymax=ub), alpha=0.1)+
  scale_fill_manual(legend.title, values=model.color, guide="none")+
  scale_color_manual(legend.title, values=model.color)+ 
  labs(x='Year', y='Total population size') +
  guides(color = guide_legend(title.position = "top"))+
  theme_classic()+
  theme(legend.position = c(0.6,0.3), text=element_text(size=20))
 
p

ggsave(filename= paste(figure_folder,"comp_pop_size_fixedK_Bci.png",sep=""),
       plot=p,   width=5, height=5)

#### Growth rate ####

rep1=500
# G1.sim=array(NA, dim=c(3,len-1))
G1.sim={}
N1.temp=array(NA, dim=c(rep1,2))
N_start.g=1000
N.grid=seq(from=N_start.g, to=130000, length=len-1)

rep1=10
for (i in 1:(len-1)){
  P0=structure(F,S,N.grid[i]*0.59)
  #for (j in 1:rep1){
  P = PopMod(t = 1, 
               A = A, 
               K = K, 
               P = P0, 
               B = F, 
               S = S, 
               R = 0, 
               H = 0, 
               H_sk = H_sk,
               E = 0, 
               E_sk = E_sk, 
               proportion = FALSE, 
               stru = FALSE)
  # N1.temp[j,]=P
  #}
  #N1.temp.q=apply(N1.temp,2,quantile, qs)
  #G1.sim[,i]=(N1.temp.q[,2]/N1.temp.q[,1]-1)*100
  G1.sim[i]=(P[2]/P[1]-1)*100
  
}

rep=100


G2.sim=array(NA, dim=c(rep,len-1))
for(i in 1:rep) {
    params.i <- lapply(samples, function(x) {
      if(ncol(as.matrix(x)) > 1) {
        x[idx[i],]} else {x[idx[i]]}
    })
    
    for (j in 1:(len-1)){
      G2.sim[i,j]<-simulate_growth(n_start=N.grid[j],fixed_K=1,
                                   params = params.i, data = data)
      
    
  }
}

G2.sim=G2.sim[1:10,]


G1=data.frame("lb"=rep(NA, length(G1.sim)), "m"=(G1.sim))
G1$ub=NA

G1$N=N.grid
G1$mod="Carrol et al., 2024"


G2=as.data.frame(t(apply(G2.sim, 2, quantile, qs)))


colnames(G2)<- c("lb", "m", "ub")
G2$N=N.grid
G2$mod="This study"


p2=ggplot(rbind(G1,G2), aes(x=N,  fill=mod))+geom_line(aes(y=m, col=mod))+
  geom_ribbon(aes( ymin=lb, ymax=ub), alpha=0.1)+
  scale_fill_manual(legend.title, values=model.color, guide="none")+
  scale_color_manual(legend.title, values=model.color, guide="none")+ 
  labs(x='Population size', y='Growth rate') +
  geom_hline(aes(yintercept=0), linetype='dashed') + 
  
#  guides(color = guide_legend(title.position = "top"))+
  theme_classic()+
  xlim(0,125000)+
  ylim(-4,9)+
  theme(legend.position = c(0.2,0.6),  text=element_text(size=20))

p2


ggsave(filename= paste(figure_folder,"comp_growth.png",sep=""),
       plot=p2,   width=5, height=5)

