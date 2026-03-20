
# Initial set up
# ======================
library(rstan)
library(extraDistr)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write=TRUE)

# Set the folders

# Update this to be the path to your project folder
project.folder = "/home/jpvanhat/Documents_extended/Projects_data/Baltic_seals/baltic_grey_seals/BalticGreySealIPM"
source("set_folders.R")

# Set model details 
year=seq(from=2003, to=2025, length=23) #years data covers
t=length(year) # number of modelled years
a=6 #number of age classes (0,1,2,3,4,5+)


##### READ IN DATA ########
# ======================

setwd(proc_data_folder)

#Aeraial count
df_count=read.csv("Aerial_count.csv")
y=df_count$Total
y[16]=0 # Stan does not accept NA

#Hunting bag and quota
df_hb=read.csv("Hunting_bag_quota_FI_SW.csv")
Q_sw=df_hb$Quota_SW
Q_fi=df_hb$Quota_FI
y_hb_sw=df_hb$Hunting_bag_SW
y_hb_fi=df_hb$Hunting_bag_FI

#Hunting samples
df_hs_sw=read.csv("Hunting_samples_SW.csv", header = TRUE)
df_hs_sw$X <- NULL
y_hs_sw=as.matrix(df_hs_sw)
y_hs_sw=y_hs_sw

df_hs_fi=read.csv("Hunting_samples_FI.csv", header = TRUE)
df_hs_fi$X <- NULL
y_hs_fi=as.matrix(df_hs_fi)
y_hs_fi=y_hs_fi

#Bycatch samples
df_bc_sw=read.csv("Bycatch_samples_SW.csv",check.names = FALSE)
df_bc_sw[,1]<- NULL
y_bc_sw=as.matrix(df_bc_sw)

df_bc_fi=read.csv("Bycatch_samples_FI.csv",check.names = FALSE)
df_bc_fi[,1] <- NULL
y_bc_fi=as.matrix(df_bc_fi)
y_bc=as.matrix(df_bc_sw)
y_bc<-y_bc[1:ncol(df_bc_fi)]+as.matrix(df_bc_fi) #merge together the two countries

#Preganancy samples
df_pr=read.csv("Pregnancy_samples_SW.csv")
y_pr=df_pr$Pregnant
y_pr_tot=df_pr$Total

#Reproductive sign samples (CA and placental scar)
df_rep_fi=read.csv("Reproductive_signs_FI.csv")
df_rep_fi$X<- NULL
df_rep_fi=df_rep_fi
Z_fi=as.matrix(df_rep_fi)

#Herring WAA (mean weight of herring over the age of five)
df_h=read.csv("Herring.csv")
h_BP_GoF=df_h$BP_GoF.norm  #Baltic Proper + Gulf of Finland
h_GoB=df_h$GoB.norm #Gulf of Bothnia

inits<- function(){ #Draw initial parameters from prior distributions
  u=matrix(nrow=3*12, ncol=t)
  for (i in 1:(3*12)){
    u[i,]=rnorm(t, 0,1)
  }

  return (list(  
    
    #Inital population size
    n0=rlnorm(1,9.8,0.1),
    
    #Natural mortality
    phi_a_sc=runif(1,0.8,1),
    phi_sc=runif(1,0.6,1), 
    c=runif(1,0,1),
    
    v0=rnorm(1,0,0.2),
    v5=rnorm(1,0.88,0.2),
    
    #Bias
    g_sw=rnorm(12,0,0.1),
    g_fi=rnorm(12,0,0.1),
    g_bc=rnorm(12,0,0.1),
    
    #Hunting
    sigma_h_sw=rhcauchy(1,0.1),
    sigma_h_fi=rhcauchy(1,0.1),
    
    #Birth
    b0max=runif(1,0.8,1),
    b0min_sc=runif(1,0,1),
    alpha_sc=rnorm(1,0,4),
    beta=rnorm(1,0,3),
    w=runif(1,0,1),
    theta0_sc=runif(1,0,1),
    
    #Carrying capacity
    K=rlnorm(1,11.3,0.3),
    
    #Aerial obs
    mu=rbeta(1,32,9),
    r=rlnorm(1,5.3,1), 

    #Reproductive sign obs  
    kappa=runif(1,0,1),
    pi_s_mean=runif(1,0,1),
    pi_c_mean=runif(1,0,1),
    
    #variance for scar observation   
    sigma_rep_c=rhnorm(1,0.1),
    sigma_rep_s=rhnorm(1,0.1),
    
    #Stanadard normals for stochasticities
    epsilon_h_sw=rhnorm(t,1),
    epsilon_h_fi=rhnorm(t,1),
    
    epsilon_rep_c=rhnorm(t,1),
    epsilon_rep_s=rhnorm(t,1),
    
    epsilon_birth=rnorm(t+1,0,1),
    epsilon_sex=rnorm(t+1,0,1),
    
    u=u
  ))
}


C=matrix(0,nrow=12, ncol=12)
C[2:5,2:5]=0.95
C[8:11,8:11]=0.95
diag(C)=1
L_bias=t(chol(C)) #Variance for multivariate Gaussian prior for biases


data=list(a=6,t=t, 
              L_bias=L_bias,  
              N0=rep(1/12,12)*20000, 
              y=y,
              h_BP_GoF=h_BP_GoF, h_GoB=h_GoB,
              z_fi=Z_fi, 
              Q_sw=Q_sw,
              Q_fi=Q_fi,
              y_pr=y_pr, y_pr_tot=y_pr_tot, 
              y_hb_sw=y_hb_sw, y_hb_fi=y_hb_fi,
              y_hs_sw=y_hs_sw, y_hs_fi=y_hs_fi, 
              y_bc=y_bc)

###### run above to get data ready for Stan ##### 


##### FIT STAN ####

setwd(code_folder)


stan_model <- stan_model("grey_seal.stan") 


fit <- sampling(stan_model, data = data, iter = 4000, warmup=2000, chains = 4, 
                         init =list(inits(),inits(),inits(),inits()),
                         control = list(adapt_delta=0.95, max_treedepth = 13))


##### SAVE RDS file ###

setwd(stan_fit_folder)

saveRDS(fit,"grey_seal_fit.rds")
