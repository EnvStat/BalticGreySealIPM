library("rstan")
library("ggplot2")
library("ggridges")
library("reshape2")
library("patchwork")
library("bayesplot")
library("reshape2")
library("pheatmap")



# Read in RDS
font_size=20
quantiles <- c(0.025, 0.5, 0.975) # quantiles for median and 95%CI 
idx <- (1:1000)*8 # indices for posterior samples utilized in the simulation. thinned out so it doesn't run too long

fit=readRDS(paste(stan_fit_folder,"grey_seal_fit.rds", sep=""))


# Sample readin
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
                                      "pi_s", "pi_c", "kappa",
                                      
                                      "y_sim",
                                      "y_pr_sim", "z_sim",
                                      "y_hb_sw_sim", "y_hb_fi_sim"
                                      
))

##### READ IN DATA ########
# ======================

setwd(proc_data_folder)

#Aerial count
df_count=read.csv("Aerial_count.csv")
y_obs=df_count$Total
#y_obs[16]=0 # Stan does not accept NA

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

setwd(code_folder)

