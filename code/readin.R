library("rstan")
library("ggplot2")
library("ggridges")
library("reshape2")
library("patchwork")
library("bayesplot")
library("reshape2")
library("pheatmap")



# Readin RDS
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



#Read in observed data

##Aerail survey count
y_obs=read.csv(paste(proc_data_folder,"Aerial_count.csv", sep=""))$Total

#Hunting bag
df_hb=read.csv((paste(proc_data_folder,"Hunting_bag_quota_FI_SW.csv", sep="")))
y_hb_sw<-df_hb$Hunting_bag_SW
y_hb_fi<-df_hb$Hunting_bag_FI

#Hunting samples SW & FI
df_hs_sw=read.csv((paste(proc_data_folder,"Hunting_samples_SW.csv", sep="")))
df_hs_sw$X<-NULL
y_hb_sw<-as.matrix(df_hs_sw)

df_hs_fi=read.csv((paste(proc_data_folder,"Hunting_samples_FI.csv", sep="")))
df_hs_fi$X<-NULL
y_hs_fi<-as.matrix(df_hs_fi)

#Bycatch samples
df_bc_sw=read.csv((paste(proc_data_folder,"Bycatch_samples_SW.csv", sep="")))
df_bc_sw$X<-NULL
y_bc_sw<-as.matrix(df_bc_sw)

df_bc_fi=read.csv((paste(proc_data_folder,"Bycatch_samples_FI.csv", sep="")))
df_bc_fi$X2024=rep(0,12)
df_bc_fi$X<-NULL

y_bc_fi<-as.matrix(df_bc_fi)
y_bc=y_bc_fi+y_bc_sw

#Pregnancy samples
df_pr=read.csv(paste(proc_data_folder,"Pregnancy_samples_SW.csv", sep=""))
y_pr=df_pr$Pregnant
y_pr_tot=df_pr$Total

#Reproductive samples(CA and placental scar)
df_rep=read.csv(paste(proc_data_folder,"Reproductive_signs_FI.csv", sep=""))
Z=as.matrix(df_rep[,2:23])


#Herring WAA (mean weight of herring over the age of five)
df_h=read.csv(paste(proc_data_folder,"Herring.csv", sep=""))
h_BP_GoF=df_h$BP_GoF.norm  #Baltic Proper + Gulf of Finland
h_GoB=df_h$GoB.norm #Gulf of Bothnia
