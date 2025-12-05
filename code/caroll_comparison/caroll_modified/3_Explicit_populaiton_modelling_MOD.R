#Daire Carroll University of Gothenburg 2022
#Grey seal PVA script 3, Fitting explicit form of the Beverton-Holt model to aerial survey data to explore long term grwoth rate under different carrying capacity scenarios

setwd(caroll_data_folder)
my.data = read.csv("Areial_Counts.csv", sep = ",", header = TRUE, fileEncoding = 'UTF-8-BOM')
attach(my.data)
head(my.data)
y = rev(Estimate)[10:length(Year)]
x = c(0:(length(y)-1))

#parametric functions and models

bev_exp = function(K,P0,lam,t){
  
  (K*P0)/(P0 + (K - P0)*lam^-t)
  
} #explicit form of the Beverton-Holt equation

f=expression((K*P0)/(P0 + (K - P0)*lam^-t))
D(f,"t")
dy.dt = function(K,P0,lam,t){
  (K * P0) * ((K - P0) * (lam^-t * log(lam)))/(P0 + (K - P0) * lam^-t)^2
} #rate of change

Rp = function(K,P0,lam,t){
  
  R = dy.dt(K,P0,lam,t)
  y = bev_exp(K,P0,lam,t)
  return(R/y)
  
} #proportional rate of change

dd_contribution = function(K, P0, lam, t){
  
  P1 = P0
  
  hold = c(P1)
  M = (1+((lam-1)/K)*P1)
  contribution = c()
  
  for(i in 2:length(t)){
    
    M = (1+((lam-1)/K)*P1)
    contribution[i] = M
    P1 = P1*lam*(M^-1)
    hold[i] = P1
    
  }
  
  return(contribution)
  
} #the effect of density dependence at t, can be thought of as the effective fraction of lam (the prpoortional rate of change) 

#two model scenarios without density dependence

nls1 = nls(y ~ m*x + c,  #line
           start = list(m = 2000, c = 16810))
nls2 = nls(y ~ a*lam^x, #exponential
           start = list(lam = 1.057, a = 16810))

#four model scenarios with density dependence (beverton-hold) with difference carrying capacities

K1 = 100000
nls3 = nls(y ~ (K1*P0)/(P0 + (K1 - P0)*(lam^-x)), 
           start = list(P0 = 15810, lam = 1.1))

K2 = 88000
nls4 = nls(y ~ (K2*P0)/(P0 + (K2 - P0)*(lam^-x)), 
           start = list(P0 = 15810, lam = 1.1))

K3 = 120000 ##Added by Vanko et al.
nls5 = nls(y ~ (K3*P0)/(P0 + (K3 - P0)*(lam^-x)), ##Added by Vanko et al.
           start = list(P0 = 15810, lam = 1.1))

est1 = c(rep(NA,length(Year)-length(x)), bev_exp(K1,coef(nls3)[1],coef(nls3)[2],x))
est2 = c(rep(NA,length(Year)-length(x)), bev_exp(K2,coef(nls4)[1],coef(nls4)[2],x))
est3 = c(rep(NA,length(Year)-length(x)), bev_exp(K3,coef(nls5)[1],coef(nls5)[2],x)) ##Added by Vanko et al.

estimates = data.frame(
  cbind(rev(Year),est1,est2,est3 ##Added by Vanko et al.
  )
) #estimates for the four carrying capacity scenarios


Rp1 = Rp(K1,coef(nls3)[1],coef(nls3)[2],x)
Rp2 = Rp(K2,coef(nls4)[1],coef(nls4)[2],x)
Rp3 = Rp(K3,coef(nls5)[1],coef(nls5)[2],x) ##Added by Vanko et al.

Prop_rates = data.frame(
  cbind(rev(Year)[10:length(Year)],Rp1,Rp2,Rp3
  )
) #estimates of proportional annual rate of growth

Q1 = dd_contribution(K1,coef(nls3)[1],coef(nls3)[2],x)
Q2 = dd_contribution(K2,coef(nls4)[1],coef(nls4)[2],x)
Q3 = dd_contribution(K3,coef(nls5)[1],coef(nls5)[2],x)  ##Added by Vanko et al.

dd_effect = data.frame(
  cbind(rev(Year)[10:length(Year)],Q1^-1,Q2^-1,Q3^-1  ##Added by Vanko et al.
  )
) #the impact of density dependence on growth for each carrying capacity scenario
names(dd_effect) = c("V1","Q1","Q2", "Q3") 


P2009_1  = coef(nls3)[1]
P2009_2  = coef(nls4)[1]
P2009_3  = coef(nls5)[1]  ##Added by Vanko et al.

lam1=   coef(nls3)[2]
lam2  = coef(nls4)[2]
lam3  = coef(nls5)[2]  ##Added by Vanko et al.


Q_m_1 = mean(na.omit(Q1)) 
Q_m_2 = mean(na.omit(Q2))
Q_m_3 = mean(na.omit(Q3))  ##Added by Vanko et al.

R1 = mean(sqrt((y - c(na.omit(estimates[,2])))^2)/na.omit(estimates[,2]))
R2 = mean(sqrt((y - c(na.omit(estimates[,3])))^2)/na.omit(estimates[,3]))
R3 = mean(sqrt((y - c(na.omit(estimates[,4])))^2)/na.omit(estimates[,4]))  ##Added by Vanko et al.

#plotting

ggplot(estimates , aes(V1)) +	
  labs(x="Year", y=expression("Population size")) +
  cleanup_grid + cleanup_text +
  geom_point(aes(y = rev(Estimate)), size = 1.3, color = "black") +
  
  geom_line(aes(y = est1), size = 1.3, color = "#88CCEE", alpha = 1)+
  geom_ribbon(aes(ymin= est1-summary(nls3)$sigma, ymax = est1+summary(nls3)$sigma), fill = "#88CCEE", alpha = 0.2)+
  
  geom_line(aes(y = est2), size = 1.3, color = "#CC6677",  alpha = 0.8)+
  geom_ribbon(aes(ymin= est2-summary(nls4)$sigma, ymax = est2+summary(nls4)$sigma), fill = "#CC6677", alpha = 0.1)+
  
  geom_hline(yintercept = K1, linetype = "dashed", color = "#88CCEE")+
  geom_hline(yintercept = K2, linetype = "dashed", color = "#CC6677")

ggplot(Prop_rates, aes(V1)) +	
  labs(x="Year", y=expression("Proportional growth rate")) +
  cleanup_grid + cleanup_text +
  
  geom_line(aes(y = Rp1), size = 1.3, color = "#88CCEE", alpha = 0.8)+
  geom_line(aes(y = Rp2), size = 1.3, color = "#CC6677", alpha = 0.8)+
  
  geom_hline(yintercept = log(coef(nls3)[2]), linetype = "dashed", color = "#88CCEE")+
  geom_hline(yintercept = log(coef(nls4)[2]), linetype = "dashed", color = "#CC6677")+
  
  ylim(0.02,0.075) #plot of proportional growth rate per year with four model scenarios (one for each K) 

ggplot(dd_effect, aes(V1)) +	
  labs(x="Year", y=expression("Density dependence modifier (Q"^-1*")",sep = "")) +
  cleanup_grid + cleanup_text +
  
  geom_line(aes(y = Q1), size = 1.3, color = "#88CCEE", alpha = 0.8)+
  geom_line(aes(y = Q2), size = 1.3, color = "#CC6677", alpha = 0.8)+
  
  ylim(0.9,1) #plot of the density dependence modifier per year with four model scenarios (one for each K) 

ggplot(dd_effect, aes(V1)) +	
  labs(x="Year", y=expression("Effective annual growth rate",sep = "")) +
  cleanup_grid + cleanup_text +
  
  geom_line(aes(y = Q1*coef(nls3)[2]), size = 1.3, color = "#88CCEE", alpha = 0.8)+
  geom_line(aes(y = Q2*coef(nls4)[2]), size = 1.3, color = "#CC6677", alpha = 0.8)+
  
  ylim(1.04,1.07) #plot of the realized annual growth rate per year with four model scenarios (one for each K) 

