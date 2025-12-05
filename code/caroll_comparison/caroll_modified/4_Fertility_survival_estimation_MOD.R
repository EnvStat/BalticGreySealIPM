#Daire Carroll University of Gothenburg 2022
#Grey seal PVA script 4, Making initial estimates of age specific fertility

setwd(paste(caroll_data_folder, "/Skold_data", sep=""))

my.data = read.csv("reproduction.csv", sep = ",", header = TRUE)
attach(my.data)
head(my.data)

age_count = my.data %>% 
  group_by(age) %>%
  summarise(no_rows = length(age)) #calculate samples per age class

fertility_count = my.data %>% 
  group_by(age) %>%
  summarise(sum(reproduce)) #calculate number reproductive per age class

fertility_levels = cbind(data.frame(age_count),data.frame(fertility_count)[,2])
names(fertility_levels) = c("Age","Count","No_Fertile")

birthrates = fertility_levels$No_Fertile/fertility_levels$Count #estimate birthrates as a proportion reproductive of sample size

fertility_levels = cbind(fertility_levels,birthrates)

ggplot(fertility_levels, aes(Age)) +	
  geom_point(aes(y = Count), size = 3)+
  labs(x="Age", y=expression("Sample size")) +
  cleanup_grid + cleanup_text #sample size

ggplot(fertility_levels, aes(Age)) +	
  geom_point(aes(y = birthrates), size = 3)+
  labs(x="Age", y=expression("Proportion fertile")) +
  cleanup_grid + cleanup_text  #fertility rates

ct = 4 #apply cupoff based on sample size

x = fertility_levels[fertility_levels$Count > ct,][,1]
y = fertility_levels[fertility_levels$Count > ct,][,4] #apply cupoff based on sample size

x = x[3:length(x)]
y = y[3:length(y)]

logistic = function(y0,K,mu,x){
  (K * y0)/(y0 + (K - y0) * exp(-mu * x))
}

zl = nls(y ~ (K * y0)/(y0 + (K - y0) * exp(-mu * x)), 
         start = list(y0 = 1e-6, K = 0.85, mu = 1.5))

bind = data.frame(cbind(fertility_levels$Age,
                        c(0,0,0,logistic(summary(zl)$coefficients[1,1],summary(zl)$coefficients[2,1],summary(zl)$coefficients[3,1],fertility_levels$Age[4:length(fertility_levels$Age)])),
                        fertility_levels$birthrates,
                        fertility_levels$birthrates,
                        fertility_levels$Count))

names(bind) = c("Age","Log_predict","Data","Data_Omitted","Count")

bind$Data[which(bind$Count <= 4)] = NA
bind$Data_Omitted[which(bind$Count > 4)] = NA

er1 = summary(zl)$sigma

ggplot(bind, aes(Age)) +	
  geom_line(aes(y = Log_predict), size = 1.3, color = "blue") +
  geom_ribbon(aes(ymin= (Log_predict-er1), ymax = (Log_predict+er1)), fill = "skyblue1", alpha = 0.4) +
  geom_point(aes(y = (Data)),size = 2)+
  geom_point(aes(y = (Data_Omitted)),size = 2, col = "grey")+
  labs(x="Age", y=expression("Reproductive rate (f)"))+
  
  cleanup_grid +
  cleanup_text #no senescence 

#considering senescence

F26plus = sum(fertility_levels$No_Fertile [which(fertility_levels$Age>25)])/sum(fertility_levels$Count[which(fertility_levels$Age>25)])
F26plus

scen_scenaro = rep(NA, nrow(bind))
scen_scenaro[27:33] =  F26plus
subad = rep(NA, nrow(bind))
subad[0:3] =  0
bind$Log_predict[27:33] =  NA
Log_predict2 = bind$Log_predict
bind$Log_predict[1:3] =  NA

bind = cbind(bind,scen_scenaro,subad,Log_predict2 )

ggplot(bind, aes(Age)) +	
  geom_line(aes(y = Log_predict2), size = 1.3, color = "blue") +
  geom_line(aes(y = scen_scenaro), size = 1.3, color = "blue") +
  geom_line(aes(y = subad), size = 1.3, color = "blue") +
  geom_ribbon(aes(ymin= (Log_predict-er1), ymax = (Log_predict+er1)), fill = "skyblue1", alpha = 0.4) +
  geom_point(aes(y = (Data)),size = 2)+
  geom_point(aes(y = (Data_Omitted)),size = 2, col = "black")+
  labs(x="Age", y=expression("Fertility rate"))+
  geom_vline(xintercept = 25.5, size = 1.3, linetype = "dashed")+
  
  
  cleanup_grid +
  cleanup_text 

#########################################################################
#to build a predictive model including density dependence, we need to consider that density may have been acting on parameters during the observation period, if we assume that we have a mean F and mean(F)/mean(Q) = F

A = 47

Fe = c(0,0,logistic(summary(zl)$coefficients[1,1],summary(zl)$coefficients[2,1],summary(zl)$coefficients[3,1],c(3:26)),rep((F26plus),(A)-26))#Each age class is assigned fertility of next cycle

F1 = Fe*Q_m_1
F2 = Fe*Q_m_2
F3 = Fe*Q_m_3  ##Added by Vanko et al. NOTE: need to run scrit 3 for this
##################################################

P2009_1 = 31163.51 
P2009_2 = 30968.65 
P2009_3 = 31358   ##Added by Vanko et al. NOTE: this is a result from script 3 opied here
lam1 = 1.099
lam2 = 1.111 
lam3 = 1.086613  ##Added by Vanko et al.

scon = function(parm){
  St = c(parm[1],rep(0.932,3),rep(0.95,(A-4)))
  ev1 = eig_val(F,St)
  return( sqrt((ev1 - goal)^2) )
} #absolute difference between "goal" eigenvalue and eigenvalue for a leslie matrix given fertility values, F, and grouped survival values  

F = F1
goal = lam1
o1 = optim(c(0.7), scon, method = "Brent", lower = 0, upper = 0.932) #minimize scon function

F = F2
goal = lam2
o2 = optim(c(0.7), scon, method = "Brent", lower = 0, upper = 0.932) #minimize scon function

F = F3   ##Added by Vanko et al.
goal = lam3
o3 = optim(c(0.7), scon, method = "Brent", lower = 0, upper = 0.932) #minimize scon function


S1 = c(o1$par[1]  ,rep(0.932 , 3), rep(0.95  , (A-4)))
S2 = c(o2$par[1]  ,rep(0.932 , 3), rep(0.95  , (A-4))) 
S3 = c(o3$par[1]  ,rep(0.932 , 3), rep(0.95  , (A-4)))   ##Added by Vanko et al.

