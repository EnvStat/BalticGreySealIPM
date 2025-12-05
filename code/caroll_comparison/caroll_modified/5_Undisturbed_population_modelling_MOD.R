#Daire Carroll University of Gothenburg 2022
#Grey seal PVA script 5, Estimating demographic rates for a theoretical undisturbed  population (no hunting)

#look at the distribution of hunting and entanglement within age classes

setwd(paste(caroll_data_folder, "/Skold_data", sep=""))

####Entanglment###

my.data = read.csv("age_structures.csv", sep = ",", header = TRUE)
attach(my.data)
head(my.data)

my.data = my.data[which(my.data$source == "bycaught"),] #first breaking down the age structure of entangled animals

age_count_e = my.data %>% 
  group_by(age) %>%
  summarise(no_rows = length(age))

age_count_e= data.frame(age_count_e)
names(age_count_e) = c("Age", "Count")
Proportion = age_count_e$Count/sum(age_count_e$Count)
age_count_e = cbind(age_count_e, Proportion)

n = data.frame(NA,NA,NA)
names(n) = names(age_count_e)
for(i in 0:46){
  if(length(which(age_count_e$Age == i)) > 0){
    n = rbind(n,age_count_e[(which(age_count_e$Age == i)),])
  }else{
    n = rbind(n,c(i,0,0))
  }
}
age_count_e = n
age_count_e = na.omit(age_count_e)

age_count_e_pl = age_count_e
age_count_e_pl[,1] = age_count_e_pl[,1]  + 1

ggplot(age_count_e_pl, aes(Age)) +	
  
  labs(x="Age class", y=expression("Proportion")) +
  cleanup_grid + cleanup_text +
  #  geom_line(data = bind, aes(y = bind[,4], x = Year), size = 1.3, color = "blue")+
  
  #geom_vline(xintercept = 4.5, linetype = "dashed")+
  geom_point(aes(y = Proportion), size = 3)+
  xlim(1,50) #plot of proportion of sample entangled per age class

###Hunting###

my.data = read.csv("hunt_structure.csv", sep = ",", header = TRUE)
attach(my.data)
head(my.data)

age_count_h = my.data %>% 
  group_by(age) %>%
  summarise(no_rows = length(age))  #breaking down the age structure of hunted animals

age_count_h= data.frame(age_count_h)
names(age_count_h) = c("Age", "Count")
Proportion = age_count_h$Count/sum(age_count_h$Count)
age_count_h = cbind(age_count_h, Proportion)

n = data.frame(NA,NA,NA)
names(n) = names(age_count_h)
for(i in 0:46){
  if(length(which(age_count_h$Age == i)) > 0){
    n = rbind(n,age_count_h[(which(age_count_h$Age == i)),])
  }else{
    n = rbind(n,c(i,0,0))
  }
}
age_count_h = n
age_count_h = na.omit(age_count_h)

age_count_h_pl = age_count_h
age_count_h_pl[,1] = age_count_h_pl[,1]  + 1

ggplot(age_count_h_pl, aes(Age)) +	
  
  labs(x="Age class", y=expression("Proportion")) +
  cleanup_grid + cleanup_text +

  #geom_vline(xintercept = 4.5, linetype = "dashed")+
  geom_point(aes(y = Proportion), size = 3)+
  xlim(1,50)

my.data = read.csv("grey_seal_hunting.csv", sep = ",", header = TRUE)
attach(my.data)
head(my.data)

ggplot(my.data, aes(year)) +	
  
  labs(x="Year", y=expression("Number of hunted seals")) +
  cleanup_grid + cleanup_text +
  geom_point(aes(y = Total), size = 3) #hunted grey seals in Baltic per year

hunter_prop = data.frame(cbind(c(rev(Estimate),NA),c(rev(Year),2021),c(NA,my.data$Total),c(NA,my.data$Total)/c(rev(Estimate),NA))) #hunting as a proportion of estimated popualiton size based on areial surveys

ggplot(hunter_prop, aes(X2)) +	
  
  labs(x="Year", y=expression("Hunted proportion of population")) +
  cleanup_grid + cleanup_text +
  geom_point(aes(y = X4), size = 3) #hunted grey seals in Baltic per year expressed as a proportion


n = data.frame(matrix(ncol = (length(year)-8),nrow = length(age_count_h[,3])))

for(i in 1:ncol(n)){
  n[,i] = Total[i+8]*age_count_h[,3]
}
names(n) = year[9:(length(year))]

n = round(n) #these are the seals (M + F) that have been removed every year thanks to hunting in each age class `(estimate)`

#What would an undisturbed population look like?

Females_hunted = n*0.43 #assume 43% of hunted animals are female

Dummy = Females_hunted*0 #empty matrix

PopMod_Hunt = function(t, A, K, P, B, S, R, Hunt){
  
  P = structure(B,S,P)
  
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
  
  P_out = matrix(0, nrow = length(P),ncol = t+1)
  P_out[,1] = P
  
  Keffect = function(K,ev1,TP){ 
    
    if(is.null(K)){
      return(1)
    }else{
      return((1+((ev1-1)/K)*TP)^-1)
    }  
  } 
  
  stocst = function(R,L){
    stoc = (1 - (rbeta(1,5,5)*2-1)*R)*L
    return(stoc)    
  } #random variation applied to L
  
  TP = sum(P)
  L = const_leslie(B,S)
  ev = eigen(L)
  ev1 = as.numeric(ev$value[1]) 
  
  Kef = Keffect(K,ev1,TP)
  LK = L*Kef 
  LK = stocst(R,LK)
  
  PI = LK%*%P + LK%*%Hunt[,1]
  
  P_out[,2] = PI
  
  for(i in 2:t){
    
    TP = sum(PI)
    L = const_leslie(B,S)
    ev = eigen(L)
    ev1 = as.numeric(ev$value[1]) #dominant eigan value (lambda), also 'long term population growth'
    
    Kef = Keffect(K,ev1,TP)
    LK = L*Kef #apply effect of carrying capacity to Leslie matrix
    LK = stocst(R,LK)
    PI = LK%*%PI + LK%*%Hunt[,i]	 #population at t = 1 
    
    P_out[,(i+1)] = PI 
    
  }
  
  return(P_out) #to get all years plus age structure
  
} #simple population model, hunt only

#and now trying to retroactively calculate the undisturbed survival

Est_pop = function(K,F,S,P,Hunt){
  
  Undis = matrix(nrow = ncol(Females_hunted[,1:13]), ncol = rep)
  
  #P = structure(F,S,P)
  
  for(i in 1:rep){
    
    Est = PopMod_Hunt(t= t,
                      A = A,
                      K = K,
                      P = P,
                      B = F,
                      S = S,
                      R = R,
                      Hunt = Hunt)
    
    Undis[,i] = colSums(Est)
    
  }
  
  ave = rowSums(Undis)/rep
  er = rowSums(sqrt((Undis  - ave)^2))/rep
  
  return(cbind(ave,er))
  
}

rep = 1
R = 0 #deterministic model
t = 12 

mod_S = function(Sfn,Red){
  S = Sfn
  S[1] = S[1] + Red[1]*(0.21)/1
  S[2:4] = S[2:4] + Red[1]*(0.27)/3
  S[5:length(S)] = S[5:length(S)] + Red[1]*(0.52)/10 #length(S[5:length(S)])
  return(S)
} #function to modify survival values

costfn = function(Red){
  
  S = mod_S(Sfn,Red)
  
  E1 =  Est_pop(K = Kfn, F = Ffn, S = S, P = Pfn, Hunt = Dummy)
  
  E2 =  Est_pop(K = Kfn, F = Ffn, S = Sfn, P = Pfn, Hunt = Females_hunted)
  
  return(sqrt(sum((E1[,1]-E2[,1])^2)))
  
} #cost function to find absolute difference between population estimates given modified S values

Kfn = K1/0.59
Ffn = F1
Sfn = S1
Pfn = P2009_1/0.59
o1 = optim(c(0.1), costfn) #0.2256055

Kfn = K2/0.59
Ffn = F2
Sfn = S2
Pfn = P2009_2/0.59
o2 =  optim(c(0.1), costfn) #0.2336328

Kfn = K3/0.59     ##Added by Vanko et al.
Ffn = F3
Sfn = S3
Pfn = P2009_3/0.59
o3 =  optim(c(0.1), costfn) #0.2336328

S_1_Hunt = mod_S(S1     , o1$par)
S_2_Hunt = mod_S(S2     , o2$par)#estimated survivals (undisturbed population)
S_3_Hunt = mod_S(S3     , o3$par) ##Added by Vanko et al.

rep = 100
R = 0.053

Undisturbed_M = matrix(nrow = length(Year)+1, ncol = rep)
Undisturbed_Re = matrix(nrow = length(Year)+1,ncol = rep)
Disturbed = matrix(nrow = length(Year)+1,ncol = rep)

for(i in 1:rep){
  
  Undisturbed_M[10:nrow(Undisturbed_M),i] = colSums( PopMod_Hunt(t = t, 
                                                                 A = A, 
                                                                 K = K1*0.59, 
                                                                 P = P2009_1*0.59, 
                                                                 B = F1, 
                                                                 S = S_1_Hunt, 
                                                                 R = R, 
                                                                 Hunt = Dummy) )
  
  Disturbed[10:nrow(Disturbed),i] = colSums( PopMod_Hunt(t = t, 
                                                         A = A, 
                                                         K = K1*0.59, 
                                                         P = P2009_1*0.59, 
                                                         B = F1, 
                                                         S = S1, 
                                                         R = R, 
                                                         Hunt = Dummy) )
  
  
  Undisturbed_Re[10:nrow(Disturbed),i] = colSums( PopMod_Hunt(t = t, 
                                                              A = A, 
                                                              K = K1*0.59, 
                                                              P = P2009_1*0.59, 
                                                              B = F1, 
                                                              S = S1, 
                                                              R = R, 
                                                              Hunt = Females_hunted) )
  
}

erUM = 1.96*( sqrt(rowSums((Undisturbed_M-(rowSums(Undisturbed_M)/rep))^2)/rep) )
erRe = 1.96*( sqrt(rowSums((Undisturbed_Re-(rowSums(Undisturbed_Re)/rep))^2)/rep) )
erD = 1.96*( sqrt(rowSums((Disturbed-(rowSums(Disturbed)/rep))^2)/rep) ) #CI 95%

outcome = data.frame(cbind(c(rev(Year),2021),c(rev(Estimate)*0.59,NA),rowSums(Undisturbed_M)/rep,rowSums(Undisturbed_Re)/rep,rowSums(Disturbed )/rep))
names(outcome) = c("Year","Estimate","M_Undis","Re_Undis","Dis")

outcome = cbind(outcome,(sqrt((outcome$Estimate - outcome$Dis)^2)/outcome$Estimate))

names(outcome) = c("Year","Estimate","M_Undis","Re_Undis","Dis","Proportional_Difference")

outcome = cbind(outcome,erUM,erRe,erD) 

mean(outcome$Proportional_Difference, na.rm = TRUE) #this is one option for R

outcome  = data.frame(outcome[which(is.na(outcome$M_Undis)==FALSE),])

ggplot(outcome  , aes(Year)) +	
  labs(x="Year", y=expression("Estimated population size")) +
  cleanup_grid + cleanup_text   +
  
  geom_ribbon(aes(ymin= M_Undis-erUM, ymax = M_Undis+erUM), fill = "red", alpha = 0.1)+
  geom_ribbon(aes(ymin= Dis-erD, ymax = Dis+erD), fill = "grey", alpha = 0.5)+
  
  geom_line(aes(y = M_Undis), size = 1.3, color = "red", alpha = 0.5, linetype = "dashed") +
  geom_line(aes(y = Dis), size = 1.3, color = "black", alpha = 0.5) +
  
  geom_point(aes(y = Estimate), size = 2, color = "black", alpha = 1) +
  
  ylim(0,90000)+
  scale_x_continuous(breaks = seq(2008,2020,by = 2))


###R-squared of model vs. observed#####################################

R2_undis = 1 -  sum((outcome$Estimate[1:12] - outcome$Dis[1:12])^2)/sum((outcome$Dis[1:12] - mean(outcome$Dis[1:12]))^2)  

R = 0
hold = c(P2009_2*0.59,rep(NA,11))
for(i in 1:11){
  
  Est = PopMod(t = 1, 
              A = A, 
              K = K2*0.59, 
              P = structure(F2,S2,hold[i]), 
              B = F2, 
              S = S_2_Hunt, 
              R = 0, 
              H = sum(Females_hunted[,i]), 
              H_sk = H_sk,
              proportion = FALSE, 
              E = 0, 
              E_sk = E_sk, 
              stru = FALSE)
  
  hold[i+1] = Est[2]
  print(Est)

}

test = data.frame(cbind(outcome$Year[1:12],outcome$Estimate[1:12],hold))

ggplot(test  , aes(V1)) +	
  labs(x="Year", y=expression("Estimated population size")) +
  cleanup_grid + cleanup_text   +
  geom_line(aes(y = hold*2), size = 1.3, color = "black", alpha = 0.5) +
  
  geom_point(aes(y = V2*2), size = 2, color = "black", alpha = 1) +
  
  ylim(0,90000)+
  scale_x_continuous(breaks = seq(2008,2020,by = 2))

R2_dis = 1 -  sum((test$V2[1:12] - test$hold[1:12])^2)/sum((test$hold[1:12] - mean(test$hold[1:12]))^2)  

#####Plotting hunting skew#####################################################

F = F1
S = S_1_Hunt

H_sk_p = age_count_h_pl$Proportion
H_sk_p = c(H_sk,(1-sum(H_sk))*(structure(F,S,1)[(length(H_sk)+1):A]/sum((structure(F,S,1)[(length(H_sk)+1):A]))))
H_sk_a = c(0,0,0,0,(1)*(structure(F,S,1)[(length(H_sk)):A]/sum((structure(F,S,1)[(length(H_sk)):A]))))
comp = data.frame(cbind(c(1:A),
                        H_sk_p, 
                        structure(F,S,1),
                        H_sk_a))
names(comp) = c("Age","Pup","Nutral","Adult")
ggplot(comp, aes(Age)) +	
  labs(x="Age class", y=expression("Proportion of hunt")) +
  cleanup_grid + cleanup_text   +
  geom_line(aes(y = Pup),linetype = "solid", size = 1, color = "#000000", alpha = 1)+
  #geom_line(aes(y = Nutral),linetype = "dashed", size = 1, color = "#000000", alpha = 1)+
  geom_line(aes(y = Adult),linetype = "dotted", size = 1, color = "#000000", alpha = 1) 



