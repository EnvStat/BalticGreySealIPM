#Daire Carroll University of Gothenburg 2022
#Grey seal PVA script 6, testing various hunting and environmental change scenarios


setwd(caroll_code_folder)
source("1_Ordered_scripts_MOD.R", echo = TRUE)
setwd(caroll_data_folder)
      
my.data = read.csv("Areial_Counts.csv", sep = ",", header = TRUE, fileEncoding = 'UTF-8-BOM')
ice.data = read.csv("Ice.csv", sep = ",", header = TRUE, fileEncoding = 'UTF-8-BOM')

setwd(caroll_code_folder)

E_Range = c(0) #increased entanglement ,0.01,0.02,0.03

H_sk_range = c("H_sk_1 = H_sk",
               #"H_sk_1 = rep(1/A,A)",
               "H_sk_1 = c(0,0,0,0,rep(1/(A-4), (A-4)))")


vals = c("K = K1", 
         "F = F1",
         "S = S1",
         "P0 = structure(F1,S1,P01)",
         "fer =c(0)",
         
         "K = K2", 
         "F = F2",
         "S = S2",
         "P0 = structure(F2,S2,P02)",
         "fer = c(0)",
         
         
         "K = K3",   ##Added by Vanko et al.
         "F = F3", 
         "S = S3",
         "P0 = structure(F3,S3,P03)",
         "fer =c(0)"
         
) #carying capacity scenarios


grp = c("eval(parse(text=vals[1:5])) ",
        "eval(parse(text=vals[6:10])) ",
        "eval(parse(text=vals[11:15])) ") ##Added by Vanko et al.

#estimating mean ice coverage and number of good and bad years

# mice = mean(ice.data$Extent) 
# fb = length(which(ice.data$Extent<mice))/length(ice.data$Extent)#frequency of bad ice years to date
# 
# x = c(1,50,100)
# yb = c(fb,((1-fb)/2 + fb),1)  
# lmb = lm(yb~x)
# P_ice = rbind(c(fb,0),lmb$coefficients) #two scenarios, first buisness as usual, constant ice coverage, next reduction of ice until there is a frequncy of 1 for bad ice years
# 
# S_G = c(0.779,S_1_Hunt[2:length(S_1_Hunt)]) #survival in a good year (high pup survival)
# S_B = c(0.579,S_1_Hunt[2:length(S_1_Hunt)]) #survival in a good year (low pup survival)
# 
# P0 = structure(F1,S_1_Hunt,P01)

#####################Stochastic (R) simulations begin#######################

H_Range_Q = c(0,2400,3600) #seals hunted, quota

H_Range_Q = H_Range_Q*0.43 #females hunted, quota

rep=100
Fer_Q = rep(NA,106)

for(i in c(1:length(H_Range_Q))){
  
  H = H_Range_Q[i]
  
  for(j in c(1:length(E_Range))){
    
    E = E_Range[j]
    
    for(m in c(1:length(grp))){
      
      eval(parse(text=grp[m])) 
      
      for(k in c(1:length(H_sk_range))){
        
        eval(parse(text=H_sk_range[k]))
        
      
          
          F_R = F 
          
          Settings = c(H,E,k,K,0) #k = 1 is skew, k = 2 is no skew
          
          for(l in c(1:rep)){
            
            P = PopMod(t = 100, 
                                  A = A, 
                                  K = K, 
                                  P = P0, 
                                  B = F_R, 
                                  S = S, 
                                  R = R, 
                                  H = H, 
                                  H_sk = H_sk_1,
                       E = 0, 
                       E_sk = E_sk, 
                                  proportion = FALSE, 
                                  stru = FALSE)
            
            Fer_Q = cbind(Fer_Q,c(Settings,P))
            
            #plot(c(rev(Estimate)/2,P)~c(rev(Year),2021:(2021 + 100)), main = Settings)
            
            print(paste(c(Settings,P[length(P)])))
            
          }
          
        
        
      }
      
    }
    
  }
  
}

write.csv(Fer_Q, "Caroll_simulation.csv") #hunting = fixed quota, no GYBY dynamics
