
#### Set up ####

t=23
years.sim <- 2025:2080
years.sim.composition <- 2024:2070
years.sim.birth=2025:2060


##### Scenarios ####

scenario.breaks <- c('Historical',"0","1311","2400","3600","4800", "3050")
bias.breaks<- c("Historical","Pup",   "Adult",  "Male", "Female")
her.breaks<-  c("Mean of historical", "+0.2g/yr", "-0.2g/yr", "+4g in 2025", "-4g in 2025")

##### Scenarios for Finnish & Swedish hunting bag ####
Q_fi.scenarios <- matrix(NA,6 , length(years.sim))
Q_fi.scenarios[,1] <- 0 # known Finnish quota for 2023
Q_fi.scenarios[1,] <- 0
Q_fi.scenarios[2,] <- y_hb_fi[22]
Q_fi.scenarios[3,] <- 1200
Q_fi.scenarios[4,] <- 1800
Q_fi.scenarios[5,] <- 2400
Q_fi.scenarios[6,] <- 1550


Q_sw.scenarios <- matrix(NA,6 , length(years.sim))
Q_sw.scenarios[,1] <- 0 
Q_sw.scenarios[1,] <- 0
Q_sw.scenarios[2,] <- y_hb_sw[23]
Q_sw.scenarios[3,] <- 1200
Q_sw.scenarios[4,] <- 1800
Q_sw.scenarios[5,] <- 2400

Q_sw.scenarios[6,] <- 1500

##### Bias scenarios ####
bias_sw=log(colMeans(samples$rho[,1,,3]/apply(samples$rho[,1,,3],1, mean))) #historical bias returned from our model fitting
bias_fi=log(colMeans(samples$rho[,1,,4]/apply(samples$rho[,1,,4],1, mean)))


bias=log(1.5) # 50% more likely to hunt seals in the "biased" groups

bias_fi.scenarios <- matrix(NA,5 , 12)
bias_fi.scenarios[1,]=bias_fi
bias_fi.scenarios[2,]=c(bias,0,0,0,0,0,bias,0,0,0,0,0) #Bias towards pups
bias_fi.scenarios[3,]=c(0,0,0,0,0,bias,0,0,0,0,0,bias) #Bias towards adults
bias_fi.scenarios[4,]=c(0,0,0,0,0,0,bias,bias,bias,bias,bias,bias)  #Bias towards males
bias_fi.scenarios[5,]=c(bias,bias,bias,bias,bias,bias,0,0,0,0,0,0)  #Bias towards females

bias_sw.scenarios <- matrix(NA,5 , 12)
bias_sw.scenarios[1,]=bias_sw
bias_sw.scenarios[2,]=c(bias,0,0,0,0,0,bias,0,0,0,0,0)
bias_sw.scenarios[3,]=c(0,0,0,0,0,bias,0,0,0,0,0,bias)
bias_sw.scenarios[4,]=c(0,0,0,0,0,0,bias,bias,bias,bias,bias,bias)
bias_sw.scenarios[5,]=c(bias,bias,bias,bias,bias,bias,0,0,0,0,0,0)


##### Herring quality scenarios ####

#Note historical mean is 0, since the data has been normalized
m_S=mean(df_h$BP_GoF)
sd_S=sd(df_h$BP_GoF)

h_S.scenarios <- matrix(NA, 5, length(years.sim))
h_S.scenarios[,1] <- 0 
h_S.scenarios[1,] <-0
h_S.scenarios[2,] <- h_S.scenarios[2,1] + c(0:(length(years.sim)-1))*0.2/sd_S
h_S.scenarios[3,] <- h_S.scenarios[3,1] - c(0:(length(years.sim)-1))*0.2/sd_S
h_S.scenarios[4,2:ncol(h_S.scenarios)] <- h_S.scenarios[4,1] + 4/sd_S
h_S.scenarios[5,2:ncol(h_S.scenarios)] <- h_S.scenarios[5,1] - 4/sd_S



m_G=mean(df_h$GoB)
sd_G=sd(df_h$GoB)

h_G.scenarios <- matrix(NA, 5, length(years.sim))
h_G.scenarios[,1] <- 0 
h_G.scenarios[1,] <-0
h_G.scenarios[2,] <- h_G.scenarios[2,1] + c(0:(length(years.sim)-1))*0.2/sd_G
h_G.scenarios[3,] <- h_G.scenarios[3,1] - c(0:(length(years.sim)-1))*0.2/sd_G
h_G.scenarios[4,2:ncol(h_G.scenarios)] <- h_G.scenarios[4,1] + 4/sd_G
h_G.scenarios[5,2:ncol(h_G.scenarios)] <- h_G.scenarios[5,1] - 4/sd_G


##### Vertical denisty plots #### 
#These are the "interesting" scenarios plotted separately
Q.list=c(1,2,3,4,4,4,4,4,5)
B.list=c(1,1,1,1,2,3,4,5,1)

H.list=c(1,2,3,4,5) #for herring quality scenarios, we only plotted under 3600 harvest size 

