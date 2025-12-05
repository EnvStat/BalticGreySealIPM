library(reshape2)
library("zoo")

setwd(paste(source_data_folder,"Herring" , sep=""))

# One Drive-> grey_seal -> Data -> Herring

df_BP_GoF=read.csv("Herring_mean_weight_25-29_32.csv", sep="")
df_BP_GoF=df_BP_GoF[df_BP_GoF$Year>=2001,]
df_BP_GoF[,2:ncol(df_BP_GoF)]=1000*df_BP_GoF[,2:ncol(df_BP_GoF)]
df_GoB=read.csv("Herring_mean_weight_GoB.csv",sep=",")
df_GoB=df_GoB[df_GoB$Year>=2001,]


df_catch_BP_GoF=read.csv("Herring_catch_in_numbers_25-29_32.csv")
df_catch_BP_GoF=df_catch_BP_GoF[df_catch_BP_GoF$Year>=2001,]


df_catch_GoB=read.csv("Herring_catch_in_numbers_GoB.csv", check.names = FALSE)
df_catch_GoB=df_catch_GoB[df_catch_GoB$Year>=2001,]

df_BP_GoF$Age5plus=colSums(t(df_BP_GoF[7:10])*colMeans(df_catch_BP_GoF, na.rm=TRUE)[6:9]/sum(colMeans(df_catch_BP_GoF, na.rm=TRUE)[6:9]))
df_GoB$Age5plus=colSums(t(df_GoB[append(paste("X",c(5:14),sep=""), "X15.")])*colMeans(df_catch_GoB)[append(c(5:14), "15+")]/sum(colMeans(df_catch_GoB)[append(c(5:14), "15+")]))


df=df_BP_GoF[c("Year")]
df$BP_GoF=df_BP_GoF$Age5plus
df$GoB=df_GoB$Age5plus

df$BP_GoF.norm=(df$BP_GoF-mean(df$BP_GoF))/sd(df$BP_GoF)
df$GoB.norm=(df$GoB-mean(df$GoB))/sd(df$GoB)

df

write.csv(df,  paste(proc_data_folder,"Herring.csv", sep=""))
