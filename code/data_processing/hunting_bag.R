setwd(paste(source_data_folder,"Hunting_bag",sep=""))

df_hb_fi=as.data.frame(t(read.csv("Hunting_bag_FI.csv", check.names = FALSE)))
colnames(df_hb_fi)<-df_hb_fi[1,]
df_hb_fi=df_hb_fi[2:nrow(df_hb_fi),]

df_hb_fi$Total_bag=as.numeric(df_hb_fi$`Total bag - mainland`)+as.numeric(df_hb_fi$`Total bag - Åland`)
df_hb_fi$Total_quota=as.numeric(df_hb_fi$`Quota -mainland`)+as.numeric(df_hb_fi$`Quota -Åland`)
df_hb_fi$Year=seq(from=1998, to=1997+nrow(df_hb_fi), length=nrow(df_hb_fi))
df_hb_fi=df_hb_fi[df_hb_fi$Year>=2002,]

df_hb_fi
df_hb_sw=as.data.frame(t(read.csv("Hunting_bag_SW.csv", check.names = FALSE)))

colnames(df_hb_sw)<-df_hb_sw[1,]
df_hb_sw=df_hb_sw[2:nrow(df_hb_sw),]
df_hb_sw$Year=df_hb_sw[,1]

df_hb=data.frame(Year=seq(from=2002, to=2024, length=23))
df_hb$Hunting_bag_FI=df_hb_fi$Total_bag
df_hb$Quota_FI=df_hb_fi$Total_quota
df_hb$Hunting_bag_SW=as.numeric(df_hb_sw$Sum)
df_hb$Quota_SW=as.numeric(df_hb_sw$Quota)

df_hb
setwd(proc_data_folder)

write.csv(df_hb, "Hunting_bag_quota_FI_SW.csv")
