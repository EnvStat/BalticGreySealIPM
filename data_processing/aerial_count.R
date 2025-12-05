setwd(paste(source_data_folder, "Aerial_count", sep=""))

df_count=as.data.frame(t(read.csv("Aerial_count_total.csv", check.names = FALSE)))

colnames(df_count)<-df_count[1,]
df_count=df_count[2:nrow(df_count),]

df_count

df_count$Year=seq(from=2003, to=2002+nrow(df_count), length=nrow(df_count))
df_count$Total_count=df_count$Total

df_count=df_count[c("Year", "Total_count")]
rownames(df_count)<-NULL
df_count
setwd(proc_data_folder)

write.csv(df_count, "Aerial_count.csv")

