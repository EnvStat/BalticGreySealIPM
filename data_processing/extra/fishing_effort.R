setwd("/home/milenava/Seals/Grey_seals/Data/Source_data/Fishing_effort")

df_fe=read.csv("Fishing_effort_SW_FI.csv", check.names = FALSE)

df_fe[is.na(df_fe)]<- 0

df_fe=aggregate(df_fe[c("FIN_gear_days_MOD", "SWE_gear_effort")], by=list(df_fe$Year, df_fe$Gear_type), FUN=sum)
colnames(df_fe) <-c("Year", "Gear_type", colnames(df_fe)[3], colnames(df_fe)[4])
df_fe=df_fe[df_fe$Gear_type %in% c("FYK", "GNS"),]
df_fe$Total=df_fe$FIN_gear_days_MOD+df_fe$SWE_gear_effort
df_fe=df_fe[df_fe$Year>=2002,]



setwd("/home/milenava/Seals/Grey_seals/Data/Processed_data")

write.csv(df_fe, "Fishing_effort_SW_FI.csv")
