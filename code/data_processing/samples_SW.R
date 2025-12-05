#### READ IN #####
setwd(paste(source_data_folder,"Hunting_and_bycatch_samples", sep=""))
df=read.csv("Samples_SW.csv")
df_24=read.csv("Samples_SW_2024.csv")
cn=intersect(colnames(df_24), colnames(df))
df=rbind(df[colnames(df) %in% cn],df_24[colnames(df_24) %in% cn])
df$Individual
df$C
#### CLEAN #####

df=df[as.numeric(df$Sam)>=2002,]
df$Group=substring(df$Individual.ID,1,1)
df
##Fix sex
df$Sex[df$Sex=="hane"] <- "Male"
df$Sex[df$Sex=="Hane"] <- "Male"

df$Sex[df$Sex=="Hona"] <- "Female"
df$Sex[df$Sex=="Hona "] <- "Female"
df$Sex[df$Sex=="hona"] <- "Female"

df$Sex[df$Sex=="Okänt"] <- "Unknown"
df$Sex[df$Sex=="okänt"] <- "Unknown"
df$Sex[df$Sex=="?"] <- "Unknown"
df$Sex[df$Sex==""] <- "Unknown"

#Fix how its found
df
df$Source<-"Bycaught"
df$Source[df$How.it.was.found=="Jakt"] <- "Hunt"
df$Source[df$How.it.was.found=="jakt"] <- "Hunt"
df$Source[df$How.it.was.found=="Jakt "] <- "Hunt"
df$Source[df$How.it.was.found=="Jakt?"] <- "Hunt"
df$Source[df$How.it.was.found=="Hunt"] <- "Hunt"


df$Source[df$How.it.was.found==""] <- "Unknown"

df$How="Unknown"
df$How[df$Source=="Hunt"]<- "Hunt"

df$How[df$How.it.was.found=="Strandad"] <- "Stranded"
df$How[df$How.it.was.found=="Fiskeredskap"] <- "Fishing gear"
df$How[df$How.it.was.found=="Fiskeredskap, trål"] <- "Fishing gear"
df$How[df$How.it.was.found=="Avlivad"] <- "Killed"


#Fix age classes
df$Age.class="Unknown"

df$Age.class[df$Age=="0"] <- 0
df$Age.class[df$Age==0] <- 0
df$Age.class[df$Age=="1"] <- 1
df$Age.class[df$Age=="2"] <- 2
df$Age.class[df$Age=="3"] <- 3
df$Age.class[df$Age=="4"] <- 4
df$Age.class[as.numeric(df$Age)>=5] <- "5+"

df$Age.class[df$Age=="16+"] <- "5+"
df$Age.class[df$Age==">25"] <- "5+"



df$Age.Sex=paste( df$Sex,df$Age.class)


df$Age.Sex[df$Status.of.the.uterus %in% c("multiparous", "mulitparous", "dräktig")]
#Make numeric numeric
df$Month=as.numeric(df$Collection.month)

df$Year=as.numeric(df$Sam)

# Fix typos

df$Weight[df$Weight %in% c(721,2530,3020)]=NA
df$Weight[df$Weight==150200]=175
df$Weight[df$Weight==200225]=212
df$Weight[df$Weight==200250]=225

#Fix pregnancy status

df$Pregnancy.status[df$Sex=="Male"] <- "Male"
df$Pregnancy.status[df$Pregnancy.status==""] <- "Unknown"

df$Status.of.the.uterus[df$Status.of.the.uterus %in% c("Juvenil", "Juvenile", "juvenil ", "juvenil","juvenil, hormonpåverkan" )]<- "Juvenile"

df$Length=as.numeric(df$Length)
df$Weight=as.numeric(df$Weight)
df$Blubber.thickness=as.numeric(df$Blubber.thickness)


#### HUNTING SAMPLES ####

df_h=df[df$Source=="Hunt",]
df_h=df_h[df_h$Sex!="Unknown"& df_h$Age.class!="Unknown",]
y_hs_sw=table(df_h$Age.Sex, df_h$Year)

y_hs_sw
write.csv(y_hs_sw,  paste(proc_data_folder,"Hunting_samples_SW.csv", sep=""))

#### BYCATCH SAMPLES ####

df_bc=df[df$Source=="Bycaught",]
df_bc=df_bc[df_bc$Sex!="Unknown"& df_bc$Age.class!="Unknown",]
y_bc_sw=table(df_bc$Age.Sex, df_bc$Year)

write.csv(y_bc_sw,   paste(proc_data_folder,"Bycatch_samples_SW.csv", sep=""))

#### PREGNANCY SAMPLES ####

df_pr=df[df$Pregnancy.status %in% c("Dräktig", "Ej dräktig") & df$Source == "Hunt",]
df_pr=df_pr[df_pr$Month>=8 & df_pr$Status.of.the.uterus!="Juvenile" & df_pr$Age.class %in% c("5+", "Unknown"),]
P=as.matrix(table(df_pr$Pregnancy.status, df_pr$Year))
P[2,]=P[1,]+P[2,]
P=t(P)
colnames(P)=(c( "Pregnant", "Total"))
write.csv(P, paste(proc_data_folder, "Pregnancy_samples_SW.csv", sep=""))
