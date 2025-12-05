#### READ IN #####

setwd(paste(source_data_folder,"Hunting_and_bycatch_samples", sep=""))
df=read.csv("Samples_FI.csv")

#### CLEAN #####

df[is.na(df)]=""
df=df[order(df$Year),]
df=df[df$Year>=2002,]
df=df[is.na(df$Year)==FALSE,]

# How it was found

df$Source[df$Source %in% c("by-catch", "By-catch", "By-caught", "Stranded")]<- "Bycaught"
df$Source[df$Mortality.reason %in% c(4,6)]<- "Bycaught"

df$Source[df$Source %in% c("hunted", "Hunted")]="Hunt"
df$Source[df$Mortality.reason %in% c(1,2,3,5)]<- "Hunt"

#Fix sex
df$Sex[df$Sex %in% c(0, "2?")]=""
df$Sex[df$Sex=="1"] <- "Male"
df$Sex[df$Sex=="2"] <- "Female"

#Fix age
df$Age.class[df$Age.class=="adult"]="Adult"
df$Age.class[df$Age.class=="juvenile"]="Juvenile"



df$Pregnent[(df$Pregnent %in% c(0,1))==FALSE]=""

df$Month[df$Month=="4.-5."]=5
df$Month[df$Month=="kevät"]=5
df$Month[df$Month==" 6-7"]=7
df$Month[df$Month==16]=""
df$Month=as.numeric(df$Month)

df[is.na(df)]=""
df$Age[df$Age..years.==0]=0
df$Age[df$Age..years.==1]=1
df$Age[df$Age..years.==2]=2
df$Age[df$Age..years.==3]=3
df$Age[df$Age..years.==4]=4
df$Age[as.numeric(df$Age..years.)>=5]="5+"

df$Age.Sex=paste( df$Sex,df$Age)


df$CA[(df$CA %in% c(0,1))==FALSE]=""
df$CA[df$Age.Sex!="Female 5+"]<-""
df$Scar[df$Age.Sex!="Female 5+"]<-""

#### HUNTING SAMPLES ####

df_h=df[df$Source=="Hunt",]
df_h=df_h[df_h$Sex%in% c("Male", "Female")& df_h$Age!="",]
y_hs_fi=table(df_h$Age.Sex, df_h$Year)
y_hs_fi=cbind(y_hs_fi,matrix(data=rep(0, 12), ncol=1, nrow=12))
colnames(y_hs_fi)[23]="2024"
y_hs_fi

write.csv(y_hs_fi,  paste(proc_data_folder,"Hunting_samples_FI.csv", sep=""))



#### BYCATCH SAMPLES ####

df_bc=df[df$Source=="Bycaught",]
df_bc=df_bc[df_bc$Sex%in% c("Male", "Female")& df_bc$Age!="",]
y_bc_fi=table(df_bc$Age.Sex, df_bc$Year)
y_bc_fi=cbind(y_bc_fi,matrix(data=rep(0, 66), ncol=6, nrow=11))
y_bc_fi=rbind(y_bc_fi,rep(0,18))
y_bc_fi
colnames(y_bc_fi)[18:23]<- c("2002","2004","2006", "2016", "2019", "2024")
rownames(y_bc_fi)[12]<- "Female 2"
df_bc_fi=as.data.frame(y_bc_fi)
df_bc_fi=df_bc_fi[order(rownames(df_bc_fi)), order(colnames(df_bc_fi))]
df_bc_fi
write.csv(df_bc_fi,  paste(proc_data_folder,"Bycatch_samples_FI.csv", sep=""))


#### REPRODUCTIVE SIGNS ####

df_rs=df[df$Month %in% c(4,5,6),] #adult, femlae, spring
df_rs$Reproduction.sign=1+as.numeric(df_rs$Scar)+2*as.numeric(df_rs$CA)

Z_fi=as.matrix(table(df_rs$Reproduction.sign, df_rs$Year))
Z_fi
write.csv(Z_fi, paste(proc_data_folder,"Reproductive_signs_FI.csv", sep=""))

