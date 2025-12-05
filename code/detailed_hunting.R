#For Appendix on hunting effort distribution

library("rstan")
library("ggplot2")
library("ggridges")
library("reshape2")
library("bayesplot")
library(zoo)


df_hb=read.csv(paste(data_folder,"Hunting_bag_quota_FI_SW.csv", sep=""))

#### SWEDEN ####


#Read in detailed hunting bags, clean and count the number of shot seals per month in each year
df=read.csv(paste(data_folder,"Detailed_hunting_bag_SW.csv", sep=""))

C=colnames(df)
n=length(C)
for (i in 1:n){
  df_temp=data.frame("Date"=df[,i])
  df_temp$County=C[i]
  
  if (i==1){ df2=df_temp}
  else{df2=rbind(df2, df_temp)}
}
df=df2
df$Date=as.Date(df$Date)
df$Year=format(as.Date(df$Date), "%Y")
df$MD=format(as.Date(df$Date), "%m/%d")


df$Month=as.numeric(format(as.Date(df$Date), "%m"))
df$D=as.numeric(format(as.Date(df$Date), "%d"))

df=na.omit(df)

df=df[df$Year>=2002,]
df=df[df$Year<=2024,]


cov=round(table(df$Year)/df_hb[df_hb$Year%in%df$Year,5]*100)
df_cov=as.data.frame(cov)
colnames(df_cov)=c("Year", "cov")
df=merge(df_cov, df)
df$Month=as.numeric(df$Month)
df$Year=as.character(df$Year)

plot={}
Y=as.numeric(unique(df[df$cov>80,]$Year))-2001 #Years with high enough coverage

cov=c(86,96,100,81,90,100,99,92,99,100,100,100,100,88,96)
unique(df[df$cov>80,]$cov)
Y=Y[is.na(Y)==FALSE]


k=0

for ( i in Y){
  k=k+1
  plot[[k]]=
    
    local({
      i<-i  
      p=ggplot()+
        geom_bar(aes(x=df[df$Year==2001+i,]$Month), fill="blue")+
        #geom_line(aes(x=seq(from=1, to=12), y=c(rep(NA, 3),h_sw.q[2,i,1:9])))+geom_ribbon(aes(x=seq(from=1, to=12), ymin=c(rep(NA, 3),h_sw.q[1,i,1:9]), ymax=c(rep(NA, 3),h_sw.q[3,i,1:9])), alpha=0.1)+
        scale_fill_continuous()+
        labs(x="Month", y="Number of seals shot")+
        scale_x_continuous(breaks=seq(from=0, to=12, length=4)) +
        # ggtitle(paste(i+2001,", ", df[df$Year==2001+i,]$cov[1],"% covergae", sep=""))
        #ggtitle(i+2001)
        ggtitle(paste(i+2001,", ", nrow(df[df$Year==2001+i,])," shot, ", cov[k],"% coverage",sep=""))
      
      print(p)})
}

p=
  
  
  plot[[1]]+ 
  plot[[2]]+
  plot[[3]]+
  plot[[4]]+
  plot[[5]]+
  plot[[6]]+
  plot[[7]]+
  plot[[8]]+
  plot[[9]]+
  plot[[10]]+
  plot[[11]]+ 
  plot[[12]]+
  plot[[13]]+
  plot[[14]]+
  plot[[15]]+
  plot_layout(widths =  rep(3,5))
p

ggsave(filename= paste(figure_folder,"detailed_h_SW_2009-2020.png",sep="" ),plot=p,   width=16, height=10)

#### FINLAND ####


df=read.csv(paste(data_folder,"Detailed_hunting_bag_FI_2001-2014.csv", sep=""))

df$Date=as.Date(df$ti_tilastopaiva)
df$Year=as.numeric(df$ti_tilastovuosi)
df=df[c("Date", "Year")]
df$Month=as.numeric(format(df$Date, "%m"))
df1=df
df=read.csv(paste(data_folder,"Detailed_hunting_bag_FI_2014-2024.csv", sep=""))

df$Date=as.Date(df$point_of_time)
df=df[c("Date")]
df$Year=format(df$Date, "%Y")
df$Month=as.numeric(format(df$Date, "%m"))

df=rbind(df1, df)

df=df[df$Year>=2002,]
cov=round(table(df$Year)/df_hb[df_hb$Year%in%df$Year,3]*100)
df_cov=as.data.frame(cov)
colnames(df_cov)=c("Year", "cov")
df=merge(df_cov, df)
df$Month=as.numeric(df$Month)
df$Year=as.character(df$Year)

i=22
plot={}
Y=as.numeric(unique(df[df$cov>60,]$Year))-2001
Y=Y[is.na(Y)==FALSE]
k=0
for ( i in Y){
  k=k+1
  plot[[k]]=
    
    local({
      i<-i  
      p=ggplot()+
        geom_bar(aes(x=df[df$Year==2001+i,]$Month), fill="blue")+
        #geom_line(aes(x=seq(from=1, to=12), y=c(rep(NA, 3),h_fi.q[2,i,1:9])))+geom_ribbon(aes(x=seq(from=1, to=12), ymin=c(rep(NA, 3),h_fi.q[1,i,1:9]), ymax=c(rep(NA, 3),h_fi.q[3,i,1:9])), alpha=0.1)+
        scale_fill_continuous()+
        labs(x="Month", y="Number of seals shot")+
        scale_x_continuous(breaks=seq(from=0, to=12, length=4)) +
        ggtitle(paste(i+2001,", ",nrow(df[df$Year==2001+i,]), " shot, " ,cov[i-2],"% coverage", sep=""))
      print(p)})
}
p=
  
  
  plot[[1]]+ 
  plot[[2]]+
  plot[[3]]+
  plot[[4]]+
  plot[[5]]+
  plot[[6]]+
  guides(fill=guide_legend(title="% of hunting bag covered"))+
  plot_layout(widths =  rep(3,3))

p
ggsave(filename=  paste(figure_folder,"detailed_h_FI.png",sep=""),plot=p,   width=8, height=5)

