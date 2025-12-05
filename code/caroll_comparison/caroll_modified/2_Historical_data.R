#Daire Carroll University of Gothenburg 2022
#Grey seal PVA script 2, plot historical reconstruction of grey seal populaiton sizes in the 1900s

setwd(caroll_data_folder)
my.data = read.csv("Historical.csv", sep = ",", header = TRUE, fileEncoding = 'UTF-8-BOM')
attach(my.data)
head(my.data)

ggplot(my.data  , aes(Year)) +	
  labs(x="Year", y=expression("Estimated population size")) +
  cleanup_grid + cleanup_text   +
  
  geom_line(aes(y = (Aerial_Surveys_Estimate*0.7)/0.6 ), linetype = "dashed", size = 1, color = "black", alpha = 1) +
  geom_line(aes(y = Aerial_Surveys_Estimate ), size = 1.3, color = "black", alpha = 1) +
  geom_line(aes(y = (Aerial_Surveys_Estimate*0.7)/0.8 ), linetype = "dashed", size = 1, color = "black", alpha = 1) +
  geom_line(aes(y = Loss_10), size = 1, color = "#CC6677", linetype = "dashed", alpha = 1)+
  geom_line(aes(y = Loss_20), size = 1.3, color = "#CC6677", alpha = 1) +
  geom_line(aes(y = Loss_30), size = 1, color = "#CC6677", linetype = "dashed", alpha = 1) +
  
  theme(plot.margin = margin(10, 10, 10, 10)) +
  
  xlim(1900,2022)  

max(na.omit(Loss_10))
max(na.omit(Loss_30))
