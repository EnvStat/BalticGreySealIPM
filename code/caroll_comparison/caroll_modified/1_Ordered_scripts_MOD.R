#Daire Carroll University of Gothenburg 2022
#Grey seal PVA script 1, load packages and base scripts for core parameters and functions

if(!require("ggplot2")){
  install.packages("ggplot2",dependencies = TRUE)
}
if(!require("dplyr")){
  install.packages("dplyr",dependencies = TRUE)
} #load in packages

cleanup_grid = theme(panel.border = element_blank(),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     panel.background = element_blank(),
                     axis.line = element_line(color = "black"),)

cleanup_text = theme(axis.text.x = element_text(color = "black", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
                     axis.text.y = element_text(color = "black", size = 12, angle = 0, hjust = 0.8, vjust = 0.4, face = "plain"),  
                     axis.title.x = element_text(color = "black", size = 20, angle = 0, hjust = .5, vjust = 0, face = "plain"),
                     axis.title.y = element_text(color = "black", size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain")) #clean up for ggplots			


setwd(caroll_code_folder) #replace with code location
source("Basic_model_MOD.R", echo = TRUE) 
source("Generation_time_calculator.R", echo = TRUE) 
source("Parameters_MOD.R", echo = TRUE) 

structure(F1,S_1_Hunt,1)
