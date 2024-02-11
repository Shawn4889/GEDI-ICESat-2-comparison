#Author: Xiaoxuan Li
library(corrplot)
library(tidyverse)
library(psych)
library(xlsx)
library(car)
library(MASS)
library(rsq)
library(ggpmisc)
library(ggrepel)
library(dplyr)
library(rlist)
library(data.table)
library(lemon)
library(lmtest)
library(reshape2)
library(ggpubr)
library(ggplot2)
library(mltools)
library(tidyr)
library(stats)
library(lidR)
library(tools) 
library(raster)
library(ggpointdensity)
library(viridis)
library(grid)
library(readxl)
library(ehaGoF)
library(Metrics)
library(rGEDI)
windowsFonts(A = windowsFont("Times New Roman"))


dir_GEDI <- "E:\\GEDI\\Result\\Result\\All_02272023_2.csv"
Data = read.csv(dir_GEDI,header=T)
Data <- na.omit(Data)
Data <- Data[Data$status == "Leaf-on",]
Data <- Data[Data$ql2a ==1,]
Data <- Data[Data$ql2b ==1,]
Data <- na.omit(Data)

Data1 <- Data[Data$RH1_98 < 2.34,]

r2 <- round(cor(Data$RH1_98, Data$p98,method = "pearson")^2,3)
RMSE<- sqrt(mean((Data$RH1_98 - Data$p98)^2))
rRMSE<- 100*sqrt(mean((Data$RH1_98 - Data$p98)^2))/mean(Data$p98)
MD <- bias(Data$p98, Data$RH1_98)
RB <- 100*bias(Data$p98, Data$RH1_98)/mean(Data$p98)

p1 <- ggplot(Data, aes(x=p98, y=RH1_98))+ 
  geom_pointdensity()+
  scale_color_viridis()+
  theme_bw()+
  annotate("text",x=1,y=13,hjust = 0,size = 10,family= "A",
           label= paste(expression(" "~R^2),": ",r2 
           ),parse=TRUE) + 
  annotate("text",x=1,y=10.3,hjust = 0,size = 10,family= "A",
           label= paste(
             "  RMSE: ", round(RMSE,3),"m",             
             "\n" , " %RMSE: ", round(rRMSE,3),"%",
             "\n" , " Bias: ", round(MD,3),"m",
             "\n" , " %Bias: ", round(RB,3),"%",
             "\n" , " Sample size: ", nrow(Data)
           )) +  
  coord_cartesian(xlim = c(0,15),ylim = c(0,15))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(text=element_text(size=40)) + 
  scale_x_continuous(minor_breaks = seq(0,15,2),breaks = seq(0,15,2))+
  scale_y_continuous(minor_breaks = seq(0,15,2),breaks = seq(0,15,2))+
  theme(text=element_text(family="A"))+
  theme(legend.title = element_blank())+
  geom_abline(intercept = 0, slope = 1,color="black", linetype="solid", size=1.5)+
  labs(x=expression("ALS P98 (m)"), 
       y=expression("GEDI RH98 (m)"))+
  theme(legend.position = "none")




dir_ICESAT <- "E:\\ICESat2\\Results\\Result\\All_ICESat2_P.csv"
Data = read.csv(dir_ICESAT,header=T)
Data <- na.omit(Data)
Data <- Data[Data$p98 >0,]

Data <- Data[Data$h_canopy_2 >0,]

Data <- Data[Data$night_flag ==1,]
Data <- Data[Data$status =="Leaf-on",]

r2 <- round(cor(Data$h_canopy_2, Data$p98,method = "pearson")^2,3)
RMSE <- sqrt(mean((Data$h_canopy_2 - Data$p98)^2))
rRMSE <- 100*sqrt(mean((Data$h_canopy_2 - Data$p98)^2))/mean(Data$p98)
MD <- bias(Data$p98, Data$h_canopy_2)
RB <- 100*bias(Data$p98, Data$h_canopy_2)/mean(Data$p98)

p2 <- ggplot(Data, aes(x=p98, y=h_canopy_2))+ 
  geom_pointdensity()+
  scale_color_viridis()+
  theme_bw()+
  coord_cartesian(xlim = c(0,15),ylim = c(0,15))+
  annotate("text",x=1,y=13,hjust = 0,size = 10,family= "A",
           label= paste(expression(" "~R^2),": ",r2 
           ),parse=TRUE) + 
  annotate("text",x=1,y=10.3,hjust = 0,size = 10,family= "A",
           label= paste(
             "  RMSE: ", round(RMSE,3),"m",             
             "\n" , " %RMSE: ", round(rRMSE,3),"%",
             "\n" , " Bias: ", round(MD,3),"m",
             "\n" , " %Bias: ", round(RB,3),"%",
             "\n" , " Sample size: ", nrow(Data)
           )) +  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(text=element_text(size=40)) + 
  scale_x_continuous(minor_breaks = seq(0,15,2),breaks = seq(0,15,2))+
  scale_y_continuous(minor_breaks = seq(0,15,2),breaks = seq(0,15,2))+
  theme(text=element_text(family="A"))+
  theme(legend.title = element_blank())+
  geom_abline(intercept = 0, slope = 1,color="black", linetype="solid", size=1.5)+
  labs(x=expression("ALS P98 (m)"), 
       y=expression("ICESat-2 canopy height (m)"))+
  theme(legend.position = "none")


ggarrange(p1,p2)
out = "E:\\ICESat2\\Results\\Figure\\composite.jpg"
ggsave(out,height=12, width=24, dpi=600)




#canopy height gedi biasplot
dir_GEDI <- "E:\\GEDI\\Result\\Result\\All_02272023_2.csv"
Data = read.csv(dir_GEDI,header=T)
Data <- na.omit(Data)
Data <- Data[Data$status == "Leaf-on",]
Data <- Data[Data$ql2a ==1,]
Data <- Data[Data$ql2b ==1,]
Data <- na.omit(Data)

Data$diff <- Data$RH1_98 - Data$p98
a1 <- tapply(Data$diff, cut(Data$p98,breaks = seq(0, 10, by=1),dig.lab=1), mean)

x <- seq(0,9,1)
l <- cbind(x,a1)
df <- data.frame(l)
data.frame(l)
#df[[2]][[1]] <- 0
#df[[2]][[2]] <- 0

p3<- ggplot(df, aes(x=x+0.5, y=a1)) + 
  geom_bar(stat='identity')+
  theme_bw()+
  coord_cartesian(ylim = c(-2.5, 2.5))+
  annotate("text",x=8,y=2,hjust = 0,size = 15, label= "(a)") + 
  scale_y_continuous(minor_breaks = seq(-2.5, 2.5, 0.5),breaks =  seq(-2.5, 2.5, 0.5))+
  scale_x_continuous(minor_breaks = round(seq(0,10,1),digits = 1),
                     breaks = round(seq(0,10,1),digits = 1))+
  labs(x="ALS P98 (m)", 
       y="Difference between GEDI and ALS canopy height (m)")+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=25))+
  theme(plot.title = element_text(hjust = 0.5))



#canopy height icesat-2 biasplot
dir_ICESAT <- "E:\\ICESat2\\Results\\Result\\All_ICESat2_P.csv"
Data = read.csv(dir_ICESAT,header=T)
Data <- na.omit(Data)
Data <- Data[Data$p98 >0,]
Data <- Data[Data$h_canopy_2 >0,]
Data <- Data[Data$night_flag ==1,]
Data <- Data[Data$status =="Leaf-on",]

Data$diff <- Data$h_canopy_2 - Data$p98
a1 <- tapply(Data$diff, cut(Data$p98,breaks = seq(0, 10, by=1),dig.lab=1), mean)

x <- seq(0,9,1)
l <- cbind(x,a1)
df <- data.frame(l)
data.frame(l)

p4<- ggplot(df, aes(x=x+0.5, y=a1)) + 
  geom_bar(stat='identity')+
  theme_bw()+
  coord_cartesian(ylim = c(-2.5, 2.5))+
  annotate("text",x=8,y=2,hjust = 0,size = 15, label= "(b)") + 
  scale_y_continuous(minor_breaks = seq(-2.5, 2.5, 0.5),breaks =  seq(-2.5, 2.5, 0.5))+
  scale_x_continuous(minor_breaks = round(seq(0,10,1),digits = 1),
                     breaks = round(seq(0,10,1),digits = 1))+
  labs(x="ALS P98 (m)", 
       y="Difference between ICESat-2 and ALS canopy height (m)")+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=25))+
  theme(plot.title = element_text(hjust = 0.5))



ggarrange(p3,p4)
out = "E:\\ICESat2\\Results\\Figure\\composite_biasplot_als_p98.jpg"
ggsave(out,height=12, width=24, dpi=600)




#ground - GEDI
dir_GEDI <- "E:\\GEDI\\Result_3\\Result\\Result_gedi_ground.csv"
Data = read.csv(dir_GEDI,header=T)
Data <- na.omit(Data)

r2 <- round(cor(Data$GEDI, Data$ALS,method = "pearson")^2,3)
RMSE<- sqrt(mean((Data$GEDI - Data$ALS)^2))
MD <- bias(Data$ALS, Data$GEDI)


#ground - ICESat-2
dir_ICESAT <- "E:\\ICESat2\\Results\\Result\\All_ICESat2_P.csv"
Data = read.csv(dir_ICESAT,header=T)
Data <- na.omit(Data)
Data <- Data[Data$p98 >0,]
Data <- Data[Data$h_canopy_2 >0,]

Data <- Data[Data$night_flag ==1,]
Data <- Data[Data$status =="Leaf-on",]

r2 <- round(cor(Data$h_te_best_, Data$dem_h,method = "pearson")^2,3)
RMSE <- sqrt(mean((Data$h_te_best_ - Data$dem_h)^2))
MD <- bias(Data$dem_h, Data$h_te_best_)


#GEDI vs. ICESat-2
filedir <- "E:\\GEDI\\Result_3\\Result\\Result_icesat_gedi_r.xlsx"
Data = read_excel(filedir, sheet = "gedi_ice")

r2 <- round(cor(Data$ICESat_RH98, Data$GEDI_RH98,method = "pearson")^2,3)
RMSE <- sqrt(mean((Data$ICESat_RH98 - Data$GEDI_RH98)^2))
rRMSE <- 100*sqrt(mean((Data$ICESat_RH98 - Data$GEDI_RH98)^2))/mean(Data$GEDI_RH98)
MD <- bias(Data$GEDI_RH98, Data$ICESat_RH98)
RB <- 100*bias(Data$GEDI_RH98, Data$ICESat_RH98)/mean(Data$GEDI_RH98)

ggplot(Data, aes(x=GEDI_RH98, y=ICESat_RH98))+ 
  geom_pointdensity()+
  scale_color_viridis()+
  theme_bw()+
  annotate("text",x=9,y=7,hjust = 0,size = 10,family= "A",
           label= paste(expression(" "~R^2),": ",r2 
           ),parse=TRUE) + 
  annotate("text",x=9,y=4.3,hjust = 0,size = 10,family= "A",
           label= paste(
             "  RMSE: ", round(RMSE,3),"m",             
             "\n" , " %RMSE: ", round(rRMSE,3),"%",
             "\n" , " MSD: ", round(MD,3),"m",
             "\n" , " RMSD: ", round(RB,3),"%",
             "\n" , " Sample size: ", nrow(Data)
           )) +  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(text=element_text(size=40)) + 
  scale_x_continuous(limits = c(0,15),minor_breaks = seq(0,15,2))+
  scale_y_continuous(limits = c(0,15),minor_breaks = seq(0,15,2))+
  theme(text=element_text(family="A"))+
  theme(legend.title = element_blank())+
  geom_abline(intercept = 0, slope = 1,color="black", linetype="solid", size=1.5)+
  labs(x=expression("GEDI RH98 (m)"), 
       y=expression("ICESat canopy height (m)"))+
  theme(legend.position = "none")

out = "E:\\ICESat2\\Results\\Figure\\gedi_icesat_2.jpg"
ggsave(out,height=12, width=12, dpi=600)
