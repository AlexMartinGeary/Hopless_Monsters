#FIGURE 1
rm(list=ls())
library(ggfortify)
library(plyr)
library(lattice)
setwd("/Users/mqbpqam6/Desktop/PhD/Last_Hurrah_Data_Gen/")
group.colours <- c("Unknown" = "deepskyblue4", "Dominant" = "black", "Recessive" = "mediumvioletred", "Both" = "darkseagreen", "None"= "grey92")
split <- read.csv("Last_Hurrah_Dataset.csv", header = T)

split<-subset(split, Ortho_Age!="0")
split <- na.omit(split) #Remove any rows with null assignments from the analysis
split$Disease_Association <- relevel(split$Disease_Association, "Recessive", "Both", "Dominant", "Unknown", "None")
ADS =count(split,c("Ortho_Age", "Disease_Association"))
ADS$Ortho_Age <- factor(ADS$Ortho_Age)
ADS$Disease_Association <- factor(ADS$Disease_Association)
ADS$Disease_Association <- relevel(ADS$Disease_Association, "None")
ADS$Ortho_Age <- factor(ADS$Ortho_Age, levels=rev(levels(ADS$Ortho_Age)))

Alldisease<-subset(split, Disease_Association=="Dominant" | Disease_Association=="Recessive" | Disease_Association=="Both" | Disease_Association=="Unknown")
domdis<-subset(split, Disease_Association=="Dominant")
recdis<-subset(split, Disease_Association=="Recessive")
rddis<-subset(split, Disease_Association=="Both")
undis<-subset(split, Disease_Association=="Unknown")
ADSdom =count(domdis,c("Ortho_Age"))
ADSrec =count(recdis,c("Ortho_Age"))
ADSrd =count(rddis,c("Ortho_Age"))
ADSun =count(undis,c("Ortho_Age"))
ADSage =count(split,c("Ortho_Age"))
DiseaseData<-subset(split, Disease_Association=="Dominant" | Disease_Association=="Recessive") 
BDS = count(na.omit(DiseaseData), c("Ortho_Age", "Disease_Association"))
cumuSET <- BDS[order(BDS$Ortho_Age),] 
BADSinv=cumuSET[order(nrow(cumuSET):1),] #invert row order
BADS1 <- subset(BADSinv, Disease_Association=="Dominant")
BADS2 <- subset(BADSinv, Disease_Association=="Recessive")
BADS3 <- subset(BADSinv, Disease_Association=="Unknown")
xout1 <- transform(BADS1, cumFreq = cumsum(freq), relative = prop.table(freq))
xout2 <- transform(BADS2, cumFreq = cumsum(freq), relative = prop.table(freq))
xout3 <- transform(BADS3, cumFreq = cumsum(freq), relative = prop.table(freq))
data<-join(xout1, xout2, type = "full", match = "all")
data<-join(data, xout3, type = "full", match = "all")
library(lattice)
tg <- ddply(data, c("Ortho_Age", "Disease_Association","cumFreq"), summarise, length=cumFreq)
library(raster)
library(gridExtra)
library(grid)
ADS$Disease_Association <- relevel(ADS$Disease_Association, "Unknown")
ADS$Disease_Association <- relevel(ADS$Disease_Association, "None")
F1 <- ggplot(data=ADS, aes(x=Ortho_Age, y=freq, fill=Disease_Association)) +
  geom_bar(stat="identity") + 
  ggtitle("a.") +
  ylab(label="Number of genes") +
  xlab(label="Gene Age MYA(MRCA)") +
  theme(plot.title = element_text(lineheight=.8, face="bold"))+ theme(panel.background = element_rect(fill = "white")) + theme(axis.text.x=element_text(angle=30,hjust=1,vjust=1),axis.text=element_text(size=10))+ 
  theme(axis.title.y=element_text(face = "bold"))+ theme(axis.title.x=element_text(face = "bold"))+
  scale_fill_manual(values=group.colours)+ 
  theme(legend.background = element_rect(fill="white", size=0.5, linetype="solid", colour ="grey86"))+
  theme(plot.title = element_text(size = 12),panel.grid.major = element_line(colour="white", size = (0.2)), panel.grid.minor = element_line(size = (0.2), colour="white"), panel.background = element_rect(fill = 'white', colour = 'transparent'), plot.margin=unit(c(0.5,0.5,0,0.2), "cm"))+
  theme(legend.position=c(0.77, 0.34))
F1 <- F1 + guides(fill=guide_legend(title="Disease Association"))                                                                
tg$Ortho_Age <- factor(tg$Ortho_Age)
tg$Ortho_Age <- factor(tg$Ortho_Age, levels=rev(levels(tg$Ortho_Age)))
F2 <- ggplot(tg, aes(x=factor(Ortho_Age), y=cumFreq, colour=Disease_Association, group=Disease_Association))+ 
  geom_point() + 
  geom_line() + 
  ggtitle("b.") +
  xlab("Gene Age MYA(Ortho_Age)") +
  ylab("Cumulative frequency") +
  scale_y_continuous(limits=c(0,2000), breaks=seq(0,2000, by = 200), expand = c(0, 0))+
  scale_x_discrete()+
  theme(plot.title = element_text(size=18,lineheight=0.8, face="bold"))+ 
  theme(panel.background = element_rect(fill = "white")) +
  theme(axis.text.x=element_text(angle=40,hjust=1,vjust=1),axis.text=element_text(size=4))+ 
  theme(axis.text.x=element_text(size=5))+ 
  theme(axis.title.y=element_text(size=7,face="bold"))+ theme(axis.title.x=element_text(size=7,face="bold"))+
  scale_colour_manual("Disease_Association",values=group.colours,breaks=c("Dominant", "Recessive","All"), labels=c( "Dominant", "Recessive", "All disease"))+ 
  theme(plot.title = element_text(size = 8),panel.grid.major = element_line(colour="snow2", size = (0.2)), panel.grid.minor = element_line(size = (0.2), colour="snow1"), panel.background = element_rect(fill="aliceblue",colour=NA), legend.position = "none")
grid.newpage()
v1<-viewport(width = 1, height = 1, x = 0.5, y = 0.5) #plot area for the main map
v2<-viewport(width = 0.4, height = 0.4, x = 0.77, y = 0.72) #plot area for the inset map
print(F1,vp=v1)  
print(F2,vp=v2)



#FIGURE 2
rm(list=ls())
library(ggfortify)
library(plyr)
library(lattice)
setwd("/Users/mqbpqam6/Desktop/PhD/Last_Hurrah_Data_Gen/")
split <- read.csv("Last_Hurrah_Dataset.csv", header = T)[c(3,5,7)]		#Insert path to DS1
split <- na.omit(split) #Remove any rows with null assignments from the analysis
split<-subset(split, Ortho_Age!="0")
split$Ortho_Age<-as.factor(split$Ortho_Age)
levels(split$Ortho_Age)[levels(split$Ortho_Age)=="73"]<-"Recent"
levels(split$Ortho_Age)[levels(split$Ortho_Age)=="67"]<-"Recent"
levels(split$Ortho_Age)[levels(split$Ortho_Age)=="43"]<-"Recent"
levels(split$Ortho_Age)[levels(split$Ortho_Age)=="29"]<-"Recent"
levels(split$Ortho_Age)[levels(split$Ortho_Age)=="20"]<-"Recent"
levels(split$Ortho_Age)[levels(split$Ortho_Age)=="15"]<-"Recent"
levels(split$Ortho_Age)[levels(split$Ortho_Age)=="9"]<-"Recent"
levels(split$Ortho_Age)[levels(split$Ortho_Age)=="0"]<-"Recent"

Alldisease<-subset(split, Disease_Association=="Dominant" | Disease_Association=="Recessive" | Disease_Association=="Both" | Disease_Association=="Unknown")
domdis<-subset(split, Disease_Association=="Dominant")
recdis<-subset(split, Disease_Association=="Recessive")
rddis<-subset(split, Disease_Association=="Both")
undis<-subset(split, Disease_Association=="Unknown")
ADSdisease =count(Alldisease,c("Ortho_Age"))
ADSdom =count(domdis,c("Ortho_Age"))
ADSrec =count(recdis,c("Ortho_Age"))
ADSrd =count(rddis,c("Ortho_Age"))
ADSun =count(undis,c("Ortho_Age"))
ADSage =count(split,c("Ortho_Age"))
finalSET <-merge(x = ADSage, y = ADSdisease, by = "Ortho_Age", all.x = TRUE)
colnames(finalSET)[2] <- "Total.genes"
colnames(finalSET)[3] <- "All.Disease"
finalSET <- transform(finalSET, Disease = All.Disease / Total.genes)
finalSET <-merge(x = finalSET, y = ADSdom, by = "Ortho_Age", all.x = TRUE)
colnames(finalSET)[5] <- "Dominant.genes"
finalSET <- transform(finalSET, Dominant = Dominant.genes / Total.genes)
finalSET <-merge(x = finalSET, y = ADSrec, by = "Ortho_Age", all.x = TRUE)
colnames(finalSET)[7] <- "Recessive.genes"
finalSET <- transform(finalSET, Recessive = Recessive.genes / Total.genes)
finalSET <-merge(x = finalSET, y = ADSrd, by = "Ortho_Age", all.x = TRUE)
colnames(finalSET)[9] <- "RecDom.genes"
finalSET <- transform(finalSET, Both = RecDom.genes / Total.genes)
finalSET <-merge(x = finalSET, y = ADSun, by = "Ortho_Age", all.x = TRUE)
colnames(finalSET)[11] <- "Unknown.genes"
finalSET <- transform(finalSET, Unknown = Unknown.genes / Total.genes)
lines <-(finalSET)[c(1,4,6,8,10,12)]
lines[is.na(lines)] <- 0
lines$Ortho_Age <- factor(lines$Ortho_Age)
lines$Ortho_Age <- factor(lines$Ortho_Age, levels=rev(levels(lines$Ortho_Age)))
library("reshape2")
lines1 <- melt(lines, id.vars="Ortho_Age", value.name="Proportion", variable.name="Disease")
F3<-ggplot(data=lines1, aes(x=Ortho_Age, y=Proportion, group = Disease, colour = Disease)) +
  geom_line() +
  ggtitle(" ") +
  ylab(label="Proportion of total genes per Age") +
  xlab(label=" ") +
  theme(axis.text.x = element_blank(), axis.ticks = element_blank())+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 32.2, b = 0, l = 0)))+
  theme(axis.title.x = element_text(margin = margin(t = 0, r = 0, b = 30.5, l = 1)))+
  theme(plot.title = element_text(lineheight=.8, face="bold"))+ theme(panel.background = element_rect(fill = "white"))+
  scale_y_continuous(limits=c(0,1), breaks=seq(0,1, by = 0.2), expand = c(0.05, 0), position = "right")+
  scale_colour_manual("Disease",values = c('red','black','mediumvioletred','darkseagreen','deepskyblue4'),breaks=c("Dominant", "Recessive", "Both", "Disease", "Unknown"), labels=c( "Dominant", "Recessive", "Both", "All Disease", "Unknown"))+
  theme(panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent"))+
  theme(legend.background = element_rect(fill="white", size=0.5, linetype="solid", colour ="grey86"))+
  theme(panel.background = element_rect(fill = 'transparent', colour = 'transparent'), legend.key.size = unit(.9, "cm"), legend.position = c(.90, .8),legend.title=element_text(size=9, face="bold"), legend.text = element_text(size = 9), plot.margin=unit(c(0.5,0.5,0,0.2), "cm"))
rm(list= ls()[!(ls() %in% c('F3'))])
library(plyr)
library(ggfortify)
library(gridExtra)
library(grid)
setwd("/Users/mqbpqam6/Desktop/PhD/Last_Hurrah_Data_Gen/")
split <- read.csv("Last_Hurrah_Old_Genes.csv", stringsAsFactors=FALSE, header = T)[c(6,8)] #Path to DS1
split <- na.omit(split)
split<-subset(split, Ortho_Age!="0")
split$Ortho_Age<-as.factor(split$Ortho_Age)
levels(split$Ortho_Age)[levels(split$Ortho_Age)=="73"]<-"Recent"
levels(split$Ortho_Age)[levels(split$Ortho_Age)=="67"]<-"Recent"
levels(split$Ortho_Age)[levels(split$Ortho_Age)=="43"]<-"Recent"
levels(split$Ortho_Age)[levels(split$Ortho_Age)=="29"]<-"Recent"
levels(split$Ortho_Age)[levels(split$Ortho_Age)=="20"]<-"Recent"
levels(split$Ortho_Age)[levels(split$Ortho_Age)=="15"]<-"Recent"
levels(split$Ortho_Age)[levels(split$Ortho_Age)=="9"]<-"Recent"
levels(split$Ortho_Age)[levels(split$Ortho_Age)=="0"]<-"Recent"

split$Haplosufficiency_Rank	<-as.numeric(split$Haplosufficiency_Rank)
# Box Plots
x1 <- split$Haplosufficiency_Rank[split$Ortho_Age=="796"]
x2 <- split$Haplosufficiency_Rank[split$Ortho_Age=="615"]
x3 <- split$Haplosufficiency_Rank[split$Ortho_Age=="435"]
x4 <- split$Haplosufficiency_Rank[split$Ortho_Age=="413"]
x5 <- split$Haplosufficiency_Rank[split$Ortho_Age=="351"]
x6 <- split$Haplosufficiency_Rank[split$Ortho_Age=="311"]
x7 <- split$Haplosufficiency_Rank[split$Ortho_Age=="176"]
x8 <- split$Haplosufficiency_Rank[split$Ortho_Age=="158"]
x9 <- split$Haplosufficiency_Rank[split$Ortho_Age=="105"]
x10 <- split$Haplosufficiency_Rank[split$Ortho_Age=="96"]
x18 <- split$Haplosufficiency_Rank[split$Ortho_Age=="Recent"]

par(cex.lab=1, cex.axis=.8) 
boxplot(x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x18, names=c("796", "615", "435", "413", "351", "311","176","158","105","96","Recent"), 
        col="aliceblue", notch=TRUE, outline=FALSE, ylim=c(100,0),frame=F) +
  theme(panel.background = element_rect(fill = 'transparent', colour = 'transparent'))
title(xlab="Gene Age MYA", ylab="Haploinsufficiency (0=HI, 100=HS)")
abline(h = 10, col = "red", lty=6, lwd=2) 
#text(11,8, "'True' Haploinsufficiency", pos = 4, col = "red", cex = 1)
v3<-viewport(width = .9, height = .89, x = 0.535, y = 0.5) #plot area for the overlay map
print(F3,vp=v3) 






#FIGURE 3
rm(list=ls())
library(ggfortify)
library(plyr)
library(tidyr)
library(gridExtra)
library(grid) 
setwd("/Users/mqbpqam6/Desktop/PhD/Last_Hurrah_Data_Gen/")
disess <- read.csv("Last_Hurrah_Dataset.csv") [c(4,5,7)]	#Path to DS1
Hap = na.omit(disess)


Hap<-subset(Hap, Disease_Association=="Recessive" | Disease_Association=="Dominant" | Disease_Association=="Unknown")
Hap <- unite(Hap, Paralog_and_Disease_Status, c(Disease_Association, Paralog_status), remove=FALSE)
Hap$decile <- Hap$Haplosufficiency_Rank


n = 0
for (i in Hap$decile){
  n = n+1
  if(i <= 10){
    Hap$decile[n]<-"10th"  
  }
  if ((i > 10) && (i <=20)){
    Hap$decile[n]<-"9th"
  }  
  if ((i > 20) && (i <=30)){
    Hap$decile[n]<-"8th"
  }  
  if ((i > 30) && (i <=40)){
    Hap$decile[n]<-"7th"
  }  
  if ((i > 40) && (i <=50)){
    Hap$decile[n]<-"6th"
  }  
  if ((i > 50) && (i <=60)){
    Hap$decile[n]<-"5th"
  }  
  if ((i > 60) && (i <=70)){
    Hap$decile[n]<-"4th"
  }  
  if ((i > 70) && (i <=80)){
    Hap$decile[n]<-"3rd"
  }  
  if ((i > 80) && (i <=90)){
    Hap$decile[n]<-"2nd"
  }  
  if ((i > 90) && (i <=100)){
    Hap$decile[n]<-"1st"
  }  
}


HaploADS = count(na.omit(Hap), c("Paralog_and_Disease_Status"))
New3 = count(na.omit(Hap), c("Paralog_and_Disease_Status","decile"))
New4 <-merge(x = HaploADS, y = New3, by = "Paralog_and_Disease_Status", all.x = TRUE)
ADS3 <- transform(New4, Proportion.of.genes = freq.y / freq.x)
ADS3$Paralog_and_Disease_Status <-factor(ADS3$Paralog_and_Disease_Status)

group.colours <- c("Unknown_Ohnolog" = "gray", "Unknown_Singleton" = "darkseagreen2", "Unknown_SSD" = "plum1", "Dominant_Singleton" = "darkgreen", "Dominant_Ohnolog" = "Black", "Dominant_SSD" = "mediumvioletred", "Recessive_Ohnolog" = "gray40", "Recessive_Singleton" = "mediumseagreen", "Recessive_UPP" = "darkcyan", "Recessive_SSD" = "violetred1")
ONE<- ggplot(data=ADS3, aes(x=decile, y=Proportion.of.genes, colour=Paralog_and_Disease_Status,linetype = Paralog_and_Disease_Status))+ 
  ggtitle("a") +
  xlab("Haplosufficiency Decile") +
  ylab("Proportion of genes") +
  geom_point(aes(shape=factor(Paralog_and_Disease_Status)), size = 1) +
  geom_smooth(aes(group=Paralog_and_Disease_Status), method= loess, se=FALSE, size=1)+
  scale_colour_manual("Paralog_and_Disease_Status", values=group.colours, breaks=c("Dominant_Singleton", "Dominant_Ohnolog", "Dominant_SSD", "Recessive_Singleton", "Recessive_Ohnolog", "Recessive_SSD","Unknown_Singleton","Unknown_Ohnolog", "Unknown_SSD"), labels=c("Dominant Singleton","Dominant Ohnolog", "Dominant SSD", "Recessive Singleton","Recessive Ohnolog", "Recessive SSD","Unknown Singleton","Unknown Ohnolog", "Unknown SSD"))+
  scale_linetype_manual(name="Paralog_and_Disease_Status", values= c("Dominant_Singleton"=1, "Dominant_Ohnolog"=1, "Dominant_UPP"=1, "Dominant_SSD"=1,  "Recessive_Singleton"=8, "Recessive_Ohnolog"=8, "Recessive_UPP"=8, "Recessive_SSD"=8, "Unknown_Singleton"=3, "Unknown_Ohnolog"=3, "Unknown_UPP"=3, "Unknown_SSD"=3), labels=c("Dominant Singleton", "Dominant Ohnolog", "Dominant SSD", "Recessive Singleton", "Recessive Ohnolog", "Recessive SSD", "Unknown Singleton", "Unknown Ohnolog", "Unknown SSD"))+
  scale_shape_manual(name="Paralog_and_Disease_Status", values= c("Dominant_Singleton"=3,"Dominant_Ohnolog"=3, "Dominant_UPP"=3, "Dominant_SSD"=3,  "Recessive_Singleton"=1,"Recessive_Ohnolog"=1, "Recessive_UPP"=1, "Recessive_SSD"=1,"Unknown_Singleton"=10,"Unknown_Ohnolog"=10, "Unknown_UPP"=10, "Unknown_SSD"=10), labels=c("Dominant Singleton", "Dominant Ohnolog", "Dominant SSD", "Recessive Singleton", "Recessive Ohnolog", "Recessive SSD", "Unknown Singleton", "Unknown Ohnolog", "Unknown SSD"))+
  #scale_y_continuous(limits=c(0,.5), breaks=seq(0,.5, by = .05), expand = c(0, 0))+
  scale_x_discrete(limits=c("1st","2nd","3rd","4th","5th","6th","7th","8th","9th","10th"))+
  theme(plot.title = element_text(lineheight=.8, face="bold"))+ 
  theme(axis.title.y=element_text(face = "bold", size=12))+
  theme(axis.title.x=element_text(face = "bold", size=12))+
  theme(panel.background = element_rect(fill = "white")) +
  theme(axis.text.x=element_text(angle=30,hjust=1,vjust=1),axis.text=element_text(size=10))+ 
  theme(plot.title = element_text(size = 14),panel.grid.major = element_line(colour="snow2", size = (0.2)), panel.grid.minor = element_line(size = (0.2), colour="snow1"), panel.background = element_rect(fill="transparent",colour=NA), legend.key.size = unit(0.7, "cm"), legend.position = c(.25, .8),legend.title=element_text(size=0), legend.text = element_text(size = 9), plot.margin=unit(c(0.5,0.5,0.5,0.5), "cm"))  + labs(linetype='custom title')
#-----------------MCA
library(ggfortify)
library(cluster)
library(plyr)
library("FactoMineR")
library("devtools")
library("factoextra")
library("corrplot")
disess2 <- read.csv("Last_Hurrah_Dataset.csv") [c(3,4,5,7)]
test = na.omit(disess2)
test<-subset(test, Ortho_Age!="0")
test$Disease_Association[test$Disease_Association == "Both"] <- "Dominant"
test$Ortho_Age <- sub("^", "Ortho_Age.", test$Ortho_Age )
test$Ortho_Age <- factor(test$Ortho_Age)
test <- na.omit(test)
test$decile2 <- test$Haplosufficiency_Rank
n = 0
for (i in test$decile2){
  n = n+1
  if(i <= 10){
    test$decile2[n]<-"10th"  
  }
  if ((i > 10) && (i <=20)){
    test$decile2[n]<-"9th"
  }  
  if ((i > 20) && (i <=30)){
    test$decile2[n]<-"8th"
  }  
  if ((i > 30) && (i <=40)){
    test$decile2[n]<-"7th"
  }  
  if ((i > 40) && (i <=50)){
    test$decile2[n]<-"6th"
  }  
  if ((i > 50) && (i <=60)){
    test$decile2[n]<-"5th"
  }  
  if ((i > 60) && (i <=70)){
    test$decile2[n]<-"4th"
  }  
  if ((i > 70) && (i <=80)){
    test$decile2[n]<-"3rd"
  }  
  if ((i > 80) && (i <=90)){
    test$decile2[n]<-"2nd"
  }  
  if ((i > 90) && (i <=100)){
    test$decile2[n]<-"1st"
  }  
}
test$decile2 <- factor(test$decile2, levels = c("1st","2nd","3rd","4th","5th","6th","7th","8th","9th","10th"))
test <- test[, c(1,2,4,5)]
test.mca <- MCA(test, graph = FALSE)
Haplo <- as.factor(test$decile2)
head(Haplo)

TWO<- fviz_mca_ind(test.mca, label ="none", habillage=Haplo, addEllipses = TRUE, ellipse.level = 0.95)+ggtitle("b")+ 
  theme(axis.title.y=element_text(face = "bold", size=14))+
  theme(axis.title.x=element_text(face = "bold", size=14))+
  theme(plot.title = element_text(size = 13, lineheight=.8, face="bold"),panel.grid.major = element_line(colour="snow2", size = (0.2)), panel.grid.minor = element_line(size = (0.2), colour="snow1"), panel.background = element_rect(fill="transparent",colour=NA),  legend.key.size = unit(0.6, "cm"), legend.position = c(.9, .8),legend.title=element_text(size=11), legend.text = element_text(size = 9), plot.margin=unit(c(0.5,0.5,0.5,0.5), "cm"))
library(gridExtra)
library(lattice)
library(grid) 
grid.arrange(arrangeGrob(ONE,TWO, ncol=7, nrow=1,layout_matrix = rbind(c(1,1,1,2,2,2,2))))


#-----------------------------------------------------------------------------------------------------------

#Family analysis

#Looking at families of 4 or more members, 
#(2nd run limit to  - in which the oldest members are at least 105MYO - corresponding with the peak that we see in gene accumulation at the time of the branching of placental mammals - Thus excluding families made entirely of "young" genes) 
#Take the "oldest" 2 genes in each family 

#What is the spread of ages across these genes (a-la figure 1)
#ascertain the proportions of Dominant, recessive and "both" disease in these genes
#ascertain the box and whisker spread of haploinsufficiency within these genes
  
 



#ascertain the proportions of Dominant, recessive and "both" disease in these genes
rm(list=ls())
library(ggfortify)
library(plyr)
library(tidyr)
library(gridExtra)
library(grid) 
setwd("/Users/mqbpqam6/Desktop/PhD/Last_Hurrah_Data_Gen/")
disess1 <- read.csv("Last_Hurrah_Old_Genes.csv") [c(10)]	#Path to DS1
disess1 = na.omit(disess1)

nd<-subset(disess1, Disease_Association=="None")
disess1<-subset(disess1, Disease_Association!="None")

nd$type<-"None"
disess1$type<-"Disease"

disess <- rbind(disess1, nd)


ADS =count(disess,c("Disease_Association","type"))
tots <- sum(ADS$freq)
ADS$Norm <- tots
ADS2 <- transform(ADS, normalized_freq = freq / tots)

group.colours <- c("Unknown" = "deepskyblue4", "Dominant" = "black", "Recessive" = "mediumvioletred", "Both" = "darkseagreen", "None"= "grey92")

F1 <- ggplot(data=ADS2, aes(x=type, y=normalized_freq, fill=Disease_Association)) +
  geom_bar(stat="identity") + 
  ggtitle("A. Old family genes") +
  ylab(label="") +
  xlab(label="") +
  scale_y_continuous(limits=c(0,0.8), breaks=seq(0,0.8, by = 0.2), expand = c(0, 0))+
  theme(plot.title = element_text(lineheight=.8, face="bold"))+ theme(panel.background = element_rect(fill = "white")) + theme(axis.text.x=element_text(angle=30,hjust=1,vjust=1),axis.text=element_text(size=10))+ 
  theme(axis.title.y=element_text(face = "bold"))+ theme(axis.title.x=element_text(face = "bold"))+
  scale_fill_manual(values=group.colours)+ 
  theme(legend.background = element_rect(fill="white", size=0.5, linetype="solid", colour ="grey86"))+
  theme(plot.title = element_text(size = 12),panel.grid.major = element_line(colour="white", size = (0.2)), panel.grid.minor = element_line(size = (0.2), colour="white"), panel.background = element_rect(fill = 'white', colour = 'transparent'), plot.margin=unit(c(0.5,0.5,0,0.2), "cm"))+
  theme(legend.position="NONE")
F1 <- F1 + guides(fill=guide_legend(title="Disease Association"))     
F1

rm(list= ls()[!(ls() %in% c('F1'))])
library(ggfortify)
library(plyr)
library(tidyr)
library(gridExtra)
library(grid) 
setwd("/Users/mqbpqam6/Desktop/PhD/Last_Hurrah_Data_Gen/")
disess1 <- read.csv("Last_Hurrah_Dataset.csv") [c(7)]	#Path to DS1
disess1 = na.omit(disess1)

nd<-subset(disess1, Disease_Association=="None")
disess1<-subset(disess1, Disease_Association!="None")

nd$type<-"None"
disess1$type<-"Disease"

disess <- rbind(disess1, nd)


ADS =count(disess,c("Disease_Association","type"))
tots <- sum(ADS$freq)
ADS$Norm <- tots
ADS2 <- transform(ADS, normalized_freq = freq / tots)

group.colours <- c("Unknown" = "deepskyblue4", "Dominant" = "black", "Recessive" = "mediumvioletred", "Both" = "darkseagreen", "None"= "grey92")

F2 <- ggplot(data=ADS2, aes(x=type, y=normalized_freq, fill=Disease_Association)) +
  geom_bar(stat="identity") + 
  ggtitle("B. All genes") +
  ylab(label="") +
  xlab(label="") +
  scale_y_continuous(limits=c(0,0.8), breaks=seq(0,0.8, by = 0.2), expand = c(0, 0))+
  theme(plot.title = element_text(lineheight=.8, face="bold"))+ theme(panel.background = element_rect(fill = "white")) + theme(axis.text.x=element_text(angle=30,hjust=1,vjust=1),axis.text=element_text(size=10))+ 
  theme(axis.title.y=element_text(face = "bold"))+ theme(axis.title.x=element_text(face = "bold"))+
  scale_fill_manual(values=group.colours)+ 
  theme(legend.background = element_rect(fill="white", size=0.5, linetype="solid", colour ="grey86"))+
  theme(plot.title = element_text(size = 12),panel.grid.major = element_line(colour="white", size = (0.2)), panel.grid.minor = element_line(size = (0.2), colour="white"), panel.background = element_rect(fill = 'white', colour = 'transparent'), plot.margin=unit(c(0.5,0.5,0,0.2), "cm"))+
  theme(legend.position=c(0.75, 0.75))
F2 <- F2 + guides(fill=guide_legend(title="Disease Association"))     
F2

grid.arrange(F1,F2, ncol=2, nrow=1, layout_matrix = rbind(c(1,2)),vp=viewport(width=0.93, height=0.93),bottom = textGrob("Disease-branch association within paralog pairs",gp=gpar(fontsize=10,font=3)),left = textGrob("Proportion of genes",rot=90,gp=gpar(fontsize=10,font=3)))




#ascertain the violin spread of haploinsufficiency within these genes

rm(list=ls())
library(ggfortify)
library(plyr)
library(tidyr)
library(gridExtra)
library(Hmisc) 
library(ggplot2)

setwd("/Users/mqbpqam6/Desktop/PhD/Last_Hurrah_Data_Gen/")
hap1 <- read.csv("Last_Hurrah_Old_Genes.csv") [c(4,8,7,10,12)]	#Path to DS1
hap1 = na.omit(hap1)
group.colours <- c("Ohnolog" = "#0072B2", "SSD" = "#CC79A7")

# Basic violin plot
o <- ggplot(hap1, aes(x=Paralog_status, y=Haplosufficiency_Rank, fill = Paralog_status)) + 
  geom_violin() + theme_classic() +
  ggtitle("A.") +
  ylab(label="") +
  xlab(label="") +
  scale_fill_manual(values=group.colours)+
  scale_y_continuous(limits=c(-15,120), breaks=seq(0,100, by = 10))+
  stat_summary(fun.data=mean_sdl, geom="pointrange", color="black") +
  theme(legend.position="NONE")+
  geom_hline(yintercept=10, color = "black", size=1, linetype="dashed")
o

hap2 <- read.csv("Last_Hurrah_Dataset.csv") [c(1,5,4,7,9)]	#Path to DS1
hap2 = na.omit(hap2)
hap2<-subset(hap2, Paralog_status!="Singleton")
hap3<-subset(hap2, !(Gene_id %in% hap1$Gene_id))
hap4<-subset(hap3, (Family_ID %in% hap1$Family_ID))

# Basic violin plot
a <- ggplot(hap4, aes(x=Paralog_status, y=Haplosufficiency_Rank, fill = Paralog_status)) + 
  geom_violin() + theme_classic()+
  ggtitle("B.") +
  ylab(label="") +
  xlab(label="") +
  scale_fill_manual(values=group.colours)+
  scale_y_continuous(limits=c(-15,120), breaks=seq(0,100, by = 10))+
  stat_summary(fun.data=mean_sdl, geom="pointrange", color="black") +
  theme(legend.position="NONE")+
  geom_hline(yintercept=10, color = "black", size=1, linetype="dashed")
a

grid.arrange(o,a, ncol=1, nrow=2, layout_matrix = rbind(c(1),c(2)),vp=viewport(width=0.93, height=0.93),bottom = textGrob("Paralog Status",gp=gpar(fontsize=10,font=3)),left = textGrob("Haplosufficiency rank",rot=90,gp=gpar(fontsize=10,font=3)))


#########################################################################################################

#ascertain the spread of disease within these genes

rm(list=ls())
library(plyr)
library(gridExtra)
library(ggplot2)

setwd("/Users/mqbpqam6/Desktop/PhD/Last_Hurrah_Data_Gen/")
hap1 <- read.csv("Last_Hurrah_Old_Genes.csv") [c(4,7,10,12)]	#Path to DS1
hap1 = na.omit(hap1)
group.colours <- c("Unknown" = "deepskyblue4", "Dominant" = "black", "Recessive" = "mediumvioletred", "Both" = "darkseagreen", "None"= "grey92")

ADS =count(hap1,c("Disease_Association","Paralog_status"))

tots2 <- sum(ADS$freq)
ADS$Norm <- tots2
ADS <- transform(ADS, normalized_freq = freq / tots2)
#levels(ADS$Paralog_status)[levels(ADS$Paralog_status)=="Ohnolog"]<-"Ohnolog-25"
#levels(ADS$Paralog_status)[levels(ADS$Paralog_status)=="SSD"]<-"SSD-25"
ADS$ON<-"Oldest ~25%"

hap2 <- read.csv("Last_Hurrah_Dataset.csv") [c(1,4,7,9)]	#Path to DS1
hap2 = na.omit(hap2)
hap2<-subset(hap2, Paralog_status!="Singleton")
hap3<-subset(hap2, !(Gene_id %in% hap1$Gene_id))
hap4<-subset(hap3, (Family_ID %in% hap1$Family_ID))

BDS =count(hap4,c("Disease_Association","Paralog_status"))
tots2 <- sum(BDS$freq)
BDS$Norm <- tots2
BDS <- transform(BDS, normalized_freq = freq / tots2)
#levels(BDS$Paralog_status)[levels(BDS$Paralog_status)=="Ohnolog"]<-"Ohnolog-75"
#levels(BDS$Paralog_status)[levels(BDS$Paralog_status)=="SSD"]<-"SSD-75"
BDS$ON<-"Younger ~75%"

CDS<-rbind(ADS, BDS)

CDS$Disease_Association <- relevel(CDS$Disease_Association, "Unknown")
CDS$Disease_Association <- relevel(CDS$Disease_Association, "None")

o <-ggplot(CDS, aes(x = Paralog_status, y = normalized_freq, fill=Disease_Association)) + 
  geom_bar(stat="identity") +
  facet_grid(. ~ ON)+
  theme_classic()+
  ggtitle("") +
  ylab(label="Proportion of total genes in group") +
  xlab(label="Paralog Status") +
  scale_fill_manual(values=group.colours)+
  labs(fill = "Disease")
o

###########################################

#Are gene families with an identifiably early duplication point more likely to have higher proportions of disease? 
  
#First of all for overall disease, then secondly for individual disease subsets

#Violin plot this? 

rm(list=ls())
library(plyr)
library(gridExtra)
library(ggplot2)

setwd("/Users/mqbpqam6/Desktop/PhD/Last_Hurrah_Data_Gen/")
fam1 <- read.csv("Last_Hurrah_Family_data.csv") [c(10:13,15)]	#Path to DS1
fam1<-subset(fam1, (Family_Root != "Singleton"))

#Combine totals for columns 10-13
fam1$tots <- rowSums(fam1[1:4])


#Are gene families with an identifiably early duplication point more likely to have higher proportions of disease? 
#Violin plot

#levels(split1$virus_classification)[levels(split1$virus_classification)=="1"]<-"VIP"

#group.colours <- c("Ohnolog" = "#0072B2", "SSD" = "#CC79A7")
fam1$Family_Root <- factor(fam1$Family_Root, levels = c("LECA","Opisthokonta","Bilateria","Chordata","Vertebrata","Gnathostomata","Euteleostomi","Sarcopterygii","Tetrapoda","Amniota","Mammalia","Theria","Eutheria","Boreoeutheria","Euarchontoglires","Primates","Haplorrhini","Simiiformes","Catarrhini","Hominoidea","Hominidae","Homininae","Homo sapiens","Singleton"))

# Basic violin plot
o <- ggplot(fam1, aes(x=Family_Root, y=tots, fill = Family_Root)) + 
  geom_violin() + theme_classic() +
  ggtitle("") +
  ylab(label="") +
  xlab(label="") +
  #scale_fill_manual(values=group.colours)+
  scale_y_continuous(limits=c(-0.5,1.7), breaks=seq(-.4,1.6, by = 0.2))+
  stat_summary(fun.data=mean_sdl, geom="pointrange", color="black") +
  theme(legend.position="NONE")+
  geom_hline(yintercept=0.75, color = "red", size=1, linetype="dashed")+
  coord_flip()
  #geom_hline(yintercept=10, color = "black", size=1, linetype="dashed")
o











#DNDS SCATTER PLOT
rm(list= ls()[!(ls() %in% c('F1'))])
library(ggfortify)
library(plyr)
library(tidyr)
library(gridExtra)
library(grid) 
setwd("/Users/mqbpqam6/Desktop/PhD/Last_Hurrah_Data_Gen/")
scaplot1 <- read.csv("Last_Hurrah_Ohnolog_DNDS_Pairs.csv") [c(3:5)]

all<-ggplot(scaplot1, aes(x=dNdS1, y=dNdS2)) +
  geom_point(size=2, shape=23)


split<-subset(scaplot1, dNdS1 < 0.01 & dNdS2 < 0.01)


min<-ggplot(split, aes(x=dNdS1, y=dNdS2)) +
  geom_point(size=2, shape=23)


grid.arrange(arrangeGrob(all,min, ncol=2, nrow=1,layout_matrix = rbind(c(1,2))))


d <- density(scaplot1$dNdS_diff)
a<-plot(d)

d2 <- density(split$dNdS_diff)
b<-plot(d2)




