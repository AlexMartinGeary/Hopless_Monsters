####################Non-Introgression analysis############################

#QUESTION 1: Do the virus and diseasome collide? 

rm(list=ls())
library(ggplot2)
library(plyr)
library(grid)
library(VennDiagram)
setwd("/Users/mqbpqam6/Desktop/NEW_VR_Data/February_2020/February_rerun/Statistical_analysis_and_interpretation")

d1<- read.csv("Masterfile_14Feb.csv", header = T, strip.white=TRUE)[c(6,7,9,24)]
d1 <- subset(d1, Disease_Association!="NA") #restrict to genes in VR and LH

d1$VIP_type<-as.factor(d1$VIP_type)
levels(d1$VIP_type)[levels(d1$VIP_type)=="DNA"]<-"VIP"
levels(d1$VIP_type)[levels(d1$VIP_type)=="RNA"]<-"VIP"
levels(d1$VIP_type)[levels(d1$VIP_type)=="Multi"]<-"VIP"
levels(d1$VIP_type)[levels(d1$VIP_type)=="NA"]<-"None"

levels(d1$Disease_Association)[levels(d1$Disease_Association)=="Recessive"]<-"Heritable_Disease"
levels(d1$Disease_Association)[levels(d1$Disease_Association)=="Dominant"]<-"Heritable_Disease"
levels(d1$Disease_Association)[levels(d1$Disease_Association)=="Both"]<-"Heritable_Disease"
levels(d1$Disease_Association)[levels(d1$Disease_Association)=="Unknown"]<-"Heritable_Disease"

vip<-nrow(subset(d1,VIP_type == "VIP"))
disease <-nrow(subset(d1, Disease_Association == "Heritable_Disease"))
VD <-nrow(subset(d1, Disease_Association == "Heritable_Disease" &VIP_type == "VIP"))

grid.newpage()
draw.pairwise.venn(area1 = vip, area2 = disease, cross.area = VD, category = c("VIP","Disease"),
                   fill = c("lightblue1", "steelblue3"), alpha = rep(0.5, 2), cat.pos = c(0,0),
                   cat.dist = rep(0.025, 2), scaled = TRUE)

###Retrieve a random sample of N genes (where N is the observed number of VIP genes) from the full data
#Get distribution of randomly sampled intersecting genes -how many genes are found

#####VIPs#####
library(dplyr)
set.seed(35)

#Observed totals in intersect
v1_obs <- VD

interout<-data.frame(matrix(NA, nrow = 10000, ncol = 3)) #make an empty dataframe with the capacity to hold 10,000 rows
sampledata2<-d1[c(3,4)] #just virus class and disease association

for(i in 1:10000){ #for 10000 replicates
  df<-sample_n(sampledata2, size = vip) #randomly sample data from "sampledata" ~ the same number of rows as there are vip genes in the real data
  uniqinter <-nrow(subset(df, Disease_Association == "Heritable_Disease" &VIP_type == "VIP" )) #Count the number of genes that were found to have a VIP association, and Heritable disease association.
  interout[i,]<-list(i, uniqinter) #add to the output data the number of VIPdisease genes found
}

vstandc<-sd(interout$X2) #standard deviation
vmnc<-mean(interout$X2) #mean
vs1c<-(vmnc-vstandc) #one standard deviation below the mean
vs2c<-(vmnc+vstandc) #one standard deviation above the mean

#plot Distribution:
a <- ggplot(interout, aes(x = X2))+
  ggtitle("VIP/Disease") +
  ylab(label="") +
  xlab(label="")+ 
  theme(plot.title = element_text(lineheight=.8, face="bold"))+
  theme(panel.grid.major = element_line(colour="white", size = (0.2)), panel.grid.minor = element_line(size = (0.2), colour="white"), panel.background = element_rect(fill = 'white', colour = 'transparent'))+
  geom_density(aes(y = ..count..), fill = "aliceblue") +
  geom_vline(aes(xintercept = mean(X2)), 
             linetype = "dashed", size = 0.6,
             color = "black")+
  geom_vline(aes(xintercept = vs1c), 
             linetype = "dashed", size = 0.6,
             color = "blue")+
  geom_vline(aes(xintercept = vs2c), 
             linetype = "dashed", size = 0.6,
             color = "blue")+
  geom_vline(aes(xintercept = v1_obs), 
             linetype = "dashed", size = 0.6,
             color = "red")

a

#Prop test to check observed proportion of VIP/Disease genes in the intersect against mean expected number
tests  <- c( vmnc, v1_obs)
totals <- c( vip, vip)
prop.test(tests, totals)

##############################################################################################

#QUESTION 1.2: Are WGD & Earlier genes in the intersect? 
d1$Ortho_Age<-as.factor(d1$Ortho_Age)
levels(d1$Ortho_Age)[levels(d1$Ortho_Age)=="796"]<-"WGD & Earlier"
levels(d1$Ortho_Age)[levels(d1$Ortho_Age)=="615"]<-"WGD & Earlier"
levels(d1$Ortho_Age)[levels(d1$Ortho_Age)=="435"]<-"WGD & Earlier"

age <-nrow(subset(d1, Ortho_Age == "WGD & Earlier")) #total "WGD & Earlier" genes
DA <-nrow(subset(d1, Disease_Association == "Heritable_Disease" & Ortho_Age == "WGD & Earlier" )) #Total disease associated WGD & Earlier genes
VA <-nrow(subset(d1, Ortho_Age == "WGD & Earlier" &VIP_type == "VIP")) #Total VIP associated WGD & Earlier genes
VAD <-nrow(subset(d1, Disease_Association == "Heritable_Disease" &VIP_type == "VIP" & Ortho_Age == "WGD & Earlier")) #Total genes that are both disease and VIP

grid.newpage()
grid.text('C)', x=.1, y=.1)
draw.triple.venn(area1 = vip, area2 = disease, area3 = age, VD,DA,VA,VAD, category = c("VIP","Disease","WGD & Earlier genes"),
                 fill = c("lightblue1", "steelblue3", "deepskyblue4"), alpha = rep(0.5, 3), cat.pos = c(-40, 40, 180),
                 cat.dist = c(0.05, 0.05, 0.025), scaled =TRUE)

###Retrieve a random samle of N genes (where N is the observed number of WGD & Earlier genes -age) from the full data
#Get distribution of randomly sampled intersecting genes -how many genes are found in each intersect
#####VIPs#####

#Observed totals
vd_obs <- VAD #Observed WGD & Earlier disease&VIP genes

interout<-data.frame(matrix(NA, nrow = 10000, ncol = 3)) #make an empty dataframe with the capacity to hold 10,000 rows
sampledata2<-d1[c(3,4)] #just virus class and disease association


set.seed(35)

for(i in 1:10000){ #for 10000 replicates
  df<-sample_n(sampledata2, size = age) #randomly sample data from "sampledata" ~ the same number of rows as there are WGD & Earlier genes in the real data
  uniqinter <-nrow(subset(df, Disease_Association == "Heritable_Disease" &VIP_type == "VIP" )) #Count the number of genes that were found to have a VIP association, and Heritable disease association.
  interout[i,]<-list(i, uniqinter) #add to the output data the number of Uniquely VIP genes found
}

#Plot VIP&disease intersect distribution from random sampling
istandc<-sd(interout$X2) #standard deviation
imnc<-mean(interout$X2) #mean
is1c<-(imnc-istandc) #one standard deviation below the mean
is2c<-(imnc+istandc) #one standard deviation above the mean

b <- ggplot(interout, aes(x = X2))+
  ggtitle("INTERSECT.") +
  ylab(label="") +
  xlab(label="")+ 
  theme(plot.title = element_text(lineheight=.8, face="bold"))+
  theme(panel.grid.major = element_line(colour="white", size = (0.2)), panel.grid.minor = element_line(size = (0.2), colour="white"), panel.background = element_rect(fill = 'white', colour = 'transparent'))+
  geom_density(aes(y = ..count..), fill = "aliceblue") +
  geom_vline(aes(xintercept = mean(X2)), 
             linetype = "dashed", size = 0.6,
             color = "black")+
  geom_vline(aes(xintercept = is1c), 
             linetype = "dashed", size = 0.6,
             color = "blue")+
  geom_vline(aes(xintercept = is2c), 
             linetype = "dashed", size = 0.6,
             color = "blue")+
  geom_vline(aes(xintercept = vd_obs), 
             linetype = "dashed", size = 0.6,
             color = "red")

b


#Prop test to check observed proportion of INTERSECTION_GENES against mean
tests  <- c( imnc, vd_obs)
totals <- c( age, age)
prop.test(tests, totals)

##############################################################################################


#QUESTION 1.3: Are Ohnolog genes in the intersect? 
ohnolog <-nrow(subset(d1, Paralog_status == "Ohnolog")) #Total number of genes that are Ohnologs

DO <-nrow(subset(d1, Disease_Association == "Heritable_Disease" & Paralog_status == "Ohnolog" )) #Total Ohnologs which are disease associated
VO <-nrow(subset(d1, Paralog_status == "Ohnolog" &VIP_type == "VIP")) #Total Ohnologs which are VIPs
VOD <-nrow(subset(d1, Disease_Association == "Heritable_Disease" &VIP_type == "VIP" & Paralog_status == "Ohnolog")) #Total Ohnologs that are both VIP & disease

grid.newpage()
grid.text('D)', x=.1, y=.1)
draw.triple.venn(area1 = vip, area2 = disease, area3 = ohnolog, VD,DO,VO,VOD, category = c("VIP","Disease","Ohnologs"),
                 fill = c("lightblue1", "steelblue3", "deepskyblue4"), alpha = rep(0.5, 3), cat.pos = c(-40, 40, 180),
                 cat.dist = c(0.05, 0.05, 0.025), scaled =TRUE)

###Retrieve a random samle of N genes (where N is the observed number of ohnolog genes -int) from the full data
#Get distribution of randomly sampled intersecting genes -how many genes are found in each intersect

#Observed totals
vd_obs <- VOD #Observed Ohnolog disease&VIP genes


#####Intersect#####
interout<-data.frame(matrix(NA, nrow = 10000, ncol = 3)) #make an empty dataframe with the capacity to hold 10,000 rows
sampledata2<-d1[c(3,4)] #just virus class and disease association


set.seed(35)
for(i in 1:10000){ #for 10000 replicates
  df<-sample_n(sampledata2, size = ohnolog) #randomly sample data from "sampledata" ~ the same number of rows as there are ohnolog genes in the real data
  uniqinter <-nrow(subset(df, Disease_Association == "Heritable_Disease" &VIP_type == "VIP" )) #Count the number of genes that were found to have a VIP association, but no Heritable disease association.
  interout[i,]<-list(i, uniqinter) #add to the output data the number of Uniquely VIP genes found
}

istandc<-sd(interout$X2) #standard deviation
imnc<-mean(interout$X2) #mean
is1c<-(imnc-istandc) #one standard deviation below the mean
is2c<-(imnc+istandc) #one standard deviation above the mean

b <- ggplot(interout, aes(x = X2))+
  ggtitle("INTERSECT.") +
  ylab(label="") +
  xlab(label="")+ 
  theme(plot.title = element_text(lineheight=.8, face="bold"))+
  theme(panel.grid.major = element_line(colour="white", size = (0.2)), panel.grid.minor = element_line(size = (0.2), colour="white"), panel.background = element_rect(fill = 'white', colour = 'transparent'))+
  geom_density(aes(y = ..count..), fill = "aliceblue") +
  geom_vline(aes(xintercept = mean(X2)), 
             linetype = "dashed", size = 0.6,
             color = "black")+
  geom_vline(aes(xintercept = is1c), 
             linetype = "dashed", size = 0.6,
             color = "blue")+
  geom_vline(aes(xintercept = is2c), 
             linetype = "dashed", size = 0.6,
             color = "blue")+
  geom_vline(aes(xintercept = vd_obs), 
             linetype = "dashed", size = 0.6,
             color = "red")

b


#Prop test to check observed proportion of INTERSECTION_GENES against mean
tests  <- c( imnc, vd_obs)
totals <- c( ohnolog, ohnolog)
prop.test(tests, totals)

##############################################################################################


#QUESTION 1.3: Are SSD genes VIP heavy? 
ssds <-nrow(subset(d1, Paralog_status == "SSD")) #Total number of SSDs

DS <-nrow(subset(d1, Disease_Association == "Heritable_Disease" & Paralog_status == "SSD" )) #Total number of Disease associated SSDs
VS <-nrow(subset(d1, Paralog_status == "SSD" &VIP_type == "VIP")) #Total number of VIP SSDs
VSD <-nrow(subset(d1, Disease_Association == "Heritable_Disease" &VIP_type == "VIP" & Paralog_status == "SSD")) #Total number of VIP and disease associated SSDs

grid.newpage()
grid.text('D)', x=.1, y=.1)
draw.triple.venn(area1 = vip, area2 = disease, area3 = ssds, VD,DS,VS,VSD, category = c("VIP","Disease","SSDs"),
                 fill = c("lightblue1", "steelblue3", "deepskyblue4"), alpha = rep(0.5, 3), cat.pos = c(-40, 40, 180),
                 cat.dist = c(0.05, 0.05, 0.025), scaled =TRUE)

###Retrieve a random samle of N genes (where N is the observed number of SSD genes) from the full data
#Get distribution of randomly sampled intersecting genes -how many genes are found in each intersect

#Observed totals
vd_obs <- VSD

#####Intersect#####
interout<-data.frame(matrix(NA, nrow = 10000, ncol = 3)) #make an empty dataframe with the capacity to hold 10,000 rows
sampledata2<-d1[c(3,4)] #just virus class and disease association

set.seed(35)
for(i in 1:10000){ #for 10000 replicates
  df<-sample_n(sampledata2, size = ssds) #randomly sample data from "sampledata" ~ the same number of rows as there are critical genes in the real data
  uniqinter <-nrow(subset(df, Disease_Association == "Heritable_Disease" &VIP_type == "VIP" )) #Count the number of genes that were found to have a VIP association, but no Heritable disease association.
  interout[i,]<-list(i, uniqinter) #add to the output data the number of Uniquely VIP genes found
}

istandc<-sd(interout$X2) #standard deviation
imnc<-mean(interout$X2) #mean
is1c<-(imnc-istandc) #one standard deviation below the mean
is2c<-(imnc+istandc) #one standard deviation above the mean

b <- ggplot(interout, aes(x = X2))+
  ggtitle("INTERSECT.") +
  ylab(label="") +
  xlab(label="")+ 
  theme(plot.title = element_text(lineheight=.8, face="bold"))+
  theme(panel.grid.major = element_line(colour="white", size = (0.2)), panel.grid.minor = element_line(size = (0.2), colour="white"), panel.background = element_rect(fill = 'white', colour = 'transparent'))+
  geom_density(aes(y = ..count..), fill = "aliceblue") +
  geom_vline(aes(xintercept = mean(X2)), 
             linetype = "dashed", size = 0.6,
             color = "black")+
  geom_vline(aes(xintercept = is1c), 
             linetype = "dashed", size = 0.6,
             color = "blue")+
  geom_vline(aes(xintercept = is2c), 
             linetype = "dashed", size = 0.6,
             color = "blue")+
  geom_vline(aes(xintercept = vd_obs), 
             linetype = "dashed", size = 0.6,
             color = "red")

b

#Prop test to check observed proportion of INTERSECTION_GENES against mean
tests  <- c( imnc, vd_obs)
totals <- c( ohnolog, ohnolog)
prop.test(tests, totals)


#Multi Venn Plot
grid.newpage()
pushViewport(viewport(layout=grid.layout(ncol=2, nrow =2)))
#Add dividing lines
grid.lines(x = unit(c(0.50, 0.50), "npc"),
           y = unit(c(0.05, 0.95), "npc"),
           default.units = "npc",
           arrow = NULL, name = NULL,
           gp=gpar(), draw = TRUE, vp = NULL)
grid.lines(y = unit(c(0.50, 0.50), "npc"),
           x = unit(c(0.05, 0.95), "npc"),
           default.units = "npc",
           arrow = NULL, name = NULL,
           gp=gpar(), draw = TRUE, vp = NULL)

pushViewport(viewport(layout.pos.row=1,layout.pos.col=1))
#grid.newpage()
grid.text('A)', x=.1, y=.1)
draw.pairwise.venn(area1 = vip, area2 = disease, cross.area = VD, category = c("VIP","Disease"),
                   fill = c("lightblue1", "steelblue3"), alpha = rep(0.5, 2), cat.pos = c(0,0),
                   cat.dist = rep(0.025, 2), scaled = TRUE)

popViewport()

pushViewport(viewport(layout.pos.row=1,layout.pos.col=2))
#grid.newpage()
grid.text('B)', x=.1, y=.1)
draw.triple.venn(area1 = vip, area2 = disease, area3 = age, VD,DA,VA,VAD, category = c("VIP","Disease","WGD & Earlier genes"),
                 fill = c("lightblue1", "steelblue3", "deepskyblue4"), alpha = rep(0.5, 3), cat.pos = c(-40, 40, 180),
                 cat.dist = c(0.05, 0.05, 0.025), scaled =TRUE)
popViewport()

pushViewport(viewport(layout.pos.row=2,layout.pos.col=1))
#grid.newpage()
grid.text('C)', x=.1, y=.1)
draw.triple.venn(area1 = vip, area2 = disease, area3 = ohnolog, VD,DO,VO,VOD, category = c("VIP","Disease","Ohnologs"),
                 fill = c("lightblue1", "steelblue3", "deepskyblue4"), alpha = rep(0.5, 3), cat.pos = c(-40, 40, 180),
                 cat.dist = c(0.05, 0.05, 0.025), scaled =TRUE)
popViewport()

pushViewport(viewport(layout.pos.row=2,layout.pos.col=2))
#grid.newpage()
grid.text('D)', x=.1, y=.1)
draw.triple.venn(area1 = vip, area2 = disease, area3 = ssds, VD,DS,VS,VSD, category = c("VIP","Disease","SSDs"),
                 fill = c("lightblue1", "steelblue3", "deepskyblue4"), alpha = rep(0.5, 3), cat.pos = c(-40, 40, 180),
                 cat.dist = c(0.05, 0.05, 0.025), scaled =TRUE)








#################Introgression######################
#QUESTION N: Are Eur genes in the intersect? 
int <-nrow(subset(d1, Introgression_Status == "EUR"))

DI <-nrow(subset(d1, Disease_Association == "Heritable_Disease" & Introgression_Status == "EUR" ))
VI <-nrow(subset(d1, Introgression_Status == "EUR" &VIP_type == "VIP"))
VID <-nrow(subset(d1, Disease_Association == "Heritable_Disease" &VIP_type == "VIP" & Introgression_Status == "EUR"))

grid.newpage()
draw.triple.venn(area1 = vip, area2 = disease, area3 = int, VD,DI,VI,VID, category = c("VIP","Disease","European Introgressed genes"),
                 fill = c("lightblue1", "steelblue3", "deepskyblue4"), alpha = rep(0.5, 3), cat.pos = c(-40, 40, 180),
                 cat.dist = c(0.05, 0.05, 0.025), scaled =TRUE)

#Find out the number of families you would expect to see within each category if the occurrence of each were random. 
#dplyr-sample_n
library(dplyr)
set.seed(35)


#Observed VIP & Disease
vd_obs <- VID


sampledata<-d1

###Retrieve a random samle of N genes (where N is the observed number of introgressed genes -int) from the full data
#Get distribution of randomly sampled intersecting genes -how many genes are found in each intersect
#####VIPs#####
interout<-data.frame(matrix(NA, nrow = 10000, ncol = 3)) #make an empty dataframe with the capacity to hold 10,000 rows
sampledata2<-sampledata[c(3,4)] #just virus class and disease association

for(i in 1:10000){ #for 10000 replicates
  df<-sample_n(sampledata2, size = int) #randomly sample data from "sampledata" ~ the same number of rows as there are critical genes in the real data
  uniqinter <-nrow(subset(df, Disease_Association == "Heritable_Disease" &VIP_type == "VIP" )) #Count the number of genes that were found to have a VIP association, but no Heritable disease association.
  interout[i,]<-list(i, uniqinter) #add to the output data the number of Uniquely VIP genes found
}

istandc<-sd(interout$X2) #standard deviation
imnc<-mean(interout$X2) #mean
is1c<-(imnc-istandc) #one standard deviation below the mean
is2c<-(imnc+istandc) #one standard deviation above the mean

b <- ggplot(interout, aes(x = X2))+
  ggtitle("INTERSECTs.") +
  ylab(label="") +
  xlab(label="")+ 
  theme(plot.title = element_text(lineheight=.8, face="bold"))+
  theme(panel.grid.major = element_line(colour="white", size = (0.2)), panel.grid.minor = element_line(size = (0.2), colour="white"), panel.background = element_rect(fill = 'white', colour = 'transparent'))+
  geom_density(aes(y = ..count..), fill = "aliceblue") +
  geom_vline(aes(xintercept = mean(X2)), 
             linetype = "dashed", size = 0.6,
             color = "black")+
  geom_vline(aes(xintercept = is1c), 
             linetype = "dashed", size = 0.6,
             color = "blue")+
  geom_vline(aes(xintercept = is2c), 
             linetype = "dashed", size = 0.6,
             color = "blue")+
  geom_vline(aes(xintercept = vd_obs), 
             linetype = "dashed", size = 0.6,
             color = "red")

b


#Prop test to check observed proportion of INTERSECTION_GENES against mean
tests  <- c( imnc, vd_obs)
totals <- c( int, int)
prop.test(tests, totals)

