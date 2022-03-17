######RESULTS PART 2.1########

#Critical nodes are often VIPs
rm(list=ls())
library(plyr)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
d1 <- read.csv("Masterfile_14Feb.csv", header = T)[c(9,16,24)]
d1 <- subset(d1, Disease_Association!="NA") #restrict to genes in VR and LH
d1<-d1[c(2:3)]
d1$VIP_type<-as.character(d1$VIP_type)
d1$VIP_type[is.na(d1$VIP_type)] <- "NONE"
d1$consol_vip<-d1$VIP_type
d1$consol_vip<-as.character(d1$consol_vip)
d1$consol_vip[d1$consol_vip != "NONE"] <- "VIP"

criticals<-nrow(subset(d1, control.category == "critical"))
intermittents<-nrow(subset(d1, control.category == "intermittent"))
redundants<-nrow(subset(d1, control.category == "redundant"))
all_VIPs<-nrow(subset(d1, VIP_type != "NONE"))

cv =plyr::count(d1,c("consol_vip", "control.category"))
cv<-subset(cv, consol_vip != "NONE")

crit_vip<-subset(cv,control.category == "critical")[c(3)]
int_vip<-subset(cv,control.category == "intermittent")[c(3)]
red_vip<-subset(cv,control.category == "redundant")[c(3)]
crit_vip<-crit_vip[1,1]
int_vip<-int_vip[1,1]
red_vip<-red_vip[1,1]

observations  <- c( crit_vip, int_vip, red_vip)
totals <- c( criticals, intermittents, redundants)
pairwise.prop.test(observations, totals, p.adjust.method="bonferroni")

#Critical VIPs are Multi
d2<-subset(d1, control.category == "critical" & consol_vip == "VIP")
control_vips<-nrow(d2)
cv2 =plyr::count(d2,c("VIP_type"))
cv2$total<-control_vips

observations2  <- cv2$freq
totals2 <- cv2$total
pairwise.prop.test(observations2, totals2, p.adjust.method="bonferroni")

#No significance for intermittents and viruses
d2<-subset(d1, control.category == "intermittent" & consol_vip == "VIP")
control_vips<-nrow(d2)
cv2 =plyr::count(d2,c("VIP_type"))
cv2$total<-control_vips

observations2  <- cv2$freq
totals2 <- cv2$total
pairwise.prop.test(observations2, totals2, p.adjust.method="bonferroni")


#--------------------------------
#VIPs and disease genes intersect

rm(list=ls())
library(plyr)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

d1 <- read.csv("Masterfile_14Feb.csv", header = T)[c(9,24)]
d1 <- subset(d1, Disease_Association!="NA") #restrict to genes in VR and LH
d1$VIP_type<-as.character(d1$VIP_type)
d1$VIP_type[is.na(d1$VIP_type)] <- "NONE"
d1$consol_vip<-d1$VIP_type
d1$consol_vip<-as.character(d1$consol_vip)
d1$consol_vip[d1$consol_vip != "NONE"] <- "VIP"
d1$consol_dis<-d1$Disease_Association
d1$consol_dis<-as.character(d1$consol_dis)
d1$consol_dis[d1$consol_dis != "None"] <- "Disease"

###Retrieve a random samle of N genes (where N is the observed number VIP genes) from the full data
#Get distribution of randomly sampled intersecting genes -how many genes are found in each intersect
#####VIPs#####

#Observed totals
vips <- nrow(subset(d1, consol_vip == "VIP" )) #VIP genes
vip_dis <- nrow(subset(d1, consol_vip == "VIP" & consol_dis == "Disease" )) #observed intersection

interout<-data.frame(matrix(NA, nrow = 10000, ncol = 1)) #make an empty dataframe with the capacity to hold 10,000 rows
sampledata2<-d1[c(3,4)] #just virus and disease association

set.seed(35)
library(dplyr)

for(i in 1:10000){ #for 10000 replicates
  df<-sample_n(sampledata2, size = vips) #randomly sample data from "sampledata" ~ the same number of rows as there are VIPs in the real data
  uniqinter <-nrow(subset(df, consol_dis == "Disease")) #Count the number of genes that were found to have a Heritable disease association.
  interout[i,1]<-(uniqinter) #add to the output data the number of Uniquely VIP genes found
}

names(interout)[1] <- "X1"
#Plot VIP&disease intersect distribution from random sampling
istandc<-sd(interout$X1) #standard deviation
imnc<-mean(interout$X1) #mean
is1c<-(imnc-istandc) #one standard deviation below the mean
is2c<-(imnc+istandc) #one standard deviation above the mean

b <- ggplot(interout, aes(x = X1))+
  ggtitle("INTERSECT.") +
  ylab(label="") +
  xlab(label="")+ 
  theme(plot.title = element_text(lineheight=.8, face="bold"))+
  theme(panel.grid.major = element_line(colour="white", size = (0.2)), panel.grid.minor = element_line(size = (0.2), colour="white"), panel.background = element_rect(fill = 'white', colour = 'transparent'))+
  geom_density(aes(y = ..count..), fill = "aliceblue") +
  geom_vline(aes(xintercept = mean(X1)), 
             linetype = "dashed", size = 0.6,
             color = "black")+
  geom_vline(aes(xintercept = is1c), 
             linetype = "dashed", size = 0.6,
             color = "blue")+
  geom_vline(aes(xintercept = is2c), 
             linetype = "dashed", size = 0.6,
             color = "blue")+
  geom_vline(aes(xintercept = vip_dis), 
             linetype = "dashed", size = 0.6,
             color = "red")

b

#Calculate significance (number >= obs / replicates)
pval<-(nrow(subset(interout, X1 >= vip_dis)))/10000
top_score<-max(interout$X1)

#Dominant specific
vips <- nrow(subset(d1, consol_vip == "VIP" )) #VIP genes
vip_dis <- nrow(subset(d1, consol_vip == "VIP" & Disease_Association == "Dominant" )) #observed intersection

interout<-data.frame(matrix(NA, nrow = 10000, ncol = 1)) #make an empty dataframe with the capacity to hold 10,000 rows
sampledata2<-d1[c(1,3)] #just virus and disease association

set.seed(35)
library(dplyr)

for(i in 1:10000){ #for 10000 replicates
  df<-sample_n(sampledata2, size = vips) #randomly sample data from "sampledata" ~ the same number of rows as there are VIPs in the real data
  uniqinter <-nrow(subset(df, Disease_Association == "Dominant")) #Count the number of genes that were found to have a Heritable disease association.
  interout[i,1]<-(uniqinter) #add to the output data the number of Uniquely VIP genes found
}
names(interout)[1] <- "X1"

#Calculate significance (number >= obs / replicates)
pval<-(nrow(subset(interout, X1 >= vip_dis)))/10000
top_score<-max(interout$X1)

#Recessive specific
vips <- nrow(subset(d1, consol_vip == "VIP" )) #VIP genes
vip_dis <- nrow(subset(d1, consol_vip == "VIP" & Disease_Association == "Recessive" )) #observed intersection

interout<-data.frame(matrix(NA, nrow = 10000, ncol = 1)) #make an empty dataframe with the capacity to hold 10,000 rows
sampledata2<-d1[c(1,3)] #just virus and disease association

set.seed(35)
library(dplyr)

for(i in 1:10000){ #for 10000 replicates
  df<-sample_n(sampledata2, size = vips) #randomly sample data from "sampledata" ~ the same number of rows as there are VIPs in the real data
  uniqinter <-nrow(subset(df, Disease_Association == "Recessive")) #Count the number of genes that were found to have a Heritable disease association.
  interout[i,1]<-(uniqinter) #add to the output data the number of Uniquely VIP genes found
}
names(interout)[1] <- "X1"

#Calculate significance (number >= obs / replicates)
pval<-(nrow(subset(interout, X1 >= vip_dis)))/10000
top_score<-max(interout$X1)

###Within the disease VIPs is there a disproportionate number of dominant?
#Pairwise prop test Dom/Rec/Both/Unknown
vd<-subset(d1, consol_vip == "VIP" & consol_dis == "Disease")
cvd <-vd %>% count(Disease_Association)
cvd$total<-nrow(vd)

observations2  <- cvd$n
totals2 <- cvd$total
pairwise.prop.test(observations2, totals2, p.adjust.method="bonferroni")


#----------------------------
#Disease Genes tend to be ancient

rm(list=ls())
library(dplyr)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
d1 <- read.csv("Masterfile_14Feb.csv", header = T)[c(6,9)]

d1 <- subset(d1, Disease_Association!="NA") #restrict to genes in VR and LH
d1 <- na.omit(d1) #Remove any rows with null assignments from the analysis
d1$Disease_Association <- relevel(d1$Disease_Association, "Recessive", "Both", "Dominant", "Unknown", "None")
d1$Disease_Association[d1$Disease_Association == "Both"] <- "Dominant"
d1$consol_dis<-d1$Disease_Association
d1$consol_dis<-as.character(d1$consol_dis)
d1$consol_dis[d1$consol_dis != "None"] <- "Disease"
d1$Age_group <- d1$Ortho_Age
d1$Age_group[d1$Age_group==435] <- "WGD & Earlier"
d1$Age_group[d1$Age_group==615] <- "WGD & Earlier"
d1$Age_group[d1$Age_group==796] <- "WGD & Earlier"
d1$Age_group[d1$Age_group!="WGD & Earlier"] <- "Since WGD"

old_dis<-nrow(subset(d1, consol_dis == "Disease" & Age_group == "WGD & Earlier"))
all_dis<-nrow(subset(d1, consol_dis == "Disease"))
old_all<-nrow(subset(d1, Age_group == "WGD & Earlier"))
all<-nrow(d1)

observations2  <- c(old_dis,old_all)
totals2 <- c(all_dis,all)
prop.test(observations2, totals2)


#-----------------------------------------
#Duplicated genes (Ohnologs) are often disease associated
rm(list=ls())
library(dplyr)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
d1 <- read.csv("Masterfile_14Feb.csv", header = T)[c(7,9)]

d1 <- subset(d1, Disease_Association!="NA") #restrict to genes in VR and LH
d1 <- na.omit(d1) #Remove any rows with null assignments from the analysis
d1$Disease_Association <- relevel(d1$Disease_Association, "Recessive", "Both", "Dominant", "Unknown", "None")

d1$consol_dis<-d1$Disease_Association
d1$consol_dis<-as.character(d1$consol_dis)
d1$consol_dis[d1$consol_dis != "None"] <- "Disease"

rec<-nrow(subset(d1, Disease_Association == "Recessive"))
both<-nrow(subset(d1, Disease_Association == "Both"))
dom<-nrow(subset(d1, Disease_Association == "Dominant"))
unk<-nrow(subset(d1, Disease_Association == "Unknown"))
none<-nrow(subset(d1, Disease_Association == "None"))
all<-nrow(d1)

####Ohnologs
od<-subset(d1, Paralog_status == "Ohnolog")
o_rec<-nrow(subset(od, Disease_Association == "Recessive"))
o_both<-nrow(subset(od, Disease_Association == "Both"))
o_dom<-nrow(subset(od, Disease_Association == "Dominant"))
o_unk<-nrow(subset(od, Disease_Association == "Unknown"))
o_none<-nrow(subset(od, Disease_Association == "None"))
o_all<-nrow(od)

#Recessive
observations  <- c(o_rec,rec)
totals <- c(o_all,all)
prop.test(observations, totals)

#Both
observations  <- c(o_both,both)
totals <- c(o_all,all)
prop.test(observations, totals)

#Dominant
observations  <- c(o_dom,dom)
totals <- c(o_all,all)
prop.test(observations, totals)

#Unknown
observations  <- c(o_unk,unk)
totals <- c(o_all,all)
prop.test(observations, totals)

#None
observations  <- c(o_none,none)
totals <- c(o_all,all)
prop.test(observations, totals)

####SSDs
od<-subset(d1, Paralog_status == "SSD")
o_rec<-nrow(subset(od, Disease_Association == "Recessive"))
o_both<-nrow(subset(od, Disease_Association == "Both"))
o_dom<-nrow(subset(od, Disease_Association == "Dominant"))
o_unk<-nrow(subset(od, Disease_Association == "Unknown"))
o_none<-nrow(subset(od, Disease_Association == "None"))
o_all<-nrow(od)

#Recessive
observations  <- c(o_rec,rec)
totals <- c(o_all,all)
prop.test(observations, totals)

#Both
observations  <- c(o_both,both)
totals <- c(o_all,all)
prop.test(observations, totals)

#Dominant
observations  <- c(o_dom,dom)
totals <- c(o_all,all)
prop.test(observations, totals)

#Unknown
observations  <- c(o_unk,unk)
totals <- c(o_all,all)
prop.test(observations, totals)

#None
observations  <- c(o_none,none)
totals <- c(o_all,all)
prop.test(observations, totals)

####Singletons
od<-subset(d1, Paralog_status == "Singleton")
o_rec<-nrow(subset(od, Disease_Association == "Recessive"))
o_both<-nrow(subset(od, Disease_Association == "Both"))
o_dom<-nrow(subset(od, Disease_Association == "Dominant"))
o_unk<-nrow(subset(od, Disease_Association == "Unknown"))
o_none<-nrow(subset(od, Disease_Association == "None"))
o_all<-nrow(od)

#Recessive
observations  <- c(o_rec,rec)
totals <- c(o_all,all)
prop.test(observations, totals)

#Both
observations  <- c(o_both,both)
totals <- c(o_all,all)
prop.test(observations, totals)

#Dominant
observations  <- c(o_dom,dom)
totals <- c(o_all,all)
prop.test(observations, totals)

#Unknown
observations  <- c(o_unk,unk)
totals <- c(o_all,all)
prop.test(observations, totals)

#None
observations  <- c(o_none,none)
totals <- c(o_all,all)
prop.test(observations, totals)


######RESULTS PART 2.2########

#Is there a link between driver node status and paralog status?
rm(list=ls())
library(dplyr)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
d1 <- read.csv("Masterfile_14Feb.csv", header = T)[c(7,9,16)]
d1 <- subset(d1, Disease_Association!="NA") #restrict to genes in VR and LH
d1<-d1[c(1,3)]
d1 <- na.omit(d1) #Remove any rows with null assignments from the analysis

ohnologs<-nrow(subset(d1, Paralog_status == "Ohnolog"))
ssds<-nrow(subset(d1, Paralog_status == "SSD"))
singletons<-nrow(subset(d1, Paralog_status == "Singleton"))
all<-nrow(d1)

###Critical
od<-subset(d1, control.category == "critical")
o_ohn<-nrow(subset(od, Paralog_status == "Ohnolog"))
o_ssd<-nrow(subset(od, Paralog_status == "SSD"))
o_sin<-nrow(subset(od, Paralog_status == "Singleton"))
o_all<-nrow(od)

#Ohnolog
observations  <- c(o_ohn,ohnologs)
totals <- c(o_all,all)
prop.test(observations, totals)

#SSD
observations  <- c(o_ssd,ssds)
totals <- c(o_all,all)
prop.test(observations, totals)

#Singleton
observations  <- c(o_sin,singletons)
totals <- c(o_all,all)
prop.test(observations, totals)

###Intermittent
od<-subset(d1, control.category == "intermittent")
o_ohn<-nrow(subset(od, Paralog_status == "Ohnolog"))
o_ssd<-nrow(subset(od, Paralog_status == "SSD"))
o_sin<-nrow(subset(od, Paralog_status == "Singleton"))
o_all<-nrow(od)

#Ohnolog
observations  <- c(o_ohn,ohnologs)
totals <- c(o_all,all)
prop.test(observations, totals)

#SSD
observations  <- c(o_ssd,ssds)
totals <- c(o_all,all)
prop.test(observations, totals)

#Singleton
observations  <- c(o_sin,singletons)
totals <- c(o_all,all)
prop.test(observations, totals)

###Redundant
od<-subset(d1, control.category == "redundant")
o_ohn<-nrow(subset(od, Paralog_status == "Ohnolog"))
o_ssd<-nrow(subset(od, Paralog_status == "SSD"))
o_sin<-nrow(subset(od, Paralog_status == "Singleton"))
o_all<-nrow(od)

#Ohnolog
observations  <- c(o_ohn,ohnologs)
totals <- c(o_all,all)
prop.test(observations, totals)

#SSD
observations  <- c(o_ssd,ssds)
totals <- c(o_all,all)
prop.test(observations, totals)

#Singleton
observations  <- c(o_sin,singletons)
totals <- c(o_all,all)
prop.test(observations, totals)


###########################################################################

#Is there a link between driver node status and paralog status?
rm(list=ls())
library(dplyr)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
d1 <- read.csv("Masterfile_14Feb.csv", header = T)[c(7,9,16)]
d1 <- subset(d1, Disease_Association!="NA") #restrict to genes in VR and LH
d1<-d1[c(1,3)]
d1 <- na.omit(d1) #Remove any rows with null assignments from the analysis

criticals<-nrow(subset(d1, control.category == "critical"))
ints<-nrow(subset(d1, control.category == "intermittent"))
reds<-nrow(subset(d1, control.category == "redundant"))
all<-nrow(d1)

###Ohnolog
od<-subset(d1, Paralog_status == "Ohnolog")
o_crit<-nrow(subset(od, control.category == "critical"))
o_int<-nrow(subset(od, control.category == "intermittent"))
o_red<-nrow(subset(od, control.category == "redundant"))
o_all<-nrow(od)

#Critical
observations  <- c(o_crit,criticals)
totals <- c(o_all,all)
prop.test(observations, totals)

#Intermittent
observations  <- c(o_int,ints)
totals <- c(o_all,all)
prop.test(observations, totals)

#Redundant
observations  <- c(o_red,reds)
totals <- c(o_all,all)
prop.test(observations, totals)

###SSD
od<-subset(d1, Paralog_status == "SSD")
o_crit<-nrow(subset(od, control.category == "critical"))
o_int<-nrow(subset(od, control.category == "intermittent"))
o_red<-nrow(subset(od, control.category == "redundant"))
o_all<-nrow(od)

#Critical
observations  <- c(o_crit,criticals)
totals <- c(o_all,all)
prop.test(observations, totals)

#Intermittent
observations  <- c(o_int,ints)
totals <- c(o_all,all)
prop.test(observations, totals)

#Redundant
observations  <- c(o_red,reds)
totals <- c(o_all,all)
prop.test(observations, totals)

###Singletons
od<-subset(d1, Paralog_status == "Singleton")
o_crit<-nrow(subset(od, control.category == "critical"))
o_int<-nrow(subset(od, control.category == "intermittent"))
o_red<-nrow(subset(od, control.category == "redundant"))
o_all<-nrow(od)

#Critical
observations  <- c(o_crit,criticals)
totals <- c(o_all,all)
prop.test(observations, totals)

#Intermittent
observations  <- c(o_int,ints)
totals <- c(o_all,all)
prop.test(observations, totals)

#Redundant
observations  <- c(o_red,reds)
totals <- c(o_all,all)
prop.test(observations, totals)

