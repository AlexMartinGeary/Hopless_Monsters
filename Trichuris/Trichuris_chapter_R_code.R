rm(list=ls())
library(plyr)
library(ggfortify)
library(gridExtra)
library(grid)

split <- read.csv("/Users/mqbpqam6/Desktop/Ogunkanbi_Complete/July_19_Paper_data_run/Last_Hurrah_Dataset.csv", header = T)[c(1,4)]		#Insert path to DS1
split <- na.omit(split) #Remove any rows with null assignments from the analysis
split <-subset(split, Paralog_status=="Ohnolog"|Paralog_status=="SSD")

split2 <- read.csv("/Users/mqbpqam6/Desktop/Ogunkanbi_Complete/July_19_Paper_data_run/TGF_For_R.csv", header = T)[c(1,4)]		#Insert path to DS1
split2 <- na.omit(split2) #Remove any rows with null assignments from the analysis

split<-split[!(split$Gene_id %in% split2$Gene_id),] #Drop TGFB genes

ADS =count(split,c("Paralog_status"))
tots <- sum(ADS$freq)
ADS$Norm <- tots
ADSN <- transform(ADS, normalized_freq = freq / tots)
ADSN$Paralog_status <- factor(ADSN$Paralog_status)
ADSN$Paralog_status <- factor(ADSN$Paralog_status, levels=rev(levels(ADSN$Paralog_status)))

ADSN <-subset(ADSN, Paralog_status=="Ohnolog"|Paralog_status=="SSD")

ADSZ<-ADSN
ADSZ$freq <- 0
ADSZ$Norm <- 0
ADSZ$normalized_freq <- 0
ADSZ <- ADSZ[c(1,2,4)]


ADS2 =count(split2,c("Paralog_status"))
tots <- sum(ADS2$freq)
ADS2$Norm <- tots
ADS2N <- transform(ADS2, normalized_freq = freq / tots)
ADS2N$Paralog_status <- factor(ADS2N$Paralog_status)

ADS2P <-subset(ADS2N, Paralog_status=="Ohnolog"|Paralog_status=="SSD")

New4 <-merge(x = ADSZ, y = ADS2P, by = c("Paralog_status"), all.x = TRUE)
New4[is.na(New4)] <- 0
New4 <- New4[c(1,4,5,6)]
New4$type<-"TGFb N=32"
ADSN$type<-"Genome N=16,036"
colnames(New4) <- c("Paralog_status", "freq", "Norm", "normalized_freq", "type")
ADS3 <- rbind(ADSN, New4)

group.colours <- c("TGFb N=32" = "#CC79A7", "Genome N=16,036" = "black")

o <-ggplot(ADS3, aes(x = Paralog_status, y = (normalized_freq*100), fill=type)) + 
  geom_bar(stat="identity",position = position_dodge()) +
  theme_classic()+
  ggtitle("") +
  ylab(label="Percentage of genes in group") +
  xlab(label="Paralog Status") +
  scale_fill_manual(values=group.colours)+
  labs(fill = "")
o

#
tests  <- c( 14, 5930)
totals <- c( 32, 16036)
prop.test(tests, totals)


#######TWO##########

rm(list=ls())
library(plyr)
library(ggfortify)
library(gridExtra)
library(grid)

#genome
split <- read.csv("/Users/mqbpqam6/Desktop/Ogunkanbi_Complete/July_19_Paper_data_run/Last_Hurrah_Dataset.csv", header = T)[c(1,3)]		#Insert path to DS1
split <- na.omit(split) #Remove any rows with null assignments from the analysis
split<-subset(split, Ortho_Age!="0")

#TGFb
split2 <- read.csv("/Users/mqbpqam6/Desktop/Ogunkanbi_Complete/July_19_Paper_data_run/TGF_For_R.csv", header = T)[c(1,3)]		#Insert path to DS1
split2 <- na.omit(split2) #Remove any rows with null assignments from the analysis

split<-split[!(split$Gene_id %in% split2$Gene_id),] #Drop TGFB genes

gn =nrow(split)
ADS =count(split,c("Ortho_Age"))
tots <- sum(ADS$freq)
ADS$Norm <- tots
ADSN <- transform(ADS, normalized_freq = freq / tots)
ADSN$Ortho_Age <- factor(ADSN$Ortho_Age)
ADSN$Ortho_Age <- factor(ADSN$Ortho_Age, levels=rev(levels(ADSN$Ortho_Age)))
ADSZ<-ADSN
ADSZ$Ortho_Age <- factor(ADSZ$Ortho_Age)
#ADSZ$Ortho_Age <- factor(ADSZ$Ortho_Age, levels=rev(levels(ADSN$Ortho_Age)))


ADS2 =count(split2,c("Ortho_Age"))
tots <- sum(ADS2$freq)
ADS2$Norm <- tots
ADS2N <- transform(ADS2, normalized_freq = freq / tots)
ADS2N$Ortho_Age <- factor(ADS2N$Ortho_Age)

ADS2N$type<-"TGFb N=32"
ADSZ$type<-"Genome N= 19,598"
ADS3 <- rbind(ADSZ, ADS2N)

group.colours <- c("TGFb N=32" = "#009E73", "Genome N= 19,598" = "black")

o <-ggplot(ADS3, aes(x = Ortho_Age, y = (normalized_freq*100), fill=type)) + 
  geom_bar(stat="identity",position = position_dodge2(width = 0.9, preserve = "single")) +
  theme_classic()+
  ggtitle("") +
  ylab(label="Percentage of genes in group") +
  xlab(label="Gene Age (MYA)") +
  scale_fill_manual(values=group.colours)+
  labs(fill = "")
o

#435
tests  <- c(11,3719)
totals <- c(32,19598)
prop.test(tests, totals)

#615
tests  <- c(12,4418)
totals <- c(32,19598)
prop.test(tests, totals)

#796
tests  <- c(1,4058)
totals <- c(32,19598)
prop.test(tests, totals)

#####THREE########
rm(list=ls())
library(plyr)
library(ggfortify)
library(gridExtra)
library(grid)
library(binr)

split <- read.csv("/Users/mqbpqam6/Desktop/Ogunkanbi_Complete/July_19_Paper_data_run/Last_Hurrah_Dataset.csv", header = T)[c(1,5)]		#Insert path to DS1
split <- na.omit(split) #Remove any rows with null assignments from the analysis

split_T <- read.csv("/Users/mqbpqam6/Desktop/Ogunkanbi_Complete/July_19_Paper_data_run/TGF_For_R.csv", header = T)[c(1,5)]		#Insert path to DS1
split_T <- na.omit(split_T) #Remove any rows with null assignments from the analysis

library(tidyverse)
split<-split %>%
  mutate(group = ntile(Haplosufficiency_Rank, 10))
detach("package:tidyverse", unload=TRUE)
detach("package:dplyr", unload=TRUE)
library(plyr)

split3<-split[(split$Gene_id %in% split_T$Gene_id),] #Keep TGFB genes
split<-split[!(split$Gene_id %in% split_T$Gene_id),] #Drop TGFB genes

split$group<- as.factor(split$group)
ADS2 =count(split,c("group"))
tots <- sum(ADS2$freq)
ADS2$Norm <- tots
ADSN <- transform(ADS2, normalized_freq = freq / tots)
ADSN<-subset(ADSN, group !="NA")
ADSZ<-ADSN
ADSZ$freq <- 0
ADSZ$Norm <- 0
ADSZ$normalized_freq <- 0

ADSN$group<-as.numeric(ADSN$group)
p<-ggplot(ADSN, aes(x=group, y=(normalized_freq*100), xaxt = "n")) +
  geom_bar(stat ="identity",color = "black", fill = "grey")+
  theme_classic()+
  ggtitle("A") +
  ylab(label="Percentage of genes in each group N=19,889") +
  xlab(label="") +
  labs(fill = "")
p


split3$group<- as.factor(split3$group)
ADS3 =count(split3,c("group"))
tots <- sum(ADS3$freq)
ADS3$Norm <- tots
ADSO <- transform(ADS3, normalized_freq = freq / tots)
ADS2N <-merge(x = ADSZ, y = ADSO, by = c("group"), all.x = TRUE)
ADS2N<-subset(ADS2N, group !="NA")
ADS2N[is.na(ADS2N)] <- 0

q<-ggplot(ADS2N, aes(x=group, y=(normalized_freq.y*100), xaxt = "n")) +
  geom_bar(stat ="identity",color = "black", fill = "#0072B2")+
  theme_classic()+
  ggtitle("B") +
  ylab(label="Percentage of genes in each group (N=32)") +
  xlab(label="Haplosufficiency decile, 1=Haploinsufficient") +
  labs(fill = "")
q

tests  <- c( 11, 1977)
totals <- c( 32, 19857)
prop.test(tests, totals)


###FOUR###

rm(list=ls())
library(plyr)
library(ggfortify)
library(gridExtra)
library(grid)

split <- read.csv("/Users/mqbpqam6/Desktop/Ogunkanbi_Complete/July_19_Paper_data_run/Last_Hurrah_Dataset.csv", header = T)[c(1,6)]		#Insert path to DS1
split <- na.omit(split) #Remove any rows with null assignments from the analysis
split$type<-"Genome"
split3 <- read.csv("/Users/mqbpqam6/Desktop/Ogunkanbi_Complete/July_19_Paper_data_run/TGF_For_R.csv", header = T)[c(1,6)]		#Insert path to DS1
split3 <- na.omit(split3) #Remove any rows with null assignments from the analysis
split3$type<-"TGFβ"
split<-split[!(split$Gene_id %in% split3$Gene_id),] #Drop TGFB genes

ADS3 <- rbind(split3, split)

library(tidyverse)
ADS3<-ADS3 %>%
  mutate(PPIdecile = ntile(PPI_connectivity, 10))

s1 <-subset(ADS3, type=="Genome")
s2 <-subset(ADS3, type=="TGFβ")

detach("package:tidyverse", unload=TRUE)
detach("package:dplyr", unload=TRUE)
library(plyr)
ADS2 =count(s1,c("PPIdecile", "type"))
tots <- sum(ADS2$freq)
ADS2$Norm <- tots
ADSN <- transform(ADS2, normalized_freq = freq / tots)
ADSN$PPIdecile <- factor(ADSN$PPIdecile)
ADSN$PPIdecile <- factor(ADSN$PPIdecile, levels=rev(levels(ADSN$PPIdecile)))

ADSZ<-ADSN
ADSZ$freq <- 0
ADSZ$Norm <- 0
ADSZ$normalized_freq <- 0
ADSZ$type<-"TGFβ"

ADSb =count(s2,c("PPIdecile", "type"))
tots <- sum(ADSb$freq)
ADSb$Norm <- tots
ADSa <- transform(ADSb, normalized_freq = freq / tots)
ADSa$PPIdecile <- factor(ADSa$PPIdecile)
ADSa$PPIdecile <- factor(ADSa$PPIdecile, levels=rev(levels(ADSa$PPIdecile)))

New4 <-merge(x = ADSZ, y = ADSa, by = c("PPIdecile"), all.x = TRUE)
New4 <-New4[c(1,2,7,8,9)]
colnames(New4) <- c("PPIdecile", "type", "freq", "Norm", "normalized_freq")
New4$Norm[is.na(New4$Norm)] <- 32
New4[is.na(New4)] <- 0

ADSF <- New4

ADSF$type<-as.factor(ADSF$type)
levels(ADSF$type)[levels(ADSF$type)=="TGFβ"]<-"TGFb N=32"


group.colours <- c("TGFb N=32" = "#E69F00")

b <-ggplot(ADSF, aes(x = PPIdecile, y = (normalized_freq*100), fill=type)) + 
  geom_bar(stat="identity",position = position_dodge()) +
  theme_classic()+
  ggtitle("") +
  ylab(label="Percentage of total genes in group") +
  xlab(label="PPI Decile") +
  scale_fill_manual(values=group.colours)+
  labs(fill = "")
b


#8
tests  <- c( 12, 2260)
totals <- c( 32, 22690)
prop.test(tests, totals)

smokers  <- c( ADSF$freq )
patients <- c( ADSF$Norm )
pairwise.prop.test(smokers, patients, p.adjust.method = "bonferroni")

