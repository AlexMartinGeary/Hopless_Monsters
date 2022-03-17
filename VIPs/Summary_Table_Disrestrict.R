##Summary Table for Paper

rm(list=ls())
library(ggplot2)
library(plyr)
library(grid)
library(VennDiagram)
setwd("/Users/mqbpqam6/Desktop/NEW_VR_Data/February_2020/February_rerun/Statistical_analysis_and_interpretation")

d1<- read.csv("Masterfile_14Feb.csv", header = T, strip.white=TRUE, stringsAsFactors=FALSE)

d1 <- subset(d1, Disease_Association!="NA") #restrict to genes in VR and LH
levels(d1$Disease_Association)[levels(d1$Disease_Association)=="Both"]<-"Dominant"

d1$VIP_type<-as.factor(d1$VIP_type)
levels(d1$VIP_type)[levels(d1$VIP_type)=="Multi"]<-"Both"
levels(d1$VIP_type)[levels(d1$VIP_type)=="NA"]<-"None"

d1$Ortho_Age<-as.factor(d1$Ortho_Age)
levels(d1$Ortho_Age)[levels(d1$Ortho_Age)=="796"]<-"WGD & Earlier"
levels(d1$Ortho_Age)[levels(d1$Ortho_Age)=="615"]<-"WGD & Earlier"
levels(d1$Ortho_Age)[levels(d1$Ortho_Age)=="435"]<-"WGD & Earlier"
levels(d1$Ortho_Age)[levels(d1$Ortho_Age)!="WGD & Earlier"]<-"Since WGD"


node.Disease =count(d1,c("control.category", "Disease_Association"))
node.VIP =count(d1,c("control.category", "VIP_type"))
node.PS =count(d1,c("control.category", "Paralog_status"))
node.Age =count(d1,c("control.category", "Ortho_Age"))
disease.VIP =count(d1,c("Disease_Association", "VIP_type"))
disease.PS =count(d1,c("Disease_Association", "Paralog_status"))
disease.Age =count(d1,c("Disease_Association", "Ortho_Age"))
VIP.PS =count(d1,c("VIP_type", "Paralog_status"))
VIP.Age =count(d1,c("VIP_type", "Ortho_Age"))
PS.Age =count(d1,c("Paralog_status", "Ortho_Age"))

nodes =count(d1,c("control.category"))
disease =count(d1,c("Disease_Association"))
VIPs =count(d1,c("VIP_type"))
