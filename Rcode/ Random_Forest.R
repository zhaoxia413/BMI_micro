library(data.table)
library(tidyverse)
library(randomForest)
library(plotROC)
library(ggsci)
library(ggthemes)
library(ggpubr)
options(stringsAsFactors = F)
expr<-fread("../database/ExMicroSy/TuMicroSy/allOTUtable_split.csv",data.table = F)
data<-fread("../database/ExMicroSy/TuMicroSy/combined_metadata.csv")[,1:10]
data<-data[!duplicated(data$sampleID),]
meta<-subset(data,body_site%in%c('stool',"oral"))
expr[1:5,1:9]
meta[1:5,1:10]
meta<-meta[which(meta$sampleID%in%colnames(expr)),]
phylum<-expr[,-c(1:2,4:9)]
genus<-expr[,-c(1:6,8:9)]
genus<-genus%>%group_by(Genus)%>%summarise_all(sum)%>%as.data.frame()
phylum<-phylum%>%group_by(Phylum)%>%summarise_all(sum)%>%as.data.frame()
genus<-genus[-1,]
phylum<-phylum[-1,]
genus<-data.frame(row.names = genus$Genus,genus[,-1])
phylum<-data.frame(row.names = phylum$Phylum,phylum[,-1])
genus[1:4,1:4]
phylum[1:4,1:4]
tg<-data.frame(t(genus))
tp<-data.frame(t(phylum))
tg$sampleID<-rownames(tg)
tp$sampleID<-rownames(tp)
tg[1:3,1:3]
tp[1:3,1:3]
meta[1:3,1:10]
pg<-merge(tp,tg,by="sampleID")
rm(genus,phylum,tp,tg)
which(is.na(meta$BMI))
microGP<-merge(meta,pg,by="sampleID")
microGP<-microGP[-which(is.na(microGP$BMI)),]
range(microGP$BMI)
microGP$BMIgroup<-microGP$BMI
microGP$BMIgroup<-ifelse(microGP$BMIgroup>=30,"BMI_high",ifelse(microGP$BMIgroup<=25,"BMI_high","Normal"))
longdata <- melt_roc(microGP, "BMIgroup", colnames(microGP)[c(13,14)])
colnames(microGP)[1:11]
ggplot(microGP, aes(Disease,BMI, fill = BMIgroup)) +
  geom_boxplot()+
  scale_fill_d3()+
  theme_classic2()
#commparisons
my_comparisons<-list(c("control", "adenoma"), 
                     c("control", "AS"), 
                     c("control", "cephalosporins"),
                     c("control", "cirrhosis"),
                     c("control", "CRC"),
                     c("control", "IBD"),
                     c("control", "T1D"),
                     c("control", "T2D"))
fun_to_plot <- function(data, group, variable) {
  p <- ggboxplot(data, x=group, y=variable,fill = group, 
                add = "jitter")+
    stat_compare_means(comparisons = my_comparisons)+
    theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))+
    guides(fill=FALSE)
    #stat_compare_means(label.y = 125)
  return(p)
}
fun_to_plot(microGP,"Disease","BMI")
CRC<-subset(microGP,Disease=="CRC")
##microbiome abundance
source("../../TuMicroAnalysor/R/OTUtable_merge.R")
OTUmerge<-OTUtablemerge(OTUtable_split = expr)
source("../TuMicroAnalysor/R/reads2abundance.R")
group<-data.frame(samples<-meta$sampleID,Group=meta$dataset_name)
micro<-reads2abundance(data,meta = group)
my_comparisons<-list(c("BMI_high", "Normal"), 
                     c("BMI_high", "Normal"), 
                     c("BMI_high", "BMI_low"),
                     c("BMI_low", "Normal"))
fun_to_plot(microGP,"BMIgroup","BMIgroup")