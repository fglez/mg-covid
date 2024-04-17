
#The following script shows an analysis of the microbiome samples from mexican people
#during the COVID pandemic

#Required libraries
#Libraries, install them if needed
#install.packages("ggplot2")
library("ggplot2")
#install.packages("readr")
library("readr")
#install.packages("BiocManager")
#libary("BiocManager")
#BiocManager::install("phyloseq")
library("phyloseq")
library("patchwork")
#install.packages("vegan")
library("vegan")
#install.packages("tidyverse")
library("tidyverse")
#library(BiocManager)
#library(devtools)
#BiocManager::install("microbiome")
library("microbiome")  
#install.packages("reshape2")
library("reshape2")
#install.packages("RColorBrewer")
library("RColorBrewer")
#install.packages("ggpubr")
library("ggpubr")


#READING FILES

#Set path to your files
pth <- "/home/mfcg/Descargas/covid/microbiome/taxonomia/taxanalysis2224" #your path to your files
setwd(pth)
#How to make this file can be read on the 01-ta-format.R script
mdnmgnm <- read_rds("mdnmgnm.rds") #Reading rds file with phyloseq object that contains metagenomes and metadata


#ALPHA DIVERSITY
alph <- estimate_richness(mdnmgnm, measure = "Observed") #Estimation of observed OTUs
mtdt <- as.data.frame(mdnmgnm@sam_data) #Retrieving metadata
mdalp <- cbind(mtdt, alph) #Joint data frame of metadata and observed OTUs


#READS AND MAGs 
#Table with the count of reads that was obtained in the script 01-ta-format.R saved as mtdtcvd5.csv
#this table was modified manually by counting the number of MAGs recovered
maremeta <- read.csv("covid_isam_table.csv") #Metadata wih the number of reads and recovered MAGs
maremeta <- maremeta [-13, ] #Removing sample with less than a 1000000 reads identified with krakenbracken
maremeta <- cbind(maremeta, alph) #Adding the observed OTUs
#Arrangin the table by symptoms and stage of COVID-19
sysafc <- as.factor(maremeta$Symptoms.Stage) #Changing symptoms and stage as factors
sysanf <- as.numeric(sysafc) #Changing factors to numeric
maremdfc <- cbind(maremeta, sysanf) #Adding numeric factors to metadata
maremdfc <- maremdfc[order(maremdfc[ , "sysanf"], decreasing = FALSE), ] #Re arrange matrix so the samples with the same symptoms and stage togertbe

#Observed OTUs in each sample ordered by symptoms and stage
#Scatter plot
obs_syst_p <- ggplot() + 
  geom_point(data = maremdfc, #Observed OTUs per sample, color by stage + symptoms
             aes(x = factor(X, levels = unique(X)), y = Observed, color = Symptoms.Stage)) +
  scale_fill_viridis_b() +
  theme_minimal() +
  theme(axis.text.x = element_blank(), #Removing labels to improve aesthetics
        axis.title.x = element_blank(),
        legend.position = "none") +
  labs(y = "Observed OTUs")
obs_syst_p
#Box plot
obso_bx <- ggplot() +
  geom_boxplot(data = maremdfc, #Observed OTUs per sample, color by stage + symptoms
               aes(x = Symptoms.Stage, y = Observed, color = Symptoms.Stage)) +
  scale_fill_viridis_b() +
  theme_minimal() +
  theme(axis.text.x = element_blank(), #Removing labels to improve aesthetics
        axis.title.x = element_blank(),
        legend.position = "none") +
  labs(y = "Observed OTUs")
obso_bx

#Number of Reads identified by Kracken Bracken in each sample
#Scatter plot showing the number of reads identified as bacterial by Kracken Bracken
ireads_p <- ggplot() + 
  geom_point(data = maremdfc, #Color show the symptom and stage when the sample was took
             aes(x = factor(X, levels = unique(X)), y = Identified.Reads, color = Symptoms.Stage)) +
  scale_fill_viridis_b() +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none") +
  labs(y = "Bacterial Reads")
ireads_p

#Number of MAGs recovered in each sample
#Scatter plot showing the number of MAGs Recovered by each sample
magsr_p <- ggplot() + 
  geom_point(data = maremdfc, #Color indicates the symptom and stage of the patient when the sample was took
             aes(x = factor(X, levels = unique(X)), y = MAGs, color = Symptoms.Stage)) +
  scale_fill_viridis_b() +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "bottom") +
  labs(y = "Recovered MAGs")
magsr_p
mags_bx <- ggplot() +
  geom_boxplot(data = maremdfc, #Color indicates the symptom and stage of the patient when the sample was took
               aes(x = Symptoms.Stage, y = MAGs, color = Symptoms.Stage)) +
  scale_fill_viridis_b() +
  theme_minimal() +
  theme(axis.text.x = element_blank(), #Removing labels to improve aesthetics
        axis.title.x = element_blank(),
        legend.position = "none") +
  labs(y = "Recovered MAGs")
mags_bx
  

scatterps <- ggarrange(obs_syst_p, ireads_p, magsr_p, nrow = 3, heights = c(6, 3, 7), align = "v")
scatterps

bxplts <- ggarrange(obso_bx, mags_bx, nrow = 2, align = "v")
bxplts

allplt <- ggarrange(scatterps, bxplts, ncol = 2, widths = c(3,1))
allplt


#Counts by stage and symptoms

table(maremeta$Symptoms.Stage)

#Obtaining means by Stage of infection and presence/absence of Symptoms

maremeta %>%
  group_by(Symptoms.Stage) %>%
  summarise_at(vars(Observed), mean)

#Number of MAGs per stage-symptom group

maremeta %>%
  group_by(Symptoms.Stage) %>%
  summarise_at(vars(MAGs), sum)

#Mean MAGs recovered by stage-symptom group
maremeta %>%
  group_by(Symptoms.Stage) %>%
  summarise_at(vars(MAGs), mean)

