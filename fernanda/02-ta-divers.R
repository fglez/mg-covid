##DIVERSITY ANALYSIS

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

#Plot of observed OTUs grouped by symptoms and color by stage
obs_symsta <- ggplot() + 
  geom_point(data = mdalp,
          aes(x = Symptoms.Stage, y = Observed, color = Symptoms.Stage)) +
  scale_fill_viridis_b() +
  theme_minimal()

obs_symsta


#Number of reads and recovered MAGs
maremeta <- read.csv("covid_isam_table.csv")
rea_symsta <- ggplot() + 
  geom_point(data = maremeta,
             aes(x = Symptoms.Stage, y = Identified.Reads, color = Symptoms.Stage)) +
  scale_fill_viridis_b() +
  theme_minimal()
rea_symsta

mag_symsta <- ggplot() + 
  geom_point(data = maremeta,
             aes(x = Symptoms.Stage, y = MAGs, color = Symptoms.Stage)) +
  scale_fill_viridis_b() +
  theme_minimal()
mag_symsta

#GENERA 

#Using microbiome function aggregate_rare to aggregate the OTUs by taxa 
#Keeping the phyloseq object
gnrr <- aggregate_rare(x = mdnmgnm, #Phyloseq object with the 31 metagenome samples
                       level = "genus", #Indication that they aggregate at genera level
                       detection = 1, #Treshold of absence presence
                       prevalence = 1/100, #Prevalence across the samples
                       include.lowest = TRUE) #Include the lower boundary of detection 1/100

View(gnrr@otu_table)
View(gnrr@tax_table)

#HEATMAP CORE TAXA

tpngnr <- names(sort(taxa_sums(gnrr), TRUE)[1:11]) #Top ten genera with more reads in the sample set
tpngnr <- tpngnr[-which(tpngnr == "Unknown")] #Removing unknown genera
tppgnr <- prune_taxa(tpngnr, gnrr) #Leaving only the top ten taxa in a phyloseq object
tplgmt <- log10(as.matrix(otu_table(tppgnr, taxa_are_rows))) #Log10 matix with the abundances of the most abundant genera in the sample set

tplglf <- melt(tplgmt) #Log10 matix with the abundances of the most abundant genera in the sample set in a long format
colnames(tplglf) <- c("Genus", "Sample", "Log10")

ggplot(data = tplglf,
       mapping = aes(Sample, Genus, fill=Log10)) +
  geom_tile() +
  scale_fill_viridis_b() +
  theme_minimal()

smpsmp <- as.numeric(as.factor(gnrr@sam_data$Symptoms.Stage))
snscol <- brewer.pal(6, "Set1")[smpsmp]
heatmap(tplgmt, ColSideColors = snscol)

#SPECIES 

#Using microbiome function aggregate_rare to aggregate the OTUs by taxa 
#Keeping the phyloseq object
spcsr <- aggregate_rare(x = mdnmgnm, #Phyloseq object with the 31 metagenome samples
                       level = "species", #Indication that they aggregate at genera level
                       detection = 1, #Treshold of absence presence
                       prevalence = 1/100, #Prevalence across the samples
                       include.lowest = TRUE) #Include the lower boundary of detection 1/50


#HEATMAP CORE TAXA

tpnspc <- names(sort(taxa_sums(spcsr), TRUE)[1:11]) #Top ten genera with more reads in the sample set
tpnspc <- tpnspc[-which(tpnspc == "Unknown")] #Removing unknown genera
tppspc <- prune_taxa(tpnspc, spcsr) #Leaving only the top ten taxa in a phyloseq object
tplsmt <- log10(as.matrix(otu_table(tppspc, taxa_are_rows))) #Log10 matrix with the abundances of the most abundant genera in the sample set

lfttmsp <- melt(tplsmt) #Log10 matix with the abundances of the most abundant genera in the sample set in a long format
colnames(lfttmsp) <- c("Species", "Sample", "Log10")

ggplot(data = lfttmsp,
       mapping = aes(Sample, Species, fill=Log10)) +
  geom_tile() +
  scale_fill_viridis_b() +
  theme_minimal()

ssmpsmp <- as.numeric(as.factor(spcsr@sam_data$Symptoms))
ssnscol <- brewer.pal(3, "Set1")[ssmpsmp]
heatmap(tplsmt, scale = "column", ColSideColors = ssnscol)



#GENERA OPTION 2

#GENERA 

#Using microbiome function aggregate_rare to aggregate the OTUs by taxa 
#Keeping the phyloseq object
gnrr <- aggregate_rare(x = mdnmgnm, #Phyloseq object with the 31 metagenome samples
                       level = "genus", #Indication that they aggregate at genera level
                       detection = 1, #Treshold of absence presence
                       prevalence = 1/100, #Prevalence across the samples
                       include.lowest = TRUE) #Include the lower boundary of detection 1/50

#HEATMAP CORE TAXA

tpngnr <- names(sort(taxa_sums(gnrr), TRUE)[1:23]) #Top ten genera with more reads in the sample set
tpngnr <- tpngnr[-which(tpngnr == "Unknown")] #Removing unknown genera
tppgnr <- prune_taxa(tpngnr, gnrr) #Leaving only the top ten taxa in a phyloseq object
tplgmt <- t(log10(as.matrix(otu_table(tppgnr, taxa_are_rows)))) #Log10 matrix with the abundances of the most abundant genera in the sample set

symsta <- as.factor(gnrr@sam_data$Symptoms.Stage) #Symptoms and stages unique values as factors
ssfac <- as.numeric(as.factor(gnrr@sam_data$Symptoms.Stage)) #Recovering as numeric factors symptoms and stage
tgnssf <- cbind(tplgmt, ssfac) #Adding symptoms and stage to the matrix
tgnssf <- tgnssf[order(tgnssf[ , "ssfac"], decreasing = FALSE), ] #Re arrange matrix so the samples with the same symptoms and stage togertbe
tgnssm <- tgnssm[ , -23] #Removing the metadata information column

tgndf <- melt(tgnssm) #Log10 matrix of top genera as dataframe
colnames(tgndf) <- c("Sample", "Genus", "Log10") #Column names for the dataframe

ggplot(data = tgndf,
       mapping = aes(Sample, Genus, fill=Log10)) +
  geom_tile() +
  scale_fill_viridis_b() +
  theme_minimal()

