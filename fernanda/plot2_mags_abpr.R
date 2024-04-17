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

#Read a csv file that recopiles all the information on the mags recovered
#that meet the requirements of >80% of completion and <5% of contamination
#this file contains the sample they came from, species, completion, contamination, L50, N50 and virulence genes
magcsv <- read.csv("goodmaglist.csv") #MAGs, taxonomy, their sampe and other measures

#Read a csv file with the metadata of the samples
maremeta <- read.csv("covid_isam_table.csv") #Metadata wih the number of reads and recovered MAGs
maremeta <- maremeta[-13, ]



#HEATMAP MAG GENERA 

#Only the metadata of the samples that have MAGs of good quality
magmeta <- maremeta[which(maremeta$MAGs > 0), ]

#Removing samples without good quality MAGs recovered, taxonomic information and other measures
jstmags <- magcsv[which(is.na(magcsv$MAG_GENERA) == FALSE), ]
gnrmags <- jstmags[ , c("SAMPLES", "MAG_GENERA")]
gmagmat <- as.matrix(table(gnrmags)) #To matrix to count how many MAGs of a genera were recovered from each sample
gmagmat <- cbind(gmagmat, SAMPLES = magmeta$SAMPLES, SYMPTOMS.STAGE = magmeta$Symptoms.Stage) #Adding Symptom and Stage to group them, as well as the sample to demonstrate each category is with its correspongin sample
gmagmat <- gmagmat[order(as.numeric(as.factor(gmagmat[ ,"SYMPTOMS.STAGE"]))), ] #Order the samples by their symptoms presence and stage of infection
gemagmt <- gmagmat[ , !colnames(gmagmat) %in% c("SAMPLES", "SYMPTOMS.STAGE")] #Only genera and samples

condtab <- as.data.frame(gmagmat[ , colnames(gmagmat) %in% c("SAMPLES", "SYMPTOMS.STAGE")])
colnames(condtab) <- c("SAMPLES", "SYMPTOMS.STAGE")
condtab$SAMPLES <- factor(condtab$SAMPLES, levels = unique(condtab$SAMPLES))
condtab$SYMPTOMS.STAGE <- as.factor(condtab$SYMPTOMS.STAGE)

gmaglng <- melt(gemagmt) #To long format
colnames(gmaglng) <- c("SAMPLES", "GENERA", "MAGS")

#Heatmap showing only the samples that have MAGs, the genera and the quantity of them recovered
qmagen_hm <- ggplot(
  data = gmaglng,
  mapping = aes(x = SAMPLES, y = GENERA, fill = MAGS)) +
  geom_tile() +
  scale_fill_viridis_d(begin = 0.35, end = 0.65, option = "B", direction = 1) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "top")

#Plot showing if the sample with MAGs is in what stage of infection and the clinical manifestations
mgsm_syst <- ggplot(data = condtab,
       mapping = aes(x = SAMPLES, y = NA, fill = SYMPTOMS.STAGE)) +
  geom_tile() +
  scale_fill_viridis_d(option = "E") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(color = "white"),
        axis.title.y = element_text(color = "white"),
        legend.position = "bottom",
        legend.title = element_blank())
  
pl2_main <- ggarrange(qmagen_hm, mgsm_syst, ncol = 1, nrow = 2, heights = c(9, 1), align = "v")




#Presence of recovered MAGs genera in reads not MAGs

#How to make this file can be read on the 01-ta-format.R script
mdnmgnm <- read_rds("mdnmgnm.rds") #Reading rds file with phyloseq object that contains metagenomes and metadata

#Using microbiome function aggregate_rare to aggregate the OTUs by taxa 
#Keeping the phyloseq object
gnrr <- aggregate_rare(x = mdnmgnm, #Phyloseq object with the 31 metagenome samples
                       level = "genus", #Indication that they aggregate at genera level
                       detection = 1, #Treshold of absence presence
                       prevalence = 1/100, #Prevalence across the samples
                       include.lowest = TRUE) #Include the lower boundary of detection 1/100

mgenera <- unique(jstmags$MAG_GENERA) #MAGs genera


mgotutb <- log10(as.matrix(otu_table(gnrr, taxa_are_rows))) #Data frame with samples and number of reads for each genera foun with Kraken Braken 2
mgotutb <- mgotutb[
  mgenera[which(mgenera %in% rownames(mgotutb) == TRUE)], #Which genera were identfied in the reads by Kraken Braken 2 and also were recovered as MAGs
]
mgotutb <- t(mgotutb) #Change orientation
mgotutb <- cbind(mgotutb, SAMPLES = maremeta$SAMPLES, SYMPTOMS.STAGE = maremeta$Symptoms.Stage)
mgotutb <- mgotutb[order(as.numeric(as.factor(mgotutb[ , "SYMPTOMS.STAGE"]))), ]
magremt <- mgotutb[ , !colnames(mgotutb) %in% c("SAMPLES", "SYMPTOMS.STAGE")] #Remove their stage and symptoms, leave only number of reads per genera

magrelg <- melt(magremt) #Number of reads identified of the MAG genera
colnames(magrelg) <- c("SAMPLES", "GENERA", "READS") 
class(magrelg$READS) <- "numeric" #Numeric the reads column

allcon <- as.data.frame(mgotutb[ , colnames(mgotutb) %in% c("SAMPLES", "SYMPTOMS.STAGE")])
colnames(allcon) <- c("SAMPLES", "SYMPTOMS.STAGE")
allcon$SAMPLES <- factor(allcon$SAMPLES, levels = unique(allcon$SAMPLES))
allcon$SYMPTOMS.STAGE <- factor(allcon$SYMPTOMS.STAGE, levels = unique(allcon$SYMPTOMS.STAGE))

#Heatmap showing all the samples with the number of reads of the MAGs genera
mgenre_hm <- ggplot(
  data = magrelg,
  mapping = aes(x = SAMPLES, y = GENERA, fill = READS)) +
  geom_tile() +
  scale_fill_viridis_c(begin = 0.35, end = 0.65, option = "B", direction = 1) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "top")

#Plot showing if the sample with MAGs is in what stage of infection and the clinical manifestations
allcon_p <- ggplot(data = allcon,
                    mapping = aes(x = SAMPLES, y = NA, fill = SYMPTOMS.STAGE)) +
  geom_tile() +
  scale_fill_viridis_d() +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(color = "white"),
        axis.title.y = element_text(color = "white"),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.box = "vertical") +
  guides(col = guide_legend(ncol = 2))

pl2_ptb <- ggarrange(mgenre_hm, allcon_p, ncol = 1, nrow = 2, heights = c(6, 1), align = "v")


#Boxplot number of virulence genes per group (Stage of infection and presencee/absence of Symptoms)

magmeta #Using previous table with only the metadata of the samples that have MAGs of good quality

magcsv #Using matrix with MAGs information<
nonamag <- magcsv[which(is.na(magcsv$MAG_GENERA) == FALSE), ] #Only information about MAGs removing the samples that doesnt have MAGs
nonamag <- inner_join(nonamag, magmeta, by = "SAMPLES") #Adding metadata to the MAG information matrix

pl2_ptc <- ggplot(
  data = nonamag,
  mapping = aes(x = Symptoms.Stage, y = VIRULENCE_GENES, fill = Symptoms.Stage)) +
  geom_boxplot() +
  scale_fill_viridis_d() +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(x = "SYMPTOMS AND STAGE",
       y = "#VIRULENCE GENES")

pl2_extra <- ggarrange(pl2_ptb, pl2_ptc, ncol = 1, nrow = 2, heights = c(4, 2))

ggarrange(pl2_main, pl2_extra, ncol = 2, nrow = 1, widths = c(5, 4))

