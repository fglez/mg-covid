###MAKING FILES THAT CAN BE ANALYZED ON R

#The following script shows how to transform to phyloseq objects, 
#the tables produced by the abundance.sh script
#The abundace.sh script can be found here https://github.com/carpentries-incubator/metagenomics/blob/gh-pages/files/abundance.sh

#This script is a modified version from the one found here https://carpentries-incubator.github.io/metagenomics/08-automating_abundance/index.html
#by Nelly Selem, Claudia Zirion and Diego Garfías
#if you wish to learn more, please visit the link above

#This script also shows how to merge the phyloseq object produced by different samples, into one.

#Finally it helps you to prepare your phyloseq object for the following analyses by leaving only bacterial reads, 
#samples with a similar read depth and normalized data

#Libraries, install them if needed
#install.packages("readr")
library("readr")
#install.packages("BiocManager")
library("BiocManager")
#BiocManager::install("phyloseq")
library("phyloseq")
library("patchwork")
#install.packages("vegan")
library("vegan")



#READING SAMPLE FILES AND TRANSFORMING THEM TO PHYLOSEQ OBJECTS

#Reading kraken bracken files
pth <- "/home/mfcg/Descargas/covid/microbiome/taxonomia/taxanalysis2224" #your path to your files
setwd(pth)

#Making the numeration in order to read the files, we have 32 files from the 32 samples
dece <- c(
  c(rep(0, 9)),
  c(rep(1, 10)),
  c(rep(2, 10)),
  c(rep(3, 3))
          )
unid <- c(
  c(rep(c(1:9, 0), 3)),
  c(1, 2)
  )
samnum <- paste(dece, unid, sep = "") #sample numbers

#Transforming from the tables of the abundance.sh script to phyloseq objects
allmtgnm <- list() #Empty list to save the results of the ab_to_ps function
#The following function automates the steps to transform from the tables done by abundace.sh to a phyloseq object
ab_to_ps <- function(smnm){
  #File names
  spth <- paste0("/home/mfcg/Descargas/covid/microbiome/taxonomia/krakenbraken/AC1ME2SS", smnm, "/results/AC1ME2SS", smnm, "_kraken.") #Path to sample files
  rnkwc <- paste0(spth,"ranked-wc") #Sample ranked-wc file name
  ltbl <- paste0(spth, "lineage_table-wc") #Sample lineage_table-wc file name
  #Reading abundance.sh files
  otus <- read_delim(rnkwc,"\t", escape_double = FALSE, trim_ws = TRUE)
  taxa <- read_delim(ltbl, "\t", escape_double = FALSE,
                     col_types = cols(subspecies = col_character(), subspecies_2 = col_character()), 
                     trim_ws = TRUE)
  #From abundance.sh files to matrixes
  abnd <- as.matrix(otus[ , -1]) #To avoid that the OTU column is taken as a sample the first column is omitted
  lngs <- as.matrix(taxa)
  row.names(abnd) <- otus$OTU #Their identity is not lost cause their names are saved as rownames before
  row.names(lngs) <- taxa$OTU
  #From matrixes to phyloseq objects
  otab <- otu_table(abnd, taxa_are_rows = TRUE)
  ttab <- tax_table(lngs)
  psmet <- phyloseq(otab, ttab) #Can be left as it is if you want untrimmed data
  psmcn <- prune_taxa(taxa_sums(psmet)>100, psmet) #Filter OTUs with less than 100 reads, I'm using this number considering the error rate and length of the sequences, avoid if you want untrimmed data
  return(psmcn)
}

#Running the function ab_to_ps for all the samples and saving the phyloseq object of all the metagenomes in a list 
for(i in 1:length(samnum)){
  ipsi <- ab_to_ps(samnum[i])
  allmtgnm <<- merge_phyloseq(allmtgnm, ipsi)
}

#Save phyloseq object of all the samples with all the taxa present
write_rds(allmtgnm, file="allmtgnm.rds")

#Read file with all the samples and taxa
allmtgnm <- read_rds("allmtgnm.rds")

#Checking how many reads each sample has
Identified.Reads <- c(sample_sums(allmtgnm))
plot(sample_sums(allmtgnm), #Visualization reads per sample
     log = "y", 
     ylim = c(1000, 100000000),
     ylab = "Number of reads",
     xlab = "Samples")

#CCLEANING DATA IN THE PHYLOSEQ OBJECT
#Removing data with less than 1000000 reads, sample 13 is removed
clnmtgnm <- prune_samples(sample_sums(allmtgnm)>=1000000, allmtgnm) 
sample_sums(clnmtgnm)
plot(sample_sums(clnmtgnm), 
     log = "y", 
     ylim = c(1000, 100000000),
     ylab = "Number of reads",
     xlab = "Samples")
#Leaving only bacterial taxa
clnmtgnm <- subset_taxa(clnmtgnm, superkingdom == "Bacteria") 
#Saving phyloseq object with clean data in a rds
write_rds(clnmtgnm, file = "clnmtgnm.rds") 
#Read the rds file object with the clean data
clnmtgnm <- read_rds("clnmtgnm.rds") 

smpdt <- read.csv("mtdtcvd4.csv", header = TRUE, row.names = 1) #Read a csv file containing a table with the sample data or metadata
smpdt <- cbind(smpdt, Identified.Reads)
write.csv(smpdt, file = "mtdtcvd5.csv")
smpdt <- smpdt[-13, ] #Remove sample 13
smpdt <- sample_data(smpdt) #Giving the csv metadata file the phyloseq format
row.names(smpdt) #Sample names in the sample data


#ADDING METADATA
#Retrieve as tables -otu table and taxa table- the components of the phyloseq object with clean data -clnmtgnm-
#Taxa table
txtbl <- tax_table(clnmtgnm)
#Otu table
ottbl <- otu_table(clnmtgnm)
colnames(ottbl) #Names of the samples in the otu table
#Giving the same sample names to ottbl and smpdt
sort(row.names((smpdt)), decreasing = FALSE) #Sorting sample names so they are in the same order as phyloseq samples
colnames(ottbl) <- sort(rownames(smpdt), decreasing = FALSE) #Changing otu table names for the ones in the sample table

#Phyloseq object with the sample slot updated
mdnmgnm <- phyloseq(txtbl, ottbl, smpdt) #Putting together otu, taxa and sample tables
write_rds(mdnmgnm, file = "mdnmgnm.rds")

mdnmgnm <- read_rds("mdnmgnm.rds") #Saving metadata and metagenome as an rds file
mdnmgnm@sam_data #Verifying the sample data was added to the file

