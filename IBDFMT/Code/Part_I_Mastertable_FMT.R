####################################
# SCRIPT 1: CREATE A MASTERTABLE
####################################
# Before starting this script
####################################
# ANALYSIS: Transfer the "Transfer" folder to my own desktop
# ------------------> You could use this script on every single "Transfer" folder
# ------------------> Go to "0. packages" to download & load all necessary packages
# ------------------> Go to "2. Update taxonomy database" to update taxonomy database
# ------------------> Go to "10.1 List eukaryotic viruses" and update this list on occasion
####################################
# 0. Packages: install packages and load them wherever needed
####################################
#install.packages("ggplot2")
#install.packages("dplyr")
#install.packages("plyr")
#install.packages("reshape")
#install.packages("reshape2")
#install.packages("scales")
#install.packages("viridis")
#install.packages("readr")
#install.packages("WriteXLS")
#install.packages("readxl")
#install.packages("devtools")
#install_github("jokergoo/ComplexHeatmap")
#install.packages("dendsort")
#install.packages("vegan")
#install.packages("ape")
#install.packages("devtools")
#install.packages("seriation")
#install.packages("taxonomizr")
#install.packages("ggThemeAssist")
#install.packages("esquisse")
#installed.packages("modeldata")
#install.packages("tidyverse")
#devtools::install_github('sherrillmix/taxonomizr')

library(taxonomizr)
library(esquisse)
library(modeldata)
library(tidyverse)
library(ggplot2)
library(ggThemeAssist)
library(ggplot2)
library(tibble)
library(dplyr)
library(plyr)
library(reshape)
library(reshape2)
library(scales)
library(grid)
library(readr)
library(viridis)
library(vegan)
library(tidyverse)
library(WriteXLS)
library(readxl)
library(devtools)
library(ComplexHeatmap)
library(circlize)
library(dendextend)
library(dendsort)
library(seriation)
library(vegan)
library(ggplot2)
library(tibble)
library(dplyr)
library(ape)
library(devtools)
library(taxonomizr)
####################################
# 1. Update taxonomy database
####################################
setwd("/Volumes/Daan_2/Desktop/DB/")
getwd()
dir()

rm(nucl_wgs.accession2taxid.gz)
rm(nucl_gb.accession2taxid.gz)
rm(nodes.dmp)
rm(names.dmp)
rm(accessionTaxa.sql-journal)
rm(accessionTaxa.sql)

getAccession2taxid(baseUrl='https://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/')
system('TMPDIR=/Volumes/Daan_2/Desktop/DB R -e "taxonomizr::prepareDatabase(\'accessionTaxa.sql\')"')
# prepareDatabase('accessionTaxa.sql') # 82 GB
####################################
# 2. Create abundance table 
####################################
## A) set directory
setwd("/Users/daan/Desktop/Transfer/Abundances/")
getwd()
dir()

## B) read the abundance text files and store in a list
list <- list()
list_txt <- dir(pattern = "*.magnitudes", full.names = FALSE)
for (k in 1: length(list_txt)){
  list[[k]] <- read.delim(list_txt[k])
}

## C) create names file
names <- sub(pattern = "*.magnitudes", "", list_txt)

## D) evaluate lists
head(list[[1]]) ## perform this command to check if the first header is an 'X', needed to merge in next command
#view(names)

## E) merge lists in a table by matching headers
table <- merge(list[[1]], list[[2]], by="X", all = TRUE) 
for (i in 1:((length(list))-2)){
  table <- merge(table, list [[i+2]], by = "X", all = TRUE)
}

## F) create a dataframe, discard first row and make contigs the row names,
table <- data.frame(table[,-1], row.names = table[,1])

## G) change all colnames by names vector
for (j in 1:length(names)){
  colnames(table) <- c(names)
}

## H) convert NA to 0 and dots with underscores
table[is.na(table)] <- 0
rownames(table) <- gsub("\\.{1}", "_", rownames(table))

## I) Remove sample with no reads mapped (in this case 7 samples found)
table <- table[, colSums(table != 0) > 0]
ncol(table)
# 313 --> 306 samples
####################################
# 3. Create vector of contigs
####################################
setwd("/Users/daan/Desktop/Transfer/FastaFiles")
getwd()
dir()

non_redundant_contig_names <- read.table("NR_contigs_names.txt", quote = "\"", comment.char = "")
colnames(non_redundant_contig_names) <- "contigs_orginal_name.txt"
rownames(non_redundant_contig_names) <- non_redundant_contig_names$contigs_orginal_name
rownames(non_redundant_contig_names) <- gsub("\\.{1}", "_", rownames(non_redundant_contig_names))
contigs_project <- non_redundant_contig_names$contigs_orginal_name 

#view(non_redundant_contig_names) # 'view' is performed with tibble
#View(contigs_project)           # 'View' is performed with utils
####################################
# 4. Insert taxonomical annotation 
####################################
# 4.1 Diamond TaxIdpath
####################################
setwd("/Users/daan/Desktop/Transfer/Taxonomy/")
getwd()
dir()

AllSamples_classificationtaxonomy_Diamond <- as.data.frame(read_delim("AllSamples.scaffolds.diamond.tab", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE))
AllSamples_classificationtaxonomy_Diamond$X3 <- NULL
colnames(AllSamples_classificationtaxonomy_Diamond) <- c('contig', 'Diamond_single_TaxIdpath')
AllSamples_classificationtaxonomy_Diamond <- AllSamples_classificationtaxonomy_Diamond[AllSamples_classificationtaxonomy_Diamond$contig %in% contigs_project,]
rownames(AllSamples_classificationtaxonomy_Diamond) <- AllSamples_classificationtaxonomy_Diamond$contig
AllSamples_classificationtaxonomy_Diamond$contig <- NULL
AllSamples_classificationtaxonomy_Diamond <- AllSamples_classificationtaxonomy_Diamond[-2] 
rownames(AllSamples_classificationtaxonomy_Diamond) <- gsub("\\.{1}", "_", rownames(AllSamples_classificationtaxonomy_Diamond))

#view(AllSamples_classificationtaxonomy_Diamond)
####################################
# 4.2 Blastn TaxIdpath
####################################
AllSamples_classificationtaxonomy_Blastn <- as.data.frame(read_delim("AllSamples.scaffolds.blastn_nocomments.out.tab", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE))[-1,]
AllSamples_classificationtaxonomy_Blastn$X3 <- NULL
colnames(AllSamples_classificationtaxonomy_Blastn) <- c('contig', 'Blastn_single_TaxIdpath')
rownames(AllSamples_classificationtaxonomy_Blastn) <- AllSamples_classificationtaxonomy_Blastn$contig
AllSamples_classificationtaxonomy_Blastn$contig <- NULL
AllSamples_classificationtaxonomy_Blastn <- AllSamples_classificationtaxonomy_Blastn[-2] # ALREADY REMOVED  
rownames(AllSamples_classificationtaxonomy_Blastn) <- gsub("\\.{1}", "_", rownames(AllSamples_classificationtaxonomy_Blastn))

#view(AllSamples_classificationtaxonomy_Blastn)
####################################
# 4.3 CAT TaxIdPath
####################################
AllSamples_classificationtaxonomy_CAT <- as.data.frame(read_delim("AllSamples.scaffolds.CAT.tab", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE))
AllSamples_classificationtaxonomy_CAT$X3 <- NULL
colnames(AllSamples_classificationtaxonomy_CAT) <- c('contig', 'CAT_single_TaxIdpath')
rownames(AllSamples_classificationtaxonomy_CAT) <- AllSamples_classificationtaxonomy_CAT$contig
AllSamples_classificationtaxonomy_CAT$contig <- NULL
AllSamples_classificationtaxonomy_CAT <- AllSamples_classificationtaxonomy_CAT[-2] 
rownames(AllSamples_classificationtaxonomy_CAT) <- gsub("\\.{1}", "_", rownames(AllSamples_classificationtaxonomy_CAT))

#view(AllSamples_classificationtaxonomy_CAT)
####################################
# 4.4 Virsorter 
####################################
AllSamples_VirsorterCategories <- as.data.frame(read_delim("Virsorter2.csv", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE))
AllSamples_VirsorterCategories$X2 <- NULL
AllSamples_VirsorterCategories$X4 <- NULL
AllSamples_VirsorterCategories$X6 <- NULL
AllSamples_VirsorterCategories$X8 <- NULL
colnames(AllSamples_VirsorterCategories) <- c('contig','dsDNA', 'ssDNA', 'RNA', 'Virsorter')
rownames(AllSamples_VirsorterCategories) <- AllSamples_VirsorterCategories$contig
AllSamples_VirsorterCategories$contig <- NULL
AllSamples_VirsorterCategories <- AllSamples_VirsorterCategories[c(4,1,2,3)]
rownames(AllSamples_VirsorterCategories) <- gsub("\\.{1}", "_", rownames(AllSamples_VirsorterCategories))
#View(AllSamples_VirsorterCategories)
####################################
# 4.5 Nayfach Genera 
####################################
setwd("/Users/daan/Desktop/Transfer/Nayfach")
getwd()
dir()

AllSamples_viral_genera <- as.data.frame(read_delim("genus_clusters.txt", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE))
colnames(AllSamples_viral_genera) <- c('contig','Genus')
rownames(AllSamples_viral_genera) <- AllSamples_viral_genera$contig
AllSamples_viral_genera$contig <- NULL
rownames(AllSamples_viral_genera) <- gsub("\\.{1}", "_", rownames(AllSamples_viral_genera))
#View(AllSamples_viral_genera)
####################################
# 5. Cenotetaker2: functions
####################################
# 5.1 Number of hallmark proteins
####################################
setwd ("/Users/daan/Desktop/Transfer/Functions") 
getwd()
dir()

AllSamples_Cenotetaker2_hallmark_proteins <- as.data.frame(read_delim("summary_total_viral_proteins.tsv", "\t", escape_double = FALSE, trim_ws = TRUE, col_names = FALSE))
AllSamples_Cenotetaker2_hallmark_proteins <- AllSamples_Cenotetaker2_hallmark_proteins[,c(1,3)]
colnames(AllSamples_Cenotetaker2_hallmark_proteins) <- c('contig', 'Hallmark proteins')
rownames(AllSamples_Cenotetaker2_hallmark_proteins) <- AllSamples_Cenotetaker2_hallmark_proteins$contig
AllSamples_Cenotetaker2_hallmark_proteins$contig <- NULL
rownames(AllSamples_Cenotetaker2_hallmark_proteins) <- gsub("\\.{1}", "_", rownames(AllSamples_Cenotetaker2_hallmark_proteins))
AllSamples_Cenotetaker2_hallmark_proteins[is.na(AllSamples_Cenotetaker2_hallmark_proteins)] <- 0
AllSamples_Cenotetaker2_hallmark_proteins <- as.data.frame(AllSamples_Cenotetaker2_hallmark_proteins)
AllSamples_Cenotetaker2_hallmark_proteins <- slice(AllSamples_Cenotetaker2_hallmark_proteins, 1:(n() - 1))
####################################
# 5.2 Lysogenic contigs
####################################
AllSamples_Cenotetaker2_lysogenic_phages <- as.data.frame(read_delim("All_predicted_protein_lysogenic_nodes.txt", "\t", escape_double = FALSE, trim_ws = TRUE, col_names = FALSE))
colnames(AllSamples_Cenotetaker2_lysogenic_phages) <- c('contig', 'lysogenic cycle')
rownames(AllSamples_Cenotetaker2_lysogenic_phages) <- AllSamples_Cenotetaker2_lysogenic_phages$contig
AllSamples_Cenotetaker2_lysogenic_phages$contig <- NULL
rownames(AllSamples_Cenotetaker2_lysogenic_phages) <- gsub("\\.{1}", "_", rownames(AllSamples_Cenotetaker2_lysogenic_phages))
AllSamples_Cenotetaker2_lysogenic_phages[is.na(AllSamples_Cenotetaker2_lysogenic_phages)] <- 0
#view(AllSamples_Cenotetaker2_lysogenic_phages)
###################################
# 5.3 Toxin-containing contigs
####################################
# re-evaluate if I don't miss any "words" for toxins as different names in orginal file
AllSamples_Cenotetaker2_toxin_phages <- as.data.frame(read_delim("All_predicted_toxin_phages_nodes.txt", "\t", escape_double = TRUE, trim_ws = TRUE, col_names = FALSE))
colnames(AllSamples_Cenotetaker2_toxin_phages) <- c('contig', 'Toxin presence', 'Toxin type')
rownames(AllSamples_Cenotetaker2_toxin_phages) <- AllSamples_Cenotetaker2_toxin_phages$contig
AllSamples_Cenotetaker2_toxin_phages$contig <- NULL
rownames(AllSamples_Cenotetaker2_toxin_phages) <- gsub("\\.{1}", "_", rownames(AllSamples_Cenotetaker2_toxin_phages))
AllSamples_Cenotetaker2_toxin_phages[is.na(AllSamples_Cenotetaker2_toxin_phages)] <- 0
#View(AllSamples_Cenotetaker2_toxin_phages)
####################################
# 5.4 Recombination-containing contigs
###################################
AllSamples_Cenotetaker2_recombination_phages <- as.data.frame(read_delim("All_predicted_recombination_phages_nodes.txt", "\t", escape_double = FALSE, trim_ws = TRUE, col_names = FALSE))
colnames(AllSamples_Cenotetaker2_recombination_phages) <- c('contig', 'recombination')
rownames(AllSamples_Cenotetaker2_recombination_phages) <- AllSamples_Cenotetaker2_recombination_phages$contig
AllSamples_Cenotetaker2_recombination_phages$contig <- NULL
rownames(AllSamples_Cenotetaker2_recombination_phages) <- gsub("\\.{1}", "_", rownames(AllSamples_Cenotetaker2_recombination_phages))
AllSamples_Cenotetaker2_recombination_phages[is.na(AllSamples_Cenotetaker2_recombination_phages)] <- 0
#view(AllSamples_Cenotetaker2_recombination_phages)
####################################
# 5.5 RT-containing contigs
####################################
AllSamples_Cenotetaker2_RT_phages <- as.data.frame(read_delim("RT_combined.txt", "\t", escape_double = FALSE, trim_ws = TRUE, col_names = FALSE))
colnames(AllSamples_Cenotetaker2_RT_phages) <- c('contig', 'RT', 'RT type')
rownames(AllSamples_Cenotetaker2_RT_phages) <- AllSamples_Cenotetaker2_RT_phages$contig
AllSamples_Cenotetaker2_RT_phages$contig <- NULL
rownames(AllSamples_Cenotetaker2_RT_phages) <- gsub("\\.{1}", "_", rownames(AllSamples_Cenotetaker2_RT_phages))
AllSamples_Cenotetaker2_RT_phages[is.na(AllSamples_Cenotetaker2_RT_phages)] <- 0
#view(AllSamples_Cenotetaker2_RT_phages)
####################################
# 5.5.6 Bacteriocin-containing contigs
####################################
AllSamples_Cenotetaker2_Bacteriocin_phages <- as.data.frame(read_delim("Bacteriocins.txt", "\t", escape_double = FALSE, trim_ws = TRUE, col_names = FALSE))
colnames(AllSamples_Cenotetaker2_Bacteriocin_phages) <- c('contig', 'bacteriocin')
rownames(AllSamples_Cenotetaker2_Bacteriocin_phages) <- AllSamples_Cenotetaker2_Bacteriocin_phages$contig
AllSamples_Cenotetaker2_Bacteriocin_phages$contig <- NULL
rownames(AllSamples_Cenotetaker2_Bacteriocin_phages) <- gsub("\\.{1}", "_", rownames(AllSamples_Cenotetaker2_Bacteriocin_phages))
AllSamples_Cenotetaker2_Bacteriocin_phages[is.na(AllSamples_Cenotetaker2_Bacteriocin_phages)] <- 0
#view(AllSamples_Cenotetaker2_Bacteriocin_phages)
####################################
# 5.7 Taxonomy Cenotetaker2
####################################
setwd("/Users/daan/Desktop/Transfer/Taxonomy/")
getwd()
dir()

AllSamples_classificationtaxonomy_Cenotaker2 <- as.data.frame(read_delim("CNoteTaker2_taxonomy2.txt", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE))[-1,]
AllSamples_classificationtaxonomy_Cenotaker2 <- AllSamples_classificationtaxonomy_Cenotaker2[,c(1,3)]
colnames(AllSamples_classificationtaxonomy_Cenotaker2) <- c('contig', 'Cenote_single_TaxIdpath')
rownames(AllSamples_classificationtaxonomy_Cenotaker2) <- AllSamples_classificationtaxonomy_Cenotaker2$contig
AllSamples_classificationtaxonomy_Cenotaker2$contig <- NULL
rownames(AllSamples_classificationtaxonomy_Cenotaker2) <- gsub("\\.{1}", "_", rownames(AllSamples_classificationtaxonomy_Cenotaker2))
AllSamples_classificationtaxonomy_Cenotaker2[is.na(AllSamples_classificationtaxonomy_Cenotaker2)] <- 0
#View(AllSamples_classificationtaxonomy_Cenotaker2)
####################################
# 6. Additional taxonomies:
####################################
# 6.1 Crass phage taxonomy
####################################
setwd("/Users/daan/Desktop/Transfer/CrAss_Ino_Micro/Crassvirales")
getwd()
dir()

AllSamples_classificationtaxonomy_Crass <- as.data.frame(read_delim("AllSamples.scaffolds.crass.blastn_reliable.nocomments.out", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE))
AllSamples_classificationtaxonomy_Crass <- AllSamples_classificationtaxonomy_Crass[,c(1,2)]
colnames(AllSamples_classificationtaxonomy_Crass) <- c('contig', 'CrAss_taxonomy')
rownames(AllSamples_classificationtaxonomy_Crass) <- AllSamples_classificationtaxonomy_Crass$contig
AllSamples_classificationtaxonomy_Crass$contig <- NULL
rownames(AllSamples_classificationtaxonomy_Crass) <- gsub("\\.{1}", "_", rownames(AllSamples_classificationtaxonomy_Crass))
#View(AllSamples_classificationtaxonomy_Crass)
####################################
# 6.2 Inoviridae phage taxonomy
####################################
setwd("/Users/daan/Desktop/Transfer/CrAss_Ino_Micro/Inoviridae")
getwd()
dir()

AllSamples_classificationtaxonomy_Ino <- as.data.frame(read_delim("AllSamples.scaffolds.ino.blastn_reliable.nocomments.out", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE))
AllSamples_classificationtaxonomy_Ino <- AllSamples_classificationtaxonomy_Ino[,c(1,2)]
colnames(AllSamples_classificationtaxonomy_Ino) <- c('contig', 'Ino_taxonomy')
rownames(AllSamples_classificationtaxonomy_Ino) <- AllSamples_classificationtaxonomy_Ino$contig
AllSamples_classificationtaxonomy_Ino$contig <- NULL
rownames(AllSamples_classificationtaxonomy_Ino) <- gsub("\\.{1}", "_", rownames(AllSamples_classificationtaxonomy_Ino))
#View(AllSamples_classificationtaxonomy_Ino)
####################################
# 6.3 Microviridae phage taxonomy
####################################
setwd("/Users/daan/Desktop/Transfer/CrAss_Ino_Micro/Microviridae")
getwd()
dir()

AllSamples_classificationtaxonomy_Micro <- as.data.frame(read_delim("AllSamples.scaffolds.micro.blastn_reliable.nocomments.out", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE))
AllSamples_classificationtaxonomy_Micro <- AllSamples_classificationtaxonomy_Micro[,c(1,2)]
colnames(AllSamples_classificationtaxonomy_Micro) <- c('contig', 'Micro_taxonomy')
rownames(AllSamples_classificationtaxonomy_Micro) <- AllSamples_classificationtaxonomy_Micro$contig
AllSamples_classificationtaxonomy_Micro$contig <- NULL
rownames(AllSamples_classificationtaxonomy_Micro) <- gsub("\\.{1}", "_", rownames(AllSamples_classificationtaxonomy_Micro))
#View(AllSamples_classificationtaxonomy_Micro)
####################################
# 7. Insert Completeness
####################################
setwd ("/Users/daan/Desktop/Transfer/Taxonomy") 
getwd() 
dir()

AllSamples_completeness_Checkv<- as.data.frame(read_delim("quality_summary_R.tsv", "\t", escape_double = FALSE, trim_ws = TRUE, col_names = FALSE))
colnames(AllSamples_completeness_Checkv) <- c('contig', 'total genes', 'viral genes', 'completeness_quality', 'completeness')
rownames(AllSamples_completeness_Checkv) <- AllSamples_completeness_Checkv$contig
AllSamples_completeness_Checkv$contig <- NULL
rownames(AllSamples_completeness_Checkv) <- gsub("\\.{1}", "_", rownames(AllSamples_completeness_Checkv))
AllSamples_completeness_Checkv[is.na(AllSamples_completeness_Checkv)] <- 0
# View(AllSamples_completeness_Checkv)
####################################
# 8. Host prediction
####################################
# 8.1 Host Phyla
####################################
setwd ("/Users/daan/Desktop/Transfer/Host") 
getwd() 
dir()

AllSamples_completeness_Host_phyla <- as.data.frame(read_delim("Host_prediction_phylum2.txt", "\t", escape_double = FALSE, trim_ws = TRUE, col_names = FALSE))
colnames(AllSamples_completeness_Host_phyla) <- c('contig', 'Host_phyla')
rownames(AllSamples_completeness_Host_phyla) <- AllSamples_completeness_Host_phyla$contig
AllSamples_completeness_Host_phyla$contig <- NULL
rownames(AllSamples_completeness_Host_phyla) <- gsub("\\.{1}", "_", rownames(AllSamples_completeness_Host_phyla))

setwd("/Volumes/Daan_2/Desktop/DB/")
getwd()
dir()

# Change taxonomy to names
contigs_host_phyla <- rownames(AllSamples_completeness_Host_phyla)
host_phyla_taxonomy <- as.data.frame(getTaxonomy(ids = AllSamples_completeness_Host_phyla$Host_phyla, 'accessionTaxa.sql', desiredTaxa = c("phylum")))
rownames(host_phyla_taxonomy) <- contigs_host_phyla
colnames(host_phyla_taxonomy) <- paste0('Phyla_', colnames(host_phyla_taxonomy))
host_phyla_taxonomy[is.na(host_phyla_taxonomy)] <- "Unannotated"
AllSamples_completeness_Host_phyla <- host_phyla_taxonomy
View(AllSamples_completeness_Host_phyla)
####################################
# 8.2 Host Family
####################################
setwd ("/Users/daan/Desktop/Transfer/Host") 
getwd() 
dir()

AllSamples_completeness_Host_family <- as.data.frame(read_delim("Host_prediction_family2.txt", "\t", escape_double = FALSE, trim_ws = TRUE, col_names = FALSE))
colnames(AllSamples_completeness_Host_family) <- c('contig', 'Host_family')
rownames(AllSamples_completeness_Host_family) <- AllSamples_completeness_Host_family$contig
AllSamples_completeness_Host_family$contig <- NULL
rownames(AllSamples_completeness_Host_family) <- gsub("\\.{1}", "_", rownames(AllSamples_completeness_Host_family))

setwd("/Volumes/Daan_2/Desktop/DB/")
getwd()
dir()

# Change taxonomy to names
contigs_host_family <- rownames(AllSamples_completeness_Host_family)
host_family_taxonomy <- as.data.frame(getTaxonomy(ids = AllSamples_completeness_Host_family$Host_family, 'accessionTaxa.sql', desiredTaxa = c("family")))
rownames(host_family_taxonomy) <- contigs_host_family
colnames(host_family_taxonomy) <- paste0('Phyla_', colnames(host_family_taxonomy))
host_family_taxonomy[is.na(host_family_taxonomy)] <- "Unannotated"
AllSamples_completeness_Host_family <- host_family_taxonomy
#View(AllSamples_completeness_Host_family)
####################################
# 8.2 Host Genus
####################################
setwd ("/Users/daan/Desktop/Transfer/Host") 
getwd() 
dir()

AllSamples_completeness_Host_genus <- as.data.frame(read_delim("Host_prediction_genus2.txt", "\t", escape_double = FALSE, trim_ws = TRUE, col_names = FALSE))
colnames(AllSamples_completeness_Host_genus) <- c('contig', 'Host_genus')
rownames(AllSamples_completeness_Host_genus) <- AllSamples_completeness_Host_genus$contig
AllSamples_completeness_Host_genus$contig <- NULL
rownames(AllSamples_completeness_Host_genus) <- gsub("\\.{1}", "_", rownames(AllSamples_completeness_Host_genus))

setwd("/Volumes/Daan_2/Desktop/DB/")
getwd()
dir()

# Change taxonomy to names
contigs_host_genus <- rownames(AllSamples_completeness_Host_genus)
host_genus_taxonomy <- as.data.frame(getTaxonomy(ids = AllSamples_completeness_Host_genus$Host_genus, 'accessionTaxa.sql', desiredTaxa = c("genus")))
rownames(host_genus_taxonomy) <- contigs_host_genus
colnames(host_genus_taxonomy) <- paste0('Phyla_', colnames(host_genus_taxonomy))
host_genus_taxonomy[is.na(host_genus_taxonomy)] <- "Unannotated"
AllSamples_completeness_Host_genus <- host_genus_taxonomy
# View(AllSamples_completeness_Host_genus)
####################################
# 9. Create mastertable: merge scaffolds, abundances and annotations
####################################
COMBO_1 <- merge(table, AllSamples_VirsorterCategories, by = 0, all = T)
rownames(COMBO_1) <- COMBO_1$Row.names
COMBO_1$Row.names <- NULL
dim(COMBO_1)

COMBO_1 <- merge(table, AllSamples_VirsorterCategories, by = 0, all = T)
rownames(COMBO_1) <- COMBO_1$Row.names
COMBO_1$Row.names <- NULL
dim(COMBO_1)

COMBO_2 <- merge(COMBO_1, AllSamples_classificationtaxonomy_Blastn, by = 0, all = T)
rownames(COMBO_2) <- COMBO_2$Row.names
COMBO_2$Row.names <- NULL
dim(COMBO_2)

COMBO_3 <- merge(COMBO_2, AllSamples_classificationtaxonomy_Diamond, by=0, all=T)
rownames(COMBO_3) <- COMBO_3$Row.names
COMBO_3$Row.names <- NULL
dim(COMBO_3)

COMBO_4 <- merge(COMBO_3,AllSamples_classificationtaxonomy_CAT, by=0, all=T)
rownames(COMBO_4) <- COMBO_4$Row.names
COMBO_4$Row.names <- NULL
dim(COMBO_4)

COMBO_5 <- merge(COMBO_4, AllSamples_Cenotetaker2_hallmark_proteins, by=0, all=T)
rownames(COMBO_5) <- COMBO_5$Row.names
COMBO_5$Row.names <- NULL
dim(COMBO_5)

COMBO_6 <- merge(COMBO_5, AllSamples_Cenotetaker2_lysogenic_phages, by=0, all=T)
rownames(COMBO_6) <- COMBO_6$Row.names
COMBO_6$Row.names <- NULL
dim(COMBO_6)

COMBO_7 <- merge(COMBO_6, AllSamples_Cenotetaker2_toxin_phages, by=0, all=T)
rownames(COMBO_7) <- COMBO_7$Row.names
COMBO_7$Row.names <- NULL
dim(COMBO_7)

COMBO_8 <- merge(COMBO_7, AllSamples_Cenotetaker2_recombination_phages, by=0, all=T)
rownames(COMBO_8) <- COMBO_8$Row.names
COMBO_8$Row.names <- NULL
dim(COMBO_8)

COMBO_9 <- merge(COMBO_8, AllSamples_Cenotetaker2_RT_phages, by=0, all=T)
rownames(COMBO_9) <- COMBO_9$Row.names
COMBO_9$Row.names <- NULL
dim(COMBO_9)

COMBO_10 <- merge(COMBO_9, AllSamples_Cenotetaker2_Bacteriocin_phages, by=0, all=T)
rownames(COMBO_10) <- COMBO_10$Row.names
COMBO_10$Row.names <- NULL
dim(COMBO_10)

COMBO_11 <- merge(COMBO_10, AllSamples_completeness_Checkv, by=0, all=T)
rownames(COMBO_11) <- COMBO_11$Row.names
COMBO_11$Row.names <- NULL
dim(COMBO_11)

COMBO_12 <- merge(COMBO_11, AllSamples_classificationtaxonomy_Crass, by=0, all=T)
rownames(COMBO_12) <- COMBO_12$Row.names
COMBO_12$Row.names <- NULL
dim(COMBO_12)

COMBO_13 <- merge(COMBO_12, AllSamples_classificationtaxonomy_Micro, by=0, all=T)
rownames(COMBO_13) <- COMBO_13$Row.names
COMBO_13$Row.names <- NULL
dim(COMBO_13)

COMBO_14 <- merge(COMBO_13, AllSamples_classificationtaxonomy_Ino, by=0, all=T)
rownames(COMBO_14) <- COMBO_14$Row.names
COMBO_14$Row.names <- NULL
dim(COMBO_14)

COMBO_15 <- merge(COMBO_14, AllSamples_completeness_Host_phyla, by=0, all=T)
rownames(COMBO_15) <- COMBO_15$Row.names
COMBO_15$Row.names <- NULL
dim(COMBO_15)

COMBO_16 <- merge(COMBO_15, AllSamples_completeness_Host_family, by=0, all=T)
rownames(COMBO_16) <- COMBO_16$Row.names
COMBO_16$Row.names <- NULL
dim(COMBO_16)

COMBO_17 <- merge(COMBO_16, AllSamples_completeness_Host_genus, by=0, all=T)
rownames(COMBO_17) <- COMBO_17$Row.names
COMBO_17$Row.names <- NULL
dim(COMBO_17)

COMBO_17 <- merge(COMBO_16, AllSamples_completeness_Host_genus, by=0, all=T)
rownames(COMBO_17) <- COMBO_17$Row.names
COMBO_17$Row.names <- NULL
dim(COMBO_17)

Mastertable <- merge(COMBO_17, AllSamples_viral_genera, by=0, all=T)
rownames(Mastertable) <- Mastertable$Row.names
Mastertable$Row.names <- NULL
dim(Mastertable)

Mastertable1 <- Mastertable
View(Mastertable1)
####################################
# 9.1 Additional information to Mastertable
####################################
# 9.1.1 Insert length of contigs
####################################
Mastertable$length <- rownames(Mastertable)
Mastertable$length <- gsub("*_cov.*","",Mastertable$length)
Mastertable$length <- gsub(".*length_","", Mastertable$length)
Mastertable$length <- as.numeric(Mastertable$length)
Mastertable$length <- Mastertable$length/1000
####################################
# 9.1.2 Insert total number of trimmed reads
####################################
Mastertable[is.na(Mastertable)] <- 0
Mastertable$Totalnumberofreads <- as.numeric(rowSums(Mastertable[,colnames(Mastertable) %in% names]))
Mastertable <- Mastertable[!(Mastertable$Totalnumberofreads == "0"),] ## remove rows with total number of zero
Mastertable$AD11_w3 <- NULL
Mastertable$P53_w4 <- NULL
Mastertable$P58_w12 < NULL
sort(colSums(Mastertable[1:304])) # Every sample has reads mapped

sum(Mastertable$Totalnumberofreads[Mastertable$Genus=="Genus15"])/sum(Mastertable$Totalnumberofreads)*100
sum(Mastertable$Totalnumberofreads[Mastertable$Blastn_single_TaxIdpath=="2219103"])/sum(Mastertable$Totalnumberofreads)*100
#View(Mastertable)
####################################
# 10. Taxonomical annotation for homology-based tools (Blastn, Diamond, CAT)
####################################
## Currently, ICTV is changing morphology-based taxonomy to genome-based taxonomies. This process will take quite a while.
## Because of a lot of changes that will be made, I will just add TaxID into Rstudio & convert them to the actual viral taxonomies. It is possible that some phages will be changed to an entirely different family, subfamily, genus and or species. To avoid this I will use TaxID, in which I can easily adapt this to the renewed taxonomies by updating taxonomy DB (~ see "1. Update taxonomy DB"). It gives me more flexibility to keep taxonomies up-to-date!

# PS: The reason why we start with these annotation is because we need them as critera for the phage-identification scheme (~ see later)
####################################
# 10.1 List eukaryotic viruses
###################################
## Some eukaryotic viruses have no phylum, or order annotation, but do have a family annotation. Here I create 3 lists that contain all annotated eukaryotic viruses in ncbi taxonomy.
## The reason of 3 list is because some eukaryotic viruses miss some hierarchical levels (eg. phylum) but are classified as eukaryotic viruses such as Spiraviridae.
## Picobirnaviruses are regarded as prokaryotic viruses, since we decided this in the labmeeting. Thus, it is not in this list.

Eukaryotic_viruses_phylum <- c("Cressdnaviricota")
Eukaryotic_viruses_order <- c("Herpesvirales", "Ortervirales", "Zurhausenvirales", "Piccovirales", "Algavirales", "Sepolyvirales", "Chitovirales","Hepelivirales","Sobelivirales","Ourlivirales","Stellavirales","Martellivirales","Picornavirales","Tolivirales","Ghabrivirales","Amarillovirales","Wolframvirales","Patatavirales","Reovirales","Tymovirales","Halopanivirales","Belfryvirales","Rowavirales","Pimascovirales","Asfuvirales","Baphyvirales","Polivirales","Cirlivirales","Geplafuvirales","Blubervirales","Priklausovirales", "Mulpavirales")
Eukaryotic_viruses_family <- c("Pithoviridae","Tolecusatellitidae","Spiraviridae","Sarthroviridae","Permutotetraviridae","Partitiviridae","Alphasatellitidae","Anelloviridae","Baculoviridae","Hytrosaviridae","Nimaviridae","Nudiviridae","Polydnaviridae","Amalgaviridae","Birnaviridae")
####################################
# 10.2 Blastn
####################################
# 10.2.1 Add taxonomies
####################################
setwd("/Volumes/Daan_2/Desktop/DB/")
getwd()
dir()

contigs_blastn <- rownames(AllSamples_classificationtaxonomy_Blastn)
Blastn_taxonomy <- as.data.frame(getTaxonomy(ids = AllSamples_classificationtaxonomy_Blastn$Blastn_single_TaxIdpath, 'accessionTaxa.sql', desiredTaxa = c("superkingdom","phylum", "class", "order", "family", "subfamily", "genus", "species")))
rownames(Blastn_taxonomy) <- contigs_blastn
colnames(Blastn_taxonomy) <- paste0('Blastn_', colnames(Blastn_taxonomy))
Blastn_taxonomy[is.na(Blastn_taxonomy)] <- "Unannotated"

#View(Blastn_taxonomy)
table(Blastn_taxonomy$Blastn_superkingdom)
nrow(Blastn_taxonomy)
View(Blastn_taxonomy)
####################################
# 10.2.2 Add additional hierarchy 
####################################
# We will use this additional hierarchy to divide the viral superkingdom in (i) eukaryotic and (ii) prokaryotic viruses
Blastn_taxonomy$Blastn_viral <- 1
Blastn_taxonomy <- Blastn_taxonomy[,colnames(Blastn_taxonomy)[c(1,9,2:8)]]
####################################
# 10.2.3 Annotate eukaryotic viruses 
####################################
Blastn_taxonomy$Blastn_viral[Blastn_taxonomy$Blastn_superkingdom == "Viruses" & grepl(paste(Eukaryotic_viruses_phylum, collapse = "|"), Blastn_taxonomy$Blastn_phylum)] <- "Blastn_eukaryotic_virus"
Blastn_taxonomy$Blastn_viral[Blastn_taxonomy$Blastn_superkingdom == "Viruses" & grepl(paste(Eukaryotic_viruses_order, collapse = "|"), Blastn_taxonomy$Blastn_order)] <- "Blastn_eukaryotic_virus"
Blastn_taxonomy$Blastn_viral[Blastn_taxonomy$Blastn_superkingdom == "Viruses" & grepl(paste(Eukaryotic_viruses_family, collapse = "|"), Blastn_taxonomy$Blastn_family)] <- "Blastn_eukaryotic_virus"
####################################
# 10.2.4 Annotate prokaryotic viruses 
####################################
Blastn_taxonomy$Blastn_viral[Blastn_taxonomy$Blastn_superkingdom == "Viruses" & Blastn_taxonomy$Blastn_viral == "1"] <- "Blastn_prokaryotic_virus"
Subset_Blastn_fages <- unique(Blastn_taxonomy[Blastn_taxonomy$Blastn_viral == "Blastn_prokaryotic_virus",])
# View(Subset_Blastn_fages)
####################################
# 10.2.5 Annotate all others
####################################
Blastn_taxonomy$Blastn_viral[!Blastn_taxonomy$Blastn_superkingdom == "Viruses"] <- Blastn_taxonomy$Blastn_superkingdom[!Blastn_taxonomy$Blastn_superkingdom == "Viruses"]
####################################
# 10.2.6 Add best hit
####################################
setwd ("/Users/daan/Desktop/Transfer/Taxonomy/")
getwd()
dir()
####################################
# 10.2.6.1 Add features
####################################
Blastn_best_hits <- as.data.frame(read_delim("Best_Blastn_hit_coverage_ANI.txt", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE))
colnames(Blastn_best_hits) <- c('NODES', 'Blastn_ANI', 'Blastn_coverage_absolute')
Blastn_best_hits$NODES <- gsub("\\.{1}", "_", Blastn_best_hits$NODES)
rownames(Blastn_best_hits) <- Blastn_best_hits$NODES
Blastn_best_hits$NODES <- NULL
####################################
# 10.2.6.2 Merge 
####################################
Blastn_taxonomy <- transform(merge(Blastn_taxonomy, Blastn_best_hits, by = 0, all = T), row.names=Row.names, Row.names = NULL)
Blastn_taxonomy$Blastn_ANI[is.na(Blastn_taxonomy$Blastn_ANI)] <- 0
Blastn_taxonomy$Blastn_coverage_absolute[is.na(Blastn_taxonomy$Blastn_coverage_absolute)] <- 0

## length
Blastn_taxonomy$length <- rownames(Blastn_taxonomy)
Blastn_taxonomy$length <- gsub("*_cov.*","",Blastn_taxonomy$length)
Blastn_taxonomy$length <- gsub(".*length_","", Blastn_taxonomy$length)
Blastn_taxonomy$length <- as.numeric(Blastn_taxonomy$length)

## calculate %
Blastn_taxonomy$Blastn_coverage_relative <- (Blastn_taxonomy$Blastn_coverage_absolute/Blastn_taxonomy$length)*100
Blastn_taxonomy$Blastn_coverage_relative[Blastn_taxonomy$Blastn_coverage_relative > 100] <- 100
Blastn_taxonomy$Blastn_coverage_absolute <- NULL
Blastn_taxonomy$length <- NULL
Blastn_taxonomy[is.na(Blastn_taxonomy)] <- 0
Blastn_taxonomy <- Blastn_taxonomy[!Blastn_taxonomy$Blastn_superkingdom==0,]
table(Blastn_taxonomy$Blastn_superkingdom)
####################################
# 10.2.7 Clean-up environment
####################################
rm(contigs_blastn)
rm(AlSamples_classificationtaxonomy_Blastn)
rm(Subset_Blastn_fages)
rm(Blastn_best_hits)
####################################
# 10.3 Diamond
####################################
# 10.3.1 Add taxonomies
####################################
setwd("/Volumes/Daan_2/Desktop/DB/")
getwd()
dir()

contigs_diamond <- rownames(AllSamples_classificationtaxonomy_Diamond)
Diamond_taxonomy <- as.data.frame(getTaxonomy(ids = AllSamples_classificationtaxonomy_Diamond$Diamond_single_TaxIdpath, 'accessionTaxa.sql', desiredTaxa = c("superkingdom", "phylum", "class", "order", "family", "subfamily", "genus", "species")))
rownames(Diamond_taxonomy) <- contigs_diamond
colnames(Diamond_taxonomy) <- paste0('Diamond_', colnames(Diamond_taxonomy))
Diamond_taxonomy[is.na(Diamond_taxonomy)] <- "Unannotated"
#View(Diamond_taxonomy)
####################################
# 10.3.2 Add additional hierarchy 
####################################
Diamond_taxonomy$Diamond_viral <- 1
Diamond_taxonomy <- Diamond_taxonomy[,colnames(Diamond_taxonomy)[c(1,9,2:8)]]
####################################
# 10.3.3 Annotate eukaryotic viruses 
####################################
Diamond_taxonomy$Diamond_viral[Diamond_taxonomy$Diamond_superkingdom == "Viruses" & grepl(paste(Eukaryotic_viruses_phylum, collapse = "|"), Diamond_taxonomy$Diamond_phylum)] <- "Diamond_eukaryotic_virus"
Diamond_taxonomy$Diamond_viral[Diamond_taxonomy$Diamond_superkingdom == "Viruses" & grepl(paste(Eukaryotic_viruses_order, collapse = "|"), Diamond_taxonomy$Diamond_order)] <- "Diamond_eukaryotic_virus"
Diamond_taxonomy$Diamond_viral[Diamond_taxonomy$Diamond_superkingdom == "Viruses" & grepl(paste(Eukaryotic_viruses_family, collapse = "|"), Diamond_taxonomy$Diamond_family)] <- "Diamond_eukaryotic_virus"
####################################
# 10.3.4 Annotate prokaryotic viruses 
####################################
Diamond_taxonomy$Diamond_viral[Diamond_taxonomy$Diamond_superkingdom == "Viruses" & Diamond_taxonomy$Diamond_viral == "1"] <- "Diamond_prokaryotic_virus"
Subset_Diamond_fages <- unique(Diamond_taxonomy[Diamond_taxonomy$Diamond_viral == "Diamond_prokaryotic_virus",])
#View(Subset_Diamond_fages)
####################################
# 10.3.5 Annotate all others
####################################
Diamond_taxonomy$Diamond_viral[!Diamond_taxonomy$Diamond_superkingdom == "Viruses"] <- Diamond_taxonomy$Diamond_superkingdom[!Diamond_taxonomy$Diamond_superkingdom == "Viruses"]
#View(Diamond_taxonomy)
####################################
# 10.3.6 Add best hit
####################################
setwd ("/Users/daan/Desktop/Transfer/Taxonomy/")
getwd()
dir()
####################################
# 10.3.6.1 Add features
####################################
Diamond_best_hits <- as.data.frame(read_delim("Best_Diamond_hit_coverage_ANI.txt", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE))
colnames(Diamond_best_hits) <- c('NODES', 'Diamond_ANI', 'Diamond_coverage_absolute')
Diamond_best_hits <- Diamond_best_hits[Diamond_best_hits$NODES %in% contigs_project,]
Diamond_best_hits$NODES <- gsub("\\.{1}", "_", Diamond_best_hits$NODES)
rownames(Diamond_best_hits) <- Diamond_best_hits$NODES
Diamond_best_hits$NODES <- NULL
####################################
# 10.3.6.2 Merge
####################################
Diamond_taxonomy <- transform(merge(Diamond_taxonomy, Diamond_best_hits, by = 0, all = T), row.names=Row.names, Row.names = NULL)
Diamond_taxonomy$Diamond_ANI[is.na(Diamond_taxonomy$Diamond_ANI)] <- 0
Diamond_taxonomy$Diamond_coverage_absolute[is.na(Diamond_taxonomy$Diamond_coverage_absolute)] <- 0
Diamond_taxonomy$Diamond_coverage_absolute <- Diamond_taxonomy$Diamond_coverage_absolute*3 # because it is on an aa level

## length
Diamond_taxonomy$length <- rownames(Diamond_taxonomy)
Diamond_taxonomy$length <- gsub("*_cov.*","",Diamond_taxonomy$length)
Diamond_taxonomy$length <- gsub(".*length_","", Diamond_taxonomy$length)
Diamond_taxonomy$length <- as.numeric(Diamond_taxonomy$length)

## calculate %
Diamond_taxonomy$Diamond_coverage_relative <- (Diamond_taxonomy$Diamond_coverage_absolute/Diamond_taxonomy$length)*100
Diamond_taxonomy$Diamond_coverage_relative[Diamond_taxonomy$Diamond_coverage_relative > 100] <- 100
Diamond_taxonomy$Diamond_coverage_absolute <- NULL
Diamond_taxonomy$length <- NULL
table(Diamond_taxonomy$Diamond_superkingdom)
table(Diamond_taxonomy$Diamond_class)
####################################
# 10.3.7 Clean-up environment
####################################
rm(contigs_diamond)
rm(Subset_Diamond_fages)
rm(Diamond_best_hits)
####################################
# 10.4 CAT
####################################
# 10.4.1 Add taxonomies
####################################
setwd("/Volumes/Daan_2/Desktop/DB/")
getwd()
dir()

contigs_CAT <- rownames(AllSamples_classificationtaxonomy_CAT)
CAT_taxonomy <- as.data.frame(getTaxonomy(ids = AllSamples_classificationtaxonomy_CAT$CAT_single_TaxIdpath, 'accessionTaxa.sql', desiredTaxa = c("superkingdom", "phylum", "class", "order", "family", "subfamily", "genus", "species")))
rownames(CAT_taxonomy) <- contigs_CAT
colnames(CAT_taxonomy) <- paste0('CAT_', colnames(CAT_taxonomy))
CAT_taxonomy[is.na(CAT_taxonomy)] <- "Unannotated"
####################################
# 10.4.2 Add additional hierarchy 
####################################
CAT_taxonomy$CAT_viral <- 1
CAT_taxonomy <- CAT_taxonomy[,colnames(CAT_taxonomy)[c(1,9,2:8)]]
####################################
# 10.4.3 Annotate eukaryotic viruses 
####################################
CAT_taxonomy$CAT_viral[CAT_taxonomy$CAT_superkingdom == "Viruses" & grepl(paste(Eukaryotic_viruses_phylum, collapse = "|"), CAT_taxonomy$CAT_phylum)] <- "CAT_eukaryotic_virus"
CAT_taxonomy$CAT_viral[CAT_taxonomy$CAT_superkingdom == "Viruses" & grepl(paste(Eukaryotic_viruses_order, collapse = "|"), CAT_taxonomy$CAT_order)] <- "CAT_eukaryotic_virus"
CAT_taxonomy$CAT_viral[CAT_taxonomy$CAT_superkingdom == "Viruses" & grepl(paste(Eukaryotic_viruses_family, collapse = "|"), CAT_taxonomy$CAT_family)] <- "CAT_eukaryotic_virus"
####################################
# 10.4.4 Annotate prokaryotic viruses 
####################################
CAT_taxonomy$CAT_viral[CAT_taxonomy$CAT_superkingdom == "Viruses" & CAT_taxonomy$CAT_viral == "1"] <- "CAT_prokaryotic_virus"
Subset_CAT_fages <- unique(CAT_taxonomy[CAT_taxonomy$CAT_viral == "CAT_prokaryotic_virus",])
#View(Subset_CAT_fages)
####################################
# 10.4.5 Annotate all others
####################################
CAT_taxonomy$CAT_viral[!CAT_taxonomy$CAT_superkingdom == "Viruses"] <- CAT_taxonomy$CAT_superkingdom[!CAT_taxonomy$CAT_superkingdom == "Viruses"]
table(CAT_taxonomy$CAT_superkingdom)
#View(CAT_taxonomy)
####################################
# 10.5 CNoteTaker2 annotation
####################################
# 10.5.1 Add taxonomies
####################################
setwd("/Volumes/Daan_2/Desktop/DB/")
getwd()
dir()

# Change older Caudoviricetes taxonomy with newer one
AllSamples_classificationtaxonomy_Cenotaker2$Cenote_single_TaxIdpath[AllSamples_classificationtaxonomy_Cenotaker2$Cenote_single_TaxIdpath=="28883"] <- "2731619" 
contigs_cenote <- rownames(AllSamples_classificationtaxonomy_Cenotaker2)
cenote_taxonomy <- as.data.frame(getTaxonomy(ids = AllSamples_classificationtaxonomy_Cenotaker2$Cenote_single_TaxIdpath, 'accessionTaxa.sql', desiredTaxa = c("superkingdom", "phylum", "class", "order", "family")))
rownames(cenote_taxonomy) <- contigs_cenote
colnames(cenote_taxonomy) <- paste0('Cenote_', colnames(cenote_taxonomy))
cenote_taxonomy[is.na(cenote_taxonomy)] <- "Unannotated"
####################################
# 10.5.2 Add additional hierarchy 
####################################
cenote_taxonomy$Cenote_viral <- 1
cenote_taxonomy <- cenote_taxonomy[,colnames(cenote_taxonomy)[c(1,6,2:5)]]
####################################
# 10.5.3 Annotate eukaryotic viruses 
####################################
cenote_taxonomy$Cenote_viral[cenote_taxonomy$Cenote_superkingdom == "Viruses" & grepl(paste(Eukaryotic_viruses_phylum, collapse = "|"), cenote_taxonomy$Cenote_phylum)] <- "Cenote_eukaryotic_virus"
cenote_taxonomy$Cenote_viral[cenote_taxonomy$Cenote_superkingdom == "Viruses" & grepl(paste(Eukaryotic_viruses_order, collapse = "|"), cenote_taxonomy$Cenote_order)] <- "Cenote_eukaryotic_virus"
cenote_taxonomy$Cenote_viral[cenote_taxonomy$Cenote_superkingdom == "Viruses" & grepl(paste(Eukaryotic_viruses_family, collapse = "|"), cenote_taxonomy$Cenote_family)] <- "Cenote_eukaryotic_virus"
####################################
# 10.5.4 Annotate prokaryotic viruses 
####################################
cenote_taxonomy$Cenote_viral[cenote_taxonomy$Cenote_superkingdom == "Viruses" & cenote_taxonomy$Cenote_viral == "1"] <- "cenote_prokaryotic_virus"
Subset_cenote_fages <- unique(cenote_taxonomy[cenote_taxonomy$Cenote_viral == "cenote_prokaryotic_virus",])
#View(Subset_cenote_fages)
####################################
# 10.5.5 Annotate all others
####################################
cenote_taxonomy$Cenote_viral[!cenote_taxonomy$Cenote_superkingdom == "Viruses"] <- cenote_taxonomy$Cenote_superkingdom[!cenote_taxonomy$Cenote_superkingdom == "Viruses"]
table(cenote_taxonomy$Cenote_superkingdom)
#View(cenote_taxonomy)
####################################
# 10.6 Merge taxonomies with mastertable
####################################
COMBO_1 <- merge(Mastertable, Blastn_taxonomy, by = 0, all = T)
rownames(COMBO_1) <- COMBO_1$Row.names
COMBO_1$Row.names <- NULL
dim(COMBO_1)

COMBO_2 <- merge(COMBO_1, Diamond_taxonomy, by = 0, all = T)
rownames(COMBO_2) <- COMBO_2$Row.names
COMBO_2$Row.names <- NULL
dim(COMBO_2)

COMBO_3 <- merge(COMBO_2, CAT_taxonomy, by = 0, all = T)
rownames(COMBO_3) <- COMBO_3$Row.names
COMBO_3$Row.names <- NULL
dim(COMBO_3)

Mastertable <- merge(COMBO_3, cenote_taxonomy, by = 0, all = T)
rownames(Mastertable) <- Mastertable$Row.names
Mastertable$Row.names <- NULL
dim(Mastertable)  
#View(Mastertable)

# Clean environment
rm(Mastertable_corrected)
rm(Blastn_taxonomy)
rm(Diamond_taxonomy)
rm(CAT_taxonomy)
####################################
# 11. Phage-identification 
####################################
Mastertable$Phage_score <- 1
Mastertable$Phage_score[Mastertable$Virsorter == "dsDNAphage" & Mastertable$completeness >= 50] <- "dsDNA phage"
Mastertable$Phage_score[Mastertable$Virsorter == "RNA" & Mastertable$completeness >= 50] <- "RNA phage"
Mastertable$Phage_score[Mastertable$Virsorter == "ssDNA" & Mastertable$length >= 3] <- "ssDNA phage"

Mastertable$Phages <- 1
Mastertable$Phages[Mastertable$Phage_score == "dsDNA phage" | Mastertable$Phage_score == "ssDNA phage" | Mastertable$Phage_score == "RNA phage"] <- "phage"
table(Mastertable$Phages) 
nrow(Mastertable)

# total: 2391 phages identified using this method
# total NR-contigs = 125.356 contigs
Mastertable$Totalnumberofreads <- as.numeric(Mastertable$Totalnumberofreads)
Mastertable$Totalnumberofreads[is.na(Mastertable$Totalnumberofreads)] <- 0
sum(Mastertable$Totalnumberofreads[Mastertable$Phages == "phage"]) # total of 3070648047 reads (~ or 2,6B PE reads)
sum(Mastertable$Totalnumberofreads)
(3091015377/4692786841)*100
####################################
# 11.5.4 Clean-up environment
####################################
rm(vector_1)
rm(vector_2)
rm(vector_3)
rm(vector_4)
####################################
# 12. Final classification
####################################
# The next step is to determine a final classification of our contigs. Now, we will handle contradicting taxonomical annotation for contigs. Because we use KtClassifyBlast for all contigs we have the LCA of all the annotation that were given with Blastn, Diamond (~ Blastx) and CAT.
# Because Diamond thrives on a non-redundant protein database and protein are more conserved than nucleotides we will prioritize diamond over blastn over CAT. CAT also uses the Diamond database but is less accurate (~ also experience showed that it doesn't add too much) then the two first tools. CAT is not used for eukaryotic viral annotation.

# Down below you can see I start with creating a final classification for eukaryotic viruses on viral level. 
# Next, I subset all these eukaryotic viruses and complete the final classification for all the other levels.
####################################
# 12.1 Final classification of eukaryotic viruses
####################################
# 12.1.1 Viral
####################################
#View(Mastertable)
Mastertable$Final_viral <- 1

Mastertable$Diamond_viral[is.na(Mastertable$Diamond_viral)] <- 0
Mastertable$Blastn_viral[is.na(Mastertable$Blastn_viral)] <- 0
Mastertable$CAT_viral[is.na(Mastertable$CAT_viral)] <- 0
Mastertable$Cenote_viral[is.na(Mastertable$Cenote_viral)] <- 0

Mastertable$Final_viral[Mastertable$Diamond_viral == "Diamond_eukaryotic_virus"] <- Mastertable$Diamond_viral[Mastertable$Diamond_viral == "Diamond_eukaryotic_virus"]
Mastertable$Final_viral[!Mastertable$Diamond_viral == "Diamond_eukaryotic_virus" & Mastertable$Blastn_viral == "Blastn_eukaryotic_virus"] <- Mastertable$Blastn_viral[!Mastertable$Diamond_viral == "Diamond_eukaryotic_virus" & Mastertable$Blastn_viral == "Blastn_eukaryotic_virus"]
Mastertable$Final_viral[Mastertable$Final_viral == "Diamond_eukaryotic_virus" | Mastertable$Final_viral == "Blastn_eukaryotic_virus"] <- "eukaryotic virus"
table(Mastertable$Final_viral) 
#View(Mastertable)
## How many contigs are classified as eukaryotic viruses in our final classification?
####################################
# 12.1.2 Superkingdom
####################################
Mastertable$Final_superkingdom <- 1
Mastertable$Final_superkingdom[Mastertable$Diamond_viral == "Diamond_eukaryotic_virus"] <- Mastertable$Diamond_superkingdom[Mastertable$Diamond_viral == "Diamond_eukaryotic_virus"]
Mastertable$Final_superkingdom[!Mastertable$Diamond_viral == "Diamond_eukaryotic_virus" & Mastertable$Blastn_viral == "Blastn_eukaryotic_virus"] <- Mastertable$Blastn_superkingdom[!Mastertable$Diamond_viral == "Diamond_eukaryotic_virus" & Mastertable$Blastn_viral == "Blastn_eukaryotic_virus"]
table(Mastertable$Final_superkingdom) 
####################################
# 12.1.3 Phylum
####################################
Mastertable$Final_phylum <- 1
Mastertable$Final_phylum[Mastertable$Diamond_viral == "Diamond_eukaryotic_virus"] <- Mastertable$Diamond_phylum[Mastertable$Diamond_viral == "Diamond_eukaryotic_virus"]
Mastertable$Final_phylum[!Mastertable$Diamond_viral == "Diamond_eukaryotic_virus" & Mastertable$Blastn_viral == "Blastn_eukaryotic_virus"] <- Mastertable$Blastn_phylum[!Mastertable$Diamond_viral == "Diamond_eukaryotic_virus" & Mastertable$Blastn_viral == "Blastn_eukaryotic_virus"]
table(Mastertable$Final_phylum) 
# How many different eukaryotic viral phyla are found?
####################################
# 12.1.4 Class
####################################
Mastertable$Final_class <- 1
Mastertable$Final_class[Mastertable$Diamond_viral == "Diamond_eukaryotic_virus"] <- Mastertable$Diamond_class[Mastertable$Diamond_viral == "Diamond_eukaryotic_virus"]
Mastertable$Final_class[!Mastertable$Diamond_viral == "Diamond_eukaryotic_virus" & Mastertable$Blastn_viral == "Blastn_eukaryotic_virus"] <- Mastertable$Blastn_class[!Mastertable$Diamond_viral == "Diamond_eukaryotic_virus" & Mastertable$Blastn_viral == "Blastn_eukaryotic_virus"]
table(Mastertable$Final_class) 
# How many different eukaryotic viral classes are found?
####################################
# 12.1.5 Order
####################################
Mastertable$Final_order <- 1
Mastertable$Final_order[Mastertable$Diamond_viral == "Diamond_eukaryotic_virus"] <- Mastertable$Diamond_order[Mastertable$Diamond_viral == "Diamond_eukaryotic_virus"]
Mastertable$Final_order[!Mastertable$Diamond_viral == "Diamond_eukaryotic_virus" & Mastertable$Blastn_viral == "Blastn_eukaryotic_virus"] <- Mastertable$Blastn_order[!Mastertable$Diamond_viral == "Diamond_eukaryotic_virus" & Mastertable$Blastn_viral == "Blastn_eukaryotic_virus"]
table(Mastertable$Final_order) 
# How many different eukaryotic viral orders are found?
####################################
# 12.1.6 Family
####################################
Mastertable$Final_family <- 1
Mastertable$Final_family[Mastertable$Diamond_viral == "Diamond_eukaryotic_virus"] <- Mastertable$Diamond_family[Mastertable$Diamond_viral == "Diamond_eukaryotic_virus"]
Mastertable$Final_family[!Mastertable$Diamond_viral == "Diamond_eukaryotic_virus" & Mastertable$Blastn_viral == "Blastn_eukaryotic_virus"] <- Mastertable$Blastn_family[!Mastertable$Diamond_viral == "Diamond_eukaryotic_virus" & Mastertable$Blastn_viral == "Blastn_eukaryotic_virus"]
table(Mastertable$Final_family) 
# How many different eukaryotic viral families are found?
####################################
# 12.1.7 Subfamily
####################################
Mastertable$Final_subfamily <- 1
Mastertable$Final_subfamily[Mastertable$Diamond_viral == "Diamond_eukaryotic_virus"] <- Mastertable$Diamond_subfamily[Mastertable$Diamond_viral == "Diamond_eukaryotic_virus"]
Mastertable$Final_subfamily[!Mastertable$Diamond_viral == "Diamond_eukaryotic_virus" & Mastertable$Blastn_viral == "Blastn_eukaryotic_virus"] <- Mastertable$Blastn_subfamily[!Mastertable$Diamond_viral == "Diamond_eukaryotic_virus" & Mastertable$Blastn_viral == "Blastn_eukaryotic_virus"]
table(Mastertable$Final_subfamily) 
# How many different eukaryotic subfamilies orders are found?
####################################
# 12.1.8 Genus
####################################
Mastertable$Final_genus <- 1
Mastertable$Final_genus[Mastertable$Diamond_viral == "Diamond_eukaryotic_virus"] <- Mastertable$Diamond_genus[Mastertable$Diamond_viral == "Diamond_eukaryotic_virus"]
Mastertable$Final_genus[!Mastertable$Diamond_viral == "Diamond_eukaryotic_virus" & Mastertable$Blastn_viral == "Blastn_eukaryotic_virus"] <- Mastertable$Blastn_genus[!Mastertable$Diamond_viral == "Diamond_eukaryotic_virus" & Mastertable$Blastn_viral == "Blastn_eukaryotic_virus"]
table(Mastertable$Final_genus) 
# How many different eukaryotic viral genera are found?
####################################
# 12.1.9 Species
####################################
Mastertable$Final_species <- 1
Mastertable$Final_species[Mastertable$Diamond_viral == "Diamond_eukaryotic_virus"] <- Mastertable$Diamond_species[Mastertable$Diamond_viral == "Diamond_eukaryotic_virus"]
Mastertable$Final_species[!Mastertable$Diamond_viral == "Diamond_eukaryotic_virus" & Mastertable$Blastn_viral == "Blastn_eukaryotic_virus"] <- Mastertable$Blastn_species[!Mastertable$Diamond_viral == "Diamond_eukaryotic_virus" & Mastertable$Blastn_viral == "Blastn_eukaryotic_virus"]
table(Mastertable$Final_species) 
# How many different eukaryotic viral species are found?
####################################
# 12.1.10 alignment quality (ANI, coverage)
####################################
Mastertable$Final_ANI <- 1
Mastertable$Final_ANI[Mastertable$Diamond_viral == "Diamond_eukaryotic_virus"] <- Mastertable$Diamond_ANI[Mastertable$Diamond_viral == "Diamond_eukaryotic_virus"]
Mastertable$Final_ANI[!Mastertable$Diamond_viral == "Diamond_eukaryotic_virus" & Mastertable$Blastn_viral == "Blastn_eukaryotic_virus"] <- Mastertable$Blastn_ANI[!Mastertable$Diamond_viral == "Diamond_eukaryotic_virus" & Mastertable$Blastn_viral == "Blastn_eukaryotic_virus"]

Mastertable$Final_coverage <- 1
Mastertable$Final_coverage[Mastertable$Diamond_viral == "Diamond_eukaryotic_virus"] <- Mastertable$Diamond_coverage_relative[Mastertable$Diamond_viral == "Diamond_eukaryotic_virus"]
Mastertable$Final_coverage[!Mastertable$Diamond_viral == "Diamond_eukaryotic_virus" & Mastertable$Blastn_viral == "Blastn_eukaryotic_virus"] <- Mastertable$Blastn_coverage_relative[!Mastertable$Diamond_viral == "Diamond_eukaryotic_virus" & Mastertable$Blastn_viral == "Blastn_eukaryotic_virus"]

#View(Mastertable)
####################################
# 12.2 Final classification of prokaryotic viruses
####################################
# 12.2.1 Viral
####################################
Mastertable$Final_viral[Mastertable$Phages == "phage" & !Mastertable$Diamond_viral == "Diamond_eukaryotic_virus" & !Mastertable$Blastn_viral == "Blastn_eukaryotic_virus"] <- Mastertable$Phages[Mastertable$Phages == "phage" & !Mastertable$Diamond_viral == "Diamond_eukaryotic_virus" & !Mastertable$Blastn_viral == "Blastn_eukaryotic_virus"]
table(Mastertable$Phage_score)
table(Mastertable$Final_viral)
####################################
# 12.2.2 Superkingdom
####################################
Mastertable$Final_superkingdom[Mastertable$Final_viral == "phage" & Mastertable$Diamond_viral == "Diamond_prokaryotic_virus"] <- Mastertable$Diamond_superkingdom[Mastertable$Final_viral == "phage" & Mastertable$Diamond_viral == "Diamond_prokaryotic_virus"]
Mastertable$Final_superkingdom[Mastertable$Final_viral == "phage" & !Mastertable$Diamond_viral == "Diamond_prokaryotic_virus" & Mastertable$Blastn_viral == "Blastn_prokaryotic_virus"] <- Mastertable$Blastn_superkingdom[Mastertable$Final_viral == "phage" & !Mastertable$Diamond_viral == "Diamond_prokaryotic_virus" & Mastertable$Blastn_viral == "Blastn_prokaryotic_virus"]
Mastertable$Final_superkingdom[Mastertable$Final_viral == "phage" & !Mastertable$Diamond_viral == "Diamond_prokaryotic_virus" & !Mastertable$Blastn_viral == "Blastn_prokaryotic_virus" & Mastertable$CAT_viral == "CAT_prokaryotic_virus"] <- Mastertable$CAT_superkingdom[Mastertable$Final_viral == "phage" & !Mastertable$Diamond_viral == "Diamond_prokaryotic_virus" & !Mastertable$Blastn_viral == "Blastn_prokaryotic_virus" & Mastertable$CAT_viral == "CAT_prokaryotic_virus"]
Mastertable$Final_superkingdom[Mastertable$Final_viral == "phage" & !Mastertable$Diamond_viral == "Diamond_prokaryotic_virus" & !Mastertable$Blastn_viral == "Blastn_prokaryotic_virus" & !Mastertable$CAT_viral == "CAT_prokaryotic_virus"] <- "Viruses"

#View(Mastertable)
table(Mastertable$Final_superkingdom[Mastertable$Final_viral == "phage" & !Mastertable$Diamond_viral == "Diamond_prokaryotic_virus" & !Mastertable$Blastn_viral == "Blastn_prokaryotic_virus" & !Mastertable$CAT_viral == "CAT_prokaryotic_virus"])
# How many phages are annotated with homology-based tools?
####################################
# 12.2.3 Phylum
####################################
Mastertable$Final_phylum[Mastertable$Final_viral == "phage" & Mastertable$Diamond_viral == "Diamond_prokaryotic_virus"] <- Mastertable$Diamond_phylum[Mastertable$Final_viral == "phage" & Mastertable$Diamond_viral == "Diamond_prokaryotic_virus"]
Mastertable$Final_phylum[Mastertable$Final_viral == "phage" & !Mastertable$Diamond_viral == "Diamond_prokaryotic_virus" & Mastertable$Blastn_viral == "Blastn_prokaryotic_virus"] <- Mastertable$Blastn_phylum[Mastertable$Final_viral == "phage" & !Mastertable$Diamond_viral == "Diamond_prokaryotic_virus" & Mastertable$Blastn_viral == "Blastn_prokaryotic_virus"]
Mastertable$Final_phylum[Mastertable$Final_viral == "phage" & !Mastertable$Diamond_viral == "Diamond_prokaryotic_virus" & !Mastertable$Blastn_viral == "Blastn_prokaryotic_virus" & Mastertable$CAT_viral == "CAT_prokaryotic_virus"] <- Mastertable$CAT_phylum[Mastertable$Final_viral == "phage" & !Mastertable$Diamond_viral == "Diamond_prokaryotic_virus" & !Mastertable$Blastn_viral == "Blastn_prokaryotic_virus" & Mastertable$CAT_viral == "CAT_prokaryotic_virus"]
Mastertable$Final_phylum[Mastertable$Final_viral == "phage" & !Mastertable$Diamond_viral == "Diamond_prokaryotic_virus" & !Mastertable$Blastn_viral == "Blastn_prokaryotic_virus" & !Mastertable$CAT_viral == "CAT_prokaryotic_virus"] <- "Unannotated"

#View(Mastertable)
table(Mastertable$Final_phylum)
# How many phages phyla are found?
####################################
# 12.2.4 Class
####################################
Mastertable$Final_class[Mastertable$Final_viral == "phage" & Mastertable$Diamond_viral == "Diamond_prokaryotic_virus"] <- Mastertable$Diamond_class[Mastertable$Final_viral == "phage" & Mastertable$Diamond_viral == "Diamond_prokaryotic_virus"]
Mastertable$Final_class[Mastertable$Final_viral == "phage" & !Mastertable$Diamond_viral == "Diamond_prokaryotic_virus" & Mastertable$Blastn_viral == "Blastn_prokaryotic_virus"] <- Mastertable$Blastn_class[Mastertable$Final_viral == "phage" & !Mastertable$Diamond_viral == "Diamond_prokaryotic_virus" & Mastertable$Blastn_viral == "Blastn_prokaryotic_virus"]
Mastertable$Final_class[Mastertable$Final_viral == "phage" & !Mastertable$Diamond_viral == "Diamond_prokaryotic_virus" & !Mastertable$Blastn_viral == "Blastn_prokaryotic_virus" & Mastertable$CAT_viral == "CAT_prokaryotic_virus"] <- Mastertable$CAT_class[Mastertable$Final_viral == "phage" & !Mastertable$Diamond_viral == "Diamond_prokaryotic_virus" & !Mastertable$Blastn_viral == "Blastn_prokaryotic_virus" & Mastertable$CAT_viral == "CAT_prokaryotic_virus"]
Mastertable$Final_class[Mastertable$Final_viral == "phage" & !Mastertable$Diamond_viral == "Diamond_prokaryotic_virus" & !Mastertable$Blastn_viral == "Blastn_prokaryotic_virus" & !Mastertable$CAT_viral == "CAT_prokaryotic_virus"] <- "Unannotated"

#View(Mastertable)
table(Mastertable$Final_class)
# How many phages classes are found?
####################################
# 12.2.5 Order
####################################
Mastertable$Final_order[Mastertable$Final_viral == "phage" & Mastertable$Diamond_viral == "Diamond_prokaryotic_virus"] <- Mastertable$Diamond_order[Mastertable$Final_viral == "phage" & Mastertable$Diamond_viral == "Diamond_prokaryotic_virus"]
Mastertable$Final_order[Mastertable$Final_viral == "phage" & !Mastertable$Diamond_viral == "Diamond_prokaryotic_virus" & Mastertable$Blastn_viral == "Blastn_prokaryotic_virus"] <- Mastertable$Blastn_order[Mastertable$Final_viral == "phage" & !Mastertable$Diamond_viral == "Diamond_prokaryotic_virus" & Mastertable$Blastn_viral == "Blastn_prokaryotic_virus"]
Mastertable$Final_order[Mastertable$Final_viral == "phage" & !Mastertable$Diamond_viral == "Diamond_prokaryotic_virus" & !Mastertable$Blastn_viral == "Blastn_prokaryotic_virus" & Mastertable$CAT_viral == "CAT_prokaryotic_virus"] <- Mastertable$CAT_order[Mastertable$Final_viral == "phage" & !Mastertable$Diamond_viral == "Diamond_prokaryotic_virus" & !Mastertable$Blastn_viral == "Blastn_prokaryotic_virus" & Mastertable$CAT_viral == "CAT_prokaryotic_virus"]
Mastertable$Final_order[Mastertable$Final_viral == "phage" & !Mastertable$Diamond_viral == "Diamond_prokaryotic_virus" & !Mastertable$Blastn_viral == "Blastn_prokaryotic_virus" & !Mastertable$CAT_viral == "CAT_prokaryotic_virus"] <- "Unannotated"

#View(Mastertable)
table(Mastertable$Final_order)
# How many phages orders are found?
####################################
# 12.2.6 Family
####################################
Mastertable$Final_family[Mastertable$Final_viral == "phage" & Mastertable$Diamond_viral == "Diamond_prokaryotic_virus"] <- Mastertable$Diamond_family[Mastertable$Final_viral == "phage" & Mastertable$Diamond_viral == "Diamond_prokaryotic_virus"]
Mastertable$Final_family[Mastertable$Final_viral == "phage" & !Mastertable$Diamond_viral == "Diamond_prokaryotic_virus" & Mastertable$Blastn_viral == "Blastn_prokaryotic_virus"] <- Mastertable$Blastn_family[Mastertable$Final_viral == "phage" & !Mastertable$Diamond_viral == "Diamond_prokaryotic_virus" & Mastertable$Blastn_viral == "Blastn_prokaryotic_virus"]
Mastertable$Final_family[Mastertable$Final_viral == "phage" & !Mastertable$Diamond_viral == "Diamond_prokaryotic_virus" & !Mastertable$Blastn_viral == "Blastn_prokaryotic_virus" & Mastertable$CAT_viral == "CAT_prokaryotic_virus"] <- Mastertable$CAT_family[Mastertable$Final_viral == "phage" & !Mastertable$Diamond_viral == "Diamond_prokaryotic_virus" & !Mastertable$Blastn_viral == "Blastn_prokaryotic_virus" & Mastertable$CAT_viral == "CAT_prokaryotic_virus"]
Mastertable$Final_family[Mastertable$Final_viral == "phage" & !Mastertable$Diamond_viral == "Diamond_prokaryotic_virus" & !Mastertable$Blastn_viral == "Blastn_prokaryotic_virus" & !Mastertable$CAT_viral == "CAT_prokaryotic_virus"] <- "Unannotated"

#View(Mastertable)
table(Mastertable$Final_family)
# How many phages families are found?
####################################
# 12.2.7 Subfamily
####################################
Mastertable$Final_subfamily[Mastertable$Final_viral == "phage" & Mastertable$Diamond_viral == "Diamond_prokaryotic_virus"] <- Mastertable$Diamond_subfamily[Mastertable$Final_viral == "phage" & Mastertable$Diamond_viral == "Diamond_prokaryotic_virus"]
Mastertable$Final_subfamily[Mastertable$Final_viral == "phage" & !Mastertable$Diamond_viral == "Diamond_prokaryotic_virus" & Mastertable$Blastn_viral == "Blastn_prokaryotic_virus"] <- Mastertable$Blastn_subfamily[Mastertable$Final_viral == "phage" & !Mastertable$Diamond_viral == "Diamond_prokaryotic_virus" & Mastertable$Blastn_viral == "Blastn_prokaryotic_virus"]
Mastertable$Final_subfamily[Mastertable$Final_viral == "phage" & !Mastertable$Diamond_viral == "Diamond_prokaryotic_virus" & !Mastertable$Blastn_viral == "Blastn_prokaryotic_virus" & Mastertable$CAT_viral == "CAT_prokaryotic_virus"] <- Mastertable$CAT_subfamily[Mastertable$Final_viral == "phage" & !Mastertable$Diamond_viral == "Diamond_prokaryotic_virus" & !Mastertable$Blastn_viral == "Blastn_prokaryotic_virus" & Mastertable$CAT_viral == "CAT_prokaryotic_virus"]
Mastertable$Final_subfamily[Mastertable$Final_viral == "phage" & !Mastertable$Diamond_viral == "Diamond_prokaryotic_virus" & !Mastertable$Blastn_viral == "Blastn_prokaryotic_virus" & !Mastertable$CAT_viral == "CAT_prokaryotic_virus"] <- "Unannotated"

#View(Mastertable)
table(Mastertable$Final_subfamily)
# How many phages subfamilies are found?
####################################
# 12.2.8 Genus
####################################
Mastertable$Final_genus[Mastertable$Final_viral == "phage" & Mastertable$Diamond_viral == "Diamond_prokaryotic_virus"] <- Mastertable$Diamond_genus[Mastertable$Final_viral == "phage" & Mastertable$Diamond_viral == "Diamond_prokaryotic_virus"]
Mastertable$Final_genus[Mastertable$Final_viral == "phage" & !Mastertable$Diamond_viral == "Diamond_prokaryotic_virus" & Mastertable$Blastn_viral == "Blastn_prokaryotic_virus"] <- Mastertable$Blastn_genus[Mastertable$Final_viral == "phage" & !Mastertable$Diamond_viral == "Diamond_prokaryotic_virus" & Mastertable$Blastn_viral == "Blastn_prokaryotic_virus"]
Mastertable$Final_genus[Mastertable$Final_viral == "phage" & !Mastertable$Diamond_viral == "Diamond_prokaryotic_virus" & !Mastertable$Blastn_viral == "Blastn_prokaryotic_virus" & Mastertable$CAT_viral == "CAT_prokaryotic_virus"] <- Mastertable$CAT_genus[Mastertable$Final_viral == "phage" & !Mastertable$Diamond_viral == "Diamond_prokaryotic_virus" & !Mastertable$Blastn_viral == "Blastn_prokaryotic_virus" & Mastertable$CAT_viral == "CAT_prokaryotic_virus"]
Mastertable$Final_genus[Mastertable$Final_viral == "phage" & !Mastertable$Diamond_viral == "Diamond_prokaryotic_virus" & !Mastertable$Blastn_viral == "Blastn_prokaryotic_virus" & !Mastertable$CAT_viral == "CAT_prokaryotic_virus"] <- "Unannotated"

#View(Mastertable)
table(Mastertable$Final_genus)
# How many phages genera are found?
####################################
# 12.2.9 Species
####################################
Mastertable$Final_species[Mastertable$Final_viral == "phage" & Mastertable$Diamond_viral == "Diamond_prokaryotic_virus"] <- Mastertable$Diamond_species[Mastertable$Final_viral == "phage" & Mastertable$Diamond_viral == "Diamond_prokaryotic_virus"]
Mastertable$Final_species[Mastertable$Final_viral == "phage" & !Mastertable$Diamond_viral == "Diamond_prokaryotic_virus" & Mastertable$Blastn_viral == "Blastn_prokaryotic_virus"] <- Mastertable$Blastn_species[Mastertable$Final_viral == "phage" & !Mastertable$Diamond_viral == "Diamond_prokaryotic_virus" & Mastertable$Blastn_viral == "Blastn_prokaryotic_virus"]
Mastertable$Final_species[Mastertable$Final_viral == "phage" & !Mastertable$Diamond_viral == "Diamond_prokaryotic_virus" & !Mastertable$Blastn_viral == "Blastn_prokaryotic_virus" & Mastertable$CAT_viral == "CAT_prokaryotic_virus"] <- Mastertable$CAT_species[Mastertable$Final_viral == "phage" & !Mastertable$Diamond_viral == "Diamond_prokaryotic_virus" & !Mastertable$Blastn_viral == "Blastn_prokaryotic_virus" & Mastertable$CAT_viral == "CAT_prokaryotic_virus"]
Mastertable$Final_species[Mastertable$Final_viral == "phage" & !Mastertable$Diamond_viral == "Diamond_prokaryotic_virus" & !Mastertable$Blastn_viral == "Blastn_prokaryotic_virus" & !Mastertable$CAT_viral == "CAT_prokaryotic_virus"] <- "Unannotated"

#View(Mastertable)
table(Mastertable$Final_species)
# How many phages genera are found?
####################################
# 12.2.10 alignment quality (ANI, coverage)
####################################
Mastertable$Final_ANI[Mastertable$Final_viral == "phage" & Mastertable$Diamond_viral == "Diamond_prokaryotic_virus"] <- Mastertable$Diamond_ANI[Mastertable$Final_viral == "phage" & Mastertable$Diamond_viral == "Diamond_prokaryotic_virus"]
Mastertable$Final_ANI[Mastertable$Final_viral == "phage" & !Mastertable$Diamond_viral == "Diamond_prokaryotic_virus" & Mastertable$Blastn_viral == "Blastn_prokaryotic_virus"] <- Mastertable$Blastn_ANI[Mastertable$Final_viral == "phage" & !Mastertable$Diamond_viral == "Diamond_prokaryotic_virus" & Mastertable$Blastn_viral == "Blastn_prokaryotic_virus"]
Mastertable$Final_ANI[Mastertable$Final_viral == "phage" & !Mastertable$Diamond_viral == "Diamond_prokaryotic_virus" & !Mastertable$Blastn_viral == "Blastn_prokaryotic_virus" & !Mastertable$CAT_viral == "CAT_prokaryotic_virus"] <- "Unannotated"

Mastertable$Final_coverage[Mastertable$Final_viral == "phage" & Mastertable$Diamond_viral == "Diamond_prokaryotic_virus"] <- Mastertable$Diamond_coverage_relative[Mastertable$Final_viral == "phage" & Mastertable$Diamond_viral == "Diamond_prokaryotic_virus"]
Mastertable$Final_coverage[Mastertable$Final_viral == "phage" & !Mastertable$Diamond_viral == "Diamond_prokaryotic_virus" & Mastertable$Blastn_viral == "Blastn_prokaryotic_virus"] <- Mastertable$Blastn_coverage_relative[Mastertable$Final_viral == "phage" & !Mastertable$Diamond_viral == "Diamond_prokaryotic_virus" & Mastertable$Blastn_viral == "Blastn_prokaryotic_virus"]
Mastertable$Final_coverage[Mastertable$Final_viral == "phage" & !Mastertable$Diamond_viral == "Diamond_prokaryotic_virus" & !Mastertable$Blastn_viral == "Blastn_prokaryotic_virus" & !Mastertable$CAT_viral == "CAT_prokaryotic_virus"] <- "Unannotated"
#View(Mastertable)
####################################
# 12.3 Final classification of eukaryotes
####################################
Mastertable$Diamond_superkingdom[is.na(Mastertable$Diamond_superkingdom)] <- 0
Mastertable$Blastn_superkingdom[is.na(Mastertable$Blastn_superkingdom)] <- 0
Mastertable$CAT_superkingdom[is.na(Mastertable$CAT_superkingdom)] <- 0
Mastertable$Final_superkingdom[is.na(Mastertable$Final_superkingdom)] <- 0

Mastertable$Final_superkingdom[!Mastertable$Final_superkingdom == "Viruses" & Mastertable$Diamond_superkingdom == "Eukaryota"] <- Mastertable$Diamond_superkingdom[!Mastertable$Final_superkingdom == "Viruses" & Mastertable$Diamond_superkingdom == "Eukaryota"]
Mastertable$Final_superkingdom[!Mastertable$Final_superkingdom == "Viruses" & !Mastertable$Diamond_superkingdom == "Eukaryota" & Mastertable$Blastn_superkingdom == "Eukaryota"] <- Mastertable$Blastn_superkingdom[!Mastertable$Final_superkingdom == "Viruses" & !Mastertable$Diamond_superkingdom == "Eukaryota" & Mastertable$Blastn_superkingdom == "Eukaryota"]
Mastertable$Final_superkingdom[!Mastertable$Final_superkingdom == "Viruses" & !Mastertable$Diamond_superkingdom == "Eukaryota" & !Mastertable$Blastn_superkingdom == "Eukaryota" & Mastertable$CAT_superkingdom == "Eukaryota"] <- Mastertable$CAT_superkingdom[!Mastertable$Final_superkingdom == "Viruses" & !Mastertable$Diamond_superkingdom == "Eukaryota" & !Mastertable$Blastn_superkingdom == "Eukaryota" & Mastertable$CAT_superkingdom == "Eukaryota"]
####################################
# 12.4 Final classification of archaea
####################################
Mastertable$Final_superkingdom[!Mastertable$Final_superkingdom == "Viruses" & !Mastertable$Final_superkingdom == "Eukaryota" & Mastertable$Diamond_superkingdom == "Archaea"] <- Mastertable$Diamond_superkingdom[!Mastertable$Final_superkingdom == "Viruses" & !Mastertable$Final_superkingdom == "Eukaryota" & Mastertable$Diamond_superkingdom == "Archaea"]
Mastertable$Final_superkingdom[!Mastertable$Final_superkingdom == "Viruses" & !Mastertable$Final_superkingdom == "Eukaryota" & !Mastertable$Diamond_superkingdom == "Archaea" & Mastertable$Blastn_superkingdom == "Archaea"] <- Mastertable$Blastn_superkingdom[!Mastertable$Final_superkingdom == "Viruses" & !Mastertable$Final_superkingdom == "Eukaryota" & !Mastertable$Diamond_superkingdom == "Archaea" & Mastertable$Blastn_superkingdom == "Archaea"]
Mastertable$Final_superkingdom[!Mastertable$Final_superkingdom == "Viruses" & !Mastertable$Final_superkingdom == "Eukaryota" & !Mastertable$Diamond_superkingdom == "Archaea" & !Mastertable$Blastn_superkingdom == "Archaea" & Mastertable$CAT_superkingdom == "Archaea"] <- Mastertable$CAT_superkingdom[!Mastertable$Final_superkingdom == "Viruses" & !Mastertable$Final_superkingdom == "Eukaryota" & !Mastertable$Diamond_superkingdom == "Archaea" & !Mastertable$Blastn_superkingdom == "Archaea" & Mastertable$CAT_superkingdom == "Archaea"]
####################################
# 12.5 Final classification of bacteria
####################################
Mastertable$Final_superkingdom[!Mastertable$Final_superkingdom == "Viruses" & !Mastertable$Final_superkingdom == "Eukaryota" & !Mastertable$Final_superkingdom == "Archaea" & Mastertable$Diamond_superkingdom == "Bacteria"] <- Mastertable$Diamond_superkingdom[!Mastertable$Final_superkingdom == "Viruses" & !Mastertable$Final_superkingdom == "Eukaryota" & !Mastertable$Final_superkingdom == "Archaea" & Mastertable$Diamond_superkingdom == "Bacteria"]
Mastertable$Final_superkingdom[!Mastertable$Final_superkingdom == "Viruses" & !Mastertable$Final_superkingdom == "Eukaryota" & !Mastertable$Final_superkingdom == "Archaea" & !Mastertable$Diamond_superkingdom == "Bacteria" & Mastertable$Blastn_superkingdom == "Bacteria"] <- Mastertable$Blastn_superkingdom[!Mastertable$Final_superkingdom == "Viruses" & !Mastertable$Final_superkingdom == "Eukaryota" & !Mastertable$Final_superkingdom == "Archaea" & !Mastertable$Diamond_superkingdom == "Bacteria" & Mastertable$Blastn_superkingdom == "Bacteria"]
Mastertable$Final_superkingdom[!Mastertable$Final_superkingdom == "Viruses" & !Mastertable$Final_superkingdom == "Eukaryota" & !Mastertable$Final_superkingdom == "Archaea" & !Mastertable$Diamond_superkingdom == "Bacteria" & !Mastertable$Blastn_superkingdom == "Bacteria" & Mastertable$CAT_superkingdom == "Bacteria"] <- Mastertable$CAT_superkingdom[!Mastertable$Final_superkingdom == "Viruses" & !Mastertable$Final_superkingdom == "Eukaryota" & !Mastertable$Final_superkingdom == "Archaea" & !Mastertable$Diamond_superkingdom == "Bacteria" & !Mastertable$Blastn_superkingdom == "Bacteria" & Mastertable$CAT_superkingdom == "Bacteria"]
####################################
# 12.6 Final classification of Other
####################################
Mastertable$Final_superkingdom[!Mastertable$Final_superkingdom == "Viruses" & !Mastertable$Final_superkingdom == "Eukaryota" & !Mastertable$Final_superkingdom == "Archaea" & !Mastertable$Final_superkingdom == "Bacteria" & Mastertable$Diamond_superkingdom == "Other"] <- Mastertable$Diamond_superkingdom[!Mastertable$Final_superkingdom == "Viruses" & !Mastertable$Final_superkingdom == "Eukaryota" & !Mastertable$Final_superkingdom == "Archaea" & !Mastertable$Final_superkingdom == "Bacteria" & Mastertable$Diamond_superkingdom == "Other"]
Mastertable$Final_superkingdom[!Mastertable$Final_superkingdom == "Viruses" & !Mastertable$Final_superkingdom == "Eukaryota" & !Mastertable$Final_superkingdom == "Archaea" & !Mastertable$Final_superkingdom == "Bacteria" & !Mastertable$Diamond_superkingdom == "Other" & Mastertable$Blastn_superkingdom == "Other"] <- Mastertable$Blastn_superkingdom[!Mastertable$Final_superkingdom == "Viruses" & !Mastertable$Final_superkingdom == "Eukaryota" & !Mastertable$Final_superkingdom == "Archaea" & !Mastertable$Final_superkingdom == "Bacteria" & !Mastertable$Diamond_superkingdom == "Other" & Mastertable$Blastn_superkingdom == "Other"]
Mastertable$Final_superkingdom[!Mastertable$Final_superkingdom == "Viruses" & !Mastertable$Final_superkingdom == "Eukaryota" & !Mastertable$Final_superkingdom == "Archaea" & !Mastertable$Final_superkingdom == "Bacteria" & !Mastertable$Diamond_superkingdom == "Other" & !Mastertable$Blastn_superkingdom == "Other" & Mastertable$CAT_superkingdom == "Other"] <- Mastertable$CAT_superkingdom[!Mastertable$Final_superkingdom == "Viruses" & !Mastertable$Final_superkingdom == "Eukaryota" & !Mastertable$Final_superkingdom == "Archaea" & !Mastertable$Final_superkingdom == "Bacteria" & !Mastertable$Diamond_superkingdom == "Other" & !Mastertable$Blastn_superkingdom == "Other" & Mastertable$CAT_superkingdom == "Other"]
####################################
# 12.7 Final classification of Unclassified
####################################
# You could, if you like, scan the "unclassified" sequences manully that have an annotation to check if you have some unclassified viral sequences. Since this we take a lot of time and makes it impossible to automize this pipeline & it doesn't add a lot I don't do this.
Mastertable$Final_superkingdom[!Mastertable$Final_superkingdom == "Viruses" & !Mastertable$Final_superkingdom == "Eukaryota" & !Mastertable$Final_superkingdom == "Archaea" & !Mastertable$Final_superkingdom == "Bacteria" & !Mastertable$Diamond_superkingdom == "Other" & Mastertable$Diamond_superkingdom == "Unannotated"] <- Mastertable$Diamond_superkingdom[!Mastertable$Final_superkingdom == "Viruses" & !Mastertable$Final_superkingdom == "Eukaryota" & !Mastertable$Final_superkingdom == "Archaea" & !Mastertable$Final_superkingdom == "Bacteria" & Mastertable$Diamond_superkingdom == "Unannotated"]
Mastertable$Final_superkingdom[!Mastertable$Final_superkingdom == "Viruses" & !Mastertable$Final_superkingdom == "Eukaryota" & !Mastertable$Final_superkingdom == "Archaea" & !Mastertable$Final_superkingdom == "Bacteria" & !Mastertable$Diamond_superkingdom == "Other" & !Mastertable$Diamond_superkingdom == "Unannotated" & Mastertable$Blastn_superkingdom == "Unannotated"] <- Mastertable$Blastn_superkingdom[!Mastertable$Final_superkingdom == "Viruses" & !Mastertable$Final_superkingdom == "Eukaryota" & !Mastertable$Final_superkingdom == "Archaea" & !Mastertable$Final_superkingdom == "Bacteria" & !Mastertable$Diamond_superkingdom == "Other" & !Mastertable$Diamond_superkingdom == "Unannotated" & Mastertable$Blastn_superkingdom == "Unannotated"]
Mastertable$Final_superkingdom[!Mastertable$Final_superkingdom == "Viruses" & !Mastertable$Final_superkingdom == "Eukaryota" & !Mastertable$Final_superkingdom == "Archaea" & !Mastertable$Final_superkingdom == "Bacteria" & !Mastertable$Diamond_superkingdom == "Other" & !Mastertable$Diamond_superkingdom == "Unannotated" & !Mastertable$Blastn_superkingdom == "Unannotated" & Mastertable$CAT_superkingdom == "Unannotated"] <- Mastertable$CAT_superkingdom[!Mastertable$Final_superkingdom == "Viruses" & !Mastertable$Final_superkingdom == "Eukaryota" & !Mastertable$Final_superkingdom == "Archaea" & !Mastertable$Final_superkingdom == "Bacteria" & !Mastertable$Diamond_superkingdom == "Other" & !Mastertable$Diamond_superkingdom == "Unannotated" & !Mastertable$Blastn_superkingdom == "Unannotated" & Mastertable$CAT_superkingdom == "Unannotated"]

#View(Mastertable[Mastertable$Final_superkingdom == "Unannotated",])
#View(Mastertable[Mastertable$Final_superkingdom == "Unannotated" & Mastertable$Diamond_superkingdom == "Unannotated",])
#View(Mastertable[Mastertable$Final_superkingdom == "Unannotated" & Mastertable$Blastn_superkingdom == "Unannotated",])
#View(Mastertable[Mastertable$Final_superkingdom == "Unannotated" & Mastertable$CAT_superkingdom == "Unannotated",])
####################################
# 12.8 Final classification of dark matter
####################################
# The only sequences left are the "unnanotated" sequences, or "dark matter"
Mastertable$Final_superkingdom[Mastertable$Final_superkingdom == "1"] <- "Dark matter"
table(Mastertable$Final_superkingdom)
table(Mastertable$Final_viral)
####################################
# 13. Additional Cenote-Taker2 Taxonomical information: 1 hallmark gene
####################################
# Only add taxonomies if hallmark gene is present
Mastertable$Final_phylum[Mastertable$Final_viral == "phage" & Mastertable$Final_phylum == "Unannotated" & Mastertable$`Hallmark proteins` >= 1 & Mastertable$Cenote_viral == "cenote_prokaryotic_virus"] <- Mastertable$Cenote_phylum[Mastertable$Final_viral == "phage" & Mastertable$Final_phylum == "Unannotated" & Mastertable$`Hallmark proteins` >= 1 & Mastertable$Cenote_viral == "cenote_prokaryotic_virus"]
Mastertable$Final_class[Mastertable$Final_viral == "phage" & Mastertable$Final_class == "Unannotated" & Mastertable$`Hallmark proteins` >= 1 & Mastertable$Cenote_viral == "cenote_prokaryotic_virus"] <- Mastertable$Cenote_class[Mastertable$Final_viral == "phage" & Mastertable$Final_class == "Unannotated" & Mastertable$`Hallmark proteins` >= 1 & Mastertable$Cenote_viral == "cenote_prokaryotic_virus"]
Mastertable$Final_order[Mastertable$Final_viral == "phage" & Mastertable$Final_order == "Unannotated" & Mastertable$`Hallmark proteins` >= 1 & Mastertable$Cenote_viral == "cenote_prokaryotic_virus"] <- Mastertable$Cenote_order[Mastertable$Final_viral == "phage" & Mastertable$Final_order == "Unannotated" & Mastertable$`Hallmark proteins` >= 1 & Mastertable$Cenote_viral == "cenote_prokaryotic_virus"]
Mastertable$Final_family[Mastertable$Final_viral == "phage" & Mastertable$Final_family == "Unannotated" & Mastertable$`Hallmark proteins` >= 1 & Mastertable$Cenote_viral == "cenote_prokaryotic_virus"] <- Mastertable$Cenote_family[Mastertable$Final_viral == "phage" & Mastertable$Final_family == "Unannotated" & Mastertable$`Hallmark proteins` >= 1 & Mastertable$Cenote_viral == "cenote_prokaryotic_virus"]
# View(Mastertable[Mastertable$Final_class == "Malgrandaviricetes",])

table(Mastertable$Final_viral)
table(Mastertable$Final_superkingdom)
table(Mastertable$Final_phylum)
table(Mastertable$Final_class)
table(Mastertable$Final_order)
table(Mastertable$Final_family)
####################################
# 14. Curate Mastertable
####################################
# 14.1 Remove NA's
####################################
Mastertable[is.na(Mastertable)] <- 0
Mastertable$Totalnumberofreads <- NULL
names <- colnames(Mastertable[1:304])
names

Mastertable$Totalnumberofreads <- as.numeric(rowSums(Mastertable[,colnames(Mastertable) %in% names]))
Mastertable <- Mastertable[!(Mastertable$Totalnumberofreads == "0"),] 
sum(Mastertable$Totalnumberofreads)
####################################
# 14.2 Identify and annotated Inoviridae and Picobirnaviridae
####################################
# If 2 out of 3 homology based tools give are identical on family level than identifiy the contig as phage & annotate as Diamond, Blastn, CAT

# Inoviridae
Mastertable$Phages[Mastertable$completeness_quality == "Not-determined" & Mastertable$length > 3 & (Mastertable$Blastn_family == "Inoviridae" & Mastertable$Diamond_family == "Inoviridae" | Mastertable$Blastn_family == "Inoviridae" & Mastertable$CAT_family == "Inoviridae" | Mastertable$CAT_family == "Inoviridae" & Mastertable$Diamond_family == "Inoviridae" | Mastertable$Blastn_family == "Inoviridae" & Mastertable$Cenote_family == "Inoviridae" | Mastertable$Diamond_family == "Inoviridae" & Mastertable$Cenote_family == "Inoviridae" | Mastertable$CAT_family == "Inoviridae" & Mastertable$Cenote_family == "Inoviridae")] <- "phage"
Mastertable$Final_viral[Mastertable$completeness_quality == "Not-determined"  & Mastertable$length > 3 & (Mastertable$Blastn_family == "Inoviridae" & Mastertable$Diamond_family == "Inoviridae" | Mastertable$Blastn_family == "Inoviridae" & Mastertable$CAT_family == "Inoviridae" | Mastertable$CAT_family == "Inoviridae" & Mastertable$Diamond_family == "Inoviridae" | Mastertable$Blastn_family == "Inoviridae" & Mastertable$Cenote_family == "Inoviridae" | Mastertable$Diamond_family == "Inoviridae" & Mastertable$Cenote_family == "Inoviridae" | Mastertable$CAT_family == "Inoviridae" & Mastertable$Cenote_family == "Inoviridae")] <- "phage"
Mastertable$Final_superkingdom[Mastertable$completeness_quality == "Not-determined"  & Mastertable$length > 3 & (Mastertable$Blastn_family == "Inoviridae" & Mastertable$Diamond_family == "Inoviridae" | Mastertable$Blastn_family == "Inoviridae" & Mastertable$CAT_family == "Inoviridae" | Mastertable$CAT_family == "Inoviridae" & Mastertable$Diamond_family == "Inoviridae" | Mastertable$Blastn_family == "Inoviridae" & Mastertable$Cenote_family == "Inoviridae" | Mastertable$Diamond_family == "Inoviridae" & Mastertable$Cenote_family == "Inoviridae" | Mastertable$CAT_family == "Inoviridae" & Mastertable$Cenote_family == "Inoviridae")] <- "Viruses"
Mastertable$Final_phylum[Mastertable$completeness_quality == "Not-determined"  & Mastertable$length > 3 & (Mastertable$Blastn_family == "Inoviridae" & Mastertable$Diamond_family == "Inoviridae" | Mastertable$Blastn_family == "Inoviridae" & Mastertable$CAT_family == "Inoviridae" | Mastertable$CAT_family == "Inoviridae" & Mastertable$Diamond_family == "Inoviridae" | Mastertable$Blastn_family == "Inoviridae" & Mastertable$Cenote_family == "Inoviridae" | Mastertable$Diamond_family == "Inoviridae" & Mastertable$Cenote_family == "Inoviridae" | Mastertable$CAT_family == "Inoviridae" & Mastertable$Cenote_family == "Inoviridae")] <- "Hofneiviricota"
Mastertable$Final_class[Mastertable$completeness_quality == "Not-determined"  & Mastertable$length > 3 & (Mastertable$Blastn_family == "Inoviridae" & Mastertable$Diamond_family == "Inoviridae" | Mastertable$Blastn_family == "Inoviridae" & Mastertable$CAT_family == "Inoviridae" | Mastertable$CAT_family == "Inoviridae" & Mastertable$Diamond_family == "Inoviridae" | Mastertable$Blastn_family == "Inoviridae" & Mastertable$Cenote_family == "Inoviridae" | Mastertable$Diamond_family == "Inoviridae" & Mastertable$Cenote_family == "Inoviridae" | Mastertable$CAT_family == "Inoviridae" & Mastertable$Cenote_family == "Inoviridae")] <- "Faserviricetes"
Mastertable$Final_order[Mastertable$completeness_quality == "Not-determined"  & Mastertable$length > 3 & (Mastertable$Blastn_family == "Inoviridae" & Mastertable$Diamond_family == "Inoviridae" | Mastertable$Blastn_family == "Inoviridae" & Mastertable$CAT_family == "Inoviridae" | Mastertable$CAT_family == "Inoviridae" & Mastertable$Diamond_family == "Inoviridae" | Mastertable$Blastn_family == "Inoviridae" & Mastertable$Cenote_family == "Inoviridae" | Mastertable$Diamond_family == "Inoviridae" & Mastertable$Cenote_family == "Inoviridae" | Mastertable$CAT_family == "Inoviridae" & Mastertable$Cenote_family == "Inoviridae")] <- "Tubulavirales"
Mastertable$Final_family[Mastertable$completeness_quality == "Not-determined"  & Mastertable$length > 3 & (Mastertable$Blastn_family == "Inoviridae" & Mastertable$Diamond_family == "Inoviridae" | Mastertable$Blastn_family == "Inoviridae" & Mastertable$CAT_family == "Inoviridae" | Mastertable$CAT_family == "Inoviridae" & Mastertable$Diamond_family == "Inoviridae" | Mastertable$Blastn_family == "Inoviridae" & Mastertable$Cenote_family == "Inoviridae" | Mastertable$Diamond_family == "Inoviridae" & Mastertable$Cenote_family == "Inoviridae" | Mastertable$CAT_family == "Inoviridae" & Mastertable$Cenote_family == "Inoviridae")] <- "Inoviridae"
nrow(Mastertable[Mastertable$Final_family=="Inoviridae",])
sum(Mastertable$Totalnumberofreads[Mastertable$Final_family=="Inoviridae"])
# 54 NR-contigs to which 118million PE reads are mapped 

# 13.2.3 Picobirnaviridae
Mastertable$Phages[Mastertable$completeness_quality == "Not-determined" & (Mastertable$Blastn_family == "Picobirnaviridae" & Mastertable$Diamond_family == "Picobirnaviridae" | Mastertable$Blastn_family == "Picobirnaviridae" & Mastertable$CAT_family == "Picobirnaviridae" | Mastertable$CAT_family == "Picobirnaviridae" & Mastertable$Diamond_family == "Picobirnaviridae")] <- "phage"
Mastertable$Final_viral[Mastertable$completeness_quality == "Not-determined" & (Mastertable$Blastn_family == "Picobirnaviridae" & Mastertable$Diamond_family == "Picobirnaviridae" | Mastertable$Blastn_family == "Picobirnaviridae" & Mastertable$CAT_family == "Picobirnaviridae" | Mastertable$CAT_family == "Picobirnaviridae" & Mastertable$Diamond_family == "Picobirnaviridae")] <- "phage"
Mastertable$Final_superkingdom[Mastertable$completeness_quality == "Not-determined" & (Mastertable$Blastn_family == "Picobirnaviridae" & Mastertable$Diamond_family == "Picobirnaviridae" | Mastertable$Blastn_family == "Picobirnaviridae" & Mastertable$CAT_family == "Picobirnaviridae" | Mastertable$CAT_family == "Picobirnaviridae" & Mastertable$Diamond_family == "Picobirnaviridae")] <- "Viruses"
Mastertable$Final_phylum[Mastertable$completeness_quality == "Not-determined" & (Mastertable$Blastn_family == "Picobirnaviridae" & Mastertable$Diamond_family == "Picobirnaviridae" | Mastertable$Blastn_family == "Picobirnaviridae" & Mastertable$CAT_family == "Picobirnaviridae" | Mastertable$CAT_family == "Picobirnaviridae" & Mastertable$Diamond_family == "Picobirnaviridae")] <- "Pisuviricota"
Mastertable$Final_class[Mastertable$completeness_quality == "Not-determined" & (Mastertable$Blastn_family == "Picobirnaviridae" & Mastertable$Diamond_family == "Picobirnaviridae" | Mastertable$Blastn_family == "Picobirnaviridae" & Mastertable$CAT_family == "Picobirnaviridae" | Mastertable$CAT_family == "Picobirnaviridae" & Mastertable$Diamond_family == "Picobirnaviridae")] <- "Duplopiviricetes"
Mastertable$Final_order[Mastertable$completeness_quality == "Not-determined" & (Mastertable$Blastn_family == "Picobirnaviridae" & Mastertable$Diamond_family == "Picobirnaviridae" | Mastertable$Blastn_family == "Picobirnaviridae" & Mastertable$CAT_family == "Picobirnaviridae" | Mastertable$CAT_family == "Picobirnaviridae" & Mastertable$Diamond_family == "Picobirnaviridae")] <- "Durnavirales"
Mastertable$Final_family[Mastertable$completeness_quality == "Not-determined" & (Mastertable$Blastn_family == "Picobirnaviridae" & Mastertable$Diamond_family == "Picobirnaviridae" | Mastertable$Blastn_family == "Picobirnaviridae" & Mastertable$CAT_family == "Picobirnaviridae" | Mastertable$CAT_family == "Picobirnaviridae" & Mastertable$Diamond_family == "Picobirnaviridae")] <- "Picobirnaviridae"
nrow(Mastertable[Mastertable$Final_family=="Picobirnaviridae",])
sum(Mastertable$Totalnumberofreads[Mastertable$Final_family=="Picobirnaviridae"])
# 69 NR-contigs to whuch 2,2M PE reads are mapped

table(Mastertable$Phages)
table(Mastertable$Final_viral) # prokaryotic versus eukaryotic virus
table(Mastertable$Final_superkingdom[Mastertable$Final_viral=="phage"]) 
table(Mastertable$Final_phylum[Mastertable$Final_viral=="phage"])
table(Mastertable$Final_class[Mastertable$Final_viral=="phage"])
table(Mastertable$Final_family[Mastertable$Final_viral=="phage"])
table(Mastertable$Final_subfamily[Mastertable$Final_viral=="phage"])
table(Mastertable$Final_genus[Mastertable$Final_viral=="phage"])

# Family
(sum(Mastertable$Totalnumberofreads[Mastertable$Final_viral=="phage" & Mastertable$Final_family=="Unannotated"])/sum(Mastertable$Totalnumberofreads[Mastertable$Final_viral=="phage"]))*100 # 90.7% unclassified reads
sum(Mastertable$Totalnumberofreads[Mastertable$Final_family=="Microviridae"])/sum(Mastertable$Totalnumberofreads) # 28.2% are Microviridae
(1517/2519)*100 # 60.2% unclassified contigs
## Subfamily: 47.6% unclassified genus reads, 60.2% contigs

# Subfamily
(sum(Mastertable$Totalnumberofreads[Mastertable$Final_viral=="phage" & Mastertable$Fin=="Unannotated"])/sum(Mastertable$Totalnumberofreads[Mastertable$Final_viral=="phage"]))*100 # 90.7% unclassified reads
(2273/2519)*100 # 90.2% unclassified contigs
## Subfamily: 90.7% unclassified genus reads, 90.2% contigs

# Genus
(sum(Mastertable$Totalnumberofreads[Mastertable$Final_viral=="phage" & Mastertable$Final_order=="Unannotated"])/sum(Mastertable$Totalnumberofreads[Mastertable$Final_viral=="phage"]))*100 # 88.5% unclassified reads
(2284/2519)*100 # 90.7% unclassified contigs
## Genus: 88.5% unclassified genus reads, 90.7% contigs

## Families
(sum(Mastertable$Totalnumberofreads[Mastertable$Final_family=="Suoliviridae"])/sum(Mastertable$Totalnumberofreads[Mastertable$Final_viral=="phage"]))*100 
#View(Mastertable[rev(order(Mastertable$Totalnumberofreads)),])
####################################
# 14.3 Run Cenote-Taker2 for protein identification on additions
####################################
## Print Contigs non-determine group, ino, picobirnaviridae in file go to HPC and find proteins
## Determine (i) temperate, (ii) toxins, (iii) other proteins

getwd()
setwd("/Users/daan/Desktop")

# Take names contigs of Inoviridae
vector_1 <- as.data.frame(rownames(Mastertable[Mastertable$completeness_quality == "Not-determined"  & Mastertable$length > 3 & (Mastertable$Blastn_family == "Inoviridae" & Mastertable$Diamond_family == "Inoviridae" | Mastertable$Blastn_family == "Inoviridae" & Mastertable$CAT_family == "Inoviridae" | Mastertable$CAT_family == "Inoviridae" & Mastertable$Diamond_family == "Inoviridae"),]))
vector_1$names <- vector_1$`rownames(Mastertable[Mastertable$completeness_quality == "Not-determined" & Mastertable$length > 3 & (Mastertable$Blastn_family == "Inoviridae" & Mastertable$Diamond_family == "Inoviridae" | Mastertable$Blastn_family == "Inoviridae" & Mastertable$CAT_family == "Inoviridae" | Mastertable$CAT_family == "Inoviridae" & Mastertable$Diamond_family == "Inoviridae"), ])`
vector_1$`rownames(Mastertable[Mastertable$completeness_quality == "Not-determined" & Mastertable$length > 3 & (Mastertable$Blastn_family == "Inoviridae" & Mastertable$Diamond_family == "Inoviridae" | Mastertable$Blastn_family == "Inoviridae" & Mastertable$CAT_family == "Inoviridae" | Mastertable$CAT_family == "Inoviridae" & Mastertable$Diamond_family == "Inoviridae"), ])` <- NULL

# Take names contigs of Picobirnaviridae
vector_2 <- as.data.frame(rownames(Mastertable[Mastertable$completeness_quality == "Not-determined" & (Mastertable$Blastn_family == "Picobirnaviridae" & Mastertable$Diamond_family == "Picobirnaviridae" | Mastertable$Blastn_family == "Picobirnaviridae" & Mastertable$CAT_family == "Picobirnaviridae" | Mastertable$CAT_family == "Picobirnaviridae" & Mastertable$Diamond_family == "Picobirnaviridae"),]))
vector_2$names <- vector_2$`rownames(Mastertable[Mastertable$completeness_quality == "Not-determined" & (Mastertable$Blastn_family == "Picobirnaviridae" & Mastertable$Diamond_family == "Picobirnaviridae" | Mastertable$Blastn_family == "Picobirnaviridae" & Mastertable$CAT_family == "Picobirnaviridae" | Mastertable$CAT_family == "Picobirnaviridae" & Mastertable$Diamond_family == "Picobirnaviridae"), ])`
vector_2$`rownames(Mastertable[Mastertable$completeness_quality == "Not-determined" & (Mastertable$Blastn_family == "Picobirnaviridae" & Mastertable$Diamond_family == "Picobirnaviridae" | Mastertable$Blastn_family == "Picobirnaviridae" & Mastertable$CAT_family == "Picobirnaviridae" | Mastertable$CAT_family == "Picobirnaviridae" & Mastertable$Diamond_family == "Picobirnaviridae"), ])` <- NULL

# Merge into one DF
vector_3 <- rbind(vector_1,vector_2)
#View(vector_3)

# Print as text file
write.table(vector_3, file = "cenote.txt", sep = "\t",row.names = FALSE)

## Now run CenoteTaker2 on those phages and merge the protein, temperateand toxin data with existing files, and run this entire script again.
####################################
# 14.4 Manual curation of phage fracion: separation with CrAss
####################################
table(Mastertable$Final_class[Mastertable$Final_viral=="phage"])
Mastertable$Final_order[Mastertable$Final_class=="Caudoviricetes" & Mastertable$Final_order=="Unannotated"] <- "Unclassified Caudoviricetes"
table(Mastertable$Final_order[Mastertable$Final_viral=="phage"])
Mastertable1 <- Mastertable
####################################
# 14. Curation viral fraction
####################################
## In case of presence of eukaryotic viruses that are clearly wrongly annotated you have to manually change them to what the correct classification should be. For example, giant viruses can't be present (filter in wetlab limits this)
####################################
# 14.1 Curate eukaryotic viruses
####################################
Mastertable2 <- Mastertable

# Add AS
Mastertable2$Blastn_AS <- ((Mastertable2$Blastn_ANI/100)*(Mastertable2$Blastn_coverage_relative/100))
Mastertable2$Blastn_AS[!Mastertable2$Blastn_viral == "Blastn_eukaryotic_virus"] <- 0

Mastertable2$Diamond_AS <- ((Mastertable2$Diamond_ANI/100)*(Mastertable2$Diamond_coverage_relative/100))
Mastertable2$Diamond_AS[!Mastertable2$Diamond_viral == "Diamond_eukaryotic_virus"] <- 0

# Annotation tools to use
Mastertable2$Best_AS[Mastertable2$Diamond_AS > Mastertable2$Blastn_AS] <- "Diamond"
Mastertable2$Best_AS[Mastertable2$Diamond_AS < Mastertable2$Blastn_AS] <- "Blastn"
Mastertable2$Best_AS[Mastertable2$Diamond_AS < 0.1 & Mastertable2$Blastn_AS < 0.1] <- "Remove"

# Re-annotate: Blastn
Mastertable2$Final_phylum[Mastertable2$Best_AS == "Blastn" & Mastertable2$Blastn_viral == "Blastn_eukaryotic_virus"] <- Mastertable2$Blastn_phylum[Mastertable2$Best_AS == "Blastn" & Mastertable2$Blastn_viral == "Blastn_eukaryotic_virus"]
Mastertable2$Final_class[Mastertable2$Best_AS == "Blastn" & Mastertable2$Blastn_viral == "Blastn_eukaryotic_virus"] <- Mastertable2$Blastn_class[Mastertable2$Best_AS == "Blastn" & Mastertable2$Blastn_viral == "Blastn_eukaryotic_virus"]
Mastertable2$Final_order[Mastertable2$Best_AS == "Blastn" & Mastertable2$Blastn_viral == "Blastn_eukaryotic_virus"] <- Mastertable2$Blastn_order[Mastertable2$Best_AS == "Blastn" & Mastertable2$Blastn_viral == "Blastn_eukaryotic_virus"]
Mastertable2$Final_family[Mastertable2$Best_AS == "Blastn" & Mastertable2$Blastn_viral == "Blastn_eukaryotic_virus"] <- Mastertable2$Blastn_family[Mastertable2$Best_AS == "Blastn" & Mastertable2$Blastn_viral == "Blastn_eukaryotic_virus"]
Mastertable2$Final_subfamily[Mastertable2$Best_AS == "Blastn" & Mastertable2$Blastn_viral == "Blastn_eukaryotic_virus"] <- Mastertable2$Blastn_subfamily[Mastertable2$Best_AS == "Blastn" & Mastertable2$Blastn_viral == "Blastn_eukaryotic_virus"]
Mastertable2$Final_genus[Mastertable2$Best_AS == "Blastn" & Mastertable2$Blastn_viral == "Blastn_eukaryotic_virus"] <- Mastertable2$Blastn_genus[Mastertable2$Best_AS == "Blastn" & Mastertable2$Blastn_viral == "Blastn_eukaryotic_virus"]
Mastertable2$Final_species[Mastertable2$Best_AS == "Blastn" & Mastertable2$Blastn_viral == "Blastn_eukaryotic_virus"] <- Mastertable2$Blastn_species[Mastertable2$Best_AS == "Blastn" & Mastertable2$Blastn_viral == "Blastn_eukaryotic_virus"]
Mastertable2$Final_ANI[Mastertable2$Best_AS == "Blastn" & Mastertable2$Blastn_viral == "Blastn_eukaryotic_virus"] <- Mastertable2$Blastn_ANI[Mastertable2$Best_AS == "Blastn" & Mastertable2$Blastn_viral == "Blastn_eukaryotic_virus"]
Mastertable2$Final_coverage[Mastertable2$Best_AS == "Blastn"& Mastertable2$Blastn_viral == "Blastn_eukaryotic_virus"] <- Mastertable2$Blastn_coverage_relative[Mastertable2$Best_AS == "Blastn" & Mastertable2$Blastn_viral == "Blastn_eukaryotic_virus"]

# Re-annotate: Diamond
Mastertable2$Final_phylum[Mastertable2$Best_AS == "Diamond" & Mastertable2$Diamond_viral == "Diamond_eukaryotic_virus"] <- Mastertable2$Diamond_phylum[Mastertable2$Best_AS == "Diamond" & Mastertable2$Diamond_viral == "Diamond_eukaryotic_virus"] 
Mastertable2$Final_class[Mastertable2$Best_AS == "Diamond" & Mastertable2$Diamond_viral == "Diamond_eukaryotic_virus"] <- Mastertable2$Diamond_class[Mastertable2$Best_AS == "Diamond" & Mastertable2$Diamond_viral == "Diamond_eukaryotic_virus"] 
Mastertable2$Final_order[Mastertable2$Best_AS == "Diamond" & Mastertable2$Diamond_viral == "Diamond_eukaryotic_virus"] <- Mastertable2$Diamond_order[Mastertable2$Best_AS == "Diamond" & Mastertable2$Diamond_viral == "Diamond_eukaryotic_virus"] 
Mastertable2$Final_family[Mastertable2$Best_AS == "Diamond" & Mastertable2$Diamond_viral == "Diamond_eukaryotic_virus"] <- Mastertable2$Diamond_family[Mastertable2$Best_AS == "Diamond" & Mastertable2$Diamond_viral == "Diamond_eukaryotic_virus"] 
Mastertable2$Final_subfamily[Mastertable2$Best_AS == "Diamond"& Mastertable2$Diamond_viral == "Diamond_eukaryotic_virus"] <- Mastertable2$Diamond_subfamily[Mastertable2$Best_AS == "Diamond" & Mastertable2$Diamond_viral == "Diamond_eukaryotic_virus"] 
Mastertable2$Final_genus[Mastertable2$Best_AS == "Diamond" & Mastertable2$Diamond_viral == "Diamond_eukaryotic_virus"] <- Mastertable2$Diamond_genus[Mastertable2$Best_AS == "Diamond" & Mastertable2$Diamond_viral == "Diamond_eukaryotic_virus"] 
Mastertable2$Final_species[Mastertable2$Best_AS == "Diamond" & Mastertable2$Diamond_viral == "Diamond_eukaryotic_virus"] <- Mastertable2$Diamond_species[Mastertable2$Best_AS == "Diamond" & Mastertable2$Diamond_viral == "Diamond_eukaryotic_virus"] 
Mastertable2$Final_ANI[Mastertable2$Best_AS == "Diamond" & Mastertable2$Diamond_viral == "Diamond_eukaryotic_virus"] <- Mastertable2$Diamond_ANI[Mastertable2$Best_AS == "Diamond" & Mastertable2$Diamond_viral == "Diamond_eukaryotic_virus"]
Mastertable2$Final_coverage[Mastertable2$Best_AS == "Diamond" & Mastertable2$Diamond_viral == "Diamond_eukaryotic_virus"] <- Mastertable2$Diamond_coverage_relative[Mastertable2$Best_AS == "Diamond" & Mastertable2$Diamond_viral == "Diamond_eukaryotic_virus"]

# Remove bad alignments from eukaryotic viral fraction
Mastertable2$Final_viral[Mastertable2$Best_AS == "Remove" & Mastertable2$Final_viral == "eukaryotic virus"] <- 1
Mastertable2$Final_superkingdom[Mastertable2$Best_AS == "Remove" & Mastertable2$Final_viral == "eukaryotic virus"] <- "Unannotated"

# Reads and contigs
(sum(Mastertable2$Totalnumberofreads[Mastertable2$Final_viral == "eukaryotic virus"])/sum(Mastertable2$Totalnumberofreads))*100
nrow(Mastertable2[Mastertable2$Final_viral == "eukaryotic virus",])
table(Mastertable2$Final_family[Mastertable2$Final_viral=="eukaryotic virus"])
table(Mastertable2$Final_viral)
table(Mastertable2$Final_superkingdom)
#View(Mastertable2[Mastertable2$Final_viral=="eukaryotic virus",])

Mastertable <- Mastertable2
####################################
# 14.1.2 Subset eukaryotic viral fraction
####################################
vector_1 <- which(colnames(Mastertable) == "Final_ANI")
vector_2 <- which(colnames(Mastertable) == "Final_coverage")
vector_3 <- which(colnames(Mastertable) == "Final_superkingdom")

vector_ANI_Blastn <- which(colnames(Mastertable) == "Blastn_species")
vector_cov_Blastn <- which(colnames(Mastertable) == "Blastn_coverage_relative")

vector_ANI_Diamond <- which(colnames(Mastertable) == "Diamond_species")
vector_cov_Diamond <- which(colnames(Mastertable) == "Diamond_coverage_relative")

length <- which(colnames(Mastertable) == "length")
reads <- which(colnames(Mastertable) == "Totalnumberofreads")

eukaryotic_viruses <- Mastertable[Mastertable$Final_viral == "eukaryotic virus", c(vector_3:vector_2, vector_ANI_Blastn:vector_cov_Blastn, vector_ANI_Diamond: vector_cov_Diamond, length, reads)]

setwd("/Users/daan/Desktop/Bioinformatics/Analysis/FMT/Data/1_Eukaryotic_viruses")
getwd()

#install.packages("writexl")
library("writexl")

write_xlsx(eukaryotic_viruses,"./eukaryotic_viruses.xlsx")
####################################
# 14.1.4 Clean-up environment
####################################
rm(vector_1)
rm(vector_2)
rm(eukaryotic_viruses)
#################################
# 15 Add Nayfach after running again with Inoviridae
#################################
# Take clusters

table(Mastertable$Final_viral)
Mastertable2 <- Mastertable[!c(Mastertable$Genus==0 & Mastertable$Final_viral=="phage"),]
(sum(Mastertable$Totalnumberofreads[Mastertable$Final_viral=="phage"])/sum(Mastertable$Totalnumberofreads))*100
(sum(Mastertable2$Totalnumberofreads[Mastertable2$Final_viral=="phage"])/sum(Mastertable2$Totalnumberofreads))*100
Mastertable <- Mastertable2
#View(Mastertable[Mastertable$Final_viral=="phage",])

DF_1 <- Mastertable[Mastertable$Final_viral=="phage",c(374:378,380,381,385,386,387,324:326)]
colnames(DF_1)

ncol(Mastertable)
Mastertable$names <- rownames(Mastertable)
DF_1$Final_ANI[DF_1$Final_ANI=="Unannotated"] <- 0
DF_1$Final_coverage[DF_1$Final_coverage=="Unannotated"] <- 0
DF_1$Final_ANI <- as.numeric(DF_1$Final_ANI)
DF_1$Final_ANI <- DF_1$Final_ANI/100
DF_1$Final_coverage <- as.numeric(DF_1$Final_coverage)
DF_1$Final_coverage <- DF_1$Final_coverage/100
DF_1$AS <- DF_1$Final_ANI*DF_1$Final_coverage
DF_1$AS[DF_1$AS<0.1] <- "Unannotated"
DF_1
View(DF_1)

# export and assign taxonomies to all aligning nodes
setwd("/Users/daan/Desktop/Bioinformatics/Analysis/FMT/Data/3_Phages")
getwd()
dir()

library(writexl)
write_xlsx(DF_1,"./DF_1.xlsx")

# insert table again and merge with Mastertable NODES
Mastertable_genus <- Mastertable

# Import table
dir()
DF_1 <- as.data.frame(read_excel("DF_1.xlsx"))
rownames(DF_1) <- DF_1$names
DF_1$names <- NULL

# Merge by rownames
Mastertable_genus1 <- merge(Mastertable_genus, DF_1, by = 0, all = T)
rownames(Mastertable_genus1) <- Mastertable_genus1$Row.names
Mastertable_genus1$Row.names <- NULL

# paste new taxonomies
Mastertable_genus1$names <- NULL
Mastertable_genus1$cluster.y <- NULL

Mastertable_genus1$Final_class.y[Mastertable_genus1$Final_viral=="eukaryotic virus"] <- Mastertable_genus1$Final_class.x [Mastertable_genus1$Final_viral=="eukaryotic virus"]
Mastertable_genus1$Final_order.y[Mastertable_genus1$Final_viral=="eukaryotic virus"] <- Mastertable_genus1$Final_order.x [Mastertable_genus1$Final_viral=="eukaryotic virus"]
Mastertable_genus1$Final_family.y[Mastertable_genus1$Final_viral=="eukaryotic virus"] <- Mastertable_genus1$Final_family.x [Mastertable_genus1$Final_viral=="eukaryotic virus"]
Mastertable_genus1$Final_subfamily.y[Mastertable_genus1$Final_viral=="eukaryotic virus"] <- Mastertable_genus1$Final_subfamily.x [Mastertable_genus1$Final_viral=="eukaryotic virus"]
Mastertable_genus1$Final_genus.y[Mastertable_genus1$Final_viral=="eukaryotic virus"] <- Mastertable_genus1$Final_genus.x [Mastertable_genus1$Final_viral=="eukaryotic virus"]

Mastertable_genus1$Final_class.x <- Mastertable_genus1$Final_class.y
Mastertable_genus1$Final_order.x <- Mastertable_genus1$Final_order.y
Mastertable_genus1$Final_family.x <- Mastertable_genus1$Final_family.y
Mastertable_genus1$Final_subfamily.x <- Mastertable_genus1$Final_subfamily.y
Mastertable_genus1$Final_genus.x <- Mastertable_genus1$Final_genus.y
Mastertable_genus1$Final_class.y <- NULL
Mastertable_genus1$Final_order.y <- NULL
Mastertable_genus1$Final_family.y <- NULL
Mastertable_genus1$Final_subfamily.y <- NULL
Mastertable_genus1$Final_genus.y <- NULL
Mastertable_genus1$cluster.x <- NULL

colnames(Mastertable_genus1)[colnames(Mastertable_genus1) == "Final_class.x"] <- "Final_class"
colnames(Mastertable_genus1)[colnames(Mastertable_genus1) == "Final_order.x"] <- "Final_order"
colnames(Mastertable_genus1)[colnames(Mastertable_genus1) == "Final_family.x"] <- "Final_family"
colnames(Mastertable_genus1)[colnames(Mastertable_genus1) == "Final_subfamily.x"] <- "Final_subfamily"
colnames(Mastertable_genus1)[colnames(Mastertable_genus1) == "Final_genus.x"] <- "Final_genus"

# Subset viral table (eukaryotic viral & phages)
Mastertable_eukaryotic <- Mastertable_genus1[Mastertable_genus1$Final_viral=="eukaryotic virus",]
Mastertable_phage <- Mastertable_genus1[Mastertable_genus1$Final_viral=="phage",]

# Count reads
table(Mastertable_phage$Final_family)

(sum(Mastertable_phage$Totalnumberofreads[Mastertable_phage$Final_family=="Intestiviridae"])/sum(Mastertable_phage$Totalnumberofreads))*100
(sum(Mastertable_phage$Totalnumberofreads[Mastertable_phage$Final_family=="Suoliviridae"])/sum(Mastertable_phage$Totalnumberofreads))*100
(sum(Mastertable_phage$Totalnumberofreads[Mastertable_phage$Final_family=="Steigviridae"])/sum(Mastertable_phage$Totalnumberofreads))*100
(sum(Mastertable_phage$Totalnumberofreads[Mastertable_phage$Final_family=="Microviridae"])/sum(Mastertable_phage$Totalnumberofreads))*100
(sum(Mastertable_phage$Totalnumberofreads[Mastertable_phage$Final_family=="Unannotated"])/sum(Mastertable_phage$Totalnumberofreads))*100
#  View(Mastertable_genus1)

table(Mastertable_phage$Final_order)
(sum(Mastertable_phage$Totalnumberofreads[Mastertable_phage$Final_order=="Crassvirales"])/sum(Mastertable_phage$Totalnumberofreads))*100
(sum(Mastertable_phage$Totalnumberofreads[Mastertable_phage$Final_order=="Tubulavirales"])/sum(Mastertable_phage$Totalnumberofreads))*100
(sum(Mastertable_phage$Totalnumberofreads[Mastertable_phage$Final_order=="Petitvirales"])/sum(Mastertable_phage$Totalnumberofreads))*100
(sum(Mastertable_phage$Totalnumberofreads[Mastertable_phage$Final_order=="Unclassified Caudoviricetes"])/sum(Mastertable_phage$Totalnumberofreads))*100

# Conclusion: mostly microviridae (43%) and unannotated (36%)
#################################
# Tables for further analysis
#################################
Mastertable # Reads & contigs figure
Mastertable_genus1 # Percentage classified on class, order, family level
Mastertable_eukaryotic # Eukaryotic viral figure
Mastertable_phage  # Phage analysis
#################################

######### !!!!!!!!!!!!!!! #########
    # SAVE Global Environment #
######### !!!!!!!!!!!!!!! #########
