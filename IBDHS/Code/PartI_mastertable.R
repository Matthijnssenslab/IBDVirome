####################################
# SCRIPT 1: CREATE A MASTERTABLE
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
## We will use "taxonomizr" for adding taxonomies later (10. Taxonomy) by converting single TaxIDs from the tab files to taxon names based on most recently downloaded taxa.

setwd("/Volumes/Daan_2/Desktop/DB/")
getwd()
dir()

prepareDatabase('accessionTaxa.sql')
####################################
# 2. Create abundance table 
####################################
## A) set directory
setwd("/Users/daan/Desktop/Transfer/Input_R/Abundances/")
getwd()
dir()

## B) read the abundance text files and store in a list
list <- list()
list_txt <- dir(pattern = "*.length.corrected.coverages.txt", full.names = FALSE)
for (k in 1: length(list_txt)){
  list[[k]] <- read.delim(list_txt[k])
}

## C) create names file
names <- sub(pattern = "*.length.corrected.coverages.txt", "", list_txt)

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

#View(table)
####################################
# 3. Create vector of contigs
####################################
setwd("/Users/daan/Desktop/Transfer/Input_R/Abundances/Multi_fasta_files")
getwd()
dir()

non_redundant_contig_names <- read.table("NR_contigs_names.txt", quote = "\"", comment.char = "")
colnames(non_redundant_contig_names) <- "contigs_orginal_name.txt"
rownames(non_redundant_contig_names) <- non_redundant_contig_names$contigs_orginal_name
rownames(non_redundant_contig_names) <- gsub("\\.{1}", "_", rownames(non_redundant_contig_names))
contigs_project <- non_redundant_contig_names$contigs_orginal_name 
####################################
# 4. Insert taxonomical annotation 
####################################
# 4.1 Diamond TaxIdpath
####################################
setwd("/Users/daan/Desktop/Transfer/Input_R/Taxonomy/")
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
#####################################
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
####################################
# 5 Insert taxonomical annotation 
####################################
# 5.5.1 Number of hallmark proteins
####################################
setwd ("/Users/daan/Desktop/Transfer/Input_R/Functions") 
getwd()
dir()

AllSamples_Cenotetaker2_hallmark_proteins <- as.data.frame(read_delim("summary_total_viral_proteins.tsv", "\t", escape_double = FALSE, trim_ws = TRUE, col_names = FALSE))
AllSamples_Cenotetaker2_hallmark_proteins <- AllSamples_Cenotetaker2_hallmark_proteins[,c(1,3)]
colnames(AllSamples_Cenotetaker2_hallmark_proteins) <- c('contig', 'Hallmark proteins')
rownames(AllSamples_Cenotetaker2_hallmark_proteins) <- AllSamples_Cenotetaker2_hallmark_proteins$contig
AllSamples_Cenotetaker2_hallmark_proteins$contig <- NULL
rownames(AllSamples_Cenotetaker2_hallmark_proteins) <- gsub("\\.{1}", "_", rownames(AllSamples_Cenotetaker2_hallmark_proteins))
AllSamples_Cenotetaker2_hallmark_proteins[is.na(AllSamples_Cenotetaker2_hallmark_proteins)] <- 0
####################################
# 5.5.2 Lysogenic contigs
####################################
setwd ("/Users/daan/Desktop/Transfer/Input_R/Functions/input") 

# re-evaluate if I don't miss any "words" in orginal file
AllSamples_Cenotetaker2_lysogenic_phages <- as.data.frame(read_delim("All_predicted_protein_lysogenic_nodes.txt", "\t", escape_double = FALSE, trim_ws = TRUE, col_names = FALSE))
colnames(AllSamples_Cenotetaker2_lysogenic_phages) <- c('contig', 'lysogenic cycle')
rownames(AllSamples_Cenotetaker2_lysogenic_phages) <- AllSamples_Cenotetaker2_lysogenic_phages$contig
AllSamples_Cenotetaker2_lysogenic_phages$contig <- NULL
rownames(AllSamples_Cenotetaker2_lysogenic_phages) <- gsub("\\.{1}", "_", rownames(AllSamples_Cenotetaker2_lysogenic_phages))
AllSamples_Cenotetaker2_lysogenic_phages[is.na(AllSamples_Cenotetaker2_lysogenic_phages)] <- 0
####################################
# 5.5.3 Toxin-containing contigs
####################################
AllSamples_Cenotetaker2_toxin_phages <- as.data.frame(read_delim("All_predicted_toxin_phages_nodes.txt", "\t", escape_double = TRUE, trim_ws = TRUE, col_names = FALSE))
colnames(AllSamples_Cenotetaker2_toxin_phages) <- c('contig', 'Toxin presence', 'Toxin type')
rownames(AllSamples_Cenotetaker2_toxin_phages) <- AllSamples_Cenotetaker2_toxin_phages$contig
AllSamples_Cenotetaker2_toxin_phages$contig <- NULL
rownames(AllSamples_Cenotetaker2_toxin_phages) <- gsub("\\.{1}", "_", rownames(AllSamples_Cenotetaker2_toxin_phages))
AllSamples_Cenotetaker2_toxin_phages[is.na(AllSamples_Cenotetaker2_toxin_phages)] <- 0
####################################
# 5.5.4 Recombination-containing contigs
####################################
AllSamples_Cenotetaker2_recombination_phages <- as.data.frame(read_delim("All_predicted_recombination_phages_nodes.txt", "\t", escape_double = FALSE, trim_ws = TRUE, col_names = FALSE))
colnames(AllSamples_Cenotetaker2_recombination_phages) <- c('contig', 'recombination')
rownames(AllSamples_Cenotetaker2_recombination_phages) <- AllSamples_Cenotetaker2_recombination_phages$contig
AllSamples_Cenotetaker2_recombination_phages$contig <- NULL
rownames(AllSamples_Cenotetaker2_recombination_phages) <- gsub("\\.{1}", "_", rownames(AllSamples_Cenotetaker2_recombination_phages))
AllSamples_Cenotetaker2_recombination_phages[is.na(AllSamples_Cenotetaker2_recombination_phages)] <- 0
####################################
# 5.5.5 RT-containing contigs
####################################
AllSamples_Cenotetaker2_RT_phages <- as.data.frame(read_delim("RT_combined.txt", "\t", escape_double = FALSE, trim_ws = TRUE, col_names = FALSE))
colnames(AllSamples_Cenotetaker2_RT_phages) <- c('contig', 'RT', 'RT type')
rownames(AllSamples_Cenotetaker2_RT_phages) <- AllSamples_Cenotetaker2_RT_phages$contig
AllSamples_Cenotetaker2_RT_phages$contig <- NULL
rownames(AllSamples_Cenotetaker2_RT_phages) <- gsub("\\.{1}", "_", rownames(AllSamples_Cenotetaker2_RT_phages))
AllSamples_Cenotetaker2_RT_phages[is.na(AllSamples_Cenotetaker2_RT_phages)] <- 0
####################################
# 5.5.6 Bacteriocin-containing contigs
####################################
AllSamples_Cenotetaker2_Bacteriocin_phages <- as.data.frame(read_delim("Bacteriocins.txt", "\t", escape_double = FALSE, trim_ws = TRUE, col_names = FALSE))
colnames(AllSamples_Cenotetaker2_Bacteriocin_phages) <- c('contig', 'bacteriocin')
rownames(AllSamples_Cenotetaker2_Bacteriocin_phages) <- AllSamples_Cenotetaker2_Bacteriocin_phages$contig
AllSamples_Cenotetaker2_Bacteriocin_phages$contig <- NULL
rownames(AllSamples_Cenotetaker2_Bacteriocin_phages) <- gsub("\\.{1}", "_", rownames(AllSamples_Cenotetaker2_Bacteriocin_phages))
AllSamples_Cenotetaker2_Bacteriocin_phages[is.na(AllSamples_Cenotetaker2_Bacteriocin_phages)] <- 0
####################################
# 5.5.7 Taxonomy Cenotetaker2
####################################
setwd("/Users/daan/Desktop/Transfer/Input_R/Taxonomy/")
getwd()
dir()
AllSamples_classificationtaxonomy_Cenotaker2 <- as.data.frame(read_delim("CNoteTaker2_taxonomy2.txt", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE))[-1,]
AllSamples_classificationtaxonomy_Cenotaker2 <- AllSamples_classificationtaxonomy_Cenotaker2[,c(1,3)]
colnames(AllSamples_classificationtaxonomy_Cenotaker2) <- c('contig', 'Cenote_single_TaxIdpath')
rownames(AllSamples_classificationtaxonomy_Cenotaker2) <- AllSamples_classificationtaxonomy_Cenotaker2$contig
AllSamples_classificationtaxonomy_Cenotaker2$contig <- NULL
rownames(AllSamples_classificationtaxonomy_Cenotaker2) <- gsub("\\.{1}", "_", rownames(AllSamples_classificationtaxonomy_Cenotaker2))
####################################
# 6. Family-like groups
####################################
AllSamples_classification_family_groups <- as.data.frame(read_delim("family_clusters5.txt", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE))[-1,]
AllSamples_classification_family_groups <- AllSamples_classification_family_groups[,c(1,2)]
colnames(AllSamples_classification_family_groups) <- c('contig', 'Family_level')
rownames(AllSamples_classification_family_groups) <- AllSamples_classification_family_groups$contig
AllSamples_classification_family_groups$contig <- NULL
rownames(AllSamples_classification_family_groups) <- gsub("\\.{1}", "_", rownames(AllSamples_classification_family_groups))
table(AllSamples_classification_family_groups$Family_level)
####################################
# 7. Crass phage taxonomy
####################################
AllSamples_classificationtaxonomy_Crass <- as.data.frame(read_delim("AllSamples.scaffolds.crass.blastn_reliable.nocomments2.out", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE))
AllSamples_classificationtaxonomy_Crass <- AllSamples_classificationtaxonomy_Crass[,c(1,2)]
colnames(AllSamples_classificationtaxonomy_Crass) <- c('contig', 'CrAss_taxonomy')
rownames(AllSamples_classificationtaxonomy_Crass) <- AllSamples_classificationtaxonomy_Crass$contig
AllSamples_classificationtaxonomy_Crass$contig <- NULL
rownames(AllSamples_classificationtaxonomy_Crass) <- gsub("\\.{1}", "_", rownames(AllSamples_classificationtaxonomy_Crass))
####################################
# 8. Insert circularity
####################################
setwd ("/Users/daan/Desktop/Transfer/Input_R/CheckV_circularity/") 
AllSamples_classification_Circularity <- as.data.frame(read_delim("CircularContigs.txt", "\t", escape_double = FALSE, trim_ws = TRUE, col_names = FALSE))
colnames(AllSamples_classification_Circularity) <- c('contig', 'circular')
rownames(AllSamples_classification_Circularity) <- AllSamples_classification_Circularity$contig
AllSamples_classification_Circularity$contig <- NULL
rownames(AllSamples_classification_Circularity) <- gsub("\\.{1}", "_", rownames(AllSamples_classification_Circularity))
####################################
# 9. Insert Completeness
####################################
setwd ("/Users/daan/Desktop/Transfer/Input_R/CheckV_circularity/") 
getwd() 
dir()

AllSamples_completeness_Checkv<- as.data.frame(read_delim("quality_summary_R.tsv", "\t", escape_double = FALSE, trim_ws = TRUE, col_names = FALSE))
colnames(AllSamples_completeness_Checkv) <- c('contig', 'total genes', 'viral genes', 'completeness_quality', 'completeness')
rownames(AllSamples_completeness_Checkv) <- AllSamples_completeness_Checkv$contig
AllSamples_completeness_Checkv$contig <- NULL
rownames(AllSamples_completeness_Checkv) <- gsub("\\.{1}", "_", rownames(AllSamples_completeness_Checkv))

#View(AllSamples_completeness_Checkv)
####################################
# 10. Create mastertable
####################################
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

COMBO_11 <- merge(COMBO_10, AllSamples_classification_Circularity, by=0, all=T)
rownames(COMBO_11) <- COMBO_11$Row.names
COMBO_11$Row.names <- NULL
dim(COMBO_11)

COMBO_12 <- merge(COMBO_11, AllSamples_completeness_Checkv, by=0, all=T)
rownames(COMBO_12) <- COMBO_12$Row.names
COMBO_12$Row.names <- NULL
dim(COMBO_12)

COMBO_13 <- merge(COMBO_12, AllSamples_classificationtaxonomy_Crass, by=0, all=T)
rownames(COMBO_13) <- COMBO_13$Row.names
COMBO_13$Row.names <- NULL
dim(COMBO_13)

Mastertable <- merge(COMBO_13, AllSamples_classification_family_groups, by=0, all=T)
rownames(Mastertable) <- Mastertable$Row.names
Mastertable$Row.names <- NULL
dim(Mastertable)
####################################
# 10.1 Additional information to Mastertable
####################################
# 10.1.1 Insert length of contigs
####################################
Mastertable$length <- rownames(Mastertable)
Mastertable$length <- gsub("*_cov.*","",Mastertable$length)
Mastertable$length <- gsub(".*length_","", Mastertable$length)
Mastertable$length <- as.numeric(Mastertable$length)
Mastertable$length <- Mastertable$length/1000
# View(Mastertable)
####################################
# 10.1.2 Total Reads
####################################
Mastertable$Totalnumberofreads <- as.numeric(rowSums(Mastertable[,colnames(Mastertable) %in% names]))
Mastertable <- Mastertable[!(Mastertable$Totalnumberofreads == "0"),] ## remove rows with total number of zero
Mastertable$Totalnumberofreads[is.na(Mastertable$Totalnumberofreads)] <- 0
####################################
# 10.1.3 Add possible phage 
####################################
Mastertable$Family_level[is.na(Mastertable$Family_level)] <- 0
sum(Mastertable$Totalnumberofreads)
# 441.0 M PE reads total mapped reads

sum(Mastertable$Totalnumberofreads[Mastertable$Family_level==0])
# 55.8M PE reads are not in any family-level (12.7% not identifyable viruses or simply bacteria)

# Add category 'possible phage' for phage identification
Mastertable$possible_phage <- 1
Mastertable$possible_phage[!Mastertable$Family_level==0] <- "possible"
(sum(Mastertable$Totalnumberofreads[Mastertable$possible_phage=="possible"])/sum(Mastertable$Totalnumberofreads))*100
# 87.3% possible phages
####################################
# 10.1.4 Clean-up environment 
####################################
rm(COMBO_1)
rm(COMBO_2)
rm(COMBO_3)
rm(COMBO_4)
rm(COMBO_5)
rm(COMBO_6)
rm(COMBO_7)
rm(COMBO_8)
rm(COMBO_9)
rm(AllSamples_VirsorterCategories)
rm(MetaPhinder2_output_filtered)
rm(DeepVirFinder_output_filtered)
rm(AllSamples_classification_numberofgenes)
rm(AllSamples_classification_numberofpVOGs)
rm(AllSamples_classification_numberofEggNOGs)
rm(AllSamples_classification_numberofPfams)
rm(AllSamples_classification_Circularity)
rm(AllSamples_completeness_Checkv)
rm(i)
rm(j)
rm(k)
rm(list)
rm(list_txt)
rm(Mastertable_kb)
rm(columns_to_merge)
rm(vector_2)
####################################
# 11. Taxonomical annotation for homology-based tools (Blastn, Diamond, CAT)
####################################
# 11.1 List eukaryotic viruses
###################################
Eukaryotic_viruses_phylum <- c("Cressdnaviricota")
Eukaryotic_viruses_order <- c("Herpesvirales", "Ortervirales", "Zurhausenvirales", "Piccovirales", "Algavirales", "Sepolyvirales", "Chitovirales","Hepelivirales","Sobelivirales","Ourlivirales","Stellavirales","Martellivirales","Picornavirales","Tolivirales","Ghabrivirales","Amarillovirales","Wolframvirales","Patatavirales","Reovirales","Tymovirales","Halopanivirales","Belfryvirales","Rowavirales","Pimascovirales","Asfuvirales","Baphyvirales","Polivirales","Cirlivirales","Geplafuvirales","Blubervirales","Priklausovirales", "Mulpavirales")
Eukaryotic_viruses_family <- c("Pithoviridae","Tolecusatellitidae","Spiraviridae","Sarthroviridae","Permutotetraviridae","Partitiviridae","Alphasatellitidae","Anelloviridae","Baculoviridae","Hytrosaviridae","Nimaviridae","Nudiviridae","Polydnaviridae","Amalgaviridae","Birnaviridae")
####################################
# 11.2 Blastn
####################################
# 11.2.1 Add taxonomies
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
####################################
# 11.2.2 Add additional hierarchy 
####################################
Blastn_taxonomy$Blastn_viral <- 1
Blastn_taxonomy <- Blastn_taxonomy[,colnames(Blastn_taxonomy)[c(1,9,2:8)]]
####################################
# 11.2.3 Annotate eukaryotic viruses 
####################################
Blastn_taxonomy$Blastn_viral[Blastn_taxonomy$Blastn_superkingdom == "Viruses" & grepl(paste(Eukaryotic_viruses_phylum, collapse = "|"), Blastn_taxonomy$Blastn_phylum)] <- "Blastn_eukaryotic_virus"
Blastn_taxonomy$Blastn_viral[Blastn_taxonomy$Blastn_superkingdom == "Viruses" & grepl(paste(Eukaryotic_viruses_order, collapse = "|"), Blastn_taxonomy$Blastn_order)] <- "Blastn_eukaryotic_virus"
Blastn_taxonomy$Blastn_viral[Blastn_taxonomy$Blastn_superkingdom == "Viruses" & grepl(paste(Eukaryotic_viruses_family, collapse = "|"), Blastn_taxonomy$Blastn_family)] <- "Blastn_eukaryotic_virus"
####################################
# 11.2.4 Annotate prokaryotic viruses 
####################################
Blastn_taxonomy$Blastn_viral[Blastn_taxonomy$Blastn_superkingdom == "Viruses" & Blastn_taxonomy$Blastn_viral == "1"] <- "Blastn_prokaryotic_virus"
Subset_Blastn_fages <- unique(Blastn_taxonomy[Blastn_taxonomy$Blastn_viral == "Blastn_prokaryotic_virus",])
####################################
# 11.2.5 Annotate all others
####################################
Blastn_taxonomy$Blastn_viral[!Blastn_taxonomy$Blastn_superkingdom == "Viruses"] <- Blastn_taxonomy$Blastn_superkingdom[!Blastn_taxonomy$Blastn_superkingdom == "Viruses"]
####################################
# 11.2.6 Add best hit
####################################
setwd ("/Users/daan/Desktop/Transfer/Input_R/Taxonomy/")
# getwd()
# dir()
####################################
# 11.2.6.1 Add features
####################################
Blastn_best_hits <- as.data.frame(read_delim("Best_Blastn_hit_coverage_ANI.txt", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE))
colnames(Blastn_best_hits) <- c('NODES', 'Blastn_ANI', 'Blastn_coverage_absolute')
Blastn_best_hits$NODES <- gsub("\\.{1}", "_", Blastn_best_hits$NODES)
rownames(Blastn_best_hits) <- Blastn_best_hits$NODES
Blastn_best_hits$NODES <- NULL
####################################
# 11.2.6.2 Merge 
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
####################################
# 11.2.7 Clean-up environment
####################################
rm(contigs_blastn)
rm(AlSamples_classificationtaxonomy_Blastn)
rm(Subset_Blastn_fages)
rm(Blastn_best_hits)
####################################
# 11.3 Diamond
####################################
# 11.3.1 Add taxonomies
####################################
setwd("/Volumes/Daan_2/Desktop/DB/")
getwd()
dir()

contigs_diamond <- rownames(AllSamples_classificationtaxonomy_Diamond)
Diamond_taxonomy <- as.data.frame(getTaxonomy(ids = AllSamples_classificationtaxonomy_Diamond$Diamond_single_TaxIdpath, 'accessionTaxa.sql', desiredTaxa = c("superkingdom", "phylum", "class", "order", "family", "subfamily", "genus", "species")))
rownames(Diamond_taxonomy) <- contigs_diamond
colnames(Diamond_taxonomy) <- paste0('Diamond_', colnames(Diamond_taxonomy))
Diamond_taxonomy[is.na(Diamond_taxonomy)] <- "Unannotated"
####################################
# 11.3.2 Add additional hierarchy 
####################################
Diamond_taxonomy$Diamond_viral <- 1
Diamond_taxonomy <- Diamond_taxonomy[,colnames(Diamond_taxonomy)[c(1,9,2:8)]]
####################################
# 11.3.3 Annotate eukaryotic viruses 
####################################
Diamond_taxonomy$Diamond_viral[Diamond_taxonomy$Diamond_superkingdom == "Viruses" & grepl(paste(Eukaryotic_viruses_phylum, collapse = "|"), Diamond_taxonomy$Diamond_phylum)] <- "Diamond_eukaryotic_virus"
Diamond_taxonomy$Diamond_viral[Diamond_taxonomy$Diamond_superkingdom == "Viruses" & grepl(paste(Eukaryotic_viruses_order, collapse = "|"), Diamond_taxonomy$Diamond_order)] <- "Diamond_eukaryotic_virus"
Diamond_taxonomy$Diamond_viral[Diamond_taxonomy$Diamond_superkingdom == "Viruses" & grepl(paste(Eukaryotic_viruses_family, collapse = "|"), Diamond_taxonomy$Diamond_family)] <- "Diamond_eukaryotic_virus"
####################################
# 11.3.4 Annotate prokaryotic viruses 
####################################
Diamond_taxonomy$Diamond_viral[Diamond_taxonomy$Diamond_superkingdom == "Viruses" & Diamond_taxonomy$Diamond_viral == "1"] <- "Diamond_prokaryotic_virus"
Subset_Diamond_fages <- unique(Diamond_taxonomy[Diamond_taxonomy$Diamond_viral == "Diamond_prokaryotic_virus",])
####################################
# 11.3.5 Annotate all others
####################################
Diamond_taxonomy$Diamond_viral[!Diamond_taxonomy$Diamond_superkingdom == "Viruses"] <- Diamond_taxonomy$Diamond_superkingdom[!Diamond_taxonomy$Diamond_superkingdom == "Viruses"]
####################################
# 11.3.6 Add best hit
####################################
setwd ("/Users/daan/Desktop/Transfer/Input_R/Taxonomy/")
getwd()
dir()
####################################
# 11.3.6.1 Add features
####################################
Diamond_best_hits <- as.data.frame(read_delim("Best_Diamond_hit_coverage_ANI.txt", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE))
colnames(Diamond_best_hits) <- c('NODES', 'Diamond_ANI', 'Diamond_coverage_absolute')
Diamond_best_hits <- Diamond_best_hits[Diamond_best_hits$NODES %in% contigs_project,]
Diamond_best_hits$NODES <- gsub("\\.{1}", "_", Diamond_best_hits$NODES)
rownames(Diamond_best_hits) <- Diamond_best_hits$NODES
Diamond_best_hits$NODES <- NULL
####################################
# 11.3.6.2 Merge
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
####################################
# 11.3.7 Clean-up environment
####################################
rm(contigs_diamond)
rm(Subset_Diamond_fages)
rm(Diamond_best_hits)
####################################
# 11.4 CAT
####################################
# 11.4.1 Add taxonomies
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
# 11.4.2 Add additional hierarchy 
####################################
CAT_taxonomy$CAT_viral <- 1
CAT_taxonomy <- CAT_taxonomy[,colnames(CAT_taxonomy)[c(1,9,2:8)]]
####################################
# 11.4.3 Annotate eukaryotic viruses 
####################################
CAT_taxonomy$CAT_viral[CAT_taxonomy$CAT_superkingdom == "Viruses" & grepl(paste(Eukaryotic_viruses_phylum, collapse = "|"), CAT_taxonomy$CAT_phylum)] <- "CAT_eukaryotic_virus"
CAT_taxonomy$CAT_viral[CAT_taxonomy$CAT_superkingdom == "Viruses" & grepl(paste(Eukaryotic_viruses_order, collapse = "|"), CAT_taxonomy$CAT_order)] <- "CAT_eukaryotic_virus"
CAT_taxonomy$CAT_viral[CAT_taxonomy$CAT_superkingdom == "Viruses" & grepl(paste(Eukaryotic_viruses_family, collapse = "|"), CAT_taxonomy$CAT_family)] <- "CAT_eukaryotic_virus"
####################################
# 11.4.4 Annotate prokaryotic viruses 
####################################
CAT_taxonomy$CAT_viral[CAT_taxonomy$CAT_superkingdom == "Viruses" & CAT_taxonomy$CAT_viral == "1"] <- "CAT_prokaryotic_virus"
Subset_CAT_fages <- unique(CAT_taxonomy[CAT_taxonomy$CAT_viral == "CAT_prokaryotic_virus",])
####################################
# 11.4.5 Annotate all others
####################################
CAT_taxonomy$CAT_viral[!CAT_taxonomy$CAT_superkingdom == "Viruses"] <- CAT_taxonomy$CAT_superkingdom[!CAT_taxonomy$CAT_superkingdom == "Viruses"]
####################################
# 11.4.6 Add ORFs
####################################
setwd ("/Users/daan/Desktop/Transfer/Input_R/Taxonomy/")
# getwd()
# dir()
####################################
# 11.4.6.1 Add features
####################################
# IN CODE, cut -f1,4 and not -f1,5
CAT_number_of_ORFs_for_annotation_LCA <- as.data.frame(read_delim("CAT_alignment_ORFs.txt", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE))
colnames(CAT_number_of_ORFs_for_annotation_LCA) <- c('NODES', 'CAT_ORFs')
CAT_number_of_ORFs_for_annotation_LCA <- CAT_number_of_ORFs_for_annotation_LCA[CAT_number_of_ORFs_for_annotation_LCA$NODES %in% contigs_project,]
CAT_number_of_ORFs_for_annotation_LCA$NODES <- gsub("\\.{1}", "_", CAT_number_of_ORFs_for_annotation_LCA$NODES)
rownames(CAT_number_of_ORFs_for_annotation_LCA) <- CAT_number_of_ORFs_for_annotation_LCA$NODES
CAT_number_of_ORFs_for_annotation_LCA$NODES <- NULL
####################################
# 11.4.6.2 Merge
####################################
CAT_taxonomy <- transform(merge(CAT_taxonomy, CAT_number_of_ORFs_for_annotation_LCA, by = 0, all = F), row.names=Row.names, Row.names = NULL)
CAT_taxonomy$CAT_ORFs[is.na(CAT_taxonomy$CAT_ORFs)] <- 0

## length
CAT_taxonomy$length <- rownames(CAT_taxonomy)
CAT_taxonomy$length <- gsub("*_cov.*","",CAT_taxonomy$length)
CAT_taxonomy$length <- gsub(".*length_","", CAT_taxonomy$length)
CAT_taxonomy$length <- as.numeric(CAT_taxonomy$length)

## calculate %
CAT_taxonomy$CAT_ORFs_per_length <- (CAT_taxonomy$CAT_ORFs/(CAT_taxonomy$length/1000))
CAT_taxonomy$length <- NULL
####################################
# 11.4.7 Clean-up environment
####################################
rm(contigs_CAT)
rm(Subset_CAT_fages)
rm(CAT_number_of_ORFs_for_annotation_LCA)
####################################
# 1.5 CNoteTaker2 functional annotation
####################################
####################################
# 11.5.1 Add taxonomies
####################################
setwd("/Volumes/Daan_2/Desktop/DB/")
getwd()
dir()

contigs_cenote <- rownames(AllSamples_classificationtaxonomy_Cenotaker2)
cenote_taxonomy <- as.data.frame(getTaxonomy(ids = AllSamples_classificationtaxonomy_Cenotaker2$Cenote_single_TaxIdpath, 'accessionTaxa.sql', desiredTaxa = c("superkingdom", "phylum", "class", "order", "family")))
rownames(cenote_taxonomy) <- contigs_cenote
colnames(cenote_taxonomy) <- paste0('Cenote_', colnames(cenote_taxonomy))
cenote_taxonomy <- na.omit(cenote_taxonomy)                          
cenote_taxonomy[is.na(cenote_taxonomy)] <- "Unannotated"
#View(cenote_taxonomy)
####################################
# 11.5.2 Add additional hierarchy 
####################################
cenote_taxonomy$Cenote_viral <- 1
cenote_taxonomy <- cenote_taxonomy[,colnames(cenote_taxonomy)[c(1,6,2:5)]]
####################################
# 11.5.3 Annotate eukaryotic viruses 
####################################
cenote_taxonomy$Cenote_viral[cenote_taxonomy$Cenote_superkingdom == "Viruses" & grepl(paste(Eukaryotic_viruses_phylum, collapse = "|"), cenote_taxonomy$Cenote_phylum)] <- "Cenote_eukaryotic_virus"
cenote_taxonomy$Cenote_viral[cenote_taxonomy$Cenote_superkingdom == "Viruses" & grepl(paste(Eukaryotic_viruses_order, collapse = "|"), cenote_taxonomy$Cenote_order)] <- "Cenote_eukaryotic_virus"
cenote_taxonomy$Cenote_viral[cenote_taxonomy$Cenote_superkingdom == "Viruses" & grepl(paste(Eukaryotic_viruses_family, collapse = "|"), cenote_taxonomy$Cenote_family)] <- "Cenote_eukaryotic_virus"
####################################
# 11.5.4 Annotate prokaryotic viruses 
####################################
cenote_taxonomy$Cenote_viral[cenote_taxonomy$Cenote_superkingdom == "Viruses" & cenote_taxonomy$Cenote_viral == "1"] <- "cenote_prokaryotic_virus"
Subset_cenote_fages <- unique(cenote_taxonomy[cenote_taxonomy$Cenote_viral == "cenote_prokaryotic_virus",])
#View(Subset_cenote_fages)
####################################
# 11.5.5 Annotate all others
####################################
cenote_taxonomy$Cenote_viral[!cenote_taxonomy$Cenote_superkingdom == "Viruses"] <- cenote_taxonomy$Cenote_superkingdom[!cenote_taxonomy$Cenote_superkingdom == "Viruses"]
####################################
# 11.6 Merge taxonomies with mastertable
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
Mastertable <- Mastertable[Mastertable$Totalnumberofreads > 0,]
# View(Mastertable)

# Clean environment
rm(Mastertable_corrected)
rm(Blastn_taxonomy)
rm(Diamond_taxonomy)
rm(CAT_taxonomy)
####################################
# 12. Phage-identification 
####################################
# 12.1 High-quality phage identification (CheckV + Vs2)
####################################
## The Phage identification is based on a combination of two tools: (i) Virsorter2 and (ii) CheckV:

## We identify phages as:
## virsorter2: dsDNAphages; CheckV > 50% complete
## virsorter2: RNA phages; CheckV > 50% complete
## virsorter2: ssDNA phages; length >= 3kb; at leas 1 viral gene predicted by CheckV
Mastertable$completeness[is.na(Mastertable$completeness)] <- 0

Mastertable$Phage_score1 <- 1
Mastertable$Phage_score1[Mastertable$Virsorter == "dsDNAphage" & Mastertable$completeness >= 50 & Mastertable$length >= 3] <- "phage"
Mastertable$Phage_score1[Mastertable$Virsorter == "RNA" & Mastertable$completeness >= 50 & Mastertable$length >= 3] <- "phage"
Mastertable$Phage_score1[Mastertable$Virsorter == "ssDNA" & Mastertable$length >= 3] <- "phage"
Mastertable$Phage_score1[!Mastertable$Diamond_viral=="Diamond_eukaryotic_virus" & !Mastertable$Blastn_viral=="Blastn_eukaryotic_virus" & !Mastertable$CAT_viral=="CAT_eukaryotic_virus" & !Mastertable$Cenote_viral=="Cenote_eukaryotic_virus"] <- 1
Mastertable$Totalnumberofreads[is.na(Mastertable$Totalnumberofreads)] <- 0
Mastertable_1 <- Mastertable[!Mastertable$Totalnumberofreads==0,]
table(Mastertable_1$Phage_score1)

Mastertable$Phage_score1[Mastertable$Virsorter == "dsDNAphage" & Mastertable$completeness >= 50 & Mastertable$length >= 3 & !Mastertable$Diamond_viral=="Diamond_eukaryotic_virus" & !Mastertable$Blastn_viral=="Blastn_eukaryotic_virus" & !Mastertable$CAT_viral=="CAT_eukaryotic_virus" & !Mastertable$Cenote_viral=="Cenote_eukaryotic_virus"] <- "phage"
Mastertable$Phage_score1[Mastertable$Virsorter == "RNA" & Mastertable$completeness >= 50 & Mastertable$length >= 3 & !Mastertable$Diamond_viral=="Diamond_eukaryotic_virus" & !Mastertable$Blastn_viral=="Blastn_eukaryotic_virus" & !Mastertable$CAT_viral=="CAT_eukaryotic_virus" & !Mastertable$Cenote_viral=="Cenote_eukaryotic_virus"] <- "phage"
Mastertable$Phage_score1[Mastertable$Virsorter == "ssDNA" & Mastertable$length >= 3 & !Mastertable$Diamond_viral=="Diamond_eukaryotic_virus" & !Mastertable$Blastn_viral=="Blastn_eukaryotic_virus" & !Mastertable$CAT_viral=="CAT_eukaryotic_virus" & !Mastertable$Cenote_viral=="Cenote_eukaryotic_virus"] <- "phage"
####################################
# 12.2 Include phages sharing same family
####################################
Mastertable$Phage_score2 <- 1
family_phages <- unique(Mastertable$Family_level[Mastertable$Phage_score1=="phage"])
Mastertable$Phage_score2[Mastertable$Family_level %in% family_phages] <- "phage_family"
Mastertable$Phage_score2[Mastertable$Blastn_viral=="Blastn_eukaryotic_virus" | Mastertable$Diamond_viral=="Diamond_eukaryotic_virus" | Mastertable$CAT_viral=="CAT_eukaryotic_virus" | Mastertable$Cenote_viral=="CAT_eukaryotic_virus"] <- 1

# Do not include fragments for phages
# dsDNA < 10kb 
Mastertable$Phage_score2[Mastertable$Virsorter == "dsDNAphage" & Mastertable$length < 3] <- 1

# ssDNA < 3kb 
Mastertable$Phage_score2[Mastertable$Virsorter == "RNA" & Mastertable$length < 3] <- 1

# RNA < 3kb 
Mastertable$Phage_score2[Mastertable$Virsorter == "ssDNA" & Mastertable$length < 3] <- 1
 
table(Mastertable$Phage_score2)
family_phages <- unique(Mastertable$Family_level[Mastertable$Phage_score2=="phage_family"])
family_phages 

# In total 122 phage-like groups
(sum(Mastertable$Totalnumberofreads[Mastertable$Phage_score2=="phage_family"])/sum(Mastertable$Totalnumberofreads))*100 

# Encompassing 58.2% of the reads
Mastertable$names <- rownames(Mastertable)
sort(colSums(Mastertable[!Mastertable$Phage_score2==1,1:127]))
Mastertable$Phages <- Mastertable$Phage_score2

Mastertable$Phages[Mastertable$Phages == "phage_family"] <- "phage"
table(Mastertable$Phage_score1)
table(Mastertable$Phages)
####################################
# 13. Final classification
####################################
# 13.1 Final classification of eukaryotic viruses
####################################
# 13.1.1 Viral
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
table(Mastertable$Final_viral) # 150 eukaryotic viruses
table(Mastertable$Phage_score2) 
####################################
# 13.1.2 Superkingdom
####################################
Mastertable$Final_superkingdom <- 1
Mastertable$Final_superkingdom[Mastertable$Diamond_viral == "Diamond_eukaryotic_virus"] <- Mastertable$Diamond_superkingdom[Mastertable$Diamond_viral == "Diamond_eukaryotic_virus"]
Mastertable$Final_superkingdom[!Mastertable$Diamond_viral == "Diamond_eukaryotic_virus" & Mastertable$Blastn_viral == "Blastn_eukaryotic_virus"] <- Mastertable$Blastn_superkingdom[!Mastertable$Diamond_viral == "Diamond_eukaryotic_virus" & Mastertable$Blastn_viral == "Blastn_eukaryotic_virus"]
table(Mastertable$Final_superkingdom) 
####################################
# 13.1.3 Phylum
####################################
Mastertable$Final_phylum <- 1
Mastertable$Final_phylum[Mastertable$Diamond_viral == "Diamond_eukaryotic_virus"] <- Mastertable$Diamond_phylum[Mastertable$Diamond_viral == "Diamond_eukaryotic_virus"]
Mastertable$Final_phylum[!Mastertable$Diamond_viral == "Diamond_eukaryotic_virus" & Mastertable$Blastn_viral == "Blastn_eukaryotic_virus"] <- Mastertable$Blastn_phylum[!Mastertable$Diamond_viral == "Diamond_eukaryotic_virus" & Mastertable$Blastn_viral == "Blastn_eukaryotic_virus"]
table(Mastertable$Final_phylum) 
# How many different eukaryotic viral phyla are found?
####################################
# 13.1.4 Class
####################################
Mastertable$Final_class <- 1
Mastertable$Final_class[Mastertable$Diamond_viral == "Diamond_eukaryotic_virus"] <- Mastertable$Diamond_class[Mastertable$Diamond_viral == "Diamond_eukaryotic_virus"]
Mastertable$Final_class[!Mastertable$Diamond_viral == "Diamond_eukaryotic_virus" & Mastertable$Blastn_viral == "Blastn_eukaryotic_virus"] <- Mastertable$Blastn_class[!Mastertable$Diamond_viral == "Diamond_eukaryotic_virus" & Mastertable$Blastn_viral == "Blastn_eukaryotic_virus"]
table(Mastertable$Final_class) 
# How many different eukaryotic viral classes are found?
####################################
# 13.1.5 Order
####################################
Mastertable$Final_order <- 1
Mastertable$Final_order[Mastertable$Diamond_viral == "Diamond_eukaryotic_virus"] <- Mastertable$Diamond_order[Mastertable$Diamond_viral == "Diamond_eukaryotic_virus"]
Mastertable$Final_order[!Mastertable$Diamond_viral == "Diamond_eukaryotic_virus" & Mastertable$Blastn_viral == "Blastn_eukaryotic_virus"] <- Mastertable$Blastn_order[!Mastertable$Diamond_viral == "Diamond_eukaryotic_virus" & Mastertable$Blastn_viral == "Blastn_eukaryotic_virus"]
table(Mastertable$Final_order) 
# How many different eukaryotic viral orders are found?
####################################
# 13.1.6 Family
####################################
Mastertable$Final_family <- 1
Mastertable$Final_family[Mastertable$Diamond_viral == "Diamond_eukaryotic_virus"] <- Mastertable$Diamond_family[Mastertable$Diamond_viral == "Diamond_eukaryotic_virus"]
Mastertable$Final_family[!Mastertable$Diamond_viral == "Diamond_eukaryotic_virus" & Mastertable$Blastn_viral == "Blastn_eukaryotic_virus"] <- Mastertable$Blastn_family[!Mastertable$Diamond_viral == "Diamond_eukaryotic_virus" & Mastertable$Blastn_viral == "Blastn_eukaryotic_virus"]
table(Mastertable$Final_family) 
# How many different eukaryotic viral families are found?
####################################
# 13.1.7 Subfamily
####################################
Mastertable$Final_subfamily <- 1
Mastertable$Final_subfamily[Mastertable$Diamond_viral == "Diamond_eukaryotic_virus"] <- Mastertable$Diamond_subfamily[Mastertable$Diamond_viral == "Diamond_eukaryotic_virus"]
Mastertable$Final_subfamily[!Mastertable$Diamond_viral == "Diamond_eukaryotic_virus" & Mastertable$Blastn_viral == "Blastn_eukaryotic_virus"] <- Mastertable$Blastn_subfamily[!Mastertable$Diamond_viral == "Diamond_eukaryotic_virus" & Mastertable$Blastn_viral == "Blastn_eukaryotic_virus"]
table(Mastertable$Final_subfamily) 
# How many different eukaryotic subfamilies orders are found?
####################################
# 13.1.8 Genus
####################################
Mastertable$Final_genus <- 1
Mastertable$Final_genus[Mastertable$Diamond_viral == "Diamond_eukaryotic_virus"] <- Mastertable$Diamond_genus[Mastertable$Diamond_viral == "Diamond_eukaryotic_virus"]
Mastertable$Final_genus[!Mastertable$Diamond_viral == "Diamond_eukaryotic_virus" & Mastertable$Blastn_viral == "Blastn_eukaryotic_virus"] <- Mastertable$Blastn_genus[!Mastertable$Diamond_viral == "Diamond_eukaryotic_virus" & Mastertable$Blastn_viral == "Blastn_eukaryotic_virus"]
table(Mastertable$Final_genus) 
# How many different eukaryotic viral genera are found?
####################################
# 13.1.9 Species
####################################
Mastertable$Final_species <- 1
Mastertable$Final_species[Mastertable$Diamond_viral == "Diamond_eukaryotic_virus"] <- Mastertable$Diamond_species[Mastertable$Diamond_viral == "Diamond_eukaryotic_virus"]
Mastertable$Final_species[!Mastertable$Diamond_viral == "Diamond_eukaryotic_virus" & Mastertable$Blastn_viral == "Blastn_eukaryotic_virus"] <- Mastertable$Blastn_species[!Mastertable$Diamond_viral == "Diamond_eukaryotic_virus" & Mastertable$Blastn_viral == "Blastn_eukaryotic_virus"]
table(Mastertable$Final_species) 
# How many different eukaryotic viral species are found?
####################################
# 13.1.10 alignment quality (ANI, coverage)
####################################
Mastertable$Final_ANI <- 1
Mastertable$Final_ANI[Mastertable$Diamond_viral == "Diamond_eukaryotic_virus"] <- Mastertable$Diamond_ANI[Mastertable$Diamond_viral == "Diamond_eukaryotic_virus"]
Mastertable$Final_ANI[!Mastertable$Diamond_viral == "Diamond_eukaryotic_virus" & Mastertable$Blastn_viral == "Blastn_eukaryotic_virus"] <- Mastertable$Blastn_ANI[!Mastertable$Diamond_viral == "Diamond_eukaryotic_virus" & Mastertable$Blastn_viral == "Blastn_eukaryotic_virus"]

Mastertable$Final_coverage <- 1
Mastertable$Final_coverage[Mastertable$Diamond_viral == "Diamond_eukaryotic_virus"] <- Mastertable$Diamond_coverage_relative[Mastertable$Diamond_viral == "Diamond_eukaryotic_virus"]
Mastertable$Final_coverage[!Mastertable$Diamond_viral == "Diamond_eukaryotic_virus" & Mastertable$Blastn_viral == "Blastn_eukaryotic_virus"] <- Mastertable$Blastn_coverage_relative[!Mastertable$Diamond_viral == "Diamond_eukaryotic_virus" & Mastertable$Blastn_viral == "Blastn_eukaryotic_virus"]

#View(Mastertable)
####################################
# 13.2 Final classification of prokaryotic viruses
####################################
# 13.2.1 Temp. Class annotation phages - Virak
####################################
# Add temp. phage class classification (vs2+comple+3kbssDNA)
Mastertable$Final_viral[Mastertable$Phages == "phage" & !Mastertable$Diamond_viral == "Diamond_eukaryotic_virus" & !Mastertable$Blastn_viral == "Blastn_eukaryotic_virus"] <- Mastertable$Phages[Mastertable$Phages == "phage" & !Mastertable$Diamond_viral == "Diamond_eukaryotic_virus" & !Mastertable$Blastn_viral == "Blastn_eukaryotic_virus"]

Mastertable$Final_class_temp <- 1
Mastertable$Final_class_temp[!Mastertable$Final_viral == "eukaryotic virus" & Mastertable$Phage_score1 == "phage" & Mastertable$Diamond_viral == "Diamond_prokaryotic_virus"] <- Mastertable$Diamond_class[!Mastertable$Final_viral == "eukaryotic virus" & Mastertable$Phage_score1 == "phage" & Mastertable$Diamond_viral == "Diamond_prokaryotic_virus"]
Mastertable$Final_class_temp[!Mastertable$Final_viral == "eukaryotic virus" & Mastertable$Phage_score1 == "phage" & !Mastertable$Diamond_viral == "Diamond_prokaryotic_virus" & Mastertable$Blastn_viral == "Blastn_prokaryotic_virus"] <- Mastertable$Blastn_class[!Mastertable$Final_viral == "eukaryotic virus" & Mastertable$Phage_score1 == "phage" & !Mastertable$Diamond_viral == "Diamond_prokaryotic_virus" & Mastertable$Blastn_viral == "Blastn_prokaryotic_virus"]
Mastertable$Final_class_temp[!Mastertable$Final_viral == "eukaryotic virus" & Mastertable$Phage_score1 == "phage" & !Mastertable$Diamond_viral == "Diamond_prokaryotic_virus" & !Mastertable$Blastn_viral == "Blastn_prokaryotic_virus" & Mastertable$CAT_viral=="CAT_prokaryotic_virus"] <- Mastertable$CAT_class[!Mastertable$Final_viral == "eukaryotic virus" & Mastertable$Phage_score1 == "phage" & !Mastertable$Diamond_viral == "Diamond_prokaryotic_virus" & !Mastertable$Blastn_viral == "Blastn_prokaryotic_virus" & Mastertable$CAT_viral=="CAT_prokaryotic_virus"]
Mastertable$Final_class_temp[!Mastertable$Final_viral == "eukaryotic virus" & Mastertable$Phage_score1 == "phage" & !Mastertable$Diamond_viral == "Diamond_prokaryotic_virus" & !Mastertable$Blastn_viral == "Blastn_prokaryotic_virus" & !Mastertable$CAT_viral=="CAT_prokaryotic_virus" & Mastertable$Cenote_viral=="cenote_prokaryotic_virus"] <- Mastertable$Cenote_class[!Mastertable$Final_viral == "eukaryotic virus" & Mastertable$Phage_score1 == "phage" & !Mastertable$Diamond_viral == "Diamond_prokaryotic_virus" & !Mastertable$Blastn_viral == "Blastn_prokaryotic_virus" & !Mastertable$CAT_viral=="CAT_prokaryotic_virus" & Mastertable$Cenote_viral=="cenote_prokaryotic_virus"]

# Adjust Classs based on length cutoff
Mastertable$Final_class_temp[Mastertable$Final_class_temp=="Caudoviricetes" & Mastertable$length<10] <- "Unannotated"
Mastertable$Final_class_temp[Mastertable$Final_class_temp=="Malgrandaviricetes" & Mastertable$length>10] <- "Unannotated"

# Convert unannotated to '1'
Mastertable$Final_class_temp[Mastertable$Final_class_temp=="Unannotated"] <- 1

# Rename unclassified ones '1' as in line with other unreliable annotations
Mastertable$Final_class_temp[Mastertable$Final_class_temp=="Unannotated"] <- 1

# Create unique vector of temp. classification for (I) Caudoviricetes and (II) Malgrandaviricetes
Phage_class_caudo <- unique(Mastertable[!Mastertable$Final_class_temp==1 & Mastertable$Final_class_temp=="Caudoviricetes",c(149,209)])
vector_family_caudo <- Phage_class_caudo$Family_level
Phage_class_malgr <- unique(Mastertable[!Mastertable$Final_class_temp==1 & Mastertable$Final_class_temp=="Malgrandaviricetes",c(149,209)])
vector_family_malg <- Phage_class_malgr$Family_level

# Classify other phages belong to same family-like group (RELIABLE PHAGE CLASS CLASSIFICATIONS)
Mastertable$Final_class_temp[Mastertable$Family_level %in% vector_family_caudo] <- "Caudoviricetes"
Mastertable$Final_class_temp[Mastertable$Family_level %in% vector_family_malg] <- "Malgrandaviricetes"
####################################
# 13.2.2 Superkingdom
####################################
Mastertable$Final_superkingdom[Mastertable$Phages=="phage" & Mastertable$Final_class_temp=="Caudoviricetes"] <- "Viruses"
Mastertable$Final_superkingdom[Mastertable$Phages=="phage" & Mastertable$Final_class_temp=="Malgrandaviricetes"] <- "Viruses"
Mastertable$Final_superkingdom[Mastertable$Phages=="phage" & Mastertable$Final_superkingdom==1] <- "Viruses"
table(Mastertable$Final_superkingdom)
####################################
# 13.2.3 Phylum
####################################
Mastertable$Final_phylum[Mastertable$Phages=="phage" & Mastertable$Final_class_temp=="Caudoviricetes"] <- "Uroviricota"
Mastertable$Final_phylum[Mastertable$Phages=="phage" & Mastertable$Final_class_temp=="Malgrandaviricetes"] <- "Phixviricota"
Mastertable$Final_phylum[Mastertable$Phages=="phage" & Mastertable$Final_phylum==1] <- "Unannotated"
table(Mastertable$Final_phylum)
####################################
# 13.2.4 Class
####################################
Mastertable$Final_class[Mastertable$Phages=="phage" & Mastertable$Final_class_temp=="Caudoviricetes"] <- "Caudoviricetes"
Mastertable$Final_class[Mastertable$Phages=="phage" & Mastertable$Final_class_temp=="Malgrandaviricetes"] <- "Malgrandaviricetes"
Mastertable$Final_class[Mastertable$Phages=="phage" & Mastertable$Final_class==1] <- "Unannotated"
table(Mastertable$Final_class)
####################################
# 13.2.5 Order
####################################
Mastertable$Final_order[Mastertable$Phages=="phage"] <- "Unannotated"
table(Mastertable$Final_order)
####################################
# 13.2.6 Family
####################################
Mastertable$Final_family[Mastertable$Phages=="phage"] <- "Unannotated"
table(Mastertable$Final_family)
####################################
# 13.2.7 Subfamily
####################################
Mastertable$Final_subfamily[Mastertable$Phages=="phage"] <- "Unannotated"
table(Mastertable$Final_subfamily)
####################################
# 13.2.8 Genus
####################################
Mastertable$Final_genus[Mastertable$Phages=="phage"] <- "Unannotated"
table(Mastertable$Final_genus)
####################################
# 13.2.9 Species
####################################
Mastertable$Final_species[Mastertable$Phages=="phage"] <- "Unannotated"
table(Mastertable$Final_species)
####################################
# 13.3 Final classification of eukaryotes
####################################
Mastertable$Diamond_superkingdom[is.na(Mastertable$Diamond_superkingdom)] <- 0
Mastertable$Blastn_superkingdom[is.na(Mastertable$Blastn_superkingdom)] <- 0
Mastertable$CAT_superkingdom[is.na(Mastertable$CAT_superkingdom)] <- 0

Mastertable$Final_superkingdom[!Mastertable$Final_superkingdom == "Viruses" & Mastertable$Diamond_superkingdom == "Eukaryota"] <- Mastertable$Diamond_superkingdom[!Mastertable$Final_superkingdom == "Viruses" & Mastertable$Diamond_superkingdom == "Eukaryota"]
Mastertable$Final_superkingdom[!Mastertable$Final_superkingdom == "Viruses" & !Mastertable$Diamond_superkingdom == "Eukaryota" & Mastertable$Blastn_superkingdom == "Eukaryota"] <- Mastertable$Blastn_superkingdom[!Mastertable$Final_superkingdom == "Viruses" & !Mastertable$Diamond_superkingdom == "Eukaryota" & Mastertable$Blastn_superkingdom == "Eukaryota"]
Mastertable$Final_superkingdom[!Mastertable$Final_superkingdom == "Viruses" & !Mastertable$Diamond_superkingdom == "Eukaryota" & !Mastertable$Blastn_superkingdom == "Eukaryota" & Mastertable$CAT_superkingdom == "Eukaryota"] <- Mastertable$CAT_superkingdom[!Mastertable$Final_superkingdom == "Viruses" & !Mastertable$Diamond_superkingdom == "Eukaryota" & !Mastertable$Blastn_superkingdom == "Eukaryota" & Mastertable$CAT_superkingdom == "Eukaryota"]
table(Mastertable$Final_superkingdom)
####################################
# 13.4 Final classification of archaea
####################################
Mastertable$Final_superkingdom[!Mastertable$Final_superkingdom == "Viruses" & !Mastertable$Final_superkingdom == "Eukaryota" & Mastertable$Diamond_superkingdom == "Archaea"] <- Mastertable$Diamond_superkingdom[!Mastertable$Final_superkingdom == "Viruses" & !Mastertable$Final_superkingdom == "Eukaryota" & Mastertable$Diamond_superkingdom == "Archaea"]
Mastertable$Final_superkingdom[!Mastertable$Final_superkingdom == "Viruses" & !Mastertable$Final_superkingdom == "Eukaryota" & !Mastertable$Diamond_superkingdom == "Archaea" & Mastertable$Blastn_superkingdom == "Archaea"] <- Mastertable$Blastn_superkingdom[!Mastertable$Final_superkingdom == "Viruses" & !Mastertable$Final_superkingdom == "Eukaryota" & !Mastertable$Diamond_superkingdom == "Archaea" & Mastertable$Blastn_superkingdom == "Archaea"]
Mastertable$Final_superkingdom[!Mastertable$Final_superkingdom == "Viruses" & !Mastertable$Final_superkingdom == "Eukaryota" & !Mastertable$Diamond_superkingdom == "Archaea" & !Mastertable$Blastn_superkingdom == "Archaea" & Mastertable$CAT_superkingdom == "Archaea"] <- Mastertable$CAT_superkingdom[!Mastertable$Final_superkingdom == "Viruses" & !Mastertable$Final_superkingdom == "Eukaryota" & !Mastertable$Diamond_superkingdom == "Archaea" & !Mastertable$Blastn_superkingdom == "Archaea" & Mastertable$CAT_superkingdom == "Archaea"]
####################################
# 13.5 Final classification of bacteria
####################################
Mastertable$Final_superkingdom[!Mastertable$Final_superkingdom == "Viruses" & !Mastertable$Final_superkingdom == "Eukaryota" & !Mastertable$Final_superkingdom == "Archaea" & Mastertable$Diamond_superkingdom == "Bacteria"] <- Mastertable$Diamond_superkingdom[!Mastertable$Final_superkingdom == "Viruses" & !Mastertable$Final_superkingdom == "Eukaryota" & !Mastertable$Final_superkingdom == "Archaea" & Mastertable$Diamond_superkingdom == "Bacteria"]
Mastertable$Final_superkingdom[!Mastertable$Final_superkingdom == "Viruses" & !Mastertable$Final_superkingdom == "Eukaryota" & !Mastertable$Final_superkingdom == "Archaea" & !Mastertable$Diamond_superkingdom == "Bacteria" & Mastertable$Blastn_superkingdom == "Bacteria"] <- Mastertable$Blastn_superkingdom[!Mastertable$Final_superkingdom == "Viruses" & !Mastertable$Final_superkingdom == "Eukaryota" & !Mastertable$Final_superkingdom == "Archaea" & !Mastertable$Diamond_superkingdom == "Bacteria" & Mastertable$Blastn_superkingdom == "Bacteria"]
Mastertable$Final_superkingdom[!Mastertable$Final_superkingdom == "Viruses" & !Mastertable$Final_superkingdom == "Eukaryota" & !Mastertable$Final_superkingdom == "Archaea" & !Mastertable$Diamond_superkingdom == "Bacteria" & !Mastertable$Blastn_superkingdom == "Bacteria" & Mastertable$CAT_superkingdom == "Bacteria"] <- Mastertable$CAT_superkingdom[!Mastertable$Final_superkingdom == "Viruses" & !Mastertable$Final_superkingdom == "Eukaryota" & !Mastertable$Final_superkingdom == "Archaea" & !Mastertable$Diamond_superkingdom == "Bacteria" & !Mastertable$Blastn_superkingdom == "Bacteria" & Mastertable$CAT_superkingdom == "Bacteria"]
####################################
# 13.6 Final classification of Other
####################################
Mastertable$Final_superkingdom[!Mastertable$Final_superkingdom == "Viruses" & !Mastertable$Final_superkingdom == "Eukaryota" & !Mastertable$Final_superkingdom == "Archaea" & !Mastertable$Final_superkingdom == "Bacteria" & Mastertable$Diamond_superkingdom == "Other"] <- Mastertable$Diamond_superkingdom[!Mastertable$Final_superkingdom == "Viruses" & !Mastertable$Final_superkingdom == "Eukaryota" & !Mastertable$Final_superkingdom == "Archaea" & !Mastertable$Final_superkingdom == "Bacteria" & Mastertable$Diamond_superkingdom == "Other"]
Mastertable$Final_superkingdom[!Mastertable$Final_superkingdom == "Viruses" & !Mastertable$Final_superkingdom == "Eukaryota" & !Mastertable$Final_superkingdom == "Archaea" & !Mastertable$Final_superkingdom == "Bacteria" & !Mastertable$Diamond_superkingdom == "Other" & Mastertable$Blastn_superkingdom == "Other"] <- Mastertable$Blastn_superkingdom[!Mastertable$Final_superkingdom == "Viruses" & !Mastertable$Final_superkingdom == "Eukaryota" & !Mastertable$Final_superkingdom == "Archaea" & !Mastertable$Final_superkingdom == "Bacteria" & !Mastertable$Diamond_superkingdom == "Other" & Mastertable$Blastn_superkingdom == "Other"]
Mastertable$Final_superkingdom[!Mastertable$Final_superkingdom == "Viruses" & !Mastertable$Final_superkingdom == "Eukaryota" & !Mastertable$Final_superkingdom == "Archaea" & !Mastertable$Final_superkingdom == "Bacteria" & !Mastertable$Diamond_superkingdom == "Other" & !Mastertable$Blastn_superkingdom == "Other" & Mastertable$CAT_superkingdom == "Other"] <- Mastertable$CAT_superkingdom[!Mastertable$Final_superkingdom == "Viruses" & !Mastertable$Final_superkingdom == "Eukaryota" & !Mastertable$Final_superkingdom == "Archaea" & !Mastertable$Final_superkingdom == "Bacteria" & !Mastertable$Diamond_superkingdom == "Other" & !Mastertable$Blastn_superkingdom == "Other" & Mastertable$CAT_superkingdom == "Other"]
####################################
# 13.7 Final classification of Unclassified
####################################
# You could, if you like, scan the "unclassified" sequences manully that have an annotation to check if you have some unclassified viral sequences. Since this we take a lot of time and makes it impossible to automize this pipeline & it doesn't add a lot I don't do this.
Mastertable$Final_superkingdom[!Mastertable$Final_superkingdom == "Viruses" & !Mastertable$Final_superkingdom == "Eukaryota" & !Mastertable$Final_superkingdom == "Archaea" & !Mastertable$Final_superkingdom == "Bacteria" & !Mastertable$Diamond_superkingdom == "Other" & Mastertable$Diamond_superkingdom == "Unannotated"] <- Mastertable$Diamond_superkingdom[!Mastertable$Final_superkingdom == "Viruses" & !Mastertable$Final_superkingdom == "Eukaryota" & !Mastertable$Final_superkingdom == "Archaea" & !Mastertable$Final_superkingdom == "Bacteria" & Mastertable$Diamond_superkingdom == "Unannotated"]
Mastertable$Final_superkingdom[!Mastertable$Final_superkingdom == "Viruses" & !Mastertable$Final_superkingdom == "Eukaryota" & !Mastertable$Final_superkingdom == "Archaea" & !Mastertable$Final_superkingdom == "Bacteria" & !Mastertable$Diamond_superkingdom == "Other" & !Mastertable$Diamond_superkingdom == "Unannotated" & Mastertable$Blastn_superkingdom == "Unannotated"] <- Mastertable$Blastn_superkingdom[!Mastertable$Final_superkingdom == "Viruses" & !Mastertable$Final_superkingdom == "Eukaryota" & !Mastertable$Final_superkingdom == "Archaea" & !Mastertable$Final_superkingdom == "Bacteria" & !Mastertable$Diamond_superkingdom == "Other" & !Mastertable$Diamond_superkingdom == "Unannotated" & Mastertable$Blastn_superkingdom == "Unannotated"]
Mastertable$Final_superkingdom[!Mastertable$Final_superkingdom == "Viruses" & !Mastertable$Final_superkingdom == "Eukaryota" & !Mastertable$Final_superkingdom == "Archaea" & !Mastertable$Final_superkingdom == "Bacteria" & !Mastertable$Diamond_superkingdom == "Other" & !Mastertable$Diamond_superkingdom == "Unannotated" & !Mastertable$Blastn_superkingdom == "Unannotated" & Mastertable$CAT_superkingdom == "Unannotated"] <- Mastertable$CAT_superkingdom[!Mastertable$Final_superkingdom == "Viruses" & !Mastertable$Final_superkingdom == "Eukaryota" & !Mastertable$Final_superkingdom == "Archaea" & !Mastertable$Final_superkingdom == "Bacteria" & !Mastertable$Diamond_superkingdom == "Other" & !Mastertable$Diamond_superkingdom == "Unannotated" & !Mastertable$Blastn_superkingdom == "Unannotated" & Mastertable$CAT_superkingdom == "Unannotated"]
####################################
# 13.8 Final classification of dark matter
####################################
# The only sequences left are the "unnanotated" sequences, or "dark matter"
Mastertable$Final_superkingdom[Mastertable$Final_superkingdom == "1"] <- "Dark matter"
table(Mastertable$Final_superkingdom)

# Total:
# Bacteria: 40.348 
# Eukaryota: 4.936
# Viruses: 4.834
# Dark matter: 1.732
# Unannotated: 545
# Archaea: 34
####################################

######### !!!!!!!!!!!!!!! #########
    # SAVE Global Environment #
######### !!!!!!!!!!!!!!! #########
