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

# ------------------> PS: I didn't insert the host prediction results in this script yet (script 4). 
####################################
# 0. Packages: install packages and load them wherever needed
####################################
install.packages("ggplot2")
install.packages("dplyr")
install.packages("plyr")
install.packages("reshape")
install.packages("reshape2")
install.packages("scales")
install.packages("viridis")
install.packages("readr")
install.packages("WriteXLS")
install.packages("readxl")
install.packages("devtools")
install_github("jokergoo/ComplexHeatmap")
install.packages("dendsort")
install.packages("vegan")
install.packages("ape")
install.packages("devtools")
install.packages("seriation")
install.packages("taxonomizr")
install.packages("ggThemeAssist")
install.packages("esquisse")
installed.packages("modeldata")
install.packages("tidyverse")

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
## Regenerate accessionTaxa.sql every single time to keep up-to-date. For that go to this directory & delete all the file.
## Go to the place where you install the DB & remove it. Installing it again will take a while (~ 1 hour)
https://cran.r-project.org/web/packages/taxonomizr/vignettes/usage.html

setwd("/Users/daan/Desktop/scripts/Rstudio/DB")
getwd()
dir()

rm(nucl_wgs.accession2taxid.gz)
rm(nucl_gb.accession2taxid.gz)
rm(nodes.dmp)
rm(names.dmp)
rm(accessionTaxa.sql-journal)
rm(accessionTaxa.sql)

prepareDatabase('accessionTaxa.sql')
####################################
# 2. Create abundance table 
####################################
## To make the dataframe, I create abundance files of reads that specifically map (~ meaning mapping of reads in a samples to clusters only comprising nr-contigs of that samples).
## A threshold of horizontal coverage of 70% is set to remove FP (streching a vertical depth of min. 1) over minimum 1 kb nr-contigs. If this threshold is not met, all reads values are converted to a zero.

## As an alternative you can use the 'read.table' function which reads the file and create directly a dataframe. Since we want to merge all the abundance files later, we opted to create lists and merge the lists for all samples by doing it down below

## Some samples have empty text files because all the reads were trimmed away. Remove these ones out of our list, if not, this code won't work.
## removed: 12 samples already that were empty
## B2, B225, B230, B75, B178, B25, B60, B43, B38, B193, B370, B67

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
## Because DIAMOND gives to many contigs, it is way easier to save a vector now that contains all the contigs (clustered files). Now take go to the HPC and adapt the .cltr file by doing the following command:
## I obtain this file by the following line on the HPC (now on Transfer folders), cat *.clstr | grep '>' | cut -f2 > non-redundant-contigs.names.txt

## check empty ones? still in NR-contigs? remove

setwd("/Users/daan/Desktop/Transfer/Multi_fasta_files/")
getwd()

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
# If it is needed, delete the header of the DIAMOND file by: cat *.txt | sed 1d > newname.txt 
# Same for rest of files

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
# 4.4.1 vs1
####################################
# Recently, it changes in which there are only three categories now for both (I) sure phages/prophages, (II) moderately sure phages/prophages and (III) less sure phages/prophages
# Duplication of NDOES are possible - remove them then
# AllSamples_VirsorterCategories <- as.data.frame(read_csv("VIRSorter_global-phage-signal.csv", col_names = FALSE)) 
# colnames(AllSamples_VirsorterCategories) <- c('contig', 'Virsorter')
# AllSamples_VirsorterCategories <- lapply(AllSamples_VirsorterCategories, gsub, pattern = "VIRSorter_", replacement = "")
# AllSamples_VirsorterCategories <- as.data.frame(lapply(AllSamples_VirsorterCategories, gsub, pattern = "-circular", replacement = ""))
# rownames(AllSamples_VirsorterCategories) <- AllSamples_VirsorterCategories$contig
# AllSamples_VirsorterCategories$contig <- NULL

#view(AllSamples_VirsorterCategories)
####################################
# 4.4.2 vs2
####################################
# vs2 create scores for the specific categories (dsDNA, ssDNA, RNA). We will not focus on the other scores (lavidaviridae, virophages NCLDV).
# We will take the category with the highest scores and assume that this is a correct score for later use in the newly-made phage-identification scheme

# Evaluate how to do this better, since I don't know actual input
# I need the group that has the highest score & give name of group

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
# 4.5 Metaphinder2
####################################
# MetaPhinder2_output_filtered <- as.data.frame(read_delim("MetaPhinder2_phage.txt", "\t", escape_double = FALSE, trim_ws = TRUE, col_names = FALSE))
# colnames(MetaPhinder2_output_filtered) <- c('contig', 'MetaPhinder2')
# rownames(MetaPhinder2_output_filtered) <- MetaPhinder2_output_filtered$contig
# MetaPhinder2_output_filtered$contig <- NULL
# rownames(MetaPhinder2_output_filtered) <- gsub("\\.{1}", "_", rownames(MetaPhinder2_output_filtered))

#view(MetaPhinder2_output_filtered)
####################################
# 4.6 DeepVirFinder
####################################
# DeepVirFinder_output_filtered <- as.data.frame(read_delim("DeepVirFinderPhage.txt", "\t", escape_double = FALSE, trim_ws = TRUE, col_names = FALSE))
# colnames(DeepVirFinder_output_filtered) <- c('contig', 'Dvf')
# rownames(DeepVirFinder_output_filtered) <- DeepVirFinder_output_filtered$contig
# DeepVirFinder_output_filtered$contig <- NULL
# rownames(DeepVirFinder_output_filtered) <- gsub("\\.{1}", "_", rownames(DeepVirFinder_output_filtered))

#view(DeepVirFinder_output_filtered)
####################################
# 5. Insert gene prediction and functional annotation
####################################
# 5.1 Genes
####################################
## After the name of the file being "NumberOfGenesPerContigOrdered.txt", the sepearation is written being a space " " in this case. previously this was a tab "\t".

# setwd ("/Users/daan/Desktop/Transfer/Input_R/Genes/") # Check up on this because we have multiple genes.txt file, if correct delete this line
# getwd()
# dir()

# AllSamples_classification_numberofgenes <- as.data.frame(read_delim("GenesPerContig.txt", " ", escape_double = FALSE, trim_ws = TRUE, col_names = FALSE))
# colnames(AllSamples_classification_numberofgenes) <- c('contig', 'Genes')
# rownames(AllSamples_classification_numberofgenes) <- AllSamples_classification_numberofgenes$contig
# AllSamples_classification_numberofgenes$contig <- NULL
# rownames(AllSamples_classification_numberofgenes) <- gsub("\\.{1}", "_", rownames(AllSamples_classification_numberofgenes))

#view(AllSamples_classification_numberofgenes)
####################################
# 5.2 pVOGs
####################################
# setwd ("/Users/daan/Desktop/Transfer/Input_R/Functions/pVOGs/") 
# getwd()
# dir()

# AllSamples_classification_numberofpVOGs <- as.data.frame(read_delim("Phage_VOGs.txt", " ", escape_double = FALSE, trim_ws = TRUE, col_names = FALSE))
# colnames(AllSamples_classification_numberofpVOGs) <- c('contig', 'pVOGs')
# rownames(AllSamples_classification_numberofpVOGs) <- AllSamples_classification_numberofpVOGs$contig
# AllSamples_classification_numberofpVOGs$contig <- NULL
# rownames(AllSamples_classification_numberofpVOGs) <- gsub("\\.{1}", "_", rownames(AllSamples_classification_numberofpVOGs))

#view(AllSamples_classification_numberofpVOGs)
####################################
# 5.3 EggNOGs
####################################
# setwd ("/Users/daan/Desktop/Transfer/Input_R/Functions/EggNOGs/") 
# getwd()
# dir()

# AllSamples_classification_numberofEggNOGs <- as.data.frame(read_delim("ViralEggNOGs.txt", " ", escape_double = FALSE, trim_ws = TRUE, col_names = FALSE))
# colnames(AllSamples_classification_numberofEggNOGs) <- c('contig', 'EggNOGs')
# rownames(AllSamples_classification_numberofEggNOGs) <- AllSamples_classification_numberofEggNOGs$contig
# AllSamples_classification_numberofEggNOGs$contig <- NULL
# rownames(AllSamples_classification_numberofEggNOGs) <- gsub("\\.{1}", "_", rownames(AllSamples_classification_numberofEggNOGs))

#View(AllSamples_classification_numberofEggNOGs)
####################################
# 5.4 InterProScan: pfams
####################################
# setwd ("/Users/daan/Desktop/Transfer/Input_R/Functions/InterProScan/") 
# getwd()
# dir()

# AllSamples_classification_numberofPfams <- as.data.frame(read_delim("Viral_Pfams.txt", " ", escape_double = FALSE, trim_ws = TRUE, col_names = FALSE))
# colnames(AllSamples_classification_numberofPfams) <- c('contig', 'Pfams')
# rownames(AllSamples_classification_numberofPfams) <- AllSamples_classification_numberofPfams$contig
# AllSamples_classification_numberofPfams$contig <- NULL
# rownames(AllSamples_classification_numberofPfams) <- gsub("\\.{1}", "_", rownames(AllSamples_classification_numberofPfams))

#view(AllSamples_classification_numberofPfams)
####################################
# 5.5 Cenotetaker2
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
#view(AllSamples_Cenotetaker2_hallmark_proteins)
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
#view(AllSamples_Cenotetaker2_lysogenic_phages)
####################################
# 5.5.3 Toxin-containing contigs
####################################
# re-evaluate if I don't miss any "words" for toxins as different names in orginal file
AllSamples_Cenotetaker2_toxin_phages <- as.data.frame(read_delim("All_predicted_toxin_phages_nodes.txt", "\t", escape_double = TRUE, trim_ws = TRUE, col_names = FALSE))
colnames(AllSamples_Cenotetaker2_toxin_phages) <- c('contig', 'Toxin presence', 'Toxin type', 'Toxin TA')
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
view(AllSamples_Cenotetaker2_recombination_phages)
####################################
# 5.5.5 RT-containing contigs
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
# view(AllSamples_Cenotetaker2_Bacteriocin_phages)
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
#View(AllSamples_classificationtaxonomy_Cenotaker2)
####################################
# 5.6 Crass phage taxonomy
####################################
# Phages that can be identified as CrAssphage.
# Comparison against ~ 1000 CrAssphage validated genomes in big public datasets (Yunnin, Guerin)
# Based on BLASTN with E-value of 1E-05 & introduce coverage length of 3kb to remove smallest crap, although if you look at data (taxonomy folder) most of the contigs will have a very high contig length (10's of kbs)
# Also tried: TBLASTN doesn't improve it.
# Also tried aligning own phages against putative CrAss halmmark protein (terminase & Polymerase), which only give data that LCA approach also gives.

AllSamples_classificationtaxonomy_Crass <- as.data.frame(read_delim("phages_crass_predictions1.txt", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE))
AllSamples_classificationtaxonomy_Crass <- AllSamples_classificationtaxonomy_Crass[,c(1,2)]
colnames(AllSamples_classificationtaxonomy_Crass) <- c('contig', 'CrAss_taxonomy')
rownames(AllSamples_classificationtaxonomy_Crass) <- AllSamples_classificationtaxonomy_Crass$contig
AllSamples_classificationtaxonomy_Crass$contig <- NULL
rownames(AllSamples_classificationtaxonomy_Crass) <- gsub("\\.{1}", "_", rownames(AllSamples_classificationtaxonomy_Crass))
#View(AllSamples_classificationtaxonomy_Crass)
####################################
# 6. Insert circularity
####################################
setwd ("/Users/daan/Desktop/Transfer/Input_R/CheckV_circularity/") 
# getwd()
# dir()

AllSamples_classification_Circularity <- as.data.frame(read_delim("CircularContigs.txt", "\t", escape_double = FALSE, trim_ws = TRUE, col_names = FALSE))
colnames(AllSamples_classification_Circularity) <- c('contig', 'circular')
rownames(AllSamples_classification_Circularity) <- AllSamples_classification_Circularity$contig
AllSamples_classification_Circularity$contig <- NULL
rownames(AllSamples_classification_Circularity) <- gsub("\\.{1}", "_", rownames(AllSamples_classification_Circularity))

#view(AllSamples_classification_Circularity)
####################################
# 7. Insert Completeness
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
# 8. Create mastertable: merge scaffolds, abundances and annotations
####################################
# Still have to insert host prediction later

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

Mastertable <- merge(COMBO_12, AllSamples_classificationtaxonomy_Crass, by=0, all=T)
rownames(Mastertable) <- Mastertable$Row.names
Mastertable$Row.names <- NULL
dim(Mastertable)

#View(Mastertable)
####################################
# 9. Additional information to Mastertable
####################################
# 9.1 Insert length of contigs
####################################
Mastertable$length <- rownames(Mastertable)
Mastertable$length <- gsub("*_cov.*","",Mastertable$length)
Mastertable$length <- gsub(".*length_","", Mastertable$length)
Mastertable$length <- as.numeric(Mastertable$length)
Mastertable$length <- Mastertable$length/1000
# View(Mastertable)
####################################
# 9.2 Insert gene length (kb) ratio of contigs
####################################
Mastertable$GeneLengthRatio <- (Mastertable$`total genes`/Mastertable$length)
# view(Mastertable)
####################################
# 9.3 Insert total gene to viral gene ratio of contigs (~ estimate of how many genes of the total are viral)
####################################
Mastertable$ViralGeneLengthRatio <- (Mastertable$`viral genes`/Mastertable$length)
####################################
# 9.4 Insert total number of trimmed reads per contig; and additional information
####################################
# 9.4.1 Total Nr. of reads
####################################
Mastertable$Totalnumberofreads <- as.numeric(rowSums(Mastertable[,colnames(Mastertable) %in% names]))
Mastertable <- Mastertable[!(Mastertable$Totalnumberofreads == "0"),] ## remove rows with total number of zero
# From here on you could determine general things like: (I) longest phage, (II) most abundant phage and (III) most present phage
# (I) 'length' column
# (II) 'Totalnumberofreads' column
# (III) create new count: see down below (~ or make pressence/absence table)

#Mastertable1 <- Mastertable
#Mastertable1 <- Mastertable1[, names] # subset all columns based on name vector
#Mastertable1[Mastertable1 > 0] <- 1
#Mastertable1$mostpresentphage <- rowSums(Mastertable1[,colnames(Mastertable1) %in% names])
#view(Mastertable1)

# (1) Repeat this one we identified all viruses with phage-identification scheme
# (2) Repeat this only for high-quality phage contigs
# (3) Put all these data quickly in one powerpoint & copy these line to the place in the script where I will actually perform these lines
####################################
# 9.4.2 Additional information
####################################
#Mastertable_contigs <- Mastertable

# This mastertable "Mastertable_contigs" is created for another script later. For now neglect.
####################################
# 9.5 Removal of low-abundant contigs (~ vertical reads coverage correction = 10X) based on reads per kb
####################################
# 9.5.1 Create vector for parts of dataframe
####################################
#vector_1 <- (which(names(Mastertable)== "Virsorter"))
#vector_1_1 <- (which(names(Mastertable)== "Virsorter")-1) 
#vector_2 <- which(names(Mastertable)== "ViralGeneLengthRatio")
#columns_to_merge <- Mastertable %>% select(vector_1:vector_2)
# head(columns_to_merge)
# head(Mastertable)

# PS: I use vector_1_1 to make this applicable for all datasets with different names yet to come.
####################################
# 9.5.2 Subset abundance table & correct for vertical depth (10X = 63 reads/kb)
####################################
#Mastertable_kb <- (Mastertable[,1:vector_1_1]/Mastertable$length)
#Mastertable_kb[Mastertable_kb < 78] <- 0
####################################
# 9.5.3 Add contig length & convert to abundance table again
####################################
#Mastertable_kb$length <- rownames(Mastertable_kb)
#Mastertable_kb$length <- gsub("*_cov.*","",Mastertable_kb$length)
#Mastertable_kb$length <- gsub(".*length_","", Mastertable_kb$length)
#Mastertable_kb$length <- as.numeric(Mastertable_kb$length)
#Mastertable_kb$length <- Mastertable_kb$length/1000
#Mastertable_corrected <- Mastertable_kb[,1:vector_1_1]*Mastertable_kb$length

# PS: the reason why we can't keep the reads per kb in the abundance table is because we know the assemblers used today create chimeras
# Therefore, some viruses have double, triple, or even more the length of the normal genome & therefore you can't correct for it because the length is often not correct.
# Also we know from experience that it is skewed towards small circular contigs by using WTA2 giving worse results after this correction.
####################################
# 9.5.4 Merge into depth correct mastertable
####################################
#Mastertable_corrected <- merge(Mastertable_corrected, columns_to_merge, by = 0, all = T)
#rownames(Mastertable_corrected) <- Mastertable_corrected$Row.names
#Mastertable_corrected$Row.names <- NULL
####################################
# 9.5.5 Add total nr. of reads & total nr. reads/kb
####################################
#Mastertable_corrected$Totalnumberofreads <- as.numeric(rowSums(Mastertable_corrected[,colnames(Mastertable_corrected) %in% names]))
#Mastertable_corrected <- Mastertable_corrected[!(Mastertable_corrected$Totalnumberofreads == "0"),] 

#Mastertable_corrected$Totalnumberofreadsperkb <- as.numeric(Mastertable_corrected$Totalnumberofreads/Mastertable_corrected$length)
####################################
# 9.5.6 Calculation of %contigs & %reads retained 
####################################
# how much low-abundant contigs with < 10X vertical depth thrown away
#(nrow(Mastertable_corrected)/nrow(Mastertable))*100 # Number of contigs retained after correction step
#(sum(Mastertable_corrected$Totalnumberofreads)/sum(Mastertable$Totalnumberofreads))*100 # Number of reads retained after correction step

## 94.66972% of contigs is kept
## 99,33058 of reads is kept

# View(Mastertable_corrected)
####################################
# 9.5.7 Clean-up environment 
####################################
## Otherwise Rstudio will start slowing down after a while.
## Do the following command to create a variable with all the stuff in your environment to evlauate what you can remove.

# Clean-up environment
# MyVariables<-objects()
# for (i in MyVariables) {
#  print(i)
# }

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

## The remove of low-abundant contigs is de facto setting a threshold for vertical read coverage. This is because we already insert bioinformatically a horizontal reads coverage of 70%.
## It means that assuming reads have a max. length of 127 bp and a min. length of 50 bp, and contigs are min. 1 kb, every contig has a minimal horizontal coverage of (1000 bp/ 127 bp)*0,7 = 5,52 reads/kb or 6 reads per kb.
## This means that minimally 5,52 reads are covering every 1 kb contig.

## Since there are also a lot of contigs that are fully covered we would like to make the threshold more conservative. Therefore, we calculate assuming the entire contig is covered with reads how many reads would minimally cover this contig.
## 1000 bp/127 = 7,78 reads are minimally covering a 1 kb contig with a depth of 1X.

## Since we want to be more conservative, I would propose a vertical depth threshold of 10X vertical read coverage.
## (1000bp contig/127bp reads)* 10X vertical read depth = 77,8 (~  78 reads per kb 100% 10X covering the contig)

## We propose that a contig covered by a minimum of 78 reads per kb is a legit threshold for vertical depth
## And at the same time to correct for low-abundant contigs.
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
## Quite some errors will occur, because data is very large. In this case press ESC so error disappaears and check output.
## Normally not made for such a huge data.
## Quite some unannotated data too.

setwd("/Users/daan/Desktop/scripts/Rstudio/DB")
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
## Convert all the other viral contigs to prokaryotic viruses, and create a list so you can manually check if contigs were miss-classified.
## Eukaryotic viruses are generally very well annotated. If it is annotated as a virus, but not an eukaryotic virus, the overwhelming rest majority will be a phage. Therefore, we will annotate it for now as such & also because the phage identification scheme later will only make a distinction between eukaryotic viruses & other viral elements, so this has no impact later.
## There are a few negligible contigs that only have a viral annotation, but are not denoted on a deeper hierarchical level. These you can check via the "Subset_Blastn_fages".

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
setwd ("/Users/daan/Desktop/Transfer/Input_R/Taxonomy/")
# getwd()
# dir()
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

#View(Blastn_taxonomy)
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
setwd("/Users/daan/Desktop/scripts/Rstudio/DB")
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
setwd ("/Users/daan/Desktop/Transfer/Input_R/Taxonomy/")
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
#View(Diamond_taxonomy)
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
setwd("/Users/daan/Desktop/scripts/Rstudio/DB")
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
#View(CAT_taxonomy)
####################################
# 10.4.6 Add ORFs
####################################
setwd ("/Users/daan/Desktop/Transfer/Input_R/Taxonomy/")
# getwd()
# dir()
####################################
# 10.4.6.1 Add features
####################################
CAT_number_of_ORFs_for_annotation_LCA <- as.data.frame(read_delim("CAT_alignment_ORFs.txt", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE))
colnames(CAT_number_of_ORFs_for_annotation_LCA) <- c('NODES', 'CAT_ORFs')
CAT_number_of_ORFs_for_annotation_LCA <- CAT_number_of_ORFs_for_annotation_LCA[CAT_number_of_ORFs_for_annotation_LCA$NODES %in% contigs_project,]
CAT_number_of_ORFs_for_annotation_LCA$NODES <- gsub("\\.{1}", "_", CAT_number_of_ORFs_for_annotation_LCA$NODES)
rownames(CAT_number_of_ORFs_for_annotation_LCA) <- CAT_number_of_ORFs_for_annotation_LCA$NODES
CAT_number_of_ORFs_for_annotation_LCA$NODES <- NULL
####################################
# 10.4.6.2 Merge
####################################
CAT_taxonomy <- transform(merge(CAT_taxonomy, CAT_number_of_ORFs_for_annotation_LCA, by = 0, all = T), row.names=Row.names, Row.names = NULL)
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
# 10.4.7 Clean-up environment
####################################
rm(contigs_CAT)
rm(Subset_CAT_fages)
rm(CAT_number_of_ORFs_for_annotation_LCA)
####################################
# 10.5 CNoteTaker2 annotation
####################################
# A lot of the phages that are annotated with blastx, blastn with not be annotated using the LCA approach
# To meet this limitations in the middle we will add annotation of identified phages (~ later stage), to which the viruses that are still Unannotated after LCA blastx; LCA blastn, will be annotated by CNoteTaker2 blastx if this is present & if a hallmark protein is present.
# In practice, I saw that often viruses were known but were Unannotated because of a crappy bacterial contig still present in LCA leading to an insane loss of annotation.
# Therefore, if: (i) phage identified, (ii) hallmark protein present, (ii) blastx annotation:
####################################
# 10.5.1 Add taxonomies
####################################
setwd("/Users/daan/Desktop/scripts/Rstudio/DB")
getwd()
dir()

contigs_cenote <- rownames(AllSamples_classificationtaxonomy_Cenotaker2)
cenote_taxonomy <- as.data.frame(getTaxonomy(ids = AllSamples_classificationtaxonomy_Cenotaker2$Cenote_single_TaxIdpath, 'accessionTaxa.sql', desiredTaxa = c("superkingdom", "phylum", "class", "order", "family")))
rownames(cenote_taxonomy) <- contigs_cenote
colnames(cenote_taxonomy) <- paste0('Cenote_', colnames(cenote_taxonomy))
cenote_taxonomy <- na.omit(cenote_taxonomy)                          
cenote_taxonomy[is.na(cenote_taxonomy)] <- "Unannotated"
# View(cenote_taxonomy)
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
#table(cenote_taxonomy$Cenote_family)
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
#View(Mastertable[Mastertable$Totalnumberofreads > 10000000,])

# Clean environment
rm(Mastertable_corrected)
rm(Blastn_taxonomy)
rm(Diamond_taxonomy)
rm(CAT_taxonomy)

# PS: hold in midn that CrassPhages are now busy with being placed in a new order called "Crassvirales".
# it might make sense to make a list later with all the Crass-like names (~ or as Crass*) and store them in a variable so you can make a distinction later with other Caudovirales fraction, or compare these as such.
####################################
# 11. Phage-identification 
####################################
## The Phage identification is based on a combination of two tools: (i) Virsorter2 and (ii) CheckV:

## We identify phages as:
## virsorter2: dsDNAphages; CheckV > 50% complete
## virsorter2: RNA phages; CheckV > 50% complete
## virsorter2: ssDNA phages; length >= 3kb; at leasr 1 viral gene predicted by CheckV

Mastertable$Phage_score <- 1
Mastertable$Phage_score[Mastertable$Virsorter == "dsDNAphage" & Mastertable$completeness >= 50] <- "dsDNA phage"
Mastertable$Phage_score[Mastertable$Virsorter == "RNA" & Mastertable$completeness >= 50] <- "RNA phage"
Mastertable$Phage_score[Mastertable$Virsorter == "ssDNA" & Mastertable$length >= 3 & Mastertable$completeness >= 50] <- "ssDNA phage"
table(Mastertable$`viral genes`[Mastertable$Phage_score == "ssDNA phage"])
#Mastertable$Phage_score[Mastertable$Phage_score == "ssDNA phage" & Mastertable$`viral genes` >= 1] <- "Real ssDNA phage"
#table(Mastertable$Phage_score)

# How many phages are identified?
table(Mastertable$Phage_score)

# dsDNA phages; 2850 near-complete phages (also lots of ssDNA phages in here)
# RNA phages; 15 near-complete phages
# ssDNA phages; 70 near-complete phages

Mastertable$Phages <- 1
Mastertable$Phages[Mastertable$Phage_score == "dsDNA phage" | Mastertable$Phage_score == "ssDNA phage" | Mastertable$Phage_score == "RNA phage"] <- "phage"
table(Mastertable$Phages) 
nrow(Mastertable)
# total: 2935 phages identified using this method
# total NR-contigs = 110.330 contigs

Mastertable$Totalnumberofreads <- as.numeric(Mastertable$Totalnumberofreads)
sum(Mastertable$Totalnumberofreads[Mastertable$Phages == "phage"]) # total of 2662662144 reads (~ or 2,6B PE reads)
####################################
# 11.5.4 Clean-up environment
####################################
# We will remove this to keep the table as least messy as possible. The only thing left is the 'Phage_score' given us the identification result. Removing data we won't need anymore will also speed things up a lot.
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
# Since we took the LCA as annotation for all homology-based methods we don't have alignment features for our taxonomy. However, I added before the alignments features (% similarity and % cov of query) of the the best alignment found by Blastn and Diamond. The program of CAT doesn't give this in output, but does allow us to use the number of ORFs the program uses in its scoring-scheme.
# I will not use those data about the quality of alignments to assign taxonomy. I will just thrive on the E-value that I put as a cutoff on the HPC & trust alignment given by the blast as output
# However, I could use this data to evaluate if the phages are unknown are can be classified at a certain hierarchy. For example genera level hierarchy according to ICTV is 70% ANI over full-length genomes.

# Thus, depending on the similarity to the sequence in the database you can or can not use the annotation that is given. For example 70% ANI over 100%coverage is given to use genus level taxonomy for phages in general. This is something you could use as a cutoff for annotation for Blastn, and then extrapolate to proteins for Diamond. If we deal with 70% ANI, the %API should be higher ofcourse since it is less variable.

# CAT is a different tools and predicts genes first, and then performs a local alignment using diamonds adjusted blastx method, then it used those alignments and takes a LCA approach to taxonomically annotate the entire contig.
# The problem is that some annotations will be performed on 3 genes & some on 10, which will depend on the entire length. We take this into account by taking a minimal of 3 ORFs predicted & a minimal of 1 gene/kb.
####################################
# 12.2.1 Viral
####################################
Mastertable$Final_viral[Mastertable$Phages == "phage" & !Mastertable$Diamond_viral == "Diamond_eukaryotic_virus" & !Mastertable$Blastn_viral == "Blastn_eukaryotic_virus"] <- Mastertable$Phages[Mastertable$Phages == "phage" & !Mastertable$Diamond_viral == "Diamond_eukaryotic_virus" & !Mastertable$Blastn_viral == "Blastn_eukaryotic_virus"]
table(Mastertable$Phage_score)
table(Mastertable$Final_viral)

sum(Mastertable$Totalnumberofreads[Mastertable$Phages == "phage"])
sum(Mastertable$Totalnumberofreads[Mastertable$Phage_score == "Real ssDNA phage"])
sum(Mastertable$Totalnumberofreads[Mastertable$Phage_score == "ssDNA phage"])
sum(Mastertable$Totalnumberofreads[Mastertable$Phage_score == "dsDNA phage"])
sum(Mastertable$Totalnumberofreads[Mastertable$Phage_score == "RNA phage"])

## Most of them are long dsDNA phages, but we know very often these are wrongly annotated as dsDNA and are ssDNA
## so check also the length ( < 10kb)
sum(Mastertable$Totalnumberofreads[Mastertable$Phages == "phage" & Mastertable$length <= 10])
sum(Mastertable$Totalnumberofreads[Mastertable$Phages == "phage" & Mastertable$length > 10])

#View(Mastertable[Mastertable$Phages == "phage",])
#View(Mastertable[Mastertable$Phage_score == "ssDNA phage",])
#View(Mastertable[Mastertable$Phage_score == "Real ssDNA phage",])

## Now, we see a more nuanced image with:
# ~ 865M phages being likely micro's (short ones); with the accurate micro's being 
# ~ 1,74B phages being likely caudo's (long ones)
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
####################################
# 12.9 Clean-up environment
####################################
#Mastertable$Blastn_superkingdom <- NULL
#Mastertable$Blastn_phylum <- NULL
#Mastertable$Blastn_class <- NULL
#Mastertable$Blastn_order <- NULL
#Mastertable$Blastn_family <- NULL
#Mastertable$Blastn_subfamily <- NULL
#Mastertable$Blastn_genus <- NULL
#Mastertable$Blastn_species<- NULL

#Mastertable$Diamond_superkingdom <- NULL
#Mastertable$Diamond_phylum <- NULL
#Mastertable$Diamond_class <- NULL
#Mastertable$Diamond_order <- NULL
#Mastertable$Diamond_family <- NULL
#Mastertable$Diamond_subfamily <- NULL
#Mastertable$Diamond_genus <- NULL
#Mastertable$Diamond_species<- NULL

#Mastertable$CAT_superkingdom <- NULL
#Mastertable$CAT_phylum <- NULL
#Mastertable$CAT_class <- NULL
#Mastertable$CAT_order <- NULL
#Mastertable$CAT_family <- NULL
#Mastertable$CAT_subfamily <- NULL
#Mastertable$CAT_genus <- NULL
#Mastertable$CAT_species<- NULL
####################################
# 13. Additional Cenote-Taker2 Taxonomical information: 1 hallmark gene
####################################
# Only add taxonomies if hallmark gene is present
Mastertable$Final_phylum[Mastertable$Final_viral == "phage" & Mastertable$Final_phylum == "Unannotated" & Mastertable$`Hallmark proteins` >= 1 & Mastertable$Cenote_viral == "cenote_prokaryotic_virus"] <- Mastertable$Cenote_phylum[Mastertable$Final_viral == "phage" & Mastertable$Final_phylum == "Unannotated" & Mastertable$`Hallmark proteins` >= 1 & Mastertable$Cenote_viral == "cenote_prokaryotic_virus"]
Mastertable$Final_class[Mastertable$Final_viral == "phage" & Mastertable$Final_class == "Unannotated" & Mastertable$`Hallmark proteins` >= 1 & Mastertable$Cenote_viral == "cenote_prokaryotic_virus"] <- Mastertable$Cenote_class[Mastertable$Final_viral == "phage" & Mastertable$Final_class == "Unannotated" & Mastertable$`Hallmark proteins` >= 1 & Mastertable$Cenote_viral == "cenote_prokaryotic_virus"]
Mastertable$Final_order[Mastertable$Final_viral == "phage" & Mastertable$Final_order == "Unannotated" & Mastertable$`Hallmark proteins` >= 1 & Mastertable$Cenote_viral == "cenote_prokaryotic_virus"] <- Mastertable$Cenote_order[Mastertable$Final_viral == "phage" & Mastertable$Final_order == "Unannotated" & Mastertable$`Hallmark proteins` >= 1 & Mastertable$Cenote_viral == "cenote_prokaryotic_virus"]
Mastertable$Final_family[Mastertable$Final_viral == "phage" & Mastertable$Final_family == "Unannotated" & Mastertable$`Hallmark proteins` >= 1 & Mastertable$Cenote_viral == "cenote_prokaryotic_virus"] <- Mastertable$Cenote_family[Mastertable$Final_viral == "phage" & Mastertable$Final_family == "Unannotated" & Mastertable$`Hallmark proteins` >= 1 & Mastertable$Cenote_viral == "cenote_prokaryotic_virus"]
# View(Mastertable[Mastertable$Final_class == "Malgrandaviricetes",])

table(Mastertable$Final_phylum)
table(Mastertable$Final_class)
table(Mastertable$Final_order)
table(Mastertable$Final_family)
table(Mastertable$Final_viral)
####################################
# 14. Additional information: Not-determined category
####################################
# 14.1 Evaluation of highly-abundant contigs
####################################
# Evaluate the non-determined fraction potentially holding known obvious viruses not identified by our method

# View(Mastertable[Mastertable$completeness_quality == "Not-determined",])
####################################
# 14.2 Identify and annotated Inoviridae and Picobirnaviridae
####################################
# If 2 out of 3 homology based tools give are identical on family level than identifiy the contig as phage & annotate as Diamond, Blastn, CAT

# 13.2.1 Inoviridae
Mastertable$Phages[Mastertable$completeness_quality == "Not-determined" & Mastertable$length > 3 & (Mastertable$Blastn_family == "Inoviridae" & Mastertable$Diamond_family == "Inoviridae" | Mastertable$Blastn_family == "Inoviridae" & Mastertable$CAT_family == "Inoviridae" | Mastertable$CAT_family == "Inoviridae" & Mastertable$Diamond_family == "Inoviridae")] <- "phage"
Mastertable$Final_viral[Mastertable$completeness_quality == "Not-determined"  & Mastertable$length > 3 & (Mastertable$Blastn_family == "Inoviridae" & Mastertable$Diamond_family == "Inoviridae" | Mastertable$Blastn_family == "Inoviridae" & Mastertable$CAT_family == "Inoviridae" | Mastertable$CAT_family == "Inoviridae" & Mastertable$Diamond_family == "Inoviridae")] <- "phage"
Mastertable$Final_superkingdom[Mastertable$completeness_quality == "Not-determined"  & Mastertable$length > 3 & (Mastertable$Blastn_family == "Inoviridae" & Mastertable$Diamond_family == "Inoviridae" | Mastertable$Blastn_family == "Inoviridae" & Mastertable$CAT_family == "Inoviridae" | Mastertable$CAT_family == "Inoviridae" & Mastertable$Diamond_family == "Inoviridae")] <- "Viruses"
Mastertable$Final_phylum[Mastertable$completeness_quality == "Not-determined"  & Mastertable$length > 3 & (Mastertable$Blastn_family == "Inoviridae" & Mastertable$Diamond_family == "Inoviridae" | Mastertable$Blastn_family == "Inoviridae" & Mastertable$CAT_family == "Inoviridae" | Mastertable$CAT_family == "Inoviridae" & Mastertable$Diamond_family == "Inoviridae")] <- "Hofneiviricota"
Mastertable$Final_class[Mastertable$completeness_quality == "Not-determined"  & Mastertable$length > 3 & (Mastertable$Blastn_family == "Inoviridae" & Mastertable$Diamond_family == "Inoviridae" | Mastertable$Blastn_family == "Inoviridae" & Mastertable$CAT_family == "Inoviridae" | Mastertable$CAT_family == "Inoviridae" & Mastertable$Diamond_family == "Inoviridae")] <- "Faserviricetes"
Mastertable$Final_order[Mastertable$completeness_quality == "Not-determined"  & Mastertable$length > 3 & (Mastertable$Blastn_family == "Inoviridae" & Mastertable$Diamond_family == "Inoviridae" | Mastertable$Blastn_family == "Inoviridae" & Mastertable$CAT_family == "Inoviridae" | Mastertable$CAT_family == "Inoviridae" & Mastertable$Diamond_family == "Inoviridae")] <- "Tubulavirales"
Mastertable$Final_family[Mastertable$completeness_quality == "Not-determined"  & Mastertable$length > 3 & (Mastertable$Blastn_family == "Inoviridae" & Mastertable$Diamond_family == "Inoviridae" | Mastertable$Blastn_family == "Inoviridae" & Mastertable$CAT_family == "Inoviridae" | Mastertable$CAT_family == "Inoviridae" & Mastertable$Diamond_family == "Inoviridae")] <- "Inoviridae"
# 21 Not-determined NR-contigs to which 29,125,703M PE reads are mapped (pre-rarefied) (total 50 NR-contig; 60M PE reads)

# 13.2.2 Adjust wrong annotations  
Mastertable$Final_subfamily[Mastertable$completeness_quality == "Not-determined"  & Mastertable$length > 3 & (Mastertable$Blastn_family == "Inoviridae" & Mastertable$Diamond_family == "Inoviridae" | Mastertable$Blastn_family == "Inoviridae" & Mastertable$CAT_family == "Inoviridae" | Mastertable$CAT_family == "Inoviridae" & Mastertable$Diamond_family == "Inoviridae")] <- "1"
Mastertable$Final_species[Mastertable$completeness_quality == "Not-determined"  & Mastertable$length > 3 & (Mastertable$Blastn_family == "Inoviridae" & Mastertable$Diamond_family == "Inoviridae" | Mastertable$Blastn_family == "Inoviridae" & Mastertable$CAT_family == "Inoviridae" | Mastertable$CAT_family == "Inoviridae" & Mastertable$Diamond_family == "Inoviridae")] <- "1"
Mastertable$Final_genus[Mastertable$completeness_quality == "Not-determined"  & Mastertable$length > 3 & (Mastertable$Blastn_family == "Inoviridae" & Mastertable$Diamond_family == "Inoviridae" | Mastertable$Blastn_family == "Inoviridae" & Mastertable$CAT_family == "Inoviridae" | Mastertable$CAT_family == "Inoviridae" & Mastertable$Diamond_family == "Inoviridae")] <- "1"
Mastertable$Final_ANI[Mastertable$completeness_quality == "Not-determined"  & Mastertable$length > 3 & (Mastertable$Blastn_family == "Inoviridae" & Mastertable$Diamond_family == "Inoviridae" | Mastertable$Blastn_family == "Inoviridae" & Mastertable$CAT_family == "Inoviridae" | Mastertable$CAT_family == "Inoviridae" & Mastertable$Diamond_family == "Inoviridae")] <- "1"
Mastertable$Final_coverage[Mastertable$completeness_quality == "Not-determined"  & Mastertable$length > 3 & (Mastertable$Blastn_family == "Inoviridae" & Mastertable$Diamond_family == "Inoviridae" | Mastertable$Blastn_family == "Inoviridae" & Mastertable$CAT_family == "Inoviridae" | Mastertable$CAT_family == "Inoviridae" & Mastertable$Diamond_family == "Inoviridae")] <- "1"

# 13.2.3 Picobirnaviridae
Mastertable$Phages[Mastertable$completeness_quality == "Not-determined" & (Mastertable$Blastn_family == "Picobirnaviridae" & Mastertable$Diamond_family == "Picobirnaviridae" | Mastertable$Blastn_family == "Picobirnaviridae" & Mastertable$CAT_family == "Picobirnaviridae" | Mastertable$CAT_family == "Picobirnaviridae" & Mastertable$Diamond_family == "Picobirnaviridae")] <- "phage"
Mastertable$Final_viral[Mastertable$completeness_quality == "Not-determined" & (Mastertable$Blastn_family == "Picobirnaviridae" & Mastertable$Diamond_family == "Picobirnaviridae" | Mastertable$Blastn_family == "Picobirnaviridae" & Mastertable$CAT_family == "Picobirnaviridae" | Mastertable$CAT_family == "Picobirnaviridae" & Mastertable$Diamond_family == "Picobirnaviridae")] <- "phage"
Mastertable$Final_superkingdom[Mastertable$completeness_quality == "Not-determined" & (Mastertable$Blastn_family == "Picobirnaviridae" & Mastertable$Diamond_family == "Picobirnaviridae" | Mastertable$Blastn_family == "Picobirnaviridae" & Mastertable$CAT_family == "Picobirnaviridae" | Mastertable$CAT_family == "Picobirnaviridae" & Mastertable$Diamond_family == "Picobirnaviridae")] <- "Viruses"
Mastertable$Final_phylum[Mastertable$completeness_quality == "Not-determined" & (Mastertable$Blastn_family == "Picobirnaviridae" & Mastertable$Diamond_family == "Picobirnaviridae" | Mastertable$Blastn_family == "Picobirnaviridae" & Mastertable$CAT_family == "Picobirnaviridae" | Mastertable$CAT_family == "Picobirnaviridae" & Mastertable$Diamond_family == "Picobirnaviridae")] <- "Pisuviricota"
Mastertable$Final_class[Mastertable$completeness_quality == "Not-determined" & (Mastertable$Blastn_family == "Picobirnaviridae" & Mastertable$Diamond_family == "Picobirnaviridae" | Mastertable$Blastn_family == "Picobirnaviridae" & Mastertable$CAT_family == "Picobirnaviridae" | Mastertable$CAT_family == "Picobirnaviridae" & Mastertable$Diamond_family == "Picobirnaviridae")] <- "Duplopiviricetes"
Mastertable$Final_order[Mastertable$completeness_quality == "Not-determined" & (Mastertable$Blastn_family == "Picobirnaviridae" & Mastertable$Diamond_family == "Picobirnaviridae" | Mastertable$Blastn_family == "Picobirnaviridae" & Mastertable$CAT_family == "Picobirnaviridae" | Mastertable$CAT_family == "Picobirnaviridae" & Mastertable$Diamond_family == "Picobirnaviridae")] <- "Durnavirales"
Mastertable$Final_family[Mastertable$completeness_quality == "Not-determined" & (Mastertable$Blastn_family == "Picobirnaviridae" & Mastertable$Diamond_family == "Picobirnaviridae" | Mastertable$Blastn_family == "Picobirnaviridae" & Mastertable$CAT_family == "Picobirnaviridae" | Mastertable$CAT_family == "Picobirnaviridae" & Mastertable$Diamond_family == "Picobirnaviridae")] <- "Picobirnaviridae"
# 66 Not-determined NR-contigs to whuch 1,211,197M PE reads are mapped (pre-rarefied).

table(Mastertable$Phages)
table(Mastertable$Final_viral) # prokaryotic versus eukaryotic virus
table(Mastertable$Final_superkingdom) # total nr. viruses identified
table(Mastertable$Final_phylum)
table(Mastertable$Final_class)
table(Mastertable$Final_family)
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

# Print as text file
write.table(vector_3, file = "cenote.txt", sep = "\t",row.names = FALSE)

## Now run CenoteTaker2 on those phages and merge the protein, temperateand toxin data with existing files, and run this entire script again.
####################################
# 15. Manual curation of phage fracion
####################################
# 15.1 Evaluate clearly wrong taxonomies
####################################
table(Mastertable$Phages)
table(Mastertable$Final_class[Mastertable$Final_viral == "phage"])

## Megaviricetes, wrongly annotated Caudoviricetes
Mastertable$Final_phylum[Mastertable$Final_viral == "phage" & Mastertable$Final_class == "Megaviricetes"] <- "Uroviricota"
Mastertable$Final_class[Mastertable$Final_viral == "phage" & Mastertable$Final_class == "Megaviricetes"] <- "Caudoviricetes"
Mastertable$Final_order[Mastertable$Final_viral == "phage" & Mastertable$Final_class == "Megaviricetes"] <- "Caudovirales"
Mastertable$Final_family[Mastertable$Final_viral == "phage" & Mastertable$Final_class == "Megaviricetes"] <- "Podoviridae"

Mastertable$Final_class[Mastertable$Final_viral == "phage" & Mastertable$Final_family == "Mimiviridae"] <- "Caudoviricetes"
Mastertable$Final_family[Mastertable$Final_viral == "phage" & Mastertable$Final_family == "Mimiviridae"] <- "Unannotated"

Mastertable$Final_genus[Mastertable$Final_viral == "phage" & Mastertable$Final_order == "Imitervirales"] <- "Unannotated"
Mastertable$Final_species[Mastertable$Final_viral == "phage" & Mastertable$Final_order == "Imitervirales"] <- "Unannotated"
Mastertable$Final_order[Mastertable$Final_viral == "phage" & Mastertable$Final_order == "Imitervirales"] <- "Caudovirales"

## Tectiliviricetes, wrongly annotated Caudoviricetes 
Mastertable$Final_phylum[Mastertable$Final_viral == "phage" & Mastertable$Final_class == "Tectiliviricetes"] <- "Uroviricota"
Mastertable$Final_order[Mastertable$Final_viral == "phage" & Mastertable$Final_family == "Corticoviridae"] <- "Caudovirales"
Mastertable$Final_class[Mastertable$Final_viral == "phage" & Mastertable$Final_family == "Corticoviridae"] <- "Caudoviricetes"
Mastertable$Final_family[Mastertable$Final_viral == "phage" & Mastertable$Final_family == "Corticoviridae"] <- "Unannotated"

table(Mastertable$Phages)
table(Mastertable$Final_class[Mastertable$Final_viral == "phage"])
####################################
# 15.2 Evaluate phages concatamers
####################################
# Phages with a length that is significantly higher then the reference of their class
table(Mastertable$Final_class[Mastertable$Final_viral == "phage"])

# 15.2.1 Malgrandaviricetes (Class level) - ssDNA
#View(Mastertable[Mastertable$Final_class == "Malgrandaviricetes" & Mastertable$length > 10,])
Mastertable$Final_phylum[Mastertable$Final_class == "Malgrandaviricetes" & Mastertable$length > 10 & Mastertable$Cenote_phylum == "Uroviricota"] <- "Uroviricota"
Mastertable$Final_class[Mastertable$Final_class == "Malgrandaviricetes" & Mastertable$length > 10 & Mastertable$Cenote_phylum == "Uroviricota"] <- "Caudoviricetes"
Mastertable$Final_order[Mastertable$Final_class == "Malgrandaviricetes" & Mastertable$length > 10 & Mastertable$Cenote_phylum == "Uroviricota"] <- "Caudovirales"
Mastertable$Final_family[Mastertable$Final_class == "Malgrandaviricetes" & Mastertable$length > 10 & Mastertable$Cenote_phylum == "Uroviricota"] <- "Podoviridae"
# 30 NR-contigs above 10kb
# 6-7kb, some are longer, but confirmed concatamers (also but Cenote-Taker2). Simply state in the table in text.

# 15.2.2 Faserviricetes - Inoviridae (Class level) - ssDNA
#View(Mastertable[Mastertable$Final_class == "Faserviricetes",])
# 7kb, some are longer, but confirmed concatamers

## After adding the Cenotaker2 data, you can cross-check the ones the weren't taxonomicallty annotated by CenoteTaker2
## The ones CenoteTaker2 cannot annotate, Cross-Check them -- alter to 2LCA & Cenote-Taker2.

# 15.2.3 Caudovirales: Crass phage (Species level only available)
# Base this on Cenote-Taker2, if Family is Podovirdae. Then accept it, irrespective of length

# 15.2.4 Duplopiviricetes - Picobirnaviridae -2 segments - ssDNA
#View(Mastertable[Mastertable$Final_class == "Duplopiviricetes",])
# 2kb, OK

# 15.2.5 Allassoviricetes - RNA phage  - ssRNA
#View(Mastertable[Mastertable$Final_class == "Allassoviricetes",])
# 4kb, OK
####################################
# 15.3 Evaluate CrAss phages
####################################
# Annotate all CrAssphages not meeting minimum criteria as "Unannotated": Take a minimum threshold of 10%coverage 
Mastertable$Final_species[Mastertable$Diamond_coverage_relative < 10 & (grepl("crAss",Mastertable$Diamond_species) | grepl("CrAssphage",Mastertable$Diamond_species))] <- "Unannotated"
Mastertable$Final_species[Mastertable$Blastn_coverage_relative < 10 & (grepl("crAss",Mastertable$Blastn_species) | grepl("CrAssphage",Mastertable$Blastn_species))] <- "Unannotated"

# Annotate properly all CrAssphages
Mastertable$Final_species[Mastertable$Diamond_coverage_relative > 10 & (grepl("crAss",Mastertable$Diamond_species) | grepl("CrAssphage",Mastertable$Diamond_species))] <- Mastertable$Diamond_species[Mastertable$Diamond_coverage_relative > 10 & (grepl("crAss",Mastertable$Diamond_species) | grepl("CrAssphage",Mastertable$Diamond_species))]
Mastertable$Final_species[Mastertable$Blastn_coverage_relative > 10 & (grepl("crAss",Mastertable$Blastn_species) | grepl("CrAssphage",Mastertable$Blastn_species))] <- Mastertable$Blastn_species[Mastertable$Blastn_coverage_relative > 10 & (grepl("crAss",Mastertable$Blastn_species) | grepl("CrAssphage",Mastertable$Blastn_species))]
Mastertable$Final_species[grepl("crAss",Mastertable$CAT_species) | grepl("CrAssphage",Mastertable$CAT_species)] <- Mastertable$CAT_species[grepl("crAss",Mastertable$CAT_species) | grepl("CrAssphage",Mastertable$CAT_species)]

# Annotate entire taxonomy of Crassphages
Mastertable$Final_phylum[grepl("crAss",Mastertable$Final_species) | grepl("CrAssphage",Mastertable$Final_species)] <- "Uroviricota"
Mastertable$Final_class[grepl("crAss",Mastertable$Final_species) | grepl("CrAssphage",Mastertable$Final_species)] <- "Caudoviricetes"
Mastertable$Final_order[grepl("crAss",Mastertable$Final_species) | grepl("CrAssphage",Mastertable$Final_species)] <- "Caudovirales"
Mastertable$Final_family[grepl("crAss",Mastertable$Final_species) | grepl("CrAssphage",Mastertable$Final_species)] <- "Podoviridae"

## Evaluate again after excluding crapy alignments
#View(Mastertable[Mastertable$Final_viral == "phage" & Mastertable$Final_family == "Inoviridae",])
#View(Mastertable[Mastertable$Final_viral == "phage" & (grepl("crAss",Mastertable$Final_species) | grepl("CrAssphage",Mastertable$Final_species)),])
sum(Mastertable$Totalnumberofreads[Mastertable$Final_viral == "phage" & (grepl("crAss",Mastertable$Final_species) | grepl("CrAssphage",Mastertable$Final_species))])

# Rename the entire taxonomy of CrAss phages by making a distinction with other Caudovirales
Mastertable$Final_class[Mastertable$Final_viral == "phage" & (grepl("crAss",Mastertable$Final_species) | grepl("CrAssphage",Mastertable$Final_species))] <- "Caudoviricetes (CrAss)"

# Rename Caudoviricetes (non-CrAss)
Mastertable$Final_class[Mastertable$Final_viral == "phage" & Mastertable$Final_class == "Caudoviricetes"] <- "Caudoviricetes (non-CrAss)"
table(Mastertable$Final_class[Mastertable$Final_viral == "phage"])

# Compare LCA Caudoviricetes (CrAss) & Crassphages annotated compared to DB
table(Mastertable$CrAss_taxonomy)
Mastertable$CrAss_taxonomy[is.na(Mastertable$CrAss_taxonomy)] <- "non-CrAss"
#View(Mastertable[Mastertable$Final_class == "Caudoviricetes (CrAss)",])
#View(Mastertable[Mastertable$CrAss_taxonomy == "Crassphage",])

# Annotate all taxonomical level
Mastertable$Final_class[Mastertable$CrAss_taxonomy == "crAssphage"] <- "Caudoviricetes (CrAss)"
Mastertable$Final_order[Mastertable$CrAss_taxonomy == "crAssphage"] <- "Caudovirales"
Mastertable$Final_family[Mastertable$CrAss_taxonomy == "crAssphage"] <- "Podoviridae"
#View(Mastertable[Mastertable$CrAss_taxonomy == "Crassphage",])
sum(Mastertable$Totalnumberofreads[Mastertable$CrAss_taxonomy == "crAssphage"])
nrow(Mastertable[Mastertable$CrAss_taxonomy == "crAssphage",])

# Final_class2 == "Caudoviricetes"
Mastertable$Final_class2[Mastertable$Final_class == "Caudoviricetes (CrAss)"] <- "Caudoviricetes"
Mastertable$Final_class2[Mastertable$Final_class == "Caudoviricetes (non-CrAss)"] <- "Caudoviricetes"
table(Mastertable$Final_class2)

# Calculate the percentage of annotated and unannotated phages (contigs/reads)
# Reads
Mastertable$Totalnumberofreads[is.na(Mastertable$Totalnumberofreads)] <- 0
table(Mastertable$Final_superkingdom)

sum(Mastertable$Totalnumberofreads)
sum(Mastertable$Totalnumberofreads[!Mastertable$Final_superkingdom == "Viruses"])
sum(Mastertable$Totalnumberofreads[Mastertable$Final_superkingdom == "Viruses"])
(sum(Mastertable$Totalnumberofreads[Mastertable$Final_superkingdom == "Viruses"])/sum(Mastertable$Totalnumberofreads))*100

sum(Mastertable$Totalnumberofreads[Mastertable$Final_superkingdom == "Viruses"])
sum(Mastertable$Totalnumberofreads[Mastertable$Final_superkingdom == "Viruses" & !Mastertable$Final_class == "Unannotated"])
sum(Mastertable$Totalnumberofreads[Mastertable$Final_superkingdom == "Viruses" & Mastertable$Final_class == "Unannotated"])
(sum(Mastertable$Totalnumberofreads[Mastertable$Final_superkingdom == "Viruses" & !Mastertable$Final_class == "Unannotated"])/sum(Mastertable$Totalnumberofreads[Mastertable$Final_superkingdom == "Viruses"]))*100

# Contigs
nrow(Mastertable[Mastertable$Final_superkingdom == "Viruses",])
nrow(Mastertable[Mastertable$Final_superkingdom == "Viruses" & !Mastertable$Final_class == "Unannotated",])
nrow(Mastertable[Mastertable$Final_superkingdom == "Viruses" & Mastertable$Final_class == "Unannotated",])
(nrow(Mastertable[Mastertable$Final_superkingdom == "Viruses" & !Mastertable$Final_class == "Unannotated",])/nrow(Mastertable[Mastertable$Final_superkingdom == "Viruses",]))*100

# We were able to identify 3067 near-complete viral genomes to which 64,53% of the PE reads are mapping
# 93% of the NR-contigs are annotated upon Class level to which 97,4% of the viral PE reads are mapping.
####################################
# 14. Curation viral fraction
####################################
## In case of presence of eukaryotic viruses that are clearly wrongly annotated you have to manually change them to what the correct classification should be. For example, giant viruses can't be present (filter in wetlab limits this)
####################################
# 14.1 Curate eukaryotic viruses
####################################
# Add AS
Mastertable$Blastn_AS <- ((Mastertable$Blastn_ANI/100)*(Mastertable$Blastn_coverage_relative/100))
Mastertable$Blastn_AS[!Mastertable$Blastn_viral == "Blastn_eukaryotic_virus"] <- 0

Mastertable$Diamond_AS <- ((Mastertable$Diamond_ANI/100)*(Mastertable$Diamond_coverage_relative/100))
Mastertable$Diamond_AS[!Mastertable$Diamond_viral == "Diamond_eukaryotic_virus"] <- 0

# Annotation tools to use
Mastertable$Best_AS[Mastertable$Diamond_AS > Mastertable$Blastn_AS] <- "Diamond"
Mastertable$Best_AS[Mastertable$Diamond_AS < Mastertable$Blastn_AS] <- "Blastn"
Mastertable$Best_AS[Mastertable$Diamond_AS < 0.1 & Mastertable$Blastn_AS < 0.1] <- "Remove"

# Re-annotate: Blastn
Mastertable$Final_phylum[Mastertable$Best_AS == "Blastn"] <- Mastertable$Blastn_phylum[Mastertable$Best_AS == "Blastn"]
Mastertable$Final_class[Mastertable$Best_AS == "Blastn"] <- Mastertable$Blastn_class[Mastertable$Best_AS == "Blastn"]
Mastertable$Final_order[Mastertable$Best_AS == "Blastn"] <- Mastertable$Blastn_order[Mastertable$Best_AS == "Blastn"]
Mastertable$Final_family[Mastertable$Best_AS == "Blastn"] <- Mastertable$Blastn_family[Mastertable$Best_AS == "Blastn"]
Mastertable$Final_subfamily[Mastertable$Best_AS == "Blastn"] <- Mastertable$Blastn_subfamily[Mastertable$Best_AS == "Blastn"]
Mastertable$Final_genus[Mastertable$Best_AS == "Blastn"] <- Mastertable$Blastn_genus[Mastertable$Best_AS == "Blastn"]
Mastertable$Final_species[Mastertable$Best_AS == "Blastn"] <- Mastertable$Blastn_species[Mastertable$Best_AS == "Blastn"]
Mastertable$Final_ANI[Mastertable$Best_AS == "Blastn"] <- Mastertable$Blastn_ANI[Mastertable$Best_AS == "Blastn"]
Mastertable$Final_coverage[Mastertable$Best_AS == "Blastn"] <- Mastertable$Blastn_coverage_relative[Mastertable$Best_AS == "Blastn"]

# Re-annotate: Diamond
Mastertable$Final_phylum[Mastertable$Best_AS == "Diamond"] <- Mastertable$Diamond_phylum[Mastertable$Best_AS == "Diamond"] 
Mastertable$Final_class[Mastertable$Best_AS == "Diamond"] <- Mastertable$Diamond_class[Mastertable$Best_AS == "Diamond"] 
Mastertable$Final_order[Mastertable$Best_AS == "Diamond"] <- Mastertable$Diamond_order[Mastertable$Best_AS == "Diamond"] 
Mastertable$Final_family[Mastertable$Best_AS == "Diamond"] <- Mastertable$Diamond_family[Mastertable$Best_AS == "Diamond"] 
Mastertable$Final_subfamily[Mastertable$Best_AS == "Diamond"] <- Mastertable$Diamond_subfamily[Mastertable$Best_AS == "Diamond"] 
Mastertable$Final_genus[Mastertable$Best_AS == "Diamond"] <- Mastertable$Diamond_genus[Mastertable$Best_AS == "Diamond"] 
Mastertable$Final_species[Mastertable$Best_AS == "Diamond"] <- Mastertable$Diamond_species[Mastertable$Best_AS == "Diamond"] 
Mastertable$Final_ANI[Mastertable$Best_AS == "Diamond"] <- Mastertable$Diamond_ANI[Mastertable$Best_AS == "Diamond"]
Mastertable$Final_coverage[Mastertable$Best_AS == "Diamond"] <- Mastertable$Diamond_coverage_relative[Mastertable$Best_AS == "Diamond"]

# Remove bad alignments from eukaryotic viral fraction
Mastertable$Final_viral[Mastertable$Best_AS == "Remove"] <- 1
Mastertable$Final_superkingdom[Mastertable$Best_AS == "Remove"] <- "Unannotated"

# Reads and contigs
sum(Mastertable$Totalnumberofreads[Mastertable$Final_viral == "eukaryotic virus"])
nrow(Mastertable[Mastertable$Final_viral == "eukaryotic virus",])

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
#View(Mastertable[Mastertable$Final_viral == "eukaryotic virus", c(vector_3:vector_2, vector_ANI_Blastn:vector_cov_Blastn, vector_ANI_Diamond:vector_cov_Diamond)])
#View(Mastertable)
table(Mastertable$Final_superkingdom)

## Blastx and Blastn taxonomical annotation gives the samples annotation on genus level & slightly different on species level.
## We prioritize Diamond over Blastn & put additional Blastn annotation if they're extra.
## We also put an absolute minimal %coverage of 10% of the contigs (although a lot of them are higher ofcourse), to curate the bad annotation.
#################################### 
# 14.1.3 Export as excel file
####################################
# You can then if you want inspect how closely or distantly related you're eukaryotic viruses are from the databases you used.
# I don'tadd here if it is nucleotide or protein sequence similairty so be carefull. Take into account that the absolute majority is derived from Diamond. I will not do this for phages because they are to different, doesn't add a lot. 

# First make this folder "eukaryotic_viruses" on you mac desktop ofcourse. This path you have to change to your own PATHs.

setwd("/Users/daan/Desktop/Bioinformatics/Analysis/Biologicals/eukaryotic_viruses")
getwd()
dir()

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

######### !!!!!!!!!!!!!!! #########
    # SAVE Global Environment #
######### !!!!!!!!!!!!!!! #########
