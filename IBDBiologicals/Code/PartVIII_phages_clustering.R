####################################
# SCRIPT 10: PROKARYOTIC VIRUSES: UNIQUENESS AND CLUSTERING
####################################
# Before starting this script
####################################
# ------------------> Load the "Global environment" output of of the script of after rarefaction.
####################################
# 0. Packages: Install packages and load them wherever needed
####################################
#install.packages("ggplot2")
library(ggplot2)
#install.packages("dplyr")
library(dplyr)
#install.packages("plyr")
library(plyr)
#install.packages("reshape")
library(reshape)
#install.packages("reshape2")
library(reshape2)
#install.packages("scales")
library(scales)
#install.packages("viridis")
library(viridis)
#install.packages("readr")
library(readr)
#install.packages("WriteXLS")
library(WriteXLS)
#install.packages("readxl")
library(readxl)
#install.packages("remotes")
library(remotes)
#install.packages("devtools") # I have the problem that Mac Big sure blocks access to all kind of system on the mac and had to disable system integratie protection
library(devtools)
#install.packages("BiocManager") # # before downlading from github you need to load devtools
library(BiocManager)
#BiocManager::install("ComplexHeatmap")
# install_github("jokergoo/ComplexHeatmap")
library(ComplexHeatmap)
#source("https://raw.githubusercontent.com/joey711/phyloseq/master/inst/scripts/installer.R",local = TRUE) ## like this it works!!
library(phyloseq)
#BiocManager::install("ALDEx2")
library(ALDEx2)
#install.packages("RSQLite")
library(RSQLite)
#install.packages("https://cran.r-project.org/src/contrib/Archive/RcppArmadillo/RcppArmadillo_0.9.900.3.0.tar.gz", repos=NULL, type="source")
library(RcppArmadillo)
#BiocManager::install("DESeq2")
library(DESeq2)
#remotes::install_github("yiluheihei/microbiomeMarker")
library(microbiomeMarker)
#install.packages("vegan")
library(vegan)
#install.packages("dendsort")
library(dendsort)
install.packages("ape")
library(ape)
#install.packages("seriation")
library(seriation)
#install.packages("taxonomizr")
library(taxonomizr)
#install.packages("esquisse")
library(esquisse)
#install.packages("tidyverse")
library(tidyverse)
#install.packages("grid")
library(grid)
#install.packages("coin")
library(coin)
#install.packages("rstatix")
library(rstatix)
#install.packages("PMCMR")
library(PMCMR)
#remotes::install_github("tapj/biotyper")
library(biotyper)
#BiocManager::install("DirichletMultinomial")
#install.packages("DirichletMultinomial")
library(DirichletMultinomial)
#install.packages("shiny")
library(shiny)
#install.packages("fpc")
library (fpc)
#install.packages("clValid")
library (clValid)
#install.packages("cluster")
library(cluster) 
#install.packages("philentropy")
library(philentropy)
####################################
# 1. Adding new class clusters based on %AAI
###################################
# 1.1 Genus-like groups input
###################################
setwd ("/Users/daan/Desktop/Transfer/Multi_fasta_files") 
getwd()
dir()

AllSamples_classification_genus<- as.data.frame(read_delim("genus_clusters.txt", "\t", escape_double = FALSE, trim_ws = TRUE, col_names = FALSE))
colnames(AllSamples_classification_genus) <- c('contig', 'Genus')
rownames(AllSamples_classification_genus) <- AllSamples_classification_genus$contig
AllSamples_classification_genus$contig <- NULL
rownames(AllSamples_classification_genus) <- gsub("\\.{1}", "_", rownames(AllSamples_classification_genus))
###################################
# 1.2 Merge with mastertable
###################################
Mastertable_viral_rarefied_1 <- merge(Mastertable_viral_rarefied, AllSamples_classification_genus,by=0, all=F)
rownames(Mastertable_viral_rarefied_1) <- Mastertable_viral_rarefied_1$Row.names
Mastertable_viral_rarefied_1$Row.names <- NULL
Mastertable_viral_rarefied_1[1:377] # good
###################################
# 1.3 Subset phages again to be sure
###################################
Mastertable_viral_rarefied_2 <- Mastertable_viral_rarefied_1[Mastertable_viral_rarefied_1$Final_viral == "phage",]
unique(Mastertable_viral_rarefied_2$Genus) # 735 genus-like phages
###################################
# 1.4. Phyloseq for clusters
####################################
#View(taxa(Mastertable_viral_rarefied_2)) 
####################################
# 1.4.1 Create phyloseq abundance table
####################################
Mastertable_viral_rarefied_2$Final_ANI <- NULL
Mastertable_viral_rarefied_2$Final_coverage <- NULL
Mastertable_viral_rarefied_2$Final_class2 <- NULL
Mastertable_viral_rarefied_2$Totalnumberofreads <- NULL
Mastertable_viral_rarefied_2$Blastn_AS <- NULL
Mastertable_viral_rarefied_2$Diamond_AS <- NULL
Mastertable_viral_rarefied_2$Best_AS <- NULL
Mastertable_viral_rarefied_2$Final_subfamily <- NULL
Mastertable_viral_rarefied_2$Final_nodes <- NULL
Mastertable_viral_rarefied_2$names <- rownames(Mastertable_viral_rarefied_2)

colnames(Mastertable_viral_rarefied_2)[which(names(Mastertable_viral_rarefied_2) == "Final_superkingdom")] <- "Final_kingdom"
vector_1 <- (which(names(Mastertable_viral_rarefied_2)== "Virsorter")-1) 

abundance_table_rarefied <- Mastertable_viral_rarefied_2[,c(1:vector_1)]
# abundance_table_rarefied <- Mastertable_viral_unrarefied_phages[,c(1:vector_1)]
sort(rowSums(abundance_table_rarefied))
####################################
# 1.4.2 Create phyloseq taxonomy table
####################################
Mastertable_viral_rarefied_2$Final_genus <- NULL
Mastertable_viral_rarefied_2 <- rename(Mastertable_viral_rarefied_2, c("Genus" = "Final_genus"))
vector_2 <- (which(names(Mastertable_viral_rarefied_2)== "Final_kingdom"))
vector_3 <- which(names(Mastertable_viral_rarefied_2)== "names")

taxonomy_table_rarefied <- Mastertable_viral_rarefied_2[,c(vector_2:vector_3)]
names(taxonomy_table_rarefied)[names(taxonomy_table_rarefied) == "Final_kingdom"] <- "Kingdom"
names(taxonomy_table_rarefied)[names(taxonomy_table_rarefied) == "Final_phylum"] <- "Phylum"
names(taxonomy_table_rarefied)[names(taxonomy_table_rarefied) == "Final_class"] <- "Class"
names(taxonomy_table_rarefied)[names(taxonomy_table_rarefied) == "Final_order"] <- "Order"
names(taxonomy_table_rarefied)[names(taxonomy_table_rarefied) == "Final_family"] <- "Family"
#names(taxonomy_table_rarefied)[names(taxonomy_table_rarefied) == "Final_subfamily"] <- "Subfamily"
names(taxonomy_table_rarefied)[names(taxonomy_table_rarefied) == "Final_genus"] <- "Genus"
taxonomy_table_rarefied$Final_species <- NULL
taxonomy_table_rarefied$Final_species <- taxonomy_table_rarefied$names
taxonomy_table_rarefied$names <- NULL
names(taxonomy_table_rarefied)[names(taxonomy_table_rarefied) == "Final_species"] <- "Species"
#View(taxonomy_table_rarefied)
####################################
# 1.4.3 Create phyloseq sample tables (contain samples + metadata)
####################################
setwd("/Users/daan/Desktop/Bioinformatics/pre-processing/Biological_study/Metadata/output_metadata")
getwd()
dir()

metadata <- as.data.frame(read_excel("Final_Metadata_subselection_Relevant_R_without-units.xlsx"))
rownames(metadata) <- metadata$`Randomized Nr.`
names(metadata)[1]<-paste("sample")

sample_table_rarefied_for_all <- (merge(t(abundance_table_rarefied),metadata, by=0, all=F)) # this is done so the metadata matches the samples
rownames(sample_table_rarefied_for_all) <- sample_table_rarefied_for_all$Row.names
sample_table_rarefied_for_all$Row.names <- NULL

vector_4 <- (which(names(sample_table_rarefied_for_all)== "sample"))
vector_5 <- which(names(sample_table_rarefied_for_all)== "Combined_remission")
sample_table_rarefied_for_all <- sample_table_rarefied_for_all[,c(vector_4:vector_5)]
sample_table_rarefied_for_all$Endoscopic_outcome_combined[sample_table_rarefied_for_all$Endoscopic_outcome_combined == "NA"] <- "unknown"
sample_table_rarefied_for_all$Combined_remission[sample_table_rarefied_for_all$Combined_remission == "NA"] <- "unknown"

# Add new columns 'timepoints'
sample_table_rarefied_for_all$timepoints <- 1
sample_table_rarefied_for_all$timepoints[sample_table_rarefied_for_all$Week == "0"] <- "baseline"
sample_table_rarefied_for_all$timepoints[sample_table_rarefied_for_all$Diagnosis == "CD" & sample_table_rarefied_for_all$Week == 14] <- "no"
sample_table_rarefied_for_all$timepoints[sample_table_rarefied_for_all$timepoints == 1] <- "primary endpoint"
table(sample_table_rarefied_for_all$timepoints)

## Remove Data (unknowns & w14)
sample_table_rarefied_for_all <- sample_table_rarefied_for_all[!c(sample_table_rarefied_for_all$Week == "14" & sample_table_rarefied_for_all$Diagnosis == "CD"),]
#sample_table_rarefied_for_all <- sample_table_rarefied_for_all[!sample_table_rarefied_for_all$Combined_remission == "unknown",]
#sample_table_rarefied_for_all <- sample_table_rarefied_for_all[sample_table_rarefied_for_all$Week == "0",]
#sample_table_rarefied_for_all <- sample_table_rarefied_for_all[sample_table_rarefied_for_all$timepoints == "baseline" & sample_table_rarefied_for_all$Diagnosis == "UC",]
#sample_table_rarefied_for_all <- sample_table_rarefied_for_all[!sample_table_rarefied_for_all$timepoints == "baseline",]
#View(sample_table_rarefied_for_all)

sample_table_rarefied_for_all$Endoscopic_outcome_combined[sample_table_rarefied_for_all$Endoscopic_outcome_combined == "1"] <- "remission"
sample_table_rarefied_for_all$Endoscopic_outcome_combined[sample_table_rarefied_for_all$Endoscopic_outcome_combined == "0"] <- "non-remission"
#sample_table_rarefied_for_all$Endoscopic_outcome_combined[sample_table_rarefied_for_all$timepoints == "baseline"] <- "non-remission"

sample_table_rarefied_for_all$Clinical_outcome_combined[sample_table_rarefied_for_all$Clinical_outcome_combined == "1"] <- "remission"
sample_table_rarefied_for_all$Clinical_outcome_combined[sample_table_rarefied_for_all$Clinical_outcome_combined == "0"] <- "non-remission"
#sample_table_rarefied_for_all$Clinical_outcome_combined[sample_table_rarefied_for_all$timepoints == "baseline"] <- "non-remission"

sample_table_rarefied_for_all$Biomarker_outcome_combined[sample_table_rarefied_for_all$Biomarker_outcome_combined == "1"] <- "remission"
sample_table_rarefied_for_all$Biomarker_outcome_combined[sample_table_rarefied_for_all$Biomarker_outcome_combined == "0"] <- "non-remission"
#sample_table_rarefied_for_all$Biomarker_outcome_combined[sample_table_rarefied_for_all$timepoints == "baseline"] <- "non-remission"

sample_table_rarefied_for_all$Combined_remission[sample_table_rarefied_for_all$Combined_remission == "1"] <- "remission"
sample_table_rarefied_for_all$Combined_remission[sample_table_rarefied_for_all$Combined_remission == "0"] <- "non-remission"
#sample_table_rarefied_for_all$Combined_remission[sample_table_rarefied_for_all$timepoints == "baseline"] <- "non-remission"

# convert Age to numeric & others
str(sample_table_rarefied_for_all)
#sample_table_rarefied_for_all$sample_ID <- NULL
#sample_table_rarefied_for_all$patient_ID <- NULL
sample_table_rarefied_for_all$Sample_date <- NULL
sample_table_rarefied_for_all$Date_diagnosis <- NULL
sample_table_rarefied_for_all$Date_of_birth <- NULL
sample_table_rarefied_for_all$`Sample dilution color` <- NULL
#sample_table_rarefied_for_all$Therapy_1 <- NULL
sample_table_rarefied_for_all$endoscopic_mayo <- NULL
sample_table_rarefied_for_all$Week <- NULL
sample_table_rarefied_for_all$Smoking_2 <- NULL
#sample_table_rarefied_for_all$Smoking <- NULL 
sample_table_rarefied_for_all$total_mayo <- NULL
sample_table_rarefied_for_all$SES_CD <- NULL
sample_table_rarefied_for_all$PRO2 <- NULL
sample_table_rarefied_for_all$sample <- NULL
sample_table_rarefied_for_all$active_disease_baseline <- NULL
#sample_table_rarefied_for_all$Clinical_outcome_combined <- NULL
#sample_table_rarefied_for_all$Biomarker_outcome_combined <- NULL
sample_table_rarefied_for_all$Age <- as.numeric(sample_table_rarefied_for_all$Age)
names(sample_table_rarefied_for_all)[names(sample_table_rarefied_for_all) == 'time_sample_to_sequencing'] <- 'Storage'
names(sample_table_rarefied_for_all)[names(sample_table_rarefied_for_all) == 'Therapy_2'] <- 'Therapy'
names(sample_table_rarefied_for_all)[names(sample_table_rarefied_for_all) == 'Disease_Duration'] <- 'Duration'
names(sample_table_rarefied_for_all)[names(sample_table_rarefied_for_all) == 'Disease_location_start'] <- 'Location'
names(sample_table_rarefied_for_all)[names(sample_table_rarefied_for_all) == 'bacterial_Cellcount'] <- 'Cellcount'
names(sample_table_rarefied_for_all)[names(sample_table_rarefied_for_all) == 'Endoscopic_outcome_combined'] <- 'Endoscopy-outcome'
names(sample_table_rarefied_for_all)[names(sample_table_rarefied_for_all) == 'biologicals_before'] <- 'Biologicals'
names(sample_table_rarefied_for_all)[names(sample_table_rarefied_for_all) == 'biological_classes_before'] <- 'Biologicals-classes'
names(sample_table_rarefied_for_all)[names(sample_table_rarefied_for_all) == 'biologicals_before'] <- 'Biologicals'
names(sample_table_rarefied_for_all)[names(sample_table_rarefied_for_all) == 'Clinical_outcome_combined'] <- 'Clinical-outcome'
names(sample_table_rarefied_for_all)[names(sample_table_rarefied_for_all) == 'Biomarker_outcome_combined'] <- 'Biomarker-outcome'
sample_table_rarefied_for_all$Gender[sample_table_rarefied_for_all$Gender == "0"] <- "female"
sample_table_rarefied_for_all$Gender[sample_table_rarefied_for_all$Gender == "1"] <- "male"
#sample_table_rarefied_for_all$Smoking[sample_table_rarefied_for_all$Smoking == "0"] <- "non-smoker"
#sample_table_rarefied_for_all$Smoking[sample_table_rarefied_for_all$Smoking == "1"] <- "smoker"
sample_table_rarefied_for_all$Smoking_2[sample_table_rarefied_for_all$Smoking_2 == "0"] <- "non-smoker"
sample_table_rarefied_for_all$Smoking_2[sample_table_rarefied_for_all$Smoking_2 == "1"] <- "smoker"
sample_table_rarefied_for_all$Smoking_2[sample_table_rarefied_for_all$Smoking_2 == "2"] <- "ex-smoker"
names(sample_table_rarefied_for_all)[names(sample_table_rarefied_for_all) == 'Smoking_2'] <- 'Smoking'
sample_table_rarefied_for_all$`Clinical-outcome`[sample_table_rarefied_for_all$`Clinical-outcome` == "0"] <- "non-remission"
sample_table_rarefied_for_all$`Clinical-outcome`[sample_table_rarefied_for_all$`Clinical-outcome` == "1"] <- "remission"
sample_table_rarefied_for_all$`Biomarker-outcome`[sample_table_rarefied_for_all$`Biomarker-outcome` == "0"] <- "non-remission"
sample_table_rarefied_for_all$`Biomarker-outcome`[sample_table_rarefied_for_all$`Biomarker-outcome` == "1"] <- "remission"
sample_table_rarefied_for_all$Gender <- as.character(sample_table_rarefied_for_all$Gender)
sample_table_rarefied_for_all$Age <- as.numeric(sample_table_rarefied_for_all$Age)
sample_table_rarefied_for_all$Week <- as.numeric(sample_table_rarefied_for_all$Week)
sample_table_rarefied_for_all$Moisture <- as.numeric(sample_table_rarefied_for_all$Moisture, replace = TRUE)
sample_table_rarefied_for_all$Cellcount <- as.numeric(sample_table_rarefied_for_all$Cellcount, replace = TRUE)
sample_table_rarefied_for_all$Calprotectin <- as.numeric(sample_table_rarefied_for_all$Calprotectin, replace = TRUE)
sample_table_rarefied_for_all$BMI <- as.numeric(sample_table_rarefied_for_all$BMI, replace = TRUE)
sample_table_rarefied_for_all$Hemoglobin <- as.numeric(sample_table_rarefied_for_all$Hemoglobin, replace = TRUE)
sample_table_rarefied_for_all$Albumin <- as.numeric(sample_table_rarefied_for_all$Albumin, replace = TRUE)
sample_table_rarefied_for_all$CRP <- as.numeric(sample_table_rarefied_for_all$CRP, replace = TRUE)
sample_table_rarefied_for_all$Hemoglobin <- as.numeric(sample_table_rarefied_for_all$Hemoglobin, replace = TRUE)
sample_table_rarefied_for_all$patient_ID <- as.character(sample_table_rarefied_for_all$patient_ID)
sample_table_rarefied_for_all[is.na(sample_table_rarefied_for_all)] <- NA
sample_table_rarefied_for_all$Therapy_1 <- NULL
sample_table_rarefied_for_all$Therapy[sample_table_rarefied_for_all$timepoints == 'baseline'] <- "untreated"
#sample_table_rarefied_for_all$Location[sample_table_rarefied_for_all$Diagnosis == "UC"] <- "L2"
str(sample_table_rarefied_for_all)
#View(sample_data(sample_table_rarefied_for_all))
####################################
# 1.4.4 Overview all unrarefied tables (abundance, taxonomy and sample tables)
####################################
## ABUNDANCE TABLE
#View(abundance_table_rarefied)

## TAXONOMY TABLE 
#View(taxonomy_table_rarefied)

## SAMPLE TABLE
#View(sample_table_rarefied_for_all)
####################################
# 1.5 Create phyloSeq objects
####################################
# 1.5.1 Create a matrix of abundance and taxonomy tables
####################################
abundance_table_rarefied_m <- as.matrix(abundance_table_rarefied)
taxonomy_table_rarefied_m <- as.matrix(taxonomy_table_rarefied)
#View(abundance_table_rarefied_m)
####################################
# 1.5.2 Transform to phylseq objects
####################################
ABUNDANCE_rarefied_family <- otu_table(abundance_table_rarefied_m, taxa_are_rows = TRUE)
TAX_rarefied_family <- tax_table(taxonomy_table_rarefied_m)
samples <- sample_data(sample_table_rarefied_for_all)
phyloseq_rarefied_family <- phyloseq(ABUNDANCE_rarefied_family,TAX_rarefied_family, samples)
####################################
# 1.5.3 Visualize data & subset only the phage fraction
####################################
sample_names(phyloseq_rarefied_family) # All sample names
rank_names(phyloseq_rarefied_family) # All taxonomies
sample_variables(phyloseq_rarefied_family) ## All metadata
# if you want to subset the dataframe into only "responders", just do "responder_data <- subset_samples(phyloseq_unrarefied_family, response == "responder")"
phyloseq_rarefied_phages <- subset_taxa(phyloseq_rarefied_family, Kingdom == "Viruses")
####################################
# 1.6 PhyloSeq format
####################################
phyloseq_rarefied_phages
####################################
# 2. Community typing
###################################
# SUMMARY:
# using DMM clustering and BIC stability indices the number is cluster are determined on 2.
# All the phages are subsequently being visualized on a PcOA plot.
###################################
# 2.1 Data trimming (aggregate common taxa)
###################################
phyloseq_rarefied_phages_visualization <- aggregate_taxa(phyloseq_rarefied_phages, 'Genus')
phyloseq_rarefied_phages_DMM <- phyloseq_rarefied_phages_visualization %>%
  microbiome::transform(transform = "identity") %>%
 aggregate_rare(level = "Genus", detection = 0.01/100, prevalence = 20/100)
###################################
# 2.2 DMM algorithm
###################################
# Pick the OTU count matrix and convert it into samples x taxa format
dat <- abundances(phyloseq_rarefied_phages_DMM)
count <- as.data.frame(t(dat))
phyloseq_rarefied_phages_DMM
## Remove 'Other' meaning all the unique viruses appearing only in less than 20% of samples (0.01 % abundance)
count$total <- rowSums(count)
count <- count[!count$total == 0,]
count$total <- NULL
count$Other <- NULL

## Remove samples composed of zero shared viruses
count$total <- rowSums(count)
count <- count[!count$total == 0,]
sort(rowSums(count))
count$total <- NULL

## How many samples and reads are left for determining DMM clusters?
nrow(count) # 363 samples left for which community-typing is possible
sum(count) # 363 samples having 120M reads
(sum(count)/sum(Mastertable_viral_rarefied$Totalnumberofreads))*100 # or 50% of the reads shared

## Median genera per sample using community typing
count_1 <- count
count_1[count_1 > 0] <- 1
count_1 <- as.data.frame(count_1)
count_1$prevalence <- rowSums(count_1)
ncol(count_1) # 18 genera shared over 20% of samples with 0.01% abundance
median(count_1$prevalence) # Median of 6 samples per samples
mean(count_1$prevalence)
nrow(count_1) # 363 samples can be used for clusters
(nrow(count_1)/377)*100 # or > 96% can be used.

# Determine DMM clusters and evaluate BIC score
count <- as.matrix(count)
fit <- lapply(1:10, dmn, count = count, verbose=TRUE) 
# lplc <- sapply(fit, laplace) # AIC / BIC / Laplace
# aic  <- sapply(fit, AIC) # AIC / BIC / Laplace
bic  <- sapply(fit, BIC) # AIC / BIC / Laplace
plot(bic, type="b", xlab="Number of Dirichlet Components", ylab="Model Fit")
#lines(lplc, type="b", lty = 2)
#lines(aic, type="b", lty = 3)

# Determine the optimal model
best <- fit[[which.min(unlist(bic))]]
mixturewt(best)
best@group
best@mixture

# Additionally, evaluate the probability of samples belonging to one of the two DMM groups. If not that interesting you could leave this out
group_probability <- as.data.frame(best@group)
mean(group_probability$V1)
mean(group_probability$V2)
median(group_probability$V1)
median(group_probability$V2)

# Subselect the clusters for all samples as determined before
cluster <- apply(mixture(best), 1, which.max)
cluster <- as.data.frame(cluster)
cluster$cluster <- as.character(cluster$cluster)
table(cluster$cluster) # cluster 1 244 samples, cluster 2 119 samples.
cluster_phylo <- sample_data(cluster)

# At last find drivers or contribution of each taxonomic group to each component
for (k in seq(ncol(fitted(best)))) {
  d <- melt(fitted(best))
  colnames(d) <- c("OTU", "cluster", "value")
  d <- subset(d, cluster == k) %>%
    # Arrange OTUs by assignment strength
    arrange(value) %>%
    mutate(OTU = factor(OTU, levels = unique(OTU))) %>%
    # Only show the most important drivers
    filter(abs(value) > quantile(abs(value), 0.99)) #0.8 default     
  
  p <- ggplot(d, aes(x = OTU, y = value)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(title = paste("Top drivers: community type", k))
  print(p)
}

# Fit genera to model clusters
# A fitted value is a statistical model’s prediction of the mean response value when you input the values of the predictors, factor levels, or components into the model. Suppose you have the following regression equation: y = 3X + 5. If you enter a value of 5 for the predictor, the fitted value is 20. Fitted values are also called predicted values.
# fitted is a generic function which extracts fitted values from objects returned by modeling functions.

# Cluster 1 
# Faecalibacterium phahes = Genus_5 # check: 
# Bacteroidees phages = Genus_69
# Unknown caudo (non-crass) = Genus_2

# Cluster 2
# Lactococcus = Genus_35 (-)

# cluster 3
# Crass = Genus 3, Genus_14, Genus_103 (+)
# Micro = Genus_1, Genus_7, Genus_8, Genus_6

contribution_viruses <- as.data.frame(fitted(best))
contribution_viruses$cluster1 <- contribution_viruses$V1
contribution_viruses$cluster2 <- contribution_viruses$V2
contribution_viruses$V1 <- NULL
contribution_viruses$V2 <- NULL
###################################
# 2.3 Add clusters to phyloseq object
###################################
phyloseq_rarefied_phages_visualization_1 <- merge_phyloseq(phyloseq_rarefied_phages_visualization, cluster_phylo)
###################################
# 2.4 Determine number of genera in samples
###################################
# How was a genus defined again?
# ----> determined as in https://github.com/snayfach/MGV/tree/master/aai_cluster
# ----> genera are binned and defined on similarity of more or equal than 20% amino acid identity (%AAI) between genomes.
# ----> and genomes with either 8 shared genes or at least 20% of shared genes (relative to both genomes)
# ----> So basically based on generic %AAI cutoff for phages

genera <- as.data.frame(t(abundances(phyloseq_rarefied_phages_visualization_1)))
genera[genera > 0] <- 1
genera$prevalence <- rowSums(genera)

# generic infiormation
# A) Min and maximum genera in samples
range(genera$prevalence) # 3 min; 74 maximum

# B) genera in total
nrow(tax_table(phyloseq_rarefied_phages_visualization_1)) # 874 genera

# C) Median genera per sample
median(genera$prevalence) # 26 median
mean(genera$prevalence) # 27.23 average
###################################
# 2.5 Bray-curtis dissimilarity and visualization
###################################
# 2.5.1 Remove samples for which no community-type can be assigned
###################################
# Remove NA's for clusters 
phyloseq_rarefied_phages_visualization_2 = subset_samples(phyloseq_rarefied_phages_visualization_1, !cluster == "NA")
# (n = 363 of the n = 377 have been assigned; or 96,3% of samples have a community-typing possible)
###################################
# 2.5.2 Transform count (log10)
###################################
# I use a log10 transformation for the simple fact that I observed before that taxa that are less abundant or often lysogenic phages
# These are more shared and less abundant than lytic taxa. I don't want to underepresent these so log10 transformation
# You could also opt for other transformation (eg. Hellinger transformation is famous one givin more weight to low-abundant taxa)
phyloseq_rarefied_phages_transformed = transform_sample_counts(phyloseq_rarefied_phages_visualization_2, function(x) log(10 + x))
###################################
# 2.5.3 Bray-curtis dissimmilarities
###################################
phyloseq_rarefied_phages_transformed_bray <- phyloseq::distance(phyloseq_rarefied_phages_transformed, method="bray")
###################################
# 2.5.4 Reduce dimensionality and plot virome composition on PcoA
###################################
phyloseq_rarefied_phages_PcoA <- ordinate(phyloseq_rarefied_phages_transformed, method="PCoA", distance=phyloseq_rarefied_phages_transformed_bray)
#View(sample_data(phyloseq_rarefied_phages_transformed))
###################################
# 2.6 Visualization of microbiome using BC dissimilarities and PcOA plot
###################################
# Phyloseq: phyloseq_rarefied_phages_transformed
# Bray-curtis data: phyloseq_rarefied_phages_transformed_bray
# Ordinationation: phyloseq_rarefied_phages_PcoA

metadata_begin <- as.data.frame(sample_data(phyloseq_rarefied_phages_transformed))
metadata_begin$timepoints <- as.character(metadata_begin$timepoints)
metadata_begin$timepoints

metadata_begin$timepoints[metadata_begin$timepoints == "baseline"] <- "pre-intervention + CD w14"
metadata_begin$timepoints[metadata_begin$timepoints == "no"] <- "pre-intervention + CD w14"
metadata_begin$timepoints[metadata_begin$timepoints == "primary endpoint"] <- "post-intervention"

metadata_begin$timepoints <- factor(metadata_begin$timepoints, levels = c("post-intervention","pre-intervention + CD w14"))
metadata_begin_ <- sample_data(metadata_begin)
phyloseq_rarefied_phages_transformed <- merge_phyloseq(metadata_begin_, phyloseq_rarefied_phages_transformed)

# PcoA 1: without elipses
plot_ordination(phyloseq_rarefied_phages_transformed, phyloseq_rarefied_phages_PcoA, color="cluster", shape = "timepoints") +
  geom_point(size = 1, alpha=1)  +
  scale_shape_manual(values = c(16,17)) +
  theme_bw() +
  scale_color_manual(values = c(`1` = "#bd7f1c", `2` = "#2a8a92")) +
  geom_hline(yintercept = 0, linetype = 3, colour = "black") + 
  geom_vline(xintercept = 0, linetype = 3, colour = "black") +
  xlab("PC1 (8.6%)") +
  ylab("PC2 (4.1%)")

  # PcoA 2: with ellipes
pcoa_plot_ <- plot_ordination(phyloseq_rarefied_phages_transformed, phyloseq_rarefied_phages_PcoA, color="cluster", shape = "Endoscopy.outcome",title="DMM clustering") +
  geom_point(size=2, alpha=0.75)  
  theme_classic() +
  scale_shape_manual(values = c(21,16, 2)) + ## depends on Therapy_2 number of variables, 2 in this case
  scale_color_manual(values = c("#bd7f1c","#2a8a92")) +
  geom_hline(yintercept = 0, linetype = 3, colour = "black") + 
  geom_vline(xintercept = 0, linetype = 3, colour = "black")  +
  xlab("PC1 (8.6%)") +
  ylab("PC2 (4.1%)") +
  theme(axis.text.y = element_text(colour = "black", size = 9, face = "bold"),
        axis.text.x = element_text(colour = "black", face = "bold", size = 9),
        legend.text = element_text(size = 9, face = "bold", colour = "black"),
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14),
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"),
        legend.title = element_blank(), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2)) +
  stat_ellipse(type='t',size =1, linetype = 2,level = 0.95,aes(group = cluster))

plot(pcoa_plot_)

# PcoA 3: Diagnosis
pcoa_plot_R <- plot_ordination(phyloseq_rarefied_phages_transformed, phyloseq_rarefied_phages_PcoA, color="Endoscopy.outcome", shape = "Diagnosis",title="DMM clustering") +
  geom_point(size=2, alpha=0.75)  +
  theme_classic() 

plot(pcoa_plot_R)

# PcoA 4: Diagnosis
pcoa_plot_IBD <- plot_ordination(phyloseq_rarefied_phages_transformed, phyloseq_rarefied_phages_PcoA, color="Diagnosis", shape = "Endoscopy.outcome",title="DMM clustering") +
  geom_point(size=2, alpha=0.75)  +
  theme_classic() +

plot(pcoa_plot_IBD)
# Some samples do not exceed the threshold of 0.01 abundance and 20% detection.
# These will not have a cluster ventoerotype assigned. We will also not show them in this PcoA plot.

phyloseq_rarefied_phages_visualization
###################################
# 2.7 Bacterial community types - viral community types 
###################################
# Insert enterotypes data
setwd("/Users/daan/Desktop/Bioinformatics/Analysis/Biologicals/Rstudio_analysis/6. Prokaryotic_viruses/6.6 uniqueness_clustering")
enterotypes <- as.data.frame(read_excel("Enterotypes_DJ.xlsx"))
rownames(enterotypes) <- enterotypes$nummer
enterotypes1 <- enterotypes[,c(3,4:7)]

# 363 samples are rarefied and have a viral community types
phyloseq_rarefied_phages_cluster <- as.data.frame(sample_data(phyloseq_rarefied_phages_transformed))
viral_community_types <- phyloseq_rarefied_phages_cluster[,c(1,4,24,25)]

rownames(viral_community_types) <- viral_community_types$sample_ID
viral_community_types$sample_ID <- NULL 
viral_community_types$cluster[viral_community_types$cluster == "1"] <- "CA"
viral_community_types$cluster[viral_community_types$cluster == "2"] <- "CrM"

# Combine data
venterotypes <- (merge(enterotypes1,viral_community_types, by=0, all=F))
venterotypes$Row.names <- NULL
venterotypes$QMP_g_Bacteroides <- as.numeric(venterotypes$QMP_g_Bacteroides)
venterotypes$QMP_g_Prevotella <- as.numeric(venterotypes$QMP_g_Prevotella)
venterotypes$RMP_g_Bacteroides <- as.numeric(venterotypes$RMP_g_Bacteroides)
venterotypes$RMP_g_Prevotella <- as.numeric(venterotypes$RMP_g_Prevotella)
str(venterotypes)
enterotypes_host <- venterotypes

table(venterotypes$cluster[venterotypes$Enterotype == "Bact1"]) # 66 samples
(40/66)*100
# 60,6%

table(venterotypes$cluster[venterotypes$Enterotype == "Bact2"]) # 178 samples
(33/178)*100
# 18.5%

table(venterotypes$cluster[venterotypes$Enterotype == "Prev"]) # 29 samples
(9/29)*100
# 31.0%

table(venterotypes$cluster[venterotypes$Enterotype == "Rum"]) # 22 samples
(17/22)*100
# 77.3%

# Correlation matrix
install.packages("corrplot")
library(corrplot)
library(Hmisc)
venterotypes2 <- venterotypes
str(venterotypes2)

ggplot(venterotypes2) +
  aes(x = Enterotype, fill = cluster) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = c(`CA` = "#BD801C", `CrM` = "#2B8A92")) +
  theme_bw()

ggplot(venterotypes2) +
 aes(x = cluster, fill = Enterotype) +
 geom_bar(position = "dodge") +
 scale_fill_hue(direction = 1) +
 theme_bw() +
  
ggplot(venterotypes2) +
  aes(x = Enterotype, fill = cluster) +
  geom_bar(position = "dodge") +
  scale_fill_manual(values = c(`CA` = "#BD801C", `CrM` = "#2B8A92")) +
  theme_bw()

ggplot(venterotypes2) +
 aes(x = cluster, fill = Enterotype) +
 geom_bar(position = "fill") +
 scale_fill_hue(direction = 1) +
 theme_classic()

# Proportion tables: IBD cohort
venterotypes2 <- venterotypes[,c(1,7)]

View(venterotypes2)
venterotypes2_table <- table(venterotypes2)
venterotypes2_table_prop <- prop.table(venterotypes2_table, margin = 2)
venterotypes2_table_prop <- as.data.frame(venterotypes2_table_prop)
venterotypes2_table_prop$perc <- venterotypes2_table_prop$Freq*100

ggplot(venterotypes2_table_prop) +
 aes(x = Enterotype, fill = cluster, weight = perc) +
 geom_bar(position = "dodge") +
 ylab("Percentage (%)") +
 scale_fill_manual(values = c(`CA` = "#BD801C", `CrM` = "#2B8A92")) +
 theme_bw()

# Statistics enterotypes
chisq.test(venterotypes2_table, correct = F)
cramer_v(venterotypes2_table)
sum(venterotypes2_table)

# pairwise_prop_test(venterotypes2_table, p.adjust.method = "BH", correct = FALSE)
row_wise_prop_test(venterotypes2_table, p.adjust.method = "BH",correct = TRUE)
row_wise_prop_test(venterotypes2_table, p.adjust.method = "BH",correct = TRUE, detailed = TRUE)
p.adjust(0.923, method = "BH", n = 4)
?row_wise_prop_test

# Proportion tables: PE
venterotypes3 <- venterotypes[,c(1,6,7)]
venterotypes3 <- venterotypes3[venterotypes3$timepoints == "post-intervention",]
venterotypes3$timepoints <- NULL
venterotypes3 <- table(venterotypes3)
venterotypes3 <- prop.table(venterotypes3, margin = 2)
venterotypes3 <- as.data.frame(venterotypes3)
venterotypes3$perc <- venterotypes3$Freq*100

ggplot(venterotypes3) +
  aes(x = Enterotype, fill = cluster, weight = perc) +
  geom_bar(position = "dodge") +
  ylab("Percentage (%)") +
  scale_fill_manual(values = c(`CA` = "#BD801C", `CrM` = "#2B8A92")) +
  theme_bw()

# within IBD: UC versus CD
venterotypesUC <- venterotypes[venterotypes$Diagnosis == "UC",c(1,8)]
venterotypesUC_table <- table(venterotypesUC)
venterotypesUC_table_prop <- prop.table(venterotypesUC_table, margin = 2)
venterotypesUC_table_prop <- as.data.frame(venterotypesUC_table_prop)
venterotypesUC_table_prop$perc <- venterotypesUC_table_prop$Freq*100
venterotypesUC_table_prop$cluster <- as.character(venterotypesUC_table_prop$cluster)
venterotypesUC_table_prop$Enterotype <- as.character(venterotypesUC_table_prop$Enterotype)

# UC
ggplot(venterotypesUC_table_prop) +
  aes(x = Enterotype, fill = cluster, weight = perc) +
  geom_bar(position = "dodge") +
  ylab("Percentage (%)") +
  scale_fill_manual(values = c(`CA` = "#BD801C", `CrM` = "#2B8A92")) +
  theme_bw()

ggplot(venterotypesUC_table_prop) +
 aes(x = cluster, fill = Enterotype, weight = perc) +
 geom_bar(position = "dodge") +
 scale_fill_manual(values = c(Bact1 = "#BD9931", Bact2 = "#9B3E1E", Prev = "#9BAD6C", Rum = "#5C8EB2")) +
 theme_bw()

chisq.test(venterotypesUC_table)
prop.test(venterotypesUC_table)
row_wise_prop_test(venterotypesUC_table, p.adjust.method = "BH",correct = TRUE)
row_wise_prop_test(venterotypesUC_table, p.adjust.method = "BH",correct = TRUE, detailed = TRUE)

xtab <- as.table(rbind(
  c(8,19), c(28,10),
  c(9,4), c(4,6)))

dimnames(xtab) <- list(
  Class = c("Bact1", "Bact2", "Prev", "Rum"),
  Clusters = c("cluster CA", "cluster CrM"))

crosstable_statistics(xtab,statistics = "phi")
# Phi(r) = 0.4042

# CD
venterotypesCD <- venterotypes[venterotypes$Diagnosis == "CD",c(1,8)]
venterotypesCD_table <- table(venterotypesCD)
venterotypesCD_table_prop <- prop.table(venterotypesCD_table, margin = 2)
venterotypesCD_table_prop <- as.data.frame(venterotypesCD_table_prop)
venterotypesCD_table_prop$perc <- venterotypesCD_table_prop$Freq*100

ggplot(venterotypesCD_table_prop) +
  aes(x = Enterotype, fill = cluster, weight = perc) +
  geom_bar(position = "dodge") +
  ylab("Percentage (%)") +
  scale_fill_manual(values = c(`CA` = "#BD801C", `CrM` = "#2B8A92")) +
  theme_bw()

ggplot(venterotypesCD_table_prop) +
  aes(x = cluster, fill = Enterotype, weight = perc) +
  geom_bar(position = "dodge") +
  scale_fill_manual(values = c(Bact1 = "#BD9931", Bact2 = "#9B3E1E", Prev = "#9BAD6C", Rum = "#5C8EB2")) +
  theme_bw()

prop.test(venterotypesCD_table)
row_wise_prop_test(venterotypesCD_table, p.adjust.method = "BH",correct = TRUE) # compare individual comparisons
row_wise_prop_test(venterotypesCD_table, p.adjust.method = "BH",correct = TRUE, detailed = TRUE)

venterotypesCD_table

xtab <- as.table(rbind(
  c(18,21), c(117,23),
  c(11,5), c(1,11)))

dimnames(xtab) <- list(
  Class = c("Bact1", "Bact2", "Prev", "Rum"),
  Clusters = c("cluster CA", "cluster CrM"))

crosstable_statistics(xtab,statistics = "phi")
# Phi(r) = 0.4042

# No clear difference for Prevotella and Bacteroides abundances
###################################
# 3. Barplots of metadata and cluster assignment: quick visualization
###################################
# 3.1 Determine cluster assignment (probability)
###################################
group_probability
names(group_probability)[names(group_probability) == "V1"] <- "cluster 1"
names(group_probability)[names(group_probability) == "V2"] <- "cluster 2"
group_probability <- merge(group_probability,cluster, by = 0, all = F)
rownames(group_probability) <- group_probability$Row.names
names(group_probability)[names(group_probability) == "Row.names"] <- "sample ID"
names(group_probability)[names(group_probability) == "cluster"] <- "selected_cluster"
group_probability$probability_selection <- "0"
group_probability$probability_selection[group_probability$selected_cluster == "1"] <- group_probability$`cluster 1`[group_probability$selected_cluster == "1"]
group_probability$probability_selection[group_probability$selected_cluster == "2"] <- group_probability$`cluster 2`[group_probability$selected_cluster == "2"]
group_probability$probability_selection <- as.numeric(group_probability$probability_selection)
mean(group_probability$probability_selection)
median(group_probability$probability_selection)
# View(group_probability)

#library("writexl")
setwd("/Users/daan/Desktop/Bioinformatics/Analysis/Biologicals/Phages/virome_composition/DMM/Determining_clusters")
write_xlsx(group_probability,"./group_probability.xlsx")
###################################
# 3.2 Associate alpha-diversity with DMM clusters
###################################
# You want to compare Richness outcome (continous) between 2 cateogrical groups (cluster 1/cluster 2)
# For comparing 2 groups (1 with continous variable: "Richness" and 1 with categorical: "Cluster") you can use wilcoxon test.
# On genus level, or if you prefer also on viral species level
phyloseq_rarefied_phages_visualization_2_richness <- phyloseq_rarefied_phages
#phyloseq_rarefied_phages_visualization_2_richness <- phyloseq_rarefied_phages_visualization_2

# richness
richness <- richness(phyloseq_rarefied_phages_visualization_2_richness, index = c("observed"))
diversity <- diversity(phyloseq_rarefied_phages_visualization_2_richness, index = "shannon", zeroes = TRUE)
alpha_diversity <- merge(richness,diversity, by = 0, all =T)
rownames(alpha_diversity) <- alpha_diversity$Row.names
alpha_diversity$Row.names <- NULL

# merge with clusters
alpha_diversity_cluster <- merge(alpha_diversity,cluster, by = 0, all =T)
rownames(alpha_diversity_cluster) <- alpha_diversity_cluster$Row.names
alpha_diversity_cluster$Row.names <- NULL

# Add to phyloseq object
alpha_diversity_cluster <- sample_data(alpha_diversity_cluster)
phyloseq_rarefied_phages_3 <- merge_phyloseq(phyloseq_rarefied_phages_visualization_2_richness, alpha_diversity_cluster)

# Remove NA's: samples without cluster (n = 351 of the n = 377 have been assigned)
phyloseq_rarefied_phages_3 = subset_samples(phyloseq_rarefied_phages_3, !cluster == "NA")

# Boxplot
## Richness & shannon diversity
phyloseq_rarefied_phages_3_stat <- as.data.frame(sample_data(phyloseq_rarefied_phages_3))
names(phyloseq_rarefied_phages_3_stat)[names(phyloseq_rarefied_phages_3_stat) == "observed"] <- "Observed richness"
names(phyloseq_rarefied_phages_3_stat)[names(phyloseq_rarefied_phages_3_stat) == "shannon"] <- "Shannon diversity"

## Make figure Diversity/richness
phyloseq_rarefied_phages_3_stat <- phyloseq_rarefied_phages_3_stat[,c(1,3,20,22,25)]
phyloseq_rarefied_phages_3_stat$patient_ID <- NULL
phyloseq_rarefied_phages_3_stat <- melt(phyloseq_rarefied_phages_3_stat, id.vars = c("cluster"))

View(phyloseq_rarefied_phages_3_stat[phyloseq_rarefied_phages_3_stat$Endoscopy.outcome == "remission",])
ggplot(phyloseq_rarefied_phages_3_stat) +
aes(x = cluster, y = value , fill = cluster) +
geom_boxplot(shape = "circle") +
geom_jitter(width = 0.15) +
scale_fill_manual(values = c(`1` = "#990011FF", `2` = "#00203FFF")) +
theme_bw() +
ylab("") +
facet_wrap(vars(variable), scales = "free_y")

## Statistics:
# To test median differences of continuous variables between two or more groups of continuous variables (alpha diversity measures, abundances, etc.), Mann-Whitney U test and Kruskal-Wallis (KW) test were performed respectively. The KW test was always followed by post hoc Dunn’s (phD) test for all pairs of comparisons between groups. Multiple testing correction was performed where appropriate using the Benjamini-Hochberg procedure (FDR adjustment set at <0.05).
# Richness statistics
unique(phyloseq_rarefied_phages_3_stat$variable)
phyloseq_rarefied_phages_3_stat_richness <- phyloseq_rarefied_phages_3_stat[phyloseq_rarefied_phages_3_stat$variable == "Observed.richness",]
wilcox.test(phyloseq_rarefied_phages_3_stat_richness$value~phyloseq_rarefied_phages_3_stat_richness$cluster, alternative = "two.sided",  paired = FALSE, conf.int = T, conf.level = 0.95, exact = FALSE, correct = FALSE)
p.adjust(2.379e-13, method = "BH", n = 2) #  P.adjsut = 4.758e-13
nrow(phyloseq_rarefied_phages_3_stat_richness[phyloseq_rarefied_phages_3_stat_richness$cluster == 1,]) # 242 sample
nrow(phyloseq_rarefied_phages_3_stat_richness[phyloseq_rarefied_phages_3_stat_richness$cluster == 2,]) # 121 samples
effect_size <- wilcox_effsize(phyloseq_rarefied_phages_3_stat_richness,value~cluster)
?wilcox_effsize
(effect_size$effsize)
(effect_size$effsize)^2
effect_size$n1
effect_size$n2
median(phyloseq_rarefied_phages_3_stat_richness$value[phyloseq_rarefied_phages_3_stat_richness$cluster == "1"])
median(phyloseq_rarefied_phages_3_stat_richness$value[phyloseq_rarefied_phages_3_stat_richness$cluster == "2"])
# viral contig level: W = 4914, p-value < 2.2e-16, P.adjust < 2.2e-16, n1 = 242; n2 = 121 samples, r = 0.542 (effect size), median cluster 1 = 27, median cluster 2 = 46
# viral genus level: W = 7945, p-value = 1.205e-12, n1 = 242; n2 = 121 samples, r = 0.373 (effect size)

## Diversity statistics
View(phyloseq_rarefied_phages_3_stat_shannon)
phyloseq_rarefied_phages_3_stat_shannon <- phyloseq_rarefied_phages_3_stat[phyloseq_rarefied_phages_3_stat$variable == "Shannon.diversity",]
wilcox.test(phyloseq_rarefied_phages_3_stat_shannon$value~phyloseq_rarefied_phages_3_stat_shannon$cluster, alternative = "two.sided",  paired = FALSE, conf.int = T, conf.level = 0.95, exact = FALSE, correct = FALSE)
p.adjust(2.379e-13, method = "BH", n = 2) # P.adjust = 4.758e-13
nrow(phyloseq_rarefied_phages_3_stat_shannon[phyloseq_rarefied_phages_3_stat_shannon$cluster == 1,]) # 242 samples
nrow(phyloseq_rarefied_phages_3_stat_shannon[phyloseq_rarefied_phages_3_stat_shannon$cluster == 2,]) # 121 samples
effect_size <- wilcox_effsize(phyloseq_rarefied_phages_3_stat_shannon,value~cluster)
(effect_size$effsize)
(effect_size$effsize)^2
effect_size$n1
effect_size$n2
# viral contig level: W = 7737, p-value = 2.379e-13, P.adjust = 4.758e-13, n1 = 242; n2 = 121 samples, r = 0.384 (effect size)
# viral genus level: W = 7945, p-value = 1.205e-12, n1 = 242; n2 = 121 samples, r = 0.373 (effect size)
###################################
# 3.3 Associate Endoscopic outcome with DMM clusters (covariate: see later)
###################################
# You want to compare clusters outcome (categorical) between 3 categorical groups (R/NR/unknown)
# For comparing 3 categorical groups (Endoscopical outcome: R, NR and unknown) with 2 categorical groups (clusters: 1 or 2) or more you can use kruskall wallis test with post-hoc dun's test and fdr test for multiple comparisons
# For comparing 2 cateogrical groyps (Endoscopical outcome: R and NR) with with categorical groups (clusters: 1 or 2) you can use chisquared test.
# Confirm later with Jurgen Vercauteren

## Extract metadata
colnames(sample_data(phyloseq_rarefied_phages_transformed))
phyloseq_rarefied_phages_transformed_meta <-as.data.frame(sample_data(phyloseq_rarefied_phages_transformed))
phyloseq_rarefied_phages_transformed_R <- phyloseq_rarefied_phages_transformed_meta[,c(4,21,25)]
phyloseq_rarefied_phages_transformed_R_UC <- phyloseq_rarefied_phages_transformed_R[phyloseq_rarefied_phages_transformed_R$Diagnosis == "UC",] 
phyloseq_rarefied_phages_transformed_R_CD <- phyloseq_rarefied_phages_transformed_R[phyloseq_rarefied_phages_transformed_R$Diagnosis == "CD",] 
#phyloseq_rarefied_phages_transformed_R_L1 <- phyloseq_rarefied_phages_transformed_R[phyloseq_rarefied_phages_transformed_R$Location == "L1",] 
#phyloseq_rarefied_phages_transformed_R_L2 <- phyloseq_rarefied_phages_transformed_R[phyloseq_rarefied_phages_transformed_R$Location == "L2",] 
#phyloseq_rarefied_phages_transformed_R_l3 <- phyloseq_rarefied_phages_transformed_R[phyloseq_rarefied_phages_transformed_R$Location == "L3",] 

## There are two questions:
# 1) How do the enterotypes look like? 
# 2) If you look at a specific variable do you see more enterotype prevalence within the value of that variables?

# 1) Create table
phyloseq_rarefied_phages_transformed_R_all <- phyloseq_rarefied_phages_transformed_R
phyloseq_rarefied_phages_transformed_R_allUC <- phyloseq_rarefied_phages_transformed_R_all[phyloseq_rarefied_phages_transformed_R_all$Diagnosis == "UC"]
phyloseq_rarefied_phages_transformed_R_allCD <- phyloseq_rarefied_phages_transformed_R_all[phyloseq_rarefied_phages_transformed_R_all$Diagnosis == "CD"]

## Remove NA's
#phyloseq_rarefied_phages_transformed_R <- phyloseq_rarefied_phages_transformed_R[!phyloseq_rarefied_phages_transformed_R$Location == "NA"]
table(phyloseq_rarefied_phages_transformed_R$Endoscopy.outcome)

## For ALL
phyloseq_rarefied_phages_transformed_R$value <- 1
phyloseq_rarefied_phages_transformed_R$freq <- 1
phyloseq_rarefied_phages_transformed_R$freq_end <- 1
phyloseq_rarefied_phages_transformed_R$freq[phyloseq_rarefied_phages_transformed_R$cluster == 1] <- phyloseq_rarefied_phages_transformed_R$value[phyloseq_rarefied_phages_transformed_R$cluster == 1]/nrow(phyloseq_rarefied_phages_transformed_R[phyloseq_rarefied_phages_transformed_R$cluster == 1])
phyloseq_rarefied_phages_transformed_R$freq[phyloseq_rarefied_phages_transformed_R$cluster == 2] <- phyloseq_rarefied_phages_transformed_R$value[phyloseq_rarefied_phages_transformed_R$cluster == 2]/nrow(phyloseq_rarefied_phages_transformed_R[phyloseq_rarefied_phages_transformed_R$cluster == 2])
phyloseq_rarefied_phages_transformed_R$freq <- phyloseq_rarefied_phages_transformed_R$freq*100
table(phyloseq_rarefied_phages_transformed_R$Endoscopy.outcome[phyloseq_rarefied_phages_transformed_R$cluster == "1"])
table(phyloseq_rarefied_phages_transformed_R$Endoscopy.outcome[phyloseq_rarefied_phages_transformed_R$cluster == "2"])

## For CD
phyloseq_rarefied_phages_transformed_R_allCD$value <- 1
phyloseq_rarefied_phages_transformed_R_allCD$freq <- 1
phyloseq_rarefied_phages_transformed_R_allCD$freq_end <- 1
phyloseq_rarefied_phages_transformed_R_allCD$freq[phyloseq_rarefied_phages_transformed_R_allCD$cluster == 1] <- phyloseq_rarefied_phages_transformed_R_allCD$value[phyloseq_rarefied_phages_transformed_R_allCD$cluster == 1]/nrow(phyloseq_rarefied_phages_transformed_R_allCD[phyloseq_rarefied_phages_transformed_R_allCD$cluster == 1])
phyloseq_rarefied_phages_transformed_R_allCD$freq[phyloseq_rarefied_phages_transformed_R_allCD$cluster == 2] <- phyloseq_rarefied_phages_transformed_R_allCD$value[phyloseq_rarefied_phages_transformed_R_allCD$cluster == 2]/nrow(phyloseq_rarefied_phages_transformed_R_allCD[phyloseq_rarefied_phages_transformed_R_allCD$cluster == 2])
phyloseq_rarefied_phages_transformed_R_allCD$freq <- phyloseq_rarefied_phages_transformed_R_allCD$freq*100
table(phyloseq_rarefied_phages_transformed_R_allCD$Endoscopy.outcome[phyloseq_rarefied_phages_transformed_R_allCD$cluster == "1"])
table(phyloseq_rarefied_phages_transformed_R_allCD$Endoscopy.outcome[phyloseq_rarefied_phages_transformed_R_allCD$cluster == "2"])

## For UC
phyloseq_rarefied_phages_transformed_R_allUC$value <- 1
phyloseq_rarefied_phages_transformed_R_allUC$freq <- 1
phyloseq_rarefied_phages_transformed_R_allUC$freq_end <- 1
phyloseq_rarefied_phages_transformed_R_allUC$freq[phyloseq_rarefied_phages_transformed_R_allUC$cluster == 1] <- phyloseq_rarefied_phages_transformed_R_allUC$value[phyloseq_rarefied_phages_transformed_R_allUC$cluster == 1]/nrow(phyloseq_rarefied_phages_transformed_R_allUC[phyloseq_rarefied_phages_transformed_R_allUC$cluster == 1])
phyloseq_rarefied_phages_transformed_R_allUC$freq[phyloseq_rarefied_phages_transformed_R_allUC$cluster == 2] <- phyloseq_rarefied_phages_transformed_R_allUC$value[phyloseq_rarefied_phages_transformed_R_allUC$cluster == 2]/nrow(phyloseq_rarefied_phages_transformed_R_allUC[phyloseq_rarefied_phages_transformed_R_allUC$cluster == 2])
phyloseq_rarefied_phages_transformed_R_allUC$freq <- phyloseq_rarefied_phages_transformed_R_allUC$freq*100
table(phyloseq_rarefied_phages_transformed_R_allUC$Endoscopy.outcome[phyloseq_rarefied_phages_transformed_R_allUC$cluster == "1"])
table(phyloseq_rarefied_phages_transformed_R_allUC$Endoscopy.outcome[phyloseq_rarefied_phages_transformed_R_allUC$cluster == "2"])

# IBD
ggplot(phyloseq_rarefied_phages_transformed_R) +
  aes(x = cluster, fill = Endoscopy.outcome, weight = freq) +
  geom_bar(position = "dodge") +
  # scale_fill_manual(values = list(`L1` = "#ad4b38", `L2` = "#499b78", `L3` = "#354b9a",)) +
  theme_classic() +
  ggtitle("IBD") +
  ylab("Endoscopic outcome (%)") +
  theme(axis.title.y = element_text(size= 16))

# CD
ggplot(phyloseq_rarefied_phages_transformed_R_allCD) +
  aes(x = cluster, fill = Endoscopy.outcome, weight = freq) +
  geom_bar(position = "dodge") +
  # scale_fill_manual(values = list(`L1` = "#ad4b38", `L2` = "#499b78", `L3` = "#354b9a",)) +
  theme_classic() +
  ggtitle("CD") +
  ylab("Endoscopic outcome (%)") +
  theme(axis.title.y = element_text(size= 16))

# UC
ggplot(phyloseq_rarefied_phages_transformed_R_allUC) +
  aes(x = cluster, fill = Endoscopy.outcome, weight = freq) +
  geom_bar(position = "dodge") +
  # scale_fill_manual(values = list(`L1` = "#ad4b38", `L2` = "#499b78", `L3` = "#354b9a",)) +
  theme_classic() +
  ggtitle("IBD") +
  ylab("Endoscopic outcome (%)") +
  theme(axis.title.y = element_text(size= 16))

# 2) Statistics
# All
phyloseq_rarefied_phages_transformed_R_all_chi <- chisq.test(table(phyloseq_rarefied_phages_transformed_R_all[,2:3])) # X-squared: 8.3016, p-value 0.01575
phyloseq_rarefied_phages_transformed_R_all_post_hoc <- chisq.posthoc.test(table(phyloseq_rarefied_phages_transformed_R_all[,2:3]), method = "BH")
phyloseq_rarefied_phages_transformed_R_all_chi
phyloseq_rarefied_phages_transformed_R_all_post_hoc #remission compared to remission in

# IBD
table(phyloseq_rarefied_phages_transformed_R_all[,2:3])
IBD_fit <- table(phyloseq_rarefied_phages_transformed_R_all[,2:3])

xtab <- as.table(rbind(
  c(14), c(42, 30),
  c(19, 9)))

dimnames(xtab) <- list(
  Class = c("non-remission", "remission", "unknown"),
  Clusters = c("cluster 1", "cluster 2")
)

# View(xtab)
row_wise_prop_test(xtab,detailed = FALSE, p.adjust.method = "BH")
row_wise_prop_test(xtab,detailed = TRUE, p.adjust.method = "BH")
pairwise_prop_test(xtab, p.adjust.method = "BH")
crosstable_statistics(xtab,statistics = "phi")
# X-squared: 8.30, p-value 0.0158, n = 363

# Non-remission: significantly lower non-remission rate in cluster 2 compared to cluster 1
# Z-test: 7.54, p-value 0.00604, P.adjust = 0.0181, n = 363
114/(114+104+24) # 47.1% non-remission cluster 1  
38/(38+69+14)    # 31.4% non-remission cluster 2

# Remission significantly higher remission rate in cluster 2 compared to cluster 1
# Z-test: 5.83, p-value 0.0157, P.adjust = 0.0314, n = 363
104/(114+104+24) # 43.0% remission cluster 1  
69/(38+69+14)    # 57.0% remission cluster 2

# unknown: no differences
# Z-test: 0.0919, p-value 0.762, P.adjust = 0.762, n = 363

# UC
phyloseq_rarefied_phages_transformed_R_all_chi <- chisq.test(table(phyloseq_rarefied_phages_transformed_R_allUC[,2:3])) # X-squared: 8.3016, p-value 0.01575
phyloseq_rarefied_phages_transformed_R_all_post_hoc <- chisq.posthoc.test(table(phyloseq_rarefied_phages_transformed_R_allUC[,2:3]), method = "BH")
phyloseq_rarefied_phages_transformed_R_all_chi
phyloseq_rarefied_phages_transformed_R_all_post_hoc #remission compared to remission in
prop.table(table(phyloseq_rarefied_phages_transformed_R_allUC[,2:3]), margin =2)
UC_fit <- table(phyloseq_rarefied_phages_transformed_R_allUC[,2:3])
CramerV(UC_fit) # r value 
phyloseq_rarefied_phages_transformed_R_allUC
# X-squared = 6.6956, df = 1,p-value = 0.009665, r = 0.287

# UC
xtab <- as.table(rbind(
  c(31, 11), c(24,29)))

dimnames(xtab) <- list(
  Class = c("non-remission", "remission"),
  Clusters = c("cluster 1", "cluster 2")
)

row_wise_prop_test(xtab, detailed = TRUE, alternative = "two.sided",p.adjust.method = "BH")
#pairwise.prop.test(xtab,p.adjust.method = "BH")
CramerV(xtab) # r value 

?prop_test()

# Two proportions
# Non-remission
prop.test(x = c(31, 11), n = c(55, 40))

# Remission
prop.test(x = c(24, 29), n = c(55, 40))

# --> Same results as before

# Non-remission: significantly lower non-remission rate in cluster 2 compared to cluster 1
# Z-test: 6.70, N = 95, p-value 0.00966, p-adjust = 0.0193
31/(31+24)    # 56.4% non-remission cluster 1  
11/(11+29)    # 27.5% non-remission cluster 2
# CramersFit: 0.286

# Remission significantly higher remission rate in cluster 2 compared to cluster 1
# Z-test: 6.70, N = 95, p-value 0.00966, p-adjust = 0.0193
24/(31+24)    # 43.6% remission cluster 1  
29/(11+29)    # 72.5% remission cluster 2
# CramersFit: 0.286

# CD
phyloseq_rarefied_phages_transformed_R_all_chi <- chisq.test(table(phyloseq_rarefied_phages_transformed_R_allCD[,2:3])) # X-squared: 8.3016, p-value 0.01575
phyloseq_rarefied_phages_transformed_R_all_post_hoc <- chisq.posthoc.test(table(phyloseq_rarefied_phages_transformed_R_allCD[,2:3]), method = "BH")
phyloseq_rarefied_phages_transformed_R_all_chi
phyloseq_rarefied_phages_transformed_R_all_post_hoc #remission compared to remission in
#View(table(phyloseq_rarefied_phages_transformed_R_all[,2:3]))
#View(phyloseq_rarefied_phages_transformed_R_allCD)

# CD
table(phyloseq_rarefied_phages_transformed_R_allCD[,2:3])
83/(80+24+80) # 45.1% NR cluster 1
27/(40+14+27) # 33.3% NR cluster 2

80/(80+24+80) # 43.5% R cluster 1
40/(40+14+27) # 49/4% R cluster 2

24/(80+24+80) # 13.0% U cluster 1
14/(40+14+27) # 17.3% U cluster 2

xtab <- as.table(rbind(
  c(37, 9), c(28,13),
  c(19, 9)))

dimnames(xtab) <- list(
  Class = c("non-remission", "remission", "unknown"),
  Clusters = c("cluster 1", "cluster 2")
)

row_wise_prop_test(xtab, detailed = TRUE, p.adjust.method = "BH")
CramerV(xtab) # r value 

# X-squared: 3.0213, p-value 0.2208, n = 268
# No P adjust not significant for CD

# https://rpkgs.datanovia.com/rstatix/reference/prop_test.html
###################################
# 3.4 Associate Location with DMM clusters (covariate: see later)
###################################
library("chisq.posthoc.test")

## Extract metadata
colnames(sample_data(phyloseq_rarefied_phages_transformed))
phyloseq_rarefied_phages_transformed_meta <- as.data.frame(sample_data(phyloseq_rarefied_phages_transformed))
phyloseq_rarefied_phages_transformed_R <- phyloseq_rarefied_phages_transformed_meta[,c(3,12,24)]

## Remove NA's
#phyloseq_rarefied_phages_transformed_R <- phyloseq_rarefied_phages_transformed_R[!phyloseq_rarefied_phages_transformed_R$Location == "NA"]
table(phyloseq_rarefied_phages_transformed_R$Location)

## For ALL
phyloseq_rarefied_phages_transformed_R$value <- 1
phyloseq_rarefied_phages_transformed_R$freq <- 1
phyloseq_rarefied_phages_transformed_R$freq_end <- 1
phyloseq_rarefied_phages_transformed_R$freq[phyloseq_rarefied_phages_transformed_R$cluster == 1] <- phyloseq_rarefied_phages_transformed_R$value[phyloseq_rarefied_phages_transformed_R$cluster == 1]/nrow(phyloseq_rarefied_phages_transformed_R[phyloseq_rarefied_phages_transformed_R$cluster == 1])
phyloseq_rarefied_phages_transformed_R$freq[phyloseq_rarefied_phages_transformed_R$cluster == 2] <- phyloseq_rarefied_phages_transformed_R$value[phyloseq_rarefied_phages_transformed_R$cluster == 2]/nrow(phyloseq_rarefied_phages_transformed_R[phyloseq_rarefied_phages_transformed_R$cluster == 2])
phyloseq_rarefied_phages_transformed_R$freq <- phyloseq_rarefied_phages_transformed_R$freq*100
table(phyloseq_rarefied_phages_transformed_R$Location[phyloseq_rarefied_phages_transformed_R$cluster == "1"])
table(phyloseq_rarefied_phages_transformed_R$Location[phyloseq_rarefied_phages_transformed_R$cluster == "2"])

ggplot(phyloseq_rarefied_phages_transformed_R) +
  aes(x = cluster, fill = Location, weight = freq) +
  geom_bar(position = "dodge") +
 # scale_fill_manual(values = list(`L1` = "#ad4b38", `L2` = "#499b78", `L3` = "#354b9a",)) +
  theme_classic() +
  ggtitle("IBD") +
  ylab("Location (%)") +
  theme(axis.title.y = element_text(size= 16))

## Statistics
# Statistical differences in the proportions of categorical variables (enterotypes) between patient groups were evaluated using pairwise Chi-square tests.

# 1) Create table
phyloseq_rarefied_phages_transformed_R_all <- phyloseq_rarefied_phages_transformed_R
table(phyloseq_rarefied_phages_transformed_R_all$Diagnosis[phyloseq_rarefied_phages_transformed_R_all$cluster == "1" & phyloseq_rarefied_phages_transformed_R_all$Location == "L2"])
table(phyloseq_rarefied_phages_transformed_R_all$Diagnosis[phyloseq_rarefied_phages_transformed_R_all$cluster == "2" & phyloseq_rarefied_phages_transformed_R_all$Location == "L2"])

# 2) Statistics
# All
#View(table(phyloseq_rarefied_phages_transformed_R_all[,2:3]))
phyloseq_rarefied_phages_transformed_R_all_chi <- chisq.test(table(phyloseq_rarefied_phages_transformed_R_all[,2:3])) # X-squared: 8.3016, p-value 0.01575
phyloseq_rarefied_phages_transformed_R_all_post_hoc <- chisq.posthoc.test(table(phyloseq_rarefied_phages_transformed_R_all[,2:3]), method = "BH")
phyloseq_rarefied_phages_transformed_R_all_chi
phyloseq_rarefied_phages_transformed_R_all_post_hoc #remission compared to remission in

# You don't need to do P;adjust since you will do a pairwise proportion test later & here you do one test not 3.
# X-squared: 1.022, p-value = 0.6, p.adjust = 1, df = 2,n = 273 
table(phyloseq_rarefied_phages_transformed_R_all[,2:3])

xtab <- as.table(rbind(
  c(67, 34), c(50, 18),
  c(73, 31), c(52, 38)))

dimnames(xtab) <- list(
  Class = c("L1", "L2", "L3", "NA"),
  Clusters = c("cluster 1", "cluster 2")
)


row_wise_prop_test(xtab)
CramerV(xtab)
#View(table(phyloseq_rarefied_phages_transformed_R_all))
###################################
# 3.5 Associate Therapy with DMM clusters (covariate: see later)
###################################
library("chisq.posthoc.test")

## Extract metadata
colnames(sample_data(phyloseq_rarefied_phages_transformed))
phyloseq_rarefied_phages_transformed_meta <- as.data.frame(sample_data(phyloseq_rarefied_phages_transformed))
phyloseq_rarefied_phages_transformed_R <- phyloseq_rarefied_phages_transformed_meta[,c(3,6,24)]

## Remove NA's
#phyloseq_rarefied_phages_transformed_R <- phyloseq_rarefied_phages_transformed_R[!phyloseq_rarefied_phages_transformed_R$Location == "NA"]
table(phyloseq_rarefied_phages_transformed_R$Therapy)

## For ALL
phyloseq_rarefied_phages_transformed_R$value <- 1
phyloseq_rarefied_phages_transformed_R$freq <- 1
phyloseq_rarefied_phages_transformed_R$freq_end <- 1
phyloseq_rarefied_phages_transformed_R$freq[phyloseq_rarefied_phages_transformed_R$cluster == 1] <- phyloseq_rarefied_phages_transformed_R$value[phyloseq_rarefied_phages_transformed_R$cluster == 1]/nrow(phyloseq_rarefied_phages_transformed_R[phyloseq_rarefied_phages_transformed_R$cluster == 1])
phyloseq_rarefied_phages_transformed_R$freq[phyloseq_rarefied_phages_transformed_R$cluster == 2] <- phyloseq_rarefied_phages_transformed_R$value[phyloseq_rarefied_phages_transformed_R$cluster == 2]/nrow(phyloseq_rarefied_phages_transformed_R[phyloseq_rarefied_phages_transformed_R$cluster == 2])
phyloseq_rarefied_phages_transformed_R$freq <- phyloseq_rarefied_phages_transformed_R$freq*100
table(phyloseq_rarefied_phages_transformed_R$Therapy[phyloseq_rarefied_phages_transformed_R$cluster == "1"])
table(phyloseq_rarefied_phages_transformed_R$Therapy[phyloseq_rarefied_phages_transformed_R$cluster == "2"])
ggplot(phyloseq_rarefied_phages_transformed_R) +
  aes(x = cluster, fill = Therapy, weight = freq) +
  geom_bar(position = "dodge") +
 # scale_fill_manual(values = list(`L1` = "#ad4b38", `L2` = "#499b78", `L3` = "#354b9a")) +
  theme_classic() +
  ggtitle("IBD") +
  ylab("Location (%)") +
  theme(axis.title.y = element_text(size= 16)) +
  facet_wrap(vars(Diagnosis))


## Statistics
# Statistical differences in the proportions of categorical variables (enterotypes) between patient groups were evaluated using pairwise Chi-square tests.

# 1) Create table
phyloseq_rarefied_phages_transformed_R_all <- phyloseq_rarefied_phages_transformed_R
table(phyloseq_rarefied_phages_transformed_R_all$Diagnosis[phyloseq_rarefied_phages_transformed_R_all$cluster == "1" & phyloseq_rarefied_phages_transformed_R_all$Location == "L2"])
table(phyloseq_rarefied_phages_transformed_R_all$Diagnosis[phyloseq_rarefied_phages_transformed_R_all$cluster == "2" & phyloseq_rarefied_phages_transformed_R_all$Location == "L2"])

# 2) Statistics
# All
#View(table(phyloseq_rarefied_phages_transformed_R_all[,2:3]))
phyloseq_rarefied_phages_transformed_R_all_chi <- chisq.test(table(phyloseq_rarefied_phages_transformed_R_all[,2:3])) # X-squared: 8.3016, p-value 0.01575
phyloseq_rarefied_phages_transformed_R_all_post_hoc <- chisq.posthoc.test(table(phyloseq_rarefied_phages_transformed_R_all[,2:3]), method = "BH")
phyloseq_rarefied_phages_transformed_R_all_chi
phyloseq_rarefied_phages_transformed_R_all_post_hoc #remission compared to remission in

# You don't need to do P;adjust since you will do a pairwise proportion test later & here you do one test not 3.
# X-squared: 1.022, p-value = 0.6, p.adjust = 1, df = 2,n = 273 
table(phyloseq_rarefied_phages_transformed_R_all[,2:3])

xtab <- as.table(rbind(
  c(67, 34), c(50, 18),
  c(73, 31), c(52, 38)))

dimnames(xtab) <- list(
  Class = c("L1", "L2", "L3", "NA"),
  Clusters = c("cluster 1", "cluster 2")
)


row_wise_prop_test(xtab)
CramerV(xtab)
#View(table(phyloseq_rarefied_phages_transformed_R_all))
###################################
# 3.5 Associate Moisture with DMM cluster (covariate: see later)
###################################
## Extract metadata
colnames(sample_data(phyloseq_rarefied_phages_transformed))
phyloseq_rarefied_phages_transformed_meta <- as.data.frame(sample_data(phyloseq_rarefied_phages_transformed))
phyloseq_rarefied_phages_transformed_R <- phyloseq_rarefied_phages_transformed_meta[,c(3,7,24)]
phyloseq_rarefied_phages_transformed_R$Moisture[is.na(phyloseq_rarefied_phages_transformed_R$Moisture)] <- "remove"
phyloseq_rarefied_phages_transformed_R <- phyloseq_rarefied_phages_transformed_R[!phyloseq_rarefied_phages_transformed_R$Moisture == "remove",]
phyloseq_rarefied_phages_transformed_R$Moisture <- as.numeric(phyloseq_rarefied_phages_transformed_R$Moisture)

# Figure
ggplot(phyloseq_rarefied_phages_transformed_R) +
  aes(x = cluster, y = Moisture, fill = cluster) +
  geom_boxplot(shape = "circle") +
  scale_fill_hue(direction = 1) +
  theme_bw() 

ggplot(phyloseq_rarefied_phages_transformed_R) +
 aes(x = cluster, y = Moisture, fill = cluster) +
 geom_boxplot(shape = "circle") +
 scale_fill_hue(direction = 1) +
 theme_bw() +
 facet_wrap(vars(Diagnosis))

# wilcoxon
str(phyloseq_rarefied_phages_transformed_R)
wilcox.test(phyloseq_rarefied_phages_transformed_R$Moisture ~ phyloseq_rarefied_phages_transformed_R$cluster, alternative = "two.sided",  paired = FALSE, exact = F)
phyloseq_rarefied_phages_transformed_R$

# UC
phyloseq_rarefied_phages_transformed_R_UC_ <- phyloseq_rarefied_phages_transformed_R[phyloseq_rarefied_phages_transformed_R$Diagnosis == "UC",] 
wilcox.test(phyloseq_rarefied_phages_transformed_R_UC_$Moisture ~ phyloseq_rarefied_phages_transformed_R_UC_$cluster, alternative = "two.sided",  paired = FALSE, exact = F)

# CD
phyloseq_rarefied_phages_transformed_R_CD_ <- phyloseq_rarefied_phages_transformed_R[phyloseq_rarefied_phages_transformed_R$Diagnosis == "CD",] 
wilcox.test(phyloseq_rarefied_phages_transformed_R_CD_$Moisture ~ phyloseq_rarefied_phages_transformed_R_CD_$cluster, alternative = "two.sided",  paired = FALSE, exact = F)

# Not significantly different
###################################
# 3.6 RA Lysogenic/lytic phages within DMM cluster
###################################
Mastertable_viral_rarefied_phages_DMM <- Mastertable_viral_rarefied_phages
vector_1_first <- (which(names(Mastertable_viral_rarefied_phages_DMM)== "Virsorter"))
vector_1_1_first <- (which(names(Mastertable_viral_rarefied_phages_DMM)== "Virsorter")-1)
vector_2_first <- (which(names(Mastertable_viral_rarefied_phages_DMM)== "Totalnumberofreads"))
vector_ab_viral_NMDS_col_first <- Mastertable_viral_rarefied_phages_DMM[,1:vector_1_1_first]
vector_ab_viral_NMDS_col_tax_first <- Mastertable_viral_rarefied_phages_DMM[,vector_1_first:vector_2_first]
names <- colnames(vector_ab_viral_NMDS_col_first)

# Lysogenic and Lytic
Mastertable_viral_rarefied_phages_DMM$lysogenic.cycle[Mastertable_viral_rarefied_phages_DMM$lysogenic.cycle == "temperate"] <- "lysogenic"
Mastertable_viral_rarefied_phages_DMM$lysogenic.cycle[Mastertable_viral_rarefied_phages_DMM$lysogenic.cycle == "0"] <- "lytic"
table(Mastertable_viral_rarefied_phages_DMM$lysogenic.cycle)

# Aggregate based on lytic/lysogenic
Mastertable_viral_rarefied_phages_DMM_lifestyle <- aggregate(. ~lysogenic.cycle, FUN = sum, data = Mastertable_viral_rarefied_phages_DMM[,colnames(Mastertable_viral_rarefied_phages_DMM) %in% names | colnames(Mastertable_viral_rarefied_phages_DMM) == 'lysogenic.cycle'])
rownames(Mastertable_viral_rarefied_phages_DMM_lifestyle) <- Mastertable_viral_rarefied_phages_DMM_lifestyle$lysogenic.cycle
Mastertable_viral_rarefied_phages_DMM_lifestyle$lysogenic.cycle <- NULL

# RA calculations
Mastertable_viral_rarefied_phages_DMM_lifestyle <- sweep(Mastertable_viral_rarefied_phages_DMM_lifestyle, 2, colSums(Mastertable_viral_rarefied_phages_DMM_lifestyle), '/')
Mastertable_viral_rarefied_phages_DMM_lifestyle[is.na(Mastertable_viral_rarefied_phages_DMM_lifestyle)] <- 0

# Merge with clusters
#Mastertable_viral_rarefied_phages_lifestyle <- merge(t(Mastertable_viral_rarefied_phages_DMM_lifestyle),cluster, by = 0, all = F)
#rownames(Mastertable_viral_rarefied_phages_lifestyle) <- Mastertable_viral_rarefied_phages_lifestyle$Row.names
#Mastertable_viral_rarefied_phages_lifestyle$Row.names <- NULL

# Mergr with other metadata
Mastertable_viral_rarefied_phages_lifestyle_ <- merge(t(Mastertable_viral_rarefied_phages_DMM_lifestyle), Metadata, by = 0, all = F)
rownames(Mastertable_viral_rarefied_phages_lifestyle_) <- Mastertable_viral_rarefied_phages_lifestyle_$Row.names
Mastertable_viral_rarefied_phages_lifestyle_$Row.names <- NULL
colnames(Mastertable_viral_rarefied_phages_lifestyle_)

# Melt
Mastertable_viral_rarefied_phages_lifestyle_ <- Mastertable_viral_rarefied_phages_lifestyle_[,c(1:2,6,27,30)]
colnames(Mastertable_viral_rarefied_phages_lifestyle_)
Mastertable_viral_rarefied_phages_lifestyle_cluster <- melt(Mastertable_viral_rarefied_phages_lifestyle_, id.vars = c("cluster", "Diagnosis","Endoscopy.outcome.1"))
Mastertable_viral_rarefied_phages_lifestyle_cluster$cluster[Mastertable_viral_rarefied_phages_lifestyle_cluster$cluster == "1"] <- "Community type CA"
Mastertable_viral_rarefied_phages_lifestyle_cluster$cluster[Mastertable_viral_rarefied_phages_lifestyle_cluster$cluster == "2"] <- "Community type CrM"
Mastertable_viral_rarefied_phages_lifestyle_cluster_CA <- Mastertable_viral_rarefied_phages_lifestyle_cluster[Mastertable_viral_rarefied_phages_lifestyle_cluster$variable == "lysogenic",]
Mastertable_viral_rarefied_phages_lifestyle_cluster_CrM <- Mastertable_viral_rarefied_phages_lifestyle_cluster[c(Mastertable_viral_rarefied_phages_lifestyle_cluster$variable == "lysogenic" & Mastertable_viral_rarefied_phages_lifestyle_cluster$cluster == "Community type CrM"),]
colnames(Mastertable_viral_rarefied_phages_lifestyle_cluster_CA)

## Boxplot
## Lysogenic/ylytic per cluster
Mastertable_viral_rarefied_phages_lifestyle_cluster_CA %>%
 filter(!(Endoscopy.outcome.1 %in% "unknown")) %>%
 ggplot() +
 aes(x = Endoscopy.outcome.1, y = value, 
 fill = cluster) +
 geom_boxplot() +
 scale_fill_hue(direction = 1) +
 theme_bw() +
 facet_wrap(vars(Diagnosis))

View(Mastertable_viral_rarefied_phages_lifestyle_cluster_CA_R)
Mastertable_viral_rarefied_phages_lifestyle_cluster_CA_R <- Mastertable_viral_rarefied_phages_lifestyle_cluster_CA[Mastertable_viral_rarefied_phages_lifestyle_cluster_CA$Endoscopy.outcome.1 == "remission",]
wilcox.test(Mastertable_viral_rarefied_phages_lifestyle_cluster_CA_R$value[Mastertable_viral_rarefied_phages_lifestyle_cluster_CA_R$Diagnosis == "UC"] ~ Mastertable_viral_rarefied_phages_lifestyle_cluster_CA_R$cluster[Mastertable_viral_rarefied_phages_lifestyle_cluster_CA_R$Diagnosis == "UC"], alternative = "two.sided",  paired = FALSE, exact = F)
wilcox.test(Mastertable_viral_rarefied_phages_lifestyle_cluster_CA_R$value[Mastertable_viral_rarefied_phages_lifestyle_cluster_CA_R$Diagnosis == "CD"] ~ Mastertable_viral_rarefied_phages_lifestyle_cluster_CA_R$cluster[Mastertable_viral_rarefied_phages_lifestyle_cluster_CA_R$Diagnosis == "CD"], alternative = "two.sided",  paired = FALSE, exact = F)
p.adjust(0.6198, method = "BH", n = 4) 
p.adjust(0.1742, method = "BH", n = 4) 
effect_size <- wilcox_effsize(Mastertable_viral_rarefied_phages_lifestyle_cluster_CA_R[Mastertable_viral_rarefied_phages_lifestyle_cluster_CA_R$Diagnosis == "UC",],value~cluster)
(effect_size$effsize)
effect_size <- wilcox_effsize(Mastertable_viral_rarefied_phages_lifestyle_cluster_CA_R[Mastertable_viral_rarefied_phages_lifestyle_cluster_CA_R$Diagnosis == "CD",],value~cluster)
(effect_size$effsize)

Mastertable_viral_rarefied_phages_lifestyle_cluster_CA_NR <- Mastertable_viral_rarefied_phages_lifestyle_cluster_CA[Mastertable_viral_rarefied_phages_lifestyle_cluster_CA$Endoscopy.outcome.1 == "non-remission",]
wilcox.test(Mastertable_viral_rarefied_phages_lifestyle_cluster_CA_NR$value[Mastertable_viral_rarefied_phages_lifestyle_cluster_CA_NR$Diagnosis == "UC"] ~ Mastertable_viral_rarefied_phages_lifestyle_cluster_CA_NR$cluster[Mastertable_viral_rarefied_phages_lifestyle_cluster_CA_NR$Diagnosis == "UC"], alternative = "two.sided",  paired = FALSE, exact = F)
wilcox.test(Mastertable_viral_rarefied_phages_lifestyle_cluster_CA_NR$value[Mastertable_viral_rarefied_phages_lifestyle_cluster_CA_NR$Diagnosis == "CD"] ~ Mastertable_viral_rarefied_phages_lifestyle_cluster_CA_NR$cluster[Mastertable_viral_rarefied_phages_lifestyle_cluster_CA_NR$Diagnosis == "CD"], alternative = "two.sided",  paired = FALSE, exact = F)
p.adjust(0.6005, method = "BH", n = 4) 
p.adjust(0.2127, method = "BH", n = 4) 
effect_size <- wilcox_effsize(Mastertable_viral_rarefied_phages_lifestyle_cluster_CA_NR[Mastertable_viral_rarefied_phages_lifestyle_cluster_CA_NR$Diagnosis == "UC",],value~cluster)
(effect_size$effsize)
effect_size <- wilcox_effsize(Mastertable_viral_rarefied_phages_lifestyle_cluster_CA_NR[Mastertable_viral_rarefied_phages_lifestyle_cluster_CA_NR$Diagnosis == "CD",],value~cluster)
(effect_size$effsize)

## Viral community types are not associated  with remission to determine RA of phages 
## Meaning it seems that viral community types don't seem to drive this outcome, but rather endoscopic outcome.

Mastertable_viral_rarefied_phages_lifestyle_cluster_CrM %>%
  filter(cluster %in% "Community type CrM") %>%
  ggplot() +
  aes(x = Endoscopy.outcome.1, y = value, 
      fill = cluster) +
  geom_boxplot() +
  scale_fill_hue(direction = 1) +
  theme_bw() +
  facet_wrap(vars(Diagnosis))

View(Mastertable_viral_rarefied_phages_lifestyle_cluster_CrM)

wilcox.test(Mastertable_viral_rarefied_phages_lifestyle_cluster_CrM$value[c(Mastertable_viral_rarefied_phages_lifestyle_cluster_CrM$Diagnosis == "UC" &! Mastertable_viral_rarefied_phages_lifestyle_cluster_CrM$Endoscopy.outcome.1 == "unknown")] ~ Mastertable_viral_rarefied_phages_lifestyle_cluster_CrM$Endoscopy.outcome.1[c(Mastertable_viral_rarefied_phages_lifestyle_cluster_CrM$Diagnosis == "UC" &! Mastertable_viral_rarefied_phages_lifestyle_cluster_CrM$Endoscopy.outcome.1 == "unknown")], alternative = "two.sided",  paired = FALSE, exact = F)
wilcox.test(Mastertable_viral_rarefied_phages_lifestyle_cluster_CrM$value[c(Mastertable_viral_rarefied_phages_lifestyle_cluster_CrM$Diagnosis == "CD" &! Mastertable_viral_rarefied_phages_lifestyle_cluster_CrM$Endoscopy.outcome.1 == "unknown")] ~ Mastertable_viral_rarefied_phages_lifestyle_cluster_CrM$Endoscopy.outcome.1[c(Mastertable_viral_rarefied_phages_lifestyle_cluster_CrM$Diagnosis == "CD" &! Mastertable_viral_rarefied_phages_lifestyle_cluster_CrM$Endoscopy.outcome.1 == "unknown")], alternative = "two.sided",  paired = FALSE, exact = F)
p.adjust(0.4821, method = "BH", n = 2) # P.adjust = 0.9642
effect_size <- wilcox_effsize(Mastertable_viral_rarefied_phages_lifestyle_cluster_CA[c(Mastertable_viral_rarefied_phages_lifestyle_cluster_CA$Diagnosis == "UC" &! Mastertable_viral_rarefied_phages_lifestyle_cluster_CA$Endoscopy.outcome.1 == "unknown"),],value~Endoscopy.outcome.1)
(effect_size$effsize)

## Conclusion: Samples are showing that viral community type CA in UC is not significantly lower in remission than non-remission.
## Intepration: This means that the lysogenic lifestyle of phages is  driven by viral community type, but not by endoscopic outcome.
## Intepration2: if they would be associated it means that community type change with endoscopic outcome and that endoscopic outcome is driving phage lifestyle. 
## Put this one line of CA, UC statistics in suppl. Tables to refer, no new data.

## CANNOT ASSOCIATE LYSOGENIC POTENTIAL OF ENDOSCOPIC OUTCOME BETWEEN VIRAL COMMUNITY TYPES. 
## Also provided CrM


## Lysogenic/lytic as such
Mastertable_viral_rarefied_phages_lifestyle_cluster %>%
 filter(!(Endoscopy.outcome %in% "unknown")) %>%
 filter(timepoints %in% "primary endpoint") %>%
 filter(variable %in% "lysogenic") %>%
 ggplot() +
 aes(x = variable, y = value, fill = Endoscopy.outcome) +
 geom_boxplot() +
 ylab("Relative abundance") +
 xlab("") +
 scale_fill_manual(values = c(`remission` = "#49ad86", `non-remission` = "#ad4b38")) +
 theme_bw() +
 theme(legend.position="none") +
 facet_wrap(vars(Diagnosis))

# Statistics
# 1) Lysogenic/lysogenic
wilcox.test(Mastertable_viral_rarefied_phages_lifestyle_cluster$value[Mastertable_viral_rarefied_phages_lifestyle_cluster$variable == "lysogenic"] ~ Mastertable_viral_rarefied_phages_lifestyle_cluster$cluster[Mastertable_viral_rarefied_phages_lifestyle_cluster$variable == "lysogenic"], alternative = "two.sided",  paired = FALSE, exact = F)
p.adjust(0.4821, method = "BH", n = 2) # P.adjust = 0.9642
effect_size <- wilcox_effsize(Mastertable_viral_rarefied_phages_lifestyle_cluster[Mastertable_viral_rarefied_phages_lifestyle_cluster$variable == "lysogenic",],value~cluster)
(effect_size$effsize)

wilcox.test(Mastertable_viral_rarefied_phages_lifestyle_cluster$value[Mastertable_viral_rarefied_phages_lifestyle_cluster$variable == "lytic"] ~ Mastertable_viral_rarefied_phages_lifestyle_cluster$cluster[Mastertable_viral_rarefied_phages_lifestyle_cluster$variable == "lytic"], alternative = "two.sided",  paired = FALSE, exact = F)
p.adjust(0.4821, method = "BH", n = 2) # P.adjust = 0.9642
effect_size <- wilcox_effsize(Mastertable_viral_rarefied_phages_lifestyle_cluster[Mastertable_viral_rarefied_phages_lifestyle_cluster$variable == "lytic",],value~cluster)
(effect_size$effsize)

# 2) post-intervention lysogenic/lytic
# UC
Mastertable_viral_rarefied_phages_lifestyle_cluster_UC <-Mastertable_viral_rarefied_phages_lifestyle_cluster[Mastertable_viral_rarefied_phages_lifestyle_cluster$timepoints == "primary endpoint" & Mastertable_viral_rarefied_phages_lifestyle_cluster$variable == "lysogenic",]
wilcox.test(Mastertable_viral_rarefied_phages_lifestyle_cluster_UC$value[Mastertable_viral_rarefied_phages_lifestyle_cluster_UC$Diagnosis == "UC"] ~ Mastertable_viral_rarefied_phages_lifestyle_cluster_UC$Endoscopy.outcome[Mastertable_viral_rarefied_phages_lifestyle_cluster_UC$Diagnosis == "UC"], alternative = "two.sided",  paired = FALSE, exact = F)
p.adjust(0.01719, method = "BH", n = 2) # P.adjust = 0.03438
effect_size <- wilcox_effsize(Mastertable_viral_rarefied_phages_lifestyle_cluster_UC[Mastertable_viral_rarefied_phages_lifestyle_cluster_UC$Diagnosis == "UC",],value~Endoscopy.outcome)
(effect_size$effsize)
effect_size$n1+effect_size$n2

# CD
View(Mastertable_viral_rarefied_phages_lifestyle_cluster_UC)
Mastertable_viral_rarefied_phages_lifestyle_cluster_UC <- Mastertable_viral_rarefied_phages_lifestyle_cluster[!Mastertable_viral_rarefied_phages_lifestyle_cluster$Endoscopy.outcome == "unknown" & Mastertable_viral_rarefied_phages_lifestyle_cluster$Diagnosis == "CD" & Mastertable_viral_rarefied_phages_lifestyle_cluster$timepoints == "primary endpoint" & Mastertable_viral_rarefied_phages_lifestyle_cluster$variable == "lysogenic",]
wilcox.test(value ~ Endoscopy.outcome,data = Mastertable_viral_rarefied_phages_lifestyle_cluster_UC, alternative = "two.sided",  paired = FALSE, exact = )
p.adjust(0.012, method = "BH", n = 2) # P.adjust = 0.03438
effect_size <- wilcox_effsize(Mastertable_viral_rarefied_phages_lifestyle_cluster_UC[Mastertable_viral_rarefied_phages_lifestyle_cluster_UC$Diagnosis == "CD" & !Mastertable_viral_rarefied_phages_lifestyle_cluster_UC$Endoscopy.outcome == "unknown",],value~Endoscopy.outcome)
(effect_size$effsize)
effect_size$n1+effect_size$n2

## Caudoviricetes only



## Make figure same nice
#Put statistcis in Suppl.
Do same for only caudoviricetes phages, subset them first..

###################################
# 3.7 Number of Lysogenic/lytic phages within DMM cluster
###################################
Mastertable_viral_rarefied_phages_DMM <- Mastertable_viral_rarefied_phages
vector_1_first <- (which(names(Mastertable_viral_rarefied_phages_DMM)== "Virsorter"))
vector_1_1_first <- (which(names(Mastertable_viral_rarefied_phages_DMM)== "Virsorter")-1)
vector_2_first <- (which(names(Mastertable_viral_rarefied_phages_DMM)== "Totalnumberofreads"))
vector_ab_viral_NMDS_col_first <- Mastertable_viral_rarefied_phages_DMM[,1:vector_1_1_first]
vector_ab_viral_NMDS_col_tax_first <- Mastertable_viral_rarefied_phages_DMM[,vector_1_first:vector_2_first]
names <- colnames(vector_ab_viral_NMDS_col_first)

# Lysogenic and Lytic
Mastertable_viral_rarefied_phages_DMM$lysogenic.cycle[Mastertable_viral_rarefied_phages_DMM$lysogenic.cycle == "temperate"] <- "lysogenic"
Mastertable_viral_rarefied_phages_DMM$lysogenic.cycle[Mastertable_viral_rarefied_phages_DMM$lysogenic.cycle == "0"] <- "lytic"
table(Mastertable_viral_rarefied_phages_DMM$lysogenic.cycle)

# Extract Class & phages
final_lyso <- (which(names(Mastertable_viral_rarefied_phages_DMM)== "lysogenic.cycle"))
final_class <- (which(names(Mastertable_viral_rarefied_phages_DMM)== "Final_class"))
Mastertable_viral_rarefied_phages_DMM_lyso <- Mastertable_viral_rarefied_phages_DMM[c(final_lyso,final_class)]

# CrAss
#esquisser(Mastertable_viral_rarefied_phages_DMM_lyso)
table(Mastertable_viral_rarefied_phages_DMM_lyso$lysogenic.cycle[Mastertable_viral_rarefied_phages_DMM_lyso$Final_class == "Caudoviricetes (non-CrAss)"])

Mastertable_viral_rarefied_phages_DMM_lyso %>%
 filter(Final_class %in% c("Caudoviricetes (CrAss)", "Malgrandaviricetes", "Caudoviricetes (non-CrAss)")) %>%
 ggplot() +
 aes(x = Final_class, fill = lysogenic.cycle) +
  theme(legend.title = element_blank()) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = c(`lytic` = "#c07572", `lysogenic` = "#a269a9")) +
  scale_color_manual(values = c(`Lytic` = "black", `Lysogenic` = "black")) +
  ylab("") +
  xlab("") +
  theme_classic()

table(Mastertable_viral_rarefied_phages_DMM_lyso$lysogenic.cycle[Mastertable_viral_rarefied_phages_DMM_lyso$Final_class == "Caudoviricetes (CrAss)"])
(8/(66+8))*100        # 10.8%

table(Mastertable_viral_rarefied_phages_DMM_lyso$lysogenic.cycle[Mastertable_viral_rarefied_phages_DMM_lyso$Final_class == "Caudoviricetes (non-CrAss)"])
(1273/(1273+562))*100 # 69.4%

table(Mastertable_viral_rarefied_phages_DMM_lyso$lysogenic.cycle[Mastertable_viral_rarefied_phages_DMM_lyso$Final_class == "Malgrandaviricetes"])
                      # 100%
  
View(Mastertable_viral_rarefied[Mastertable_viral_rarefied$Final_class == "Caudoviricetes (CrAss)" & Mastertable_viral_rarefied$lysogenic.cycle == "temperate",])

# Extract Lysogenic phages & convert to presence/absence table
#Mastertable_viral_rarefied_phages_DMM <- Mastertable_viral_rarefied_phages_DMM[Mastertable_viral_rarefied_phages_DMM$Final_class == "Caudoviricetes (non-CrAss)",]
final_lyso <- (which(names(Mastertable_viral_rarefied_phages_DMM)== "lysogenic.cycle"))
Mastertable_viral_rarefied_phages_DMM_lyso <- Mastertable_viral_rarefied_phages_DMM[c(final_lyso)]
Mastertable_viral_rarefied_phages_DMM_ <-  Mastertable_viral_rarefied_phages_DMM[,c(1:vector_1_1_first)] 
Mastertable_viral_rarefied_phages_DMM_[Mastertable_viral_rarefied_phages_DMM_ > 0] <- 1

# Merge
Mastertable_viral_rarefied_phages_DMM_lyso_ <- merge(Mastertable_viral_rarefied_phages_DMM_,Mastertable_viral_rarefied_phages_DMM_lyso, by = 0, all =F)
rownames(Mastertable_viral_rarefied_phages_DMM_lyso_) <- Mastertable_viral_rarefied_phages_DMM_lyso_$Row.names
Mastertable_viral_rarefied_phages_DMM_lyso_$Row.names <- NULL

# Aggregate
Mastertable_viral_rarefied_phages_DMM_lyso__lifestyle <- aggregate(. ~lysogenic.cycle, FUN = sum, data = Mastertable_viral_rarefied_phages_DMM_lyso_[,colnames(Mastertable_viral_rarefied_phages_DMM_lyso_) %in% names | colnames(Mastertable_viral_rarefied_phages_DMM_lyso_) == 'lysogenic.cycle'])
rownames(Mastertable_viral_rarefied_phages_DMM_lyso__lifestyle) <- Mastertable_viral_rarefied_phages_DMM_lyso__lifestyle$lysogenic.cycle
Mastertable_viral_rarefied_phages_DMM_lyso__lifestyle$lysogenic.cycle <- NULL

# RA calculations
#Mastertable_viral_rarefied_phages_DMM_lyso__lifestyle <- sweep(Mastertable_viral_rarefied_phages_DMM_lyso__lifestyle, 2, colSums(Mastertable_viral_rarefied_phages_DMM_lyso__lifestyle), '/')
#Mastertable_viral_rarefied_phages_DMM_lyso__lifestyle[is.na(Mastertable_viral_rarefied_phages_DMM_lyso__lifestyle)] <- 0

# Merge with clusters
Mastertable_viral_rarefied_phages_lifestyle <- merge(t(Mastertable_viral_rarefied_phages_DMM_lyso__lifestyle),cluster, by = 0, all = F)
rownames(Mastertable_viral_rarefied_phages_lifestyle) <- Mastertable_viral_rarefied_phages_lifestyle$Row.names
Mastertable_viral_rarefied_phages_lifestyle$Row.names <- NULL

# Melt
Mastertable_viral_rarefied_phages_lifestyle_cluster <- melt(Mastertable_viral_rarefied_phages_lifestyle, id.vars = c("cluster"))
View(Mastertable_viral_rarefied_phages_lifestyle_cluster)

# Boxplot
ggplot(Mastertable_viral_rarefied_phages_lifestyle_cluster) +
 aes(x = variable, y = value, fill = variable, color = variable) +
  geom_violin(shape = "circle", trim = FALSE,show.legend = TRUE, draw_quantiles = c(0.5)) +
  scale_fill_manual(values = c(`lytic` = "#9a1919", `lysogenic` = "#15243c")) +
  scale_color_manual(values = c(`lytic` = "black", `lysogenic` = "black")) +
  theme_bw() +
  xlab("") +
  ylab("") +
  facet_wrap(vars(cluster))

# Barplot lytic/lysogneic per class phages
esquisser(Mastertable_viral_rarefied_phages_lifestyle_cluster)
  
ggplot(Mastertable_viral_rarefied_phages_lifestyle_cluster) +
 aes(x = cluster, fill = variable, weight = value) +
 geom_bar(position = "dodge") +
 scale_fill_manual(values = c(lysogenic = "#9a1919", lytic = "#15243c")) +
 coord_flip() +
 theme_minimal()

ggplot(Mastertable_viral_rarefied_phages_lifestyle_cluster) +
 aes(x = cluster, fill = variable, weight = value) +
 geom_bar(position = "fill") +
 scale_fill_manual(values = c(lysogenic = "#9a1919", lytic = "#15243c")) +
 coord_flip() +
 theme_minimal()

# Statistics
# Lysogenic
Mastertable_viral_rarefied_phages_lifestyle_cluster_lyso <- Mastertable_viral_rarefied_phages_lifestyle_cluster[Mastertable_viral_rarefied_phages_lifestyle_cluster$cluster == "1",]
wilcox.test(Mastertable_viral_rarefied_phages_lifestyle_cluster_lyso$value ~ Mastertable_viral_rarefied_phages_lifestyle_cluster_lyso$variable, alternative = "two.sided",  paired = FALSE, exact = F)
p.adjust(0.06165, method = "BH", n = 4) # P.adjust = 0.2466
effect_size <- wilcox_effsize(Mastertable_viral_rarefied_phages_lifestyle_cluster_lyso,value~variable)
(effect_size$effsize)
effect_size$n1
effect_size$n2

Mastertable_viral_rarefied_phages_lifestyle_cluster_lyt <- Mastertable_viral_rarefied_phages_lifestyle_cluster[Mastertable_viral_rarefied_phages_lifestyle_cluster$cluster == "2",]
wilcox.test(Mastertable_viral_rarefied_phages_lifestyle_cluster_lyt$value ~ Mastertable_viral_rarefied_phages_lifestyle_cluster_lyt$variable, alternative = "two.sided",  paired = FALSE, exact = F)
p.adjust(1.715e-06, method = "BH", n = 4) # P.adjust = 6.86e-06
effect_size <- wilcox_effsize(Mastertable_viral_rarefied_phages_lifestyle_cluster_lyt,value~variable)
(effect_size$effsize)
effect_size$n1
effect_size$n2

Mastertable_viral_rarefied_phages_lifestyle_cluster_lyt <- Mastertable_viral_rarefied_phages_lifestyle_cluster[Mastertable_viral_rarefied_phages_lifestyle_cluster$variable == "lysogenic",]
wilcox.test(Mastertable_viral_rarefied_phages_lifestyle_cluster_lyt$value ~ Mastertable_viral_rarefied_phages_lifestyle_cluster_lyt$cluster, alternative = "two.sided",  paired = FALSE, exact = F)
p.adjust(2.958e-08, method = "BH", n = 4) # P.adjust = 1.1832e-07
effect_size <- wilcox_effsize(Mastertable_viral_rarefied_phages_lifestyle_cluster_lyt,value~cluster)
(effect_size$effsize)
effect_size$n1
effect_size$n2

Mastertable_viral_rarefied_phages_lifestyle_cluster_lyt <- Mastertable_viral_rarefied_phages_lifestyle_cluster[Mastertable_viral_rarefied_phages_lifestyle_cluster$variable == "lytic",]
wilcox.test(Mastertable_viral_rarefied_phages_lifestyle_cluster_lyt$value ~ Mastertable_viral_rarefied_phages_lifestyle_cluster_lyt$cluster, alternative = "two.sided",  paired = FALSE, exact = F)
p.adjust(2.2e-16, method = "BH", n = 4) # P.adjust = <  2.2e-16
effect_size <- wilcox_effsize(Mastertable_viral_rarefied_phages_lifestyle_cluster_lyt,value~cluster)
(effect_size$effsize)
effect_size$n1
effect_size$n2

# Scatterplot
#esquisser(Mastertable_viral_rarefied_phages_lifestyle)
ggplot(Mastertable_viral_rarefied_phages_lifestyle) +
  aes(x = lysogenic, y = lytic, colour = cluster) +
  geom_point(shape = "circle", size = 1.5) +
  ylab("Lytic phages") +
  xlab("Lysogenic phages") +
  geom_smooth(se = FALSE, method = lm) +
  scale_color_manual(values = c(`1` = "#990011FF", `2` = "#00203FFF")) +
  theme_classic() +
  geom_smooth(aes(shape = NULL), method = lm, se = FALSE)

ggplot(Mastertable_viral_rarefied_phages_lifestyle) +
  aes(x = Mastertable_viral_rarefied_phages_lifestyle_cluster$value[Mastertable_viral_rarefied_phages_lifestyle_cluster$variable == "lysogenic"], y = Mastertable_viral_rarefied_phages_lifestyle_cluster$value[Mastertable_viral_rarefied_phages_lifestyle_cluster$variable == "lytic"], colour = cluster) +
  geom_point(shape = "circle", size = 1.5) +
  ylab("Lytic phages") +
  xlab("Lysogenic phages") +
  geom_smooth(se = FALSE, method = lm) +
  scale_color_manual(values = c(`1` = "#990011FF", `2` = "#00203FFF")) +
  theme_classic() +
  geom_smooth(aes(shape = NULL), method = lm, se = FALSE)

ggplot(Mastertable_viral_rarefied_phages_lifestyle, aes(x = lysogenic, y= lytic, color = cluster, fill = cluster, shape = Diagnosis)) +
  geom_point() +
  scale_shape_manual(values = c(1,16)) +
   ylab("Lytic phages") +
  xlab("Lysogenic phages") +
  scale_color_manual(values = c(`1` = "#990011FF", `2` = "#00203FFF")) +
  scale_fill_manual(values = c(`1` = "#990011FF", `2` = "#00203FFF")) +
  theme_classic() +
  geom_smooth(aes(shape = NULL), method = lm, se = FALSE)

ggplot(Mastertable_viral_rarefied_phages_lifestyle, aes(x = lysogenic, y= lytic, color = cluster, fill = cluster, shape = timepoints)) +
  geom_point() +
  scale_shape_manual(values = c(2, 1,16)) +
  ylab("Lytic phages") +
  xlab("Lysogenic phages") +
  scale_color_manual(values = c(`1` = "#990011FF", `2` = "#00203FFF")) +
  scale_fill_manual(values = c(`1` = "#990011FF", `2` = "#00203FFF")) +
  theme_classic() +
  geom_smooth(aes(shape = NULL), method = lm, se = FALSE)

# Statistics: linear regression
# Add trendline to figure
model <- lm(lytic ~ lysogenic , data = Mastertable_viral_rarefied_phages_lifestyle)
summary(model)
# Multiple R-squared: 0.2452, Adjusted R-squared: 0.2431, F-statistic: 117.3, p-value = < 2.2e-16, Residual standard error: 8.986

model <- lm(lytic ~ lysogenic , data = Mastertable_viral_rarefied_phages_lifestyle[Mastertable_viral_rarefied_phages_lifestyle$cluster == "1",])
summary(model)
# Multiple R-squared: 0.1949, Adjusted R-squared: 0.1916, F-statistic: 58.11, p-value =  5.741e-13, Residual standard error: 6.331

model <- lm(lytic ~ lysogenic , data = Mastertable_viral_rarefied_phages_lifestyle[Mastertable_viral_rarefied_phages_lifestyle$cluster == "2",])
summary(model)
# Multiple R-squared: 0.1757, Adjusted R-squared: 0.1688, F-statistic: 25.37, p-value =  1.703e-06, Residual standard error: 9.348
###################################
# 3.8 Host prediction within DMM cluster: RA at phyla
###################################
Mastertable_viral_rarefied_phages_DMM <- Mastertable_viral_rarefied_phages
vector_1_first <- (which(names(Mastertable_viral_rarefied_phages_DMM)== "Virsorter"))
vector_1_1_first <- (which(names(Mastertable_viral_rarefied_phages_DMM)== "Virsorter")-1)
vector_2_first <- (which(names(Mastertable_viral_rarefied_phages_DMM)== "Totalnumberofreads"))
vector_ab_viral_NMDS_col_first <- Mastertable_viral_rarefied_phages_DMM[,1:vector_1_1_first]
vector_ab_viral_NMDS_col_tax_first <- Mastertable_viral_rarefied_phages_DMM[,vector_1_first:vector_2_first]
names <- colnames(vector_ab_viral_NMDS_col_first)

# Merge with Host
Mastertable_viral_rarefied_phages_DMM <- merge(Mastertable_viral_rarefied_phages_DMM,Host_phylum_taxonomy, by = 0, all = F)
rownames(Mastertable_viral_rarefied_phages_DMM) <- Mastertable_viral_rarefied_phages_DMM$Row.names
Mastertable_viral_rarefied_phages_DMM$Row.names <- NULL
table(Mastertable_viral_rarefied_phages_DMM$phylum)

# Aggregate based on lytic/lysogenic
Mastertable_viral_rarefied_phages_DMM_lifestyle <- aggregate(. ~phylum, FUN = sum, data = Mastertable_viral_rarefied_phages_DMM[,colnames(Mastertable_viral_rarefied_phages_DMM) %in% names | colnames(Mastertable_viral_rarefied_phages_DMM) == 'phylum'])
rownames(Mastertable_viral_rarefied_phages_DMM_lifestyle) <- Mastertable_viral_rarefied_phages_DMM_lifestyle$phylum
Mastertable_viral_rarefied_phages_DMM_lifestyle$v <- NULL
View(Mastertable_viral_rarefied_phages_DMM_lifestyle)

# RA calculations
Mastertable_viral_rarefied_phages_DMM_lifestyle <- sweep(Mastertable_viral_rarefied_phages_DMM_lifestyle, 2, colSums(Mastertable_viral_rarefied_phages_DMM_lifestyle), '/')
Mastertable_viral_rarefied_phages_DMM_lifestyle[is.na(Mastertable_viral_rarefied_phages_DMM_lifestyle)] <- 0

# Merge with clusters
Mastertable_viral_rarefied_phages_lifestyle <- merge(t(Mastertable_viral_rarefied_phages_DMM_lifestyle),cluster, by = 0, all = F)
rownames(Mastertable_viral_rarefied_phages_lifestyle) <- Mastertable_viral_rarefied_phages_lifestyle$Row.names
Mastertable_viral_rarefied_phages_lifestyle$Row.names <- NULL

# Melt
Mastertable_viral_rarefied_phages_lifestyle_cluster <- melt(Mastertable_viral_rarefied_phages_lifestyle, id.vars = c("cluster"))
str(Mastertable_viral_rarefied_phages_lifestyle_cluster)

## Boxplot
Mastertable_viral_rarefied_phages_lifestyle_cluster %>%
 filter(variable %in% c("Actinobacteria", "Bacteroidetes", 
"Chlamydiae", "Firmicutes", "Proteobacteria")) %>%
 ggplot() +
 aes(x = cluster, y = value, fill = cluster) +
 geom_boxplot(shape = "circle") +
 scale_fill_manual(values = c(`1` = "#990011FF", `2` = "#00203FFF")) +
 theme_bw() +
 facet_wrap(vars(variable), scales = "free_y")

# Statistics
wilcox.test(Mastertable_viral_rarefied_phages_lifestyle_cluster$value[Mastertable_viral_rarefied_phages_lifestyle_cluster$variable == "Proteobacteria"] ~ Mastertable_viral_rarefied_phages_lifestyle_cluster$cluster[Mastertable_viral_rarefied_phages_lifestyle_cluster$variable == "Proteobacteria"], alternative = "two.sided",  paired = FALSE, exact = F)
wilcox.test(Mastertable_viral_rarefied_phages_lifestyle_cluster$value[Mastertable_viral_rarefied_phages_lifestyle_cluster$variable == "Bacteroidetes"] ~ Mastertable_viral_rarefied_phages_lifestyle_cluster$cluster[Mastertable_viral_rarefied_phages_lifestyle_cluster$variable == "Proteobacteria"], alternative = "two.sided",  paired = FALSE, exact = F)
wilcox.test(Mastertable_viral_rarefied_phages_lifestyle_cluster$value[Mastertable_viral_rarefied_phages_lifestyle_cluster$variable == "Firmicutes"] ~ Mastertable_viral_rarefied_phages_lifestyle_cluster$cluster[Mastertable_viral_rarefied_phages_lifestyle_cluster$variable == "Proteobacteria"], alternative = "two.sided",  paired = FALSE, exact = F)
wilcox.test(Mastertable_viral_rarefied_phages_lifestyle_cluster$value[Mastertable_viral_rarefied_phages_lifestyle_cluster$variable == "Chlamydiae"] ~ Mastertable_viral_rarefied_phages_lifestyle_cluster$cluster[Mastertable_viral_rarefied_phages_lifestyle_cluster$variable == "Proteobacteria"], alternative = "two.sided",  paired = FALSE, exact = F)
wilcox.test(Mastertable_viral_rarefied_phages_lifestyle_cluster$value[Mastertable_viral_rarefied_phages_lifestyle_cluster$variable == "Actinobacteria"] ~ Mastertable_viral_rarefied_phages_lifestyle_cluster$cluster[Mastertable_viral_rarefied_phages_lifestyle_cluster$variable == "Proteobacteria"], alternative = "two.sided",  paired = FALSE, exact = F)
p.adjust(0.1877, method = "BH", n = 5) 
p.adjust(0.1127, method = "BH", n = 5) 
p.adjust(0.8815, method = "BH", n = 5) 
p.adjust(3.779e-13, method = "BH", n = 5)
p.adjust(1.702e-05, method = "BH", n = 5)
# Note that I remove the Host with < 0.01 RA
###################################
# 3.9 Host prediction within DMM cluster: RA at genus level
###################################
Mastertable_viral_rarefied_phages_DMM <- Mastertable_viral_rarefied_phages
vector_1_first <- (which(names(Mastertable_viral_rarefied_phages_DMM)== "Virsorter"))
vector_1_1_first <- (which(names(Mastertable_viral_rarefied_phages_DMM)== "Virsorter")-1)
vector_2_first <- (which(names(Mastertable_viral_rarefied_phages_DMM)== "Totalnumberofreads"))
vector_ab_viral_NMDS_col_first <- Mastertable_viral_rarefied_phages_DMM[,1:vector_1_1_first]
vector_ab_viral_NMDS_col_tax_first <- Mastertable_viral_rarefied_phages_DMM[,vector_1_first:vector_2_first]
names <- colnames(vector_ab_viral_NMDS_col_first)

# Merge with Host genus
Mastertable_viral_rarefied_phages_DMM <- merge(Mastertable_viral_rarefied_phages_DMM,Host_phylum_genus, by = 0, all = F)
rownames(Mastertable_viral_rarefied_phages_DMM) <- Mastertable_viral_rarefied_phages_DMM$Row.names
Mastertable_viral_rarefied_phages_DMM$Row.names <- NULL
table(Mastertable_viral_rarefied_phages_DMM$genus)

# Aggregate based on lytic/lysogenic
Mastertable_viral_rarefied_phages_DMM_lifestyle <- aggregate(. ~genus, FUN = sum, data = Mastertable_viral_rarefied_phages_DMM[,colnames(Mastertable_viral_rarefied_phages_DMM) %in% names | colnames(Mastertable_viral_rarefied_phages_DMM) == 'genus'])
rownames(Mastertable_viral_rarefied_phages_DMM_lifestyle) <- Mastertable_viral_rarefied_phages_DMM_lifestyle$genus
Mastertable_viral_rarefied_phages_DMM_lifestyle$genus <- NULL

# RA calculations
Mastertable_viral_rarefied_phages_DMM_lifestyle <- sweep(Mastertable_viral_rarefied_phages_DMM_lifestyle, 2, colSums(Mastertable_viral_rarefied_phages_DMM_lifestyle), '/')
Mastertable_viral_rarefied_phages_DMM_lifestyle[is.na(Mastertable_viral_rarefied_phages_DMM_lifestyle)] <- 0

# Merge with clusters
Mastertable_viral_rarefied_phages_lifestyle <- merge(t(Mastertable_viral_rarefied_phages_DMM_lifestyle),cluster, by = 0, all = F)
rownames(Mastertable_viral_rarefied_phages_lifestyle) <- Mastertable_viral_rarefied_phages_lifestyle$Row.names
Mastertable_viral_rarefied_phages_lifestyle$Row.names <- NULL

# Melt
Mastertable_viral_rarefied_phages_lifestyle_cluster <- melt(Mastertable_viral_rarefied_phages_lifestyle, id.vars = c("cluster"))
str(Mastertable_viral_rarefied_phages_lifestyle_cluster)

## Boxplot
#esquisser(Mastertable_viral_rarefied_phages_lifestyle_cluster)

Mastertable_viral_rarefied_phages_lifestyle_cluster %>%
 ggplot() +
 aes(x = cluster, y = value, fill = cluster) +
  geom_boxplot(shape = "circle") +
  geom_jitter(width = 0.15) +
# geom_violin(shape = "circle", trim = TRUE) +
 # geom_boxplot(width=0.05, color = "black") +
 #  stat_summary(fun.data = "median_hilow", geom = "crossbar",
              # colour = "white", width = 0.03) +
  scale_fill_manual(values = c(`1` = "#990011FF", `2` = "#00203FFF")) +
  theme_bw() +
 facet_wrap(vars(variable))

# Cutoff to show, total RA > 1%
table(Mastertable_viral_rarefied_phages_lifestyle_cluster$variable)
median(Mastertable_viral_rarefied_phages_lifestyle_cluster$value[Mastertable_viral_rarefied_phages_lifestyle_cluster$variable == "Bacteroides" & Mastertable_viral_rarefied_phages_lifestyle_cluster$cluster == "1"])*100 # or 7,39%
median(Mastertable_viral_rarefied_phages_lifestyle_cluster$value[Mastertable_viral_rarefied_phages_lifestyle_cluster$variable == "Bacteroides" & Mastertable_viral_rarefied_phages_lifestyle_cluster$cluster == "2"])*100 # or 7.74%
median(Mastertable_viral_rarefied_phages_lifestyle_cluster$value[Mastertable_viral_rarefied_phages_lifestyle_cluster$variable == "Prevotella" & Mastertable_viral_rarefied_phages_lifestyle_cluster$cluster == "1"])*100 # or 0.20%
median(Mastertable_viral_rarefied_phages_lifestyle_cluster$value[Mastertable_viral_rarefied_phages_lifestyle_cluster$variable == "Prevotella" & Mastertable_viral_rarefied_phages_lifestyle_cluster$cluster == "2"])*100 # or 10.14%

# Statistics
# Bacteroides
wilcox.test(Mastertable_viral_rarefied_phages_lifestyle_cluster$value[Mastertable_viral_rarefied_phages_lifestyle_cluster$variable == "Bacteroides"] ~ Mastertable_viral_rarefied_phages_lifestyle_cluster$cluster[Mastertable_viral_rarefied_phages_lifestyle_cluster$variable == "Bacteroides"], alternative = "two.sided",  paired = FALSE, exact = F)
p.adjust(0.92, method = "BH", n = 2) # P.adjust = 0.9642
effect_size <- wilcox_effsize(Mastertable_viral_rarefied_phages_lifestyle_cluster[Mastertable_viral_rarefied_phages_lifestyle_cluster$variable == "Bacteroides",],value~cluster)
(effect_size$effsize)

# Prevotella
wilcox.test(Mastertable_viral_rarefied_phages_lifestyle_cluster$value[Mastertable_viral_rarefied_phages_lifestyle_cluster$variable == "Prevotella"] ~ Mastertable_viral_rarefied_phages_lifestyle_cluster$cluster[Mastertable_viral_rarefied_phages_lifestyle_cluster$variable == "Prevotella"], alternative = "two.sided",  paired = FALSE, exact = F)
p.adjust(7.546e-06, method = "BH", n = 2) # P.adjust = 1.5092e-05
effect_size <- wilcox_effsize(Mastertable_viral_rarefied_phages_lifestyle_cluster[Mastertable_viral_rarefied_phages_lifestyle_cluster$variable == "Prevotella",],value~cluster)
(effect_size$effsize)
###################################
# 3.10 Host prediction: percentage for cluster 1 (donut plot)
###################################
## Donut plot
Mastertable_viral_rarefied_phages_DMM <- Mastertable_viral_rarefied_phages
vector_1_first <- (which(names(Mastertable_viral_rarefied_phages_DMM)== "Virsorter"))
vector_1_1_first <- (which(names(Mastertable_viral_rarefied_phages_DMM)== "Virsorter")-1)
vector_2_first <- (which(names(Mastertable_viral_rarefied_phages_DMM)== "Totalnumberofreads"))
vector_ab_viral_NMDS_col_first <- Mastertable_viral_rarefied_phages_DMM[,1:vector_1_1_first]
vector_ab_viral_NMDS_col_tax_first <- Mastertable_viral_rarefied_phages_DMM[,vector_1_first:vector_2_first]
names <- colnames(vector_ab_viral_NMDS_col_first)

# Extract cluster
Mastertable_viral_rarefied_phages_DMM_cluster <- merge(t(vector_ab_viral_NMDS_col_first),cluster,  by = 0, all = F)
rownames(Mastertable_viral_rarefied_phages_DMM_cluster) <- Mastertable_viral_rarefied_phages_DMM_cluster$Row.names
Mastertable_viral_rarefied_phages_DMM_cluster$Row.names <- NULL
Mastertable_viral_rarefied_phages_DMM_cluster_1 <- Mastertable_viral_rarefied_phages_DMM_cluster[Mastertable_viral_rarefied_phages_DMM_cluster$cluster == "1",]
Mastertable_viral_rarefied_phages_DMM_cluster_1$cluster <- NULL
Mastertable_viral_rarefied_phages_DMM_cluster_1 <- as.data.frame(t(Mastertable_viral_rarefied_phages_DMM_cluster_1))

# Merge with Host genus
Mastertable_viral_rarefied_phages_DMM_cluster_1$Totalnumberofreads <- rowSums(Mastertable_viral_rarefied_phages_DMM_cluster_1)
Mastertable_viral_rarefied_phages_DMM <- merge(Mastertable_viral_rarefied_phages_DMM_cluster_1,Host_phylum_genus, by = 0, all = T)
rownames(Mastertable_viral_rarefied_phages_DMM) <- Mastertable_viral_rarefied_phages_DMM$Row.names
Mastertable_viral_rarefied_phages_DMM$Row.names <- NULL
Mastertable_viral_rarefied_phages_DMM$genus[is.na(Mastertable_viral_rarefied_phages_DMM$genus)] <- "unknown"
table(Mastertable_viral_rarefied_phages_DMM$genus)

# Extract Host & reads
vector_1_1_first <- (which(names(Mastertable_viral_rarefied_phages_DMM)== "genus"))
vector_2_first <- (which(names(Mastertable_viral_rarefied_phages_DMM)== "Totalnumberofreads"))
Mastertable_viral_rarefied_phages_DMM_donut_plot <- Mastertable_viral_rarefied_phages_DMM[,vector_1_1_first:vector_2_first]

# Aggregate
names <- (which(names(Mastertable_viral_rarefied_phages_DMM_donut_plot)== "Totalnumberofreads"))
Mastertable_viral_rarefied_phages_DMM_donut_plot_ <- aggregate(. ~genus, FUN = sum, data = Mastertable_viral_rarefied_phages_DMM_donut_plot)
rownames(Mastertable_viral_rarefied_phages_DMM_donut_plot_) <- Mastertable_viral_rarefied_phages_DMM_donut_plot_$genus
str(Mastertable_viral_rarefied_phages_DMM_donut_plot_)

# Proportions
Mastertable_viral_rarefied_phages_DMM_donut_plot_$percentage <- (Mastertable_viral_rarefied_phages_DMM_donut_plot_$Totalnumberofreads/sum(Mastertable_viral_rarefied_phages_DMM_donut_plot_$Totalnumberofreads))*100
str(Mastertable_viral_rarefied_phages_DMM_donut_plot_)

# Remove insanely small percentages
Mastertable_viral_rarefied_phages_DMM_donut_plot_$genus[Mastertable_viral_rarefied_phages_DMM_donut_plot_$percentage < 1] <- "Other"
Mastertable_viral_rarefied_phages_DMM_donut_plot_ <- aggregate(. ~genus, FUN = sum, data = Mastertable_viral_rarefied_phages_DMM_donut_plot_)
rownames(Mastertable_viral_rarefied_phages_DMM_donut_plot_) <- Mastertable_viral_rarefied_phages_DMM_donut_plot_$genus
Mastertable_viral_rarefied_phages_DMM_donut_plot_1 <- Mastertable_viral_rarefied_phages_DMM_donut_plot_
str(Mastertable_viral_rarefied_phages_DMM_donut_plot_2)
View(Mastertable_viral_rarefied_phages_DMM_donut_plot_2)

# Order
x <- c("unknown","Prevotella","Bacteroides","Other","Escherichia","Chlamydia","Veillonella","Odoribacter","Lactococcus","Clostridium","Parabacteroides","Blautia")
Mastertable_viral_rarefied_phages_DMM_donut_plot_2$genus <- factor(Mastertable_viral_rarefied_phages_DMM_donut_plot_2$genus,levels=x)
str(Mastertable_viral_rarefied_phages_DMM_donut_plot_2)

# Donut plot 
ggplot(Mastertable_viral_rarefied_phages_DMM_donut_plot_2,aes(x = 2, y = percentage, fill = genus)) +
  geom_bar(stat = "identity", color = "black") +
  coord_polar("y", start = 1) +
  # geom_text(aes(label = paste0(round(values), "%")),
  #  position = position_stack(vjust = 0.5)) +
  theme_void() +
  xlim(.2,2.5) +
  # scale_fill_brewer(palette = "Greys") +
  theme(legend.title = element_blank()) 
#theme(legend.position = "none") 
###################################
# 3.11 Host prediction: percentage for cluster 2 (donut plot)
###################################
## Donut plot
Mastertable_viral_rarefied_phages_DMM <- Mastertable_viral_rarefied_phages
vector_1_first <- (which(names(Mastertable_viral_rarefied_phages_DMM)== "Virsorter"))
vector_1_1_first <- (which(names(Mastertable_viral_rarefied_phages_DMM)== "Virsorter")-1)
vector_2_first <- (which(names(Mastertable_viral_rarefied_phages_DMM)== "Totalnumberofreads"))
vector_ab_viral_NMDS_col_first <- Mastertable_viral_rarefied_phages_DMM[,1:vector_1_1_first]
vector_ab_viral_NMDS_col_tax_first <- Mastertable_viral_rarefied_phages_DMM[,vector_1_first:vector_2_first]
names <- colnames(vector_ab_viral_NMDS_col_first)

# Extract cluster
Mastertable_viral_rarefied_phages_DMM_cluster <- merge(t(vector_ab_viral_NMDS_col_first),cluster,  by = 0, all = F)
rownames(Mastertable_viral_rarefied_phages_DMM_cluster) <- Mastertable_viral_rarefied_phages_DMM_cluster$Row.names
Mastertable_viral_rarefied_phages_DMM_cluster$Row.names <- NULL
Mastertable_viral_rarefied_phages_DMM_cluster_1 <- Mastertable_viral_rarefied_phages_DMM_cluster[Mastertable_viral_rarefied_phages_DMM_cluster$cluster == "2",]
Mastertable_viral_rarefied_phages_DMM_cluster_1$cluster <- NULL
Mastertable_viral_rarefied_phages_DMM_cluster_1 <- as.data.frame(t(Mastertable_viral_rarefied_phages_DMM_cluster_1))

# Merge with Host genus
Mastertable_viral_rarefied_phages_DMM_cluster_1$Totalnumberofreads <- rowSums(Mastertable_viral_rarefied_phages_DMM_cluster_1)
Mastertable_viral_rarefied_phages_DMM <- merge(Mastertable_viral_rarefied_phages_DMM_cluster_1,Host_phylum_genus, by = 0, all = T)
rownames(Mastertable_viral_rarefied_phages_DMM) <- Mastertable_viral_rarefied_phages_DMM$Row.names
Mastertable_viral_rarefied_phages_DMM$Row.names <- NULL
Mastertable_viral_rarefied_phages_DMM$genus[is.na(Mastertable_viral_rarefied_phages_DMM$genus)] <- "unknown"
table(Mastertable_viral_rarefied_phages_DMM$genus)

# Extract Host & reads
vector_1_1_first <- (which(names(Mastertable_viral_rarefied_phages_DMM)== "genus"))
vector_2_first <- (which(names(Mastertable_viral_rarefied_phages_DMM)== "Totalnumberofreads"))
Mastertable_viral_rarefied_phages_DMM_donut_plot <- Mastertable_viral_rarefied_phages_DMM[,vector_1_1_first:vector_2_first]

# Aggregate
names <- (which(names(Mastertable_viral_rarefied_phages_DMM_donut_plot)== "Totalnumberofreads"))
Mastertable_viral_rarefied_phages_DMM_donut_plot_ <- aggregate(. ~genus, FUN = sum, data = Mastertable_viral_rarefied_phages_DMM_donut_plot)
rownames(Mastertable_viral_rarefied_phages_DMM_donut_plot_) <- Mastertable_viral_rarefied_phages_DMM_donut_plot_$genus
str(Mastertable_viral_rarefied_phages_DMM_donut_plot_)

# Proportions
Mastertable_viral_rarefied_phages_DMM_donut_plot_$percentage <- (Mastertable_viral_rarefied_phages_DMM_donut_plot_$Totalnumberofreads/sum(Mastertable_viral_rarefied_phages_DMM_donut_plot_$Totalnumberofreads))*100
str(Mastertable_viral_rarefied_phages_DMM_donut_plot_)

# Remove insanely small percentages
Mastertable_viral_rarefied_phages_DMM_donut_plot_$genus[Mastertable_viral_rarefied_phages_DMM_donut_plot_$percentage < 1] <- "Other"
Mastertable_viral_rarefied_phages_DMM_donut_plot_ <- aggregate(. ~genus, FUN = sum, data = Mastertable_viral_rarefied_phages_DMM_donut_plot_)
rownames(Mastertable_viral_rarefied_phages_DMM_donut_plot_) <- Mastertable_viral_rarefied_phages_DMM_donut_plot_$genus
str(Mastertable_viral_rarefied_phages_DMM_donut_plot_)
Mastertable_viral_rarefied_phages_DMM_donut_plot_1 <- Mastertable_viral_rarefied_phages_DMM_donut_plot_
View(Mastertable_viral_rarefied_phages_DMM_donut_plot_1)

# Order
x <- c("unknown","Prevotella","Bacteroides","Chlamydia","Other","Alistipes","Parabacteroides","Escherichia","Odoribacter","Blautia","Akkermansia","Lactococcus","Veillonella")
Mastertable_viral_rarefied_phages_DMM_donut_plot_$genus <- factor(Mastertable_viral_rarefied_phages_DMM_donut_plot_$genus,levels=x)
str(Mastertable_viral_rarefied_phages_DMM_donut_plot_)

# Donut plot 
ggplot(Mastertable_viral_rarefied_phages_DMM_donut_plot_,aes(x = 2, y = percentage, fill = genus)) +
  geom_bar(stat = "identity", color = "black") +
  coord_polar("y", start = 1) +
  # geom_text(aes(label = paste0(round(values), "%")),
  #  position = position_stack(vjust = 0.5)) +
  theme_void() +
  xlim(.2,2.5) +
  # scale_fill_brewer(palette = "Greys") +
  theme(legend.title = element_blank()) 
#theme(legend.position = "none") 

## Statistics
### chi-squared what? <-> LB thesis same there..

###################################
# 3.12 Host prediction: percentage for all clusters
###################################
## Donut plot
Mastertable_viral_rarefied_phages_DMM <- Mastertable_viral_rarefied_phages
vector_1_first <- (which(names(Mastertable_viral_rarefied_phages_DMM)== "Virsorter"))
vector_1_1_first <- (which(names(Mastertable_viral_rarefied_phages_DMM)== "Virsorter")-1)
vector_2_first <- (which(names(Mastertable_viral_rarefied_phages_DMM)== "Totalnumberofreads"))
vector_ab_viral_NMDS_col_first <- Mastertable_viral_rarefied_phages_DMM[,1:vector_1_1_first]
vector_ab_viral_NMDS_col_tax_first <- Mastertable_viral_rarefied_phages_DMM[,vector_1_first:vector_2_first]
names <- colnames(vector_ab_viral_NMDS_col_first)

# Merge with Host genus
vector_ab_viral_NMDS_col_first$Totalnumberofreads <- rowSums(vector_ab_viral_NMDS_col_first)
Mastertable_viral_rarefied_phages_DMM <- merge(vector_ab_viral_NMDS_col_first,Host_phylum_genus, by = 0, all = T)
rownames(Mastertable_viral_rarefied_phages_DMM) <- Mastertable_viral_rarefied_phages_DMM$Row.names
Mastertable_viral_rarefied_phages_DMM$Row.names <- NULL
Mastertable_viral_rarefied_phages_DMM$genus[is.na(Mastertable_viral_rarefied_phages_DMM$genus)] <- "unknown"
table(Mastertable_viral_rarefied_phages_DMM$genus)

# Extract Host & reads
View(Mastertable_viral_rarefied_phages_DMM)
Mastertable_viral_rarefied_phages_DMM$total <- 1 # to obtain contigs instead of RA
vector_1_1_first <- (which(names(Mastertable_viral_rarefied_phages_DMM)== "genus"))
#vector_2_first <- (which(names(Mastertable_viral_rarefied_phages_DMM)== "Totalnumberofreads"))
vector_2_first <- (which(names(Mastertable_viral_rarefied_phages_DMM)== "total"))
Mastertable_viral_rarefied_phages_DMM_donut_plot <- Mastertable_viral_rarefied_phages_DMM[,vector_1_1_first:vector_2_first]
View(Mastertable_viral_rarefied_phages_DMM_donut_plot)

# Aggregate
#names <- (which(names(Mastertable_viral_rarefied_phages_DMM_donut_plot)== "Totalnumberofreads"))
names <- (which(names(Mastertable_viral_rarefied_phages_DMM_donut_plot)== "total"))
Mastertable_viral_rarefied_phages_DMM_donut_plot_ <- aggregate(. ~genus, FUN = sum, data = Mastertable_viral_rarefied_phages_DMM_donut_plot)
rownames(Mastertable_viral_rarefied_phages_DMM_donut_plot_) <- Mastertable_viral_rarefied_phages_DMM_donut_plot_$genus
str(Mastertable_viral_rarefied_phages_DMM_donut_plot_)
#View(Mastertable_viral_rarefied_phages_DMM_donut_plot_)

# Proportions
#Mastertable_viral_rarefied_phages_DMM_donut_plot_$percentage <- (Mastertable_viral_rarefied_phages_DMM_donut_plot_$Totalnumberofreads/sum(Mastertable_viral_rarefied_phages_DMM_donut_plot_$Totalnumberofreads))*100
Mastertable_viral_rarefied_phages_DMM_donut_plot_$percentage <- (Mastertable_viral_rarefied_phages_DMM_donut_plot_$total/sum(Mastertable_viral_rarefied_phages_DMM_donut_plot_$total))*100
#View(Mastertable_viral_rarefied_phages_DMM_donut_plot_)

# Remove insanely small percentages
#Mastertable_viral_rarefied_phages_DMM_donut_plot_$genus[Mastertable_viral_rarefied_phages_DMM_donut_plot_$percentage < 1] <- "Other"
Mastertable_viral_rarefied_phages_DMM_donut_plot_$genus[!c(Mastertable_viral_rarefied_phages_DMM_donut_plot_$genus == "Prevotella" | Mastertable_viral_rarefied_phages_DMM_donut_plot_$genus == "Bacteroides" | Mastertable_viral_rarefied_phages_DMM_donut_plot_$genus == "unknown")] <- "Other"
View(Mastertable_viral_rarefied_phages_DMM_donut_plot_)
Mastertable_viral_rarefied_phages_DMM_donut_plot_ <- aggregate(. ~genus, FUN = sum, data = Mastertable_viral_rarefied_phages_DMM_donut_plot_)
rownames(Mastertable_viral_rarefied_phages_DMM_donut_plot_) <- Mastertable_viral_rarefied_phages_DMM_donut_plot_$genus
str(Mastertable_viral_rarefied_phages_DMM_donut_plot_)
Mastertable_viral_rarefied_phages_DMM_donut_plot_1 <- Mastertable_viral_rarefied_phages_DMM_donut_plot_
View(Mastertable_viral_rarefied_phages_DMM_donut_plot_1)

# Order
x <- c("unknown","Prevotella","Bacteroides","Other","Escherichia","Chlamydia","Veillonella","Parabacteroides","Odoribacter","Alistipes","Lactococcus","Clostridium","Blautia")
Mastertable_viral_rarefied_phages_DMM_donut_plot_$genus <- factor(Mastertable_viral_rarefied_phages_DMM_donut_plot_$genus,levels=x)
str(Mastertable_viral_rarefied_phages_DMM_donut_plot_)
View(Mastertable_viral_rarefied_phages_DMM_donut_plot_)
100-58.093567-11.247811-8.115700

# Donut plot 
ggplot(Mastertable_viral_rarefied_phages_DMM_donut_plot_,aes(x = 2, y = percentage, fill = genus)) +
  geom_bar(stat = "identity", color = "black") +
  coord_polar("y", start = 1) +
  # geom_text(aes(label = paste0(round(values), "%")),
  #  position = position_stack(vjust = 0.5)) +
  theme_void() +
  xlim(-1,2.5) +
  # scale_fill_brewer(palette = "Greys") +
  theme(legend.title = element_blank()) 
#theme(legend.position = "none") 

# Change to others
Mastertable_viral_rarefied_phages_DMM_donut_plot_1 = Mastertable_viral_rarefied_phages_DMM_donut_plot_[-c(1,3:8,10,13),]
100-58.093567-11.247811-8.115700
Mastertable_viral_rarefied_phages_DMM_donut_plot_1$percentage[Mastertable_viral_rarefied_phages_DMM_donut_plot_1$genus == "Other"] <- 22.54292

# Donut plot 
ggplot(Mastertable_viral_rarefied_phages_DMM_donut_plot_1,aes(x = 2, y = percentage, fill = genus)) +
  geom_bar(stat = "identity", color = "black") +
  coord_polar("y", start = 1) +
  # geom_text(aes(label = paste0(round(values), "%")),
  #  position = position_stack(vjust = 0.5)) +
  theme_void() +
  xlim(-0.5,2.5) 
  # scale_fill_brewer(palette = "Greys") 
  theme(legend.title = element_blank)
###################################
# 3.13 PcoA visualization: diversity and clusters
###################################
# Input metadata
alpha_diversity_cluster

# Create PC axis
phyloseq_rarefied_phages_PcoA_vegan <- as.data.frame(phyloseq_rarefied_phages_PcoA$vectors[,1:2])
phyloseq_rarefied_phages_vegan <- merge(phyloseq_rarefied_phages_PcoA_vegan,alpha_diversity_cluster, by = 0,  all = F)
rownames(phyloseq_rarefied_phages_vegan) <- phyloseq_rarefied_phages_vegan$Row.names
phyloseq_rarefied_phages_vegan$Row.names <- NULL

## Modify richness
phyloseq_rarefied_phages_vegan$richness_categorical <- phyloseq_rarefied_phages_vegan$observed
median(phyloseq_rarefied_phages_vegan$observed) # 31
phyloseq_rarefied_phages_vegan$richness_categorical[phyloseq_rarefied_phages_vegan$observed >= 31] <- "High (> 31)"
phyloseq_rarefied_phages_vegan$richness_categorical[phyloseq_rarefied_phages_vegan$observed < 31] <- "Low (< 31)"
names(phyloseq_rarefied_phages_vegan)[names(phyloseq_rarefied_phages_vegan) == "observed"] <- "Observed richness"
names(phyloseq_rarefied_phages_vegan)[names(phyloseq_rarefied_phages_vegan) == "shannon"] <- "Shannon diversity"

# PCOA: richness
# Richness >= 30 is high & < 30 is low
 ggplot(phyloseq_rarefied_phages_vegan, aes(x = Axis.1, y= Axis.2, color = cluster, shape = richness_categorical)) +
  geom_point() +
  scale_shape_manual(values = c(16,1)) +
 #  tat_ellipse(type='t',size =1, linetype = 2,level = 0.95, aes(group = cluster)) +
  theme_bw() +
  scale_color_manual(values = c("#ae4c39","#39ae85")) +
  geom_hline(yintercept = 0, linetype = 3, colour = "black") + 
  geom_vline(xintercept = 0, linetype = 3, colour = "black")  +
  xlab("PC1 (8.6%)") +
  ylab("PC2 (4.1%)") 

## Adonis2 test: add to phyloseq & calculate significance
Plot_pcoa_richness_phylo <- as.data.frame(phyloseq_rarefied_phages_vegan[,c(3:4,6)])
Plot_pcoa_richness_phylo <- sample_data(Plot_pcoa_richness_phylo)
phyloseq_rarefied_phages_pcoa_1 <- merge_phyloseq(phyloseq_rarefied_phages_visualization_2_richness, Plot_pcoa_richness_phylo)

# Categorical richness
colnames(sample_data(phyloseq_rarefied_phages_pcoa_1))
adonis2(phyloseq_rarefied_phages_transformed_bray ~ sample_data(phyloseq_rarefied_phages_pcoa_1)$`Observed richness`)
# P-value (Pr(>F) = 0.001; Df= 1, R2 = 0.01928, F = 7.0953

# Observed Richness
adonis2(phyloseq_rarefied_phages_transformed_bray ~ sample_data(phyloseq_rarefied_phages_pcoa_1)$`Observed richness`)
# P-value (Pr(>F) = 0.001; Df= 1, R2 = 0.01358, F = 4.9696

# Shannon diversity
adonis2(phyloseq_rarefied_phages_transformed_bray ~ sample_data(phyloseq_rarefied_phages_pcoa_1)$`Shannon diversity`)
# P-value (Pr(>F) = 1; Df= 1, R2 = -0.0003, F = -0.1069

## Thus, we can associate both Observed as well ass categorical richness to cluster 1.
## Shannon diversity cannot be assigned to either of them. Non-significant results.
###################################
# 3.14 Statistical test
###################################
# If you want you could already do some exploratory analysis using the unconstrained distance-based PcoA plot.
# Here you can try to associate all metadata one-by-one to the distance metrix to assess significance & degree of association.
# This is purely exploratory and later one we will do more detailed analysis about this
# Remember this is not only done on PC1 and PC2, but on all variation.

# PERMANOVA (Permutation Based Analysis of Variance) statistcal test
# PERMANOVA tests if the centroids, similar to means, of each group are significantly different from each other. Likewise, an R2 statistic is calculated, showing the percentage of the variance explained by the groups.
###################################
# 3.15.1 Adonis
###################################
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)

#Tables needed
phyloseq_rarefied_phages_transformed #phyloseq_rarefied_phages_transformed # log10 transformed table
phyloseq_rarefied_phages_transformed_bray # Bray-curtis distance metrics
phyloseq_rarefied_phages_PcoA # Ordinationation plot

# R2 is percentage of variation that is explained & Pr > F is significance.
# Simple adonis2 test (<= 2 groups): categorical variables
colnames(sample_data(phyloseq_rarefied_phages_transformed))
adonis2(phyloseq_rarefied_phages_transformed_bray ~ sample_data(phyloseq_rarefied_phages_transformed)$Storage)
adonis2(phyloseq_rarefied_phages_transformed_bray ~ sample_data(phyloseq_rarefied_phages_transformed)$timepoints) # timepoints not significant 
adonis2(phyloseq_rarefied_phages_transformed_bray ~ sample_data(phyloseq_rarefied_phages_transformed)$Diagnosis)
adonis2(phyloseq_rarefied_phages_transformed_bray ~ sample_data(phyloseq_rarefied_phages_transformed)$Gender)
adonis2(phyloseq_rarefied_phages_transformed_bray ~ sample_data(phyloseq_rarefied_phages_transformed)$Smoking)
adonis2(phyloseq_rarefied_phages_transformed_bray ~ sample_data(phyloseq_rarefied_phages_transformed)$Biologicals)
adonis2(phyloseq_rarefied_phages_transformed_bray ~ sample_data(phyloseq_rarefied_phages_transformed)$Biologicals.classes)
adonis2(phyloseq_rarefied_phages_transformed_bray ~ sample_data(phyloseq_rarefied_phages_transformed)$Endoscopy.outcome)
adonis2(phyloseq_rarefied_phages_transformed_bray ~ sample_data(phyloseq_rarefied_phages_transformed))
View(sample_data(phyloseq_rarefied_phages_transformed))


#Simple adonis2 test (>= 2 groups): continous variables (can take long time)
adonis2(phyloseq_rarefied_phages_transformed_bray ~ sample_data(phyloseq_rarefied_phages_transformed)$Moisture) # not significant
adonis2(phyloseq_rarefied_phages_transformed_bray ~ sample_data(phyloseq_rarefied_phages_transformed)$Cellcount) # not significant
adonis2(phyloseq_rarefied_phages_transformed_bray ~ sample_data(phyloseq_rarefied_phages_transformed)$Calprotectin) # not significant
adonis2(phyloseq_rarefied_phages_transformed_bray ~ sample_data(phyloseq_rarefied_phages_transformed)$BMI)
adonis2(phyloseq_rarefied_phages_transformed_bray ~ sample_data(phyloseq_rarefied_phages_transformed)$Duration)
adonis2(phyloseq_rarefied_phages_transformed_bray ~ sample_data(phyloseq_rarefied_phages_transformed)$Hemoglobin) # not significant
adonis2(phyloseq_rarefied_phages_transformed_bray ~ sample_data(phyloseq_rarefied_phages_transformed)$Albumin) #not significant
adonis2(phyloseq_rarefied_phages_transformed_bray ~ sample_data(phyloseq_rarefied_phages_transformed)$CRP) # not significant

# Adonis2 test followed by pairwise comparisons (>= 3 groups)
unique(sample_data(phyloseq_rarefied_phages_transformed)$Therapy)
adonis2(phyloseq_rarefied_phages_transformed_bray ~ sample_data(phyloseq_rarefied_phages_transformed)$Therapy)
pairwise.adonis(phyloseq_rarefied_phages_transformed_bray,sample_data(phyloseq_rarefied_phages_transformed)$Therapy, p.adjust.m = "BH")

unique(sample_data(phyloseq_rarefied_phages_transformed)$Location)
adonis2(phyloseq_rarefied_phages_transformed_bray ~ sample_data(phyloseq_rarefied_phages_transformed)$Location)
pairwise.adonis(phyloseq_rarefied_phages_transformed_bray,sample_data(phyloseq_rarefied_phages_transformed)$Location, p.adjust.m = "BH")

## Conclusion:
## Only a few variables are not significant to virome composition being Weeks and timepoints.
## Paves way as usage of baseline for predictive value
## Perhaps decide to remove clusters from metadata associations because we added it, just like diversity.
####################################
# 4. DMM clusters: composition
####################################
# 4.1 RA: Class
####################################
# Create Dataframe
phyloseq_rarefied_phages_visualization_2_RA <- phyloseq_rarefied_phages_visualization_1
phyloseq_rarefied_phages_tax <- as.data.frame(tax_table(phyloseq_rarefied_phages_visualization_2_RA))
phyloseq_rarefied_phages_OTU <- as.data.frame(otu_table(phyloseq_rarefied_phages_visualization_2_RA))
phyloseq_rarefied_phages_metadata <- as.data.frame(sample_data(phyloseq_rarefied_phages_visualization_2_RA))

# Merge
phyloseq_rarefied_phages_ggplot <- merge(phyloseq_rarefied_phages_OTU,phyloseq_rarefied_phages_tax, by = 0, all = TRUE)
rownames(phyloseq_rarefied_phages_ggplot) <- phyloseq_rarefied_phages_ggplot$Row.names
phyloseq_rarefied_phages_ggplot$Row.names <- NULL

# Aggregate
vector_1 <- (which(names(phyloseq_rarefied_phages_ggplot)== "Kingdom")-1) 
ggplot_names <- colnames(phyloseq_rarefied_phages_ggplot[,c(1:vector_1)])

phyloseq_rarefied_phages_ggplot_class_agg <- aggregate(. ~Class, FUN = sum, data = phyloseq_rarefied_phages_ggplot[,colnames(phyloseq_rarefied_phages_ggplot) %in% ggplot_names | colnames(phyloseq_rarefied_phages_ggplot) == "Class"])
rownames(phyloseq_rarefied_phages_ggplot_class_agg) <- phyloseq_rarefied_phages_ggplot_class_agg$Class
phyloseq_rarefied_phages_ggplot_class_agg$Class <- NULL

# Calculate RA
phyloseq_rarefied_phages_ggplot_class <- sweep(phyloseq_rarefied_phages_ggplot_class_agg, 2, colSums(phyloseq_rarefied_phages_ggplot_class_agg), '/')
phyloseq_rarefied_phages_ggplot_class[is.na(phyloseq_rarefied_phages_ggplot_class)] <- 0

# Combine with metadata
phyloseq_rarefied_phages_ggplot_class <- merge(t(phyloseq_rarefied_phages_ggplot_class), phyloseq_rarefied_phages_metadata, by = 0, all = F)
rownames(phyloseq_rarefied_phages_ggplot_class) <- phyloseq_rarefied_phages_ggplot_class$Row.names
phyloseq_rarefied_phages_ggplot_class$Row.names <- NULL

# 7.6 Melt
vector_1 <- (which(names(phyloseq_rarefied_phages_ggplot_class)== "Unannotated")) 
vector_1_1 <- (which(names(phyloseq_rarefied_phages_ggplot_class)== "cluster")) 

phyloseq_rarefied_phages_ggplot_no_crass <- phyloseq_rarefied_phages_ggplot_class[,c(1:vector_1,vector_1_1)]
phyloseq_rarefied_phages_ggplot_class <- melt(phyloseq_rarefied_phages_ggplot_no_crass, id.vars = c("cluster"))

phyloseq_rarefied_phages_ggplot_class$cluster[is.na(phyloseq_rarefied_phages_ggplot_class$cluster)] <- "remove"
phyloseq_rarefied_phages_ggplot_class <- phyloseq_rarefied_phages_ggplot_class[!phyloseq_rarefied_phages_ggplot_class$cluster == "remove",]

# 7.7 Relative abundance boxplots
phyloseq_rarefied_phages_ggplot_class %>%
 filter(!is.na(cluster)) %>%
 filter(variable %in% c("Caudoviricetes (CrAss)", 
"Caudoviricetes (non-CrAss)", "Malgrandaviricetes")) %>%
 ggplot() +
 aes(x = cluster, y = value, fill = cluster) +
  scale_fill_manual(values = c(`1` = "#990011", `2` = "#00203f")) +
  geom_boxplot(shape = "circle") +
  geom_jitter(width = 0.15) +
 ylab("Relative abundance") +
 theme_bw() +
 facet_wrap(vars(variable), scales = "free_y")

# Statistics
# CrAss
phyloseq_rarefied_phages_3_stat_richness <- phyloseq_rarefied_phages_ggplot_class[phyloseq_rarefied_phages_ggplot_class$variable == "Caudoviricetes (CrAss)",]
wilcox.test(phyloseq_rarefied_phages_3_stat_richness$value~phyloseq_rarefied_phages_3_stat_richness$cluster, alternative = "two.sided",  paired = FALSE, conf.int = T, conf.level = 0.95, exact = FALSE, correct = FALSE)
p.adjust(1.386e-06, method = "BH", n = 3) #  P.adjust = 4.758e-13
nrow(phyloseq_rarefied_phages_3_stat_richness[phyloseq_rarefied_phages_3_stat_richness$cluster == 1,]) # 242 sample
nrow(phyloseq_rarefied_phages_3_stat_richness[phyloseq_rarefied_phages_3_stat_richness$cluster == 2,]) # 121 samples
effect_size <- wilcox_effsize(phyloseq_rarefied_phages_ggplot_class[phyloseq_rarefied_phages_ggplot_class$variable == "Caudoviricetes (CrAss)",],value~cluster)
(effect_size$effsize)

# non-CrAss
phyloseq_rarefied_phages_3_stat_richness <- phyloseq_rarefied_phages_ggplot_class[phyloseq_rarefied_phages_ggplot_class$variable == "Caudoviricetes (non-CrAss)",]
wilcox.test(phyloseq_rarefied_phages_3_stat_richness$value~phyloseq_rarefied_phages_3_stat_richness$cluster, alternative = "two.sided",  paired = FALSE, conf.int = T, conf.level = 0.95, exact = FALSE, correct = FALSE)
p.adjust(0.008882, method = "BH", n = 3) #  P.adjust = 4.758e-13
nrow(phyloseq_rarefied_phages_3_stat_richness[phyloseq_rarefied_phages_3_stat_richness$cluster == 1,]) # 242 sample
nrow(phyloseq_rarefied_phages_3_stat_richness[phyloseq_rarefied_phages_3_stat_richness$cluster == 2,]) # 121 samples
effect_size <- wilcox_effsize(phyloseq_rarefied_phages_ggplot_class[phyloseq_rarefied_phages_ggplot_class$variable == "Caudoviricetes (non-CrAss)",],value~cluster)
(effect_size$effsize)

# Malgrandaviricetes
phyloseq_rarefied_phages_3_stat_richness <- phyloseq_rarefied_phages_ggplot_class[phyloseq_rarefied_phages_ggplot_class$variable == "Malgrandaviricetes",]
wilcox.test(phyloseq_rarefied_phages_3_stat_richness$value~phyloseq_rarefied_phages_3_stat_richness$cluster, alternative = "two.sided",  paired = FALSE, conf.int = T, conf.level = 0.95, exact = FALSE, correct = FALSE)
p.adjust(0.0001078, method = "BH", n = 3) #  P.adjust = 0.0003234
nrow(phyloseq_rarefied_phages_3_stat_richness[phyloseq_rarefied_phages_3_stat_richness$cluster == 1,]) # 242 sample
nrow(phyloseq_rarefied_phages_3_stat_richness[phyloseq_rarefied_phages_3_stat_richness$cluster == 2,]) # 121 samples
effect_size <- wilcox_effsize(phyloseq_rarefied_phages_ggplot_class[phyloseq_rarefied_phages_ggplot_class$variable == "Caudoviricetes (non-CrAss)",],value~cluster)
(effect_size$effsize)

# > 0.01% abundance
median(phyloseq_rarefied_phages_ggplot_class$value[phyloseq_rarefied_phages_ggplot_class$variable == "Malgrandaviricetes"])
median(phyloseq_rarefied_phages_ggplot_class$value[phyloseq_rarefied_phages_ggplot_class$variable == "Caudoviricetes (non-CrAss)"])
median(phyloseq_rarefied_phages_ggplot_class$value[phyloseq_rarefied_phages_ggplot_class$variable == "Caudoviricetes (CrAss)"])

# < 0.01% abundance
median(phyloseq_rarefied_phages_ggplot_class$value[phyloseq_rarefied_phages_ggplot_class$variable == "Allassoviricetes"])
median(phyloseq_rarefied_phages_ggplot_class$value[phyloseq_rarefied_phages_ggplot_class$variable == "Faserviricetes"])
median(phyloseq_rarefied_phages_ggplot_class$value[phyloseq_rarefied_phages_ggplot_class$variable == "Unannotated"])
# check median > 0.01 or 1% otherwise leave
####################################
# 4.2 RA: genus
####################################
phyloseq_rarefied_phages_visualization_2_RA <- phyloseq_rarefied_phages_visualization_1

# Determine most abundant taxa
OTU10 = names(sort(taxa_sums(phyloseq_rarefied_phages_visualization_2_RA), TRUE)[1:10])
OTU10
taxa_extract <- c("genus_1", "genus_37", "genus_14", "genus_103", "genus_7", "genus_6","genus_5", "genus_95", "genus_11", "genus_4")

# Create Dataframe
phyloseq_rarefied_phages_tax <- as.data.frame(tax_table(phyloseq_rarefied_phages_visualization_2_RA))
phyloseq_rarefied_phages_OTU <- as.data.frame(otu_table(phyloseq_rarefied_phages_visualization_2_RA))
phyloseq_rarefied_phages_metadata <- as.data.frame(sample_data(phyloseq_rarefied_phages_visualization_2_RA))

# Merge
phyloseq_rarefied_phages_ggplot <- merge(phyloseq_rarefied_phages_OTU,phyloseq_rarefied_phages_tax, by = 0, all = TRUE)
rownames(phyloseq_rarefied_phages_ggplot) <- phyloseq_rarefied_phages_ggplot$Row.names
phyloseq_rarefied_phages_ggplot$Row.names <- NULL

# Aggregate
vector_1 <- (which(names(phyloseq_rarefied_phages_ggplot)== "Kingdom")-1) 
ggplot_names <- colnames(phyloseq_rarefied_phages_ggplot[,c(1:vector_1)])

phyloseq_rarefied_phages_ggplot_class_agg <- aggregate(. ~Genus, FUN = sum, data = phyloseq_rarefied_phages_ggplot[,colnames(phyloseq_rarefied_phages_ggplot) %in% ggplot_names | colnames(phyloseq_rarefied_phages_ggplot) == "Genus"])
rownames(phyloseq_rarefied_phages_ggplot_class_agg) <- phyloseq_rarefied_phages_ggplot_class_agg$Genus
phyloseq_rarefied_phages_ggplot_class_agg$Genus <- NULL

# Calculate RA
phyloseq_rarefied_phages_ggplot_class <- sweep(phyloseq_rarefied_phages_ggplot_class_agg, 2, colSums(phyloseq_rarefied_phages_ggplot_class_agg), '/')
phyloseq_rarefied_phages_ggplot_class[is.na(phyloseq_rarefied_phages_ggplot_class)] <- 0

# Extract most abundant taxa (~ as described above)
phyloseq_rarefied_phages_ggplot_class_tax <- as.data.frame(t(phyloseq_rarefied_phages_ggplot_class))
phyloseq_rarefied_phages_ggplot_class_tax <- phyloseq_rarefied_phages_ggplot_class_tax[, names(phyloseq_rarefied_phages_ggplot_class_tax) %in% taxa_extract]

# Remove empty samples
phyloseq_rarefied_phages_ggplot_class_tax <- as.data.frame(phyloseq_rarefied_phages_ggplot_class_tax)
phyloseq_rarefied_phages_ggplot_class_tax$total <- 0
phyloseq_rarefied_phages_ggplot_class_tax$total <- rowSums(phyloseq_rarefied_phages_ggplot_class_tax)
phyloseq_rarefied_phages_ggplot_class_tax <- phyloseq_rarefied_phages_ggplot_class_tax[!phyloseq_rarefied_phages_ggplot_class_tax$total == 0,]
phyloseq_rarefied_phages_ggplot_class_tax$total <- NULL

# Combine with metadata
phyloseq_rarefied_phages_ggplot_class <- merge(phyloseq_rarefied_phages_ggplot_class_tax, phyloseq_rarefied_phages_metadata, by = 0, all = F)
rownames(phyloseq_rarefied_phages_ggplot_class) <- phyloseq_rarefied_phages_ggplot_class$Row.names
phyloseq_rarefied_phages_ggplot_class$Row.names <- NULL

# 7.6 Melt
vector_1 <- (which(names(phyloseq_rarefied_phages_ggplot_class)== "genus_95")) 
vector_1_1 <- (which(names(phyloseq_rarefied_phages_ggplot_class)== "cluster")) 

phyloseq_rarefied_phages_ggplot_no_crass <- phyloseq_rarefied_phages_ggplot_class[,c(1:vector_1,vector_1_1)]
phyloseq_rarefied_phages_ggplot_class <- melt(phyloseq_rarefied_phages_ggplot_no_crass, id.vars = c("cluster"))

# 7.7 Relative abundance boxplots
phyloseq_rarefied_phages_ggplot_class %>%
  filter(!is.na(cluster)) %>%
  ggplot() +
  aes(x = cluster, y = value, fill = cluster) +
  scale_fill_manual(values = c(`1` = "#990011", `2` = "#00203f")) +
  geom_boxplot(shape = "circle") +
  ylab("Relative abundance") +
  theme_bw() +
  facet_wrap(vars(variable), scales = "free_y")

# 7.8 Statistics
# genus_1
phyloseq_rarefied_phages_ggplot_class_ <- phyloseq_rarefied_phages_ggplot_class[phyloseq_rarefied_phages_ggplot_class$variable == "genus_1",]
wilcox.test(value ~ cluster, data = phyloseq_rarefied_phages_ggplot_class_, alternative = "two.sided",  paired = FALSE, conf.int = T, conf.level = 0.95, exact = FALSE, correct = FALSE)

# genus_103
phyloseq_rarefied_phages_ggplot_class_ <- phyloseq_rarefied_phages_ggplot_class[phyloseq_rarefied_phages_ggplot_class$variable == "genus_103",]
wilcox.test(value ~ cluster, data = phyloseq_rarefied_phages_ggplot_class_, alternative = "two.sided",  paired = FALSE, conf.int = T, conf.level = 0.95, exact = FALSE, correct = FALSE)

# genus_11
phyloseq_rarefied_phages_ggplot_class_ <- phyloseq_rarefied_phages_ggplot_class[phyloseq_rarefied_phages_ggplot_class$variable == "genus_11",]
wilcox.test(value ~ cluster, data = phyloseq_rarefied_phages_ggplot_class_, alternative = "two.sided",  paired = FALSE, conf.int = T, conf.level = 0.95, exact = FALSE, correct = FALSE)

# genus_14
phyloseq_rarefied_phages_ggplot_class_ <- phyloseq_rarefied_phages_ggplot_class[phyloseq_rarefied_phages_ggplot_class$variable == "genus_14",]
wilcox.test(value ~ cluster, data = phyloseq_rarefied_phages_ggplot_class_, alternative = "two.sided",  paired = FALSE, conf.int = T, conf.level = 0.95, exact = FALSE, correct = FALSE)

# genus_37
phyloseq_rarefied_phages_ggplot_class_ <- phyloseq_rarefied_phages_ggplot_class[phyloseq_rarefied_phages_ggplot_class$variable == "genus_37",]
wilcox.test(value ~ cluster, data = phyloseq_rarefied_phages_ggplot_class_, alternative = "two.sided",  paired = FALSE, conf.int = T, conf.level = 0.95, exact = FALSE, correct = FALSE)

# genus_4
phyloseq_rarefied_phages_ggplot_class_ <- phyloseq_rarefied_phages_ggplot_class[phyloseq_rarefied_phages_ggplot_class$variable == "genus_4",]
wilcox.test(value ~ cluster, data = phyloseq_rarefied_phages_ggplot_class_, alternative = "two.sided",  paired = FALSE, conf.int = T, conf.level = 0.95, exact = FALSE, correct = FALSE)

# genus_5
phyloseq_rarefied_phages_ggplot_class_ <- phyloseq_rarefied_phages_ggplot_class[phyloseq_rarefied_phages_ggplot_class$variable == "genus_5",]
wilcox.test(value ~ cluster, data = phyloseq_rarefied_phages_ggplot_class_, alternative = "two.sided",  paired = FALSE, conf.int = T, conf.level = 0.95, exact = FALSE, correct = FALSE)

# genus_6, * not significant
phyloseq_rarefied_phages_ggplot_class_ <- phyloseq_rarefied_phages_ggplot_class[phyloseq_rarefied_phages_ggplot_class$variable == "genus_6",]
wilcox.test(value ~ cluster, data = phyloseq_rarefied_phages_ggplot_class_, alternative = "two.sided",  paired = FALSE, conf.int = T, conf.level = 0.95, exact = FALSE, correct = FALSE)

# genus_7
phyloseq_rarefied_phages_ggplot_class_ <- phyloseq_rarefied_phages_ggplot_class[phyloseq_rarefied_phages_ggplot_class$variable == "genus_7",]
wilcox.test(value ~ cluster, data = phyloseq_rarefied_phages_ggplot_class_, alternative = "two.sided",  paired = FALSE, conf.int = T, conf.level = 0.95, exact = FALSE, correct = FALSE)

# genus_95
phyloseq_rarefied_phages_ggplot_class_ <- phyloseq_rarefied_phages_ggplot_class[phyloseq_rarefied_phages_ggplot_class$variable == "genus_95",]
wilcox.test(value ~ cluster, data = phyloseq_rarefied_phages_ggplot_class_, alternative = "two.sided",  paired = FALSE, conf.int = T, conf.level = 0.95, exact = FALSE, correct = FALSE)
####################################
# 4.3 RA: genus & class summary
####################################
phyloseq_rarefied_phages_visualization_2_RA <- phyloseq_rarefied_phages_visualization_1
top20 <- names(sort(taxa_sums(phyloseq_rarefied_phages_visualization_2_RA), decreasing=TRUE)[1:10])
top20 #shows 20 results
dat.aglo = tax_glom(phyloseq_rarefied_phages_visualization_2_RA, taxrank = "Class")
dat.trans = transform_sample_counts(dat.aglo, function(x) x/sum(x))
prune.dat.two = prune_taxa(top20, dat.trans)
dat.dataframe = psmelt(prune.dat.two)
dat.agr = aggregate(Abundance~cluster+Class, data=dat.dataframe, FUN=mean)

ggplot(dat.agr, aes(x=cluster, y=Abundance, fill=Class)) + 
  geom_bar(stat="identity") +
  theme_classic()

phyloseq_rarefied_phages_visualization_2_RA <- phyloseq_rarefied_phages_visualization_1
top20 <- names(sort(taxa_sums(phyloseq_rarefied_phages_visualization_2_RA), decreasing=TRUE)[1:10])
top20 #shows 20 results
dat.aglo = tax_glom(phyloseq_rarefied_phages_visualization_2_RA, taxrank = "Genus")
dat.trans = transform_sample_counts(dat.aglo, function(x) x/sum(x))
prune.dat.two = prune_taxa(top20, dat.trans)
dat.dataframe = psmelt(prune.dat.two)
dat.agr = aggregate(Abundance~cluster+Genus, data=dat.dataframe, FUN=mean)

ggplot(dat.agr, aes(x=cluster, y=Abundance, fill=Genus)) + 
  geom_bar(stat="identity") +
  theme_classic()
####################################
# 4.4 Summary
####################################
# We can create a relative abundance figure on class level because we will aggregate all the reads per cluster (1 and 2), giving us a realistic picture.
# If you would do this on genus level you'll have to chose what to show, usually the most abundant genera (eg. top 10 or top 20 most abundant genera).
# Also samples are not shared equally well over samples, so you will simply gather samples with more reads & more biased viruses it (biased towards micro's eg).
# Upon that we cannot in any reliably way add a taxonomy to the genus level

# For that reason, we use genus level to explain cluster direction.
# Cluster 2 seems to have an increase of the class (I) Caudoviricetes (CrAss) and (II) Malgrandaviricetes
# Cluster 1 seems to have an increase of the class (I) Caudoviricetes (non-CrAss)
####################################
# 5 Metadata variables correlating to virome composition on all samples (baseline+PE+w14CD)
####################################
# 5.1 Input 
####################################
#Tables needed
phyloseq_rarefied_phages_transformed # log10 transformed table
phyloseq_rarefied_phages_transformed_bray # Bray-curtis distance metrics
phyloseq_rarefied_phages_PcoA # Ordinationation plot
####################################
# 5.2 Forward selection  & Distance-based redundancy analysis (db-RDA)
####################################
# Variable selection is a procedure for selecting a subset of explanatory variables from the set of all variables available for constrained ordination (RDA, CCA, db-RDA). The goal is to reduce the number of explanatory variables entering the analysis while keeping the variation explained by them to the maximum.
# The standard method is forward selection, which is adding explanatory variables one by one
# information: https://www.davidzeleny.net/anadat-r/doku.php/en:forward_sel

# Steps:
# 1) first, test the significance of the global test with all explanatory variables included; if it is significant, you may proceed to forward selection, while if it is not, it is better not to
# 2) use each variable one by one as explanatory in constrained ordination, and record the explained variation
# 3) sort variables according to variation explained by them with the highest values at the top
# 4) check whether the variation explained by the best variable is significant using Monte Carlo permutation test - if yes, include it to the model, if not, stop the selection
# 5) use each of remaining explanatory variables and check how much variation they (each separately) explain if put as explanatory (with the already selected variable acting as covariable)
# 6) sort again the variables according to the decreasing variation explained by them (now this variation represents the partial effect of this variable) and choose the one explaining the most; test whether the variation is significant, and if yes, select it into the model; if not, stop the selection;
# 7) continue by step 5 until the variation explained by the best variable is not significant

# Distance-based redundancy analysis (db-RDA) is a method for carrying out constrained ordinations on data using non-Euclidean distance measures.
# (1) A distance matrix is calculated using the distance measure of choice, (2) A principle coordinates analysis (PCoA) is done on the matrix, (3) The eigenvectors (PC's) obtained in the PCoA are plugged into an RDA.
# (2) In order to calculate multivariate cumulative R2 Metadata is selected using a foward manner using automatic stepwise model building for constrained ordinatio methods (cca, rda, capscale). 
# (3) In order to calculate univariate R2 you add only the variable you want to the model.

#More information: https://archetypalecology.wordpress.com/2018/02/21/distance-based-redundancy-analysis-db-rda-in-r/
####################################
# 5.2.1 Multivariate db-RDA
####################################
# Convert missing values or "NAs" to median value to solve the NA problem for na
# ordiR2step will be upgraded at some point and then it should be possible to deal with NA's better.
# Montral classification: name UC L2 for colonic affect, then you also have CD L1 (ileal presentation) and L3 (ieo-colonic presentations)
Metadata <- data.frame(sample_data(phyloseq_rarefied_phages_transformed))

# univariate first (see next header) --> variables?

# RDA data
# Do NA exclude on whole dataset (also the bray-curtis metric to keep data matched)
dbRDA_0 <- capscale(phyloseq_rarefied_phages_transformed_bray ~ 1, data = Metadata)  # scale=TRUE
#dbRDA_1 <- capscale(phyloseq_rarefied_phages_transformed_bray ~ patient_ID + Location + Age + BMI + Smoking + Endoscopy.outcome + Diagnosis + Hemoglobin + Moisture + Biologicals.classes + Albumin + Biologicals + Cellcount + CRP + Gender, data = Metadata) # All timepoints
#dbRDA_1 <- capscale(phyloseq_rarefied_phages_transformed_bray ~ Diagnosis + Age + BMI + Hemoglobin + Smoking, data = Metadata) # BASELINE
dbRDA_1 <- capscale(phyloseq_rarefied_phages_transformed_bray ~ Location + Diagnosis + Age + BMI + Endoscopy.outcome, data = Metadata) # PE
dbRDA_0
dbRDA_1

# Ordistep (R2) 
# This is a package in Vegan for automatic stepwise model building for constrained ordination methods (eg. CCA, RDA, capscale)
# Take into account that in forward selection you need to add all variables successively 
# If you just give all metadata as such you actually perform backward selection.
step.res_cum_ALL <- ordiR2step(dbRDA_0, scope = formula(dbRDA_1), direction="forward", R2scope = TRUE, permutations = 999)
RsquareAdj(step.res_cum_ALL)$r.squared
RsquareAdj(step.res_cum_ALL)$adj.r.squared
(step.res_cum_ALL)$anova

# Multivariate model: P-value and Adjusted R2 of individual covariates
step.res_cum_ALL$anova  # Summary table 
anova(step.res_cum_ALL)
anova(step.res_cum_ALL, by="terms", permu=200) # test for sign. environ. variables
anova(step.res_cum_ALL, by="axis", perm.max=500) # test axes for significance

# Multivariate model: P-value & R2
anova(step.res_cum_ALL)
model_anova <- anova(step.res_cum_ALL)$'Pr(>F)'[1] # P-value of 0.001
model_anova

(RsquareAdj(step.res_cum_ALL)$adj.r.squared)*100 # models explain 63.7% of microbiome variation
nrow(Metadata) # samples size of 363 samples

table(Metadata$Endoscopy.outcome)
101+68+104
# Multivariate model explains 78.1% of microbiome variation with P-value of 0.001 with 363 samples

# Both timepoints (adapted therapy as 3 categories): Individuality, Location + Therapy (CC) + Moisture
# Baseline (active IBD): BMI + Hemoglobin + Age
# PE: Location + Endoscopic outcome
####################################
# 5.2.2 Univariate db-RDA
####################################
# Couldnt find easy package for anova, so I simply redo ordiR2step but only on one variable
# This provides all data I need for statistics.
# I have 22 metadata variable, so I use 22 for fdr calculations ofr univarate analysis

# Patient ID*
dbRDA_0 <- capscale(phyloseq_rarefied_phages_transformed_bray ~ 1, data = Metadata, na.action = na.exclude)
dbRDA_1 <- capscale(phyloseq_rarefied_phages_transformed_bray ~  patient_ID, data = Metadata, na.action = na.exclude)
step.res_cum_patient_ID <- ordiR2step(dbRDA_0, dbRDA_1, direction="forward", R2scope = TRUE, na.action = na.exclude)
RsquareAdj(step.res_cum_patient_ID)$r.squared
RsquareAdj(step.res_cum_patient_ID)$adj.r.squared
step.res_cum_patient_ID$anova
step.res_cum_patient_ID <- step.res_cum_patient_ID$anova[1,1:5]

# Location*
dbRDA_0 <- capscale(phyloseq_rarefied_phages_transformed_bray ~ 1, data = Metadata, na.action = na.exclude)
dbRDA_1 <- capscale(phyloseq_rarefied_phages_transformed_bray ~ Location, data = Metadata, na.action = na.exclude)
step.res_cum_location <- ordiR2step(dbRDA_0, dbRDA_1, direction="forward", R2scope = TRUE, na.action = na.exclude)
RsquareAdj(step.res_cum_location)$r.squared
RsquareAdj(step.res_cum_location)$adj.r.squared
step.res_cum_location <- step.res_cum_location$anova[1,1:5]

# endoscopical outcome*
dbRDA_0 <- capscale(phyloseq_rarefied_phages_transformed_bray ~ 1, data = Metadata, na.action = na.exclude)
dbRDA_1 <- capscale(phyloseq_rarefied_phages_transformed_bray ~ Endoscopy.outcome, data = Metadata, na.action = na.exclude)
step.res_cum_endo <- ordiR2step(dbRDA_0, dbRDA_1, direction="forward", R2scope = TRUE, na.action = na.exclude)
RsquareAdj(step.res_cum_endo)$r.squared
RsquareAdj(step.res_cum_endo)$adj.r.squared
step.res_cum_endo <- step.res_cum_endo$anova[1,1:5]

# Moisture*
dbRDA_0 <- capscale(phyloseq_rarefied_phages_transformed_bray ~ 1, data = Metadata, na.action = na.exclude)
dbRDA_1 <- capscale(phyloseq_rarefied_phages_transformed_bray ~ Moisture, data = Metadata, na.action = na.exclude)
step.res_cum_mois <- ordiR2step(dbRDA_0, dbRDA_1, direction="forward", R2scope = TRUE, na.action = na.exclude)
RsquareAdj(step.res_cum_mois)$r.squared
RsquareAdj(step.res_cum_mois)$adj.r.squared
step.res_cum_mois <- step.res_cum_mois$anova[1,1:5]

# Other variables
# Storage
dbRDA_0 <- capscale(phyloseq_rarefied_phages_transformed_bray ~ 1, data = Metadata, na.action = na.exclude)
dbRDA_1 <- capscale(phyloseq_rarefied_phages_transformed_bray ~ Storage, data = Metadata, na.action = na.exclude)
step.res_cum_Storage <- ordiR2step(dbRDA_0, dbRDA_1, direction="forward", R2scope = TRUE, na.action = na.exclude)
RsquareAdj(step.res_cum_Storage)$r.squared
RsquareAdj(step.res_cum_Storage)$adj.r.squared
step.res_cum_Storage <- step.res_cum_Storage$anova[1,1:5]

# Diagnosis*
dbRDA_0 <- capscale(phyloseq_rarefied_phages_transformed_bray ~ 1, data = Metadata, na.action = na.exclude)
dbRDA_1 <- capscale(phyloseq_rarefied_phages_transformed_bray ~ Diagnosis, data = Metadata, na.action = na.exclude)
step.res_cum_Diagnosis <- ordiR2step(dbRDA_0, dbRDA_1, direction="forward", R2scope = TRUE, na.action = na.exclude)
RsquareAdj(step.res_cum_Diagnosis)$r.squared
RsquareAdj(step.res_cum_Diagnosis)$adj.r.squared
step.res_cum_Diagnosis <- step.res_cum_Diagnosis$anova[1,1:5]
 
# Age*
dbRDA_0 <- capscale(phyloseq_rarefied_phages_transformed_bray ~ 1, data = Metadata, na.action = na.exclude)
dbRDA_1 <- capscale(phyloseq_rarefied_phages_transformed_bray ~ Age, data = Metadata, na.action = na.exclude)
step.res_cum_Age <- ordiR2step(dbRDA_0, dbRDA_1, direction="forward", R2scope = TRUE, na.action = na.exclude)
RsquareAdj(step.res_cum_Age)$r.squared
RsquareAdj(step.res_cum_Age)$adj.r.squared
step.res_cum_Age <- step.res_cum_Age$anova[1,1:5]

# Gender
dbRDA_0 <- capscale(phyloseq_rarefied_phages_transformed_bray ~ 1, data = Metadata, na.action = na.exclude)
dbRDA_1 <- capscale(phyloseq_rarefied_phages_transformed_bray ~ Gender, data = Metadata, na.action = na.exclude)
step.res_cum_Gender <- ordiR2step(dbRDA_0, dbRDA_1, direction="forward", R2scope = TRUE, na.action = na.exclude)
RsquareAdj(step.res_cum_Gender)$r.squared
RsquareAdj(step.res_cum_Gender)$adj.r.squared
step.res_cum_Gender <- step.res_cum_Gender$anova[1,1:5]

# Therapy*
dbRDA_0 <- capscale(phyloseq_rarefied_phages_transformed_bray ~ 1, data = Metadata, na.action = na.exclude)
dbRDA_1 <- capscale(phyloseq_rarefied_phages_transformed_bray ~ Therapy, data = Metadata, na.action = na.exclude)
step.res_cum_Therapy <- ordiR2step(dbRDA_0, dbRDA_1, direction="forward", R2scope = TRUE, na.action = na.exclude)
RsquareAdj(step.res_cum_Therapy)$r.squared
RsquareAdj(step.res_cum_Therapy)$adj.r.squared
step.res_cum_Therapy <- step.res_cum_Therapy$anova[1,1:5]

# Cellcount
dbRDA_0 <- capscale(phyloseq_rarefied_phages_transformed_bray ~ 1, data = Metadata, na.action = na.exclude)
dbRDA_1 <- capscale(phyloseq_rarefied_phages_transformed_bray ~ Cellcount, data = Metadata, na.action = na.exclude)
step.res_cum_Cellcount <- ordiR2step(dbRDA_0, dbRDA_1, direction="forward", R2scope = TRUE, na.action = na.exclude)
RsquareAdj(step.res_cum_Cellcount)$r.squared
RsquareAdj(step.res_cum_Cellcount)$adj.r.squared
step.res_cum_Cellcount <- step.res_cum_Cellcount$anova[1,1:5]

# Calprotectin
dbRDA_0 <- capscale(phyloseq_rarefied_phages_transformed_bray ~ 1, data = Metadata, na.action = na.exclude)
dbRDA_1 <- capscale(phyloseq_rarefied_phages_transformed_bray ~ Calprotectin, data = Metadata, na.action = na.exclude)
step.res_cum_Calprotectin <- ordiR2step(dbRDA_0, dbRDA_1, direction="forward", R2scope = TRUE, na.action = na.exclude)
RsquareAdj(step.res_cum_Calprotectin)$r.squared
RsquareAdj(step.res_cum_Calprotectin)$adj.r.squared
step.res_cum_Calprotectin <- step.res_cum_Calprotectin$anova[1,1:5]

# BMI*
dbRDA_0 <- capscale(phyloseq_rarefied_phages_transformed_bray ~ 1, data = Metadata, na.action = na.exclude)
dbRDA_1 <- capscale(phyloseq_rarefied_phages_transformed_bray ~ BMI, data = Metadata, na.action = na.exclude)
step.res_cum_BMI <- ordiR2step(dbRDA_0, dbRDA_1, direction="forward", R2scope = TRUE, na.action = na.exclude)
RsquareAdj(step.res_cum_BMI)$r.squared
RsquareAdj(step.res_cum_BMI)$adj.r.squared
step.res_cum_BMI <- step.res_cum_BMI$anova[1,1:5]

# Duration
dbRDA_0 <- capscale(phyloseq_rarefied_phages_transformed_bray ~ 1, data = Metadata, na.action = na.exclude)
dbRDA_1 <- capscale(phyloseq_rarefied_phages_transformed_bray ~ Duration, data = Metadata, na.action = na.exclude)
step.res_cum_Duration <- ordiR2step(dbRDA_0, dbRDA_1, direction="forward", R2scope = TRUE, na.action = na.exclude)
RsquareAdj(step.res_cum_Duration)$r.squared
RsquareAdj(step.res_cum_Duration)$adj.r.squared
step.res_cum_Duration <- step.res_cum_Duration$anova[1,1:5]

# Smoking
dbRDA_0 <- capscale(phyloseq_rarefied_phages_transformed_bray ~ 1, data = Metadata, na.action = na.exclude)
dbRDA_1 <- capscale(phyloseq_rarefied_phages_transformed_bray ~ Smoking, data = Metadata, na.action = na.exclude)
step.res_cum_Smoking <- ordiR2step(dbRDA_0, dbRDA_1, direction="forward", R2scope = TRUE, na.action = na.exclude)
RsquareAdj(step.res_cum_Smoking)$r.squared
RsquareAdj(step.res_cum_Smoking)$adj.r.squared
step.res_cum_Smoking <- step.res_cum_Smoking$anova[1,1:5]

# Hemoglobin
dbRDA_0 <- capscale(phyloseq_rarefied_phages_transformed_bray ~ 1, data = Metadata, na.action = na.exclude)
dbRDA_1 <- capscale(phyloseq_rarefied_phages_transformed_bray ~ Hemoglobin, data = Metadata, na.action = na.exclude)
step.res_cum_Hemoglobin <- ordiR2step(dbRDA_0, dbRDA_1, direction="forward", R2scope = TRUE, na.action = na.exclude)
RsquareAdj(step.res_cum_Hemoglobin)$r.squared
RsquareAdj(step.res_cum_Hemoglobin)$adj.r.squared
step.res_cum_Hemoglobin <- step.res_cum_Hemoglobin$anova[1,1:5]
 
# Albumin
dbRDA_0 <- capscale(phyloseq_rarefied_phages_transformed_bray ~ 1, data = Metadata, na.action = na.exclude)
dbRDA_1 <- capscale(phyloseq_rarefied_phages_transformed_bray ~ Albumin, data = Metadata, na.action = na.exclude)
step.res_cum_Albumin <- ordiR2step(dbRDA_0, dbRDA_1, direction="forward", R2scope = TRUE, na.action = na.exclude)
RsquareAdj(step.res_cum_Albumin)$r.squared
RsquareAdj(step.res_cum_Albumin)$adj.r.squared
step.res_cum_Albumin <- step.res_cum_Albumin$anova[1,1:5]

# CRP
dbRDA_0 <- capscale(phyloseq_rarefied_phages_transformed_bray ~ 1, data = Metadata, na.action = na.exclude)
dbRDA_1 <- capscale(phyloseq_rarefied_phages_transformed_bray ~ CRP, data = Metadata, na.action = na.exclude)
step.res_cum_CRP <- ordiR2step(dbRDA_0, dbRDA_1, direction="forward", R2scope = TRUE, na.action = na.exclude)
RsquareAdj(step.res_cum_CRP)$r.squared
RsquareAdj(step.res_cum_CRP)$adj.r.squared
step.res_cum_CRP <- step.res_cum_CRP$anova[1,1:5]

# Biologicals
dbRDA_0 <- capscale(phyloseq_rarefied_phages_transformed_bray ~ 1, data = Metadata, na.action = na.exclude)
dbRDA_1 <- capscale(phyloseq_rarefied_phages_transformed_bray ~ Biologicals, data = Metadata, na.action = na.exclude)
step.res_cum_Biologicals <- ordiR2step(dbRDA_0, dbRDA_1, direction="forward", R2scope = TRUE, na.action = na.exclude)
RsquareAdj(step.res_cum_Biologicals)$r.squared
RsquareAdj(step.res_cum_Biologicals)$adj.r.squared
step.res_cum_Biologicals <- step.res_cum_Biologicals$anova[1,1:5]

# Biologicals.classes
dbRDA_0 <- capscale(phyloseq_rarefied_phages_transformed_bray ~ 1, data = Metadata, na.action = na.exclude)
dbRDA_1 <- capscale(phyloseq_rarefied_phages_transformed_bray ~ Biologicals.classes, data = Metadata, na.action = na.exclude)
step.res_cum_Biologicals.classes <- ordiR2step(dbRDA_0, dbRDA_1, direction="forward", R2scope = TRUE, na.action = na.exclude)
RsquareAdj(step.res_cum_Biologicals.classes)$r.squared
RsquareAdj(step.res_cum_Biologicals.classes)$adj.r.squared
step.res_cum_Biologicals.classes <- step.res_cum_Biologicals.classes$anova[1,1:5]

# Clinical.outcome
dbRDA_0 <- capscale(phyloseq_rarefied_phages_transformed_bray ~ 1, data = Metadata, na.action = na.exclude)
dbRDA_1 <- capscale(phyloseq_rarefied_phages_transformed_bray ~ Clinical.outcome, data = Metadata, na.action = na.exclude)
step.res_cum_clinical <- ordiR2step(dbRDA_0, dbRDA_1, direction="forward", R2scope = TRUE, na.action = na.exclude)
RsquareAdj(step.res_cum_clinical)$r.squared
RsquareAdj(step.res_cum_clinical)$adj.r.squared
step.res_cum_clinical <- step.res_cum_clinical$anova[1,1:5]

# Biomarker.outcome
dbRDA_0 <- capscale(phyloseq_rarefied_phages_transformed_bray ~ 1, data = Metadata, na.action = na.exclude)
dbRDA_1 <- capscale(phyloseq_rarefied_phages_transformed_bray ~ Biomarker.outcome, data = Metadata, na.action = na.exclude)
step.res_cum_Biomarker <- ordiR2step(dbRDA_0, dbRDA_1, direction="forward", R2scope = TRUE, na.action = na.exclude)
RsquareAdj(step.res_cum_Biomarker)$r.squared
RsquareAdj(step.res_cum_Biomarker)$adj.r.squared
step.res_cum_Biomarker <- step.res_cum_Biomarker$anova[1,1:5]

# Combined outcome
dbRDA_0 <- capscale(phyloseq_rarefied_phages_transformed_bray ~ 1, data = Metadata, na.action = na.exclude)
dbRDA_1 <- capscale(phyloseq_rarefied_phages_transformed_bray ~ Combined_remission, data = Metadata, na.action = na.exclude)
step.res_cum_Biomarker <- ordiR2step(dbRDA_0, dbRDA_1, direction="forward", R2scope = TRUE, na.action = na.exclude)
RsquareAdj(step.res_cum_Biomarker)$r.squared
RsquareAdj(step.res_cum_Biomarker)$adj.r.squared
step.res_cum_Biomarker <- step.res_cum_Biomarker$anova[1,1:5]

# timepoints
dbRDA_0 <- capscale(phyloseq_rarefied_phages_transformed_bray ~ 1, data = Metadata, na.action = na.exclude)
dbRDA_1 <- capscale(phyloseq_rarefied_phages_transformed_bray ~ timepoints, data = Metadata, na.action = na.exclude)
step.res_cum_timepoints <- ordiR2step(dbRDA_0, dbRDA_1, direction="forward", R2scope = TRUE, na.action = na.exclude)
RsquareAdj(step.res_cum_timepoints)$r.squared
RsquareAdj(step.res_cum_timepoints)$adj.r.squared
step.res_cum_timepoints <- step.res_cum_timepoints$anova[1,1:5]

# Univariate modelling: combine all data
step.res_cum_patient_ID
step.res_cum_location
step.res_cum_Age
step.res_cum_mois
step.res_cum_endo
step.res_cum_Storage
step.res_cum_Diagnosis
step.res_cum_Gender
step.res_cum_Cellcount
step.res_cum_BMI
step.res_cum_Therapy
step.res_cum_Duration
step.res_cum_Smoking
step.res_cum_Hemoglobin
step.res_cum_Albumin
step.res_cum_CRP
step.res_cum_Biologicals
step.res_cum_Biologicals.classes

# Merge DF: all IBD
DF_1 <- bind_rows(step.res_cum_patient_ID,anti_join(step.res_cum_location,step.res_cum_patient_ID,by="R2.adj"))
DF_2 <- bind_rows(DF_1,anti_join(step.res_cum_endo,DF_1,by="R2.adj"))
DF_3 <- bind_rows(DF_2,anti_join(step.res_cum_mois,DF_2,by="R2.adj"))
DF_4 <- bind_rows(DF_3,anti_join(step.res_cum_Storage,DF_3,by="R2.adj"))
DF_5 <- bind_rows(DF_4,anti_join(step.res_cum_Diagnosis,DF_4,by="R2.adj"))
DF_6 <- bind_rows(DF_5,anti_join(step.res_cum_Age,DF_5,by="R2.adj"))
DF_7 <- bind_rows(DF_6,anti_join(step.res_cum_Gender,DF_6,by="R2.adj"))
DF_8 <- bind_rows(DF_7,anti_join(step.res_cum_Cellcount,DF_7,by="R2.adj"))
DF_9 <- bind_rows(DF_8,anti_join(step.res_cum_BMI,DF_8,by="R2.adj"))
DF_10 <- bind_rows(DF_9,anti_join(step.res_cum_Smoking,DF_9,by="R2.adj"))
DF_11 <- bind_rows(DF_10,anti_join(step.res_cum_Hemoglobin,DF_10,by="R2.adj"))
DF_12 <- bind_rows(DF_11,anti_join(step.res_cum_Albumin,DF_11,by="R2.adj"))
DF_13 <- bind_rows(DF_12,anti_join(step.res_cum_CRP,DF_12,by="R2.adj"))
DF_14 <- bind_rows(DF_13,anti_join(step.res_cum_Biologicals,DF_13,by="R2.adj"))
DF_15 <- bind_rows(DF_14,anti_join(step.res_cum_Biologicals.classes,DF_14,by="R2.adj"))
View(DF_15)

# Baseline
#DF_1 <- bind_rows(step.res_cum_Age,anti_join(step.res_cum_Diagnosis,step.res_cum_Age,by="R2.adj"))
#DF_2 <- bind_rows(DF_1,anti_join(step.res_cum_BMI,DF_1,by="R2.adj"))
#DF_3 <- bind_rows(DF_2,anti_join(step.res_cum_Smoking,DF_2,by="R2.adj"))
#DF_4 <- bind_rows(DF_3,anti_join(step.res_cum_Hemoglobin,DF_3,by="R2.adj"))
#View(DF_4)

# PE
#DF_1 <- bind_rows(step.res_cum_location,anti_join(step.res_cum_endo,step.res_cum_location,by="R2.adj"))
#DF_2 <- bind_rows(DF_1,anti_join(step.res_cum_Diagnosis,DF_1,by="R2.adj"))
#DF_3 <- bind_rows(DF_2,anti_join(step.res_cum_Age,DF_2,by="R2.adj"))
#DF_4 <- bind_rows(DF_3,anti_join(step.res_cum_endo,DF_3,by="R2.adj"))
#DF_5 <- bind_rows(DF_4,anti_join(step.res_cum_BMI,DF_3,by="R2.adj"))
#View(DF_5)

# Calling df, F-value, p_value, P_adjusted values, R2.adj, AIC
DF_all <- DF_15
DF_all$p_value <- DF_all$`Pr(>F)`
DF_all$`Pr(>F)` <- NULL
n.tests <- nrow (DF_all)  # number of tests equals to number of all variables from which is being selected
pval.adj <- p.adjust(DF_all$p_value, method = 'BH', n = n.tests)
DF_all$pval.adj <- pval.adj
View(DF_all)

# Stacked Barplot
DF_1 <- as.data.frame(c("Patient ID", "Location", "Age", "Moisture", "Unexplained"))
DF_1$names <- DF_1$`c("Patient ID", "Location", "Age", "Moisture", "Unexplained")`
DF_1$`c("Patient ID", "Location", "Age", "Moisture", "Unexplained")`<- NULL
DF_1$percentage <- c(75.8, 1.40, 0.5, 0.3, 22.0)
DF_1$general <- "general"
DF_1$general[DF_1$general == "general" & DF_1$names == "Unexplained"] <- "unexplained"
DF_1$names <- factor(DF_1$names,levels=c("Moisture","Age","Location","Patient ID", "Unexplained"))

ggplot(DF_1) +
 aes(x = general, fill = names, weight = percentage) +
 geom_bar() +
 scale_fill_hue(direction = 1) +
 coord_flip() +
  ylab("Multivariate effect size (%)") +
  theme_bw() +
  ylab("") +
  xlab("")


# ALl data in here for table later!
####################################
# 5.2.3 Visualization effect size variables: donut plot
####################################
# Output Multivariate data
step.res_cum$anova*100

# Output Univariate data
step.res_cum_patient_ID$anova*100
step.res_cum_location$anova*100
step.res_cum_endo$anova*100
step.res_cum_mois$anova*100
step.res_cum_Storage$anova*100
step.res_cum_Diagnosis$anova*100
step.res_cum_Age$anova*100
step.res_cum_Gender$anova*100
step.res_cum_Therapy$anova*100
step.res_cum_Cellcount$anova*100
step.res_cum_Calprotectin$anova # not significant
step.res_cum_BMI$anova*100
step.res_cum_Duration$anova*100
step.res_cum_Smoking$anova*100
step.res_cum_Hemoglobin$anova*100
step.res_cum_Albumin$anova*100
step.res_cum_CRP$anova*100
step.res_cum_Biologicals$anova*100
step.res_cum_Biologicals.classes$anova*100
step.res_cum_clinical$anova*100 # not significant
step.res_cum_Biomarker$anova*100
step.res_cum_timepoints$anova # not significant
colnames(Metadata)

## Effect size figure_2:IBD
DF_variables <- as.data.frame(c("Patient_ID","Location","Age","Moisture","BMI","Smoking","Endoscopic outcome","Diagnosis","Hemoglobin","Biologicals classes","Albumin","Biologicals","Cellcount","CRP","Gender","Storage","Patient_ID","Location","Age","Moisture","BMI","Smoking","Endoscopic outcome","Diagnosis","Hemoglobin","Biologicals classes","Albumin","Biologicals","Cellcount","CRP","Gender","Storage"))
DF_variables$variables <- DF_variables$`c("Patient_ID", "Location", "Age", "Moisture", "BMI", "Smoking", "Endoscopic outcome", "Diagnosis", "Hemoglobin", "Biologicals classes", "Albumin", "Biologicals", "Cellcount", "CRP", "Gender", "Storage", "Patient_ID", "Location", "Age", "Moisture", "BMI", "Smoking", "Endoscopic outcome", "Diagnosis", "Hemoglobin", "Biologicals classes", "Albumin", "Biologicals", "Cellcount", "CRP", "Gender", "Storage")`
DF_variables$`c("Patient_ID", "Location", "Age", "Moisture", "BMI", "Smoking", "Endoscopic outcome", "Diagnosis", "Hemoglobin", "Biologicals classes", "Albumin", "Biologicals", "Cellcount", "CRP", "Gender", "Storage", "Patient_ID", "Location", "Age", "Moisture", "BMI", "Smoking", "Endoscopic outcome", "Diagnosis", "Hemoglobin", "Biologicals classes", "Albumin", "Biologicals", "Cellcount", "CRP", "Gender", "Storage")` <- NULL
DF_variables$variable_type <- c("multivariate","multivariate","multivariate","multivariate","multivariate","multivariate","multivariate","multivariate","multivariate","multivariate","multivariate","multivariate","multivariate","multivariate","multivariate","multivariate","univariate","univariate","univariate","univariate","univariate","univariate","univariate","univariate","univariate","univariate","univariate","univariate","univariate","univariate","univariate","univariate")
DF_variables$values <- c("75.8","77.2","77.7","78.0","78.0","78.0","78.0","78.0","78.0","78.0","78.0","78.0","78.0","78.0","78.0","78.0","75.8","1.34","0.812","0.390","0.773","0.596","0.506","0.473","0.421","0.292","0.270","0.237","0.211","0.209","0.207","0.172")
ord <- c("Patient_ID","Location","Age","Moisture","BMI","Smoking","Endoscopic outcome","Diagnosis","Hemoglobin","Biologicals classes","Albumin","Biologicals","Cellcount","CRP","Gender","Storage")
DF_variables$variables <- factor(DF_variables$variables, levels=rev(ord))
DF_variables$values <- as.numeric(DF_variables$values)

ggplot(DF_variables) +
  aes(x = variables, fill = variable_type, weight = values) +
  geom_bar(position = "dodge") +
  scale_fill_manual(values = c(multivariate = "#8B8A8B", univariate = "#000000")) +
  coord_flip() +
  theme_bw() +
  labs(x= "covariates of virome composition", y = "effect size (%)") +
  theme(legend.title = element_blank()) +
  geom_hline(yintercept = 78.0, linetype = 5, colour = "#49006a")

## Effect size figure_2: baseline
DF_variables <- as.data.frame(c("BMI", "Hemoglobin", "Age","Smoking","Diagnosis", "BMI", "Hemoglobin", "Age","Smoking","Diagnosis"))
DF_variables$variables <- DF_variables$`c("BMI", "Hemoglobin", "Age", "Smoking", "Diagnosis", "BMI", "Hemoglobin", "Age", "Smoking", "Diagnosis")`
DF_variables$`c("BMI", "Hemoglobin", "Age", "Smoking", "Diagnosis", "BMI", "Hemoglobin", "Age", "Smoking", "Diagnosis")` <- NULL
DF_variables$variable_type <- c("multivariate","multivariate","multivariate","multivariate","multivariate","univariate","univariate","univariate","univariate","univariate")
DF_variables$values <- c("0.725","1.21", "1.62", "1.94", "1.94","0.725","0.717","0.553","0.322","0.360")
ord <- c("BMI", "Hemoglobin", "Age","Smoking", "Diagnosis")
DF_variables$variables <- factor(DF_variables$variables, levels=rev(ord))
DF_variables$values <- as.numeric(DF_variables$values)

ggplot(DF_variables) +
  aes(x = variables, fill = variable_type, weight = values) +
  geom_bar(position = "dodge") +
  scale_fill_manual(values = c(multivariate = "#8B8A8B", univariate = "#000000")) +
  coord_flip() +
  theme_bw() +
  labs(x= "covariates of virome composition", y = "effect size (%)") +
  theme(legend.title = element_blank()) +
  geom_hline(yintercept = 1.94, linetype = 5, colour = "#49006a")

## Effect size figure_2: PE
DF_variables <- as.data.frame(c("Location","Endoscopic outcome", "Age", "Diagnosis", "BMI","Location","Endoscopic outcome","Age","Diagnosis", "BMI"))
DF_variables$variables <- DF_variables$`c("Location", "Endoscopic outcome", "Age", "Diagnosis", "BMI", "Location", "Endoscopic outcome", "Age", "Diagnosis", "BMI")`
DF_variables$`c("Location", "Endoscopic outcome", "Age", "Diagnosis", "BMI", "Location", "Endoscopic outcome", "Age", "Diagnosis", "BMI")` <- NULL
DF_variables$variable_type <- c("multivariate","multivariate","multivariate","multivariate","multivariate","univariate","univariate","univariate","univariate","univariate")
DF_variables$values <- c("1.01","1.47","1.84","1.84", "1.84", "1.01", "0.707", "0.427","0.474","0.346")
ord <- c("Location","Endoscopic outcome", "Age","Diagnosis","BMI")
DF_variables$variables <- factor(DF_variables$variables, levels=rev(ord))
DF_variables$values <- as.numeric(DF_variables$values)

ggplot(DF_variables) +
  aes(x = variables, fill = variable_type, weight = values) +
  geom_bar(position = "dodge") +
  scale_fill_manual(values = c(multivariate = "#8B8A8B", univariate = "#000000")) +
  coord_flip() +
  theme_bw() +
  labs(x= "covariates of virome composition", y = "effect size (%)") +
  theme(legend.title = element_blank()) +
  geom_hline(yintercept = 1.84, linetype = 5, colour = "#49006a")

# Multivariate model for all timepoints
# Total Multivariate model explains 78.1% of microbiome variation with P-value of 0.001 with 363 samples
# Significant variables in multivariate model, but here univariate contribution is given.
# (1) patient_ID = 75.8%; p-value = 0.002, n = 363
# (2) Location = 1.40%; p-value = 0.002, n = 363
# (3) Endoscopy.outcome = 0.90%; p-value = 0.002, n = 363

## Summary: DMM clustering + beta-diversity on all virome composition
## DMM clustering: original abundance table is aggregated on phage genus level. This is done in accordance with 20%AAI (see paper above). Taxa are removed if lower than 0.01% abundance and 20% detection limit. DMM clustering is performed on the rarefied count and visualizated later on a PcOA plot.
## Beta-diversity: original abundance table on genus level is log10 transformed and bray-curtis distances are calculated. This distance metric is subsequently visualized using a PcoA plot using PC1 and PC2. If preferred community-level differences in between metadata groups can be tested using Adonis non-parametric test of vegan package? If more than two groups are compared, a post-hoc Adonis test was used in a pairwise way, corrrecting for multiple testing.
## Metadata: to investigate which metadata covariates contribute to variation in the complete virome community, dbRDA was performed on genus level (bray-curtis distances). Covariates found to significantly contribute to the ordination outcome were furthur implemented in forward modem selection on dbRDA using ordiR2step to determine nonredundant comulative contribution of metadata variables to the variation stepwise dbRDA).
## This works and we also found no significance in the Weeks (baseline or primary endpoint), meaning that perhaps we can use baseline samples as predictive markers for endoscopic outcome (R/NR).
## To find that result we will now repeat this analysis but only with baseline samples
## The most import variable we will plot upon the PcOA plot over all samples to find which way we can explain result.
####################################
# 5.2.4 Visualization effect size variables: donut plot
####################################
# A) Subset dataframe
DF_variables_multi_ <- DF_variables[DF_variables$variable_type == "multivariate",]
DF_variables_multi <- DF_variables_multi_[1:4,]
View(DF_variables)

# B) Fix multivariate data per contributing factor
(DF_variables_multi$variables)
DF_variables_multi$values[DF_variables_multi$variables == "BMI ID"] <- "0.725"
DF_variables_multi$values[DF_variables_multi$variables == "Hemoglobin"] <- "0.717"
DF_variables_multi$values[DF_variables_multi$variables == "Age"] <- "0.553"
DF_variables_multi$values[DF_variables_multi$variables == "Smoking"] <- "0.322"
DF_variables_multi$variable_type <- NULL
View(DF_variables_multi)

# C) Add row with unexplained fraction
# 78.1% explained  variation
100-1.94

# 21.9% unexplained variation
colnames(DF_variables_multi)
df <- data.frame("unexplained","98.06")
names(df)<-c("variables", "values")

newdf <- rbind(DF_variables_multi, df)

# D) Create frequency for donut plot
Effect_size_donut_plot <- as.data.frame(newdf)
Effect_size_donut_plot$values <- as.numeric(Effect_size_donut_plot$values)
#View(Effect_size_donut_plot)

# Make the donut plot
ggplot(Effect_size_donut_plot, aes(x = 2, y = values, fill = variables)) +
  geom_bar(stat = "identity") +
  coord_polar("y", start = 1) +
  # geom_text(aes(label = paste0(round(values), "%")),
        #  position = position_stack(vjust = 0.5)) +
  theme_void() +
  xlim(.2,2.5) +
  scale_fill_brewer(palette = "Greys") +
  theme(legend.title = element_blank())  +
  theme(legend.position = "none") 
   # theme(legend.key.size = unit(2,"line"), legend.text = element_text(size = 12))

# Create dataframe explained versus non-explained
DF_variables_multi_ex <- as.data.frame(c(78.1, 21.9))
DF_variables_multi_ex$values <- DF_variables_multi_ex$`c(78.1, 21.9)`
DF_variables_multi_ex$`c(78.1, 21.9)` <- NULL
DF_variables_multi_ex$variables <- c("explained", "unexplained")

ggplot(DF_variables_multi_ex, aes(x = 2, y = values, fill = variables)) +
  geom_bar(stat = "identity") +
  coord_polar("y", start = 1) +
  geom_text(aes(label = paste0(round(values), "%")),
            position = position_stack(vjust = 0.5)) +
  theme_void() +
  scale_fill_brewer(palette = "Greys") +
  xlim(.2,2.5) +
  theme(legend.title = element_blank())  +
  theme(legend.position = "none")
  theme(legend.key.size = unit(2,"line"), legend.text = element_text(size = 12))
  
# Make both figures with an without legend so it's easier to make in illustrator.
####################################
# 6.Associate covariates with community types: multivariable logistic regression
####################################
# 6.1  Both timepoints
####################################
## Input metadata in multivariate model:
# 1. Patient ID (can drop this now)
# 2. Location (L1, L2, L3)
# 3. Therapy (aTNF, VDZ, UST)
# 4. Moisture (%)

############ 6.2.1 Univariate logistic regression: Location ############ 
Metadata <- data.frame(sample_data(phyloseq_rarefied_phages_transformed))
colnames(Metadata)
Metadata_GLM <- Metadata[,c(12,23)]
Metadata_GLM$timepoints <- NULL

# Add community types
Metadata_GLM <- merge(Metadata_GLM,cluster, by = 0, all = F)
rownames(Metadata_GLM) <- Metadata_GLM$Row.names   
Metadata_GLM$Row.names <- NULL 
#View(Metadata_GLM)

# Convert data to detect presence community type 2 (comm 2 = 1, comm 1 = 0)
names(Metadata_GLM)[names(Metadata_GLM) == "cluster"] <- "Community_type_2"
Metadata_GLM$Community_type_2[Metadata_GLM$Community_type_2 == "1"] <- "0"
Metadata_GLM$Community_type_2[Metadata_GLM$Community_type_2 == "2"] <- "1"
Metadata_GLM$Community_type_2 <- as.numeric(Metadata_GLM$Community_type_2)
str(Metadata_GLM)

# Create multvariable logistic regression model (LRM)
glm.fit <- glm(formula = Community_type_2 ~ Location, family = binomial("logit"), data = Metadata_GLM, na.action(na.omit))
summary(glm.fit)

# Analysis of variance for individual terms
Anova(glm.fit, type="II", test="Wald")

# Pseudo-R-squared
library(pscl)
library(lmtest)
pR2(glm.fit) # McFadden output*100 is best measure for R2
lrtest(glm.fit) # Significance of model

# Calculate relative risk to visualize.
coef(glm.fit)
RR <- exp(coef(glm.fit))[-1]
RR
exp(confint(glm.fit))[-1,]

## Simply show the significant ones out of the GLR model. Therefore in this case only show remission
glm.fit$coefficients

# Conclusion: Location is not significant
############ 6.2.2 Univariate logistic regression: Age ############ 
Metadata <- data.frame(sample_data(phyloseq_rarefied_phages_transformed))
colnames(Metadata)
Metadata_GLM <- Metadata[,c(4,23)]
Metadata_GLM$timepoints <- NULL

# Add community types
Metadata_GLM <- merge(Metadata_GLM,cluster, by = 0, all = F)
rownames(Metadata_GLM) <- Metadata_GLM$Row.names   
Metadata_GLM$Row.names <- NULL 

# Convert data to detect presence community type 2 (comm 2 = 1, comm 1 = 0)
names(Metadata_GLM)[names(Metadata_GLM) == "cluster"] <- "Community_type_2"
Metadata_GLM$Community_type_2[Metadata_GLM$Community_type_2 == "1"] <- "0"
Metadata_GLM$Community_type_2[Metadata_GLM$Community_type_2 == "2"] <- "1"
Metadata_GLM$Community_type_2 <- as.numeric(Metadata_GLM$Community_type_2)
str(Metadata_GLM)

# Create multvariable logistic regression model (LRM)
glm.fit <- glm(formula = Community_type_2 ~ Age, family = binomial("logit"), data = Metadata_GLM, na.action(na.omit))
summary(glm.fit)

# Analysis of variance for individual terms
Anova(glm.fit, type="II", test="Wald")
# Conclusion: Location is not significant

# Pseudo-R-squared
pR2(glm.fit) # McFadden output*100 is best measure for R2
lrtest(glm.fit) # Significance of model

# Calculate relative risk to visualize.
RR <- exp(coef(glm.fit))[-1]
RR
exp(confint(glm.fit))[-1,]

## Simply show the significant ones out of the GLR model. Therefore in this case only show remission
glm.fit$coefficients

## RR figure
variable <- c("Age")
RR_value <- c(1.03)
lowerCI <- c(1.01)
upperCI <- c(1.04)
RR_df <- data.frame(variable, RR_value,lowerCI,upperCI)
#View(RR_df)

ggplot(RR_df) +
  geom_point(aes(x = variable, y = RR_value), size = 5,  color = "black") +
  geom_errorbar(aes(x = variable, y = RR_value, ymin = lowerCI, ymax = upperCI), width = 0, color = "black", ) +
  geom_text(aes(x = variable, y = RR_value, label = round(RR_value, 2)), vjust = - 1.5, size = 3) +
  ylim(0,2) +
  geom_hline(yintercept=1, linetype='dotted', col = 'black') +
  coord_flip() +
  xlab("") +
  ylab("Community typ 2 Relative risk") +
  theme_bw()

## Age 
# Create figure logistic regression
ggplot(Metadata_GLM, aes(x=Age, y=Community_type_2)) + 
  stat_smooth(method="glm", se=TRUE, method.args = list(family=binomial), col="#b87a4d", fill="#b87a4d", lty=1) +
  ylab("Community type 2 prevalence") +
  theme_bw() +
  ylim(0,1) 
############ 6.2.3 Univariate logistic regression: Moisture ############ 
Metadata <- data.frame(sample_data(phyloseq_rarefied_phages_transformed))
colnames(Metadata)
Metadata_GLM <- Metadata[,c(7,23)]
Metadata_GLM$timepoints <- NULL

# Add community types
Metadata_GLM <- merge(Metadata_GLM,cluster, by = 0, all = F)
rownames(Metadata_GLM) <- Metadata_GLM$Row.names   
Metadata_GLM$Row.names <- NULL 

# Convert data to detect presence community type 2 (comm 2 = 1, comm 1 = 0)
names(Metadata_GLM)[names(Metadata_GLM) == "cluster"] <- "Community_type_2"
Metadata_GLM$Community_type_2[Metadata_GLM$Community_type_2 == "1"] <- "0"
Metadata_GLM$Community_type_2[Metadata_GLM$Community_type_2 == "2"] <- "1"
Metadata_GLM$Community_type_2 <- as.numeric(Metadata_GLM$Community_type_2)
str(Metadata_GLM)

# Create multvariable logistic regression model (LRM)
glm.fit <- glm(formula = Community_type_2 ~ Moisture, family = binomial("logit"), data = Metadata_GLM, na.action(na.omit))
summary(glm.fit)

# Analysis of variance for individual terms
Anova(glm.fit, type="II", test="Wald")

# Pseudo-R-squared
pR2(glm.fit) # McFadden output*100 is best measure for R2
lrtest(glm.fit) # Significance of model

# Calculate relative risk to visualize.
RR <- exp(coef(glm.fit))[-1]
RR
exp(confint(glm.fit))[-1,]

## Simply show the significant ones out of the GLR model. Therefore in this case only show remission
glm.fit$coefficients

# Conclusion: Location is not significant
############ 6.2.4 Create Multiple logistic regression model (LRM) ############ 
# Select variables from multivariate model & appropriate samples
Metadata <- data.frame(sample_data(phyloseq_rarefied_phages_transformed))
Metadata_GLM <- Metadata
colnames(Metadata_GLM)
Metadata_GLM <- Metadata_GLM[,c(7,23)]

# Add community types
Metadata_GLM <- merge(Metadata_GLM,cluster, by = 0, all = F)
rownames(Metadata_GLM) <- Metadata_GLM$Row.names   
Metadata_GLM$Row.names <- NULL 

# Convert data to detect presence community type 2 (comm 2 = 1, comm 1 = 0)
names(Metadata_GLM)[names(Metadata_GLM) == "cluster"] <- "Community_type_2"
Metadata_GLM$Community_type_2[Metadata_GLM$Community_type_2 == "1"] <- "0"
Metadata_GLM$Community_type_2[Metadata_GLM$Community_type_2 == "2"] <- "1"
Metadata_GLM$Community_type_2 <- as.numeric(Metadata_GLM$Community_type_2)
str(Metadata_GLM)

# Create multvariable logistic regression model (LRM)
glm.fit <- glm(formula = Community_type_2 ~ Moisture,family = binomial("logit"), data = Metadata_GLM)
summary(glm.fit)

# Analysis of variance for individual terms
Anova(glm.fit, type="II", test="Wald")

# Pseudo-R-squared
pR2(glm.fit)*100 # McFadden output*100 is best measure for R2
lrtest(glm.fit) # Significance of model

# Calculate relative risk to visualize.
coef(glm.fit)[-1]
exp(coef(glm.fit))[-1]
exp(confint(glm.fit))[-1,]

## RR figure

# MAKE FIGURE & FINISH WITH BASELINE
variable <- c("Location [L2CD/UC]", "Location [L3]","Age","Moisture")
RR_value <- c(0.71,0.84,1.03,0.98)
lowerCI <- c(0.355,0.463,1.01,0.934)
upperCI <- c(1.38,1.51,1.04,1.02)
RR_df <- data.frame(variable, RR_value,lowerCI,upperCI)

ggplot(RR_df) +
  geom_point(aes(x = variable, y = RR_value), size = 5,  color = "black") +
  geom_errorbar(aes(x = variable, y = RR_value, ymin = lowerCI, ymax = upperCI), width = 0, color = "black", ) +
  geom_text(aes(x = variable, y = RR_value, label = round(RR_value, 2)), vjust = - 1.5, size = 3) +
  ylim(0,1.5) +
  geom_hline(yintercept=1, linetype='dotted', col = 'black') +
  coord_flip() +
  xlab("") +
  ylab("Relative risk: community type M") +
  theme_bw()

# Figure Therapy
Metadata_GLM_ <- Metadata_GLM
Metadata_GLM_$Hemoglobin <- NULL
Metadata_GLM_$Community_type_2[Metadata_GLM_$Community_type_2 == "0"] <- "community type 1"
Metadata_GLM_$Community_type_2[Metadata_GLM_$Community_type_2 == "1"] <- "community type 2"

# IBD patients
# Create frequency and proportion table for IBD patients + Statistics
Metadata_GLM_table <- table(Metadata_GLM_)
Metadata_GLM_table_prop <-  as.data.frame(prop.table(Metadata_GLM_table, margin = 1))

ggplot(Metadata_GLM_table_prop) +
  aes(x = Therapy, fill = Community_type_2, weight = Freq) +
  geom_bar(position = "dodge", alpha=0.9) +
  ylab("Community type prevalence") +
  scale_fill_manual(values = list(`community type 1` = "#4d74b8", `community type 2` = "#b87a4d")) +
  theme_classic()

chisq.test(Metadata_GLM_table, correct = F)
cramer_v(Metadata_GLM_table)

# Community type 1
prop.test(x = c(104,58), n = c(158,73), correct = F) #  TNF vs UST***
prop.test(x = c(104,80), n = c(158,132), correct = F) #  TNF vs VDZ
prop.test(x = c(58,80), n = c(73,132), correct = F) #  UST vs VDZ***

# Community type 2
prop.test(x = c(54,15), n = c(158,73), correct = F) #  TNF vs UST***
prop.test(x = c(54,52), n = c(158,132), correct = F) #  TNF vs VDZ
prop.test(x = c(15,52), n = c(73,132), correct = F) #  UST vs VDZ***
####################################
# 6.2  Baseline
####################################
# Run lines before again but now specifiy to subselect baseline samples in creating phyloseq object

## Input metadata in multivariate model:
# 1. BMI
# 2. Hemoglobin 
# 3. Age 

############ 6.2.1 Univariate logistic regression: BMI ############ 
Metadata <- data.frame(sample_data(phyloseq_rarefied_phages_transformed))
colnames(Metadata)
Metadata_GLM <- Metadata[,c(9,22)]
Metadata_GLM <- Metadata_GLM[Metadata_GLM$timepoints == "pre-intervention + CD w14",]
Metadata_GLM$timepoints <- NULL

# Add community types
Metadata_GLM <- merge(Metadata_GLM,cluster, by = 0, all = F)
rownames(Metadata_GLM) <- Metadata_GLM$Row.names   
Metadata_GLM$Row.names <- NULL 

# Convert data to detect presence community type 2 (comm 2 = 1, comm 1 = 0)
names(Metadata_GLM)[names(Metadata_GLM) == "cluster"] <- "Community_type_2"
Metadata_GLM$Community_type_2[Metadata_GLM$Community_type_2 == "1"] <- "0"
Metadata_GLM$Community_type_2[Metadata_GLM$Community_type_2 == "2"] <- "1"
Metadata_GLM$Community_type_2 <- as.numeric(Metadata_GLM$Community_type_2)
str(Metadata_GLM)

# Create multvariable logistic regression model (LRM)
glm.fit <- glm(formula = Community_type_2 ~ BMI, family = binomial("logit"), data = Metadata_GLM, na.action(na.omit))
summary(glm.fit)

# Analysis of variance for individual terms
Anova(glm.fit, type="II", test="Wald")
# Conclusion: BMI is not significant

# Pseudo-R-squared
pR2(glm.fit) # McFadden output*100 is best measure for R2
lrtest(glm.fit) # Significance of model

# Calculate relative risk to visualize.
coef(glm.fit)[-1]
exp(coef(glm.fit))[-1]
exp(confint(glm.fit))[-1,]

############ 6.2.2 Univariate logistic regression: Hemoglobin ############ 
Metadata <- data.frame(sample_data(phyloseq_rarefied_phages_transformed))
colnames(Metadata)
Metadata_GLM <- Metadata[,c(13,22)]
Metadata_GLM <- Metadata_GLM[Metadata_GLM$timepoints == "pre-intervention + CD w14",]
Metadata_GLM$timepoints <- NULL

# Add community types
Metadata_GLM <- merge(Metadata_GLM,cluster, by = 0, all = F)
rownames(Metadata_GLM) <- Metadata_GLM$Row.names   
Metadata_GLM$Row.names <- NULL 

# Convert data to detect presence community type 2 (comm 2 = 1, comm 1 = 0)
names(Metadata_GLM)[names(Metadata_GLM) == "cluster"] <- "Community_type_2"
Metadata_GLM$Community_type_2[Metadata_GLM$Community_type_2 == "1"] <- "0"
Metadata_GLM$Community_type_2[Metadata_GLM$Community_type_2 == "2"] <- "1"
Metadata_GLM$Community_type_2 <- as.numeric(Metadata_GLM$Community_type_2)
str(Metadata_GLM)

# Create multvariable logistic regression model (LRM)
glm.fit <- glm(formula = Community_type_2 ~ Hemoglobin, family = binomial("logit"), data = Metadata_GLM, na.action(na.omit))
summary(glm.fit)

# Analysis of variance for individual terms
Anova(glm.fit, type="II", test="Wald")
# Conclusion: Hemoglobin is not significant

# Pseudo-R-squared
pR2(glm.fit) # McFadden output*100 is best measure for R2
lrtest(glm.fit) # Significance of model

# Calculate relative risk to visualize.
coef(glm.fit)[-1]
exp(coef(glm.fit))[-1]
exp(confint(glm.fit))[-1,]

############ 6.2.3 Univariate logistic regression: Age ############ 
Metadata <- data.frame(sample_data(phyloseq_rarefied_phages_transformed))
colnames(Metadata)
Metadata_GLM <- Metadata[,c(3,22)]
Metadata_GLM <- Metadata_GLM[Metadata_GLM$timepoints == "pre-intervention + CD w14",]
Metadata_GLM$timepoints <- NULL

# Add community types
Metadata_GLM <- merge(Metadata_GLM,cluster, by = 0, all = F)
rownames(Metadata_GLM) <- Metadata_GLM$Row.names   
Metadata_GLM$Row.names <- NULL 

# Convert data to detect presence community type 2 (comm 2 = 1, comm 1 = 0)
names(Metadata_GLM)[names(Metadata_GLM) == "cluster"] <- "Community_type_2"
Metadata_GLM$Community_type_2[Metadata_GLM$Community_type_2 == "1"] <- "0"
Metadata_GLM$Community_type_2[Metadata_GLM$Community_type_2 == "2"] <- "1"
Metadata_GLM$Community_type_2 <- as.numeric(Metadata_GLM$Community_type_2)
str(Metadata_GLM)

# Create multvariable logistic regression model (LRM)
glm.fit <- glm(formula = Community_type_2 ~ Age, family = binomial("logit"), data = Metadata_GLM, na.action(na.omit))
summary(glm.fit)

# Analysis of variance for individual terms
Anova(glm.fit, type="II", test="Wald")
# Conclusion: Hemoglobin is not significant

# Pseudo-R-squared
pR2(glm.fit) # McFadden output*100 is best measure for R2
lrtest(glm.fit) # Significance of model

# Calculate relative risk to visualize.
coef(glm.fit)[-1]
exp(coef(glm.fit))[-1]
exp(confint(glm.fit))[-1,]

############ 6.2.4 Univariate logistic regression: Smoking ############ 
Metadata <- data.frame(sample_data(phyloseq_rarefied_phages_transformed))
colnames(Metadata)
Metadata_GLM <- Metadata[,c(12,22)]
Metadata_GLM <- Metadata_GLM[Metadata_GLM$timepoints == "pre-intervention + CD w14",]
Metadata_GLM$timepoints <- NULL

# Add community types
Metadata_GLM <- merge(Metadata_GLM,cluster, by = 0, all = F)
rownames(Metadata_GLM) <- Metadata_GLM$Row.names   
Metadata_GLM$Row.names <- NULL 

# Convert data to detect presence community type 2 (comm 2 = 1, comm 1 = 0)
names(Metadata_GLM)[names(Metadata_GLM) == "cluster"] <- "Community_type_2"
Metadata_GLM$Community_type_2[Metadata_GLM$Community_type_2 == "1"] <- "0"
Metadata_GLM$Community_type_2[Metadata_GLM$Community_type_2 == "2"] <- "1"
Metadata_GLM$Community_type_2 <- as.numeric(Metadata_GLM$Community_type_2)
str(Metadata_GLM)

# Create multvariable logistic regression model (LRM)
glm.fit <- glm(formula = Community_type_2 ~ Smoking, family = binomial("logit"), data = Metadata_GLM, na.action(na.omit))
summary(glm.fit)

# Analysis of variance for individual terms
Anova(glm.fit, type="II", test="Wald")
# Conclusion: Hemoglobin is not significant

# Pseudo-R-squared
pR2(glm.fit) # McFadden output*100 is best measure for R2
lrtest(glm.fit) # Significance of model
#p.adjust(p = 0.402, method ="BH", n = 4)

# Calculate relative risk to visualize.
coef(glm.fit)[-1]
exp(coef(glm.fit))[-1]
exp(confint(glm.fit))[-1,]

############ 6.2.5 Create multvariable logistic regression model (LRM) ############ 
# Select variables from multivariate model & appropriate samples
Metadata <- data.frame(sample_data(phyloseq_rarefied_phages_transformed))
Metadata_GLM <- Metadata[Metadata$timepoints == "baseline",]
colnames(Metadata_GLM)
Metadata_GLM <- Metadata_GLM[,c(4,14,4)]

# Add community types
Metadata_GLM <- merge(Metadata_GLM,cluster, by = 0, all = F)
rownames(Metadata_GLM) <- Metadata_GLM$Row.names   
Metadata_GLM$Row.names <- NULL 

# Convert data to detect presence community type 2 (comm 2 = 1, comm 1 = 0)
names(Metadata_GLM)[names(Metadata_GLM) == "cluster"] <- "Community_type_2"
Metadata_GLM$Community_type_2[Metadata_GLM$Community_type_2 == "1"] <- "0"
Metadata_GLM$Community_type_2[Metadata_GLM$Community_type_2 == "2"] <- "1"
Metadata_GLM$Community_type_2 <- as.numeric(Metadata_GLM$Community_type_2)
str(Metadata_GLM)

# Create multvariable logistic regression model (LRM)
glm.fit <- glm(formula = Community_type_2 ~ Hemoglobin + Age,family = binomial("logit"), data = Metadata_GLM)
summary(glm.fit)

# Analysis of variance for individual terms
Anova(glm.fit, type="II", test="Wald")

# Pseudo-R-squared
pR2(glm.fit)*100 # McFadden output*100 is best measure for R2
lrtest(glm.fit) # Significance of model

# Calculate relative risk to visualize.
exp(coef(glm.fit))[-1]
exp(confint(glm.fit))[-1,]

## RR figure
variable <- c("BMI","Hemoglobin","Age","Smoking")
RR_value <- c(1.07,1.41,1.03,1.41)
lowerCI <- c(1.01,1.13,1.01,0.62)
upperCI <- c(1.14,1.79,1.05,3.14)
RR_df <- data.frame(variable, RR_value,lowerCI,upperCI)
View(RR_df)

ggplot(RR_df) +
  geom_point(aes(x = variable, y = RR_value), size = 5,  color = "black") +
  geom_errorbar(aes(x = variable, y = RR_value, ymin = lowerCI, ymax = upperCI), width = 0, color = "black", ) +
  geom_text(aes(x = variable, y = RR_value, label = round(RR_value, 2)), vjust = - 1.5, size = 3) +
  ylim(0.9,3.2) +
  geom_hline(yintercept=1, linetype='dotted', col = 'black') +
  coord_flip() +
  xlab("") +
  ylab("Relative risk: community type M") +
  theme_bw()

## Hemoglobin
# Create figure logistic regression
ggplot(Metadata_GLM, aes(x=Hemoglobin, y=Community_type_2)) + 
  stat_smooth(method="glm", se=TRUE, method.args = list(family=binomial), col="#b87a4d", fill="#b87a4d", lty=1) +
  ylab("Community type 2 prevalence") +
  theme_bw() +
  ylim(0,1)

## Age
# Create figure logistic regression
ggplot(Metadata_GLM, aes(x=Age, y=Community_type_2)) + 
  stat_smooth(method="glm", se=TRUE, method.args = list(family=binomial), col="#b87a4d", fill="#b87a4d", lty=1) +
  ylab("Community type 2 prevalence") +
  theme_bw() +
  ylim(0,1)
############ 6.2.6 Association between BMI & HG: explains Frailty ####
Metadata <- data.frame(sample_data(phyloseq_rarefied_phages_transformed))
colnames(Metadata)
Metadata_GLM <- Metadata[,c(11,15)]

# Add community types
Metadata_GLM <- merge(Metadata_GLM,cluster, by = 0, all = F)
rownames(Metadata_GLM) <- Metadata_GLM$Row.names   
Metadata_GLM$Row.names <- NULL 
#View(Metadata_GLM)

# Convert data to detect presence community type 2 (comm 2 = 1, comm 1 = 0)
names(Metadata_GLM)[names(Metadata_GLM) == "cluster"] <- "Community_type_2"
Metadata_GLM$Community_type_2[Metadata_GLM$Community_type_2 == "1"] <- "Viral community type CA"
Metadata_GLM$Community_type_2[Metadata_GLM$Community_type_2 == "2"] <- "Viral community type CrM"
#View(Metadata_GLM)
#str(Metadata_GLM)

# Correlate BMI and Hemoglobin
ggplot(Metadata_GLM) +
 aes(x = Hemoglobin, y = BMI, colour = Community_type_2) +
 geom_point(shape = "circle", 
 size = 1.5) +
 geom_smooth(span = 0.75) +
 scale_color_manual(values = c(`Viral community type CA` = "#BD801C", `Viral community type CrM` = "#2B8A92")) +
 theme_bw()

cor.test(Metadata_GLM$Hemoglobin, Metadata_GLM$BMI, method = c("spearman"))
# Positive correlation meaning the higher the BMI the higher the hemoglobin level
# Or reverselt, the lower the BMI the lower the hemoglobin levels
# Thereby confirming frailty of IBD pateints as a potential factor.
####################################
# 6.3  PE
####################################
# Run lines before again but now specifiy to subselect PE samples in creating phyloseq object

## Input metadata in multivariate model:
# 1. Location (L1, L2, L3)
# 2. Endoscopic outcome (R/NR)
############  6.3.1 Univariate logistic regression: Location ############ 

# Select variables from multivariate model & appropriate samples
Metadata <- data.frame(sample_data(phyloseq_rarefied_phages_transformed))
colnames(Metadata)
Metadata_GLM <- Metadata[,c(11,22)]
colnames(Metadata_GLM)
Metadata_GLM <- Metadata_GLM[Metadata_GLM$timepoints == "post-intervention",]
Metadata_GLM$timepoints <- NULL

# Add community types
Metadata_GLM <- merge(Metadata_GLM,cluster, by = 0, all = F)
rownames(Metadata_GLM) <- Metadata_GLM$Row.names   
Metadata_GLM$Row.names <- NULL 

# Convert data to detect presence community type 2 (comm 2 = 1, comm 1 = 0)
names(Metadata_GLM)[names(Metadata_GLM) == "cluster"] <- "Community_type_2"
Metadata_GLM$Community_type_2[Metadata_GLM$Community_type_2 == "1"] <- "0"
Metadata_GLM$Community_type_2[Metadata_GLM$Community_type_2 == "2"] <- "1"
Metadata_GLM$Community_type_2 <- as.numeric(Metadata_GLM$Community_type_2)
str(Metadata_GLM)
#View(Metadata_GLM)

# Create multvariable logistic regression model (LRM)
#install.packages("pscl")
library(pscl)
glm.fit <- glm(formula = Community_type_2 ~ Location, family = binomial("logit"), data = Metadata_GLM, na.action(na.omit))
summary(glm.fit)

# Analysis of variance for individual terms
Anova(glm.fit, type="II", test="Wald")

# Pseudo-R-squared
pR2(glm.fit) # McFadden output*100 is best measure for R2
lrtest(glm.fit) # Significance of model

# Calculate relative risk to visualize.
RR <- exp(coef(glm.fit))[-1]
exp(confint(glm.fit))[-1,]

# Coefficients
glm.fit$coefficients

# Conclusion: Location is not significant
############  6.3.2 Univariate logistic regression: Endoscopic outcome ############  

# Select variables from multivariate model & appropriate samples
Metadata <- data.frame(sample_data(phyloseq_rarefied_phages_transformed))
colnames(Metadata)
Metadata_GLM <- Metadata[,c(21,25)]
colnames(Metadata_GLM)
Metadata_GLM <- Metadata_GLM[Metadata_GLM$timepoints == "post-intervention",]
Metadata_GLM$timepoints <- NULL
is.na(Metadata_GLM$Endoscopy.outcome) <- Metadata_GLM$Endoscopy.outcome == "unknown"

# Add community types
Metadata_GLM <- merge(Metadata_GLM,cluster, by = 0, all = F)
rownames(Metadata_GLM) <- Metadata_GLM$Row.names   
Metadata_GLM$Row.names <- NULL 
View(Metadata_GLM)

# Convert data to detect presence community type 2 (comm 2 = 1, comm 1 = 0)
names(Metadata_GLM)[names(Metadata_GLM) == "cluster"] <- "Community_type_2"
Metadata_GLM$Community_type_2[Metadata_GLM$Community_type_2 == "1"] <- "0"
Metadata_GLM$Community_type_2[Metadata_GLM$Community_type_2 == "2"] <- "1"
Metadata_GLM$Community_type_2 <- as.numeric(Metadata_GLM$Community_type_2)
str(Metadata_GLM)
View(Metadata_GLM)

# Create multvariable logistic regression model (LRM)
glm.fit <- glm(formula = Community_type_2 ~ Endoscopy.outcome, family = binomial("logit"), data = Metadata_GLM, na.action(na.exclude))
summary(glm.fit)

# Analysis of variance for individual terms
Anova(glm.fit, type="II", test="Wald")

# Pseudo-R-squared
pR2(glm.fit) # McFadden output*100 is best measure for R2
lrtest(glm.fit) # Significance of model

# Calculate relative risk to visualize.
exp(coef(glm.fit))[-1]
exp(confint(glm.fit))[-1,]
coef(glm.fit)

## Relative risk figure
variable <- c("Endoscopic remission")
RR_value <- c(2.65)
lowerCI <- c(1.27)
upperCI <- c(5.76)
RR_df <- data.frame(variable, RR_value,lowerCI,upperCI)
#View(RR_df)

ggplot(RR_df) +
  geom_point(aes(x = variable, y = RR_value), size = 5,  color = "black") +
  geom_errorbar(aes(x = variable, y = RR_value, ymin = lowerCI, ymax = upperCI), width = 0, color = "black", ) +
  geom_text(aes(x = variable, y = RR_value, label = round(RR_value, 2)), vjust = - 1.5, size = 3) +
  ylim(0,6) +
  geom_hline(yintercept=1, linetype='dotted', col = 'black') +
  coord_flip() +
  xlab("") +
  ylab("Community typ 2 Relative risk") +
  theme_bw()
# You can still remove grid lines & add number 1, but is easier in illustrator
## RR of 2.65 for endoscopic outcome to predict Comm type 2.

## Endoscopy.outcome
Metadata_GLM_ <- Metadata_GLM
Metadata_GLM_$Endoscopy.outcome[is.na(Metadata_GLM_$Endoscopy.outcome)] <- "unknown"
Metadata_GLM_$Community_type_2[Metadata_GLM_$Community_type_2 == "0"] <- "community type 1"
Metadata_GLM_$Community_type_2[Metadata_GLM_$Community_type_2 == "1"] <- "community type 2"

# IBD patients
# Create frequency and proportion table for IBD patients + Statistics
Metadata_GLM_table <- table(Metadata_GLM_)
Metadata_GLM_table_prop <-  as.data.frame(prop.table(Metadata_GLM_table, margin = 2))

ggplot(Metadata_GLM_table_prop) +
  aes(x = Community_type_2, fill = Endoscopy.outcome, weight = Freq) +
  geom_bar(position = "dodge") +
  scale_fill_manual(values = list(`non-remission` = "#ad4b38", `remission` = "#499b78", `unknown` = "#354b9a")) +
  theme_classic() +
  ggtitle("IBD") +
  ylab("Endoscopic outcome (%)") +
  theme(axis.title.y = element_text(size= 16))

# poster
ggplot(Metadata_GLM_table_prop) +
  aes(x = Community_type_2, fill = Endoscopy.outcome, weight = Freq) +
  geom_bar(position = "dodge") +
  scale_fill_manual(values = list(`non-remission` = "#a7d3b2", `remission` = "#66996b", `unknown` = "#2a4a20")) +
  theme_classic() +
  ggtitle("IBD") +
  xlab("") +
  ylab("Endoscopic outcome (%)") +
  theme(axis.title.y = element_text(size= 16)) +
  coord_flip()

# Statistics (this is immediately corrected by BH)
# We do not need a multiple correction test since we are comparing only 2 values
# These comparisons are done by proportion test between 2 values
# I only need to add effect size (using CramersPhi) if I compare two or more & subsequently compare post-hoc comparisons

# Comparing 2 variables only
Metadata_GLM_table
prop.test(x = c(30,42), n = c(113,65),correct =F) # R
prop.test(x = c(52,14), n = c(113,65),correct =F) # NR
prop.test(x = c(19,9), n = c(113,65),correct =F) # unknown

# CD patients
colnames(Metadata)
Metadata_GLM_CD <- Metadata[,c(3,20,23,24)]
Metadata_GLM_CD <- Metadata_GLM_CD[Metadata_GLM_CD$timepoints == "primary endpoint",]
Metadata_GLM_CD <- Metadata_GLM_CD[Metadata_GLM_CD$Diagnosis == "CD",]
Metadata_GLM_CD$Diagnosis <- NULL
Metadata_GLM_CD$timepoints <- NULL
Metadata_GLM_CD$cluster[Metadata_GLM_CD$cluster == "1"] <- "community type 1"
Metadata_GLM_CD$cluster[Metadata_GLM_CD$cluster == "2"] <- "community type 2"
nrow(Metadata_GLM_CD)

# Create tables
Metadata_GLM_table <- table(Metadata_GLM_CD)
Metadata_GLM_table_prop <-  as.data.frame(prop.table(Metadata_GLM_table, margin = 2))

ggplot(Metadata_GLM_table_prop) +
  aes(x = cluster, fill = Endoscopy.outcome, weight = Freq) +
  geom_bar(position = "dodge") +
  scale_fill_manual(values = list(`non-remission` = "#ad4b38", `remission` = "#499b78", `unknown` = "#354b9a")) +
  theme_classic() +
  ggtitle("CD") +
  ylab("Endoscopic outcome (%)") +
  theme(axis.title.y = element_text(size= 16))

# poster
ggplot(Metadata_GLM_table_prop) +
  aes(x = cluster, fill = Endoscopy.outcome, weight = Freq) +
  geom_bar(position = "dodge") +
  scale_fill_manual(values = list(`non-remission` = "#a7d3b2", `remission` = "#66996b", `unknown` = "#2a4a20")) +
  theme_classic() +
  ggtitle("CD") +
  xlab("") +
  ylab("Endoscopic outcome (%)") +
  theme(axis.title.y = element_text(size= 16)) +
  coord_flip()
  
# Statistics: Comparing 2 variables only
Metadata_GLM_table
prop.test(x = c(37,9), n = c(84,31), correct =F) # R
prop.test(x = c(28,13), n = c(84,31), correct =F) # NR
prop.test(x = c(19,9), n = c(84,31), correct =F) # unknown

# UC patients
colnames(Metadata)
colnames(Metadata_GLM_CD)
Metadata_GLM_CD <- Metadata[,c(4,21,24,25)]
Metadata_GLM_CD <- Metadata_GLM_CD[Metadata_GLM_CD$timepoints == "post-intervention",]
Metadata_GLM_CD <- Metadata_GLM_CD[Metadata_GLM_CD$Diagnosis == "UC",]
Metadata_GLM_CD$Diagnosis <- NULL
Metadata_GLM_CD$timepoints <- NULL
Metadata_GLM_CD$cluster[Metadata_GLM_CD$cluster == "1"] <- "community type 1"
Metadata_GLM_CD$cluster[Metadata_GLM_CD$cluster == "2"] <- "community type 2"
nrow(Metadata_GLM_CD)
Metadata_GLM_CD

# Create tables
Metadata_GLM_table <- table(Metadata_GLM_CD)
Metadata_GLM_table_prop <-  as.data.frame(prop.table(Metadata_GLM_table, margin = 2))
Metadata_GLM_table_prop
ggplot(Metadata_GLM_table_prop) +
  aes(x = cluster, fill = Endoscopy.outcome, weight = Freq) +
  geom_bar(position = "dodge") +
  scale_fill_manual(values = list(`non-remission` = "#ad4b38", `remission` = "#499b78")) +
  theme_classic() +
  ggtitle("UC") +
  ylab("Endoscopic outcome (%)") +
  theme(axis.title.y = element_text(size= 16))

# poster
ggplot(Metadata_GLM_table_prop) +
  aes(x = cluster, fill = Endoscopy.outcome, weight = Freq) +
  geom_bar(position = "dodge") +
  scale_fill_manual(values = list(`non-remission` = "#a7d3b2", `remission` = "#66996b", `unknown` = "#2a4a20")) +
  theme_classic() +
  ggtitle("UC") +
  xlab("") +
  ylab("Endoscopic outcome (%)") +
  theme(axis.title.y = element_text(size= 16)) +
  coord_flip()

# Statistics= Comparing 2 variables only
Metadata_GLM_table
prop.test(x = c(15,5), n = c(29,22), correct = F) # R
prop.test(x = c(14,17), n = c(29,22), correct = F) # NR

#### ALTERNATIVELY SEEK COMMUNITY TYPE PREVALENCE OVER VARIABLE
# IBD patients
Metadata_GLM_ <- Metadata_GLM
Metadata_GLM_$Endoscopy.outcome[is.na(Metadata_GLM_$Endoscopy.outcome)] <- "unknown"
Metadata_GLM_$Community_type_2[Metadata_GLM_$Community_type_2 == "0"] <- "community type 1"
Metadata_GLM_$Community_type_2[Metadata_GLM_$Community_type_2 == "1"] <- "community type 2"

# IBD patients
# Create frequency and proportion table for IBD patients + Statistics
Metadata_GLM_table <- table(Metadata_GLM_)
Metadata_GLM_table_prop <-  as.data.frame(prop.table(Metadata_GLM_table, margin = 1))
Metadata_GLM_table_prop$Freq <- Metadata_GLM_table_prop$Freq*100

# Barplot next to eachother
ggplot(Metadata_GLM_table_prop) +
 aes(x = Endoscopy.outcome, fill = Community_type_2, weight = Freq) +
 geom_bar(position = "dodge", alpha=0.9) +
 ylab("Community type prevalence") +
 scale_fill_manual(values = list(`community type 1` = "#bd801c", `community type 2` = "#2b8a92")) +
 theme_classic()

# Stacked barplot
ggplot(Metadata_GLM_table_prop) +
 aes(x = Endoscopy.outcome, fill = Community_type_2, weight = Freq) +
 geom_bar() +
 ylab("Community type prevalence") +
 scale_fill_manual(values = list(`community type 1` = "#bd801c", `community type 2` = "#2b8a92")) +
 theme_classic()

# Statistics (this is immediately corrected by BH)
# Because we now compare three groups you can do chi-squared with post-hoc proportion test
# And can also calculate cramersPhi for effect size
Metadata_GLM_table_prop
chisq.test(Metadata_GLM_table, correct = F)
cramer_v(Metadata_GLM_table)
str(Metadata_GLM_table)
sum(Metadata_GLM_table)

pairwise_prop_test(Metadata_GLM_table, p.adjust.method = "BH", correct = FALSE)

# Community type 2
Metadata_GLM_table
prop.test(x = c(14,30), n = c(66,72), correct = F) #  R vs NR****
prop.test(x = c(30,9), n = c(72,28), correct = F) #  R vs unknown
prop.test(x = c(14,9), n = c(66,28), correct = F) #  NR vs unknown

# Community type 1
prop.test(x = c(52,42,19), n = c(66,72,28), correct = F) #  similar as previous line
prop.test(x = c(52,42), n = c(66,72), correct = F) #  R vs NR****
prop.test(x = c(42,19), n = c(72,28), correct = F) #  R vs unknown
prop.test(x = c(52,19), n = c(66,28), correct = F) #  NR vs unknown
p.adjust(0.0100, method = "BH", n = 2)

# Prevalence of community type 2 is higher in remission compared to non-remission

# CD patients
colnames(Metadata)
Metadata_GLM_CD <- Metadata[,c(4,5,21,24,25)]
Metadata_GLM_CD <- Metadata_GLM_CD[Metadata_GLM_CD$timepoints == "post-intervention",]
Metadata_GLM_CD <- Metadata_GLM_CD[Metadata_GLM_CD$Diagnosis == "CD",]
Metadata_GLM_CD$Diagnosis <- NULL
Metadata_GLM_CD$timepoints <- NULL
Metadata_GLM_CD$Age <- NULL
Metadata_GLM_CD$cluster[Metadata_GLM_CD$cluster == "1"] <- "community type 1"
Metadata_GLM_CD$cluster[Metadata_GLM_CD$cluster == "2"] <- "community type 2"
#View(Metadata_GLM_CD)

# Proportion tables
table(Metadata_GLM_CD)
Metadata_GLM_table_CD <- table(Metadata_GLM_CD)
Metadata_GLM_table_prop_CD <-  as.data.frame(prop.table(Metadata_GLM_table_CD, margin = 1))
Metadata_GLM_table_prop_CD
sum(Metadata_GLM_table_CD)

ggplot(Metadata_GLM_table_prop_CD) +
  aes(x = Endoscopy.outcome, fill = cluster, weight = Freq) +
  geom_bar(position = "dodge", alpha = 0.9) +
  ylab("Community type prevalence") +
  scale_fill_manual(values = list(`community type 1` = "#4d74b8", `community type 2` = "#b87a4d")) +
  theme_classic()

# Stacked barplot
ggplot(Metadata_GLM_table_prop_CD) +
  aes(x = Endoscopy.outcome, fill = cluster, weight = Freq) +
  geom_bar() +
  ylab("Community type prevalence") +
  scale_fill_manual(values = list(`community type 1` = "#bd801c", `community type 2` = "#2b8a92")) +
  theme_classic()

# Statistics (this is immediately corrected by BH)
chisq.test(Metadata_GLM_table_CD, correct = F) # NS
cramer_v(Metadata_GLM_table_CD)
pairwise_prop_test(Metadata_GLM_table_CD, p.adjust.method = "BH", correct = FALSE)
sum(Metadata_GLM_table_CD)

# UC patients
colnames(Metadata)
Metadata_GLM_CD <- Metadata[,c(4,21,24,25)]
Metadata_GLM_CD <- Metadata_GLM_CD[Metadata_GLM_CD$timepoints == "post-intervention",]
Metadata_GLM_CD <- Metadata_GLM_CD[Metadata_GLM_CD$Diagnosis == "UC",]
Metadata_GLM_CD$Diagnosis <- NULL
Metadata_GLM_CD$timepoints <- NULL
Metadata_GLM_CD$cluster[Metadata_GLM_CD$cluster == "1"] <- "community type 1"
Metadata_GLM_CD$cluster[Metadata_GLM_CD$cluster == "2"] <- "community type 2"
nrow(Metadata_GLM_CD)

# Proportion tables
Metadata_GLM_table_UC <- table(Metadata_GLM_CD)
Metadata_GLM_table_prop_UC <-  as.data.frame(prop.table(Metadata_GLM_table_UC, margin = 1))
Metadata_GLM_table_prop_UC

ggplot(Metadata_GLM_table_prop_UC) +
  aes(x = Endoscopy.outcome, fill = cluster, weight = Freq) +
  geom_bar(position = "dodge", alpha = 0.9) +
  ylab("Community type prevalence") +
  scale_fill_manual(values = list(`community type 1` = "#4d74b8", `community type 2` = "#b87a4d")) +
  theme_classic()

# Stacked barplot
ggplot(Metadata_GLM_table_prop_UC) +
  aes(x = Endoscopy.outcome, fill = cluster, weight = Freq) +
  geom_bar() +
  ylab("Community type prevalence") +
  scale_fill_manual(values = list(`community type 1` = "#bd801c", `community type 2` = "#2b8a92")) +
  theme_classic()

# Statistics (this is immediately corrected by BH)
Metadata_GLM_table_UC
Metadata_GLM_table_prop_UC
chisq_test(Metadata_GLM_table_UC, correct = F)
cramer_v(Metadata_GLM_table_UC)
pairwise_prop_test(Metadata_GLM_table_UC, p.adjust.method = "BH", correct = FALSE)

# Community type 1
prop.test(x = c(15,14), n = c(20,31), correct = F) # **

# Community type 2
prop.test(x = c(5,17), n = c(20,31), correct = F) #  **

# R and NR patients hosting community type 2 in respectively, 54.8% and 25.0% of samples. 
# R and NR patients hosting community type 1 in respectively, 45.2% and 75.0% of samples.

# Conclusion: Endoscopic outcome is significant and has a predictive power of 3.91% (McFadden)
# Endoscopic remission is significantly higher in patients harboring community type 2, majorly driven by UC although non-significant.

# Made figures for poster and chose colors for later
# I can also chose community type 2 as main color for logistic regression figures

## Compare community type 2 remission between UC and CD
Metadata_GLM_table_UC
Metadata_GLM_table_CD
17/(14+17)*100 # remission UC
15/(15+33)*100 # remission CD
prop.test(x = c(15,17), n = c(48,31), correct = F)

## Compare community type 2 non-remission between UC and CD
Metadata_GLM_table_UC
Metadata_GLM_table_CD
5/(5+15)*100 # remission UC
10/(40+10)*100 # remission CD
prop.test(x = c(10,5), n = c(50,20), correct = F)

############  6.3.3 Univariate logistic regression: Age ############ 
Metadata <- data.frame(sample_data(phyloseq_rarefied_phages_transformed))
colnames(Metadata)
Metadata_GLM <- Metadata[,c(3,22)]
Metadata_GLM <- Metadata_GLM[Metadata_GLM$timepoints == "post-intervention",]
Metadata_GLM$timepoints <- NULL

# Add community types
Metadata_GLM <- merge(Metadata_GLM,cluster, by = 0, all = F)
rownames(Metadata_GLM) <- Metadata_GLM$Row.names   
Metadata_GLM$Row.names <- NULL 

# Convert data to detect presence community type 2 (comm 2 = 1, comm 1 = 0)
names(Metadata_GLM)[names(Metadata_GLM) == "cluster"] <- "Community_type_2"
Metadata_GLM$Community_type_2[Metadata_GLM$Community_type_2 == "1"] <- "0"
Metadata_GLM$Community_type_2[Metadata_GLM$Community_type_2 == "2"] <- "1"
Metadata_GLM$Community_type_2 <- as.numeric(Metadata_GLM$Community_type_2)
str(Metadata_GLM)

# Create multvariable logistic regression model (LRM)
glm.fit <- glm(formula = Community_type_2 ~ Age, family = binomial("logit"), data = Metadata_GLM, na.action(na.omit))
summary(glm.fit)

# Analysis of variance for individual terms
Anova(glm.fit, type="II", test="Wald")
# Conclusion: Age is not significant

# Pseudo-R-squared
pR2(glm.fit) # McFadden output*100 is best measure for R2
lrtest(glm.fit) # Significance of model

# Calculate relative risk to visualize.
RR <- exp(coef(glm.fit))[-1]
RR
exp(confint(glm.fit))[-1,]

## Simply show the significant ones out of the GLR model. Therefore in this case only show remission
glm.fit$coefficients

############  6.3.4 Multiple logistic regression ############  
# Since only endoscopy outcome came out as being singificantly associate using univariate logistic regression
# We cannot perform multiple logistic regression. I simply show here how it should be done.

# Select variables from multivariate model & appropriate samples
Metadata <- data.frame(sample_data(phyloseq_rarefied_phages_transformed))
colnames(Metadata)
Metadata_GLM <- Metadata[,c(3,11,19,22)]
Metadata_GLM <- Metadata_GLM[Metadata_GLM$timepoints == "post-intervention",]
Metadata_GLM$timepoints <- NULL
# Metadata_GLM <- Metadata_GLM[!Metadata_GLM$Location == "NA",]
#is.na(Metadata_GLM$Location) <- Metadata_GLM$Location == "NA"
is.na(Metadata_GLM$Endoscopy.outcome) <- Metadata_GLM$Endoscopy.outcome == "unknown"

# Add community types
Metadata_GLM <- merge(Metadata_GLM,cluster, by = 0, all = F)
rownames(Metadata_GLM) <- Metadata_GLM$Row.names   
Metadata_GLM$Row.names <- NULL 

# Convert data to detect presence community type 2 (comm 2 = 1, comm 1 = 0)
names(Metadata_GLM)[names(Metadata_GLM) == "cluster"] <- "Community_type_2"
Metadata_GLM$Community_type_2[Metadata_GLM$Community_type_2 == "1"] <- "0"
Metadata_GLM$Community_type_2[Metadata_GLM$Community_type_2 == "2"] <- "1"
Metadata_GLM$Community_type_2 <- as.numeric(Metadata_GLM$Community_type_2)
str(Metadata_GLM)
#View(Metadata_GLM)

# Create multvariable logistic regression model (LRM)
#install.packages("pscl")
library(pscl)
library(lmtest)

glm.fit <- glm(formula = Community_type_2 ~ Location + Endoscopy.outcome + Age, family = binomial("logit"), data = Metadata_GLM, method = "glm.fit", na.action(na.omit))
summary(glm.fit)

# Analysis of variance for individual terms
Anova(glm.fit, type="II", test="Wald")

# Pseudo-R-squared
pR2(glm.fit)*100 # McFadden output*100 is best measure for R2
lrtest(glm.fit) # Significance of model

# Calculate relative risk to visualize.
RR <- exp(coef(glm.fit))[-1]
RR
exp(confint(glm.fit))[-1,]

## Simply show the significant ones out of the GLR model. Therefore in this case only show remission
glm.fit$coefficients

## RR figure
variable <- c("Location [L2CD/UC]","Location [L3]","Endoscopic outcome [remission]", "Age")
RR_value <- c(0.660,0.618, 2.65, 1.02)
lowerCI <- c(0.259,0.228,1.27, 1.00)
upperCI <- c(1.93,1.56,5.76, 1.04)
RR_df <- data.frame(variable, RR_value,lowerCI,upperCI)
View(RR_df)

ggplot(RR_df) +
  geom_point(aes(x = variable, y = RR_value), size = 4,  color = "black") +
  geom_errorbar(aes(x = variable, y = RR_value, ymin = lowerCI, ymax = upperCI), width = 0, color = "black", ) +
  geom_text(aes(x = variable, y = RR_value, label = round(RR_value, 2)), vjust = - 1.5, size = 2.5) +
  ylim(0,6) +
  geom_hline(yintercept=1, linetype='dotted', col = 'black') +
  coord_flip() +
  xlab("") +
  ylab("Community typ M Relative risk") +
  theme_bw()

## Age 
# Create figure logistic regression
ggplot(Metadata_GLM, aes(x=Age, y=Community_type_2)) + 
  stat_smooth(method="glm", se=TRUE, method.args = list(family=binomial), col="#b87a4d", fill="#b87a4d", lty=1) +
  ylab("Community type 2 prevalence") +
  theme_bw() +
  ylim(0,1)

#p.adjust(0.2193, method = "BH", n = 3)
#p.adjust(0.00933, method = "BH", n = 3)
#p.adjust(0.0540, method = "BH", n = 3)

## Enterotype prevalence
## (1) Endoscopy.outcome: IBD
Metadata_GLM_ <- Metadata_GLM
Metadata_GLM_$Location <- NULL
Metadata_GLM_$Age <- NULL
Metadata_GLM_$Endoscopy.outcome[is.na(Metadata_GLM_$Endoscopy.outcome)] <- "unknown"
Metadata_GLM_$Community_type_2[Metadata_GLM_$Community_type_2 == "0"] <- "community type 1"
Metadata_GLM_$Community_type_2[Metadata_GLM_$Community_type_2 == "1"] <- "community type 2"

# IBD patients
# Create frequency and proportion table for IBD patients + Statistics
Metadata_GLM_table <- table(Metadata_GLM_)
Metadata_GLM_table_prop <-  as.data.frame(prop.table(Metadata_GLM_table, margin = 1))

ggplot(Metadata_GLM_table_prop) +
  aes(x = Endoscopy.outcome, fill = Community_type_2, weight = Freq) +
  geom_bar(position = "dodge", alpha=0.9) +
  ylab("Community type prevalence") +
  scale_fill_manual(values = list(`community type 1` = "#4d74b8", `community type 2` = "#b87a4d")) +
  theme_classic()

# Statistics (this is immediately corrected by BH)
# Because we now compare three groups you can do chi-squared with post-hoc proportion test
# And can also calculate cramersPhi for effect size
Metadata_GLM_table
chisq.test(Metadata_GLM_table, correct = F)
cramer_v(Metadata_GLM_table)
?cramer_v

# Community type 1
prop.test(x = c(14,30), n = c(66,72), correct = F) #  R vs NR****
prop.test(x = c(30,9), n = c(72,28), correct = F) #  R vs unknown
prop.test(x = c(14,9), n = c(66,28), correct = F) #  NR vs unknown

# Community type 2
prop.test(x = c(52,42,19), n = c(66,72,28), correct = F) #  similar as previous line
prop.test(x = c(52,42), n = c(66,72), correct = F) #  R vs NR****
prop.test(x = c(42,19), n = c(72,28), correct = F) #  R vs unknown
prop.test(x = c(52,19), n = c(66,28), correct = F) #  NR vs unknown

# Prevalence of community type 2 is higher in remission compared to non-remission

# Conclusion
# Multivariable logistic regression model explains 5.87% of community-type 2 prevalence with p-value of 0.038.
# Only endoscopic remission seems to be significant.
# The relative risk of remission is 2.37 with confidence intervals of [1.10-5.27]

values <- c("1.01", "0.460", "0.370")
variables <- c("Location", "Endoscopic outcome", "Age")
total <- c("total", "total", "total")
figure_ <- as.data.frame(cbind(values, variables, total))
figure_$values <- as.numeric(figure_$values)
esquisser(figure_)

ggplot(figure_) +
 aes(x = total, fill = variables, weight = values) +
 geom_bar() +
 scale_fill_manual(values = c(Age = "#4C3778", 
`Endoscopic outcome` = "#61914E", Location = "#786237")) +
 theme_classic()
############  6.3.5 Association between Location & endoscopic outcome: explains location is related to endoscopic outcome ######
Metadata <- data.frame(sample_data(phyloseq_rarefied_phages_transformed))
colnames(Metadata)
Metadata_GLM <- Metadata[,c(4,13,21,25)]
#View(Metadata_GLM)

# Convert data to detect presence community type 2 (comm 2 = 1, comm 1 = 0)
names(Metadata_GLM)[names(Metadata_GLM) == "cluster"] <- "Community_type_2"
Metadata_GLM$Community_type_2[Metadata_GLM$Community_type_2 == "1"] <- "Viral community type CA"
Metadata_GLM$Community_type_2[Metadata_GLM$Community_type_2 == "2"] <- "Viral community type CrM"
View(Metadata_GLM)
#str(Metadata_GLM)

ggplot(Metadata_GLM) +
 aes(x = Endoscopy.outcome, fill = Location) +
 geom_bar(position = "dodge") +
 scale_fill_hue(direction = 1) +
 theme_bw() +
  facet_grid(~Diagnosis)

ggplot(Metadata_GLM) +
  aes(x = Location, fill = Endoscopy.outcome) +
  geom_bar(position = "dodge") +
  scale_fill_hue(direction = 1) +
  theme_bw()

esquisser(Metadata_GLM)

Metadata_GLM$Community_type_2 <- NULL
#Metadata_GLM_CD <- Metadata_GLM[Metadata_GLM$Diagnosis == "CD",]
Metadata_GLM_CD <- Metadata_GLM[!Metadata_GLM$Endoscopy.outcome == "unknown",]
#Metadata_GLM_CD <- Metadata_GLM_CD[!Metadata_GLM_CD$Location == "NA",]
Metadata_GLM_table <- table(Metadata_GLM_CD$Location,Metadata_GLM_CD$Endoscopy.outcome)
prop.table(Metadata_GLM_table, margin = 1)
chisq.test(Metadata_GLM_table)
row_wise_prop_test(Metadata_GLM_table, p.adjust.method = "BH",correct = TRUE) # compare individual comparisons
row_wise_prop_test(Metadata_GLM_table, p.adjust.method = "BH",correct = TRUE, detailed = TRUE)
pairwise_prop_test(Metadata_GLM_table, p.adjust.method = "BH",correct = TRUE)

xtab <- as.table(rbind(c(22,15), c(15,14),c(15,12), c(14,31)))

dimnames(xtab) <- list(
  Class = c("L1", "L2", "L3", "unknown"),
  Clusters = c("non-remission", "remission"))

xtab

crosstable_statistics(xtab,statistics = "phi")
###################################
# 6.4 Sankey Diagram: pre-versus-post intervention
####################################
Metadata <- data.frame(sample_data(phyloseq_rarefied_phages_transformed))
colnames(Metadata)
Metadata_Sankey <- Metadata[,c(2,4,21,24,25)]
Metadata_Sankey$timepoints <- as.character(Metadata_Sankey$timepoints)
Metadata_Sankey$timepoints[Metadata_Sankey$timepoints == "pre-intervention + CD w14"] <- "pre-intervention"
table(Metadata_Sankey$Endoscopy.outcome)
Metadata_Sankey <- Metadata_Sankey[Metadata_Sankey$Endoscopy.outcome == "remission",]
Metadata_Sankey$ID <- rownames(Metadata_Sankey)

# remove rows withouth match
samples_remove <- c("B421","B84","B264","B186","B152", "B71","B242","B73","B280","B135","B80","B291","B93","B82","B202","B92","B184","B284","B400","B232","B279","B256","B356","B281","B318","B257","B351","B365","B376","B417","B430","B124","B425","B13","B382","B161","B130","B371","B251","B77","B315","B416","B49","B45","B50","B17","B163","B268","B252","B62","B68","B368","B401","B398","B427","B369","B377","B404","B443","B196","B261","B243","B343","B270","B297","B447","B241")
Metadata_Sankey <- Metadata_Sankey[!Metadata_Sankey$ID %in% samples_remove,]
Metadata_Sankey$ID <- NULL
View(Metadata_Sankey)

# Save to excel file
install.packages("writexl")
library("writexl")
setwd("/Users/daan/Desktop")
getwd()
write_xlsx(Metadata_Sankey,".Sankey.xlsx")

# Load remade excel
library("readxl")
Metadata_Sankey_remission <-  as.data.frame(read_excel("remission_comm.xlsx"))
Metadata_Sankey_remission$`Patient ID` <- NULL
Metadata_Sankey_remission_UC <- Metadata_Sankey_remission[Metadata_Sankey_remission$Diagnosis == "UC",]
Metadata_Sankey_remission_CD <- Metadata_Sankey_remission[Metadata_Sankey_remission$Diagnosis == "CD",]
#View(Metadata_Sankey_remission)

# Sankey Plot: IBD remission samples
#install.packages("ggalluvial")
library(ggalluvial)
Metadata_Sankey_remission$`Pre-intervention comm type` <- as.character(Metadata_Sankey_remission$`Pre-intervention comm type`)
Metadata_Sankey_remission$`Post-intervention comm type` <- as.character(Metadata_Sankey_remission$`Post-intervention comm type`)
Metadata_Sankey_remission$`Pre-intervention comm type`[Metadata_Sankey_remission$`Pre-intervention comm type` == "CA"] <- "comm_C_pre"
Metadata_Sankey_remission$`Post-intervention comm type`[Metadata_Sankey_remission$`Post-intervention comm type` == "CA"] <- "comm_C_post"
Metadata_Sankey_remission$`Pre-intervention comm type`[Metadata_Sankey_remission$`Pre-intervention comm type` == "CrM"] <- "comm_M_pre"
Metadata_Sankey_remission$`Post-intervention comm type`[Metadata_Sankey_remission$`Post-intervention comm type` == "CrM"] <- "comm_M_post"
#Metadata_Sankey_remission$freq <- 1
#Metadata_Sankey_remission$freq <- NULL

# Summarize and count
#Metadata_Sankey_remission$Diagnosis <- NULL
Metadata_Sankey_remission_count <- as.data.frame(table(Metadata_Sankey_remission))

View(list(Metadata_Sankey_remission))
Metadata_Sankey_remission_count_prop <- prop.table(Metadata_Sankey_remission_count, margin =2)
#Metadata_Sankey_remission_count$perc <- (Metadata_Sankey_remission_count$/sum(Metadata_Sankey_remission_count$freq))*100

Metadata_Sankey_remission_UC <- Metadata_Sankey_remission_count[Metadata_Sankey_remission_count$Diagnosis == "UC",]
Metadata_Sankey_remission_UC$Diagnosis <- NULL
Metadata_Sankey_remission_UC$perc <- (Metadata_Sankey_remission_UC$Freq/sum(Metadata_Sankey_remission_UC$Freq))*100
View(Metadata_Sankey_remission_UC)

Metadata_Sankey_remission_CD <- Metadata_Sankey_remission_count[Metadata_Sankey_remission_count$Diagnosis == "CD",]
Metadata_Sankey_remission_CD$Diagnosis <- NULL
Metadata_Sankey_remission_CD$perc <- (Metadata_Sankey_remission_CD$Freq/sum(Metadata_Sankey_remission_CD$Freq))*100
View(Metadata_Sankey_remission_CD)

ggplot(data = Metadata_Sankey_remission_count,aes(axis2 = Pre.intervention.comm.type, axis1 = Diagnosis, axis3 = Post.intervention.comm.type, y = Freq)) +
  geom_alluvium(aes(fill = Diagnosis)) +
  geom_stratum() +
  geom_text(stat = "stratum",aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Survey", "Response"),expand = c(0.15, 0.05)) +
  theme_void()

# UC
ggplot(data = Metadata_Sankey_remission_UC,aes(axis1 = Pre.intervention.comm.type, axis2 = Post.intervention.comm.type, y = Freq)) +
  geom_alluvium(aes(fill = Pre.intervention.comm.type)) +
  geom_stratum() +
  geom_text(stat = "stratum",aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Survey", "Response"),expand = c(0.15, 0.05)) +
  theme_void()

# CD
ggplot(data = Metadata_Sankey_remission_CD,aes(axis1 = Pre.intervention.comm.type, axis2 = Post.intervention.comm.type, y = Freq)) +
  geom_alluvium(aes(fill = Pre.intervention.comm.type)) +
  geom_stratum() +
  geom_text(stat = "stratum",aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Survey", "Response"),expand = c(0.15, 0.05)) +
  theme_void()

# Logistic regression
Metadata <- data.frame(sample_data(phyloseq_rarefied_phages_transformed))
colnames(Metadata)
Metadata_GLM <- Metadata[,c(2,4,21,24)]
colnames(Metadata_GLM)
Metadata_GLM <- Metadata_GLM[Metadata_GLM$timepoints == "post-intervention",]
Metadata_GLM$timepoints <- NULL

## Logistic regression IBD
## You can control for a confounder by adding a variable to the equation, in our case 'patient_ID'. By doing this you can see that P-value drops a bit for endoscopic timepoints but 
## Here we wanted to evaluate if patients going into endoscopic remission will have a shift in community type from baseline to endpoint.
## We can evaluate this by modeling logistic regression on community type M carrier status, confounding for patient_ID, associated with timepoints
## If there is a significant shift from community type C towards community type M at endpoint, we should see a signifiant association.
library(pscl)
Logistic_remission <-  as.data.frame(read_excel("remission_comm_all.xlsx"))
names(Logistic_remission)[names(Logistic_remission) == "cluster"] <- "Community_type_2"
Logistic_remission$Community_type_2[Logistic_remission$Community_type_2 == "1"] <- "0"
Logistic_remission$Community_type_2[Logistic_remission$Community_type_2 == "2"] <- "1"
Logistic_remission$Community_type_2 <- as.numeric(Logistic_remission$Community_type_2)
str(Logistic_remission)
View(Logistic_remission)

glm.fit <- glm(Community_type_2 ~ timepoints, family = binomial("logit"), data = Logistic_remission)
summary(glm.fit)
glm.fit <- glm(Community_type_2 ~ timepoints+patient_ID, family = binomial("logit"), data = Logistic_remission)
summary(glm.fit)

# Analysis of variance for individual terms
Anova(glm.fit, type="II", test="Wald")

# Pseudo-R-squared
pR2(glm.fit) # McFadden output*100 is best measure for R2
lrtest(glm.fit) # Significance of model

## Logistic regression UC
Logistic_remission_UC <- Logistic_remission[Logistic_remission$Diagnosis == "UC",]
glm.fit <- glm(Community_type_2 ~ timepoints, family = binomial("logit"), data = Logistic_remission_UC)
summary(glm.fit)
glm.fit <- glm(Community_type_2 ~ timepoints+patient_ID, family = binomial("logit"), data = Logistic_remission_UC)
summary(glm.fit)
Anova(glm.fit, type="II", test="Wald")

## Logistic regression CD
Logistic_remission_CD <- Logistic_remission[Logistic_remission$Diagnosis == "CD",]
glm.fit <- glm(Community_type_2 ~ timepoints, family = binomial("logit"), data = Logistic_remission_CD)
summary(glm.fit)
glm.fit <- glm(Community_type_2 ~ timepoints+patient_ID, family = binomial("logit"), data = Logistic_remission_CD)
summary(glm.fit) 
Anova(glm.fit, type="II", test="Wald")
####################################
# 6.3 Visualization of virome composition fit for adding metadata 
####################################
PcOA_plot_ %>%
  ggplot() +
  geom_point(aes(x = Axis.1, y= Axis.2, color = cluster, shape = Endoscopy.outcome), alpha = 1) + # size=2, alpha=1
  scale_shape_manual(values = c(21,16, 2)) +
  theme_bw() +
  scale_color_manual(values = c("#990011FF","#00203FFF")) +
  geom_hline(yintercept = 0, linetype = 3, colour = "black") + 
  geom_vline(xintercept = 0, linetype = 3, colour = "black")  +
  xlab("PC1 (8.6%)") +
  ylab("PC2 (4.1%)")
###################################

######### !!!!!!!!!!!!!!! #########
# SAVE Global Environment #