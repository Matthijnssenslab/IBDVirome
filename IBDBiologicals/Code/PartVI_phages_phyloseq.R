####################################
# SCRIPT 7: PROKARYOTIC VIRUSES: PHYLOSEQ
####################################
# Before starting this script
####################################
# ------------------> Load the "Global environment" output of script 3 "Rarefaction.R"
# ------------------> You can show presence/absence table without the need for rarefaction based on "Mastertable"
####################################
# 0. Packages: Install packages and load them wherever needed
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
library(devtools) # before downlading from github you need to load devtools
#install.packages("BiocManager")
library(BiocManager)
#BiocManager::install("ComplexHeatmap")
#install_github("jokergoo/ComplexHeatmap")
library(ComplexHeatmap)
#source("https://raw.githubusercontent.com/joey711/phyloseq/master/inst/scripts/installer.R",local = TRUE) ## like this it works!!
library(phyloseq)
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
#install.packages("ggpubr")
#install.packages("ggforce")
#install.packages("rJava")
#install.packages("UpSetR")
#install.packages("tidyverse")
#install.packages("venneuler")
#install.packages("grid")
#install.packages("phyloseq")
#install.packages("coin")
#install.packages("rstatix")
library(rJava)
library(phyloseq)
library(tidyverse)
library(venneuler)
library(grid)
library(ggplot2)
library(ggforce)
library(ggpubr)
library(dplyr)
library(esquisse)
library(modeldata)
library(tidyverse)
library(ggThemeAssist)
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
library(ape)
library(taxonomizr)
library(readxl)
library(BiocManager)
library(gridExtra)
library(RColorBrewer)
library(coin)
library(rstatix)
library(robustrank)
####################################
# 1. Recapitulation of mastertables
####################################
## Take into account that most of diversity analysis already correct for sequencing depth (alpha-diversity), so you have to feed it the unrarefied completeness data

# global mastertable
Mastertable
Mastertable_unrarefied
Mastertable_rarefied

# viral mastertable
# unrarefied
Mastertable_viral_unrarefied  # PRESENCE/ABSENCE, ALPHA-DIVERSITY (Shannon diversity)

# rarefied
Mastertable_viral_rarefied # RELATIVE ABUNDANCES, RICHNESS, EVENESS 
####################################
# 2. subset phages
####################################
Mastertable_viral_unrarefied_phages <- Mastertable_viral_unrarefied[Mastertable_viral_unrarefied$Final_viral == "phage",]
Mastertable_viral_rarefied_phages <- Mastertable_viral_rarefied[Mastertable_viral_rarefied$Final_viral == "phage",]
####################################
# Rarefied
####################################
# 3. Input tables 
####################################
# Raw
Mastertable_viral_unrarefied_phages # unrarefied 
Mastertable_viral_rarefied_phages # rarefied
####################################
# 4. Create phyloseq tables (abundances table, taxonomy table and sample table)
####################################
# 4.1 unrarefied data 
####################################
# check: https://vaulot.github.io/tutorials/Phyloseq_tutorial.html
# In our case we have a few additional nr-contigs in the unrarefied fraction. Hold this in mind
####################################
# 4.1.1 Create phyloseq abundance table
####################################
vector_1 <- (which(names(Mastertable_viral_unrarefied_phages)== "Virsorter")-1) 
vector_2 <- (which(names(Mastertable_viral_unrarefied_phages)== "Final_superkingdom"))
vector_3 <- which(names(Mastertable_viral_unrarefied_phages)== "Final_species")

abundance_table_unrarefied <- Mastertable_viral_unrarefied_phages[,c(1:vector_1)]
####################################
# 4.1.2 Create phyloseq taxonomy table
####################################
taxonomy_table_unrarefied <- Mastertable_viral_unrarefied_phages[,c(vector_2:vector_3)]
names(taxonomy_table_unrarefied)[names(taxonomy_table_unrarefied) == "Final_superkingdom"] <- "Superkingdom"
names(taxonomy_table_unrarefied)[names(taxonomy_table_unrarefied) == "Final_phylum"] <- "Phylum"
names(taxonomy_table_unrarefied)[names(taxonomy_table_unrarefied) == "Final_class"] <- "Class"
names(taxonomy_table_unrarefied)[names(taxonomy_table_unrarefied) == "Final_order"] <- "Order"
names(taxonomy_table_unrarefied)[names(taxonomy_table_unrarefied) == "Final_family"] <- "Family"
names(taxonomy_table_unrarefied)[names(taxonomy_table_unrarefied) == "Final_subfamily"] <- "Subfamily"
names(taxonomy_table_unrarefied)[names(taxonomy_table_unrarefied) == "Final_genus"] <- "Genus"
names(taxonomy_table_unrarefied)[names(taxonomy_table_unrarefied) == "Final_species"] <- "Species"
####################################
# 4.1.3 Create phyloseq sample tables (contain samples + metadata)
####################################
setwd("/Users/daan/Desktop/Bioinformatics/pre-processing/Biological_study/Metadata/output_metadata")
getwd()
dir()

metadata <- read_excel("Final_Metadata_subselection_All.xlsx")
rownames(metadata) <- metadata$`Randomized Nr.`
names(metadata)[1]<-paste("sample")

sample_table_unrarefied_for_all <- (merge(t(abundance_table_unrarefied),metadata, by=0, all=T)) # this is done so the metadata matches the samples
rownames(sample_table_unrarefied_for_all) <- sample_table_unrarefied_for_all$Row.names
sample_table_unrarefied_for_all$Row.names <- NULL

vector_4 <- (which(names(sample_table_unrarefied_for_all)== "sample"))
vector_5 <- which(names(sample_table_unrarefied_for_all)== "Biomarker_outcome_combined_2 (NR/R/RE)")

sample_table_unrarefied_for_all <- sample_table_unrarefied_for_all[,c(vector_4:vector_5)]
####################################
## 4.1.4 Overview all unrarefied tables (abundance, taxonomy and sample tables)
####################################
## ABUNDANCE TABLE
View(abundance_table_unrarefied)

## TAXONOMY TABLE 
View(taxonomy_table_unrarefied)

## SAMPLE TABLE
View(sample_table_unrarefied_for_all)
####################################
# 4.2 rarefied data 
####################################
# check: https://vaulot.github.io/tutorials/Phyloseq_tutorial.html
####################################
# 4.2.1 Create phyloseq abundance table
####################################
vector_1 <- (which(names(Mastertable_viral_rarefied_phages)== "Virsorter")-1) 
vector_2 <- (which(names(Mastertable_viral_rarefied_phages)== "Final_superkingdom"))
vector_3 <- which(names(Mastertable_viral_rarefied_phages)== "Final_species")

abundance_table_rarefied <- Mastertable_viral_rarefied_phages[,c(1:vector_1)]
####################################
# 4.2.2 Create phyloseq taxonomy table
####################################
taxonomy_table_rarefied <- Mastertable_viral_rarefied_phages[,c(vector_2:vector_3)]
names(taxonomy_table_rarefied)[names(taxonomy_table_rarefied) == "Final_superkingdom"] <- "Superkingdom"
names(taxonomy_table_rarefied)[names(taxonomy_table_rarefied) == "Final_phylum"] <- "Phylum"
names(taxonomy_table_rarefied)[names(taxonomy_table_rarefied) == "Final_class"] <- "Class"
names(taxonomy_table_rarefied)[names(taxonomy_table_rarefied) == "Final_order"] <- "Order"
names(taxonomy_table_rarefied)[names(taxonomy_table_rarefied) == "Final_family"] <- "Family"
names(taxonomy_table_rarefied)[names(taxonomy_table_rarefied) == "Final_subfamily"] <- "Subfamily"
names(taxonomy_table_rarefied)[names(taxonomy_table_rarefied) == "Final_genus"] <- "Genus"
names(taxonomy_table_rarefied)[names(taxonomy_table_rarefied) == "Final_species"] <- "Species"
####################################
# 4.2.3 Create phyloseq sample tables (contain samples + metadata)
####################################
setwd("/Users/daan/Desktop/Bioinformatics/pre-processing/Biological_study/Metadata/output_metadata")
getwd()
dir()

metadata <- read_excel("Final_Metadata_subselection_All.xlsx")
rownames(metadata) <- metadata$`Randomized Nr.`
names(metadata)[1]<-paste("sample")

sample_table_rarefied_for_all <- (merge(t(abundance_table_rarefied),metadata, by=0, all=T)) # this is done so the metadata matches the samples
rownames(sample_table_rarefied_for_all) <- sample_table_rarefied_for_all$Row.names
sample_table_rarefied_for_all$Row.names <- NULL

vector_4 <- (which(names(sample_table_rarefied_for_all)== "sample"))
vector_5 <- which(names(sample_table_rarefied_for_all)== "Biomarker_outcome_combined_2 (NR/R/RE)")

sample_table_rarefied_for_all <- sample_table_rarefied_for_all[,c(vector_4:vector_5)]
####################################
## 4.2.4 Overview all unrarefied tables (abundance, taxonomy and sample tables)
####################################
## ABUNDANCE TABLE
View(abundance_table_rarefied)

## TAXONOMY TABLE 
View(taxonomy_table_rarefied)

## SAMPLE TABLE
View(sample_table_rarefied_for_all)
####################################
# 5. Create phyloSeq objects
####################################
# 5.1 Unrarefied data
####################################
# 5.1.1 Create a matrix of abundance and taxonomy tables
####################################
abundance_table_unrarefied_m <- as.matrix(abundance_table_unrarefied)
taxonomy_table_unrarefied_m <- as.matrix(taxonomy_table_unrarefied)
####################################
# 5.1.2 Transform to phylseq objects
####################################
ABUNDANCE_unrarefied_family <- otu_table(abundance_table_unrarefied_m, taxa_are_rows = TRUE)
TAX_unrarefied_family <- tax_table(taxonomy_table_unrarefied_m)
samples <- sample_data(sample_table_unrarefied_for_all)
phyloseq_unrarefied_family <- phyloseq(ABUNDANCE_unrarefied_family,TAX_unrarefied_family, samples)
####################################
# 5.1.3 Visualize data & subset only the phage fraction
####################################
sample_names(phyloseq_unrarefied_family) # All sample names
rank_names(phyloseq_unrarefied_family) # All taxonomies
sample_variables(phyloseq_unrarefied_family) ## All metadata
# if you want to subset the dataframe into only "responders", just do "responder_data <- subset_samples(phyloseq_unrarefied_family, response == "responder")"
phyloseq_unrarefied_phages <- subset_taxa(phyloseq_unrarefied_family, Superkingdom == "Viruses")
####################################
# 5.2 rarefied data
####################################
# 5.2.1 Create a matrix of abundance and taxonomy tables
####################################
abundance_table_rarefied_m <- as.matrix(abundance_table_rarefied)
taxonomy_table_rarefied_m <- as.matrix(taxonomy_table_rarefied)
####################################
# 5.2.2 Transform to phylseq objects
####################################
ABUNDANCE_rarefied_family <- otu_table(abundance_table_rarefied_m, taxa_are_rows = TRUE)
TAX_rarefied_family <- tax_table(taxonomy_table_rarefied_m)
samples <- sample_data(sample_table_rarefied_for_all)
phyloseq_rarefied_family <- phyloseq(ABUNDANCE_rarefied_family,TAX_rarefied_family, samples)
####################################
# 5.2.3 Visualize data & subset only the phage fraction
####################################
sample_names(phyloseq_rarefied_family) # All sample names
rank_names(phyloseq_rarefied_family) # All taxonomies
sample_variables(phyloseq_rarefied_family) ## All metadata
# if you want to subset the dataframe into only "responders", just do "responder_data <- subset_samples(phyloseq_unrarefied_family, response == "responder")"
phyloseq_rarefied_phages <- subset_taxa(phyloseq_rarefied_family, Superkingdom == "Viruses")
####################################
# 5.3 PhyloSeq format
####################################
## Unrarefied
phyloseq_unrarefied_phages

## Rarefied
phyloseq_rarefied_phages
####################################
# 6. Relative abundances (use rarefied data): convert to normal DF
####################################
# 6.1 Create dataframe (Tax, OTU)
####################################
phyloseq_rarefied_phages_tax <- as.data.frame(tax_table(phyloseq_rarefied_phages))
phyloseq_rarefied_phages_OTU <- as.data.frame(otu_table(phyloseq_rarefied_phages))
####################################
# 6.2 Merge 
####################################
phyloseq_rarefied_phages_ggplot <- merge(phyloseq_rarefied_phages_OTU,phyloseq_rarefied_phages_tax, by = 0, all = TRUE)
rownames(phyloseq_rarefied_phages_ggplot) <- phyloseq_rarefied_phages_ggplot$Row.names
phyloseq_rarefied_phages_ggplot$Row.names <- NULL
####################################
# 6.3 Aggregate
####################################
vector_1 <- (which(names(phyloseq_rarefied_phages_ggplot)== "Kingdom")-1) 
ggplot_names <- colnames(phyloseq_rarefied_phages_ggplot[,c(1:vector_1)])

phyloseq_rarefied_phages_ggplot_class_agg <- aggregate(. ~Class, FUN = sum, data = phyloseq_rarefied_phages_ggplot[,colnames(phyloseq_rarefied_phages_ggplot) %in% ggplot_names | colnames(phyloseq_rarefied_phages_ggplot) == "Class"])
rownames(phyloseq_rarefied_phages_ggplot_class_agg) <- phyloseq_rarefied_phages_ggplot_class_agg$Class
phyloseq_rarefied_phages_ggplot_class_agg$Class <- NULL
####################################
# 6.4 Calculate RA
####################################
phyloseq_rarefied_phages_ggplot_class <- sweep(phyloseq_rarefied_phages_ggplot_class_agg, 2, colSums(phyloseq_rarefied_phages_ggplot_class_agg), '/')
phyloseq_rarefied_phages_ggplot_class[is.na(phyloseq_rarefied_phages_ggplot_class)] <- 0
####################################
# 6.5 Combine with metadata
####################################
phyloseq_rarefied_phages_ggplot_class <- merge(t(phyloseq_rarefied_phages_ggplot_class), sample_table_rarefied_for_all, by = 0, all = F)
rownames(phyloseq_rarefied_phages_ggplot_class) <- phyloseq_rarefied_phages_ggplot_class$Row.names
phyloseq_rarefied_phages_ggplot_class$Row.names <- NULL
####################################
# 6.6 Melt
####################################
vector_1 <- (which(names(phyloseq_rarefied_phages_ggplot_class)== "Unannotated")) 
vector_2 <- (which(names(phyloseq_rarefied_phages_ggplot_class)== "Diagnosis")) 
vector_3 <- (which(names(phyloseq_rarefied_phages_ggplot_class)== "Gender")) 
vector_4 <- (which(names(phyloseq_rarefied_phages_ggplot_class)== "timepoints")) 
vector_5 <- (which(names(phyloseq_rarefied_phages_ggplot_class)== "Endoscopy-outcome")) 
vector_6 <- (which(names(phyloseq_rarefied_phages_ggplot_class)é== "Therapy")) 

phyloseq_rarefied_phages_ggplot_no_crass <- phyloseq_rarefied_phages_ggplot_class[,c(1:vector_1,vector_2,vector_3,vector_5,vector_6,vector_4)]
View(phyloseq_rarefied_phages_ggplot_no_crass)

phyloseq_rarefied_phages_ggplot_no_crass$`Endoscopy-outcome`[phyloseq_rarefied_phages_ggplot_no_crass$`Endoscopy-outcome` == 1] <- "remission"
phyloseq_rarefied_phages_ggplot_no_crass$`Endoscopy-outcome`[phyloseq_rarefied_phages_ggplot_no_crass$`Endoscopy-outcome` == 0] <- "non-remission"
phyloseq_rarefied_phages_ggplot_no_crass$`Endoscopy-outcome`[phyloseq_rarefied_phages_ggplot_no_crass$`Endoscopy-outcome` == "NA"] <- "unknown"
phyloseq_rarefied_phages_ggplot_class <- melt(phyloseq_rarefied_phages_ggplot_no_crass, id.vars = c("Diagnosis", "Gender", "timepoints", "Endoscopy-outcome", "Therapy"))
####################################
# 6.7 Relative abundance boxplots
####################################
# 6.7.1 UC
####################################
phyloseq_rarefied_phages_ggplot_class_UC <- phyloseq_rarefied_phages_ggplot_class[phyloseq_rarefied_phages_ggplot_class$Diagnosis == "UC",]
phyloseq_rarefied_phages_ggplot_class_UC$timepoints[phyloseq_rarefied_phages_ggplot_class_UC$timepoints == "w0"] <- "baseline"
phyloseq_rarefied_phages_ggplot_class_UC$timepoints[phyloseq_rarefied_phages_ggplot_class_UC$timepoints == "w14"] <- "primary endpoint"
#View(phyloseq_rarefied_phages_ggplot_class_UC)

# figure:
phyloseq_rarefied_phages_ggplot_class_UC %>%
  filter(variable %in% c("Caudoviricetes (CrAss)", "Caudoviricetes (non-CrAss)", "Malgrandaviricetes", "phages")) %>%
  filter(!`Endoscopy-outcome` %in% "unknown") %>%
   ggplot() +
  aes(x = variable, y = value, fill = `Endoscopy-outcome`) +
  geom_boxplot() +
  scale_fill_manual(values = c(`non-remission` = "#AD4B38", remission = "#49AD86")) +
  theme_bw() +
  xlab("") +
  ylab("Relative abundance") +
  theme(legend.position = "none")  

# print long way 
# use different statistic within & between RA.
####################################
# 6.7.2 Statistics within timepoints
####################################
# Statistics within baseline
# Caudo (non-crAss)
table(phyloseq_rarefied_phages_ggplot_class_UC$variable)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- phyloseq_rarefied_phages_ggplot_class_UC[c(phyloseq_rarefied_phages_ggplot_class_UC$Week == "baseline" & phyloseq_rarefied_phages_ggplot_class_UC$Diagnosis == "UC" & phyloseq_rarefied_phages_ggplot_class_UC$variable == "Caudoviricetes (non-CrAss)"),]
table(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$`Endoscopic_outcome_combined (NR/R)`) ## 22 R patients/ 22 NR patients
wilcox.test(value~`Endoscopic_outcome_combined (NR/R)`, data = total_virome_R_NR_total_phages_w0_UC_caudo_non_crass, alternative = "two.sided",  paired = FALSE, exact = F)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$remission <- total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$`Endoscopic_outcome_combined (NR/R)`
effect_size <- wilcox_effsize(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass,value~remission)
(effect_size$effsize)^2
effect_size$n1
effect_size$n2
# P-value = 0.7361
# R-squared: 0.002817431
#  n1: 22
#  n2: 22

# Caudo (crAss)
table(phyloseq_rarefied_phages_ggplot_class_UC$variable)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- phyloseq_rarefied_phages_ggplot_class_UC[c(phyloseq_rarefied_phages_ggplot_class_UC$Week == "baseline" & phyloseq_rarefied_phages_ggplot_class_UC$Diagnosis == "UC" & phyloseq_rarefied_phages_ggplot_class_UC$variable == "Caudoviricetes (CrAss)"),]
table(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$`Endoscopic_outcome_combined (NR/R)`) ## 22 R patients/ 22 NR patients
wilcox.test(value~`Endoscopic_outcome_combined (NR/R)`, data = total_virome_R_NR_total_phages_w0_UC_caudo_non_crass, alternative = "two.sided",  paired = FALSE, exact = F)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$remission <- total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$`Endoscopic_outcome_combined (NR/R)`
effect_size <- wilcox_effsize(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass,value~remission, alternative = "two.sided", paired = FALSE)
(effect_size$effsize)^2 
effect_size$n1
effect_size$n2
# P-value = 0.2818
# R-squared: 0.02691363
#  n1: 22
#  n2: 22

# Caudo (Micro)
table(phyloseq_rarefied_phages_ggplot_class_UC$variable)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- phyloseq_rarefied_phages_ggplot_class_UC[c(phyloseq_rarefied_phages_ggplot_class_UC$Week == "baseline" & phyloseq_rarefied_phages_ggplot_class_UC$Diagnosis == "UC" & phyloseq_rarefied_phages_ggplot_class_UC$variable == "Malgrandaviricetes"),]
table(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$`Endoscopic_outcome_combined (NR/R)`) ## 22 R patients/ 22 NR patients
wilcox.test(value~`Endoscopic_outcome_combined (NR/R)`, data = total_virome_R_NR_total_phages_w0_UC_caudo_non_crass, alternative = "two.sided",  paired = FALSE, exact = F)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$remission <- total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$`Endoscopic_outcome_combined (NR/R)`
effect_size <- wilcox_effsize(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass,value~remission, alternative = "two.sided", paired = FALSE)
(effect_size$effsize)^2 
effect_size$n1
effect_size$n2
# P-value = 0.6899
# R-squared: 0.003835106
#  n1: 22
#  n2: 22

# Statistics within primary endpoint
# Caudo (non-crAss)
table(phyloseq_rarefied_phages_ggplot_class_UC$variable)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- phyloseq_rarefied_phages_ggplot_class_UC[c(phyloseq_rarefied_phages_ggplot_class_UC$Week == "primary endpoint" & phyloseq_rarefied_phages_ggplot_class_UC$Diagnosis == "UC" & phyloseq_rarefied_phages_ggplot_class_UC$variable == "Caudoviricetes (non-CrAss)"),]
table(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$`Endoscopic_outcome_combined (NR/R)`) ## 31 R patients/ 20 NR patients
wilcox.test(value~`Endoscopic_outcome_combined (NR/R)`, data = total_virome_R_NR_total_phages_w0_UC_caudo_non_crass, alternative = "two.sided",  paired = FALSE, exact = F)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$remission <- total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$`Endoscopic_outcome_combined (NR/R)`
effect_size <- wilcox_effsize(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass,value~remission, alternative = "two.sided", paired = FALSE)
unique(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$remission)
(effect_size$effsize)^2
effect_size$n1
effect_size$n2
# P-value = 0.03147*****
# R-squared: 0.09154868
#  n1 (non-remission): 20
#  n2 (remission): 31

# Caudo (crAss)
table(phyloseq_rarefied_phages_ggplot_class_UC$variable)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- phyloseq_rarefied_phages_ggplot_class_UC[c(phyloseq_rarefied_phages_ggplot_class_UC$Week == "primary endpoint" & phyloseq_rarefied_phages_ggplot_class_UC$Diagnosis == "UC" & phyloseq_rarefied_phages_ggplot_class_UC$variable == "Caudoviricetes (CrAss)"),]
table(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$`Endoscopic_outcome_combined (NR/R)`)  ## 31 R patients/ 20 NR patients
wilcox.test(value~`Endoscopic_outcome_combined (NR/R)`, data = total_virome_R_NR_total_phages_w0_UC_caudo_non_crass, alternative = "two.sided",  paired = FALSE, exact = F)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$remission <- total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$`Endoscopic_outcome_combined (NR/R)`
effect_size <- wilcox_effsize(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass,value~remission, alternative = "two.sided", paired = FALSE)
unique(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$remission)
(effect_size$effsize)^2 
effect_size$n1
effect_size$n2
# P-value = 0.1112
# R-squared: 0.05034894
#  n1 (non-remission): 20
#  n2 (remission): 31

# Caudo (Micro)
table(phyloseq_rarefied_phages_ggplot_class_UC$variable)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- phyloseq_rarefied_phages_ggplot_class_UC[c(phyloseq_rarefied_phages_ggplot_class_UC$Week == "primary endpoint" & phyloseq_rarefied_phages_ggplot_class_UC$Diagnosis == "UC" & phyloseq_rarefied_phages_ggplot_class_UC$variable == "Malgrandaviricetes"),]
table(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$`Endoscopic_outcome_combined (NR/R)`) ## 22 R patients/ 22 NR patients
wilcox.test(value~`Endoscopic_outcome_combined (NR/R)`, data = total_virome_R_NR_total_phages_w0_UC_caudo_non_crass, alternative = "two.sided",  paired = FALSE, exact = F)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$remission <- total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$`Endoscopic_outcome_combined (NR/R)`
effect_size <- wilcox_effsize(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass,value~remission, alternative = "two.sided", paired = FALSE)
unique(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$remission)
(effect_size$effsize)^2 
effect_size$n1
effect_size$n2
# P-value = 0.2063
# R-squared: 0.03179098
#  n1 (non-remission): 20
#  n2 (remission): 31
####################################
# 6.7.3 Statistics between timepoints
####################################
# Statistics between baseline and primary endpoint
# Statistics within Remission

# caudoviricetes (non-CrAss)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- phyloseq_rarefied_phages_ggplot_class_UC[c(phyloseq_rarefied_phages_ggplot_class_UC$`Endoscopic_outcome_combined (NR/R)` == "remission" & phyloseq_rarefied_phages_ggplot_class_UC$Diagnosis == "UC" & phyloseq_rarefied_phages_ggplot_class_UC$variable == "Caudoviricetes (non-CrAss)"),]
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired <- total_virome_R_NR_total_phages_w0_UC_caudo_non_crass
#View(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired)

colnames_de <- colnames(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired)
de1 <-data.frame("3678",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("5520",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("7903",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("11057",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("11133",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("12255",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("12494",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("13419",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("10213",NA, NA,"baseline",NA,NA,NA,NA)

#View(de1)
names(de1)<-c(colnames_de)
unpaired_paired_test <- rbind(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired, de1)
unpaired_paired_test_baseline <- unpaired_paired_test[unpaired_paired_test$Week == "baseline",]
unpaired_paired_test_baseline <- unpaired_paired_test_baseline[order(unpaired_paired_test_baseline$patient_ID) , ]
unpaired_paired_test_baseline$value <- as.numeric(unpaired_paired_test_baseline$value)
unpaired_paired_test_PE <- unpaired_paired_test[!unpaired_paired_test$Week == "baseline",]
unpaired_paired_test_PE <- unpaired_paired_test_PE[order(unpaired_paired_test_PE$patient_ID) , ]
unpaired_paired_test_PE$value <- as.numeric(unpaired_paired_test_PE$value)

#View(unpaired_paired_test)
#View(unpaired_paired_test_baseline)
#View(unpaired_paired_test_PE)

pm.wilcox.test(unpaired_paired_test_PE$value,unpaired_paired_test_baseline$value,alternative = "two.sided", correct = F, method = "MW-MW")
?pm.wilcox.test
# P-value = 0.5359 (x: baseline, y: PE)
# Xpaired: 22
# Ypaired: 22
# Xextra: 0
# Y-extra: 9

# Caudo (CrAss)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- phyloseq_rarefied_phages_ggplot_class_UC[c(phyloseq_rarefied_phages_ggplot_class_UC$`Endoscopic_outcome_combined (NR/R)` == "remission" & phyloseq_rarefied_phages_ggplot_class_UC$Diagnosis == "UC" & phyloseq_rarefied_phages_ggplot_class_UC$variable == "Caudoviricetes (CrAss)"),]
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired <- total_virome_R_NR_total_phages_w0_UC_caudo_non_crass
#View(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired)

colnames_de <- colnames(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired)
de1 <-data.frame("3678",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("5520",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("7903",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("11057",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("11133",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("12255",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("12494",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("13419",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("10213",NA, NA,"baseline",NA,NA,NA,NA)

#View(de1)
names(de1)<-c(colnames_de)
unpaired_paired_test <- rbind(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired, de1)
unpaired_paired_test_baseline <- unpaired_paired_test[unpaired_paired_test$Week == "baseline",]
unpaired_paired_test_baseline <- unpaired_paired_test_baseline[order(unpaired_paired_test_baseline$patient_ID) , ]
unpaired_paired_test_baseline$value <- as.numeric(unpaired_paired_test_baseline$value)
unpaired_paired_test_PE <- unpaired_paired_test[!unpaired_paired_test$Week == "baseline",]
unpaired_paired_test_PE <- unpaired_paired_test_PE[order(unpaired_paired_test_PE$patient_ID) , ]
unpaired_paired_test_PE$value <- as.numeric(unpaired_paired_test_PE$value)

#View(unpaired_paired_test)
#View(unpaired_paired_test_baseline)
#View(unpaired_paired_test_PE)

pm.wilcox.test(unpaired_paired_test_PE$value,unpaired_paired_test_baseline$value,alternative = "two.sided", correct = F, method = "MW-MW")

# P-value = 0.7781 (x: baseline, y: PE)
# Xpaired: 22
# Ypaired: 22
# Xextra: 0
# Y-extra: 9

# Caudo (Micro)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- phyloseq_rarefied_phages_ggplot_class_UC[c(phyloseq_rarefied_phages_ggplot_class_UC$`Endoscopic_outcome_combined (NR/R)` == "remission" & phyloseq_rarefied_phages_ggplot_class_UC$Diagnosis == "UC" & phyloseq_rarefied_phages_ggplot_class_UC$variable == "Malgrandaviricetes"),]
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired <- total_virome_R_NR_total_phages_w0_UC_caudo_non_crass
#View(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired)

colnames_de <- colnames(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired)
de1 <-data.frame("3678",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("5520",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("7903",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("11057",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("11133",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("12255",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("12494",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("13419",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("10213",NA, NA,"baseline",NA,NA,NA,NA)

#View(de1)
names(de1)<-c(colnames_de)
unpaired_paired_test <- rbind(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired, de1)
unpaired_paired_test_baseline <- unpaired_paired_test[unpaired_paired_test$Week == "baseline",]
unpaired_paired_test_baseline <- unpaired_paired_test_baseline[order(unpaired_paired_test_baseline$patient_ID) , ]
unpaired_paired_test_baseline$value <- as.numeric(unpaired_paired_test_baseline$value)
unpaired_paired_test_PE <- unpaired_paired_test[!unpaired_paired_test$Week == "baseline",]
unpaired_paired_test_PE <- unpaired_paired_test_PE[order(unpaired_paired_test_PE$patient_ID) , ]
unpaired_paired_test_PE$value <- as.numeric(unpaired_paired_test_PE$value)

#View(unpaired_paired_test)
#View(unpaired_paired_test_baseline)
#View(unpaired_paired_test_PE)

pm.wilcox.test(unpaired_paired_test_PE$value,unpaired_paired_test_baseline$value,alternative = "two.sided", correct = F, method = "MW-MW")

# P-value = 0.9298 (x: baseline, y: PE)
# Xpaired: 22
# Ypaired: 22
# Xextra: 0
# Y-extra: 9

# Statistics within non-Remission
# caudoviricetes (non-CrAss)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- phyloseq_rarefied_phages_ggplot_class_UC[c(phyloseq_rarefied_phages_ggplot_class_UC$`Endoscopic_outcome_combined (NR/R)` == "non-remission" & phyloseq_rarefied_phages_ggplot_class_UC$Diagnosis == "UC" & phyloseq_rarefied_phages_ggplot_class_UC$variable == "Caudoviricetes (non-CrAss)"),]
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired <- total_virome_R_NR_total_phages_w0_UC_caudo_non_crass
#View(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired)

colnames_de <- colnames(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired)
de1 <-data.frame("7164",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("7270",NA, NA,"primary endpoint",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("10449",NA, NA,"primary endpoint",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("12515",NA, NA,"primary endpoint",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("13453",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("9516",NA, NA,"primary endpoint",NA,NA,NA,NA)

#View(de1)
names(de1)<-c(colnames_de)
unpaired_paired_test <- rbind(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired, de1)
unpaired_paired_test$patient_ID[unpaired_paired_test$patient_ID == "12770" & unpaired_paired_test$Therapy_2 == "TNF"] <- "1"
unpaired_paired_test$patient_ID[unpaired_paired_test$patient_ID == "13093" & unpaired_paired_test$Therapy_2 == "TNF"] <- "2"
unpaired_paired_test$patient_ID[unpaired_paired_test$patient_ID == "9516" & unpaired_paired_test$Therapy_2 == "TNF"] <- "3"

unpaired_paired_test_baseline <- unpaired_paired_test[unpaired_paired_test$Week == "baseline",]
unpaired_paired_test_baseline <- unpaired_paired_test_baseline[order(unpaired_paired_test_baseline$patient_ID) , ]
unpaired_paired_test_baseline$value <- as.numeric(unpaired_paired_test_baseline$value)
unpaired_paired_test_PE <- unpaired_paired_test[!unpaired_paired_test$Week == "baseline",]
unpaired_paired_test_PE <- unpaired_paired_test_PE[order(unpaired_paired_test_PE$patient_ID) , ]
unpaired_paired_test_PE$value <- as.numeric(unpaired_paired_test_PE$value)

#View(unpaired_paired_test)
#View(unpaired_paired_test_baseline)
#View(unpaired_paired_test_PE)

pm.wilcox.test(unpaired_paired_test_PE$value,unpaired_paired_test_baseline$value,alternative = "two.sided", correct = F, method = "MW-MW")

# P-value = 0.008545***** (x: baseline, y: PE)
# Xpaired: 18
# Ypaired: 18
# Xextra: 2
# Y-extra: 4

# Caudo (CrAss)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- phyloseq_rarefied_phages_ggplot_class_UC[c(phyloseq_rarefied_phages_ggplot_class_UC$`Endoscopic_outcome_combined (NR/R)` == "non-remission" & phyloseq_rarefied_phages_ggplot_class_UC$Diagnosis == "UC" & phyloseq_rarefied_phages_ggplot_class_UC$variable == "Caudoviricetes (CrAss)"),]
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired <- total_virome_R_NR_total_phages_w0_UC_caudo_non_crass
#View(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired)

colnames_de <- colnames(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired)
de1 <-data.frame("7164",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("7270",NA, NA,"primary endpoint",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("10449",NA, NA,"primary endpoint",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("12515",NA, NA,"primary endpoint",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("13453",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("9516",NA, NA,"primary endpoint",NA,NA,NA,NA)

#View(de1)
names(de1)<-c(colnames_de)
unpaired_paired_test <- rbind(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired, de1)
unpaired_paired_test$patient_ID[unpaired_paired_test$patient_ID == "12770" & unpaired_paired_test$Therapy_2 == "TNF"] <- "1"
unpaired_paired_test$patient_ID[unpaired_paired_test$patient_ID == "13093" & unpaired_paired_test$Therapy_2 == "TNF"] <- "2"
unpaired_paired_test$patient_ID[unpaired_paired_test$patient_ID == "9516" & unpaired_paired_test$Therapy_2 == "TNF"] <- "3"

unpaired_paired_test_baseline <- unpaired_paired_test[unpaired_paired_test$Week == "baseline",]
unpaired_paired_test_baseline <- unpaired_paired_test_baseline[order(unpaired_paired_test_baseline$patient_ID) , ]
unpaired_paired_test_baseline$value <- as.numeric(unpaired_paired_test_baseline$value)
unpaired_paired_test_PE <- unpaired_paired_test[!unpaired_paired_test$Week == "baseline",]
unpaired_paired_test_PE <- unpaired_paired_test_PE[order(unpaired_paired_test_PE$patient_ID) , ]
unpaired_paired_test_PE$value <- as.numeric(unpaired_paired_test_PE$value)

#View(unpaired_paired_test)
#View(unpaired_paired_test_baseline)
#View(unpaired_paired_test_PE)

pm.wilcox.test(unpaired_paired_test_PE$value,unpaired_paired_test_baseline$value,alternative = "two.sided", correct = F, method = "MW-MW")

# P-value: 0.5634: baseline, y: PE)
# Xpaired: 18
# Ypaired: 19
# Xextra: 2
# Y-extra: 4

# Caudo (Micro)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- phyloseq_rarefied_phages_ggplot_class_UC[c(phyloseq_rarefied_phages_ggplot_class_UC$`Endoscopic_outcome_combined (NR/R)` == "non-remission" & phyloseq_rarefied_phages_ggplot_class_UC$Diagnosis == "UC" & phyloseq_rarefied_phages_ggplot_class_UC$variable == "Malgrandaviricetes"),]
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired <- total_virome_R_NR_total_phages_w0_UC_caudo_non_crass
#View(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired)

colnames_de <- colnames(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired)
de1 <-data.frame("7164",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("7270",NA, NA,"primary endpoint",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("10449",NA, NA,"primary endpoint",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("12515",NA, NA,"primary endpoint",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("13453",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("9516",NA, NA,"primary endpoint",NA,NA,NA,NA)

#View(de1)
names(de1)<-c(colnames_de)
unpaired_paired_test <- rbind(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired, de1)
unpaired_paired_test$patient_ID[unpaired_paired_test$patient_ID == "12770" & unpaired_paired_test$Therapy_2 == "TNF"] <- "1"
unpaired_paired_test$patient_ID[unpaired_paired_test$patient_ID == "13093" & unpaired_paired_test$Therapy_2 == "TNF"] <- "2"
unpaired_paired_test$patient_ID[unpaired_paired_test$patient_ID == "9516" & unpaired_paired_test$Therapy_2 == "TNF"] <- "3"

unpaired_paired_test_baseline <- unpaired_paired_test[unpaired_paired_test$Week == "baseline",]
unpaired_paired_test_baseline <- unpaired_paired_test_baseline[order(unpaired_paired_test_baseline$patient_ID) , ]
unpaired_paired_test_baseline$value <- as.numeric(unpaired_paired_test_baseline$value)
unpaired_paired_test_PE <- unpaired_paired_test[!unpaired_paired_test$Week == "baseline",]
unpaired_paired_test_PE <- unpaired_paired_test_PE[order(unpaired_paired_test_PE$patient_ID) , ]
unpaired_paired_test_PE$value <- as.numeric(unpaired_paired_test_PE$value)

#View(unpaired_paired_test)
#View(unpaired_paired_test_baseline)
#View(unpaired_paired_test_PE)

pm.wilcox.test(unpaired_paired_test_PE$value,unpaired_paired_test_baseline$value,alternative = "two.sided", correct = F, method = "MW-MW")

# P-value: 0.07408 baseline, y: PE)
# Xpaired: 18
# Ypaired: 19
# Xextra: 2
# Y-extra: 4
####################################
# 6.7.4 CD
####################################
phyloseq_rarefied_phages_ggplot_class_UC <- phyloseq_rarefied_phages_ggplot_class[phyloseq_rarefied_phages_ggplot_class$Diagnosis == "CD",]
phyloseq_rarefied_phages_ggplot_class_UC$Week[phyloseq_rarefied_phages_ggplot_class_UC$Week == "w0"] <- "baseline"
phyloseq_rarefied_phages_ggplot_class_UC <- phyloseq_rarefied_phages_ggplot_class_UC[!phyloseq_rarefied_phages_ggplot_class_UC$Week == "w14",]
phyloseq_rarefied_phages_ggplot_class_UC$Week[phyloseq_rarefied_phages_ggplot_class_UC$Week == "w24"] <- "primary endpoint"
phyloseq_rarefied_phages_ggplot_class_UC <- phyloseq_rarefied_phages_ggplot_class_UC[!phyloseq_rarefied_phages_ggplot_class_UC$`Endoscopic_outcome_combined (NR/R)` == "unknown",]
#View(phyloseq_rarefied_phages_ggplot_class_UC)

phyloseq_rarefied_phages_ggplot_class_UC %>%
  filter(variable %in% c("Caudoviricetes (CrAss)", "Caudoviricetes (non-CrAss)", 
                         "Malgrandaviricetes")) %>%
  ggplot() +
  aes(x = `Endoscopic_outcome_combined (NR/R)`, y = value, fill = variable) +
  geom_boxplot(shape = "circle") +
  theme(axis.text = element_text(colour="black")) +
  facet_wrap(vars(Week)) + 
  xlab("") +
  ylab("Relative Abundance") +
  scale_fill_manual(values=c("#59b031", "#ee766f", "#c49b05"),name="Class") +
  theme_bw() 
# geom_jitter() don't use this because boxplots to close together.

# print long way 
# use different statistic within & between RA.
####################################
# 6.7.5 Statistics within timepoints
####################################
# Statistics within baseline
# Caudo (non-crAss)
table(phyloseq_rarefied_phages_ggplot_class_UC$variable)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- phyloseq_rarefied_phages_ggplot_class_UC[c(phyloseq_rarefied_phages_ggplot_class_UC$Week == "baseline" & phyloseq_rarefied_phages_ggplot_class_UC$Diagnosis == "CD" & phyloseq_rarefied_phages_ggplot_class_UC$variable == "Caudoviricetes (non-CrAss)"),]
table(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$`Endoscopic_outcome_combined (NR/R)`) ## 22 R patients/ 22 NR patients
wilcox.test(value~`Endoscopic_outcome_combined (NR/R)`, data = total_virome_R_NR_total_phages_w0_UC_caudo_non_crass, alternative = "two.sided",  paired = FALSE, exact = F)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$remission <- total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$`Endoscopic_outcome_combined (NR/R)`
effect_size <- wilcox_effsize(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass,value~remission)
(effect_size$effsize)^2
effect_size$n1
effect_size$n2
# P-value = 0.05061
# R-squared: 0.04263385
#  n1: 47
#  n2: 43

# Caudo (crAss)
table(phyloseq_rarefied_phages_ggplot_class_UC$variable)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- phyloseq_rarefied_phages_ggplot_class_UC[c(phyloseq_rarefied_phages_ggplot_class_UC$Week == "baseline" & phyloseq_rarefied_phages_ggplot_class_UC$Diagnosis == "CD" & phyloseq_rarefied_phages_ggplot_class_UC$variable == "Caudoviricetes (CrAss)"),]
table(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$`Endoscopic_outcome_combined (NR/R)`) ## 22 R patients/ 22 NR patients
wilcox.test(value~`Endoscopic_outcome_combined (NR/R)`, data = total_virome_R_NR_total_phages_w0_UC_caudo_non_crass, alternative = "two.sided",  paired = FALSE, exact = F)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$remission <- total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$`Endoscopic_outcome_combined (NR/R)`
effect_size <- wilcox_effsize(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass,value~remission)
(effect_size$effsize)^2
effect_size$n1
effect_size$n2
# P-value = 0.5216
# R-squared: 0.004622027
#  n1: 47
#  n2: 43

# Caudo (Micro)
table(phyloseq_rarefied_phages_ggplot_class_UC$variable)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- phyloseq_rarefied_phages_ggplot_class_UC[c(phyloseq_rarefied_phages_ggplot_class_UC$Week == "baseline" & phyloseq_rarefied_phages_ggplot_class_UC$Diagnosis == "CD" & phyloseq_rarefied_phages_ggplot_class_UC$variable == "Malgrandaviricetes"),]
table(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$`Endoscopic_outcome_combined (NR/R)`) ## 22 R patients/ 22 NR patients
wilcox.test(value~`Endoscopic_outcome_combined (NR/R)`, data = total_virome_R_NR_total_phages_w0_UC_caudo_non_crass, alternative = "two.sided",  paired = FALSE, exact = F)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$remission <- total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$`Endoscopic_outcome_combined (NR/R)`
effect_size <- wilcox_effsize(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass,value~remission)
(effect_size$effsize)^2
effect_size$n1
effect_size$n2
# P-value = 0.2303
# R-squared: 0.0161001
#  n1: 47
#  n2: 43

# Statistics within primary endpoint
# Caudo (non-crAss)
table(phyloseq_rarefied_phages_ggplot_class_UC$variable)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- phyloseq_rarefied_phages_ggplot_class_UC[c(phyloseq_rarefied_phages_ggplot_class_UC$Week == "primary endpoint" & phyloseq_rarefied_phages_ggplot_class_UC$Diagnosis == "CD" & phyloseq_rarefied_phages_ggplot_class_UC$variable == "Caudoviricetes (non-CrAss)"),]
table(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$`Endoscopic_outcome_combined (NR/R)`) ## 22 R patients/ 22 NR patients
wilcox.test(value~`Endoscopic_outcome_combined (NR/R)`, data = total_virome_R_NR_total_phages_w0_UC_caudo_non_crass, alternative = "two.sided",  paired = FALSE, exact = F)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$remission <- total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$`Endoscopic_outcome_combined (NR/R)`
effect_size <- wilcox_effsize(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass,value~remission)
(effect_size$effsize)^2
effect_size$n1
effect_size$n2
# P-value = 0.1794
# R-squared: 0.01928561
#  n1: 48
#  n2: 46

# Caudo (crAss)
table(phyloseq_rarefied_phages_ggplot_class_UC$variable)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- phyloseq_rarefied_phages_ggplot_class_UC[c(phyloseq_rarefied_phages_ggplot_class_UC$Week == "primary endpoint" & phyloseq_rarefied_phages_ggplot_class_UC$Diagnosis == "CD" & phyloseq_rarefied_phages_ggplot_class_UC$variable == "Caudoviricetes (CrAss)"),]
table(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$`Endoscopic_outcome_combined (NR/R)`) ## 22 R patients/ 22 NR patients
wilcox.test(value~`Endoscopic_outcome_combined (NR/R)`, data = total_virome_R_NR_total_phages_w0_UC_caudo_non_crass, alternative = "two.sided",  paired = FALSE, exact = F)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$remission <- total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$`Endoscopic_outcome_combined (NR/R)`
effect_size <- wilcox_effsize(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass,value~remission)
(effect_size$effsize)^2
effect_size$n1
effect_size$n2
# P-value = 0.03898*****
# R-squared: 0.04550569
#  n1: 48
#  n2: 46

# Caudo (Micro)
table(phyloseq_rarefied_phages_ggplot_class_UC$variable)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- phyloseq_rarefied_phages_ggplot_class_UC[c(phyloseq_rarefied_phages_ggplot_class_UC$Week == "primary endpoint" & phyloseq_rarefied_phages_ggplot_class_UC$Diagnosis == "CD" & phyloseq_rarefied_phages_ggplot_class_UC$variable == "Malgrandaviricetes"),]
table(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$`Endoscopic_outcome_combined (NR/R)`) ## 22 R patients/ 22 NR patients
wilcox.test(value~`Endoscopic_outcome_combined (NR/R)`, data = total_virome_R_NR_total_phages_w0_UC_caudo_non_crass, alternative = "two.sided",  paired = FALSE, exact = F)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$remission <- total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$`Endoscopic_outcome_combined (NR/R)`
effect_size <- wilcox_effsize(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass,value~remission)
(effect_size$effsize)^2
effect_size$n1
effect_size$n2
# P-value = .4148
# R-squared: 0.007140787
#  n1: 47
#  n2: 43
####################################
# 6.7.6 Statistics between timepoints
####################################
# Statistics between baseline and primary endpoint
# Statistics within Remission

# caudoviricetes (non-CrAss)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- phyloseq_rarefied_phages_ggplot_class_UC[c(phyloseq_rarefied_phages_ggplot_class_UC$`Endoscopic_outcome_combined (NR/R)` == "remission" & phyloseq_rarefied_phages_ggplot_class_UC$Diagnosis == "CD" & phyloseq_rarefied_phages_ggplot_class_UC$variable == "Caudoviricetes (non-CrAss)"),]
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired <- total_virome_R_NR_total_phages_w0_UC_caudo_non_crass
#View(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass)

colnames_de <- colnames(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired)
de1 <-data.frame("6365",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("12943",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("12868",NA, NA,"baseline",NA,NA,NA,NA)

#View(de1)
names(de1)<-c(colnames_de)
unpaired_paired_test <- rbind(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired, de1)
unpaired_paired_test_baseline <- unpaired_paired_test[unpaired_paired_test$Week == "baseline",]
unpaired_paired_test_baseline <- unpaired_paired_test_baseline[order(unpaired_paired_test_baseline$patient_ID) , ]
unpaired_paired_test_baseline$value <- as.numeric(unpaired_paired_test_baseline$value)
unpaired_paired_test_PE <- unpaired_paired_test[!unpaired_paired_test$Week == "baseline",]
unpaired_paired_test_PE <- unpaired_paired_test_PE[order(unpaired_paired_test_PE$patient_ID) , ]
unpaired_paired_test_PE$value <- as.numeric(unpaired_paired_test_PE$value)

#View(unpaired_paired_test)
#View(unpaired_paired_test_baseline)
#View(unpaired_paired_test_PE)

pm.wilcox.test(unpaired_paired_test_PE$value,unpaired_paired_test_baseline$value,alternative = "two.sided", correct = F, method = "MW-MW")

# P-value = 0.1703 (x: baseline, y: PE)
# Xpaired: 43
# Ypaired: 43
# Xextra: 0
# Y-extra: 3

# Caudo CrAss)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- phyloseq_rarefied_phages_ggplot_class_UC[c(phyloseq_rarefied_phages_ggplot_class_UC$`Endoscopic_outcome_combined (NR/R)` == "remission" & phyloseq_rarefied_phages_ggplot_class_UC$Diagnosis == "CD" & phyloseq_rarefied_phages_ggplot_class_UC$variable == "Caudoviricetes (CrAss)"),]
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired <- total_virome_R_NR_total_phages_w0_UC_caudo_non_crass
#View(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass)

colnames_de <- colnames(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired)
de1 <-data.frame("6365",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("12943",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("12868",NA, NA,"baseline",NA,NA,NA,NA)

#View(de1)
names(de1)<-c(colnames_de)
unpaired_paired_test <- rbind(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired, de1)
unpaired_paired_test_baseline <- unpaired_paired_test[unpaired_paired_test$Week == "baseline",]
unpaired_paired_test_baseline <- unpaired_paired_test_baseline[order(unpaired_paired_test_baseline$patient_ID) , ]
unpaired_paired_test_baseline$value <- as.numeric(unpaired_paired_test_baseline$value)
unpaired_paired_test_PE <- unpaired_paired_test[!unpaired_paired_test$Week == "baseline",]
unpaired_paired_test_PE <- unpaired_paired_test_PE[order(unpaired_paired_test_PE$patient_ID) , ]
unpaired_paired_test_PE$value <- as.numeric(unpaired_paired_test_PE$value)

#View(unpaired_paired_test)
#View(unpaired_paired_test_baseline)
#View(unpaired_paired_test_PE)

pm.wilcox.test(unpaired_paired_test_PE$value,unpaired_paired_test_baseline$value,alternative = "two.sided", correct = F, method = "MW-MW")

# P-value = 0.6691 (x: baseline, y: PE)
# Xpaired: 43
# Ypaired: 43
# Xextra: 0
# Y-extra: 3

# Caudo (Micro)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- phyloseq_rarefied_phages_ggplot_class_UC[c(phyloseq_rarefied_phages_ggplot_class_UC$`Endoscopic_outcome_combined (NR/R)` == "remission" & phyloseq_rarefied_phages_ggplot_class_UC$Diagnosis == "CD" & phyloseq_rarefied_phages_ggplot_class_UC$variable == "Malgrandaviricetes"),]
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired <- total_virome_R_NR_total_phages_w0_UC_caudo_non_crass
#View(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass)

colnames_de <- colnames(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired)
de1 <-data.frame("6365",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("12943",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("12868",NA, NA,"baseline",NA,NA,NA,NA)

#View(de1)
names(de1)<-c(colnames_de)
unpaired_paired_test <- rbind(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired, de1)
unpaired_paired_test_baseline <- unpaired_paired_test[unpaired_paired_test$Week == "baseline",]
unpaired_paired_test_baseline <- unpaired_paired_test_baseline[order(unpaired_paired_test_baseline$patient_ID) , ]
unpaired_paired_test_baseline$value <- as.numeric(unpaired_paired_test_baseline$value)
unpaired_paired_test_PE <- unpaired_paired_test[!unpaired_paired_test$Week == "baseline",]
unpaired_paired_test_PE <- unpaired_paired_test_PE[order(unpaired_paired_test_PE$patient_ID) , ]
unpaired_paired_test_PE$value <- as.numeric(unpaired_paired_test_PE$value)

#View(unpaired_paired_test)
#View(unpaired_paired_test_baseline)
#View(unpaired_paired_test_PE)

pm.wilcox.test(unpaired_paired_test_PE$value,unpaired_paired_test_baseline$value,alternative = "two.sided", correct = F, method = "MW-MW")

# P-value = 0.2508 (x: baseline, y: PE)
# Xpaired: 43
# Ypaired: 43
# Xextra: 0
# Y-extra: 3

# Statistics within non-Remission
# caudoviricetes (non-CrAss)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- phyloseq_rarefied_phages_ggplot_class_UC[c(phyloseq_rarefied_phages_ggplot_class_UC$`Endoscopic_outcome_combined (NR/R)` == "non-remission" & phyloseq_rarefied_phages_ggplot_class_UC$Diagnosis == "CD" & phyloseq_rarefied_phages_ggplot_class_UC$variable == "Caudoviricetes (non-CrAss)"),]
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired <- total_virome_R_NR_total_phages_w0_UC_caudo_non_crass
#View(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired)

colnames_de <- colnames(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired)
de1 <-data.frame("3709",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("5364",NA, NA,"primary endpoint",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("5479",NA, NA,"primary endpoint",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("7265",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("9172",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("10409",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("10777",NA, NA,"primary endpoint",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("11040",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("11169",NA, NA,"primary endpoint",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("11180",NA, NA,"primary endpoint",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("12308",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("12405",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("12946",NA, NA,"primary endpoint",NA,NA,NA,NA)

#View(de1)
names(de1)<-c(colnames_de)
unpaired_paired_test <- rbind(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired, de1)
unpaired_paired_test$patient_ID[unpaired_paired_test$patient_ID == "6716" & unpaired_paired_test$Therapy_2 == "TNF"] <- "1"
unpaired_paired_test_baseline <- unpaired_paired_test[unpaired_paired_test$Week == "baseline",]
unpaired_paired_test_baseline <- unpaired_paired_test_baseline[order(unpaired_paired_test_baseline$patient_ID) , ]
unpaired_paired_test_baseline$value <- as.numeric(unpaired_paired_test_baseline$value)
unpaired_paired_test_PE <- unpaired_paired_test[!unpaired_paired_test$Week == "baseline",]
unpaired_paired_test_PE <- unpaired_paired_test_PE[order(unpaired_paired_test_PE$patient_ID) , ]
unpaired_paired_test_PE$value <- as.numeric(unpaired_paired_test_PE$value)

#View(unpaired_paired_test)
#View(unpaired_paired_test_baseline)
#View(unpaired_paired_test_PE)

pm.wilcox.test(unpaired_paired_test_PE$value,unpaired_paired_test_baseline$value,alternative = "two.sided", correct = F, method = "MW-MW")

# P-value = 0.2222: baseline, y: PE)
# Xpaired: 41
# Ypaired: 41
# Xextra: 7
# Y-extra: 6

# Caudo (CrAss)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- phyloseq_rarefied_phages_ggplot_class_UC[c(phyloseq_rarefied_phages_ggplot_class_UC$`Endoscopic_outcome_combined (NR/R)` == "non-remission" & phyloseq_rarefied_phages_ggplot_class_UC$Diagnosis == "CD" & phyloseq_rarefied_phages_ggplot_class_UC$variable == "Caudoviricetes (CrAss)"),]
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired <- total_virome_R_NR_total_phages_w0_UC_caudo_non_crass
#View(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired)

colnames_de <- colnames(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired)
de1 <-data.frame("3709",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("5364",NA, NA,"primary endpoint",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("5479",NA, NA,"primary endpoint",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("7265",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("9172",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("10409",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("10777",NA, NA,"primary endpoint",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("11040",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("11169",NA, NA,"primary endpoint",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("11180",NA, NA,"primary endpoint",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("12308",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("12405",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("12946",NA, NA,"primary endpoint",NA,NA,NA,NA)

#View(de1)
names(de1)<-c(colnames_de)
unpaired_paired_test <- rbind(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired, de1)
unpaired_paired_test$patient_ID[unpaired_paired_test$patient_ID == "6716" & unpaired_paired_test$Therapy_2 == "TNF"] <- "1"
unpaired_paired_test_baseline <- unpaired_paired_test[unpaired_paired_test$Week == "baseline",]
unpaired_paired_test_baseline <- unpaired_paired_test_baseline[order(unpaired_paired_test_baseline$patient_ID) , ]
unpaired_paired_test_baseline$value <- as.numeric(unpaired_paired_test_baseline$value)
unpaired_paired_test_PE <- unpaired_paired_test[!unpaired_paired_test$Week == "baseline",]
unpaired_paired_test_PE <- unpaired_paired_test_PE[order(unpaired_paired_test_PE$patient_ID) , ]
unpaired_paired_test_PE$value <- as.numeric(unpaired_paired_test_PE$value)

#View(unpaired_paired_test)
#View(unpaired_paired_test_baseline)
#View(unpaired_paired_test_PE)

pm.wilcox.test(unpaired_paired_test_PE$value,unpaired_paired_test_baseline$value,alternative = "two.sided", correct = F, method = "MW-MW")

# P-value = 0.1532
# Xpaired: 41
# Ypaired: 41
# Xextra: 7
# Y-extra: 6

# Caudo (Micro)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- phyloseq_rarefied_phages_ggplot_class_UC[c(phyloseq_rarefied_phages_ggplot_class_UC$`Endoscopic_outcome_combined (NR/R)` == "non-remission" & phyloseq_rarefied_phages_ggplot_class_UC$Diagnosis == "CD" & phyloseq_rarefied_phages_ggplot_class_UC$variable == "Malgrandaviricetes"),]
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired <- total_virome_R_NR_total_phages_w0_UC_caudo_non_crass
#View(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired)

colnames_de <- colnames(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired)
de1 <-data.frame("3709",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("5364",NA, NA,"primary endpoint",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("5479",NA, NA,"primary endpoint",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("7265",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("9172",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("10409",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("10777",NA, NA,"primary endpoint",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("11040",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("11169",NA, NA,"primary endpoint",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("11180",NA, NA,"primary endpoint",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("12308",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("12405",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("12946",NA, NA,"primary endpoint",NA,NA,NA,NA)

#View(de1)
names(de1)<-c(colnames_de)
unpaired_paired_test <- rbind(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired, de1)
unpaired_paired_test$patient_ID[unpaired_paired_test$patient_ID == "6716" & unpaired_paired_test$Therapy_2 == "TNF"] <- "1"
unpaired_paired_test_baseline <- unpaired_paired_test[unpaired_paired_test$Week == "baseline",]
unpaired_paired_test_baseline <- unpaired_paired_test_baseline[order(unpaired_paired_test_baseline$patient_ID) , ]
unpaired_paired_test_baseline$value <- as.numeric(unpaired_paired_test_baseline$value)
unpaired_paired_test_PE <- unpaired_paired_test[!unpaired_paired_test$Week == "baseline",]
unpaired_paired_test_PE <- unpaired_paired_test_PE[order(unpaired_paired_test_PE$patient_ID) , ]
unpaired_paired_test_PE$value <- as.numeric(unpaired_paired_test_PE$value)

#View(unpaired_paired_test)
#View(unpaired_paired_test_baseline)
#View(unpaired_paired_test_PE)

pm.wilcox.test(unpaired_paired_test_PE$value,unpaired_paired_test_baseline$value,alternative = "two.sided", correct = F, method = "MW-MW")

# P-value = 0.3848: baseline, y: PE)
# Xpaired: 41
# Ypaired: 41
# Xextra: 7
# Y-extra: 6
####################################
# 6.7.7 RA boxplots between disease (UC/CD) at primary endpoint to compare inflammed and non-inflammed gut
####################################
phyloseq_rarefied_phages_ggplot_class_UC_CD <- phyloseq_rarefied_phages_ggplot_class
phyloseq_rarefied_phages_ggplot_class_UC_CD$Week[c(phyloseq_rarefied_phages_ggplot_class_UC_CD$Diagnosis == "CD" & phyloseq_rarefied_phages_ggplot_class_UC_CD$Week == "w14")] <- "discard"
phyloseq_rarefied_phages_ggplot_class_UC_CD <- phyloseq_rarefied_phages_ggplot_class_UC_CD[!phyloseq_rarefied_phages_ggplot_class_UC_CD$Week == "discard",]
phyloseq_rarefied_phages_ggplot_class_UC_CD <- phyloseq_rarefied_phages_ggplot_class_UC_CD[!phyloseq_rarefied_phages_ggplot_class_UC_CD$`Endoscopic_outcome_combined (NR/R)` == "unknown",]
phyloseq_rarefied_phages_ggplot_class_UC_CD$Week[phyloseq_rarefied_phages_ggplot_class_UC_CD$Week == "w0"] <- "baseline"
phyloseq_rarefied_phages_ggplot_class_UC_CD$Week[c(phyloseq_rarefied_phages_ggplot_class_UC_CD$Week == "w14" | phyloseq_rarefied_phages_ggplot_class_UC_CD$Week == "w24")] <- "primary endpoint"
phyloseq_rarefied_phages_ggplot_class_UC_CD <- phyloseq_rarefied_phages_ggplot_class_UC_CD[!phyloseq_rarefied_phages_ggplot_class_UC_CD$Week == "baseline",]
phyloseq_rarefied_phages_ggplot_class_UC_CD$Diagnosis <- factor(phyloseq_rarefied_phages_ggplot_class_UC_CD$Diagnosis, levels=c("UC", "CD"))
#View(phyloseq_rarefied_phages_ggplot_class_UC_CD)

# esquisser()
phyloseq_rarefied_phages_ggplot_class_UC_CD %>%
 filter(variable %in% c("Caudoviricetes (CrAss)", "Caudoviricetes (non-CrAss)", 
"Malgrandaviricetes")) %>%
 ggplot() +
 aes(x = `Endoscopic_outcome_combined (NR/R)`, y = value, 
 fill = variable) +
 geom_boxplot(shape = "circle") +
 theme_bw() +
 xlab("") +
 ylab("Relative Abundance") +
 scale_fill_manual(values=c("#59b031", "#ee766f", "#c49b05"),name="Class") +
 facet_wrap(vars(Diagnosis))
####################################
# 7. Distinction Crass versus non-Crass
####################################
# Very often CrAss phages are being related to situations of health, or in our case possibly remission. Evluate if this is correct.
####################################
# 7.1 Merge
####################################
phyloseq_rarefied_phages_ggplot_crass <- phyloseq_rarefied_phages_ggplot
phyloseq_rarefied_phages_ggplot_crass$Class[grepl("crAss",phyloseq_rarefied_phages_ggplot_crass$Species) | grepl("CrAssphage",phyloseq_rarefied_phages_ggplot_crass$Species)] <- "Caudoviricetes (CrAss)"
phyloseq_rarefied_phages_ggplot_crass$Class[phyloseq_rarefied_phages_ggplot_crass$Class == "Caudoviricetes"] <- "Caudoviricetes (non-CrAss)"
table(phyloseq_rarefied_phages_ggplot_crass$Class)
####################################
# 7.2 Aggregate
####################################
vector_1 <- (which(names(phyloseq_rarefied_phages_ggplot_crass)== "Superkingdom")-1) 
ggplot_names <- colnames(phyloseq_rarefied_phages_ggplot_crass[,c(1:vector_1)])

phyloseq_rarefied_phages_ggplot_crass_agg <- aggregate(. ~Class, FUN = sum, data = phyloseq_rarefied_phages_ggplot_crass[,colnames(phyloseq_rarefied_phages_ggplot_crass) %in% ggplot_names | colnames(phyloseq_rarefied_phages_ggplot_crass) == "Class"])
rownames(phyloseq_rarefied_phages_ggplot_crass_agg) <- phyloseq_rarefied_phages_ggplot_crass_agg$Class
phyloseq_rarefied_phages_ggplot_crass_agg$Class <- NULL
####################################
# 7.4 Calculate RA
####################################
phyloseq_rarefied_phages_ggplot_RA <- sweep(phyloseq_rarefied_phages_ggplot_crass_agg, 2, colSums(phyloseq_rarefied_phages_ggplot_crass_agg), '/')
phyloseq_rarefied_phages_ggplot_RA[is.na(phyloseq_rarefied_phages_ggplot_RA)] <- 0
####################################
# 7.5 Combine with metadata
####################################
phyloseq_rarefied_phages_ggplot_RA <- merge(t(phyloseq_rarefied_phages_ggplot_RA), sample_table_rarefied_for_all, by = 0, all = F)
rownames(phyloseq_rarefied_phages_ggplot_RA) <- phyloseq_rarefied_phages_ggplot_RA$Row.names
phyloseq_rarefied_phages_ggplot_RA$Row.names <- NULL
####################################
# 7.6 Melt
####################################
vector_1 <- (which(names(phyloseq_rarefied_phages_ggplot_RA)== "Unannotated")) 
vector_2 <- (which(names(phyloseq_rarefied_phages_ggplot_RA)== "Diagnosis")) 
vector_3 <- (which(names(phyloseq_rarefied_phages_ggplot_RA)== "Gender")) 
vector_4 <- (which(names(phyloseq_rarefied_phages_ggplot_RA)== "Week")) 
vector_5 <- (which(names(phyloseq_rarefied_phages_ggplot_RA)== "Endoscopic_outcome_combined (NR/R)")) 
vector_6 <- (which(names(phyloseq_rarefied_phages_ggplot_RA)== "Therapy_2")) 

phyloseq_rarefied_phages_ggplot_RA <- phyloseq_rarefied_phages_ggplot_RA[,c(1:vector_1,vector_2,vector_3,vector_4,vector_5,vector_6)]
phyloseq_rarefied_phages_ggplot_RA$`Endoscopic_outcome_combined (NR/R)`[phyloseq_rarefied_phages_ggplot_RA$`Endoscopic_outcome_combined (NR/R)` == 1] <- "remission"
phyloseq_rarefied_phages_ggplot_RA$`Endoscopic_outcome_combined (NR/R)`[phyloseq_rarefied_phages_ggplot_RA$`Endoscopic_outcome_combined (NR/R)` == 0] <- "non-remission"
phyloseq_rarefied_phages_ggplot_RA$`Endoscopic_outcome_combined (NR/R)`[phyloseq_rarefied_phages_ggplot_RA$`Endoscopic_outcome_combined (NR/R)` == "NA"] <- "unknown"

phyloseq_rarefied_phages_RA_ggplot <- melt(phyloseq_rarefied_phages_ggplot_RA, id.vars = c("Diagnosis", "Gender", "Week", "Endoscopic_outcome_combined (NR/R)", "Therapy_2"))
####################################
# 7.7 Boxplot 
####################################
# 7.7.1 Create DF
####################################
# No Crass 
phyloseq_rarefied_phages_ggplot_class
phyloseq_rarefied_phages_ggplot_class_UC_no_crass <- phyloseq_rarefied_phages_ggplot_class[phyloseq_rarefied_phages_ggplot_class$Diagnosis == "UC",]
phyloseq_rarefied_phages_ggplot_class_CD_no_crass <- phyloseq_rarefied_phages_ggplot_class[phyloseq_rarefied_phages_ggplot_class$Diagnosis == "CD",]

# Crass
phyloseq_rarefied_phages_RA_ggplot
phyloseq_rarefied_phages_ggplot_class_UC <- phyloseq_rarefied_phages_RA_ggplot[phyloseq_rarefied_phages_RA_ggplot$Diagnosis == "UC",]
phyloseq_rarefied_phages_ggplot_class_CD <- phyloseq_rarefied_phages_RA_ggplot[phyloseq_rarefied_phages_RA_ggplot$Diagnosis == "CD",]

# weeks
phyloseq_rarefied_phages_ggplot_class_UC_CD_w0 <- phyloseq_rarefied_phages_RA_ggplot[phyloseq_rarefied_phages_RA_ggplot$Week == "w0",]
phyloseq_rarefied_phages_ggplot_class_UC_CD_w14 <- phyloseq_rarefied_phages_RA_ggplot[phyloseq_rarefied_phages_RA_ggplot$Week == "w14",]

# weeks and diagnosis
phyloseq_rarefied_phages_ggplot_class_UC_w14 <- phyloseq_rarefied_phages_RA_ggplot[phyloseq_rarefied_phages_RA_ggplot$Diagnosis == "UC" & phyloseq_rarefied_phages_RA_ggplot$Week == "w14",]
phyloseq_rarefied_phages_ggplot_class_CD_w14 <- phyloseq_rarefied_phages_RA_ggplot[phyloseq_rarefied_phages_RA_ggplot$Diagnosis == "CD" & phyloseq_rarefied_phages_RA_ggplot$Week == "w14",]

phyloseq_rarefied_phages_ggplot_class_UC_w0 <- phyloseq_rarefied_phages_RA_ggplot[phyloseq_rarefied_phages_RA_ggplot$Diagnosis == "UC" & phyloseq_rarefied_phages_RA_ggplot$Week == "w0",]
phyloseq_rarefied_phages_ggplot_class_CD_w0 <- phyloseq_rarefied_phages_RA_ggplot[phyloseq_rarefied_phages_RA_ggplot$Diagnosis == "CD" & phyloseq_rarefied_phages_RA_ggplot$Week == "w0",]
####################################
# 7.7.1  IBD
####################################
# 1. IBD patients: 
# 1.1 Non-seperated Classes
# 1.1.1 example 1
ggplot(phyloseq_rarefied_phages_ggplot_class) +
  aes(x = variable, y = value, fill = `Endoscopic_outcome_combined (NR/R)`) +
  geom_boxplot(shape = "circle") +
  scale_fill_hue(direction = 1) +
  theme_bw() +
  facet_wrap(vars(Diagnosis)) +
  theme(axis.text.x = element_text(vjust = 0.5, angle = 70))

# 1.1.2 example 2
phyloseq_rarefied_phages_ggplot_class %>%
  filter(!(`Endoscopic_outcome_combined (NR/R)` %in% "unknown")) %>%
 filter(variable %in% c("Caudoviricetes (CrAss)", "Caudoviricetes (non-CrAss","Malgrandaviricetes"
)) %>%
 ggplot() +
 aes(x = `Endoscopic_outcome_combined (NR/R)`, y = value, 
 fill = `Endoscopic_outcome_combined (NR/R)`) +
 geom_boxplot(shape = "circle") +
 theme_bw() +
  scale_fill_manual(values = list(
    `non-remission` = "#DB4137", remission = "#2CA228")) +
 facet_wrap(vars(variable)) +
 theme(axis.text.x = element_text(vjust = 0.5, angle = 70))

# 1.1.3 example 3
ggplot(phyloseq_rarefied_phages_ggplot_class) +
  aes(x = variable, y = value, fill = variable) +
  geom_boxplot(shape = "circle") +
  scale_fill_hue(direction = 1) +
  theme_minimal() +
  facet_wrap(vars(`Endoscopic_outcome_combined (NR/R)`)) +
  theme_bw() +
  theme(axis.text.x = element_text(vjust = 0.5, angle = 70))

# 1.2 Seperated Classes: CrAss & Non-Crass
# 1.2.1 example 1
ggplot(phyloseq_rarefied_phages_RA_ggplot) +
 aes(x = variable, y = value, fill = `Endoscopic_outcome_combined (NR/R)`) +
 geom_boxplot(shape = "circle") +
 scale_fill_hue(direction = 1) +
 theme_bw() +
 facet_wrap(vars(Diagnosis)) +
theme(axis.text.x = element_text(vjust = 0.5, angle = 70))

# 1.2.2 example 2
phyloseq_rarefied_phages_RA_ggplot %>%
  filter(!(`Endoscopic_outcome_combined (NR/R)` %in% "unknown")) %>%
  filter(variable %in% c("Caudoviricetes (CrAss)", "Caudoviricetes (non-CrAss)", 
                         "Malgrandaviricetes")) %>%
  ggplot() +
  aes(x = `Endoscopic_outcome_combined (NR/R)`, y = value, 
      fill = `Endoscopic_outcome_combined (NR/R)`) +
  geom_boxplot(shape = "circle") +
  scale_fill_manual(values = list(
    `non-remission` = "#DB4137", remission = "#2CA228")) +
  theme_bw() +
  facet_wrap(vars(variable))

# 1.2.3 example 3
ggplot(phyloseq_rarefied_phages_RA_ggplot) +
  aes(x = variable, y = value, fill = variable) +
  geom_boxplot(shape = "circle") +
  scale_fill_hue(direction = 1) +
  theme_minimal() +
  facet_wrap(vars(`Endoscopic_outcome_combined (NR/R)`)) +
  theme_bw() +
  theme(axis.text.x = element_text(vjust = 0.5, angle = 70))

# 1.3 Statistics
# 1.3.1 Caudoviricetes (non-CrAss)
table(total_virome_R_NR_total_phages$`Endoscopic_outcome_combined (NR/R)`)
total_virome_R_NR_total_phages <- phyloseq_rarefied_phages_RA_ggplot[phyloseq_rarefied_phages_RA_ggplot$variable == "Caudoviricetes (non-CrAss)",]
total_virome_R_NR_total_phages <- total_virome_R_NR_total_phages[!total_virome_R_NR_total_phages$`Endoscopic_outcome_combined (NR/R)` == "unknown",]
wilcox.test(value~`Endoscopic_outcome_combined (NR/R)`, data =total_virome_R_NR_total_phages, alternative = "two.sided",  paired = FALSE, conf.int = T, conf.level = 0.95)
wilcox.test(value~`Endoscopic_outcome_combined (NR/R)`, data =total_virome_R_NR_total_phages, alternative = "greater",  paired = FALSE, conf.int = T, conf.level = 0.95)
wilcox.test(value~`Endoscopic_outcome_combined (NR/R)`, data =total_virome_R_NR_total_phages, alternative = "less",  paired = FALSE, conf.int = T, conf.level = 0.95)
# p-value: 0.003576; The RA of caudoviricetes (non-Crass) is significantly different between the remission and non-remission group
# p-value: 0.001788; The RA of caudoviricetes (non-Crass) is significantly greater than in the non-remission group compared to the remission group
# p-value: 0.9982; The RA of caudoviricetes (non-Crass) is NOT significantly less than in the non-remission group compared to the remission group

# 1.3.2 Caudoviricetes (CrAss)
total_virome_R_NR_total_phages <- phyloseq_rarefied_phages_RA_ggplot[phyloseq_rarefied_phages_RA_ggplot$variable == "Caudoviricetes (CrAss)",]
total_virome_R_NR_total_phages <- total_virome_R_NR_total_phages[!total_virome_R_NR_total_phages$`Endoscopic_outcome_combined (NR/R)` == "unknown",]
wilcox.test(value~`Endoscopic_outcome_combined (NR/R)`, data =total_virome_R_NR_total_phages, alternative = "two.sided",  paired = FALSE, conf.int = T, conf.level = 0.95)
wilcox.test(value~`Endoscopic_outcome_combined (NR/R)`, data =total_virome_R_NR_total_phages, alternative = "greater",  paired = FALSE, conf.int = T, conf.level = 0.95)
wilcox.test(value~`Endoscopic_outcome_combined (NR/R)`, data =total_virome_R_NR_total_phages, alternative = "less",  paired = FALSE, conf.int = T, conf.level = 0.95)
# p-value: 0.0002565; The RA of caudoviricetes (Crass) is significantly different between the remission and non-remission group
# p-value: 0.9999; The RA of caudoviricetes (Crass) is NOT significantly greater than in the non-remission group compared to the remission group
# p-value: 0.0001282; The RA of caudoviricetes (Crass) is significantly less than in the non-remission group compared to the remission group

# 1.3.3 Malgrandaviricetes
total_virome_R_NR_total_phages <- phyloseq_rarefied_phages_RA_ggplot[phyloseq_rarefied_phages_RA_ggplot$variable == "Malgrandaviricetes",]
total_virome_R_NR_total_phages <- total_virome_R_NR_total_phages[!total_virome_R_NR_total_phages$`Endoscopic_outcome_combined (NR/R)` == "unknown",]
wilcox.test(value~`Endoscopic_outcome_combined (NR/R)`, data =total_virome_R_NR_total_phages, alternative = "two.sided",  paired = FALSE, conf.int = T, conf.level = 0.95)
wilcox.test(value~`Endoscopic_outcome_combined (NR/R)`, data =total_virome_R_NR_total_phages, alternative = "greater",  paired = FALSE, conf.int = T, conf.level = 0.95)
wilcox.test(value~`Endoscopic_outcome_combined (NR/R)`, data =total_virome_R_NR_total_phages, alternative = "less",  paired = FALSE, conf.int = T, conf.level = 0.95)
# p-value: 0.09188; The RA of Malgrandaviricetes is NOT significantly different between the remission and non-remission group
# p-value: 0.9542; The RA of Malgrandaviricetes is NOT significantly greater than in the non-remission group compared to the remission group
# p-value: 0.04594; The RA of Malgrandaviricetes is significantly less than in the non-remission group compared to the remission group
####################################
# 7.7.2 UC or CD
####################################
# 3. Split by diagnosis
# 3.1 UC
# 3.1.1 No Crass
phyloseq_rarefied_phages_ggplot_class_UC_no_crass %>%
  filter(variable %in% c("Caudoviricetes", "Malgrandaviricetes"
  )) %>%
  ggplot() +
  aes(x = `Endoscopic_outcome_combined (NR/R)`, y = value, 
      fill = `Endoscopic_outcome_combined (NR/R)`) +
  geom_boxplot(shape = "circle") +
  scale_fill_hue(direction = 1) +
  theme_bw() +
  facet_wrap(vars(variable)) +
  theme(axis.text.x = element_text(vjust = 0.5, angle = 70))

# 3.1.2 Crass
phyloseq_rarefied_phages_ggplot_class_UC %>%
  filter(!(`Endoscopic_outcome_combined (NR/R)` %in% "unknown")) %>%
  filter(variable %in% c("Caudoviricetes (non-CrAss)","Caudoviricetes (CrAss)" ,"Malgrandaviricetes"
  )) %>%
  ggplot() +
  aes(x = `Endoscopic_outcome_combined (NR/R)`, y = value, 
      fill = `Endoscopic_outcome_combined (NR/R)`) +
  geom_boxplot(shape = "circle") +
  scale_fill_manual(values = list(
    `non-remission` = "#DB4137", remission = "#2CA228")) +
  theme_bw() +
  facet_wrap(vars(variable)) +
  theme(axis.text.x = element_text(vjust = 0.5, angle = 70))

# 3.1.3 Statistics
# 3.1.3.1 Caudoviricetes (non-CrAss)
table(phyloseq_rarefied_phages_ggplot_class_UC$`Endoscopic_outcome_combined (NR/R)`)
phyloseq_rarefied_phages_ggplot_class_UC_s <- phyloseq_rarefied_phages_ggplot_class_UC[phyloseq_rarefied_phages_ggplot_class_UC$variable == "Caudoviricetes (non-CrAss)",]
phyloseq_rarefied_phages_ggplot_class_UC_s <- phyloseq_rarefied_phages_ggplot_class_UC_s[!phyloseq_rarefied_phages_ggplot_class_UC_s$`Endoscopic_outcome_combined (NR/R)` == "unknown",]
wilcox.test(value~`Endoscopic_outcome_combined (NR/R)`, data = phyloseq_rarefied_phages_ggplot_class_UC_s, alternative = "two.sided",  paired = FALSE, conf.int = T, conf.level = 0.95)
# p-value: 0.1513

# 3.1.3.2 Caudoviricetes (CrAss)
phyloseq_rarefied_phages_ggplot_class_UC_s <- phyloseq_rarefied_phages_ggplot_class_UC[phyloseq_rarefied_phages_ggplot_class_UC$variable == "Caudoviricetes (CrAss)",]
phyloseq_rarefied_phages_ggplot_class_UC_s <- phyloseq_rarefied_phages_ggplot_class_UC_s[!phyloseq_rarefied_phages_ggplot_class_UC_s$`Endoscopic_outcome_combined (NR/R)` == "unknown",]
wilcox.test(value~`Endoscopic_outcome_combined (NR/R)`, data = phyloseq_rarefied_phages_ggplot_class_UC_s, alternative = "two.sided",  paired = FALSE, conf.int = T, conf.level = 0.95)
# p-value: 0.02419

# 3.1.3.3 Malgrandaviricetes
phyloseq_rarefied_phages_ggplot_class_UC_s <- phyloseq_rarefied_phages_ggplot_class_UC[phyloseq_rarefied_phages_ggplot_class_UC$variable == "Malgrandaviricetes",]
phyloseq_rarefied_phages_ggplot_class_UC_s <- phyloseq_rarefied_phages_ggplot_class_UC_s[!phyloseq_rarefied_phages_ggplot_class_UC_s$`Endoscopic_outcome_combined (NR/R)` == "unknown",]
wilcox.test(value~`Endoscopic_outcome_combined (NR/R)`, data = phyloseq_rarefied_phages_ggplot_class_UC_s, alternative = "two.sided",  paired = FALSE, conf.int = T, conf.level = 0.95)
# p-value: 0.6885

# 3.2 CD
# 3.2.1 No Crass
phyloseq_rarefied_phages_ggplot_class_CD_no_crass %>%
  filter(variable %in% c("Caudoviricetes", "Malgrandaviricetes"
  )) %>%
  ggplot() +
  aes(x = `Endoscopic_outcome_combined (NR/R)`, y = value, 
      fill = `Endoscopic_outcome_combined (NR/R)`) +
  geom_boxplot(shape = "circle") +
  scale_fill_hue(direction = 1) +
  theme_bw() +
  facet_wrap(vars(variable)) +
  theme(axis.text.x = element_text(vjust = 0.5, angle = 70))

# 3.2.2 Crass
phyloseq_rarefied_phages_ggplot_class_CD %>%
  filter(!(`Endoscopic_outcome_combined (NR/R)` %in% "unknown")) %>%
  filter(variable %in% c("Caudoviricetes (non-CrAss)","Caudoviricetes (CrAss)" ,"Malgrandaviricetes"
  )) %>%
  ggplot() +
  aes(x = `Endoscopic_outcome_combined (NR/R)`, y = value, 
  fill = `Endoscopic_outcome_combined (NR/R)`) +
  geom_boxplot(shape = "circle") +
  scale_fill_manual(values = list(
    `non-remission` = "#DB4137", remission = "#2CA228")) +
    theme_bw() +
  facet_wrap(vars(variable)) +
  theme(axis.text.x = element_text(vjust = 0.5, angle = 70))

# 3.2.3 Statistics
# 3.2.3.1 Caudoviricetes (non-CrAss)
table(phyloseq_rarefied_phages_ggplot_class_CD$`Endoscopic_outcome_combined (NR/R)`)
phyloseq_rarefied_phages_ggplot_class_CD_s <- phyloseq_rarefied_phages_ggplot_class_CD[phyloseq_rarefied_phages_ggplot_class_CD$variable == "Caudoviricetes (non-CrAss)",]
phyloseq_rarefied_phages_ggplot_class_CD_s <- phyloseq_rarefied_phages_ggplot_class_CD_s[!phyloseq_rarefied_phages_ggplot_class_CD_s$`Endoscopic_outcome_combined (NR/R)` == "unknown",]
wilcox.test(value~`Endoscopic_outcome_combined (NR/R)`, data = phyloseq_rarefied_phages_ggplot_class_CD_s, alternative = "two.sided",  paired = FALSE, conf.int = T, conf.level = 0.95)
# p-value: 0.009414

# 3.2.3.2 Caudoviricetes (CrAss)
phyloseq_rarefied_phages_ggplot_class_CD_s <- phyloseq_rarefied_phages_ggplot_class_CD[phyloseq_rarefied_phages_ggplot_class_UC$variable == "Caudoviricetes (CrAss)",]
phyloseq_rarefied_phages_ggplot_class_CD_s <- phyloseq_rarefied_phages_ggplot_class_CD_s[!phyloseq_rarefied_phages_ggplot_class_CD_s$`Endoscopic_outcome_combined (NR/R)` == "unknown",]
wilcox.test(value~`Endoscopic_outcome_combined (NR/R)`, data = phyloseq_rarefied_phages_ggplot_class_CD_s, alternative = "two.sided",  paired = FALSE, conf.int = T, conf.level = 0.95)
# p-value: 0.0004492

# 3.2.3.3 Malgrandaviricetes
phyloseq_rarefied_phages_ggplot_class_CD_s <- phyloseq_rarefied_phages_ggplot_class_CD[phyloseq_rarefied_phages_ggplot_class_UC$variable == "Malgrandaviricetes",]
phyloseq_rarefied_phages_ggplot_class_CD_s <- phyloseq_rarefied_phages_ggplot_class_CD_s[!phyloseq_rarefied_phages_ggplot_class_CD_s$`Endoscopic_outcome_combined (NR/R)` == "unknown",]
wilcox.test(value~`Endoscopic_outcome_combined (NR/R)`, data = phyloseq_rarefied_phages_ggplot_class_CD_s, alternative = "two.sided",  paired = FALSE, conf.int = T, conf.level = 0.95)
# p-value: 0.8842
####################################
# 7.7.3 Baseline versus timepint
####################################
# 4. Split additionaly in timepoints
# 4.1 UC 
# 4.1.1 w0
phyloseq_rarefied_phages_ggplot_class_UC_w0 %>%
 filter(variable %in% c("Caudoviricetes (CrAss)", "Caudoviricetes (non-CrAss)", 
"Malgrandaviricetes")) %>%
 ggplot() +
 aes(x = `Endoscopic_outcome_combined (NR/R)`, y = value, 
 fill = `Endoscopic_outcome_combined (NR/R)`) +
 geom_boxplot(shape = "circle") +
 scale_fill_manual(values = list(
 `non-remission` = "#DB4137", remission = "#2CA228")) +
 theme_bw() +
 facet_wrap(vars(variable))

# 4.1.2 w14
phyloseq_rarefied_phages_ggplot_class_UC_w14 %>%
  filter(variable %in% c("Caudoviricetes (CrAss)", "Caudoviricetes (non-CrAss)", 
                         "Malgrandaviricetes")) %>%
  ggplot() +
  aes(x = `Endoscopic_outcome_combined (NR/R)`, y = value, 
      fill = `Endoscopic_outcome_combined (NR/R)`) +
  geom_boxplot(shape = "circle") +
  theme_bw() +
  facet_wrap(vars(variable))

# 4.1.3 Statistics
# 4.1.3.1 Caudoviricetes (non-CrAss)
phyloseq_rarefied_phages_ggplot_class_UC_w14_s <- phyloseq_rarefied_phages_ggplot_class_UC_w14[phyloseq_rarefied_phages_ggplot_class_UC_w14$variable == "Caudoviricetes (non-CrAss)",]
phyloseq_rarefied_phages_ggplot_class_UC_w14_s <- phyloseq_rarefied_phages_ggplot_class_UC_w14_s[!phyloseq_rarefied_phages_ggplot_class_UC_w14_s$`Endoscopic_outcome_combined (NR/R)` == "unknown",]
table(phyloseq_rarefied_phages_ggplot_class_UC_w14_s$`Endoscopic_outcome_combined (NR/R)`)
wilcox.test(value~`Endoscopic_outcome_combined (NR/R)`, data = phyloseq_rarefied_phages_ggplot_class_UC_w14_s, alternative = "two.sided",  paired = FALSE, conf.int = T, conf.level = 0.95)
# p-value: 0.0276

# 4.1.3.2 Caudoviricetes (CrAss)
phyloseq_rarefied_phages_ggplot_class_UC_w14_s <- phyloseq_rarefied_phages_ggplot_class_UC_w14[phyloseq_rarefied_phages_ggplot_class_UC_w14$variable == "Caudoviricetes (CrAss)",]
phyloseq_rarefied_phages_ggplot_class_UC_w14_s <- phyloseq_rarefied_phages_ggplot_class_UC_w14_s[!phyloseq_rarefied_phages_ggplot_class_UC_w14_s$`Endoscopic_outcome_combined (NR/R)` == "unknown",]
table(phyloseq_rarefied_phages_ggplot_class_UC_w14_s$`Endoscopic_outcome_combined (NR/R)`)
wilcox.test(value~`Endoscopic_outcome_combined (NR/R)`, data = phyloseq_rarefied_phages_ggplot_class_UC_w14_s, alternative = "two.sided",  paired = FALSE, conf.int = T, conf.level = 0.95)
# p-value: 0.03237

# 4.1.3.3 Malgrandaviricetes
phyloseq_rarefied_phages_ggplot_class_UC_w14_s <- phyloseq_rarefied_phages_ggplot_class_UC_w14[phyloseq_rarefied_phages_ggplot_class_UC_w14$variable == "Malgrandaviricetes",]
phyloseq_rarefied_phages_ggplot_class_UC_w14_s <- phyloseq_rarefied_phages_ggplot_class_UC_w14_s[!phyloseq_rarefied_phages_ggplot_class_UC_w14_s$`Endoscopic_outcome_combined (NR/R)` == "unknown",]
table(phyloseq_rarefied_phages_ggplot_class_UC_w14_s$`Endoscopic_outcome_combined (NR/R)`)
wilcox.test(value~`Endoscopic_outcome_combined (NR/R)`, data = phyloseq_rarefied_phages_ggplot_class_UC_w14_s, alternative = "two.sided",  paired = FALSE, conf.int = T, conf.level = 0.95)
# p-value: 0.2015

# 4.2 CD
# 4.2.1 w0
   phyloseq_rarefied_phages_ggplot_class_CD_w0 %>%
  filter(!(`Endoscopic_outcome_combined (NR/R)` %in% "unknown")) %>%
    filter(variable %in% c("Caudoviricetes (CrAss)", "Caudoviricetes (non-CrAss)", 
                           "Malgrandaviricetes")) %>%
    ggplot() +
    aes(x = `Endoscopic_outcome_combined (NR/R)`, y = value, 
        fill = `Endoscopic_outcome_combined (NR/R)`) +
    geom_boxplot(shape = "circle") +
    scale_fill_manual(values = list(
      `non-remission` = "#DB4137", remission = "#2CA228")) +
    theme_bw() +
    facet_wrap(vars(variable))
   
# Statistics
# Caudoviricetes (non-CrAss)
phyloseq_rarefied_phages_ggplot_class_UC_w14_s <- phyloseq_rarefied_phages_ggplot_class_CD_w0[phyloseq_rarefied_phages_ggplot_class_CD_w0$variable == "Caudoviricetes (non-CrAss)",]
phyloseq_rarefied_phages_ggplot_class_UC_w14_s <- phyloseq_rarefied_phages_ggplot_class_UC_w14_s[!phyloseq_rarefied_phages_ggplot_class_UC_w14_s$`Endoscopic_outcome_combined (NR/R)` == "unknown",]
table(phyloseq_rarefied_phages_ggplot_class_UC_w14_s$`Endoscopic_outcome_combined (NR/R)`)
wilcox.test(value~`Endoscopic_outcome_combined (NR/R)`, data = phyloseq_rarefied_phages_ggplot_class_UC_w14_s, alternative = "two.sided",  paired = FALSE, conf.int = T, conf.level = 0.95)
# p-value: 0.0325
   
# Caudoviricetes (CrAss)
phyloseq_rarefied_phages_ggplot_class_UC_w14_s <- phyloseq_rarefied_phages_ggplot_class_CD_w0[phyloseq_rarefied_phages_ggplot_class_CD_w0$variable == "Caudoviricetes (CrAss)",]
phyloseq_rarefied_phages_ggplot_class_UC_w14_s <- phyloseq_rarefied_phages_ggplot_class_UC_w14_s[!phyloseq_rarefied_phages_ggplot_class_UC_w14_s$`Endoscopic_outcome_combined (NR/R)` == "unknown",]
table(phyloseq_rarefied_phages_ggplot_class_UC_w14_s$`Endoscopic_outcome_combined (NR/R)`)
wilcox.test(value~`Endoscopic_outcome_combined (NR/R)`, data = phyloseq_rarefied_phages_ggplot_class_UC_w14_s, alternative = "two.sided",  paired = FALSE, conf.int = T, conf.level = 0.95)
# p-value: 0.1823
   
# Malgrandaviricetes
phyloseq_rarefied_phages_ggplot_class_UC_w14_s <- phyloseq_rarefied_phages_ggplot_class_CD_w0[phyloseq_rarefied_phages_ggplot_class_CD_w0$variable == "Malgrandaviricetes",]
phyloseq_rarefied_phages_ggplot_class_UC_w14_s <- phyloseq_rarefied_phages_ggplot_class_UC_w14_s[!phyloseq_rarefied_phages_ggplot_class_UC_w14_s$`Endoscopic_outcome_combined (NR/R)` == "unknown",]
table(phyloseq_rarefied_phages_ggplot_class_UC_w14_s$`Endoscopic_outcome_combined (NR/R)`)
wilcox.test(value~`Endoscopic_outcome_combined (NR/R)`, data = phyloseq_rarefied_phages_ggplot_class_UC_w14_s, alternative = "two.sided",  paired = FALSE, conf.int = T, conf.level = 0.95)
# p-value: 0.2164
   
# 4.2.2 w14
phyloseq_rarefied_phages_ggplot_class_CD_w14 %>%
     filter(!(`Endoscopic_outcome_combined (NR/R)` %in% "unknown")) %>%
  filter(!(`Therapy_2` %in% "UST")) %>%
  filter(variable %in% c("Caudoviricetes (CrAss)", "Caudoviricetes (non-CrAss)", 
                            "Malgrandaviricetes")) %>%
     ggplot() +
     aes(x = `Endoscopic_outcome_combined (NR/R)`, y = value, 
         fill = `Endoscopic_outcome_combined (NR/R)`) +
     geom_boxplot(shape = "circle") +
     scale_fill_manual(values = list(
       `non-remission` = "#DB4137", remission = "#2CA228")) +
     theme_bw() +
     facet_wrap(vars(variable))

# Statistics: w24
# Caudoviricetes (non-CrAss)
phyloseq_rarefied_phages_ggplot_class_UC_w14_s <- phyloseq_rarefied_phages_ggplot_class_CD_w24[phyloseq_rarefied_phages_ggplot_class_CD_w24$variable == "Caudoviricetes (non-CrAss)",]
phyloseq_rarefied_phages_ggplot_class_UC_w14_s <- phyloseq_rarefied_phages_ggplot_class_UC_w14_s[!phyloseq_rarefied_phages_ggplot_class_UC_w14_s$`Endoscopic_outcome_combined (NR/R)` == "unknown",]
table(phyloseq_rarefied_phages_ggplot_class_UC_w14_s$`Endoscopic_outcome_combined (NR/R)`)
wilcox.test(value~`Endoscopic_outcome_combined (NR/R)`, data = phyloseq_rarefied_phages_ggplot_class_UC_w14_s, alternative = "two.sided",  paired = FALSE, conf.int = T, conf.level = 0.95)
# p-value: 0.2487

# Caudoviricetes (CrAss)
phyloseq_rarefied_phages_ggplot_class_UC_w14_s <- phyloseq_rarefied_phages_ggplot_class_CD_w24[phyloseq_rarefied_phages_ggplot_class_CD_w24$variable == "Caudoviricetes (CrAss)",]
phyloseq_rarefied_phages_ggplot_class_UC_w14_s <- phyloseq_rarefied_phages_ggplot_class_UC_w14_s[!phyloseq_rarefied_phages_ggplot_class_UC_w14_s$`Endoscopic_outcome_combined (NR/R)` == "unknown",]
table(phyloseq_rarefied_phages_ggplot_class_UC_w14_s$`Endoscopic_outcome_combined (NR/R)`)
wilcox.test(value~`Endoscopic_outcome_combined (NR/R)`, data = phyloseq_rarefied_phages_ggplot_class_UC_w14_s, alternative = "two.sided",  paired = FALSE, conf.int = T, conf.level = 0.95)
# p-value: 0.02448

# Malgrandaviricetes
phyloseq_rarefied_phages_ggplot_class_UC_w14_s <- phyloseq_rarefied_phages_ggplot_class_CD_w24[phyloseq_rarefied_phages_ggplot_class_CD_w24$variable == "Malgrandaviricetes",]
phyloseq_rarefied_phages_ggplot_class_UC_w14_s <- phyloseq_rarefied_phages_ggplot_class_UC_w14_s[!phyloseq_rarefied_phages_ggplot_class_UC_w14_s$`Endoscopic_outcome_combined (NR/R)` == "unknown",]
table(phyloseq_rarefied_phages_ggplot_class_UC_w14_s$`Endoscopic_outcome_combined (NR/R)`)
wilcox.test(value~`Endoscopic_outcome_combined (NR/R)`, data = phyloseq_rarefied_phages_ggplot_class_UC_w14_s, alternative = "two.sided",  paired = FALSE, conf.int = T, conf.level = 0.95)
# p-value: 0.4734
####################################
# 7.7.4 Split by treatment
####################################
# 5. Split additionaly in timepoints
# 5.1 IBD
phyloseq_rarefied_phages_RA_ggplot %>%
  filter(!(`Endoscopic_outcome_combined (NR/R)` %in% "unknown")) %>%
  filter(variable %in% c("Caudoviricetes (CrAss)", "Caudoviricetes (non-CrAss)", 
                         "Malgrandaviricetes")) %>%
  ggplot() +
  aes(x = `Endoscopic_outcome_combined (NR/R)`, y = value, 
      fill = Therapy_2) +
  geom_boxplot(shape = "circle") +
  scale_fill_hue(direction = 1) +
  theme_bw() +
  facet_wrap(vars(variable))

# 5.2 UC 
phyloseq_rarefied_phages_ggplot_class_UC %>%
  filter(!(`Endoscopic_outcome_combined (NR/R)` %in% "unknown")) %>%
 filter(variable %in% c("Caudoviricetes (CrAss)", "Caudoviricetes (non-CrAss)", 
"Malgrandaviricetes")) %>%
 ggplot() +
 aes(x = `Endoscopic_outcome_combined (NR/R)`, y = value, 
 fill = Therapy_2) +
 geom_boxplot(shape = "circle") +
 scale_fill_hue(direction = 1) +
 theme_bw() +
 facet_wrap(vars(variable))

# 5.3 CD
phyloseq_rarefied_phages_ggplot_class_CD %>%
  filter(!(`Endoscopic_outcome_combined (NR/R)` %in% "unknown")) %>%
  filter(variable %in% c("Caudoviricetes (CrAss)", "Caudoviricetes (non-CrAss)", 
                         "Malgrandaviricetes")) %>%
  ggplot() +
  aes(x = `Endoscopic_outcome_combined (NR/R)`, y = value, 
      fill = Therapy_2) +
  geom_boxplot(shape = "circle") +
  scale_fill_hue(direction = 1) +
  theme_bw() +
  facet_wrap(vars(variable))
####################################

######### !!!!!!!!!!!!!!! #########
# SAVE Global Environment #
######### !!!!!!!!!!!!!!! #########