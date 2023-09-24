####################################
# SCRIPT 4: PROKARYOTIC VIRUSES
####################################
# 0. Packages: Install packages and load them wherever needed
####################################
#install.packages("FSA")
library(lme4)
library (microbiome)
library(FSA)
library(metafor)
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
#install.packages("ape")
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
# 1. Datasets
####################################
Mastertable # Reads & contigs figure
Mastertable_genus1 # Percentage classified on class, order, family level
Mastertable_eukaryotic # Eukaryotic viral figure
Mastertable_phage  # Phage analysis
####################################
# 2. Upload Metadata
####################################
# 2.1  Upload Metadata per sample
#################################### 
setwd("/Users/daan/Desktop/Bioinformatics/pre-processing/FMT_study/Metadata")
dir()

metadata <- read_excel("samples_FMT_project_Final_DJ_R.xlsx")
metadata <- as.data.frame(metadata)
metadata$Bacterial_cellcount <- as.numeric(metadata$Bacterial_cellcount)
metadata$CRP <- as.numeric(metadata$CRP)
metadata$Hemoglobin <- as.numeric(metadata$Hemoglobin)
metadata$Age <- as.numeric(metadata$Age)
metadata$BMI <- as.numeric(metadata$BMI)
metadata$Disease_duration <- as.numeric(metadata$Disease_duration)
metadata$Smoking <- as.character(metadata$Smoking)
metadata$Calprotectin <- as.numeric(metadata$Calprotectin)
metadata$Moisture <- as.numeric(metadata$Moisture)
rownames(metadata) <- metadata$Sample_ID
colnames(metadata)
metadata[metadata=="NA"] <- NA
str(metadata)

samples <- metadata
#################################### 
# 2.2  Upload Patient metadata
#################################### 
setwd("/Users/daan/Desktop/Bioinformatics/pre-processing/FMT_study/Metadata")
dir()

metadata <- read_excel("samples_FMT_project_Final_DJ_Patient_R.xlsx")
metadata <- as.data.frame(metadata)
metadata[metadata=="NA"] <- NA
rownames(metadata) <- metadata$Patient_identifier
metadata$Age <- as.numeric(metadata$Age)
metadata$BMI <- as.numeric(metadata$BMI)
metadata$Disease_duration <- as.numeric(metadata$Disease_duration)
metadata$Biologicals <- as.character(metadata$Biologicals)
metadata$Steroids <- as.character(metadata$Steroids)
metadata$`5_ASA` <- as.character(metadata$`5_ASA`)
metadata$Smoking <- as.character(metadata$Smoking)
metadata$IBD_treatment_number <- as.numeric(metadata$IBD_treatment_number)
metadata$Treatment_number <- as.numeric(metadata$Treatment_number)
metadata_P <- metadata
str(metadata_P)
View(metadata_P)

## All data
metadata_P
Taxonomy_table
#################################### 
# 3. Subset prokaryotic viruses
####################################
# 3.1 Samples
####################################
Mastertable_viral_rarefied_2 <- Mastertable_phage

# Select phages
table(Mastertable_viral_rarefied_2$Final_viral)
sort(colSums(Mastertable_viral_rarefied_2[1:304]))
names <- colnames(Mastertable_viral_rarefied_2[1:304][,colSums(Mastertable_viral_rarefied_2[1:304])!=0])
names # 302 samples with phages lefr

# Select samples with phages present and create abundance table
abundance_table <- Mastertable_viral_rarefied_2
abundance_table$names <- rownames(abundance_table)
abundance_tabe1 <- aggregate(. ~names, FUN = sum, data = abundance_table[,colnames(abundance_table) %in% names | colnames(abundance_table) == 'names'])
rownames(abundance_tabe1) <- abundance_tabe1$names
abundance_tabe1$names <- NULL
colnames(abundance_tabe1)

# Create taxonomy table
colnames(Mastertable_viral_rarefied_2)
Mastertable_viral_rarefied_2$Final_ANI <- NULL
Mastertable_viral_rarefied_2$Final_coverage <- NULL
Mastertable_viral_rarefied_2$Totalnumberofreads <- NULL
Mastertable_viral_rarefied_2$Blastn_AS <- NULL
Mastertable_viral_rarefied_2$Diamond_AS <- NULL
Mastertable_viral_rarefied_2$Best_AS <- NULL
Mastertable_viral_rarefied_2$Final_genus <- Mastertable_viral_rarefied_2$Genus
vector_1 <- (which(names(Mastertable_viral_rarefied_2)== "Final_superkingdom"))
vector_2 <- which(names(Mastertable_viral_rarefied_2)== "Final_species")
Taxonomy_table <- Mastertable_viral_rarefied_2[,c(vector_1:vector_2)]

# Rename some Taxonomy elements
colnames(Taxonomy_table)
names(Taxonomy_table)[names(Taxonomy_table) == "Final_superkingdom"] <- "Kingdom"
names(Taxonomy_table)[names(Taxonomy_table) == "Final_phylum"] <- "Phylum"
names(Taxonomy_table)[names(Taxonomy_table) == "Final_class"] <- "Class"
names(Taxonomy_table)[names(Taxonomy_table) == "Final_order"] <- "Order"
names(Taxonomy_table)[names(Taxonomy_table) == "Final_family"] <- "Family"
names(Taxonomy_table)[names(Taxonomy_table) == "Final_subfamily"] <- "Subfamily"
names(Taxonomy_table)[names(Taxonomy_table) == "Final_genus"] <- "Genus"
names(Taxonomy_table)[names(Taxonomy_table) == "Final_species"] <- "Species"
Taxonomy_table$Subfamily <- NULL

Taxonomy_table <- Taxonomy_table[,c(7,1, 2, 3,4,5,6)]
#View(Taxonomy_table)
####################################
# 3.2 Patient
####################################
Mastertable_phage1 <- Mastertable_phage[1:304]

# Select phages
sort(colSums(Mastertable_phage1[1:304]))
names <- colnames(Mastertable_phage1[1:304][,colSums(Mastertable_phage1[1:304])!=0])
names 

# Select samples with phages present and create abundance table
abundance_table_P <- Mastertable_phage1
abundance_table_P$names <- rownames(abundance_table_P)
abundance_table_Pa <- aggregate(. ~names, FUN = sum, data = abundance_table_P[,colnames(abundance_table_P) %in% names | colnames(abundance_table_P) == 'names'])
rownames(abundance_table_Pa) <- abundance_table_Pa$names
abundance_table_Pa$names <- NULL
colnames(abundance_table_Pa)

# # Aggregated by patient number
Mastertable_phage1_t <- t(abundance_table_Pa)
Mastertable_phage1_t_met <- merge(Mastertable_phage1_t, samples, by=0, all =F)
rownames(Mastertable_phage1_t_met) <- Mastertable_phage1_t_met$Row.names
Mastertable_phage1_t_met$Row.names <- NULL

vectorA <- which(colnames(Mastertable_phage1_t_met) == "Sample_ID")-1
vector_patient <- colnames(Mastertable_phage1_t_met[, c(1:vectorA)])

# Aggregated by patient number
Mastertable_phage1_t_met$Patient_ID_2
Mastertable_patients_1 <- aggregate(. ~Patient_ID_2, FUN = sum, data = Mastertable_phage1_t_met[,colnames(Mastertable_phage1_t_met) %in% vector_patient | colnames(Mastertable_phage1_t_met) == 'Patient_ID_2'])
rownames(Mastertable_patients_1) <- Mastertable_patients_1$Patient_ID_2
Mastertable_patients_1$Patient_ID_2 <- NULL

# Remove column with sum zero 
Mastertable_patients_DF <-  as.data.frame(t(Mastertable_patients_1))
abundance_table1 <- Mastertable_patients_DF[, which(colSums(Mastertable_patients_DF) != 0)]
abundance_table_P <- abundance_table1
#View(abundance_table_P)
####################################
# 4. Create phyloSeq objects
####################################
# 4.1 Samples
####################################
abundance_table_rarefied_m_samples <- as.matrix(abundance_tabe1)
taxonomy_table_rarefied_m <- as.matrix(Taxonomy_table)
colnames(abundance_table_rarefied_m_samples)
colnames(taxonomy_table_rarefied_m)
####################################
# 4.2 Patient
####################################
abundance_table_rarefied_m_patient <- as.matrix(abundance_table_P)
taxonomy_table_rarefied_m_P <- as.matrix(Taxonomy_table)
####################################
# 5.Phyloseq object
####################################
# 5.1 Samples: Transform to phylseq objects
####################################
ABUNDANCE_rarefied_family <- otu_table(abundance_table_rarefied_m_samples, taxa_are_rows = T)
TAX_rarefied_family <- tax_table(taxonomy_table_rarefied_m)
samples <- sample_data(metadata)
phyloseq_rarefied_samples <- phyloseq(ABUNDANCE_rarefied_family,TAX_rarefied_family, samples)
phyloseq_rarefied_samples
####################################
# 5.1 Patient: Transform to phylseq objects
####################################
ABUNDANCE_rarefied_family <- otu_table(abundance_table_rarefied_m_patient, taxa_are_rows = T)
TAX_rarefied_family <- tax_table(taxonomy_table_rarefied_m_P)
metadata_P <- sample_data(metadata_P)
phyloseq_rarefied_patient <- phyloseq(ABUNDANCE_rarefied_family,TAX_rarefied_family, metadata_P)
rownames(sample_data(phyloseq_rarefied_patient))
####################################
# 5.2 Visualize data & subset only the phage fraction
###################################
sample_names(phyloseq_rarefied_samples) # All sample names
rank_names(phyloseq_rarefied_samples) # All taxonomies
sample_variables(phyloseq_rarefied_samples) ## All metadata
phyloseq_rarefied_phages <- phyloseq_rarefied_samples
phyloseq_rarefied_phages_P <-phyloseq_rarefied_patient
####################################
# 5.3 Transform data with Hellinger
####################################
# 5.3.1 Samples
####################################
phyloseq_phages <- phyloseq_rarefied_phages
colnames(sample_data(phyloseq_phages))

#phyloseq_phages_P = subset_samples(phyloseq_phages, Disease_status == "patient")
#phyloseq_phages = subset_samples(phyloseq_phages_P, Disease_status == "patient")

phyloseq_rarefied_phages_transformed <- phyloseq_phages %>% 
  microbiome::transform(transform = "hellinger") %>%
  aggregate_rare(level = "Genus", detection = 0.0005/100, prevalence = 0.0001/100)

colnames(sample_data(phyloseq_rarefied_phages_transformed))

# Add clusters
# alpha_diversity_cluster <- sample_data(alpha_diversity_cluster)
alpha_diversity_cluster <- sample_data(cluster)
phyloseq_rarefied_phages_transformed <- merge_phyloseq(phyloseq_rarefied_phages_transformed, alpha_diversity_cluster)

# Distances
phyloseq_rarefied_phages_transformed_bray <- phyloseq::distance(phyloseq_rarefied_phages_transformed, method="bray")
phyloseq_rarefied_phages_PcoA <- ordinate(phyloseq_rarefied_phages_transformed, method="PCoA", distance=phyloseq_rarefied_phages_transformed_bray)

# order metada variables
metadata_Pcoa <- as.data.frame(sample_data(phyloseq_rarefied_phages_transformed))
metadata_Pcoa$Disease_status <- factor(metadata_Pcoa$Disease_status, levels = c("patient", "donor"))
phyloseq_rarefied_phages_transformed <- merge_phyloseq(metadata_Pcoa, phyloseq_rarefied_phages_transformed)
colnames(sample_data(phyloseq_rarefied_phages_transformed))
####################################
# 5.3.2 Patient
####################################
phyloseq_phages_P <- phyloseq_rarefied_patient
phyloseq_phages_P

phyloseq_rarefied_phages_transformed <- phyloseq_phages_P %>% 
  microbiome::transform(transform = "hellinger") %>%
  aggregate_rare(level = "Genus", detection = 0.000001/100, prevalence = 0.00001/100)

phyloseq_rarefied_phages_transformed_bray <- phyloseq::distance(phyloseq_rarefied_phages_transformed, method="bray")
phyloseq_rarefied_phages_PcoA <- ordinate(phyloseq_rarefied_phages_transformed, method="PCoA", distance=phyloseq_rarefied_phages_transformed_bray)

# order metada variables
metadata_Pcoa <- as.data.frame(sample_data(phyloseq_rarefied_phages_transformed))
metadata_Pcoa$Disease_status <- factor(metadata_Pcoa$Disease_status, levels = c("patient", "donor"))
phyloseq_rarefied_phages_transformed <- merge_phyloseq(metadata_Pcoa, phyloseq_rarefied_phages_transformed)
colnames(metadata_Pcoa)
####################################
# 6. PcoA Visualization
####################################
# 6.1 Exploratory PcoA: samples
####################################
# Visualization PcoA 1: Disease status
sample_data(phyloseq_rarefied_phages_transformed)

plot_ordination(phyloseq_rarefied_phages_transformed, phyloseq_rarefied_phages_PcoA, color="cluster") +
  geom_point(alpha=1, size=1) +
  theme_bw() +
  scale_color_manual(values = c("1" = "#bd7f1c", "2" = "#2a8a92")) +
  geom_hline(yintercept = 0, linetype = 3, colour = "black") + 
  geom_vline(xintercept = 0, linetype = 3, colour = "black") +
  xlab("PC2 (7.2)") +
  ylab("PC1 (8.6)")

#  scale_color_manual(values = c("Patient (-autologous)" = "#AD6518", "Donor (-autologous)" = "#cf9079", "Patient (-untreated)" = "#7b8fd1", "Donor (-superdonor)" = "#a1cf79","Patient (-superdonor)" = "#648f47")) +
#  geom_hline(yintercept = 0, linetype = 3, colour = "black") + 
#  geom_vline(xintercept = 0, linetype = 3, colour = "black") +
  
# PcoA 2: Disease_extent 
  plot_ordination(phyloseq_rarefied_phages_transformed, phyloseq_rarefied_phages_PcoA, color="Disease_extent") +
    geom_point(size = 1, alpha=1) +
    theme_bw() 

# PcoA 3: FMT_Treatment_type
plot_ordination(phyloseq_rarefied_phages_transformed, phyloseq_rarefied_phages_PcoA, color="FMT_treatment_type") +
    geom_point(size = 1, alpha=1) +
    theme_bw() 

# PcoA 3: FMT_Treatment_type2
plot_ordination(phyloseq_rarefied_phages_transformed, phyloseq_rarefied_phages_PcoA, color="FMT_treatment_type2") +
  geom_point(size = 1, alpha=1) +
  scale_color_manual(values = c("Patient (-autologous)" = "#1DA0AA", "Donor (-autologous)" = "#1ECBE1", "Patient (-untreated)" = "#683fc0", "Donor (-superdonor)" = "#81C837","Patient (-superdonor)" = "#6EA52F")) +
  geom_hline(yintercept = 0, linetype = 3, colour = "black") + 
  geom_vline(xintercept = 0, linetype = 3, colour = "black") +
  xlab("PC2 (9.1)") +
  ylab("PC1 (7.3)") +
  theme_bw() 

# PcoA 3: FMT_Treatment_type2
plot_ordination(phyloseq_rarefied_phages_transformed, phyloseq_rarefied_phages_PcoA, color="FMT_treatment_type2") +
  geom_point(size = 1, alpha=1) +
  #scale_color_manual(values = c("Patient (-autologous)" = "#D78371", "Donor (-autologous)" = "#993C29", "Patient (-untreated)" = "#993C29", "Donor (-superdonor)" = "#1D8F6C","Patient (-superdonor)" = "#95DAD5")) +
 # scale_color_manual(values = c("Patient (-autologous)" = "#D1A1A2", "Donor (-autologous)" = "#AA726A", "Patient (-untreated)" = "#EBB5B6", "Donor (-superdonor)" = "#1A2D38","Patient (-superdonor)" = "#A1BFD1")) +
  scale_color_manual(values = c("Patient (-autologous)" = "#AD6518", "Donor (-autologous)" = "#cf9079", "Patient (-untreated)" = "#7b8fd1", "Donor (-superdonor)" = "#a1cf79","Patient (-superdonor)" = "#648f47")) +
  geom_hline(yintercept = 0, linetype = 3, colour = "black") + 
  geom_vline(xintercept = 0, linetype = 3, colour = "black") +
  xlab("PC2 (9.1)") +
  ylab("PC1 (7.3)") +
  theme_bw()

# PcoA 4: Mesalamine   
plot_ordination(phyloseq_rarefied_phages_transformed, phyloseq_rarefied_phages_PcoA, color="Concomitant_treatment") +
  geom_point(size = 1, alpha=1) +
  theme_bw() 

# PcoA 5: Steroids   
plot_ordination(phyloseq_rarefied_phages_transformed, phyloseq_rarefied_phages_PcoA, color="Steroid_concomitant") +
  geom_point(size = 1, alpha=1) +
  theme_bw() 

# PcoA 6: Biologicals   
plot_ordination(phyloseq_rarefied_phages_transformed, phyloseq_rarefied_phages_PcoA, color="Biologicals_concomitant") +
  geom_point(size = 1, alpha=1) +
  theme_bw() 
####################################
# 6.2 Exploratory PcoA: patient
####################################
# Visualization PcoA 1: Disease status
plot_ordination(phyloseq_rarefied_phages_transformed, phyloseq_rarefied_phages_PcoA, color="Disease_extent") +
  geom_point(size = 1, alpha=1) 
  theme_bw() 

# PcoA : Type of donors 
plot_ordination(phyloseq_rarefied_phages_transformed, phyloseq_rarefied_phages_PcoA, color="Disease_status") +
  geom_point(size = 1, alpha=1) +
  theme_bw() +
  scale_color_manual(values = c("patient" = "#bd7f1c", "donor" = "#2a8a92")) +
  geom_hline(yintercept = 0, linetype = 3, colour = "black") + 
  geom_vline(xintercept = 0, linetype = 3, colour = "black") +
  xlab("PC2 (9.1)") 
  ylab("PC1 (7.3)")

# PcoA 3: patient SE
plot_ordination(phyloseq_rarefied_phages_transformed, phyloseq_rarefied_phages_PcoA, color="FMT_Treatment_type") +
  geom_point(size = 1, alpha=1) +
  theme_bw() 

# PcoA 4: Disease_extent 
plot_ordination(phyloseq_rarefied_phages_transformed, phyloseq_rarefied_phages_PcoA, color="Disease_extent") +
  geom_point(size = 1, alpha=1) +
  theme_bw() 

# PcoA 5: X5_ASA 
plot_ordination(phyloseq_rarefied_phages_transformed, phyloseq_rarefied_phages_PcoA, color="X5_ASA") +
  geom_point(size = 1, alpha=1) +
  theme_bw() 

# PcoA 6: X5_ASA 
plot_ordination(phyloseq_rarefied_phages_transformed, phyloseq_rarefied_phages_PcoA, color="Biologicals") +
  geom_point(size = 1, alpha=1) +
  theme_bw() 
####################################
# 5.5.3 Envfit
####################################
library(vegan)

# Select metadata
metadata <- data.frame(sample_data(phyloseq_rarefied_phages_transformed))
colnames(metadata)
metadata_patient_1 <- metadata[7]

# Create DF 
scrs <- as.data.frame(phyloseq_rarefied_phages_PcoA$vectors[,1:2])

# Combine with metadata
scrs <- merge(scrs,metadata_patient_1,by=0,all=F)
rownames(scrs) <- scrs$Row.names
scrs$Row.names <- NULL

# Envfit 
set.seed(123)
phyloseq_rarefied_phages_transformed_1 <- data.frame(sample_data(phyloseq_rarefied_phages_transformed))
phyloseq_rarefied_phages_transformed_2 <- phyloseq_rarefied_phages_transformed_1[,7:8] # At least two variables
vf <- envfit(phyloseq_rarefied_phages_PcoA$vectors, phyloseq_rarefied_phages_transformed_2, perm = 999)
vf 

# Select factors you want to plot
spp.scrs <- as.data.frame(scores(vf, display = "factors"))
spp.scrs$Group <- rownames(spp.scrs)

spp.scrs$Group
spp.scrs <- spp.scrs[spp.scrs$Group=="FMT_treatment_type2Donor (-autologous)" |
                 spp.scrs$Group=="FMT_treatment_type2Patient (-autologous)" |
                 spp.scrs$Group=="FMT_treatment_type2Patient (-untreated)" |
                 spp.scrs$Group=="FMT_treatment_type2Donor (-superdonor)" | 
                 spp.scrs$Group=="FMT_treatment_type2Patient (-superdonor)"
                ,]

# Figure with arrows
ggplot(scrs) +
  geom_point(mapping = aes(x = Axis.1, y = Axis.2, colour = FMT_treatment_type2), alpha=1, size=1.1) +
  coord_fixed() +
  theme_bw() +
  scale_color_manual(values = c("Patient (-autologous)" = "#AD6518", "Donor (-autologous)" = "#cf9079", "Patient (-untreated)" = "#7b8fd1", "Donor (-superdonor)" = "#a1cf79","Patient (-superdonor)" = "#648f47")) +
  geom_hline(yintercept = 0, linetype = 3, colour = "black") + 
  geom_vline(xintercept = 0, linetype = 3, colour = "black") +
    xlab("PC2 (9.1)") +
    ylab("PC1 (7.3)") +
  geom_segment(data = spp.scrs,
              aes(x = 0, xend = Axis.1, y = 0, yend = Axis.2),
               arrow = arrow(length = unit(0.25, "cm")), colour = "black") + 
    geom_text(data = spp.scrs, aes(x = Axis.1, y = Axis.2, label = Group, colour="grey"),
              size = 3) 

# Boxplot adding to PcoA plot
ggplot(scrs) +
 aes(x = FMT_treatment_type2, y = Axis.1, color = FMT_treatment_type2) +
 geom_boxplot(lwd=0.8) +
  scale_color_manual(values = c("Patient (-autologous)" = "#AD6518", "Donor (-autologous)" = "#cf9079", "Patient (-untreated)" = "#7b8fd1", "Donor (-superdonor)" = "#a1cf79","Patient (-superdonor)" = "#648f47")) +
  theme_minimal()

ggplot(scrs) +
  aes(x = FMT_treatment_type2, y = Axis.2, color = FMT_treatment_type2) +
  geom_boxplot(lwd=0.8) +
  scale_color_manual(values = c("Patient (-autologous)" = "#AD6518", "Donor (-autologous)" = "#cf9079", "Patient (-untreated)" = "#7b8fd1", "Donor (-superdonor)" = "#a1cf79","Patient (-superdonor)" = "#648f47")) +
  theme_minimal()


# Save PcoA plot colors arrows accordingly..
# Save barplot as insert and put underneath, legend on left.

# PS: you can adapt size en shape of arrows best in illustrator.
####################################
# 6. Linear regression analysis
####################################
# 6.1 Regression analysis
####################################
# 6.1.1 Code variables
####################################
distance = "bray"
fdr = 0.05
####################################
# 6.1.2 Univariate anlysis: Baseline patients
####################################
# A) Input
####################################
phyloseq_rarefied_phages_transformed1 <- subset_samples(phyloseq_rarefied_phages_transformed, Timepoints=="Baseline" & grepl("^patient", sample_data(phyloseq_rarefied_phages_transformed)$Patient_ID_2))
phyloseq_rarefied_phages_transformed1 <- subset_samples(phyloseq_rarefied_phages_transformed, !Timepoints=="Baseline" & !FMT_treatment_type2=="Donor (-superdonor)" & grepl("^patient", sample_data(phyloseq_rarefied_phages_transformed)$Patient_ID_2))
phyloseq_rarefied_phages_transformed1 <- subset_samples(phyloseq_rarefied_phages_transformed)
GENUS <-as.data.frame(t(otu_table(phyloseq_rarefied_phages_transformed1))) #genus table
m <- data.frame(sample_data(phyloseq_rarefied_phages_transformed1))

# All
m$Sample_ID <- NULL
#m$Patient_ID <- NULL
m$Patient_ID_2 <- NULL
m$Patient_identifier2 <- NULL
m$Patient_identifier <- NULL
m$FMT_treatment_sucess <- NULL
m$FMT_treatment_donor <- NULL
m$FMT_treatment_batch <- NULL
m$FMT_treatment_type <- NULL
m$Enterotype <- NULL
m[m == "NA"] <- NA
m$FMT_treatment_type3 <- NULL
m$Mesalamine_concomitant <- NULL
m$Thiopurine_concomitant <- NULL
m$Biologicals_concomitant <- NULL
m$Steroid_concomitant <- NULL
View(m)

# Patients baseline
m$Sample_ID <- NULL
#m$Patient_ID <- NULL
m$Patient_ID_2 <- NULL
m$Patient_identifier2 <- NULL
m$Patient_identifier <- NULL
m$FMT_treatment_sucess <- NULL
m$FMT_treatment_donor <- NULL
m$FMT_treatment_batch <- NULL
m$FMT_treatment_type <- NULL
m$Enterotype <- NULL
m$FMT_treatment_type3 <- NULL
m$Disease_status <- NULL
#m$FMT_treatment_type2 <- NULL
#m$Timepoints <- NULL
m$Mesalamine_concomitant <- NULL
m$Thiopurine_concomitant <- NULL
m$Biologicals_concomitant <- NULL
m$Steroid_concomitant <- NULL
colnames(m)
View(m)
####################################
# B) Rows should be ordered the same
####################################
meta

data_rows <- rownames(m)
genus_rows <- rownames(GENUS)

# Check if the row names are ordered the same
if (!identical(metadata_rows, genus_rows)) {
  # Reorder the metadata data frame based on the row order of the genus data frame
  m <- m[genus_rows, ]
  
  # Update the row names of the metadata data frame
  rownames(m) <- genus_rows
  
  print("The rows of metadata have been reordered to match the genus table.")
} else {
  print("The rows of metadata are already ordered the same as the genus table.")
}
####################################
# C) Univariate regression analyses
####################################
all <- c()
for (ii in 1:ncol(m)) { 
  capsc <- capscale(GENUS ~ m[,ii], distance = distance, na.action=na.exclude,strata=m$Patient_ID) 
  an <- anova.cca(capsc) 
  pval <- an["Pr(>F)"][[1]][[1]] 
  Fa <- an["F"][[1]][[1]] 
  r2 <- RsquareAdj(capsc)[[1]] 
  r2adj <- RsquareAdj(capsc)[[2]] 
  all <- rbind(all,cbind(Fa,r2,r2adj,pval))
}

FDR = p.adjust(all[,"pval"],method="BH")
all = cbind(all,FDR)
colnames(all) <- c("F","r2","r2adj","p-value","FDR") #generate table
row.names(all) <- colnames(m)
all <- as.data.frame(all)
all_DF <- all[rev(order(all$r2adj)),]
all_DF
# 23 of the 24 variables are significant
####################################
# D) Multivariate analysis
####################################
# All
####################################
# Univariate and multivariate linear mixed-effects model
# These models still assume normality of underlying distribution which is often not the case for microbiome data. Because of this reason we perform first Hellinger transformation to give more weight to the more rare taxa and create a more normal distribution.
sig.vars = row.names(all[all[,"FDR"] < fdr & all[,"p-value"] < 0.05,])
m0 <- get_variable(phyloseq_rarefied_phages_transformed1,sig.vars) 
sig.vars_without_NA <- c("Patient_ID","Disease_extent","FMT_treatment_type2", "Disease_status","Smoking","Concomitant_treatment")
m0 <- get_variable(phyloseq_rarefied_phages_transformed1,sig.vars_without_NA) 
m <- m0

GENUS0 <- t(otu_table(phyloseq_rarefied_phages_transformed1))
GENUS=GENUS0[row.names(GENUS0) %in% row.names(m),] 

# 1) All samples: Linear mixed effect model constrained on strate patient_id to account for multiple samples
attach(m)
mod0=capscale(GENUS ~ 1, distance=distance,strata = m$Patient_ID,  na.action=na.exclude) # Model with intercept only
mod1=capscale(GENUS ~ . , data=m, distance=distance,strata = m$Patient_ID, na.action=na.exclude) # Model with all expl. variables
set.seed(1)
step.res<-ordiR2step(mod0, scope=formula(mod1), data=m, direction="forward", Pin = 0.05, R2scope = FALSE, pstep = 1000, trace = F)
step.res$anova
anova(step.res) # Summary table of stepwise ordination

ordiR2step.tab = data.matrix(step.res$anova[,1:5])
row.names(ordiR2step.tab) = gsub("\\+ ","",row.names(ordiR2step.tab))
ordiR2step.tab
all_DF
detach(m)
####################################
# Patients (-before treatment) 
####################################
sig.vars = row.names(all[all[,"FDR"] < fdr & all[,"p-value"] < 0.05 ,])
m0 <- get_variable(phyloseq_rarefied_phages_transformed1,sig.vars) 
m <- na.exclude(m0) 
GENUS0 <- t(otu_table(phyloseq_rarefied_phages_transformed1))
GENUS=GENUS0[row.names(GENUS0) %in% row.names(m),] 
View(m)

# 2) Patient after treatment: only calprotectin
# --> linear mixed effect model accoutn for dependent samples
attach(m)
distance="bray"
mod0=capscale(GENUS ~ 1, distance=distance,strata = m$Patient_ID) # Model with intercept only
mod1=capscale(GENUS ~ . - Patient_ID, data=m, distance=distance,strata = m$Patient_ID) # Model with all expl. variables
set.seed(1)
step.res<-ordiR2step(mod0, scope=formula(mod1), data=m, direction="forward", Pin = 0.05, R2scope = FALSE, pstep = 1000, trace = T)
step.res$anova
anova(step.res) # Summary table of stepwise ordination

ordiR2step.tab = data.matrix(step.res$anova[,1:5])
row.names(ordiR2step.tab) = gsub("\\+ ","",row.names(ordiR2step.tab))
ordiR2step.tab
all_DF
detach(m)
####################################
# Patients (-after treatment)
####################################
sig.vars = row.names(all[all[,"FDR"] < fdr & all[,"p-value"] < 0.05,])
m0 <- get_variable(phyloseq_rarefied_phages_transformed1,sig.vars) 
sig.vars_without_NA <- c("Patient_ID","Disease_extent","FMT_treatment_type2", "Age","Smoking","Concomitant_treatment","Endoscopic_outcome")
m0 <- get_variable(phyloseq_rarefied_phages_transformed1,sig.vars_without_NA) 
m <- m0

GENUS0 <- t(otu_table(phyloseq_rarefied_phages_transformed1))
GENUS=GENUS0[row.names(GENUS0) %in% row.names(m),] 

# 2) Patient after treatment: 
# --> linear mixed effect model accoutn for dependent samples
attach(m)
distance="bray"
mod0=capscale(GENUS ~ 1, distance=distance,strata = m$Patient_ID, na.action=na.exclude) # Model with intercept only
mod1=capscale(GENUS ~ . -Patient_ID, data=m, distance=distance,strata = m$Patient_ID, na.action=na.exclude) # Model with all expl. variables
set.seed(1)
step.res<-ordiR2step(mod0, scope=formula(mod1), data=m, direction="forward", Pin = 0.05, R2scope = FALSE, pstep = 1000, trace = T)
step.res$anova
anova(step.res) # Summary table of stepwise ordination

ordiR2step.tab = data.matrix(step.res$anova[,1:5])
row.names(ordiR2step.tab) = gsub("\\+ ","",row.names(ordiR2step.tab))
ordiR2step.tab
all_DF
detach(m)
####################################
# 6.1.3 Effect size figures
####################################
# All individuals
ordiR2step.tab
all_DF

DF_variables <- data.frame(
  variables = c("Patient ID", "Treatment type", "Treatment (-concomitant)" ,"Patient ID", "Treatment type", "Treatment (-concomitant)"),
  variable_type = c("multivariate", "multivariate","multivariate", "univariate", "univariate","univariate"),
  values = c(54.9, 55.6, 55.7,54.9, 3.13,1.30))

ord <- c("Patient ID", "Treatment type","Treatment (-concomitant)")
DF_variables$variables <- factor(DF_variables$variables, levels = rev(ord))

ggplot(DF_variables, aes(x = variables, fill = variable_type, weight = values)) +
  geom_bar(position = "dodge") +
  scale_fill_manual(values = c(multivariate = "#8B8A8B", univariate = "#000000")) +
  coord_flip() +
  theme_bw() +
  labs(x = "covariates of virome composition", y = "effect size (%)") +
  theme(legend.title = element_blank()) +
  geom_hline(yintercept = 55.6, linetype = 5, colour = "#49006a")


# Baseline (-treated patients)
ordiR2step.tab
all_DF

DF_variables <- data.frame(
  variables = c("Disease extent","Treatment type","Endoscopic outcome","Age","Treatment (-concomitant)", "Smoking","Disease extent","Treatment type","Endoscopic outcome","Age","Treatment (-concomitant)", "Smoking"),
  variable_type = c("multivariate","multivariate","multivariate","multivariate","multivariate","multivariate","univariate","univariate", "univariate","univariate","univariate","univariate"),
  values = c(2.1, 3.4, 5.0, 5.9,6.7,7.6,
             2.1,1.4,1.5,0.9,0.9,0.8))

ord <- c("Disease extent","Treatment type","Endoscopic outcome","Age","Treatment (-concomitant)", "Smoking")
DF_variables$variables <- factor(DF_variables$variables, levels = rev(ord))

ggplot(DF_variables, aes(x = variables, fill = variable_type, weight = values)) +
  geom_bar(position = "dodge") +
  scale_fill_manual(values = c(multivariate = "#8B8A8B", univariate = "#000000")) +
  coord_flip() +
  theme_bw() +
  labs(x = "covariates of virome composition", y = "effect size (%)") +
  theme(legend.title = element_blank()) +
  geom_hline(yintercept = 7.6, linetype = 5, colour = "#49006a")
####################################
# 6.1.4 PcoA coloured by significant variables
####################################
ordiR2step.tab

# Obtain Baseline patient samples
phyloseq_rarefied_phages_transformed1 <- subset_samples(phyloseq_rarefied_phages_transformed)
phyloseq_rarefied_phages_transformed1_bray <- phyloseq::distance(phyloseq_rarefied_phages_transformed1, method="bray")
phyloseq_rarefied_phages_PcoA_B <- ordinate(phyloseq_rarefied_phages_transformed1, method="PCoA", distance=phyloseq_rarefied_phages_transformed1_bray)

# PcoA 1: Disease.status 
plot_ordination(phyloseq_rarefied_phages_transformed1, phyloseq_rarefied_phages_PcoA_B, color="Calprotectin") +
  geom_point(size = 1, alpha=1) +
  theme_bw() +
  geom_hline(yintercept = 0, linetype = 3, colour = "black") + 
  geom_vline(xintercept = 0, linetype = 3, colour = "black") 

# PcoA 2: Endoscopic.outcome 
plot_ordination(phyloseq_rarefied_phages_transformed1, phyloseq_rarefied_phages_PcoA_B, color="FMT_treatment_type2", shape="Concomitant_treatment") +
  geom_point(size = 1, alpha=1) +
  theme_bw() +
  scale_color_manual(values = c("Patient (-autologous)" = "#AD6518", "Donor (-autologous)" = "#cf9079", "Patient (-untreated)" = "#7b8fd1", "Donor (-superdonor)" = "#a1cf79","Patient (-superdonor)" = "#648f47")) +
  geom_hline(yintercept = 0, linetype = 3, colour = "black") + 
  geom_vline(xintercept = 0, linetype = 3, colour = "black") 
  xlab("PC2 (8.6)") +
  ylab("PC1 (7.2)") +
  theme_bw()

# PcoA 3: Endoscopic.outcome 
plot_ordination(phyloseq_rarefied_phages_transformed1, phyloseq_rarefied_phages_PcoA_B, color="Disease_extent") +
  geom_point(size = 1, alpha=1) +
  theme_bw() +
  geom_hline(yintercept = 0, linetype = 3, colour = "black") + 
  geom_vline(xintercept = 0, linetype = 3, colour = "black") 

# PcoA 4: Moisture 
plot_ordination(phyloseq_rarefied_phages_transformed1, phyloseq_rarefied_phages_PcoA_B, color="Moisture") +
  geom_point(size = 1, alpha=1) +
  theme_bw() +
  geom_hline(yintercept = 0, linetype = 3, colour = "black") + 
  geom_vline(xintercept = 0, linetype = 3, colour = "black") 

# PcoA 5: Steroids 
plot_ordination(phyloseq_rarefied_phages_transformed1, phyloseq_rarefied_phages_PcoA_B, color="Steroids") +
  geom_point(size = 1, alpha=1) +
  theme_bw() +
  geom_hline(yintercept = 0, linetype = 3, colour = "black") + 
  geom_vline(xintercept = 0, linetype = 3, colour = "black") 

# PcoA 6: steroids 
plot_ordination(phyloseq_rarefied_phages_transformed1, phyloseq_rarefied_phages_PcoA_B, color="X5_ASA") +
  geom_point(size = 1, alpha=1) +
  theme_bw() +
  geom_hline(yintercept = 0, linetype = 3, colour = "black") + 
  geom_vline(xintercept = 0, linetype = 3, colour = "black") 
####################################
# 6.1.5 Associations between important variables 
####################################
phyloseq_rarefied_phages_transformed1 <- subset_samples(phyloseq_rarefied_phages_transformed, Timepoints=="Baseline") 
phyloseq_rarefied_phages_transformed1 <- subset_samples(phyloseq_rarefied_phages_transformed)
GENUS <-as.data.frame(t(otu_table(phyloseq_rarefied_phages_transformed1))) #genus table
m <- data.frame(sample_data(phyloseq_rarefied_phages_transformed1))
#View(sample_data(phyloseq_rarefied_phages_transformed1))

chisq.test(m$Calprotectin, m$Moisture)
m[is.na(m)] <- "NA"
m1 <- m[!m$Moisture=="NA" & !m$Hemoglobin=="NA",]
m1$Calprotectin <- as.numeric(m1$Calprotectin)
m1$Moisture <- as.numeric(m1$Moisture)
m1$Hemoglobin <- as.numeric(m1$Hemoglobin)
cor.test(m1$Moisture,m1$Calprotectin,method = "spearman")
cor.test(m1$Hemoglobin,m1$Calprotectin,method = "spearman")
cor.test(m1$Hemoglobin,m1$Moisture,method = "spearman")
wilcox.test(Calprotectin~Endoscopic_outcome,data=m1)
wilcox.test(Moisture~Endoscopic_outcome,data=m1)
chisq.test(m1$Disease_extent,m1$Endoscopic_outcome)
####################################
# 7. Barplot of strong associations between covariates (> 1%) and FMT treatment 
####################################
# Results of generalised linear mixed models are communicated in a similar way to results for linear models. In your results section you should mention that you are using mixed models with R package lme4, and list your random and fixed effects
# To assess if the effect of either of these variables is associated to FMT_treatment while accounting for sample dependency we used generalized linear mixed effect models
# "Patient_ID" is included as a random effect in the model. It represents the individual patients in the dataset, and the random intercept specified by "(1 | Patient_ID)" allows for patient-specific variation in the outcome. The random effect accounts for the potential correlation or clustering of observations within the same patient, capturing any unobserved patient-level factors that may influence the outcome.
# "FMT_treatment" and "Concomitant_treatment" are considered fixed effects, while "Patient_ID" is included as a random effect.
# You have univariate and multivariate generalized mixed effect models: we use univariate!!
# The different between linear mixed effect models and generalized linear mixed effect models is that the latter does not assume normality and can be used for non-parametric data
# Show numbers and percentages in Excel
####################################
# 7.1 FMT_ treatment ~ Disease extent (generalized linear mixed-effect model)
####################################
# Subset the relevant samples from the phyloseq object
phyloseq_subset <- subset_samples(phyloseq_rarefied_phages_transformed, !Timepoints=="Baseline" & !FMT_treatment_type2=="Donor (-superdonor)" & grepl("^patient", sample_data(phyloseq_rarefied_phages_transformed)$Patient_ID_2))

# Extract the relevant variables from the sample data
metadata <- sample_data(phyloseq_subset)
FMT_treatment <- metadata$FMT_treatment_type2
Disease_extent <- metadata$Disease_extent
Patient_ID <- metadata$Patient_ID
combined_table <- data.frame(FMT_treatment,Disease_extent,Patient_ID)
nrow(combined_table)

# Create a contingency table
contingency_table <- table(combined_table$FMT_treatment, combined_table$Disease_extent)
chisq.test(contingency_table) # This is not the test we will use, but gives me a quick indication
contingency_table

# Count N per disease extent
# E1:N=16
# E2:N=87
# E3: N=36
# N=16+87+36=139

# Calculate proportions from the contingency table
proportions <- prop.table(contingency_table, margin = 1)
percentages <- round(proportions * 100, 1)
percentages

# Convert the percentages to a data frame
percentages_df <- as.data.frame(percentages)
percentages_df

# Rename the column names for better visualization
colnames(percentages_df) <- c("FMT Treatment", "Disease extent", "Percentage")

# Create the bar plot using ggplot2
plot <- ggplot(percentages_df) +
  aes(x = `Disease extent`, fill = `FMT Treatment`, weight = Percentage) +
  geom_bar(position = "dodge") +
  scale_fill_manual(values = c(`Patient (-autologous)` = "#AE6618", `Patient (-superdonor)` = "#658F46")) +
  labs(y = "Percentage (%)") +
  geom_text(aes(label = paste0(Percentage, "%"), y = Percentage), position = position_dodge(width = 1), vjust = -0.5) +
  theme_bw()

# Calculate statistics: generalized linear mixed-effect model to account for sample dependencies (multiple timepoints a cross same patients)
# Convert the 'FMT_treatment' variable to numeric (1 or 0)
combined_table$FMT_treatment <- as.numeric(combined_table$FMT_treatment == "Patient (-superdonor)")
str(combined_table)

# Fit the generalized linear mixed-effects model (GLMM)
# If problems increase nAGQ
fit <- glmer(FMT_treatment ~ Disease_extent + (1 | Patient_ID), data = combined_table, family = binomial,nAGQ = 20)

# random = ~intercept + fixed effect | random effect
summary(fit)
Anova(fit)

# Extract p-value from model summary
anova_result <- Anova(fit)
p_value <- round(anova_result$"Pr(>Chisq)", 3)
p_value

# Add p-value as a text annotation to the plot
plot + annotate("text", x = 2.5, y = 50, label = paste("p =", p_value), size = 5)
####################################
# 7.2 FMT_ treatment ~ Endoscopic outcome
####################################
# Subset the relevant samples from the phyloseq object
phyloseq_subset <- subset_samples(phyloseq_rarefied_phages_transformed, !Timepoints=="Baseline" & !FMT_treatment_type2=="Donor (-superdonor)" & grepl("^patient", sample_data(phyloseq_rarefied_phages_transformed)$Patient_ID_2))

# Extract the relevant variables from the sample data
metadata <- sample_data(phyloseq_subset)
metadata$Endoscopic_outcome[metadata$Endoscopic_outcome=="1"] <- "Endoscopic remission"
metadata$Endoscopic_outcome[metadata$Endoscopic_outcome=="0"] <- "Endoscopic non-remission"
FMT_treatment <- metadata$FMT_treatment_type2
Endoscopic_outcome <- metadata$Endoscopic_outcome
Patient_ID <- metadata$Patient_ID
combined_table <- data.frame(FMT_treatment,Endoscopic_outcome,Patient_ID)
nrow(combined_table)

# Create a contingency table
contingency_table <- table(combined_table$FMT_treatment, combined_table$Endoscopic_outcome)
chisq.test(contingency_table) # This is not the test we will use, but gives me a quick indication
contingency_table

# Count N per endoscopic outcome
# Endoscopic remission:N=93
# Endoscopic non-remission:N=46
# N=93+46=139

# Calculate proportions from the contingency table
proportions <- prop.table(contingency_table, margin = 1)
percentages <- round(proportions * 100, 1)
percentages

# Convert the percentages to a data frame
percentages_df <- as.data.frame(percentages)

# Rename the column names for better visualization
colnames(percentages_df) <- c("FMT Treatment", "Endoscopic outcome", "Percentage")
percentages_df

# Create the bar plot using ggplot2
plot <- ggplot(percentages_df) +
  aes(x = `Endoscopic outcome`, fill = `FMT Treatment`, weight = Percentage) +
  geom_bar(position = "dodge") +
  scale_fill_manual(values = c(`Patient (-autologous)` = "#AE6618", `Patient (-superdonor)` = "#658F46")) +
  labs(y = "Percentage (%)") +
  geom_text(aes(label = paste0(Percentage, "%"), y = Percentage), position = position_dodge(width = 1), vjust = -0.5) +
  theme_bw()

# Calculate statistics: generalized linear mixed-effect model to account for sample dependencies (multiple timepoints a cross same patients)
# Convert the 'FMT_treatment' variable to numeric (1 or 0)
combined_table$FMT_treatment <- as.numeric(combined_table$FMT_treatment == "Patient (-superdonor)")
str(combined_table)

# Fit the generalized linear mixed-effects model (GLMM)
# If problems increase nAGQ
fit <- glmer(FMT_treatment ~ Endoscopic_outcome + (1 | Patient_ID), data = combined_table, family = binomial,nAGQ = 20)

# random = ~intercept + fixed effect | random effect
summary(fit)
Anova(fit)

# Extract p-value from model summary
anova_result <- Anova(fit)
p_value <- round(anova_result$"Pr(>Chisq)", 3)
p_value

# Add p-value as a text annotation to the plot
plot + annotate("text", x = 2.5, y = 50, label = paste("p =", p_value), size = 5)
####################################
# 7.3 FMT_treatment ~ Age
####################################
# Subset the relevant samples from the phyloseq object
phyloseq_subset <- subset_samples(phyloseq_rarefied_phages_transformed, !Timepoints=="Baseline" & !FMT_treatment_type2=="Donor (-superdonor)" & grepl("^patient", sample_data(phyloseq_rarefied_phages_transformed)$Patient_ID_2))

# Extract the relevant variables from the sample data
metadata <- sample_data(phyloseq_subset)
FMT_treatment <- metadata$FMT_treatment_type2
Age <- metadata$Age
Patient_ID <- metadata$Patient_ID
combined_table <- data.frame(FMT_treatment,Age,Patient_ID)
nrow(combined_table)

# Convert the 'FMT_treatment' variable to numeric (1 or 0)
combined_table$FMT_treatment <- as.numeric(combined_table$FMT_treatment == "Patient (-superdonor)")

# Fit the generalized linear mixed-effects model (GLMM)
# If problems increase nAGQ
fit <- glmer(FMT_treatment ~ Age + (1 | Patient_ID), data = combined_table, family = binomial,nAGQ = 20)

# random = ~intercept + fixed effect | random effect
summary(fit)
Anova(fit)

# Extract p-value from model summary
anova_result <- Anova(fit)
p_value <- round(anova_result$"Pr(>Chisq)", 3)
p_value
####################################
# 7.4  FMT_treatment ~ Treatment (-concomitant)
####################################
# Subset the relevant samples from the phyloseq object
phyloseq_subset <- subset_samples(phyloseq_rarefied_phages_transformed, !Timepoints=="Baseline" & !FMT_treatment_type2=="Donor (-superdonor)" & grepl("^patient", sample_data(phyloseq_rarefied_phages_transformed)$Patient_ID_2))

# Extract the relevant variables from the sample data
metadata <- sample_data(phyloseq_subset)
FMT_treatment <- metadata$FMT_treatment_type2
Concomitant_treatment <- metadata$Concomitant_treatment
Patient_ID <- metadata$Patient_ID
combined_table <- data.frame(FMT_treatment,Concomitant_treatment,Patient_ID)
nrow(combined_table)

# Create a contingency table
contingency_table <- table(combined_table$FMT_treatment, combined_table$Concomitant_treatment)
chisq.test(contingency_table) # This is not the test we will use, but gives me a quick indication
contingency_table

# Count N per Concomitant_treatment
# Concomitant_treatment N=64
# No Concomitant_treatmentN=75
# N=64+75=139

# Calculate proportions from the contingency table
proportions <- prop.table(contingency_table, margin = 1)
percentages <- round(proportions * 100, 1)
percentages

# Convert the percentages to a data frame
percentages_df <- as.data.frame(percentages)

# Rename the column names for better visualization
colnames(percentages_df) <- c("FMT Treatment", "Concomitant treatment", "Percentage")
percentages_df

# Create the bar plot using ggplot2
plot <- ggplot(percentages_df) +
  aes(x = `Concomitant treatment`, fill = `FMT Treatment`, weight = Percentage) +
  geom_bar(position = "dodge") +
  scale_fill_manual(values = c(`Patient (-autologous)` = "#AE6618", `Patient (-superdonor)` = "#658F46")) +
  labs(y = "Percentage (%)") +
  geom_text(aes(label = paste0(Percentage, "%"), y = Percentage), position = position_dodge(width = 1), vjust = -0.5) +
  theme_bw()

# Calculate statistics: generalized linear mixed-effect model to account for sample dependencies (multiple timepoints a cross same patients)
# Convert the 'FMT_treatment' variable to numeric (1 or 0)
combined_table$FMT_treatment <- as.numeric(combined_table$FMT_treatment == "Patient (-superdonor)")
str(combined_table)

# Fit the generalized linear mixed-effects model (GLMM)
# If problems increase nAGQ
fit <- glmer(FMT_treatment ~ Concomitant_treatment + (1 | Patient_ID), data = combined_table, family = binomial,nAGQ = 20)

# random = ~intercept + fixed effect | random effect
summary(fit)
Anova(fit)

# Extract p-value from model summary
anova_result <- Anova(fit)
p_value <- round(anova_result$"Pr(>Chisq)", 3)
p_value

# Add p-value as a text annotation to the plot
plot + annotate("text", x = 2.5, y = 50, label = paste("p =", p_value), size = 5)
####################################
# 7.5  FMT_treatment ~ Smoking
####################################
# Subset the relevant samples from the phyloseq object
phyloseq_subset <- subset_samples(phyloseq_rarefied_phages_transformed, !Timepoints=="Baseline" & !FMT_treatment_type2=="Donor (-superdonor)" & grepl("^patient", sample_data(phyloseq_rarefied_phages_transformed)$Patient_ID_2))

# Extract the relevant variables from the sample data
metadata <- sample_data(phyloseq_subset)
metadata$Smoking[metadata$Smoking=="1"] <- "Smoker"
metadata$Smoking[metadata$Smoking=="0"] <- "Non-smoker"
FMT_treatment <- metadata$FMT_treatment_type2
Smoking <- metadata$Smoking
Patient_ID <- metadata$Patient_ID
combined_table <- data.frame(FMT_treatment,Smoking,Patient_ID)
nrow(combined_table)

# Create a contingency table
contingency_table <- table(combined_table$FMT_treatment, combined_table$Smoking)
chisq.test(contingency_table) # This is not the test we will use, but gives me a quick indication
contingency_table

# Count N per Concomitant_treatment
# Smoking N=10
# No Smoking N=129
# N=10+129=139

# Calculate proportions from the contingency table
proportions <- prop.table(contingency_table, margin = 1)
percentages <- round(proportions * 100, 1)
percentages

# Convert the percentages to a data frame
percentages_df <- as.data.frame(percentages)

# Rename the column names for better visualization
colnames(percentages_df) <- c("FMT Treatment", "Smoking", "Percentage")
percentages_df

# Create the bar plot using ggplot2
plot <- ggplot(percentages_df) +
  aes(x = `Smoking`, fill = `FMT Treatment`, weight = Percentage) +
  geom_bar(position = "dodge") +
  scale_fill_manual(values = c(`Patient (-autologous)` = "#AE6618", `Patient (-superdonor)` = "#658F46")) +
  labs(y = "Percentage (%)") +
  geom_text(aes(label = paste0(Percentage, "%"), y = Percentage), position = position_dodge(width = 1), vjust = -0.5) +
  theme_bw()

# Calculate statistics: generalized linear mixed-effect model to account for sample dependencies (multiple timepoints a cross same patients)
# Convert the 'FMT_treatment' variable to numeric (1 or 0)
combined_table$FMT_treatment <- as.numeric(combined_table$FMT_treatment == "Patient (-superdonor)")
str(combined_table)

# Fit the generalized linear mixed-effects model (GLMM)
# If problems increase nAGQ
fit <- glmer(FMT_treatment ~ Smoking + (1 | Patient_ID), data = combined_table, family = binomial,nAGQ = 20)

# random = ~intercept + fixed effect | random effect
summary(fit)
Anova(fit)

# Extract p-value from model summary
anova_result <- Anova(fit)
p_value <- round(anova_result$"Pr(>Chisq)", 3)
p_value

# Add p-value as a text annotation to the plot
plot + annotate("text", x = 2.5, y = 50, label = paste("p =", p_value), size = 5)
####################################

###################################