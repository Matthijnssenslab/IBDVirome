####################################################
# SCRIPT 7: PROKARYOTIC VIRUSES: link in silico predicted host with 16S data
####################################################
# 1. Load packages
####################################################
#install.packages("FSA")
library(FSA)
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
library(xfun)

#devtools::install_github("jbisanz/qiime2R", force = TRUE) 

library(qiime2R)
library(phyloseq)
####################################################
# 2. Load previous data
####################################################
# Go to "/Users/daan/Desktop/Bioinformatics/Analysis/HS/Rstudio_analysis/4. Prokaryotic_viruses"
# And load the "Visualization_effect size Rdata".
####################################################
# 3. load your Qiime2 artificacts in R  
##################################################
#https://forum.qiime2.org/t/tutorial-integrating-qiime2-and-r-for-data-visualization-and-analysis-using-qiime2r/4121
#metadatafile needs to have an additional line with types of data of each columns  (categorical or numeric): 
# id             patient         gender
# #q2:types      categorical     categorical

setwd("/Users/daan/Desktop/Bioinformatics/Analysis/HS/Rstudio_analysis/5. Host_Link")
dir()

HS2210_ps <- qza_to_phyloseq(
  features="featuretable-dada2-HS2210-nomock-SRS-12225.qza",
  taxonomy="taxonomy-sequence-HS2210.qza",
  metadata="HS_2210_metadata-nomock.tsv",
  tree = "rootedtree-iqtree-HS2210.qza")
####################################################
# 4. Collapse at taxonomic level and obtain RA of some bacterial genera
##################################################
# Collapse at bacterial genus
HS2210_ps_genus <- aggregate_taxa(HS2210_ps, "Genus", verbose = FALSE)
#View(otu_table(HS2210_ps_genus))

# Extract table
HS2210_ps_genus_matrix <- as.data.frame(otu_table(HS2210_ps_genus))
HS2210_ps_genus_metadata <- as.data.frame(sample_data(HS2210_ps_genus))
Heatmap_genus <- merge(as.data.frame(t(HS2210_ps_genus_matrix)), HS2210_ps_genus_metadata, by=0, all=F)
rownames(Heatmap_genus) <- Heatmap_genus$Row.names
Heatmap_genus$Row.names <- NULL

# Find Vector
vectorA <- which(colnames(Heatmap_genus) == "patient")-1
vector_patient <- colnames(Heatmap_genus[, c(1:vectorA)])

# Aggregated by patient number
Heatmap_genus <- aggregate(. ~patient_short, FUN = sum, data = Heatmap_genus[,colnames(Heatmap_genus) %in% vector_patient | colnames(Heatmap_genus) == 'patient_short'])
rownames(Heatmap_genus) <- Heatmap_genus$patient_short
Heatmap_genus$patient_short <- NULL
colSums(Heatmap_genus)
sum(Heatmap_genus$Corynebacterium)
sum(Heatmap_genus$Staphylococcus)
sum(Heatmap_genus$Pseudomonas)
sum(Heatmap_genus$`Escherichia-Shigella`)

# Calculate RA bacterial genus
Heatmap_genus <- sweep(Heatmap_genus, 2, colSums(Heatmap_genus), '/') 
Heatmap_genus[is.na(Heatmap_genus)] <- 0
Heatmap_genus <- as.data.frame(t(Heatmap_genus[colSums(Heatmap_genus) > 0]))

# Create new patient data matching Lene's and my data
setwd("/Users/daan/Desktop/Bioinformatics/pre-processing/HS project/Metadata/Final")
dir()
metadata_patient <- read_excel("HS_Metadata_Final2_Patient.xlsx")
metadata_patient <- as.data.frame(metadata_patient)
rownames(metadata_patient) <- metadata_patient$Patient_link

# Merge with Patient metadata
Heatmap_genus_ <- t(Heatmap_genus)
Mastertable_genus_bac <- merge(Heatmap_genus_,metadata_patient,by=0,all=F)
rownames(Mastertable_genus_bac) <- Mastertable_genus_bac$Patient_ID
Mastertable_genus_bac$Row.names <- NULL

# Extract Corynebacterium, Staphyloccocus, Finegoldia and Peptoniphilus
Mastertable_genus_bac_4 <- Mastertable_genus_bac[,c("Corynebacterium","Staphylococcus","Finegoldia","Peptoniphilus", "Escherichia-Shigella")]
#View(Mastertable_genus_bac_4)
####################################################
# 5. Obtain RA of cVC (cVC1, cVC3 and cVC4)
##################################################
# Input phyloseq_rarefied_phages cVC

Mastertable_genus_vir <- microbiome::transform(Mastertable_family_like_approach_, "compositional")
Mastertable_genus_vir_abu <- as.data.frame(t(otu_table(Mastertable_genus_vir)))
Mastertable_genus_vir_3 <- Mastertable_genus_vir_abu[,c("family1","family3","family4")]
#View(Mastertable_genus_vir_3)
####################################################
# 6. Merge viral and bacterial taxa
####################################################
Mastertable_genus_microbes <- merge(Mastertable_genus_bac_4, Mastertable_genus_vir_3, by=0, all=F)
rownames(Mastertable_genus_microbes) <- Mastertable_genus_microbes$Row.names
Mastertable_genus_microbes$Row.names <- NULL
View(Mastertable_genus_microbes)
#rownames(metadata_patient) <- metadata_patient$Patient_link
#Mastertable_genus_microbes1 <- merge(Mastertable_genus_microbes, metadata_patient, by=0, all=F)
#rownames(Mastertable_genus_microbes1) <- Mastertable_genus_microbes1$Row.names
#Mastertable_genus_microbes1$Row.names <- NULL

#Mastertable_genus_microbes <- merge(Mastertable_genus_bac_4, Heatmap_genus, by=0, all=F)
#rownames(Mastertable_genus_microbes) <- Mastertable_genus_microbes$Row.names
#Mastertable_genus_microbes$Row.names <- NULL
####################################################
# 7. Figures
####################################################
# esquisser(Mastertable_genus_microbes)

nrow(Mastertable_genus_microbes_cor)
nrow(Mastertable_genus_microbes)

# 1) Corynebacterium
Mastertable_genus_microbes_cor <- Mastertable_genus_microbes[!Mastertable_genus_microbes$family1==0,]
ggplot(Mastertable_genus_microbes_cor) +
 aes(x = Corynebacterium, y = family1) +
 geom_point(shape = "circle", 
 size = 1.5, colour = "#112446") +
 theme_bw() +
  geom_smooth(method='lm')
cor.test(Mastertable_genus_microbes_cor$Corynebacterium, Mastertable_genus_microbes_cor$family1, method="spearman")

# 2) Staphylococcus
Mastertable_genus_microbes_stap <- Mastertable_genus_microbes[!Mastertable_genus_microbes$family3==0,]

ggplot(Mastertable_genus_microbes_stap) +
  aes(x = Staphylococcus, y = family3) +
  geom_point(shape = "circle", 
             size = 1.5, colour = "#112446") +
  theme_bw() +
  geom_smooth(method='lm')
cor.test(Mastertable_genus_microbes_stap$Staphylococcus, Mastertable_genus_microbes_stap$family3, method="spearman")

# 3) Peptoniphilus
Mastertable_genus_microbes_fin_pep <- Mastertable_genus_microbes[!Mastertable_genus_microbes$family4==0,]

ggplot(Mastertable_genus_microbes_fin_pep) +
  aes(x = Peptoniphilus, y = family4) +
  geom_point(shape = "circle", 
             size = 1.5, colour = "#112446") +
  theme_bw() +
  geom_smooth(method='lm')
cor.test(Mastertable_genus_microbes_fin_pep$Peptoniphilus, Mastertable_genus_microbes_fin_pep$family4, method="pearson")

# 4) Finegoldia
ggplot(Mastertable_genus_microbes_fin_pep) +
  aes(x = Finegoldia, y = family4) +
  geom_point(shape = "circle", 
             size = 1.5, colour = "#112446") +
  theme_bw() +
  geom_smooth(method='lm')
cor.test(Mastertable_genus_microbes_fin_pep$Finegoldia, Mastertable_genus_microbes_fin_pep$family4, method="pearson")
####################################################
# 8. Correlate different ones
####################################################
# 8.1 Correlate cV1 ot all bacteria and Escherichia (most abundant one)
####################################################
# 1) Corynebacterium: 
ggplot(Mastertable_genus_microbes) +
  aes(x = Corynebacterium, y = family1) +
  geom_point(shape = "circle", 
             size = 1.5, colour = "#112446") +
  theme_bw() +
  geom_smooth(method='lm', se = FALSE)
cor.test(Mastertable_genus_microbes$Corynebacterium, Mastertable_genus_microbes$family1, method="pearson")

# 2) Corynebacterium
ggplot(Mastertable_genus_microbes) +
  aes(x = Corynebacterium, y = family3) +
  geom_point(shape = "circle", 
             size = 1.5, colour = "#112446") +
  theme_bw() +
  geom_smooth(method='lm', se = FALSE)
cor.test(Mastertable_genus_microbes$Corynebacterium, Mastertable_genus_microbes$family3, method="pearson")

# 3) Corynebacterium
ggplot(Mastertable_genus_microbes) +
  aes(x = Corynebacterium, y = family4) +
  geom_point(shape = "circle", 
             size = 1.5, colour = "#112446") +
  theme_bw() +
  geom_smooth(method='lm', se = FALSE)
cor.test(Mastertable_genus_microbes$Corynebacterium, Mastertable_genus_microbes$family4, method="pearson")

# Multiple testing correction 
P_values <-  c(0.1335,0.1656,0.0159)
P_values
p.adjust(P_values, method = "BH", 3)
####################################################
# 8.2 Correlate Staphylococucus associations 
####################################################
# 1) Staphylococucus 
ggplot(Mastertable_genus_microbes) +
  aes(x = Staphylococcus, y = family3) +
  geom_point(shape = "circle", size = 1.5, colour = "#112446") +
  theme_bw() +
  geom_smooth(method='lm', se = FALSE)+
xlim(0,0.075)
cor.test(Mastertable_genus_microbes$Staphylococcus, Mastertable_genus_microbes$family3, method="pearson")

# 2) Staphylococucus
ggplot(Mastertable_genus_microbes) +
  aes(x = Staphylococcus, y = family1) +
  geom_point(shape = "circle", 
             size = 1.5, colour = "#112446") +
  theme_bw() +
  geom_smooth(method='lm', se = FALSE) +
xlim(0,0.075)
cor.test(Mastertable_genus_microbes$Staphylococcus, Mastertable_genus_microbes$family1, method="pearson")

# 3) Staphylococucus
ggplot(Mastertable_genus_microbes) +
  aes(x = Staphylococcus, y = family4) +
  geom_point(shape = "circle", 
             size = 1.5, colour = "#112446") +
  theme_bw() +
  geom_smooth(method='lm', se = FALSE) +
xlim(0,0.075)
cor.test(Mastertable_genus_microbes$Staphylococcus, Mastertable_genus_microbes$family4, method="pearson")
####################################################
# 8.3 Correlate Finegoldia associations 
####################################################
# 1) Finegoldia 
ggplot(Mastertable_genus_microbes) +
  aes(x = Finegoldia, y = family4) +
  geom_point(shape = "circle", size = 1.5, colour = "#112446") +
  theme_bw() +
  geom_smooth(method='lm', se = FALSE)
cor.test(Mastertable_genus_microbes$Finegoldia, Mastertable_genus_microbes$family4, method="pearson")

# 2) Finegoldia 
ggplot(Mastertable_genus_microbes) +
  aes(x = Finegoldia, y = family1) +
  geom_point(shape = "circle", size = 1.5, colour = "#112446") +
  theme_bw() +
  geom_smooth(method='lm', se = FALSE)
cor.test(Mastertable_genus_microbes$Finegoldia, Mastertable_genus_microbes$family1, method="pearson")

# 3) Finegoldia 
ggplot(Mastertable_genus_microbes) +
  aes(x = Finegoldia, y = family3) +
  geom_point(shape = "circle", size = 1.5, colour = "#112446") +
  theme_bw() +
  geom_smooth(method='lm', se = FALSE)
cor.test(Mastertable_genus_microbes$Finegoldia, Mastertable_genus_microbes$family3, method="pearson")
####################################################
# 8.4 Correlate Peptoniphilus associations 
####################################################
# 1) Finegoldia 
ggplot(Mastertable_genus_microbes) +
  aes(x = Peptoniphilus, y = family4) +
  geom_point(shape = "circle", size = 1.5, colour = "#112446") +
  theme_bw() +
  geom_smooth(method='lm', se = FALSE)
cor.test(Mastertable_genus_microbes$Peptoniphilus, Mastertable_genus_microbes$family4, method="pearson")

# 2) Finegoldia 
ggplot(Mastertable_genus_microbes) +
  aes(x = Peptoniphilus, y = family1) +
  geom_point(shape = "circle", size = 1.5, colour = "#112446") +
  theme_bw() +
  geom_smooth(method='lm', se = FALSE)
cor.test(Mastertable_genus_microbes$Peptoniphilus, Mastertable_genus_microbes$family1, method="pearson")

# 3) Finegoldia 
ggplot(Mastertable_genus_microbes) +
  aes(x = Peptoniphilus, y = family3) +
  geom_point(shape = "circle", size = 1.5, colour = "#112446") +
  theme_bw() +
  geom_smooth(method='lm', se = FALSE)
cor.test(Mastertable_genus_microbes$Peptoniphilus, Mastertable_genus_microbes$family3, method="pearson")
####################################################
# 8.5 Correlate with most abundant bacteria: 'Escherichia-Shigella'
####################################################
# 1) Escherichia-Shigella 

ggplot(Mastertable_genus_microbes) +
  aes(x = `Escherichia-Shigella`, y = family1) +
  geom_point(shape = "circle", size = 1.5, colour = "#112446") +
  theme_bw() +
  geom_smooth(method='lm', se = FALSE)
cor.test(Mastertable_genus_microbes$`Escherichia-Shigella`, Mastertable_genus_microbes$family1, method="pearson")

# 2) Escherichia-Shigella 
ggplot(Mastertable_genus_microbes) +
  aes(x = `Escherichia-Shigella`, y = family3) +
  geom_point(shape = "circle", size = 1.5, colour = "#112446") +
  theme_bw() +
  geom_smooth(method='lm', se = FALSE)
cor.test(Mastertable_genus_microbes$`Escherichia-Shigella`, Mastertable_genus_microbes$family3, method="pearson")

# 3) Escherichia-Shigella 
View(Mastertable_genus_microbes)
ggplot(Mastertable_genus_microbes) +
  aes(x = `Escherichia-Shigella`, y = family4) +
  geom_point(shape = "circle", size = 1.5, colour = "#112446") +
  theme_bw() +
  geom_smooth(method='lm', se = FALSE)
cor.test(Mastertable_genus_microbes$`Escherichia-Shigella`, Mastertable_genus_microbes$family4, method="pearson")
####################################################







