####################################
# SCRIPT 5: PROKARYOTIC VIRUSES
####################################
# 0. Packages: Install packages and load them wherever needed
####################################
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
####################################
# Datasets
####################################
# family-like approach
Mastertable
Mastertable_viral1 <- Mastertable[!Mastertable$Final_viral==1,]
sum(Mastertable_viral1$Totalnumberofreads)
Mastertable_viral1$total<- rowSums(Mastertable_viral1[,1:127])
sort(colSums(Mastertable_viral1[,1:127])) # 127-18= 109 samples kept
sum(Mastertable_viral1$total)
sum(Mastertable_viral1$total)/109

Mastertable_family_like_approach__
Mastertable_family_like_approach_

# Classical approach
Mastertable_classical_approach

View(Mastertable)
####################################
# 0. Subset Viruses
####################################
# Load Mastertable first script
Mastertable$Final_viral[Mastertable$Phages=="phage"] <- "prokaryotic virus"
Mastertable$Final_superkingdom[Mastertable$Phages=="phage"] <- "Viruses"
table(Mastertable$Final_viral)
table(Mastertable$Final_superkingdom)

# Subset phages and eukaryotic viruses
#Mastertable_virome <- Mastertable[Mastertable$Final_superkingdom=="Viruses",]
Mastertable_phages <- Mastertable[Mastertable$Final_viral=="prokaryotic virus",]
#Mastertable_euk_virus <- Mastertable[Mastertable$Final_viral=="eukaryotic virus",]
#table(Mastertable_phages$Family_level)
#View(Mastertable_phages)

## Complete virome
#table(Mastertable_virome$Final_viral)
#Mastertable_virome$Family_level[Mastertable_virome$Final_viral=="eukaryotic virus"] <- Mastertable_virome$Final_family[Mastertable_virome$Final_viral=="eukaryotic virus"]
#table(Mastertable_virome$Family_level)

#Mastertable_virome$Family_level[!Mastertable$Phages=="phage"] <-  Mastertable_virome$Final_family[!Mastertable$Phages=="phage"]
#table(Mastertable_virome$Family_level)
#Mastertable_virome[!Mastertable$Phages=="phage",]
#View(Mastertable_virome[Mastertable$Final_viral=="eukaryotic virus"])

#View(Mastertable_phages)
#View(Mastertable_virome)
View(Mastertable_phages)
#View(Mastertable_euk_virus)
####################################
# 1. Adding higher taxonomies to mastertable (genus/family)
###################################
# 1.1 Family-like groups 
###################################
#setwd("/Users/daan/Desktop/Transfer/Input_R/Taxonomy/")
#getwd()
#dir()

#AllSamples_classification_family_groups <- as.data.frame(read_delim("family_clusters5.txt", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE))[-1,]
#colnames(AllSamples_classification_family_groups) <- c('contig', 'Family')
#rownames(AllSamples_classification_family_groups) <- AllSamples_classification_family_groups$contig
#AllSamples_classification_family_groups$contig <- NULL
#rownames(AllSamples_classification_family_groups) <- gsub("\\.{1}", "_", rownames(AllSamples_classification_family_groups))
###################################
# 1.2 Merge 
###################################
#Mastertable_viral_rarefied_1 <- merge(Mastertable_virome, AllSamples_classification_family_groups,by=0, all=F)
#rownames(Mastertable_viral_rarefied_1) <- Mastertable_viral_rarefied_1$Row.names
#Mastertable_viral_rarefied_1$Row.names <- NULL
#View(Mastertable_viral_rarefied_1)
#table(Mastertable_viral_rarefied_1$Final_viral)
###################################
# 1.3 Select prokaryotic viruses
###################################
#table(Mastertable_viral_rarefied_1$Final_viral)
#Mastertable_viral_rarefied_2 <- Mastertable_viral_rarefied_1[Mastertable_viral_rarefied_1$Final_viral == "prokaryotic virus",]
#unique(Mastertable_viral_rarefied_2$Family) # 126 family-like phages
###################################
# 2. Upload Metadata
####################################
# 2.1  Upload Metadata per sample
#################################### 
setwd("/Users/daan/Desktop/Bioinformatics/pre-processing/HS project/Metadata/Final")
dir()
metadata <- read_excel("HS_Metadata_Final2_R.xlsx")
metadata <- as.data.frame(metadata)
rownames(metadata) <- metadata$Patient_ID
metadata$Age <- as.numeric(metadata$Age)
metadata$Hurley <- as.character(metadata$Hurley)
metadata$`Genetics_family HS` <- as.character(metadata$`Genetics_family HS`)
metadata$`Genetics_family skin diseases` <- as.character(metadata$`Genetics_family skin diseases`)
metadata$Roken <- as.character(metadata$Roken)
#metadata <- metadata[,c(1,2,3,12)]
colnames(metadata)
#################################### 
# 2.2 Upload Metadata per patient
#################################### 
setwd("/Users/daan/Desktop/Bioinformatics/pre-processing/HS project/Metadata/Final")
dir()

metadata_patient <- read_excel("HS_Metadata_Final2_Patient.xlsx")
metadata_patient <- as.data.frame(metadata_patient)
rownames(metadata_patient) <- metadata_patient$Patient_ID
metadata_patient$Age <- as.numeric(metadata_patient$Age)
metadata_patient$Hurley <- as.character(metadata_patient$Hurley)
metadata_patient$`Genetics_family HS` <- as.character(metadata_patient$`Genetics_family HS`)
metadata_patient$`Genetics_family skin diseases` <- as.character(metadata_patient$`Genetics_family skin diseases`)
metadata_patient$Roken <- as.character(metadata_patient$Roken)
#metadata_patient <- metadata_patient[!metadata_patient$Obese=="Underweight",]
#metadata_patient <- metadata_patient[,c(1,2,8)]

#View(metadata_patient)

ggplot(metadata_patient) +
 aes(x = Hurley, fill = Roken_2) +
 geom_bar() +
 scale_fill_hue(direction = 1) +
 theme_minimal()

ggplot(metadata_patient) +
 aes(x = Hurley, fill = Roken) +
 geom_bar(position = "dodge") +
 scale_fill_hue(direction = 1) +
 theme_minimal()

ggplot(metadata_patient) +
  aes(x = Hurley, fill = Obese) +
  geom_bar(position = "dodge") +
  scale_fill_hue(direction = 1) +
  theme_minimal()


# Statistics summary table
str(Metadata)
sum(Mastertable$Totalnumberofreads)

# 1) # Age = 0.03686
wilcox.test(Age~Disease, data =metadata_patient,alternative = "two.sided",  paired = FALSE, exact = T)
# Padj = 0.0829350

# 2) Ethnicity = 0.07931
chisq.test(metadata_patient$Disease, metadata_patient$Ethnicity)
# Padj = 0.1189650

# 3) Smoking = 0.0002405
metadata_patient$Roken[metadata_patient$Roken=="2"] <- "0"
chisq.test(metadata_patient$Disease, metadata_patient$Roken)
# Padj = 0.0021645

# 4) IBD = 0.07685
metadata_patient$IBD <- as.character(metadata_patient$IBD)
chisq.test(metadata_patient$Disease, metadata_patient$IBD)
# Padj = 0.1189650

# 5) Genetics_family.HS = 0.3946
chisq.test(metadata_patient$Disease, metadata_patient$`Genetics_family HS`)
# Padj = 0.4439250

# 6) Genetics_family.skin.diseases = 0.09645
chisq.test(metadata_patient$Disease, metadata_patient$`Genetics_family skin diseases`)
# Padj = 0.1240071

# 7) Sex = 1
chisq.test(metadata_patient$Disease, metadata_patient$Geslacht)
# Padj = 1.0000000

# 8) Categorical BMI = 0.003896 
metadata_patient_Obese <- metadata_patient[!metadata_patient$Obese=="NA",]
chisq.test(metadata_patient_Obese$Disease, metadata_patient_Obese$Obese)
# Padj = 0.0175320

# 9) Disease location = 0.02834
# Normally you need the mcnemar test to compare paired data, like here, which is more powerfull
# but since cases and controls are perfectaly matched, you will have partially matched data leading to lots of exclusions to match it.
# Therefore, simply use unpaired chisquared data.
chisq.test(metadata$Disease, metadata$Place)
# Padj = 0.0829350

# Multiple testing correction 
P_values <-  c(0.03686,0.07931,0.0002405,0.07685,0.3946,0.09645,1,0.003896,0.02834)
P_values
p.adjust(P_values, method = "BH", 9)

View(Mastertable[Mastertable$Family_level=="family1",])
####################################
# 3. Subset prokaryotic viruses
####################################
Mastertable_viral_rarefied_2 <- Mastertable_family_like_approach_
colnames(Mastertable_viral_rarefied_2)

Mastertable_viral_rarefied_2 <- Mastertable_phage
table(Mastertable_viral_rarefied_2$Final_viral)
#Mastertable_viral_rarefied_2 <- Mastertable_contigs_viral[Mastertable_contigs_viral$Final_viral=="prokaryotic virus",]
Mastertable_viral_rarefied_2

# Select phages
table(Mastertable_viral_rarefied_2$Final_viral)
sort(colSums(Mastertable_viral_rarefied_2[1:127]))

names <- colnames(Mastertable_viral_rarefied_2[1:99][,colSums(Mastertable_viral_rarefied_2[1:99])!=0])
names <- colnames(Mastertable_viral_rarefied_2[1:127][,colSums(Mastertable_viral_rarefied_2[1:127])!=0])

# Keep ones with more than 10000 reads
names # 117 samples left with entire virome (only 10 discarded)
names2 <- names

# Select samples with phages present and create abundance table
#Mastertable_viral_rarefied_2$NODES <- rownames(Mastertable_viral_rarefied_2)
abundance_table <- aggregate(. ~Family_level, FUN = sum, data = Mastertable_viral_rarefied_2[,colnames(Mastertable_viral_rarefied_2) %in% names | colnames(Mastertable_viral_rarefied_2) == 'Family_level'])
rownames(abundance_table) <- abundance_table$Family_level
abundance_table$Family_level <- NULL
abundance_table_samples <- abundance_table

# Create taxonomy table
Mastertable_viral_rarefied_2$Final_ANI <- NULL
Mastertable_viral_rarefied_2$Final_coverage <- NULL
Mastertable_viral_rarefied_2$Blastn_AS <- NULL
Mastertable_viral_rarefied_2$Diamond_AS <- NULL
Mastertable_viral_rarefied_2$Best_AS <- NULL
Mastertable_viral_rarefied_2$Final_subfamily <- NULL
Mastertable_phages$Final_species <- NULL
Mastertable_viral_rarefied_2$Final_family <- NULL
Mastertable_viral_rarefied_2$Final_class_temp <- NULL
Mastertable_viral_rarefied_2$Phages <- NULL
Mastertable_viral_rarefied_2$Family <- Mastertable_viral_rarefied_2$Family_level
vector_1 <- (which(names(Mastertable_viral_rarefied_2)== "Final_superkingdom"))
vector_2 <- which(names(Mastertable_viral_rarefied_2)== "Family")
Taxonomy_table <- Mastertable_viral_rarefied_2[,c(vector_1:vector_2)]
#Taxonomy_table$NODES <- NULL
View(Taxonomy_table)

# Rename some Taxonomy elements
names(Taxonomy_table)[names(Taxonomy_table) == "Final_superkingdom"] <- "Kingdom"
names(Taxonomy_table)[names(Taxonomy_table) == "Fina_phylum"] <- "Phylum"
names(Taxonomy_table)[names(Taxonomy_table) == "Final_class"] <- "Class"
names(Taxonomy_table)[names(Taxonomy_table) == "Final_order"] <- "Order"
names(Taxonomy_table)[names(Taxonomy_table) == "Final_genus"] <- "Genus"
names(Taxonomy_table)[names(Taxonomy_table) == "Family"] <- "Family"
names(Taxonomy_table)[names(Taxonomy_table) == "NODES"] <- "Species"
Taxonomy_table$NODES <- NULL
# Taxonomy_table_VC
Taxonomy_table_VC

names(Taxonomy_table)[names(Taxonomy_table) == "Final_species"] <- "Species"
Taxonomy_table <- Taxonomy_table[,c(7,1, 2, 3,4,5,6)]
View(Taxonomy_table)
####################################
# (3) Alternative 3
####################################
# (3.1) Alternative Aggregate by patient
####################################
Mastertable_family_like_patient <- Mastertable_family_like_approach__
Mastertable_family_like_patient <- Mastertable_family_like_approach_agg
Mastertable_family_like_approach_agg$totalreads <- NULL
#View(Mastertable_family_like_approach__)

# Merge with metadata samples
Mastertable_family_like_patient_pcoa <- t(Mastertable_family_like_patient)
Mastertable_family_like_patient_pcoa <- merge(Mastertable_family_like_patient_pcoa, metadata, by=0, all =F)
rownames(Mastertable_family_like_patient_pcoa) <- Mastertable_family_like_patient_pcoa$Row.names
Mastertable_family_like_patient_pcoa$Row.names <- NULL

# Find Vector
vectorA <- which(colnames(Mastertable_family_like_patient_pcoa) == "Patient_ID")-1
vector_patient <- colnames(Mastertable_family_like_patient_pcoa[, c(1:vectorA)])

# Aggregated by patient number
#View(Mastertable_family_like_patient_pcoa)
Mastertable_family_like_patient_pcoa_ <- aggregate(. ~HS, FUN = sum, data = Mastertable_family_like_patient_pcoa[,colnames(Mastertable_family_like_patient_pcoa) %in% vector_patient | colnames(Mastertable_family_like_patient_pcoa) == 'HS'])
rownames(Mastertable_family_like_patient_pcoa_) <- Mastertable_family_like_patient_pcoa_$HS
Mastertable_family_like_patient_pcoa_$HS <- NULL

# Remove column with sum zero 
Mastertable_family_like_patient_pcoa_ <-  as.data.frame(t(Mastertable_family_like_patient_pcoa_))
Mastertable_family_like_patient_pcoa_ <- Mastertable_family_like_patient_pcoa_[, which(colSums(Mastertable_family_like_patient_pcoa_) != 0)]
abundance_table <- Mastertable_family_like_patient_pcoa_
#View(abundance_table)
####################################
# (3.2) Prepare for Phyloseq
####################################
#names <- colnames(Mastertable_family_like_patient_pcoa_[1:51][,colSums(Mastertable_family_like_patient_pcoa_[1:51])!=0])
names # 51 patients left with family-like approach (adapt with viral identification > 3kb)
names <- colnames(Mastertable_viral_rarefied_2[1:127][,colSums(Mastertable_viral_rarefied_2[1:127])!=0])
names

# Create taxonomy table
Mastertable_family_like_patient_pcoa_$Family <- rownames(Mastertable_family_like_patient_pcoa_)
Taxonomy_table <- Mastertable_family_like_patient_pcoa_[52]
Taxonomy_table$Kingdom <- "Unannotated"
Taxonomy_table$Phylum <- "Unannotated"
Taxonomy_table$Class <- "Unannotated"
Taxonomy_table$Order <- "Unannotated"
Taxonomy_table$Genus <- "Unannotated"
Taxonomy_table$Species <- "Unannotated"
Taxonomy_table_VC <- Taxonomy_table
View(Taxonomy_table)
####################################
# 4. Phyloseq tables
####################################
# 4.1 Samples
metadata
abundance_table_samples
Taxonomy_table_VC

# 4.2 Patients
metadata_patient
abundance_table
Taxonomy_table
####################################
# 5. Create phyloSeq objects
####################################
# 5.1 Create a matrix of abundance and taxonomy tables
####################################
# Samples
abundance_table_rarefied_m_samples <- as.matrix(abundance_table_samples)
taxonomy_table_rarefied_m_samples <- as.matrix(Taxonomy_table_VC)

# patients
abundance_table_rarefied_m <- as.matrix(abundance_table)
taxonomy_table_rarefied_m <- as.matrix(Taxonomy_table)
#View(abundance_table_rarefied_m)
####################################
# 5.2 Transform to phylseq objects
####################################
# Samples
ABUNDANCE_rarefied_family <- otu_table(abundance_table_rarefied_m_samples, taxa_are_rows = TRUE)
TAX_rarefied_family <- tax_table(taxonomy_table_rarefied_m_samples)
samples <- sample_data(metadata)
phyloseq_rarefied_family <- phyloseq(ABUNDANCE_rarefied_family,TAX_rarefied_family, samples)

# Patient
ABUNDANCE_rarefied_family <- otu_table(abundance_table_rarefied_m, taxa_are_rows = TRUE)
TAX_rarefied_family <- tax_table(taxonomy_table_rarefied_m)
samples <- sample_data(metadata_patient)
phyloseq_rarefied_family <- phyloseq(ABUNDANCE_rarefied_family,TAX_rarefied_family, samples)
####################################
# 5.3 Visualize data & subset only the phage fraction
####################################
sample_names(phyloseq_rarefied_family) # All sample names
rank_names(phyloseq_rarefied_family) # All taxonomies
sample_variables(phyloseq_rarefied_family) ## All metadata
# if you want to subset the dataframe into only "responders", just do "responder_data <- subset_samples(phyloseq_unrarefied_family, response == "responder")"
#phyloseq_rarefied_phages <- subset_taxa(phyloseq_rarefied_family, Kingdom == "Viruses")

phyloseq_rarefied_phages <- phyloseq_rarefied_family
####################################
# 5.4 PhyloSeq format
####################################
phyloseq_rarefied_phages
#View(sample_data(phyloseq_rarefied_phages))
####################################
# 5.5 PcoA Visualization
####################################
# 5.5.1 Plot compositional figure
####################################
#install.packages("hrbrthemes")
library(hrbrthemes)
#install.packages("gcookbook")
library(gcookbook)
library(microbiome)

# A) All samples: Disease and obesity
ps1.rel <- microbiome::transform(phyloseq_rarefied_phages, "compositional")
ps1.fam.rel <-aggregate_rare(ps1.rel, level = "Family", detection = 0.005, prevalence = 0.25) # 20% presence
sample_data(ps1.fam.rel)$Obese <- factor(sample_data(ps1.fam.rel)$Obese, levels = c("Underweight","Normal","Overweight","Obese"))
ps1.fam.rel = subset_samples(ps1.fam.rel, Obese= "Norma")

plot_composition(ps1.fam.rel,average_by = "Disease") + 
  guides(fill = guide_legend(ncol = 1)) +
  labs(x = "Samples", 
       y = "Relative abundance",
       title = "Relative abundance data", 
       subtitle = "Subtitle",
       caption = "Caption text.")  +
       scale_fill_brewer("Family", palette = "Paired") +
        theme_bw() +

plot_composition(ps1.fam.rel,average_by = "Obese") + 
  guides(fill = guide_legend(ncol = 1)) +
  labs(x = "Samples", 
       y = "Relative abundance",
       title = "Relative abundance data", 
       subtitle = "Subtitle",
       caption = "Caption text.")  + 
  scale_fill_brewer("Family", palette = "Paired") + 
  theme_bw() 

plot_composition(ps1.fam.rel,average_by = "Hurley") + 
  guides(fill = guide_legend(ncol = 1)) +
  labs(x = "Samples", 
       y = "Relative abundance",
       title = "Relative abundance data", 
       subtitle = "Subtitle",
       caption = "Caption text.")  + 
  scale_fill_brewer("Family", palette = "Paired") + 
  theme_bw() 

ps1.fam.rel %>%
  group_by(Disease) %>%
  summarise(mean_run = mean(Disease))

# B) HS patients: HS stages and obesity.
ps1.rel <- microbiome::transform(phyloseq_rarefied_phages, "compositional")
ps1.fam.rel <-aggregate_rare(ps1.rel, level = "Family", detection = 0.005, prevalence = 0.25) # 20% presence
ps1.fam.rel_HS <- subset_samples(ps1.fam.rel, Disease == "HS")
sample_data(ps1.fam.rel_HS)$Obese <- factor(sample_data(ps1.fam.rel_HS)$Obese, levels = c("Underweight","Normal","Overweight","Obese"))

plot_composition(ps1.fam.rel_HS,
                 average_by = "Disease") + 
  guides(fill = guide_legend(ncol = 1)) +
  labs(x = "Samples", 
       y = "Relative abundance",
       title = "Relative abundance data", 
       subtitle = "Subtitle",
       caption = "Caption text.")  +
  scale_fill_brewer("Family", palette = "Paired") +
  theme_bw() +
  
  plot_composition(ps1.fam.rel_HS,average_by = "Obese") + 
  guides(fill = guide_legend(ncol = 1)) +
  labs(x = "Samples", 
       y = "Relative abundance",
       title = "Relative abundance data", 
       subtitle = "Subtitle",
       caption = "Caption text.")  + 
  scale_fill_brewer("Family", palette = "Paired") + 
  theme_bw()

# C) Normal BMI
library("microbiome")
ps1.fam.rel = subset_samples(phyloseq_rarefied_phages, Obese=="Normal")
ps1.rel <- microbiome::transform(ps1.fam.rel, "compositional")
ps1.fam.rel <-aggregate_rare(ps1.rel, level = "Family", detection = 0.01, prevalence = 0.25) # 20% presence

plot_composition(ps1.fam.rel,average_by = "Disease") + 
  guides(fill = guide_legend(ncol = 1)) +
  labs(x = "Samples", 
       y = "Relative abundance",
       title = "Relative abundance data", 
       subtitle = "Subtitle",
       caption = "Caption text.")  +
  scale_fill_brewer("Family", palette = "Paired") +
  theme_bw() 

plot_composition(ps1.fam.rel,
                 average_by = "Hurley") + 
  guides(fill = guide_legend(ncol = 1)) +
  labs(x = "Samples", 
       y = "Relative abundance",
       title = "Relative abundance data", 
       subtitle = "Subtitle",
       caption = "Caption text.")  +
  scale_fill_brewer("Family", palette = "Paired") +
  theme_bw() 

# Conclusion: family1  and family 3 show in the least-confounded dataset also a lower abundance of famil1 and family3.
  
####################################
# 5.5.2 Plot PcoA
####################################
# Rarefy
#phyloseq_rarefied_phages = rarefy_even_depth(phyloseq_rarefied_phages, rngseed=1, sample.size=10000, replace=F)

# Remove "Others
#phyloseq_rarefied_phages_visualization <- subset_taxa(phyloseq_rarefied_phages_visualization, !Family == "Other")

# Remove samples (colsums = 0)

# 1) Aggregate families: cutoff 0.01% abundance, 25% prevalence (top 21 shared phages)
phyloseq_rarefied_phages_visualization <- phyloseq_rarefied_phages %>% 
  microbiome::transform(transform = "identity") %>%
  aggregate_rare(level = "Family", detection = 0.005/100, prevalence = 20/100)

phyloseq_rarefied_phages_visualization <- phyloseq_rarefied_phages
#phyloseq_rarefied_phages_visualization <- subset_samples(phyloseq_rarefied_phages_visualization, Disease == "HS")
#phyloseq_rarefied_phages_visualization <- subset_taxa(phyloseq_rarefied_phages_visualization, !Family == "Other")

#phyloseq_rarefied_phages_visualization <- prune_samples(sample_sums(phyloseq_rarefied_phages_visualization) > 10000, phyloseq_rarefied_phages_visualization)
#phyloseq_rarefied_phages_visualization <- subset_samples(phyloseq_rarefied_phages_visualization, Disease == "HS")
#View(otu_table(phyloseq_rarefied_phages_visualization))
#View(sample_data(phyloseq_rarefied_phages_visualization))

# 2) Aggregate families
# Perform (1) RA calculation & (2) Hellinger transformation
#phyloseq_rarefied_phages_transformed = transform_sample_counts(phyloseq_rarefied_phages_visualization, function(x) log(10 + x))
#phyloseq_rarefied_phages_transformed <- phyloseq_rarefied_phages_visualization %>%
#microbiome::transform(transform = "compositional") 

phyloseq_rarefied_phages_transformed <- phyloseq_rarefied_phages_visualization
phyloseq_rarefied_phages_transformed_ <- phyloseq_rarefied_phages_transformed %>%
microbiome::transform(transform = "hellinger")
phyloseq_rarefied_phages_transformed <- phyloseq_rarefied_phages_transformed_
#View(otu_table(phyloseq_rarefied_phages_transformed))

phyloseq_rarefied_phages_transformed_bray <- phyloseq::distance(phyloseq_rarefied_phages_transformed, method="bray")
phyloseq_rarefied_phages_PcoA <- ordinate(phyloseq_rarefied_phages_transformed, method="PCoA", distance=phyloseq_rarefied_phages_transformed_bray)

# order metada variables
metadata_Pcoa <- as.data.frame(sample_data(phyloseq_rarefied_phages_transformed))
metadata_Pcoa$Disease <- factor(metadata_Pcoa$Disease, levels = c("HS", "HE"))
phyloseq_rarefied_phages_transformed <- merge_phyloseq(metadata_Pcoa, phyloseq_rarefied_phages_transformed)
colnames(metadata_Pcoa)

# Visualization PcoA 1:
plot_ordination(phyloseq_rarefied_phages_transformed, phyloseq_rarefied_phages_PcoA, color="Disease", shape = "Obese") +
  geom_point(size = 1, alpha=1) +
  theme_bw() +
  scale_color_manual(values = c(`HS` = "#bd7f1c", `HE` = "#2a8a92")) +
  scale_shape_manual(values = c(16, 3, 2, 0)) +
  geom_hline(yintercept = 0, linetype = 3, colour = "black") + 
  geom_vline(xintercept = 0, linetype = 3, colour = "black") +
  xlab("PC2 (19.6)") +
  ylab("PC1 (10.5)")

# PcoA 2: with ellipes: Disease
plot_ordination(phyloseq_rarefied_phages_transformed, phyloseq_rarefied_phages_PcoA, color="Disease") + #,shape = "Endoscopy.outcome"#
  geom_point(size=1, alpha=0.75)  +
  theme_bw() +
  scale_shape_manual(values = c(21,16, 2)) + ## depends on Therapy_2 number of variables, 2 in this case
  scale_color_manual(values = c("#bd7f1c","#2a8a92")) +
  geom_hline(yintercept = 0, linetype = 3, colour = "black") + 
  geom_vline(xintercept = 0, linetype = 3, colour = "black") +
  ylab("PC1 (10.5)") +
  xlab("PC2 (19.6)") +
  stat_ellipse(type='t',size =0.8, linetype = 2,level = 0.80,aes(color = Disease)) # 80% ellipses.

# Hurley stages
plot_ordination(phyloseq_rarefied_phages_transformed, phyloseq_rarefied_phages_PcoA, color="Hurley") + #,shape = "Endoscopy.outcome"#
  geom_point(size=2, alpha=0.75)  +
  theme_classic() +
  scale_shape_manual(values = c(21,16, 2)) + ## depends on Therapy_2 number of variables, 2 in this case
 # scale_color_manual(values = c(`0` = "brown", `1` = "grey",`2` = "grey",`3` = "#2a8a92")) +
  geom_hline(yintercept = 0, linetype = 3, colour = "black") + 
  geom_vline(xintercept = 0, linetype = 3, colour = "black") +
  xlab("PC1 (10.6)") +
  ylab("PC2 (17.5)") +
  theme(axis.text.y = element_text(colour = "black", size = 9, face = "bold"),
        axis.text.x = element_text(colour = "black", face = "bold", size = 9),
        legend.text = element_text(size = 9, face = "bold", colour = "black"),
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14),axis.title.x = element_text(face = "bold", size = 14, colour = "black"),legend.title = element_blank(), 
    panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2)) +
  stat_ellipse(type='t',size =1, linetype = 2,level = 0.95, data = . %>% filter(Hurley != 1) %>% filter(Hurley != 2), aes(color = Hurley))

# Roken 
plot_ordination(phyloseq_rarefied_phages_transformed, phyloseq_rarefied_phages_PcoA, color="Roken") + #,shape = "Endoscopy.outcome"#
  geom_point(size=2, alpha=0.75)  +
  theme_classic() +
  scale_shape_manual(values = c(21,16, 2)) + ## depends on Therapy_2 number of variables, 2 in this case
  # scale_color_manual(values = c(`0` = "brown", `1` = "grey",`2` = "grey",`3` = "#2a8a92")) +
  geom_hline(yintercept = 0, linetype = 3, colour = "black") + 
  geom_vline(xintercept = 0, linetype = 3, colour = "black") 
  xlab("PC1 (10.6)") +
  ylab("PC2 (17.5)") +
  theme(axis.text.y = element_text(colour = "black", size = 9, face = "bold"),
        axis.text.x = element_text(colour = "black", face = "bold", size = 9),
        legend.text = element_text(size = 9, face = "bold", colour = "black"),
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14),axis.title.x = element_text(face = "bold", size = 14, colour = "black"),legend.title = element_blank(), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2)) +
  stat_ellipse(type='t',size =1, linetype = 2,level = 0.95, data = . %>% filter(Hurley != 1) %>% filter(Hurley != 2), aes(color = Hurley))


# PcoA 2: with ellipes: Place
plot_ordination(phyloseq_rarefied_phages_transformed, phyloseq_rarefied_phages_PcoA, color="Place") + 
  geom_point(size=2, alpha=0.75)  +
  theme_classic() +
  scale_color_manual(values = c("Chest" = "white", "Abdominal" = "white","Axilliary" = "#2a8a92","Gluteaal" = "brown","Groin and genitals" = "#bd7f1c")) +
  geom_hline(yintercept = 0, linetype = 3, colour = "black") + 
  geom_vline(xintercept = 0, linetype = 3, colour = "black") +
  xlab("PC1 (17.2)") +
  ylab("PC2 (22.7)") +
  theme(axis.text.y = element_text(colour = "black", size = 9, face = "bold"),
        axis.text.x = element_text(colour = "black", face = "bold", size = 9),
        legend.text = element_text(size = 9, face = "bold", colour = "black"),
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14),
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"),
        legend.title = element_blank(), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2)) 
#stat_ellipse(type='t',size =1, linetype = 2,level = 0.95, data = . %>% filter(Place != "Chest") %>% filter(Place != "Abdominal"),aes(color = Place)) 

#  Obese
plot_ordination(phyloseq_rarefied_phages_transformed, phyloseq_rarefied_phages_PcoA, color="Province") + #,shape = "Endoscopy.outcome"#
  geom_point(size=2, alpha=0.75)  +
  theme_classic() +
  scale_color_manual(values = c("Normal" = "brown", "Underweight" = "grey","Overweight" = "#bd7f1c","Obese" = "#2a8a92")) +
  geom_hline(yintercept = 0, linetype = 3, colour = "black") + 
  geom_vline(xintercept = 0, linetype = 3, colour = "black") +
  xlab("PC1 (9.1)") +
  ylab("PC2 (13.3)") +
  theme(axis.text.y = element_text(colour = "black", size = 9, face = "bold"),
        axis.text.x = element_text(colour = "black", face = "bold", size = 9),
        legend.text = element_text(size = 9, face = "bold", colour = "black"),
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14),axis.title.x = element_text(face = "bold", size = 14, colour = "black"),legend.title = element_blank(), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2)) 
  #stat_ellipse(type='t',size =1, linetype = 2,level = 0.95, data = . %>% filter(Obese != "Underweight"),aes(color = Obese)) 

plot_ordination(phyloseq_rarefied_phages_transformed, phyloseq_rarefied_phages_PcoA, color="Obese") +
  geom_point(size = 1, alpha=1) +
  scale_shape_manual(values = c(16,17)) +
  theme_bw() +
  geom_hline(yintercept = 0, linetype = 3, colour = "black") + 
  geom_vline(xintercept = 0, linetype = 3, colour = "black") +
  xlab("PC1 (10.5)") +
  ylab("PC2 (19.6)") 

# Environment
plot_ordination(phyloseq_rarefied_phages_transformed, phyloseq_rarefied_phages_PcoA, color="Environment") + #,shape = "Endoscopy.outcome"#
  geom_point(size=2, alpha=0.75)  +
  theme_classic() +
 # scale_color_manual(values = c("Normal" = "brown", "Underweight" = "grey","Overweight" = "#bd7f1c","Obese" = "#2a8a92")) +
  geom_hline(yintercept = 0, linetype = 3, colour = "black") + 
  geom_vline(xintercept = 0, linetype = 3, colour = "black") +  
  xlab("PC1 (9.1)") +
  ylab("PC2 (13.3)") +
  theme(axis.text.y = element_text(colour = "black", size = 9, face = "bold"),
        axis.text.x = element_text(colour = "black", face = "bold", size = 9),
        legend.text = element_text(size = 9, face = "bold", colour = "black"),
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14),axis.title.x = element_text(face = "bold", size = 14, colour = "black"),legend.title = element_blank(), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2)) 
 stat_ellipse(type='t',size =1, linetype = 2,level = 0.95,aes(group = Environment))

 plot_ordination(phyloseq_rarefied_phages_transformed, phyloseq_rarefied_phages_PcoA, color="Roken_2") + #,shape = "Endoscopy.outcome"#
   geom_point(size=2, alpha=0.75)  +
   theme_classic() +
   # scale_color_manual(values = c("Normal" = "brown", "Underweight" = "grey","Overweight" = "#bd7f1c","Obese" = "#2a8a92")) +
   geom_hline(yintercept = 0, linetype = 3, colour = "black") + 
   geom_vline(xintercept = 0, linetype = 3, colour = "black") +  
   xlab("PC1 (9.1)") +
   ylab("PC2 (13.3)") +
   theme(axis.text.y = element_text(colour = "black", size = 9, face = "bold"),
         axis.text.x = element_text(colour = "black", face = "bold", size = 9),
         legend.text = element_text(size = 9, face = "bold", colour = "black"),
         legend.position = "right", axis.title.y = element_text(face = "bold", size = 14),axis.title.x = element_text(face = "bold", size = 14, colour = "black"),legend.title = element_blank(), 
         panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2)) +
 stat_ellipse(type='t',size =1, linetype = 2,level = 0.95,aes(group = Roken_2))
####################################
# 5.5.3 Envfit
####################################
library(vegan)

# Select metadata
colnames(metadata_patient)
metadata_patient_1 <- metadata_patient[c(2,10)]

# Create DF 
scrs <- as.data.frame(phyloseq_rarefied_phages_PcoA$vectors[,1:2])

# Combine with metadata
scrs <- merge(scrs,metadata_patient_1,by=0,all=F)
rownames(scrs) <- scrs$Row.names
scrs$Row.names <- NULL

# Envfit 
set.seed(123)
phyloseq_rarefied_phages_transformed_1 <- sample_data(phyloseq_rarefied_phages_transformed) 
vf <- envfit(phyloseq_rarefied_phages_PcoA$vectors, phyloseq_rarefied_phages_transformed_1, perm = 999)
vf # Several signifcant association variables to ordination matrix (first two PC scores)

# Select factors you want to plot
spp.scrs <- as.data.frame(scores(vf, display = "factors"))
spp.scrs$Group <- rownames(spp.scrs)
spp.scrs <- spp.scrs[spp.scrs$Group=="ObeseNormal" |
                 spp.scrs$Group=="ObeseObese" |
                 spp.scrs$Group=="ObeseOverweight" |
                 spp.scrs$Group=="ObeseUnderweight" | 
                 spp.scrs$Group=="DiseaseHS" |
                 spp.scrs$Group=="DiseaseHE",]

# Figure with arrows
ggplot(scrs) +
  geom_point(mapping = aes(x = Axis.1, y = Axis.2, colour = Disease, shape=Obese)) +
  coord_fixed() +
  theme_bw() +
  scale_color_manual(values = c(`HS` = "#bd7f1c", `HE` = "#2a8a92")) +
  scale_shape_manual(values = c(16, 3, 2, 0)) +
  geom_hline(yintercept = 0, linetype = 3, colour = "black") + 
  geom_vline(xintercept = 0, linetype = 3, colour = "black") +
  xlab("PC2 (19.6)") +
  ylab("PC1 (10.5)") +
 geom_segment(data = spp.scrs,
               aes(x = 0, xend = Axis.1, y = 0, yend = Axis.2),
               arrow = arrow(length = unit(0.25, "cm")), colour = "black") + 
    geom_text(data = spp.scrs, aes(x = Axis.1, y = Axis.2, label = Group, colour="grey"),
              size = 3)
# PS: you can adapt size en shape of arrows best in illustrator.
####################################
# 6. Linear regression analysis
####################################
# 6.1 Univariate regression analysis
####################################
Metadata <- data.frame(sample_data(phyloseq_rarefied_phages_transformed))
Metadata$Roken_1 <- "1"
Metadata$Roken_1[Metadata$Roken=="0" | Metadata$Roken=="2"] <- "0"
Metadata$Hurley
#Metadata$Patient_ID <- NULL
Metadata$IBD <- as.character(Metadata$IBD)
colnames(Metadata)
str(Metadata)
#View(Metadata)

# Disease
dbRDA_0 <- capscale(phyloseq_rarefied_phages_transformed_bray ~ 1, data = Metadata)
dbRDA_1 <- capscale(phyloseq_rarefied_phages_transformed_bray ~  Disease, data = Metadata)
step.res_cum_Disease <- ordiR2step(dbRDA_0, dbRDA_1, direction="forward", R2scope = TRUE)
RsquareAdj(step.res_cum_Disease)$r.squared
RsquareAdj(step.res_cum_Disease)$adj.r.squared
step.res_cum_Disease$anova
step.res_cum_Disease <- step.res_cum_Disease$anova[1,1:5]

# Time
dbRDA_0 <- capscale(phyloseq_rarefied_phages_transformed_bray ~ 1, data = Metadata)
dbRDA_1 <- capscale(phyloseq_rarefied_phages_transformed_bray ~  Time, data = Metadata)
step.res_cum_time <- ordiR2step(dbRDA_0, dbRDA_1, direction="forward", R2scope = TRUE)
RsquareAdj(step.res_cum_time)$r.squared
RsquareAdj(step.res_cum_time)$adj.r.squared
step.res_cum_time$anova
step.res_cum_time <- step.res_cum_time$anova[1,1:5]

# Place
dbRDA_0 <- capscale(phyloseq_rarefied_phages_transformed_bray ~ 1, data = Metadata)
dbRDA_1 <- capscale(phyloseq_rarefied_phages_transformed_bray ~  Place, data = Metadata)
step.res_cum_place <- ordiR2step(dbRDA_0, dbRDA_1, direction="forward", R2scope = TRUE)
RsquareAdj(step.res_cum_place)$r.squared
RsquareAdj(step.res_cum_place)$adj.r.squared
step.res_cum_place$anova
step.res_cum_place <- step.res_cum_place$anova[1,1:5]

# Ethnicity
dbRDA_0 <- capscale(phyloseq_rarefied_phages_transformed_bray ~ 1, data = Metadata)
dbRDA_1 <- capscale(phyloseq_rarefied_phages_transformed_bray ~  Ethnicity, data = Metadata)
step.res_cum_ethncity <- ordiR2step(dbRDA_0, dbRDA_1, direction="forward", R2scope = TRUE)
RsquareAdj(step.res_cum_ethncity)$r.squared
RsquareAdj(step.res_cum_ethncity)$adj.r.squared
step.res_cum_ethncity$anova
step.res_cum_ethncity <- step.res_cum_ethncity$anova[1,1:5]

# Gender
dbRDA_0 <- capscale(phyloseq_rarefied_phages_transformed_bray ~ 1, data = Metadata)
dbRDA_1 <- capscale(phyloseq_rarefied_phages_transformed_bray ~  Geslacht, data = Metadata)
step.res_cum_Gender <- ordiR2step(dbRDA_0, dbRDA_1, direction="forward", R2scope = TRUE)
RsquareAdj(step.res_cum_Gender)$r.squared
RsquareAdj(step.res_cum_Gender)$adj.r.squared
step.res_cum_Gender$anova
step.res_cum_Gender <- step.res_cum_Gender$anova[1,1:5]

# Environment1
dbRDA_0 <- capscale(phyloseq_rarefied_phages_transformed_bray ~ 1, data = Metadata)
dbRDA_1 <- capscale(phyloseq_rarefied_phages_transformed_bray ~  Environment1, data = Metadata)
step.res_cum_environment <- ordiR2step(dbRDA_0, dbRDA_1, direction="forward", R2scope = TRUE)
RsquareAdj(step.res_cum_environment)$r.squared
RsquareAdj(step.res_cum_environment)$adj.r.squared
step.res_cum_environment$anova
step.res_cum_environment <- step.res_cum_environment$anova[1,1:5]

# Province
dbRDA_0 <- capscale(phyloseq_rarefied_phages_transformed_bray ~ 1, data = Metadata)
dbRDA_1 <- capscale(phyloseq_rarefied_phages_transformed_bray ~  Province, data = Metadata)
step.res_cum_province <- ordiR2step(dbRDA_0, dbRDA_1, direction="forward", R2scope = TRUE)
RsquareAdj(step.res_cum_province)$r.squared
RsquareAdj(step.res_cum_province)$adj.r.squared
step.res_cum_province$anova
step.res_cum_province <- step.res_cum_province$anova[1,1:5]

# Hurley
dbRDA_0 <- capscale(phyloseq_rarefied_phages_transformed_bray ~ 1, data = Metadata, na.action = na.exclude)
dbRDA_1 <- capscale(phyloseq_rarefied_phages_transformed_bray ~  Hurley, data = Metadata, na.action = na.exclude)
step.res_cum_hurley <- ordiR2step(dbRDA_0, dbRDA_1, direction="forward", add = TRUE, R2scope = FALSE)
RsquareAdj(step.res_cum_hurley)$r.squared
RsquareAdj(step.res_cum_hurley)$adj.r.squared
step.res_cum_hurley$anova
step.res_cum_hurley <- step.res_cum_hurley$anova[1,1:5]

# Age
dbRDA_0 <- capscale(phyloseq_rarefied_phages_transformed_bray ~ 1, data = Metadata)
dbRDA_1 <- capscale(phyloseq_rarefied_phages_transformed_bray ~  Age, data = Metadata)
step.res_cum_Age <- ordiR2step(dbRDA_0, dbRDA_1, direction="forward", R2scope = TRUE)
RsquareAdj(step.res_cum_Age)$r.squared
RsquareAdj(step.res_cum_Age)$adj.r.squared
step.res_cum_Age$anova
step.res_cum_Age <- step.res_cum_Age$anova[1,1:5]

# Obesity: re-run if not significant, mildly signifcant for patient
dbRDA_0 <- capscale(phyloseq_rarefied_phages_transformed_bray ~ 1, data = Metadata)
dbRDA_1 <- capscale(phyloseq_rarefied_phages_transformed_bray ~  Obese, data = Metadata)
step.res_cum_obesity <- ordiR2step(dbRDA_0, dbRDA_1, direction="forward", R2scope = TRUE)
RsquareAdj(step.res_cum_obesity)$r.squared
RsquareAdj(step.res_cum_obesity)$adj.r.squared
step.res_cum_obesity$anova
step.res_cum_obesity <- step.res_cum_obesity$anova[1,1:5]

# Smoking
dbRDA_0 <- capscale(phyloseq_rarefied_phages_transformed_bray ~ 1, data = Metadata)
dbRDA_1 <- capscale(phyloseq_rarefied_phages_transformed_bray ~  Roken_1, data = Metadata)
step.res_cum_smoking <- ordiR2step(dbRDA_0, dbRDA_1, direction="forward", R2scope = TRUE)
RsquareAdj(step.res_cum_smoking)$r.squared
RsquareAdj(step.res_cum_smoking)$adj.r.squared
step.res_cum_smoking$anova
step.res_cum_smoking <- step.res_cum_smoking$anova[1,1:5]

# IBD
dbRDA_0 <- capscale(phyloseq_rarefied_phages_transformed_bray ~ 1, data = Metadata)
dbRDA_1 <- capscale(phyloseq_rarefied_phages_transformed_bray ~  IBD, data = Metadata)
step.res_cum_IBD <- ordiR2step(dbRDA_0, dbRDA_1, direction="forward", R2scope = TRUE)
RsquareAdj(step.res_cum_IBD)$r.squared
RsquareAdj(step.res_cum_IBD)$adj.r.squared
step.res_cum_IBD$anova
step.res_cum_IBD <- step.res_cum_IBD$anova[1,1:5]

# Genetics_family.HS
dbRDA_0 <- capscale(phyloseq_rarefied_phages_transformed_bray ~ 1, data = Metadata)
dbRDA_1 <- capscale(phyloseq_rarefied_phages_transformed_bray ~  Genetics_family.HS, data = Metadata)
step.res_cum_GFH <- ordiR2step(dbRDA_0, dbRDA_1, direction="forward", R2scope = TRUE)
RsquareAdj(step.res_cum_GFH)$r.squared
RsquareAdj(step.res_cum_GFH)$adj.r.squared
step.res_cum_GFH$anova
step.res_cum_GFH <- step.res_cum_GFH$anova[1,1:5]

# Genetics_family.skin.diseases
dbRDA_0 <- capscale(phyloseq_rarefied_phages_transformed_bray ~ 1, data = Metadata)
dbRDA_1 <- capscale(phyloseq_rarefied_phages_transformed_bray ~  Genetics_family.skin.diseases, data = Metadata)
step.res_cum_GFS <- ordiR2step(dbRDA_0, dbRDA_1, direction="forward", R2scope = TRUE)
RsquareAdj(step.res_cum_GFS)$r.squared
RsquareAdj(step.res_cum_GFS)$adj.r.squared
step.res_cum_GFS$anova
step.res_cum_GFS <- step.res_cum_GFS$anova[1,1:5]

# Select only significant ones for complete dataset
# s variables
step.res_cum_Disease
step.res_cum_place
step.res_cum_Age
step.res_cum_smoking
step.res_cum_obesity
step.res_cum_Gender
step.res_cum_ethncity
step.res_cum_IBD
step.res_cum_GFH
step.res_cum_GFS
#step.res_cum_hurley

# ns variables
step.res_cum_time
step.res_cum_province
step.res_cum_environment

# Merge all data
# All
DF_1 <- bind_rows(step.res_cum_Disease, anti_join(step.res_cum_province,step.res_cum_Disease, by ="R2.adj"))
DF_2 <- bind_rows(DF_1,anti_join(step.res_cum_Age,DF_1,by="R2.adj"))
DF_3 <- bind_rows(DF_2,anti_join(step.res_cum_obesity,DF_2,by="R2.adj"))

DF_all <- DF_3
DF_all$p_value <-DF_all$`Pr(>F)`
DF_all$`Pr(>F)` <- NULL
n.tests <- nrow (DF_all) 
pval.adj <- p.adjust(DF_all$p_value, method = 'BH', n = n.tests)
DF_all$pval.adj <- pval.adj
DF_all <- DF_all[rev(order(DF_all$R2.adj)),]
DF_all

# HS
DF_1 <- bind_rows(step.res_cum_hurley, anti_join(step.res_cum_obesity,step.res_cum_hurley, by ="R2.adj"))
DF_all <- DF_1
DF_all$p_value <-DF_all$`Pr(>F)`
DF_all$`Pr(>F)` <- NULL
n.tests <- nrow (DF_all) 
pval.adj <- p.adjust(DF_all$p_value, method = 'BH', n = 2)
DF_all$pval.adj <- pval.adj
DF_all <- DF_all[rev(order(DF_all$R2.adj)),]
DF_all

# All: Obese + Environment + Disease + Province
# Disease: Roken + Hurley

# All samples
DF_1 <- bind_rows(step.res_cum_Disease, anti_join(step.res_cum_place,step.res_cum_Disease, by ="R2.adj"))
DF_2 <- bind_rows(DF_1,anti_join(step.res_cum_Age,DF_1,by="R2.adj"))
DF_3 <- bind_rows(DF_2,anti_join(step.res_cum_obesity,DF_2,by="R2.adj"))
DF_4 <- bind_rows(DF_3,anti_join(step.res_cum_smoking,DF_2,by="R2.adj"))
DF_5 <- bind_rows(DF_4,anti_join(step.res_cum_Gender,DF_2,by="R2.adj"))
DF_6 <- bind_rows(DF_5,anti_join(step.res_cum_ethncity,DF_2,by="R2.adj"))
DF_7 <- bind_rows(DF_6,anti_join(step.res_cum_IBD,DF_2,by="R2.adj"))
DF_8 <- bind_rows(DF_7,anti_join(step.res_cum_GFH,DF_2,by="R2.adj"))
DF_9 <- bind_rows(DF_8,anti_join(step.res_cum_GFS,DF_2,by="R2.adj"))

DF_all <- DF_9
DF_all$p_value <-DF_all$`Pr(>F)`
DF_all$`Pr(>F)` <- NULL
n.tests <- nrow (DF_all) 
pval.adj <- p.adjust(DF_all$p_value, method = 'BH', n = n.tests)
DF_all$pval.adj <- pval.adj
DF_all <- DF_all[rev(order(DF_all$R2.adj)),]
DF_all
####################################
# 6.3 Multivariate analysis
####################################
Metadata <- data.frame(sample_data(phyloseq_rarefied_phages_transformed))
#Metadata$Patient_ID <- NULL
unique(Metadata$Environment)
unique(Metadata$Obese)
unique(Metadata$Disease)
# View(Metadata)

# All
dbRDA_0 <- capscale(phyloseq_rarefied_phages_transformed_bray ~ 1, data = Metadata, na.action = na.exclude)  
dbRDA_1 <- capscale(phyloseq_rarefied_phages_transformed_bray ~ Obese + Disease  + Age, data = Metadata, na.action = na.exclude) # HS
#dbRDA_0 <- capscale(phyloseq_rarefied_phages_transformed_bray ~ 1 + Condition(Obese), data = Metadata, na.action = na.exclude)  
#dbRDA_1 <- capscale(phyloseq_rarefied_phages_transformed_bray ~ Disease + Condition(Obese), data = Metadata, na.action = na.exclude) # HS
step.res_cum_ALL <- ordiR2step(dbRDA_0, scope = formula(dbRDA_1), direction="forward", R2scope = TRUE)

RsquareAdj(step.res_cum_ALL)$r.squared
RsquareAdj(step.res_cum_ALL)$adj.r.squared
(step.res_cum_ALL)V
View(step.res_cum_ALL)

# Variables are Hurley, Obesity and Environment

# Hurley
dbRDA_0 <- capscale(phyloseq_rarefied_phages_transformed_bray ~ 1, data = Metadata)  
dbRDA_1 <- capscale(phyloseq_rarefied_phages_transformed_bray ~ Hurley + Obese, data = Metadata) 
step.res_cum_ALL <- ordiR2step(dbRDA_0, scope = formula(dbRDA_1), direction="forward", R2scope = TRUE)
RsquareAdj(step.res_cum_ALL)$r.squared
RsquareAdj(step.res_cum_ALL)$adj.r.squared
(step.res_cum_ALL)$anova
# Variables are Roken and Hurley
View(Metadata)

# All samples
colnames(Metadata)

dbRDA_0 <- capscale(phyloseq_rarefied_phages_transformed_bray ~ 1, data = Metadata)  
dbRDA_1 <- capscale(phyloseq_rarefied_phages_transformed_bray ~ Place+Obese+Disease+Roken+Genetics_family.HS+Ethnicity+Geslacht+Genetics_family.skin.diseases+Age+IBD, data = Metadata) 
step.res_cum_ALL <- ordiR2step(dbRDA_0, scope = formula(dbRDA_1), direction="forward", R2scope = TRUE)
RsquareAdj(step.res_cum_ALL)$r.squared
RsquareAdj(step.res_cum_ALL)$adj.r.squared
(step.res_cum_ALL)$anova
# Variables are Roken and Hurley
View(Metadata)
####################################
# 6.4 Figure Effect size
####################################
# All individuals
DF_all # univariate
(step.res_cum_ALL)$anova # multivariate

DF_variables <- as.data.frame(c("Obese","Disease","Province","Age","Obese","Disease","Province","Age"))
DF_variables$variables <- DF_variables$`c("Obese", "Disease", "Province", "Age", "Obese", "Disease", "Province", "Age")`
DF_variables$`c("Obese", "Disease", "Province", "Age", "Obese", "Disease", "Province", "Age")` <- NULL  
DF_variables$variable_type <- c("multivariate","multivariate","multivariate","multivariate","univariate","univariate","univariate","univariate")
DF_variables$values <- c("5.16","6.79","6.79","6.79","5.16","3.32","2.69","1.56")
ord <- c("Obese","Disease","Province","Age")
DF_variables$variables <- factor(DF_variables$variables, levels = rev(ord))
DF_variables$values <- as.numeric(DF_variables$values)

ggplot(DF_variables) +
  aes(x = variables, fill = variable_type, weight = values) +
  geom_bar(position = "dodge") +
  scale_fill_manual(values = c(multivariate = "#000000", univariate = "#8B8A8B")) +
  coord_flip() +
  theme_bw() +
  labs(x= "covariates of virome composition", y = "effect size (%)") +
  theme(legend.title = element_blank()) +
  geom_hline(yintercept = 6.79, linetype = 5, colour = "#49006a")

# All individuals; explaines versus non-explained
DF_variables <- as.data.frame(c("Obese","Disease","Unexplained"))
DF_variables$variables <- DF_variables$`c("Obese", "Disease", "Unexplained")`
DF_variables$`c("Obese", "Disease", "Unexplained")` <- NULL
DF_variables$variable_type <- c("multivariate","multivariate","multivariate")
DF_variables$values <- c("5.16","1.63","93.21")
ord <- c("Unexplained","Obese","Disease")
DF_variables$variables <- factor(DF_variables$variables, levels = rev(ord))
DF_variables$values <- as.numeric(DF_variables$values)

ggplot(DF_variables) +
 aes(x = variable_type, fill = variables, weight = values) +
 geom_bar() +
 scale_fill_manual(values = c(Disease = "#44341d", Obese = "#fdc32e", Unexplained = "#aec3e6")) +
 theme_bw()

# HS
DF_all # univariate
(step.res_cum_ALL)$anova # multivariate

DF_variables <- as.data.frame(c("Hurley","Obese","Hurley","Obese"))
DF_variables$variables <- DF_variables$`c("Hurley", "Obese", "Hurley", "Obese")`
DF_variables$`c("Hurley", "Obese", "Hurley", "Obese")` <- NULL  
DF_variables$variable_type <- c("multivariate","multivariate","univariate","univariate")
DF_variables$values <- c("3.44","6.25","3.44","2.56")
ord <- c("Hurley","Obese")
DF_variables$variables <- factor(DF_variables$variables, levels = rev(ord))
DF_variables$values <- as.numeric(DF_variables$values)

ggplot(DF_variables) +
  aes(x = variables, fill = variable_type, weight = values) +
  geom_bar(position = "dodge") +
  scale_fill_manual(values = c(multivariate = "#000000", univariate = "#8B8A8B")) +
  coord_flip() +
  theme_bw() +
  labs(x= "covariates of virome composition", y = "effect size (%)") +
  theme(legend.title = element_blank()) +
  geom_hline(yintercept = 6.25, linetype = 5, colour = "#49006a")
####################################
# 6.5 Association variables: Disease and Obesity
####################################
library(chisq.posthoc.test)

# 1) Disease versus Obesity
Metadata <- data.frame(sample_data(phyloseq_rarefied_phages_transformed))
Metadata <- Metadata[,c(2,9)]
Metadata_P_HS <- table(Metadata$Disease, Metadata$Obese)
Metadata_GLM_table_prop <-  as.data.frame(prop.table(Metadata_P_HS, margin = 1))
Metadata_GLM_table_prop$Freq <- Metadata_GLM_table_prop$Freq*100
Metadata_GLM_table_prop$Var2 <- factor(Metadata_GLM_table_prop$Var2, levels = c("Underweight","Normal","Overweight","Obese"))
Metadata_GLM_table_prop$Var1 <- factor(Metadata_GLM_table_prop$Var1, levels = c("HE","HS"))
Metadata_P_HS_DF <- as.data.frame(Metadata_P_HS)
Metadata_P_HS_DF$Var2 <- factor(Metadata_P_HS_DF$Var2, levels = c("Underweight","Normal","Overweight","Obese"))
Metadata_P_HS_DF$Var1 <- factor(Metadata_P_HS_DF$Var1, levels = c("HE","HS"))

ggplot(Metadata_P_HS_DF) +
  aes(x = Var2, fill = Var1, weight = Freq) +
  geom_bar(position = "dodge") +
  ylab("Patient (#)") +
  scale_fill_manual(values = list("HS"="#bd801c","HE" = "#2b8a92")) +
  theme_bw()

ggplot(Metadata_GLM_table_prop) +
  aes(x = Var1, fill = Var2, weight = Freq) +
  geom_bar(position = "dodge") +
  ylab("Prevalence (%)") +
  xlab("Disease status") +
  scale_fill_manual(values = list("Underweight"="#e7298a","Normal" = "#1b9e77", "Overweight" = "#7570b3", "Obese" = "#d95f02")) +
  theme_classic()

## Make proportion table like in figure underneath
Metadata_GLM_table_prop2 <-  as.data.frame(prop.table(Metadata_P_HS, margin = 2))
Metadata_GLM_table_prop2$Freq <- Metadata_GLM_table_prop2$Freq*100
Metadata_GLM_table_prop2$Var1 <- factor(Metadata_GLM_table_prop2$Var1, levels = c("HE","HS"))
Metadata_GLM_table_prop2$Var2 <- factor(Metadata_GLM_table_prop2$Var2, levels = c("Underweight","Normal","Overweight","Obese"))

ggplot(Metadata_GLM_table_prop2) +
  aes(x = Var2, fill = Var1, weight = Freq) +
  geom_bar() +
  ylab("Prevalence (%)") +
  xlab("") +
  scale_fill_manual(values = list("HS"="#bd801c","HE" = "#2b8a92")) +
  theme_classic()

ggplot(Metadata_GLM_table_prop) +
  aes(x = Var1, fill = Var2, weight = Freq) +
  geom_bar() +
  ylab("Prevalence (%)") +
  xlab("Disease status") +
  scale_fill_manual(values = list("Underweight"="#e7298a","Normal" = "#1b9e77", "Overweight" = "#7570b3", "Obese" = "#d95f02")) +
  theme_classic()

chisq.test(Metadata$Disease, Metadata$Obese)
cramer_v(Metadata$Disease, Metadata$Obese)

# If significant, do pairwise proportion test
row_wise_prop_test(t(Metadata_P_HS), detailed=T) ## Normal and obese people differ between groups

# 2) Hurley versus Obesity
Metadata <- data.frame(sample_data(phyloseq_rarefied_phages_transformed))
Metadata_P <- Metadata
Metadata_P <- Metadata[!Metadata$Hurley=="0",]
Metadata_P <- Metadata_P[,c(7,9)]
Metadata_P$Obese <- factor(Metadata_P$Obese, levels = c("Underweight","Normal","Overweight","Obese"))
#View(Metadata_P)

ggplot(Metadata_P) +
 aes(x = Hurley, fill = Obese) +
 geom_bar(position = "dodge") +
  ylab("Patient (#)") +
  xlab("Disease status") +
  scale_fill_manual(values = list("Underweight"="#e7298a","Normal" = "#1b9e77", "Overweight" = "#7570b3", "Obese" = "#d95f02")) +
 theme_bw()

ggplot(Metadata_P) +
  aes(x = Obese, fill = Hurley) +
  geom_bar(position = "dodge") +
  ylab("Patient (#)") +
  xlab("") +
  scale_fill_manual(values = list("2" = "#1b9e77", "1" = "#7570b3", "3" = "#d95f02")) +
  theme_bw()

## Make proportion table like in figure underneath
Metadata_P_HS <- table(Metadata_P$Hurley, Metadata_P$Obese)
Metadata_GLM_table_prop2 <-  as.data.frame(prop.table(Metadata_P_HS, margin = 2))
Metadata_GLM_table_prop2$Freq <- Metadata_GLM_table_prop2$Freq*100
Metadata_GLM_table_prop2$Var1 <- factor(Metadata_GLM_table_prop2$Var1, levels = c("0","1", "2", "3"))
Metadata_GLM_table_prop2$Var2 <- factor(Metadata_GLM_table_prop2$Var2, levels = c("Underweight","Normal","Overweight","Obese"))

ggplot(Metadata_GLM_table_prop2) +
  aes(x = Var2, fill = Var1, weight = Freq) +
  geom_bar() +
  ylab("Prevalence (%)") +
  xlab("") +
  scale_fill_manual(values = list("0"="#2b8a92","2" = "#1b9e77", "1" = "#7570b3", "3" = "#d95f02")) +
  theme_classic()
  
# Pairwise proportion test
Metadata_GLM_table_prop <-  as.data.frame(prop.table(Metadata_P_HS, margin = 1))
Metadata_GLM_table_prop$Freq <- Metadata_GLM_table_prop$Freq*100
Metadata_GLM_table_prop$Var2 <- factor(Metadata_GLM_table_prop$Var2, levels = c("Normal","Overweight","Obese"))

# Barplot next to eachother
ggplot(Metadata_GLM_table_prop) +
  aes(x = Var1, fill = Var2, weight = Freq) +
  geom_bar(position = "dodge", alpha=0.9) +
  ylab("Prevalence (%)") +
  scale_fill_manual(values = list("Normal" = "#1b9e77", "Overweight" = "#7570b3", "Obese" = "#d95f02")) +
  theme_classic()

ggplot(Metadata_GLM_table_prop) +
  aes(x = Var1, fill = Var2, weight = Freq) +
  geom_bar() +
  ylab("Prevalence (%)") +
  scale_fill_manual(values = list("Normal" = "#1b9e77", "Overweight" = "#7570b3", "Obese" = "#d95f02")) +
  theme_classic()

chisq.test(Metadata_P$Hurley, Metadata_P$Obese)
?chisq.test
row_wise_prop_test(t(Metadata_P_HS)) ## Normal and obese people differ between groups
?row_wise_prop_test
cramer_v(Metadata_P$Hurley, Metadata_P$Obese)
####################################
# 6.6 RA barplot
####################################
# 6.6.1 RA barplot: Disease
####################################
# Disease status
names <- c("family1","family2","family3","family4","family5","family6","family10","family11","family13","family17","family30","family34","family60")
#View(otu_table(phyloseq_rarefied_phages))

phyloseq_rarefied_phages_20 <- phyloseq_rarefied_phages %>% 
  microbiome::transform(transform = "compositional")

#phyloseq_rarefied_phages_20_ <- subset_samples(phyloseq_rarefied_phages_20,Obese=="Normal")
phyloseq_rarefied_phages_20_ <- phyloseq_rarefied_phages_20
phyloseq_rarefied_phages_20_DF <- as.data.frame(otu_table(phyloseq_rarefied_phages_20_))
phyloseq_rarefied_phages_20_DF$x1 <- rownames(phyloseq_rarefied_phages_20_DF)
phyloseq_rarefied_phages_20_DF_disease <- filter(phyloseq_rarefied_phages_20_DF, x1 %in% names)   

# Merge metadata
Heatmap_family_RA_m <- merge(t(phyloseq_rarefied_phages_20_DF_disease), metadata_patient, by = 0, all = F)
rownames(Heatmap_family_RA_m) <- Heatmap_family_RA_m$Row.names
Heatmap_family_RA_m$Row.names <- NULL

# Select columns
colnames(Heatmap_family_RA_m)
vector_1 <- (which(names(Heatmap_family_RA_m)== "Patient_ID"))
vector_2 <- (which(names(Heatmap_family_RA_m)== "Disease")) 
Heatmap_family_RA_met <- Heatmap_family_RA_m[,c(1:vector_1,vector_2)]

# Melt
Heatmap_family_RA_melt <- melt(Heatmap_family_RA_met, id.vars = c("Patient_ID","Disease"))
Heatmap_family_RA_melt$value <- as.numeric(Heatmap_family_RA_melt$value)
#Heatmap_family_RA_melt$variable <- factor(Heatmap_family_RA_melt$variable, levels=rev(c("family1","family2","family3","family4","family5","family6","family10","family11","family13","family17","family30","family34","family60")))
Heatmap_family_RA_melt$variable <- factor(Heatmap_family_RA_melt$variable, levels=rev(c("family60","family34","family30","family17","family13","family11","family10","family6","family5","family4","family3","family2","family1")))
#Heatmap_family_RA_melt$Disease <- factor(Heatmap_family_RA_melt$Disease, levels=rev(c("HE","HS")))
Heatmap_family_RA_melt$Hurley <- factor(Heatmap_family_RA_melt$Disease, levels=rev(c("HS","HE")))

#plot_ef_bar(mm_lefse_table, label_level = 0) +
#  scale_fill_manual(values = c(HE = "#2d8b92", HS = "#bd801c")) 
Heatmap_family_RA_melt

ggplot(Heatmap_family_RA_melt) +
  aes(x = variable, y = value, fill = Disease) +
  geom_boxplot() +
  scale_fill_manual(values = c(HE = "#2d8b92", HS = "#bd801c")) +
  theme_bw() +
  xlab("Relative abundance")+
  ylim(0,0.5)

ggplot(Heatmap_family_RA_melt) +
 aes(y = variable, x = value, fill = Disease) +
 geom_boxplot() +
 scale_fill_manual(values = c(HE = "#2d8b92", HS = "#bd801c")) +
 theme_bw() +
  xlab("Relative abundance")+
 xlim(0,0.5)

Statistics_f <- Heatmap_family_RA_melt[c(Heatmap_family_RA_melt$variable == "family1"),]
wilcox.test(value~Disease, data =Statistics_f,alternative = "two.sided",  paired = FALSE, exact = F)
effect_size <- wilcox_effsize(Statistics_f,value~Disease)
(effect_size$effsize)^2
effect_size$n1 
effect_size$n2 

Statistics_f <- Heatmap_family_RA_melt[c(Heatmap_family_RA_melt$variable == "family2"),]
wilcox.test(value~Disease, data =Statistics_f,alternative = "two.sided",  paired = FALSE, exact = F)
effect_size <- wilcox_effsize(Statistics_f,value~Disease)
(effect_size$effsize)^2
effect_size$n1 
effect_size$n2 

Statistics_f <- Heatmap_family_RA_melt[c(Heatmap_family_RA_melt$variable == "family4"),]
wilcox.test(value~Disease, data =Statistics_f,alternative = "two.sided",  paired = FALSE, exact = F)
effect_size <- wilcox_effsize(Statistics_f,value~Disease)
(effect_size$effsize)^2
effect_size$n1 
effect_size$n2 

Statistics_f <- Heatmap_family_RA_melt[c(Heatmap_family_RA_melt$variable == "family4"),]
wilcox.test(value~Disease, data =Statistics_f,alternative = "two.sided",  paired = FALSE, exact = F)
effect_size <- wilcox_effsize(Statistics_f,value~Disease)
(effect_size$effsize)^2
effect_size$n1 
effect_size$n2 

Statistics_f <- Heatmap_family_RA_melt[c(Heatmap_family_RA_melt$variable == "family5"),]
wilcox.test(value~Disease, data =Statistics_f,alternative = "two.sided",  paired = FALSE, exact = F)
effect_size <- wilcox_effsize(Statistics_f,value~Disease)
(effect_size$effsize)^2
effect_size$n1 
effect_size$n2 
 
Statistics_f <- Heatmap_family_RA_melt[c(Heatmap_family_RA_melt$variable == "family6"),]
wilcox.test(value~Disease, data =Statistics_f,alternative = "two.sided",  paired = FALSE, exact = F)
effect_size <- wilcox_effsize(Statistics_f,value~Disease)
(effect_size$effsize)^2
effect_size$n1 +effect_size$n2 

# Multiple testing correction 
P_values <-  c(0.00353,0.660,0.00321,0.0659,0.992,0.0383,0.0259,0.0814,0.135,0.965,0.104,0.154,0.534)
p.adjust(P_values, method = "BH")

# Association Disease status and location
colnames(Mastertable_phage)

# Find all samples with phages inside
names <- colnames(Mastertable_viral1[colSums(Mastertable_viral1[,1:127]) > 0])
names <- names[1:109]
names2
metadata$names <- rownames(metadata)
metadata_samples <- metadata[is.element(metadata$names, names2), ]
metadata_samples$Place[!metadata_samples$Place=="Axilliary" & !metadata_samples$Place == "Groin and genitals"] <- "Abd/Glu/Chest"
metadata_samples <- metadata_samples[!metadata_samples$Place=="Abd/Glu/Chest",]
table(metadata_samples$Place)

# Preliminary
metadata_samples_HE <- metadata_samples[metadata_samples$Disease=="HE",]
table(metadata_samples_HE$Place)

metadata_samples_HS <- metadata_samples[metadata_samples$Disease=="HS",]
table(metadata_samples_HS$Place)

# Association
nrow(metadata_samples)
Metadata_P_HS <- table(metadata_samples$Disease, metadata_samples$Place)
Metadata_GLM_table_prop <-  as.data.frame(prop.table(Metadata_P_HS, margin = 1))
Metadata_GLM_table_prop$Freq <- Metadata_GLM_table_prop$Freq*100
Metadata_GLM_table_prop$Var2 <- factor(Metadata_GLM_table_prop$Var2, levels = c("Axilliary","Groin and genitals"))
Metadata_GLM_table_prop$Var1 <- factor(Metadata_GLM_table_prop$Var1, levels = c("HE","HS"))
Metadata_P_HS_DF <- as.data.frame(Metadata_P_HS)
Metadata_P_HS_DF$Var2 <- factor(Metadata_P_HS_DF$Var2, levels = c("Axilliary","Groin and genitals"))
Metadata_P_HS_DF$Var1 <- factor(Metadata_P_HS_DF$Var1, levels = c("HE","HS"))
Metadata_GLM_table_prop
Metadata_P_HS_DF
Metadata_P_HS

# Figure
ggplot(Metadata_GLM_table_prop) +
 aes(x = Var2, fill = Var1, weight = Freq) +
 geom_bar(position = "dodge") +
  ylab("Prevalence (%)") +
  xlab("Disease location") +
  scale_fill_manual(values = list("HS"="#bd801c","HE" = "#2b8a92")) +
 theme_bw()

chisq.test(Metadata_P_HS)
?chisq.test
cramer_v(Metadata_P_HS)
pairwise_prop_test(Metadata_P_HS, alternative = c("two.sided"), p.adjust.method = "BH")
row_wise_prop_test(t(Metadata_P_HS), alternative = c("two.sided"), p.adjust.method = "BH")
row_wise_prop_test(t(Metadata_P_HS), detailed=T, alternative = c("two.sided"), p.adjust.method = "BH")
?row_wise_prop_test
####################################
# 6.6.2 RA barplot: Hurley
####################################
# Hurley status
names <- c("family1","family2","family3","family4","family5","family6","family10","family11","family13","family17","family30", "family34","family60")

phyloseq_rarefied_phages_20 <- phyloseq_rarefied_phages %>% 
  microbiome::transform(transform = "compositional") 

View(Mastertable[c(Mastertable$Family_level=="family83" & Mastertable$Phage_score1=='phage'),])
View(Mastertable[Mastertable$names=="NODE_A1_length_40163_cov_80_378162_HS2_T2_L2",])
#phyloseq_rarefied_phages_20_ <- subset_samples(phyloseq_rarefied_phages_20,Place=="Groin and genitals" && Place=="Axilliary")
phyloseq_rarefied_phages_20_ <- phyloseq_rarefied_phages_20

#View(sample_data(phyloseq_rarefied_phages_20_))
phyloseq_rarefied_phages_20_DF <- as.data.frame(otu_table(phyloseq_rarefied_phages_20_))
phyloseq_rarefied_phages_20_DF$x1 <- rownames(phyloseq_rarefied_phages_20_DF)
phyloseq_rarefied_phages_20_DF_disease <- filter(phyloseq_rarefied_phages_20_DF, x1 %in% names)   
#View(sample_data(phyloseq_rarefied_phages_20_DF_disease))

# Merge metadata
Heatmap_family_RA_m <- merge(t(phyloseq_rarefied_phages_20_DF_disease), metadata_patient, by = 0, all = F)
rownames(Heatmap_family_RA_m) <- Heatmap_family_RA_m$Row.names
Heatmap_family_RA_m$Row.names <- NULL
#View(Heatmap_family_RA_m)

# Select columns
colnames(Heatmap_family_RA_m)
vector_1 <- (which(names(Heatmap_family_RA_m)== "Patient_ID"))
vector_2 <- (which(names(Heatmap_family_RA_m)== "Hurley")) 
Heatmap_family_RA_met <- Heatmap_family_RA_m[,c(1:vector_1,vector_2)]

# Melt
#View(Heatmap_family_RA_met)
Heatmap_family_RA_melt <- melt(Heatmap_family_RA_met, id.vars = c("Patient_ID","Hurley"))
Heatmap_family_RA_melt$value <- as.numeric(Heatmap_family_RA_melt$value)
#Heatmap_family_RA_melt$variable <- factor(Heatmap_family_RA_melt$variable, levels=rev(c("family1","family2","family3","family4","family5","family6","family10","family11","family13","family17","family30", "family34","family60")))
Heatmap_family_RA_melt$variable <- factor(Heatmap_family_RA_melt$variable, levels=rev(c("family60","family34","family30","family17","family13","family11","family10","family6","family5","family4","family3","family2","family1")))
#Heatmap_family_RA_melt$Hurley <- factor(Heatmap_family_RA_melt$Hurley, levels=rev(c("0","1","2","3")))
Heatmap_family_RA_melt$Hurley <- factor(Heatmap_family_RA_melt$Hurley, levels=rev(c("3","2","1","0")))

#View(Heatmap_family_RA_melt)
#plot_ef_bar(mm_lefse_table_HS, label_level = 0) +
#  scale_fill_manual(values = c("0" = "#499B78", "1"="blue", "2"="orange","3"="yellow"))

#equisser(Heatmap_family_RA_melt)

ggplot(Heatmap_family_RA_melt) +
 aes(x = variable, y = value, fill = Hurley) +
 geom_boxplot() +
  scale_fill_manual(values = c("0" = "#499B78", "1"="blue", "2"="orange","3"="yellow")) +
  theme_bw() +
  ylim(0.0,0.5)

ggplot(Heatmap_family_RA_melt) +
  aes(x = variable, y = value, fill = Hurley) +
  geom_boxplot() +
  scale_fill_manual(values = c("0" = "#499B78", "1"="blue", "2"="orange","3"="yellow")) +
  theme_bw() +
  xlim(0.0,0.5)

# Hurley 3
Statistics_family <- Heatmap_family_RA_melt[c(Heatmap_family_RA_melt$variable == "family60"),]
kruskal.test(value ~ Hurley, data = Statistics_family)
effect_size <- kruskal_effsize(Statistics_family,value~Hurley)
effect_size$effsize
effect_size$n
dunnTest(value ~ Hurley, data = Statistics_family, method="bh")

Statistics_family <- Heatmap_family_RA_melt[c(Heatmap_family_RA_melt$variable == "family2"),]
kruskal.test(value ~ Hurley, data = Statistics_family)
dunnTest(value ~ Hurley, data = Statistics_family, method="bh")

# Hurley 2
Statistics_family <- Heatmap_family_RA_melt[c(Heatmap_family_RA_melt$variable == "family3"),]
kruskal.test(value ~ Hurley, data = Statistics_family)
kruskal_effsize(Statistics_f,value~Disease)
effect_size <- wilcox_effsize(Statistics_f,value~Disease)
dunnTest(value ~ Hurley, data = Statistics_family, method="bh")

Statistics_family <- Heatmap_family_RA_melt[c(Heatmap_family_RA_melt$variable == "family4"),]
kruskal.test(value ~ Hurley, data = Statistics_family)
dunnTest(value ~ Hurley, data = Statistics_family, method="bh")

# Hurley 1
Statistics_family <- Heatmap_family_RA_melt[c(Heatmap_family_RA_melt$variable == "family5"),]
kruskal.test(value ~ Hurley, data = Statistics_family)
dunnTest(value ~ Hurley, data = Statistics_family, method="bh")

# Hurley 0
Statistics_family <- Heatmap_family_RA_melt[c(Heatmap_family_RA_melt$variable == "family6"),]
kruskal.test(value ~ Hurley, data = Statistics_family)
dunnTest(value ~ Hurley, data = Statistics_family, method="bh")

Statistics_family <- Heatmap_family_RA_melt[c(Heatmap_family_RA_melt$variable == "family60"),]
kruskal.test(value ~ Hurley, data = Statistics_family)
dunnTest(value ~ Hurley, data = Statistics_family, method="bh")

# Multiple testing correction 
P_values <-  c(0.235,
               0.659,
               0.143,
               0.00132,
               0.0242,
               0.0108,
               0.767,
               0.128,
               0.228,
               0.393,
               0.547,
               0.997,
               0.452)

p.adjust(P_values, method = "BH")

# Association Disease seveirty versus location
phyloseq_rarefied_phages_20 <- phyloseq_rarefied_family
sample_data(phyloseq_rarefied_phages_20)$Time <- as.character(sample_data(phyloseq_rarefied_phages_20)$Time)
phyloseq_rarefied_phages_20_ <- subset_samples(phyloseq_rarefied_phages_20, Time=="0"| c(HS=="HS6" & Time=="5")| c(HS=="HS3" & Time=="8")| c(HS=="HS10" & Time=="5")| c(HS=="HS15" & Time=="3")| c(HS=="HS18" & Time=="8") | c(HS=="HS22" & Time=="6")| c(HS=="HS24" & Time=="9") | c(HS=="HS34" & Time=="10"))
phyloseq_rarefied_phages_20_1 <- subset_samples(phyloseq_rarefied_phages_20_, !Patient_ID=="HS11_T0_L2")
phyloseq_rarefied_phages_20_2 <- subset_samples(phyloseq_rarefied_phages_20_1, Place =="Axilliary" | Place =="Groin and genitals")
metadata_samples <- sample_data(phyloseq_rarefied_phages_20_2)
table(metadata_samples$Place)

#names <- colnames(Mastertable_viral1[colSums(Mastertable_viral1[,1:127]) > 0])
#names <- names[1:109]
#names
#metadata$names <- rownames(metadata)
#metadata_samples <- metadata[is.element(metadata$names, names2), ]
#metadata_samples <- metadata_samples[metadata_samples$Place=="Axilliary" | metadata_samples$Place == "Groin and genitals",]
#metadata_samples <- metadata_samples[metadata_samples$Time=="0",]
table(metadata_samples$Place)

# Preliminary
metadata_samples_HE <- metadata_samples[metadata_samples$Disease=="HE",]
table(metadata_samples_HE$Place)

metadata_samples_HS <- metadata_samples[metadata_samples$Disease=="HS",]
table(metadata_samples_HS$Place)

# Association
nrow(metadata_samples)
Metadata_P_HS <- table(metadata_samples$Hurley, metadata_samples$Place)
Metadata_GLM_table_prop <-  as.data.frame(prop.table(Metadata_P_HS, margin = 2))
Metadata_GLM_table_prop$Freq <- Metadata_GLM_table_prop$Freq*100
Metadata_GLM_table_prop$Var2 <- factor(Metadata_GLM_table_prop$Var2, levels = c("Axilliary","Groin and genitals","Gluteaal","Chest","Abdominal"))

Metadata_GLM_table_prop$Var1
Metadata_GLM_table_prop$Var1 <- factor(Metadata_GLM_table_prop$Var1, levels = c("0","1","2","3"))
Metadata_P_HS_DF <- as.data.frame(Metadata_P_HS)
Metadata_P_HS_DF$Var2 <- factor(Metadata_P_HS_DF$Var2, levels = c("Axilliary","Groin and genitals","Gluteaal","Chest","Abdominal"))
Metadata_P_HS_DF$Var1 <- factor(Metadata_P_HS_DF$Var1, levels = c("0","1","2","3"))
Metadata_GLM_table_prop
Metadata_P_HS_DF
Metadata_P_HS
sum(Metadata_P_HS_DF$Freq)

# Figure
Metadata_GLM_table_prop

ggplot(Metadata_GLM_table_prop) +
  aes(x = Var2, fill = Var1, weight = Freq) +
  geom_bar(position = "dodge") +
  ylab("Prevalence (%)") +
  xlab("Disease location") +
  scale_fill_manual(values = list("2" = "#1b9e77", "1" = "#7570b3", "3" = "#d95f02", "0"="2b8a92")) +
  theme_bw()

ggplot(Metadata_GLM_table_prop) +
 aes(x = Var1, fill = Var2, weight = Freq) +
 geom_bar(position = "dodge") +
 scale_fill_hue(direction = 1) +
 theme_bw()

chisq.test(Metadata_P_HS) # But
cramer_v(Metadata_P_HS)
Metadata_P_HS

# Axillary
pairwise_prop_test(Metadata_P_HS, alternative = c("two.sided"), p.adjust.method = "BH")
row_wise_prop_test((Metadata_P_HS), alternative = c("two.sided"), p.adjust.method = "BH")
row_wise_prop_test(t(Metadata_P_HS), detailed=T, alternative = c("two.sided"), p.adjust.method = "BH")
pairwise_prop_test(Metadata_P_HS)
nrow(Metadata_P_HS)
# --> Gr-Axil/Hurley = 55  X-squared = 0.14836, df = 3, p-value = 0.985, cramper-phi= 0.05193726
#
####################################
# 6.6.3 RA barplot: Obesity
####################################
# Obesity 
names <- c("family1","family2","family3","family4","family5","family6","family10","family11","family13","family17","family30", "family34","family60")
phyloseq_rarefied_phages_20_DF <- as.data.frame(otu_table(phyloseq_rarefied_phages_20))
phyloseq_rarefied_phages_20_DF$x1 <- rownames(phyloseq_rarefied_phages_20_DF)
phyloseq_rarefied_phages_20_DF_disease <- filter(phyloseq_rarefied_phages_20_DF, x1 %in% names)   

# Merge metadata
Heatmap_family_RA_m <- merge(t(phyloseq_rarefied_phages_20_DF_disease), metadata_patient, by = 0, all = F)
rownames(Heatmap_family_RA_m) <- Heatmap_family_RA_m$Row.names
Heatmap_family_RA_m$Row.names <- NULL

# Select columns
colnames(Heatmap_family_RA_m)
vector_1 <- (which(names(Heatmap_family_RA_m)== "Patient_ID"))
vector_2 <- (which(names(Heatmap_family_RA_m)== "Obese")) 
Heatmap_family_RA_met <- Heatmap_family_RA_m[,c(1:vector_1,vector_2)]

# Melt
Heatmap_family_RA_melt <- melt(Heatmap_family_RA_met, id.vars = c("Patient_ID","Obese"))
Heatmap_family_RA_melt$value <- as.numeric(Heatmap_family_RA_melt$value)
Heatmap_family_RA_melt$variable <- factor(Heatmap_family_RA_melt$variable, levels=rev(c("family4","family2","family7")))

plot_ef_bar(mm_lefse_table_Obese, label_level = 0) +
  scale_fill_manual(values = c("Obese" = "orange", "Overweight"="brown")) +

  ggplot(Heatmap_family_RA_melt) +
  aes(x = value, y = variable, fill = Obese) +
  geom_boxplot() +
  theme_bw() 

# Statistics
Statistics_family <- Heatmap_family_RA_melt[c(Heatmap_family_RA_melt$variable == "family4"),]
kruskal.test(value ~ Obese, data = Statistics_family)
dunnTest(value ~ Obese, data = Statistics_family,method="bh")

Statistics_family <- Heatmap_family_RA_melt[c(Heatmap_family_RA_melt$variable == "family2"),]
kruskal.test(value ~ Obese, data = Statistics_family)
dunnTest(value ~ Obese, data = Statistics_family,method="bh")

Statistics_family <- Heatmap_family_RA_melt[c(Heatmap_family_RA_melt$variable == "family7"),]
kruskal.test(value ~ Obese, data = Statistics_family)
dunnTest(value ~ Obese, data = Statistics_family,method="bh")
####################################
# 6.7 Phage diversity changes between disease status
####################################
library(microbiome)
library(knitr)

#ps1.fam.rel <- phyloseq_rarefied_phages
ps1.fam.rel <- subset_samples(phyloseq_rarefied_phages, Obese=="Normal")
tab <- microbiome::alpha(ps1.fam.rel, index = c("observed","Shannon", "Simpson"))
#tab <- microbiome::alpha(phyloseq_rarefied_phages, index = c("observed","Shannon", "Simpson"))
kable(head(tab))
View(tab)

# Metadata
metadata_patient_2 <- merge(tab,metadata_patient,by=0,all=F)
metadata_patient_2$Row.names <- rownames(metadata_patient_2)
metadata_patient_2$Row.names <- NULL
colnames(metadata_patient_2) 
metadata_patient_2 <- metadata_patient_2[,c(1:4,6,12,14)]

wilcox.test(diversity_shannon~Disease, data =metadata_patient_2,alternative = "two.sided",  paired = FALSE, exact = F)
effect_size <- wilcox_effsize(diversity_shannon~Disease, data =metadata_patient_2)
(effect_size$effsize)^2
effect_size$n1 
effect_size$n2 

  # Diversity
ggplot(metadata_patient_2) +
 aes(x = "", y = diversity_shannon, fill = Disease) +
 geom_boxplot() +
  xlab("") +
 scale_fill_hue(direction = 1) +
 theme_bw()

ggplot(metadata_patient_2) +
  aes(x = "", y = diversity_shannon, fill = Obese) +
  geom_violin() +
  xlab("") +
  scale_fill_hue(direction = 1) +
  theme_bw()

ggplot(metadata_patient_2) +
  aes(x = "", y = diversity_shannon, fill = Hurley) +
  geom_violin() +
  xlab("") +
  scale_fill_hue(direction = 1) +
  theme_bw()

# Simpson
ggplot(metadata_patient_2) +
  aes(x = "", y = evenness_simpson, fill = Disease) +
  geom_boxplot() +
  xlab("") +
  scale_fill_hue(direction = 1) +
  theme_bw()

ggplot(metadata_patient_2) +
  aes(x = "", y = evenness_simpson, fill = Obese) +
  geom_boxplot() +
  xlab("") +
  scale_fill_hue(direction = 1) +
  theme_bw()

ggplot(metadata_patient_2) +
  aes(x = "", y = evenness_simpson, fill = Hurley) +
  geom_boxplot() +
  xlab("") +
  scale_fill_hue(direction = 1) +
  theme_bw()

  # Observed richness
ggplot(metadata_patient_2) +
 aes(x = "", y = observed, fill = Disease) +
  geom_violin() +
  xlab("") +
 scale_fill_hue(direction = 1) +
  theme_bw()

ggplot(metadata_patient_2) +
  aes(x = "", y = observed, fill = Obese) +
  geom_boxplot() +
  xlab("") +
  scale_fill_hue(direction = 1) +
  theme_bw()

ggplot(metadata_patient_2) +
  aes(x = "", y = observed, fill = Hurley) +
  geom_boxplot() +
  xlab("") +
  scale_fill_hue(direction = 1) +
  theme_bw()

# Compare diversity between Hurley stages
# extract Shannon diversity and Disease variable
shannon = diversity(phyloseq_rarefied_phages, index = "shannon")
disease = sample_data(phyloseq_rarefied_phages)$Hurley

# create boxplot
ggplot(data = data.frame(shannon, disease), aes(x = disease, y = shannon)) +
  geom_boxplot() +
  xlab("Disease") +
  ylab("Shannon Diversity")

kruskal.test(shannon ~ disease, data = data.frame(shannon, disease))

# HS vs HE
shannon = diversity(ps1.fam.rel, index = "shannon")
disease = sample_data(ps1.fam.rel)$Disease

ggplot(data = data.frame(shannon, disease), 
       aes(x = disease, y = shannon), fill = disease) +
  geom_boxplot() +
  xlab("Disease") +
  ylab("Shannon Diversity") + 
  xlab("") +
  scale_fill_hue(direction = 1) +
  theme_bw()

wilcox_test(shannon ~ disease, data = data.frame(shannon, disease))

####################################
# 6.8 Beta diversity
####################################
# 6.8.1 Disease status
####################################
library("phyloseq")
library("ggplot2")
library("dplyr")
library("ggpubr")
library("Matrix")
library("reshape2")
library("vegan")

#relab_genera = transform_sample_counts(ps1.fam.rel, function(x) x / sum(x) * 100) 
relab_genera = transform_sample_counts(ps1.fam.rel, function(x) x / sum(x) * 100) 
head(otu_table(relab_genera)[,1:6])

abrel_bray <- phyloseq::distance(relab_genera, method = "bray")
abrel_bray <- as.matrix(abrel_bray)
head(abrel_bray)[,1:6]

sub_dist <- list()
groups_all <- sample_data(relab_genera)$Disease
groups_all

unique(groups_all)

for (group in unique(groups_all)) { 
  row_group <- which(groups_all == group)
  sample_group <- sample_names(relab_genera)[row_group]
  sub_dist[[group]] <- abrel_bray[ sample_group, sample_group]
  sub_dist[[group]][!lower.tri(sub_dist[[group]])] <- NA
}

braygroups<- melt(sub_dist)
df.bray <- braygroups[complete.cases(braygroups), ]
df.bray$L1 <- factor(df.bray$L1, levels=names(sub_dist))
head(df.bray)

ggplot(df.bray, aes(x=L1, y=value, colour=L1)) +
  geom_jitter() + 
  geom_violin() +  
  theme(legend.position="none") +
  ylab("Bray-Curtis diversity") +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=12), axis.text.y=element_text(size=12)) +
  theme_bw()

ggplot(df.bray) +
 aes(x = L1, y = value, fill = L1) +
 # geom_jitter(width=0.2) + 
  geom_violin(draw_quantiles = c(0.5)) +
  ylab("Bray Curtis dissimilarity") +
  xlab("") +
 scale_fill_manual(values = c(HE = "#2d8b92", HS = "#bd801c")) +
 theme_bw()

ggplot(df.bray) +
  aes(x = L1, y = value, fill = L1) +
  # geom_jitter(width=0.2) + 
  geom_boxplot(draw_quantiles = c(0.5)) +
  ylab("Bray Curtis dissimilarity") +
  xlab("") +
  scale_fill_manual(values = c(HE = "#2d8b92", HS = "#bd801c")) +
  theme_bw()

wilcox.test(value~L1, data =df.bray,alternative = "two.sided",  paired = FALSE, exact = F)
effect_size <- wilcox_effsize(df.bray,value~L1)
(effect_size$effsize)^2
effect_size$n1 
effect_size$n2 

# Multiple testing correction 
P_values <-  c(0.715,2.2e-16)
P_values
p.adjust(P_values, method = "BH", 2)

# Comes from this paper: https://www.frontiersin.org/articles/10.3389/fmicb.2019.01743/full
####################################
# 6.8.2 Hurley
####################################
relab_genera = transform_sample_counts(phyloseq_rarefied_phages, function(x) x / sum(x) * 100) 
head(otu_table(relab_genera)[,1:6])

abrel_bray <- phyloseq::distance(relab_genera, method = "bray")
abrel_bray <- as.matrix(abrel_bray)
head(abrel_bray)[,1:6]

sub_dist <- list()
groups_all <- sample_data(relab_genera)$Hurley
groups_all

unique(groups_all)

for (group in unique(groups_all)) { 
  row_group <- which(groups_all == group)
  sample_group <- sample_names(relab_genera)[row_group]
  sub_dist[[group]] <- abrel_bray[ sample_group, sample_group]
  sub_dist[[group]][!lower.tri(sub_dist[[group]])] <- NA
}

braygroups<- melt(sub_dist)
df.bray <- braygroups[complete.cases(braygroups), ]
df.bray$L1 <- factor(df.bray$L1, levels=names(sub_dist))
head(df.bray)

ggplot(df.bray, aes(x=L1, y=value, colour=L1)) +
  geom_jitter() + 
  geom_boxplot(alpha=0.6) +  
  theme(legend.position="none") +
  ylab("Bray-Curtis diversity") +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=12), axis.text.y=element_text(size=12)) +
  theme_bw()

ggplot(df.bray) +
  aes(x = L1, y = value, fill = L1) +
  geom_boxplot() +
  ylab("Bray Curtis dissimilarity") +
  xlab("") +
  scale_fill_manual(values = c("0" = "#499B78", "1"="blue", "2"="orange","3"="yellow")) +
  theme_bw()

kruskal.test(value ~ L1, data = df.bray)
dunnTest(value ~ L1, data = df.bray, method="bh")
####################################
# 6.8.3 BMI
####################################
relab_genera = transform_sample_counts(phyloseq_rarefied_phages, function(x) x / sum(x) * 100) 
head(otu_table(relab_genera)[,1:6])

abrel_bray <- phyloseq::distance(relab_genera, method = "bray")
abrel_bray <- as.matrix(abrel_bray)
head(abrel_bray)[,1:6]

sub_dist <- list()
groups_all <- sample_data(relab_genera)$Obese
groups_all

unique(groups_all)

for (group in unique(groups_all)) { 
  row_group <- which(groups_all == group)
  sample_group <- sample_names(relab_genera)[row_group]
  sub_dist[[group]] <- abrel_bray[ sample_group, sample_group]
  sub_dist[[group]][!lower.tri(sub_dist[[group]])] <- NA
}

braygroups<- melt(sub_dist)
df.bray <- braygroups[complete.cases(braygroups), ]
df.bray$L1 <- factor(df.bray$L1, levels=names(sub_dist))
head(df.bray)

ggplot(df.bray, aes(x=L1, y=value, colour=L1)) +
  geom_jitter() + 
  geom_boxplot(alpha=0.6) +  
  theme(legend.position="none") +
  ylab("Bray-Curtis diversity") +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=12), axis.text.y=element_text(size=12)) +
  theme_bw()

ggplot(df.bray) +
  aes(x = L1, y = value, fill = L1) +
  geom_boxplot() +
  ylab("Bray Curtis dissimilarity") +
  xlab("") +
  scale_fill_manual(values = c("Obese" = "orange", "Overweight"="brown", "Normal"="darkgreen")) +
  theme_bw()

kruskal.test(value ~ L1, data = df.bray)
dunnTest(value ~ L1, data = df.bray, method="bh")
####################################
# 6.9 Sharing of viral clusters
####################################
# 6.9.1 All individuals
####################################
# Input phyloseq_rarefied_phages

# DF
#RA_matrix_p <- as.data.frame(otu_table(ps1.fam.rel))
RA_matrix_p <- as.data.frame(otu_table(phyloseq_rarefied_phages))
presenceabsencep<-RA_matrix_p
presenceabsencep[presenceabsencep>0]<-1

# Add metadata
presenceabsence_agg_metp <- merge(t(presenceabsencep), metadata_patient, by = 0, all =F)
rownames(presenceabsence_agg_metp) <-  presenceabsence_agg_metp$Row.names
presenceabsence_agg_metp$Row.names <- NULL
View(presenceabsence_agg_metp)

# Calculate shared phage families
pre <-aggregate(presenceabsence_agg_metp[, 1:nrow(presenceabsencep)], by=list(presenceabsence_agg_metp$Patient_ID),FUN=max)
pre[,-1]<-lapply(pre[,-1],as.numeric)
rowSums(unique(pre[,-1]))
pre$totalcontigs <-rowSums(unique(pre[,-1]))

#find sharing of contigs between infants
sharedcontigs <-as.data.frame(colSums(pre[,-1]))
sharedcontigs <-  as.data.frame(t(sharedcontigs))
sharedcontigs$totalcontigs <- NULL
sharedcontigs <-  as.data.frame(t(sharedcontigs))
sharedcontigs_table <- as.data.frame(table(sharedcontigs$`colSums(pre[, -1])`))
colnames(sharedcontigs_table) <- c("Individuals", "Shared")
head(sharedcontigs_table)

# Figure
ggplot(sharedcontigs_table) +
  aes(x = Individuals, y = Shared) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values="cadetblue") +
  xlab("Number of individuals") +
  ylab("Number of viral clusters shared") +
  theme_classic()
###################################
# 6.9.2 Healthy individuals
####################################
# Input phyloseq_rarefied_phages

# DF
phyloseq_rarefied_phages_HE <- subset_samples(phyloseq_rarefied_phages, Disease == "HE")
#phyloseq_rarefied_phages_HE <- subset_samples(phyloseq_rarefied_phages, Disease == "HE" & Obese=="Normal") # 15 individuals
#View(sample_data(phyloseq_rarefied_phages_HE))

phyloseq_rarefied_phages_HE <- prune_samples(sample_sums(phyloseq_rarefied_phages_HE) > 0, phyloseq_rarefied_phages_HE)
RA_matrix_p <- as.data.frame(otu_table(phyloseq_rarefied_phages_HE))
presenceabsencep<-RA_matrix_p
presenceabsencep[presenceabsencep>0]<-1
presenceabsencep <- presenceabsencep[rowSums(presenceabsencep[, -1])>0, ]

# Which VC are shared highest?
presenceabsencep_Sharing <- presenceabsencep
presenceabsencep_Sharing$sharing <- rowSums(presenceabsencep_Sharing) # which viral 
presenceabsencep_Sharing <- presenceabsencep_Sharing[rev(order(presenceabsencep_Sharing$sharing)),]
rownames(presenceabsencep_Sharing)
#View(presenceabsencep_Sharing) # Most shared viral cluster 1, viral cluster 3; makes sense

# 7 VC are not shared with other HE, 39 VC are shared with 2 or more other HE.
(7/46)*100  # 15.2%
(39/46)*100 # 84.8%

(8/46)*100  # 17.4%
(38/46)*100 # 82.6%

# Merge with metadata
presenceabsence_agg_metp <- merge(t(presenceabsencep), metadata_patient, by = 0, all =F)
rownames(presenceabsence_agg_metp) <-  presenceabsence_agg_metp$Row.names
presenceabsence_agg_metp$Row.names <- NULL

# Calculate shared phage families
pre <-aggregate(presenceabsence_agg_metp[, 1:nrow(presenceabsencep)], by=list(presenceabsence_agg_metp$Patient_ID),FUN=max)
pre[,-1]<-lapply(pre[,-1],as.numeric)
rowSums(unique(pre[,-1]))
pre$totalcontigs <-rowSums(unique(pre[,-1]))

#find sharing of contigs between infants
sharedcontigs <-as.data.frame(colSums(pre[,-1]))
sharedcontigs <-  as.data.frame(t(sharedcontigs))
sharedcontigs$totalcontigs <- NULL
sharedcontigs <-  as.data.frame(t(sharedcontigs))
sharedcontigs_table <- as.data.frame(table(sharedcontigs$`colSums(pre[, -1])`))
colnames(sharedcontigs_table) <- c("Individuals", "Shared")
sharedcontigs_table
pre$totalcontig
pre$totalcontigs <- NULL
pre$Group.1 <- NULL
pre_DF <- as.data.frame(colSums(pre))
pre_DF$percentage <- (pre_DF$`colSums(pre)`/18)*100
View(pre_DF)

# Figure
ggplot(sharedcontigs_table) +
  aes(x = Individuals, y = Shared) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values="cadetblue") +
  xlab("Number of healthy individuals") +
  ylab("Number of viral clusters shared") +
  theme_classic() +
  geom_text(aes(label=Shared), position=position_dodge(width=0.9), vjust=-0.75, size=3.5)

# Total: 18 healthy individuals
(16/18)*100 # 88.9%  most common VC shared

# In BMI normalized 5 VC >= 50% of healthy individuals, which ones?
# shared by 8:f4,f10,f11,f30,f60
# shared by 9: f13,f2,f34
# shared by 10: f5
# shared by 11: f3
# shared by 14: f1
####################################
# 6.9.3 HS patients
####################################
# Input phyloseq_rarefied_phages

# DF
#phyloseq_rarefied_phages_HS <- subset_samples(phyloseq_rarefied_phages, Disease == "HS")
phyloseq_rarefied_phages_HS <- subset_samples(phyloseq_rarefied_family, Disease == "HS") # 12 patients
#phyloseq_rarefied_phages_HS_1 <- subset_samples(phyloseq_rarefied_phages_HS, Place=="Groin and genitals" | Place=="Axilliary")
#phyloseq_rarefied_phages_HS <- phyloseq_rarefied_phages_HS_1
View(sample_data(phyloseq_rarefied_phages_HS))
14+7+8

phyloseq_rarefied_phages_HS <- prune_samples(sample_sums(phyloseq_rarefied_phages_HS) > 0, phyloseq_rarefied_phages_HS)
RA_matrix_p <- as.data.frame(otu_table(phyloseq_rarefied_phages_HS))
presenceabsencep<-RA_matrix_p
presenceabsencep[presenceabsencep>0]<-1
presenceabsencep <- presenceabsencep[rowSums(presenceabsencep[, -1])>0, ]
#View(sample_data(presenceabsencep))

# Which VC are shared highest?
presenceabsencep_Sharing <- presenceabsencep
presenceabsencep_Sharing$sharing <- rowSums(presenceabsencep_Sharing) # which viral 
presenceabsencep_Sharing <- presenceabsencep_Sharing[rev(order(presenceabsencep_Sharing$sharing)),]
rownames(presenceabsencep_Sharing)
#View(presenceabsencep_Sharing)
#View(Metadata)

# 7 VC are not shared with other HE, 39 VC are shared with 2 or more other HE.
(48/123)*100  # 39.0%
(75/123)*100 # 61.0%

# 7 VC are not shared with other HE, 39 VC are shared with 2 or more other HE.
table(presenceabsencep_Sharing$sharing)
(33/71)*100  # 46.5%
(38/71)*100 # 53.5%

# Add metadata
presenceabsence_agg_metp <- merge(t(presenceabsencep), metadata, by = 0, all =F)
presenceabsence_agg_metp <- merge(t(presenceabsencep), metadata_patient, by = 0, all =F)
rownames(presenceabsence_agg_metp) <-  presenceabsence_agg_metp$Row.names
presenceabsence_agg_metp$Row.names <- NULL

# Calculate shared phage families
pre <-aggregate(presenceabsence_agg_metp[, 1:nrow(presenceabsencep)], by=list(presenceabsence_agg_metp$Patient_ID),FUN=max)
pre[,-1]<-lapply(pre[,-1],as.numeric)
rowSums(unique(pre[,-1]))
pre$totalcontigs <-rowSums(unique(pre[,-1]))
colSums(pre[-1])
pre$totalcontigs <- NULL
pre$Group.1 <- NULL
pre_DF <- as.data.frame(colSums(pre))
pre_DF$percentage <- (pre_DF$`colSums(pre)`/51)*100
View(pre_DF[!pre_DF$`colSums(pre)`==1,])

#find sharing of contigs between infants
sharedcontigs <-as.data.frame(colSums(pre[,-1]))
sharedcontigs <-  as.data.frame(t(sharedcontigs))
sharedcontigs$totalcontigs <- NULL
sharedcontigs <-  as.data.frame(t(sharedcontigs))
sharedcontigs_table <- as.data.frame(table(sharedcontigs$`colSums(pre[, -1])`))
colnames(sharedcontigs_table) <- c("Individuals", "Shared")
sharedcontigs_table

# Figure
ggplot(sharedcontigs_table) +
  aes(x = Individuals, y = Shared) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values="cadetblue") +
  xlab("Number of HS patients") +
  ylab("Number of viral clusters shared") +
  theme_classic() +
  geom_text(aes(label=Shared), position=position_dodge(width=0.9), vjust=-0.75, size=3.5)

# Total: 39 patients
(24/39)*100 # 61%  most common VC shared

# In BMI normalized 5 VC >= 50% of patients, which ones?
# shared by 9: f2
# shared by 7: f6,f5
# shared by 6: f1,f17,f3
####################################