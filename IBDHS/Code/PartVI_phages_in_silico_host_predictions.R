####################################
# SCRIPT 6: PROKARYOTIC VIRUSES
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
#install.packages("coin")
#install.packages("rstatix")
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
library(coin)
library(rstatix)
####################################
# 1. Host prediction (RaFah): Insert table
####################################
setwd("/Users/daan/Desktop/Transfer/Input_R/Host")
getwd()
dir()

Host_genus$names <- rownames(Host_genus)

Mastertable_phage_HQ <- Mastertable_phage[Mastertable_phage$Phage_score1=="phage",]
nrow(Mastertable_phage_HQ)

## Phylum
Host_phylum <- as.data.frame(read_delim("Host_prediction_phyla_final.txt", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE))
colnames(Host_phylum) <- c('contig', 'phylum')
rownames(Host_phylum) <- Host_phylum$contig
Host_phylum$contig <- NULL
rownames(Host_phylum) <- gsub("\\.{1}", "_", rownames(Host_phylum))

## Family
Host_family <- as.data.frame(read_delim("Host_prediction_family_final.txt", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE))
colnames(Host_family) <- c('contig', 'family')
rownames(Host_family) <- Host_family$contig
Host_family$contig <- NULL
rownames(Host_family) <- gsub("\\.{1}", "_", rownames(Host_family))

## Genus
Host_genus <- as.data.frame(read_delim("Host_prediction_genus_final.txt", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE))
colnames(Host_genus) <- c('contig', 'genus')
rownames(Host_genus) <- Host_genus$contig
Host_genus$contig <- NULL
rownames(Host_genus) <- gsub("\\.{1}", "_", rownames(Host_genus))
####################################
# 2 Merge with Mastertable
####################################
# Merge
Merged_tax <- transform(merge(Mastertable_phage_HQ, Host_phylum, by = 0, all = T), row.names=Row.names, Row.names = NULL)
Merged_tax <- transform(merge(Merged_tax, Host_family, by = 0, all = T), row.names=Row.names, Row.names = NULL)
Mastertable_viral_rarefied_phages_host <- transform(merge(Merged_tax, Host_genus, by = 0, all = T), row.names=Row.names, Row.names = NULL)
Mastertable_viral_rarefied_phages_host$phylum[is.na(Mastertable_viral_rarefied_phages_host$phylum)] <- "Not predicted"
Mastertable_viral_rarefied_phages_host$family[is.na(Mastertable_viral_rarefied_phages_host$family)] <- "Not predicted"
Mastertable_viral_rarefied_phages_host$genus[is.na(Mastertable_viral_rarefied_phages_host$genus)] <- "Not predicted"

# Curate (No Chlamydiae on skin)
Mastertable_viral_rarefied_phages_host$phylum[Mastertable_viral_rarefied_phages_host$phylum=="Chlamydiae"] <- "Not predicted"
Mastertable_viral_rarefied_phages_host$family[Mastertable_viral_rarefied_phages_host$family=="Chlamydiaceae"] <- "Not predicted"
Mastertable_viral_rarefied_phages_host$genus[Mastertable_viral_rarefied_phages_host$genus=="Chlamydia"] <- "Not predicted"

table(Mastertable_viral_rarefied_phages_host$phylum)
table(Mastertable_viral_rarefied_phages_host$family)
table(Mastertable_viral_rarefied_phages_host$genus)

# Percentage predicted
(nrow(Mastertable_viral_rarefied_phages_host[!Mastertable_viral_rarefied_phages_host$phylum=="Not predicted",])/nrow(Mastertable_viral_rarefied_phages_host[Mastertable_viral_rarefied_phages_host$Phage_score1=="phage",]))*100
(nrow(Mastertable_viral_rarefied_phages_host[!Mastertable_viral_rarefied_phages_host$family=="Not predicted",])/nrow(Mastertable_viral_rarefied_phages_host[Mastertable_viral_rarefied_phages_host$Phage_score1=="phage",]))*100
(nrow(Mastertable_viral_rarefied_phages_host[!Mastertable_viral_rarefied_phages_host$genus=="Not predicted",])/nrow(Mastertable_viral_rarefied_phages_host[Mastertable_viral_rarefied_phages_host$Phage_score1=="phage",]))*100
#View(Mastertable_viral_rarefied_phages_host)

# Mastertable_viral_rarefied_phages_host
ncol(Mastertable_viral_rarefied_phages_host)
Tax_table <- Mastertable_viral_rarefied_phages_host[tail(names(Mastertable_viral_rarefied_phages_host), 3)]
#View(Tax_table)
####################################
# 3. Load metadata
####################################
# Load sample metadata
setwd("/Users/daan/Desktop/Bioinformatics/pre-processing/HS project/Metadata/Final")
dir()
metadata <- read_excel("HS_Metadata_Final2_R.xlsx")
metadata <- as.data.frame(metadata)
rownames(metadata) <- metadata$Patient_ID

# Load patient metadata
setwd("/Users/daan/Desktop/Bioinformatics/pre-processing/HS project/Metadata/Final")
metadata_patient <- read_excel("HS_Metadata_Final2_Patient.xlsx")
metadata_patient <- as.data.frame(metadata_patient)
rownames(metadata_patient) <- metadata_patient$Patient_ID
metadata_patient$Hurley <- as.character(metadata_patient$Hurley)
####################################
# 4. Aggregate by patient_ID
####################################
# Aggregate by patient
Mastertable_viral_rarefied_phages_host_p <- Mastertable_viral_rarefied_phages_host[1:127] 
Mastertable_family_like_patient_pcoa <- t(Mastertable_viral_rarefied_phages_host_p)
Mastertable_family_like_patient_pcoa <- merge(Mastertable_family_like_patient_pcoa, metadata, by=0, all =F)
rownames(Mastertable_family_like_patient_pcoa) <- Mastertable_family_like_patient_pcoa$Row.names
Mastertable_family_like_patient_pcoa$Row.names <- NULL

# Find Vector
vectorA <- which(colnames(Mastertable_family_like_patient_pcoa) == "Patient_ID")-1
vector_patient <- colnames(Mastertable_family_like_patient_pcoa[, c(1:vectorA)])

# Aggregated by patient number
Mastertable_family_like_patient_pcoa_ <- aggregate(. ~HS, FUN = sum, data = Mastertable_family_like_patient_pcoa[,colnames(Mastertable_family_like_patient_pcoa) %in% vector_patient | colnames(Mastertable_family_like_patient_pcoa) == 'HS'])
rownames(Mastertable_family_like_patient_pcoa_) <- Mastertable_family_like_patient_pcoa_$HS
Mastertable_family_like_patient_pcoa_$HS <- NULL

# Remove column with sum zero 
Mastertable_family_like_patient_pcoa_ <-  as.data.frame(t(Mastertable_family_like_patient_pcoa_))
Mastertable_family_like_patient_pcoa_ <- Mastertable_family_like_patient_pcoa_[, which(colSums(Mastertable_family_like_patient_pcoa_) != 0)]
abundance_table <- Mastertable_family_like_patient_pcoa_
nrow(abundance_table) # 247 HQ-phages
ncol(abundance_table) # 43 patients left
####################################
# 5. Make phyloseq object with host
####################################
# 5.1 Abundance table
####################################
Taxonomy_table <- Tax_table

# 5.1 Create taxonomy table
Taxonomy_table$Kingdom <- "Unannotated"
Taxonomy_table$Class <- "Unannotated"
Taxonomy_table$Order <- "Unannotated"
Taxonomy_table$Species <- "Unannotated"
names(Taxonomy_table)[names(Taxonomy_table) == "phylum"] <- "Phylum"
names(Taxonomy_table)[names(Taxonomy_table) == "family"] <- "Family"
names(Taxonomy_table)[names(Taxonomy_table) == "genus"] <- "Genus"
Taxonomy_table <- Taxonomy_table[,c(4,1,5,6,2,3,7)]
####################################
# 5.2 Create matrix
####################################
abundance_table_rarefied_m <- as.matrix(abundance_table)
taxonomy_table_rarefied_m <- as.matrix(Taxonomy_table)
####################################
# 5.3 Transform to phyloseq objects
####################################
library(phyloseq)

ABUNDANCE_rarefied_family <- otu_table(abundance_table_rarefied_m, taxa_are_rows = TRUE)
TAX_rarefied_family <- tax_table(taxonomy_table_rarefied_m)
samples <- sample_data(metadata_patient)

phyloseq_rarefied_family <- phyloseq(ABUNDANCE_rarefied_family,TAX_rarefied_family, samples)
####################################
# 6. Relative abundance
####################################
# 6.1 Barplot
####################################
table(phyloseq::tax_table(phyloseq_rarefied_family)[, "Phylum"])
table(phyloseq::tax_table(phyloseq_rarefied_family)[, "Family"])
table(phyloseq::tax_table(phyloseq_rarefied_family)[, "Genus"])

ps_rel_abund = phyloseq::transform_sample_counts(phyloseq_rarefied_family, function(x){x / sum(x)})
phyloseq::otu_table(phyloseq_rarefied_family)[1:5, 1:5]
phyloseq::otu_table(ps_rel_abund)[1:5, 1:5]

phyloseq::plot_bar(ps_rel_abund, fill = "Phylum") +
  geom_bar(aes(color = Phylum, fill = Phylum), stat = "identity", position = "stack") +
  labs(x = "", y = "Relative Abundance\n") +
  facet_wrap(~ Disease, scales = "free") +
  theme(panel.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

phyloseq::plot_bar(ps_rel_abund, fill = "Phylum") +
  geom_bar(aes(color = Phylum, fill = Phylum), stat = "identity", position = "stack") +
  labs(x = "", y = "Relative Abundance\n") +
  facet_wrap(~ Hurley, scales = "free") +
  theme(panel.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

# Barplot no RA
phyloseq::plot_bar(phyloseq_rarefied_family, fill = "Phylum") +
  geom_bar(aes(color = Phylum, fill = Phylum), stat = "identity") +
  labs(x = "", y = "Relative Abundance\n") +
  facet_wrap(~ Disease, scales = "free") +
  theme(panel.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
####################################
# 6.2 Disease boxplot
####################################
# 6.2.1 Phyla
####################################
# Phyla
ps_phylum <- phyloseq::tax_glom(phyloseq_rarefied_family, "Phylum")
phyloseq::taxa_names(ps_phylum) <- phyloseq::tax_table(ps_phylum)[, "Phylum"]
phyloseq::otu_table(ps_phylum)[1:5, 1:5]
phyloseq_rarefied_phages <- ps_phylum

# Compositional
phyloseq_rarefied_phages_20 <- phyloseq_rarefied_phages %>% 
  microbiome::transform(transform = "compositional")

phyloseq_rarefied_phages_20_ <- phyloseq_rarefied_phages_20
phyloseq_rarefied_phages_20_DF <- as.data.frame(otu_table(phyloseq_rarefied_phages_20_))
phyloseq_rarefied_phages_20_DF$x1 <- rownames(phyloseq_rarefied_phages_20_DF)
phyloseq_rarefied_phages_20_DF$x1 <- NULL
#phyloseq_rarefied_phages_20_DF_disease <- filter(phyloseq_rarefied_phages_20_DF, x1 %in% names)   

# Merge metadata
Heatmap_family_RA_m <- merge(t(phyloseq_rarefied_phages_20_DF), metadata_patient, by = 0, all = F)
rownames(Heatmap_family_RA_m) <- Heatmap_family_RA_m$Row.names
Heatmap_family_RA_m$Row.names <- NULL
View(Heatmap_family_RA_m)

# Select columns
colnames(Heatmap_family_RA_m)
vector_1 <- (which(names(Heatmap_family_RA_m)== "Patient_ID"))
vector_2 <- (which(names(Heatmap_family_RA_m)== "Disease")) 
vector_3 <- (which(names(Heatmap_family_RA_m)== "Hurley")) 
Heatmap_family_RA_met <- Heatmap_family_RA_m[,c(1:vector_1,vector_2,vector_3)]

# Melt
Heatmap_family_RA_melt <- melt(Heatmap_family_RA_met, id.vars = c("Patient_ID","Disease", "Hurley"))
Heatmap_family_RA_melt$value <- as.numeric(Heatmap_family_RA_melt$value)

# Figure
ggplot(Heatmap_family_RA_melt) +
 aes(x = variable, y = value, fill = Disease) +
 geom_boxplot() +
 scale_fill_hue(direction = 1) +
 theme_bw()

# Statistics
# Proteobacteria
Statistics_f <- Heatmap_family_RA_melt[c(Heatmap_family_RA_melt$variable == "Proteobacteria"),]
wilcox.test(value~Disease, data =Statistics_f,alternative = "two.sided",  paired = FALSE, exact = F)
effect_size <- wilcox_effsize(Statistics_f,value~Disease)
(effect_size$effsize)^2
effect_size$n1 
effect_size$n2 

# Bacteroidetes
Statistics_f <- Heatmap_family_RA_melt[c(Heatmap_family_RA_melt$variable == "Bacteroidetes"),]
wilcox.test(value~Disease, data =Statistics_f,alternative = "two.sided",  paired = FALSE, exact = F)
effect_size <- wilcox_effsize(Statistics_f,value~Disease)
(effect_size$effsize)^2
effect_size$n1 
effect_size$n2 

# Actinobacteria
Statistics_f <- Heatmap_family_RA_melt[c(Heatmap_family_RA_melt$variable == "Actinobacteria"),]
wilcox.test(value~Disease, data =Statistics_f,alternative = "two.sided",  paired = FALSE, exact = F)
effect_size <- wilcox_effsize(Statistics_f,value~Disease)
(effect_size$effsize)^2
effect_size$n1 
effect_size$n2 

# Bacteroidetes
Statistics_f <- Heatmap_family_RA_melt[c(Heatmap_family_RA_melt$variable == "Firmicutes"),]
wilcox.test(value~Disease, data =Statistics_f,alternative = "two.sided",  paired = FALSE, exact = F)
effect_size <- wilcox_effsize(Statistics_f,value~Disease)
(effect_size$effsize)^2
effect_size$n1 
effect_size$n2 

# Not predicted
Statistics_f <- Heatmap_family_RA_melt[c(Heatmap_family_RA_melt$variable == "Not predicted"),]
wilcox.test(value~Disease, data =Statistics_f,alternative = "two.sided",  paired = FALSE, exact = F)
effect_size <- wilcox_effsize(Statistics_f,value~Disease)
(effect_size$effsize)^2
effect_size$n1 
effect_size$n2 

# Multiple testing correction 
P_values <-  c(0.0003781,0.01154,0.2042,0.1432,0.8946)
p.adjust(P_values, method = "BH")
####################################
# 6.2.2 Genera
####################################
# Genus
ps1.rel <- microbiome::transform(phyloseq_rarefied_family, "compositional")
ps1.fam.rel_gen <- aggregate_rare(ps1.rel, level = "Genus", detection = 0.1, prevalence = 0.01) 
phyloseq_rarefied_phages <- ps1.fam.rel_gen

phyloseq_rarefied_phages_20_ <- phyloseq_rarefied_phages
phyloseq_rarefied_phages_20_DF <- as.data.frame(otu_table(phyloseq_rarefied_phages_20_))
phyloseq_rarefied_phages_20_DF$x1 <- rownames(phyloseq_rarefied_phages_20_DF)
phyloseq_rarefied_phages_20_DF$x1 <- NULL
#phyloseq_rarefied_phages_20_DF_disease <- filter(phyloseq_rarefied_phages_20_DF, x1 %in% names)   

# Merge metadata
Heatmap_family_RA_m <- merge(t(phyloseq_rarefied_phages_20_DF), metadata_patient, by = 0, all = F)
rownames(Heatmap_family_RA_m) <- Heatmap_family_RA_m$Row.names
Heatmap_family_RA_m$Row.names <- NULL
#View(Heatmap_family_RA_m)

# Select columns
colnames(Heatmap_family_RA_m)
vector_1 <- (which(names(Heatmap_family_RA_m)== "Patient_ID"))
vector_2 <- (which(names(Heatmap_family_RA_m)== "Disease")) 
vector_3 <- (which(names(Heatmap_family_RA_m)== "Hurley")) 
Heatmap_family_RA_met <- Heatmap_family_RA_m[,c(1:vector_1,vector_2,vector_3)]

# Melt
Heatmap_family_RA_melt <- melt(Heatmap_family_RA_met, id.vars = c("Patient_ID","Disease", "Hurley"))
Heatmap_family_RA_melt$value <- as.numeric(Heatmap_family_RA_melt$value)

# Figure
ggplot(Heatmap_family_RA_melt) +
  aes(x = variable, y = value, fill = Disease) +
  geom_boxplot() +
  scale_fill_hue(direction = 1) +
  theme_bw()

# Statistics
# Not predicted
Statistics_f <- Heatmap_family_RA_melt[c(Heatmap_family_RA_melt$variable == "Not predicted"),]
wilcox.test(value~Disease, data =Statistics_f,alternative = "two.sided",  paired = FALSE, exact = F)
effect_size <- wilcox_effsize(Statistics_f,value~Disease)
(effect_size$effsize)^2
effect_size$n1 
effect_size$n2 

# Staphylococcus
Statistics_f <- Heatmap_family_RA_melt[c(Heatmap_family_RA_melt$variable == "Staphylococcus"),]
wilcox.test(value~Disease, data =Statistics_f,alternative = "two.sided",  paired = FALSE, exact = F)
effect_size <- wilcox_effsize(Statistics_f,value~Disease)
(effect_size$effsize)^2
effect_size$n1 
effect_size$n2 

# Streptococcus
Statistics_f <- Heatmap_family_RA_melt[c(Heatmap_family_RA_melt$variable == "Streptococcus"),]
wilcox.test(value~Disease, data =Statistics_f,alternative = "two.sided",  paired = FALSE, exact = F)
effect_size <- wilcox_effsize(Statistics_f,value~Disease)
(effect_size$effsize)^2
effect_size$n1 
effect_size$n2 

# Gordonia
Statistics_f <- Heatmap_family_RA_melt[c(Heatmap_family_RA_melt$variable == "Gordonia"),]
wilcox.test(value~Disease, data =Statistics_f,alternative = "two.sided",  paired = FALSE, exact = F)
effect_size <- wilcox_effsize(Statistics_f,value~Disease)
(effect_size$effsize)^2
effect_size$n1 
effect_size$n2 

# Enterococcus
Statistics_f <- Heatmap_family_RA_melt[c(Heatmap_family_RA_melt$variable == "Enterococcus"),]
wilcox.test(value~Disease, data =Statistics_f,alternative = "two.sided",  paired = FALSE, exact = F)
effect_size <- wilcox_effsize(Statistics_f,value~Disease)
(effect_size$effsize)^2
effect_size$n1 
effect_size$n2 

# Corynebacterium
Statistics_f <- Heatmap_family_RA_melt[c(Heatmap_family_RA_melt$variable == "Corynebacterium"),]
wilcox.test(value~Disease, data =Statistics_f,alternative = "two.sided",  paired = FALSE, exact = F)
effect_size <- wilcox_effsize(Statistics_f,value~Disease)
(effect_size$effsize)^2
effect_size$n1 
effect_size$n2 

# Cutibacterium
Statistics_f <- Heatmap_family_RA_melt[c(Heatmap_family_RA_melt$variable == "Cutibacterium"),]
wilcox.test(value~Disease, data =Statistics_f,alternative = "two.sided",  paired = FALSE, exact = F)
effect_size <- wilcox_effsize(Statistics_f,value~Disease)
(effect_size$effsize)^2
effect_size$n1 
effect_size$n2 

# Bacteroides
Statistics_f <- Heatmap_family_RA_melt[c(Heatmap_family_RA_melt$variable == "Bacteroides"),]
wilcox.test(value~Disease, data =Statistics_f,alternative = "two.sided",  paired = FALSE, exact = F)
effect_size <- wilcox_effsize(Statistics_f,value~Disease)
(effect_size$effsize)^2
effect_size$n1 
effect_size$n2 

# Cutibacterium
Statistics_f <- Heatmap_family_RA_melt[c(Heatmap_family_RA_melt$variable == "Clostridium"),]
wilcox.test(value~Disease, data =Statistics_f,alternative = "two.sided",  paired = FALSE, exact = F)
effect_size <- wilcox_effsize(Statistics_f,value~Disease)
(effect_size$effsize)^2
effect_size$n1 
effect_size$n2 

# Cutibacterium
Statistics_f <- Heatmap_family_RA_melt[c(Heatmap_family_RA_melt$variable == "Actinomyces"),]
wilcox.test(value~Disease, data =Statistics_f,alternative = "two.sided",  paired = FALSE, exact = F)
effect_size <- wilcox_effsize(Statistics_f,value~Disease)
(effect_size$effsize)^2
effect_size$n1 
effect_size$n2 

# Multiple testing correction 
P_values <-  c(0.0003781,0.01154,0.2042,0.1432,0.8946)
p.adjust(P_values, method = "BH")
####################################
# 6.2 Hurley boxplot
####################################
# 6.2.1 Phyla
####################################
# Phyla
ps_phylum <- phyloseq::tax_glom(phyloseq_rarefied_family, "Phylum")
phyloseq::taxa_names(ps_phylum) <- phyloseq::tax_table(ps_phylum)[, "Phylum"]
phyloseq::otu_table(ps_phylum)[1:5, 1:5]
phyloseq_rarefied_phages <- ps_phylum

# Compositional
phyloseq_rarefied_phages_20 <- phyloseq_rarefied_phages %>% 
  microbiome::transform(transform = "compositional")

phyloseq_rarefied_phages_20_ <- phyloseq_rarefied_phages_20
phyloseq_rarefied_phages_20_DF <- as.data.frame(otu_table(phyloseq_rarefied_phages_20_))
phyloseq_rarefied_phages_20_DF$x1 <- rownames(phyloseq_rarefied_phages_20_DF)
phyloseq_rarefied_phages_20_DF$x1 <- NULL
#phyloseq_rarefied_phages_20_DF_disease <- filter(phyloseq_rarefied_phages_20_DF, x1 %in% names)   

# Merge metadata
Heatmap_family_RA_m <- merge(t(phyloseq_rarefied_phages_20_DF), metadata_patient, by = 0, all = F)
rownames(Heatmap_family_RA_m) <- Heatmap_family_RA_m$Row.names
Heatmap_family_RA_m$Row.names <- NULL
#View(Heatmap_family_RA_m)

# Select columns
colnames(Heatmap_family_RA_m)
vector_1 <- (which(names(Heatmap_family_RA_m)== "Patient_ID"))
vector_2 <- (which(names(Heatmap_family_RA_m)== "Disease")) 
vector_3 <- (which(names(Heatmap_family_RA_m)== "Hurley")) 
Heatmap_family_RA_met <- Heatmap_family_RA_m[,c(1:vector_1,vector_2,vector_3)]

# Melt
Heatmap_family_RA_melt <- melt(Heatmap_family_RA_met, id.vars = c("Patient_ID","Disease", "Hurley"))
Heatmap_family_RA_melt$value <- as.numeric(Heatmap_family_RA_melt$value)

# Figure
ggplot(Heatmap_family_RA_melt) +
  aes(x = variable, y = value, fill = Hurley) +
  geom_boxplot() +
  scale_fill_hue(direction = 1) +
  theme_bw()

# Statistics
# Proteobacteria
library("FSA")
Statistics_f <- Heatmap_family_RA_melt[c(Heatmap_family_RA_melt$variable == "Proteobacteria"),]
kruskal.test(value ~ Hurley, data = Statistics_f)
effect_size <- kruskal_effsize(Statistics_f,value~Hurley)
effect_size$effsize
effect_size$n
dunnTest(value ~ Hurley, data = Statistics_f, method="bh")

# Bacteroidetes
Statistics_f <- Heatmap_family_RA_melt[c(Heatmap_family_RA_melt$variable == "Bacteroidetes"),]
kruskal.test(value ~ Hurley, data = Statistics_f)
effect_size <- kruskal_effsize(Statistics_f,value~Hurley)
effect_size$effsize
effect_size$n
dunnTest(value ~ Hurley, data = Statistics_f, method="bh")

# Actinobacteria
Statistics_f <- Heatmap_family_RA_melt[c(Heatmap_family_RA_melt$variable == "Actinobacteria"),]
kruskal.test(value ~ Hurley, data = Statistics_f)
effect_size <- kruskal_effsize(Statistics_f,value~Hurley)
effect_size$effsize
effect_size$n
dunnTest(value ~ Hurley, data = Statistics_f, method="bh")

# Firmicutes
Statistics_f <- Heatmap_family_RA_melt[c(Heatmap_family_RA_melt$variable == "Firmicutes"),]
kruskal.test(value ~ Hurley, data = Statistics_f)
effect_size <- kruskal_effsize(Statistics_f,value~Hurley)
effect_size$effsize
effect_size$n
dunnTest(value ~ Hurley, data = Statistics_f, method="bh")

# Not predicted
Statistics_f <- Heatmap_family_RA_melt[c(Heatmap_family_RA_melt$variable == "Not predicted"),]
kruskal.test(value ~ Hurley, data = Statistics_f)
effect_size <- kruskal_effsize(Statistics_f,value~Hurley)
effect_size$effsize
effect_size$n
dunnTest(value ~ Hurley, data = Statistics_f, method="bh")

# Multiple testing correction 
P_values <-  c(0.003315,0.02243,0.3761,0.0586,0.2353)
p.adjust(P_values, method = "BH")
####################################
# 6.2.2 Genera
####################################
# Genus
ps1.rel <- microbiome::transform(phyloseq_rarefied_family, "compositional")
ps1.fam.rel_gen <- aggregate_rare(ps1.rel, level = "Genus", detection = 0.1, prevalence = 0.01) 
phyloseq_rarefied_phages <- ps1.fam.rel_gen

phyloseq_rarefied_phages_20_ <- phyloseq_rarefied_phages
phyloseq_rarefied_phages_20_DF <- as.data.frame(otu_table(phyloseq_rarefied_phages_20_))
phyloseq_rarefied_phages_20_DF$x1 <- rownames(phyloseq_rarefied_phages_20_DF)
phyloseq_rarefied_phages_20_DF$x1 <- NULL
#phyloseq_rarefied_phages_20_DF_disease <- filter(phyloseq_rarefied_phages_20_DF, x1 %in% names)   

# Merge metadata
Heatmap_family_RA_m <- merge(t(phyloseq_rarefied_phages_20_DF), metadata_patient, by = 0, all = F)
rownames(Heatmap_family_RA_m) <- Heatmap_family_RA_m$Row.names
Heatmap_family_RA_m$Row.names <- NULL
#View(Heatmap_family_RA_m)

# Select columns
colnames(Heatmap_family_RA_m)
vector_1 <- (which(names(Heatmap_family_RA_m)== "Patient_ID"))
vector_2 <- (which(names(Heatmap_family_RA_m)== "Disease")) 
vector_3 <- (which(names(Heatmap_family_RA_m)== "Hurley")) 
Heatmap_family_RA_met <- Heatmap_family_RA_m[,c(1:vector_1,vector_2,vector_3)]

# Melt
Heatmap_family_RA_melt <- melt(Heatmap_family_RA_met, id.vars = c("Patient_ID","Disease", "Hurley"))
Heatmap_family_RA_melt$value <- as.numeric(Heatmap_family_RA_melt$value)

# Figure
ggplot(Heatmap_family_RA_melt) +
  aes(x = variable, y = value, fill = Hurley) +
  geom_boxplot() +
  scale_fill_hue(direction = 1) +
  theme_bw()

# Statistics
# Not predicted
Statistics_f <- Heatmap_family_RA_melt[c(Heatmap_family_RA_melt$variable == "Not predicted"),]
kruskal.test(value ~ Hurley, data = Statistics_f)
effect_size <- kruskal_effsize(Statistics_f,value~Hurley)
effect_size$effsize
effect_size$n
dunnTest(value ~ Hurley, data = Statistics_f, method="bh")

# Staphylococcus
Statistics_f <- Heatmap_family_RA_melt[c(Heatmap_family_RA_melt$variable == "Staphylococcus"),]
kruskal.test(value ~ Hurley, data = Statistics_f)
effect_size <- kruskal_effsize(Statistics_f,value~Hurley)
effect_size$effsize
effect_size$n
dunnTest(value ~ Hurley, data = Statistics_f, method="bh")

# Streptococcus
Statistics_f <- Heatmap_family_RA_melt[c(Heatmap_family_RA_melt$variable == "Streptococcus"),]
kruskal.test(value ~ Hurley, data = Statistics_f)
effect_size <- kruskal_effsize(Statistics_f,value~Hurley)
effect_size$effsize
effect_size$n
dunnTest(value ~ Hurley, data = Statistics_f, method="bh")

# Gordonia
Statistics_f <- Heatmap_family_RA_melt[c(Heatmap_family_RA_melt$variable == "Gordonia"),]
kruskal.test(value ~ Hurley, data = Statistics_f)
effect_size <- kruskal_effsize(Statistics_f,value~Hurley)
effect_size$effsize
effect_size$n
dunnTest(value ~ Hurley, data = Statistics_f, method="bh")

# Enterococcus
Statistics_f <- Heatmap_family_RA_melt[c(Heatmap_family_RA_melt$variable == "Enterococcus"),]
kruskal.test(value ~ Hurley, data = Statistics_f)
effect_size <- kruskal_effsize(Statistics_f,value~Hurley)
effect_size$effsize
effect_size$n
dunnTest(value ~ Hurley, data = Statistics_f, method="bh")

# Corynebacterium
Statistics_f <- Heatmap_family_RA_melt[c(Heatmap_family_RA_melt$variable == "Corynebacterium"),]
kruskal.test(value ~ Hurley, data = Statistics_f)
effect_size <- kruskal_effsize(Statistics_f,value~Hurley)
effect_size$effsize
effect_size$n
dunnTest(value ~ Hurley, data = Statistics_f, method="bh")

# Cutibacterium
Statistics_f <- Heatmap_family_RA_melt[c(Heatmap_family_RA_melt$variable == "Cutibacterium"),]
kruskal.test(value ~ Hurley, data = Statistics_f)
effect_size <- kruskal_effsize(Statistics_f,value~Hurley)
effect_size$effsize
effect_size$n
dunnTest(value ~ Hurley, data = Statistics_f, method="bh")

# Bacteroides
Statistics_f <- Heatmap_family_RA_melt[c(Heatmap_family_RA_melt$variable == "Bacteroides"),]
kruskal.test(value ~ Hurley, data = Statistics_f)
effect_size <- kruskal_effsize(Statistics_f,value~Hurley)
effect_size$effsize
effect_size$n
dunnTest(value ~ Hurley, data = Statistics_f, method="bh")

# Clostridium
Statistics_f <- Heatmap_family_RA_melt[c(Heatmap_family_RA_melt$variable == "Clostridium"),]
kruskal.test(value ~ Hurley, data = Statistics_f)
effect_size <- kruskal_effsize(Statistics_f,value~Hurley)
effect_size$effsize
effect_size$n
dunnTest(value ~ Hurley, data = Statistics_f, method="bh")

# Actinomyces
Statistics_f <- Heatmap_family_RA_melt[c(Heatmap_family_RA_melt$variable == "Actinomyces"),]
kruskal.test(value ~ Hurley, data = Statistics_f)
effect_size <- kruskal_effsize(Statistics_f,value~Hurley)
effect_size$effsize
effect_size$n
dunnTest(value ~ Hurley, data = Statistics_f, method="bh")

# Multiple testing correction 
P_values <-  c(0.01175,0.005951,0.01034,0.2813,0.6005,0.3266,0.6397,0.05042,0.02985,0.2859)
p.adjust(P_values, method = "BH")
####################################
# 6.3 Compositional barplot 
####################################
# 6.3.1 Phyla
####################################
library("microbiome")
phyloseq_rarefied_family

ps1.rel <- microbiome::transform(phyloseq_rarefied_family, "compositional")
ps1.fam.rel_phyl <- aggregate_rare(ps1.rel, level = "Phylum", detection = 0.1, prevalence = 0.01) 

plot_composition(ps1.fam.rel,average_by = "Disease") + 
  guides(fill = guide_legend(ncol = 1)) +
  labs(x = "Samples", 
       y = "Relative abundance",
       title = "Relative abundance data", 
       subtitle = "Subtitle",
       caption = "Caption text.")  +
  scale_fill_brewer("Family", palette = "Paired") +
  scale_y_continuous(labels = percent) +
  theme_bw() 

plot_composition(ps1.fam.rel,group_by = "Disease") + 
  guides(fill = guide_legend(ncol = 1)) +
  labs(x = "Samples", 
       y = "Relative abundance",
       title = "Relative abundance data", 
       subtitle = "Subtitle",
       caption = "Caption text.",)  +
  scale_fill_brewer("Family", palette = "Paired") +
  theme_bw() 

plot_composition(ps1.fam.rel_phyl,average_by = "Hurley") + 
  guides(fill = guide_legend(ncol = 1)) +
  labs(x = "Samples", 
       y = "Relative abundance",
       title = "Relative abundance data", 
       subtitle = "Subtitle",
       caption = "Caption text.")  +
  scale_fill_brewer("Family", palette = "Paired") +
  scale_y_continuous(labels = percent) +
  theme_bw() 
####################################
# 6.3.2 Family
####################################
ps1.rel <- microbiome::transform(phyloseq_rarefied_family, "compositional")
ps1.fam.rel_fami <- aggregate_rare(ps1.rel, level = "Family", detection = 0.05, prevalence = 0.01) 

plot_composition(ps1.fam.rel,average_by = "Disease") + 
  guides(fill = guide_legend(ncol = 1)) +
  labs(x = "Samples", 
       y = "Relative abundance",
       title = "Relative abundance data", 
       subtitle = "Subtitle",
       caption = "Caption text.") 
  scale_fill_brewer("Family", palette = "Paired") +
  scale_y_continuous(labels = percent) +
  theme_bw() 
  
plot_composition(ps1.fam.rel,group_by = "Disease") + 
    guides(fill = guide_legend(ncol = 1)) +
    labs(x = "Samples", 
         y = "Relative abundance",
         title = "Relative abundance data", 
         subtitle = "Subtitle",
         caption = "Caption text.") 
  scale_fill_brewer("Family", palette = "Paired") +
    scale_y_continuous(labels = percent) +
    theme_bw() 
  
plot_composition(ps1.fam.rel_fami,average_by = "Hurley") + 
    guides(fill = guide_legend(ncol = 1)) +
    labs(x = "Samples", 
         y = "Relative abundance",
         title = "Relative abundance data", 
         subtitle = "Subtitle",
         caption = "Caption text.") 
# scale_fill_brewer("Family", palette = "Paired") +
#  scale_y_continuous(labels = percent) +
#  theme_bw() 

plot_composition(ps1.fam.rel,group_by = "Hurley") + 
  guides(fill = guide_legend(ncol = 1)) +
  labs(x = "Samples", 
       y = "Relative abundance",
       title = "Relative abundance data", 
       subtitle = "Subtitle",
       caption = "Caption text.") 
# scale_fill_brewer("Family", palette = "Paired") +
#  scale_y_continuous(labels = percent) +
#  theme_bw() 
####################################
# 6.3.3 Genus
####################################
ps1.rel <- microbiome::transform(phyloseq_rarefied_family, "compositional")
ps1.fam.rel_gen <- aggregate_rare(ps1.rel, level = "Genus", detection = 0.1, prevalence = 0.01) 

plot_composition(ps1.fam.rel_gen,average_by = "Disease") +
  guides(fill = guide_legend(ncol = 1)) +
  labs(x = "Samples", 
       y = "Relative abundance",
       title = "Relative abundance data", 
       subtitle = "Subtitle",
       caption = "Caption text.") +
scale_fill_brewer("Family", palette = "Paired") +
scale_y_continuous(labels = percent) +
theme_bw() 

plot_composition(ps1.fam.rel_gen,average_by = "Hurley") + 
  guides(fill = guide_legend(ncol = 1)) +
  labs(x = "Samples", 
       y = "Relative abundance",
       title = "Relative abundance data", 
       subtitle = "Subtitle",
       caption = "Caption text.") +
scale_fill_brewer("Family", palette = "Paired") +
scale_y_continuous(labels = percent) +
theme_bw() 

plot_composition(ps1.fam.rel_gen,average_by = "Hurley") + 
  guides(fill = guide_legend(ncol = 1)) +
  labs(x = "Samples", 
       y = "Relative abundance",
       title = "Relative abundance data", 
       subtitle = "Subtitle",
       caption = "Caption text.") +
  scale_fill_brewer("Family", palette = "Paired") +
  scale_y_continuous(labels = percent) +
  theme_bw() 


plot_composition(ps1.fam.rel,group_by = "Disease") + 
  guides(fill = guide_legend(ncol = 1)) +
  labs(x = "Samples", 
       y = "Relative abundance",
       title = "Relative abundance data", 
       subtitle = "Subtitle",
       caption = "Caption text.") 
scale_fill_brewer("Family", palette = "Paired") +
  scale_y_continuous(labels = percent) +
  theme_bw() 

plot_composition(ps1.fam.rel,group_by = "Hurley") + 
  guides(fill = guide_legend(ncol = 1)) +
  labs(x = "Samples", 
       y = "Relative abundance",
       title = "Relative abundance data", 
       subtitle = "Subtitle",
       caption = "Caption text.") 
# scale_fill_brewer("Family", palette = "Paired") +
#  scale_y_continuous(labels = percent) +
#  theme_bw() 
####################################
# 6.3.4 Combined figure
####################################
# Disease
plot_composition(ps1.fam.rel_phyl,average_by = "Disease") + 
  guides(fill = guide_legend(ncol = 1)) +
  labs(x = "Samples", 
       y = "Relative abundance",
       title = "Relative abundance data", 
       subtitle = "Subtitle",
       caption = "Caption text.")  +
  scale_fill_brewer("Family", palette = "Paired") +
  scale_y_continuous(labels = percent) +
  theme_bw() +
  
  plot_composition(ps1.fam.rel_gen,average_by = "Disease") + 
  guides(fill = guide_legend(ncol = 1)) +
  labs(x = "Samples", 
       y = "Relative abundance",
       title = "Relative abundance data", 
       subtitle = "Subtitle",
       caption = "Caption text.") +
  #  scale_fill_brewer("Family", palette = "Paired") +
  #  scale_y_continuous(labels = percent) +
  theme_bw() 

# Hurley
plot_composition(ps1.fam.rel_phyl,average_by = "Hurley") + 
  guides(fill = guide_legend(ncol = 1)) +
  labs(x = "Samples", 
       y = "Relative abundance",
       title = "Relative abundance data", 
       subtitle = "Subtitle",
       caption = "Caption text.")  +
  scale_fill_brewer("Family", palette = "Paired") +
  scale_y_continuous(labels = percent) +
  theme_bw() +

plot_composition(ps1.fam.rel_gen,average_by = "Hurley") + 
  guides(fill = guide_legend(ncol = 1)) +
  labs(x = "Samples", 
       y = "Relative abundance",
       title = "Relative abundance data", 
       subtitle = "Subtitle",
       caption = "Caption text.") +
#  scale_fill_brewer("Family", palette = "Paired") +
#  scale_y_continuous(labels = percent) +
  theme_bw() 

# 5% detection cutoff
####################################
