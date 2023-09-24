####################################
# SCRIPT 3: EUKARYOTIC VIRUSES
####################################
# Before starting this script
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
# install_github("jokergoo/ComplexHeatmap")
library(ComplexHeatmap)
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
library(rJava)
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
library(ggplot2)
library(tibble)
library(dplyr)
library(ape)
library(taxonomizr)
library(rJava)
library(UpSetR)
library(tidyverse)
library(venneuler)
library(grid)
####################################
# 1. Dataset
####################################
Mastertable # Reads & contigs figure
Mastertable_genus1 # Percentage classified on class, order, family level
Mastertable_eukaryotic # Eukaryotic viral figure
Mastertable_phage  # Phage analysis
####################################
# 2. Manual curation
####################################
Mastertable_viral_unrarefied_euk_virus <- Mastertable_eukaryotic
table(Mastertable_viral_unrarefied_euk_virus$Final_phylum)

# Family
table(Mastertable_viral_unrarefied_euk_virus$Final_family)
Mastertable_viral_unrarefied_euk_virus$Final_family[Mastertable_viral_unrarefied_euk_virus$Final_phylum == "Cressdnaviricota" & Mastertable_viral_unrarefied_euk_virus$Final_family=="Unannotated"] <- "Unclassified CRESS"

# Genus
table(Mastertable_viral_unrarefied_euk_virus$Final_genus)
Mastertable_viral_unrarefied_euk_virus$Final_genus[Mastertable_viral_unrarefied_euk_virus$Final_family == "Anelloviridae" & Mastertable_viral_unrarefied_euk_virus$Final_genus == "Unannotated"] <- "Unclassified Anelloviridae"
Mastertable_viral_unrarefied_euk_virus$Final_genus[Mastertable_viral_unrarefied_euk_virus$Final_family == "Genomoviridae" & Mastertable_viral_unrarefied_euk_virus$Final_genus == "Unannotated"] <- "Unclassified Genomoviridae"
Mastertable_viral_unrarefied_euk_virus$Final_genus[Mastertable_viral_unrarefied_euk_virus$Final_family == "Unclassified CRESS" & Mastertable_viral_unrarefied_euk_virus$Final_genus == "Unannotated"] <- "Unclassified CRESS"
####################################
# 3. Create a list of all eukaryotic virus
####################################
# 3.1 List of all eukaryotic viruses (~ family level)
####################################
Mastertable_viral_unrarefied_euk_virus
View(Mastertable_viral_unrarefied_euk_virus[Mastertable_viral_unrarefied_euk_virus$Final_family=="Anelloviridae",])
unique(Mastertable_viral_unrarefied_euk_virus$Final_family)
# ---------> Nr. of 13 unique eukaryotic viral families

unique(Mastertable_viral_unrarefied_euk_virus$Final_genus)
# ---------> Nr. of 17 unique eukaryotic viral genera
####################################
# 3.2 write down hosts
####################################
## We divide the above list into 3 sub-list: (1) animal-infecting viruses, (2) plant/fungi-infecting virus, (3) small circular viruses

table(Mastertable_viral_unrarefied_euk_virus$Final_family)

# 1) Human-infecting viruses
# Family
Adenoviridae

# 2) Plant/fungi-infecting viruses
# Family
Alphaflexiviridae
Chrysoviridae
Partitiviridae
Tombusviridae
Totiviridae
Tymoviridae
Virgaviridae

# 3) Small circular viruses
# Family
Anelloviridae
Genomoviridae
Geminiviridae
Smacoviridae
Unclassified CRESS
####################################
# 3.3 Add hosts to mastertable
####################################
# 3.3.1 Unrarefied mastertable 
####################################
## Adapt the lines here according to what you saw qua hosts before.

# (A) Create a new table to add 'group host annotation'
Mastertable_viral_unrarefied_euk_virus_host <- Mastertable_viral_unrarefied_euk_virus
unique(Mastertable_viral_unrarefied_euk_virus$Final_family)

# (B) Add 'group host annotation' in a new column for plant and fungal viruses
Mastertable_viral_unrarefied_euk_virus$Host <- Mastertable_viral_unrarefied_euk_virus$Final_family
Mastertable_viral_unrarefied_euk_virus$Host[Mastertable_viral_unrarefied_euk_virus$Final_family == "Alphaflexiviridae"] <- "plant & fungal viruses"
Mastertable_viral_unrarefied_euk_virus$Host[Mastertable_viral_unrarefied_euk_virus$Final_family == "Chrysoviridae"] <- "plant & fungal viruses"
Mastertable_viral_unrarefied_euk_virus$Host[Mastertable_viral_unrarefied_euk_virus$Final_family == "Partitiviridae"] <- "plant & fungal viruses"
Mastertable_viral_unrarefied_euk_virus$Host[Mastertable_viral_unrarefied_euk_virus$Final_family == "Tombusviridae"] <- "plant & fungal viruses"
Mastertable_viral_unrarefied_euk_virus$Host[Mastertable_viral_unrarefied_euk_virus$Final_family == "Totiviridae"] <- "plant & fungal viruses"
Mastertable_viral_unrarefied_euk_virus$Host[Mastertable_viral_unrarefied_euk_virus$Final_family == "Tymoviridae"] <- "plant & fungal viruses"
Mastertable_viral_unrarefied_euk_virus$Host[Mastertable_viral_unrarefied_euk_virus$Final_family == "Virgaviridae"] <- "plant & fungal viruses"

# (C) Add 'group host annotation' in a new column for animal viruses
Mastertable_viral_unrarefied_euk_virus$Host[Mastertable_viral_unrarefied_euk_virus$Final_family == "Adenoviridae"] <- "human viruses"

# (D) Add 'group host annotation' in a new column for small circular and unclassified viruses
Mastertable_viral_unrarefied_euk_virus$Host[Mastertable_viral_unrarefied_euk_virus$Final_family == "Anelloviridae"] <- "small circular viruses"
Mastertable_viral_unrarefied_euk_virus$Host[Mastertable_viral_unrarefied_euk_virus$Final_family == "Genomoviridae"] <- "small circular viruses"
Mastertable_viral_unrarefied_euk_virus$Host[Mastertable_viral_unrarefied_euk_virus$Final_family == "Geminiviridae"] <- "small circular viruses"
Mastertable_viral_unrarefied_euk_virus$Host[Mastertable_viral_unrarefied_euk_virus$Final_family == "Smacoviridae"] <- "small circular viruses"
Mastertable_viral_unrarefied_euk_virus$Host[Mastertable_viral_unrarefied_euk_virus$Final_family == "Unclassified CRESS"] <- "small circular viruses"

table(Mastertable_viral_unrarefied_euk_virus$Host)
nrow(Mastertable_viral_unrarefied_euk_virus)
unique(Mastertable_viral_unrarefied_euk_virus$Host)
#view(Mastertable_viral_unrarefied_euk_virus)
Mastertable_viral_unrarefied_euk_virus_host <- Mastertable_viral_unrarefied_euk_virus
####################################
# 3.3.2 Curated mastertables
#################################### 
Mastertable_viral_unrarefied_euk_virus_host 
####################################
# 4. Download the the curated eukaryotic viral list and save in excel format 
####################################
# Make sure here you download Family, genera level
# From the Unrarefied table & from rarefied table

setwd("/Users/daan/Desktop/Bioinformatics/Analysis/FMT/Data/1_Eukaryotic_viruses")
getwd()
dir()

vector_1 <- which(colnames(Mastertable_viral_unrarefied_euk_virus_host) == "Final_class")
vector_2 <- which(colnames(Mastertable_viral_unrarefied_euk_virus_host) == "Final_species")
vector_3 <- which(colnames(Mastertable_viral_unrarefied_euk_virus_host) == "Host")
vector_4 <- which(colnames(Mastertable_viral_unrarefied_euk_virus_host) == "length")
vector_5 <- which(colnames(Mastertable_viral_unrarefied_euk_virus_host) == "Totalnumberofreads")
vector_6 <- which(colnames(Mastertable_viral_unrarefied_euk_virus_host) == "Final_ANI")
vector_7 <- which(colnames(Mastertable_viral_unrarefied_euk_virus_host) == "Final_coverage")
vector_8 <- which(colnames(Mastertable_viral_unrarefied_euk_virus_host) == "Blastn_AS")
vector_9 <- which(colnames(Mastertable_viral_unrarefied_euk_virus_host) == "Diamond_AS")
vector_10 <- which(colnames(Mastertable_viral_unrarefied_euk_virus_host) == "Best_AS")

eukaryotic_viruses <- Mastertable_viral_unrarefied_euk_virus_host[, c(vector_1:vector_2, vector_3, vector_4, vector_5, vector_6, vector_7, vector_8,vector_9, vector_10)]
eukaryotic_viruses <- with(eukaryotic_viruses,  eukaryotic_viruses[order(Final_family) , ])
table(eukaryotic_viruses$Host)

# Vector
vector_euk_ <- which(colnames(Mastertable_viral_unrarefied_euk_virus_host) == "Virsorter")-1
vector_euk <- colnames(Mastertable_viral_unrarefied_euk_virus_host[, c(1:vector_euk_)])

library("writexl")
write_xlsx(eukaryotic_viruses,"./eukaryotic_viruses.xlsx")
####################################
# 5. Complex heatmap: hierarchical clustering based on euclidean/pearson
####################################
# 5.1 Preparation of matrix
####################################
# info: https://jokergoo.github.io/ComplexHeatmap-reference/book/introduction.html

# (A) Create an abundance table: groups as rownames for y-axis & sample ID as colnames for x-axis
table(Mastertable_viral_unrarefied_euk_virus_host$Final_species)

Mastertable_global_rarefied_fam_complex_heatm <- Mastertable_viral_unrarefied_euk_virus_host
Mastertable_global_rarefied_fam_complex_heatm <- aggregate(. ~Final_species, FUN = sum, data = Mastertable_global_rarefied_fam_complex_heatm[,colnames(Mastertable_global_rarefied_fam_complex_heatm) %in% vector_euk | colnames(Mastertable_global_rarefied_fam_complex_heatm) == 'Final_species'])
rownames(Mastertable_global_rarefied_fam_complex_heatm) <- Mastertable_global_rarefied_fam_complex_heatm$Final_species
Mastertable_global_rarefied_fam_complex_heatm$Final_species <- NULL

# (B) convert reads to relative abundances + add sample ID
Mastertable_global_rarefied_fam_complex_heatm <- sweep(Mastertable_global_rarefied_fam_complex_heatm, 2, colSums(Mastertable_global_rarefied_fam_complex_heatm), '/') 
Mastertable_global_rarefied_fam_complex_heatm[is.na(Mastertable_global_rarefied_fam_complex_heatm)] <- 0
Mastertable_global_rarefied_fam_complex_heatm <- Mastertable_global_rarefied_fam_complex_heatm[colSums(Mastertable_global_rarefied_fam_complex_heatm) > 0]
matrix_family <- as.matrix(Mastertable_global_rarefied_fam_complex_heatm)
#View(matrix_family)
####################################
# 5.2 Eukaryotic viral families
##################################
## (i) Create a percentage of the presence of a specific eukaryotic viral family over all the samples - by this measure you can see which families is the most prevalent in the samples. This is obviously based on all the column & not only on this subset.
Mastertable_global_rarefied_fam_complex_heatm <- Mastertable_viral_unrarefied_euk_virus_host
used_for_heatmap_later <- aggregate(. ~Final_species, FUN = sum, data = Mastertable_global_rarefied_fam_complex_heatm[,colnames(Mastertable_global_rarefied_fam_complex_heatm) %in% vector_euk | colnames(Mastertable_global_rarefied_fam_complex_heatm) == 'Final_species'])
Mastertable_global_rarefied_fam_complex_heatm <- aggregate(. ~Final_species, FUN = sum, data = Mastertable_global_rarefied_fam_complex_heatm[,colnames(Mastertable_global_rarefied_fam_complex_heatm) %in% vector_euk | colnames(Mastertable_global_rarefied_fam_complex_heatm) == 'Final_species'])
rownames(Mastertable_global_rarefied_fam_complex_heatm) <- Mastertable_global_rarefied_fam_complex_heatm$Final_species
Mastertable_global_rarefied_fam_complex_heatm$Final_species <- NULL
rownames(used_for_heatmap_later) <- used_for_heatmap_later$Final_species
used_for_heatmap_later$Final_species <- NULL
vector_eukfamily_names_row <- rownames(used_for_heatmap_later)
vector_eukfamily_names_row
ncol(used_for_heatmap_later) # 93 samples with colsums zero
colnames(used_for_heatmap_later)

# IBD
setwd("/Users/daan/Desktop/Bioinformatics/pre-processing/FMT_study/Metadata")
dir()

metadata <- read_excel("samples_FMT_project_Final_DJ_R.xlsx")
metadata <- as.data.frame(metadata)
rownames(metadata) <- metadata$Sample_ID
metadata$Sample_ID
#View(metadata)

Meta <- merge(t(used_for_heatmap_later),metadata,by=0,all=F)
rownames(Meta) <- Meta$Row.names
Meta$Row.names <- NULL
Meta$Disease_status[Meta$FMT_Treatment_type=="autologous donor"] <- "patient"
nrow(Meta)
rownames(Meta)
View(Meta)

percentage_families_row <- vector()  
for (l in 1: length(vector_eukfamily_names_row)){
  percentage_families_row[l] <- ((length(which(used_for_heatmap_later[l,] > 0)))/ncol(Mastertable[1:304]))*100
}

percentage_families_row 

View(Meta)

sum(Meta$`Pepper mild mottle virus`)


## (ii) Create a percentage of the number of different eukaryotic viral families over all families found within one sample
family_names_column <- colnames(matrix_family)
Mastertable_global_rarefied_fam_complex_heatm_dendo_reordered <- Mastertable_global_rarefied_fam_complex_heatm

richness_families_column <- vector()  
for (l in 1: length(family_names_column)){
  richness_families_column[l] <- (length(which(used_for_heatmap_later[,l] > 0)))
}

richness_families_column

## (iii) complex family heatmap: WITHOUT ROW DENDROGRAM (ALPHABETICALLY)
haDisease <- HeatmapAnnotation(df = data.frame(disease = Meta$`Disease status`),
                               col = list(disease = c("donor" = "lightgreen", "patient" = "orange")), show_annotation_name = TRUE,border = TRUE)


col_fun = colorRamp2(c(0,0.001,0.3,0.5,0.7,1), c("gray98","#ffdb00", "#ffa904","#ee7b06","#a12424","#400b0b"))
order_fam <- c("Adenoviridae","Alphaflexiviridae","Chrysoviridae","Partitiviridae","Tombusviridae","Totiviridae","Tymoviridae","Virgaviridae","Anelloviridae","Genomoviridae","Geminiviridae","Smacoviridae","Unclassified CRESS")
group_annotations <- c("human virus","plant and fungal virus","small circular virus","plant and fungal virus","small circular virus","small circular virus","plant and fungal virus","small circular virus","plant and fungal virus","plant and fungal virus","plant and fungal virus","small circular virus","plant and fungal virus")
col = c("plant and fungal virus" = "lightgreen", "human virus" = "orange", "small circular virus" = "#74a9cf")

Heatmap(matrix_family, name = "RA", col = col_fun,  heatmap_width = unit(20, "cm"), heatmap_height = unit(10, "cm"),
        column_title = "", column_title_gp = gpar(fontsize = 20, fontface = "bold"),
        cluster_rows = FALSE, cluster_columns = TRUE, clustering_distance_columns = "euclidean", column_dend_height = unit(1, "cm"),border = TRUE,
        row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 6),  column_names_rot = 45,
        right_annotation = rowAnnotation(groups = anno_simple(group_annotations, col = col, gp = gpar(col = "black")), bar2 = anno_points(percentage_families_row), gp = gpar(col = "black")),
        row_split = group_annotations, row_title = NULL,row_gap = unit(1, "mm"),
        top_annotation = haDisease)

ncol(matrix_family) # 93 samples have at least one euk. virus present
####################################
# 5.3. Aggregate by patients
####################################
Heatmap_family <- t(used_for_heatmap_later)
Heatmap_family <- merge(Heatmap_family, metadata, by=0, all =F)
rownames(Heatmap_family) <- Heatmap_family$Row.names
Heatmap_family$Row.names <- NULL
nrow(Heatmap_family)

# Find Vector
vectorA <- which(colnames(Heatmap_family) == "Sample_ID")-1
vector_patient <- colnames(Heatmap_family[, c(1:vectorA)])

# Aggregated by patient number
Heatmap_family <- aggregate(. ~Patient_ID, FUN = sum, data = Heatmap_family[,colnames(Heatmap_family) %in% vector_patient | colnames(Heatmap_family) == 'Patient_ID'])
rownames(Heatmap_family) <- Heatmap_family$Patient_ID
Heatmap_family$Patient_ID <- NULL
View(Heatmap_family)

# Remove column with sum zero 
Heatmap_family <-  as.data.frame(t(Heatmap_family))
Heatmap_family <- Heatmap_family[, which(colSums(Heatmap_family) != 0)]
View(Heatmap_family)

# Merge with Patient metadata
Heatmap_family_ <- t(Heatmap_family)

setwd("/Users/daan/Desktop/Bioinformatics/pre-processing/FMT_study/Metadata")
dir()

metadata_patient <- read_excel("samples_FMT_project_Final_DJ_Patient_R.xlsx")
metadata_patient <- as.data.frame(metadata_patient)
rownames(metadata_patient) <- metadata_patient$Patient_identifier
metadata_patient$`Disease status`[metadata_patient$FMT_Treatment_type=="autologous donor"] <- "patient"

Meta <- merge(Heatmap_family_,metadata_patient,by=0,all=F)
rownames(Meta) <- Meta$Row.names
Meta$Row.names <- NULL
View(Meta)

####################################
# 5.4 Heatmap patients
####################################
Complex_Heatmap_family <- t(Heatmap_family_)
Complex_Heatmap_family <- sweep(Complex_Heatmap_family, 2, colSums(Complex_Heatmap_family), '/') 
Complex_Heatmap_family[is.na(Complex_Heatmap_family)] <- 0
matrix_family <- as.matrix(Complex_Heatmap_family)
vector_eukfamily_names_row <- rownames(Complex_Heatmap_family)
ncol(matrix_family)
nrow(matrix_family)

View(Complex_Heatmap_family)
# Percentages eukaryotic viral families
percentage_families_row <- vector()  
for (l in 1: length(vector_eukfamily_names_row)){
  percentage_families_row[l] <- ((length(which(Complex_Heatmap_family[l,] > 0)))/ncol(Mastertable[1:57]))*100
}

percentage_families_row # There are only 57 patients in total and thus we subset 57 random columns from mastertable to be used as a total, and make prevalence calculations.

# Richness eukaryotic viral families
family_names_column <- colnames(matrix_family)
Mastertable_global_rarefied_fam_complex_heatm_dendo_reordered <- Heatmap_family

richness_families_column <- vector()  
for (l in 1: length(family_names_column)){
  richness_families_column[l] <- (length(which(Mastertable_global_rarefied_fam_complex_heatm_dendo_reordered[,l] > 0)))
}

richness_families_column

## Group
haDisease <- HeatmapAnnotation(df = data.frame(disease = Meta$`Disease status`),
                               col = list(disease = c("donor" = "lightgreen", "patient" = "orange")), show_annotation_name = TRUE,border = TRUE)


## (iii) complex family heatmap: WITHOUT ROW DENDROGRAM (ALPHABETICALLY)
col_fun = colorRamp2(c(0,0.001,0.3,0.5,0.7,1), c("gray98","#ffdb00", "#ffa904","#ee7b06","#a12424","#400b0b"))
order_fam <- c("Adenoviridae","Alphaflexiviridae","Chrysoviridae","Partitiviridae","Tombusviridae","Totiviridae","Tymoviridae","Virgaviridae","Anelloviridae","Genomoviridae","Geminiviridae","Smacoviridae","Unclassified CRESS")
group_annotations <- c("human virus","plant and fungal virus","small circular virus","plant and fungal virus","small circular virus","small circular virus","plant and fungal virus","small circular virus","plant and fungal virus","plant and fungal virus","plant and fungal virus","small circular virus","plant and fungal virus")
col = c("plant and fungal virus" = "lightgreen", "human virus" = "orange", "small circular virus" = "#74a9cf", "protozoan virus" = "red", "arthropod virus" = "yellow")

Heatmap(matrix_family, name = "RA", col = col_fun,  heatmap_width = unit(20, "cm"), heatmap_height = unit(10, "cm"),
        column_title = "", column_title_gp = gpar(fontsize = 20, fontface = "bold"),
        cluster_rows = FALSE, cluster_columns = TRUE, clustering_distance_columns = "euclidean", column_dend_height = unit(1, "cm"),border = TRUE,
        row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 6),  column_names_rot = 45,
        right_annotation = rowAnnotation(groups = anno_simple(group_annotations, col = col, gp = gpar(col = "black")), bar2 = anno_points(percentage_families_row), gp = gpar(col = "black")),
        row_split = group_annotations, row_title = NULL,row_gap = unit(1, "mm"),
        top_annotation = haDisease)

nrow(metadata_patient[metadata_patient$`Disease status`=="patient",])
nrow(metadata_patient[metadata_patient$`Disease status`=="donor",])

ncol(matrix_family)
View(Mastertable_eukaryotic)
####################################

######### !!!!!!!!!!!!!!! #########