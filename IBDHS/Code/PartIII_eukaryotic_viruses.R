####################################
# SCRIPT 3: EUKARYOTIC VIRUSES
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
# 1. Upload Metadata
#################################### 
# 1.1  Upload Metadata per sample
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
metadata <- metadata[,c(1,2,3,12)]
colnames(metadata)
#################################### 
# 1.2 Upload Metadata per patient
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
metadata_patient <- metadata_patient[,c(1,2,8)]
#View(metadata_patient)
#################################### 
# 2. Subset viral approaches
####################################
Mastertable_classical_approach <- Mastertable[Mastertable$Phage_score1=="phage" & !Mastertable$Final_viral=="eukaryotic virus",]
Mastertable_classical_approach$Final_class[Mastertable_classical_approach$Final_class=="Chrymotiviricetes" | Mastertable_classical_approach$Final_class==1] <- "Unannotated"
table(Mastertable_classical_approach$Final_class)
nrow(Mastertable_classical_approach)
((76+23)/(151+76+23))*100 
# 39.6% classified contigs of 250 contigs
(sum(Mastertable_classical_approach$Totalnumberofreads[!Mastertable_classical_approach$Final_class=="Unannotated"])/sum(Mastertable_classical_approach$Totalnumberofreads))*100
# 91.95% classified reads

Mastertable_family_like_approach <- Mastertable[Mastertable$Phage_score2=="phage_family" & !Mastertable$Final_viral=="eukaryotic virus",]
table(Mastertable_family_like_approach$Final_class)
unique(Mastertable_family_like_approach$Family_level)
((1962+60)/(1962+60+2662))*100 
# 43.2% classified contigs of 4684 contigs or 126 families
(sum(Mastertable_family_like_approach$Totalnumberofreads[!Mastertable_family_like_approach$Final_class=="Unannotated"])/sum(Mastertable_family_like_approach$Totalnumberofreads))*100
# 80.7% classified reads
####################################
# 3. Subset the 100 most abundant contigs or families
####################################
# Classical approach
Mastertable_classical_approach_100 <- Mastertable_classical_approach[rev(order(Mastertable_classical_approach$Totalnumberofreads)),]
Mastertable_classical_approach_100 <- Mastertable_classical_approach_100[1:100,]
(sum(Mastertable_classical_approach_100$Totalnumberofreads)/sum(Mastertable_classical_approach$Totalnumberofreads))*100 
# Most abundant 100 contigs cover 99.5% of reads
View(Mastertable_classical_approach_100)

# Family-like approach
Mastertable_family_like_approach_ <- Mastertable_family_like_approach
names <- colnames(Mastertable_family_like_approach_[1:127][,colSums(Mastertable_family_like_approach_[1:127])!=0])
names # 114 samples
Mastertable_family_like_approach_agg <- aggregate(. ~Family_level, FUN = sum, data = Mastertable_family_like_approach_[,colnames(Mastertable_family_like_approach_) %in% names | colnames(Mastertable_family_like_approach_) == 'Family_level'])
rownames(Mastertable_family_like_approach_agg) <- Mastertable_family_like_approach_agg$Family_level
Mastertable_family_like_approach_agg$Family_level <- NULL
Mastertable_family_like_approach_agg$totalreads <- rowSums(Mastertable_family_like_approach_agg)
Mastertable_family_like_approach_100 <- Mastertable_family_like_approach_agg[rev(order(Mastertable_family_like_approach_agg$totalreads)),]
Mastertable_family_like_approach_100 <- Mastertable_family_like_approach_100[1:100,]
Mastertable_family_like_approach_ <- Mastertable_family_like_approach_agg
(sum(Mastertable_family_like_approach_100$totalreads)/sum(Mastertable_family_like_approach_$totalreads))*100 
# Most abundant 100 family-like approach cover 99.9% of reads
Mastertable_family_like_approach_$totalreads <- NULL
Mastertable_family_like_approach_100$totalreads <- NULL
####################################
# 4. Tables for complex heatmap
####################################
# Classical approach
Mastertable_classical_approach
Mastertable_classical_approach_100

# Family-like approach:
Mastertable_family_like_approach_
Mastertable_family_like_approach_100
####################################
# 5. Complex heatmap of prokaryotic viruses: classical approach
####################################
Mastertable_viral_unrarefied_euk_virus <- Mastertable[Mastertable$Final_viral == "eukaryotic virus" & Mastertable$Final_coverage > 10,]
table(Mastertable_viral_unrarefied_euk_virus$Final_viral)
table(Mastertable_viral_unrarefied_euk_virus$Final_superkingdom)
####################################
# 6. Create a list of all eukaryotic virus
####################################
# 6.1 List of all eukaryotic viruses (~ family level)
####################################
Mastertable_viral_unrarefied_euk_virus

unique(Mastertable_viral_unrarefied_euk_virus$Final_family)
# ---------> Nr. of 12 unique eukaryotic viral families

unique(Mastertable_viral_unrarefied_euk_virus$Final_genus)
# ---------> Nr. of 24 unique eukaryotic viral genera
####################################
# 6.2 write down hosts
####################################
## We divide the above list into 3 sub-list: (1) animal-infecting viruses, (2) plant/fungi-infecting virus, (3) small circular viruses

Mastertable_viral_unrarefied_euk_virus <- Mastertable_viral_unrarefied_euk_virus[!Mastertable_viral_unrarefied_euk_virus$Final_family=="Retroviridae",]
table(Mastertable_viral_unrarefied_euk_virus$Final_family)
table(Mastertable_viral_unrarefied_euk_virus$Final_genus)

# --> Housemites

# Totiviridae
# --> Some infect Thrichomonas; Can be Protozoan-infecting virues

# 1) Human-infecting viruses
# Family
Papillomaviridae
Parvoviridae
Polyomaviridae
Reoviridae

# Genus
Alphapapillomavirus
Alphapolyomavirus
Betapapillomavirus
Bocaparvovirus
Deltapolyomavirus
Gammapapillomavirus
Rotavirus
Trichomonasvirus
unclassified Papillomaviridae

# 2) Arthopod-infecting viruses
# Family
Dicistroviridae

# Genus
unclassified Dicistroviridae

# 3) Protozoan-infecting viruses
# Family
Totiviridae

# Genus
Trichomonasvirus

# 4) Plant/fungi-infecting viruses
# Family
Alphaflexiviridae
Virgaviridae
Totiviridae

# Genus
Potexvirus
Tobamovirus
Totivirus
unclassified Totiviridae

# 5) Small circular viruses
# Family
Anelloviridae
Genomoviridae
Secoviridae

# Genus
Alphatorquevirus
Betatorquevirus
Gammatorquevirus
Gemyduguivirus
Gemykibivirus
Gemykrogvirus
Torradovirus
unclassified Anelloviridae
####################################
# 6.3 Add hosts to mastertable
####################################
## Adapt the lines here according to what you saw qua hosts before.

# (A) Create a new table to add 'group host annotation'
Mastertable_viral_unrarefied_euk_virus_host <- Mastertable_viral_unrarefied_euk_virus
unique(Mastertable_viral_unrarefied_euk_virus$Final_family)

# (B) Add 'group host annotation' in a new column for plant and fungal viruses
Mastertable_viral_unrarefied_euk_virus$Host <- Mastertable_viral_unrarefied_euk_virus$Final_family
Mastertable_viral_unrarefied_euk_virus$Host[Mastertable_viral_unrarefied_euk_virus$Final_family == "Alphaflexiviridae"] <- "plant & fungal viruses"
Mastertable_viral_unrarefied_euk_virus$Host[Mastertable_viral_unrarefied_euk_virus$Final_family == "Virgaviridae"] <- "plant & fungal viruses"
Mastertable_viral_unrarefied_euk_virus$Host[Mastertable_viral_unrarefied_euk_virus$Final_family == "Totiviridae"] <- "plant & fungal viruses"

# (C) Protozoan viruses
Mastertable_viral_unrarefied_euk_virus$Host[Mastertable_viral_unrarefied_euk_virus$Final_genus == "Trichomonasvirus"] <- "protozoan viruses"

# (D) 
Mastertable_viral_unrarefied_euk_virus$Host[Mastertable_viral_unrarefied_euk_virus$Final_family == "Dicistroviridae"] <- "arthropod viruses"

# (E) Add 'group host annotation' in a new column for animal viruses
Mastertable_viral_unrarefied_euk_virus$Host[Mastertable_viral_unrarefied_euk_virus$Final_family == "Papillomaviridae"] <- "human viruses"
Mastertable_viral_unrarefied_euk_virus$Host[Mastertable_viral_unrarefied_euk_virus$Final_family == "Parvoviridae"] <- "human viruses"
Mastertable_viral_unrarefied_euk_virus$Host[Mastertable_viral_unrarefied_euk_virus$Final_family == "Polyomaviridae"] <- "human viruses"
Mastertable_viral_unrarefied_euk_virus$Host[Mastertable_viral_unrarefied_euk_virus$Final_family == "Reoviridae"] <- "human viruses"

# (F) Add 'group host annotation' in a new column for small circular and unclassified viruses
Mastertable_viral_unrarefied_euk_virus$Host[Mastertable_viral_unrarefied_euk_virus$Final_family == "Anelloviridae"] <- "small circular viruses"
Mastertable_viral_unrarefied_euk_virus$Host[Mastertable_viral_unrarefied_euk_virus$Final_family == "Genomoviridae"] <- "small circular viruses"
Mastertable_viral_unrarefied_euk_virus$Host[Mastertable_viral_unrarefied_euk_virus$Final_family == "Secoviridae"] <- "small circular viruses"

unique(Mastertable_viral_unrarefied_euk_virus$Host)
#view(Mastertable_viral_unrarefied_euk_virus)
####################################
# 6.4 Download the the curated eukaryotic viral list 
####################################
setwd("/Users/daan/Desktop/Bioinformatics/Analysis/HS/HS_project/eukaryotic_viruses")
getwd()
dir()

#View(Mastertable_viral_unrarefied_euk_virus_host)
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
# 7. Complex heatmap: hierarchical clustering based on euclidean/pearson
####################################
# 7.1 Preparation of matrix
####################################
# (A) Create an abundance table: groups as rownames for y-axis & sample ID as colnames for x-axis
# (A.1) families
table(Mastertable_viral_unrarefied_euk_virus_host$Final_family)
table(Mastertable_viral_unrarefied_euk_virus_host$Final_genus)
table(Mastertable_viral_unrarefied_euk_virus_host$Host)

Mastertable_global_rarefied_fam_complex_heatm <- Mastertable_viral_unrarefied_euk_virus_host
Mastertable_global_rarefied_fam_complex_heatm <- aggregate(. ~Final_family, FUN = sum, data = Mastertable_global_rarefied_fam_complex_heatm[,colnames(Mastertable_global_rarefied_fam_complex_heatm) %in% vector_euk | colnames(Mastertable_global_rarefied_fam_complex_heatm) == 'Final_family'])
rownames(Mastertable_global_rarefied_fam_complex_heatm) <- Mastertable_global_rarefied_fam_complex_heatm$Final_family
Mastertable_global_rarefied_fam_complex_heatm$Final_family <- NULL

# (A.2) genera
Mastertable_global_rarefied_gen_complex_heatm <- Mastertable_viral_unrarefied_euk_virus_host
Mastertable_global_rarefied_gen_complex_heatm <- aggregate(. ~Final_genus, FUN = sum, data = Mastertable_global_rarefied_gen_complex_heatm[,colnames(Mastertable_global_rarefied_gen_complex_heatm) %in% vector_euk | colnames(Mastertable_global_rarefied_gen_complex_heatm) == 'Final_genus'])
rownames(Mastertable_global_rarefied_gen_complex_heatm) <- Mastertable_global_rarefied_gen_complex_heatm$Final_genus
Mastertable_global_rarefied_gen_complex_heatm$Final_genus <- NULL
#View(Mastertable_global_rarefied_gen_complex_heatm)

# (B) convert reads to relative abundances + add sample ID
# (B.1) families
Mastertable_global_rarefied_fam_complex_heatm <- sweep(Mastertable_global_rarefied_fam_complex_heatm, 2, colSums(Mastertable_global_rarefied_fam_complex_heatm), '/') 
Mastertable_global_rarefied_fam_complex_heatm[is.na(Mastertable_global_rarefied_fam_complex_heatm)] <- 0
Mastertable_global_rarefied_fam_complex_heatm <- Mastertable_global_rarefied_fam_complex_heatm[colSums(Mastertable_global_rarefied_fam_complex_heatm) > 0]
matrix_family <- as.matrix(Mastertable_global_rarefied_fam_complex_heatm)
colnames(matrix_family)

# (B.2) genera
Mastertable_global_rarefied_gen_complex_heatm <- sweep(Mastertable_global_rarefied_gen_complex_heatm, 2, colSums(Mastertable_global_rarefied_gen_complex_heatm), '/') 
Mastertable_global_rarefied_gen_complex_heatm[is.na(Mastertable_global_rarefied_gen_complex_heatm)] <- 0
Mastertable_global_rarefied_fam_complex_heatm <- Mastertable_global_rarefied_gen_complex_heatm[colSums(Mastertable_global_rarefied_gen_complex_heatm) > 0]
matrix_genus <- as.matrix(Mastertable_global_rarefied_gen_complex_heatm)
View(matrix_genus)

# (D) Create list for side annotations and create complex heatmap
####################################
# 7.2 Families
####################################
## (i) Create a percentage of the presence of a specific eukaryotic viral family over all the samples - by this measure you can see which families is the most prevalent in the samples. This is obviously based on all the column & not only on this subset.
#           --> I did a few extra step because I thought I needed the order for the sample_ID's (x-axis) after dendrogram, but this is not needed

Mastertable_global_rarefied_fam_complex_heatm <- Mastertable_viral_unrarefied_euk_virus_host
used_for_heatmap_later <- aggregate(. ~Final_family, FUN = sum, data = Mastertable_global_rarefied_fam_complex_heatm[,colnames(Mastertable_global_rarefied_fam_complex_heatm) %in% vector_euk | colnames(Mastertable_global_rarefied_fam_complex_heatm) == 'Final_family'])
Mastertable_global_rarefied_fam_complex_heatm <- aggregate(. ~Final_family, FUN = sum, data = Mastertable_global_rarefied_fam_complex_heatm[,colnames(Mastertable_global_rarefied_fam_complex_heatm) %in% vector_euk | colnames(Mastertable_global_rarefied_fam_complex_heatm) == 'Final_family'])
rownames(Mastertable_global_rarefied_fam_complex_heatm) <- Mastertable_global_rarefied_fam_complex_heatm$Final_family
Mastertable_global_rarefied_fam_complex_heatm$Final_family <- NULL
rownames(used_for_heatmap_later) <- used_for_heatmap_later$Final_family
used_for_heatmap_later$Final_family <- NULL
vector_eukfamily_names_row <- rownames(used_for_heatmap_later)

percentage_families_row <- vector()  
for (l in 1: length(vector_eukfamily_names_row)){
  percentage_families_row[l] <- ((length(which(Mastertable_global_rarefied_fam_complex_heatm[l,] > 0)))/ncol(Mastertable[1:127]))*100
}

percentage_families_row
# Put here total samples on 127 samples being 100% adapt later if needed

## (ii) Create a percentage of the number of different eukaryotic viral families over all families found within one sample
family_names_column <- colnames(matrix_family)
Mastertable_global_rarefied_fam_complex_heatm_dendo_reordered <- Mastertable_global_rarefied_fam_complex_heatm

richness_families_column <- vector()  
for (l in 1: length(family_names_column)){
  richness_families_column[l] <- (length(which(Mastertable_global_rarefied_fam_complex_heatm_dendo_reordered[,l] > 0)))
}

richness_families_column

## (iv) Add host group
Mastertable_global_rarefied_fam_complex_heatm <- Mastertable_viral_unrarefied_euk_virus_host
used_for_heatmap_later <- aggregate(. ~Final_family, FUN = sum, data = Mastertable_global_rarefied_fam_complex_heatm[,colnames(Mastertable_global_rarefied_fam_complex_heatm) %in% vector_euk | colnames(Mastertable_global_rarefied_fam_complex_heatm) == 'Final_family'])
Mastertable_global_rarefied_fam_complex_heatm <- aggregate(. ~Final_family, FUN = sum, data = Mastertable_global_rarefied_fam_complex_heatm[,colnames(Mastertable_global_rarefied_fam_complex_heatm) %in% vector_euk | colnames(Mastertable_global_rarefied_fam_complex_heatm) == 'Final_family'])
rownames(Mastertable_global_rarefied_fam_complex_heatm) <- Mastertable_global_rarefied_fam_complex_heatm$Final_family
Mastertable_global_rarefied_fam_complex_heatm$Final_family <- NULL
rownames(used_for_heatmap_later) <- used_for_heatmap_later$Final_family
used_for_heatmap_later$Final_family <- NULL

##  Make vector of host by subsetting instead of writing everything!
## Look for an easy way
## Also add disease annotation on top instead of richness/decide which one to go I say numbers
rownames(used_for_heatmap_later)

## (iii) complex family heatmap: WITHOUT ROW DENDROGRAM (ALPHABETICALLY)
col_fun = colorRamp2(c(0,0.001,0.3,0.5,0.7,1), c("gray98","#ffdb00", "#ffa904","#ee7b06","#a12424","#400b0b"))
order_fam <- c("Alphaflexiviridae","Virgaviridae","Totiviridae","Papillomaviridae","Parvoviridae","Polyomaviridae","Reoviridae","Anelloviridae","Genomoviridae","Secoviridae","Dicistroviridae","Trichomonasvirus")
group_annotations <- c("plant and fungal virus","small circular virus","arthropod virus","small circular virus","human virus","human virus","human virus","human virus","small circular virus","plant and fungal virus","protozoan virus","plant and fungal virus")
col = c("plant and fungal virus" = "lightgreen", "human virus" = "orange", "small circular virus" = "#74a9cf", "protozoan virus" = "red", "arthropod virus" = "yellow")

Heatmap(matrix_family, name = "RA", col = col_fun,  heatmap_width = unit(20, "cm"), heatmap_height = unit(10, "cm"),
        column_title = "", column_title_gp = gpar(fontsize = 20, fontface = "bold"),
        cluster_rows = FALSE, cluster_columns = TRUE, clustering_distance_columns = "euclidean", column_dend_height = unit(1, "cm"),border = TRUE,
        row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 6),  column_names_rot = 45,
        right_annotation = rowAnnotation(presence = percentage_families_row, groups = anno_simple(group_annotations, col = col, gp = gpar(col = "black")), bar2 = anno_points(percentage_families_row), gp = gpar(col = "black")),
        row_split = group_annotations, row_title = NULL,row_gap = unit(1, "mm"))

ncol(matrix_family) # 127 samples have at least one euk. virus present

# Finish in illustrators
#       --> Increase annotation size: make sure round circiles fall in borders.
#       --> remove titles next to bar & re-order legend next to figures
#       --> Make sure everything is in arial + Family level in Italic
####################################
# 7.3 Family: aggregated by patient_ID
####################################
# Merge with metadata samples
Heatmap_family <- t(used_for_heatmap_later)
Heatmap_family <- merge(Heatmap_family, metadata, by=0, all =F)
rownames(Heatmap_family) <- Heatmap_family$Row.names
Heatmap_family$Row.names <- NULL

# Find Vector
vectorA <- which(colnames(Heatmap_family) == "Patient_ID")-1
vector_patient <- colnames(Heatmap_family[, c(1:vectorA)])

# Aggregated by patient number
Heatmap_family <- aggregate(. ~HS, FUN = sum, data = Heatmap_family[,colnames(Heatmap_family) %in% vector_patient | colnames(Heatmap_family) == 'HS'])
rownames(Heatmap_family) <- Heatmap_family$HS
Heatmap_family$HS <- NULL
View(Heatmap_family)

# Remove column with sum zero 
Heatmap_family <-  as.data.frame(t(Heatmap_family))
Heatmap_family <- Heatmap_family[, which(colSums(Heatmap_family) != 0)]

# Merge with Patient metadata
Heatmap_family_ <- t(Heatmap_family)
Meta <- merge(Heatmap_family_,metadata_patient,by=0,all=F)
####################################
# 7.4 Complex heatmap
####################################
Complex_Heatmap_family <- Heatmap_family
Complex_Heatmap_family <- sweep(Complex_Heatmap_family, 2, colSums(Complex_Heatmap_family), '/') 
Complex_Heatmap_family[is.na(Complex_Heatmap_family)] <- 0
matrix_family <- as.matrix(Complex_Heatmap_family)
vector_eukfamily_names_row <- rownames(Complex_Heatmap_family)
ncol(matrix_family)
vector_eukfamily_names_row

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
Meta$Disease

order_fam_disease <- c("HE","HE","HE","HE","HE","HE","HE","HE","HE","HE",
                       "HS","HS","HS","HS","HS","HS","HS","HS","HS","HS",
                       "HS","HS","HS","HS","HS","HS","HS","HS","HS","HS",
                       "HS","HS","HS","HS","HS","HS","HS","HS")

order_fam_disease

haDisease <- HeatmapAnnotation(df = data.frame(disease = Meta$Disease),
                               col = list(disease = c("HE" = "lightgreen", "HS" = "orange")), show_annotation_name = TRUE,border = TRUE)


## (iii) complex family heatmap: WITHOUT ROW DENDROGRAM (ALPHABETICALLY)
col_fun = colorRamp2(c(0,0.001,0.3,0.5,0.7,1), c("gray98","#ffdb00", "#ffa904","#ee7b06","#a12424","#400b0b"))
order_fam <- c("Alphaflexiviridae","Virgaviridae","Totiviridae","Papillomaviridae","Parvoviridae","Polyomaviridae","Reoviridae","Anelloviridae","Genomoviridae","Secoviridae","Dicistroviridae","Trichomonasvirus")
group_annotations <- c("plant and fungal virus","small circular virus","arthropod virus","small circular virus","human virus","human virus","human virus","human virus","small circular virus","plant and fungal virus","protozoan virus","plant and fungal virus")
col = c("plant and fungal virus" = "lightgreen", "human virus" = "orange", "small circular virus" = "#74a9cf", "protozoan virus" = "red", "arthropod virus" = "yellow")

Heatmap(matrix_family, name = "RA", col = col_fun,  heatmap_width = unit(20, "cm"), heatmap_height = unit(10, "cm"),
        column_title = "", column_title_gp = gpar(fontsize = 20, fontface = "bold"),
        cluster_rows = FALSE, cluster_columns = TRUE, clustering_distance_columns = "euclidean", column_dend_height = unit(1, "cm"),border = TRUE,
        row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 6),  column_names_rot = 45,
        right_annotation = rowAnnotation(presence = percentage_families_row, groups = anno_simple(group_annotations, col = col, gp = gpar(col = "black")), bar2 = anno_points(percentage_families_row), gp = gpar(col = "black")),
        top_annotation = haDisease,
        column_split = order_fam_disease,
        row_split = group_annotations, row_title = NULL,row_gap = unit(1, "mm"))

# Boxplot
RA_matrix_family <- merge(t(matrix_family),metadata_patient,by=0,all=F)
rownames(RA_matrix_family) <- RA_matrix_family$Row.names
RA_matrix_family$Row.names <- NULL
colnames(RA_matrix_family)
vector_1 <- (which(names(RA_matrix_family)== "Patient_ID"))
vector_2 <- (which(names(RA_matrix_family)== "Disease")) 
vector_3 <- (which(names(RA_matrix_family)== "Hurley")) 
RA_matrix_family_m <- RA_matrix_family[,c(1:vector_1,vector_2,vector_3)]
Heatmap_family_RA_melt <- melt(RA_matrix_family_m, id.vars = c("Patient_ID","Disease", "Hurley"))
View(Heatmap_family_RA_melt)

## ggplot
Heatmap_family_RA_melt %>%
 filter(variable %in% c("Alphaflexiviridae", "Anelloviridae", "Totiviridae", "Virgaviridae")) %>%
 ggplot() +
 aes(x = variable, y = value, fill = Disease) +
 # geom_boxplot() +
  geom_violin() +
  scale_fill_manual(values = c(HE = "#12841B", HS = "#DD641F")) +
 theme_bw()

# Statistics
library(rstatix)
Statistics_f <- Heatmap_family_RA_melt[c(Heatmap_family_RA_melt$variable == "Anelloviridae"),]
wilcox.test(value~Disease, data =Statistics_f,alternative = "two.sided",  paired = FALSE, exact = F)
effect_size <- wilcox_effsize(Statistics_f,value~Disease)
(effect_size$effsize)^2
effect_size$n1 # n = HE
effect_size$n2 # n = HS

Statistics_f <- Heatmap_family_RA_melt[c(Heatmap_family_RA_melt$variable == "Alphaflexiviridae"),]
wilcox.test(value~Disease, data =Statistics_f,alternative = "two.sided",  paired = FALSE, exact = F)
effect_size <- wilcox_effsize(Statistics_f,value~Disease)
(effect_size$effsize)^2
effect_size$n1 # n = HE
effect_size$n2 # n = HS

Statistics_f <- Heatmap_family_RA_melt[c(Heatmap_family_RA_melt$variable == "Totiviridae"),]
wilcox.test(value~Disease, data =Statistics_f,alternative = "two.sided",  paired = FALSE, exact = F)
effect_size <- wilcox_effsize(Statistics_f,value~Disease)
(effect_size$effsize)^2
effect_size$n1 # n = HE
effect_size$n2 # n = HS

Statistics_f <- Heatmap_family_RA_melt[c(Heatmap_family_RA_melt$variable == "Virgaviridae"),]
wilcox.test(value~Disease, data =Statistics_f,alternative = "two.sided",  paired = FALSE, exact = F)
effect_size <- wilcox_effsize(Statistics_f,value~Disease)
(effect_size$effsize)^2
effect_size$n1 # n = HE
effect_size$n2 # n = HS

# Multiple testing correction 
P_values <-  c(0.04784,0.7809,0.9571,0.1829)
P_values
p.adjust(P_values, method = "BH", 4)
####################################

######### !!!!!!!!!!!!!!! #########