####################################
# SCRIPT 4: PROKARYOTIC VIRUSES
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
# Metadata
####################################
setwd("/Users/daan/Desktop/Bioinformatics/pre-processing/HS project/Metadata/Final")
dir()

metadata <- read_excel("HS_Metadata_Final2_R.xlsx")
metadata <- as.data.frame(metadata)
rownames(metadata) <- metadata$Patient_ID
####################################
# 1. Barplot: Final class 
####################################
table(Mastertable_family_like_approach$Final_viral)
Mastertable_phages <- Mastertable_family_like_approach[Mastertable_family_like_approach$Final_viral=="phage",]
unique(Mastertable_phages$Family_level)
sort(colSums(Mastertable_phages[1:127]))
names <- colnames(Mastertable_phages[1:127][,colSums(Mastertable_phages[1:127])!=0])
names # 99 samples of 127 samples or 77.95%
(99/127)*100

Mastertable_phages_categories <- aggregate(. ~Final_class, FUN = sum, data = Mastertable_phages[,colnames(Mastertable_phages) %in% names | colnames(Mastertable_phages) == 'Final_class'])
rownames(Mastertable_phages_categories) <- Mastertable_phages_categories$Final_class
Mastertable_phages_categories$Final_class <- NULL

Mastertable_phages_categories <- sweep(Mastertable_phages_categories, 2, colSums(Mastertable_phages_categories), '/')
Mastertable_phages_categories[is.na(Mastertable_phages_categories)] <- 0

Mastertable_phages_categories <- merge(t(Mastertable_phages_categories), metadata, by = 0, all = F)
rownames(Mastertable_phages_categories) <- Mastertable_phages_categories$Row.names
Mastertable_phages_categories$Row.names <- NULL

vector_1 <- (which(names(Mastertable_phages_categories)== "Patient_ID"))
vector_2 <- (which(names(Mastertable_phages_categories)== "Disease")) 
vector_3 <- (which(names(Mastertable_phages_categories)== "HS")) 
vector_4 <- (which(names(Mastertable_phages_categories)== "Hurley")) 

Mastertable_phages_categories_m <- Mastertable_phages_categories[,c(1:vector_1,vector_2,vector_3,vector_4)]
Mastertable_phages_categories_melt <- melt(Mastertable_phages_categories_m, id.vars = c("Patient_ID","Disease","HS","Hurley"))

ggplot(Mastertable_phages_categories_melt) +
 aes(x = variable, y = value, fill = Disease) +
 geom_boxplot() +
 scale_fill_manual(values = c(HE = "#609462", HS = "#b23434")) +
 theme_bw()

ggplot(Mastertable_phages_categories_melt) +
 aes(x = variable, y = value, fill = Disease) +
 geom_boxplot() +
 scale_fill_manual(values = c(HE = "#609462", HS = "#b23434")) +
 theme_bw() +
 facet_wrap(vars(Disease))

ggplot(Mastertable_phages_categories_melt) +
        aes(x = variable, y = value, fill = variable) +
        geom_col() +
        scale_fill_manual(values = c(Caudoviricetes = "#AB5049", Malgrandaviricetes = "#5F8E85", Unannotated = "#A23278")) +
        theme_classic()

esquisser(Mastertable_phages_categories_melt)

ggplot(Mastertable_phages_categories_melt) +
 aes(x = variable, fill = variable, weight = value) +
 geom_bar() +
  scale_fill_manual(values = c(Caudoviricetes = "#AB5049", Malgrandaviricetes = "#5F8E85", Unannotated = "#A23278")) +
theme_classic()

sum(Mastertable_phages_categories_melt$value[Mastertable_phages_categories_melt$variable=="Malgrandaviricetes"])
sum(Mastertable_phages_categories_melt$value[Mastertable_phages_categories_melt$variable=="Caudoviricetes"])
sum(Mastertable_phages_categories_melt$value[Mastertable_phages_categories_melt$variable=="Unannotated"])
6.06+53.98+38.96
####################################
# 2. Barplot: Final subfamily 
####################################
Mastertable_phages <- Mastertable_viral_family[Mastertable_viral_family$Final_viral=="phage",]
unique(Mastertable_phages$Family_level)
sort(colSums(Mastertable_phages[1:127]))
names <- colnames(Mastertable_phages[1:127][,colSums(Mastertable_phages[1:127])!=0])
names # 99 samples of 127 samples or 77.95%
(99/127)*100

Mastertable_phages_categories <- aggregate(. ~Final_class, FUN = sum, data = Mastertable_phages[,colnames(Mastertable_phages) %in% names | colnames(Mastertable_phages$) == 'Final_class'])
rownames(Mastertable_phages_categories) <- Mastertable_phages_categories$Final_class
Mastertable_phages_categories$Final_class <- NULL

Mastertable_phages_categories <- sweep(Mastertable_phages_categories, 2, colSums(Mastertable_phages_categories), '/')
Mastertable_phages_categories[is.na(Mastertable_phages_categories)] <- 0

Mastertable_phages_categories <- merge(t(Mastertable_phages_categories), metadata, by = 0, all = F)
rownames(Mastertable_phages_categories) <- Mastertable_phages_categories$Row.names
Mastertable_phages_categories$Row.names <- NULL

vector_1 <- (which(names(Mastertable_phages_categories)== "Patient_ID"))
vector_2 <- (which(names(Mastertable_phages_categories)== "Disease")) 
vector_3 <- (which(names(Mastertable_phages_categories)== "HS")) 
vector_4 <- (which(names(Mastertable_phages_categories)== "Hurley")) 

Mastertable_phages_categories_m <- Mastertable_phages_categories[,c(1:vector_1,vector_2,vector_3,vector_4)]
Mastertable_phages_categories_melt <- melt(Mastertable_phages_categories_m, id.vars = c("Patient_ID","Disease","HS","Hurley"))

ggplot(Mastertable_phages_categories_melt) +
        aes(x = variable, y = value, fill = Disease) +
        geom_boxplot() +
        scale_fill_manual(values = c(HE = "#609462", HS = "#b23434")) +
        theme_bw()

ggplot(Mastertable_phages_categories_melt) +
        aes(x = variable, y = value, fill = Disease) +
        geom_boxplot() +
        scale_fill_manual(values = c(HE = "#609462", HS = "#b23434")) +
        theme_bw() +
        facet_wrap(vars(Disease))

ggplot(Mastertable_phages_categories_melt) +
        aes(x = variable, y = value, fill = variable) +
        geom_col() +
        scale_fill_manual(values = c(Caudoviricetes = "#AB5049", Malgrandaviricetes = "#5F8E85", Unannotated = "#A23278"
        )) +
        theme_bw()

sum(Mastertable_phages_categories_melt$value[Mastertable_phages_categories_melt$variable=="Malgrandaviricetes"])
sum(Mastertable_phages_categories_melt$value[Mastertable_phages_categories_melt$variable=="Caudoviricetes"])
sum(Mastertable_phages_categories_melt$value[Mastertable_phages_categories_melt$variable=="Unannotated"])
6.06+53.98+38.96
####################################
# 3. Barplot: prokaryotic viruses - lifestyle
####################################
Mastertable_phages$`lysogenic cycle`[is.na(Mastertable_phages$`lysogenic cycle`)] <- "non-lysogenic"
Mastertable_phages$`lysogenic cycle`[Mastertable_phages$`lysogenic cycle` == "temperate"] <- "lysogenic"
table(Mastertable_phages$`lysogenic cycle`)

Mastertable_phages_categories <- aggregate(. ~`lysogenic cycle`, FUN = sum, data = Mastertable_phages[,colnames(Mastertable_phages) %in% names | colnames(Mastertable_phages) == 'lysogenic cycle'])
rownames(Mastertable_phages_categories) <- Mastertable_phages_categories$`lysogenic cycle`
Mastertable_phages_categories$`lysogenic cycle` <- NULL

Mastertable_phages_categories <- sweep(Mastertable_phages_categories, 2, colSums(Mastertable_phages_categories), '/')
Mastertable_phages_categories[is.na(Mastertable_phages_categories)] <- 0

Mastertable_phages_categories <- merge(t(Mastertable_phages_categories), metadata, by = 0, all = F)
rownames(Mastertable_phages_categories) <- Mastertable_phages_categories$Row.names
Mastertable_phages_categories$Row.names <- NULL

Mastertable_phages_categories_m <- melt(Mastertable_phages_categories, id.vars = c("Patient_ID","Disease","HS","Hurley"))

ggplot(Mastertable_phages_categories_m) +
 aes(x = variable, y = value, fill = Disease) +
 geom_boxplot() +
 scale_fill_hue(direction = 1) +
 theme_bw()
####################################
# 4. Median completeness
####################################
# In addition, to determine completeness is hard because it is based on current databases that often give ‘not-determined’ value to the phage contig identfiied with the viral cluster (VC) approach.
mean(Mastertable_viral_family$completeness[Mastertable_viral_family$Final_class=="Unannotated"])
mean(Mastertable_viral_family$completeness[Mastertable_viral_family$Final_class=="Caudoviricetes"])
mean(Mastertable_viral_family$completeness[Mastertable_viral_family$Final_class=="Malgrandaviricetes"])
####################################