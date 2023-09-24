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
# 1. Dataset
####################################
Mastertable # Reads & contigs figure
Mastertable_genus1 # Percentage classified on class, order, family level
Mastertable_eukaryotic # Eukaryotic viral figure
Mastertable_phage  # Phage analysis
####################################
# 2. Barplot
####################################
# 2.1 Class
####################################
Mastertable_phages <- Mastertable_phage
sort(colSums(Mastertable_phages[1:304]))
names <- colnames(Mastertable_phages[1:304][,colSums(Mastertable_phages[1:304])!=0])
names

Mastertable_phages_categories <- aggregate(. ~Final_class, FUN = sum, data = Mastertable_phages[,colnames(Mastertable_phages) %in% names | colnames(Mastertable_phages) == 'Final_class'])
rownames(Mastertable_phages_categories) <- Mastertable_phages_categories$Final_class
Mastertable_phages_categories$Final_class <- NULL
Mastertable_phages_categories$Total <- rowSums(Mastertable_phages_categories)
Mastertable_phages_categories$Percentage <- (Mastertable_phages_categories$Total/sum(Mastertable_phages_categories$Total))*100
ncol(Mastertable_phages_categories)
Mastertable_phages_categories <- Mastertable_phages_categories[303:304]
Mastertable_phages_categories$names <- rownames(Mastertable_phages_categories)

ggplot(Mastertable_phages_categories) +
 aes(x = reorder(names,-Percentage), weight = Percentage) +
 geom_bar(fill = "#112446") +
 theme_bw()

# Caudoviricetes: 52.5%
# Malgrandaviricetes: 43.2%
# Faserviricetes: 3.74%
# Unannotated 0.47%
# Duplopiviricetes: 0.07%
####################################
# 2.2 Order
####################################
Mastertable_phages_order <- aggregate(. ~Final_order, FUN = sum, data = Mastertable_phages[,colnames(Mastertable_phages) %in% names | colnames(Mastertable_phages) == 'Final_order'])
rownames(Mastertable_phages_order) <- Mastertable_phages_order$Final_order
Mastertable_phages_order$Final_order <- NULL
Mastertable_phages_order$Total <- rowSums(Mastertable_phages_order)
Mastertable_phages_order$Percentage <- (Mastertable_phages_order$Total/sum(Mastertable_phages_order$Total))*100
Mastertable_phages_order <- Mastertable_phages_order[303:304]
Mastertable_phages_order$names <- rownames(Mastertable_phages_order)
View(Mastertable_phages_order)

ggplot(Mastertable_phages_order) +
  aes(x = reorder(names,-Percentage), weight = Percentage) +
  geom_bar(fill = "#112446") +
  theme_bw()

# Unclassified Caudoviricetes: 43.7%
# Petitvirales: 43.2%
# Crassvirales: 8.83%
# Tubulavirales: 3.74%
# Unannotated 0.48%
# Durnavirales: 0.07%
####################################
# 2.3 Family
####################################
Mastertable_phages_family <- aggregate(. ~Final_family, FUN = sum, data = Mastertable_phages[,colnames(Mastertable_phages) %in% names | colnames(Mastertable_phages) == 'Final_family'])
rownames(Mastertable_phages_family) <- Mastertable_phages_family$Final_family
Mastertable_phages_family$Final_family <- NULL
Mastertable_phages_family$Total <- rowSums(Mastertable_phages_family)
Mastertable_phages_family$Percentage <- (Mastertable_phages_family$Total/sum(Mastertable_phages_family$Total))*100
Mastertable_phages_family <- Mastertable_phages_family[303:304]
Mastertable_phages_family$names <- rownames(Mastertable_phages_family)

ggplot(Mastertable_phages_family) +
  aes(x = reorder(names,-Percentage), weight = Percentage) +
  geom_bar(fill = "#112446") +
  theme_bw()
####################################
# 2.4 Genus
####################################
Mastertable_phages_genus <- aggregate(. ~Final_genus, FUN = sum, data = Mastertable_phages[,colnames(Mastertable_phages) %in% names | colnames(Mastertable_phages) == 'Final_genus'])
rownames(Mastertable_phages_genus) <- Mastertable_phages_genus$Final_genus
Mastertable_phages_genus$Final_genus <- NULL
Mastertable_phages_genus$Total <- rowSums(Mastertable_phages_genus)
Mastertable_phages_genus$Percentage <- (Mastertable_phages_genus$Total/sum(Mastertable_phages_genus$Total))*100
Mastertable_phages_genus <- Mastertable_phages_genus[303:304]
Mastertable_phages_genus$names <- rownames(Mastertable_phages_genus)

ggplot(Mastertable_phages_genus) +
  aes(x = reorder(names,-Percentage), weight = Percentage) +
  geom_bar(fill = "#112446") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
####################################