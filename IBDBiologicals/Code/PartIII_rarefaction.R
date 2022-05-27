####################################
# SCRIPT 3: CREATE RAREFIED TABLES AND ADJUST FOR COMPLETENESS
####################################
# Before starting this script
####################################
# ------------------> Load the "Global environment" output of script 1
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
#install.packages("modeldata")
#install.packages("tidyverse")
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
####################################
# 1. Rarefaction
####################################
# In the previous script "Reads_and_contigs" we determine the threshold to use for viral rarefaction. This is as given & that number we will use.
####################################
# 1.1 Viral rarefaction
####################################
# 1.1.1 Create a viral mastertable
####################################
#Mastertable_viral_unrarefied <- Mastertable
#Mastertable_viral_unrarefied <- Mastertable_viral_unrarefied[Mastertable_viral_unrarefied$Final_viral == "eukaryotic virus" | Mastertable_viral_unrarefied$Final_viral == "phage",]
Mastertable_unrarefied <- Mastertable
####################################
# 1.1.2 Create vector for samples name
####################################
# Check if "Final_coverage" is the last colnames, otherwise replace second line here to last column name
vector_1 <- (which(names(Mastertable)== "Virsorter"))
vector_1_1 <- (which(names(Mastertable)== "Virsorter")-1)
vector_2 <- (which(names(Mastertable)== "Final_class2"))
vector_ab_viral <- colnames(Mastertable[,1:vector_1_1])
vector_tax_viral <- colnames(Mastertable[,vector_1:vector_2])
####################################
# 1.1.3 Subset data to viral data if need
####################################
Mastertable_unrarefied_ab <- cbind(Mastertable[,vector_ab_viral])
Mastertable_unrarefied_tax <- cbind(Mastertable[,vector_tax_viral])
Mastertable_unrarefied_ab[is.na(Mastertable_unrarefied_ab)] <- 0
####################################
# 1.1.4 Discard samples based on threshold
####################################
colSums(Mastertable_unrarefied_ab)
min(colSums(Mastertable_unrarefied_ab))
head(Mastertable_unrarefied_ab[,colSums(Mastertable_unrarefied_ab) < 987944])

# Adopt the rarefaction number here to your preference
# How many and which samples will be discarded be rarefying to the number of reads?
# ------------------> 12 samples were already not included in this mastertable
# ------------------> 43 additional samples will be discard by doing this

Mastertable_unrarefied_ab <- Mastertable_unrarefied_ab[,colSums(Mastertable_unrarefied_ab) >= 979620]
ncol(Mastertable_unrarefied_ab)
View(Mastertable_unrarefied_ab)
# How many samples are left?
# ------------------> 377 samples are left (432-12-43)
# ------------------> SEE excel table: in Analysis/Biologicals/Reads_and_contigs/Reads

# Samples Discarded?
# B2, B225, B230, B75, B178, B25, B60, B43, B38, B193, B370, B67
# B113, B119, B122, B134, B138, B140, B141, B168, B170, B182
# B188, B192, B206, B210, B219, B220, B238, B24, B244, B246
# B247, B29, B30, B300, B321, B322, B327, B330, B338, B349, B363,
# B374, B378,B381, B418, B419, B426, B458, B53, B76, B79, B83, B94

# Check other excel file (/Users/daan/Desktop/Bioinformatics/Analysis/Biologicals/Reads_contigs/Reads/Rarefaction)
# How many are still paired?
  
# ------------------> total patients: 56 UC; 154 CD; 210 patients
# ------------------> 39 paired UC; 17 unpaired UC
# ------------------> 128 paired CD; 26 unpaired CD
####################################
# 1.1.5 Actual rarefaction
####################################
min(colSums(Mastertable_unrarefied_ab))
# The number to which we want to rarefy should be given as an output here. What is this number?
# ------------------> 979620 PE reads

Mastertable_rarefied_ab <- t(rrarefy(t(Mastertable_unrarefied_ab), sample = 987944))
Mastertable_rarefied_ab <- as.data.frame(Mastertable_rarefied_ab)
####################################
# 1.1.6 Evaluate rarefaction
####################################
sum(colSums(Mastertable_unrarefied_ab)) 
sum(colSums(Mastertable_rarefied_ab)) 
# ------------------> 4165524116 PE unrarefied reads
# ------------------> 372454888 PE rarefied reads (~ or 8,96%)
####################################
# 1.1.7 Combine rarefied abundance tables with taxonomy table
####################################
Mastertable_rarefied <- transform(merge(Mastertable_rarefied_ab,Mastertable_unrarefied_tax, by = 0), row.names=Row.names, Row.names = NULL)
Mastertable_rarefied$Totalnumberofreads <- NULL

vector_1 <- (which(names(Mastertable_rarefied)== "Virsorter"))
vector_1_1 <- (which(names(Mastertable_rarefied)== "Virsorter")-1)

Mastertable_rarefied$Totalnumberofreads <- as.numeric(rowSums(Mastertable_rarefied[,1:vector_1_1]))
####################################
# 1.1.8 Discard rows with 0 reads & calculate differences
####################################
# Total 
sum(Mastertable$Totalnumberofreads) # 4127511423
sum(Mastertable$Totalnumberofreads[Mastertable$Final_superkingdom == "Viruses"]) # 2670560838

# Viruses
Mastertable_rarefied <- Mastertable_rarefied[!(Mastertable_rarefied$Totalnumberofreads == "0"),]
Mastertable_viral_rarefied <- Mastertable_rarefied[Mastertable_rarefied$Final_superkingdom == "Viruses",]
Mastertable_viral_unrarefied <- Mastertable_unrarefied[Mastertable_unrarefied$Final_superkingdom == "Viruses",]

# ____________RAREFIED_VS_UNRAREFIED_VIRAL_DATA_______________________________________________________________________________________________________________________________________

# Evaluation RAREFIED viral reads
sum(Mastertable_viral_rarefied$Totalnumberofreads)
sum(Mastertable_viral_unrarefied$Totalnumberofreads)
(sum(Mastertable_rarefied$Totalnumberofreads)/sum(Mastertable_viral_unrarefied$Totalnumberofreads))*100
# unrarefied viral reads ------------------------> 2693456553 PE reads
# rarefied viral reads --------------------------> 237655536 PE reads
# % rarefied viral reads left  ------------------> 8.94% viral reads
# % rarefied viral reads discarded  -------------> 91.06% viral reads

# Evaluation of RAREFIED viral contigs
nrow(Mastertable_viral_unrarefied)
nrow(Mastertable_viral_rarefied)
(nrow(Mastertable_viral_rarefied)/nrow(Mastertable_viral_rarefied))*100
# unrarefied viral contigs ------------------------> 3058 NR-contigs
# rarefied viral contigs --------------------------> 3033 NR-contigs
# % rarefied viral contigs left  ------------------> 99.18% viral contigs left
# % rarefied viral contigs discard  ---------------> 0.82% viral contigs discarded

# ____________RAREFIED_DATA_VIRAL_PORTION_________________________________________________________________________________________________________________________________

# Evaluation reads of rarefied data
sum(Mastertable_rarefied$Totalnumberofreads)
sum(Mastertable_viral_rarefied$Totalnumberofreads)
(sum(Mastertable_viral_rarefied$Totalnumberofreads)/sum(Mastertable_rarefied$Totalnumberofreads))*100

# rarefied mapped reads -------------------------> 372454888 PE reads
# rarefied viral reads --------------------------> 237655536 PE reads
# % rarefied viral reads  -----------------------> 63.81% viral reads  
# % rarefied reads-------------------------------> 36.19% reads 

# Evaluation contigs of rarefied data
nrow(Mastertable_rarefied)
nrow(Mastertable_viral_rarefied)
(nrow(Mastertable_viral_rarefied)/nrow(Mastertable_rarefied))*100

# rarefied mapped contigs ------------------------> 109.452 NR-contigs
# rarefied viral contigs -------------------------> 3033 NR-contigs 
# % rarefied viral contigs -----------------------> 2.77% NR-viral-contigs 
# % rarefied contigs  ----------------------------> 97.23% NR-contigs

## Conclusion: We use rarefaction threshold of 987.944 PE reads. We rarefied to this threshold.
## Total Nr. of PE reads that were left were 372.454.888 PE reads mapping to 109.452 NR-contigs.
## Of the total, the majority 237.655.536 PE reads (~ 64,01%) were classified as viral.
## Of the total, the minority 3.033 NR-contigs (~ 3,46%) were classified as viral.
## You can make a nice donut plot out of this (rarefied reads & rarefied NR-contigs)

# ____________RAREFIED_DATA_VIRAL_CATEGORY_________________________________________________________________________________________________________________________________
sum(Mastertable_viral_rarefied$Totalnumberofreads[Mastertable_viral_rarefied$Final_viral == "eukaryotic virus"])
sum(Mastertable_viral_rarefied$Totalnumberofreads[Mastertable_viral_rarefied$Final_viral == "phage"])
(580802/(237074734+580802))*100
100-(580802/(237074734+580802))*100

# eukaryotic viral reads ------------------> 580.802 PE reads 
# phage reads -----------------------------> 237.074.734 PE reads
# % eukaryotic viral reads ----------------> 0,24% PE reads 
# % prokaryotic viral reads ---------------> 99.76 PE reads 

nrow(Mastertable_viral_rarefied[Mastertable_viral_rarefied$Final_viral == "eukaryotic virus",])
nrow(Mastertable_viral_rarefied[Mastertable_viral_rarefied$Final_viral == "phage",])
(61/(2872+61))*100

# eukaryotic viral contigs ----------------> 61 NR-contigs
# phage contigs ---------------------------> 2972 NR-contugs
# % eukaryotic viral contigs ----------------> 2.08 % NR-contigs 
# % prokaryotic viral contigs ---------------> 97.02 % NR-contigs 
####################################
# 2. Figures: pie-charts (~ perhaps also try donut plots here)
####################################
# 2.1. Pre-rarefaction: global mastertable
####################################
# 2.1.1 read level: total
####################################
# View(Mastertable)
# Because the kingdoms "Archaea, Unclassified, dark matter & Eukaryota" are very low. We can merge these and name them "other"
sum(Mastertable_rarefied$Totalnumberofreads)
nrow(Mastertable_rarefied)
Mastertable_global_rarefied_categories <- Mastertable_rarefied
table(Mastertable_global_rarefied_categories$Final_superkingdom)

Mastertable_global_rarefied_categories$Final_superkingdom[Mastertable_global_rarefied_categories$Final_superkingdom == "Archaea" | Mastertable_global_rarefied_categories$Final_superkingdom == "Dark matter" | Mastertable_global_rarefied_categories$Final_superkingdom == "Eukaryota" ] <- "Other"
unique(Mastertable_global_rarefied_categories$Final_superkingdom)

# View(Mastertable)
Mastertable_global_rarefied_categories <- aggregate(. ~Final_superkingdom, FUN = sum, data = Mastertable_global_rarefied_categories[,colnames(Mastertable_global_rarefied_categories) %in% names | colnames(Mastertable_global_rarefied_categories) == 'Final_superkingdom'])
rownames(Mastertable_global_rarefied_categories) <- Mastertable_global_rarefied_categories$Final_superkingdom
Mastertable_global_rarefied_categories$Final_superkingdom <- NULL

Mastertable_global_rarefied_categories$absolute <- rowSums(Mastertable_global_rarefied_categories)
Mastertable_global_rarefied_categories$relative <- as.numeric(((Mastertable_global_rarefied_categories$absolute/sum(Mastertable_global_rarefied_categories$absolute))*100))
Mastertable_global_rarefied_categories$relative <- round(Mastertable_global_rarefied_categories$relative, digits = 1)
Mastertable_global_rarefied_categories$names <- rownames(Mastertable_global_rarefied_categories)
table(Mastertable_global_rarefied_categories$relative)
#View(Mastertable_global_rarefied_categories)

## This you'll have to adopt quickly & do manually for some part
Mastertable_global_rarefied_categories <- Mastertable_global_rarefied_categories %>%
  select(absolute:names)

#Mastertable_global_rarefied_categories$names <- factor(Mastertable_global_rarefied_categories$names, levels = rev(c("Unannotated", "Other", "Bacteria", "Viruses")))
Mastertable_global_rarefied_categories$names <- factor(Mastertable_global_rarefied_categories$names, levels = rev(c("eukaryotic virus", "no eukaryotic virus")))
Mastertable_global_rarefied_categories$relative_adapted <- Mastertable_global_rarefied_categories$relative
Mastertable_global_rarefied_categories$relative_adapted <- c("32.8","67.2")

View(Mastertable_global_rarefied_categories)
100-32.8
#Mastertable_global_rarefied_categories$relative_adapted <- c("31.2","3.8", "1.0", "64.0")
Mastertable_global_rarefied_categories$relative_adapted <- as.numeric(Mastertable_global_rarefied_categories$relative_adapted)
view(Mastertable_global_rarefied_categories)

# (ii) Add label position (~ perhaps usefull)
Mastertable_global_rarefied_categories <- Mastertable_global_rarefied_categories %>%
  arrange(desc(names)) %>%
  mutate(lab.ypos = cumsum(relative) - 0.5*relative)

# (iii) Pie-chart (everything < 1% is not shown)
mycols <- c("#e41a1c", "#377eb8", "#4daf4a","#984ea3")
ggplot(Mastertable_global_rarefied_categories, aes(x = "", relative, fill = names)) +
  geom_bar(width = 1, stat = "identity", color = "black") +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = mycols) +
  geom_label(aes(label = scales::percent(relative_adapted, scale = 1, accuracy = 0.1), x = 1.1),
             position = position_stack(vjust = 0.2),
             color = "black",
             fill="#c6dbef") +
  labs(fill = "Category",x = NULL,y = NULL, title = "") +
  theme_void()

# (D) Pie-chart: viral category
Mastertable_global_rarefied_categories$relative <- NULL
Mastertable_global_rarefied_categories$relative_adapted <- NULL
Mastertable_global_rarefied_categories$lab.ypos <- NULL
Mastertable_global_rarefied_categories$absolute <- c("235800674", "585930", "","")
Mastertable_global_rarefied_categories$names <- c("phage", "eukaryotic virus", "","")
Mastertable_global_rarefied_categories <- Mastertable_global_rarefied_categories[1:2,]

Mastertable_global_rarefied_categories$absolute <- as.numeric(Mastertable_global_rarefied_categories$absolute)
Mastertable_global_rarefied_categories$relative <- (Mastertable_global_rarefied_categories$absolute/sum(Mastertable_global_rarefied_categories$absolute))*100
Mastertable_global_rarefied_categories$relative <- round(Mastertable_global_rarefied_categories$relative, digits = 1)

Mastertable_global_rarefied_categories <- Mastertable_global_rarefied_categories %>%
  mutate(lab.ypos = cumsum(relative) - 0.5*relative)
#view(Mastertable_global_rarefied_categories)

mycols <- c("yellow","darkorange")
ggplot(Mastertable_global_rarefied_categories, aes(x = "", y = relative, fill = names)) +
  geom_bar(width = 1, stat = "identity", color = "black") +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = mycols) +
  theme_void() +
  labs(fill = "Viruses",x = NULL,y = NULL, title = "") +
  geom_label(aes(y = lab.ypos, label = scales::percent(relative, scale = 1, accuracy = 0.1), x = 1.1), color ="black", fill="#c6dbef")
####################################
# 2.1.2 contig level: total
####################################
## Still have to rewrite this script, wrong names

# (i) determine values
Mastertable_global_rarefied_categories <- Mastertable_rarefied
Mastertable_global_rarefied_categories$Final_superkingdom[Mastertable_global_rarefied_categories$Final_superkingdom == "Archaea" | Mastertable_global_rarefied_categories$Final_superkingdom == "Dark matter" | Mastertable_global_rarefied_categories$Final_superkingdom == "Eukaryota" ] <- "Other"
table(Mastertable_global_rarefied_categories$Final_superkingdom)
contig_total <- as.data.frame(unique(Mastertable_global_rarefied_categories$Final_superkingdom))
#View(contig_total)
#View(Mastertable_global_rarefied_categories)

print(NROW(Mastertable_global_rarefied_categories[Mastertable_global_rarefied_categories$Final_superkingdom == "Viruses",])) # 3033
print(NROW(Mastertable_global_rarefied_categories[Mastertable_global_rarefied_categories$Final_superkingdom == "Bacteria",])) # 101765
print(NROW(Mastertable_global_rarefied_categories[Mastertable_global_rarefied_categories$Final_superkingdom == "Other",])) # 3301
print(NROW(Mastertable_global_rarefied_categories[Mastertable_global_rarefied_categories$Final_superkingdom == "Unannotated",])) # 1353

contig_total$absolute <- c("3033","101765","3301","1353")
contig_total$absolute <- as.numeric(contig_total$absolute)
colnames(contig_total) <- c("names", "absolute")
contig_total$relative <- as.numeric(contig_total$absolute/sum(contig_total$absolute)*100)
contig_total$relative <- round(contig_total$relative, digits = 1)
View(contig_total)

# contig_total$relative_adapted <- c("12.0", "84.7", "3.1", "")
contig_total$names <- factor(contig_total$names, levels = rev(c("Bacteria", "Viruses","Other","Unannotated")))

# (ii) Add label position (~ perhaps usefull)
contig_total <- contig_total %>%
  arrange(desc(absolute)) %>%
  mutate(lab.ypos = cumsum(relative) - 0.5*relative)

# (iii) Pie-chart (everything < 1% is not shown)
mycols <- c("#984ea3","#4daf4a","#e41a1c","#377eb8")
ggplot(contig_total, aes(x = "", relative, fill = names)) +
  geom_bar(width = 1, stat = "identity", color = "black") +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = mycols) +
  geom_label(aes(label = scales::percent(relative, scale = 1, accuracy = 0.1), x = 1.3),
             position = position_stack(vjust = 0.5),
             color = "black",
             fill="#c6dbef") +
  labs(fill = "Category",x = NULL,y = NULL, title = "") +
  theme_void()

# (B) Viral categories
print(NROW(Mastertable_global_rarefied_categories[Mastertable_global_rarefied_categories$Final_viral == "phage",])) # 3523
print(NROW(Mastertable_global_rarefied_categories[Mastertable_global_rarefied_categories$Final_viral == "eukaryotic virus",])) # 61
view(contig_total)

contig_total$relative <- NULL
contig_total$relative_adapted <- NULL
contig_total$lab.ypos <- NULL
contig_total <- contig_total[1:2,]
contig_total$names <- c("phage", "eukaryotic virus")
contig_total$absolute <- c("2972", "61")

contig_total$absolute <- as.numeric(contig_total$absolute)
contig_total <- na.omit(contig_total)
contig_total$relative <- (contig_total$absolute/sum(contig_total$absolute))*100
contig_total$relative <- round(contig_total$relative, digits = 2)
View(contig_total)

contig_total <- contig_total %>%
  mutate(lab.ypos = cumsum(relative) - 0.5*relative)

mycols <- c("yellow","orange")
ggplot(contig_total, aes(x = "", y = relative, fill = names)) +
  geom_bar(width = 1, stat = "identity", color = "black") +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = mycols) +
  theme_void() +
  labs(fill = "viral category",x = NULL,y = NULL, title = "") +
  geom_label(aes(y = lab.ypos, label = scales::percent(relative, scale = 1, accuracy = 0.1), x = 1.1), color ="black", fill="#c6dbef")
####################################
# 2.2 Post-rarefaction: viral mastertable
####################################
# 2.2.1 read level: total
####################################
## still have to we-write this entire piece of code in "2.2.1", vector already made before see how I do it
view(Mastertable_viral_rarefied)

Mastertable_viral_rarefied_categories <- Mastertable_viral_rarefied

Mastertable_viral_rarefied_categories_vector <- colnames(Mastertable_viral_rarefied_categories %>% select(starts_with("TIPS")))

# (A) aggregate data
Mastertable_viral_rarefied_categories <- aggregate(. ~FinalClassificationPlots, FUN = sum, data = Mastertable_viral_rarefied_categories[,colnames(Mastertable_viral_rarefied_categories) %in% Mastertable_viral_rarefied_categories_vector | colnames(Mastertable_viral_rarefied_categories) == 'FinalClassificationPlots'])
rownames(Mastertable_viral_rarefied_categories) <- Mastertable_viral_rarefied_categories$FinalClassificationPlots
Mastertable_viral_rarefied_categories$FinalClassificationPlots <- NULL

# (B) Create an absolute and relative column
Mastertable_viral_rarefied_categories$absolute <- rowSums(Mastertable_viral_rarefied_categories)
Mastertable_viral_rarefied_categories$relative <- as.numeric(((Mastertable_viral_rarefied_categories$absolute/sum(Mastertable_viral_rarefied_categories$absolute))*100))
Mastertable_viral_rarefied_categories$relative <- round(Mastertable_viral_rarefied_categories$relative, digits = 0.01)
Mastertable_viral_rarefied_categories$names <- rownames(Mastertable_viral_rarefied_categories)
Mastertable_viral_rarefied_categories$relative_adapted <- c("2","98")
Mastertable_global_rarefied_categories$relative_adapted <- as.numeric(Mastertable_global_rarefied_categories$relative_adapted)

# (C) Pie-chart: categories on read level
# (i) subset cumulative data and Add label position

Mastertable_viral_rarefied_categories <- Mastertable_viral_rarefied_categories %>%
  select(absolute:relative_adapted)

Mastertable_viral_rarefied_categories <- Mastertable_viral_rarefied_categories %>%
  arrange(desc(names)) %>%
  mutate(lab.ypos = cumsum(relative) - 0.5*relative)

#view(Mastertable_viral_rarefied_categories)

# (ii) figure
mycols <- c("#6baed6","#084594")
ggplot(Mastertable_viral_rarefied_categories, aes(x = "", y = relative, fill = names)) +
  geom_bar(width = 1, stat = "identity", color = "black") +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = mycols) +
  theme_void() +
  labs(fill = "viral category",x = NULL,y = NULL, title = "") +
  geom_label(aes(y = lab.ypos, label = scales::percent(relative, scale = 1, accuracy = 0.1), x = 1), color ="black", fill="#c6dbef")
####################################
# 2.2.2 contig level: total
###################################
## still have to we-write this entire piece of code in "2.2.1", vector already made before see how I do it

Mastertable_viral_rarefied_categories <- Mastertable_viral_rarefied
contig_total <- as.data.frame(unique(Mastertable_viral_rarefied_categories$FinalClassificationPlots))

unique(Mastertable_viral_rarefied_categories$FinalClassificationPlots)
print(NROW(Mastertable_viral_rarefied_categories[Mastertable_viral_rarefied_categories$FinalClassificationPlots == "prokaryotic viruses",])) # 4455
print(NROW(Mastertable_viral_rarefied_categories[Mastertable_viral_rarefied_categories$FinalClassificationPlots == "eukaryotic viruses",])) # 73

contig_total$absolute <- c("4455","73")
contig_total$absolute <- as.numeric(contig_total$absolute)
colnames(contig_total) <- c("names", "absolute")
contig_total$relative <- (contig_total$absolute/sum(contig_total$absolute)*100)
contig_total$relative <- round(contig_total$relative, digits = 1)
contig_total$relative_adapted <- c("98.4", "1.6")
contig_total$relative_adapted <- as.numeric(contig_total$relative_adapted)
contig_total$names <- factor(contig_total$names, levels = rev(c("prokaryotic viruses","eukaryotic viruses")))


contig_total <- contig_total %>%
  mutate(lab.ypos = cumsum(relative) - 0.5*relative)

mycols <- c("#6baed6","#084594")
ggplot(contig_total, aes(x = "", y = relative, fill = names)) +
  geom_bar(width = 1, stat = "identity", color = "black") +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = mycols) +
  theme_void() +
  labs(fill = "viral category",x = NULL,y = NULL, title = "") +
  geom_label(aes(y = lab.ypos, label = scales::percent(relative, scale = 1, accuracy = 0.1), x = 1.1), color ="black", fill="#c6dbef")
####################################
# 3. Completeness correction: obtain medium- to-high-quality genomes
####################################
# CheckV allows us to put a threshold for the completeness of the identified viral genomes.
# We choose for the sake of diversity analysis to include a completeness cutoff of 50% or more
# Also look to the error rate ofcourse, a low error rate meaning a thrustworthy completeness prediction.

# The question remains what we do with undetermined contigs? here the authors say that because (i) it is a novel virus, (ii) it is no virus, and (iii) it is too short you can best use a large vast length cutoff such as > 30 kb to still include it in your diversity analysis.
####################################
# 3.1 Viral unrarefied
####################################
Mastertable_viral_unrarefied_completeness <- Mastertable_viral_unrarefied
Mastertable_viral_unrarefied_completeness[, Mastertable_viral_unrarefied_completeness$Final_viral == "phage" & Mastertable_viral_unrarefied_completeness$completeness_quality == "Not-determined" & Mastertable_viral_unrarefied_completeness$length < 30000] <- NULL
table(Mastertable_viral_unrarefied$completeness_quality)
table(Mastertable_viral_unrarefied_completeness$completeness_quality)

# Did the 30kb cutoff onder "Not-determined" group work? What is the number of viral contigs into the different quality tiers?
# Complete ------------------------------>
# High-quality -------------------------->
# Medium-quality  ----------------------->
# Low-quality  -------------------------->
# Not-determined (> 30 kb) --------------->

Mastertable_viral_unrarefied_completeness <- Mastertable_viral_unrarefied_completeness[Mastertable_viral_unrarefied_completeness$Final_viral == "phage" & Mastertable_viral_unrarefied_completeness$completeness_quality == "Complete" | Mastertable_viral_unrarefied_completeness$Final_viral == "phage" & Mastertable_viral_unrarefied_completeness$completeness_quality == "High-quality" | Mastertable_viral_unrarefied_completeness$Final_viral == "phage" & Mastertable_viral_unrarefied_completeness$completeness_quality == "Medium_quality" | Mastertable_viral_unrarefied_completeness$Final_viral == "phage" & Mastertable_viral_unrarefied_completeness$completeness_quality == "Not-determined",]
table(Mastertable_viral_unrarefied_completeness$Final_order)
table(Mastertable_viral_unrarefied_completeness$Final_family)
nrow(Mastertable_viral_unrarefied_completeness)
(nrow(Mastertable_viral_unrarefied_completeness)/nrow(Mastertable_viral_unrarefied))*100

# How many different order do you find with high-to medium-quality contigs?
# -------------------------->

# How many different families do you find with high-to medium-quality contigs?
# -------------------------->

# How many different contigs do you find with high-to medium-quality contigs?
# -------------------------->

# What is percentages of contigs left for diversity analysis?
# -------------------------->
####################################
# 3.2 Viral rarefied: completeness value AND table creation
####################################
#View(Mastertable_viral_rarefied)
nrow(Mastertable_viral_rarefied) # 3033 NR-contigs
table(Mastertable_viral_rarefied$completeness_quality)

# Distribution quality-tiers?
# --------------------------> Complete: 645 (21,3%)
# --------------------------> High-quality: 1101 (36,3)   FC8D59
# --------------------------> Low-quality: 15 (0,5)     FEE08B
# --------------------------> Medium-quality: 1184 (39,1) E6F598
# --------------------------> Not determined: 87 (2,87) 99D594
# --------------------------> Total: 3153 NR-contigs

contig_comlpl <- as.data.frame(c("36.3","39.1","21.3","2.87", "0.5"))
contig_comlpl$relative <- contig_comlpl$`c("36.3", "39.1", "21.3", "2.87", "0.5")`
contig_comlpl$relative <- as.numeric(contig_comlpl$relative)
contig_comlpl$`c("36.3", "39.1", "21.3", "2.87", "0.5")` <- NULL
contig_comlpl$names <- c("not determined", "medium-quality","complete","high-quality", "low-quality")
#contig_comlpl$names <- factor(contig_comlpl$names, levels = rev(c("complete","high-quality", "low-quality", "medium-quality", "not determined")))

contig_comlpl <- contig_comlpl %>%
  mutate(lab.ypos = cumsum(relative) - 0.5*relative)

mycols <- c("#d53e4f","#fc8d59","#fee08b","#e6f598","#99d594")
ggplot(contig_comlpl, aes(x = "", y = relative, fill = names)) +
  geom_bar(width = 1, stat = "identity", color = "black") +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = mycols) +
  theme_void() +
  labs(fill = "quality tiers",x = NULL,y = NULL, title = "") +
  geom_label(aes(y = lab.ypos, label = scales::percent(relative, scale = 1, accuracy = 0.1), x = 1.1), color ="black", fill="#c6dbef")

#____________________________________________________________________________________________________________________________

Mastertable_viral_rarefied_completeness_est$lysogenic.cycle[is.na(Mastertable_viral_rarefied_completeness_est$lysogenic.cycle)] <- 0
Mastertable_viral_rarefied_completeness_est$Toxin.presence[is.na(Mastertable_viral_rarefied_completeness_est$Toxin.presence)] <- 0
Mastertable_viral_rarefied_completeness_est$Final_class2[!Mastertable_viral_rarefied_completeness_est$Final_class == "Caudoviricetes (CrAss)" & !Mastertable_viral_rarefied_completeness_est$Final_class == "Caudoviricetes (non-CrAss)"] <- 0
range(Mastertable_viral_rarefied_completeness_est$length[Mastertable_viral_rarefied_completeness_est$Final_class2 == "Caudoviricetes"])

## for average completess, subset without Not-determined fraction in calculations..
table(Mastertable_viral_rarefied_completeness_est$Final_viral)
(Mastertable_viral_rarefied_completeness_est[Mastertable_viral_rarefied_completeness_est$Final_viral == "eukaryotic virus",])

Mastertable_viral_rarefied$completeness[Mastertable_viral_rarefied$completeness_quality == "Not-determined"] <- 0
table(Mastertable_viral_rarefied$completeness == 0)
Mastertable_viral_rarefied_completeness_est <- Mastertable_viral_rarefied
sum(Mastertable_viral_rarefied_completeness_est$Totalnumberofreads) # 237655536
nrow(Mastertable_viral_rarefied_completeness_est) # 3033 total viruses
mean(Mastertable_viral_rarefied_completeness_est$length)
range(Mastertable_viral_rarefied_completeness_est$length)
table(Mastertable_viral_rarefied_completeness_est$lysogenic.cycle)
table(Mastertable_viral_rarefied_completeness_est$Toxin.presence)

nrow(Mastertable_viral_rarefied_completeness_est[Mastertable_viral_rarefied_completeness_est$Final_class == "Unannotated",]) # 457 Unannotated; 2696 annotated
table(Mastertable_viral_rarefied_completeness_est$completeness)
sum(Mastertable_viral_rarefied_completeness_est$completeness)/3033
# --------------------------> average completeness: 83.4%; Medium-quality-tier
# --------------------------> (219/3033)*100 = 7,22% of viruses annotated
(219/3033)*100
#____________________________________________________________________________________________________________________________

# Average completeness for phages?
table(Mastertable_viral_rarefied_completeness_est$Final_viral)
table(Mastertable_viral_rarefied_completeness_est$Final_class[Mastertable_viral_rarefied_completeness_est$Final_viral == "phage"])
sum(Mastertable_viral_rarefied_completeness_est$Totalnumberofreads[Mastertable_viral_rarefied_completeness_est$Final_viral == "phage"]) # 235,779,256 phage reads
nrow(Mastertable_viral_rarefied_completeness_est[Mastertable_viral_rarefied_completeness_est$Final_viral == "phage",]) # 3102
sum(Mastertable_viral_rarefied_completeness_est$completeness[Mastertable_viral_rarefied_completeness_est$Final_viral == "phage"])/3033
median(Mastertable_viral_rarefied_completeness_est$length[Mastertable_viral_rarefied_completeness_est$Final_viral == "phage"]) # 
range(Mastertable_viral_rarefied_completeness_est$length[Mastertable_viral_rarefied_completeness_est$Final_viral == "phage"]) # 
table(Mastertable_viral_rarefied_completeness_est$lysogenic.cycle[Mastertable_viral_rarefied_completeness_est$Final_viral == "phage"]) # 1364 temperate phages
table(Mastertable_viral_rarefied_completeness_est$Toxin.presence[Mastertable_viral_rarefied_completeness_est$Final_viral == "phage"]) #

# 3108 contigs; 235777278 reads; 1420 temperate; 24,5kb, 2.071 278.848; 1420 temperate; 220 toxins
# --------------------------> 80,2%; medium-quality-tier

# Average completeness for eukaryotic viruses?
table(Mastertable_viral_rarefied_completeness_est$Final_viral)
sum(Mastertable_viral_rarefied_completeness_est$Totalnumberofreads[Mastertable_viral_rarefied_completeness_est$Final_viral == "eukaryotic virus"]) # 585498 reads
nrow(Mastertable_viral_rarefied_completeness_est[Mastertable_viral_rarefied_completeness_est$Final_viral == "eukaryotic virus",]) # 45
median(Mastertable_viral_rarefied_completeness_est$length[Mastertable_viral_rarefied_completeness_est$Final_viral == "eukaryotic virus"]) # 3.4
sum(Mastertable_viral_rarefied_completeness_est$completeness[Mastertable_viral_rarefied_completeness_est$Final_viral == "eukaryotic virus"])/61
range(Mastertable_viral_rarefied_completeness_est$length[Mastertable_viral_rarefied_completeness_est$Final_viral == "eukaryotic virus"]) # 
table(Mastertable_viral_rarefied_completeness_est$lysogenic.cycle[Mastertable_viral_rarefied_completeness_est$Final_viral == "eukaryotic virus"]) # 1364 temperate phages
table(Mastertable_viral_rarefied_completeness_est$Toxin.presence[Mastertable_viral_rarefied_completeness_est$Final_viral == "eukaryotic virus"]) #
# 45 contigs; 585498 reads; 0 temperate; 0 toxin; 3.4kb, 2.1-8.3; 73.3%complete
# --------------------------> 80,2%; medium-quality-tier

#____________________________________________________________________________________________________________________________
# Phage Classes:
# prokaryotic viral families reads: 
table(Mastertable_viral_rarefied_completeness_est$Final_superkingdom)
table(Mastertable_viral_rarefied_completeness_est$Final_class[Mastertable_viral_rarefied_completeness_est$Final_viral == "phage"])
sum(Mastertable_viral_rarefied_completeness_est$Totalnumberofreads[Mastertable_viral_rarefied_completeness_est$Final_viral == "phage"])

sum(Mastertable_viral_rarefied_completeness_est$Totalnumberofreads[Mastertable_viral_rarefied_completeness_est$Final_class2 == "Caudoviricetes"])
(sum(Mastertable_viral_rarefied_completeness_est$Totalnumberofreads[Mastertable_viral_rarefied_completeness_est$Final_class2 == "Caudoviricetes"])/sum(Mastertable_viral_rarefied_completeness_est$Totalnumberofreads[Mastertable_viral_rarefied_completeness_est$Final_superkingdom == "Viruses"]))*100

sum(Mastertable_viral_rarefied_completeness_est$Totalnumberofreads[Mastertable_viral_rarefied_completeness_est$Final_class == "Caudoviricetes (non-CrAss)"])
(sum(Mastertable_viral_rarefied_completeness_est$Totalnumberofreads[Mastertable_viral_rarefied_completeness_est$Final_class == "Caudoviricetes (non-CrAss)"])/sum(Mastertable_viral_rarefied_completeness_est$Totalnumberofreads[Mastertable_viral_rarefied_completeness_est$Final_superkingdom == "Viruses"]))*100

sum(Mastertable_viral_rarefied_completeness_est$Totalnumberofreads[Mastertable_viral_rarefied_completeness_est$Final_class == "Caudoviricetes (CrAss)"])
(sum(Mastertable_viral_rarefied_completeness_est$Totalnumberofreads[Mastertable_viral_rarefied_completeness_est$Final_class == "Caudoviricetes (CrAss)"])/sum(Mastertable_viral_rarefied_completeness_est$Totalnumberofreads[Mastertable_viral_rarefied_completeness_est$Final_superkingdom == "Viruses"]))*100 

sum(Mastertable_viral_rarefied_completeness_est$Totalnumberofreads[Mastertable_viral_rarefied_completeness_est$Final_class == "Malgrandaviricetes"])
(sum(Mastertable_viral_rarefied_completeness_est$Totalnumberofreads[Mastertable_viral_rarefied_completeness_est$Final_class == "Malgrandaviricetes"])/sum(Mastertable_viral_rarefied_completeness_est$Totalnumberofreads[Mastertable_viral_rarefied_completeness_est$Final_superkingdom == "Viruses"]))*100

sum(Mastertable_viral_rarefied_completeness_est$Totalnumberofreads[Mastertable_viral_rarefied_completeness_est$Final_class == "Unannotated" & Mastertable_viral_rarefied_completeness_est$Final_viral == "phage"])
(sum(Mastertable_viral_rarefied_completeness_est$Totalnumberofreads[Mastertable_viral_rarefied_completeness_est$Final_class == "Unannotated" & Mastertable_viral_rarefied_completeness_est$Final_viral == "phage"])/sum(Mastertable_viral_rarefied_completeness_est$Totalnumberofreads[Mastertable_viral_rarefied_completeness_est$Final_superkingdom == "Viruses"]))*100

sum(Mastertable_viral_rarefied_completeness_est$Totalnumberofreads[Mastertable_viral_rarefied_completeness_est$Final_class == "Faserviricetes"])
(sum(Mastertable_viral_rarefied_completeness_est$Totalnumberofreads[Mastertable_viral_rarefied_completeness_est$Final_class == "Faserviricetes"])/sum(Mastertable_viral_rarefied_completeness_est$Totalnumberofreads[Mastertable_viral_rarefied_completeness_est$Final_superkingdom == "Viruses"]))*100

sum(Mastertable_viral_rarefied_completeness_est$Totalnumberofreads[Mastertable_viral_rarefied_completeness_est$Final_class == "Duplopiviricetes"])
(sum(Mastertable_viral_rarefied_completeness_est$Totalnumberofreads[Mastertable_viral_rarefied_completeness_est$Final_class == "Duplopiviricetes"])/sum(Mastertable_viral_rarefied_completeness_est$Totalnumberofreads[Mastertable_viral_rarefied_completeness_est$Final_superkingdom == "Viruses"]))*100

sum(Mastertable_viral_rarefied_completeness_est$Totalnumberofreads[Mastertable_viral_rarefied_completeness_est$Final_class == "Allassoviricetes"])
(sum(Mastertable_viral_rarefied_completeness_est$Totalnumberofreads[Mastertable_viral_rarefied_completeness_est$Final_class == "Allassoviricetes"])/sum(Mastertable_viral_rarefied_completeness_est$Totalnumberofreads[Mastertable_viral_rarefied_completeness_est$Final_superkingdom == "Viruses"]))*100

# prokaryotic viral class completeness (average length): 
sum(Mastertable_viral_rarefied_completeness_est$completeness[Mastertable_viral_rarefied_completeness_est$Final_class2 == "Caudoviricetes"])/1909
sum(Mastertable_viral_rarefied_completeness_est$completeness[Mastertable_viral_rarefied_completeness_est$Final_class == "Caudoviricetes (non-CrAss)"])/1835
sum(Mastertable_viral_rarefied_completeness_est$completeness[Mastertable_viral_rarefied_completeness_est$Final_class == "Caudoviricetes (CrAss)"])/74
sum(Mastertable_viral_rarefied_completeness_est$completeness[Mastertable_viral_rarefied_completeness_est$Final_class == "Malgrandaviricetes"])/757
sum(Mastertable_viral_rarefied_completeness_est$completeness[Mastertable_viral_rarefied_completeness_est$Final_class == "Unannotated" & Mastertable_viral_rarefied_completeness_est$Final_viral == "phage"])/200
sum(Mastertable_viral_rarefied_completeness_est$completeness[Mastertable_viral_rarefied_completeness_est$Final_class == "Duplopiviricetes"])/65
sum(Mastertable_viral_rarefied_completeness_est$completeness[Mastertable_viral_rarefied_completeness_est$Final_class == "Faserviricetes"])/39
sum(Mastertable_viral_rarefied_completeness_est$completeness[Mastertable_viral_rarefied_completeness_est$Final_class == "Allassoviricetes"])/2
 
# prokaryotic viral families median length: 
mean(Mastertable_viral_rarefied_completeness_est$length[Mastertable_viral_rarefied_completeness_est$Final_class2 == "Caudoviricetes"])
mean(Mastertable_viral_rarefied_completeness_est$length[Mastertable_viral_rarefied_completeness_est$Final_class == "Caudoviricetes (non-CrAss)"])
mean(Mastertable_viral_rarefied_completeness_est$length[Mastertable_viral_rarefied_completeness_est$Final_class == "Caudoviricetes (CrAss)"])
mean(Mastertable_viral_rarefied_completeness_est$length[Mastertable_viral_rarefied_completeness_est$Final_class == "Malgrandaviricetes"])
mean(Mastertable_viral_rarefied_completeness_est$length[Mastertable_viral_rarefied_completeness_est$Final_class == "Unannotated" & Mastertable_viral_rarefied_completeness_est$Final_viral == "phage"])
mean(Mastertable_viral_rarefied_completeness_est$length[Mastertable_viral_rarefied_completeness_est$Final_class == "Duplopiviricetes"])
mean(Mastertable_viral_rarefied_completeness_est$length[Mastertable_viral_rarefied_completeness_est$Final_class == "Faserviricetes"])
mean(Mastertable_viral_rarefied_completeness_est$length[Mastertable_viral_rarefied_completeness_est$Final_class == "Allassoviricetes"])

# prokaryotic viral families length range: 
range(Mastertable_viral_rarefied_completeness_est$length[Mastertable_viral_rarefied_completeness_est$Final_class2 == "Caudoviricetes"])
range(Mastertable_viral_rarefied_completeness_est$length[Mastertable_viral_rarefied_completeness_est$Final_class == "Caudoviricetes (non-CrAss)"])
range(Mastertable_viral_rarefied_completeness_est$length[Mastertable_viral_rarefied_completeness_est$Final_class == "Caudoviricetes (CrAss)"])
range(Mastertable_viral_rarefied_completeness_est$length[Mastertable_viral_rarefied_completeness_est$Final_class == "Malgrandaviricetes"])
range(Mastertable_viral_rarefied_completeness_est$length[Mastertable_viral_rarefied_completeness_est$Final_class == "Unannotated" & Mastertable_viral_rarefied_completeness_est$Final_viral == "phage"])
range(Mastertable_viral_rarefied_completeness_est$length[Mastertable_viral_rarefied_completeness_est$Final_class == "Duplopiviricetes"])
range(Mastertable_viral_rarefied_completeness_est$length[Mastertable_viral_rarefied_completeness_est$Final_class == "Faserviricetes"])
range(Mastertable_viral_rarefied_completeness_est$length[Mastertable_viral_rarefied_completeness_est$Final_class == "Allassoviricetes"])
#View(Mastertable_viral_rarefied_completeness_est[Mastertable_viral_rarefied_completeness_est$Final_class2 == "Caudoviricetes" & Mastertable_viral_rarefied_completeness_est$length < 10,])

# lysogenic phages
Mastertable_viral_rarefied_completeness_est$lysogenic.cycle[is.na(Mastertable_viral_rarefied_completeness_est$lysogenic.cycle)] <- 0
table(Mastertable_viral_rarefied_completeness_est$lysogenic.cycle[Mastertable_viral_rarefied_completeness_est$Final_viral == "phage"])
table(Mastertable_viral_rarefied_completeness_est$lysogenic.cycle[Mastertable_viral_rarefied_completeness_est$Final_class2 == "Caudoviricetes"])
table(Mastertable_viral_rarefied_completeness_est$lysogenic.cycle[Mastertable_viral_rarefied_completeness_est$Final_class == "Caudoviricetes (non-CrAss)"])
table(Mastertable_viral_rarefied_completeness_est$lysogenic.cycle[Mastertable_viral_rarefied_completeness_est$Final_class == "Caudoviricetes (CrAss)"])
table(Mastertable_viral_rarefied_completeness_est$lysogenic.cycle[Mastertable_viral_rarefied_completeness_est$Final_class == "Malgrandaviricetes"])
table(Mastertable_viral_rarefied_completeness_est$lysogenic.cycle[Mastertable_viral_rarefied_completeness_est$Final_class == "Unannotated" & Mastertable_viral_rarefied_completeness_est$Final_viral == "phage"])
table(Mastertable_viral_rarefied_completeness_est$lysogenic.cycle[Mastertable_viral_rarefied_completeness_est$Final_class == "Duplopiviricetes"])
table(Mastertable_viral_rarefied_completeness_est$lysogenic.cycle[Mastertable_viral_rarefied_completeness_est$Final_class == "Faserviricetes"])
table(Mastertable_viral_rarefied_completeness_est$lysogenic.cycle[Mastertable_viral_rarefied_completeness_est$Final_class == "Allassoviricetes"])

# CrAssPhages:
table(Mastertable_viral_rarefied_completeness_est$CrAss_taxonomy[Mastertable_viral_rarefied_completeness_est$Final_class2 == "Caudoviricetes"])
sum(Mastertable_viral_rarefied_completeness_est$Totalnumberofreads[Mastertable_viral_rarefied_completeness_est$Final_class == "Caudoviricetes" & Mastertable_viral_rarefied_completeness_est$CrassPhage == "crAssphage"])
sum(Mastertable_viral_rarefied_completeness_est$completeness[Mastertable_viral_rarefied_completeness_est$Final_class == "Caudoviricetes" & Mastertable_viral_rarefied_completeness_est$CrassPhage == "crAssphage"])/73
median(Mastertable_viral_rarefied_completeness_est$length[Mastertable_viral_rarefied_completeness_est$Final_class == "Caudoviricetes" & Mastertable_viral_rarefied_completeness_est$CrassPhage == "crAssphage"])
range(Mastertable_viral_rarefied_completeness_est$length[Mastertable_viral_rarefied_completeness_est$Final_class == "Caudoviricetes" & Mastertable_viral_rarefied_completeness_est$CrassPhage == "crAssphage"])
table(Mastertable_viral_rarefied_completeness_est$lysogenic.cycle[Mastertable_viral_rarefied_completeness_est$Final_class == "Caudoviricetes" & Mastertable_viral_rarefied_completeness_est$CrassPhage == "crAssphage"])
table(Mastertable_viral_rarefied_completeness_est$Toxin.presence[Mastertable_viral_rarefied_completeness_est$Final_class == "Caudoviricetes" & Mastertable_viral_rarefied_completeness_est$CrassPhage == "crAssphage"])

# Non-CrAssPhages
table(Mastertable_viral_rarefied_completeness_est$CrassPhage[Mastertable_viral_rarefied_completeness_est$Final_class == "Caudoviricetes"])
sum(Mastertable_viral_rarefied_completeness_est$Totalnumberofreads[Mastertable_viral_rarefied_completeness_est$Final_class == "Caudoviricetes" & !Mastertable_viral_rarefied_completeness_est$CrassPhage == "crAssphage"])
sum(Mastertable_viral_rarefied_completeness_est$completeness[Mastertable_viral_rarefied_completeness_est$Final_class == "Caudoviricetes" & !Mastertable_viral_rarefied_completeness_est$CrassPhage == "crAssphage"])/1826
median(Mastertable_viral_rarefied_completeness_est$length[Mastertable_viral_rarefied_completeness_est$Final_class == "Caudoviricetes" & !Mastertable_viral_rarefied_completeness_est$CrassPhage == "crAssphage"])
range(Mastertable_viral_rarefied_completeness_est$length[Mastertable_viral_rarefied_completeness_est$Final_class == "Caudoviricetes" & !Mastertable_viral_rarefied_completeness_est$CrassPhage == "crAssphage"])
table(Mastertable_viral_rarefied_completeness_est$lysogenic.cycle[Mastertable_viral_rarefied_completeness_est$Final_class == "Caudoviricetes" & !Mastertable_viral_rarefied_completeness_est$CrassPhage == "crAssphage"])
table(Mastertable_viral_rarefied_completeness_est$Toxin.presence[Mastertable_viral_rarefied_completeness_est$Final_class == "Caudoviricetes" & !Mastertable_viral_rarefied_completeness_est$CrassPhage == "crAssphage"])

# Toxins
Mastertable_viral_rarefied_completeness_est$Toxin.presence[is.na(Mastertable_viral_rarefied_completeness_est$Toxin.presence)] <- 0
table(Mastertable_viral_rarefied_completeness_est$Toxin.presence[Mastertable_viral_rarefied_completeness_est$Final_viral == "phage"])
table(Mastertable_viral_rarefied_completeness_est$Toxin.presence[Mastertable_viral_rarefied_completeness_est$Final_class2 == "Caudoviricetes"])
table(Mastertable_viral_rarefied_completeness_est$Toxin.presence[Mastertable_viral_rarefied_completeness_est$Final_class == "Caudoviricetes (non-CrAss)"])
table(Mastertable_viral_rarefied_completeness_est$Toxin.presence[Mastertable_viral_rarefied_completeness_est$Final_class == "Caudoviricetes (CrAss)"])
table(Mastertable_viral_rarefied_completeness_est$Toxin.presence[Mastertable_viral_rarefied_completeness_est$Final_class == "Malgrandaviricetes"])
table(Mastertable_viral_rarefied_completeness_est$Toxin.presence[Mastertable_viral_rarefied_completeness_est$Final_class == "Unannotated" & Mastertable_viral_rarefied_completeness_est$Final_viral == "phage"])
table(Mastertable_viral_rarefied_completeness_est$Toxin.presence[Mastertable_viral_rarefied_completeness_est$Final_class == "Duplopiviricetes"])
table(Mastertable_viral_rarefied_completeness_est$Toxin.presence[Mastertable_viral_rarefied_completeness_est$Final_class == "Faserviricetes"])
table(Mastertable_viral_rarefied_completeness_est$Toxin.presence[Mastertable_viral_rarefied_completeness_est$Final_class == "Allassoviricetes"])

# SUMMARY:
# -----------------> Caudoviricetes (caudo); 1898 NR-contigs; 136803739 reads; 84,4% complete; medium-quality-tier, 36,7 kb length; range: 10.9kb-278.8kb; temperate: 1222; CrAssphage; 73; toxins: 178
# -----------------subcategory ----------------> CrAssphages; 73 NR-contigs; 33865807 reads; 92,91% complete, high-quality-tier; 98,7kb length, range: 41,4kb-190,92kb; temperate: 8 (~ 11,0); toxin: 5
# -----------------subcategory ----------------> non-CrAssphages; 1825 NR-contigs; 102937932 reads; 84,1% complete, Medium-quality-tier; 36,1kb length, range: 10,9kb-278,8kb; temperate: 1214 (~ 66,5%); toxin: 173
# -----------------> Malgrandaviricetes (micro); 743 NR-contigs; 90218005 reads; 94,5% complete; high-quality-tier, 5,9 kb length; range: 2.7kb-18,1kb; temperate: 0; toxin: 0
# -----------------> Unannotated; 440 NR-contigs; 6586276 reads; 36,9% complete; low-quality-tier, 6,1 kb length; range: 2,4 - 208,6kb; temperate: 140 (~ 31,8%); toxin: 24
# -----------------> Faserviricetes (ino); 18 NR-contigs; 2058821 reads; 93,4% complete; high-quality-tier,7,3 kb length; range: 3,3kb - 11,2 kb; temperate: 1; toxin: 0
# -----------------> Allassoviricetes (levi); 2 NR-contigs; 111516 reads; 100% complete; high-quality tier, 3,89 kb length; range: 3,5kb -4,2kb; temperate: 0; toxin: 0
# -----------------> Tectiliviricetes; 1 NR-contigs; 899 reads; 98,22% complete; high-quality tier, 39,5 kb length; range: 39,5kb; temperate: 1; toxin: 0

# Bacterial and archaeal DNA viruses are often capable of integrating their genomes into the host chromosome thereby becoming proviruses. Even though proviruses related to Microviridae have not been previously reported, we set out to verify this possibility by performing searches against genomic sequences available in public databases. 
# Microviridae majorly not temperate phages

## do eukaryotic viral part
## make barplot (see pie chart) & put next to this table
## check circlator test2

#____________________________________________________________________________________________________________________________
# Average completeness for eukaryotic viruses?
## !!!! Take into account the CheckV wasn't made for eukaryotic viruses, often low completeness% is given but %cov & %ANI is very high.
table(Mastertable_viral_rarefied_completeness_est$Final_viral)
table(Mastertable_viral_rarefied_completeness_est$Final_family[Mastertable_viral_rarefied_completeness_est$Final_viral == "eukaryotic virus"])
nrow(Mastertable_viral_rarefied_completeness_est[Mastertable_viral_rarefied_completeness_est$Final_viral == "eukaryotic virus",])

#View(Mastertable_viral_rarefied_completeness_est[Mastertable_viral_rarefied_completeness_est$Final_viral == "eukaryotic virus",])
sum(Mastertable_viral_rarefied_completeness_est$Totalnumberofreads[Mastertable_viral_rarefied_completeness_est$Final_viral == "eukaryotic virus"]) # 562215 eukaryotic viral reads
nrow(Mastertable_viral_rarefied_completeness_est[Mastertable_viral_rarefied_completeness_est$Final_viral == "eukaryotic virus",]) # 45
sum(Mastertable_viral_rarefied_completeness_est$completeness[Mastertable_viral_rarefied_completeness_est$Final_viral == "eukaryotic virus"])/61
median(Mastertable_viral_rarefied_completeness_est$length[Mastertable_viral_rarefied_completeness_est$Final_viral == "eukaryotic virus"])
range(Mastertable_viral_rarefied_completeness_est$length[Mastertable_viral_rarefied_completeness_est$Final_viral == "eukaryotic virus"])

# --------------------------> 75,84%; medium-quality-tier; 44 NR-contigs; 562215 reads

#____________________________________________________________________________________________________________________________
# eukaryotic viral families reads: 
table(Mastertable_viral_rarefied_completeness_est$Final_family[Mastertable_viral_rarefied_completeness_est$Final_viral == "eukaryotic virus"])
sum(Mastertable_viral_rarefied_completeness_est$Totalnumberofreads[Mastertable_viral_rarefied_completeness_est$Final_family == "Anelloviridae"])
(sum(Mastertable_viral_rarefied_completeness_est$Totalnumberofreads[Mastertable_viral_rarefied_completeness_est$Final_family == "Anelloviridae"])/sum(Mastertable_viral_rarefied_completeness_est$Totalnumberofreads[Mastertable_viral_rarefied_completeness_est$Final_superkingdom == "Viruses"]))*100

sum(Mastertable_viral_rarefied_completeness_est$Totalnumberofreads[Mastertable_viral_rarefied_completeness_est$Final_family == "Virgaviridae"])
(sum(Mastertable_viral_rarefied_completeness_est$Totalnumberofreads[Mastertable_viral_rarefied_completeness_est$Final_family == "Virgaviridae"])/sum(Mastertable_viral_rarefied_completeness_est$Totalnumberofreads[Mastertable_viral_rarefied_completeness_est$Final_superkingdom == "Viruses"]))*100

sum(Mastertable_viral_rarefied_completeness_est$Totalnumberofreads[Mastertable_viral_rarefied_completeness_est$Final_family == "Unannotated" & Mastertable_viral_rarefied_completeness_est$Final_viral == "eukaryotic virus"])
(sum(Mastertable_viral_rarefied_completeness_est$Totalnumberofreads[Mastertable_viral_rarefied_completeness_est$Final_family == "Unannotated" & Mastertable_viral_rarefied_completeness_est$Final_viral == "eukaryotic virus"])/sum(Mastertable_viral_rarefied_completeness_est$Totalnumberofreads[Mastertable_viral_rarefied_completeness_est$Final_superkingdom == "Viruses"]))*100

sum(Mastertable_viral_rarefied_completeness_est$Totalnumberofreads[Mastertable_viral_rarefied_completeness_est$Final_family == "Caliciviridae"])
sum(Mastertable_viral_rarefied_completeness_est$Totalnumberofreads[Mastertable_viral_rarefied_completeness_est$Final_family == "Circoviridae"])
sum(Mastertable_viral_rarefied_completeness_est$Totalnumberofreads[Mastertable_viral_rarefied_completeness_est$Final_family == "Chrysoviridae"])
sum(Mastertable_viral_rarefied_completeness_est$Totalnumberofreads[Mastertable_viral_rarefied_completeness_est$Final_family == "Genomoviridae"])
sum(Mastertable_viral_rarefied_completeness_est$Totalnumberofreads[Mastertable_viral_rarefied_completeness_est$Final_family == "Picornaviridae"])
sum(Mastertable_viral_rarefied_completeness_est$Totalnumberofreads[Mastertable_viral_rarefied_completeness_est$Final_family == "Secoviridae"])
sum(Mastertable_viral_rarefied_completeness_est$Totalnumberofreads[Mastertable_viral_rarefied_completeness_est$Final_family == "Alphaflexiviridae"])
sum(Mastertable_viral_rarefied_completeness_est$Totalnumberofreads[Mastertable_viral_rarefied_completeness_est$Final_family == "Betaflexiviridae"])
sum(Mastertable_viral_rarefied_completeness_est$Totalnumberofreads[Mastertable_viral_rarefied_completeness_est$Final_family == "Papillomaviridae"])
sum(Mastertable_viral_rarefied_completeness_est$Totalnumberofreads[Mastertable_viral_rarefied_completeness_est$Final_family == "Smacoviridae"])
sum(Mastertable_viral_rarefied_completeness_est$Totalnumberofreads[Mastertable_viral_rarefied_completeness_est$Final_family == "Astroviridae"])
sum(Mastertable_viral_rarefied_completeness_est$Totalnumberofreads[Mastertable_viral_rarefied_completeness_est$Final_family == "Closteroviridae"])
sum(Mastertable_viral_rarefied_completeness_est$Totalnumberofreads[Mastertable_viral_rarefied_completeness_est$Final_family == "Endornaviridae"])
sum(Mastertable_viral_rarefied_completeness_est$Totalnumberofreads[Mastertable_viral_rarefied_completeness_est$Final_family == "Geminiviridae"])
sum(Mastertable_viral_rarefied_completeness_est$Totalnumberofreads[Mastertable_viral_rarefied_completeness_est$Final_family == "Parvoviridae"])

# eukaryotic viral families completeness (average length):
mean(Mastertable_viral_rarefied_completeness_est$completeness[Mastertable_viral_rarefied_completeness_est$Final_family == "Anelloviridae"])
mean(Mastertable_viral_rarefied_completeness_est$completeness[Mastertable_viral_rarefied_completeness_est$Final_family == "Virgaviridae"])
mean(Mastertable_viral_rarefied_completeness_est$completeness[Mastertable_viral_rarefied_completeness_est$Final_family == "Unannotated" & Mastertable_viral_rarefied_completeness_est$Final_viral == "eukaryotic virus"])
mean(Mastertable_viral_rarefied_completeness_est$completeness[Mastertable_viral_rarefied_completeness_est$Final_family == "Caliciviridae"])
mean(Mastertable_viral_rarefied_completeness_est$completeness[Mastertable_viral_rarefied_completeness_est$Final_family == "Circoviridae"])
mean(Mastertable_viral_rarefied_completeness_est$completeness[Mastertable_viral_rarefied_completeness_est$Final_family == "Chrysoviridae"])
mean(Mastertable_viral_rarefied_completeness_est$completeness[Mastertable_viral_rarefied_completeness_est$Final_family == "Genomoviridae"])
mean(Mastertable_viral_rarefied_completeness_est$completeness[Mastertable_viral_rarefied_completeness_est$Final_family == "Picornaviridae"])
mean(Mastertable_viral_rarefied_completeness_est$completeness[Mastertable_viral_rarefied_completeness_est$Final_family == "Secoviridae"])
mean(Mastertable_viral_rarefied_completeness_est$completeness[Mastertable_viral_rarefied_completeness_est$Final_family == "Alphaflexiviridae"])
mean(Mastertable_viral_rarefied_completeness_est$completeness[Mastertable_viral_rarefied_completeness_est$Final_family == "Betaflexiviridae"])
mean(Mastertable_viral_rarefied_completeness_est$completeness[Mastertable_viral_rarefied_completeness_est$Final_family == "Papillomaviridae"])
mean(Mastertable_viral_rarefied_completeness_est$completeness[Mastertable_viral_rarefied_completeness_est$Final_family == "Smacoviridae"])
mean(Mastertable_viral_rarefied_completeness_est$completeness[Mastertable_viral_rarefied_completeness_est$Final_family == "Astroviridae"])
mean(Mastertable_viral_rarefied_completeness_est$completeness[Mastertable_viral_rarefied_completeness_est$Final_family == "Closteroviridae"])
mean(Mastertable_viral_rarefied_completeness_est$completeness[Mastertable_viral_rarefied_completeness_est$Final_family == "Endornaviridae"])
mean(Mastertable_viral_rarefied_completeness_est$completeness[Mastertable_viral_rarefied_completeness_est$Final_family == "Geminiviridae"])
mean(Mastertable_viral_rarefied_completeness_est$completeness[Mastertable_viral_rarefied_completeness_est$Final_family == "Parvoviridae"])
View(Mastertable_viral_rarefied_completeness_est[Mastertable_viral_rarefied_completeness_est$Final_viral == "eukaryotic virus",])

# eukaryotic viral families length: 
mean(Mastertable_viral_rarefied_completeness_est$length[Mastertable_viral_rarefied_completeness_est$Final_family == "Anelloviridae"])
mean(Mastertable_viral_rarefied_completeness_est$length[Mastertable_viral_rarefied_completeness_est$Final_family == "Virgaviridae"])
mean(Mastertable_viral_rarefied_completeness_est$length[Mastertable_viral_rarefied_completeness_est$Final_family == "Unannotated" & Mastertable_viral_rarefied_completeness_est$Final_viral == "eukaryotic virus"])
mean(Mastertable_viral_rarefied_completeness_est$length[Mastertable_viral_rarefied_completeness_est$Final_family == "Caliciviridae"])
mean(Mastertable_viral_rarefied_completeness_est$length[Mastertable_viral_rarefied_completeness_est$Final_family == "Circoviridae"])
mean(Mastertable_viral_rarefied_completeness_est$length[Mastertable_viral_rarefied_completeness_est$Final_family == "Chrysoviridae"])
mean(Mastertable_viral_rarefied_completeness_est$length[Mastertable_viral_rarefied_completeness_est$Final_family == "Genomoviridae"])
mean(Mastertable_viral_rarefied_completeness_est$length[Mastertable_viral_rarefied_completeness_est$Final_family == "Picornaviridae"])
mean(Mastertable_viral_rarefied_completeness_est$length[Mastertable_viral_rarefied_completeness_est$Final_family == "Secoviridae"])
mean(Mastertable_viral_rarefied_completeness_est$length[Mastertable_viral_rarefied_completeness_est$Final_family == "Alphaflexiviridae"])
mean(Mastertable_viral_rarefied_completeness_est$length[Mastertable_viral_rarefied_completeness_est$Final_family == "Betaflexiviridae"])
mean(Mastertable_viral_rarefied_completeness_est$length[Mastertable_viral_rarefied_completeness_est$Final_family == "Papillomaviridae"])
mean(Mastertable_viral_rarefied_completeness_est$length[Mastertable_viral_rarefied_completeness_est$Final_family == "Smacoviridae"])
mean(Mastertable_viral_rarefied_completeness_est$length[Mastertable_viral_rarefied_completeness_est$Final_family == "Astroviridae"])
mean(Mastertable_viral_rarefied_completeness_est$length[Mastertable_viral_rarefied_completeness_est$Final_family == "Closteroviridae"])
mean(Mastertable_viral_rarefied_completeness_est$length[Mastertable_viral_rarefied_completeness_est$Final_family == "Endornaviridae"])
mean(Mastertable_viral_rarefied_completeness_est$length[Mastertable_viral_rarefied_completeness_est$Final_family == "Geminiviridae"])
mean(Mastertable_viral_rarefied_completeness_est$length[Mastertable_viral_rarefied_completeness_est$Final_family == "Parvoviridae"])

# eukaryotic viral families length range: 
range(Mastertable_viral_rarefied_completeness_est$length[Mastertable_viral_rarefied_completeness_est$Final_family == "Anelloviridae"])
range(Mastertable_viral_rarefied_completeness_est$length[Mastertable_viral_rarefied_completeness_est$Final_family == "Virgaviridae"])
range(Mastertable_viral_rarefied_completeness_est$length[Mastertable_viral_rarefied_completeness_est$Final_family == "Unannotated" & Mastertable_viral_rarefied_completeness_est$Final_viral == "eukaryotic virus"])
range(Mastertable_viral_rarefied_completeness_est$length[Mastertable_viral_rarefied_completeness_est$Final_family == "Caliciviridae"])
range(Mastertable_viral_rarefied_completeness_est$length[Mastertable_viral_rarefied_completeness_est$Final_family == "Circoviridae"])
range(Mastertable_viral_rarefied_completeness_est$length[Mastertable_viral_rarefied_completeness_est$Final_family == "Chrysoviridae"])
range(Mastertable_viral_rarefied_completeness_est$length[Mastertable_viral_rarefied_completeness_est$Final_family == "Genomoviridae"])
range(Mastertable_viral_rarefied_completeness_est$length[Mastertable_viral_rarefied_completeness_est$Final_family == "Picornaviridae"])
range(Mastertable_viral_rarefied_completeness_est$length[Mastertable_viral_rarefied_completeness_est$Final_family == "Secoviridae"])
range(Mastertable_viral_rarefied_completeness_est$length[Mastertable_viral_rarefied_completeness_est$Final_family == "Alphaflexiviridae"])
range(Mastertable_viral_rarefied_completeness_est$length[Mastertable_viral_rarefied_completeness_est$Final_family == "Betaflexiviridae"])
range(Mastertable_viral_rarefied_completeness_est$length[Mastertable_viral_rarefied_completeness_est$Final_family == "Papillomaviridae"])
range(Mastertable_viral_rarefied_completeness_est$length[Mastertable_viral_rarefied_completeness_est$Final_family == "Smacoviridae"])
range(Mastertable_viral_rarefied_completeness_est$length[Mastertable_viral_rarefied_completeness_est$Final_family == "Astroviridae"])
range(Mastertable_viral_rarefied_completeness_est$length[Mastertable_viral_rarefied_completeness_est$Final_family == "Closteroviridae"])
range(Mastertable_viral_rarefied_completeness_est$length[Mastertable_viral_rarefied_completeness_est$Final_family == "Endornaviridae"])
range(Mastertable_viral_rarefied_completeness_est$length[Mastertable_viral_rarefied_completeness_est$Final_family == "Geminiviridae"])
range(Mastertable_viral_rarefied_completeness_est$length[Mastertable_viral_rarefied_completeness_est$Final_family == "Parvoviridae"])
View(Mastertable_viral_rarefied_completeness_est[Mastertable_viral_rarefied_completeness_est$Final_viral == "eukaryotic virus",])

# -----------------> Virgaviridae; 7 NR-contigs; 234080 reads; 61% complete; medium-quality-tier; 3,4kb length; range: 2.42-6.5; host group: plants/fungi
# -----------------> Anelloviridae; 6 NR-contigs; 5807 reads; 90,88% complete; High-quality-tier; 2,6kb length; range: 2.1-3.5; host group: small circular 
# -----------------> Circoviridae; 4 NR-contigs; 47230 reads; 92,59% complete; High-quality-tier; 3,1kb length; range: 2.8-4.1; host group: small circular 
# -----------------> Unannotated (~ Cressdnaviricota); 4 NR-contigs; 118362 reads; 89,81% complete; medium-quality-tier; 5,4kb length; range: 3.6-6.2; host group: small circular 
# -----------------> Caliciviridae; 3 NR-contigs; 4324 reads; 98,39% complete; High-quality-tier; 7,5kb length; range: 7.2-7.7; host group: human
# -----------------> Genomoviridae; 3 NR-contigs; 1875 reads; 100% complete; Complete; 2,3kb length; range: 2.2-2.3; host group: small circular 
# -----------------> Alphaflexiviridae; 2 NR-contigs; 82513 reads; 95,22% complete; High-quality-tier; 6,4kb length; range: 6.4-6.4; host group: plants/fungi
# -----------------> Betaflexiviridae; 2 NR-contigs; 4152 reads; 76,01% complete; medium-quality-tier;  6,3kb length; range: 4.2-8.3; host group: plants/fungi
# -----------------> Chrysoviridae; 2 NR-contigs; 6564 reads; 26,59% complete; Low-quality-tier; 3,1kb length; range: 2.8-3.4; host group: plants/fungi
# -----------------> Papillomaviridae; 2 NR-contigs; 551 reads; 87,76% complete; medium-quality-tier; 7,0kb length; range: 6.1-7.9; host group: human
# -----------------> Smacoviridae; 2 NR-contigs; 51183 reads; 100% complete; Complete-tier; 2,5kb length; range: NA; host group: small circular 
# -----------------> Astroviridae; 1 NR-contigs; 257 reads; 8,89% complete; Low-quality-tier;  7,2kb length; range: NA; host group: human
# -----------------> Closteroviridae; 1 NR-contigs; 1132 reads; 40,25% complete; Low-quality-tier; 6,9kb length; range: NA; host group: plants/fungi
# -----------------> Endornaviridae; 1 NR-contigs; 1881 reads; 1,2% complete; Low-quality-tier; 2,1kb length; range: NA; host group: plants/fungi
# -----------------> Geminiviridae; 1 NR-contigs; 253 reads; 56,92% complete; medium-quality-tier; 2,9kb length; range: NA; host group: small circular 
# -----------------> Picornaviridae; 1 NR-contigs; 707 reads; 100% complete; High-quality-tier; 7,4kb length; range: NA; host group: human 

#____________________________________________________________________________________________________________________________

median(Mastertable_viral_rarefied_completeness_est$length[Mastertable_viral_rarefied_completeness_est$Final_superkingdom == "Viruses"])
table(Mastertable_viral_rarefied_completeness_est$Final_superkingdom)
## Add host group, being small circular viruses, animal viruses, plant/fungal viruses (just name plant if no fungi present)
## Add eukaryotic viruses first herein because this will be brief

## Conclusions: 
# 1) Small eukaryotic virome; mostly passenger viruses
# 2) Phageome annotated on class (monophyletic groups) dominated with Caudoviricetes phages
# 3) Range and median of phage length fits with what is expected for these phage classes
# 4) Most viruses are near-complete, reducing bias for diversity analysis yet to come
# 5) Little to no RNA phages detected
# 6) small ssDNA phages (~ micro's) second most abundant group
# 7) NO lysogenic phages in micro group, most of them present in Caudo group
# 8) Most of the phages are annotated (~ 85,6%)
#   -> Toxins in caudo's: HicA_toxin_of_bacterial_toxin-antitoxin OR HicA_antitoxin_of_bacterial_toxin-antitoxin; Toxin-antitoxin (TA) systems are little genetic units generally composed of two genes encoding antitoxin and toxin. These systems are known to be involved in many functions that can lead to growth arrest and cell death. Among the different types of TA systems, the type II gathers together systems where the antitoxin directly binds and inhibits the toxin. Among these type II TA systems, the HicAB module is widely distributed in free-living Bacteria and Archaea and the toxin HicA functions via RNA binding and cleavage. The genome of the symbiotic Sinorhizobium meliloti encodes numerous TA systems and only a few of them are functional.
#   -> Toxins in caudo's: do we find only toxin, antitoxin, or both?

#   -> Toxins in micro's; absent after rarefaction; Zonular_occludens_toxin_(Zot); Zonula occludens toxin (Zot) is an enterotoxin elaborated by Vibrio cholerae that increases intestinal permeability by interacting with a mammalian cell receptor with subsequent activation of intracellular signaling leading to the disassembly of the intercellular tight junctions.
###################################
# 4. Summary: all mastertables with normalization for high-quality genomes (~ phages)
###################################
## Take into account that most of diversity analysis already correct for sequencing depth (alpha-diversity), so you have to feed it the unrarefied completeness data

# global mastertable
Mastertable
Mastertable_unrarefied
Mastertable_rarefied

# viral mastertable
# unrarefied
Mastertable_viral_unrarefied

# rarefied
Mastertable_viral_rarefied
####################################

######### !!!!!!!!!!!!!!! #########
# SAVE Global Environment #
######### !!!!!!!!!!!!!!! #########
