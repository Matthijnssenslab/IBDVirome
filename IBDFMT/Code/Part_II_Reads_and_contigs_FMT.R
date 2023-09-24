####################################
# SCRIPT 2: CREATE READ AND CONTIG FIGURE + GUIDANCE FIGURES FOR RAREFACTION
####################################
# All of the data is stored in "/Users/daan/Desktop/Transfer/read_count"

# -------------------> Load global mastertable of  "1. Create Mastertable"
# -------------------> Create directory folder called "Reads"  in /Users/daan/Desktop/Bioinformatics/Analysis/Biologicals"

# Read files:
# ------------------> I+T_RawPairedReads.txt
# ------------------> I+T_TrimmedReads.txt
# ------------------> I+T_NCOutReads.txt
# ------------------> I+T_QualityControlledReads.txt
# ------------------> I+T_Mapped_bwa2_NR_Reads_threshold.txt (mapped)

# Contigs files:
# ------------------> NumberOfAssembledContigs.txt
# ------------------> NumberOfLengthFilteredContigs.txt (1 kb cutoff)
# ------------------> I+T_NumberOfNonRedundantContigs.txt
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
#install.packages("WriteXLS")
#install.packages("circlize")
#install.packages("esquisse")
#installed.packages("modeldata")
#install.packages("tidyverse")
#install.packages("ggThemeAssist")
library(esquisse)
library(modeldata)
library(tidyverse)
library(ggplot2)
library(circlize)
library(WriteXLS)
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
library(ggThemeAssist)
####################################
# 1. Reads
####################################
# 1.1 Insert files
####################################
setwd("/Users/daan/Desktop//Transfer/Read_Contig_count")
dir()
####################################
# 1.1.1 Raw reads
####################################
raw_reads <- as.data.frame(read_delim("I+T_RawPairedReads.txt", " ",escape_double = FALSE, col_names = FALSE, trim_ws = TRUE))
colnames(raw_reads) <- c("NODES", "Raw reads")
rownames(raw_reads) <- raw_reads$NODES
rownames(raw_reads)[rownames(raw_reads) == "TotalPairedReads"] <- "Total"
raw_reads$NODES <- NULL
raw_reads
####################################
# 1.1.2 Trimmed reads
####################################
trimmed_reads <- as.data.frame(read_delim("I+T_PairedTrimmedReads.txt", " ",escape_double = FALSE, col_names = FALSE, trim_ws = TRUE))
colnames(trimmed_reads) <- c("NODES", "Trimmed reads")
rownames(trimmed_reads) <- trimmed_reads$NODES
rownames(trimmed_reads)[rownames(trimmed_reads) == "TotalTrimmedPaired"] <- "Total"
trimmed_reads$NODES <- NULL
trimmed_reads
####################################
# 1.1.3 Contaminome trimmed reads
####################################
contaminome_removed_reads <- as.data.frame(read_delim("I+T_RawPairedNCoutReads.txt", " ",escape_double = FALSE, col_names = FALSE, trim_ws = TRUE))
colnames(contaminome_removed_reads) <- c("NODES", "Trimmed reads (−contaminome)")
rownames(contaminome_removed_reads) <- contaminome_removed_reads$NODES
rownames(contaminome_removed_reads)[rownames(contaminome_removed_reads) == "TotalNCoutPaired"] <- "Total"
contaminome_removed_reads$NODES <- NULL
contaminome_removed_reads
####################################
# 1.1.4 Quality-controlled reads
####################################
contaminome_human_removed_reads <- as.data.frame(read_delim("I+T_RawPairedHostoutReads.txt", " ",escape_double = FALSE, col_names = FALSE, trim_ws = TRUE))
colnames(contaminome_human_removed_reads) <- c("NODES", "Trimmed reads (−human genome)")
rownames(contaminome_human_removed_reads) <- contaminome_human_removed_reads$NODES
contaminome_human_removed_reads$NODES <- NULL
rownames(contaminome_human_removed_reads)[rownames(contaminome_human_removed_reads) == "TotalPairedHostoutReads"] <- "Total"
contaminome_human_removed_reads
####################################
# 1.1.5 Mapped quality-controlled reads (70% horizontal coverage)
####################################
reads_mapped <- as.data.frame(read_delim("I+T_MappedReads.txt", " ",escape_double = FALSE, col_names = FALSE, trim_ws = TRUE))
colnames(reads_mapped) <- c("NODES", "Trimmed reads (−mapped)")
rownames(reads_mapped) <- reads_mapped$NODES
reads_mapped$NODES <- NULL
rownames(reads_mapped)[rownames(reads_mapped) == "TotalPairedReads"] <- "Total"
reads_mapped
####################################
# 1.1.6 Viral reads
####################################
table(Mastertable$Final_viral)
Mastertable_viral <- Mastertable[Mastertable$Final_viral == "eukaryotic virus" | Mastertable$Final_viral == "phage",]
vector_1_1 <- (which(names(Mastertable_viral)== "Virsorter")-1)
Viral_Reads <- as.data.frame(colSums(Mastertable_viral[,1:vector_1_1]))
colnames(Viral_Reads) <- "Trimmed reads (-viral)`"
Viral_Reads <- Viral_Reads %>% rbind(sum(Viral_Reads$`Trimmed reads (-viral)`))
rownames(Viral_Reads)[rownames(Viral_Reads) == 305] <- "Total"
Viral_Reads
####################################
# 1.2 Merge 
####################################
COMBO_READ_1 <- transform(merge(raw_reads,trimmed_reads, by = 0, all = T), row.names = Row.names, Row.names = NULL)
COMBO_READ_2 <- transform(merge(COMBO_READ_1,contaminome_removed_reads, by = 0, all = T), row.names = Row.names, Row.names = NULL)
COMBO_READ_3 <- transform(merge(COMBO_READ_2,contaminome_human_removed_reads, by = 0, all = T), row.names = Row.names, Row.names = NULL)
COMBO_READ_4 <- transform(merge(COMBO_READ_3,reads_mapped, by = 0, all = T), row.names = Row.names, Row.names = NULL)
Mastertable_reads <- transform(merge(COMBO_READ_4,Viral_Reads, by = 0, all = T), row.names = Row.names, Row.names = NULL)
colnames(Mastertable_reads) <- c("Reads (-raw)", "Reads (-trimmed)", "Reads (-contaminome)", "Reads (-human genome)", "Reads (-mapped)","Reads (-viral)")
Mastertable_reads[is.na(Mastertable_reads)] <- 0
# View(Mastertable_reads)
####################################
# 1.3 Create Excel file 
####################################
setwd("/Users/daan/Desktop/Bioinformatics/Analysis/FMT/Data/2_Phages")
getwd()
dir()

Mastertable_reads[is.na(Mastertable_reads)] <- 0
WriteXLS("Mastertable_reads", "reads_summary.xlsx", row.names = TRUE)
#View(Mastertable_reads)
####################################
# 1.4 Extract reads per sample and total reads
####################################
# 1.4.1 Extract
####################################
Mastertable_reads_total <- Mastertable_reads[grep ("Total", rownames(Mastertable_reads)),]
Mastertable_reads_individual <- Mastertable_reads[grep ("Total", rownames(Mastertable_reads), invert=TRUE),]
#View(Mastertable_reads_total)
#View(Mastertable_reads_individual)
####################################
# 1.4.2 Create tables
####################################
Mastertable_reads_total <-  t(Mastertable_reads_total)
Mastertable_reads_total <- melt(Mastertable_reads_total)
colnames(Mastertable_reads_total) <- c("groups", "sample", "reads")
#View(Mastertable_reads_total)
Mastertable_reads_individual <-  t(Mastertable_reads_individual)
Mastertable_reads_individual <- melt(Mastertable_reads_individual)
colnames(Mastertable_reads_individual) <- c("groups", "sample", "reads")
####################################
# 1.4.3 Additional data
####################################
# Add reads per million and percentages that are left
Mastertable_reads_total$M_reads <- (Mastertable_reads_total$reads/10^6)
Mastertable_reads_total$percentage_left <- ((Mastertable_reads_total$M_reads/Mastertable_reads_total$M_reads[1])*100)

Mastertable_reads_individual$M_reads <- (Mastertable_reads_individual$reads/10^6)
# View(Mastertable_reads_total)
#View(Mastertable_reads_individual)
####################################
# 1.5 Figures
####################################
# 1.5.1 Barplot
####################################
# The colors are now manually done and too much. Adapt this for 5 groups I thought it was.
# Use colorpanels online

Mastertable_reads_total$Reads <- "Reads"
View(Mastertable_reads_total)

3155.751/4715.545
Reads_total <- ggplot(Mastertable_reads_total) +
 aes(x = groups, fill = groups, weight = M_reads) +
 geom_bar(position = "dodge", colour = "black") +
 labs(y = "Reads (M)") +
 scale_fill_manual(values = c("#d53e4f","#fc8d59","#fee08b","#ffffbf","#e6f598","#99d594","#3288bd")) +
 theme_bw() +
 theme(axis.title.y = element_text(size = 12L, 
 face = "bold")) +
 facet_wrap(vars(Reads), scales = "free_y") +
 theme(panel.grid.major = element_line(colour = "gray96", linetype = "dashed"), panel.grid.minor = element_line(colour = "gray96", linetype = "dashed"), axis.text.x = element_text(size = 8), axis.ticks.x=element_blank()) +labs(x = NULL, fill = NULL)+labs(y = "reads (M)") +
 scale_y_continuous(expand = c(0.03,0)) +
 theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
          
Reads_total
#################################### 
# 1.5.2 Pie chart: reads
#################################### 
table(Mastertable$Final_superkingdom)
Mastertable_unrarefied <- Mastertable
Mastertable_global_rarefied_categories <- Mastertable_unrarefied
Mastertable_global_rarefied_categories$Final_superkingdom[Mastertable_global_rarefied_categories$Final_superkingdom=="Archaea" | Mastertable_global_rarefied_categories$Final_superkingdom=="Dark matter" | Mastertable_global_rarefied_categories$Final_superkingdom=="Eukaryota"] <- "Other"
table(Mastertable_global_rarefied_categories$Final_superkingdom)
sort(colSums(Mastertable_global_rarefied_categories[1:304]))
names <- colnames(Mastertable_global_rarefied_categories[1:304][,colSums(Mastertable_global_rarefied_categories[1:304])!=0])
names # 100 samples with viral reads found of 144 samples or 69.4%

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

Mastertable_global_rarefied_categories$names <- rownames(Mastertable_global_rarefied_categories)
Mastertable_global_rarefied_categories$names <- factor(Mastertable_global_rarefied_categories$names, levels = rev(c("Unannotated", "Other", "Bacteria", "Viruses")))
#Mastertable_global_rarefied_categories$names <- factor(Mastertable_global_rarefied_categories$names, levels = rev(c("eukaryotic virus", "phage")))
Mastertable_global_rarefied_categories$relative_adapted <- Mastertable_global_rarefied_categories$relative
View(Mastertable_global_rarefied_categories)
Mastertable_global_rarefied_categories$relative_adapted <- c("0.7","99.3")
Mastertable_global_rarefied_categories$relative_adapted <- as.numeric(Mastertable_global_rarefied_categories$relative_adapted)
#view(Mastertable_global_rarefied_categories)

# Contigs
table(Mastertable_unrarefied$Final_superkingdom) # 39+24165+358+1130+810+2224
Mastertable_global_rarefied_categories$relative_adapted[Mastertable_global_rarefied_categories$relative_adapted==2.9] <- 2.8
Mastertable_global_rarefied_categories$relative_adapted[Mastertable_global_rarefied_categories$relative_adapted==4.6] <- 5.3
Mastertable_global_rarefied_categories$relative_adapted[Mastertable_global_rarefied_categories$relative_adapted==24.8] <- 84.1
Mastertable_global_rarefied_categories$relative_adapted[Mastertable_global_rarefied_categories$relative_adapted==67.8] <- 7.7
Mastertable_global_rarefied_categories$relative <- Mastertable_global_rarefied_categories$relative_adapted
View(Mastertable_global_rarefied_categories)

# (ii) Add label position (~ perhaps usefull)
Mastertable_global_rarefied_categories <- Mastertable_global_rarefied_categories %>%
  arrange(desc(names)) %>%
  mutate(lab.ypos = cumsum(relative) - 0.5*relative)

# (iii) Pie-chart (ever
mycols <- c("#377fb9", "#e41d1f","#4faf4b","#94509a")
#View(Mastertable_global_rarefied_categories)

Mastertable_global_rarefied_categories$relative_adapted <- Mastertable_global_rarefied_categories$relative
Mastertable_global_rarefied_categories$relative <- Mastertable_global_rarefied_categories$relative_adapted

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

# Contigs
# Total: 28726 NR-contigs
((39+358+1130)/(28726))*100 # 5.3% Other
(24165/(28726))*100 # 84.1% bacterial
(2224/(28726))*100 # 7.7% Viral
(810/(28726))*100 # 2.8% Unannoated

# Reads 
# Total: 4.667.752.674 billion reads
(213604973/(4667752674))*100 # 4.6% Other
(1156132135/(4667752674))*100 # 24.8% bacterial
(3164628531/(4667752674))*100 # 67.8% Viral
(133387035/(4667752674))*100 # 2.9% Unannoated

View(Mastertable_phage)
#################################### 
# 2. Contigs
####################################
# 2.1 Insert files
####################################
setwd("/Users/daan/Desktop//Transfer/Read_Contig_count")
dir()
####################################
# 2.1.1 MetaSPAdes contigs
####################################
Metaspades_contigs <- as.data.frame(read_delim("NumberOfAssembledContigs.txt", " ",escape_double = FALSE, col_names = FALSE, trim_ws = TRUE))
colnames(Metaspades_contigs) <- c("NODES", "contigs (-MetaSPAdes)")
rownames(Metaspades_contigs) <- Metaspades_contigs$NODES
rownames(Metaspades_contigs)[rownames(Metaspades_contigs) == "TotalOfMetaSpadesContigs"] <- "Total"
Metaspades_contigs$NODES <- NULL
####################################
# 2.1.2 Filtered 1kb contigs
####################################
filtered_contigs <- as.data.frame(read_delim("NumberOf1kbContigs.txt", " ",escape_double = FALSE, col_names = FALSE, trim_ws = TRUE))
colnames(filtered_contigs) <- c("NODES", "contigs (-1kb)")
rownames(filtered_contigs) <- filtered_contigs$NODES
rownames(filtered_contigs)[rownames(filtered_contigs) == "TotalOf1kbContigs"] <- "Total"
filtered_contigs$NODES <- NULL
filtered_contigs
####################################
# 2.1.3 Clustered abundant contigs (70% HC and 10X VC)
####################################
# 2.1.3.1 Total NR-contigs
####################################
NR_contigs <- Mastertable
NR_contigs <- cbind(NR_contigs[,1:vector_1_1])
NR_contigs <- as.data.frame(nrow(NR_contigs))
colnames(NR_contigs) <- "contigs (-NR)"
rownames(NR_contigs)[rownames(NR_contigs) == "1"] <- "Total"
####################################
# 2.1.3.2 Individual abundant NR-contigs
####################################
NR_contigs_individual <- Mastertable
NR_contigs_individual <- cbind(NR_contigs_individual[,1:vector_1_1])
NR_ones_individual <- vector()  

for (l in 1: length(names)){
  NR_ones_individual[l] <- print(NROW(NR_contigs_individual[NR_contigs_individual[,l] > 0, l])) 
}

NR_ones_individual <- as.data.frame(NR_ones_individual, row.names = names)
colnames(NR_ones_individual) <- "contigs (-NR)"
#View(NR_ones_individual)
####################################
# 2.1.4 Viral 
####################################
# 2.1.4.1 Total viral NR-contigs
####################################
Mastertable_viral_contigs <- Mastertable[Mastertable$Final_viral == "eukaryotic virus" | Mastertable$Final_viral == "phage",]
Mastertable_viral_contigs_ab_table <- cbind(Mastertable_viral_contigs[,1:vector_1_1])
viral_contigs <- as.data.frame(nrow(Mastertable_viral_contigs_ab_table))
colnames(viral_contigs) <- "contigs (-viral)"
rownames(viral_contigs)[rownames(viral_contigs) == "1"] <- "Total"
#View(viral_contigs)
####################################
# 2.1.4.2 Individual viral NR-contigs
####################################
Mastertable_viral_contigs <- Mastertable[Mastertable$Final_viral == "eukaryotic virus" | Mastertable$Final_viral == "phage",]
Mastertable_viral_contigs_ab_table <- cbind(Mastertable_viral_contigs[,1:vector_1_1])
viral_contigs_individual <- vector()  

for (l in 1: length(names)){
  viral_contigs_individual[l] <- print(NROW(Mastertable_viral_contigs_ab_table[Mastertable_viral_contigs_ab_table[,l] > 0, l])) 
}

viral_contigs_individual <- as.data.frame(viral_contigs_individual, row.names = names)
colnames(viral_contigs_individual) <- "contigs (-viral)"
#View(viral_contigs_individual)
####################################
# 2.2 Merge 
####################################
# 2.2.1 Merge dataframes
####################################
# 2.2.1.1 Total samples
####################################
COMBO_CONTIG_1 <- transform(merge(Metaspades_contigs,filtered_contigs, by = 0, all = T), row.names = Row.names, Row.names = NULL)
COMBO_CONTIG_2 <- transform(merge(COMBO_CONTIG_1,NR_contigs, by = 0, all = T), row.names = Row.names, Row.names = NULL)
Mastertable_contigs <- transform(merge(COMBO_CONTIG_2,viral_contigs, by = 0, all = T), row.names = Row.names, Row.names = NULL)
colnames(Mastertable_contigs) <- c("contigs (-MetaSPAdes)", "contigs (-1kb)", "contigs (-NR)", "contigs (-viral)")
####################################
# 2.2.1.2 Individual samples
####################################
#COMBO_CONTIG_1 <- transform(merge(Metaspades_contigs,filtered_contigs, by = 0, all = T), row.names = Row.names, Row.names = NULL)
##COMBO_CONTIG_2 <- transform(merge(COMBO_CONTIG_1,NR_ones_individual, by = 0, all = T), row.names = Row.names, Row.names = NULL)
#COMBO_CONTIG_3 <- transform(merge(COMBO_CONTIG_2,abundant_ones_individual, by = 0, all = T), row.names = Row.names, Row.names = NULL)
#Mastertable_CONTIGS <- transform(merge(COMBO_CONTIG_3,viral_contigs_individual, by = 0, all = T), row.names = Row.names, Row.names = NULL)
#colnames(Mastertable_CONTIGS) <- c("contigs (-MetaSPAdes)", "contigs (-1kb)", "contigs (-NR)", "contigs (-abundant)", "contigs (-viral)")
#Mastertable_CONTIGS[is.na(Mastertable_CONTIGS)] <- 0
#View(Mastertable_CONTIGS)
####################################
# 2.2.2 Extract group totals
####################################
# 2.2.2.1 Total samples
####################################
Mastertable_contigs_total <- Mastertable_contigs[grep ("Total", rownames(Mastertable_contigs)),]
#View(Mastertable_contigs_total)
####################################
# 2.2.2.2 Individual samples
####################################
#Mastertable_CONTIGS <- as.data.frame(Mastertable_CONTIGS[grep("Total", rownames(Mastertable_CONTIGS), invert=TRUE),])
#View(Mastertable_CONTIGS)
####################################
# 2.2.3 Transverse rows & melt
####################################
# 2.2.3.1 Total samples
####################################
Mastertable_contigs_total <-  t(Mastertable_contigs_total)
Mastertable_contigs_total <- melt(Mastertable_contigs_total)
colnames(Mastertable_contigs_total) <- c("groups", "sample", "contigs")
View(Mastertable_contigs_total)
####################################
# 2.2.3.2 Individual samples
####################################
#Mastertable_CONTIGS <-  t(Mastertable_CONTIGS)
#Mastertable_CONTIGS <- melt(Mastertable_CONTIGS)
#colnames(Mastertable_CONTIGS) <- c("groups", "sample", "contigs")
#View(Mastertable_CONTIGS)
####################################
# 2.2.4 Additional data
####################################
# 2.2.4.1 Total samples
####################################
Mastertable_contigs_total$contigs_thousand <- Mastertable_contigs_total$contigs/1000
Mastertable_contigs_total$log_10 <- log10(Mastertable_contigs_total$contigs)
Mastertable_contigs_total$relative <- ((Mastertable_contigs_total$contigs/Mastertable_contigs_total$contigs[1])*100)
Mastertable_contigs_total$contigs_name <- "contigs"
#view(Mastertable_contigs_total)
####################################
# 2.2.4.2 Individual samples
####################################
#Mastertable_CONTIGS$thousand_contigs <- (Mastertable_CONTIGS$contigs/1000)
#Mastertable_CONTIGS$log_10 <- log10(Mastertable_CONTIGS$contigs)
#Mastertable_CONTIGS$log_10[Mastertable_CONTIGS$log_10 == -Inf] <- 0
#View(Mastertable_CONTIGS)
####################################
# 2.3 Figures
####################################
# 2.3.1 Barplot
####################################
# 2.3.1.1 All samples (log10)
####################################
Mastertable_contigs_total$groups <- factor(Mastertable_contigs_total$groups, levels = c("contigs (-MetaSPAdes)","contigs (-1kb)","contigs (-NR)","contigs (-viral)"))
View(Mastertable_contigs_total)

ggplot(Mastertable_contigs_total) +
 aes(x = groups, fill = groups, weight = log_10) +
  geom_bar(position = "dodge", colour = "black") +
 labs(y = "log10(contigs)") +
 theme_bw() +
 scale_fill_manual(values = c("#d53e4f","#fc8d59","#fee08b","#99d594")) +
#scale_fill_manual(values = c("#d53e4f","#fc8d59","#fee08b","#e6f598","#99d594")) +
 theme(axis.title.y = element_text(size = 12L, face = "bold")) +
 facet_wrap(vars(contigs_name), scales = "free_y") +
 theme(panel.grid.major = element_line(colour = "gray96", linetype = "dashed"), panel.grid.minor = element_line(colour = "gray96", linetype = "dashed"), axis.text.x = element_text(size = 8), axis.ticks.x=element_blank()) +
  labs(x = NULL, fill = NULL) +
 scale_y_continuous(expand = c(0.03,0)) +
 theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

Mastertable_contigs_total_plot

(2216/22552127)*100
####################################
# 2.3.1.3 Individual samples
####################################
#ggplot(Mastertable_CONTIGS, aes(x = sample, y = thousand_contigs, fill = groups)) +
  geom_bar(mapping = aes(x = sample, y = thousand_contigs, fill = groups), stat = "identity",  position = position_dodge(width = 0), width = 3, colour = "black") +
  xlab("") +
  ylab("thousand contigs") +
  ggtitle("") +
  theme(legend.title = element_blank(), legend.text = element_text(size = 7)) + 
  theme(plot.title = element_text(hjust = 0.5, vjust = -2)) +
  theme(axis.title.y = element_text(size=9)) +
  scale_fill_manual(values = c("#d53e4f","#fc8d59","#fee08b","#99d594","#3288bd")) + 
  theme(axis.text.x = element_text(size = 5, angle = 45, vjust = 1 , hjust = 1), axis.title = element_text(size = 5)) + 
  theme(panel.background = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(color = "black", linetype = "solid", size = 0.3)) +
  theme(axis.title.x = element_text(vjust = -1)) +
  theme(axis.text = element_text(color = "black")) +
  scale_y_continuous(expand = c(0,0)) 
# Save in 4_14 format or even larger
####################################
# 2.3.1.4 Individual samples log 10
####################################
esquisser(Mastertable_CONTIGS)

ggplot(Mastertable_CONTIGS) +
  aes(x = sample, fill = groups, weight = log_10) +
  geom_bar() +
  # scale_fill_hue(direction = 1) +
  scale_fill_manual(values = c("#d53e4f","#fc8d59","#fee08b","#ffffbf","#e6f598","#99d594","#3288bd")) +
  theme_classic() +
  labs(y = "log10(contigs)", subtitle = "contig per sample") +
  theme(plot.subtitle = element_text(size = 12L, face = "bold"), axis.title.y = element_text(size = 10L)) + 
  theme(axis.text.x = element_blank(), axis.ticks.x=element_blank()) +
  coord_cartesian(ylim=c(0,18), clip="off") 
  

    theme(legend.position = "none") +labs(x = NULL, fill = NULL) +
  theme(axis.text.x = element_blank(), axis.ticks.x=element_blank()) +
  geom_hline(yintercept=1, color="darkgreen", size=0.5) +
  scale_y_continuous(breaks = c(0,1,5,10,15,20,25,30)) +
  expand_limits(x= c(0, 440))
####################################
# 3. Completeness
####################################
setwd("/Users/daan/Desktop/Bioinformatics/pre-processing/HS project/virome/Data/Skin_test_protocol/Transfer/Contigs")
dir()
completeness <- as.data.frame(read_delim("CheckV_completeness_contigs.txt", "\t",escape_double = FALSE, col_names = FALSE, trim_ws = TRUE))
colnames(completeness) <- c("Completeness", "value")
rownames(completeness) <- completeness$Completeness
rownames(completeness)[rownames(completeness) == "Total non-redundant contigs"] <- "Total"
completeness$Completeness <- NULL

ggplot(completeness) +
  aes(x = Completeness, fill = Completeness, weight = value) +
  geom_bar(position = "dodge", colour = "black") +
  labs(y = "total contig count") +
  scale_fill_manual(values = c("#d53e4f","#fc8d59","#fee08b","#ffffbf","#e6f598","#99d594","#3288bd")) +
  theme_bw() +
  theme(axis.title.y = element_text(size = 12L, face = "bold")) +
  theme(panel.grid.major = element_line(colour = "gray96", linetype = "dashed"), panel.grid.minor = element_line(colour = "gray96", linetype = "dashed"), axis.text.x = element_text(size = 8), axis.ticks.x=element_blank()) +labs(x = NULL, fill = NULL)+labs(y = "total contig count") +
  scale_y_continuous(expand = c(0.03,0)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
######### !!!!!!!!!!!!!!! #########
# SAVE Global Environment #
######### !!!!!!!!!!!!!!! #########