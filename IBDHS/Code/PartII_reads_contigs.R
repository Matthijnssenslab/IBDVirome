####################################
# SCRIPT 2: CREATE READ AND CONTIG FIGURE + GUIDANCE FIGURES FOR RAREFACTION
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
setwd("/Users/daan/Desktop/HS/Input_R/Reads")
dir()
####################################
# 1.1.1 Raw reads
####################################
raw_reads <- as.data.frame(read_delim("I+T_RawPairedReads.txt", " ",escape_double = FALSE, col_names = FALSE, trim_ws = TRUE))
colnames(raw_reads) <- c("NODES", "Raw reads")
rownames(raw_reads) <- raw_reads$NODES
#rownames(raw_reads)[rownames(raw_reads) == "TotalPairedReads"] <- "Total"
raw_reads$NODES <- NULL
####################################
# 1.1.2 Trimmed reads
####################################
trimmed_reads <- as.data.frame(read_delim("I+T_PairedTrimmedReads.txt", " ",escape_double = FALSE, col_names = FALSE, trim_ws = TRUE))
colnames(trimmed_reads) <- c("NODES", "Trimmed reads")
rownames(trimmed_reads) <- trimmed_reads$NODES
#rownames(trimmed_reads)[rownames(trimmed_reads) == "TotalTrimmed"] <- "Total"
trimmed_reads$NODES <- NULL
####################################
# 1.1.3 Contaminome trimmed reads
####################################
contaminome_removed_reads <- as.data.frame(read_delim("I+T_NCoutPairedReads.txt", " ",escape_double = FALSE, col_names = FALSE, trim_ws = TRUE))
colnames(contaminome_removed_reads) <- c("NODES", "Trimmed reads (−contaminome)")
rownames(contaminome_removed_reads) <- contaminome_removed_reads$NODES
#rownames(contaminome_removed_reads)[rownames(contaminome_removed_reads) == "TotalNCout"] <- "Total"
contaminome_removed_reads$NODES <- NULL
####################################
# 1.1.4 Quality-controlled reads
####################################
contaminome_human_removed_reads <- as.data.frame(read_delim("I+T_ContaminomeHumanOutPairedReads.txt", " ",escape_double = FALSE, col_names = FALSE, trim_ws = TRUE))
colnames(contaminome_human_removed_reads) <- c("NODES", "Trimmed reads (−human genome)")
rownames(contaminome_human_removed_reads) <- contaminome_human_removed_reads$NODES
contaminome_human_removed_reads$NODES <- NULL
####################################
# 1.1.5 Mapped quality-controlled reads (70% horizontal coverage)
####################################
reads_mapped <- as.data.frame(read_delim("I+T_Mapped_bwa2_NR_Reads_threshold.txt", "\t",escape_double = FALSE, col_names = FALSE, trim_ws = TRUE))
colnames(reads_mapped) <- c("NODES", "Trimmed reads (−mapped)")
rownames(reads_mapped) <- reads_mapped$NODES
reads_mapped$NODES <- NULL
rownames(reads_mapped)[rownames(reads_mapped) == "TotalAbundances"] <- "Total"
reads_mapped[is.na(reads_mapped)] <- 0
#View(reads_mapped)
####################################
# 1.1.7 Viral reads
####################################
# 1.1.7.1 All Contigs
Mastertable_contigs <- Mastertable
Mastertable_contigs$Final_viral[Mastertable_contigs$Phage_score1=="phage" & !Mastertable_contigs$Final_viral=="eukaryotic virus"] <- "prokaryotic virus"
table(Mastertable_contigs$Final_viral)
Mastertable_contigs_viral <- Mastertable_contigs[Mastertable_contigs$Final_viral == "eukaryotic virus" | Mastertable_contigs$Final_viral == "prokaryotic virus",]
vector_1_1 <- (which(names(Mastertable_contigs_viral)== "Virsorter")-1)
Viral_Reads_contigs <- as.data.frame(colSums(Mastertable_contigs_viral[,1:vector_1_1]))
colnames(Viral_Reads_contigs) <- "Trimmed reads (-viral)`"
Viral_Reads_contigs <- Viral_Reads_contigs %>% rbind(sum(Viral_Reads_contigs$`Trimmed reads (-viral)`))
rownames(Viral_Reads_contigs)[rownames(Viral_Reads_contigs) == 128] <- "Total"
table(Mastertable_contigs_viral$Final_viral)
#View(Viral_Reads_contigs)

# 1.1.7.2 All Families
Mastertable$Final_viral[Mastertable$Phages=="phage" & !Mastertable$Final_viral=="eukaryotic virus"] <- "prokaryotic virus"
table(Mastertable$Final_viral)
Mastertable_viral <- Mastertable[Mastertable$Final_viral == "eukaryotic virus" | Mastertable$Final_viral == "prokaryotic virus",]
vector_1_1 <- (which(names(Mastertable_viral)== "Virsorter")-1)
Viral_Reads <- as.data.frame(colSums(Mastertable_viral[,1:vector_1_1]))
colnames(Viral_Reads) <- "Trimmed reads (-viral)`"
Viral_Reads <- Viral_Reads %>% rbind(sum(Viral_Reads$`Trimmed reads (-viral)`))
rownames(Viral_Reads)[rownames(Viral_Reads) == 128] <- "Total"
#View(Viral_Reads)
####################################
# 1.2 Merge 
####################################
# 1.2.1. All contigs
COMBO_READ_1 <- transform(merge(raw_reads,trimmed_reads, by = 0, all = T), row.names = Row.names, Row.names = NULL)
COMBO_READ_2 <- transform(merge(COMBO_READ_1,contaminome_removed_reads, by = 0, all = T), row.names = Row.names, Row.names = NULL)
COMBO_READ_3 <- transform(merge(COMBO_READ_2,contaminome_human_removed_reads, by = 0, all = T), row.names = Row.names, Row.names = NULL)
COMBO_READ_4 <- transform(merge(COMBO_READ_3,reads_mapped, by = 0, all = T), row.names = Row.names, Row.names = NULL)
Mastertable_reads_contigs <- transform(merge(COMBO_READ_4,Viral_Reads_contigs, by = 0, all = T), row.names = Row.names, Row.names = NULL)
colnames(Mastertable_reads_contigs) <- c("Reads (-raw)", "Reads (-trimmed)", "Reads (-contaminome)", "Reads (-human genome)", "Reads (-mapped)","Reads (-viral)")
Mastertable_reads_contigs[is.na(Mastertable_reads_contigs)] <- 0
View(Mastertable_reads_contigs)

# 1.2.1. All Families
COMBO_READ_1 <- transform(merge(raw_reads,trimmed_reads, by = 0, all = T), row.names = Row.names, Row.names = NULL)
COMBO_READ_2 <- transform(merge(COMBO_READ_1,contaminome_removed_reads, by = 0, all = T), row.names = Row.names, Row.names = NULL)
COMBO_READ_3 <- transform(merge(COMBO_READ_2,contaminome_human_removed_reads, by = 0, all = T), row.names = Row.names, Row.names = NULL)
COMBO_READ_4 <- transform(merge(COMBO_READ_3,reads_mapped, by = 0, all = T), row.names = Row.names, Row.names = NULL)
Mastertable_reads <- transform(merge(COMBO_READ_4,Viral_Reads, by = 0, all = T), row.names = Row.names, Row.names = NULL)
colnames(Mastertable_reads) <- c("Reads (-raw)", "Reads (-trimmed)", "Reads (-contaminome)", "Reads (-human genome)", "Reads (-mapped)","Reads (-viral)")
Mastertable_reads[is.na(Mastertable_reads)] <- 0
#View(Mastertable_reads)
####################################
# 1.3 Create Excel file 
####################################
View(Mastertable)
# 1.3.1 All contigs
Mastertable_reads_contigs[is.na(Mastertable_reads_contigs)] <- 0
WriteXLS("Mastertable_reads_contigs", "reads_summary.xlsx", row.names = TRUE)

# 1.3.2 All families
Mastertable_reads[is.na(Mastertable_reads)] <- 0
WriteXLS("Mastertable_reads", "reads_summary.xlsx", row.names = TRUE)
####################################
# 1.4 Extract reads per sample and total reads
####################################
# 1.4.1 Extract
####################################
# 1.4.1 All contigs
Mastertable_reads_total_contigs  <- Mastertable_reads_contigs[grep ("Total", rownames(Mastertable_reads_contigs)),]
Mastertable_reads_individual_contigs <- Mastertable_reads_contigs[grep ("Total", rownames(Mastertable_reads_contigs), invert=TRUE),]

# 1.4.2 All families
Mastertable_reads_total <- Mastertable_reads[grep ("Total", rownames(Mastertable_reads)),]
Mastertable_reads_individual <- Mastertable_reads[grep ("Total", rownames(Mastertable_reads), invert=TRUE),]
####################################
# 1.4.2 Create tables
####################################
# 1.4.2.1 All contigs
Mastertable_reads_total_contigs <-  t(Mastertable_reads_total_contigs)
Mastertable_reads_total_contigs <- melt(Mastertable_reads_total_contigs)
colnames(Mastertable_reads_total_contigs) <- c("groups", "sample", "reads")

Mastertable_reads_individual_contigs <-  t(Mastertable_reads_individual_contigs)
Mastertable_reads_individual_contigs <- melt(Mastertable_reads_individual_contigs)
colnames(Mastertable_reads_individual_contigs) <- c("groups", "sample", "reads")

# 1.4.2.2 All families
Mastertable_reads_total <-  t(Mastertable_reads_total)
Mastertable_reads_total <- melt(Mastertable_reads_total)
colnames(Mastertable_reads_total) <- c("groups", "sample", "reads")

Mastertable_reads_individual <-  t(Mastertable_reads_individual)
Mastertable_reads_individual <- melt(Mastertable_reads_individual)
colnames(Mastertable_reads_individual) <- c("groups", "sample", "reads")
####################################
# 1.4.3 Additional data
####################################
# 1.4.3.1 All contigs
Mastertable_reads_total_contigs$M_reads <- (Mastertable_reads_total_contigs$reads/10^6)
Mastertable_reads_total_contigs$percentage_left <- ((Mastertable_reads_total_contigs$M_reads/Mastertable_reads_total_contigs$M_reads[1])*100)
Mastertable_reads_individual_contigs$M_reads <- (Mastertable_reads_individual_contigs$reads/10^6)

# 1.4.3.2 All families
Mastertable_reads_total$M_reads <- (Mastertable_reads_total$reads/10^6)
Mastertable_reads_total$percentage_left <- ((Mastertable_reads_total$M_reads/Mastertable_reads_total$M_reads[1])*100)
Mastertable_reads_individual$M_reads <- (Mastertable_reads_individual$reads/10^6)
####################################
# 1.5 Figures
####################################
# 1.5.1 Barplot
####################################
# 1.5.1.1 All samples contigs
####################################
Mastertable_reads_total_contigs$Reads <- "Reads"

Reads_total <- ggplot(Mastertable_reads_total_contigs) +
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

View(Mastertable_reads_total_contigs)

(176726685/440099377)*100
#################################### 
# 1.5.1.1 All samples family
####################################
Mastertable_reads_total$Reads <- "Reads"

Reads_total_family <- ggplot(Mastertable_reads_total) +
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

Reads_total_family

View(Mastertable_reads_total)
Mastertable_reds_total_contigs$M_reads
(294336627/440099377)*100
#################################### 
# 2. Pie chart viruses (phages/eukaryotic viruses)
#################################### 
# 2.1 All Samples contigs
#################################### 
table(Mastertable_contigs_viral$Final_viral)
Mastertable_unrarefied <- Mastertable_contigs_viral
Mastertable_global_rarefied_categories <- Mastertable_unrarefied
table(Mastertable_global_rarefied_categories$Final_viral)
sort(colSums(Mastertable_global_rarefied_categories[1:127]))
names <- colnames(Mastertable_global_rarefied_categories[1:127][,colSums(Mastertable_global_rarefied_categories[1:127])!=0])
names # 100 samples with viral reads found of 144 samples or 69.4%

# View(Mastertable)
Mastertable_global_rarefied_categories <- aggregate(. ~Final_viral, FUN = sum, data = Mastertable_global_rarefied_categories[,colnames(Mastertable_global_rarefied_categories) %in% names | colnames(Mastertable_global_rarefied_categories) == 'Final_viral'])
rownames(Mastertable_global_rarefied_categories) <- Mastertable_global_rarefied_categories$Final_viral
Mastertable_global_rarefied_categories$Final_viral <- NULL
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
Mastertable_global_rarefied_categories$names <- factor(Mastertable_global_rarefied_categories$names, levels = rev(c("eukaryotic virus", "prokaryotic virus")))
Mastertable_global_rarefied_categories$relative_adapted <- Mastertable_global_rarefied_categories$relative
Mastertable_global_rarefied_categories$relative_adapted <- c("26.4","73.6")
Mastertable_global_rarefied_categories$relative_adapted <- as.numeric(Mastertable_global_rarefied_categories$relative_adapted)
#view(Mastertable_global_rarefied_categories)

# (ii) Add label position (~ perhaps usefull)
Mastertable_global_rarefied_categories <- Mastertable_global_rarefied_categories %>%
  arrange(desc(names)) %>%
  mutate(lab.ypos = cumsum(relative) - 0.5*relative)

# (iii) Pie-chart (everything < 1% is not shown)
mycols <- c("#ef8a62", "#67a9cf")
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
#################################### 
# 2.2. All samples family
#################################### 
table(Mastertable$Final_viral)
Mastertable_unrarefied <- Mastertable[Mastertable$Final_viral=="eukaryotic virus" | Mastertable$Final_viral=="prokaryotic virus",]
Mastertable_global_rarefied_categories <- Mastertable_unrarefied
table(Mastertable_global_rarefied_categories$Final_viral)
sort(colSums(Mastertable_global_rarefied_categories[1:127]))
names <- colnames(Mastertable_global_rarefied_categories[1:127][,colSums(Mastertable_global_rarefied_categories[1:127])!=0])
names # 117 samples of 144 samples or 81.3%

# View(Mastertable)
Mastertable_global_rarefied_categories <- aggregate(. ~Final_viral, FUN = sum, data = Mastertable_global_rarefied_categories[,colnames(Mastertable_global_rarefied_categories) %in% names | colnames(Mastertable_global_rarefied_categories) == 'Final_viral'])
rownames(Mastertable_global_rarefied_categories) <- Mastertable_global_rarefied_categories$Final_viral
Mastertable_global_rarefied_categories$Final_viral <- NULL

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
Mastertable_global_rarefied_categories$names <- factor(Mastertable_global_rarefied_categories$names, levels = rev(c("eukaryotic virus", "prokaryotic virus")))
Mastertable_global_rarefied_categories$relative_adapted <- Mastertable_global_rarefied_categories$relative
Mastertable_global_rarefied_categories$relative_adapted <- c("15.9","84.1")
Mastertable_global_rarefied_categories$relative_adapted <- as.numeric(Mastertable_global_rarefied_categories$relative_adapted)
#view(Mastertable_global_rarefied_categories)

# (ii) Add label position (~ perhaps usefull)
Mastertable_global_rarefied_categories <- Mastertable_global_rarefied_categories %>%
  arrange(desc(names)) %>%
  mutate(lab.ypos = cumsum(relative) - 0.5*relative)

# (iii) Pie-chart (everything < 1% is not shown)
mycols <- c("#ef8a62", "#67a9cf")
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
#################################### 

#################################### 

######### !!!!!!!!!!!!!!! #########
# SAVE Global Environment #
######### !!!!!!!!!!!!!!! #########