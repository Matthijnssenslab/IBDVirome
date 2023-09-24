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
# 2. Basic characteristics
####################################
# 2.1 Plot scatterplot
####################################
Mastertable_phage_figure <- Mastertable_phage[c(320,321,327:330,331,374:375,382)]

# Total figure: colored by class
Mastertable_phage_figure %>%
 filter(!(Phyla_phylum %in% c("Euryarchaeota", "Thermodesulfobacteriota","Mycoplasmota", "Cyanobacteriota", "Mycoplasmatota"))) %>%
 filter(length >= 0L & length <= 200) %>%
 filter(!(Final_order %in% "Durnavirales")) %>%
 ggplot() +
 aes(x = `viral genes`, y = length, colour = Phyla_phylum, size = Totalnumberofreads) +
 geom_point(shape = "circle") +
 scale_color_hue(direction = 1) +
 theme_bw() +
 xlim(0, 64)

# Petitvirales
Mastertable_phage_figure %>%
 filter(!(Phyla_phylum %in% c("Euryarchaeota", "Thermodesulfobacteriota", "Mycoplasmota", "Cyanobacteriota", "Mycoplasmatota"))) %>%
 filter(length >= 0L & length <= 10) %>%
  filter((Final_order %in% "Petitvirales")) %>%
  ggplot() +
 aes(x = `viral genes`, y = length, colour = Phyla_phylum, size = Totalnumberofreads) +
 geom_point(shape = "circle") +
 scale_color_hue(direction = 1) +
 theme_bw() +
 facet_wrap(vars(Final_order), 
 scales = "free") +
  xlim(0,10)

# Crassvirales
Mastertable_phage_figure %>%
  filter(!(Phyla_phylum %in% c("Euryarchaeota", "Thermodesulfobacteriota", "Mycoplasmota", "Cyanobacteriota", "Mycoplasmatota"))) %>%
  filter(length >= 0L & length <= 200) %>%
  filter((Final_order %in% "Crassvirales")) %>%
  ggplot() +
  aes(x = `viral genes`, y = length, colour = Phyla_phylum, size = Totalnumberofreads) +
  geom_point(shape = "circle") +
  scale_color_hue(direction = 1) +
  theme_bw() +
  facet_wrap(vars(Final_order),scales = "free") +
  ylim(40,200)

# Unclassified Caudoviricetes
Mastertable_phage_figure %>%
  filter(!(Phyla_phylum %in% c("Euryarchaeota", "Thermodesulfobacteriota", "Mycoplasmota", "Cyanobacteriota", "Mycoplasmatota"))) %>%
  filter(length >= 0L & length <= 200) %>%
  filter((Final_order %in% "Unclassified Caudoviricetes")) %>%
  ggplot() +
  aes(x = `viral genes`, y = length, colour = Phyla_phylum, size = Totalnumberofreads) +
  geom_point(shape = "circle") +
  scale_color_hue(direction = 1) +
  theme_bw() +
  facet_wrap(vars(Final_order),scales = "free") +
  xlim(0,60) +
  ylim(10,200)

# Tubulavirales
Mastertable_phage_figure %>%
  filter(!(Phyla_phylum %in% c("Euryarchaeota", "Thermodesulfobacteriota", "Mycoplasmota", "Cyanobacteriota", "Mycoplasmatota"))) %>%
  filter(length >= 0L & length <= 10) %>%
  filter((Final_order %in% "Tubulavirales")) %>%
  ggplot() +
  aes(x = `viral genes`, y = length, colour = Phyla_phylum, size = Totalnumberofreads) +
  geom_point(shape = "circle") +
  scale_color_hue(direction = 1) +
  theme_bw() +
  facet_wrap(vars(Final_order),scales = "free") 

# Durnavirales
Mastertable_phage_figure %>%
  filter(!(Phyla_phylum %in% c("Euryarchaeota", "Thermodesulfobacteriota", "Mycoplasmota", "Cyanobacteriota", "Mycoplasmatota"))) %>%
  filter(length >= 0L & length <= 10) %>%
  filter((Final_order %in% "Durnavirales")) %>%
  ggplot() +
  aes(x = `viral genes`, y = length, colour = Phyla_phylum, size = Totalnumberofreads) +
  geom_point(shape = "circle") +
  scale_color_hue(direction = 1) +
  theme_bw() +
  facet_wrap(vars(Final_order),scales = "free") 
####################################
# 2.2 Plot barplot per viral order
####################################
Mastertable_phage_figure %>%
 filter(Phyla_phylum %in% c("Chlamydiota", "Bacillota", "Bacteroidota", "0", 
"Pseudomonadota", "Actinomycetota", "Bdellovibrionota")) %>%
 ggplot() +
 aes(x = Phyla_phylum) +
 geom_bar(fill = "#112446") +
 theme_bw()

Mastertable_phage_figure %>%
 filter(!(Phyla_phylum %in% c("Mycoplasmatota", "Cyanobacteriota", "Euryarchaeota", 
"Thermodesulfobacteriota"))) %>%
 filter(Final_order %in% "Unclassified Caudoviricetes") %>%
 ggplot() +
 aes(x = Final_order, fill = Phyla_phylum) +
 geom_bar(position = "dodge") +
 scale_fill_hue(direction = 1) +
 theme_bw()

Mastertable_phage_figure %>%
  filter(!(Phyla_phylum %in% c("Mycoplasmatota", "Cyanobacteriota", "Euryarchaeota", 
                               "Thermodesulfobacteriota"))) %>%
  filter(Final_order %in% "Petitvirales") %>%
  ggplot() +
  aes(x = Final_order, fill = Phyla_phylum) +
  geom_bar(position = "dodge") +
  scale_fill_hue(direction = 1) +
  theme_bw()

Mastertable_phage_figure %>%
  filter(!(Phyla_phylum %in% c("Mycoplasmatota", "Cyanobacteriota", "Euryarchaeota", 
                               "Thermodesulfobacteriota"))) %>%
  filter(Final_order %in% "Crassvirales") %>%
  ggplot() +
  aes(x = Final_order, fill = Phyla_phylum) +
  geom_bar(position = "dodge") +
  scale_fill_hue(direction = 1) +
  theme_bw()

Mastertable_phage_figure %>%
  filter(!(Phyla_phylum %in% c("Mycoplasmatota", "Cyanobacteriota", "Euryarchaeota", 
                               "Thermodesulfobacteriota"))) %>%
  filter(Final_order %in% "Tubulavirales") %>%
  ggplot() +
  aes(x = Final_order, fill = Phyla_phylum) +
  geom_bar(position = "dodge") +
  scale_fill_hue(direction = 1) +
  theme_bw()

Mastertable_phage_figure %>%
  filter(!(Phyla_phylum %in% c("Mycoplasmatota", "Cyanobacteriota", "Euryarchaeota", 
                               "Thermodesulfobacteriota"))) %>%
  filter(Final_order %in% "Durnavirales") %>%
  ggplot() +
  aes(x = Final_order, fill = Phyla_phylum) +
  geom_bar(position = "dodge") +
  scale_fill_hue(direction = 1) +
  theme_bw()

View(Mastertable_phage[Mastertable_phage$Final_order=="Unclassified Caudoviricetes" & Mastertable_phage$Phyla_family=="Enterobacteriaceae",])
####################################
# 2.3 Split by donor/patient
####################################
# 2.3.1  Upload metadata
#################################### 
setwd("/Users/daan/Desktop/Bioinformatics/pre-processing/FMT_study/Metadata")
dir()

metadata <- read_excel("samples_FMT_project_Final_DJ_R.xlsx")
metadata <- as.data.frame(metadata)
metadata$`cellcounts cell/gram` <- as.numeric(metadata$`cellcounts cell/gram`)
metadata$`cellcounts cells/mL` <- as.numeric(metadata$`cellcounts cells/mL`)
metadata$CRP <- as.numeric(metadata$CRP)
metadata$Hemoglobin <- as.numeric(metadata$Hemoglobin)
metadata$Age <- as.numeric(metadata$Age)
metadata$Weight <- as.numeric(metadata$Weight)
metadata$Length <- as.numeric(metadata$Length)
metadata$`Donation time` <- as.numeric(metadata$`Donation time`)
metadata$`Disease duration (y)` <- as.numeric(metadata$`Disease duration (y)`)
metadata$`Type of donor`[metadata$`Type of donor`=="open-label superdonor"] <- "superdonor"
metadata$`Disease status`[metadata$`Disease status`=="donor" & metadata$`Type of donor`=="autologous donor"] <- "patient"
rownames(metadata) <- metadata$unique_identifier
colnames(metadata)
samples <- metadata
#################################### 
# 2.3.2 Donors
#################################### 
table(metadata$`Disease status`)
samples_donor <- rownames(metadata[metadata$`Disease status`=="donor",])
Mastertable_phage_met <- Mastertable_phage[,305:385]
Mastertable_phage_ab <- Mastertable_phage[,1:304]
Mastertable_phage_donors <- Mastertable_phage_ab %>% select(samples_donor,)
Mastertable_phage_donors$reads <- rowSums(Mastertable_phage_donors)
Mastertable_phage_donors1 <- Mastertable_phage_donors[!Mastertable_phage_donors$reads==0,]
Mastertable_phage_donors_complete <- as.data.frame(merge(Mastertable_phage_donors1,Mastertable_phage_met,by=0,all=F))
rownames(Mastertable_phage_donors_complete) <- Mastertable_phage_donors_complete$Row.names
Mastertable_phage_donors_complete$Row.names <- NULL

colnames(Mastertable_phage_donors_complete[320])
Mastertable_phage_figure <- Mastertable_phage_donors_complete[c(72,71,78:81,82,125:126,133)]

ggplot(Mastertable_phage_figure) +
 aes(x = `viral genes`, y = length, colour = Phyla_phylum, size = Totalnumberofreads) +
 geom_point(shape = "circle") +
 scale_color_hue(direction = 1) +
 theme_bw() +
 ylim(0, 200)

# Petitvirales
Mastertable_phage_figure %>%
  filter(!(Phyla_phylum %in% c("Euryarchaeota", "Thermodesulfobacteriota", "Mycoplasmota", "Cyanobacteriota", "Mycoplasmatota"))) %>%
  filter(length >= 0L & length <= 10) %>%
  filter((Final_order %in% "Petitvirales")) %>%
  ggplot() +
  aes(x = `viral genes`, y = length, colour = Phyla_phylum, size = Totalnumberofreads) +
  geom_point(shape = "circle") +
  scale_color_hue(direction = 1) +
  theme_bw() +
  facet_wrap(vars(Final_order), 
             scales = "free") +
  xlim(0,10)

# Unclassified Caudoviricetes
Mastertable_phage_figure %>%
  filter(!(Phyla_phylum %in% c("Euryarchaeota", "Thermodesulfobacteriota", "Mycoplasmota", "Cyanobacteriota", "Mycoplasmatota"))) %>%
  filter(length >= 0L & length <= 200) %>%
  filter((Final_order %in% "Unclassified Caudoviricetes")) %>%
  ggplot() +
  aes(x = `viral genes`, y = length, colour = Phyla_phylum, size = Totalnumberofreads) +
  geom_point(shape = "circle") +
  scale_color_hue(direction = 1) +
  theme_bw() +
  facet_wrap(vars(Final_order),scales = "free") +
  xlim(0,60) +
  ylim(10,100)

#  Crassvirales
Mastertable_phage_figure %>%
  filter(!(Phyla_phylum %in% c("Euryarchaeota", "Thermodesulfobacteriota", "Mycoplasmota", "Cyanobacteriota", "Mycoplasmatota"))) %>%
  filter(length >= 0L & length <= 200) %>%
  filter((Final_order %in% "Crassvirales")) %>%
  ggplot() +
  aes(x = `viral genes`, y = length, colour = Phyla_phylum, size = Totalnumberofreads) +
  geom_point(shape = "circle") +
  scale_color_hue(direction = 1) +
  theme_bw() +
  facet_wrap(vars(Final_order),scales = "free") +
  xlim(0,60) +
  ylim(10,200)
####################################
# 2.3.3 Patients
####################################
metadata_P <- metadata[!metadata$unique_identifier=="P57_w4" & !metadata$unique_identifier=="P57_w8" & !metadata$unique_identifier=="P73_w12" & !metadata$unique_identifier=="P22_w0" & !metadata$unique_identifier=="P52_w0",]
metadata_P1 <- metadata_P[!metadata_P$unique_identifier=="P53_w4" & !metadata_P$unique_identifier=="P83_w8" & !metadata_P$unique_identifier=="P23_w0OL" & !metadata_P$unique_identifier=="P23_w8OL" & !metadata_P$unique_identifier=="P52_w0OL",]
metadata_P2 <- metadata_P1[!metadata_P1$unique_identifier=="P53_w0OL" & !metadata_P1$unique_identifier=="P61_w0OL" & !metadata_P1$unique_identifier=="P70_w0OL" & !metadata_P1$unique_identifier=="P83_w0OL" & !metadata_P1$unique_identifier=="AD11_w1",]
metadata_P3 <- metadata_P2[!metadata_P2$unique_identifier=="AD22_w1" & !metadata_P2$unique_identifier=="AD53_w1" & !metadata_P2$unique_identifier=="AD70_w2",]

samples_patient <- rownames(metadata_P3[metadata_P3$`Disease status`=="patient",])
Mastertable_phage_met <- Mastertable_phage[,305:385]
colnames(Mastertable_phage_ab)
Mastertable_phage_ab <- Mastertable_phage[,1:304]
Mastertable_phage_patient <- Mastertable_phage_ab %>% select(samples_patient,)

Mastertable_phage_patient$reads <- rowSums(Mastertable_phage_patient)
Mastertable_phage_donors1 <- Mastertable_phage_patient[!Mastertable_phage_patient$reads==0,]
Mastertable_phage_donors_complete <- merge(Mastertable_phage_donors1,Mastertable_phage_met,by=0,all=F)
rownames(Mastertable_phage_donors_complete) <- Mastertable_phage_donors_complete$Row.names
Mastertable_phage_donors_complete$Row.names <- NULL

Mastertable_phage_figure <- Mastertable_phage_donors_complete[c(265,266,272:275,276,319:320,327)]

ggplot(Mastertable_phage_figure) +
  aes(x = `viral genes`, y = length, colour = Phyla_phylum, size = Totalnumberofreads) +
  geom_point(shape = "circle") +
  scale_color_hue(direction = 1) +
  theme_bw() +
  ylim(0, 200) +
  xlim(0,65)

# Petitvirales
Mastertable_phage_figure %>%
  filter(!(Phyla_phylum %in% c("Euryarchaeota", "Thermodesulfobacteriota", "Mycoplasmota", "Cyanobacteriota", "Mycoplasmatota"))) %>%
  filter(length >= 0L & length <= 10) %>%
  filter((Final_order %in% "Petitvirales")) %>%
  ggplot() +
  aes(x = `viral genes`, y = length, colour = Phyla_phylum, size = Totalnumberofreads) +
  geom_point(shape = "circle") +
  scale_color_hue(direction = 1) +
  theme_bw() +
  facet_wrap(vars(Final_order), scales = "free") +
  xlim(0,10)

# Unclassified Caudoviricetes
Mastertable_phage_figure %>%
  filter(!(Phyla_phylum %in% c("Euryarchaeota", "Thermodesulfobacteriota", "Mycoplasmota", "Cyanobacteriota", "Mycoplasmatota"))) %>%
  filter(length >= 0L & length <= 200) %>%
  filter((Final_order %in% "Unclassified Caudoviricetes")) %>%
  ggplot() +
  aes(x = `viral genes`, y = length, colour = Phyla_phylum, size = Totalnumberofreads) +
  geom_point(shape = "circle") +
  scale_color_hue(direction = 1) +
  theme_bw() +
  facet_wrap(vars(Final_order),scales = "free") +
  xlim(0,60) +
  ylim(10,100)

#  Crassvirales
Mastertable_phage_figure %>%
  filter(!(Phyla_phylum %in% c("Euryarchaeota", "Thermodesulfobacteriota", "Mycoplasmota", "Cyanobacteriota", "Mycoplasmatota"))) %>%
  filter(length >= 0L & length <= 200) %>%
  filter((Final_order %in% "Crassvirales")) %>%
  ggplot() +
  aes(x = `viral genes`, y = length, colour = Phyla_phylum, size = Totalnumberofreads) +
  geom_point(shape = "circle") +
  scale_color_hue(direction = 1) +
  theme_bw() +
  facet_wrap(vars(Final_order),scales = "free") +
  xlim(0,60) +
  ylim(10,200)
####################################