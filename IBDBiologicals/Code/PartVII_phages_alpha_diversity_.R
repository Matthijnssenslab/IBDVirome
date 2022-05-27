####################################
# SCRIPT 8: PROKARYOTIC VIRUSES: DIVERSITY ANALYSES 
####################################
# Before starting this script
####################################
# ------------------> Load the "Global environment" output of script 3 "Rarefaction.R"
# ------------------> You can show presence/absence table without the need for rarefaction based on "Mastertable"
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
#install_github("jokergoo/ComplexHeatmap")
library(ComplexHeatmap)
#if (!require(remotes)) install.packages("remotes")
#remotes::install_github("yiluheihei/microbiomeMarker")
library(microbiomeMarker)
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
#install.packages("coin")
#install.packages("rstatix")
#install.packages("robustrank")
library(phyloseq)
library(rJava)
library(UpSetR)
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
library(rstatix)
library(grid)
library(robustrank)
####################################
# 1. Recapitulation of mastertables
####################################
## Take into account that most of diversity analysis already correct for sequencing depth (alpha-diversity), so you have to feed it the unrarefied completeness data

# global mastertable
Mastertable
Mastertable_unrarefied
Mastertable_rarefied

# viral mastertable
# unrarefied
Mastertable_viral_unrarefied  # PRESENCE/ABSENCE, ALPHA-DIVERSITY (Shannon diversity)

# rarefied
Mastertable_viral_rarefied # RELATIVE ABUNDANCES, RICHNESS, EVENESS 
####################################
# 2. Subset phages
####################################
Mastertable_viral_unrarefied_phages <- Mastertable_viral_unrarefied[Mastertable_viral_unrarefied$Final_viral == "phage",]
Mastertable_viral_rarefied_phages <- Mastertable_viral_rarefied[Mastertable_viral_rarefied$Final_viral == "phage",]
####################################
# 3. Input tables 
####################################
# Default mastertables
Mastertable_viral_unrarefied_phages # unrarefied 
Mastertable_viral_rarefied_phages # rarefied

# Phyloseq mastertables
phyloseq_rarefied_phages
phyloseq_unrarefied_phages

## Remake this Phyloseq tables with classification non-crass, crass
## Need it here
####################################
# 4. Richness and Shannon diversity
####################################
# 4.1 All phages
####################################
# Use rarefied data: observed
alpha_meas <- c("Evennes")
phyloseq_rarefied_phages_alpha <- estimate_richness(phyloseq_rarefied_phages, measures = alpha_meas)
names(phyloseq_rarefied_phages_alpha)[names(phyloseq_rarefied_phages_alpha) == "Observed"] <- "Observed_total_phages"
# names(phyloseq_unrarefied_phages_alpha)[names(phyloseq_unrarefied_phages_alpha) == "Chao1"] <- "Chao1_total_phages"
# names(phyloseq_unrarefied_phages_alpha)[names(phyloseq_unrarefied_phages_alpha) == "Shannon"] <- "Shannon_total_phages"
# names(phyloseq_unrarefied_phages_alpha)[names(phyloseq_unrarefied_phages_alpha) == "Simpson"] <- "Simpson_total_phages"
# names(phyloseq_unrarefied_phages_alpha)[names(phyloseq_unrarefied_phages_alpha) == "InvSimpson"] <- "InvSimpson_total_phages"
iew(phyloseq_rarefied_phages_alpha)

# Use rarefied data: shannon
alpha_meas <- c("Shannon")
phyloseq_unrarefied_phages_alpha <- estimate_richness(phyloseq_unrarefied_phages, measures = alpha_meas)
#names(phyloseq_unrarefied_phages_alpha)[names(phyloseq_unrarefied_phages_alpha) == "Chao1"] <- "Chao1_total_phages"
names(phyloseq_unrarefied_phages_alpha)[names(phyloseq_unrarefied_phages_alpha) == "Shannon"] <- "Shannon_total_phages"
#names(phyloseq_unrarefied_phages_alpha)[names(phyloseq_unrarefied_phages_alpha) == "Simpson"] <- "Simpson_total_phages"
#names(phyloseq_unrarefied_phages_alpha)[names(phyloseq_unrarefied_phages_alpha) == "InvSimpson"] <- "InvSimpson_total_phages"
#View(phyloseq_unrarefied_phages_alpha)
####################################
# 4.2 Caudoviricetes
####################################
# In some samples there is no caudo & therefore 0 richness infinite
# After measures, you can even go deeper and specify the family you want to investigate.
####################################
# 4.2.1 Caudoviricetes (Crass)
####################################
# Use rarefied data: observed
alpha_meas <- c("Observed")
Alpha_diversity_measures_caudo_CrAss_rarefied <- estimate_richness(subset_taxa(phyloseq_rarefied_phages, Class == "Caudoviricetes (CrAss)"), measures = alpha_meas)
names(Alpha_diversity_measures_caudo_CrAss_rarefied)[names(Alpha_diversity_measures_caudo_CrAss_rarefied) == "Observed"] <- "Observed_caudovirales_CrAss"
#View(Alpha_diversity_measures_caudo_CrAss_rarefied)
 
# Use unrarefied data: shannon
alpha_meas <- c("Shannon")
Alpha_diversity_measures_caudo_CrAss_unrarefied <- estimate_richness(subset_taxa(phyloseq_unrarefied_phages, Class == "Caudoviricetes (CrAss)"), measures = alpha_meas)
names(Alpha_diversity_measures_caudo_CrAss_unrarefied)[names(Alpha_diversity_measures_caudo_CrAss_unrarefied) == "Shannon"] <- "Shannon_caudovirales_CrAss"
#View(Alpha_diversity_measures_caudo_CrAss_unrarefied)
####################################
# 4.2.2 Caudoviricetes (non-CrAss)
####################################
# Use rarefied data: observed
alpha_meas <- c("Observed")
Alpha_diversity_measures_caudo_non_CrAss_rarefied <- estimate_richness(subset_taxa(phyloseq_rarefied_phages, Class == "Caudoviricetes (non-CrAss)"), measures = alpha_meas)
names(Alpha_diversity_measures_caudo_non_CrAss_rarefied)[names(Alpha_diversity_measures_caudo_non_CrAss_rarefied) == "Observed"] <- "Observed_caudovirales_non_CrAss"
#View(Alpha_diversity_measures_caudo_non_CrAss_rarefied)

# Use unrarefied data: shannon
alpha_meas <- c("Shannon")
Alpha_diversity_measures_caudo_non_CrAss_unrarefied <- estimate_richness(subset_taxa(phyloseq_unrarefied_phages, Class == "Caudoviricetes (non-CrAss)"), measures = alpha_meas)
names(Alpha_diversity_measures_caudo_non_CrAss_unrarefied)[names(Alpha_diversity_measures_caudo_non_CrAss_unrarefied) == "Shannon"] <- "Shannon_caudovirales_non_CrAss"
#View(Alpha_diversity_measures_caudo_non_CrAss_unrarefied)
####################################
# 4.3 Malgrandaviricetes
####################################
# Use rarefied data: observed
alpha_meas <- c("Observed")
Alpha_diversity_measures_caudo_Malgr_rarefied <- estimate_richness(subset_taxa(phyloseq_rarefied_phages, Class == "Malgrandaviricetes"), measures = alpha_meas)
names(Alpha_diversity_measures_caudo_Malgr_rarefied)[names(Alpha_diversity_measures_caudo_Malgr_rarefied) == "Observed"] <- "Observed_Malgrandaviricetes"
#View(Alpha_diversity_measures_caudo_Malgr_rarefied)

# Use unrarefied data: shannon
alpha_meas <- c("Shannon")
Alpha_diversity_measures_caudo_Malgr_unrarefied <- estimate_richness(subset_taxa(phyloseq_unrarefied_phages, Class == "Malgrandaviricetes"), measures = alpha_meas)
names(Alpha_diversity_measures_caudo_Malgr_unrarefied)[names(Alpha_diversity_measures_caudo_Malgr_unrarefied) == "Shannon"] <- "Shannon_Malgrandaviricetes"
#View(Alpha_diversity_measures_caudo_Malgr_unrarefied)
####################################
# 5. Create DF
####################################
# 5.1 DF
####################################
# DF: rarefied
phyloseq_rarefied_phages_alpha
Alpha_diversity_measures_caudo_CrAss_rarefied
Alpha_diversity_measures_caudo_non_CrAss_rarefied
Alpha_diversity_measures_caudo_Malgr_rarefied

# DF: unrarefied
phyloseq_unrarefied_phages_alpha
Alpha_diversity_measures_caudo_CrAss_unrarefied
Alpha_diversity_measures_caudo_non_CrAss_unrarefied
Alpha_diversity_measures_caudo_Malgr_unrarefied

# sample object
phyloseq_unrarefied_sample_data_ggplot <- as.data.frame(sample_table_unrarefied_for_all)
####################################
# 5.2 Merge DF
####################################
# 5.2.1 Observed Richness
####################################
Mastertable_unrarefied_observed_richness <-merge(phyloseq_unrarefied_sample_data_ggplot,phyloseq_rarefied_phages_alpha[1],by=0,all=F)
rownames(Mastertable_unrarefied_observed_richness) <- Mastertable_unrarefied_observed_richness$Row.names
Mastertable_unrarefied_observed_richness$Row.names <- NULL
Mastertable_unrarefied_observed_richness <-merge(Mastertable_unrarefied_observed_richness,Alpha_diversity_measures_caudo_CrAss_rarefied[1],by=0,all=F)
rownames(Mastertable_unrarefied_observed_richness) <- Mastertable_unrarefied_observed_richness$Row.names
Mastertable_unrarefied_observed_richness$Row.names <- NULL

Mastertable_unrarefied_observed_richness <-merge(Mastertable_unrarefied_observed_richness,Alpha_diversity_measures_caudo_non_CrAss_rarefied[1],by=0,all=F)
rownames(Mastertable_unrarefied_observed_richness) <- Mastertable_unrarefied_observed_richness$Row.names
Mastertable_unrarefied_observed_richness$Row.names <- NULL

Mastertable_unrarefied_observed_richness <-merge(Mastertable_unrarefied_observed_richness,Alpha_diversity_measures_caudo_Malgr_rarefied[1],by=0,all=F)
rownames(Mastertable_unrarefied_observed_richness) <- Mastertable_unrarefied_observed_richness$Row.names
Mastertable_unrarefied_observed_richness$Row.names <- NULL
#View(Mastertable_unrarefied_observed_richness)

## Melt dataframe
vector_2 <- (which(names(Mastertable_unrarefied_observed_richness)== "Diagnosis")) 
vector_2_2 <- (which(names(Mastertable_unrarefied_observed_richness)== "patient_ID")) 
vector_3 <- (which(names(Mastertable_unrarefied_observed_richness)== "Gender")) 
vector_4 <- (which(names(Mastertable_unrarefied_observed_richness)== "Week")) 
vector_5 <- (which(names(Mastertable_unrarefied_observed_richness)== "Endoscopic_outcome_combined (NR/R)")) 
vector_6 <- (which(names(Mastertable_unrarefied_observed_richness)== "Therapy_2")) 
vector_7 <- (which(names(Mastertable_unrarefied_observed_richness)== "Observed_total_phages")) 
vector_8 <- (which(names(Mastertable_unrarefied_observed_richness)== "Observed_Malgrandaviricetes")) 

Mastertable_unrarefied_observed_richness_DF <- Mastertable_unrarefied_observed_richness[,c(vector_2,vector_2_2,vector_3,vector_4,vector_5,vector_6, vector_7:vector_8)]
Mastertable_unrarefied_observed_richness_DF$`Endoscopic_outcome_combined (NR/R)`[Mastertable_unrarefied_observed_richness_DF$`Endoscopic_outcome_combined (NR/R)` == 1] <- "remission"
Mastertable_unrarefied_observed_richness_DF$`Endoscopic_outcome_combined (NR/R)`[Mastertable_unrarefied_observed_richness_DF$`Endoscopic_outcome_combined (NR/R)` == 0] <- "non-remission"
Mastertable_unrarefied_observed_richness_DF$`Endoscopic_outcome_combined (NR/R)`[Mastertable_unrarefied_observed_richness_DF$`Endoscopic_outcome_combined (NR/R)` == "NA"] <- "unknown"
phyloseq_rarefied_phages_ggplot_class_melt <- melt(Mastertable_unrarefied_observed_richness_DF, id.vars = c("patient_ID","Diagnosis", "Gender", "Week", "Endoscopic_outcome_combined (NR/R)", "Therapy_2"))
phyloseq_rarefied_phages_ggplot_class_melt$variable <- as.character(phyloseq_rarefied_phages_ggplot_class_melt$variable)
#View(phyloseq_rarefied_phages_ggplot_class_melt)

## Rename variable back to phage fraction's name
phyloseq_rarefied_phages_ggplot_class_melt$variable[phyloseq_rarefied_phages_ggplot_class_melt$variable == "Observed_total_phages"] <- "phages"
phyloseq_rarefied_phages_ggplot_class_melt$variable[phyloseq_rarefied_phages_ggplot_class_melt$variable == "Observed_Malgrandaviricetes"] <- "Malgrandaviricetes"
phyloseq_rarefied_phages_ggplot_class_melt$variable[phyloseq_rarefied_phages_ggplot_class_melt$variable == "Observed_caudovirales_CrAss"] <- "Caudoviricetes (CrAss)"
phyloseq_rarefied_phages_ggplot_class_melt$variable[phyloseq_rarefied_phages_ggplot_class_melt$variable == "Observed_caudovirales_non_CrAss"] <- "Caudoviricetes (non-CrAss)"

## Change back to a factor for using in a plot
phyloseq_rarefied_phages_ggplot_class_melt$variable <- as.factor(phyloseq_rarefied_phages_ggplot_class_melt$variable)
names(phyloseq_rarefied_phages_ggplot_class_melt)[names(phyloseq_rarefied_phages_ggplot_class_melt) == "variable"] <- "Observed richness"
table(phyloseq_rarefied_phages_ggplot_class_melt$`Observed richness`)
#View(phyloseq_rarefied_phages_ggplot_class_melt)
####################################
# 5.2.2 Shannon diversity
####################################
Mastertable_unrarefied_observed_shannon <- merge(phyloseq_unrarefied_sample_data_ggplot,phyloseq_unrarefied_phages_alpha[1],by=0,all=F)
rownames(Mastertable_unrarefied_observed_shannon) <- Mastertable_unrarefied_observed_shannon$Row.names
Mastertable_unrarefied_observed_shannon$Row.names <- NULL

Mastertable_unrarefied_observed_shannon <-merge(Mastertable_unrarefied_observed_shannon,Alpha_diversity_measures_caudo_CrAss_unrarefied[1],by=0,all=F)
rownames(Mastertable_unrarefied_observed_shannon) <- Mastertable_unrarefied_observed_shannon$Row.names
Mastertable_unrarefied_observed_shannon$Row.names <- NULL

Mastertable_unrarefied_observed_shannon <-merge(Mastertable_unrarefied_observed_shannon,Alpha_diversity_measures_caudo_non_CrAss_unrarefied[1],by=0,all=F)
rownames(Mastertable_unrarefied_observed_shannon) <- Mastertable_unrarefied_observed_shannon$Row.names
Mastertable_unrarefied_observed_shannon$Row.names <- NULL

Mastertable_unrarefied_observed_shannon <-merge(Mastertable_unrarefied_observed_shannon,Alpha_diversity_measures_caudo_Malgr_unrarefied[1],by=0,all=F)
rownames(Mastertable_unrarefied_observed_shannon) <- Mastertable_unrarefied_observed_shannon$Row.names
Mastertable_unrarefied_observed_shannon$Row.names <- NULL

## Melt dataframe
vector_2 <- (which(names(Mastertable_unrarefied_observed_shannon)== "Diagnosis")) 
vector_2_2 <- (which(names(Mastertable_unrarefied_observed_shannon)== "patient_ID")) 
vector_3 <- (which(names(Mastertable_unrarefied_observed_shannon)== "Gender")) 
vector_4 <- (which(names(Mastertable_unrarefied_observed_shannon)== "Week")) 
vector_5 <- (which(names(Mastertable_unrarefied_observed_shannon)== "Endoscopic_outcome_combined (NR/R)")) 
vector_6 <- (which(names(Mastertable_unrarefied_observed_shannon)== "Therapy_2")) 
vector_7 <- (which(names(Mastertable_unrarefied_observed_shannon)== "Shannon_total_phages")) 
vector_8 <- (which(names(Mastertable_unrarefied_observed_shannon)== "Shannon_Malgrandaviricetes")) 

Mastertable_unrarefied_observed_shannon_DF <- Mastertable_unrarefied_observed_shannon[,c(vector_2,vector_2_2,vector_3,vector_4,vector_5,vector_6, vector_7:vector_8)]
Mastertable_unrarefied_observed_shannon_DF$`Endoscopic_outcome_combined (NR/R)`[Mastertable_unrarefied_observed_shannon_DF$`Endoscopic_outcome_combined (NR/R)` == 1] <- "remission"
Mastertable_unrarefied_observed_shannon_DF$`Endoscopic_outcome_combined (NR/R)`[Mastertable_unrarefied_observed_shannon_DF$`Endoscopic_outcome_combined (NR/R)` == 0] <- "non-remission"
Mastertable_unrarefied_observed_shannon_DF$`Endoscopic_outcome_combined (NR/R)`[Mastertable_unrarefied_observed_shannon_DF$`Endoscopic_outcome_combined (NR/R)` == "NA"] <- "unknown"

Mastertable_unrarefied_observed_shannon_DF_melt <- melt(Mastertable_unrarefied_observed_shannon_DF, id.vars = c("patient_ID","Diagnosis", "Gender", "Week", "Endoscopic_outcome_combined (NR/R)", "Therapy_2"))
Mastertable_unrarefied_observed_shannon_DF_melt$variable <- as.character(Mastertable_unrarefied_observed_shannon_DF_melt$variable)
table(Mastertable_unrarefied_observed_shannon_DF_melt$variable)
#View(Mastertable_unrarefied_observed_shannon_DF_melt)

## Rename variable back to phage fraction's name
Mastertable_unrarefied_observed_shannon_DF_melt$variable[Mastertable_unrarefied_observed_shannon_DF_melt$variable == "Shannon_total_phages"] <- "phages"
Mastertable_unrarefied_observed_shannon_DF_melt$variable[Mastertable_unrarefied_observed_shannon_DF_melt$variable == "Shannon_Malgrandaviricetes"] <- "Malgrandaviricetes"
Mastertable_unrarefied_observed_shannon_DF_melt$variable[Mastertable_unrarefied_observed_shannon_DF_melt$variable == "Shannon_caudovirales_CrAss"] <- "Caudoviricetes (CrAss)"
Mastertable_unrarefied_observed_shannon_DF_melt$variable[Mastertable_unrarefied_observed_shannon_DF_melt$variable == "Shannon_caudovirales_non_CrAss"] <- "Caudoviricetes (non-CrAss)"

## Change back to a factor for using in a plot
Mastertable_unrarefied_observed_shannon_DF_melt$variable <- as.factor(Mastertable_unrarefied_observed_shannon_DF_melt$variable)
names(Mastertable_unrarefied_observed_shannon_DF_melt)[names(Mastertable_unrarefied_observed_shannon_DF_melt) == "variable"] <- "Shannon diversity"
table(Mastertable_unrarefied_observed_shannon_DF_melt$`Shannon diversity`)
#View(Mastertable_unrarefied_observed_shannon_DF_melt)
####################################
# 5.2.3 Pielou's eveness
####################################
# H/log(specnumber(BCI))
# A calculated value of Pielou's evenness ranges from 0 (no evenness) to 1 (complete evenness).

# Richness
Mastertable_unrarefied_observed_richness_DF

# Diversity
shannon <- Mastertable_unrarefied_observed_shannon_DF[,7:10]

# Pielou's eveness
Pilous_evennes <- merge(Mastertable_unrarefied_observed_richness_DF,shannon, by = 0, all =F)
rownames(Pilous_evennes) <- Pilous_evennes$Row.names
Pilous_evennes$Row.names <- NULL
Pilous_evennes$Eveness_phages <- Pilous_evennes$Shannon_total_phages/(log(Pilous_evennes$Observed_total_phages))
Pilous_evennes$Eveness_CrAss <- Pilous_evennes$Shannon_caudovirales_CrAss/(log(Pilous_evennes$Observed_caudovirales_CrAss))
Pilous_evennes$Eveness_non_CrAss <- Pilous_evennes$Shannon_caudovirales_non_CrAss/(log(Pilous_evennes$Observed_caudovirales_non_CrAss))
Pilous_evennes$Eveness_Malgranda <- Pilous_evennes$Shannon_Malgrandaviricetes/(log(Pilous_evennes$Observed_Malgrandaviricetes))
Pilous_evennes$Observed_total_phages <- NULL
Pilous_evennes$Observed_caudovirales_CrAss <- NULL
Pilous_evennes$Observed_caudovirales_non_CrAss <- NULL
Pilous_evennes$Observed_Malgrandaviricetes <- NULL
Pilous_evennes$Shannon_caudovirales_CrAss <- NULL
Pilous_evennes$Shannon_caudovirales_non_CrAss <- NULL
Pilous_evennes$Shannon_Malgrandaviricetes <- NULL
Pilous_evennes$Shannon_total_phages <- NULL
####################################
# 4.5 Boxplot
####################################
# 4.6.1 Richness boxplots between timepoints
####################################
# 4.6.1.1 UC
####################################
phyloseq_rarefied_phages_ggplot_class_melt_UC <- phyloseq_rarefied_phages_ggplot_class_melt[phyloseq_rarefied_phages_ggplot_class_melt$Diagnosis == "UC",]
phyloseq_rarefied_phages_ggplot_class_melt_UC$Week[phyloseq_rarefied_phages_ggplot_class_melt_UC$Week == "w0"] <- "baseline"
phyloseq_rarefied_phages_ggplot_class_melt_UC$Week[phyloseq_rarefied_phages_ggplot_class_melt_UC$Week == "w14"] <- "primary endpoint"
phyloseq_rarefied_phages_ggplot_class_melt_UC$`Observed richness` <- factor(phyloseq_rarefied_phages_ggplot_class_melt_UC$`Observed richness`, levels=c("phages", "Caudoviricetes (CrAss)","Caudoviricetes (non-CrAss)", "Malgrandaviricetes"))
phyloseq_rarefied_phages_ggplot_class_melt_UC <- phyloseq_rarefied_phages_ggplot_class_melt_UC[phyloseq_rarefied_phages_ggplot_class_melt_UC$Week == "primary endpoint",]
table(phyloseq_rarefied_phages_ggplot_class_melt_UC$`Observed richness`)

# Figure:
phyloseq_rarefied_phages_ggplot_class_melt_UC %>%
  filter(`Observed richness` %in% c("Caudoviricetes (CrAss)", "Caudoviricetes (non-CrAss)", "Malgrandaviricetes")) %>%
  filter(!`Endoscopic_outcome_combined (NR/R)` %in% c("unknown")) %>%
    ggplot() +
  aes(x = `Observed richness`, y = value, fill = `Endoscopic_outcome_combined (NR/R)`) +
  geom_boxplot() +
  scale_fill_manual(values = c(`non-remission` = "#AD4B38", remission = "#49AD86")) +
  theme_bw() +
  xlab("") +
  ylab("Observed richness") +
  theme(legend.position = "none")  

phyloseq_rarefied_phages_ggplot_class_melt_UC$value[phyloseq_rarefied_phages_ggplot_class_melt_UC$`Observed richness` == "Caudoviricetes (CrAss)" & phyloseq_rarefied_phages_ggplot_class_melt_UC$`Endoscopic_outcome_combined (NR/R)` == "remission"] > 0
# R in UC:  24/31 or 77%
phyloseq_rarefied_phages_ggplot_class_melt_UC$value[phyloseq_rarefied_phages_ggplot_class_melt_UC$`Observed richness` == "Caudoviricetes (CrAss)" & phyloseq_rarefied_phages_ggplot_class_melt_UC$`Endoscopic_outcome_combined (NR/R)` == "non-remission"] > 0
# NR in UC: 12/20 or 60%

prop.test(x = c(24, 12), n = c(31, 20),alternative = "two.sided",conf.level = 0.95, correct = TRUE)
p.adjust(0.3086, method = "BH", n = 2)

# Fig2
phyloseq_rarefied_phages_ggplot_class_melt_UC %>%
 filter(value >= 0L & value <= 80L | is.na(value)) %>%
 ggplot() +
 aes(x = `Endoscopic_outcome_combined (NR/R)`, 
 y = value, fill = `Observed richness`) +
 geom_boxplot(shape = "circle") +
 theme_bw() +
 facet_wrap(vars(Week)) +
 xlab("") +
 ylab("Observed Richness") +
 scale_fill_manual(values=c("#669bbc","#59b031", "#ee766f", "#c49b05"),name="Class") 
####################################
# 4.6.1.2 Statistics within timepoints
####################################
# Statistics within baseline
# Caudo (non-crAss)
table(phyloseq_rarefied_phages_ggplot_class_melt_UC$`Observed richness`)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- phyloseq_rarefied_phages_ggplot_class_melt_UC[c(phyloseq_rarefied_phages_ggplot_class_melt_UC$Week == "baseline" & phyloseq_rarefied_phages_ggplot_class_melt_UC$Diagnosis == "UC" & phyloseq_rarefied_phages_ggplot_class_melt_UC$`Observed richness` == "Caudoviricetes (non-CrAss)"),]
wilcox.test(value~`Endoscopic_outcome_combined (NR/R)`, data = total_virome_R_NR_total_phages_w0_UC_caudo_non_crass, alternative = "two.sided",  paired = FALSE, exact = F)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$remission <- total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$`Endoscopic_outcome_combined (NR/R)`
effect_size <- wilcox_effsize(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass,value~remission)
(effect_size$effsize)^2
effect_size$n1
effect_size$n2
# P-value = 0.6217
# R-squared: 0.00579888
#  n1: 22
#  n2: 22

# Caudo (crAss)
table(phyloseq_rarefied_phages_ggplot_class_melt_UC$`Observed richness`)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- phyloseq_rarefied_phages_ggplot_class_melt_UC[c(phyloseq_rarefied_phages_ggplot_class_melt_UC$Week == "baseline" & phyloseq_rarefied_phages_ggplot_class_melt_UC$Diagnosis == "UC" & phyloseq_rarefied_phages_ggplot_class_melt_UC$`Observed richness` == "Caudoviricetes (CrAss)"),]
wilcox.test(value~`Endoscopic_outcome_combined (NR/R)`, data = total_virome_R_NR_total_phages_w0_UC_caudo_non_crass, alternative = "two.sided",  paired = FALSE, exact = F)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$remission <- total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$`Endoscopic_outcome_combined (NR/R)`
effect_size <- wilcox_effsize(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass,value~remission)
(effect_size$effsize)^2
effect_size$n1
effect_size$n2
# P-value = 0.1014
# R-squared: 0.06188915
#  n1: 22
#  n2: 22

# Caudo (Micro)
table(phyloseq_rarefied_phages_ggplot_class_melt_UC$`Observed richness`)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- phyloseq_rarefied_phages_ggplot_class_melt_UC[c(phyloseq_rarefied_phages_ggplot_class_melt_UC$Week == "baseline" & phyloseq_rarefied_phages_ggplot_class_melt_UC$Diagnosis == "UC" & phyloseq_rarefied_phages_ggplot_class_melt_UC$`Observed richness` == "Malgrandaviricetes"),]
wilcox.test(value~`Endoscopic_outcome_combined (NR/R)`, data = total_virome_R_NR_total_phages_w0_UC_caudo_non_crass, alternative = "two.sided",  paired = FALSE, exact = F)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$remission <- total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$`Endoscopic_outcome_combined (NR/R)`
effect_size <- wilcox_effsize(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass,value~remission)
(effect_size$effsize)^2
effect_size$n1
effect_size$n2
# P-value = 0.06478
# R-squared: 0.07850333
#  n1: 22
#  n2: 22

# Total phages
table(phyloseq_rarefied_phages_ggplot_class_melt_UC$`Observed richness`)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- phyloseq_rarefied_phages_ggplot_class_melt_UC[c(phyloseq_rarefied_phages_ggplot_class_melt_UC$Week == "baseline" & phyloseq_rarefied_phages_ggplot_class_melt_UC$Diagnosis == "UC" & phyloseq_rarefied_phages_ggplot_class_melt_UC$`Observed richness` == "phages"),]
wilcox.test(value~`Endoscopic_outcome_combined (NR/R)`, data = total_virome_R_NR_total_phages_w0_UC_caudo_non_crass, alternative = "two.sided",  paired = FALSE, exact = F)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$remission <- total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$`Endoscopic_outcome_combined (NR/R)`
effect_size <- wilcox_effsize(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass,value~remission)
(effect_size$effsize)^2
effect_size$n1
effect_size$n2
# P-value = 0.2355
# R-squared: 0.03261547
#  n1: 22
#  n2: 22

# Statistics within primary endpoint
# Caudo (non-crAss)
phyloseq_rarefied_phages_ggplot_class_melt_UC$`Endoscopic_outcome_combined (NR/R)`
table(phyloseq_rarefied_phages_ggplot_class_melt_UC$`Observed richness`)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- phyloseq_rarefied_phages_ggplot_class_melt_UC[c(phyloseq_rarefied_phages_ggplot_class_melt_UC$Week == "primary endpoint" & phyloseq_rarefied_phages_ggplot_class_melt_UC$Diagnosis == "UC" & phyloseq_rarefied_phages_ggplot_class_melt_UC$`Observed richness` == "Caudoviricetes (non-CrAss)"),]
wilcox.test(value~`Endoscopic_outcome_combined (NR/R)`, data = total_virome_R_NR_total_phages_w0_UC_caudo_non_crass, alternative = "two.sided",  paired = FALSE, exact = F)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$remission <- total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$`Endoscopic_outcome_combined (NR/R)`
p.adjust(0.401, method = "BH", n=3)
effect_size <- wilcox_effsize(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass,value~remission)
(effect_size$effsize)
effect_size$n1
effect_size$n2
# P-value = 0.401
# R-squared: 0.01414725
#  n1: 20
#  n2: 31

# Caudo (crAss)
table(phyloseq_rarefied_phages_ggplot_class_melt_UC$`Observed richness`)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- phyloseq_rarefied_phages_ggplot_class_melt_UC[c(phyloseq_rarefied_phages_ggplot_class_melt_UC$Week == "primary endpoint" & phyloseq_rarefied_phages_ggplot_class_melt_UC$Diagnosis == "UC" & phyloseq_rarefied_phages_ggplot_class_melt_UC$`Observed richness` == "Caudoviricetes (CrAss)"),]
wilcox.test(value~`Endoscopic_outcome_combined (NR/R)`, data = total_virome_R_NR_total_phages_w0_UC_caudo_non_crass, alternative = "two.sided",  paired = FALSE, exact = F)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$remission <- total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$`Endoscopic_outcome_combined (NR/R)`
p.adjust(0.143, method = "BH", n=3)
effect_size <- wilcox_effsize(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass,value~remission)
(effect_size$effsize)
effect_size$n1
effect_size$n2
# P-value = 0.1427
# R-squared: 0.04268847
#  n1: 20
#  n2: 31

# Caudo (Micro)
table(phyloseq_rarefied_phages_ggplot_class_melt_UC$`Observed richness`)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- phyloseq_rarefied_phages_ggplot_class_melt_UC[c(phyloseq_rarefied_phages_ggplot_class_melt_UC$Week == "primary endpoint" & phyloseq_rarefied_phages_ggplot_class_melt_UC$Diagnosis == "UC" & phyloseq_rarefied_phages_ggplot_class_melt_UC$`Observed richness` == "Malgrandaviricetes"),]
wilcox.test(value~`Endoscopic_outcome_combined (NR/R)`, data = total_virome_R_NR_total_phages_w0_UC_caudo_non_crass, alternative = "two.sided",  paired = FALSE, exact = F, correct = F)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$remission <- total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$`Endoscopic_outcome_combined (NR/R)`
p.adjust(0.02522, method = "BH", n=3)
effect_size <- wilcox_effsize(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass,value~remission)
(effect_size$effsize)
effect_size$n1
effect_size$n2
# P-value = 0.02586****
# R-squared: 0.09820416
#  n1: 20
#  n2: 31

# Total phages
table(phyloseq_rarefied_phages_ggplot_class_melt_UC$`Observed richness`)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- phyloseq_rarefied_phages_ggplot_class_melt_UC[c(phyloseq_rarefied_phages_ggplot_class_melt_UC$Week == "primary endpoint" & phyloseq_rarefied_phages_ggplot_class_melt_UC$Diagnosis == "UC" & phyloseq_rarefied_phages_ggplot_class_melt_UC$`Observed richness` == "phages"),]
wilcox.test(value~`Endoscopic_outcome_combined (NR/R)`, data = total_virome_R_NR_total_phages_w0_UC_caudo_non_crass, alternative = "two.sided",  paired = FALSE, exact = F)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$remission <- total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$`Endoscopic_outcome_combined (NR/R)`
effect_size <- wilcox_effsize(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass,value~remission)
(effect_size$effsize)^2
effect_size$n1
effect_size$n2
# P-value = 0.1178
# R-squared: 0.04854899
#  n1: 20
#  n2: 31
####################################
# 4.6.1.3 Statistics between timepoints
####################################
# Statistics between baseline and primary endpoint
# Statistics within Remission
# Caudo (non-crAss)
table(phyloseq_rarefied_phages_ggplot_class_melt_UC$`Observed richness`)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- phyloseq_rarefied_phages_ggplot_class_melt_UC[c(phyloseq_rarefied_phages_ggplot_class_melt_UC$`Endoscopic_outcome_combined (NR/R)` == "remission" & phyloseq_rarefied_phages_ggplot_class_melt_UC$Diagnosis == "UC" & phyloseq_rarefied_phages_ggplot_class_melt_UC$`Observed richness` == "Caudoviricetes (non-CrAss)"),]
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired <- total_virome_R_NR_total_phages_w0_UC_caudo_non_crass
# View(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired)

colnames_de <- colnames(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired)
de1 <-data.frame("5520",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("7903",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("11133",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("12255",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("12494",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("13419",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("13523",NA, NA,"primary endpoint",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("13525",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("10213",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("3678",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("11057",NA, NA,"baseline",NA,NA,NA,NA)

#View(de1)
names(de1)<-c(colnames_de)
unpaired_paired_test <- rbind(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired, de1)
unpaired_paired_test_baseline <- unpaired_paired_test[unpaired_paired_test$Week == "baseline",]
unpaired_paired_test_baseline <- unpaired_paired_test_baseline[order(unpaired_paired_test_baseline$patient_ID) , ]
unpaired_paired_test_baseline$value <- as.numeric(unpaired_paired_test_baseline$value)
unpaired_paired_test_PE <- unpaired_paired_test[!unpaired_paired_test$Week == "baseline",]
unpaired_paired_test_PE <- unpaired_paired_test_PE[order(unpaired_paired_test_PE$patient_ID) , ]
unpaired_paired_test_PE$value <- as.numeric(unpaired_paired_test_PE$value)

View(unpaired_paired_test)
#View(unpaired_paired_test_baseline)
#View(unpaired_paired_test_PE)

pm.wilcox.test(unpaired_paired_test_PE$value,unpaired_paired_test_baseline$value,alternative = "two.sided", correct = F, method = "MW-MW")

# P-value = 0.6347 (x: baseline, y: PE)
# Xpaired: 21
# Ypaired: 21
# Xextra: 10
# Y-extra: 1

# Caudo (crAss)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- phyloseq_rarefied_phages_ggplot_class_melt_UC[c(phyloseq_rarefied_phages_ggplot_class_melt_UC$`Endoscopic_outcome_combined (NR/R)` == "remission" & phyloseq_rarefied_phages_ggplot_class_melt_UC$Diagnosis == "UC" & phyloseq_rarefied_phages_ggplot_class_melt_UC$`Observed richness` == "Caudoviricetes (CrAss)"),]
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired <- total_virome_R_NR_total_phages_w0_UC_caudo_non_crass
#View(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired)

colnames_de <- colnames(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired)
de1 <-data.frame("5520",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("7903",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("11133",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("12255",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("12494",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("13419",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("13523",NA, NA,"primary endpoint",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("13525",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("10213",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("3678",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("11057",NA, NA,"baseline",NA,NA,NA,NA)

#View(de1)
names(de1)<-c(colnames_de)
unpaired_paired_test <- rbind(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired, de1)
unpaired_paired_test_baseline <- unpaired_paired_test[unpaired_paired_test$Week == "baseline",]
unpaired_paired_test_baseline <- unpaired_paired_test_baseline[order(unpaired_paired_test_baseline$patient_ID) , ]
unpaired_paired_test_baseline$value <- as.numeric(unpaired_paired_test_baseline$value)
unpaired_paired_test_PE <- unpaired_paired_test[!unpaired_paired_test$Week == "baseline",]
unpaired_paired_test_PE <- unpaired_paired_test_PE[order(unpaired_paired_test_PE$patient_ID) , ]
unpaired_paired_test_PE$value <- as.numeric(unpaired_paired_test_PE$value)

#View(unpaired_paired_test)
#View(unpaired_paired_test_baseline)
#View(unpaired_paired_test_PE)

pm.wilcox.test(unpaired_paired_test_PE$value,unpaired_paired_test_baseline$value,alternative = "two.sided", correct = F, method = "MW-MW")

# P-value = 0.5965 (x: baseline, y: PE)
# Xpaired: 21
# Ypaired: 21
# Xextra: 10
# Y-extra: 1

# Micro
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- phyloseq_rarefied_phages_ggplot_class_melt_UC[c(phyloseq_rarefied_phages_ggplot_class_melt_UC$`Endoscopic_outcome_combined (NR/R)` == "remission" & phyloseq_rarefied_phages_ggplot_class_melt_UC$Diagnosis == "UC" & phyloseq_rarefied_phages_ggplot_class_melt_UC$`Observed richness` == "Malgrandaviricetes"),]
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired <- total_virome_R_NR_total_phages_w0_UC_caudo_non_crass
#View(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired)

colnames_de <- colnames(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired)
de1 <-data.frame("5520",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("7903",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("11133",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("12255",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("12494",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("13419",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("13523",NA, NA,"primary endpoint",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("13525",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("10213",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("3678",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("11057",NA, NA,"baseline",NA,NA,NA,NA)

#View(de1)
names(de1)<-c(colnames_de)
unpaired_paired_test <- rbind(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired, de1)
unpaired_paired_test_baseline <- unpaired_paired_test[unpaired_paired_test$Week == "baseline",]
unpaired_paired_test_baseline <- unpaired_paired_test_baseline[order(unpaired_paired_test_baseline$patient_ID) , ]
unpaired_paired_test_baseline$value <- as.numeric(unpaired_paired_test_baseline$value)
unpaired_paired_test_PE <- unpaired_paired_test[!unpaired_paired_test$Week == "baseline",]
unpaired_paired_test_PE <- unpaired_paired_test_PE[order(unpaired_paired_test_PE$patient_ID) , ]
unpaired_paired_test_PE$value <- as.numeric(unpaired_paired_test_PE$value)

#View(unpaired_paired_test)
#View(unpaired_paired_test_baseline)
#View(unpaired_paired_test_PE)

pm.wilcox.test(unpaired_paired_test_PE$value,unpaired_paired_test_baseline$value,alternative = "two.sided", correct = F, method = "MW-MW")

# P-value = 0.8072 (x: baseline, y: PE)
# Xpaired: 21
# Ypaired: 21
# Xextra: 10
# Y-extra: 1

# Phages
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- phyloseq_rarefied_phages_ggplot_class_melt_UC[c(phyloseq_rarefied_phages_ggplot_class_melt_UC$`Endoscopic_outcome_combined (NR/R)` == "remission" & phyloseq_rarefied_phages_ggplot_class_melt_UC$Diagnosis == "UC" & phyloseq_rarefied_phages_ggplot_class_melt_UC$`Observed richness` == "phages"),]
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired <- total_virome_R_NR_total_phages_w0_UC_caudo_non_crass
#View(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired)

colnames_de <- colnames(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired)
de1 <-data.frame("5520",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("7903",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("11133",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("12255",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("12494",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("13419",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("13523",NA, NA,"primary endpoint",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("13525",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("10213",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("3678",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("11057",NA, NA,"baseline",NA,NA,NA,NA)

#View(de1)
names(de1)<-c(colnames_de)
unpaired_paired_test <- rbind(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired, de1)
unpaired_paired_test_baseline <- unpaired_paired_test[unpaired_paired_test$Week == "baseline",]
unpaired_paired_test_baseline <- unpaired_paired_test_baseline[order(unpaired_paired_test_baseline$patient_ID) , ]
unpaired_paired_test_baseline$value <- as.numeric(unpaired_paired_test_baseline$value)
unpaired_paired_test_PE <- unpaired_paired_test[!unpaired_paired_test$Week == "baseline",]
unpaired_paired_test_PE <- unpaired_paired_test_PE[order(unpaired_paired_test_PE$patient_ID) , ]
unpaired_paired_test_PE$value <- as.numeric(unpaired_paired_test_PE$value)

#View(unpaired_paired_test)
#View(unpaired_paired_test_baseline)
#View(unpaired_paired_test_PE)

pm.wilcox.test(unpaired_paired_test_PE$value,unpaired_paired_test_baseline$value,alternative = "two.sided", correct = F, method = "MW-MW")

# P-value = 0.5403 (x: baseline, y: PE)
# Xpaired: 21
# Ypaired: 21
# Xextra: 10
# Y-extra: 1

# Statistics within non-Remission
# Caudo (non-crAss)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- phyloseq_rarefied_phages_ggplot_class_melt_UC[c(phyloseq_rarefied_phages_ggplot_class_melt_UC$`Endoscopic_outcome_combined (NR/R)` == "non-remission" & phyloseq_rarefied_phages_ggplot_class_melt_UC$Diagnosis == "UC" & phyloseq_rarefied_phages_ggplot_class_melt_UC$`Observed richness` == "Caudoviricetes (non-CrAss)"),]
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired <- total_virome_R_NR_total_phages_w0_UC_caudo_non_crass
#View(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired)

colnames_de <- colnames(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired)
de1 <-data.frame("7164",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("7270",NA, NA,"primary endpoint",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("10449",NA, NA,"primary endpoint",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("12515",NA, NA,"primary endpoint",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("13453",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("9516",NA, NA,"primary endpoint",NA,"VDZ",NA,NA)

#View(de1)
names(de1)<-c(colnames_de)
unpaired_paired_test <- rbind(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired, de1)
unpaired_paired_test$patient_ID[c(unpaired_paired_test$patient_ID == "12770" & unpaired_paired_test$Therapy_2 == "TNF")] <- "1"
unpaired_paired_test$patient_ID[c(unpaired_paired_test$patient_ID == "13093" & unpaired_paired_test$Therapy_2 == "TNF")] <- "2"
unpaired_paired_test$patient_ID[c(unpaired_paired_test$patient_ID == "9516" & unpaired_paired_test$Therapy_2 == "TNF")] <- "3"

unpaired_paired_test_baseline <- unpaired_paired_test[unpaired_paired_test$Week == "baseline",]
unpaired_paired_test_baseline <- unpaired_paired_test_baseline[order(unpaired_paired_test_baseline$patient_ID) , ]
unpaired_paired_test_baseline$value <- as.numeric(unpaired_paired_test_baseline$value)
unpaired_paired_test_PE <- unpaired_paired_test[!unpaired_paired_test$Week == "baseline",]
unpaired_paired_test_PE <- unpaired_paired_test_PE[order(unpaired_paired_test_PE$patient_ID) , ]
unpaired_paired_test_PE$value <- as.numeric(unpaired_paired_test_PE$value)

#View(unpaired_paired_test)
#View(unpaired_paired_test_baseline)
#View(unpaired_paired_test_PE)

pm.wilcox.test(unpaired_paired_test_PE$value,unpaired_paired_test_baseline$value,alternative = "two.sided", correct = F, method = "MW-MW")

# P-value = 0.834 (x: baseline, y: PE)
# Xpaired: 18
# Ypaired: 18
# Xextra: 2
# Y-extra: 4

# Caudo (crAss)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- phyloseq_rarefied_phages_ggplot_class_melt_UC[c(phyloseq_rarefied_phages_ggplot_class_melt_UC$`Endoscopic_outcome_combined (NR/R)` == "non-remission" & phyloseq_rarefied_phages_ggplot_class_melt_UC$Diagnosis == "UC" & phyloseq_rarefied_phages_ggplot_class_melt_UC$`Observed richness` == "Caudoviricetes (CrAss)"),]
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired <- total_virome_R_NR_total_phages_w0_UC_caudo_non_crass
#View(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired)

colnames_de <- colnames(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired)
de1 <-data.frame("7164",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("7270",NA, NA,"primary endpoint",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("10449",NA, NA,"primary endpoint",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("12515",NA, NA,"primary endpoint",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("13453",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("9516",NA, NA,"primary endpoint",NA,"VDZ",NA,NA)

#View(de1)
names(de1)<-c(colnames_de)
unpaired_paired_test <- rbind(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired, de1)
unpaired_paired_test$patient_ID[c(unpaired_paired_test$patient_ID == "12770" & unpaired_paired_test$Therapy_2 == "TNF")] <- "1"
unpaired_paired_test$patient_ID[c(unpaired_paired_test$patient_ID == "13093" & unpaired_paired_test$Therapy_2 == "TNF")] <- "2"
unpaired_paired_test$patient_ID[c(unpaired_paired_test$patient_ID == "9516" & unpaired_paired_test$Therapy_2 == "TNF")] <- "3"

unpaired_paired_test_baseline <- unpaired_paired_test[unpaired_paired_test$Week == "baseline",]
unpaired_paired_test_baseline <- unpaired_paired_test_baseline[order(unpaired_paired_test_baseline$patient_ID) , ]
unpaired_paired_test_baseline$value <- as.numeric(unpaired_paired_test_baseline$value)
unpaired_paired_test_PE <- unpaired_paired_test[!unpaired_paired_test$Week == "baseline",]
unpaired_paired_test_PE <- unpaired_paired_test_PE[order(unpaired_paired_test_PE$patient_ID) , ]
unpaired_paired_test_PE$value <- as.numeric(unpaired_paired_test_PE$value)

#View(unpaired_paired_test)
#View(unpaired_paired_test_baseline)
#View(unpaired_paired_test_PE)

pm.wilcox.test(unpaired_paired_test_PE$value,unpaired_paired_test_baseline$value,alternative = "two.sided", correct = F, method = "MW-MW")

# P-value = 0.834 (x: baseline, y: PE)
# Xpaired: 18
# Ypaired: 18
# Xextra: 2
# Y-extra: 4

# Micro
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- phyloseq_rarefied_phages_ggplot_class_melt_UC[c(phyloseq_rarefied_phages_ggplot_class_melt_UC$`Endoscopic_outcome_combined (NR/R)` == "non-remission" & phyloseq_rarefied_phages_ggplot_class_melt_UC$Diagnosis == "UC" & phyloseq_rarefied_phages_ggplot_class_melt_UC$`Observed richness` == "Malgrandaviricetes"),]
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired <- total_virome_R_NR_total_phages_w0_UC_caudo_non_crass
#View(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired)

colnames_de <- colnames(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired)
de1 <-data.frame("7164",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("7270",NA, NA,"primary endpoint",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("10449",NA, NA,"primary endpoint",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("12515",NA, NA,"primary endpoint",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("13453",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("9516",NA, NA,"primary endpoint",NA,"VDZ",NA,NA)

#View(de1)
names(de1)<-c(colnames_de)
unpaired_paired_test <- rbind(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired, de1)
unpaired_paired_test$patient_ID[c(unpaired_paired_test$patient_ID == "12770" & unpaired_paired_test$Therapy_2 == "TNF")] <- "1"
unpaired_paired_test$patient_ID[c(unpaired_paired_test$patient_ID == "13093" & unpaired_paired_test$Therapy_2 == "TNF")] <- "2"
unpaired_paired_test$patient_ID[c(unpaired_paired_test$patient_ID == "9516" & unpaired_paired_test$Therapy_2 == "TNF")] <- "3"

unpaired_paired_test_baseline <- unpaired_paired_test[unpaired_paired_test$Week == "baseline",]
unpaired_paired_test_baseline <- unpaired_paired_test_baseline[order(unpaired_paired_test_baseline$patient_ID) , ]
unpaired_paired_test_baseline$value <- as.numeric(unpaired_paired_test_baseline$value)
unpaired_paired_test_PE <- unpaired_paired_test[!unpaired_paired_test$Week == "baseline",]
unpaired_paired_test_PE <- unpaired_paired_test_PE[order(unpaired_paired_test_PE$patient_ID) , ]
unpaired_paired_test_PE$value <- as.numeric(unpaired_paired_test_PE$value)

#View(unpaired_paired_test)
#View(unpaired_paired_test_baseline)
#View(unpaired_paired_test_PE)

pm.wilcox.test(unpaired_paired_test_PE$value,unpaired_paired_test_baseline$value,alternative = "two.sided", correct = F, method = "MW-MW")

# P-value = 0.4508 (x: baseline, y: PE)
# Xpaired: 18
# Ypaired: 18
# Xextra: 2
# Y-extra: 4

# Phages
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- phyloseq_rarefied_phages_ggplot_class_melt_UC[c(phyloseq_rarefied_phages_ggplot_class_melt_UC$`Endoscopic_outcome_combined (NR/R)` == "non-remission" & phyloseq_rarefied_phages_ggplot_class_melt_UC$Diagnosis == "UC" & phyloseq_rarefied_phages_ggplot_class_melt_UC$`Observed richness` == "phages"),]
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired <- total_virome_R_NR_total_phages_w0_UC_caudo_non_crass
#View(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired)

colnames_de <- colnames(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired)
de1 <-data.frame("7164",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("7270",NA, NA,"primary endpoint",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("10449",NA, NA,"primary endpoint",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("12515",NA, NA,"primary endpoint",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("13453",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("9516",NA, NA,"primary endpoint",NA,"VDZ",NA,NA)

#View(de1)
names(de1)<-c(colnames_de)
unpaired_paired_test <- rbind(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired, de1)
unpaired_paired_test$patient_ID[c(unpaired_paired_test$patient_ID == "12770" & unpaired_paired_test$Therapy_2 == "TNF")] <- "1"
unpaired_paired_test$patient_ID[c(unpaired_paired_test$patient_ID == "13093" & unpaired_paired_test$Therapy_2 == "TNF")] <- "2"
unpaired_paired_test$patient_ID[c(unpaired_paired_test$patient_ID == "9516" & unpaired_paired_test$Therapy_2 == "TNF")] <- "3"

unpaired_paired_test_baseline <- unpaired_paired_test[unpaired_paired_test$Week == "baseline",]
unpaired_paired_test_baseline <- unpaired_paired_test_baseline[order(unpaired_paired_test_baseline$patient_ID) , ]
unpaired_paired_test_baseline$value <- as.numeric(unpaired_paired_test_baseline$value)
unpaired_paired_test_PE <- unpaired_paired_test[!unpaired_paired_test$Week == "baseline",]
unpaired_paired_test_PE <- unpaired_paired_test_PE[order(unpaired_paired_test_PE$patient_ID) , ]
unpaired_paired_test_PE$value <- as.numeric(unpaired_paired_test_PE$value)

View(unpaired_paired_test)
#View(unpaired_paired_test_baseline)
#View(unpaired_paired_test_PE)

pm.wilcox.test(unpaired_paired_test_PE$value,unpaired_paired_test_baseline$value,alternative = "two.sided", correct = F, method = "MW-MW")

# P-value = 0.9803 (x: baseline, y: PE)
# Xpaired: 18
# Ypaired: 18
# Xextra: 2
# Y-extra: 4
####################################
# 4.6.1.4 CD
####################################
phyloseq_rarefied_phages_ggplot_class_melt_CD <- phyloseq_rarefied_phages_ggplot_class_melt[phyloseq_rarefied_phages_ggplot_class_melt$Diagnosis == "CD",]
phyloseq_rarefied_phages_ggplot_class_melt_CD$Week[phyloseq_rarefied_phages_ggplot_class_melt_CD$Week == "w0"] <- "baseline"
phyloseq_rarefied_phages_ggplot_class_melt_CD$Week[phyloseq_rarefied_phages_ggplot_class_melt_CD$Week == "w24"] <- "primary endpoint"
phyloseq_rarefied_phages_ggplot_class_melt_CD <- phyloseq_rarefied_phages_ggplot_class_melt_CD[!phyloseq_rarefied_phages_ggplot_class_melt_CD$Week == "w14",]
phyloseq_rarefied_phages_ggplot_class_melt_CD <- phyloseq_rarefied_phages_ggplot_class_melt_CD[!phyloseq_rarefied_phages_ggplot_class_melt_CD$`Endoscopic_outcome_combined (NR/R)` == "unknown",]
phyloseq_rarefied_phages_ggplot_class_melt_CD$`Observed richness` <- factor(phyloseq_rarefied_phages_ggplot_class_melt_CD$`Observed richness`, levels=c("phages", "Caudoviricetes (CrAss)","Caudoviricetes (non-CrAss)", "Malgrandaviricetes"))
table(phyloseq_rarefied_phages_ggplot_class_melt_CD$`Observed richness`)

phyloseq_rarefied_phages_ggplot_class_melt_CD %>%
  ggplot() +
  filter(value >= 0L & value <= 80L | is.na(value)) %>%
  aes(x = `Endoscopic_outcome_combined (NR/R)`, 
      y = value, fill = `Observed richness`) +
  geom_boxplot(shape = "circle") +
  theme_bw() +
  facet_wrap(vars(Week)) +
  scale_y_continuous(breaks=c(0,20,40,60,80), limits = c(0,80)) +
  xlab("") +
  ylab("Observed Richness") +
  scale_fill_manual(values=c("#669bbc","#59b031", "#ee766f", "#c49b05"),name="Class") 

phyloseq_rarefied_phages_ggplot_class_melt_CD$value[phyloseq_rarefied_phages_ggplot_class_melt_CD$`Observed richness` == "Caudoviricetes (CrAss)" & phyloseq_rarefied_phages_ggplot_class_melt_CD$`Endoscopic_outcome_combined (NR/R)` == "remission"]
table(phyloseq_rarefied_phages_ggplot_class_melt_CD$value[phyloseq_rarefied_phages_ggplot_class_melt_CD$`Observed richness` == "Caudoviricetes (CrAss)" & phyloseq_rarefied_phages_ggplot_class_melt_CD$`Endoscopic_outcome_combined (NR/R)` == "remission"]) 
# R in CD:  64/89 or 57.9%
89-25
64/89

phyloseq_rarefied_phages_ggplot_class_melt_CD$value[phyloseq_rarefied_phages_ggplot_class_melt_CD$`Observed richness` == "Caudoviricetes (CrAss)" & phyloseq_rarefied_phages_ggplot_class_melt_CD$`Endoscopic_outcome_combined (NR/R)` == "non-remission"]
table(phyloseq_rarefied_phages_ggplot_class_melt_CD$value[phyloseq_rarefied_phages_ggplot_class_melt_CD$`Observed richness` == "Caudoviricetes (CrAss)" & phyloseq_rarefied_phages_ggplot_class_melt_CD$`Endoscopic_outcome_combined (NR/R)` == "non-remission"])
# NR in CD: 59/95 or 71.9%
95-36

89+95
prop.test(x = c(64, 59), n = c(89, 95),alternative = "two.sided",conf.level = 0.95, correct = TRUE)
p.adjust(0.2094, method = "BH", n = 2)
####################################
# 4.6.1.5 Statistics within timepoints
####################################
# Statistics within baseline
# Caudo (non-crAss)
table(phyloseq_rarefied_phages_ggplot_class_melt_CD$`Observed richness`)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- phyloseq_rarefied_phages_ggplot_class_melt_CD[c(phyloseq_rarefied_phages_ggplot_class_melt_CD$Week == "baseline" & phyloseq_rarefied_phages_ggplot_class_melt_CD$Diagnosis == "CD" & phyloseq_rarefied_phages_ggplot_class_melt_CD$`Observed richness` == "Caudoviricetes (non-CrAss)"),]
table(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$`Endoscopic_outcome_combined (NR/R)`) ## 22 R patients/ 22 NR patients
wilcox.test(value~`Endoscopic_outcome_combined (NR/R)`, data = total_virome_R_NR_total_phages_w0_UC_caudo_non_crass, alternative = "two.sided",  paired = FALSE, exact = F)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$remission <- total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$`Endoscopic_outcome_combined (NR/R)`
effect_size <- wilcox_effsize(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass,value~remission)
(effect_size$effsize)^2
effect_size$n1
effect_size$n2
# P-value = 0.6537
# R-squared: 0.002276187
#  n1: 47
#  n2: 43

# Caudo (crAss)
table(phyloseq_rarefied_phages_ggplot_class_melt_CD$`Observed richness`)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- phyloseq_rarefied_phages_ggplot_class_melt_CD[c(phyloseq_rarefied_phages_ggplot_class_melt_CD$Week == "baseline" & phyloseq_rarefied_phages_ggplot_class_melt_CD$Diagnosis == "CD" & phyloseq_rarefied_phages_ggplot_class_melt_CD$`Observed richness` == "Caudoviricetes (CrAss)"),]
table(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$`Endoscopic_outcome_combined (NR/R)`) ## 22 R patients/ 22 NR patients
wilcox.test(value~`Endoscopic_outcome_combined (NR/R)`, data = total_virome_R_NR_total_phages_w0_UC_caudo_non_crass, alternative = "two.sided",  paired = FALSE, exact = F)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$remission <- total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$`Endoscopic_outcome_combined (NR/R)`
effect_size <- wilcox_effsize(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass,value~remission)
(effect_size$effsize)^2
effect_size$n1
effect_size$n2
# P-value = 0.05924
# R-squared: 0.03971436
#  n1: 47
#  n2: 43

# Malgrandaviricetes
table(phyloseq_rarefied_phages_ggplot_class_melt_CD$`Observed richness`)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- phyloseq_rarefied_phages_ggplot_class_melt_CD[c(phyloseq_rarefied_phages_ggplot_class_melt_CD$Week == "baseline" & phyloseq_rarefied_phages_ggplot_class_melt_CD$Diagnosis == "CD" & phyloseq_rarefied_phages_ggplot_class_melt_CD$`Observed richness` == "Malgrandaviricetes"),]
table(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$`Endoscopic_outcome_combined (NR/R)`) ## 22 R patients/ 22 NR patients
wilcox.test(value~`Endoscopic_outcome_combined (NR/R)`, data = total_virome_R_NR_total_phages_w0_UC_caudo_non_crass, alternative = "two.sided",  paired = FALSE, exact = F)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$remission <- total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$`Endoscopic_outcome_combined (NR/R)`
effect_size <- wilcox_effsize(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass,value~remission)
(effect_size$effsize)^2
effect_size$n1
effect_size$n2
# P-value = 0.08634
# R-squared: 0.0328366
#  n1: 47
#  n2: 43

# Total phages
table(phyloseq_rarefied_phages_ggplot_class_melt_CD$`Observed richness`)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- phyloseq_rarefied_phages_ggplot_class_melt_CD[c(phyloseq_rarefied_phages_ggplot_class_melt_CD$Week == "baseline" & phyloseq_rarefied_phages_ggplot_class_melt_CD$Diagnosis == "CD" & phyloseq_rarefied_phages_ggplot_class_melt_CD$`Observed richness` == "phages"),]
table(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$`Endoscopic_outcome_combined (NR/R)`) ## 22 R patients/ 22 NR patients
wilcox.test(value~`Endoscopic_outcome_combined (NR/R)`, data = total_virome_R_NR_total_phages_w0_UC_caudo_non_crass, alternative = "two.sided",  paired = FALSE, exact = F)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$remission <- total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$`Endoscopic_outcome_combined (NR/R)`
effect_size <- wilcox_effsize(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass,value~remission)
(effect_size$effsize)^2
effect_size$n1
effect_size$n2
# P-value = 0.605
# R-squared: 0.003019488
#  n1: 47
#  n2: 43

# Statistics within primary endpoint
# Caudo (non-crAss)
table(phyloseq_rarefied_phages_ggplot_class_melt_CD$`Observed richness`)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- phyloseq_rarefied_phages_ggplot_class_melt_CD[c(phyloseq_rarefied_phages_ggplot_class_melt_CD$Week == "primary endpoint" & phyloseq_rarefied_phages_ggplot_class_melt_CD$Diagnosis == "CD" & phyloseq_rarefied_phages_ggplot_class_melt_CD$`Observed richness` == "Caudoviricetes (non-CrAss)"),]
table(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$`Endoscopic_outcome_combined (NR/R)`) ## 22 R patients/ 22 NR patients
wilcox.test(value~`Endoscopic_outcome_combined (NR/R)`, data = total_virome_R_NR_total_phages_w0_UC_caudo_non_crass, alternative = "two.sided",  paired = FALSE, exact = F)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$remission <- total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$`Endoscopic_outcome_combined (NR/R)`
effect_size <- wilcox_effsize(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass,value~remission)
(effect_size$effsize)
effect_size$n1+effect_size$n2
# P-value = 0.9096
# R-squared: 0.0001463437
#  n1: 48
#  n2: 46

# Caudo (crAss)
table(phyloseq_rarefied_phages_ggplot_class_melt_CD$`Observed richness`)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- phyloseq_rarefied_phages_ggplot_class_melt_CD[c(phyloseq_rarefied_phages_ggplot_class_melt_CD$Week == "primary endpoint" & phyloseq_rarefied_phages_ggplot_class_melt_CD$Diagnosis == "CD" & phyloseq_rarefied_phages_ggplot_class_melt_CD$`Observed richness` == "Caudoviricetes (CrAss)"),]
table(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$`Endoscopic_outcome_combined (NR/R)`) ## 22 R patients/ 22 NR patients
wilcox.test(value~`Endoscopic_outcome_combined (NR/R)`, data = total_virome_R_NR_total_phages_w0_UC_caudo_non_crass, alternative = "two.sided",  paired = FALSE, exact = F)
p.adjust(0.2638, method = "BH", n = 3)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$remission <- total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$`Endoscopic_outcome_combined (NR/R)`
effect_size <- wilcox_effsize(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass,value~remission)
(effect_size$effsize)
effect_size$n1
effect_size$n2
# P-value = 0.2638
# R-squared: 0.01337725
#  n1: 48
#  n2: 46

# Malgrandaviricetes
table(phyloseq_rarefied_phages_ggplot_class_melt_CD$`Observed richness`)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- phyloseq_rarefied_phages_ggplot_class_melt_CD[c(phyloseq_rarefied_phages_ggplot_class_melt_CD$Week == "primary endpoint" & phyloseq_rarefied_phages_ggplot_class_melt_CD$Diagnosis == "CD" & phyloseq_rarefied_phages_ggplot_class_melt_CD$`Observed richness` == "Malgrandaviricetes"),]
table(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$`Endoscopic_outcome_combined (NR/R)`) ## 22 R patients/ 22 NR patients
wilcox.test(value~`Endoscopic_outcome_combined (NR/R)`, data = total_virome_R_NR_total_phages_w0_UC_caudo_non_crass, alternative = "two.sided",  paired = FALSE, exact = F)
p.adjust(0.2015, method = "BH", n = 3)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$remission <- total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$`Endoscopic_outcome_combined (NR/R)`
effect_size <- wilcox_effsize(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass,value~remission)
(effect_size$effsize)
effect_size$n1+effect_size$n2
# P-value = 0.2015
# R-squared: 0.0328366
#  n1: 48
#  n2: 46

# Total phages
table(phyloseq_rarefied_phages_ggplot_class_melt_CD$`Observed richness`)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- phyloseq_rarefied_phages_ggplot_class_melt_CD[c(phyloseq_rarefied_phages_ggplot_class_melt_CD$Week == "primary endpoint" & phyloseq_rarefied_phages_ggplot_class_melt_CD$Diagnosis == "CD" & phyloseq_rarefied_phages_ggplot_class_melt_CD$`Observed richness` == "phages"),]
table(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$`Endoscopic_outcome_combined (NR/R)`) ## 22 R patients/ 22 NR patients
wilcox.test(value~`Endoscopic_outcome_combined (NR/R)`, data = total_virome_R_NR_total_phages_w0_UC_caudo_non_crass, alternative = "two.sided",  paired = FALSE, exact = F)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$remission <- total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$`Endoscopic_outcome_combined (NR/R)`
effect_size <- wilcox_effsize(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass,value~remission)
(effect_size$effsize)^2
effect_size$n1
effect_size$n2
# P-value = 0.5525
# R-squared: 0.003800871
#  n1: 48
#  n2: 46
####################################
# 4.6.1.6 Statistics between timepoints
####################################
# Statistics between baseline and primary endpoint
# Statistics within Remission
# Caudo (non-crAss)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- phyloseq_rarefied_phages_ggplot_class_melt_CD[c(phyloseq_rarefied_phages_ggplot_class_melt_CD$`Endoscopic_outcome_combined (NR/R)` == "remission" & phyloseq_rarefied_phages_ggplot_class_melt_CD$Diagnosis == "CD" & phyloseq_rarefied_phages_ggplot_class_melt_CD$`Observed richness` == "Caudoviricetes (non-CrAss)"),]
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired <- total_virome_R_NR_total_phages_w0_UC_caudo_non_crass
# View(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired)

colnames_de <- colnames(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired)
de1 <-data.frame("6365",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("12868",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("12943",NA, NA,"baseline",NA,NA,NA,NA)

#View(de1)
names(de1)<-c(colnames_de)
unpaired_paired_test <- rbind(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired, de1)
unpaired_paired_test_baseline <- unpaired_paired_test[unpaired_paired_test$Week == "baseline",]
unpaired_paired_test_baseline <- unpaired_paired_test_baseline[order(unpaired_paired_test_baseline$patient_ID) , ]
unpaired_paired_test_baseline$value <- as.numeric(unpaired_paired_test_baseline$value)
unpaired_paired_test_PE <- unpaired_paired_test[!unpaired_paired_test$Week == "baseline",]
unpaired_paired_test_PE <- unpaired_paired_test_PE[order(unpaired_paired_test_PE$patient_ID) , ]
unpaired_paired_test_PE$value <- as.numeric(unpaired_paired_test_PE$value)

#View(unpaired_paired_test)
#View(unpaired_paired_test_baseline)
#View(unpaired_paired_test_PE)

pm.wilcox.test(unpaired_paired_test_PE$value,unpaired_paired_test_baseline$value,alternative = "two.sided", correct = F, method = "MW-MW")

# P-value = 0.8285 (x: baseline, y: PE)
# Xpaired: 43
# Ypaired: 43
# Xextra: 0
# Y-extra: 3

# Caudo (crAss)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- phyloseq_rarefied_phages_ggplot_class_melt_CD[c(phyloseq_rarefied_phages_ggplot_class_melt_CD$`Endoscopic_outcome_combined (NR/R)` == "remission" & phyloseq_rarefied_phages_ggplot_class_melt_CD$Diagnosis == "CD" & phyloseq_rarefied_phages_ggplot_class_melt_CD$`Observed richness` == "Caudoviricetes (CrAss)"),]
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired <- total_virome_R_NR_total_phages_w0_UC_caudo_non_crass
# View(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired)

colnames_de <- colnames(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired)
de1 <-data.frame("6365",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("12868",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("12943",NA, NA,"baseline",NA,NA,NA,NA)

#View(de1)
names(de1)<-c(colnames_de)
unpaired_paired_test <- rbind(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired, de1)
unpaired_paired_test_baseline <- unpaired_paired_test[unpaired_paired_test$Week == "baseline",]
unpaired_paired_test_baseline <- unpaired_paired_test_baseline[order(unpaired_paired_test_baseline$patient_ID) , ]
unpaired_paired_test_baseline$value <- as.numeric(unpaired_paired_test_baseline$value)
unpaired_paired_test_PE <- unpaired_paired_test[!unpaired_paired_test$Week == "baseline",]
unpaired_paired_test_PE <- unpaired_paired_test_PE[order(unpaired_paired_test_PE$patient_ID) , ]
unpaired_paired_test_PE$value <- as.numeric(unpaired_paired_test_PE$value)

#View(unpaired_paired_test)
#View(unpaired_paired_test_baseline)
#View(unpaired_paired_test_PE)

pm.wilcox.test(unpaired_paired_test_PE$value,unpaired_paired_test_baseline$value,alternative = "two.sided", correct = F, method = "MW-MW")

# P-value = 0.615 (x: baseline, y: PE)
# Xpaired: 43
# Ypaired: 43
# Xextra: 0
# Y-extra: 3

# Caudo (Micro)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- phyloseq_rarefied_phages_ggplot_class_melt_CD[c(phyloseq_rarefied_phages_ggplot_class_melt_CD$`Endoscopic_outcome_combined (NR/R)` == "remission" & phyloseq_rarefied_phages_ggplot_class_melt_CD$Diagnosis == "CD" & phyloseq_rarefied_phages_ggplot_class_melt_CD$`Observed richness` == "Malgrandaviricetes"),]
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired <- total_virome_R_NR_total_phages_w0_UC_caudo_non_crass
# View(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired)

colnames_de <- colnames(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired)
de1 <-data.frame("6365",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("12868",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("12943",NA, NA,"baseline",NA,NA,NA,NA)

#View(de1)
names(de1)<-c(colnames_de)
unpaired_paired_test <- rbind(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired, de1)
unpaired_paired_test_baseline <- unpaired_paired_test[unpaired_paired_test$Week == "baseline",]
unpaired_paired_test_baseline <- unpaired_paired_test_baseline[order(unpaired_paired_test_baseline$patient_ID) , ]
unpaired_paired_test_baseline$value <- as.numeric(unpaired_paired_test_baseline$value)
unpaired_paired_test_PE <- unpaired_paired_test[!unpaired_paired_test$Week == "baseline",]
unpaired_paired_test_PE <- unpaired_paired_test_PE[order(unpaired_paired_test_PE$patient_ID) , ]
unpaired_paired_test_PE$value <- as.numeric(unpaired_paired_test_PE$value)

#View(unpaired_paired_test)
#View(unpaired_paired_test_baseline)
#View(unpaired_paired_test_PE)

pm.wilcox.test(unpaired_paired_test_PE$value,unpaired_paired_test_baseline$value,alternative = "two.sided", correct = F, method = "MW-MW")

# P-value = 0.2629 (x: baseline, y: PE)
# Xpaired: 43
# Ypaired: 43
# Xextra: 0
# Y-extra: 3

# phages
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- phyloseq_rarefied_phages_ggplot_class_melt_CD[c(phyloseq_rarefied_phages_ggplot_class_melt_CD$`Endoscopic_outcome_combined (NR/R)` == "remission" & phyloseq_rarefied_phages_ggplot_class_melt_CD$Diagnosis == "CD" & phyloseq_rarefied_phages_ggplot_class_melt_CD$`Observed richness` == "phages"),]
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired <- total_virome_R_NR_total_phages_w0_UC_caudo_non_crass
# View(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired)

colnames_de <- colnames(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired)
de1 <-data.frame("6365",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("12868",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("12943",NA, NA,"baseline",NA,NA,NA,NA)

#View(de1)
names(de1)<-c(colnames_de)
unpaired_paired_test <- rbind(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired, de1)
unpaired_paired_test_baseline <- unpaired_paired_test[unpaired_paired_test$Week == "baseline",]
unpaired_paired_test_baseline <- unpaired_paired_test_baseline[order(unpaired_paired_test_baseline$patient_ID) , ]
unpaired_paired_test_baseline$value <- as.numeric(unpaired_paired_test_baseline$value)
unpaired_paired_test_PE <- unpaired_paired_test[!unpaired_paired_test$Week == "baseline",]
unpaired_paired_test_PE <- unpaired_paired_test_PE[order(unpaired_paired_test_PE$patient_ID) , ]
unpaired_paired_test_PE$value <- as.numeric(unpaired_paired_test_PE$value)

#View(unpaired_paired_test)
#View(unpaired_paired_test_baseline)
#View(unpaired_paired_test_PE)

pm.wilcox.test(unpaired_paired_test_PE$value,unpaired_paired_test_baseline$value,alternative = "two.sided", correct = F, method = "MW-MW")

# P-value = 0.7833 (x: baseline, y: PE)
# Xpaired: 43
# Ypaired: 43
# Xextra: 0
# Y-extra: 3

# Statistics within non-Remission
# Caudo (non-crAss)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- phyloseq_rarefied_phages_ggplot_class_melt_CD[c(phyloseq_rarefied_phages_ggplot_class_melt_CD$`Endoscopic_outcome_combined (NR/R)` == "non-remission" & phyloseq_rarefied_phages_ggplot_class_melt_CD$Diagnosis == "CD" & phyloseq_rarefied_phages_ggplot_class_melt_CD$`Observed richness` == "Caudoviricetes (non-CrAss)"),]
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired <- total_virome_R_NR_total_phages_w0_UC_caudo_non_crass
#View(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired)

colnames_de <- colnames(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired)
de1 <-data.frame("3709",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("5364",NA, NA,"primary endpoint",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("5479",NA, NA,"primary endpoint",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("7265",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("9172",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("10409",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("10777",NA, NA,"primary endpoint",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("11040",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("11169",NA, NA,"primary endpoint",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("11180",NA, NA,"primary endpoint",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("12308",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("12405",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("12946",NA, NA,"primary endpoint",NA,NA,NA,NA)

#View(de1)
names(de1)<-c(colnames_de)
unpaired_paired_test <- rbind(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired, de1)
unpaired_paired_test$patient_ID[unpaired_paired_test$patient_ID == "6716" & unpaired_paired_test$Therapy_2 == "TNF"] <- "1"
unpaired_paired_test_baseline <- unpaired_paired_test[unpaired_paired_test$Week == "baseline",]
unpaired_paired_test_baseline <- unpaired_paired_test_baseline[order(unpaired_paired_test_baseline$patient_ID) , ]
unpaired_paired_test_baseline$value <- as.numeric(unpaired_paired_test_baseline$value)
unpaired_paired_test_PE <- unpaired_paired_test[!unpaired_paired_test$Week == "baseline",]
unpaired_paired_test_PE <- unpaired_paired_test_PE[order(unpaired_paired_test_PE$patient_ID) , ]
unpaired_paired_test_PE$value <- as.numeric(unpaired_paired_test_PE$value)

#View(unpaired_paired_test)
#View(unpaired_paired_test_baseline)
#View(unpaired_paired_test_PE)

pm.wilcox.test(unpaired_paired_test_PE$value,unpaired_paired_test_baseline$value,alternative = "two.sided", correct = F, method = "MW-MW")

# P-value = 0.797 (x: baseline, y: PE)
# Xpaired: 41
# Ypaired: 41
# Xextra: 7
# Y-extra: 6

# Caudo (crAss)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- phyloseq_rarefied_phages_ggplot_class_melt_CD[c(phyloseq_rarefied_phages_ggplot_class_melt_CD$`Endoscopic_outcome_combined (NR/R)` == "non-remission" & phyloseq_rarefied_phages_ggplot_class_melt_CD$Diagnosis == "CD" & phyloseq_rarefied_phages_ggplot_class_melt_CD$`Observed richness` == "Caudoviricetes (CrAss)"),]
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired <- total_virome_R_NR_total_phages_w0_UC_caudo_non_crass
# View(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired)

colnames_de <- colnames(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired)
de1 <-data.frame("3709",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("5364",NA, NA,"primary endpoint",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("5479",NA, NA,"primary endpoint",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("7265",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("9172",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("10409",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("10777",NA, NA,"primary endpoint",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("11040",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("11169",NA, NA,"primary endpoint",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("11180",NA, NA,"primary endpoint",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("12308",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("12405",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("12946",NA, NA,"primary endpoint",NA,NA,NA,NA)

#View(de1)
names(de1)<-c(colnames_de)
unpaired_paired_test <- rbind(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired, de1)
unpaired_paired_test$patient_ID[unpaired_paired_test$patient_ID == "6716" & unpaired_paired_test$Therapy_2 == "TNF"] <- "1"
unpaired_paired_test_baseline <- unpaired_paired_test[unpaired_paired_test$Week == "baseline",]
unpaired_paired_test_baseline <- unpaired_paired_test_baseline[order(unpaired_paired_test_baseline$patient_ID) , ]
unpaired_paired_test_baseline$value <- as.numeric(unpaired_paired_test_baseline$value)
unpaired_paired_test_PE <- unpaired_paired_test[!unpaired_paired_test$Week == "baseline",]
unpaired_paired_test_PE <- unpaired_paired_test_PE[order(unpaired_paired_test_PE$patient_ID) , ]
unpaired_paired_test_PE$value <- as.numeric(unpaired_paired_test_PE$value)
#View(unpaired_paired_test)
#View(unpaired_paired_test_baseline)
#View(unpaired_paired_test_PE)

pm.wilcox.test(unpaired_paired_test_PE$value,unpaired_paired_test_baseline$value,alternative = "two.sided", correct = F, method = "MW-MW")

# P-value = 0.7835 (x: baseline, y: PE)
# Xpaired: 43
# Ypaired: 43
# Xextra: 0
# Y-extra: 3

# Caudo (Micro)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- phyloseq_rarefied_phages_ggplot_class_melt_CD[c(phyloseq_rarefied_phages_ggplot_class_melt_CD$`Endoscopic_outcome_combined (NR/R)` == "non-remission" & phyloseq_rarefied_phages_ggplot_class_melt_CD$Diagnosis == "CD" & phyloseq_rarefied_phages_ggplot_class_melt_CD$`Observed richness` == "Malgrandaviricetes"),]
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired <- total_virome_R_NR_total_phages_w0_UC_caudo_non_crass
# View(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired)

colnames_de <- colnames(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired)
de1 <-data.frame("3709",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("5364",NA, NA,"primary endpoint",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("5479",NA, NA,"primary endpoint",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("7265",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("9172",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("10409",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("10777",NA, NA,"primary endpoint",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("11040",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("11169",NA, NA,"primary endpoint",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("11180",NA, NA,"primary endpoint",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("12308",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("12405",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("12946",NA, NA,"primary endpoint",NA,NA,NA,NA)

#View(de1)
names(de1)<-c(colnames_de)
unpaired_paired_test <- rbind(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired, de1)
unpaired_paired_test$patient_ID[unpaired_paired_test$patient_ID == "6716" & unpaired_paired_test$Therapy_2 == "TNF"] <- "1"
unpaired_paired_test_baseline <- unpaired_paired_test[unpaired_paired_test$Week == "baseline",]
unpaired_paired_test_baseline <- unpaired_paired_test_baseline[order(unpaired_paired_test_baseline$patient_ID) , ]
unpaired_paired_test_baseline$value <- as.numeric(unpaired_paired_test_baseline$value)
unpaired_paired_test_PE <- unpaired_paired_test[!unpaired_paired_test$Week == "baseline",]
unpaired_paired_test_PE <- unpaired_paired_test_PE[order(unpaired_paired_test_PE$patient_ID) , ]
unpaired_paired_test_PE$value <- as.numeric(unpaired_paired_test_PE$value)

#View(unpaired_paired_test)
#View(unpaired_paired_test_baseline)
#View(unpaired_paired_test_PE)

pm.wilcox.test(unpaired_paired_test_PE$value,unpaired_paired_test_baseline$value,alternative = "two.sided", correct = F, method = "MW-MW")

# P-value = 0.8268 (x: baseline, y: PE)
# Xpaired: 43
# Ypaired: 43
# Xextra: 0
# Y-extra: 3

# phages
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- phyloseq_rarefied_phages_ggplot_class_melt_CD[c(phyloseq_rarefied_phages_ggplot_class_melt_CD$`Endoscopic_outcome_combined (NR/R)` == "non-remission" & phyloseq_rarefied_phages_ggplot_class_melt_CD$Diagnosis == "CD" & phyloseq_rarefied_phages_ggplot_class_melt_CD$`Observed richness` == "phages"),]
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired <- total_virome_R_NR_total_phages_w0_UC_caudo_non_crass
# View(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired)

colnames_de <- colnames(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired)
de1 <-data.frame("3709",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("5364",NA, NA,"primary endpoint",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("5479",NA, NA,"primary endpoint",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("7265",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("9172",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("10409",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("10777",NA, NA,"primary endpoint",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("11040",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("11169",NA, NA,"primary endpoint",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("11180",NA, NA,"primary endpoint",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("12308",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("12405",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("12946",NA, NA,"primary endpoint",NA,NA,NA,NA)

#View(de1)
names(de1)<-c(colnames_de)
unpaired_paired_test <- rbind(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired, de1)
unpaired_paired_test$patient_ID[unpaired_paired_test$patient_ID == "6716" & unpaired_paired_test$Therapy_2 == "TNF"] <- "1"
unpaired_paired_test_baseline <- unpaired_paired_test[unpaired_paired_test$Week == "baseline",]
unpaired_paired_test_baseline <- unpaired_paired_test_baseline[order(unpaired_paired_test_baseline$patient_ID) , ]
unpaired_paired_test_baseline$value <- as.numeric(unpaired_paired_test_baseline$value)
unpaired_paired_test_PE <- unpaired_paired_test[!unpaired_paired_test$Week == "baseline",]
unpaired_paired_test_PE <- unpaired_paired_test_PE[order(unpaired_paired_test_PE$patient_ID) , ]
unpaired_paired_test_PE$value <- as.numeric(unpaired_paired_test_PE$value)

#View(unpaired_paired_test)
#View(unpaired_paired_test_baseline)
#View(unpaired_paired_test_PE)

pm.wilcox.test(unpaired_paired_test_PE$value,unpaired_paired_test_baseline$value,alternative = "two.sided", correct = F, method = "MW-MW")

# P-value = 0.8992 (x: baseline, y: PE)
# Xpaired: 43
# Ypaired: 43
# Xextra: 0
# Y-extra: 3
####################################
# 4.6.2 Richness boxplots between disease (UC/CD) at primary endpoint to compare inflammed and non-inflammed gut
####################################
# 4.6.2.1 Richness within Primary endpoint between disease
####################################
phyloseq_rarefied_phages_ggplot_class_melt_CD <- phyloseq_rarefied_phages_ggplot_class_melt
phyloseq_rarefied_phages_ggplot_class_melt_CD$Week[phyloseq_rarefied_phages_ggplot_class_melt_CD$Week == "w0"] <- "baseline"
phyloseq_rarefied_phages_ggplot_class_melt_CD$Week[c(phyloseq_rarefied_phages_ggplot_class_melt_CD$Week == "w24" & phyloseq_rarefied_phages_ggplot_class_melt_CD$Diagnosis == "CD") | c(phyloseq_rarefied_phages_ggplot_class_melt_CD$Week == "w14" & phyloseq_rarefied_phages_ggplot_class_melt_CD$Diagnosis == "UC")] <- "primary endpoint"
phyloseq_rarefied_phages_ggplot_class_melt_CD <- phyloseq_rarefied_phages_ggplot_class_melt_CD[!phyloseq_rarefied_phages_ggplot_class_melt_CD$Week == "w14",]
phyloseq_rarefied_phages_ggplot_class_melt_CD <- phyloseq_rarefied_phages_ggplot_class_melt_CD[!phyloseq_rarefied_phages_ggplot_class_melt_CD$`Endoscopic_outcome_combined (NR/R)` == "unknown",]
phyloseq_rarefied_phages_ggplot_class_melt_CD <- phyloseq_rarefied_phages_ggplot_class_melt_CD[!phyloseq_rarefied_phages_ggplot_class_melt_CD$Week == "baseline",]
phyloseq_rarefied_phages_ggplot_class_melt_CD$`Observed richness` <- factor(phyloseq_rarefied_phages_ggplot_class_melt_CD$`Observed richness`, levels=c("phages", "Caudoviricetes (CrAss)","Caudoviricetes (non-CrAss)", "Malgrandaviricetes"))
phyloseq_rarefied_phages_ggplot_class_melt_CD$Diagnosis <- factor(phyloseq_rarefied_phages_ggplot_class_melt_CD$Diagnosis, levels=c("UC","CD"))
table(phyloseq_rarefied_phages_ggplot_class_melt_CD$`Observed richness`)
#View(phyloseq_rarefied_phages_ggplot_class_melt_CD)

phyloseq_rarefied_phages_ggplot_class_melt_CD %>%
 filter(value >= 0L & value <= 80L) %>%
 ggplot() +
 aes(x = `Endoscopic_outcome_combined (NR/R)`, 
 y = value, fill = `Observed richness`) +
 geom_boxplot(shape = "circle") +
 theme_bw() +
 facet_wrap(vars(Diagnosis)) +
 xlab("") +
 ylab("Observed Richness") +
 scale_fill_manual(values=c("#669bbc","#59b031", "#ee766f", "#c49b05"),name="Class") +
 scale_y_continuous(breaks=c(0,20,40,60,80), limits = c(0,80)) 
####################################
# 4.6.2.2 between remission between CD -UC
####################################
phyloseq_rarefied_phages_ggplot_class_melt_CD

# Statistics within UC
# Caudo (non-crAss)
table(phyloseq_rarefied_phages_ggplot_class_melt_CD$`Observed richness`)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- phyloseq_rarefied_phages_ggplot_class_melt_CD[c(phyloseq_rarefied_phages_ggplot_class_melt_CD$`Endoscopic_outcome_combined (NR/R)` == "non-remission" & phyloseq_rarefied_phages_ggplot_class_melt_CD$`Observed richness` == "Caudoviricetes (non-CrAss)"),]
table(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$`Endoscopic_outcome_combined (NR/R)`) ## 33 R patients/ 27 NR patients
wilcox.test(value~Diagnosis, data = total_virome_R_NR_total_phages_w0_UC_caudo_non_crass, alternative = "two.sided",  paired = FALSE, exact = F)
p.adjust(0.2466, method = "BH", n =3)
effect_size <- wilcox_effsize(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass,value~Diagnosis)
(effect_size$effsize)
effect_size$n1+effect_size$n2
# P-value = 0.2361

# Caudo (crAss)
table(phyloseq_rarefied_phages_ggplot_class_melt_CD$`Observed richness`)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- phyloseq_rarefied_phages_ggplot_class_melt_CD[c(phyloseq_rarefied_phages_ggplot_class_melt_CD$`Endoscopic_outcome_combined (NR/R)` == "non-remission" & phyloseq_rarefied_phages_ggplot_class_melt_CD$`Observed richness` == "Caudoviricetes (CrAss)"),]
table(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$`Endoscopic_outcome_combined (NR/R)`) ## 33 R patients/ 27 NR patients
wilcox.test(value~Diagnosis, data = total_virome_R_NR_total_phages_w0_UC_caudo_non_crass, alternative = "two.sided",  paired = FALSE, exact = F)
p.adjust(0.6883, method = "BH", n =3)
effect_size <- wilcox_effsize(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass,value~Diagnosis)
(effect_size$effsize)
effect_size$n1+effect_size$n2
# P-value = 0.6883

# Caudo (Micro)
table(phyloseq_rarefied_phages_ggplot_class_melt_CD$`Observed richness`)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- phyloseq_rarefied_phages_ggplot_class_melt_CD[c(phyloseq_rarefied_phages_ggplot_class_melt_CD$`Endoscopic_outcome_combined (NR/R)` == "non-remission" & phyloseq_rarefied_phages_ggplot_class_melt_CD$`Observed richness` == "Malgrandaviricetes"),]
table(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$`Endoscopic_outcome_combined (NR/R)`) ## 33 R patients/ 27 NR patients
wilcox.test(value~Diagnosis, data = total_virome_R_NR_total_phages_w0_UC_caudo_non_crass, alternative = "two.sided",  paired = FALSE, exact = F)
p.adjust(0.248, method = "BH", n =3)
effect_size <- wilcox_effsize(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass,value~Diagnosis)
(effect_size$effsize)
effect_size$n1+effect_size$n2
# P-value = 0.008126*****
# R-squared: 0.09133822
#  n1: 31
#  n2: 46

# Total phages
table(phyloseq_rarefied_phages_ggplot_class_melt_CD$`Observed richness`)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- phyloseq_rarefied_phages_ggplot_class_melt_CD[c(phyloseq_rarefied_phages_ggplot_class_melt_CD$`Endoscopic_outcome_combined (NR/R)` == "remission" & phyloseq_rarefied_phages_ggplot_class_melt_CD$`Observed richness` == "phages"),]
table(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$`Endoscopic_outcome_combined (NR/R)`) ## 33 R patients/ 27 NR patients
wilcox.test(value~Diagnosis, data = total_virome_R_NR_total_phages_w0_UC_caudo_non_crass, alternative = "two.sided",  paired = FALSE, exact = F)
effect_size <- wilcox_effsize(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass,value~Diagnosis, alternative = "two.sided",  paired = F, exact = F)
(effect_size$effsize)^2
effect_size$n1
effect_size$n2

# P-value = 0.0472*****
# R-squared: 0.05141276
#  n1: 31
#  n2: 46
####################################
# 4.6.2.3 between non-remission between CD -UC
####################################
phyloseq_rarefied_phages_ggplot_class_melt_CD

# Statistics within UC
# Caudo (non-crAss)
table(phyloseq_rarefied_phages_ggplot_class_melt_CD$`Observed richness`)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- phyloseq_rarefied_phages_ggplot_class_melt_CD[c(phyloseq_rarefied_phages_ggplot_class_melt_CD$`Endoscopic_outcome_combined (NR/R)` == "non-remission" & phyloseq_rarefied_phages_ggplot_class_melt_CD$`Observed richness` == "Caudoviricetes (non-CrAss)"),]
table(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$`Endoscopic_outcome_combined (NR/R)`) ## 33 R patients/ 27 NR patients
wilcox.test(value~Diagnosis, data = total_virome_R_NR_total_phages_w0_UC_caudo_non_crass, alternative = "two.sided",  paired = FALSE, exact = F)
# P-value = 0.7414

# Caudo (crAss)
table(phyloseq_rarefied_phages_ggplot_class_melt_CD$`Observed richness`)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- phyloseq_rarefied_phages_ggplot_class_melt_CD[c(phyloseq_rarefied_phages_ggplot_class_melt_CD$`Endoscopic_outcome_combined (NR/R)` == "non-remission" & phyloseq_rarefied_phages_ggplot_class_melt_CD$`Observed richness` == "Caudoviricetes (CrAss)"),]
table(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$`Endoscopic_outcome_combined (NR/R)`) ## 33 R patients/ 27 NR patients
wilcox.test(value~Diagnosis, data = total_virome_R_NR_total_phages_w0_UC_caudo_non_crass, alternative = "two.sided",  paired = FALSE, exact = F)
# P-value = 0.77

# Caudo (Micro)
table(phyloseq_rarefied_phages_ggplot_class_melt_CD$`Observed richness`)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- phyloseq_rarefied_phages_ggplot_class_melt_CD[c(phyloseq_rarefied_phages_ggplot_class_melt_CD$`Endoscopic_outcome_combined (NR/R)` == "non-remission" & phyloseq_rarefied_phages_ggplot_class_melt_CD$`Observed richness` == "Malgrandaviricetes"),]
table(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$`Endoscopic_outcome_combined (NR/R)`) ## 33 R patients/ 27 NR patients
wilcox.test(value~Diagnosis, data = total_virome_R_NR_total_phages_w0_UC_caudo_non_crass, alternative = "two.sided",  paired = FALSE, exact = F)
effect_size <- wilcox_effsize(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass,value~Diagnosis, alternative = "two.sided",  paired = F, exact = F)
(effect_size$effsize)^2
effect_size$n1
effect_size$n2
# P-value = 0.2476
# R-squared: 0.09133822
#  n1: 31
#  n2: 46

# Total phages
table(phyloseq_rarefied_phages_ggplot_class_melt_CD$`Observed richness`)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- phyloseq_rarefied_phages_ggplot_class_melt_CD[c(phyloseq_rarefied_phages_ggplot_class_melt_CD$`Endoscopic_outcome_combined (NR/R)` == "non-remission" & phyloseq_rarefied_phages_ggplot_class_melt_CD$`Observed richness` == "phages"),]
table(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$`Endoscopic_outcome_combined (NR/R)`) ## 33 R patients/ 27 NR patients
wilcox.test(value~Diagnosis, data = total_virome_R_NR_total_phages_w0_UC_caudo_non_crass, alternative = "two.sided",  paired = FALSE, exact = F)
effect_size <- wilcox_effsize(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass,value~Diagnosis, alternative = "two.sided",  paired = F, exact = F)
(effect_size$effsize)^2
effect_size$n1
effect_size$n2

# P-value = 0.3669
# R-squared: 0.05141276
#  n1: 31
#  n2: 46
####################################
# 4.6.3 Shannon-diversity boxplots
####################################
# 4.6.3.1 UC
####################################
phyloseq_rarefied_phages_ggplot_class_melt_UC <- Mastertable_unrarefied_observed_shannon_DF_melt[Mastertable_unrarefied_observed_shannon_DF_melt$Diagnosis == "UC",]
phyloseq_rarefied_phages_ggplot_class_melt_UC$Week[phyloseq_rarefied_phages_ggplot_class_melt_UC$Week == "w0"] <- "baseline"
phyloseq_rarefied_phages_ggplot_class_melt_UC$Week[phyloseq_rarefied_phages_ggplot_class_melt_UC$Week == "w14"] <- "primary endpoint"
phyloseq_rarefied_phages_ggplot_class_melt_UC$`Shannon diversity` <- factor(phyloseq_rarefied_phages_ggplot_class_melt_UC$`Shannon diversity`, levels=c("phages", "Caudoviricetes (CrAss)","Caudoviricetes (non-CrAss)", "Malgrandaviricetes"))
phyloseq_rarefied_phages_ggplot_class_melt_UC <- phyloseq_rarefied_phages_ggplot_class_melt_UC[phyloseq_rarefied_phages_ggplot_class_melt_UC$Week == "primary endpoint",]
table(phyloseq_rarefied_phages_ggplot_class_melt_UC$`Shannon diversity`)

View(phyloseq_rarefied_phages_ggplot_class_melt_UC)

# Figure
phyloseq_rarefied_phages_ggplot_class_melt_UC %>%
  filter(`Shannon diversity` %in% c("Caudoviricetes (CrAss)", "Caudoviricetes (non-CrAss)", "Malgrandaviricetes")) %>%
  ggplot() +
  aes(x = `Shannon diversity`, y = value, fill = `Endoscopic_outcome_combined (NR/R)`) +
  geom_boxplot() +
  scale_fill_manual(values = c(`non-remission` = "#AD4B38", remission = "#49AD86")) +
 theme_bw() +
 xlab("") +
 ylab("Shannon diversity") +
 theme(legend.position = "none")  
####################################
# 4.6.3.2 Statistics within timepoints
####################################
# Statistics within baseline
# Caudo (non-crAss)
table(phyloseq_rarefied_phages_ggplot_class_melt_UC$`Shannon diversity`)
View(phyloseq_rarefied_phages_ggplot_class_melt_UC)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- phyloseq_rarefied_phages_ggplot_class_melt_UC[c(phyloseq_rarefied_phages_ggplot_class_melt_UC$Week == "primary endpoint" & phyloseq_rarefied_phages_ggplot_class_melt_UC$Diagnosis == "UC" & phyloseq_rarefied_phages_ggplot_class_melt_UC$`Shannon diversity` == "Caudoviricetes (non-CrAss)"),]
wilcox.test(value~`Endoscopic_outcome_combined (NR/R)`, data = total_virome_R_NR_total_phages_w0_UC_caudo_non_crass, alternative = "two.sided",  paired = FALSE, exact = F)
p.adjust(0.00477, method = "BH", n = 8)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$remission <- total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$`Endoscopic_outcome_combined (NR/R)`
effect_size <- wilcox_effsize(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass,value~remission)
(effect_size$effsize)
effect_size$n1+effect_size$n2
# P-value = 0.7861
# R-squared: 0.001395576
#  n1: 26
#  n2: 30

# Caudo (crAss)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- phyloseq_rarefied_phages_ggplot_class_melt_UC[c(phyloseq_rarefied_phages_ggplot_class_melt_UC$Week == "primary endpoint" & phyloseq_rarefied_phages_ggplot_class_melt_UC$Diagnosis == "UC" & phyloseq_rarefied_phages_ggplot_class_melt_UC$`Shannon diversity` == "Caudoviricetes (CrAss)"),]
wilcox.test(value~`Endoscopic_outcome_combined (NR/R)`, data = total_virome_R_NR_total_phages_w0_UC_caudo_non_crass, alternative = "two.sided",  paired = FALSE, exact = F)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$remission <- total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$`Endoscopic_outcome_combined (NR/R)`
effect_size <- wilcox_effsize(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass,value~remission)
(effect_size$effsize)
effect_size$n1+effect_size$n2
p.adjust(0.2895, method = "BH", n = 8)
# P-value = 0.2895
# R-squared: 0.02038195
#  n1: 26
#  n2: 30

# Caudo (Micro)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- phyloseq_rarefied_phages_ggplot_class_melt_UC[c(phyloseq_rarefied_phages_ggplot_class_melt_UC$Week == "primary endpoint" & phyloseq_rarefied_phages_ggplot_class_melt_UC$Diagnosis == "UC" & phyloseq_rarefied_phages_ggplot_class_melt_UC$`Shannon diversity` == "Malgrandaviricetes"),]
wilcox.test(value~`Endoscopic_outcome_combined (NR/R)`, data = total_virome_R_NR_total_phages_w0_UC_caudo_non_crass, alternative = "two.sided",  paired = FALSE, exact = F)
p.adjust(0.8496, method = "BH", n = 8)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$remission <- total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$`Endoscopic_outcome_combined (NR/R)`
effect_size <- wilcox_effsize(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass,value~remission)
(effect_size$effsize)
effect_size$n1+effect_size$n2
# P-value = 0.8496
# R-squared: 0.0006993007
#  n1: 26
#  n2: 30

# Total phages
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- phyloseq_rarefied_phages_ggplot_class_melt_UC[c(phyloseq_rarefied_phages_ggplot_class_melt_UC$Week == "baseline" & phyloseq_rarefied_phages_ggplot_class_melt_UC$Diagnosis == "UC" & phyloseq_rarefied_phages_ggplot_class_melt_UC$`Shannon diversity` == "phages"),]
wilcox.test(value~`Endoscopic_outcome_combined (NR/R)`, data = total_virome_R_NR_total_phages_w0_UC_caudo_non_crass, alternative = "two.sided",  paired = FALSE, exact = F)
p.adjust(.9672, method = "BH", n = 8)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$remission <- total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$`Endoscopic_outcome_combined (NR/R)`
effect_size <- wilcox_effsize(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass,value~remission)
(effect_size$effsize)
effect_size$n1+effect_size$n2
# P-value = 0.9672
# R-squared: 4.33925e-05 
#  n1: 26
#  n2: 30

# Statistics within primary endpoint
# Caudo (non-crAss)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- phyloseq_rarefied_phages_ggplot_class_melt_UC[c(phyloseq_rarefied_phages_ggplot_class_melt_UC$Week == "primary endpoint" & phyloseq_rarefied_phages_ggplot_class_melt_UC$Diagnosis == "UC" & phyloseq_rarefied_phages_ggplot_class_melt_UC$`Shannon diversity` == "Caudoviricetes (non-CrAss)"),]
View(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass)
wilcox.test(value~`Endoscopic_outcome_combined (NR/R)`, data = total_virome_R_NR_total_phages_w0_UC_caudo_non_crass, alternative = "two.sided",  paired = FALSE, exact = F)
p.adjust(0.00477, method = "BH", n = 3)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$remission <- total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$`Endoscopic_outcome_combined (NR/R)`
effect_size <- wilcox_effsize(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass,value~remission)
(effect_size$effsize)
effect_size$n1+effect_size$n2
# P-value = 0.00477*****
# R-squared: 0.1380856
#  n1: 26
#  n2: 32

# Caudo (crAss)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- phyloseq_rarefied_phages_ggplot_class_melt_UC[c(phyloseq_rarefied_phages_ggplot_class_melt_UC$Week == "primary endpoint" & phyloseq_rarefied_phages_ggplot_class_melt_UC$Diagnosis == "UC" & phyloseq_rarefied_phages_ggplot_class_melt_UC$`Shannon diversity` == "Caudoviricetes (CrAss)"),]
wilcox.test(value~`Endoscopic_outcome_combined (NR/R)`, data = total_virome_R_NR_total_phages_w0_UC_caudo_non_crass, alternative = "two.sided",  paired = FALSE, exact = F)
p.adjust(0.2587, method = "BH", n = 3)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$remission <- total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$`Endoscopic_outcome_combined (NR/R)`
effect_size <- wilcox_effsize(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass,value~remission)
(effect_size$effsize)
effect_size$n1+effect_size$n2
# P-value = 0.2587
# R-squared: 0.02232674
#  n1: 26
#  n2: 32

# Caudo (Micro)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- phyloseq_rarefied_phages_ggplot_class_melt_UC[c(phyloseq_rarefied_phages_ggplot_class_melt_UC$Week == "primary endpoint" & phyloseq_rarefied_phages_ggplot_class_melt_UC$Diagnosis == "UC" & phyloseq_rarefied_phages_ggplot_class_melt_UC$`Shannon diversity` == "Malgrandaviricetes"),]
wilcox.test(value~`Endoscopic_outcome_combined (NR/R)`, data = total_virome_R_NR_total_phages_w0_UC_caudo_non_crass, alternative = "two.sided",  paired = FALSE, exact = F)
p.adjust(0.01261, method = "BH", n = 3)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$remission <- total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$`Endoscopic_outcome_combined (NR/R)`
effect_size <- wilcox_effsize(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass,value~remission)
(effect_size$effsize)
effect_size$n1+effect_size$n2
# P-value = 0.01261****
# R-squared: 0.1079655
#  n1: 26
#  n2: 32

# Total phages
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- phyloseq_rarefied_phages_ggplot_class_melt_UC[c(phyloseq_rarefied_phages_ggplot_class_melt_UC$Week == "primary endpoint" & phyloseq_rarefied_phages_ggplot_class_melt_UC$Diagnosis == "UC" & phyloseq_rarefied_phages_ggplot_class_melt_UC$`Shannon diversity` == "phages"),]
wilcox.test(value~`Endoscopic_outcome_combined (NR/R)`, data = total_virome_R_NR_total_phages_w0_UC_caudo_non_crass, alternative = "two.sided",  paired = FALSE, exact = F)
p.adjust(0.1235, method = "BH", n = 4)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$remission <- total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$`Endoscopic_outcome_combined (NR/R)`
effect_size <- wilcox_effsize(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass,value~remission)
(effect_size$effsize)
effect_size$n1+effect_size$n2
# P-value = 0.1235
# R-squared: 0.04130934
#  n1: 26
#  n2: 32
####################################
# 4.6.3.3 Statistics between timepoints
####################################
# Statistics between baseline and primary endpoint
# Statistics within Remission
# Caudo (non-crAss)
View(phyloseq_rarefied_phages_ggplot_class_melt_UC)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- phyloseq_rarefied_phages_ggplot_class_melt_UC[c(phyloseq_rarefied_phages_ggplot_class_melt_UC$`Endoscopic_outcome_combined (NR/R)` == "remission" & phyloseq_rarefied_phages_ggplot_class_melt_UC$Diagnosis == "UC" & phyloseq_rarefied_phages_ggplot_class_melt_UC$`Shannon diversity` == "Caudoviricetes (non-CrAss)"),]
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired <- total_virome_R_NR_total_phages_w0_UC_caudo_non_crass
View(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired)

colnames_de <- colnames(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired)
de1 <-data.frame("10213",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("12255",NA, NA,"baseline",NA,NA,NA,NA)

#View(de1)
names(de1)<-c(colnames_de)
unpaired_paired_test <- rbind(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired, de1)
unpaired_paired_test_baseline <- unpaired_paired_test[unpaired_paired_test$Week == "baseline",]
unpaired_paired_test_baseline <- unpaired_paired_test_baseline[order(unpaired_paired_test_baseline$patient_ID) , ]
unpaired_paired_test_baseline$value <- as.numeric(unpaired_paired_test_baseline$value)
unpaired_paired_test_PE <- unpaired_paired_test[!unpaired_paired_test$Week == "baseline",]
unpaired_paired_test_PE <- unpaired_paired_test_PE[order(unpaired_paired_test_PE$patient_ID) , ]
unpaired_paired_test_PE$value <- as.numeric(unpaired_paired_test_PE$value)

pm.wilcox.test(unpaired_paired_test_PE$value,unpaired_paired_test_baseline$value,alternative = "two.sided", correct = F, method = "MW-MW")
p.adjust(0.03691, method = "BH", n = 8)

unpaired_paired_test_PE_baseline <- rbind.fill(unpaired_paired_test_PE,unpaired_paired_test_baseline)
unpaired_paired_test_PE_baseline <- unpaired_paired_test_PE_baseline[order(unpaired_paired_test_PE_baseline$patient_ID),]
effect_size <- wilcox_effsize(unpaired_paired_test_PE_baseline,value~Week,paired = TRUE, alternative = "two.sided")
(effect_size$effsize)
effect_size$n1+effect_size$n2
# P-value = 0.03691***** (x: baseline, y: PE)
# Xpaired: 30
# Ypaired: 30
# Xextra: 0
# Y-extra: 2

# Caudo (crAss)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- phyloseq_rarefied_phages_ggplot_class_melt_UC[c(phyloseq_rarefied_phages_ggplot_class_melt_UC$`Endoscopic_outcome_combined (NR/R)` == "remission" & phyloseq_rarefied_phages_ggplot_class_melt_UC$Diagnosis == "UC" & phyloseq_rarefied_phages_ggplot_class_melt_UC$`Shannon diversity` == "Caudoviricetes (CrAss)"),]
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired <- total_virome_R_NR_total_phages_w0_UC_caudo_non_crass
#View(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired)

colnames_de <- colnames(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired)
de1 <-data.frame("10213",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("12255",NA, NA,"baseline",NA,NA,NA,NA)

#View(de1)
names(de1)<-c(colnames_de)
unpaired_paired_test <- rbind(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired, de1)
unpaired_paired_test_baseline <- unpaired_paired_test[unpaired_paired_test$Week == "baseline",]
unpaired_paired_test_baseline <- unpaired_paired_test_baseline[order(unpaired_paired_test_baseline$patient_ID) , ]
unpaired_paired_test_baseline$value <- as.numeric(unpaired_paired_test_baseline$value)
unpaired_paired_test_PE <- unpaired_paired_test[!unpaired_paired_test$Week == "baseline",]
unpaired_paired_test_PE <- unpaired_paired_test_PE[order(unpaired_paired_test_PE$patient_ID) , ]
unpaired_paired_test_PE$value <- as.numeric(unpaired_paired_test_PE$value)

pm.wilcox.test(unpaired_paired_test_PE$value,unpaired_paired_test_baseline$value,alternative = "two.sided", correct = F, method = "MW-MW")
p.adjust(0.71, method = "BH", n = 8)

unpaired_paired_test_PE_baseline <- rbind.fill(unpaired_paired_test_PE,unpaired_paired_test_baseline)
unpaired_paired_test_PE_baseline <- unpaired_paired_test_PE_baseline[order(unpaired_paired_test_PE_baseline$patient_ID),]
effect_size <- wilcox_effsize(unpaired_paired_test_PE_baseline,value~Week,paired = TRUE, alternative = "two.sided")
(effect_size$effsize)
effect_size$n1+effect_size$n2

# P-value = 0.71 (x: baseline, y: PE)
# Xpaired: 30
# Ypaired: 30
# Xextra: 0
# Y-extra: 2

# Micro
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- phyloseq_rarefied_phages_ggplot_class_melt_UC[c(phyloseq_rarefied_phages_ggplot_class_melt_UC$`Endoscopic_outcome_combined (NR/R)` == "remission" & phyloseq_rarefied_phages_ggplot_class_melt_UC$Diagnosis == "UC" & phyloseq_rarefied_phages_ggplot_class_melt_UC$`Shannon diversity` == "Malgrandaviricetes"),]
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired <- total_virome_R_NR_total_phages_w0_UC_caudo_non_crass
#View(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired)

colnames_de <- colnames(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired)
de1 <-data.frame("10213",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("12255",NA, NA,"baseline",NA,NA,NA,NA)

#View(de1)
names(de1)<-c(colnames_de)
unpaired_paired_test <- rbind(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired, de1)
unpaired_paired_test_baseline <- unpaired_paired_test[unpaired_paired_test$Week == "baseline",]
unpaired_paired_test_baseline <- unpaired_paired_test_baseline[order(unpaired_paired_test_baseline$patient_ID) , ]
unpaired_paired_test_baseline$value <- as.numeric(unpaired_paired_test_baseline$value)
unpaired_paired_test_PE <- unpaired_paired_test[!unpaired_paired_test$Week == "baseline",]
unpaired_paired_test_PE <- unpaired_paired_test_PE[order(unpaired_paired_test_PE$patient_ID) , ]
unpaired_paired_test_PE$value <- as.numeric(unpaired_paired_test_PE$value)

pm.wilcox.test(unpaired_paired_test_PE$value,unpaired_paired_test_baseline$value,alternative = "two.sided", correct = F, method = "MW-MW")
p.adjust(0.004983, method = "BH", n = 8)

unpaired_paired_test_PE_baseline <- rbind.fill(unpaired_paired_test_PE,unpaired_paired_test_baseline)
unpaired_paired_test_PE_baseline <- unpaired_paired_test_PE_baseline[order(unpaired_paired_test_PE_baseline$patient_ID),]
effect_size <- wilcox_effsize(unpaired_paired_test_PE_baseline,value~Week,paired = TRUE, alternative = "two.sided")
(effect_size$effsize)
effect_size$n1+effect_size$n2

# P-value = 0.004983***** (x: baseline, y: PE)
# Xpaired: 30
# Ypaired: 30
# Xextra: 0
# Y-extra: 2

# Phages
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- phyloseq_rarefied_phages_ggplot_class_melt_UC[c(phyloseq_rarefied_phages_ggplot_class_melt_UC$`Endoscopic_outcome_combined (NR/R)` == "remission" & phyloseq_rarefied_phages_ggplot_class_melt_UC$Diagnosis == "UC" & phyloseq_rarefied_phages_ggplot_class_melt_UC$`Shannon diversity` == "phages"),]
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired <- total_virome_R_NR_total_phages_w0_UC_caudo_non_crass
#View(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired)

colnames_de <- colnames(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired)
de1 <-data.frame("10213",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("12255",NA, NA,"baseline",NA,NA,NA,NA)

#View(de1)
names(de1)<-c(colnames_de)
unpaired_paired_test <- rbind(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired, de1)
unpaired_paired_test_baseline <- unpaired_paired_test[unpaired_paired_test$Week == "baseline",]
unpaired_paired_test_baseline <- unpaired_paired_test_baseline[order(unpaired_paired_test_baseline$patient_ID) , ]
unpaired_paired_test_baseline$value <- as.numeric(unpaired_paired_test_baseline$value)
unpaired_paired_test_PE <- unpaired_paired_test[!unpaired_paired_test$Week == "baseline",]
unpaired_paired_test_PE <- unpaired_paired_test_PE[order(unpaired_paired_test_PE$patient_ID) , ]
unpaired_paired_test_PE$value <- as.numeric(unpaired_paired_test_PE$value)

pm.wilcox.test(unpaired_paired_test_PE$value,unpaired_paired_test_baseline$value,alternative = "two.sided", correct = F, method = "MW-MW")
p.adjust(0.005092, method = "BH", n = 8)

unpaired_paired_test_PE_baseline <- rbind.fill(unpaired_paired_test_PE,unpaired_paired_test_baseline)
unpaired_paired_test_PE_baseline <- unpaired_paired_test_PE_baseline[order(unpaired_paired_test_PE_baseline$patient_ID),]
effect_size <- wilcox_effsize(unpaired_paired_test_PE_baseline,value~Week,paired = TRUE, alternative = "two.sided")
(effect_size$effsize)
effect_size$n1+effect_size$n2

# P-value = 0.005092***** (x: baseline, y: PE)
# Xpaired: 30
# Ypaired: 30
# Xextra: 0
# Y-extra: 2

# Statistics within non-Remission
# Caudo (non-crAss)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- phyloseq_rarefied_phages_ggplot_class_melt_UC[c(phyloseq_rarefied_phages_ggplot_class_melt_UC$`Endoscopic_outcome_combined (NR/R)` == "non-remission" & phyloseq_rarefied_phages_ggplot_class_melt_UC$Diagnosis == "UC" & phyloseq_rarefied_phages_ggplot_class_melt_UC$`Shannon diversity` == "Caudoviricetes (non-CrAss)"),]
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired <- total_virome_R_NR_total_phages_w0_UC_caudo_non_crass
#View(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired)

colnames_de <- colnames(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired)
de1 <-data.frame("8798",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("10449",NA, NA,"primary endpoint",NA,NA,NA,NA)

#View(de1)
names(de1)<-c(colnames_de)
unpaired_paired_test <- rbind(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired, de1)
unpaired_paired_test$patient_ID[c(unpaired_paired_test$patient_ID == "13093" & unpaired_paired_test$Therapy_2 == "TNF")] <- "2"
unpaired_paired_test$patient_ID[c(unpaired_paired_test$patient_ID == "12770" & unpaired_paired_test$Therapy_2 == "TNF")] <- "1"
unpaired_paired_test$patient_ID[c(unpaired_paired_test$patient_ID == "9516" & unpaired_paired_test$Therapy_2 == "TNF")] <- "3"

unpaired_paired_test_baseline <- unpaired_paired_test[unpaired_paired_test$Week == "baseline",]
unpaired_paired_test_baseline <- unpaired_paired_test_baseline[order(unpaired_paired_test_baseline$patient_ID) , ]
unpaired_paired_test_baseline$value <- as.numeric(unpaired_paired_test_baseline$value)
unpaired_paired_test_PE <- unpaired_paired_test[!unpaired_paired_test$Week == "baseline",]
unpaired_paired_test_PE <- unpaired_paired_test_PE[order(unpaired_paired_test_PE$patient_ID) , ]
unpaired_paired_test_PE$value <- as.numeric(unpaired_paired_test_PE$value)

pm.wilcox.test(unpaired_paired_test_PE$value,unpaired_paired_test_baseline$value,alternative = "two.sided", correct = F, method = "MW-MW")
p.adjust(0.11, method = "BH", n = 8)

unpaired_paired_test_PE_baseline <- rbind.fill(unpaired_paired_test_PE,unpaired_paired_test_baseline)
unpaired_paired_test_PE_baseline <- unpaired_paired_test_PE_baseline[order(unpaired_paired_test_PE_baseline$patient_ID),]
effect_size <- wilcox_effsize(unpaired_paired_test_PE_baseline,value~Week,paired = TRUE, alternative = "two.sided")
(effect_size$effsize)
effect_size$n1+effect_size$n2

id <-unique(unpaired_paired_test_PE_baseline$patient_ID)
id
View(unpaired_paired_test_PE_baseline)
# P-value = 0.11 (x: baseline, y: PE)
# Xpaired: 25
# Ypaired: 25
# Xextra: 1
# Y-extra: 1

# Caudo (crAss)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- phyloseq_rarefied_phages_ggplot_class_melt_UC[c(phyloseq_rarefied_phages_ggplot_class_melt_UC$`Endoscopic_outcome_combined (NR/R)` == "non-remission" & phyloseq_rarefied_phages_ggplot_class_melt_UC$Diagnosis == "UC" & phyloseq_rarefied_phages_ggplot_class_melt_UC$`Shannon diversity` == "Caudoviricetes (CrAss)"),]
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired <- total_virome_R_NR_total_phages_w0_UC_caudo_non_crass
#View(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired)

colnames_de <- colnames(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired)
de1 <-data.frame("8798",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("10449",NA, NA,"primary endpoint",NA,NA,NA,NA)

#View(de1)
names(de1)<-c(colnames_de)
unpaired_paired_test <- rbind(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired, de1)
unpaired_paired_test$patient_ID[c(unpaired_paired_test$patient_ID == "12770" & unpaired_paired_test$Therapy_2 == "TNF")] <- "1"
unpaired_paired_test$patient_ID[c(unpaired_paired_test$patient_ID == "13093" & unpaired_paired_test$Therapy_2 == "TNF")] <- "2"
unpaired_paired_test$patient_ID[c(unpaired_paired_test$patient_ID == "9516" & unpaired_paired_test$Therapy_2 == "TNF")] <- "3"

unpaired_paired_test_baseline <- unpaired_paired_test[unpaired_paired_test$Week == "baseline",]
unpaired_paired_test_baseline <- unpaired_paired_test_baseline[order(unpaired_paired_test_baseline$patient_ID) , ]
unpaired_paired_test_baseline$value <- as.numeric(unpaired_paired_test_baseline$value)
unpaired_paired_test_PE <- unpaired_paired_test[!unpaired_paired_test$Week == "baseline",]
unpaired_paired_test_PE <- unpaired_paired_test_PE[order(unpaired_paired_test_PE$patient_ID) , ]
unpaired_paired_test_PE$value <- as.numeric(unpaired_paired_test_PE$value)

pm.wilcox.test(unpaired_paired_test_PE$value,unpaired_paired_test_baseline$value,alternative = "two.sided", correct = F, method = "MW-MW")
p.adjust(0.8051, method = "BH", n = 8)

unpaired_paired_test_PE_baseline <- rbind.fill(unpaired_paired_test_PE,unpaired_paired_test_baseline)
unpaired_paired_test_PE_baseline <- unpaired_paired_test_PE_baseline[order(unpaired_paired_test_PE_baseline$patient_ID),]
effect_size <- wilcox_effsize(unpaired_paired_test_PE_baseline,value~Week,paired = TRUE, alternative = "two.sided")
(effect_size$effsize)
effect_size$n1+effect_size$n2

# P-value = 0.8051 (x: baseline, y: PE)
# Xpaired: 25
# Ypaired: 25
# Xextra: 1
# Y-extra: 1

# Micro
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- phyloseq_rarefied_phages_ggplot_class_melt_UC[c(phyloseq_rarefied_phages_ggplot_class_melt_UC$`Endoscopic_outcome_combined (NR/R)` == "non-remission" & phyloseq_rarefied_phages_ggplot_class_melt_UC$Diagnosis == "UC" & phyloseq_rarefied_phages_ggplot_class_melt_UC$`Shannon diversity` == "Malgrandaviricetes"),]
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired <- total_virome_R_NR_total_phages_w0_UC_caudo_non_crass
#View(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired)

colnames_de <- colnames(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired)
de1 <-data.frame("8798",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("10449",NA, NA,"primary endpoint",NA,NA,NA,NA)

#View(de1)
names(de1)<-c(colnames_de)
unpaired_paired_test <- rbind(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired, de1)
unpaired_paired_test$patient_ID[c(unpaired_paired_test$patient_ID == "12770" & unpaired_paired_test$Therapy_2 == "TNF")] <- "1"
unpaired_paired_test$patient_ID[c(unpaired_paired_test$patient_ID == "13093" & unpaired_paired_test$Therapy_2 == "TNF")] <- "2"
unpaired_paired_test$patient_ID[c(unpaired_paired_test$patient_ID == "9516" & unpaired_paired_test$Therapy_2 == "TNF")] <- "3"

unpaired_paired_test_baseline <- unpaired_paired_test[unpaired_paired_test$Week == "baseline",]
unpaired_paired_test_baseline <- unpaired_paired_test_baseline[order(unpaired_paired_test_baseline$patient_ID) , ]
unpaired_paired_test_baseline$value <- as.numeric(unpaired_paired_test_baseline$value)
unpaired_paired_test_PE <- unpaired_paired_test[!unpaired_paired_test$Week == "baseline",]
unpaired_paired_test_PE <- unpaired_paired_test_PE[order(unpaired_paired_test_PE$patient_ID) , ]
unpaired_paired_test_PE$value <- as.numeric(unpaired_paired_test_PE$value)

pm.wilcox.test(unpaired_paired_test_PE$value,unpaired_paired_test_baseline$value,alternative = "two.sided", correct = F, method = "MW-MW")
p.adjust(0.8595, method = "BH", n = 8)

unpaired_paired_test_PE_baseline <- rbind.fill(unpaired_paired_test_PE,unpaired_paired_test_baseline)
unpaired_paired_test_PE_baseline <- unpaired_paired_test_PE_baseline[order(unpaired_paired_test_PE_baseline$patient_ID),]
effect_size <- wilcox_effsize(unpaired_paired_test_PE_baseline,value~Week,paired = TRUE, alternative = "two.sided")
(effect_size$effsize)
effect_size$n1+effect_size$n2
View(unpaired_paired_test_PE_baseline)
# P-value = 0.8595 (x: baseline, y: PE)
# Xpaired: 25
# Ypaired: 25
# Xextra: 1
# Y-extra: 1

# Phages
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- phyloseq_rarefied_phages_ggplot_class_melt_UC[c(phyloseq_rarefied_phages_ggplot_class_melt_UC$`Endoscopic_outcome_combined (NR/R)` == "non-remission" & phyloseq_rarefied_phages_ggplot_class_melt_UC$Diagnosis == "UC" & phyloseq_rarefied_phages_ggplot_class_melt_UC$`Shannon diversity` == "phages"),]
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired <- total_virome_R_NR_total_phages_w0_UC_caudo_non_crass
#View(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired)

colnames_de <- colnames(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired)
de1 <-data.frame("8798",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("10449",NA, NA,"primary endpoint",NA,NA,NA,NA)

#View(de1)
names(de1)<-c(colnames_de)
unpaired_paired_test <- rbind(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired, de1)
unpaired_paired_test$patient_ID[c(unpaired_paired_test$patient_ID == "12770" & unpaired_paired_test$Therapy_2 == "TNF")] <- "1"
unpaired_paired_test$patient_ID[c(unpaired_paired_test$patient_ID == "13093" & unpaired_paired_test$Therapy_2 == "TNF")] <- "2"
unpaired_paired_test$patient_ID[c(unpaired_paired_test$patient_ID == "9516" & unpaired_paired_test$Therapy_2 == "TNF")] <- "3"

unpaired_paired_test_baseline <- unpaired_paired_test[unpaired_paired_test$Week == "baseline",]
unpaired_paired_test_baseline <- unpaired_paired_test_baseline[order(unpaired_paired_test_baseline$patient_ID) , ]
unpaired_paired_test_baseline$value <- as.numeric(unpaired_paired_test_baseline$value)
unpaired_paired_test_PE <- unpaired_paired_test[!unpaired_paired_test$Week == "baseline",]
unpaired_paired_test_PE <- unpaired_paired_test_PE[order(unpaired_paired_test_PE$patient_ID) , ]
unpaired_paired_test_PE$value <- as.numeric(unpaired_paired_test_PE$value)

pm.wilcox.test(unpaired_paired_test_PE$value,unpaired_paired_test_baseline$value,alternative = "two.sided", correct = F, method = "MW-MW")
p.adjust(0.736, method = "BH", n = 8)

unpaired_paired_test_PE_baseline <- rbind.fill(unpaired_paired_test_PE,unpaired_paired_test_baseline)
unpaired_paired_test_PE_baseline <- unpaired_paired_test_PE_baseline[order(unpaired_paired_test_PE_baseline$patient_ID),]
effect_size <- wilcox_effsize(unpaired_paired_test_PE_baseline,value~Week,paired = TRUE, alternative = "two.sided")
(effect_size$effsize)
effect_size$n1+effect_size$n2

# P-value = 0.736 (x: baseline, y: PE)
# Xpaired: 25
# Ypaired: 25
# Xextra: 1
# Y-extra: 1
####################################
# 4.6.3.4 CD
####################################
phyloseq_rarefied_phages_ggplot_class_melt_UC <- Mastertable_unrarefied_observed_shannon_DF_melt[Mastertable_unrarefied_observed_shannon_DF_melt$Diagnosis == "CD",]
phyloseq_rarefied_phages_ggplot_class_melt_UC$Week[phyloseq_rarefied_phages_ggplot_class_melt_UC$Week == "w0"] <- "baseline"
phyloseq_rarefied_phages_ggplot_class_melt_UC$Week[phyloseq_rarefied_phages_ggplot_class_melt_UC$Week == "w24"] <- "primary endpoint"
phyloseq_rarefied_phages_ggplot_class_melt_UC <- phyloseq_rarefied_phages_ggplot_class_melt_UC[!phyloseq_rarefied_phages_ggplot_class_melt_UC$Week == "w14",] 
phyloseq_rarefied_phages_ggplot_class_melt_UC <- phyloseq_rarefied_phages_ggplot_class_melt_UC[!phyloseq_rarefied_phages_ggplot_class_melt_UC$`Endoscopic_outcome_combined (NR/R)` == "unknown",]
phyloseq_rarefied_phages_ggplot_class_melt_UC$`Shannon diversity` <- factor(phyloseq_rarefied_phages_ggplot_class_melt_UC$`Shannon diversity`, levels=c("phages", "Caudoviricetes (CrAss)","Caudoviricetes (non-CrAss)", "Malgrandaviricetes"))
table(phyloseq_rarefied_phages_ggplot_class_melt_UC$`Shannon diversity`)
#View(phyloseq_rarefied_phages_ggplot_class_melt_UC)

# Figure:
phyloseq_rarefied_phages_ggplot_class_melt_UC %>%
  filter(`Shannon diversity` %in% c("Caudoviricetes (CrAss)", "Caudoviricetes (non-CrAss)", "Malgrandaviricetes")) %>%
  ggplot() +
  aes(x = `Shannon diversity`, y = value, fill = `Endoscopic_outcome_combined (NR/R)`) +
  geom_boxplot() +
  scale_fill_manual(values = c(`non-remission` = "#AD4B38", remission = "#49AD86")) +
  theme_bw() +
  xlab("") +
  ylab("Shannon diversity") +
  theme(legend.position = "none")  

## RA
## Richness
## Pilou's evenesss lastly then put udnerneath each other & claculate significance
## Shannon 

####################################
# 4.6.3.5 Statistics within timepoints
####################################
# Statistics within baseline
# Caudo (non-crAss)
table(phyloseq_rarefied_phages_ggplot_class_melt_UC$`Shannon diversity`)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- phyloseq_rarefied_phages_ggplot_class_melt_UC[c(phyloseq_rarefied_phages_ggplot_class_melt_UC$Week == "baseline" & phyloseq_rarefied_phages_ggplot_class_melt_UC$Diagnosis == "CD" & phyloseq_rarefied_phages_ggplot_class_melt_UC$`Shannon diversity` == "Caudoviricetes (non-CrAss)"),]
wilcox.test(value~`Endoscopic_outcome_combined (NR/R)`, data = total_virome_R_NR_total_phages_w0_UC_caudo_non_crass, alternative = "two.sided",  paired = FALSE, exact = F)
p.adjust(0.7693, method = "BH", n = 8)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$remission <- total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$`Endoscopic_outcome_combined (NR/R)`
effect_size <- wilcox_effsize(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass,value~remission)
(effect_size$effsize)
effect_size$n1+effect_size$n2
# P-value = 0.7693
# R-squared: 0.001395576
#  n1: 26
#  n2: 30

# Caudo (crAss)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- phyloseq_rarefied_phages_ggplot_class_melt_UC[c(phyloseq_rarefied_phages_ggplot_class_melt_UC$Week == "baseline" & phyloseq_rarefied_phages_ggplot_class_melt_UC$Diagnosis == "CD" & phyloseq_rarefied_phages_ggplot_class_melt_UC$`Shannon diversity` == "Caudoviricetes (CrAss)"),]
wilcox.test(value~`Endoscopic_outcome_combined (NR/R)`, data = total_virome_R_NR_total_phages_w0_UC_caudo_non_crass, alternative = "two.sided",  paired = FALSE, exact = F)
p.adjust(0.08787, method = "BH", n = 8)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$remission <- total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$`Endoscopic_outcome_combined (NR/R)`
effect_size <- wilcox_effsize(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass,value~remission)
(effect_size$effsize)
effect_size$n1+effect_size$n2
# P-value = 0.08787
# R-squared: 0.02038195
#  n1: 26
#  n2: 30

# Caudo (Micro)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- phyloseq_rarefied_phages_ggplot_class_melt_UC[c(phyloseq_rarefied_phages_ggplot_class_melt_UC$Week == "baseline" & phyloseq_rarefied_phages_ggplot_class_melt_UC$Diagnosis == "CD" & phyloseq_rarefied_phages_ggplot_class_melt_UC$`Shannon diversity` == "Malgrandaviricetes"),]
wilcox.test(value~`Endoscopic_outcome_combined (NR/R)`, data = total_virome_R_NR_total_phages_w0_UC_caudo_non_crass, alternative = "two.sided",  paired = FALSE, exact = F)
p.adjust(0.04074, method = "BH", n = 8)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$remission <- total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$`Endoscopic_outcome_combined (NR/R)`
effect_size <- wilcox_effsize(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass,value~remission)
(effect_size$effsize)
effect_size$n1+effect_size$n2
# P-value = 0.04074*****
# R-squared: 0.04159534
#  n1: 55
#  n2: 46

# Total phages
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- phyloseq_rarefied_phages_ggplot_class_melt_UC[c(phyloseq_rarefied_phages_ggplot_class_melt_UC$Week == "baseline" & phyloseq_rarefied_phages_ggplot_class_melt_UC$Diagnosis == "CD" & phyloseq_rarefied_phages_ggplot_class_melt_UC$`Shannon diversity` == "phages"),]
wilcox.test(value~`Endoscopic_outcome_combined (NR/R)`, data = total_virome_R_NR_total_phages_w0_UC_caudo_non_crass, alternative = "two.sided",  paired = FALSE, exact = F)
p.adjust(0.6874, method = "BH", n = 8)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$remission <- total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$`Endoscopic_outcome_combined (NR/R)`
effect_size <- wilcox_effsize(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass,value~remission)
(effect_size$effsize)
effect_size$n1+effect_size$n2
# P-value = 0.6874
# R-squared: 4.33925e-05 
#  n1: 26
#  n2: 30

# Statistics within primary endpoint
# Caudo (non-crAss)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- phyloseq_rarefied_phages_ggplot_class_melt_UC[c(phyloseq_rarefied_phages_ggplot_class_melt_UC$Week == "primary endpoint" & phyloseq_rarefied_phages_ggplot_class_melt_UC$Diagnosis == "CD" & phyloseq_rarefied_phages_ggplot_class_melt_UC$`Shannon diversity` == "Caudoviricetes (non-CrAss)"),]
wilcox.test(value~`Endoscopic_outcome_combined (NR/R)`, data = total_virome_R_NR_total_phages_w0_UC_caudo_non_crass, alternative = "two.sided",  paired = FALSE, exact = F)
p.adjust(0.6866, method = "BH", n = 3)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$remission <- total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$`Endoscopic_outcome_combined (NR/R)`
effect_size <- wilcox_effsize(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass,value~remission)
(effect_size$effsize)
effect_size$n1+effect_size$n2
# P-value: 0.6866
# R-squared: 0.1380856
#  n1: 26
#  n2: 32

# Caudo (crAss)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- phyloseq_rarefied_phages_ggplot_class_melt_UC[c(phyloseq_rarefied_phages_ggplot_class_melt_UC$Week == "primary endpoint" & phyloseq_rarefied_phages_ggplot_class_melt_UC$Diagnosis == "CD" & phyloseq_rarefied_phages_ggplot_class_melt_UC$`Shannon diversity` == "Caudoviricetes (CrAss)"),]
wilcox.test(value~`Endoscopic_outcome_combined (NR/R)`, data = total_virome_R_NR_total_phages_w0_UC_caudo_non_crass, alternative = "two.sided",  paired = FALSE, exact = F)
p.adjust(0.1708, method = "BH", n = 3)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$remission <- total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$`Endoscopic_outcome_combined (NR/R)`
effect_size <- wilcox_effsize(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass,value~remission)
(effect_size$effsize)
effect_size$n1+effect_size$n2
# P-value = 0.1708
# R-squared: 0.02232674
#  n1: 26
#  n2: 32

# Caudo (Micro)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- phyloseq_rarefied_phages_ggplot_class_melt_UC[c(phyloseq_rarefied_phages_ggplot_class_melt_UC$Week == "primary endpoint" & phyloseq_rarefied_phages_ggplot_class_melt_UC$Diagnosis == "CD" & phyloseq_rarefied_phages_ggplot_class_melt_UC$`Shannon diversity` == "Malgrandaviricetes"),]
wilcox.test(value~`Endoscopic_outcome_combined (NR/R)`, data = total_virome_R_NR_total_phages_w0_UC_caudo_non_crass, alternative = "two.sided",  paired = FALSE, exact = F)
p.adjust(0.1057, method = "BH", n = 3)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$remission <- total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$`Endoscopic_outcome_combined (NR/R)`
effect_size <- wilcox_effsize(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass,value~remission)
(effect_size$effsize)
effect_size$n1+effect_size$n2
# P-value = 0.1057
# R-squared: 0.1079655
#  n1: 26
#  n2: 32

# Total phages
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- phyloseq_rarefied_phages_ggplot_class_melt_UC[c(phyloseq_rarefied_phages_ggplot_class_melt_UC$Week == "primary endpoint" & phyloseq_rarefied_phages_ggplot_class_melt_UC$Diagnosis == "CD" & phyloseq_rarefied_phages_ggplot_class_melt_UC$`Shannon diversity` == "phages"),]
wilcox.test(value~`Endoscopic_outcome_combined (NR/R)`, data = total_virome_R_NR_total_phages_w0_UC_caudo_non_crass, alternative = "two.sided",  paired = FALSE, exact = F)
p.adjust(0.4984, method = "BH", n = 8)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$remission <- total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$`Endoscopic_outcome_combined (NR/R)`
effect_size <- wilcox_effsize(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass,value~remission)
(effect_size$effsize)
effect_size$n1+effect_size$n2
# P-value = 0.4984
# R-squared: 0.04130934
#  n1: 26
#  n2: 32
####################################
# 4.6.3.6 Statistics between timepoints
####################################
# Statistics between baseline and primary endpoint
# Statistics within Remission
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- phyloseq_rarefied_phages_ggplot_class_melt_UC[c(phyloseq_rarefied_phages_ggplot_class_melt_UC$`Endoscopic_outcome_combined (NR/R)` == "remission" & phyloseq_rarefied_phages_ggplot_class_melt_UC$Diagnosis == "CD" & phyloseq_rarefied_phages_ggplot_class_melt_UC$`Shannon diversity` == "Caudoviricetes (non-CrAss)"),]
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired <- total_virome_R_NR_total_phages_w0_UC_caudo_non_crass
#View(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired)

unpaired_paired_test_baseline <- total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired[total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired$Week == "baseline",]
unpaired_paired_test_baseline <- unpaired_paired_test_baseline[order(unpaired_paired_test_baseline$patient_ID) , ]
unpaired_paired_test_baseline$value <- as.numeric(unpaired_paired_test_baseline$value)
unpaired_paired_test_PE <- total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired[!total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired$Week == "baseline",]
unpaired_paired_test_PE <- unpaired_paired_test_PE[order(unpaired_paired_test_PE$patient_ID) , ]
unpaired_paired_test_PE$value <- as.numeric(unpaired_paired_test_PE$value)

pm.wilcox.test(unpaired_paired_test_PE$value,unpaired_paired_test_baseline$value,alternative = "two.sided", correct = F, method = "MW-MW")
p.adjust(0.9741, method = "BH", n = 8)

unpaired_paired_test_PE_baseline <- rbind.fill(unpaired_paired_test_PE,unpaired_paired_test_baseline)
unpaired_paired_test_PE_baseline <- unpaired_paired_test_PE_baseline[order(unpaired_paired_test_PE_baseline$patient_ID),]
effect_size <- wilcox_effsize(unpaired_paired_test_PE_baseline,value~Week,paired = TRUE, alternative = "two.sided")
(effect_size$effsize)
effect_size$n1+effect_size$n2
# P-value = 0.9741 (x: baseline, y: PE)
# Xpaired: 46
# Ypaired: 46
# Xextra: 0
# Y-extra: 0

# Caudo (crAss)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- phyloseq_rarefied_phages_ggplot_class_melt_UC[c(phyloseq_rarefied_phages_ggplot_class_melt_UC$`Endoscopic_outcome_combined (NR/R)` == "remission" & phyloseq_rarefied_phages_ggplot_class_melt_UC$Diagnosis == "CD" & phyloseq_rarefied_phages_ggplot_class_melt_UC$`Shannon diversity` == "Caudoviricetes (CrAss)"),]
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired <- total_virome_R_NR_total_phages_w0_UC_caudo_non_crass
#View(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired)

unpaired_paired_test_baseline <- total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired[total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired$Week == "baseline",]
unpaired_paired_test_baseline <- unpaired_paired_test_baseline[order(unpaired_paired_test_baseline$patient_ID) , ]
unpaired_paired_test_baseline$value <- as.numeric(unpaired_paired_test_baseline$value)
unpaired_paired_test_PE <- total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired[!total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired$Week == "baseline",]
unpaired_paired_test_PE <- unpaired_paired_test_PE[order(unpaired_paired_test_PE$patient_ID) , ]
unpaired_paired_test_PE$value <- as.numeric(unpaired_paired_test_PE$value)

pm.wilcox.test(unpaired_paired_test_PE$value,unpaired_paired_test_baseline$value,alternative = "two.sided", correct = F, method = "MW-MW")
p.adjust(0.416, method = "BH", n = 8)

unpaired_paired_test_PE_baseline <- rbind.fill(unpaired_paired_test_PE,unpaired_paired_test_baseline)
unpaired_paired_test_PE_baseline <- unpaired_paired_test_PE_baseline[order(unpaired_paired_test_PE_baseline$patient_ID),]
effect_size <- wilcox_effsize(unpaired_paired_test_PE_baseline,value~Week,paired = TRUE, alternative = "two.sided")
(effect_size$effsize)
effect_size$n1+effect_size$n2

# P-value = 0.416 (x: baseline, y: PE)
# Xpaired: 30
# Ypaired: 30
# Xextra: 0
# Y-extra: 2

# Micro
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- phyloseq_rarefied_phages_ggplot_class_melt_UC[c(phyloseq_rarefied_phages_ggplot_class_melt_UC$`Endoscopic_outcome_combined (NR/R)` == "remission" & phyloseq_rarefied_phages_ggplot_class_melt_UC$Diagnosis == "CD" & phyloseq_rarefied_phages_ggplot_class_melt_UC$`Shannon diversity` == "Malgrandaviricetes"),]
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired <- total_virome_R_NR_total_phages_w0_UC_caudo_non_crass
#View(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired)

unpaired_paired_test_baseline <- total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired[total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired$Week == "baseline",]
unpaired_paired_test_baseline <- unpaired_paired_test_baseline[order(unpaired_paired_test_baseline$patient_ID) , ]
unpaired_paired_test_baseline$value <- as.numeric(unpaired_paired_test_baseline$value)
unpaired_paired_test_PE <- total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired[!total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired$Week == "baseline",]
unpaired_paired_test_PE <- unpaired_paired_test_PE[order(unpaired_paired_test_PE$patient_ID) , ]
unpaired_paired_test_PE$value <- as.numeric(unpaired_paired_test_PE$value)

pm.wilcox.test(unpaired_paired_test_PE$value,unpaired_paired_test_baseline$value,alternative = "two.sided", correct = F, method = "MW-MW")
p.adjust(0.7264, method = "BH", n = 8)

unpaired_paired_test_PE_baseline <- rbind.fill(unpaired_paired_test_PE,unpaired_paired_test_baseline)
unpaired_paired_test_PE_baseline <- unpaired_paired_test_PE_baseline[order(unpaired_paired_test_PE_baseline$patient_ID),]
effect_size <- wilcox_effsize(unpaired_paired_test_PE_baseline,value~Week,paired = TRUE, alternative = "two.sided")
(effect_size$effsize)
effect_size$n1+effect_size$n2
# P-value = 0.7264 (x: baseline, y: PE)
# Xpaired: 46
# Ypaired: 46
# Xextra: 0
# Y-extra: 0

# Phages
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- phyloseq_rarefied_phages_ggplot_class_melt_UC[c(phyloseq_rarefied_phages_ggplot_class_melt_UC$`Endoscopic_outcome_combined (NR/R)` == "remission" & phyloseq_rarefied_phages_ggplot_class_melt_UC$Diagnosis == "CD" & phyloseq_rarefied_phages_ggplot_class_melt_UC$`Shannon diversity` == "phages"),]
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired <- total_virome_R_NR_total_phages_w0_UC_caudo_non_crass
#View(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired)

unpaired_paired_test_baseline <- total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired[total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired$Week == "baseline",]
unpaired_paired_test_baseline <- unpaired_paired_test_baseline[order(unpaired_paired_test_baseline$patient_ID) , ]
unpaired_paired_test_baseline$value <- as.numeric(unpaired_paired_test_baseline$value)
unpaired_paired_test_PE <- total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired[!total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired$Week == "baseline",]
unpaired_paired_test_PE <- unpaired_paired_test_PE[order(unpaired_paired_test_PE$patient_ID) , ]
unpaired_paired_test_PE$value <- as.numeric(unpaired_paired_test_PE$value)

pm.wilcox.test(unpaired_paired_test_PE$value,unpaired_paired_test_baseline$value,alternative = "two.sided", correct = F, method = "MW-MW")
p.adjust(0.3127, method = "BH", n = 8)

unpaired_paired_test_PE_baseline <- rbind.fill(unpaired_paired_test_PE,unpaired_paired_test_baseline)
unpaired_paired_test_PE_baseline <- unpaired_paired_test_PE_baseline[order(unpaired_paired_test_PE_baseline$patient_ID),]
effect_size <- wilcox_effsize(unpaired_paired_test_PE_baseline,value~Week,paired = TRUE, alternative = "two.sided")
(effect_size$effsize)
effect_size$n1+effect_size$n2
# P-value = 0.3127 (x: baseline, y: PE)
# Xpaired: 46
# Ypaired: 46
# Xextra: 0
# Y-extra: 0

# Statistics within non-Remission
# Caudo (non-crAss)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- phyloseq_rarefied_phages_ggplot_class_melt_UC[c(phyloseq_rarefied_phages_ggplot_class_melt_UC$`Endoscopic_outcome_combined (NR/R)` == "non-remission" & phyloseq_rarefied_phages_ggplot_class_melt_UC$Diagnosis == "CD" & phyloseq_rarefied_phages_ggplot_class_melt_UC$`Shannon diversity` == "Caudoviricetes (non-CrAss)"),]
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired <- total_virome_R_NR_total_phages_w0_UC_caudo_non_crass
#View(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired)

colnames_de <- colnames(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired)
de1 <-data.frame("10777",NA, NA,"primary endpoint",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("12011",NA, NA,"primary endpoint",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("12405",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("5479",NA, NA,"primary endpoint",NA,NA,NA,NA)

#View(de1)
names(de1)<-c(colnames_de)
unpaired_paired_test <- rbind(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired, de1)
unpaired_paired_test$patient_ID[c(unpaired_paired_test$patient_ID == "6716" & unpaired_paired_test$Therapy_2 == "TNF")] <- "1"
unpaired_paired_test_baseline <- unpaired_paired_test[unpaired_paired_test$Week == "baseline",]
unpaired_paired_test_baseline <- unpaired_paired_test_baseline[order(unpaired_paired_test_baseline$patient_ID) , ]
unpaired_paired_test_baseline$value <- as.numeric(unpaired_paired_test_baseline$value)
unpaired_paired_test_PE <- unpaired_paired_test[!unpaired_paired_test$Week == "baseline",]
unpaired_paired_test_PE <- unpaired_paired_test_PE[order(unpaired_paired_test_PE$patient_ID) , ]
unpaired_paired_test_PE$value <- as.numeric(unpaired_paired_test_PE$value)

pm.wilcox.test(unpaired_paired_test_PE$value,unpaired_paired_test_baseline$value,alternative = "two.sided", correct = F, method = "MW-MW")
p.adjust(0.6572, method = "BH", n = 8)

unpaired_paired_test_PE_baseline <- rbind.fill(unpaired_paired_test_PE,unpaired_paired_test_baseline)
unpaired_paired_test_PE_baseline <- unpaired_paired_test_PE_baseline[order(unpaired_paired_test_PE_baseline$patient_ID),]
effect_size <- wilcox_effsize(unpaired_paired_test_PE_baseline,value~Week,paired = TRUE, alternative = "two.sided")
(effect_size$effsize)
effect_size$n1+effect_size$n2

# P-value = 0.6572 (x: baseline, y: PE)
# Xpaired: 52
# Ypaired: 52
# Xextra: 1
# Y-extra: 3

# Caudo (crAss)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- phyloseq_rarefied_phages_ggplot_class_melt_UC[c(phyloseq_rarefied_phages_ggplot_class_melt_UC$`Endoscopic_outcome_combined (NR/R)` == "non-remission" & phyloseq_rarefied_phages_ggplot_class_melt_UC$Diagnosis == "CD" & phyloseq_rarefied_phages_ggplot_class_melt_UC$`Shannon diversity` == "Caudoviricetes (CrAss)"),]
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired <- total_virome_R_NR_total_phages_w0_UC_caudo_non_crass
#View(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired)

colnames_de <- colnames(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired)
de1 <-data.frame("10777",NA, NA,"primary endpoint",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("12011",NA, NA,"primary endpoint",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("12405",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("5479",NA, NA,"primary endpoint",NA,NA,NA,NA)

#View(de1)
names(de1)<-c(colnames_de)
unpaired_paired_test <- rbind(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired, de1)
unpaired_paired_test$patient_ID[c(unpaired_paired_test$patient_ID == "6716" & unpaired_paired_test$Therapy_2 == "TNF")] <- "1"
unpaired_paired_test_baseline <- unpaired_paired_test[unpaired_paired_test$Week == "baseline",]
unpaired_paired_test_baseline <- unpaired_paired_test_baseline[order(unpaired_paired_test_baseline$patient_ID) , ]
unpaired_paired_test_baseline$value <- as.numeric(unpaired_paired_test_baseline$value)
unpaired_paired_test_PE <- unpaired_paired_test[!unpaired_paired_test$Week == "baseline",]
unpaired_paired_test_PE <- unpaired_paired_test_PE[order(unpaired_paired_test_PE$patient_ID) , ]
unpaired_paired_test_PE$value <- as.numeric(unpaired_paired_test_PE$value)

pm.wilcox.test(unpaired_paired_test_PE$value,unpaired_paired_test_baseline$value,alternative = "two.sided", correct = F, method = "MW-MW")
p.adjust(0.2606, method = "BH", n = 8)

unpaired_paired_test_PE_baseline <- rbind.fill(unpaired_paired_test_PE,unpaired_paired_test_baseline)
unpaired_paired_test_PE_baseline <- unpaired_paired_test_PE_baseline[order(unpaired_paired_test_PE_baseline$patient_ID),]
effect_size <- wilcox_effsize(unpaired_paired_test_PE_baseline,value~Week,paired = TRUE, alternative = "two.sided")
(effect_size$effsize)
effect_size$n1+effect_size$n2
# P-value = 0.2606 (x: baseline, y: PE)*****
# Xpaired: 52
# Ypaired: 52
# Xextra: 1
# Y-extra: 3

# Micro
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- phyloseq_rarefied_phages_ggplot_class_melt_UC[c(phyloseq_rarefied_phages_ggplot_class_melt_UC$`Endoscopic_outcome_combined (NR/R)` == "non-remission" & phyloseq_rarefied_phages_ggplot_class_melt_UC$Diagnosis == "CD" & phyloseq_rarefied_phages_ggplot_class_melt_UC$`Shannon diversity` == "Malgrandaviricetes"),]
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired <- total_virome_R_NR_total_phages_w0_UC_caudo_non_crass
#View(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired)

colnames_de <- colnames(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired)
de1 <-data.frame("10777",NA, NA,"primary endpoint",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("12011",NA, NA,"primary endpoint",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("12405",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("5479",NA, NA,"primary endpoint",NA,NA,NA,NA)

#View(de1)
names(de1)<-c(colnames_de)
unpaired_paired_test <- rbind(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired, de1)
unpaired_paired_test$patient_ID[c(unpaired_paired_test$patient_ID == "6716" & unpaired_paired_test$Therapy_2 == "TNF")] <- "1"
unpaired_paired_test_baseline <- unpaired_paired_test[unpaired_paired_test$Week == "baseline",]
unpaired_paired_test_baseline <- unpaired_paired_test_baseline[order(unpaired_paired_test_baseline$patient_ID) , ]
unpaired_paired_test_baseline$value <- as.numeric(unpaired_paired_test_baseline$value)
unpaired_paired_test_PE <- unpaired_paired_test[!unpaired_paired_test$Week == "baseline",]
unpaired_paired_test_PE <- unpaired_paired_test_PE[order(unpaired_paired_test_PE$patient_ID) , ]
unpaired_paired_test_PE$value <- as.numeric(unpaired_paired_test_PE$value)

pm.wilcox.test(unpaired_paired_test_PE$value,unpaired_paired_test_baseline$value,alternative = "two.sided", correct = F, method = "MW-MW")
p.adjust(0.9117, method = "BH", n = 8)

unpaired_paired_test_PE_baseline <- rbind.fill(unpaired_paired_test_PE,unpaired_paired_test_baseline)
unpaired_paired_test_PE_baseline <- unpaired_paired_test_PE_baseline[order(unpaired_paired_test_PE_baseline$patient_ID),]
effect_size <- wilcox_effsize(unpaired_paired_test_PE_baseline,value~Week,paired = TRUE, alternative = "two.sided")
(effect_size$effsize)
effect_size$n1+effect_size$n2

# P-value = 0.9117 (x: baseline, y: PE)
# Xpaired: 52
# Ypaired: 52
# Xextra: 1
# Y-extra: 3

# Phages
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- phyloseq_rarefied_phages_ggplot_class_melt_UC[c(phyloseq_rarefied_phages_ggplot_class_melt_UC$`Endoscopic_outcome_combined (NR/R)` == "non-remission" & phyloseq_rarefied_phages_ggplot_class_melt_UC$Diagnosis == "CD" & phyloseq_rarefied_phages_ggplot_class_melt_UC$`Shannon diversity` == "phages"),]
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired <- total_virome_R_NR_total_phages_w0_UC_caudo_non_crass
#View(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired)

colnames_de <- colnames(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired)
de1 <-data.frame("10777",NA, NA,"primary endpoint",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("12011",NA, NA,"primary endpoint",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("12405",NA, NA,"baseline",NA,NA,NA,NA)
de1[nrow(de1) + 1,] = c("5479",NA, NA,"primary endpoint",NA,NA,NA,NA)

#View(de1)
names(de1)<-c(colnames_de)
unpaired_paired_test <- rbind(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass_paired, de1)
unpaired_paired_test$patient_ID[c(unpaired_paired_test$patient_ID == "6716" & unpaired_paired_test$Therapy_2 == "TNF")] <- "1"
unpaired_paired_test_baseline <- unpaired_paired_test[unpaired_paired_test$Week == "baseline",]
unpaired_paired_test_baseline <- unpaired_paired_test_baseline[order(unpaired_paired_test_baseline$patient_ID) , ]
unpaired_paired_test_baseline$value <- as.numeric(unpaired_paired_test_baseline$value)
unpaired_paired_test_PE <- unpaired_paired_test[!unpaired_paired_test$Week == "baseline",]
unpaired_paired_test_PE <- unpaired_paired_test_PE[order(unpaired_paired_test_PE$patient_ID) , ]
unpaired_paired_test_PE$value <- as.numeric(unpaired_paired_test_PE$value)

pm.wilcox.test(unpaired_paired_test_PE$value,unpaired_paired_test_baseline$value,alternative = "two.sided", correct = F, method = "MW-MW")
p.adjust(0.7516, method = "BH", n = 8)

unpaired_paired_test_PE_baseline <- rbind.fill(unpaired_paired_test_PE,unpaired_paired_test_baseline)
unpaired_paired_test_PE_baseline <- unpaired_paired_test_PE_baseline[order(unpaired_paired_test_PE_baseline$patient_ID),]
effect_size <- wilcox_effsize(unpaired_paired_test_PE_baseline,value~Week,paired = TRUE, alternative = "two.sided")
(effect_size$effsize)
effect_size$n1+effect_size$n2

# P-value = 0.7516 (x: baseline, y: PE)
# Xpaired: 52
# Ypaired: 52
# Xextra: 1
# Y-extra: 3
####################################
# 4.6.4 Shannon diversity boxplots between disease (UC/CD) at primary endpoint to compare inflammed and non-inflammed gut
####################################
# 4.6.4.1 Shannon diversity within Primary endpoint between disease
####################################
phyloseq_rarefied_phages_ggplot_class_melt_CD <- Mastertable_unrarefied_observed_shannon_DF_melt
phyloseq_rarefied_phages_ggplot_class_melt_CD$Week[phyloseq_rarefied_phages_ggplot_class_melt_CD$Week == "w0"] <- "baseline"
phyloseq_rarefied_phages_ggplot_class_melt_CD$Week[c(phyloseq_rarefied_phages_ggplot_class_melt_CD$Week == "w24" & phyloseq_rarefied_phages_ggplot_class_melt_CD$Diagnosis == "CD") | c(phyloseq_rarefied_phages_ggplot_class_melt_CD$Week == "w14" & phyloseq_rarefied_phages_ggplot_class_melt_CD$Diagnosis == "UC")] <- "primary endpoint"
phyloseq_rarefied_phages_ggplot_class_melt_CD <- phyloseq_rarefied_phages_ggplot_class_melt_CD[!phyloseq_rarefied_phages_ggplot_class_melt_CD$Week == "w14",]
phyloseq_rarefied_phages_ggplot_class_melt_CD <- phyloseq_rarefied_phages_ggplot_class_melt_CD[!phyloseq_rarefied_phages_ggplot_class_melt_CD$`Endoscopic_outcome_combined (NR/R)` == "unknown",]
phyloseq_rarefied_phages_ggplot_class_melt_CD <- phyloseq_rarefied_phages_ggplot_class_melt_CD[!phyloseq_rarefied_phages_ggplot_class_melt_CD$Week == "baseline",]
phyloseq_rarefied_phages_ggplot_class_melt_CD$`Shannon diversity` <- factor(phyloseq_rarefied_phages_ggplot_class_melt_CD$`Shannon diversity`, levels=c("phages", "Caudoviricetes (CrAss)","Caudoviricetes (non-CrAss)", "Malgrandaviricetes"))
phyloseq_rarefied_phages_ggplot_class_melt_CD$Diagnosis <- factor(phyloseq_rarefied_phages_ggplot_class_melt_CD$Diagnosis, levels=c("UC","CD"))
table(phyloseq_rarefied_phages_ggplot_class_melt_CD$`Shannon diversity`)
#View(phyloseq_rarefied_phages_ggplot_class_melt_CD)

phyloseq_rarefied_phages_ggplot_class_melt_CD %>%
  filter(value >= 0L & value <= 80L) %>%
  ggplot() +
  aes(x = `Endoscopic_outcome_combined (NR/R)`, 
      y = value, fill = `Shannon diversity`) +
  geom_boxplot(shape = "circle") +
  theme_bw() +
  facet_wrap(vars(Diagnosis)) +
  xlab("") +
  ylab("Shannon diversity") +
  scale_fill_manual(values=c("#669bbc","#59b031", "#ee766f", "#c49b05"),name="Class") 
####################################
# 4.6.4.2 between remission between CD -UC
####################################
phyloseq_rarefied_phages_ggplot_class_melt_CD

# Statistics within remission
# Caudo (non-crAss)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- phyloseq_rarefied_phages_ggplot_class_melt_CD[c(phyloseq_rarefied_phages_ggplot_class_melt_CD$`Endoscopic_outcome_combined (NR/R)` == "non-remission" & phyloseq_rarefied_phages_ggplot_class_melt_CD$`Shannon diversity` == "Caudoviricetes (non-CrAss)"),]
table(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$`Endoscopic_outcome_combined (NR/R)`) ## 33 R patients/ 27 NR patients
wilcox.test(value~Diagnosis, data = total_virome_R_NR_total_phages_w0_UC_caudo_non_crass, alternative = "two.sided",  paired = FALSE, exact = F)
p.adjust(0.280, method = "BH", n = 3)
effect_size <- wilcox_effsize(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass,value~Diagnosis, alternative = "two.sided",  paired = F, exact = F)
(effect_size$effsize)
effect_size$n1+effect_size$n2
# P-value = 0.04704*****
# R-squared: 0.05082342
#  n1: 32
#  n2: 46

# Caudo (crAss)
table(phyloseq_rarefied_phages_ggplot_class_melt_CD$`Observed richness`)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- phyloseq_rarefied_phages_ggplot_class_melt_CD[c(phyloseq_rarefied_phages_ggplot_class_melt_CD$`Endoscopic_outcome_combined (NR/R)` == "non-remission" & phyloseq_rarefied_phages_ggplot_class_melt_CD$`Shannon diversity` == "Caudoviricetes (CrAss)"),]
table(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$`Endoscopic_outcome_combined (NR/R)`) ## 33 R patients/ 27 NR patients
wilcox.test(value~Diagnosis, data = total_virome_R_NR_total_phages_w0_UC_caudo_non_crass, alternative = "two.sided",  paired = FALSE, exact = F)
p.adjust(0.6987, method = "BH", n = 3)
effect_size <- wilcox_effsize(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass,value~Diagnosis, alternative = "two.sided",  paired = F, exact = F)
(effect_size$effsize)
effect_size$n1+effect_size$n2
# P-value = 0.8533

# Caudo (Micro)
table(phyloseq_rarefied_phages_ggplot_class_melt_CD$`Observed richness`)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- phyloseq_rarefied_phages_ggplot_class_melt_CD[c(phyloseq_rarefied_phages_ggplot_class_melt_CD$`Endoscopic_outcome_combined (NR/R)` == "non-remission" & phyloseq_rarefied_phages_ggplot_class_melt_CD$`Shannon diversity` == "Malgrandaviricetes"),]
wilcox.test(value~Diagnosis, data = total_virome_R_NR_total_phages_w0_UC_caudo_non_crass, alternative = "two.sided",  paired = FALSE, exact = F)
p.adjust(0.00813, method = "BH", n = 6)
effect_size <- wilcox_effsize(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass,value~Diagnosis, alternative = "two.sided",  paired = F, exact = F)
(effect_size$effsize)
effect_size$n1+effect_size$n2
# P-value = 0.002334*****
# R-squared: 0.1191944
#  n1: 32
#  n2: 46

# Total phages
table(phyloseq_rarefied_phages_ggplot_class_melt_CD$`Observed richness`)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- phyloseq_rarefied_phages_ggplot_class_melt_CD[c(phyloseq_rarefied_phages_ggplot_class_melt_CD$`Endoscopic_outcome_combined (NR/R)` == "remission" & phyloseq_rarefied_phages_ggplot_class_melt_CD$`Shannon diversity` == "phages"),]
table(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$`Endoscopic_outcome_combined (NR/R)`) ## 33 R patients/ 27 NR patients
wilcox.test(value~Diagnosis, data = total_virome_R_NR_total_phages_w0_UC_caudo_non_crass, alternative = "two.sided",  paired = FALSE, exact = F)
p.adjust(0.8028, method = "BH", n = 8)
effect_size <- wilcox_effsize(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass,value~Diagnosis, alternative = "two.sided",  paired = F, exact = F)
(effect_size$effsize)
effect_size$n1+effect_size$n2

# P-value = 0.07799
# R-squared: 0.0400544
#  n1: 32
#  n2: 46
####################################
# 4.6.4.3 between non-remission between CD -UC
####################################
phyloseq_rarefied_phages_ggplot_class_melt_CD

# Statistics within non-remissioon
# Caudo (non-crAss)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- phyloseq_rarefied_phages_ggplot_class_melt_CD[c(phyloseq_rarefied_phages_ggplot_class_melt_CD$`Endoscopic_outcome_combined (NR/R)` == "non-remission" & phyloseq_rarefied_phages_ggplot_class_melt_CD$`Shannon diversity` == "Caudoviricetes (non-CrAss)"),]
table(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$`Endoscopic_outcome_combined (NR/R)`) ## 33 R patients/ 27 NR patients
wilcox.test(value~Diagnosis, data = total_virome_R_NR_total_phages_w0_UC_caudo_non_crass, alternative = "two.sided",  paired = FALSE, exact = F)
p.adjust(0.9879, method = "BH", n = 8)
effect_size <- wilcox_effsize(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass,value~Diagnosis, alternative = "two.sided",  paired = F, exact = F)
(effect_size$effsize)
effect_size$n1+effect_size$n2
# P-value = 0.2802

# Caudo (crAss)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- phyloseq_rarefied_phages_ggplot_class_melt_CD[c(phyloseq_rarefied_phages_ggplot_class_melt_CD$`Endoscopic_outcome_combined (NR/R)` == "non-remission" & phyloseq_rarefied_phages_ggplot_class_melt_CD$`Shannon diversity` == "Caudoviricetes (CrAss)"),]
table(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$`Endoscopic_outcome_combined (NR/R)`) ## 33 R patients/ 27 NR patients
wilcox.test(value~Diagnosis, data = total_virome_R_NR_total_phages_w0_UC_caudo_non_crass, alternative = "two.sided",  paired = FALSE, exact = F)
p.adjust(0.8951, method = "BH", n = 8)
effect_size <- wilcox_effsize(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass,value~Diagnosis, alternative = "two.sided",  paired = F, exact = F)
(effect_size$effsize)
effect_size$n1+effect_size$n2
# P-value = 0.6987

# Caudo (Micro)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- phyloseq_rarefied_phages_ggplot_class_melt_CD[c(phyloseq_rarefied_phages_ggplot_class_melt_CD$`Endoscopic_outcome_combined (NR/R)` == "non-remission" & phyloseq_rarefied_phages_ggplot_class_melt_CD$`Shannon diversity` == "Malgrandaviricetes"),]
table(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$`Endoscopic_outcome_combined (NR/R)`) ## 33 R patients/ 27 NR patients
wilcox.test(value~Diagnosis, data = total_virome_R_NR_total_phages_w0_UC_caudo_non_crass, alternative = "two.sided",  paired = FALSE, exact = F)
p.adjust(0.1624, method = "BH", n = 8)
effect_size <- wilcox_effsize(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass,value~Diagnosis, alternative = "two.sided",  paired = F, exact = F)
(effect_size$effsize)
effect_size$n1+effect_size$n2
# P-value = 0.1584
# R-squared: 0.09133822
#  n1: 31
#  n2: 46

# Total phages
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- phyloseq_rarefied_phages_ggplot_class_melt_CD[c(phyloseq_rarefied_phages_ggplot_class_melt_CD$`Endoscopic_outcome_combined (NR/R)` == "non-remission" & phyloseq_rarefied_phages_ggplot_class_melt_CD$`Shannon diversity` == "phages"),]
table(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$`Endoscopic_outcome_combined (NR/R)`) ## 33 R patients/ 27 NR patients
wilcox.test(value~Diagnosis, data = total_virome_R_NR_total_phages_w0_UC_caudo_non_crass, alternative = "two.sided",  paired = FALSE, exact = F)
p.adjust(0.6381, method = "BH", n = 8)
effect_size <- wilcox_effsize(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass,value~Diagnosis, alternative = "two.sided",  paired = F, exact = F)
(effect_size$effsize)
effect_size$n1+effect_size$n2

# P-value = 0.6499
# R-squared: 0.05141276
#  n1: 31
#  n2: 46
####################################
# 4.7 Pielous eveness
####################################
# UC
Pilous_evennes_melt <- melt(Pilous_evennes, id.vars = c("patient_ID","Diagnosis", "Gender", "Week", "Endoscopic_outcome_combined (NR/R)", "Therapy_2"))
Pilous_evennes_melt$variable <- as.character(Pilous_evennes_melt$variable)
Pilous_evennes_melt <- Pilous_evennes_melt[Pilous_evennes_melt$Diagnosis == "UC",]
Pilous_evennes_melt$Week[Pilous_evennes_melt$Week == "w0"] <- "baseline"
Pilous_evennes_melt$Week[Pilous_evennes_melt$Week == "w14"] <- "primary endpoint"
Pilous_evennes_melt$variable[Pilous_evennes_melt$variable == "Eveness_CrAss"] <- "Caudoviricetes (CrAss)"
Pilous_evennes_melt$variable[Pilous_evennes_melt$variable == "Eveness_non_CrAss"] <- "Caudoviricetes (non-CrAss)"
Pilous_evennes_melt$variable[Pilous_evennes_melt$variable == "Eveness_Malgranda"] <- "Malgrandaviricetes"
Pilous_evennes_melt$variable[Pilous_evennes_melt$variable == "Eveness_phages"] <- "Phages"
Pilous_evennes_melt <- Pilous_evennes_melt[!Pilous_evennes_melt$variable == "Phages",]
Pilous_evennes_melt$variable <- factor(Pilous_evennes_melt$variable, levels=c("Caudoviricetes (CrAss)","Caudoviricetes (non-CrAss)", "Malgrandaviricetes"))
Pilous_evennes_melt <- Pilous_evennes_melt[Pilous_evennes_melt$Week == "primary endpoint",]
table(Pilous_evennes_melt$variable)
View(Pilous_evennes_melt)

# Figure
Pilous_evennes_melt %>%
  filter(variable %in% c("Caudoviricetes (CrAss)", "Caudoviricetes (non-CrAss)", "Malgrandaviricetes")) %>%
 filter(!`Endoscopic_outcome_combined (NR/R)`== "unknown") %>%  
  ggplot() +
  aes(x = variable, y = value, fill = `Endoscopic_outcome_combined (NR/R)`) +
  geom_boxplot() +
  scale_fill_manual(values = c(`non-remission` = "#AD4B38", remission = "#49AD86")) +
  theme_bw() +
  xlab("") +
  ylab("Pielou's evenness") +
  theme(legend.position = "none")  

# CD
View(Pilous_evennes)
Pilous_evennes_melt_ <- melt(Pilous_evennes, id.vars = c("patient_ID","Diagnosis", "Gender", "Week", "Endoscopic_outcome_combined (NR/R)", "Therapy_2"))
Pilous_evennes_melt_$variable <- as.character(Pilous_evennes_melt_$variable)
Pilous_evennes_melt_CD <- Pilous_evennes_melt_[Pilous_evennes_melt_$Diagnosis == "CD",]
Pilous_evennes_melt_CD$Week[Pilous_evennes_melt_CD$Week == "w0"] <- "baseline"
Pilous_evennes_melt_CD$Week[Pilous_evennes_melt_CD$Week == "w24"] <- "primary endpoint"
Pilous_evennes_melt_CD$variable[Pilous_evennes_melt_CD$variable == "Eveness_CrAss"] <- "Caudoviricetes (CrAss)"
Pilous_evennes_melt_CD$variable[Pilous_evennes_melt_CD$variable == "Eveness_non_CrAss"] <- "Caudoviricetes (non-CrAss)"
Pilous_evennes_melt_CD$variable[Pilous_evennes_melt_CD$variable == "Eveness_Malgranda"] <- "Malgrandaviricetes"
Pilous_evennes_melt_CD$variable[Pilous_evennes_melt_CD$variable == "Eveness_phages"] <- "Phages"
Pilous_evennes_melt_CD <- Pilous_evennes_melt_CD[!Pilous_evennes_melt_CD$variable == "Phages",]
Pilous_evennes_melt_CD$variable <- factor(Pilous_evennes_melt_CD$variable, levels=c("Caudoviricetes (CrAss)","Caudoviricetes (non-CrAss)", "Malgrandaviricetes"))
Pilous_evennes_melt_CD <- Pilous_evennes_melt_CD[Pilous_evennes_melt_CD$Week == "primary endpoint",]
table(Pilous_evennes_melt_CD$variable)
View(Pilous_evennes_melt_CD)

# Figure
Pilous_evennes_melt_CD %>%
  filter(variable %in% c("Caudoviricetes (CrAss)", "Caudoviricetes (non-CrAss)", "Malgrandaviricetes")) %>%
  filter(!`Endoscopic_outcome_combined (NR/R)`== "unknown") %>%  
  ggplot() +
  aes(x = variable, y = value, fill = `Endoscopic_outcome_combined (NR/R)`) +
  geom_boxplot() +
  scale_fill_manual(values = c(`non-remission` = "#AD4B38", remission = "#49AD86")) +
  theme_bw() +
  xlab("") +
  ylab("Pielou's evenness") +
  theme(legend.position = "none")  

# Statistics
#UC
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- Pilous_evennes_melt[c(Pilous_evennes_melt$variable == "Malgrandaviricetes"),]
table(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$`Endoscopic_outcome_combined (NR/R)`) ## 33 R patients/ 27 NR patients
wilcox.test(value~`Endoscopic_outcome_combined (NR/R)`, data = total_virome_R_NR_total_phages_w0_UC_caudo_non_crass, alternative = "two.sided",  paired = FALSE, exact = F)
p.adjust(0.549, method = "BH", n = 3)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$remission <- total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$`Endoscopic_outcome_combined (NR/R)`
effect_size <- wilcox_effsize(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass,value~remission)
(effect_size$effsize)
effect_size$n1+effect_size$n2

# CD
Pilous_evennes_melt_CD <- Pilous_evennes_melt_CD[!Pilous_evennes_melt_CD$`Endoscopic_outcome_combined (NR/R)` == "unknown",]
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass <- Pilous_evennes_melt_CD[c(Pilous_evennes_melt_CD$variable == "Caudoviricetes (non-CrAss)"),]
wilcox.test(value~`Endoscopic_outcome_combined (NR/R)`, data = total_virome_R_NR_total_phages_w0_UC_caudo_non_crass, alternative = "two.sided",  paired = FALSE, exact = F)
p.adjust(0.8411, method = "BH", n = 3)
total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$remission <- total_virome_R_NR_total_phages_w0_UC_caudo_non_crass$`Endoscopic_outcome_combined (NR/R)`
effect_size <- wilcox_effsize(total_virome_R_NR_total_phages_w0_UC_caudo_non_crass,value~remission)
(effect_size$effsize)
effect_size$n1+effect_size$n2

# between disease
table(Pilous_evennes_melt$variable)
Pilous_evennes_melt$Week[Pilous_evennes_melt$Week == "w24"] <- "primary endpoint"
Pilous_evennes_melt$Week[c(Pilous_evennes_melt$Week == "w14" & Pilous_evennes_melt$Diagnosis == "UC")] <- "primary endpoint"
Pilous_evennes_melt <- Pilous_evennes_melt[Pilous_evennes_melt$Week == "primary endpoint",]
Pilous_evennes_melt_ <- Pilous_evennes_melt[c(Pilous_evennes_melt$`Endoscopic_outcome_combined (NR/R)` == "non-remission" & Pilous_evennes_melt$variable == "Eveness_non_CrAss"),]
table(Pilous_evennes_melt_$`Endoscopic_outcome_combined (NR/R)`) ## 33 R patients/ 27 NR patients
wilcox.test(value~Diagnosis, data = Pilous_evennes_melt_, alternative = "two.sided",  paired = FALSE, exact = F)
p.adjust(0.5673, method = "BH", n = 3)
effect_size <- wilcox_effsize(Pilous_evennes_melt_,value~Diagnosis, alternative = "two.sided",  paired = F, exact = F)
(effect_size$effsize)
effect_size$n1+effect_size$n2
####################################
# 5. LefSe Analysis on primary endpoint: between UC & CD at remission
####################################
### LefSe: Differential abundance analysis ###
####################################
# 5.1. Create phyloseq Objects (needed for LefSe)
####################################
Mastertable_viral_rarefied_phages
####################################
# 5.1.1 Create phyloseq abundance table
####################################
Mastertable_viral_rarefied_phages_ <- Mastertable_viral_rarefied_phages
Mastertable_viral_rarefied_phages_$Final_ANI <- NULL
Mastertable_viral_rarefied_phages_$Final_coverage <- NULL
Mastertable_viral_rarefied_phages_$Final_class2 <- NULL
Mastertable_viral_rarefied_phages_$Totalnumberofreads <- NULL
Mastertable_viral_rarefied_phages_$Blastn_AS <- NULL
Mastertable_viral_rarefied_phages_$Diamond_AS <- NULL
Mastertable_viral_rarefied_phages_$Best_AS <- NULL
Mastertable_viral_rarefied_phages_$Final_subfamily <- NULL
Mastertable_viral_rarefied_phages_$Final_nodes <- NULL

colnames(Mastertable_viral_rarefied_phages_)[which(names(Mastertable_viral_rarefied_phages_) == "Final_superkingdom")] <- "Final_kingdom"
vector_1 <- (which(names(Mastertable_viral_rarefied_phages_)== "Virsorter")-1) 
vector_2 <- (which(names(Mastertable_viral_rarefied_phages_)== "Final_kingdom"))
vector_3 <- which(names(Mastertable_viral_rarefied_phages_)== "Final_species")

abundance_table_rarefied <- Mastertable_viral_rarefied_phages_[,c(1:vector_1)]
####################################
# 5.1.2 Create phyloseq taxonomy table
####################################
taxonomy_table_rarefied <- Mastertable_viral_rarefied_phages_[,c(vector_2:vector_3)]
names(taxonomy_table_rarefied)[names(taxonomy_table_rarefied) == "Final_kingdom"] <- "Kingdom"
names(taxonomy_table_rarefied)[names(taxonomy_table_rarefied) == "Final_phylum"] <- "Phylum"
names(taxonomy_table_rarefied)[names(taxonomy_table_rarefied) == "Final_class"] <- "Class"
names(taxonomy_table_rarefied)[names(taxonomy_table_rarefied) == "Final_order"] <- "Order"
names(taxonomy_table_rarefied)[names(taxonomy_table_rarefied) == "Final_family"] <- "Family"
#names(taxonomy_table_rarefied)[names(taxonomy_table_rarefied) == "Final_subfamily"] <- "Subfamily"
names(taxonomy_table_rarefied)[names(taxonomy_table_rarefied) == "Final_genus"] <- "Genus"
names(taxonomy_table_rarefied)[names(taxonomy_table_rarefied) == "Final_species"] <- "Species"
####################################
# 5.1.3 Create phyloseq sample tables (contain samples + metadata)
####################################
setwd("/Users/daan/Desktop/Bioinformatics/pre-processing/Biological_study/Metadata/output_metadata")
getwd()
dir()

metadata <- read_excel("Final_Metadata_subselection_Relevant_R_without-units.xlsx")
rownames(metadata) <- metadata$`Randomized Nr.`
names(metadata)[1]<-paste("sample")

sample_table_rarefied_for_all <- (merge(t(abundance_table_rarefied),metadata, by=0, all=T)) # this is done so the metadata matches the samples
rownames(sample_table_rarefied_for_all) <- sample_table_rarefied_for_all$Row.names
sample_table_rarefied_for_all$Row.names <- NULL

vector_4 <- (which(names(sample_table_rarefied_for_all)== "sample"))
vector_5 <- which(names(sample_table_rarefied_for_all)== "Biomarker_outcome_combined")
sample_table_rarefied_for_all <- sample_table_rarefied_for_all[,c(vector_4:vector_5)]
sample_table_rarefied_for_all$Endoscopic_outcome_combined[sample_table_rarefied_for_all$Endoscopic_outcome_combined == "NA"] <- "unknown"

# Add new columns 'timepoints'
sample_table_rarefied_for_all$timepoints[sample_table_rarefied_for_all$Week == "0"] <- "baseline"
sample_table_rarefied_for_all$timepoints[!sample_table_rarefied_for_all$Week == "0"] <- "primary endpoint"

## Remove Data (unknowns & w14)
sample_table_rarefied_for_all <- sample_table_rarefied_for_all[!c(sample_table_rarefied_for_all$Week == "14" & sample_table_rarefied_for_all$Diagnosis == "CD"),]
sample_table_rarefied_for_all <- sample_table_rarefied_for_all[!sample_table_rarefied_for_all$Endoscopic_outcome_combined == "unknown",]
sample_table_rarefied_for_all <- sample_table_rarefied_for_all[c(sample_table_rarefied_for_all$timepoints == "primary endpoint" & sample_table_rarefied_for_all$Endoscopic_outcome_combined == "1"),]
#View(sample_table_rarefied_for_all)

sample_table_rarefied_for_all$Endoscopic_outcome_combined[sample_table_rarefied_for_all$Endoscopic_outcome_combined == "1"] <- "remission"
sample_table_rarefied_for_all$Endoscopic_outcome_combined[sample_table_rarefied_for_all$Endoscopic_outcome_combined == "0"] <- "non-remission"
####################################
# 5.1.4 Overview all unrarefied tables (abundance, taxonomy and sample tables)
####################################
## ABUNDANCE TABLE
#View(abundance_table_rarefied)

## TAXONOMY TABLE 
#View(taxonomy_table_rarefied)

## SAMPLE TABLE
#View(sample_table_rarefied_for_all)
####################################
# 5.2 Create phyloSeq objects
####################################
# 5.2.1 Create a matrix of abundance and taxonomy tables
####################################
abundance_table_rarefied_m <- as.matrix(abundance_table_rarefied)
taxonomy_table_rarefied_m <- as.matrix(taxonomy_table_rarefied)
####################################
# 5.2.2 Transform to phylseq objects
####################################
ABUNDANCE_rarefied_family <- otu_table(abundance_table_rarefied_m, taxa_are_rows = TRUE)
TAX_rarefied_family <- tax_table(taxonomy_table_rarefied_m)
samples <- sample_data(sample_table_rarefied_for_all)
phyloseq_rarefied_family <- phyloseq(ABUNDANCE_rarefied_family,TAX_rarefied_family, samples)
####################################
# 5.2.3 Visualize data & subset only the phage fraction
####################################
sample_names(phyloseq_rarefied_family) # All sample names
rank_names(phyloseq_rarefied_family) # All taxonomies
sample_variables(phyloseq_rarefied_family) ## All metadata
# if you want to subset the dataframe into only "responders", just do "responder_data <- subset_samples(phyloseq_unrarefied_family, response == "responder")"
phyloseq_rarefied_phages <- subset_taxa(phyloseq_rarefied_family, Kingdom == "Viruses")
####################################
# 5.2.4 PhyloSeq format
####################################
phyloseq_rarefied_phages
####################################
# 5.3 Actual LefSe analysis
####################################
# 5.3.1 LDA plot: enriched group
####################################
#?run_lefse()

#install.packages("xlsx")
#library(xlsx)

mm_lefse <- run_lefse(
  phyloseq_rarefied_phages,
  class = "Diagnosis",
  transform = c("identity"),
  taxa_rank = "none",
  kw_cutoff = 0.05,
  wilcoxon_cutoff = 0.05,
  norm = "none",
  multicls_strat = FALSE,
  lda_cutoff = 3,
  correct = "1")

# lefse return a microbioMarker class inherits from phyloseq
mm_lefse

# Evaluate LEfSe
setwd("/Users/daan/Desktop/Bioinformatics/Analysis/Biologicals/Phages/Complex_heatmap_clustering/primary_endpoint/Lefse/Diagnosis")
View(marker_table(mm_lefse))
mm_lefse_table <- marker_table(mm_lefse)

Mastertable_viral_rarefied_phages$names <- rownames(Mastertable_viral_rarefied_phages)
View(Mastertable_viral_rarefied_phages[Mastertable_viral_rarefied_phages$names == "NODE_4_length_96496_cov_192_866458_B354",1:377])

UC_colnames <- c("B107","B110","B112","B124","B130","B176","B184","B205","B21","B211","B229","B256","B264","B297","B305","B310","B340","B347","B36","B364","B365","B383","B400","B42","B430","B433","B445","B45","B57","B6","B71")
sample_names_t

sample_names_t <- sample_names(phyloseq_rarefied_family) # All sample names
table(Mastertable_viral_rarefied_phages[Mastertable_viral_rarefied_phages$names == "NODE_4_length_57512_cov_289_045477_B331",UC_colnames] > 1)
table(Mastertable_viral_rarefied_phages[Mastertable_viral_rarefied_phages$names == "NODE_4_length_57512_cov_289_045477_B331",sample_names_t] > 1)

# total= 84 IBD (53 CD + 31 UC)

# ENRICHED IN UC
# NODE_4_length_96496_cov_192_866458_B354
(8/31)*100 # 25,81% in UC 
(4/53)*100 # 7,55% in CD

# NODE_12_length_26561_cov_46_724475_B333
(10/31)*100 # 32,26% in UC 
(5/53)*100 # 9,44% in CD

# NODE_2_length_84718_cov_385_032029_B386
(6/31)*100 # 19,35% in UC 
(3/53)*100 # 5,66% in CD

# NODE_369_length_3367_cov_219_504559_B102
(3/31)*100 # 9,68% in UC 
(0/53)*100 # 0% in CD

# NODE_4_length_57512_cov_289_045477_B331
(5/31)*100 # 16,13% in UC 
(2/53)*100 # 3,77% in CD

write.xlsx(mm_lefse_table, file = "LefSe_primary_endpoint_remission_enriched.xlsx",
           sheetName = "Lefse", append = FALSE)
####################################
# 5.3.2 Visualization
####################################
plot_ef_bar(mm_lefse, label_level = 0) +
  scale_fill_manual(values = c("UC" = "#d9846b", "CD" = "#8080b4"))

# CrAss phage: 6 in UC & 3 in CD. 
####################################