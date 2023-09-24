####################################
# SCRIPT 6: VIRAL TRANSPLANT
####################################
# 0. Packages: Install packages and load them wherever needed
####################################
#install.packages("ggplot2")
library(ggplot2)
#install.packages("dplyr")
library(dplyr)
#install.packages("plyr")
library(plyr)
#install.packages("reshape")
library(reshape)
#install.packages("reshape2")
library(reshape2)
#install.packages("scales")
library(scales)
#install.packages("viridis")
library(viridis)
#install.packages("readr")
library(readr)
#install.packages("WriteXLS")
library(WriteXLS)
#install.packages("readxl")
library(readxl)
#install.packages("remotes")
library(remotes)
#install.packages("devtools") # I have the problem that Mac Big sure blocks access to all kind of system on the mac and had to disable system integratie protection
library(devtools)
#install.packages("BiocManager") # # before downlading from github you need to load devtools
library(BiocManager)
#BiocManager::install("ComplexHeatmap")
#install_github("jokergoo/ComplexHeatmap")
library(ComplexHeatmap)
#source("https://raw.githubusercontent.com/joey711/phyloseq/master/inst/scripts/installer.R",local = TRUE) ## like this it works!!
library(phyloseq)
#BiocManager::install("ALDEx2")
library(ALDEx2)
#install.packages("RSQLite")
library(RSQLite)
#install.packages("https://cran.r-project.org/src/contrib/Archive/RcppArmadillo/RcppArmadillo_0.9.900.3.0.tar.gz", repos=NULL, type="source")
library(RcppArmadillo)
#BiocManager::install("DESeq2")
library(DESeq2)
#remotes::install_github("yiluheihei/microbiomeMarker")
#install.packages("microbiomeMarker")
library(microbiomeMarker)
#install.packages("vegan")
library(vegan)
#install.packages("dendsort")
library(dendsort)
#install.packages("ape")
library(ape)
#install.packages("seriation")
library(seriation)
#install.packages("taxonomizr")
library(taxonomizr)
#install.packages("esquisse")
library(esquisse)
#install.packages("tidyverse")
library(tidyverse)
#install.packages("grid")
library(grid)
#install.packages("coin")
library(coin)
#install.packages("rstatix")
library(rstatix)
#install.packages("PMCMR")
library(PMCMR)
#remotes::install_github("tapj/biotyper")
library(biotyper)
#BiocManager::install("DirichletMultinomial")
#install.packages("DirichletMultinomial")
library(DirichletMultinomial)
#install.packages("shiny")
library(shiny)
#install.packages("fpc")
library (fpc)
#install.packages("clValid")
library (clValid)
#install.packages("cluster")
library(cluster) 
#install.packages("philentropy")
library(philentropy)
#source("https://raw.githubusercontent.com/joey711/phyloseq/master/inst/scripts/installer.R",local = TRUE) ## like this it works!!
library(phyloseq)
#BiocManager::install("microbiome")
library(microbiome)
####################################
# 1. Datasets
####################################
Mastertable # Reads & contigs figure
Mastertable_genus1 # Percentage classified on class, order, family level
Mastertable_eukaryotic # Eukaryotic viral figure
Mastertable_phage  # Phage analysis
#################################### 
# 2. Phyloseq object
####################################
# 2.1 Upload Metadata
####################################
setwd("/Users/daan/Desktop/Bioinformatics/pre-processing/FMT_study/Metadata/viral_transplant")
dir()

metadata <- read_excel("samples_FMT_project_Final_DJ_R.xlsx")
metadata <- as.data.frame(metadata)
metadata$Bacterial_cellcount <- as.numeric(metadata$Bacterial_cellcount)
metadata$CRP <- as.numeric(metadata$CRP)
metadata$Hemoglobin <- as.numeric(metadata$Hemoglobin)
metadata$Age <- as.numeric(metadata$Age)
metadata$BMI <- as.numeric(metadata$BMI)
metadata$Disease_duration <- as.numeric(metadata$Disease_duration)
metadata$Smoking <- as.character(metadata$Smoking)
metadata$Calprotectin <- as.numeric(metadata$Calprotectin)
metadata$Moisture <- as.numeric(metadata$Moisture)
rownames(metadata) <- metadata$Sample_ID
colnames(metadata)
metadata[metadata=="NA"] <- NA
str(metadata)

samples <- metadata

nrow(Mastertable_phage)
table(Mastertable_phage$Phyla_phylum[Mastertable_phage$Genus=="Genus5"])
table(Mastertable_phage$Phyla_family[Mastertable_phage$Genus=="Genus5"])
table(Mastertable_phage$Phyla_genus[Mastertable_phage$Genus=="Genus5"])
(2/(39+9+11+6))*100
(7/(14+9))*100
####################################
# 2.2 Abundance table
####################################
Mastertable_viral_rarefied_2 <- Mastertable_phage

# Select phages
table(Mastertable_viral_rarefied_2$Final_viral)
sort(colSums(Mastertable_viral_rarefied_2[1:304]))
names <- colnames(Mastertable_viral_rarefied_2[1:304][,colSums(Mastertable_viral_rarefied_2[1:304])!=0])
names # 302 samples with phages lefr

# Select samples with phages present and create abundance table
abundance_table <- Mastertable_viral_rarefied_2
abundance_table$names <- rownames(abundance_table)
abundance_tabe1 <- aggregate(. ~names, FUN = sum, data = abundance_table[,colnames(abundance_table) %in% names | colnames(abundance_table) == 'names'])
rownames(abundance_tabe1) <- abundance_tabe1$names
abundance_tabe1$names <- NULL
abundance_table_rarefied_m_samples <- as.matrix(abundance_tabe1)

table(Mastertable_viral_rarefied_2$Final_order[Mastertable_viral_rarefied_2$Genus=='Genus13'])
table(Mastertable_viral_rarefied_2$Phyla_phylum[Mastertable_viral_rarefied_2$Genus=='Genus13'])
table(Mastertable_viral_rarefied_2$Phyla_genus[Mastertable_viral_rarefied_2$Genus=='Genus13'])
table(Mastertable_viral_rarefied_2$`lysogenic cycle`[Mastertable_viral_rarefied_2$Genus=='Genus13'])
####################################
# 2.3 Taxonomy table
####################################
# Create taxonomy table
colnames(Mastertable_viral_rarefied_2)
Mastertable_viral_rarefied_2$Final_ANI <- NULL
Mastertable_viral_rarefied_2$Final_coverage <- NULL
Mastertable_viral_rarefied_2$Totalnumberofreads <- NULL
Mastertable_viral_rarefied_2$Blastn_AS <- NULL
Mastertable_viral_rarefied_2$Diamond_AS <- NULL
Mastertable_viral_rarefied_2$Best_AS <- NULL
Mastertable_viral_rarefied_2$Final_genus <- Mastertable_viral_rarefied_2$Genus
vector_1 <- (which(names(Mastertable_viral_rarefied_2)== "Final_superkingdom"))
vector_2 <- which(names(Mastertable_viral_rarefied_2)== "Final_species")
Taxonomy_table <- Mastertable_viral_rarefied_2[,c(vector_1:vector_2)]

# Extract phyla-infect phages
#phyla_column <- Mastertable_viral_rarefied_2$Final_order
#rownames_column <- rownames(Mastertable_viral_rarefied_2)
#phyla_column[phyla_column == 0] <- "Unknown"
#modified_data <- data.frame(RowNames = rownames_column, Phyla_phylum = phyla_column)
#rownames(modified_data) <- modified_data$RowNames
#modified_data$RowNames <- NULL
#View(modified_data)

# Rename some Taxonomy elements
colnames(Taxonomy_table)
names(Taxonomy_table)[names(Taxonomy_table) == "Final_superkingdom"] <- "Kingdom"
names(Taxonomy_table)[names(Taxonomy_table) == "Final_phylum"] <- "Phylum"
names(Taxonomy_table)[names(Taxonomy_table) == "Final_class"] <- "Class"
names(Taxonomy_table)[names(Taxonomy_table) == "Final_order"] <- "Order"
names(Taxonomy_table)[names(Taxonomy_table) == "Final_family"] <- "Family"
names(Taxonomy_table)[names(Taxonomy_table) == "Final_subfamily"] <- "Subfamily"
names(Taxonomy_table)[names(Taxonomy_table) == "Final_genus"] <- "Genus"
names(Taxonomy_table)[names(Taxonomy_table) == "Final_species"] <- "Species"
Taxonomy_table$Subfamily <- NULL

# Add table
#Taxonomy_table <- merge(Taxonomy_table,modified_data,by=0,all=F)
#rownames(Taxonomy_table) <- Taxonomy_table$Row.names
#Taxonomy_table$Row.names <- NULL
#Taxonomy_table$Final_order <- NULL
#names(Taxonomy_table)[names(Taxonomy_table) == "Phyla_phylum"] <- "Order"

#Taxonomy_table <- Taxonomy_table[,c(7,1, 2, 3,5,6,4)]
Taxonomy_table <- Taxonomy_table[,c(7,1, 2, 3,4,5,6)]
taxonomy_table_rarefied_m <- as.matrix(Taxonomy_table)
#View(taxonomy_table_rarefied_m)

# Species level
#Taxonomy_table$Genus <- rownames(Taxonomy_table)
#taxonomy_table_rarefied_m <- as.matrix(Taxonomy_table)
#View(taxonomy_table_rarefied_m)
####################################
# 2.4 Insert Host Iphop
####################################
# 2.4.1 Enter original file
####################################
setwd ("/Users/daan/Desktop/Transfer/Host") 
getwd() 
dir()

# Take only those rows with column 1 duplicate name & highest column 5 confidence
AllSamples_completeness_Host_phyla <- as.data.frame(read_delim("Host1.txt", "\t", escape_double = FALSE, trim_ws = TRUE, col_names = FALSE))
modified_data <- AllSamples_completeness_Host_phyla %>%
  group_by(X1) %>%
  slice_max(order_by = X5) %>%
  ungroup()

duplicated_column1 <- duplicated(modified_data$X1)
duplicated_strings <- modified_data$X1[duplicated_column1]
print(duplicated_strings)

write.table(modified_data, file = "./modified_data.txt", sep = "\t",row.names = FALSE)
####################################
# 2.4.2 Enter new file
####################################
# Correct by same name with same confidence
AllSamples_completeness_Host_phyla1 <- as.data.frame(read_delim("modified_data1.txt", "\t", escape_double = FALSE, trim_ws = TRUE, col_names = TRUE))
AllSamples_completeness_Host_phyla1$X5 <- NULL
colnames(AllSamples_completeness_Host_phyla1) <- c('contig', 'Host_phyla1', "Host_family1","Host_genus1")
rownames(AllSamples_completeness_Host_phyla1) <- AllSamples_completeness_Host_phyla1$contig
AllSamples_completeness_Host_phyla1$contig <- NULL
rownames(AllSamples_completeness_Host_phyla1) <- gsub("\\.{1}", "_", rownames(AllSamples_completeness_Host_phyla1))

Mastertable_viral_rarefied_2_H <- merge(Mastertable_viral_rarefied_2,AllSamples_completeness_Host_phyla1,by=0,all=T)
rownames(Mastertable_viral_rarefied_2_H) <- Mastertable_viral_rarefied_2_H$Row.names
Mastertable_viral_rarefied_2_H$Row.names <- NULL

Mastertable_viral_rarefied_2_H$Host_phyla1[is.na(Mastertable_viral_rarefied_2_H$Host_phyla1)] <- "Unknown"
Mastertable_viral_rarefied_2_H$Host_family1[is.na(Mastertable_viral_rarefied_2_H$Host_family1)] <- "Unknown"
Mastertable_viral_rarefied_2_H$Host_genus1[is.na(Mastertable_viral_rarefied_2_H$Host_genus1)] <- "Unknown"
table(Mastertable_viral_rarefied_2_H$Host_phyla1)
table(Mastertable_viral_rarefied_2_H$Host_family1)
table(Mastertable_viral_rarefied_2_H$Host_genus1)
####################################  
# 2.5 Create phyloseq object
####################################
ABUNDANCE_rarefied_family <- otu_table(abundance_table_rarefied_m_samples, taxa_are_rows = T)
TAX_rarefied_family <- tax_table(taxonomy_table_rarefied_m)
samples <- sample_data(metadata)
phyloseq_rarefied_samples <- phyloseq(ABUNDANCE_rarefied_family,TAX_rarefied_family, samples)

# Remove crappy sample with less than 100 k reads
# phyloseq_rarefied_samples <- rarefy_even_depth(phyloseq_rarefied_samples, sample.size = 500000)
#sample_reads <- sample_sums(phyloseq_rarefied_samples)
#samples_to_remove <- sample_reads < 100000
#phyloseq_filtered <- subset_samples(phyloseq_rarefied_samples, !samples_to_remove)
#phyloseq_rarefied_samples <- phyloseq_filtered
#phyloseq_rarefied_samples
library(microbiome)

phyloseq_rarefied_phages <- phyloseq_rarefied_samples %>%
  transform(transform = "identity") %>%
  aggregate_rare(level = "Genus", detection = 0.01/100, prevalence = 1/100)

#View(tax_table(phyloseq_rarefied_phages))

#phyloseq_phages_DMM <- phyloseq_rarefied_samples1 %>%
#  microbiome::transform(transform = "identity") %>%
#  aggregate_rare(level = "Genus", detection = 0.00001/100, prevalence = 1/100)
####################################
# 3. Keep only donor-patient pairs
####################################
# Autologous donor + patient pairs
list_autologous_keep <- c("AP83", "AP78","AP70","AP69","AP68","AP66","AP61","AP58","AP53","AP52","AP50","AP46","AP41","AP23","AP22","AP20","AP12","AP11","AP1")

# Superdonor + patients pairs
list_superdonor_keep <- c("99","98","90","88","83","80","76","75","74","73","70","67","65","50","40","37","36","31","27","23","2","15","11","109","1")

# Combine and get the sample names to keep from the list
combined_list_pair <- c(list_autologous_keep, list_superdonor_keep)

# Extract the samples that match the "FMT_treatment_batch" in the metadata
phyloseq_rarefied_phages_subset <- subset_samples(phyloseq_rarefied_phages, FMT_treatment_batch %in% combined_list_pair)
####################################
# 4. Extract donor samples
####################################
# Create donor/patient list
list_donor <- c("donor")
list_patient <- c("patient")

# Subset donor samples
phyloseq_rarefied_phages_subset
phyloseq_phages_donor <- subset_samples(phyloseq_rarefied_phages_subset, Disease_status %in% list_donor)

# Subset patient samples
phyloseq_phages_patient <- subset_samples(phyloseq_rarefied_phages_subset, Disease_status %in% list_patient)
####################################
# 5. Aggregate donors samples by specific donor batch 
####################################
Donor_samples <- merge_samples(phyloseq_phages_donor, "FMT_treatment_batch",fun=mean)

# Transpose the otu_table
transposed_otu_table <- t(otu_table(Donor_samples))

# Create a new phyloseq object with the transposed otu_table
phyloseq_transposed <- phyloseq(tax_table(Donor_samples),
                                sample_data(Donor_samples),
                                transposed_otu_table)

# Add FMT_treatment_batch names by rownames
sample_data(phyloseq_transposed)$FMT_treatment_batch <- rownames(sample_data(phyloseq_transposed))

# Only keep this column
phyloseq_transposed <- subset_samples(phyloseq_transposed, select = FMT_treatment_batch)

# Add column with new names
donor_numbering <- 1:length(sample_data(phyloseq_transposed)$FMT_treatment_batch)
sample_data(phyloseq_transposed)$Sample_ID <- paste0("Donor_Batch", donor_numbering)

# Add new sample names to phyloseq object
sample_names(phyloseq_transposed) <- sample_data(phyloseq_transposed)$Sample_ID
####################################
# 6. Combine patients and donors
####################################
# Add column FMT_treatment_type2
sample_data(phyloseq_transposed)$FMT_treatment_type2 <- rep('Donor (-superdonor)', nrow(sample_data(phyloseq_transposed)))

# Merge the two phyloseq objects
combined_phyloseq <- merge_phyloseq(phyloseq_transposed, phyloseq_phages_patient)

# Calculate RA matrix
phyloseq_phages_donor_patient <- combined_phyloseq %>% 
  microbiome::transform(transform = "hellinger") %>%
  aggregate_rare(level = "Genus", detection = 0.01/100, prevalence = 0.001/100)
####################################
# 7 Create a list of phyloseq table for each donor batch - patient pair
####################################
# Get unique values of 'FMT_treatment_batch'
unique_batches <- unique(sample_data(phyloseq_phages_donor_patient)$FMT_treatment_batch)

# Create an empty list to store the subsetted phyloseq objects
subsetted_phyloseqs <- list()

# Loop through each unique value of 'FMT_treatment_batch'
for (batch in unique_batches) {
  # Subset samples with the current 'FMT_treatment_batch' value
  subset_phyloseq <- subset_samples(phyloseq_phages_donor_patient, FMT_treatment_batch == batch)
  
  # Store the subsetted phyloseq object in the list
  subsetted_phyloseqs[[as.character(batch)]] <- subset_phyloseq
}
####################################
# 8. Remove Taxa with zero abundance in each list
####################################
# Loop through each phyloseq object and filter rows with zero row sums in OTU table
for (i in seq_along(subsetted_phyloseqs)) {
  otu_table <- otu_table(subsetted_phyloseqs[[i]])
  
  # Calculate row sums of the OTU table
  row_sums <- rowSums(otu_table)
  
  # Keep only rows with non-zero row sums
  filtered_otu_table <- otu_table[row_sums > 0, ]
  
  # Update the OTU table in the current phyloseq object
  subsetted_phyloseqs[[i]] <- phyloseq(sample_data(subsetted_phyloseqs[[i]]), filtered_otu_table, tax_table(subsetted_phyloseqs[[i]]))
}

subsetted_phyloseqs_trimmed <- subsetted_phyloseqs
# The physeq_list now contains the filtered phyloseq objects with rows removed from OTU tables where the row sums were zero.

#View(otu_table(subsetted_phyloseqs_trimmed[[1]]))
####################################
# 9. Obtain dataframes out of phyloseq objects
####################################
# Create an empty list to store the otu_tables
otu_tables_list <- list()

# Loop through each phyloseq object in the original list
for (i in seq_along(subsetted_phyloseqs_trimmed)) {
  # Extract the otu_table from the ith phyloseq object
  otu_table_i <- as.data.frame(otu_table(subsetted_phyloseqs_trimmed[[i]]))
  
  # Add the otu_table to the new list
  otu_tables_list[[i]] <- otu_table_i
}

# Function to remove taxa not present in the first column
#remove_non_first_column_taxa <- function(df) {
  # Identify taxa present in the first column
#  first_column_taxa <- rownames(df)[df[, 1] != 0]
  
  # Remove the identified taxa not present in the first column from the dataframe
#  df <- df[rownames(df) %in% first_column_taxa, ]
  
#  return(df)
# }

# Apply the function to each dataframe in the list
#otu_tables_list <- lapply(otu_tables_list, remove_non_first_column_taxa)
#View(otu_tables_list[[3]])
####################################
# 10. Calculate BC dissimilarity
####################################
#calculate_spearman_correlation = function(df) {
#  cor_matrix = as.data.frame(cor(df, method = 'pearson'))
#  cor_matrix = cor_matrix[1, -1, drop = FALSE]
#  return(cor_matrix)
#}

library("vegan")

get_bray_curtis_dist <- function(df) {
  bc_dist <- vegdist(t(df), method = "bray")
  bc_dist <- as.data.frame(as.matrix(bc_dist))
  return(bc_dist[1, -1, drop = FALSE])
}

# Loop through the list of dataframes and calculate correlations
bc_df = data.frame()

for (i in seq_along(otu_tables_list)) {

# Calculate the correlations
  bc_dist_result <- get_bray_curtis_dist(otu_tables_list[[i]])
  colnames(bc_dist_result) <- gsub("^P\\d+_(\\w+)", "\\1",
                                   gsub("w8.*", "w8",
                                    colnames(bc_dist_result)))
  # append
  bc_df <- bind_rows(bc_df, bc_dist_result)
}

cor_df <- bc_df
# Replace negative correlations with zero
#cor_df[cor_df < 0] <- NA
#cor_df[cor_df < 0] <- 0
####################################
# 11. Organize dataframe accordingly
####################################
cor_df1 <- cor_df[, c('w0', 'w4', 'w8', 'w12', 'm6', 'm12')]
####################################
# 12. Merge with metadata
####################################
# Upload metadata
setwd("/Users/daan/Desktop/Bioinformatics/pre-processing/FMT_study/Metadata/viral_transplant")
dir()

metadata_D <- read_excel("samples_FMT_project2.xlsx")
metadata_D <- as.data.frame(metadata_D)
rownames(metadata_D) <- metadata_D$Sample_ID
metadata_D$Sample_ID <- NULL

# merge with data
cor_df_m <- merge(cor_df,metadata_D,by=0,all=F)
rownames(cor_df_m) <- cor_df_m$Row.names
cor_df_m$Row.names <- NULL
#View(cor_df_m)
####################################
# 13. BC dissimilarity visualization
####################################
# Melt
cor_df_melt <- melt(cor_df_m, id.vars = c("Patient_ID","FMT_treatment_type","FMT_treatment_sucess"))
cor_df_melt$value <- as.numeric(cor_df_melt$value)
#View(cor_df_melt)

# Define the desired order of levels
desired_order <- c('w0', 'w4', 'w8', 'w12', 'm6', 'm12')

# Convert the 'variable' column to a factor with the desired order
cor_df_melt$variable <- factor(cor_df_melt$variable, levels = desired_order)

# Figure 3
ggplot(cor_df_melt) +
  aes(x = variable, y = value, fill = FMT_treatment_sucess) +
  geom_boxplot() +
  scale_fill_manual(values = c(`Non-remission` = "#F8766D", `Remission` = "#0DAC86")) +
  theme_bw() +
  facet_wrap(vars(FMT_treatment_type))

# Good figure
cor_df_melt %>%
  ggplot() +
  aes(x = variable, y = value, fill = FMT_treatment_type, line = "subject") +
  geom_boxplot() +  # You can adjust the alpha value for transparency
  # geom_point(size = 1.4, alpha = 0.3, position = position_jitterdodge(jitter.width = 0.0)) +
  scale_fill_manual(values = c(`Autologous FMT` = "#F8766D", `Superdonor FMT` = "#0DAC86")) +
  theme_bw() +
  facet_wrap(vars(FMT_treatment_type), scales = "free_y")

cor_df_melt_sup <- cor_df_melt[cor_df_melt$FMT_treatment_type=="Superdonor FMT",]
cor_df_melt_aut <- cor_df_melt[cor_df_melt$FMT_treatment_type=="Autologous FMT",]

# Statistics
# Convert response to factor for logistic regression
library(lmerTest)
library(lme4)

# reference is Vieautomatically chosen to be w0)
# Superdonor FMT
model_sup <- lmer(value ~ variable + (1 | Patient_ID),data = cor_df_melt_sup) # anova(model_sup)
summary(model_sup)
anova(model_sup)
coeff_summary <- coef(summary(model_sup))
p_values_sup <- coeff_summary[, "Pr(>|t|)"]
p_values_sup_adj <- p.adjust(p_values_sup, method = "BH")
coefficients_super <- coeff_summary[, "Estimate"]
combined_superdonor <- cbind(p_values_sup, p_values_sup_adj, coefficients_super)
count_point <- cor_df_melt_sup
count_point$variable[is.na(count_point$value)] <- NA
n_per_timepoint <- as.data.frame(table(count_point$variable, useNA = "no"))
rownames(n_per_timepoint) <- n_per_timepoint$Var1
n_per_timepoint$Var1 <- NULL
rownames(combined_superdonor) <- gsub("variable", "", rownames(combined_superdonor))
combined_superdonor_final <- merge(combined_superdonor,n_per_timepoint,by=0, all=F)
combined_superdonor_final

# Autologous FMT
model_aut <- lmer(value ~ variable + (1 | Patient_ID),data = cor_df_melt_aut)
summary(model_aut)
anova(model_aut)
coeff_summary <- coef(summary(model_aut))
p_values_aut <- coeff_summary[, "Pr(>|t|)"]
p_values_aut_adj <- p.adjust(p_values_aut, method = "BH")
coefficients_aut <- coeff_summary[, "Estimate"]
combined_auto <- cbind(p_values_aut, p_values_aut_adj, coefficients_aut)
count_point <- cor_df_melt_aut
count_point$variable[is.na(count_point$value)] <- NA
n_per_timepoint <- as.data.frame(table(count_point$variable, useNA = "no"))
rownames(n_per_timepoint) <- n_per_timepoint$Var1
n_per_timepoint$Var1 <- NULL
rownames(combined_auto) <- gsub("variable", "", rownames(combined_auto))
combined_auto_final <- merge(combined_auto,n_per_timepoint,by=0, all=F)

# Paper tables
# P-values, coefficients and N
combined_superdonor_final
combined_auto_final

# Statistics endoscopic outcome R/NR
cor_df_melt_sup_w8 <- cor_df_melt[cor_df_melt$FMT_treatment_type=="Superdonor FMT" & cor_df_melt$variable=='m12',]
wilcox.test(value ~ FMT_treatment_sucess, data = cor_df_melt_sup_w8, exact=F)
effect_size <- wilcox_effsize(cor_df_melt_sup_w8,value~FMT_treatment_sucess)
(effect_size$effsize)^2
effect_size$n1 
effect_size$n2 

P_adst <- c(0.1078,0.2782,0.8155,0.5101,0.594,0.817)
p.adjust(P_adst, method = "BH")
View(cor_df_melt_sup_w8)

# Wilcoxon test.
####################################
# 14. Calculate overlapping viral genera between donor and patient
####################################
# 14.1 Calculate pre/abs matrix and removal taxa not in donor
####################################
# Loop through the list to find data frames with one column

# Remove columns with one column
#indices_to_remove <- c(30,11)
#dataframes_list <- otu_tables_list[-indices_to_remove]

dataframes_list <- otu_tables_list

# Do actual function
convert_to_presence_absence <- function(df) {
  # Check if the DataFrame has more than one column
  if (ncol(df) <= 1) {
    message("The DataFrame has only one column and cannot be converted to presence/absence.")
    return(df)
  }
  # Identify taxa present in the first column
  first_column_taxa <- rownames(df)[df[, 1] != 0]
  # Remove the identified taxa not present in the first column from the dataframe
  df <- df[rownames(df) %in% first_column_taxa, ]
  # Convert abundance values to binary presence/absence
  presence_absence_df <- apply(df, 2, function(x) ifelse(x > 0, 1, 0))
  # Set row and column names
  rownames(presence_absence_df) <- rownames(df)
  colnames(presence_absence_df) <- colnames(df)
  return(presence_absence_df)
}
# Apply the function to each dataframe in the list
presence_absence_list <- lapply(dataframes_list, convert_to_presence_absence)
presence_absence_list
####################################
# 14.2 Calculate overlapping taxa
####################################
# Convert the list of tables to a list of data frames
list_of_dataframes <- lapply(presence_absence_list, as.data.frame)

get_viral_overlap <- function(df) {
  numeric_df <- df
  donor_column <- df[, 1]
  matching_cells <- sapply(numeric_df, function(col) sum(donor_column == col))
  matching_cells <- as.data.frame(t(matching_cells))
  rownames(matching_cells) <- colnames(matching_cells[1])
  matching_cells <- (matching_cells[-1]/matching_cells[,1])*100 # remove this when not calculating %
  #  matching_cells <- matching_cells[, -1, drop = FALSE] ##. TO CALCULATE absolute number of transfrred genera
  return(matching_cells)
}

# Loop through the list of dataframes and calculate correlations
overlap_df = data.frame()
for (i in seq_along(list_of_dataframes)) {
  # Calculate the correlations
  overlap_result <- get_viral_overlap(list_of_dataframes[[i]])
  colnames(overlap_result) <- gsub("^P\\d+_(\\w+)", "\\1", colnames(overlap_result))
  colnames(overlap_result) <- gsub("w8.*", "w8", colnames(overlap_result))
  # append
  overlap_df <- bind_rows(overlap_df, overlap_result)
}

overlap_DF <- overlap_df
####################################
# 14.3 Merge metadata
####################################
overlap_df_m <- merge(overlap_DF,metadata_D,by=0,all=F)
rownames(overlap_df_m) <- overlap_df_m$Row.names
overlap_df_m$Row.names <- NULL
####################################
# 14.4 Overlap visualization
####################################
# Melt
overlap_df_melt <- melt(overlap_df_m, id.vars = c("Patient_ID","FMT_treatment_type","FMT_treatment_sucess"))
overlap_df_melt$value <- as.numeric(overlap_df_melt$value)

# Define the desired order of levels
desired_order <- c('w0', 'w4', 'w8', 'w12', 'm6', 'm12')

# Convert the 'variable' column to a factor with the desired order
overlap_df_melt$variable <- factor(overlap_df_melt$variable, levels = desired_order)

# Figure 1
ggplot(overlap_df_melt) +
  aes(x = variable, y = value, fill = FMT_treatment_type) +
  geom_boxplot() +
  scale_fill_manual(values = c(`Autologous FMT` = "#F8766D", `Superdonor FMT` = "#0DAC86")) +
  theme_bw() +
  facet_wrap(vars(FMT_treatment_type), scales = "free_y")

ggplot(overlap_df_melt) +
 aes(x = variable, y = value, fill = FMT_treatment_sucess) +
 geom_boxplot() +
 scale_fill_hue(direction = 1) +
 theme_bw() +
 facet_wrap(vars(FMT_treatment_type))

# Statistics
cor_df_melt_sup <- overlap_df_melt[overlap_df_melt$FMT_treatment_type=="Superdonor FMT",]
cor_df_melt_aut <- overlap_df_melt[overlap_df_melt$FMT_treatment_type=="Autologous FMT",]
cor_df_melt_sup$variable <- factor(cor_df_melt_sup$variable)
cor_df_melt_aut$variable <- factor(cor_df_melt_aut$variable)
# reference is automatically chosen to be w0)

# Superdonor FMT
model_sup <- lmer(value ~ variable + (1 | Patient_ID),data = cor_df_melt_sup) # anova(model_sup)
coeff_summary <- coef(summary(model_sup))
p_values_sup <- coeff_summary[, "Pr(>|t|)"]
p_values_sup_adj <- p.adjust(p_values_sup, method = "BH")
coefficients_super <- coeff_summary[, "Estimate"]
combined_superdonor <- cbind(p_values_sup, p_values_sup_adj, coefficients_super)
count_point <- cor_df_melt_sup
count_point$variable[is.na(count_point$value)] <- NA
n_per_timepoint <- as.data.frame(table(count_point$variable, useNA = "no"))
rownames(n_per_timepoint) <- n_per_timepoint$Var1
n_per_timepoint$Var1 <- NULL
rownames(combined_superdonor) <- gsub("variable", "", rownames(combined_superdonor))
combined_superdonor_final <- merge(combined_superdonor,n_per_timepoint,by=0, all=F)
combined_superdonor_final

# Autologous FMT
model_aut <- lmer(value ~ variable + (1 | Patient_ID),data = cor_df_melt_aut)
summary(model_aut)
anova(model_aut)
coeff_summary <- coef(summary(model_aut))
p_values_aut <- coeff_summary[, "Pr(>|t|)"]
p_values_aut_adj <- p.adjust(p_values_aut, method = "BH")
coefficients_aut <- coeff_summary[, "Estimate"]
combined_auto <- cbind(p_values_aut, p_values_aut_adj, coefficients_aut)
count_point <- cor_df_melt_aut
count_point$variable[is.na(count_point$value)] <- NA
n_per_timepoint <- as.data.frame(table(count_point$variable, useNA = "no"))
rownames(n_per_timepoint) <- n_per_timepoint$Var1
n_per_timepoint$Var1 <- NULL
rownames(combined_auto) <- gsub("variable", "", rownames(combined_auto))
combined_auto_final <- merge(combined_auto,n_per_timepoint,by=0, all=F)

# Paper tables
# P-values, coefficients and N
combined_superdonor_final
combined_auto_final

# Statatistcs endoscopic outcome R/NR
cor_df_melt_sup_w0 <- cor_df_melt_sup[cor_df_melt_sup$FMT_treatment_type=="Superdonor FMT" & cor_df_melt_sup$variable=='m12',]
wilcox.test(value ~ FMT_treatment_sucess, data = cor_df_melt_sup_w0, exact=F)
effect_size <- wilcox_effsize(cor_df_melt_sup_w0,value~FMT_treatment_sucess)
(effect_size$effsize)^2
effect_size$n1 
effect_size$n2 
P_adst <- c(0.3181,0.1148,0.7139,0.4208,0.04283,0.1052)
p.adjust(P_adst, method = "BH")
# wilcoxon test R/NR
####################################
# 15. Compositional analyses
####################################
## I showed that the transplantation is partially sucessfully transplanted but is not associated to endoscopic outcome.
## So how do some patient respond and others not? A different donor composition!
####################################
# 15.1 Keep only donor-patient pairs
####################################
phyloseq_rarefied_samples_D <- phyloseq_rarefied_samples

# Autologous donor + patient pairs
list_autologous_keep <- c("AP83", "AP78","AP70","AP69","AP68","AP66","AP61","AP58","AP53","AP52","AP50","AP46","AP41","AP23","AP22","AP20","AP12","AP11","AP1")

# Superdonor + patients pairs
list_superdonor_keep <- c("99","98","90","88","83","80","76","75","74","73","70","67","65","50","40","37","36","31","27","23","2","15","11","109","1")

# Combine and get the sample names to keep from the list
combined_list_pair <- c(list_autologous_keep, list_superdonor_keep)

# Extract the samples that match the "FMT_treatment_batch" in the metadata
phyloseq_rarefied_phages_subset_D <- subset_samples(phyloseq_rarefied_samples_D, FMT_treatment_batch %in% combined_list_pair)
####################################
# 15.2. Subset donor
####################################
# Create donor/patient list
list_donor <- c("donor")

# Subset donor samples
phyloseq_phages_donor1 <- subset_samples(phyloseq_rarefied_phages_subset_D, Disease_status %in% list_donor)
####################################
# 15.3 Aggregate donor in donor batch
####################################
Donor_samples <- merge_samples(phyloseq_phages_donor1, "FMT_treatment_batch",fun=mean)

# Transpose the otu_table
transposed_otu_table <- t(otu_table(Donor_samples))

# Create a new phyloseq object with the transposed otu_table
phyloseq_transposed <- phyloseq(tax_table(Donor_samples),
                                sample_data(Donor_samples),
                                transposed_otu_table)

# Add FMT_treatment_batch names by rownames
sample_data(phyloseq_transposed)$FMT_treatment_batch <- rownames(sample_data(phyloseq_transposed))

# Only keep this column
phyloseq_transposed <- subset_samples(phyloseq_transposed, select = FMT_treatment_batch)

# Add column with new names
donor_numbering <- 1:length(sample_data(phyloseq_transposed)$FMT_treatment_batch)
sample_data(phyloseq_transposed)$Sample_ID <- paste0("Donor_Batch", donor_numbering)

# Add new sample names to phyloseq object
sample_names(phyloseq_transposed) <- sample_data(phyloseq_transposed)$Sample_ID
####################################
# 15.4 Add endoscopic outcome for matching patients
####################################
phyloseq_transposed_D <- phyloseq_transposed

# Get the current sample_data from the phyloseq object
sample_data <- sample_data(phyloseq_transposed_D)

# Add a new column with value '0'
sample_data$endoscopic_outcome <- "Non-remission"
remission_list <- c("Donor_Batch1", "Donor_Batch6","Donor_Batch8","Donor_Batch13","Donor_Batch19","Donor_Batch24","Donor_Batch27","Donor_Batch33","Donor_Batch37","Donor_Batch41","Donor_Batch43")
sample_data$endoscopic_outcome[sample_data$Sample_ID %in% remission_list] <- "Remission"

# Add a new column with
labels <- c(rep("Superdonor FMT", 25), rep("Autologous FMT", nrow(sample_data) - 25))
sample_data$FMT_treatment_type <- labels

# Update the sample_data in the phyloseq object
sample_data(phyloseq_transposed_D) <- sample_data
####################################
# 15.5 Calculate RA matrix
####################################
#phyloseq_transposed_D <- phyloseq_transposed # use this if not combined batch
phyloseq_transposed_D

# Extract tables
otu_table <- as.data.frame(otu_table(phyloseq_transposed_D))
tax_table <- as.data.frame(tax_table(phyloseq_transposed_D))
metadata_table <- as.data.frame(sample_data(phyloseq_transposed_D))

# Select data
#metadata_table <- metadata_table[, c("Patient_ID", "FMT_treatment_batch", "FMT_treatment_sucess", "FMT_treatment_type", 'Timepoints')]
metadata_table <- metadata_table[, c("FMT_treatment_batch", "Sample_ID", "endoscopic_outcome", "FMT_treatment_type")]

# Merge table
DF_1 <- merge(otu_table,tax_table, by = 0, all =F)
rownames(DF_1) <- DF_1$Row.names
DF_1$Row.names <- NULL

# aggregate by order
names <- colnames(DF_1[1:44])
names
DF_2 <- aggregate(. ~Order, FUN = sum, data = DF_1[,colnames(DF_1) %in% names | colnames(DF_1) == 'Order'])
rownames(DF_2) <- DF_2$Order
DF_2$Order <- NULL

# Convert DF to RA matrix
DF_2_RA <- sweep(DF_2, 2, colSums(DF_2), '/')
DF_2_RA[is.na(DF_2_RA)] <- 0

# Merge with metadata
DF_3 <- merge(t(DF_2_RA),metadata_table, by =0,all=F)
rownames(DF_3) <- DF_3$Row.names
DF_3$Row.names <- NULL

# Melt
colnames(DF_3)
DF_4 <- melt(DF_3, id.vars = c("FMT_treatment_batch", "Sample_ID", "endoscopic_outcome", "FMT_treatment_type"))
DF_4$value <- as.numeric(DF_4$value)
#DF_4$FMT_treatment_sucess[DF_4$FMT_treatment_sucess==0] <- "Non-remission"
#DF_4$FMT_treatment_sucess[DF_4$FMT_treatment_sucess==1] <- "Remission"
DF_4$endoscopic_outcome[DF_4$endoscopic_outcome==0] <- "Non-remission"
DF_4$endoscopic_outcome[DF_4$endoscopic_outcome==1] <- "Remission"

ggplot(DF_4) +
  aes(x = variable, y = value, fill = FMT_treatment_type) +
  geom_boxplot() +
  scale_fill_hue(direction = 1) +
  theme_bw() +
  facet_wrap(vars(FMT_treatment_type))

ggplot(DF_4) +
  aes(x = variable, y = value, fill = endoscopic_outcome) +
  geom_boxplot() +
  scale_fill_hue(direction = 1) +
  theme_bw() +
  facet_wrap(vars(FMT_treatment_type))

DF_4_stat <- DF_4[c(DF_4$variable=="Unclassified Caudoviricetes" & DF_4$FMT_treatment_type=="Superdonor FMT"),]
wilcox.test(value ~ endoscopic_outcome, data = DF_4_stat)

DF_4_stat <- DF_4[c(DF_4$variable=="Petitvirales" & DF_4$FMT_treatment_type=="Superdonor FMT"),]
wilcox.test(value ~ endoscopic_outcome, data = DF_4_stat)
####################################
# 16. Host phyla
####################################
#phyloseq_rarefied_samples_D <- phyloseq_rarefied_samples
phyloseq_transposed_Host <- phyloseq_transposed_D

# Extract tables
otu_table_H <- as.data.frame(otu_table(phyloseq_transposed_Host))
metadata_table_H <- as.data.frame(sample_data(phyloseq_transposed_Host))
metadata_table_H <- metadata_table_H[, c("Patient_ID", "FMT_treatment_batch", "FMT_treatment_sucess", "FMT_treatment_type","FMT_treatment_type2","Timepoints")]

# Extract Thost from original host tables
tax_table_H <- Mastertable_viral_rarefied_2[, "Phyla_phylum", drop = FALSE]
#tax_table_H <- Mastertable_viral_rarefied_2_H[, "Host_family1", drop = FALSE]
table(Mastertable_viral_rarefied_2_H$Phyla_phylum)

# Merge with taxonomies
DF_1_H <- merge(otu_table_H,tax_table_H, by = 0, all =F)
rownames(DF_1_H) <- DF_1_H$Row.names
DF_1_H$Row.names <- NULL

# aggregate by order
names <- colnames(otu_table_H[1:44])
DF_2_H <- aggregate(. ~Phyla_phylum, FUN = sum, data = DF_1_H[,colnames(DF_1_H) %in% names | colnames(DF_1_H) == 'Phyla_phylum'])
rownames(DF_2_H) <- DF_2_H$Phyla_phylum
DF_2_H$Phyla_phylum <- NULL

# Check highly abundant genera
DF_2_H$total <- rowSums(DF_2_H)
sorted_DF_2_H <- DF_2_H[rev(order(DF_2_H$total)), ]
print(sorted_DF_2_H)

# Convert DF to RA matrix
DF_2_H_RA <- sweep(DF_2_H, 2, colSums(DF_2_H), '/')
DF_2_H_RA[is.na(DF_2_H_RA)] <- 0

# Merge with metadata
DF_3_H <- merge(t(DF_2_H_RA),metadata_table_H, by =0,all=F)
rownames(DF_3_H) <- DF_3_H$Row.names
DF_3_H$Row.names <- NULL

# Melt
colnames(DF_3_H)
DF_4_H <- melt(DF_3_H, id.vars = c("Patient_ID", "FMT_treatment_batch", "FMT_treatment_sucess", "FMT_treatment_type","FMT_treatment_type2","Timepoints"))
DF_4_H <- melt(DF_3_H, id.vars = c("FMT_treatment_batch", "Sample_ID", "endoscopic_outcome", "FMT_treatment_type"))
DF_4_H$value <- as.numeric(DF_4_H$value)
DF_4_H$FMT_treatment_sucess[DF_4_H$FMT_treatment_sucess==0] <- "Non-remission"
DF_4_H$FMT_treatment_sucess[DF_4_H$FMT_treatment_sucess==1] <- "Remission"
DF_4_H$endoscopic_outcome[DF_4_H$endoscopic_outcome==0] <- "Non-remission"
DF_4_H$endoscopic_outcome[DF_4_H$endoscopic_outcome==1] <- "Remission"

# Figure host-infecting phages
DF_4_H %>%
  filter(variable %in% c("Actinomycetota", "Bacillota", "Bacteroidota", "Verrucomicrobiota", "Pseudomonadota", "Unknown")) %>%
  ggplot() +
  aes(x = variable, y = value, fill = endoscopic_outcome) +
  geom_boxplot() +
  scale_fill_hue(direction = 1) +
  theme_bw() +
  facet_wrap(vars(FMT_treatment_type))

DF_4_H %>%
  filter(variable %in% c("Akkermansiaceae","Bacteroidaceae","Bifidobacteriaceae","Burkholderiaceae","Enterobacteriaceae","Lachnospiraceae","Oscillospiraceae","Lachnospiraceae","Marinifilaceae","Rikenellaceae","Ruminococcaceae")) %>%
  ggplot() +
  aes(x = variable, y = value, fill = endoscopic_outcome) +
  geom_boxplot() +
  scale_fill_hue(direction = 1) +
  theme_bw() +
  facet_wrap(vars(FMT_treatment_type))

DF_4_H %>%
  filter(variable %in% c("Prevotella","Bacteroides","Alistipes","Escherichia","Bifidobacterium","Faecalibacterium","Parabacteroides","Blautia","Phocaeicola","Collinsella","Odoribacter")) %>%
  ggplot() +
  aes(x = variable, y = value, fill = endoscopic_outcome) +
  geom_boxplot() +
  scale_fill_hue(direction = 1) +
  theme_bw() +
  facet_wrap(vars(FMT_treatment_type))


DF_4_H_Stat <- DF_4_H[DF_4_H$FMT_treatment_type=="Superdonor FMT" & DF_4_H$variable=="Pseudomonadota",]
wilcox.test(value ~ endoscopic_outcome,data=DF_4_H)
esquisser(DF_4_H)
####################################
# 17. What viral genera are transplanted?
####################################
# 17.1 Keep only donor-patient pairs
####################################
# Autologous donor + patient pairs
list_autologous_keep <- c("AP83", "AP78","AP70","AP69","AP68","AP66","AP61","AP58","AP53","AP52","AP50","AP46","AP41","AP23","AP22","AP20","AP12","AP11","AP1")

# Superdonor + patients pairs
list_superdonor_keep <- c("99","98","90","88","83","80","76","75","74","73","70","67","65","50","40","37","36","31","27","23","2","15","11","109","1")

# Combine and get the sample names to keep from the list
combined_list_pair <- c(list_autologous_keep, list_superdonor_keep)

# REMOVE THIS WHEN GENUS: Extract the samples that match the "FMT_treatment_batch" in the metadata
#tax_table_original <- tax_table(phyloseq_rarefied_samples)
#tax_table_swapped <- tax_table_original
#tax_table_swapped[, c('Order', 'Genus')] <- tax_table_swapped[, c('Genus', 'Order')]
#physeq_swapped <- phyloseq(otu_table(phyloseq_rarefied_samples), tax_table_swapped, sample_data(phyloseq_rarefied_samples))
#tax_table(physeq_swapped)[tax_table(physeq_swapped)[,"Genus"]== "","Genus"] <- "Unmatched_genus"

#phyloseq_rarefied_phages_1 <- physeq_swapped %>%
#  aggregate_rare(level = "Genus", detection = 0.000001/100, prevalence = 0.00001/100)

# GENUS
phyloseq_rarefied_phages <- phyloseq_rarefied_phages
phyloseq_rarefied_phages_subset <- subset_samples(phyloseq_rarefied_phages, FMT_treatment_batch %in% combined_list_pair)

#View(sample_data(phyloseq_rarefied_phages_subset))
####################################
# 17.2. Extract donor samples
####################################
# Create donor/patient list
list_donor <- c("donor")
list_patient <- c("patient")

# Subset donor samples
phyloseq_phages_donor <- subset_samples(phyloseq_rarefied_phages_subset, Disease_status %in% list_donor)

# Subset patient samples
phyloseq_phages_patient <- subset_samples(phyloseq_rarefied_phages_subset, Disease_status %in% list_patient)
####################################
# 17.3. Aggregate donors samples by specific donor batch 
####################################
Donor_samples <- merge_samples(phyloseq_phages_donor, "FMT_treatment_batch",fun=mean)

# Transpose the otu_table
transposed_otu_table <- t(otu_table(Donor_samples))

# Create a new phyloseq object with the transposed otu_table
phyloseq_transposed <- phyloseq(tax_table(Donor_samples),
                                sample_data(Donor_samples),
                                transposed_otu_table)

# Add FMT_treatment_batch names by rownames
sample_data(phyloseq_transposed)$FMT_treatment_batch <- rownames(sample_data(phyloseq_transposed))

# Only keep this column
phyloseq_transposed <- subset_samples(phyloseq_transposed, select = FMT_treatment_batch)

# Add column with new names
donor_numbering <- 1:length(sample_data(phyloseq_transposed)$FMT_treatment_batch)
sample_data(phyloseq_transposed)$Sample_ID <- paste0("Donor_Batch", donor_numbering)

# Add new sample names to phyloseq object
sample_names(phyloseq_transposed) <- sample_data(phyloseq_transposed)$Sample_ID
####################################
# 17.4. Combine patients and donors
####################################
# Add column FMT_treatment_type2
sample_data(phyloseq_transposed)$FMT_treatment_type2 <- rep('Donor (-superdonor)', nrow(sample_data(phyloseq_transposed)))

# Merge the two phyloseq objects
combined_phyloseq <- merge_phyloseq(phyloseq_transposed, phyloseq_phages_patient)

# Calculate RA matrix
phyloseq_phages_donor_patient <- combined_phyloseq %>% 
  microbiome::transform(transform = "identity") %>%
  aggregate_rare(level = "Genus", detection = 0.0001/100, prevalence = 0.001/100)

phyloseq_without_other <- subset_taxa(phyloseq_phages_donor_patient, !taxa_names(phyloseq_phages_donor_patient) %in% c("Other"))
phyloseq_phages_donor_patient <- phyloseq_without_other
####################################
# 17.5 Keep patient (-superdonor) samples
####################################
# Select patient (-superdonor) sample
phyloseq_phages_superdonor <- subset_samples(phyloseq_phages_donor_patient, FMT_treatment_batch %in% list_superdonor_keep)

# Assuming phyloseq_phages_superdonor is your phyloseq object
sample_data_subset <- subset_samples(phyloseq_phages_superdonor, grepl("_w4|_w0|Donor", Sample_ID))

# Discard patient with no baseline and week 4 samples
list_superdonor_discard <- c("15","31","40","67","70","99","36")
#list_superdonor_discard <- c("15","31","40","67","70","99")
phyloseq_phages_superdonor1 <- subset_samples(sample_data_subset, !FMT_treatment_batch %in% list_superdonor_discard)
####################################
# 17.6 Rarefy samples
####################################
colSums(otu_table(phyloseq_phages_superdonor1))
read_counts <- colSums(otu_table(phyloseq_phages_superdonor1))
ordered_samples <- names(read_counts)[order(read_counts)]

# Print the ordered sample names and their corresponding read counts
for (sample in ordered_samples) {
  print(paste("Sample:", sample, "- Read Count:", read_counts[sample]))
}

# Assuming phyloseq_phages_superdonor1 is your phyloseq object
lowest_read_count <- min(colSums(otu_table(phyloseq_phages_superdonor1)))
phyloseq_rarefied_superdonor <- rarefy_even_depth(phyloseq_phages_superdonor1, sample.size = lowest_read_count)

# Print the rarefied phyloseq object
print(phyloseq_rarefied_superdonor)
####################################
# 17.7 Create a list of phyloseq table for each donor batch - patient pair
####################################
# Get unique values of 'FMT_treatment_batch'
unique_batches <- unique(sample_data(phyloseq_rarefied_superdonor)$FMT_treatment_batch)

# Create an empty list to store the subsetted phyloseq objects
subsetted_phyloseqs <- list()

# Loop through each unique value of 'FMT_treatment_batch'
for (batch in unique_batches) {
  # Subset samples with the current 'FMT_treatment_batch' value
  subset_phyloseq <- subset_samples(phyloseq_rarefied_superdonor, FMT_treatment_batch == batch)
  
  # Store the subsetted phyloseq object in the list
  subsetted_phyloseqs[[as.character(batch)]] <- subset_phyloseq
}
####################################
# 17.8 Remove Taxa with zero abundance in each list
####################################
# Loop through each phyloseq object and filter rows with zero row sums in OTU table
for (i in seq_along(subsetted_phyloseqs)) {
  otu_table <- otu_table(subsetted_phyloseqs[[i]])
  
  # Calculate row sums of the OTU table
  row_sums <- rowSums(otu_table)
  
  # Keep only rows with non-zero row sums
  filtered_otu_table <- otu_table[row_sums > 0, ]
  
  # Update the OTU table in the current phyloseq object
  subsetted_phyloseqs[[i]] <- phyloseq(sample_data(subsetted_phyloseqs[[i]]), filtered_otu_table, tax_table(subsetted_phyloseqs[[i]]))
}

subsetted_phyloseqs_trimmed <- subsetted_phyloseqs
# The physeq_list now contains the filtered phyloseq objects with rows removed from OTU tables where the row sums were zero.
####################################
# 17.9 Obtain dataframes out of phyloseq objects 
####################################
# 1) Create an empty list to store the otu_tables
otu_tables_list <- list()

# 1) Loop through each phyloseq object in the original list
for (i in seq_along(subsetted_phyloseqs_trimmed)) {
  # Extract the otu_table from the ith phyloseq object
  otu_table_i <- as.data.frame(otu_table(subsetted_phyloseqs_trimmed[[i]]))
  
  # Add the otu_table to the new list
  otu_tables_list[[i]] <- otu_table_i
}

# 2) Function to remove taxa not present in the first column
remove_non_first_column_taxa <- function(df) {
# Identify taxa present in the first column
first_column_taxa <- rownames(df)[df[, 1] != 0]

# Remove the identified taxa not present in the first column from the dataframe
df <- df[rownames(df) %in% first_column_taxa, ]

return(df)
 }

# 2) Apply the function to each dataframe in the list
otu_tables_list <- lapply(otu_tables_list, remove_non_first_column_taxa)

# 3) Remove first column
# 3) Define a function to remove the first column from a data frame
remove_first_column <- function(df) {
  df <- df[, -1]
  return(df)
}

# 3) Use lapply to remove the first column from each data frame in the list
otu_tables_list_removed <- lapply(otu_tables_list, remove_first_column)

# 4) Again remove taxa not present in second column
remove_non_first_column_taxa <- function(df) {
  # Identify taxa present in the first column
  first_column_taxa <- rownames(df)[df[, 2] != 0]
  
  # Remove the identified taxa not present in the first column from the dataframe
  df <- df[rownames(df) %in% first_column_taxa, ]
  
  return(df)
}

# 4) Apply the function to each dataframe in the list
otu_tables_list <- lapply(otu_tables_list_removed, remove_non_first_column_taxa)
otu_tables_list
# Conclusion: every phage in the list is present in the donor and also present at w4
# And thus might have been transplanted for which we can calculate a LogFoldChange
####################################
# 17.10 Add of pseudocount
####################################
# Define a function to add a pseudocount of 1 to a data frame
add_pseudocount <- function(df) {
  df <- df + 1
  return(df)
}

# Use lapply to add a pseudocount of 1 to each data frame in the list
otu_tables_list_pseudocount <- lapply(otu_tables_list, add_pseudocount)

# Print the list of data frames with pseudocount added
print(otu_tables_list_pseudocount)
####################################
# 17.11 Calculate LogFoldChange
####################################
otu_tables_list_pseudocount

# Define a function to calculate log-fold change and add a new column
calculate_logfoldchange <- function(df) {
  df$LogFoldChange <- log2(df[, 2] / df[, 1])
  return(df)
}

# Use lapply to calculate log-fold change and add the new column to each data frame in the list
otu_tables_list_logfoldchange <- lapply(otu_tables_list_pseudocount, calculate_logfoldchange)
####################################
# 17.13 Rename LogFoldChange column
####################################
# Define a function to process each data frame
process_dataframe <- function(df) {
  extracted_col <- df[, 3, drop = FALSE]  # Extract the third column
  new_col_name <- colnames(df)[2]  # Name of the second column
  colnames(extracted_col) <- new_col_name  # Rename the extracted column
  return(extracted_col)
}

# Use lapply to process each data frame and create a list of data frames
extracted_columns_list <- lapply(otu_tables_list_logfoldchange, process_dataframe)

# Print the combined data frame
print(extracted_columns_list)
####################################
# 17.14 Combine in one large dataframe
####################################
# Merge all data frames in the list
merged_data <- extracted_columns_list %>%
  reduce(function(x, y) {
    merged <- merge(x, y, by = 0, all = TRUE)
    rownames(merged) <- merged$Row.names
    merged$Row.names <- NULL
    return(merged)
  })

# Extract some taxa
merged_data[merged_data < 0] <- NA
merged_data$Prevalence <- rowSums(!is.na(merged_data[, -1])) / ncol(merged_data)
merged_data1 <- merged_data

# Extract rows with prevalence above 10%
selected_rows <- merged_data[merged_data$Prevalence > 0.2, ]

# Remove the "Prevalence" column
selected_rows$Prevalence <- NULL

# Convert NA to zero
selected_rows[is.na(selected_rows)] <- 0
selected_rows[selected_rows < 0] <- 0
####################################
# 17.15 Merge with metadata
####################################
# Upload metadata
setwd("/Users/daan/Desktop/Bioinformatics/pre-processing/FMT_study/Metadata/viral_transplant")
dir()

metadata_D <- read_excel("samples_FMT_project_Final_DJ_R.xlsx")
metadata_D <- as.data.frame(metadata_D)
rownames(metadata_D) <- metadata_D$Sample_ID
metadata_D$Sample_ID <- NULL
colnames(metadata_D)
selected_columns <- metadata_D[, c("Patient_ID", "FMT_treatment_type", "FMT_treatment_sucess")]

# merge with data
cor_df_m <- merge(t(selected_rows),selected_columns,by=0,all=F)
rownames(cor_df_m) <- cor_df_m$Row.names
cor_df_m$Row.names <- NULL
####################################
# 17.16 Viral transplant boxplot figure
####################################
# Melt
cor_df_melt <- melt(cor_df_m, id.vars = c("Patient_ID","FMT_treatment_type","FMT_treatment_sucess"))
cor_df_melt$value <- as.numeric(cor_df_melt$value)

# Endoscopic remission
cor_df_melt$FMT_treatment_sucess[cor_df_melt$FMT_treatment_sucess==1] <- "Remission"
cor_df_melt$FMT_treatment_sucess[cor_df_melt$FMT_treatment_sucess==0] <- "Non-remission"
unique(cor_df_melt$variable)

# Figure
ggplot(cor_df_melt) +
 aes(x = variable, y = value, fill = FMT_treatment_sucess) +
 geom_boxplot() +
 scale_fill_hue(direction = 1) +
 theme_bw()

# Statistics
cor_df_melt_15 <- cor_df_melt[cor_df_melt$variable=="Petitvirales",]
wilcox.test(value ~ FMT_treatment_sucess,data=cor_df_melt_15,exact=FALSE)

cor_df_melt_15 <- cor_df_melt[cor_df_melt$variable=="Crassvirales",]
wilcox.test(value ~ FMT_treatment_sucess,data=cor_df_melt_15,exact=FALSE)

cor_df_melt_15 <- cor_df_melt[cor_df_melt$variable=="Unclassified Caudoviricetes",]
wilcox.test(value ~ FMT_treatment_sucess,data=cor_df_melt_15, exact=FALSE)

p_values_sup <- c(0.4293,0.7791,0.02714)
p.adjust(p_values_sup, method = "BH")
####################################
# 17.17 Heatmap figure
####################################
cor_df_m$Patient_ID <- NULL
cor_df_m$FMT_treatment_type <- NULL
cor_df_m$FMT_treatment_sucess[cor_df_m$FMT_treatment_sucess==1] <- "Remission"
cor_df_m$FMT_treatment_sucess[cor_df_m$FMT_treatment_sucess==0] <- "Non-remission"
merged_data1[merged_data1 < 1] <- 0

fmt_treatment_success <- cor_df_m$FMT_treatment_sucess

# Remove the 'FMT_treatment_success' column from the original dataframe
cor_df_m1 <- cor_df_m[, -ncol(cor_df_m)]
normalized_data <- cor_df_m1 

haDisease <- HeatmapAnnotation(df = data.frame(disease = fmt_treatment_success),
                               col = list(disease = c("Remission" = "lightgreen", "Non-remission" = "orange")), show_annotation_name = FALSE,border = TRUE)

order_samples <- c("P21_w4","P55_w4","P65_w4","P76_w4","P88_w4",
                   "P24_w4","P25_w4", "P38_w4","P47_w4","P49_w4","P56_w4",
                   "P60_w4","P62_w4","P73_w4","P74_w4","P75_w4",
                   "P81_w4","P82_w4")

order_samples <- c("P21_w8","P55_w8","P65_w8","P76_w8","P88_w8","P16_w8",
                   "P24_w8","P25_w8", "P38_w8","P47_w8","P49_w8","P56_w8",
                   "P60_w8","P62_w8","P73_w8","P74_w8","P75_w8",
                   "P81_w8","P82_w8")

fmt_treatment_success
order_samples

#order_superdonor_samples <- c(
#  "superdonor 1","superdonor 2","superdonor 4","superdonor 11","superdonor 1",
#  "superdonor 5","superdonor 12", "superdonor 13","superdonor 10","superdonor 7",
#  "superdonor 9","superdonor 9","superdonor 9","superdonor 9","superdonor 9",
#  "superdonor 10","superdonor 11","superdonor 3")

# haDisease <- HeatmapAnnotation(df = data.frame(superdonors = order_superdonor_samples), show_annotation_name = TRUE,border = TRUE)
library("circlize")
col_fun = colorRamp2(c(0,0.001,4,6,8,10), c("gray98","#ffdb00", "#ffa904","#ee7b06","#a12424","#400b0b"))

# Create a Heatmap object with clustering based on 'FMT_treatment_success_numeric'
Heatmap(
  t(normalized_data),  # Transpose the data to have samples as rows
  col = col_fun,
  name = "Log2FoldChange",
  row_names_side = "left",
  show_row_names = TRUE,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  border = TRUE,
  column_order = order_samples,
  top_annotation = haDisease)
####################################
# 17.18 Prevalence figure
####################################
merged_data1 <- merged_data
merged_data1$Prevalence <- NULL

# Insert threshold of LogFoldIncrease 
# A Log2FoldIncrease of 100% is a doubling of the number of reads for a genera in that samples
# This will be Log2(2)=1 (set a threshold of 1)
# merged_data1[merged_data1 > 0] <- 1
merged_data1[merged_data1 >= 1] <- 1
merged_data1[merged_data1 < 1] <- 0
merged_data1[is.na(merged_data1)] <- 0
merged_data1$Prevalence <- rowSums(merged_data1)
merged_data1$number <- 1
merged_data1$rownames <- rownames(merged_data1) 
#View(merged_data1)

# Create a data frame with counts of shared taxa for each number of individuals
taxa_counts <- as.data.frame(table(merged_data1$Prevalence))
taxa_counts <- taxa_counts[-1,]
taxa_counts$Var1

# Create the plot using ggplot2
ggplot(taxa_counts, aes(x = Var1, y = Freq)) +
  geom_bar(stat = "identity") +
  labs(x = "Number of Individuals", y = "Number of Transplanted Taxa",
       title = "Distribution of Transplanted Taxa Between Individuals")

ggplot(taxa_counts) +
 aes(x = Var1, y = Freq) +
 geom_col(fill = "cadetblue") +
 theme_minimal() +
  labs(x = "Number of Individuals", y = "Number of Transplanted Taxa",
       title = "Distribution of Transplanted Taxa Between Individuals")

ggplot(taxa_counts) +
 aes(x = Var1, y = Freq) +
 geom_col(fill = "cadetblue") +
 theme_bw() +
  labs(x = "Number of Individuals", y = "Number of Transplanted Taxa",
  title = "Distribution of Transplanted Taxa Between Individuals") +
  scale_y_continuous(breaks = seq(0, max(taxa_counts$Freq), by = 2))

# Add N on x-axis total of w4 18 samples (19 for w8)
# The idea is that Genus 2 is a driver of VCT2 & mostly transmitted..
# Good link to VCT..

table(Mastertable_viral_rarefied_2_)
table(Mastertable_viral_rarefied_2_H$Phyla_phylum)
table(Mastertable_viral_rarefied_2_H$Host_phyla1[Mastertable_viral_rarefied_2_H$Genus=='Genus5'])
table(Mastertable_viral_rarefied_2_H$Host_genus1[Mastertable_viral_rarefied_2_H$Genus=='Genus15'])

table(Mastertable_viral_rarefied_2_H$`lysogenic cycle`[Mastertable_viral_rarefied_2_H$Genus=='Genus3'])
(3/14)*100

####################################
# 18. What viral species (NR-contigs) are transplanted?
####################################
# 18.1 Keep only donor-patient pairs
####################################
# Autologous donor + patient pairs
list_autologous_keep <- c("AP83", "AP78","AP70","AP69","AP68","AP66","AP61","AP58","AP53","AP52","AP50","AP46","AP41","AP23","AP22","AP20","AP12","AP11","AP1")

# Superdonor + patients pairs
list_superdonor_keep <- c("99","98","90","88","83","80","76","75","74","73","70","67","65","50","40","37","36","31","27","23","2","15","11","109","1")

# Combine and get the sample names to keep from the list
combined_list_pair <- c(list_autologous_keep, list_superdonor_keep)
# REMOVE THIS WHEN GENUS: Extract the samples that match the "FMT_treatment_batch" in the metadata
#tax_table_original <- tax_table(phyloseq_rarefied_samples)
#tax_table_swapped <- tax_table_original
#tax_table_swapped[, c('Order', 'Genus')] <- tax_table_swapped[, c('Genus', 'Order')]
#physeq_swapped <- phyloseq(otu_table(phyloseq_rarefied_samples), tax_table_swapped, sample_data(phyloseq_rarefied_samples))
#tax_table(physeq_swapped)[tax_table(physeq_swapped)[,"Genus"]== "","Genus"] <- "Unmatched_genus"

#phyloseq_rarefied_phages_1 <- physeq_swapped %>%
#  aggregate_rare(level = "Genus", detection = 0.000001/100, prevalence = 0.00001/100)

# NODES/Species
phyloseq_rarefied_phages <- phyloseq_rarefied_samples
phyloseq_rarefied_phages_subset <- subset_samples(phyloseq_rarefied_phages, FMT_treatment_batch %in% combined_list_pair)
####################################
# 18.2. Extract donor samples
####################################
# Create donor/patient list
list_donor <- c("donor")
list_patient <- c("patient")

# Subset donor samples
phyloseq_phages_donor <- subset_samples(phyloseq_rarefied_phages_subset, Disease_status %in% list_donor)

# Subset patient samples
phyloseq_phages_patient <- subset_samples(phyloseq_rarefied_phages_subset, Disease_status %in% list_patient)
####################################
# 18.3. Aggregate donors samples by specific donor batch 
####################################
Donor_samples <- merge_samples(phyloseq_phages_donor, "FMT_treatment_batch",fun=mean)

# Transpose the otu_table
transposed_otu_table <- t(otu_table(Donor_samples))

# Create a new phyloseq object with the transposed otu_table
phyloseq_transposed <- phyloseq(tax_table(Donor_samples),
                                sample_data(Donor_samples),
                                transposed_otu_table)

# Add FMT_treatment_batch names by rownames
sample_data(phyloseq_transposed)$FMT_treatment_batch <- rownames(sample_data(phyloseq_transposed))

# Only keep this column
phyloseq_transposed <- subset_samples(phyloseq_transposed, select = FMT_treatment_batch)

# Add column with new names
donor_numbering <- 1:length(sample_data(phyloseq_transposed)$FMT_treatment_batch)
sample_data(phyloseq_transposed)$Sample_ID <- paste0("Donor_Batch", donor_numbering)

# Add new sample names to phyloseq object
sample_names(phyloseq_transposed) <- sample_data(phyloseq_transposed)$Sample_ID
####################################
# 18.4. Combine patients and donors
####################################
# Add column FMT_treatment_type2
sample_data(phyloseq_transposed)$FMT_treatment_type2 <- rep('Donor (-superdonor)', nrow(sample_data(phyloseq_transposed)))

# Merge the two phyloseq objects
combined_phyloseq <- merge_phyloseq(phyloseq_transposed, phyloseq_phages_patient)

# Calculate RA matrix
#phyloseq_phages_donor_patient <- combined_phyloseq %>% 
#  microbiome::transform(transform = "identity") %>%
#  aggregate_rare(level = "Genus", detection = 0.0001/100, prevalence = 0.001/100)
phyloseq_phages_donor_patient <- combined_phyloseq

phyloseq_without_other <- subset_taxa(phyloseq_phages_donor_patient, !taxa_names(phyloseq_phages_donor_patient) %in% c("Other"))
phyloseq_phages_donor_patient <- phyloseq_without_other

#View(otu_table(phyloseq_phages_donor_patient))
####################################
# 18.5 Keep patient (-superdonor) samples
####################################
# Select patient (-superdonor) sample
phyloseq_phages_superdonor <- subset_samples(phyloseq_phages_donor_patient, FMT_treatment_batch %in% list_superdonor_keep)

# Assuming phyloseq_phages_superdonor is your phyloseq object
sample_data_subset <- subset_samples(phyloseq_phages_superdonor, grepl("_w8|_w0|Donor", Sample_ID))

# Discard patient with no baseline and week 4 samples
#list_superdonor_discard <- c("15","31","40","67","70","99","36")
list_superdonor_discard <- c("15","31","40","67","70","99")
phyloseq_phages_superdonor1 <- subset_samples(sample_data_subset, !FMT_treatment_batch %in% list_superdonor_discard)
####################################
# 18.6 Rarefy samples
####################################
colSums(otu_table(phyloseq_phages_superdonor1))
read_counts <- colSums(otu_table(phyloseq_phages_superdonor1))
ordered_samples <- names(read_counts)[order(read_counts)]

# Print the ordered sample names and their corresponding read counts
for (sample in ordered_samples) {
  print(paste("Sample:", sample, "- Read Count:", read_counts[sample]))
}

# Assuming phyloseq_phages_superdonor1 is your phyloseq object
lowest_read_count <- min(colSums(otu_table(phyloseq_phages_superdonor1)))
phyloseq_rarefied_superdonor <- rarefy_even_depth(phyloseq_phages_superdonor1, sample.size = lowest_read_count)

# Print the rarefied phyloseq object
print(phyloseq_rarefied_superdonor)
####################################
# 18.7 Create a list of phyloseq table for each donor batch - patient pair
####################################
# Get unique values of 'FMT_treatment_batch'
unique_batches <- unique(sample_data(phyloseq_rarefied_superdonor)$FMT_treatment_batch)

# Create an empty list to store the subsetted phyloseq objects
subsetted_phyloseqs <- list()

# Loop through each unique value of 'FMT_treatment_batch'
for (batch in unique_batches) {
  # Subset samples with the current 'FMT_treatment_batch' value
  subset_phyloseq <- subset_samples(phyloseq_rarefied_superdonor, FMT_treatment_batch == batch)
  
  # Store the subsetted phyloseq object in the list
  subsetted_phyloseqs[[as.character(batch)]] <- subset_phyloseq
}
####################################
# 18.8 Remove Taxa with zero abundance in each list
####################################
# Loop through each phyloseq object and filter rows with zero row sums in OTU table
for (i in seq_along(subsetted_phyloseqs)) {
  otu_table <- otu_table(subsetted_phyloseqs[[i]])
  
  # Calculate row sums of the OTU table
  row_sums <- rowSums(otu_table)
  
  # Keep only rows with non-zero row sums
  filtered_otu_table <- otu_table[row_sums > 0, ]
  
  # Update the OTU table in the current phyloseq object
  subsetted_phyloseqs[[i]] <- phyloseq(sample_data(subsetted_phyloseqs[[i]]), filtered_otu_table, tax_table(subsetted_phyloseqs[[i]]))
}

subsetted_phyloseqs_trimmed <- subsetted_phyloseqs

# The physeq_list now contains the filtered phyloseq objects with rows removed from OTU tables where the row sums were zero.
####################################
# 18.9 Obtain dataframes out of phyloseq objects 
####################################
# 1) Create an empty list to store the otu_tables
otu_tables_list <- list()

# 1) Loop through each phyloseq object in the original list
for (i in seq_along(subsetted_phyloseqs_trimmed)) {
  # Extract the otu_table from the ith phyloseq object
  otu_table_i <- as.data.frame(otu_table(subsetted_phyloseqs_trimmed[[i]]))
  
  # Add the otu_table to the new list
  otu_tables_list[[i]] <- otu_table_i
}

# 2) Function to remove taxa not present in the first column
remove_non_first_column_taxa <- function(df) {
  # Identify taxa present in the first column
  first_column_taxa <- rownames(df)[df[, 1] != 0]
  
  # Remove the identified taxa not present in the first column from the dataframe
  df <- df[rownames(df) %in% first_column_taxa, ]
  
  return(df)
}

# 2) Apply the function to each dataframe in the list
otu_tables_list <- lapply(otu_tables_list, remove_non_first_column_taxa)

# 3) Remove first column
# Define a function to remove the first column from a data frame
remove_first_column <- function(df) {
  df <- df[, -1]
  return(df)
}

# 3) Use lapply to remove the first column from each data frame in the list
otu_tables_list_removed <- lapply(otu_tables_list, remove_first_column)
otu_tables_list_removed

## Remove list 7
otu_tables_list_removed[[7]] <- NULL

# 4) Again remove taxa not present in second column
remove_non_first_column_taxa <- function(df) {
  # Identify taxa present in the first column
  first_column_taxa <- rownames(df)[df[, 2] != 0]
  
  # Remove the identified taxa not present in the first column from the dataframe
  df <- df[rownames(df) %in% first_column_taxa, ]
  
  return(df)
}

# 4) Apply the function to each dataframe in the list
otu_tables_list <- lapply(otu_tables_list_removed, remove_non_first_column_taxa)
otu_tables_list

# Conclusion: every phage in the list is present in the donor and also present at w4
# And thus might have been transplanted for which we can calculate a LogFoldChange
####################################
# 18.10 Add of pseudocount
####################################
# Define a function to add a pseudocount of 1 to a data frame
add_pseudocount <- function(df) {
  df <- df + 1
  return(df)
}

# Use lapply to add a pseudocount of 1 to each data frame in the list
otu_tables_list_pseudocount <- lapply(otu_tables_list, add_pseudocount)

# Print the list of data frames with pseudocount added
print(otu_tables_list_pseudocount)
####################################
# 18.11 Calculate LogFoldChange
####################################
otu_tables_list_pseudocount

# Define a function to calculate log-fold change and add a new column
calculate_logfoldchange <- function(df) {
  df$LogFoldChange <- log2(df[, 2] / df[, 1])
  return(df)
}

# Use lapply to calculate log-fold change and add the new column to each data frame in the list
otu_tables_list_logfoldchange <- lapply(otu_tables_list_pseudocount, calculate_logfoldchange)
otu_tables_list_logfoldchange
####################################
# 18.12 Rename LogFoldChange column
####################################
# Define a function to process each data frame
process_dataframe <- function(df) {
  extracted_col <- df[, 3, drop = FALSE]  # Extract the third column
  new_col_name <- colnames(df)[2]  # Name of the second column
  colnames(extracted_col) <- new_col_name  # Rename the extracted column
  return(extracted_col)
}

# Use lapply to process each data frame and create a list of data frames
extracted_columns_list <- lapply(otu_tables_list_logfoldchange, process_dataframe)

# Print the combined data frame
print(extracted_columns_list)
####################################
# 18.13 Combine in one large dataframe
####################################
# Merge all data frames in the list
merged_data <- extracted_columns_list %>%
  reduce(function(x, y) {
    merged <- merge(x, y, by = 0, all = TRUE)
    rownames(merged) <- merged$Row.names
    merged$Row.names <- NULL
    return(merged)
  })

# Extract some taxa
merged_data[merged_data < 1] <- NA
merged_data$Prevalence <- rowSums(!is.na(merged_data[, -1])) / ncol(merged_data)
####################################
# 18.14 Combine in one large dataframe
####################################
merged_data1 <- merged_data
merged_data1$Prevalence <- NULL

# Insert threshold of LogFoldIncrease 
# A Log2FoldIncrease of 100% is a doubling of the number of reads for a genera in that samples
# This will be Log2(2)=1 (set a threshold of 1)
# merged_data1[merged_data1 > 0] <- 1
merged_data1[merged_data1 >= 1] <- 1
merged_data1[merged_data1 < 1] <- 0
merged_data1[is.na(merged_data1)] <- 0
merged_data1$Prevalence <- rowSums(merged_data1)
merged_data1$number <- 1
merged_data1$rownames <- rownames(merged_data1) 

# Create a data frame with counts of shared taxa for each number of individuals
taxa_counts <- as.data.frame(table(merged_data1$Prevalence))
taxa_counts <- taxa_counts[-1,]
taxa_counts$Var1

# Create the plot using ggplot2
ggplot(taxa_counts, aes(x = Var1, y = Freq)) +
  geom_bar(stat = "identity") +
  labs(x = "Number of Individuals", y = "Number of Transplanted Taxa",
       title = "Distribution of Transplanted Taxa Between Individuals")

ggplot(taxa_counts) +
  aes(x = Var1, y = Freq) +
  geom_col(fill = "cadetblue") +
  theme_minimal() +
  labs(x = "Number of Individuals", y = "Number of Transplanted Taxa",
       title = "Distribution of Transplanted Taxa Between Individuals")

ggplot(taxa_counts) +
  aes(x = Var1, y = Freq) +
  geom_col(fill = "cadetblue") +
  theme_bw() +
  labs(x = "Number of Individuals", y = "Number of Transplanted Taxa",
       title = "Distribution of Transplanted Taxa Between Individuals") +
  scale_y_continuous(breaks = seq(0, max(taxa_counts$Freq), by = 2))
####################################
# 18.15 Add Genera
####################################
# Extract Genus column
genus_df <- data.frame(Genus = Taxonomy_table$Genus, row.names = rownames(Taxonomy_table))

# Merge with Output table
merged_data_DF <- merge(merged_data1,genus_df,by=0,all=F)
rownames(merged_data_DF) <- merged_data_DF$Row.names
merged_data_DF$Row.names <- NULL

# Remove columns with only NA's
merged_data_DF <- merged_data_DF[, colSums(is.na(merged_data_DF)) != nrow(merged_data_DF)]

# Add prevalence
merged_data_DF$Prevalence <- NULL

# Select genera
merged_data_DF_Genus2 <- merged_data_DF[merged_data_DF$Genus=="Genus1",]
View(merged_data_DF)



#
# 125 NR-contigs transplanted

####################################
# 18.16 Select Genus of interest
####################################
# Week 4
# Define a vector of Genus values you want to extract
selected_genera <- c("Genus15","Genus2", "Genus3", "Genus4", "Genus5", "Genus22")

# Create an empty data frame to store the combined data
combined_data <- data.frame()

# Loop through each selected Genus value and combine the data
for (genus_value in selected_genera) {
  genus_data <- merged_data_DF %>%
    filter(Genus == genus_value) %>%
    select(-number, -rownames, -Genus) %>%
    as.data.frame()
  
  # Add a 'Genus' column to identify the genus in each data frame
  genus_data$Genus <- genus_value
  
  combined_data <- bind_rows(combined_data, genus_data)
}

# Select names of samples
names <- colnames(combined_data %>%
  select(-ncol(combined_data)))

# Aggregate by genus
table1 <- aggregate(. ~Genus, FUN = sum, data = combined_data[,colnames(combined_data) %in% names | colnames(combined_data) == 'Genus'])
rownames(table1) <- table1$Genus
table1$Genus <- NULL

# change to present (value of 1)
table1[table1>0] <- 1
table1$prevalence <- rowSums(table1)
table1$Genus <- rownames(table1)

# Reorder 'Genus' by 'prevalence' in descending order
table1 <- table1 %>%
  arrange(desc(prevalence)) %>%
  mutate(Genus = factor(Genus, levels = Genus))

# Figure
ggplot(table1) +
 aes(x = Genus, weight = prevalence) +
 geom_bar(fill = "cadetblue") +
 coord_flip() +
 theme_bw() +
  labs(y = "Number of patients (n=18)", x = "Transplanted contigs (w4)",title = "") +
  theme(axis.text.y = element_text(angle = 0, hjust = 1))  # Rotate y-axis labels if needed

# CONCLUSION:
# A Log2FoldIncrease of 100% is a doubling of the number of reads for a genera in that samples
# This will be Log2(2)=1 (set a threshold of 1).
####################################
# Week 8
####################################
# Define a vector of Genus values you want to extract
selected_genera <- c("Genus2","Genus1", "Genus4", "Genus8", "Genus9", "Genus15", "Genus3","Genus13")

# Create an empty data frame to store the combined data
combined_data <- data.frame()

# Loop through each selected Genus value and combine the data
for (genus_value in selected_genera) {
  genus_data <- merged_data_DF %>%
    filter(Genus == genus_value) %>%
    select(-number, -rownames, -Genus) %>%
    as.data.frame()
  
  # Add a 'Genus' column to identify the genus in each data frame
  genus_data$Genus <- genus_value
  
  combined_data <- bind_rows(combined_data, genus_data)
}

# Select names of samples
names <- colnames(combined_data %>%
                    select(-ncol(combined_data)))

# Aggregate by genus
table1 <- aggregate(. ~Genus, FUN = sum, data = combined_data[,colnames(combined_data) %in% names | colnames(combined_data) == 'Genus'])
rownames(table1) <- table1$Genus
table1$Genus <- NULL

# change to present (value of 1)
table1[table1>0] <- 1
table1$prevalence <- rowSums(table1)
table1$Genus <- rownames(table1)

# Reorder 'Genus' by 'prevalence' in descending order
table1 <- table1 %>%
  arrange(desc(prevalence)) %>%
  mutate(Genus = factor(Genus, levels = Genus))

# Figure
ggplot(table1) +
  aes(x = Genus, weight = prevalence) +
  geom_bar(fill = "cadetblue") +
  coord_flip() +
  theme_bw() +
  labs(y = "Number of patients (n=19)", x = "Transplanted contigs (w8)",title = "") +
  theme(axis.text.y = element_text(angle = 0, hjust = 1))  # Rotate y-axis labels if needed
####################################
# 19 Individual metrics (vTax, hTax and VCT)
####################################
# 19.1 Insert VCT tables
####################################
library(utils)

setwd("/Users/daan/Desktop/Transfer")
VCT_1 <- read.table("cluster_data_frame.txt", quote = "\"", comment.char = "",header = TRUE)
rownames(VCT_1) <- VCT_1$Samples
VCT_1$Samples <- NULL

# Convert VCT2 to 
VCT_1$Cluster[VCT_1$Cluster==1] <- 0
VCT_1$Cluster[VCT_1$Cluster==2] <- 1

# Obtain colnames of  samples
names <- colnames(Mastertable_viral_rarefied_2[1:304])

# Download excel with link of samples and donor batch
library(readxl)
setwd("/Users/daan/Desktop/Bioinformatics/pre-processing/FMT_study/Metadata/viral_transplant")
metadata <- as.data.frame(read_excel("samples_FMT_project_Final_DJ_R.xlsx"))
rownames(metadata) <- metadata$Sample_ID
metadata_Batch <- metadata[, c('FMT_treatment_batch', 'Disease_status'), drop = FALSE]
metadata_Batch_donor <- metadata_Batch[metadata_Batch$Disease_status=="donor",]
metadata_Batch_donor$Disease_status <- NULL
metadata_Batch_patient <- metadata_Batch[metadata_Batch$Disease_status=="patient",]
metadata_Batch_patient$Disease_status <- NULL

# Add this excel to otu table
Mastertable_phages_Agg_merged <- merge(VCT_1,metadata_Batch_donor,by=0,all=F)
rownames(Mastertable_phages_Agg_merged) <- Mastertable_phages_Agg_merged$Row.names
Mastertable_phages_Agg_merged$Row.names <- NULL

Mastertable_phages_Agg_merged_P <- merge(VCT_1,metadata_Batch_patient,by=0,all=F)
rownames(Mastertable_phages_Agg_merged_P) <- Mastertable_phages_Agg_merged_P$Row.names
Mastertable_phages_Agg_merged_P$Row.names <- NULL

# FOR DONORS: Calculate the mean for each group within 'FMT_treatment_batch' column, excluding NAs
# Percentage of VCT2
library(dplyr)
grouped_means <- Mastertable_phages_Agg_merged %>%
  group_by(FMT_treatment_batch) %>%
  summarize_all(mean)

# Add column with new names
donor_numbering <- 1:length(grouped_means$FMT_treatment_batch)
grouped_means$Sample_ID <- paste0("Donor_Batch", donor_numbering)
grouped_means <- as.data.frame(grouped_means)
rownames(grouped_means) <- grouped_means$Sample_ID
grouped_means$Sample_ID <- NULL 

# Merge based on similar colummn
merged_data <- bind_rows(grouped_means,Mastertable_phages_Agg_merged_P)

# Keep patients with w8 and baseline present
list_superdonor_discard <- c("15","31","40","67","70","99","4","5","13","16","68","12", "AP67","AP72","AP77", "AP51","AP22","AP52","AP83")
merged_data_vTax <- subset(merged_data, !FMT_treatment_batch %in% list_superdonor_discard)
batch_column <- data.frame(FMT_treatment_batch = merged_data_vTax$FMT_treatment_batch, row.names = rownames(merged_data_vTax))
# 19 patients undergoing superdonor treatment (donor, w0 and w8)
# 16 patients undergoing autologous treatment (donor, w0 and w8)

# Select cluster2
Clusters1 <- merged_data_vTax
Clusters1$FMT_treatment_batch <- NULL
Clusters2 <- Clusters1 %>%
  rename(VCT2 = Cluster)

Clusters2$VCT1 <- 1-Clusters2$VCT2

# Merge with batch number link for (figure later on)
Mastertable_phages_Agg1_figure <- merge(Clusters2,batch_column,by=0,all=F)
rownames(Mastertable_phages_Agg1_figure) <- Mastertable_phages_Agg1_figure$Row.names
colnames(Mastertable_phages_Agg1_figure)[colnames(Mastertable_phages_Agg1_figure) == 'Row.names'] <- 'Sample_ID'

# Add timescale
Mastertable_phages_Agg1_figure$Timepoint <- Mastertable_phages_Agg1_figure$Sample_ID
Mastertable_phages_Agg1_figure$Timepoint <- sub(".*_", "", Mastertable_phages_Agg1_figure$Timepoint)
Mastertable_phages_Agg1_figure$Timepoint <- gsub("^Batch", "Donor", Mastertable_phages_Agg1_figure$Timepoint)
Mastertable_phages_Agg1_figure$Timepoint <- sub("Donor*.", "Donor", Mastertable_phages_Agg1_figure$Timepoint)
Mastertable_phages_Agg1_figure$Timepoint <- sub("Donor*.", "Donor", Mastertable_phages_Agg1_figure$Timepoint)

# Add Donor type
Mastertable_phages_Agg1_figure$Donor_type <- Mastertable_phages_Agg1_figure$FMT_treatment_batch
Mastertable_phages_Agg1_figure$Donor_type <- gsub("^AP", "Autologous", Mastertable_phages_Agg1_figure$Donor_type)
Mastertable_phages_Agg1_figure$Donor_type <- sub("Autologous*.", "Autologous", Mastertable_phages_Agg1_figure$Donor_type)
Mastertable_phages_Agg1_figure$Donor_type <- sub("Autologous*.", "Autologous", Mastertable_phages_Agg1_figure$Donor_type)
Mastertable_phages_Agg1_figure$Donor_type[!Mastertable_phages_Agg1_figure$Donor_type=="Autologous"] <- "Superdonor"

# FIGURE
# Melt
library(reshape2)
Mastertable_phages_Agg1_figure1 <- melt(Mastertable_phages_Agg1_figure, id.vars = c("FMT_treatment_batch","Timepoint","Donor_type","Sample_ID"))
Mastertable_phages_Agg1_figure1$value <- as.numeric(Mastertable_phages_Agg1_figure1$value)

# Define the desired order of levels
desired_order <- c('Donor','w0', 'w4', 'w8', 'w12', 'm6', 'm12')

# Convert the 'variable' column to a factor with the desired order
Mastertable_phages_Agg1_figure1$Timepoint <- factor(Mastertable_phages_Agg1_figure1$Timepoint, levels = desired_order)

# Esquisse figure
library(ggplot2)

ggplot(Mastertable_phages_Agg1_figure1) +
  aes(x = Timepoint, y = value, fill = variable) +
  geom_col() +
 scale_fill_manual(values = c(VCT2 = "#2b8a92", VCT1 = "#bd801a")) +
  theme_bw() +
  facet_wrap(vars(FMT_treatment_batch))
####################################
# 19.2 Insert vTax
####################################
# Obtain colnames of  samples
names <- colnames(Mastertable_viral_rarefied_2[1:304])

# Extract orders
Mastertable_phages_Agg <- aggregate(. ~Final_order, FUN = sum, data = Mastertable_viral_rarefied_2[,colnames(Mastertable_viral_rarefied_2) %in% names | colnames(Mastertable_viral_rarefied_2) == 'Final_order'])
rownames(Mastertable_phages_Agg) <- Mastertable_phages_Agg$Final_order
Mastertable_phages_Agg$Final_order <- NULL

# Select the columns you want to aggregate
Mastertable_phages_Agg <- as.data.frame(t(Mastertable_phages_Agg))
columns_to_aggregate <- c('Durnavirales', 'Tubulavirales', 'Unannotated')
Mastertable_phages_Agg$Other <- rowSums(Mastertable_phages_Agg[columns_to_aggregate], na.rm = TRUE)
cols_to_remove <- which(colnames(Mastertable_phages_Agg) %in% columns_to_aggregate)
Mastertable_phages_Agg <- Mastertable_phages_Agg[, -cols_to_remove]
#Mastertable_phages_Agg <- as.data.frame(t(Mastertable_phages_Agg))

# Download excel with link of samples and donor batch
library(readxl)
setwd("/Users/daan/Desktop/Bioinformatics/pre-processing/FMT_study/Metadata/viral_transplant")
metadata <- as.data.frame(read_excel("samples_FMT_project_Final_DJ_R.xlsx"))
rownames(metadata) <- metadata$Sample_ID
metadata_Batch <- metadata[, c('FMT_treatment_batch', 'Disease_status'), drop = FALSE]
metadata_Batch_donor <- metadata_Batch[metadata_Batch$Disease_status=="donor",]
metadata_Batch_donor$Disease_status <- NULL
metadata_Batch_patient <- metadata_Batch[metadata_Batch$Disease_status=="patient",]
metadata_Batch_patient$Disease_status <- NULL

# Add this excel to otu table
Mastertable_phages_Agg_merged <- merge(Mastertable_phages_Agg,metadata_Batch_donor,by=0,all=F)
rownames(Mastertable_phages_Agg_merged) <- Mastertable_phages_Agg_merged$Row.names
Mastertable_phages_Agg_merged$Row.names <- NULL

Mastertable_phages_Agg_merged_P <- merge(Mastertable_phages_Agg,metadata_Batch_patient,by=0,all=F)
rownames(Mastertable_phages_Agg_merged_P) <- Mastertable_phages_Agg_merged_P$Row.names
Mastertable_phages_Agg_merged_P$Row.names <- NULL

# FOR DONORS: Calculate the mean for each group within 'FMT_treatment_batch' column, excluding NAs
library(dplyr)
grouped_means <- Mastertable_phages_Agg_merged %>%
  group_by(FMT_treatment_batch) %>%
  summarize_all(mean)

# Add column with new names
donor_numbering <- 1:length(grouped_means$FMT_treatment_batch)
grouped_means$Sample_ID <- paste0("Donor_Batch", donor_numbering)
grouped_means <- as.data.frame(grouped_means)
rownames(grouped_means) <- grouped_means$Sample_ID
grouped_means$Sample_ID <- NULL 

# Merge based on similar colummn
merged_data <- bind_rows(grouped_means,Mastertable_phages_Agg_merged_P)

# Keep patients with w8 and baseline present
list_superdonor_discard <- c("15","31","40","67","70","99","4","5","13","16","68","12", "AP67","AP72","AP77", "AP51","AP22","AP52","AP83")
merged_data_vTax <- subset(merged_data, !FMT_treatment_batch %in% list_superdonor_discard)
# 19 patients undergoing superdonor treatment (donor, w0 and w8)
# 16 patients undergoing autologous treatment (donor, w0 and w8)

# Rarefaction
batch_column <- data.frame(FMT_treatment_batch = merged_data_vTax$FMT_treatment_batch, row.names = rownames(merged_data_vTax))
merged_data_vTax1 <- as.data.frame(t(merged_data_vTax[-1]))

read_counts <- colSums(merged_data_vTax1)
ordered_samples <- names(read_counts)[order(read_counts)]

# Print the ordered sample names and their corresponding read counts
for (sample in ordered_samples) {
  print(paste("Sample:", sample, "- Read Count:", read_counts[sample]))
}

# Set the target rarefaction depth (e.g., 7000)
target_depth <- 7299

# Convert your otu_table dataframe to a phyloseq object
library(phyloseq)
otu_table_phyloseq <- otu_table(merged_data_vTax1, taxa_are_rows = TRUE)

# Rarefy the phyloseq object to the specified depth
rarefied_phyloseq <- rarefy_even_depth(otu_table_phyloseq, sample.size = target_depth)

# Convert the rarefied phyloseq object back to a dataframe
rarefied_dataframe <- as.data.frame(rarefied_phyloseq)

# Calculate RA
Mastertable_phages_Agg1 <- sweep(rarefied_dataframe, 2, colSums(rarefied_dataframe), '/')
Mastertable_phages_Agg1[is.na(Mastertable_phages_Agg1)] <- 0

# Transverse Dataframe
Mastertable_phages_Agg1 <- as.data.frame(t(Mastertable_phages_Agg1))

# Merge with batch number link for (figure later on)
Mastertable_phages_Agg1_figure <- merge(Mastertable_phages_Agg1,batch_column,by=0,all=F)
rownames(Mastertable_phages_Agg1_figure) <- Mastertable_phages_Agg1_figure$Row.names
colnames(Mastertable_phages_Agg1_figure)[colnames(Mastertable_phages_Agg1_figure) == 'Row.names'] <- 'Sample_ID'

# Add timescale
Mastertable_phages_Agg1_figure$Timepoint <- Mastertable_phages_Agg1_figure$Sample_ID
Mastertable_phages_Agg1_figure$Timepoint <- sub(".*_", "", Mastertable_phages_Agg1_figure$Timepoint)
Mastertable_phages_Agg1_figure$Timepoint <- gsub("^Batch", "Donor", Mastertable_phages_Agg1_figure$Timepoint)
Mastertable_phages_Agg1_figure$Timepoint <- sub("Donor*.", "Donor", Mastertable_phages_Agg1_figure$Timepoint)
Mastertable_phages_Agg1_figure$Timepoint <- sub("Donor*.", "Donor", Mastertable_phages_Agg1_figure$Timepoint)

# Add Donor type
Mastertable_phages_Agg1_figure$Donor_type <- Mastertable_phages_Agg1_figure$FMT_treatment_batch
Mastertable_phages_Agg1_figure$Donor_type <- gsub("^AP", "Autologous", Mastertable_phages_Agg1_figure$Donor_type)
Mastertable_phages_Agg1_figure$Donor_type <- sub("Autologous*.", "Autologous", Mastertable_phages_Agg1_figure$Donor_type)
Mastertable_phages_Agg1_figure$Donor_type <- sub("Autologous*.", "Autologous", Mastertable_phages_Agg1_figure$Donor_type)
Mastertable_phages_Agg1_figure$Donor_type[!Mastertable_phages_Agg1_figure$Donor_type=="Autologous"] <- "Superdonor"
View(Mastertable_phages_Agg1_figure)

# FIGURE
# Melt
library(reshape2)
Mastertable_phages_Agg1_figure1 <- melt(Mastertable_phages_Agg1_figure, id.vars = c("FMT_treatment_batch","Timepoint","Donor_type","Sample_ID"))
Mastertable_phages_Agg1_figure1$value <- as.numeric(Mastertable_phages_Agg1_figure1$value)

# Define the desired order of levels
desired_order <- c('Donor','w0', 'w4', 'w8', 'w12', 'm6', 'm12')

# Convert the 'variable' column to a factor with the desired order
Mastertable_phages_Agg1_figure1$Timepoint <- factor(Mastertable_phages_Agg1_figure1$Timepoint, levels = desired_order)

# Esquisse figure
library(esquisse)
#esquisser(Mastertable_phages_Agg1_figure1)

library(ggplot2)
ggplot(Mastertable_phages_Agg1_figure1) +
 aes(x = Timepoint, y = value, fill = variable) +
 geom_col() +
 scale_fill_manual(values = c(Crassvirales = "#599b1b", Petitvirales = "#c49a00", `Unclassified Caudoviricetes` = "#c15c59", Other = "gray85")) +
 theme_bw() +
 facet_wrap(vars(FMT_treatment_batch))
####################################
# 19.3 Insert hTax
####################################
# Obtain colnames of  samples
names <- colnames(Mastertable_viral_rarefied_2[1:304])

# Extract orders
Mastertable_phages_Agg <- aggregate(. ~Phyla_phylum, FUN = sum, data = Mastertable_viral_rarefied_2[,colnames(Mastertable_viral_rarefied_2) %in% names | colnames(Mastertable_viral_rarefied_2) == 'Phyla_phylum'])
rownames(Mastertable_phages_Agg) <- Mastertable_phages_Agg$Phyla_phylum
Mastertable_phages_Agg$Phyla_phylum <- NULL

# Select the columns you want to aggregate
Mastertable_phages_Agg <- as.data.frame(t(Mastertable_phages_Agg))
columns_to_aggregate <- c('Campylobacterota', 'Cyanobacteriota', 'Euryarchaeota','Fusobacteriota','Mycoplasmatota','Thermodesulfobacteriota','Verrucomicrobiota', '0', 'Bdellovibrionota')
Mastertable_phages_Agg$Other <- rowSums(Mastertable_phages_Agg[columns_to_aggregate], na.rm = TRUE)
cols_to_remove <- which(colnames(Mastertable_phages_Agg) %in% columns_to_aggregate)
Mastertable_phages_Agg <- Mastertable_phages_Agg[, -cols_to_remove]
View(Mastertable_phages_Agg)

# Download excel with link of samples and donor batch
library(readxl)
setwd("/Users/daan/Desktop/Bioinformatics/pre-processing/FMT_study/Metadata/viral_transplant")
metadata <- as.data.frame(read_excel("samples_FMT_project_Final_DJ_R.xlsx"))
rownames(metadata) <- metadata$Sample_ID
metadata_Batch <- metadata[, c('FMT_treatment_batch', 'Disease_status'), drop = FALSE]
metadata_Batch_donor <- metadata_Batch[metadata_Batch$Disease_status=="donor",]
metadata_Batch_donor$Disease_status <- NULL
metadata_Batch_patient <- metadata_Batch[metadata_Batch$Disease_status=="patient",]
metadata_Batch_patient$Disease_status <- NULL

# Add this excel to otu table
Mastertable_phages_Agg_merged <- merge(Mastertable_phages_Agg,metadata_Batch_donor,by=0,all=F)
rownames(Mastertable_phages_Agg_merged) <- Mastertable_phages_Agg_merged$Row.names
Mastertable_phages_Agg_merged$Row.names <- NULL

Mastertable_phages_Agg_merged_P <- merge(Mastertable_phages_Agg,metadata_Batch_patient,by=0,all=F)
rownames(Mastertable_phages_Agg_merged_P) <- Mastertable_phages_Agg_merged_P$Row.names
Mastertable_phages_Agg_merged_P$Row.names <- NULL

# FOR DONORS: Calculate the mean for each group within 'FMT_treatment_batch' column, excluding NAs
library(dplyr)
grouped_means <- Mastertable_phages_Agg_merged %>%
  group_by(FMT_treatment_batch) %>%
  summarize_all(mean)

# Add column with new names
donor_numbering <- 1:length(grouped_means$FMT_treatment_batch)
grouped_means$Sample_ID <- paste0("Donor_Batch", donor_numbering)
grouped_means <- as.data.frame(grouped_means)
rownames(grouped_means) <- grouped_means$Sample_ID
grouped_means$Sample_ID <- NULL 

# Merge based on similar colummn
merged_data <- bind_rows(grouped_means,Mastertable_phages_Agg_merged_P)

# Keep patients with w8 and baseline present
list_superdonor_discard <- c("15","31","40","67","70","99","4","5","13","16","68","12", "AP67","AP72","AP77", "AP51","AP22","AP52","AP83")
merged_data_vTax <- subset(merged_data, !FMT_treatment_batch %in% list_superdonor_discard)
# 19 patients undergoing superdonor treatment (donor, w0 and w8)
# 16 patients undergoing autologous treatment (donor, w0 and w8)

# Rarefaction
batch_column <- data.frame(FMT_treatment_batch = merged_data_vTax$FMT_treatment_batch, row.names = rownames(merged_data_vTax))
merged_data_vTax1 <- as.data.frame(t(merged_data_vTax[-1]))

read_counts <- colSums(merged_data_vTax1)
ordered_samples <- names(read_counts)[order(read_counts)]

# Print the ordered sample names and their corresponding read counts
for (sample in ordered_samples) {
  print(paste("Sample:", sample, "- Read Count:", read_counts[sample]))
}

# Set the target rarefaction depth (e.g., 7000)
target_depth <- 7299

# Convert your otu_table dataframe to a phyloseq object
library(phyloseq)
otu_table_phyloseq <- otu_table(merged_data_vTax1, taxa_are_rows = TRUE)

# Rarefy the phyloseq object to the specified depth
rarefied_phyloseq <- rarefy_even_depth(otu_table_phyloseq, sample.size = target_depth)

# Convert the rarefied phyloseq object back to a dataframe
rarefied_dataframe <- as.data.frame(rarefied_phyloseq)

# Calculate RA
Mastertable_phages_Agg1 <- sweep(rarefied_dataframe, 2, colSums(rarefied_dataframe), '/')
Mastertable_phages_Agg1[is.na(Mastertable_phages_Agg1)] <- 0

# Transverse Dataframe
Mastertable_phages_Agg1 <- as.data.frame(t(Mastertable_phages_Agg1))

# Merge with batch number link for (figure later on)
Mastertable_phages_Agg1_figure <- merge(Mastertable_phages_Agg1,batch_column,by=0,all=F)
rownames(Mastertable_phages_Agg1_figure) <- Mastertable_phages_Agg1_figure$Row.names
colnames(Mastertable_phages_Agg1_figure)[colnames(Mastertable_phages_Agg1_figure) == 'Row.names'] <- 'Sample_ID'

# Add timescale
Mastertable_phages_Agg1_figure$Timepoint <- Mastertable_phages_Agg1_figure$Sample_ID
Mastertable_phages_Agg1_figure$Timepoint <- sub(".*_", "", Mastertable_phages_Agg1_figure$Timepoint)
Mastertable_phages_Agg1_figure$Timepoint <- gsub("^Batch", "Donor", Mastertable_phages_Agg1_figure$Timepoint)
Mastertable_phages_Agg1_figure$Timepoint <- sub("Donor*.", "Donor", Mastertable_phages_Agg1_figure$Timepoint)
Mastertable_phages_Agg1_figure$Timepoint <- sub("Donor*.", "Donor", Mastertable_phages_Agg1_figure$Timepoint)

# Add Donor type
Mastertable_phages_Agg1_figure$Donor_type <- Mastertable_phages_Agg1_figure$FMT_treatment_batch
Mastertable_phages_Agg1_figure$Donor_type <- gsub("^AP", "Autologous", Mastertable_phages_Agg1_figure$Donor_type)
Mastertable_phages_Agg1_figure$Donor_type <- sub("Autologous*.", "Autologous", Mastertable_phages_Agg1_figure$Donor_type)
Mastertable_phages_Agg1_figure$Donor_type <- sub("Autologous*.", "Autologous", Mastertable_phages_Agg1_figure$Donor_type)
Mastertable_phages_Agg1_figure$Donor_type[!Mastertable_phages_Agg1_figure$Donor_type=="Autologous"] <- "Superdonor"

# FIGURE
# Melt
library(reshape2)
Mastertable_phages_Agg1_figure1 <- melt(Mastertable_phages_Agg1_figure, id.vars = c("FMT_treatment_batch","Timepoint","Donor_type","Sample_ID"))
Mastertable_phages_Agg1_figure1$value <- as.numeric(Mastertable_phages_Agg1_figure1$value)

# Define the desired order of levels
desired_order <- c('Donor','w0', 'w4', 'w8', 'w12', 'm6', 'm12')

# Convert the 'variable' column to a factor with the desired order
Mastertable_phages_Agg1_figure1$Timepoint <- factor(Mastertable_phages_Agg1_figure1$Timepoint, levels = desired_order)

# Esquisse figure
library(ggplot2)
ggplot(Mastertable_phages_Agg1_figure1) +
  aes(x = Timepoint, y = value, fill = variable) +
  geom_col() +
  scale_fill_manual(values = c(Actinomycetota = "#d89000", Bacillota = "#a3a500", `Bacteroidota` = "#39b600", Chlamydiota = "gray85",Pseudomonadota = "#e76bf3", Other = "gray85")) +
  theme_bw() +
  facet_wrap(vars(FMT_treatment_batch))
####################################

####################################

