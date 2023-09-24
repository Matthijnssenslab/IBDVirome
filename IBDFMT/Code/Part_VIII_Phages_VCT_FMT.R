####################################
# SCRIPT 5: PROKARYOTIC VIRUSES_2
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
library(lme4)
library(multcomp)
####################################
# 1. Datasets
####################################
phyloseq_rarefied_phages
####################################
# 2. Community typing
####################################
# 2.1 Data trimming 
####################################
phyloseq_rarefied_samples1 <- phyloseq_rarefied_phages

# Filter taxa based on abundance threshold
phyloseq_phages_DMM <- phyloseq_rarefied_samples1 %>%
  aggregate_rare(level = "Genus", detection = 0.001, prevalence = 1/100) 


# Prevalence threshold at 1% to increase accuracy because multiple samples are from the same individual (1.75% to be exact).
###################################
# 2.2 DMM algorithm
###################################
# Pick the OTU count matrix and convert it into samples x taxa format
dat <- abundances(phyloseq_phages_DMM)
count <- as.data.frame(t(dat))

## Remove 'Other' meaning all the unique viruses appearing only in less than X% of samples (X% abundance)
count$total <- rowSums(count)
count <- count[!count$total == 0,]
count$total <- NULL
count$Other <- NULL

## Remove samples composed of zero shared viruses
count$total <- rowSums(count)
count <- count[!count$total == 0,]
sort(rowSums(count))
count$total <- NULL

## How many samples and reads are left for determining DMM clusters?
nrow(count) # 296 samples left for which community-typing is possible
sum(count) # 296 samples having 2.96 billion reads
(sum(count)/sum(Mastertable_phage$Totalnumberofreads))*100 # or 94.6% of the reads were shared

## Median genera per sample using community typing
count_1 <- count
count_1 <- (count_1/100)-0.01
count_1[count_1 > 0] <- 1
count_1 <- as.data.frame(count_1)
count_1$prevalence <- rowSums(count_1)
ncol(count_1) # 257 genera shared over 1% of samples with 0.01% abundance
count_1$prevalence
median(count_1$prevalence) # Median of 17 genera  per sample
mean(count_1$prevalence) # Average 17.1 genera per sample

# Determine DMM clusters and evaluate BIC score
# You're matrix cannot have decimals after the komma
count <- as.matrix(count)
#fit <- lapply(1:5, dmn, count = count, verbose=TRUE) 
bic  <- sapply(fit, BIC) # BIC 
aic  <- sapply(fit, AIC) # AIC 
lplc  <- sapply(fit, laplace) # laplace 
plot(bic, type="b", xlab="Number of Dirichlet Components", ylab="Model Fit")
plot(aic, type="b", xlab="Number of Dirichlet Components", ylab="Model Fit")
plot(lplc, type="b", xlab="Number of Dirichlet Components", ylab="Model Fit")
plot(lplc, type="b", xlab="Number of Dirichlet Components", ylab="Model Fit")

# Determine the optimal model
best <- fit[[which.min(unlist(bic))]]
mixturewt(best)
best@group
best@mixture

# Additionally, evaluate the probability of samples belonging to one of the two DMM groups. If not that interesting you could leave this out
group_probability <- as.data.frame(best@group)
mean(group_probability$V1)
mean(group_probability$V2)
median(group_probability$V1)
median(group_probability$V2)

# Subselect the clusters for all samples as determined before
cluster <- apply(mixture(best), 1, which.max)
cluster <- as.data.frame(cluster)
cluster$cluster <- as.character(cluster$cluster)
table(cluster$cluster) # cluster 1; 216 samples, cluster 2 80 samples.
cluster_phylo <- sample_data(cluster)

# At last find drivers or contribution of each taxonomic group to each component
for (k in seq(ncol(fitted(best)))) {
  d <- melt(fitted(best))
  colnames(d) <- c("OTU", "cluster", "value")
  d <- subset(d, cluster == k) %>%
    # Arrange OTUs by assignment strength
    arrange(value) %>%
    mutate(OTU = factor(OTU, levels = unique(OTU))) %>%
    # Only show the most important drivers
    filter(abs(value) > quantile(abs(value), 0.90)) #0.99 default     
  
  p <- ggplot(d, aes(x = OTU, y = value)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(title = paste("Top drivers: community type", k))
  print(p)
}

# Fit genera to model clusters
# A fitted value is a statistical model’s prediction of the mean response value when you input the values of the predictors, factor levels, or components into the model. Suppose you have the following regression equation: y = 3X + 5. If you enter a value of 5 for the predictor, the fitted value is 20. Fitted values are also called predicted values.
# fitted is a generic function which extracts fitted values from objects returned by modeling functions.

contribution_viruses <- as.data.frame(fitted(best))
contribution_viruses$cluster1 <- contribution_viruses$V1
contribution_viruses$cluster2 <- contribution_viruses$V2
contribution_viruses$V1 <- NULL
contribution_viruses$V2 <- NULL
###################################
# 2.3 Add clusters to phyloseq object
###################################
phyloseq_rarefied_phages_visualization_1 <- merge_phyloseq(phyloseq_phages_DMM, cluster_phylo)
###################################
# 2.4 Determine number of genera in samples
###################################
genera <- as.data.frame(t(abundances(phyloseq_rarefied_phages_visualization_1)))
genera[genera > 0] <- 1
genera$prevalence <- rowSums(genera)

# generic infiormation
# A) Min and maximum genera in samples
range(genera$prevalence) # 3 min; 74 maximum

# B) genera in total
nrow(tax_table(phyloseq_rarefied_phages_visualization_1)) # 874 genera

# C) Median genera per sample
median(genera$prevalence) # 26 median
mean(genera$prevalence) # 27.23 average
###################################
# 2.5 Bray-curtis dissimilarity and visualization
###################################
# 2.5.1 Remove samples for which no community-type can be assigned
###################################
# Remove NA's for clusters 
phyloseq_rarefied_phages_visualization_2 = subset_samples(phyloseq_rarefied_phages_visualization_1, !cluster == "NA")
###################################
# 2.5.2 Hellinger transformation
###################################
phyloseq_rarefied_phages_transformed <- phyloseq_rarefied_phages_visualization_2 %>%
  microbiome::transform(transform = "hellinger") %>%
  aggregate_rare(level = "Genus", detection = 0.000001/100, prevalence = 0.000001/100)
###################################
# 2.5.3 Bray-Curtis dissimmilarities
###################################
phyloseq_rarefied_phages_transformed_bray <- phyloseq::distance(phyloseq_rarefied_phages_transformed, method="bray")
###################################
# 2.5.4 Reduce dimensionality and plot virome composition on PcoA
###################################
phyloseq_rarefied_phages_PcoA <- ordinate(phyloseq_rarefied_phages_transformed, method="PCoA", distance=phyloseq_rarefied_phages_transformed_bray)
#View(sample_data(phyloseq_rarefied_phages_transformed))
###################################
# 2.6 Visualization of microbiome using BC dissimilarities and PcOA plot
###################################
# Phyloseq: phyloseq_rarefied_phages_transformed
# Bray-curtis data: phyloseq_rarefied_phages_transformed_bray
# Ordinationation: phyloseq_rarefied_phages_PcoA

metadata_begin <- as.data.frame(sample_data(phyloseq_rarefied_phages_transformed))
#View(metadata_begin)

# PcoA 1: without elipses
plot_ordination(phyloseq_rarefied_phages_transformed, phyloseq_rarefied_phages_PcoA, color="cluster") +
  geom_point(size = 1, alpha=1) +
  scale_shape_manual(values = c(16,17,18)) +
  theme_bw() +
  scale_color_manual(values = c(`1` = "#bd7f1c", `2` = "#2a8a92",`3` = "green")) +
  geom_hline(yintercept = 0, linetype = 3, colour = "black") + 
  geom_vline(xintercept = 0, linetype = 3, colour = "black") +
  xlab("PC1 (8.6%)") +
  ylab("PC2 (4.1%)")
###################################
# 3. Barplots of metadata and cluster assignment
###################################
# 3.1 Determine VCT assignment (probability)
###################################
group_probability
names(group_probability)[names(group_probability) == "V1"] <- "cluster 1"
names(group_probability)[names(group_probability) == "V2"] <- "cluster 2"
names(group_probability)[names(group_probability) == "V3"] <- "cluster 3"
names(group_probability)[names(group_probability) == "V4"] <- "cluster 4"
group_probability <- merge(group_probability,cluster, by = 0, all = F)
rownames(group_probability) <- group_probability$Row.names
names(group_probability)[names(group_probability) == "Row.names"] <- "unique_identifier"
names(group_probability)[names(group_probability) == "cluster"] <- "selected_cluster"
group_probability$probability_selection <- "0"
group_probability$probability_selection[group_probability$selected_cluster == "1"] <- group_probability$`cluster 1`[group_probability$selected_cluster == "1"]
group_probability$probability_selection[group_probability$selected_cluster == "2"] <- group_probability$`cluster 2`[group_probability$selected_cluster == "2"]
group_probability$probability_selection[group_probability$selected_cluster == "3"] <- group_probability$`cluster 2`[group_probability$selected_cluster == "3"]
group_probability$probability_selection[group_probability$selected_cluster == "4"] <- group_probability$`cluster 2`[group_probability$selected_cluster == "4"]
group_probability$probability_selection <- as.numeric(group_probability$probability_selection)
mean(group_probability$probability_selection)
median(group_probability$probability_selection)
#View(group_probability)

#library("writexl")
#dir()
#setwd("/Users/daan/Desktop/Bioinformatics/Analysis/FMT/Figures/Figure3")
#write_xlsx(group_probability,"./group_probability.xlsx")
###################################
# 3.2 Associate FMT_treatment_type with VCT
###################################
# Input:
phyloseq_rarefied_samples1

# Merge clusters
phyloseq_rarefied_phages_3 <- merge_phyloseq(phyloseq_rarefied_samples1, cluster_phylo)

# Remove NA's: samples without cluster (n = 301)
nrow(cluster) # All samples (n=301) have a cluster
phyloseq_rarefied_phages_3 # 301 samples
phyloseq_rarefied_phages_3 = subset_samples(phyloseq_rarefied_phages_3, !cluster == "NA")

# Subset metadata
phyloseq_rarefied_phages_3_stat <- as.data.frame(sample_data(phyloseq_rarefied_phages_3))
phyloseq_rarefied_phages_3_stat <- phyloseq_rarefied_phages_3_stat[,c(2,4,6,7,11,25,32)]
colnames(phyloseq_rarefied_phages_3_stat)

# Level FMT Treatment Type2
phyloseq_rarefied_phages_3_stat$FMT_treatment_type2 <- factor(phyloseq_rarefied_phages_3_stat$FMT_treatment_type2,levels = c('Donor (-superdonor)', 'Patient (-superdonor)', "Donor (-autologous)" ,'Patient (-autologous)', 'Patient (-untreated)'))

# Figure 1: FMT Treatment Type 2 versus cluster
ggplot(phyloseq_rarefied_phages_3_stat) +
  aes(x = FMT_treatment_type2, fill = cluster) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = c(`1` = "#bd7f1b", `2` = "#2a8a92")) +
  theme_bw()

# Statistics
# Extract the response variable 'cluster' from the 'phyloseq_rarefied_phages_3_stat' phyloseq object
response_variable <- sample_data(phyloseq_rarefied_phages_3_stat)$cluster

# Convert the response variable to binary: replace "1" with "0" and "2" with "1"
response_variable <- ifelse(response_variable == "1", "0", "1")

# Convert the response variable to numeric type (0 and 1)
response_variable <- as.numeric(response_variable)

# Extract the confounding variable 'Patient_ID' from the 'phyloseq_rarefied_phages_3_stat' phyloseq object
confounding_variable <- sample_data(phyloseq_rarefied_phages_3_stat)$Patient_ID

# Extract the predictor variable 'FMT_treatment_type2' from the 'phyloseq_rarefied_phages_3_stat' phyloseq object
predictor_variable <- sample_data(phyloseq_rarefied_phages_3_stat)$FMT_treatment_type2

# Combine the response, predictor, and confounding variables into a new data frame 'data_df'
data_df <- data.frame(cluster = response_variable,
                      FMT_treatment_type2 = predictor_variable,
                      Patient_ID = confounding_variable)

# Convert FMT_treatment_type2 to a factor and reorder the levels; Donor (-superdonor)' will be the reference
data_df$FMT_treatment_type2 <- factor(data_df$FMT_treatment_type2, levels = c('Donor (-superdonor)', 'Donor (-autologous)', 'Patient (-autologous)', 'Patient (-superdonor)', 'Patient (-untreated)'))

# Full model
null_model <- glmer(cluster ~ 1 + (1 | Patient_ID), data = data_df, family = binomial("logit"), nAGQ = 100) # Fit the null model with only the intercept
full_model <- glmer(cluster ~ FMT_treatment_type2 + (1 | Patient_ID), data = data_df, family = binomial("logit"), nAGQ = 100) # Fit the mixed-effects logistic regression model
lrt_result <- anova(null_model, full_model) # Perform the Likelihood Ratio Test
lrt_result
full_model

# Comparison with reference variable (# compared to reference = Donor (-superdonor))
model <- glmer(cluster ~ FMT_treatment_type2 + (1 | Patient_ID), data = data_df, family = binomial("logit"),nAGQ = 100)
summary(model) 

# Perform Tukey HSD post hoc test (# all versus all comparison)
tukey_result <- glht(model, linfct = mcp(FMT_treatment_type2 = "Tukey"))
summary(tukey_result)

# Obtain _values
coeff_summary <- coef(summary(tukey_result))
p_values_adjusted <- summary(tukey_result)$test$pvalues
combined_superdonor <- cbind(coeff_summary,p_values_adjusted)
count_point <- data_df
n_per_timepoint <- as.data.frame(table(count_point$FMT_treatment_type2, useNA = "no"))
rownames(n_per_timepoint) <- n_per_timepoint$Var1
n_per_timepoint$Var1 <- NULL
n_per_timepoint
combined_superdonor

# Conclusion:
# Statistically significant differences between These groups driven by Donor (-superdonors)

# The glht() function from the multcomp library performs the Tukey HSD post hoc test. The linfct argument specifies that we want to perform the post hoc test on the FMT_treatment_type2 variable using the Tukey method.
# The summary of tukey_result will show the pairwise comparisons between the levels of FMT_treatment_type2, along with the adjusted p-values (or confidence intervals) for each comparison. Statistically significant comparisons indicate that the log-odds of the response variable cluster significantly differ between the corresponding levels.
# Keep in mind that with multiple comparisons, the risk of Type I error increases. Therefore, it's essential to interpret the results cautiously and consider adjusting the significance level (e.g., using Bonferroni correction) to control for the family-wise error rate if you are making a large number of comparisons.
###################################
# 3.3 Associate patients (-superdonor) with VCT over time
################################### 
# Figure 2: FMT Treatment Type 2 versus cluster divided by Treatment success for superdonors
df_percent <- phyloseq_rarefied_phages_3_stat %>%
  group_by(cluster) %>%
  mutate(Percentage = n() / sum(n()) * 100)

df_percent <- df_percent[df_percent$FMT_treatment_type2=="Patient (-superdonor)" | df_percent$FMT_treatment_type2=="Patient (-untreated)",]
list_patient <- c("patient 1","patient 11","patient 12","patient 20","patient 41","patient 46","patient 50","patient 51","patient 52","patient 58","patient 67","patient 66","patient 68","patient 69","patient 72","patient 77","patient 78","patient 83")

filtered_df_percent <- df_percent[!(df_percent$Patient_ID %in% list_patient), ]
filtered_df_percent$Timepoints <- factor(filtered_df_percent$Timepoints,levels = c("Baseline","w4","w8","w12","m6","m12"))
filtered_df_percent$FMT_treatment_type2 <- factor(filtered_df_percent$FMT_treatment_type2,levels = c("Patient (-untreated)",'Patient (-superdonor)'))

# Figure 1: comparison VCT over time (only Donor (-superdonors))
ggplot(filtered_df_percent) +
 aes(x = Timepoints, fill = cluster) +
 geom_bar(position = "fill") +
 scale_fill_manual(values = c(`1` = "#bd7f1b", `2` = "#2a8a92")) +
 theme_bw()

filtered_df_percent1 <- filtered_df_percent
response_variable <- filtered_df_percent1$cluster
response_variable <- ifelse(response_variable == "1", "0", "1")
response_variable <- as.numeric(response_variable)
confounding_variable <- filtered_df_percent1$Patient_ID
predictor_variable <- filtered_df_percent1$Timepoints
data_df <- data.frame(cluster = response_variable, Timepoints = predictor_variable,Patient_ID = confounding_variable)

# Full model
view(data_df)
null_model <- glmer(cluster ~ 1 + (1 | Patient_ID), data = data_df, family = binomial("logit")) # Fit the null model with only the intercept
full_model <- glmer(cluster ~ Timepoints + (1 | Patient_ID), data = data_df, family = binomial("logit"), nAGQ = 100) # Fit the mixed-effects logistic regression model
lrt_result <- anova(null_model, full_model) # Perform the Likelihood Ratio Test
lrt_result

# Conclusion:
# No difference in VCT over time
################################### 
# 3.4 Associate patients (-superdonor) with VCT divided by endoscopic outcome - Fig 1
################################### 
# Figure 2: comparison VCT time before/after treatment between response and non-response
ggplot(filtered_df_percent) +
 aes(x = FMT_treatment_type2, fill = cluster) +
 geom_bar(position = "fill") +
 scale_fill_manual(values = c(`1` = "#bd7f1b", `2` = "#2a8a92")) +
 theme_bw() +
 facet_wrap(vars(Endoscopic_outcome))

# Statistics
# Remission
filtered_df_percent_R <- filtered_df_percent[filtered_df_percent$Endoscopic_outcome=="1",]
response_variable <- filtered_df_percent_R$cluster
response_variable <- ifelse(response_variable == "1", "0", "1")
response_variable <- as.numeric(response_variable)
confounding_variable <- filtered_df_percent_R$Patient_ID
predictor_variable <- filtered_df_percent_R$FMT_treatment_type2
data_df <- data.frame(cluster = response_variable, FMT_treatment_type2 = predictor_variable,Patient_ID = confounding_variable)
View(data_df)
# Full model
null_model <- glmer(cluster ~ 1 + (1 | Patient_ID), data = data_df, family = binomial) # Fit the null model with only the intercept
full_model <- glmer(cluster ~ FMT_treatment_type2 + (1 | Patient_ID), data = data_df, family = binomial, nAGQ = 100) # Fit the mixed-effects logistic regression model
lrt_result <- anova(null_model, full_model) # Perform the Likelihood Ratio Test
lrt_result

# Non-remission
filtered_df_percent_NR <- filtered_df_percent[filtered_df_percent$Endoscopic_outcome=="0",]
response_variable <- filtered_df_percent_NR$cluster
response_variable <- ifelse(response_variable == "1", "0", "1")
response_variable <- as.numeric(response_variable)
confounding_variable <- filtered_df_percent_NR$Patient_ID
predictor_variable <- filtered_df_percent_NR$FMT_treatment_type2
data_df <- data.frame(cluster = response_variable, FMT_treatment_type2 = predictor_variable,Patient_ID = confounding_variable)

# Full model
null_model <- glmer(cluster ~ 1 + (1 | Patient_ID), data = data_df, family = binomial) # Fit the null model with only the intercept
full_model <- glmer(cluster ~ FMT_treatment_type2 + (1 | Patient_ID), data = data_df, family = binomial, nAGQ = 100) # Fit the mixed-effects logistic regression model
lrt_result <- anova(null_model, full_model) # Perform the Likelihood Ratio Test
lrt_result

# Both remission and non-remission are not significantly distinct, if controlled for biologically depencence of samples
# Thus use this as extended Data Figures but leave
###################################
# 3.5 Associate patients (-superdonor) with VCT divided by endoscopic outcome - Fig 2
###################################
# Figure 2: FMT Treatment Type 2 versus cluster divided by Treatment success for superdonors
df_percent <- phyloseq_rarefied_phages_3_stat %>%
  group_by(cluster) %>%
  mutate(Percentage = n() / sum(n()) * 100)

df_percent <- df_percent[df_percent$FMT_treatment_type2=="Patient (-superdonor)" | df_percent$FMT_treatment_type2=="Patient (-untreated)" | df_percent$FMT_treatment_type2=="Donor (-superdonor)",]
list_patient <- c("patient 1","patient 11","patient 12","patient 20","patient 41","patient 46","patient 50","patient 51","patient 52","patient 58","patient 67","patient 66","patient 68","patient 69","patient 72","patient 77","patient 78","patient 83")
filtered_df_percent <- df_percent[!(df_percent$Patient_ID %in% list_patient), ]

filtered_df_percent$Timepoints[filtered_df_percent$FMT_treatment_type2=="Donor (-superdonor)"] <- "Donor (-superdonor)"
filtered_df_percent$Timepoints <- factor(filtered_df_percent$Timepoints,levels = c("Donor (-superdonor)","Baseline","w4","w8","w12","m6","m12"))
filtered_df_percent$FMT_treatment_type2 <- factor(filtered_df_percent$FMT_treatment_type2,levels = c("Patient (-untreated)",'Patient (-superdonor)'))

# Figure 1: comparison VCT over time (only Donor (-superdonors))
ggplot(filtered_df_percent) +
  aes(x = Timepoints, fill = cluster) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = c(`1` = "#bd7f1b", `2` = "#2a8a92")) +
  theme_bw()

filtered_df_percent1 <- filtered_df_percent
response_variable <- filtered_df_percent1$cluster
response_variable <- ifelse(response_variable == "1", "0", "1")
response_variable <- as.numeric(response_variable)
confounding_variable <- filtered_df_percent1$Patient_ID
predictor_variable <- filtered_df_percent1$Timepoints
data_df <- data.frame(cluster = response_variable, Timepoints = predictor_variable,Patient_ID = confounding_variable)
#data_df <- data_df[data_df$Timepoints=="w4" | data_df$Timepoints=="w8" | data_df$Timepoints=="Donor (-superdonor)",]
#table(data_df$Timepoints)

# Full model
control <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))
null_model <- glmer(cluster ~ 1 + (1 | Patient_ID), data = data_df, family = binomial("logit"),control = control) # Fit the null model with only the intercept
full_model <- glmer(cluster ~ Timepoints + (1 | Patient_ID), data = data_df, family = binomial("logit"), nAGQ = 100) # Fit the mixed-effects logistic regression model
lrt_result <- anova(null_model, full_model) # Perform the Likelihood Ratio Test
lrt_result

# Comparison with reference variable (# compared to reference = Donor (-superdonor))
model <- glmer(cluster ~ Timepoints + (1 | Patient_ID), data = data_df, family = binomial("logit"),nAGQ = 100)
summary(model) 

# Conclusion:
# No difference in VCT over time when added to Donor (-superdonor)
###################################
# 4. Compositional analyses of clusters
###################################
# 4.1 Genus level drivers
###################################
# Abundances per clusters
ps1.rel <- microbiome::transform(phyloseq_rarefied_phages_3, "compositional")
ps1.fam.rel <-aggregate_rare(ps1.rel, level = "Genus", detection = 0.01, prevalence = 15/100) 

plot_composition(ps1.fam.rel,average_by = "cluster", otu.sort="abundance") + 
  guides(fill = guide_legend(ncol = 1)) +
  labs(x = "Samples", 
       y = "Relative abundance",
       title = "Relative abundance data", 
       subtitle = "Subtitle",
       caption = "Caption text.")  +
  scale_fill_brewer("Family", palette = "Paired") +
  theme_bw() 

# Extract the relative abundance data
abundance_data <- as.data.frame(t(otu_table(ps1.fam.rel)))
cluster_column <- sample_data(ps1.rel)$cluster
cluster_data <- data.frame(Cluster = sample_data(ps1.fam.rel)$cluster, row.names = rownames(sample_data(ps1.fam.rel)))

# Extract 'Patient_ID' data from the phyloseq object
patient_data <- sample_data(ps1.rel)$Patient_ID
patient_data <- data.frame(Patient_ID = patient_data, row.names = rownames(sample_data(ps1.rel)))

# Merge abundance_data, cluster_data, and patient_data
DF_1 <- merge(abundance_data, cluster_data, by = 0, all = FALSE)
rownames(DF_1) <- DF_1$Row.names
DF_1$Row.names <- NULL

# Merge with patient_data based on row names
DF_1 <- merge(DF_1, patient_data, by = 0, all = FALSE)
rownames(DF_1) <- DF_1$Row.names
DF_1$Row.names <- NULL

# Convert DF_1 to long format (melt)
DF_1_m <- melt(DF_1, id.vars = c("Cluster", "Patient_ID"))

# Figure Genus
DF_1_m %>%
  filter(!(variable %in% "Other")) %>%
  ggplot() +
  aes(x = variable, y = value, fill = Cluster) +
  geom_boxplot() +
  scale_fill_manual(values = c(`1` = "#bd7f1b", `2` = "#2a8a92")) +
  theme_bw()

# Statistics
# Convert 'Cluster' and 'Patient_ID' to factors
DF_1_m$Cluster <- as.factor(DF_1_m$Cluster)
DF_1_m$Patient_ID <- as.factor(DF_1_m$Patient_ID)

# List of genera to fit the linear mixed-effects model for
genera_of_interest <- c("Genus2", "Genus10", "Genus15", "Genus3", "Genus4", "Genus5", "Genus6")

# Initialize an empty data frame to store results
results_df <- data.frame(Genus = character(),
                         Estimate = numeric(),
                         p_value = numeric(),
                         stringsAsFactors = FALSE)

# Loop over each genus of interest
for (genus in genera_of_interest) {
  cat("Genus:", genus, "\n")
  
  # Subset the data for the current genus
  DF_1_m_sub <- subset(DF_1_m, variable == genus)
  
  # Convert 'variable' to a factor
  DF_1_m_sub$variable <- factor(DF_1_m_sub$variable, levels = genus)
  
  # Fit the null model with only intercept
  null_model <- lmer(value ~ 1 + (1 | Patient_ID), data = DF_1_m_sub)
  
  # Fit the full model with Cluster as a fixed effect
  full_model <- lmer(value ~ Cluster + (1 | Patient_ID), data = DF_1_m_sub,)
  
  # Perform the Likelihood Ratio Test
  lrt_result <- anova(null_model, full_model)
  
  # Store the p-value and estimate for the full model (Cluster effect) in the results data frame
  p_value <- lrt_result$`Pr(>Chisq)`[2]  # Extract the p-value for the full model
  chisq <- lrt_result$Chisq[2] # Extract chi-squared
  estimate <- summary(full_model)$coefficients[2, 1]  # Extract the estimate for the full model
  
  # Add results to the data frame
  results_df <- rbind(results_df, c(genus, estimate, chisq, p_value))
}

# Rename the columns
colnames(results_df) <- c("Genus", "Estimate","Chisq","p_value")

# Add a column for adjusted p-values using BH method
results_df$adjusted_p_value <- p.adjust(results_df$p_value, method = "BH")

# View the results data frame
print(results_df)

  # Conclusion:
# Genus5 has a significantly higher RA in VCT1 compared to VCT2
# Genus2, Genus15 and Genus 10 has a significantly higher RA in VCT2 compared to VCT1.
###################################
# 4.2 Order level drivers
###################################
# Assuming 'phyloseq_rarefied_phages_3' is already defined
phyloseq_rarefied_phages_VCT <- phyloseq_rarefied_phages_3
tax_table_df <- as.data.frame(tax_table(phyloseq_rarefied_phages_VCT))

# Set custom taxonomy "Viruses" for 'Species' and 'Phylum' taxa
tax_table_df$Species <- "Viruses"
tax_table_df$Phylum <- "Viruses"

# Convert the updated tax_table_df back to tax_table object
tax_table(phyloseq_rarefied_phages_VCT) <- as.matrix(tax_table_df)

# Aggregate at Order level and transform to compositional data
ps1.order.rel <- aggregate_rare(phyloseq_rarefied_phages_VCT, level = "Order", detection = 0.01, prevalence = 1/100)
ps1.rel <- microbiome::transform(ps1.order.rel, "compositional")

# Extract the relative abundance data
abundance_data <- as.data.frame(t(otu_table(ps1.rel)))
cluster_column <- sample_data(ps1.rel)$cluster
cluster_data <- data.frame(Cluster = sample_data(ps1.rel)$cluster, row.names = rownames(sample_data(ps1.rel)))

# Extract 'Patient_ID' data from the phyloseq object
patient_data <- sample_data(ps1.rel)$Patient_ID
patient_data <- data.frame(Patient_ID = patient_data, row.names = rownames(sample_data(ps1.rel)))

# Merge abundance_data, cluster_data, and patient_data
DF_1 <- merge(abundance_data, cluster_data, by = 0, all = FALSE)
rownames(DF_1) <- DF_1$Row.names
DF_1$Row.names <- NULL

# Merge with patient_data based on row names
DF_1 <- merge(DF_1, patient_data, by = 0, all = FALSE)
rownames(DF_1) <- DF_1$Row.names
DF_1$Row.names <- NULL

# Convert DF_1 to long format (melt)
DF_1_m <- melt(DF_1, id.vars = c("Cluster", "Patient_ID"))

# Figure Genus
Boxplot_Order1 <- DF_1_m %>%
 filter(variable %in% c("Crassvirales", "Unclassified Caudoviricetes", "Petitvirales")) %>%
 ggplot() +
 aes(x = variable, y = value, fill = Cluster) +
 geom_boxplot() +
 scale_fill_manual(values = c(`1` = "#bd7f1b", `2` = "#2a8a92")) +
 theme_bw() +
 theme(legend.position = "bottom") 

#ggsave("test.pdf",Boxplot_Order1, width = 8, height = 6)

# Statistics
# Convert 'Cluster' and 'Patient_ID' to factors
DF_1_m$Cluster <- as.factor(DF_1_m$Cluster)
DF_1_m$Patient_ID <- as.factor(DF_1_m$Patient_ID)

# List of genera to fit the linear mixed-effects model for
genera_of_interest <- c("Crassvirales", "Unclassified Caudoviricetes", "Petitvirales")

# Initialize an empty data frame to store results
results_df <- data.frame(Genus = character(),
                         Estimate = numeric(),
                         p_value = numeric(),
                         stringsAsFactors = FALSE)

# Loop over each genus of interest
for (genus in genera_of_interest) {
  cat("Genus:", genus, "\n")
  
  # Subset the data for the current genus (performing partial matching)
  DF_1_m_sub <- subset(DF_1_m, grepl(genus, variable))
  
  # Convert 'variable' to a factor
  DF_1_m_sub$variable <- factor(DF_1_m_sub$variable, levels = genus)
  
  # Fit the null model with only intercept
  null_model <- lmer(value ~ 1 + (1 | Patient_ID), data = DF_1_m_sub)
  
  # Fit the full model with Cluster as a fixed effect
  full_model <- lmer(value ~ Cluster + (1 | Patient_ID), data = DF_1_m_sub)
  
  # Perform the Likelihood Ratio Test
  lrt_result <- anova(null_model, full_model)
  
  # Store the p-value and estimate for the full model (Cluster effect) in the results data frame
  p_value <- lrt_result$`Pr(>Chisq)`[2]  # Extract the p-value for the full model
  chisq <- lrt_result$Chisq[2] # Extract chi-squared
  estimate <- summary(full_model)$coefficients[2, 1]  # Extract the estimate for the full model

  # Add results to the data frame
  results_df <- rbind(results_df, c(genus, estimate, chisq, p_value))
}

# Rename the columns
colnames(results_df) <- c("Order", "Estimate","Chisq","p_value")

# Add a column for adjusted p-values using BH method
results_df$p_value[results_df$Order=="Unclassified Caudoviricetes"] <- "0.03011"
results_df$adjusted_p_value <- p.adjust(results_df$p_value, method = "BH")

# View the results data frame
print(results_df)

# Conclusion:
# VCT1: Higher RA of Unclassified Caudoviricetes
# VCT2: Higher RA of Petitvirales
# Significant difference in RA Crassvirales
###################################
# 4.3 In Silico predicted host phyla
###################################
# Hypothesis:
# Since The majorly driving genera being Genus 15, Genus 10 and Genus 2 majorly infected Pseudomonadota
# We expected that VCT2 is driving an infection of this bacteria as well

table(Mastertable$Final_order[Mastertable$Genus=='Genus15'])
table(Mastertable$Phyla_phylum[Mastertable$Genus=='Genus5'])
table(Mastertable$Phyla_family[Mastertable$Genus=='Genus5'])
table(Mastertable$Phyla_genus[Mastertable$Genus=='Genus5'])
table(Mastertable$`lysogenic cycle`[Mastertable$Genus=='Genus5'])

table(Mastertable$Final_order[Mastertable$Genus=='Genus'])
table(Mastertable$Phyla_phylum[Mastertable$Genus=='Genus10'])
table(Mastertable$Phyla_genus[Mastertable$Genus=='Genus10'])
table(Mastertable$`lysogenic cycle`[Mastertable$Genus=='Genus5'])

table(Mastertable$Final_order[Mastertable$Genus=='Genus2'])
table(Mastertable$Phyla_phylum[Mastertable$Genus=='Genus2'])
table(Mastertable$Phyla_genus[Mastertable$Genus=='Genus2'])
table(Mastertable$`lysogenic cycle`[Mastertable$Genus=='Genus2'])

table(Mastertable$`lysogenic cycle`[Mastertable$Phyla_phylum=='Chlamydiota'])
table(Mastertable$`lysogenic cycle`[Mastertable$Phyla_phylum=='Bacillota'])
table(Mastertable$`lysogenic cycle`[Mastertable$Phyla_phylum=='Pseudomonadota'])

# Add this information in a table of In Silico prediction of host
# Extract columns "Patient_ID", "cluster", and "FMT_treatment_type2" and keep the rownames as "Sample_ID"
selected_data <- sample_data(phyloseq_rarefied_phages_3)[, c("Patient_ID", "cluster", "FMT_treatment_type2")]
rownames(selected_data) <- sample_data(phyloseq_rarefied_phages_3)$Sample_ID

# Extract Abundance table
otu_table_df <- as.data.frame(otu_table(phyloseq_rarefied_phages_3))

# Extract In silico host
Phyla_phylum_phage <- Mastertable_viral_rarefied_2[, "Phyla_phylum", drop = FALSE]

# Merge otu table with phyla
otu_table_df_phyla <- merge(otu_table_df,Phyla_phylum_phage,by=0,all=F)
rownames(otu_table_df_phyla) <- otu_table_df_phyla$Row.names
otu_table_df_phyla$Row.names <- NULL
otu_table_df_phyla$Phyla_phylum[otu_table_df_phyla$Phyla_phylum==0] <- "Unknown"

# Aggregate by host
names <- colnames(otu_table_df[-ncol(otu_table_df)])
otu_table_df_phyla <- aggregate(. ~Phyla_phylum, FUN = sum, data = otu_table_df_phyla[,colnames(otu_table_df_phyla) %in% names | colnames(otu_table_df_phyla) == 'Phyla_phylum'])
rownames(otu_table_df_phyla) <- otu_table_df_phyla$Phyla_phylum
otu_table_df_phyla$Phyla_phylum <- NULL

# Print most abundant host
otu_table_df_phyla_sum <- otu_table_df_phyla
otu_table_df_phyla_sum$total <- rowSums(otu_table_df_phyla_sum)
otu_table_df_phyla_sum <- otu_table_df_phyla_sum[rev(order(otu_table_df_phyla_sum$total)),]
otu_table_df_phyla_sum$percentage <- (otu_table_df_phyla_sum$total/sum(otu_table_df_phyla_sum$total))*100
selected_columns <- otu_table_df_phyla_sum[, c("total", "percentage")]
print(selected_columns)

# Calculate RA
otu_table_df_phyla_RA <- sweep(otu_table_df_phyla, 2, colSums(otu_table_df_phyla), '/')
otu_table_df_phyla[is.na(otu_table_df_phyla)] <- 0

# Merge with metadata
otu_meta_phyla <- merge(t(otu_table_df_phyla_RA), selected_data, by = 0, all = F)
rownames(otu_meta_phyla) <- otu_meta_phyla$Row.names
otu_meta_phyla$Row.names <- NULL

# Melt
otu_meta_phyla_m <- melt(otu_meta_phyla, id.vars = c("Patient_ID", "cluster", "FMT_treatment_type2"))
otu_meta_phyla_m$variable <- factor(otu_meta_phyla_m$variable, levels = c("Actinomycetota","Bacillota","Bacteroidota","Chlamydiota","Pseudomonadota","Unknown"))

# Figure
otu_meta_phyla_m %>%
 filter(variable %in% c("Bacillota", "Bacteroidota", "Pseudomonadota", "Unknown","Chlamydiota","Actinomycetota")) %>%
 ggplot() +
 aes(x = variable, y = value, fill = cluster) +
 geom_boxplot() +
 scale_fill_manual(values = c(`1` = "#bd7f1b", `2` = "#2a8a92")) +
 theme_bw()

#ggsave("test.pdf",Phyla1, width = 8, height = 4)

# Statistics
# Convert 'Cluster' and 'Patient_ID' to factors
otu_meta_phyla_m$Cluster <- as.factor(otu_meta_phyla_m$cluster)
otu_meta_phyla_m$Patient_ID <- as.factor(otu_meta_phyla_m$Patient_ID)

# List of host phyla to fit the linear mixed-effects model for
Host_phyla_to_extract <- c("Actinomycetota","Bacillota","Bacteroidota","Chlamydiota","Pseudomonadota","Unknown")

# Initialize an empty data frame to store results
results_df <- data.frame(Genus = character(),
                         Estimate = numeric(),
                         p_value = numeric(),
                         stringsAsFactors = FALSE)

# Loop over each genus of interest
for (genus in Host_phyla_to_extract) {
  cat("Genus:", genus, "\n")
  
  # Subset the data for the current genus (performing partial matching)
  DF_1_m_sub <- subset(otu_meta_phyla_m, grepl(genus, variable))
  
  # Convert 'variable' to a factor
  DF_1_m_sub$variable <- factor(DF_1_m_sub$variable, levels = genus)
  
  # Fit the null model with only intercept
  null_model <- lmer(value ~ 1 + (1 | Patient_ID), data = DF_1_m_sub)
  
  # Fit the full model with Cluster as a fixed effect
  full_model <- lmer(value ~ Cluster + (1 | Patient_ID), data = DF_1_m_sub)
  
  # Perform the Likelihood Ratio Test
  lrt_result <- anova(null_model, full_model)
  
  # Store the p-value and estimate for the full model (Cluster effect) in the results data frame
  p_value <- lrt_result$`Pr(>Chisq)`[2]  # Extract the p-value for the full model
  chisq <- lrt_result$Chisq[2] # Extract chi-squared
  estimate <- summary(full_model)$coefficients[2, 1]  # Extract the estimate for the full model

  # Add results to the data frame
  results_df <- rbind(results_df, c(genus, estimate, chisq, p_value))
}

# Rename the columns
colnames(results_df) <- c("Host phyla", "Estimate","Chisq","p_value")

# Add a column for adjusted p-values using BH method
results_df$adjusted_p_value <- p.adjust(results_df$p_value, method = "BH")

# View the results data frame
print(results_df)

# Conclusion:
# VCT2 higher in RA of Pseudomonadota-infecting phages and lower in Bacillota-infecting phages (Firmicutes)
# Host shown > 1% of reads
###################################
# 4.4 In Silico predicted host family
###################################
# Extract columns "Patient_ID", "cluster", and "FMT_treatment_type2" and keep the rownames as "Sample_ID"
selected_data <- sample_data(phyloseq_rarefied_phages_3)[, c("Patient_ID", "cluster", "FMT_treatment_type2")]
rownames(selected_data) <- sample_data(phyloseq_rarefied_phages_3)$Sample_ID

# Extract Abundance table
otu_table_df <- as.data.frame(otu_table(phyloseq_rarefied_phages_3))

# Extract In silico host
Phyla_phylum_phage <- Mastertable_viral_rarefied_2[, "Phyla_family", drop = FALSE]

# Merge otu table with phyla
otu_table_df_phyla <- merge(otu_table_df,Phyla_phylum_phage,by=0,all=F)
rownames(otu_table_df_phyla) <- otu_table_df_phyla$Row.names
otu_table_df_phyla$Row.names <- NULL
otu_table_df_phyla$Phyla_family[otu_table_df_phyla$Phyla_family==0] <- "Unknown"
table(otu_table_df_phyla$Phyla_family)

# Aggregate by host
names <- colnames(otu_table_df[-ncol(otu_table_df)])
otu_table_df_phyla <- aggregate(. ~Phyla_family, FUN = sum, data = otu_table_df_phyla[,colnames(otu_table_df_phyla) %in% names | colnames(otu_table_df_phyla) == 'Phyla_family'])
rownames(otu_table_df_phyla) <- otu_table_df_phyla$Phyla_family
otu_table_df_phyla$Phyla_family <- NULL

# Print most abundant host
otu_table_df_phyla_sum <- otu_table_df_phyla
otu_table_df_phyla_sum$total <- rowSums(otu_table_df_phyla_sum)
otu_table_df_phyla_sum <- otu_table_df_phyla_sum[rev(order(otu_table_df_phyla_sum$total)),]
otu_table_df_phyla_sum$percentage <- (otu_table_df_phyla_sum$total/sum(otu_table_df_phyla_sum$total))*100
selected_columns <- otu_table_df_phyla_sum[, c("total", "percentage")]
print(selected_columns)

# Calculate RA
otu_table_df_phyla_RA <- sweep(otu_table_df_phyla, 2, colSums(otu_table_df_phyla), '/')
otu_table_df_phyla[is.na(otu_table_df_phyla)] <- 0

# Merge with metadata
otu_meta_phyla <- merge(t(otu_table_df_phyla_RA), selected_data, by = 0, all = F)
rownames(otu_meta_phyla) <- otu_meta_phyla$Row.names
otu_meta_phyla$Row.names <- NULL

# Melt
otu_meta_phyla_m <- melt(otu_meta_phyla, id.vars = c("Patient_ID", "cluster", "FMT_treatment_type2"))
otu_meta_phyla_m$variable <- factor(otu_meta_phyla_m$variable, levels = c("Bacteroidaceae","Chlamydiaceae","Enterobacteriaceae","Lachnospiraceae","Oscillospiraceae","Prevotellaceae","Streptococcaceae","Tannerellaceae","Unknown","Veillonellaceae"))

# Figure (> 1%)
otu_meta_phyla_m %>%
  filter(variable %in% c("Unknown", "Enterobacteriaceae", "Prevotellaceae", "Chlamydiaceae","Bacteroidaceae","Lachnospiraceae","Oscillospiraceae","Streptococcaceae","Tannerellaceae","Veillonellaceae")) %>%
  ggplot() +
  aes(x = variable, y = value, fill = cluster) +
  geom_boxplot() +
  scale_fill_manual(values = c(`1` = "#bd7f1b", `2` = "#2a8a92")) +
  theme_bw()

# Statistics
# Convert 'Cluster' and 'Patient_ID' to factors
otu_meta_phyla_m$Cluster <- as.factor(otu_meta_phyla_m$cluster)
otu_meta_phyla_m$Patient_ID <- as.factor(otu_meta_phyla_m$Patient_ID)

# List of host phyla to fit the linear mixed-effects model for
Host_phyla_to_extract <- c("Bacteroidaceae","Chlamydiaceae","Enterobacteriaceae","Lachnospiraceae","Oscillospiraceae","Prevotellaceae","Streptococcaceae","Tannerellaceae","Unknown","Veillonellaceae")

# Initialize an empty data frame to store results
results_df <- data.frame(Genus = character(),
                         Estimate = numeric(),
                         p_value = numeric(),
                         stringsAsFactors = FALSE)

# Loop over each genus of interest
for (genus in Host_phyla_to_extract) {
  cat("Genus:", genus, "\n")
  
  # Subset the data for the current genus (performing partial matching)
  DF_1_m_sub <- subset(otu_meta_phyla_m, grepl(genus, variable))
  
  # Convert 'variable' to a factor
  DF_1_m_sub$variable <- factor(DF_1_m_sub$variable, levels = genus)
  
  # Fit the null model with only intercept
  null_model <- lmer(value ~ 1 + (1 | Patient_ID), data = DF_1_m_sub)
  
  # Fit the full model with Cluster as a fixed effect
  full_model <- lmer(value ~ Cluster + (1 | Patient_ID), data = DF_1_m_sub)
  
  # Perform the Likelihood Ratio Test
  lrt_result <- anova(null_model, full_model)
  
  # Store the p-value and estimate for the full model (Cluster effect) in the results data frame
  p_value <- lrt_result$`Pr(>Chisq)`[2]  # Extract the p-value for the full model
  chisq <- lrt_result$Chisq[2] # Extract chi-squared
  estimate <- summary(full_model)$coefficients[2, 1]  # Extract the estimate for the full model

  # Add results to the data frame
  results_df <- rbind(results_df, c(genus, estimate, chisq, p_value))
}

# Rename the columns
colnames(results_df) <- c("Host family", "Estimate","Chisq","p_value")
# Add a column for adjusted p-values using BH method
results_df$adjusted_p_value <- p.adjust(results_df$p_value, method = "BH")

# View the results data frame
print(results_df)

# Conclusion:
# VCT2 higher in RA of Enterobactericeae-infecting phages compared to VC1.
###################################
# 4.5 In Silico predicted host genus
###################################
# Extract columns "Patient_ID", "cluster", and "FMT_treatment_type2" and keep the rownames as "Sample_ID"
selected_data <- sample_data(phyloseq_rarefied_phages_3)[, c("Patient_ID", "cluster", "FMT_treatment_type2")]
rownames(selected_data) <- sample_data(phyloseq_rarefied_phages_3)$Sample_ID

# Extract Abundance table
otu_table_df <- as.data.frame(otu_table(phyloseq_rarefied_phages_3))

# Extract In silico host
Phyla_phylum_phage <- Mastertable_viral_rarefied_2[, "Phyla_genus", drop = FALSE]

# Merge otu table with phyla
otu_table_df_phyla <- merge(otu_table_df,Phyla_phylum_phage,by=0,all=F)
rownames(otu_table_df_phyla) <- otu_table_df_phyla$Row.names
otu_table_df_phyla$Row.names <- NULL
otu_table_df_phyla$Phyla_genus[otu_table_df_phyla$Phyla_genus==0] <- "Unknown"
table(otu_table_df_phyla$Phyla_genus)

# Aggregate by host
names <- colnames(otu_table_df[-ncol(otu_table_df)])
otu_table_df_phyla <- aggregate(. ~Phyla_genus, FUN = sum, data = otu_table_df_phyla[,colnames(otu_table_df_phyla) %in% names | colnames(otu_table_df_phyla) == 'Phyla_genus'])
rownames(otu_table_df_phyla) <- otu_table_df_phyla$Phyla_genus
otu_table_df_phyla$Phyla_genus <- NULL

# Print most abundant host
otu_table_df_phyla_sum <- otu_table_df_phyla
otu_table_df_phyla_sum$total <- rowSums(otu_table_df_phyla_sum)
otu_table_df_phyla_sum <- otu_table_df_phyla_sum[rev(order(otu_table_df_phyla_sum$total)),]
otu_table_df_phyla_sum$percentage <- (otu_table_df_phyla_sum$total/sum(otu_table_df_phyla_sum$total))*100
selected_columns <- otu_table_df_phyla_sum[, c("total", "percentage")]
print(selected_columns)

# Calculate RA
otu_table_df_phyla_RA <- sweep(otu_table_df_phyla, 2, colSums(otu_table_df_phyla), '/')
otu_table_df_phyla[is.na(otu_table_df_phyla)] <- 0

# Merge with metadata
otu_meta_phyla <- merge(t(otu_table_df_phyla_RA), selected_data, by = 0, all = F)
rownames(otu_meta_phyla) <- otu_meta_phyla$Row.names
otu_meta_phyla$Row.names <- NULL

# Melt
otu_meta_phyla_m <- melt(otu_meta_phyla, id.vars = c("Patient_ID", "cluster", "FMT_treatment_type2"))
otu_meta_phyla_m$variable <- factor(otu_meta_phyla_m$variable, levels = c("Bacteroides","Blautia","Chlamydia","Escherichia","Faecalibacterium","Lactococcus","Parabacteroides","Prevotella","Unknown","Veillonella"))

# Figure (> 1%)
otu_meta_phyla_m %>%
  filter(variable %in% c("Bacteroides","Blautia","Chlamydia","Escherichia","Faecalibacterium","Lactococcus","Parabacteroides","Prevotella","Unknown","Veillonella")) %>%
  ggplot() +
  aes(x = variable, y = value, fill = cluster) +
  geom_boxplot() +
  scale_fill_manual(values = c(`1` = "#bd7f1b", `2` = "#2a8a92")) +
  theme_bw()

# Statistics
# Convert 'Cluster' and 'Patient_ID' to factors
otu_meta_phyla_m$Cluster <- as.factor(otu_meta_phyla_m$cluster)
otu_meta_phyla_m$Patient_ID <- as.factor(otu_meta_phyla_m$Patient_ID)

# List of host phyla to fit the linear mixed-effects model for
Host_phyla_to_extract <- c("Bacteroides","Blautia","Chlamydia","Escherichia","Faecalibacterium","Lactococcus","Parabacteroides","Prevotella","Unknown","Veillonella")

# Initialize an empty data frame to store results
results_df <- data.frame(Genus = character(),
                         Estimate = numeric(),
                         p_value = numeric(),
                         stringsAsFactors = FALSE)

# Loop over each genus of interest
for (genus in Host_phyla_to_extract) {
  cat("Genus:", genus, "\n")
  
  # Subset the data for the current genus (performing partial matching)
  DF_1_m_sub <- subset(otu_meta_phyla_m, grepl(genus, variable))
  
  # Convert 'variable' to a factor
  DF_1_m_sub$variable <- factor(DF_1_m_sub$variable, levels = genus)
  
  # Fit the null model with only intercept
  null_model <- lmer(value ~ 1 + (1 | Patient_ID), data = DF_1_m_sub)
  
  # Fit the full model with Cluster as a fixed effect
  full_model <- lmer(value ~ Cluster + (1 | Patient_ID), data = DF_1_m_sub)
  
  # Perform the Likelihood Ratio Test
  lrt_result <- anova(null_model, full_model)
  
  # Store the p-value and estimate for the full model (Cluster effect) in the results data frame
  p_value <- lrt_result$`Pr(>Chisq)`[2]  # Extract the p-value for the full model
  chisq <- lrt_result$Chisq[2] # Extract chi-squared
  estimate <- summary(full_model)$coefficients[2, 1]  # Extract the estimate for the full model

  # Add results to the data frame
  results_df <- rbind(results_df, c(genus, estimate, chisq, p_value))
}

# Rename the columns
colnames(results_df) <- c("Host genus", "Estimate","Chisq","p_value")

# Add a column for adjusted p-values using BH method
results_df$adjusted_p_value <- p.adjust(results_df$p_value, method = "BH")
print(results_df)

# Conclusion:
# VCT2 higher in RA of Escherichia-infecting phages compared to VC1.
# This is quite clear.
# Mostly unknowns though.
###################################
# 5. Viral lifestyle
###################################
# 5.1 All phages
###################################
# Extract columns "Patient_ID", "cluster", and "FMT_treatment_type2" and keep the rownames as "Sample_ID"
selected_data <- sample_data(phyloseq_rarefied_phages_3)[, c("Patient_ID", "cluster", "FMT_treatment_type2")]
rownames(selected_data) <- sample_data(phyloseq_rarefied_phages_3)$Sample_ID

# Extract Abundance table
otu_table_df <- as.data.frame(otu_table(phyloseq_rarefied_phages_3))

# Extract In silico host
Phyla_phylum_phage <- Mastertable_viral_rarefied_2[, "lysogenic cycle", drop = FALSE]
table(Phyla_phylum_phage$`lysogenic cycle`)

# Merge otu table with phyla
otu_table_df_phyla <- merge(otu_table_df,Phyla_phylum_phage,by=0,all=F)
rownames(otu_table_df_phyla) <- otu_table_df_phyla$Row.names
otu_table_df_phyla$Row.names <- NULL
otu_table_df_phyla$`lysogenic cycle`[otu_table_df_phyla$`lysogenic cycle`==0] <- "Lytic"
otu_table_df_phyla$`lysogenic cycle`[otu_table_df_phyla$`lysogenic cycle`=="temperate"] <- "Lysogenic"
table(otu_table_df_phyla$`lysogenic cycle`)

# Aggregate by host
names <- colnames(otu_table_df[-ncol(otu_table_df)])
otu_table_df_phyla <- aggregate(. ~`lysogenic cycle`, FUN = sum, data = otu_table_df_phyla[,colnames(otu_table_df_phyla) %in% names | colnames(otu_table_df_phyla) == 'lysogenic cycle'])
rownames(otu_table_df_phyla) <- otu_table_df_phyla$`lysogenic cycle`
otu_table_df_phyla$`lysogenic cycle` <- NULL

# Calculate RA
otu_table_df_phyla_RA <- sweep(otu_table_df_phyla, 2, colSums(otu_table_df_phyla), '/')
otu_table_df_phyla_RA[is.na(otu_table_df_phyla_RA)] <- 0

# Merge with metadata
otu_meta_phyla <- merge(t(otu_table_df_phyla_RA), selected_data, by = 0, all = F)
rownames(otu_meta_phyla) <- otu_meta_phyla$Row.names
otu_meta_phyla$Row.names <- NULL

# Melt
otu_meta_phyla_m <- melt(otu_meta_phyla, id.vars = c("Patient_ID", "cluster", "FMT_treatment_type2"))

# Figure 
otu_meta_phyla_m %>%
 filter(variable %in% "Lysogenic") %>%
 ggplot() +
 aes(x = variable, y = value, fill = cluster) +
 geom_boxplot() +
 scale_fill_manual(values = c(`1` = "#bd7f1b", `2` = "#2a8a92")) +
 theme_bw() +
 theme(legend.position = "none")

# Statistics
# Convert 'Cluster' and 'Patient_ID' to factors
otu_meta_phyla_m$Cluster <- as.factor(otu_meta_phyla_m$cluster)
otu_meta_phyla_m$Patient_ID <- as.factor(otu_meta_phyla_m$Patient_ID)

# List of host phyla to fit the linear mixed-effects model for
Host_phyla_to_extract <- c("Lysogenic")

# Initialize an empty data frame to store results
results_df <- data.frame(Genus = character(),
                         Estimate = numeric(),
                         p_value = numeric(),
                         stringsAsFactors = FALSE)

# Loop over each genus of interest
for (genus in Host_phyla_to_extract) {
  cat("Genus:", genus, "\n")
  
  # Subset the data for the current genus (performing partial matching)
  DF_1_m_sub <- subset(otu_meta_phyla_m, grepl(genus, variable))
  
  # Convert 'variable' to a factor
  DF_1_m_sub$variable <- factor(DF_1_m_sub$variable, levels = genus)
  
  # Fit the null model with only intercept
  null_model <- lmer(value ~ 1 + (1 | Patient_ID), data = DF_1_m_sub)
  
  # Fit the full model with Cluster as a fixed effect
  full_model <- lmer(value ~ Cluster + (1 | Patient_ID), data = DF_1_m_sub)
  
  # Perform the Likelihood Ratio Test
  lrt_result <- anova(null_model, full_model)
  
  # Store the p-value and estimate for the full model (Cluster effect) in the results data frame
  p_value <- lrt_result$`Pr(>Chisq)`[2]  # Extract the p-value for the full model
  chisq <- lrt_result$Chisq[2] # Extract chi-squared
  estimate <- summary(full_model)$coefficients[2, 1]  # Extract the estimate for the full model
  
  # Add results to the data frame
  results_df <- rbind(results_df, c(genus, estimate, chisq, p_value))
}

# Rename the columns
colnames(results_df) <- c("Lifestyle", "Estimate","Chisq","p_value")

# Add a column for adjusted p-values using BH method
results_df$adjusted_p_value <- p.adjust(results_df$p_value, method = "BH")

# View the results data frame
print(results_df)

# Conclusion:
# VCT2 lower lysogenic potential of phages or more lytic phages
###################################
# 5.2 Unclassified Caudoviricetes
###################################
# Extract columns "Patient_ID", "cluster", and "FMT_treatment_type2" and keep the rownames as "Sample_ID"
selected_data <- sample_data(phyloseq_rarefied_phages_3)[, c("Patient_ID", "cluster", "FMT_treatment_type2")]
rownames(selected_data) <- sample_data(phyloseq_rarefied_phages_3)$Sample_ID

# Extract Abundance table
otu_table_df <- as.data.frame(otu_table(phyloseq_rarefied_phages_3))

# Extract In silico host
Mastertable_viral_rarefied_2_CA <- Mastertable_viral_rarefied_2[Mastertable_viral_rarefied_2$Final_order=="Unclassified Caudoviricetes",]
Phyla_phylum_phage <- Mastertable_viral_rarefied_2_CA[, "lysogenic cycle", drop = FALSE]
table(Mastertable_viral_rarefied_2_CA$`lysogenic cycle`)

# Merge otu table with phyla
otu_table_df_phyla <- merge(otu_table_df,Phyla_phylum_phage,by=0,all=F)
rownames(otu_table_df_phyla) <- otu_table_df_phyla$Row.names
otu_table_df_phyla$Row.names <- NULL
otu_table_df_phyla$`lysogenic cycle`[otu_table_df_phyla$`lysogenic cycle`==0] <- "Lytic"
otu_table_df_phyla$`lysogenic cycle`[otu_table_df_phyla$`lysogenic cycle`=="temperate"] <- "Lysogenic"
table(otu_table_df_phyla$`lysogenic cycle`)

# Aggregate by host
names <- colnames(otu_table_df[-ncol(otu_table_df)])
otu_table_df_phyla <- aggregate(. ~`lysogenic cycle`, FUN = sum, data = otu_table_df_phyla[,colnames(otu_table_df_phyla) %in% names | colnames(otu_table_df_phyla) == 'lysogenic cycle'])
rownames(otu_table_df_phyla) <- otu_table_df_phyla$`lysogenic cycle`
otu_table_df_phyla$`lysogenic cycle` <- NULL

# Calculate RA
otu_table_df_phyla_RA <- sweep(otu_table_df_phyla, 2, colSums(otu_table_df_phyla), '/')
otu_table_df_phyla_RA[is.na(otu_table_df_phyla_RA)] <- 0

# Merge with metadata
otu_meta_phyla <- merge(t(otu_table_df_phyla_RA), selected_data, by = 0, all = F)
rownames(otu_meta_phyla) <- otu_meta_phyla$Row.names
otu_meta_phyla$Row.names <- NULL

# Melt
otu_meta_phyla_m <- melt(otu_meta_phyla, id.vars = c("Patient_ID", "cluster", "FMT_treatment_type2"))

# Figure 
otu_meta_phyla_m %>%
  filter(variable %in% "Lysogenic") %>%
  ggplot() +
  aes(x = variable, y = value, fill = cluster) +
  geom_boxplot() +
  scale_fill_manual(values = c(`1` = "#bd7f1b", `2` = "#2a8a92")) +
  theme_bw() +
  theme(legend.position = "none")

# Statistics
# Convert 'Cluster' and 'Patient_ID' to factors
otu_meta_phyla_m$Cluster <- as.factor(otu_meta_phyla_m$cluster)
otu_meta_phyla_m$Patient_ID <- as.factor(otu_meta_phyla_m$Patient_ID)

# List of host phyla to fit the linear mixed-effects model for
Host_phyla_to_extract <- c("Lysogenic")

# Initialize an empty data frame to store results
results_df <- data.frame(Genus = character(),
                         Estimate = numeric(),
                         p_value = numeric(),
                         stringsAsFactors = FALSE)

# Loop over each genus of interest
for (genus in Host_phyla_to_extract) {
  cat("Genus:", genus, "\n")
  
  # Subset the data for the current genus (performing partial matching)
  DF_1_m_sub <- subset(otu_meta_phyla_m, grepl(genus, variable))
  
  # Convert 'variable' to a factor
  DF_1_m_sub$variable <- factor(DF_1_m_sub$variable, levels = genus)
  
  # Fit the null model with only intercept
  null_model <- lmer(value ~ 1 + (1 | Patient_ID), data = DF_1_m_sub)
  
  # Fit the full model with Cluster as a fixed effect
  full_model <- lmer(value ~ Cluster + (1 | Patient_ID), data = DF_1_m_sub)
  
  # Perform the Likelihood Ratio Test
  lrt_result <- anova(null_model, full_model)
  
  # Store the p-value and estimate for the full model (Cluster  effect) in the results data frame
  p_value <- lrt_result$`Pr(>Chisq)`[2]  # Extract the p-value for the full model
  chisq <- lrt_result$Chisq[2] # Extract chi-squared
  estimate <- summary(full_model)$coefficients[2, 1]  # Extract the estimate for the full model
  
  # Add results to the data frame
  results_df <- rbind(results_df, c(genus, estimate, chisq, p_value))
}

# Rename the columns
colnames(results_df) <- c("Lifestyle Caudoviricetes", "Estimate","Chisq","p_value")

# Add a column for adjusted p-values using BH method
results_df$adjusted_p_value <- p.adjust(results_df$p_value, method = "BH")

# View the results data frame
print(results_df)

# Conclusion:
# VCT2 lower Caudoviricetes lysogenic potential of phages or more lytic Caudoviricetes phages
###################################


###################################
# Save samples and cluster
selected_data <- as.data.frame(sample_data(phyloseq_rarefied_phages_3)[, c("cluster")])

# Extract the 'cluster' column
cluster_column <- selected_data$cluster

# Create a new data frame with row names preserved
cluster_data_frame <- data.frame(cluster = cluster_column)
rownames(cluster_data_frame) <- rownames(selected_data)

# Save VCT
setwd ("/Users/daan/Desktop/Transfer/") 
write.table(cluster_data_frame, file = "./cluster_data_frame.txt", sep = "\t",row.names = FALSE)



###################################