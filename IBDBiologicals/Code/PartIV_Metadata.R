####################################
# SCRIPT 4: METADATA
####################################
# Before starting this script
####################################
# ------------------> We will load metadata in this script and assess it in the IBD cohort
# ------------------> This is a very individual and not automatic script
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
library(devtools) # before downlading from github you need to load devtools
#install.packages("BiocManager")
library(BiocManager)
#BiocManager::install("ComplexHeatmap")
#install_github("jokergoo/ComplexHeatmap")
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
####################################
# 1. Metadata (individual)
####################################
# Take and insert your own PATH to your Metadata. Make sure it is save in a comma-delimited file.
####################################
# 1.1 Insert data
####################################
setwd("/Users/daan/Desktop/Bioinformatics/pre-processing/Biological_study/Metadata/output_metadata")
getwd()
dir()
Metadata <- read.csv("Final_Metadata_subselection1_Relevant_R_without-units.csv", sep = ";")
rownames(Metadata) <- Metadata$Randomized.Nr.
View(Metadata)
# Females are 0 and males are 1
# Moisture is in %
# Bacterial cell count is in counts
# Disease duration is in years
# Smoking: 0 is non-roker & 1 is roker
# Smoking_2: 0 is non-roker, 1 is ex-roker, 2 is roker
####################################
# 1.2 Navigate through it using new package 'esquisse' - manual here of how to download
####################################
# 1.2.1 Install Homebrew
####################################
# Go to mac terminal & install homebrew. This you can use to download a bunch a programs on your mac.
# /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
####################################
# 1.2.2 Install Git & setup
####################################
# 1.2.2.1 Install Git
####################################
# Go the terminal in mac and download Git
# command: "brew install git"
# command: "which git" & put in "R/Tools/Global_options/git & paste it in the first path.
# see in : https://www.geo.uzh.ch/microsite/reproducible_research/post/rr-rstudio-git/
####################################
# 1.2.2.2 Configure Github
####################################
# username: DaanJansen94
# account: daan.jansen@kuleuven.be
# Paswoord: Heroesneverdie1

# In mac terminal do the following command
# git config --global user.name 'DaanJansen94'
# git config --global user.email 'daan.jansen@kuleuven.be'
# git config --global credential.helper osxkeychain
####################################
# 1.3.3 Open Git Project
####################################
# Go to "File/Git/"
# File in the repository we want to load to have all those packages in our "3. Rarefaction folder"
# Navigate here to use the package if we would liek that.
####################################
# 1.4 Time series
####################################
# 1.4.1 Install and load libraries
####################################
# check: TIPS_16 in "R/Rarefaction/TIPS/TIPS_16"
# in youtube: https://www.youtube.com/watch?v=Nf8FwFCJz2c&list=PLo32uKohmrXtCExTRr-e8SxyQK6yT_Xri&index=17
install.packages("timetk")
#library(plotly)
#library(tidyverse)
#library(timetk)
####################################
# 1.4.2 Quick visualization
####################################
# Make sure we have our dates in the USA style of writting. You can adopt this the easiest in excel before inputting the metadata.
####################################
# 1.4.2.1 Create tibble DF
####################################
Metadata <- tibble(Metadata)
Metadata$Sample_date <- as.Date(Metadata$Sample_date)
Metadata$Date_diagnosis <- as.Date(Metadata$Date_diagnosis)
Metadata$Date_of_birth <- as.Date(Metadata$Date_of_birth)
####################################
# 1.4.2.2 Plot time-series for sample collection of dataset
####################################
## In these figures you can easily see how long and when every single sample was collected from this dataset.
## This is divided for both Crohn's disease and Ulcerative Colitis

# ABSOLUTE
Metadata_sample_time_series <- Metadata[order(as.Date(Metadata$Sample_date, format="%Y/%m/%d")),]
Metadata_sample_time_series$cum_value <- 1
Metadata_sample_time_series$cum_value[Metadata_sample_time_series$Diagnosis == "CD"] <- seq.int(nrow(Metadata_sample_time_series[Metadata_sample_time_series$Diagnosis == "CD",]))
Metadata_sample_time_series$cum_value[Metadata_sample_time_series$Diagnosis == "UC"] <- seq.int(nrow(Metadata_sample_time_series[Metadata_sample_time_series$Diagnosis == "UC",]))

# RELATIVE
Metadata_sample_time_series$Percentages_cum_value <- 1
nrow(Metadata_sample_time_series[Metadata_sample_time_series$Diagnosis == "CD",]) # 312 samples of patients with CD
Metadata_sample_time_series$Percentages_cum_value[Metadata_sample_time_series$Diagnosis == "CD"] <- (Metadata_sample_time_series$cum_value[Metadata_sample_time_series$Diagnosis == "CD"]/312)*100
nrow(Metadata_sample_time_series[Metadata_sample_time_series$Diagnosis == "UC",]) # 120 samples of patients with UC
Metadata_sample_time_series$Percentages_cum_value[Metadata_sample_time_series$Diagnosis == "UC"] <- (Metadata_sample_time_series$cum_value[Metadata_sample_time_series$Diagnosis == "UC"]/120)*100
# View(Metadata_sample_time_series)

# Sample Date Absolute (put interactive to false if you want to adapt figure later)
Metadata_sample_time_series_p_absolute <- Metadata_sample_time_series %>%
  group_by(Diagnosis) %>%
  plot_time_series(
    .date_var    = Sample_date,
    .value       = cum_value,
    .smooth      = TRUE,
    .facet_ncol  = 3,
    .interactive = FALSE
  )

Metadata_sample_time_series_p_absolute

# Sample Date relative
Metadata_sample_time_series_p_relative <- Metadata_sample_time_series %>%
  group_by(Diagnosis) %>%
  plot_time_series(
    .date_var    = Sample_date,
    .value       = Percentages_cum_value,
    .smooth      = TRUE,
    .facet_ncol  = 3,
    .interactive = FALSE
  )

Metadata_sample_time_series_p_relative

# PS: use Illustrator to put a marked line on when 50% of data was collected ,as you can see this is < 2 y ago for both disease.
# PS1: you could use SPSS for simple and fast creation of simple barplots in which you can easily see the number of samples into w0, w14 and w24.
# PS2: Same for pie-chart of treatments used for every patient. Easiest in SPSS.
# PS3: Also a distribution of age & BMI is something you could easily put in a barplot and show.
# PS4 you could also use these time-series to show numeric variables in y-axis (in video) & eg. plot CRP, BMI or something else.
####################################
# 1.4.2.3 Plot time-series from scratch in ggplot
####################################
## ABSOLUTE
Metadata_sample_time_series %>%
  ggplot(aes(Sample_date, cum_value, color = Diagnosis, group = Diagnosis)) +
  geom_line() +
  geom_smooth(
    aes(Sample_date, cum_value),
    inherit.aes = FALSE,
    se = FALSE
  ) +
  facet_wrap(
    facets = ~ Diagnosis,
    scales = "free_y",
    ncol   = 2
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")
####################################
# 1.4.2.4 Cyberpunk style
####################################
## RELATIVE

# Setup a Color Palette and show what they look like
clrs <- colorRampPalette(c("#00ff9f", "#00b8ff", "#001eff", "#bd00ff", "#d600ff"))(7)
clrs2 <- colorRampPalette(c("#9f9fdf", "#FFB099"))(2)

scales::show_col(clrs)

# Setup colors for background, gridlines and text
clr_bg   <- "black"
clr_bg2  <- "gray10"
clr_grid <- "gray30"
clr_text <- "#d600ff"

# Setup a ggplot theme
# This is very important. This will create they frame of your ggplot schemes that we can use later.
# This is somethih we could add to any ggplot function
theme_cyberpunk <- function() {
  theme(
    # Plot / Panel
    plot.background = element_rect(fill = clr_bg, colour = clr_bg),
    # plot.margin = margin(1.5, 2, 1.5, 1.5, "cm"),
    panel.background = element_rect(fill = clr_bg, color = clr_bg),
    # Grid
    panel.grid = element_line(colour = clr_grid, size = 1),
    panel.grid.major = element_line(colour = clr_grid, size = 1),
    panel.grid.minor = element_line(colour = clr_grid, size = 1),
    axis.ticks.x = element_line(colour = clr_grid, size = 1),
    axis.line.y = element_line(colour = clr_grid, size = 0.5),
    axis.line.x = element_line(colour = clr_grid, size = 0.5),
    # Text
    plot.title = element_text(colour = clr_text),
    plot.subtitle = element_text(colour = clr_text),
    axis.text = element_text(colour = clr_text),
    axis.title = element_text(colour = clr_text),
    # Legend
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.title = element_text(colour = clr_text),
    legend.text = element_text(colour = "gray80", size = 12, face = "bold"),
    # Strip
    strip.background = element_rect(fill = clr_bg2, color = clr_bg2)
  )
}

# Make Custom ggplot
# We will use our previous quick visualization as input now

Metadata_sample_time_series_p_relative_esthetics_improved_creepy <- Metadata_sample_time_series_p_relative +
  geom_area(position = "identity", alpha = 0.) +
  geom_line(aes(color = Diagnosis), size = 1) +
  scale_color_manual(values = clrs) +
  geom_smooth(
    aes(Sample_date, Percentages_cum_value),
    inherit.aes = FALSE,
    se = FALSE
  ) +
  # facet_wrap(facets = NULL) +
  labs(title = "Sample collection IBD patients", subtitle = "Time Series Plot") +
  theme_cyberpunk()

Metadata_sample_time_series_p_relative_esthetics_improved <- Metadata_sample_time_series_p_relative +
  geom_area(position = "identity", alpha = 0.15) +
  geom_line(aes(color = Diagnosis), size = 1) +
  scale_color_manual(values = clrs2) +
  #geom_smooth(
  # aes(Sample_date, Percentages_cum_value),
  #  inherit.aes = FALSE,
  #  se = FALSE
  #  ) +
  # facet_wrap(facets = NULL) +
  labs(title = "Sample collection IBD patients", subtitle = "Time Series Plot") +
  coord_cartesian(ylim=c(3,100)) +
  theme()

Metadata_sample_time_series_p_relative_esthetics_improved_bw <- Metadata_sample_time_series_p_relative +
  geom_area(position = "identity", alpha = 0.15) +
  geom_line(aes(color = Diagnosis), size = 1) +
  scale_color_manual(values = clrs2) +
  #geom_smooth(
  # aes(Sample_date, Percentages_cum_value),
  #  inherit.aes = FALSE,
  #  se = FALSE
  #  ) +
  # facet_wrap(facets = NULL) +
  labs(title = "Sample collection IBD patients", subtitle = "Time Series Plot") +
  coord_cartesian(ylim=c(3,100)) +
  theme_bw() 

Metadata_sample_time_series_p_relative
Metadata_sample_time_series_p_relative_esthetics_improved
Metadata_sample_time_series_p_relative_esthetics_improved_bw
Metadata_sample_time_series_p_relative_esthetics_improved_creepy

# * Make Interactive plotly
ggplotly(g_cyberpunk)

## ABSOLUTE
Metadata_sample_time_series_p_absolute_esthetics_improve <- Metadata_sample_time_series_p_absolute +
  geom_area(position = "identity", alpha = 0.15) +
  geom_line(aes(color = Diagnosis), size = 1) +
  scale_color_manual(values = clrs2) +
  #geom_smooth(
  # aes(Sample_date, Percentages_cum_value),
  #  inherit.aes = FALSE,
  #  se = FALSE
  #  ) +
  # facet_wrap(facets = NULL) +
  labs(title = "Sample collection IBD patients", subtitle = "Time Series Plot") +
  coord_cartesian(ylim=c(10,310)) +
  theme()

Metadata_sample_time_series_p_absolute_esthetics_improve_bw <- Metadata_sample_time_series_p_absolute +
  geom_area(position = "identity", alpha = 0.15) +
  geom_line(aes(color = Diagnosis), size = 1) +
  scale_color_manual(values = clrs2) +
  #geom_smooth(
  # aes(Sample_date, Percentages_cum_value),
  #  inherit.aes = FALSE,
  #  se = FALSE
  #  ) +
  # facet_wrap(facets = NULL) +
  labs(title = "Sample collection IBD patients", subtitle = "Time Series Plot") +
  coord_cartesian(ylim=c(3,310)) +
  theme_bw() 

Metadata_sample_time_series_p_absolute
Metadata_sample_time_series_p_absolute_esthetics_improve
Metadata_sample_time_series_p_absolute_esthetics_improve_bw
####################################
# 2. Esquisse: create Rcode based on interactive graph building
####################################
# 2.1 Create barplot of Sample collection
####################################
esquisser()

ggplot(Metadata_sample_time_series) +
  aes(x = Sample_date, fill = Diagnosis, group = Diagnosis) +
  geom_histogram(bins = 40L) +
  scale_fill_manual(values = list(CD = "#8080B4", UC = "#D9846B")) +
  labs(x = "date", y = "count", 
       title = "Samples collection IBD patients", subtitle = "Time series") +
  theme_bw() +
  theme(plot.title = element_text(size = 16L), 
        plot.subtitle = element_text(face = "bold"), plot.caption = element_text(face = "bold"), axis.title.y = element_text(face = "bold"), 
        axis.title.x = element_text(face = "bold")) +
  facet_wrap(vars(Diagnosis), scales = "free", nrow = 2L)
####################################
# 2.2 Create stacked plot of Date of Diagnosis
####################################
esquisser()

patient_date_diagnosis <- Metadata_sample_time_series %>%
  filter(Week >= 0 & Week <= 9.4) %>%
  ggplot() +
  aes(x = Date_diagnosis, fill = Diagnosis) +
  geom_histogram(bins = 30L) +
  scale_fill_manual(values = list(CD = "#8080B4", UC = "#D9846B")) +
  labs(x = "date", y = "count", 
       title = "IBD patients", subtitle = "date of diagnosis") +
  theme_classic() +
  theme(plot.title = element_text(size = 16L), 
        plot.subtitle = element_text(size = 12L, face = "bold"), axis.title.y = element_text(face = "bold"), 
        axis.title.x = element_text(face = "bold")) +
  coord_cartesian(ylim=c(1,23))

patient_date_diagnosis
####################################
# 2.3 Date of birth
####################################
# 2.3.1 Create stacked barplot of date of birth
####################################
esquisser()

stacked_barplot <- Metadata_sample_time_series %>%
  filter(Week >= 0 & Week <= 4.6) %>%
  ggplot() +
  aes(x = Date_of_birth, fill = Diagnosis) +
  geom_histogram(bins = 30L) +
  scale_fill_manual(values = list(CD = "#8080B4", UC = "#D9846B")) +
  labs(x = "date", y = "patients", 
       title = "IBD patients", subtitle = "date of birth") +
  theme_classic() +
  theme(plot.title = element_text(size = 16L), 
        plot.subtitle = element_text(face = "bold"), axis.title.y = element_text(face = "bold"), axis.title.x = element_text(face = "bold")) +
  coord_cartesian(ylim=c(0.61,13)) 
geom_density()

stacked_barplot
####################################
# 2.3.2 Evaluate normality of age with density plot
####################################
# For this I will need the variable "Age" and not written as a date.
# It is important that this has point and no comma's, as you cannot make density plots out of this
# You could also use esquisser for density plots ofcourse.

esquisser()

density_plot <- Metadata_sample_time_series %>%
  filter(Week >= 0 & Week <= 9.2) %>%
  ggplot() +
  aes(x = Age, fill = Diagnosis) +
  geom_density(adjust = 1L) +
  scale_fill_manual(values = list(
    CD = "#8080B4", UC = "#D9846B")) +
  labs(x = "age", y = "density", title = "IBD patients", subtitle = "density plot") +
  theme_bw() +
  theme(plot.title = element_text(size = 16L), plot.subtitle = element_text(size = 12L, 
                                                                            face = "bold"), axis.title.y = element_text(face = "bold"), axis.title.x = element_text(face = "bold")) +
  facet_wrap(vars(Diagnosis)) +
  coord_cartesian(ylim=c(0.001,0.023)) 

density_plot
####################################
# 2.3.3 Normality assumption
####################################
# IBD
ggqqplot(Metadata_sample_time_series$Age)
shapiro.test(Metadata_sample_time_series$Age)

# CD
ggqqplot(Metadata_sample_time_series$Age[Metadata_sample_time_series$Diagnosis == "CD"])
shapiro.test(Metadata_sample_time_series$Age[Metadata_sample_time_series$Diagnosis == "CD"])

# UC
ggqqplot(Metadata_sample_time_series$Age[Metadata_sample_time_series$Diagnosis == "UC"])
shapiro.test(Metadata_sample_time_series$Age[Metadata_sample_time_series$Diagnosis == "UC"])

# From the output, the p-value > 0.05 implying that the distribution of the data are not significantly different from normal distribution. 
# In other words, in none of our cases we can assume the normality. 
# therefore, our Age is a confounding variable (also based on previous theory).
####################################
# 3. Remission
####################################
# 3.1 Endoscopic remission
####################################
Metadata_1 <- Metadata
Metadata_1$Randomized.Nr. <- NULL
Metadata_1$Sample_date <- NULL
Metadata_1$Date_of_birth <- NULL
Metadata_1$Sample.dilution.color <- NULL
Metadata_1$time_sample_to_sequencing <- NULL
Metadata_1$Clinical_outcome_combined[is.na(Metadata_1$Clinical_outcome_combined)] <- "unknown"
Metadata_1$Biomarker_outcome_combined[is.na(Metadata_1$Biomarker_outcome_combined)] <- "unknown"
Metadata_1$Endoscopic_outcome_combined[is.na(Metadata_1$Endoscopic_outcome_combined)] <- "unknown"

Metadata_1$Endoscopic_outcome_combined[Metadata_1$Endoscopic_outcome_combined == "1"] <- "Remission"
Metadata_1$Endoscopic_outcome_combined[Metadata_1$Endoscopic_outcome_combined == "0"] <- "Non-Remission"

esquisser()
ggThemeAssist

View(Metadata_1)
# Nr. patients
# All input samples are paired and therefore by looking at all w0 or baseline sample you can count the number of patients present.
# 
nrow(Metadata_1[Metadata_1$Week == 0,]) #191
nrow(Metadata_1[Metadata_1$Therapy_1 == "IFX" & Metadata_1$Week == 0,]) # 51
nrow(Metadata_1[Metadata_1$Therapy_1 == "ADA" & Metadata_1$Week == 0,]) # 28
nrow(Metadata_1[Metadata_1$Therapy_1 == "VDZ" & Metadata_1$Week == 0,]) # 63
nrow(Metadata_1[Metadata_1$Therapy_1 == "UST" & Metadata_1$Week == 0,]) # 49
51+28+63+49

nrow(Metadata_1[Metadata_1$Week == 0 & Metadata_1$Diagnosis == "UC",]) # 60
nrow(Metadata_1[Metadata_1$Week == 14 & Metadata_1$Diagnosis == "CD",]) # 50

## Add R+NR based on combination clinical & biomarker:
Metadata_1$Endoscopic_outcome_combined[c(Metadata_1$Endoscopic_outcome_combined == "unknown" & Metadata_1$Clinical_outcome_combined == "1" & Metadata_1$Biomarker_outcome_combined == "1")] <- "Remission"
Metadata_1$Endoscopic_outcome_combined[c(Metadata_1$Endoscopic_outcome_combined == "unknown" & Metadata_1$Clinical_outcome_combined == "0" & Metadata_1$Biomarker_outcome_combined == "0")] <- "Non-Remission"

View(Metadata_1)

Treatment <- Metadata_1 %>%
  filter(Week < 14L) %>%
  ggplot() +
  aes(x = Endoscopic_outcome_combined, fill = Diagnosis) +
  geom_bar() +
  scale_fill_manual(values = list(
    CD = "#8080B4", UC = "#D9846B")) +
  labs(y = "patients (#)", subtitle = "treatment") +
  theme_bw() +
  theme(plot.subtitle = element_text(size = 12L, face = "bold"), axis.title.y = element_text(face = "bold")) +
  facet_wrap(vars(Therapy_1)) + theme(panel.grid.major = element_line(colour = "gray96"), 
        panel.grid.minor = element_line(colour = "gray97", 
                                        linetype = "blank"), axis.title = element_text(family = "serif"), 
        plot.title = element_text(family = "serif"), 
        legend.title = element_text(face = "bold")) +labs(x = NULL)

Treatment

Treatment_endoscopy <- Treatment + theme(plot.subtitle = element_text(size = 10, 
                                                                      colour = "gray0"), plot.title = element_text(size = 16)) +labs(title = "Treatment", x = NULL, subtitle = "endoscopic remission")
Treatment_endoscopy
####################################
# 3.2 Biomarker remission
####################################
Metadata_1$Biomarker_outcome_combined[Metadata_1$Biomarker_outcome_combined == "1"] <- "Remission"
Metadata_1$Biomarker_outcome_combined[Metadata_1$Biomarker_outcome_combined == "0"] <- "Non-Remission"
# View(Metadata_1)

esquisser()
ggThemeAssist

Treatment <- Metadata_1 %>%
  filter(Week < 14L) %>%
  ggplot() +
  aes(x = Biomarker_outcome_combined, fill = Diagnosis) +
  geom_bar() +
  scale_fill_manual(values = list(
    CD = "#8080B4", UC = "#D9846B")) +
  labs(y = "patients (#)", subtitle = "treatment") +
  theme_bw() +
  theme(plot.subtitle = element_text(size = 12L, face = "bold"), axis.title.y = element_text(face = "bold")) +
  facet_wrap(vars(Therapy_1))

Treatment

Treatment_biomarker <- Treatment + theme(plot.subtitle = element_text(size = 10, 
                                                                      colour = "gray0"), plot.title = element_text(size = 16)) +labs(title = "Treatment", x = NULL, subtitle = "biomarker remission")
Treatment_biomarker
###################################
# 3.3 Clinical remission
####################################
Metadata_1$Clinical_outcome_combined[Metadata_1$Clinical_outcome_combined == "1"] <- "Remission"
Metadata_1$Clinical_outcome_combined[Metadata_1$Clinical_outcome_combined == "0"] <- "Non-Remission"
# View(Metadata_1)

Treatment <- Metadata_1 %>%
  filter(Week < 14L) %>%
  ggplot() +
  aes(x = Clinical_outcome_combined, fill = Diagnosis) +
  geom_bar() +
  scale_fill_manual(values = list(
    CD = "#8080B4", UC = "#D9846B")) +
  labs(y = "patients (#)", subtitle = "treatment") +
  theme_bw() +
  theme(plot.subtitle = element_text(size = 12L, face = "bold"), axis.title.y = element_text(face = "bold")) +
  facet_wrap(vars(Therapy_1))

Treatment

Treatment_clinical <- Treatment + theme(plot.subtitle = element_text(size = 10, 
                                                                     colour = "gray0"), plot.title = element_text(size = 16)) +labs(title = "Treatment", x = NULL, subtitle = "clinical remission")
Treatment_clinical
####################################
# 4. Bacterial cell count
####################################
## If you want to show only patients then selected week >= 0 & week < 14.
## If not you will take all samples, but by doing this you subset based on baseline samples, which are 100% present in all of the samples and are unique.

Metadata_1$bacterial_Cellcount <- as.numeric(gsub(",", ".", Metadata_1$bacterial_Cellcount))
Metadata_1$Week[Metadata_1$Week == "0"] <- "Baseline"
Metadata_1$Week[Metadata_1$Week == "14"] <- "week 14"
Metadata_1$Week[Metadata_1$Week == "24"] <- "week 24"
Metadata_1$Week[Metadata_1$Week == "NA"] <- "unknown"

esquisser()

bacterial_cellcount_remission_endo <- Metadata_1 %>%
 filter(bacterial_Cellcount >= 5394089318 & bacterial_Cellcount <= 230098998194 | is.na(bacterial_Cellcount)) %>%
 filter(!is.na(Endoscopic_outcome_combined)) %>%
 ggplot() +
 aes(x = Week, y = bacterial_Cellcount, fill = Diagnosis) +
 geom_boxplot(shape = "circle") +
 scale_fill_manual(values = list(CD = "#8080B4", UC = "#D9846B")) +
 labs(x = "week", y = "bacterial cellcount", 
 title = "IBD patients", subtitle = "endoscopical remission") +
 theme_bw() +
 theme(legend.position = "bottom", 
 plot.title = element_text(size = 16L), plot.subtitle = element_text(size = 12L, face = "bold"), axis.title.y = element_text(face = "bold"), 
 axis.title.x = element_text(face = "bold")) +
 facet_wrap(vars(Endoscopic_outcome_combined))

bacterial_cellcount_remission_endo

## Nr. patients in each group
# UC
nrow(Metadata_1[Metadata_1$Diagnosis == "UC",]) # 120 patients
nrow(Metadata_1[Metadata_1$Endoscopic_outcome_combined == "Remission" & Metadata_1$Diagnosis == "UC",]) # 66 patients
nrow(Metadata_1[Metadata_1$Endoscopic_outcome_combined == "Remission" & Metadata_1$Diagnosis == "UC" & Metadata_1$Week == "Baseline",]) # 33 data points
nrow(Metadata_1[Metadata_1$Endoscopic_outcome_combined == "Remission" & Metadata_1$Diagnosis == "UC" & Metadata_1$Week == "week 14",]) # 33 data points

nrow(Metadata_1[Metadata_1$Endoscopic_outcome_combined == "Non-Remission" & Metadata_1$Diagnosis == "UC",]) ## 54 patients
nrow(Metadata_1[Metadata_1$Endoscopic_outcome_combined == "Non-Remission" & Metadata_1$Diagnosis == "UC" & Metadata_1$Week == "Baseline",]) ## 27 patients
nrow(Metadata_1[Metadata_1$Endoscopic_outcome_combined == "Non-Remission" & Metadata_1$Diagnosis == "UC" & Metadata_1$Week == "week 14",]) ## 27 patients

# CD 
nrow(Metadata_1[Metadata_1$Diagnosis == "CD",]) # 312 patients
nrow(Metadata_1[Metadata_1$Endoscopic_outcome_combined == "Remission" & Metadata_1$Diagnosis == "CD",]) # 187 patients
nrow(Metadata_1[Metadata_1$Endoscopic_outcome_combined == "Remission" & Metadata_1$Diagnosis == "CD" & Metadata_1$Week == "Baseline",]) # 75 data points
nrow(Metadata_1[Metadata_1$Endoscopic_outcome_combined == "Remission" & Metadata_1$Diagnosis == "CD" & Metadata_1$Week == "week 14",]) # 37 data points
nrow(Metadata_1[Metadata_1$Endoscopic_outcome_combined == "Remission" & Metadata_1$Diagnosis == "CD" & Metadata_1$Week == "week 24",]) # 75 data points

nrow(Metadata_1[Metadata_1$Endoscopic_outcome_combined == "Non-Remission" & Metadata_1$Diagnosis == "CD",]) ## 194 patients
nrow(Metadata_1[Metadata_1$Endoscopic_outcome_combined == "Non-Remission" & Metadata_1$Diagnosis == "CD" & Metadata_1$Week == "Baseline",]) ## 84 patients
nrow(Metadata_1[Metadata_1$Endoscopic_outcome_combined == "Non-Remission" & Metadata_1$Diagnosis == "CD" & Metadata_1$Week == "week 14",]) ## 24 patients
nrow(Metadata_1[Metadata_1$Endoscopic_outcome_combined == "Non-Remission" & Metadata_1$Diagnosis == "CD" & Metadata_1$Week == "week 24",]) ## 85 patients

## Finding: decrease bacterial cell count at baseline in patients that are receptive for treatment
## 1. Baseline endoscopical bacterial count: density plot
Metadata_1$active_disease_baseline <- NULL
esquisser()

baseline_bacterial_count_density_w0_w14 <- Metadata_1 %>%
 filter(bacterial_Cellcount >= 12415541742 & bacterial_Cellcount <= 5.25e+11 | is.na(bacterial_Cellcount)) %>%
 filter(!is.na(Clinical_outcome_combined)) %>%
 filter(Endoscopic_outcome_combined %in% "Remission") %>%
 filter(!is.na(Biomarker_outcome_combined)) %>%
 ggplot() +
 aes(x = bacterial_Cellcount, fill = Week) +
 geom_density(adjust = 1L) +
 scale_fill_brewer(palette = "Dark2", 
 direction = 1) +
 labs(x = "bacterial cellcount", y = "density", title = "shift in bacterial cell count", 
 subtitle = "remission") +
 theme_bw() +
 theme(plot.title = element_text(size = 16L, face = "bold"), 
 plot.subtitle = element_text(size = 12L, face = "bold"), axis.title.y = element_text(face = "bold"), 
 axis.title.x = element_text(face = "bold")) +
 facet_wrap(vars(Diagnosis), scales = "free")

baseline_bacterial_count_density_w0_w14 + theme(axis.ticks = element_line(colour = "gray13"), 
    axis.text.y = element_text(size = 0),
    axis.ticks.y = element_blank(),
    panel.background = element_rect(fill = NA), 
    plot.background = element_rect(linetype = "solid"))



baseline_bacterial_count_density <- Metadata_1 %>%
 filter(Week %in% "Baseline") %>%
 filter(bacterial_Cellcount >= 14170904848 & bacterial_Cellcount <= 
 5.25e+11 | is.na(bacterial_Cellcount)) %>%
 filter(!is.na(Endoscopic_outcome_combined)) %>%
 ggplot() +
 aes(x = bacterial_Cellcount, fill = Diagnosis) +
 geom_density(adjust = 1L) +
 scale_fill_manual(values = list(
 CD = "#8080B4", UC = "#D9846B")) +
 labs(x = "bacterial cell count", y = "density", title = "Endoscopical remission", 
 subtitle = "baseline") +
 theme_bw() +
 theme(plot.title = element_text(size = 16L, face = "bold"), 
 plot.subtitle = element_text(size = 12L, face = "bold"), axis.title.y = element_text(face = "bold"), 
 axis.title.x = element_text(face = "bold")) +
 facet_wrap(vars(Endoscopic_outcome_combined))

baseline_bacterial_count_density

## 2. Reproduce with biomarker & clinical remission which shows the same thing.
## you could do this later to make claim more prominent! now it doesn't really matter.
####################################
# 5. Moisture(%)
####################################
Metadata_1$Moisture <- as.numeric(gsub(",", ".", Metadata_1$Moisture))

moisture_endoscopical_remission_box <- Metadata_1 %>%
 filter(Moisture >= 63L & Moisture <= 99L | is.na(Moisture)) %>%
 filter(!is.na(Endoscopic_outcome_combined)) %>%
 ggplot() +
 aes(x = Week, y = Moisture, fill = Diagnosis) +
 geom_boxplot(shape = "circle") +
 scale_fill_manual(values = list(CD = "#8080B4", UC = "#D9846B")) +
 theme_bw() +
 theme(legend.position = "bottom") +
 facet_wrap(vars(Endoscopic_outcome_combined))

moisture_endoscopical_remission_box

moisture_endoscopical_remission <- Metadata_1 %>%
 filter(Moisture >= 63L & Moisture <= 99L | is.na(Moisture)) %>%
 filter(!is.na(Endoscopic_outcome_combined)) %>%
 ggplot() +
 aes(x = Week, y = Moisture, fill = Diagnosis) +
 geom_violin(adjust = 1L, scale = "area") +
 scale_fill_manual(values = list(CD = "#8080B4", UC = "#D9846B")) +
 theme_bw() +
 theme(legend.position = "bottom") +
 facet_wrap(vars(Endoscopic_outcome_combined))

moisture_endoscopical_remission
####################################
# 6. Donut plot: remission/non-remission
####################################
# Create two donut plots to show the patients in remission/non-remssion in percentages
# Put in the middle the absolute number of patients (CD) & (UC).
# Put in excel next to eachother.

# https://www.r-graph-gallery.com/128-ring-or-donut-plot.html

Metadata_2 <- Metadata_1
unique(Metadata_2$Week)
View(Metadata_2)

unique(Metadata$patient_ID[Metadata$Endoscopic_outcome_combined == "unknown"])
unique(Metadata_2$patient_ID[Metadata_2$Endoscopic_outcome_combined == "unknown"])
# In total we have data of 176 patients, in which we use endoscopical remission to determine response of the patient at primary endpoint.
# In 29 patients the endoscopical response data was not available
# For these patients we used the combination of clinical and biomarker response, but only if those results were consistent.
# After implementing that we don't have remission response data for only 18 patients.

# Add count 
Metadata_2$count <- "1"

# Check if correct
nrow(Metadata_2[Metadata_2$Diagnosis == "CD" & Metadata_2$Week== "Baseline",]) # 131 CD patients
nrow(Metadata_2[Metadata_2$Diagnosis == "CD" & Metadata_2$Week == "Baseline" & Metadata_2$Endoscopic_outcome_combined == "Remission",]) # 54 patients
nrow(Metadata_2[Metadata_2$Diagnosis == "CD" & Metadata_2$Week == "Baseline" & Metadata_2$Endoscopic_outcome_combined == "Non-Remission",]) #  60 patients
nrow(Metadata_2[Metadata_2$Diagnosis == "CD" & Metadata_2$Week == "Baseline" & Metadata_2$Endoscopic_outcome_combined == "unknown",]) # 17 patients
54+60+17
unique(Metadata_2$Endoscopic_outcome_combined)
str(Metadata_2$Endoscopic_outcome_combined)

nrow(Metadata_2[Metadata_2$Diagnosis == "UC" & Metadata_2$Week == "Baseline",]) # 60 UC patients
nrow(Metadata_2[Metadata_2$Diagnosis == "UC" & Metadata_2$Week == "Baseline" & Metadata_2$Endoscopic_outcome_combined == "Remission",]) ## 33
nrow(Metadata_2[Metadata_2$Diagnosis == "UC" & Metadata_2$Week == "Baseline" & Metadata_2$Endoscopic_outcome_combined == "Non-Remission",]) # 27
nrow(Metadata_2[Metadata_2$Diagnosis == "UC" & Metadata_2$Week == "Baseline" & Metadata_2$Endoscopic_outcome_combined == "unknown",]) # 0
33+27+0

# Subset baseline (~ unique patient 1 sample ~ 1 patient)
Metadata_donut_plot <- Metadata_2[Metadata_2$Week == "Baseline",]
Metadata_donut_plot_subset <- Metadata_donut_plot[,c("Endoscopic_outcome_combined","Diagnosis", "count")]
Metadata_donut_plot_subset_CD <- Metadata_donut_plot_subset[Metadata_donut_plot_subset$Diagnosis == "CD",]
Metadata_donut_plot_subset_UC <- Metadata_donut_plot_subset[Metadata_donut_plot_subset$Diagnosis == "UC",]

# Create frequency for donut plot
Metadata_donut_plot_subset_CD <- data.frame(with(Metadata_donut_plot_subset_CD, table(Endoscopic_outcome_combined)))
Metadata_donut_plot_subset_UC <- data.frame(with(Metadata_donut_plot_subset_UC, table(Endoscopic_outcome_combined)))
View(Metadata_donut_plot_subset_UC)

## CD (131 patients)
str(Metadata_donut_plot_subset_CD)
Metadata_donut_plot_subset_CD$Freq <- as.numeric(Metadata_donut_plot_subset_CD$Freq)
Metadata_donut_plot_subset_CD$fraction <- "1"
Metadata_donut_plot_subset_CD$fraction <- Metadata_donut_plot_subset_CD$Freq/131
Metadata_donut_plot_subset_CD

## UC (60 patients)
str(Metadata_donut_plot_subset_UC)
Metadata_donut_plot_subset_UC$Freq <- as.numeric(Metadata_donut_plot_subset_UC$Freq)
Metadata_donut_plot_subset_UC$fraction <- "1"
Metadata_donut_plot_subset_UC$fraction <- Metadata_donut_plot_subset_UC$Freq/60
Metadata_donut_plot_subset_UC

# Percentages
Metadata_donut_plot_subset_CD$percentage <- Metadata_donut_plot_subset_CD$fraction*100
Metadata_donut_plot_subset_UC$percentage <- Metadata_donut_plot_subset_UC$fraction*100

# Labels
#Metadata_donut_plot_subset_CD$percentage <- formatC(Metadata_donut_plot_subset_CD$percentage, digits = 3)
#Metadata_donut_plot_subset_UC$percentage <- formatC(Metadata_donut_plot_subset_UC$percentage, digits = 3)
Metadata_donut_plot_subset_CD
str(Metadata_donut_plot_subset_UC)

## Save colors
colors <- c("#b80000","#b8a100","#b8a1ae")

# Make the donut plot
# CD
donut_plot_CD_remission <- ggplot(Metadata_donut_plot_subset_CD, aes(x = 2, y = percentage, fill = Endoscopic_outcome_combined)) +
  geom_bar(stat = "identity") +
  coord_polar("y", start = 1) +
  geom_text(aes(label = paste0(round(percentage), "%")),
            position = position_stack(vjust = 0.5)) +
  theme_void() +
  scale_fill_brewer(palette = "Dark2") +
  xlim(.2,2.5) +
  theme(legend.title = element_blank())  +
  # theme(legend.position = "none")
  theme(legend.key.size = unit(2,"line"), legend.text = element_text(size = 12))

donut_plot_CD_remission

# UC
donut_plot_UC_remission <- ggplot(Metadata_donut_plot_subset_UC, aes(x = 2, y = percentage, fill = Endoscopic_outcome_combined)) +
  geom_bar(stat = "identity") +
  coord_polar("y", start = 1) +
  geom_text(aes(label = paste0(round(percentage), "%")),
            position = position_stack(vjust = 0.5)) +
  theme_void() +
  scale_fill_brewer(palette = "Dark2") +
  xlim(.2,2.5) +
  theme(legend.title = element_blank())  +
  theme(legend.position = "none")
  # theme(legend.key.size = unit(2,"line"), legend.text = element_text(size = 12))

donut_plot_UC_remission

donut_plot_UC_remission + donut_plot_CD_remission
####################################

######### !!!!!!!!!!!!!!! #########
# SAVE Global Environment #
######### !!!!!!!!!!!!!!! #########