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
install.packages("ggplot2")
install.packages("dplyr")
install.packages("plyr")
install.packages("reshape")
install.packages("reshape2")
install.packages("scales")
install.packages("viridis")
install.packages("readr")
install.packages("WriteXLS")
install.packages("readxl")
install.packages("devtools")
library(devtools) # before downlading from github you need to load devtools
install.packages("BiocManager")
library(BiocManager)
BiocManager::install("ComplexHeatmap")
# install_github("jokergoo/ComplexHeatmap")
library(ComplexHeatmap)

install.packages("dendsort")
install.packages("vegan")
install.packages("ape")
install.packages("devtools")
install.packages("seriation")
install.packages("taxonomizr")
install.packages("ggThemeAssist")
install.packages("esquisse")
installed.packages("modeldata")
install.packages("tidyverse")
install.packages("ggpubr")
install.packages("ggforce")
install.packages("rJava")
install.packages("UpSetR")
install.packages("tidyverse")
install.packages("venneuler")
install.packages("grid")

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
library(grid)
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
library(plotly)
library(tidyverse)
library(timetk)
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
# 2.1 Work with Esquise to digg into metadata
####################################
# Esquisser is a program written in which you can get the Rcode for ggplot via interactive browsing of your dataframe. 
# This is very easy and user-friendly, especially good for beginners.

# You just have to press the next command and it opens all the interactive options for you.
# Next, I obtained the following code.

esquisser()
####################################
# 2.1.1 Create barplot of Sample collection
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
# 2.1.2 Create stacked plot of Date of Diagnosis
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
# 2.1.3 Date of birth
####################################
# 2.1.3.1 Create stacked barplot of date of birth
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
# 2.1.3.2 Evaluate normality of age with density plot
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
# 2.1.3.3 Normality assumption
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
# 2.1.4.1 Endoscopic remission
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
# 2.1.4.2 Biomarker remission
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
# 2.1.4.3 Clinical remission
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

esquisser()


library(dplyr)
library(ggplot2)

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
# 7. Show similarities and differences between groups
####################################
# A Venn diagram is an illustration that uses circles to show the relationships among things or finite groups of things. Circles that overlap have a commonality while circles that do not overlap do not share those traits. Venn diagrams help to visually represent the similarities and differences between two concepts
####################################
# 7.1 Venn Diagram: overlapping remission indicators (UC/CD)
####################################
# https://www.littlemissdata.com/blog/set-analysis
# Venn Diagram is used for 4 or less groups within the matrix.
####################################
# 7.1.1 Install and load packages
####################################
install.packages("rJava")
install.packages("UpSetR")
install.packages("tidyverse")
install.packages("venneuler")
install.packages("grid")

library(rJava)
library(UpSetR)
library(tidyverse)
library(venneuler)
library(grid)
####################################
# 7.1.2 Rename columns to make it more readable
####################################
colnames(Metadata_2)
names(Metadata_2)[names(Metadata_2) == "Endoscopic_outcome_combined"] <- "endoscopical outcome"
names(Metadata_2)[names(Metadata_2) == "Clinical_outcome_combined"] <- "clinical outcome"
names(Metadata_2)[names(Metadata_2) == "Biomarker_outcome_combined"] <- "biomarker outcome"
names(Metadata_2)[names(Metadata_2) == "Therapy_1"] <- "therapy"
names(Metadata_2)[names(Metadata_2) == "Week"] <- "week"
names(Metadata_2)[names(Metadata_2) == "Diagnosis"] <- "diagnosis"
####################################
# 7.1.3 Subset columns you want to make a Venn Diagram
####################################
Metadata_2_venn_set <- Metadata_2[,c("diagnosis", "week", "therapy", "endoscopical outcome", "clinical outcome",  "biomarker outcome")]
colnames(Metadata_2_venn_set)
View(Metadata_2_venn_set)
####################################
# 7.1.4 Venn diagram:
####################################
# 7.1.4.1 UC (n = 60)
####################################
# Calculate patients in remission in the specific category
nrow(Metadata_2_venn_set[Metadata_2_venn_set$diagnosis == "UC" & Metadata_2_venn_set$week == "Baseline" & Metadata_2_venn_set$`endoscopical outcome` == "Remission",]) # 33 patients
nrow(Metadata_2_venn_set[Metadata_2_venn_set$diagnosis == "UC" & Metadata_2_venn_set$week == "Baseline" & Metadata_2_venn_set$`clinical outcome` == "Remission",]) # 37
nrow(Metadata_2_venn_set[Metadata_2_venn_set$diagnosis == "UC" & Metadata_2_venn_set$week == "Baseline" & Metadata_2_venn_set$`biomarker outcome` == "Remission",]) # 26

# Calculate patients in remission in combined categories
nrow(Metadata_2_venn_set[Metadata_2_venn_set$diagnosis == "UC" & Metadata_2_venn_set$week == "Baseline" & Metadata_2_venn_set$`endoscopical outcome` == "Remission" & Metadata_2_venn_set$`clinical outcome` == "Remission",]) # 30 patients
nrow(Metadata_2_venn_set[Metadata_2_venn_set$diagnosis == "UC" & Metadata_2_venn_set$week == "Baseline" & Metadata_2_venn_set$`endoscopical outcome` == "Remission" & Metadata_2_venn_set$`biomarker outcome` == "Remission",]) # 20 patients
nrow(Metadata_2_venn_set[Metadata_2_venn_set$diagnosis == "UC" & Metadata_2_venn_set$week == "Baseline" & Metadata_2_venn_set$`clinical outcome` == "Remission" & Metadata_2_venn_set$`biomarker outcome` == "Remission",]) # 23 patients
nrow(Metadata_2_venn_set[Metadata_2_venn_set$diagnosis == "UC" & Metadata_2_venn_set$week == "Baseline" & Metadata_2_venn_set$`clinical outcome` == "Remission" & Metadata_2_venn_set$`biomarker outcome` == "Remission" & Metadata_2_venn_set$`endoscopical outcome`=="Remission",]) # 19 patients

# Set input data Venn Diagram for 3 categories (having remission/non-remission as result)
# The problem you'd have now is that if you combine two categories such as 'endoscpy&biomarker' it will take a aprt of each circle & therefore you have to multiply the value by 2.
expressionInput <- c(`endoscopy` = 33, `clinical` = 37, `biomarker` = 26, `endoscopy&clinical` = 30, `endoscopy&biomarker` = 20, `clinical&biomarker` = 23, `endoscopy&clinical&biomarker` = 19)
expressionInput2 <- c(`endoscopy` = 33, `clinical` = 37, `biomarker` = 26, `endoscopy&clinical` = 60, `endoscopy&biomarker` = 40, `clinical&biomarker` = 46, `endoscopy&clinical&biomarker` = 57)
expressionInput2

myExpVenn_with_title <- venneuler(expressionInput2) 
par(cex=0.95)
plot(myExpVenn_with_title, main = "remission indicator")
grid.text(
  "UC",
  x = 0.52,
  y = 0.15,
  gp = gpar(
    fontsize = 15,
    fontface = 2
  )
)

myExpVenn_without_title <- venneuler(expressionInput2) 
par(cex=0.95)
plot(myExpVenn_without_title, main = "")
grid.text(
  "UC",
  x = 0.52,
  y = 0.15,
  gp = gpar(
    fontsize = 15,
    fontface = 2
  )
)
####################################
# 7.1.4.2 CD (n = 131)
####################################
View(Metadata_2_venn_set) # also quite some unknowns here

# Calculate patients in remission in the specific category
nrow(Metadata_2_venn_set[Metadata_2_venn_set$diagnosis == "CD" & Metadata_2_venn_set$week == "Baseline" & Metadata_2_venn_set$`endoscopical outcome` == "Remission",]) # 46 patients
nrow(Metadata_2_venn_set[Metadata_2_venn_set$diagnosis == "CD" & Metadata_2_venn_set$week == "Baseline" & Metadata_2_venn_set$`endoscopical outcome` == "unknown",]) # 29 patients

nrow(Metadata_2_venn_set[Metadata_2_venn_set$diagnosis == "CD" & Metadata_2_venn_set$week == "Baseline" & Metadata_2_venn_set$`clinical outcome` == "Remission",]) # 86 patients
nrow(Metadata_2_venn_set[Metadata_2_venn_set$diagnosis == "CD" & Metadata_2_venn_set$week == "Baseline" & Metadata_2_venn_set$`clinical outcome` == "unknown",]) # 0 patients

nrow(Metadata_2_venn_set[Metadata_2_venn_set$diagnosis == "CD" & Metadata_2_venn_set$week == "Baseline" & Metadata_2_venn_set$`biomarker outcome` == "Remission",]) # 43 patients
nrow(Metadata_2_venn_set[Metadata_2_venn_set$diagnosis == "CD" & Metadata_2_venn_set$week == "Baseline" & Metadata_2_venn_set$`biomarker outcome` == "unknown",]) # 4 patients

# Calculate patients in remission in combined categories
nrow(Metadata_2_venn_set[Metadata_2_venn_set$diagnosis == "CD" & Metadata_2_venn_set$week == "Baseline" & Metadata_2_venn_set$`endoscopical outcome` == "Remission" & Metadata_2_venn_set$`clinical outcome` == "Remission",]) # 39 patients
nrow(Metadata_2_venn_set[Metadata_2_venn_set$diagnosis == "CD" & Metadata_2_venn_set$week == "Baseline" & Metadata_2_venn_set$`endoscopical outcome` == "Remission" & Metadata_2_venn_set$`biomarker outcome` == "Remission",]) # 26 patients
nrow(Metadata_2_venn_set[Metadata_2_venn_set$diagnosis == "CD" & Metadata_2_venn_set$week == "Baseline" & Metadata_2_venn_set$`clinical outcome` == "Remission" & Metadata_2_venn_set$`biomarker outcome` == "Remission",]) # 30 patients
nrow(Metadata_2_venn_set[Metadata_2_venn_set$diagnosis == "CD" & Metadata_2_venn_set$week == "Baseline" & Metadata_2_venn_set$`clinical outcome` == "Remission" & Metadata_2_venn_set$`biomarker outcome` == "Remission" & Metadata_2_venn_set$`endoscopical outcome`=="Remission",]) # 21 patients

# Set input data Venn Diagram for 3 categories (having remission/non-remission as result)
# The problem you'd have now is that if you combine two categories such as 'endoscpy&biomarker' it will take a aprt of each circle & therefore you have to multiply the value by 2.
expressionInput <- c(`endoscopy` = 46, `clinical` = 86, `biomarker` = 43, `endoscopy&clinical` = 29, `endoscopy&biomarker` = 26, `clinical&biomarker` = 30, `endoscopy&clinical&biomarker` = 21)
expressionInput3 <- c(`endoscopy` = 46, `clinical` = 86, `biomarker` = 43, `endoscopy&clinical` = 58, `endoscopy&biomarker` = 52, `clinical&biomarker` = 60, `endoscopy&clinical&biomarker` = 42)

myExpVenn_CD <- venneuler(expressionInput3) 
par(cex=0.95)
plot(myExpVenn_CD, main = "")
grid.text(
  "CD",
  x = 0.52,
  y = 0.15,
  gp = gpar(
    fontsize = 15,
    fontface = 2
  )
)

myExpVenn_without_title
####################################
# 7.2 UpSetR plot
####################################
# For showing more than 4 groups we like to use UpSetR plot because Venn Diagrams become to complex.
# The great thing is that we can also create an UpSet plot using the same basic expression list.  You simply pass the fromExpression() function into the upset() function.  The remaining code is to format the labels and font size.
# How to read an UpSet plot:  UpSet plots offer a straight forward way for us to view set data by frequency.  On the bottom left hand side horizontal bar chart, we show the entire size of each set.  In this case, each set is of size 8.  The vertical bar chart on the upper right hand side shows the sizes of isolated set participation.  In the example, 5 values only belong to the #rstats set or only belong to the memes set.  3 values belong to both sets.  
####################################
# 7.2.1 Simple UpsetR plot (~ simple input)
####################################
# Here you can put yourself numbers as input to make a simple UpsetR plot
####################################
# 7.2.1.1 UC
####################################
expressionInput_UC <- c(`endoscopy` = 33, `clinical` = 37, `biomarker` = 26, `endoscopy&clinical` = 30, `endoscopy&biomarker` = 20, `clinical&biomarker` = 23, `endoscopy&clinical&biomarker` = 19)

# Create an UpsetR Plot for UC
upset(fromExpression(expressionInput_UC), order.by = "freq")
grid.text(
  "UC patients",
  x = 0.80,
  y = 0.05,
  gp = gpar(
    fontsize = 13,
    fontface = 2
  )
)
####################################
# 7.2.1.2 CD
####################################
expressionInput_CD <- c(`endoscopy` = 46, `clinical` = 86, `biomarker` = 43, `endoscopy&clinical` = 29, `endoscopy&biomarker` = 26, `clinical&biomarker` = 30, `endoscopy&clinical&biomarker` = 21)

# Create an UpsetR Plot for CD
upset(fromExpression(expressionInput_CD), order.by = "freq")
grid.text(
  "CD patients",
  x = 0.80,
  y = 0.05,
  gp = gpar(
    fontsize = 13,
    fontface = 2
  )
)
####################################
# 7.3 Complicated UpsetR plot (~ matrix input with integers)
####################################
# https://www.littlemissdata.com/blog/set-analysis
# https://cran.r-project.org/web/packages/UpSetR/vignettes/basic.usage.html
# https://jokergoo.github.io/ComplexHeatmap-reference/book/upset-plot.html

# Or you make an UpSetR plot from this package, or you can integrate it in the ComplexHeatmap package.
# The latter one offers a lot of annotations.
####################################
# 7.3.1 UpSetR package
####################################
# 7.3.1.1 Create matrix
####################################
# The point of Set analysis is that you have to create sets. This is often time-consuming.
View(Metadata_2)

Metadata_2$`endoscopical outcome`[is.na(Metadata_2$`endoscopical outcome`)] <- "unknown"
Metadata_2$`endoscopical outcome`
View(Metadata)

# Remove "unknowns"
Metadata_2 <- Metadata_2[!c(Metadata_2$`endoscopical outcome` == "unknown"),]
Metadata_2 <- Metadata_2[!c(Metadata_2$`biomarker outcome` == "unknown"),]

Metadata_2_setUp <- Metadata_2[Metadata_2$week == "Baseline",]
unique(Metadata_2$patient_ID[Metadata_2$diagnosis == "CD"]) # 100 CD patients with 3R results
unique(Metadata_2$patient_ID[Metadata_2$diagnosis == "UC"]) # 55 UC patients with 3 R results
View(Metadata_2)

Metadata_2_setUp <- Metadata_2_setUp[, c("diagnosis", "therapy","endoscopical outcome","clinical outcome","biomarker outcome")]
nrow(Metadata_2_setUp) # 162 patients unique with data available (102 CD + 60 UC)
nrow(Metadata_2_setUp[Metadata_2_setUp$diagnosis == "UC",]) # 60 UC patients with known endoscopical outcome
nrow(Metadata_2_setUp[Metadata_2_setUp$diagnosis == "CD",]) # 114 CD patients with known endoscopical outcome
#View(Metadata_2_setUp)

# Convert diagnosis to matrix format
Metadata_2_setUp$`crohn's disease` <- 1
Metadata_2_setUp$`ulcerative colitis` <- 1
Metadata_2_setUp$`crohn's disease`[Metadata_2_setUp$diagnosis == "CD"] <- "1"
Metadata_2_setUp$`crohn's disease`[Metadata_2_setUp$diagnosis == "UC"] <- "0"
Metadata_2_setUp$`ulcerative colitis`[Metadata_2_setUp$diagnosis == "UC"] <- "1"
Metadata_2_setUp$`ulcerative colitis`[Metadata_2_setUp$diagnosis == "CD"] <- "0"

# Convert treatment to matrix format
Metadata_2_setUp$infliximab <- 1
Metadata_2_setUp$infliximab[Metadata_2_setUp$therapy == "IFX"] <- "1"
Metadata_2_setUp$infliximab[!Metadata_2_setUp$therapy == "IFX"] <- "0"

Metadata_2_setUp$adalimumab <- 1
Metadata_2_setUp$adalimumab[Metadata_2_setUp$therapy == "ADA"] <- "1"
Metadata_2_setUp$adalimumab[!Metadata_2_setUp$therapy == "ADA"] <- "0"

Metadata_2_setUp$vedolizumab <- 1
Metadata_2_setUp$vedolizumab[Metadata_2_setUp$therapy == "VDZ"] <- "1"
Metadata_2_setUp$vedolizumab[!Metadata_2_setUp$therapy == "VDZ"] <- "0"

Metadata_2_setUp$ustekinumab <- 1
Metadata_2_setUp$ustekinumab[Metadata_2_setUp$therapy == "UST"] <- "1"
Metadata_2_setUp$ustekinumab[!Metadata_2_setUp$therapy == "UST"] <- "0"

# Convert ALL response to matrix format
# endoscopical
Metadata_2_setUp$`remission (endoscopical)` <- 1
Metadata_2_setUp$`remission (endoscopical)`[Metadata_2_setUp$`endoscopical outcome` == "Remission"] <- "1"
Metadata_2_setUp$`remission (endoscopical)`[!Metadata_2_setUp$`endoscopical outcome` == "Remission"] <- "0"

Metadata_2_setUp$`non-remission (endoscopical)` <- 1
Metadata_2_setUp$`non-remission (endoscopical)`[Metadata_2_setUp$`endoscopical outcome` == "Non-Remission"] <- "1"
Metadata_2_setUp$`non-remission (endoscopical)`[!Metadata_2_setUp$`endoscopical outcome` == "Non-Remission"] <- "0"

# biomarker
Metadata_2_setUp$`remission (biomarker)` <- 1
Metadata_2_setUp$`remission (biomarker)`[Metadata_2_setUp$`biomarker outcome` == "Remission"] <- "1"
Metadata_2_setUp$`remission (biomarker)`[!Metadata_2_setUp$`biomarker outcome` == "Remission"] <- "0"

Metadata_2_setUp$`non-remission (biomarker)` <- 1
Metadata_2_setUp$`non-remission (biomarker)`[Metadata_2_setUp$`biomarker outcome` == "Non-Remission"] <- "1"
Metadata_2_setUp$`non-remission (biomarker)`[!Metadata_2_setUp$`biomarker outcome` == "Non-Remission"] <- "0"

# clinical
Metadata_2_setUp$`remission (clinical)` <- 1
Metadata_2_setUp$`remission (clinical)`[Metadata_2_setUp$`clinical outcome` == "Remission"] <- "1"
Metadata_2_setUp$`remission (clinical)`[!Metadata_2_setUp$`clinical outcome` == "Remission"] <- "0"

Metadata_2_setUp$`non-remission (clinical)` <- 1
Metadata_2_setUp$`non-remission (clinical)`[Metadata_2_setUp$`clinical outcome` == "Non-Remission"] <- "1"
Metadata_2_setUp$`non-remission (clinical)`[!Metadata_2_setUp$`clinical outcome` == "Non-Remission"] <- "0"
#View(Metadata_2_setUp)
####################################
# 7.3.1.2 Integers
####################################
# Make sure that this is not a tibble DF, but simply a DF
# and make sure the numbers are integers
Metadata_2_set_analysis <- as.data.frame(Metadata_2_setUp)
Metadata_2_set_analysis$`crohn's disease` <- as.integer(Metadata_2_set_analysis$`crohn's disease`)
Metadata_2_set_analysis$`ulcerative colitis` <- as.integer(Metadata_2_set_analysis$`ulcerative colitis`)
Metadata_2_set_analysis$infliximab <- as.integer(Metadata_2_set_analysis$infliximab)
Metadata_2_set_analysis$adalimumab <- as.integer(Metadata_2_set_analysis$adalimumab)
Metadata_2_set_analysis$vedolizumab <- as.integer(Metadata_2_set_analysis$vedolizumab)
Metadata_2_set_analysis$ustekinumab <- as.integer(Metadata_2_set_analysis$ustekinumab)

Metadata_2_set_analysis$`remission (endoscopical)` <- as.integer(Metadata_2_set_analysis$`remission (endoscopical)`)
Metadata_2_set_analysis$`non-remission (endoscopical)` <- as.integer(Metadata_2_set_analysis$`non-remission (endoscopical)`)

Metadata_2_set_analysis$`remission (biomarker)` <- as.integer(Metadata_2_set_analysis$`remission (biomarker)`)
Metadata_2_set_analysis$`non-remission (biomarker)` <- as.integer(Metadata_2_set_analysis$`non-remission (biomarker)`)

Metadata_2_set_analysis$`remission (clinical)` <- as.integer(Metadata_2_set_analysis$`remission (clinical)`)
Metadata_2_set_analysis$`non-remission (clinical)` <- as.integer(Metadata_2_set_analysis$`non-remission (clinical)`)
colnames(Metadata_2_set_analysis)
str(Metadata_2_set_analysis)

#View(Metadata_2_set_analysis)
####################################
# 7.3.1.3 UpsetR plot
####################################
# Choose Sets: select all of them & orden by frequency
upset(Metadata_2_set_analysis, sets = c("crohn's disease", "ulcerative colitis", "infliximab", "adalimumab", "vedolizumab", "ustekinumab","endoscopical-remission","endoscopical-non-remission"), mb.ratio = c(0.55, 0.45), order.by = "freq")

# orden by degree
upset(Metadata_2_set_analysis, sets = c("crohn's disease", "ulcerative colitis", "infliximab", "adalimumab", "vedolizumab", "ustekinumab","remission","non-remission"), mb.ratio = c(0.55, 0.45), order.by = "degree")

# orden by degree and frequency
upset(Metadata_2_set_analysis, sets = c("crohn's disease", "ulcerative colitis", "infliximab", "adalimumab", "vedolizumab", "ustekinumab","remission","non-remission"), mb.ratio = c(0.55, 0.45), order.by = c("degree", "freq"))

# group by sets
upset(Metadata_2_set_analysis, sets = c("crohn's disease", "ulcerative colitis", "infliximab", "adalimumab", "vedolizumab", "ustekinumab","remission","non-remission"), group.by = "sets")

# Orden by degree (diagnosis) & make more esthethic
upset(Metadata_2_set_analysis, sets = c("crohn's disease", "ulcerative colitis", "infliximab", "adalimumab", "vedolizumab", "ustekinumab","remission","non-remission"), mb.ratio = c(0.7, 0.3), order.by = "degree",
      number.angles = 30, point.size = 2, line.size = 1,
      mainbar.y.label = "endoscopical outcome", sets.x.label = "patients (#)"
)
grid.text(
  "IBD patients",
  x = 0.90,
  y = 0.05,
  gp = gpar(
    fontsize = 10,
    fontface = 2
  )
)


make_comb_mat(Metadata_2_set_analysis)
####################################
# 7.3.2 Complex heatmap package
####################################
# 7.3.2.1 Create combination matrix as alternative for ComplexHeatmaps
####################################
# https://jokergoo.github.io/ComplexHeatmap-reference/book/upset-plot.html
# A combination matrix is much more flexible
# Here we re-implemented UpSet plots with the ComplexHeatmap package with some improvements.

# session_info()
# install.packages("BiocManager")
# library(BiocManager)
# BiocManager::install("ComplexHeatmap")
# library(ComplexHeatmap)

# lt <- list(diagnosis = c("crohn's disease", "ulcerative colitis"),
#     therapy = c("infliximab", "adalimumab", "vedolizumab", "ustekinumab"),
#    `endoscopical outcome` = c("remission", "non-remission"))

# list_to_matrix(lt)
colnames(Metadata_2_set_analysis)
Metadata_2_set_analysis$diagnosis <- NULL
Metadata_2_set_analysis$therapy <- NULL
Metadata_2_set_analysis$`endoscopical outcome` <- NULL
Metadata_2_set_analysis$`clinical outcome` <- NULL
Metadata_2_set_analysis$`biomarker outcome` <- NULL
colnames(Metadata_2_set_analysis)
#View(Metadata_2_set_analysis)

Metadata_2_set_analysis_complex_heatmap_plot <- Metadata_2_set_analysis
Upset_Complex_heatmap_plot <- make_comb_mat(Metadata_2_set_analysis_complex_heatmap_plot) # default mode = "distinct"
Upset_Complex_heatmap_plot
str(Upset_Complex_heatmap_plot)
####################################
# 7.3.2.2 UpsetR Complex Heatmap plot
####################################
UpSet(Upset_Complex_heatmap_plot)
UpSet(t(Upset_Complex_heatmap_plot))

UpSet(Upset_Complex_heatmap_plot, set_order = c("crohn's disease", "ulcerative colitis","remission" , "non-remission" ,"infliximab", "adalimumab", "vedolizumab", "ustekinumab"), comb_order = order(comb_size(test_)))

UpSet(Upset_Complex_heatmap_plot, set_order = c("crohn's disease", "ulcerative colitis","remission" , "non-remission" ,"infliximab", "adalimumab", "vedolizumab", "ustekinumab"), 
      pt_size = unit(3, "mm"),
      lwd = 2,
      row_names_side = "left",
      comb_col = c(rep(1,7), rep(8,12),2))

UpSet(Upset_Complex_heatmap_plot, 
      pt_size = unit(3, "mm"),
      lwd = 2,
      row_names_side = "left",
      comb_col = c(rep(1,7), rep(8,12),2))

# In this case you can see that CD patients are th biggest groups with most connections to ustekinumab treatment under non-remission.
# Second place are CD that undergo vedolizumab treatment and are in under remission.
# Third, UC patients in remission undergoing vedolizumab treatment.

# 162 patients with data available for patients
# 60 UC 
# 102 CD
####################################
# 7.3.2.3 UpsetR Complex Heatmap plot per diagnosis
####################################
# https://jokergoo.github.io/ComplexHeatmap-reference/book/upset-plot.html

## What is the choice between: intersect, distinct, union
## (i) distinct mode:  mode: 1 means in that set and 0 means not in that set.
## (ii )intersect mode: mode: 1 means in that set and 0 is not taken into account.
## (iii) union mode: 1 means in that set and 0 is not taken into account.
## WE WILL USE DISTINCT MODE ALWAYS!

# Drop one diagnosis (UC or CD)
Metadata_2_set_analysis_complex_heatmap_plot_CD <- Metadata_2_set_analysis_complex_heatmap_plot[Metadata_2_set_analysis_complex_heatmap_plot$`crohn's disease` == "1",]
Metadata_2_set_analysis_complex_heatmap_plot_UC <- Metadata_2_set_analysis_complex_heatmap_plot[Metadata_2_set_analysis_complex_heatmap_plot$`ulcerative colitis` == "1",]
colnames(Metadata_2_set_analysis_complex_heatmap_plot_CD)

# Take remissions only
Metadata_2_set_analysis_complex_heatmap_plot_CD$`crohn's disease` <- NULL
Metadata_2_set_analysis_complex_heatmap_plot_CD$`ulcerative colitis` <- NULL
Metadata_2_set_analysis_complex_heatmap_plot_CD$infliximab <- NULL
Metadata_2_set_analysis_complex_heatmap_plot_CD$adalimumab <- NULL
Metadata_2_set_analysis_complex_heatmap_plot_CD$vedolizumab <- NULL
Metadata_2_set_analysis_complex_heatmap_plot_CD$ustekinumab <- NULL

Metadata_2_set_analysis_complex_heatmap_plot_UC$`ulcerative colitis` <- NULL
Metadata_2_set_analysis_complex_heatmap_plot_UC$`crohn's disease` <- NULL
Metadata_2_set_analysis_complex_heatmap_plot_UC$infliximab <- NULL
Metadata_2_set_analysis_complex_heatmap_plot_UC$adalimumab <- NULL
Metadata_2_set_analysis_complex_heatmap_plot_UC$vedolizumab <- NULL
Metadata_2_set_analysis_complex_heatmap_plot_UC$ustekinumab <- NULL
colnames(Metadata_2_set_analysis_complex_heatmap_plot_UC)

# Create a combination matrix
Upset_Complex_heatmap_plot_UC <- make_comb_mat(Metadata_2_set_analysis_complex_heatmap_plot_UC) 
Upset_Complex_heatmap_plot_CD <- make_comb_mat(Metadata_2_set_analysis_complex_heatmap_plot_CD) 
Upset_Complex_heatmap_plot_UC
Upset_Complex_heatmap_plot_CD

# Create multiple UpsetPlot
## UC
Upset_Complex_heatmap_UC <- UpSet((Upset_Complex_heatmap_plot_UC), 
      column_title = "Colitis", 
      comb_order = order(comb_size(Upset_Complex_heatmap_plot_UC)),
      set_order = c("remission (endoscopical)","remission (clinical)","remission (biomarker)" ,"non-remission (endoscopical)", "non-remission (clinical)","non-remission (biomarker)"),
      top_annotation = upset_top_annotation(Upset_Complex_heatmap_plot_UC, gp = gpar(fill = "#D9846B"), height = unit(5, "cm")))
     
Upset_Complex_heatmap_UC


## CD
Upset_Complex_heatmap_CD <- UpSet((Upset_Complex_heatmap_plot_CD), 
      column_title = "Crohn's", 
      comb_order = order(comb_size(Upset_Complex_heatmap_plot_CD)),
      set_order = c("remission (endoscopical)","remission (clinical)","remission (biomarker)" ,"non-remission (endoscopical)", "non-remission (clinical)","non-remission (biomarker)"),
      top_annotation = upset_top_annotation(Upset_Complex_heatmap_plot_CD, gp = gpar(fill = "#8080B4"), height = unit(5, "cm")
      ))

Upset_Complex_heatmap_CD

# Combine
Upset_Complex_heatmap_combi <- Upset_Complex_heatmap_UC + Upset_Complex_heatmap_CD
Upset_Complex_heatmap_combi # finish details in illustrator

# With complex heatmap lot more options to attach to figures.
####################################
# X.Number of treamtent per biological group
####################################
view(Metadata)
Metadata$biological_classes_before <- as.character(Metadata$biological_classes_before)
Metadata$biologicals_before <- as.character(Metadata$biologicals_before)

str(Metadata$biological_classes_before)

esquisser()

Metadata %>%
 filter(Week >= 0L & Week <= 0L) %>%
 ggplot() +
 aes(x = biologicals_before, fill = Therapy_1) +
 geom_bar() +
scale_shape_manual(values = c(0, 16, 3, 4)) +
 theme_classic()

Metadata %>%
 filter(Week >= 0L & Week <= 0L) %>%
 ggplot() +
 aes(x = biological_classes_before, fill = Therapy_1) +
 geom_bar() +
 scale_fill_hue(direction = 1) +
 theme_classic()

table(Metadata$biologicals_before[Metadata$Week == "0"])
table(Metadata$biological_classes_before[Metadata$Week == "0"])
######### !!!!!!!!!!!!!!! #########
# SAVE Global Environment #
######### !!!!!!!!!!!!!!! #########