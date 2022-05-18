####################################
# SCRIPT 5: EUKARYOTIC VIRUSES
####################################
# Before starting this script
####################################
# ------------------> Load the "Global environment" output of script 3 "Rarefaction.R"

# ------------------> You can show presence/absence table without the need for rarefaction based on "Mastertable"
# ------------------> If you want to show difference in abundances, richness or diversity you need to work with Rarefied tables.

# ------------------> Since eukaryotic viruses are an absolute minority in the total viral fraction in our IBD sample we take a subset of it to have a better look at them.
# ------------------> Therefore, we will only represent abundances & presence/absence here.

# PS: if not done yet, make a pie-chart (~ Metadata script) to show the difference between the relative abundance of phages and eukaryotic viruses (have to be done on rarefied table)
####################################
# 0. Packages: Install packages and load them wherever needed
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
# 1. Recapitulation of mastertables
####################################
## Take into account that most of diversity analysis already correct for sequencing depth (alpha-diversity), so you have to feed it the unrarefied completeness data

# global mastertable
Mastertable
Mastertable_unrarefied
Mastertable_rarefied

# viral mastertable
# unrarefied
Mastertable_viral_unrarefied # PRESENCE/ABSENCE, ALPHA-DIVERSITY (Shannon diversity)

# rarefied
Mastertable_viral_rarefied # RELATIVE ABUNDANCES, RICHNESS,
####################################
# 2. Evaluate Taxonomical annotation
####################################
# Since eukaryotic viruses are thought to be quite well annotated we will manually check the %ANI, %coverage of Blastn, Diamond.
# We will prioritize always Diamond over Blastn over CAT. Because Diamond thrives a NR-protein database.
# We will subset all eukaryotic viruses out of the tables above.
# Next you have to manually check all the eukaryotic viruses with < 70% ANI; 70% coverage


View(Mastertable_rarefied)

####################################
# 2.2 Rarefied viral mastertable
####################################
Mastertable_viral_rarefied_euk_virus <- Mastertable_viral_rarefied[Mastertable_viral_rarefied$Final_viral == "eukaryotic virus" & Mastertable_viral_rarefied$Final_coverage > 10,]
####################################
# 2.3 Mastertable: evaluate (~ use for figures)
####################################
# Now open every table en look at that contigs lower than 70%ANI; 70% coverage to see if the taxonomical alignments are thrusthworthy or not.
# If not remove them from our table, it could be that you have to look for the specific family, genus what the cutoff are.
# It could be that the one with reasonable but crappy alignment are novel viruses, so you know that.
# If the %ANI, %cov is too crap you can best remove this, because it would be wrong.

# unrarefied mastertable (> 0.1 AS)
#View(Mastertable_viral_unrarefied_euk_virus)        # PRESENCE/ABSENCE, ALPHA-DIVERSITY
table(Mastertable_viral_unrarefied_euk_virus$Final_family)
sum(Mastertable_viral_unrarefied_euk_virus$Totalnumberofreads)

# rarefied mastertable (> 0.1 AS)
#View(Mastertable_viral_rarefied_euk_virus)              # RELATIVE ABUNDANCES, RICHNESS
table(Mastertable_viral_rarefied_euk_virus$Final_family)
sum(Mastertable_viral_rarefied_euk_virus$Totalnumberofreads)
####################################
# 2.4 Calculate AS (AF*AAI/ANI) (if you did this before skip it)
####################################
# 2.4.1 Unrarefied
####################################
# Add AS
Mastertable_viral_unrarefied_euk_virus$Blastn_AS <- ((Mastertable_viral_unrarefied_euk_virus$Blastn_ANI/100)*(Mastertable_viral_unrarefied_euk_virus$Blastn_coverage_relative/100))
Mastertable_viral_unrarefied_euk_virus$Blastn_AS[!Mastertable_viral_unrarefied_euk_virus$Blastn_viral == "Blastn_eukaryotic_virus"] <- 0

Mastertable_viral_unrarefied_euk_virus$Diamond_AS <- ((Mastertable_viral_unrarefied_euk_virus$Diamond_ANI/100)*(Mastertable_viral_unrarefied_euk_virus$Diamond_coverage_relative/100))
Mastertable_viral_unrarefied_euk_virus$Diamond_AS[!Mastertable_viral_unrarefied_euk_virus$Diamond_viral == "Diamond_eukaryotic_virus"] <- 0

# Annotation tools to use
Mastertable_viral_unrarefied_euk_virus$Best_AS[Mastertable_viral_unrarefied_euk_virus$Diamond_AS > Mastertable_viral_unrarefied_euk_virus$Blastn_AS] <- "Diamond"
Mastertable_viral_unrarefied_euk_virus$Best_AS[Mastertable_viral_unrarefied_euk_virus$Diamond_AS < Mastertable_viral_unrarefied_euk_virus$Blastn_AS] <- "Blastn"
Mastertable_viral_unrarefied_euk_virus$Best_AS[Mastertable_viral_unrarefied_euk_virus$Diamond_AS < 0.1 & Mastertable_viral_unrarefied_euk_virus$Blastn_AS < 0.1] <- "Remove"

# Re-annotate: Blastn
Mastertable_viral_unrarefied_euk_virus$Final_phylum[Mastertable_viral_unrarefied_euk_virus$Best_AS == "Blastn"] <- Mastertable_viral_unrarefied_euk_virus$Blastn_phylum[Mastertable_viral_unrarefied_euk_virus$Best_AS == "Blastn"]
Mastertable_viral_unrarefied_euk_virus$Final_class[Mastertable_viral_unrarefied_euk_virus$Best_AS == "Blastn"] <- Mastertable_viral_unrarefied_euk_virus$Blastn_class[Mastertable_viral_unrarefied_euk_virus$Best_AS == "Blastn"]
Mastertable_viral_unrarefied_euk_virus$Final_order[Mastertable_viral_unrarefied_euk_virus$Best_AS == "Blastn"] <- Mastertable_viral_unrarefied_euk_virus$Blastn_order[Mastertable_viral_unrarefied_euk_virus$Best_AS == "Blastn"]
Mastertable_viral_unrarefied_euk_virus$Final_family[Mastertable_viral_unrarefied_euk_virus$Best_AS == "Blastn"] <- Mastertable_viral_unrarefied_euk_virus$Blastn_family[Mastertable_viral_unrarefied_euk_virus$Best_AS == "Blastn"]
Mastertable_viral_unrarefied_euk_virus$Final_subfamily[Mastertable_viral_unrarefied_euk_virus$Best_AS == "Blastn"] <- Mastertable_viral_unrarefied_euk_virus$Blastn_subfamily[Mastertable_viral_unrarefied_euk_virus$Best_AS == "Blastn"]
Mastertable_viral_unrarefied_euk_virus$Final_genus[Mastertable_viral_unrarefied_euk_virus$Best_AS == "Blastn"] <- Mastertable_viral_unrarefied_euk_virus$Blastn_genus[Mastertable_viral_unrarefied_euk_virus$Best_AS == "Blastn"]
Mastertable_viral_unrarefied_euk_virus$Final_species[Mastertable_viral_unrarefied_euk_virus$Best_AS == "Blastn"] <- Mastertable_viral_unrarefied_euk_virus$Blastn_species[Mastertable_viral_unrarefied_euk_virus$Best_AS == "Blastn"]

# Re-annotate: Diamond
Mastertable_viral_unrarefied_euk_virus$Final_phylum[Mastertable_viral_unrarefied_euk_virus$Best_AS == "Diamond"] <- Mastertable_viral_unrarefied_euk_virus$Diamond_phylum[Mastertable_viral_unrarefied_euk_virus$Best_AS == "Diamond"] 
Mastertable_viral_unrarefied_euk_virus$Final_class[Mastertable_viral_unrarefied_euk_virus$Best_AS == "Diamond"] <- Mastertable_viral_unrarefied_euk_virus$Diamond_class[Mastertable_viral_unrarefied_euk_virus$Best_AS == "Diamond"] 
Mastertable_viral_unrarefied_euk_virus$Final_order[Mastertable_viral_unrarefied_euk_virus$Best_AS == "Diamond"] <- Mastertable_viral_unrarefied_euk_virus$Diamond_order[Mastertable_viral_unrarefied_euk_virus$Best_AS == "Diamond"] 
Mastertable_viral_unrarefied_euk_virus$Final_family[Mastertable_viral_unrarefied_euk_virus$Best_AS == "Diamond"] <- Mastertable_viral_unrarefied_euk_virus$Diamond_family[Mastertable_viral_unrarefied_euk_virus$Best_AS == "Diamond"] 
Mastertable_viral_unrarefied_euk_virus$Final_subfamily[Mastertable_viral_unrarefied_euk_virus$Best_AS == "Diamond"] <- Mastertable_viral_unrarefied_euk_virus$Diamond_subfamily[Mastertable_viral_unrarefied_euk_virus$Best_AS == "Diamond"] 
Mastertable_viral_unrarefied_euk_virus$Final_genus[Mastertable_viral_unrarefied_euk_virus$Best_AS == "Diamond"] <- Mastertable_viral_unrarefied_euk_virus$Diamond_genus[Mastertable_viral_unrarefied_euk_virus$Best_AS == "Diamond"] 
Mastertable_viral_unrarefied_euk_virus$Final_species[Mastertable_viral_unrarefied_euk_virus$Best_AS == "Diamond"] <- Mastertable_viral_unrarefied_euk_virus$Diamond_species[Mastertable_viral_unrarefied_euk_virus$Best_AS == "Diamond"] 

# Remove bad alignments
Mastertable_viral_unrarefied_euk_virus$Final_viral[Mastertable_viral_unrarefied_euk_virus$Best_AS == "Remove"] <- 1
Mastertable_viral_unrarefied_euk_virus$Final_superkingdom[Mastertable_viral_unrarefied_euk_virus$Best_AS == "Remove"] <- "Unannotated"
Mastertable_viral_unrarefied_euk_virus <- Mastertable_viral_unrarefied_euk_virus[Mastertable_viral_unrarefied_euk_virus$Final_viral == "eukaryotic virus",]

# Reads and contigs
table(Mastertable_viral_unrarefied_euk_virus$Final_viral)
nrow(Mastertable_viral_unrarefied_euk_virus[Mastertable_viral_unrarefied_euk_virus$Final_viral == "eukaryotic virus"])
sum(Mastertable_viral_unrarefied_euk_virus$Totalnumberofreads[Mastertable_viral_unrarefied_euk_virus$Final_viral == "eukaryotic virus"])
####################################
# 2.4.2 Rarefied
####################################
# Add AS
Mastertable_viral_rarefied_euk_virus$Blastn_AS <- ((Mastertable_viral_rarefied_euk_virus$Blastn_ANI/100)*(Mastertable_viral_rarefied_euk_virus$Blastn_coverage_relative/100))
Mastertable_viral_rarefied_euk_virus$Blastn_AS[!Mastertable_viral_rarefied_euk_virus$Blastn_viral == "Blastn_eukaryotic_virus"] <- 0

Mastertable_viral_rarefied_euk_virus$Diamond_AS <- ((Mastertable_viral_rarefied_euk_virus$Diamond_ANI/100)*(Mastertable_viral_rarefied_euk_virus$Diamond_coverage_relative/100))
Mastertable_viral_rarefied_euk_virus$Diamond_AS[!Mastertable_viral_rarefied_euk_virus$Diamond_viral == "Diamond_eukaryotic_virus"] <- 0

# Annotation tools to use
Mastertable_viral_rarefied_euk_virus$Best_AS[Mastertable_viral_rarefied_euk_virus$Diamond_AS > Mastertable_viral_rarefied_euk_virus$Blastn_AS] <- "Diamond"
Mastertable_viral_rarefied_euk_virus$Best_AS[Mastertable_viral_rarefied_euk_virus$Diamond_AS < Mastertable_viral_rarefied_euk_virus$Blastn_AS] <- "Blastn"
Mastertable_viral_rarefied_euk_virus$Best_AS[Mastertable_viral_rarefied_euk_virus$Diamond_AS < 0.1 & Mastertable_viral_rarefied_euk_virus$Blastn_AS < 0.1] <- "Remove"

# Re-annotate: Blastn
Mastertable_viral_rarefied_euk_virus$Final_phylum[Mastertable_viral_rarefied_euk_virus$Best_AS == "Blastn"] <- Mastertable_viral_rarefied_euk_virus$Blastn_phylum[Mastertable_viral_rarefied_euk_virus$Best_AS == "Blastn"]
Mastertable_viral_rarefied_euk_virus$Final_class[Mastertable_viral_rarefied_euk_virus$Best_AS == "Blastn"] <- Mastertable_viral_rarefied_euk_virus$Blastn_class[Mastertable_viral_rarefied_euk_virus$Best_AS == "Blastn"]
Mastertable_viral_rarefied_euk_virus$Final_order[Mastertable_viral_rarefied_euk_virus$Best_AS == "Blastn"] <- Mastertable_viral_rarefied_euk_virus$Blastn_order[Mastertable_viral_rarefied_euk_virus$Best_AS == "Blastn"]
Mastertable_viral_rarefied_euk_virus$Final_family[Mastertable_viral_rarefied_euk_virus$Best_AS == "Blastn"] <- Mastertable_viral_rarefied_euk_virus$Blastn_family[Mastertable_viral_rarefied_euk_virus$Best_AS == "Blastn"]
Mastertable_viral_rarefied_euk_virus$Final_subfamily[Mastertable_viral_rarefied_euk_virus$Best_AS == "Blastn"] <- Mastertable_viral_rarefied_euk_virus$Blastn_subfamily[Mastertable_viral_rarefied_euk_virus$Best_AS == "Blastn"]
Mastertable_viral_rarefied_euk_virus$Final_genus[Mastertable_viral_rarefied_euk_virus$Best_AS == "Blastn"] <- Mastertable_viral_rarefied_euk_virus$Blastn_genus[Mastertable_viral_rarefied_euk_virus$Best_AS == "Blastn"]
Mastertable_viral_rarefied_euk_virus$Final_species[Mastertable_viral_rarefied_euk_virus$Best_AS == "Blastn"] <- Mastertable_viral_rarefied_euk_virus$Blastn_species[Mastertable_viral_rarefied_euk_virus$Best_AS == "Blastn"]
Mastertable_viral_rarefied_euk_virus$Final_ANI[Mastertable_viral_rarefied_euk_virus$Best_AS == "Blastn"] <- Mastertable_viral_rarefied_euk_virus$Blastn_ANI[Mastertable_viral_rarefied_euk_virus$Best_AS == "Blastn"]
Mastertable_viral_rarefied_euk_virus$Final_coverage[Mastertable_viral_rarefied_euk_virus$Best_AS == "Blastn"] <- Mastertable_viral_rarefied_euk_virus$Blastn_coverage_relative[Mastertable_viral_rarefied_euk_virus$Best_AS == "Blastn"]

# Re-annotate: Diamond
Mastertable_viral_rarefied_euk_virus$Final_phylum[Mastertable_viral_rarefied_euk_virus$Best_AS == "Diamond"] <- Mastertable_viral_rarefied_euk_virus$Diamond_phylum[Mastertable_viral_rarefied_euk_virus$Best_AS == "Diamond"] 
Mastertable_viral_rarefied_euk_virus$Final_class[Mastertable_viral_rarefied_euk_virus$Best_AS == "Diamond"] <- Mastertable_viral_rarefied_euk_virus$Diamond_class[Mastertable_viral_rarefied_euk_virus$Best_AS == "Diamond"] 
Mastertable_viral_rarefied_euk_virus$Final_order[Mastertable_viral_rarefied_euk_virus$Best_AS == "Diamond"] <- Mastertable_viral_rarefied_euk_virus$Diamond_order[Mastertable_viral_rarefied_euk_virus$Best_AS == "Diamond"] 
Mastertable_viral_rarefied_euk_virus$Final_family[Mastertable_viral_rarefied_euk_virus$Best_AS == "Diamond"] <- Mastertable_viral_rarefied_euk_virus$Diamond_family[Mastertable_viral_rarefied_euk_virus$Best_AS == "Diamond"] 
Mastertable_viral_rarefied_euk_virus$Final_subfamily[Mastertable_viral_rarefied_euk_virus$Best_AS == "Diamond"] <- Mastertable_viral_rarefied_euk_virus$Diamond_subfamily[Mastertable_viral_rarefied_euk_virus$Best_AS == "Diamond"] 
Mastertable_viral_rarefied_euk_virus$Final_genus[Mastertable_viral_rarefied_euk_virus$Best_AS == "Diamond"] <- Mastertable_viral_rarefied_euk_virus$Diamond_genus[Mastertable_viral_rarefied_euk_virus$Best_AS == "Diamond"] 
Mastertable_viral_rarefied_euk_virus$Final_species[Mastertable_viral_rarefied_euk_virus$Best_AS == "Diamond"] <- Mastertable_viral_rarefied_euk_virus$Diamond_species[Mastertable_viral_rarefied_euk_virus$Best_AS == "Diamond"] 
Mastertable_viral_rarefied_euk_virus$Final_ANI[Mastertable_viral_rarefied_euk_virus$Best_AS == "Diamond"] <- Mastertable_viral_rarefied_euk_virus$Diamond_ANI[Mastertable_viral_rarefied_euk_virus$Best_AS == "Diamond"]
Mastertable_viral_rarefied_euk_virus$Final_coverage[Mastertable_viral_rarefied_euk_virus$Best_AS == "Diamond"] <- Mastertable_viral_rarefied_euk_virus$Diamond_coverage_relative[Mastertable_viral_rarefied_euk_virus$Best_AS == "Diamond"]

# Remove bad alignments
Mastertable_viral_rarefied_euk_virus$Final_viral[Mastertable_viral_rarefied_euk_virus$Best_AS == "Remove"] <- 1
Mastertable_viral_rarefied_euk_virus$Final_superkingdom[Mastertable_viral_rarefied_euk_virus$Best_AS == "Remove"] <- "Unannotated"
Mastertable_viral_rarefied_euk_virus <- Mastertable_viral_rarefied_euk_virus[Mastertable_viral_rarefied_euk_virus$Final_viral == "eukaryotic virus",]

# Reads and contigs
table(Mastertable_viral_rarefied_euk_virus$Final_viral)
nrow(Mastertable_viral_rarefied_euk_virus[Mastertable_viral_rarefied_euk_virus$Final_viral == "eukaryotic virus"])
sum(Mastertable_viral_rarefied_euk_virus$Totalnumberofreads[Mastertable_viral_rarefied_euk_virus$Final_viral == "eukaryotic virus"])
mean(Mastertable_viral_rarefied_euk_virus$completeness[Mastertable_viral_rarefied_euk_virus$Final_viral == "eukaryotic virus"])
mean(Mastertable_viral_rarefied_euk_virus$length[Mastertable_viral_rarefied_euk_virus$Final_viral == "eukaryotic virus"])
range(Mastertable_viral_rarefied_euk_virus$length[Mastertable_viral_rarefied_euk_virus$Final_viral == "eukaryotic virus"])
table(Mastertable_viral_rarefied_euk_virus$Final_family)
####################################
# 2.5 Manual curation
####################################
# Unrarefied
#view(Mastertable_viral_unrarefied_euk_virus)

# Genus
Mastertable_viral_unrarefied_euk_virus$Final_genus[Mastertable_viral_unrarefied_euk_virus$Final_species == "Shallot virus S"] <- "Carlavirus"
Mastertable_viral_unrarefied_euk_virus$Final_genus[Mastertable_viral_unrarefied_euk_virus$Final_phylum == "Cressdnaviricota" & Mastertable_viral_unrarefied_euk_virus$Final_family == "Circoviridae"] <- "unclassified Circoviridae"
Mastertable_viral_unrarefied_euk_virus$Final_genus[Mastertable_viral_unrarefied_euk_virus$Final_phylum == "Cressdnaviricota" & Mastertable_viral_unrarefied_euk_virus$Final_family == "Unannotated"] <- "unclassified CRESS"
Mastertable_viral_unrarefied_euk_virus$Final_genus[Mastertable_viral_unrarefied_euk_virus$Final_species == "Circovirus VS6600032"] <- "Circovirus"
Mastertable_viral_unrarefied_euk_virus$Final_genus[Mastertable_viral_unrarefied_euk_virus$Final_species == "uncultured densovirus"] <- "unclassified Parvoviridae"
Mastertable_viral_unrarefied_euk_virus$Final_genus[Mastertable_viral_unrarefied_euk_virus$Final_species == "Anelloviridae sp."] <- "unclassified Anelloviridae"
Mastertable_viral_unrarefied_euk_virus$Final_genus[Mastertable_viral_unrarefied_euk_virus$Final_species == "Paradise bay virus"] <- "unclassified Tombusviridae"

# Family
Mastertable_viral_unrarefied_euk_virus$Final_family[Mastertable_viral_unrarefied_euk_virus$Final_family == "Unannotated"] <- "unclassified CRESS"

# Rarefied 
#View(Mastertable_viral_rarefied_euk_virus)

# Genus
Mastertable_viral_rarefied_euk_virus$Final_genus[Mastertable_viral_rarefied_euk_virus$Final_species == "Shallot virus S"] <- "Carlavirus"
Mastertable_viral_rarefied_euk_virus$Final_genus[Mastertable_viral_rarefied_euk_virus$Final_phylum == "Cressdnaviricota" & Mastertable_viral_rarefied_euk_virus$Final_family == "Circoviridae"] <- "unclassified Circoviridae"
Mastertable_viral_rarefied_euk_virus$Final_genus[Mastertable_viral_rarefied_euk_virus$Final_phylum == "Cressdnaviricota" & Mastertable_viral_rarefied_euk_virus$Final_family == "Unannotated"] <- "unclassified CRESS"
Mastertable_viral_rarefied_euk_virus$Final_genus[Mastertable_viral_rarefied_euk_virus$Final_species == "Circovirus VS6600032"] <- "Circovirus"
Mastertable_viral_rarefied_euk_virus$Final_genus[Mastertable_viral_rarefied_euk_virus$Final_species == "uncultured densovirus"] <- "unclassified Parvoviridae"
Mastertable_viral_rarefied_euk_virus$Final_genus[Mastertable_viral_rarefied_euk_virus$Final_species == "Anelloviridae sp."] <- "unclassified Anelloviridae"
Mastertable_viral_rarefied_euk_virus$Final_genus[Mastertable_viral_rarefied_euk_virus$Final_species == "Paradise bay virus"] <- "unclassified Tombusviridae"

# Family
Mastertable_viral_rarefied_euk_virus$Final_family[Mastertable_viral_rarefied_euk_virus$Final_family == "Unannotated"] <- "unclassified CRESS"

#View(Mastertable_viral_rarefied_euk_virus)
####################################
# 3. Create a list of all eukaryotic virus
####################################
# 3.1 List of all eukaryotic viruses (~ family level)
####################################
Mastertable_viral_rarefied_euk_virus
Mastertable_viral_unrarefied_euk_virus

unique(Mastertable_viral_unrarefied_euk_virus$Final_family)
# ---------> Nr. of 18 unique eukaryotic viral families

unique(Mastertable_viral_unrarefied_euk_virus$Final_genus)
# ---------> Nr. of 27 unique eukaryotic viral genera

####################################
# 3.2 write down hosts
####################################
## We divide the above list into 3 sub-list: (1) animal-infecting viruses, (2) plant/fungi-infecting virus, (3) small circular viruses

#View(Mastertable_viral_unrarefied_euk_virus)
#View(Mastertable_viral_rarefied_euk_virus)

# 1) Animal-infecting viruses
# Family
Caliciviridae
Papillomaviridae
Picornaviridae
Parvoviridae

# Genus
Alphapapillomavirus
Cosavirus
Enterovirus
Norovirus
Sapovirus
unclassified Astroviridae
unclassified Parvoviridae

# 2) Plant/fungi-infecting viruses
# Family
Alphaflexiviridae
Betaflexiviridae
Chrysoviridae
Closteroviridae
Endornaviridae
Secoviridae
Virgaviridae

# Genus
Alphachrysovirus
Alphaendornavirus
Carlavirus
Crinivirus
Nepovirus
Potexvirus
Tobamovirus
Trichovirus

# 3) Small circular viruses
# Family
Anelloviridae
Circoviridae
Cressdnaviricota
Genomoviridae
Geminiviridae
Smacoviridae

# Genus
Alphatorquevirus
Begomovirus
Betatorquevirus
Circovirus
Gammatorquevirus
Gemycircularvirus
Gemykrogvirus
Gyrovirus
Huchismacovirus
unclassified Anelloviridae
unclassified Circoviridae
unclassified CRESS
####################################
# 3.3 Add hosts to mastertable
####################################
# 3.3.1 Unrarefied mastertable 
####################################
## Adapt the lines here according to what you saw qua hosts before.

# (A) Create a new table to add 'group host annotation'
Mastertable_viral_unrarefied_euk_virus_host <- Mastertable_viral_unrarefied_euk_virus
unique(Mastertable_viral_unrarefied_euk_virus$Final_family)

# (B) Add 'group host annotation' in a new column for plant and fungal viruses
Mastertable_viral_unrarefied_euk_virus$Host <- Mastertable_viral_unrarefied_euk_virus$Final_family
Mastertable_viral_unrarefied_euk_virus$Host[Mastertable_viral_unrarefied_euk_virus$Final_family == "Alphaflexiviridae"] <- "plant & fungal viruses"
Mastertable_viral_unrarefied_euk_virus$Host[Mastertable_viral_unrarefied_euk_virus$Final_family == "Closteroviridae"] <- "plant & fungal viruses"
Mastertable_viral_unrarefied_euk_virus$Host[Mastertable_viral_unrarefied_euk_virus$Final_family == "Betaflexiviridae"] <- "plant & fungal viruses"
Mastertable_viral_unrarefied_euk_virus$Host[Mastertable_viral_unrarefied_euk_virus$Final_family == "Chrysoviridae"] <- "plant & fungal viruses"
Mastertable_viral_unrarefied_euk_virus$Host[Mastertable_viral_unrarefied_euk_virus$Final_family == "Virgaviridae"] <- "plant & fungal viruses"
Mastertable_viral_unrarefied_euk_virus$Host[Mastertable_viral_unrarefied_euk_virus$Final_family == "Endornaviridae"] <- "plant & fungal viruses"
Mastertable_viral_unrarefied_euk_virus$Host[Mastertable_viral_unrarefied_euk_virus$Final_family == "Secoviridae"] <- "plant & fungal viruses"

# (C) Add 'group host annotation' in a new column for animal viruses
Mastertable_viral_unrarefied_euk_virus$Host[Mastertable_viral_unrarefied_euk_virus$Final_family == "Caliciviridae"] <- "animal viruses"
Mastertable_viral_unrarefied_euk_virus$Host[Mastertable_viral_unrarefied_euk_virus$Final_family == "Parvoviridae"] <- "animal viruses"
Mastertable_viral_unrarefied_euk_virus$Host[Mastertable_viral_unrarefied_euk_virus$Final_family == "Picornaviridae"] <- "animal viruses"
Mastertable_viral_unrarefied_euk_virus$Host[Mastertable_viral_unrarefied_euk_virus$Final_family == "Papillomaviridae"] <- "animal viruses"

# (D) Add 'group host annotation' in a new column for small circular and unclassified viruses
Mastertable_viral_unrarefied_euk_virus$Host[Mastertable_viral_unrarefied_euk_virus$Final_family == "unclassified CRESS"] <- "small circular viruses"
Mastertable_viral_unrarefied_euk_virus$Host[Mastertable_viral_unrarefied_euk_virus$Final_family == "Anelloviridae"] <- "small circular viruses"
Mastertable_viral_unrarefied_euk_virus$Host[Mastertable_viral_unrarefied_euk_virus$Final_family == "Genomoviridae"] <- "small circular viruses"
Mastertable_viral_unrarefied_euk_virus$Host[Mastertable_viral_unrarefied_euk_virus$Final_family == "Circoviridae"] <- "small circular viruses"
Mastertable_viral_unrarefied_euk_virus$Host[Mastertable_viral_unrarefied_euk_virus$Final_family == "Geminiviridae"] <- "small circular viruses"
Mastertable_viral_unrarefied_euk_virus$Host[Mastertable_viral_unrarefied_euk_virus$Final_family == "Smacoviridae"] <- "small circular viruses"

unique(Mastertable_viral_unrarefied_euk_virus$Host)
#view(Mastertable_viral_unrarefied_euk_virus)
####################################
# 3.3.2 Rarefied mastertable 
####################################
#View(Mastertable_viral_rarefied_euk_virus)

# (A) Create a new table to add 'group host annotation'
Mastertable_viral_rarefied_euk_virus_host <- Mastertable_viral_rarefied_euk_virus
unique(Mastertable_viral_rarefied_euk_virus_host$Final_family)

# (B) Add 'group host annotation' in a new column for plant and fungal viruses
Mastertable_viral_rarefied_euk_virus_host$Host <- Mastertable_viral_rarefied_euk_virus_host$Final_family
Mastertable_viral_rarefied_euk_virus_host$Host[Mastertable_viral_rarefied_euk_virus_host$Final_family == "Alphaflexiviridae"] <- "plant & fungal viruses"
Mastertable_viral_rarefied_euk_virus_host$Host[Mastertable_viral_rarefied_euk_virus_host$Final_family == "Closteroviridae"] <- "plant & fungal viruses"
Mastertable_viral_rarefied_euk_virus_host$Host[Mastertable_viral_rarefied_euk_virus_host$Final_family == "Betaflexiviridae"] <- "plant & fungal viruses"
Mastertable_viral_rarefied_euk_virus_host$Host[Mastertable_viral_rarefied_euk_virus_host$Final_family == "Chrysoviridae"] <- "plant & fungal viruses"
Mastertable_viral_rarefied_euk_virus_host$Host[Mastertable_viral_rarefied_euk_virus_host$Final_family == "Virgaviridae"] <- "plant & fungal viruses"
Mastertable_viral_rarefied_euk_virus_host$Host[Mastertable_viral_rarefied_euk_virus_host$Final_family == "Endornaviridae"] <- "plant & fungal viruses"
Mastertable_viral_rarefied_euk_virus_host$Host[Mastertable_viral_rarefied_euk_virus_host$Final_family == "Secoviridae"] <- "plant & fungal viruses"

# (C) Add 'group host annotation' in a new column for animal viruses
Mastertable_viral_rarefied_euk_virus_host$Host[Mastertable_viral_rarefied_euk_virus_host$Final_family == "Caliciviridae"] <- "animal viruses"
Mastertable_viral_rarefied_euk_virus_host$Host[Mastertable_viral_rarefied_euk_virus_host$Final_family == "Parvoviridae"] <- "animal viruses"
Mastertable_viral_rarefied_euk_virus_host$Host[Mastertable_viral_rarefied_euk_virus_host$Final_family == "Picornaviridae"] <- "animal viruses"
Mastertable_viral_rarefied_euk_virus_host$Host[Mastertable_viral_rarefied_euk_virus_host$Final_family == "Papillomaviridae"] <- "animal viruses"

# (D) Add 'group host annotation' in a new column for small circular and unclassified viruses
Mastertable_viral_rarefied_euk_virus_host$Host[Mastertable_viral_rarefied_euk_virus_host$Final_family == "unclassified CRESS"] <- "small circular viruses"
Mastertable_viral_rarefied_euk_virus_host$Host[Mastertable_viral_rarefied_euk_virus_host$Final_family == "Anelloviridae"] <- "small circular viruses"
Mastertable_viral_rarefied_euk_virus_host$Host[Mastertable_viral_rarefied_euk_virus_host$Final_family == "Genomoviridae"] <- "small circular viruses"
Mastertable_viral_rarefied_euk_virus_host$Host[Mastertable_viral_rarefied_euk_virus_host$Final_family == "Circoviridae"] <- "small circular viruses"
Mastertable_viral_rarefied_euk_virus_host$Host[Mastertable_viral_rarefied_euk_virus_host$Final_family == "Geminiviridae"] <- "small circular viruses"
Mastertable_viral_rarefied_euk_virus_host$Host[Mastertable_viral_rarefied_euk_virus_host$Final_family == "Smacoviridae"] <- "small circular viruses"

unique(Mastertable_viral_rarefied_euk_virus_host$Host)
#view(Mastertable_viral_rarefied_euk_virus_host)
####################################
# 3.3 Curated mastertables
#################################### 
Mastertable_viral_unrarefied_euk_virus_host # PRESENCE/ABSENCE
Mastertable_viral_rarefied_euk_virus_host # RA,ALPHA-DIVERSITY
####################################
# 3.4 Download the the curated eukaryotic viral list and save in excel format 
####################################
# Make sure here you download Family, genera level
# From the Unrarefied table & from rarefied table

setwd("/Users/daan/Desktop/Bioinformatics/Analysis/Biologicals/eukaryotic_viruses")
getwd()
dir()

# View(Mastertable_viral_rarefied_euk_virus_host)
vector_1 <- which(colnames(Mastertable_viral_rarefied_euk_virus_host) == "Final_class")
vector_2 <- which(colnames(Mastertable_viral_rarefied_euk_virus_host) == "Final_species")
vector_3 <- which(colnames(Mastertable_viral_rarefied_euk_virus_host) == "Host")
vector_4 <- which(colnames(Mastertable_viral_rarefied_euk_virus_host) == "length")
vector_5 <- which(colnames(Mastertable_viral_rarefied_euk_virus_host) == "Totalnumberofreads")
vector_6 <- which(colnames(Mastertable_viral_rarefied_euk_virus_host) == "Final_ANI")
vector_7 <- which(colnames(Mastertable_viral_rarefied_euk_virus_host) == "Final_coverage")
vector_8 <- which(colnames(Mastertable_viral_rarefied_euk_virus_host) == "Blastn_AS")
vector_9 <- which(colnames(Mastertable_viral_rarefied_euk_virus_host) == "Diamond_AS")
vector_10 <- which(colnames(Mastertable_viral_rarefied_euk_virus_host) == "Best_AS")

eukaryotic_viruses <- Mastertable_viral_rarefied_euk_virus_host[, c(vector_1:vector_2, vector_3, vector_4, vector_5, vector_6, vector_7, vector_8,vector_9, vector_10)]
eukaryotic_viruses <- with(eukaryotic_viruses,  eukaryotic_viruses[order(Final_family) , ])
#View(eukaryotic_viruses)

library("writexl")
write_xlsx(eukaryotic_viruses,"./eukaryotic_viruses.xlsx")
####################################
# 4. Relative abundances
####################################
## Use rarefied mastertable with hosts.
Mastertable_viral_rarefied_euk_virus_host
####################################
# 4.1 stacked barplot and others
####################################
# use esquisser
# create area plot or stacked barplot (have code for that in script)
## Repeat this for family level, genus level, hosts, ...
## create other smple graphs here
## could be that you have to change characters to integers or others for some plots.
## Optimze it using ggThemeAssist.

## Normalized barplot or stacked barplot.

esquisser(Mastertable_viral_rarefied_euk_virus_host)

####################################
# 4.2 stacked barplot: absolute number of reads
####################################
nrow(Mastertable_viral_rarefied_euk_virus_host)
vector_1_1 <- (which(names(Mastertable_viral_rarefied_euk_virus_host)== "Virsorter")-1) 
Mastertable_viral_rarefied_categories_vector <- colnames(Mastertable_viral_rarefied_euk_virus_host[, c(1:vector_1_1)])

# (A) Aggregate
# (A.1) families
Mastertable_viral_rarefied_categories_euk_fam <- aggregate(. ~Final_family, FUN = sum, data = Mastertable_viral_rarefied_euk_virus_host[,colnames(Mastertable_viral_rarefied_euk_virus_host) %in% Mastertable_viral_rarefied_categories_vector | colnames(Mastertable_viral_rarefied_euk_virus_host) == 'Final_family'])
rownames(Mastertable_viral_rarefied_categories_euk_fam) <- Mastertable_viral_rarefied_categories_euk_fam$Final_family
Mastertable_viral_rarefied_categories_euk_fam$Final_family <- NULL

# (A.2) genera
Mastertable_viral_rarefied_categories_euk_gen <- aggregate(. ~Final_genus, FUN = sum, data = Mastertable_viral_rarefied_euk_virus_host[,colnames(Mastertable_viral_rarefied_euk_virus_host) %in% Mastertable_viral_rarefied_categories_vector | colnames(Mastertable_viral_rarefied_euk_virus_host) == 'Final_genus'])
rownames(Mastertable_viral_rarefied_categories_euk_gen) <- Mastertable_viral_rarefied_categories_euk_gen$Final_genus
Mastertable_viral_rarefied_categories_euk_gen$Final_genus <- NULL
View(Mastertable_viral_rarefied_categories_euk_gen)

# (B) melt 
# (B.1) families
Mastertable_viral_rarefied_categories_euk_fam <-  t(Mastertable_viral_rarefied_categories_euk_fam)
Mastertable_viral_rarefied_categories_euk_fam <- melt(Mastertable_viral_rarefied_categories_euk_fam)
colnames(Mastertable_viral_rarefied_categories_euk_fam) <- c("sample", "groups", "reads")
#View(Mastertable_viral_rarefied_categories_euk_fam)

# (B.2) genera
Mastertable_viral_rarefied_categories_euk_gen <-  t(Mastertable_viral_rarefied_categories_euk_gen)
Mastertable_viral_rarefied_categories_euk_gen <- melt(Mastertable_viral_rarefied_categories_euk_gen)
colnames(Mastertable_viral_rarefied_categories_euk_gen) <- c("sample", "groups", "reads")
#view(Mastertable_viral_rarefied_categories_euk_gen)

# (C) stacked barplot: absolute number of reads
# (C.1) families
ggplot(Mastertable_viral_rarefied_categories_euk_fam, aes(x=sample, y= reads)) +
  geom_area(aes(y= reads, group = groups, fill = groups), position = 'stack', show.legend = T) +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_fill_manual(values = rev(c(inferno(length(unique(Mastertable_viral_rarefied_categories_euk_fam$groups))/2,1, begin = 0.2, end = 0.9), rev(plasma(length(unique(Mastertable_viral_rarefied_categories_euk_fam$groups))/2,1, begin = 0.3, end = 1))))) +
  theme(legend.title = element_blank(), legend.text = element_text(size = 9)) + 
  xlab("") +
  ylab("reads (#)") +
  theme(axis.title.y = element_text(size=9)) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  theme(axis.text.x = element_text(size = 5, angle = 45, vjust = 1 , hjust = 1, color = "black"), axis.title = element_text(size = 9)) +
  theme(axis.text.y = element_text(color = "black"))

# (C.2) genera
ggplot(Mastertable_viral_rarefied_categories_euk_gen, aes(x=sample, y= reads)) +
  geom_area(aes(y= reads, group = groups, fill = groups), position = 'stack', show.legend = T) +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_fill_manual(values = rev(c(inferno(length(unique(Mastertable_viral_rarefied_categories_euk_gen$groups))/2,1, begin = 0.2, end = 0.9), rev(plasma(length(unique(Mastertable_viral_rarefied_categories_euk_gen$groups))/2,1, begin = 0.3, end = 1))))) +
  theme(legend.title = element_blank(), legend.text = element_text(size = 9)) + 
  xlab("") +
  ylab("reads (#)") +
  theme(axis.title.y = element_text(size=9)) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  theme(axis.text.x = element_text(size = 5, angle = 45, vjust = 1 , hjust = 1, color = "black"), axis.title = element_text(size = 9)) +
  theme(axis.text.y = element_text(color = "black"))
####################################
# 4.3 Heatmap relative abundances
####################################
# 4.3.1 All samples
####################################
# (A) aggregate 
# (A.1) families
Mastertable_viral_rarefied_categories_euk_fam <- aggregate(. ~Final_family, FUN = sum, data = Mastertable_viral_rarefied_euk_virus_host[,colnames(Mastertable_viral_rarefied_euk_virus_host) %in% Mastertable_viral_rarefied_categories_vector | colnames(Mastertable_viral_rarefied_euk_virus_host) == 'Final_family'])
rownames(Mastertable_viral_rarefied_categories_euk_fam) <- Mastertable_viral_rarefied_categories_euk_fam$Final_family
Mastertable_viral_rarefied_categories_euk_fam$Final_family <- NULL
view(Mastertable_viral_rarefied_categories_euk_fam)

# (A.2) genera
Mastertable_viral_rarefied_categories_euk_gen <- aggregate(. ~Final_genus, FUN = sum, data = Mastertable_viral_rarefied_euk_virus_host[,colnames(Mastertable_viral_rarefied_euk_virus_host) %in% Mastertable_viral_rarefied_categories_vector | colnames(Mastertable_viral_rarefied_euk_virus_host) == 'Final_genus'])
rownames(Mastertable_viral_rarefied_categories_euk_gen) <- Mastertable_viral_rarefied_categories_euk_gen$Final_genus
Mastertable_viral_rarefied_categories_euk_gen$Final_genus <- NULL
view(Mastertable_viral_rarefied_categories_euk_gen)

# (B) normalization: relative abundances
# (B.1) families
Mastertable_viral_rarefied_categories_euk_fam_normalized <- sweep(Mastertable_viral_rarefied_categories_euk_fam, 2, colSums(Mastertable_viral_rarefied_categories_euk_fam), '/') 
Mastertable_viral_rarefied_categories_euk_fam_normalized[is.na(Mastertable_viral_rarefied_categories_euk_fam_normalized)] <- 0
Mastertable_viral_rarefied_categories_euk_fam_normalized <- t(Mastertable_viral_rarefied_categories_euk_fam_normalized)

# (B.2) genera
Mastertable_viral_rarefied_categories_euk_gen_normalized <- sweep(Mastertable_viral_rarefied_categories_euk_gen, 2, colSums(Mastertable_viral_rarefied_categories_euk_gen), '/') 
Mastertable_viral_rarefied_categories_euk_gen_normalized[is.na(Mastertable_viral_rarefied_categories_euk_gen_normalized)] <- 0
Mastertable_viral_rarefied_categories_euk_gen_normalized <- t(Mastertable_viral_rarefied_categories_euk_gen_normalized)

# (C) Add samples ID in dataframe
# (C.1) families
Mastertable_viral_rarefied_categories_euk_fam_normalized <- rownames_to_column(as.data.frame(Mastertable_viral_rarefied_categories_euk_fam_normalized))
names(Mastertable_viral_rarefied_categories_euk_fam_normalized)[names(Mastertable_viral_rarefied_categories_euk_fam_normalized) == "rowname"] <- "sample"
rownames(Mastertable_viral_rarefied_categories_euk_fam_normalized) <- Mastertable_viral_rarefied_categories_euk_fam_normalized$sample

# (C.2) genera
Mastertable_viral_rarefied_categories_euk_gen_normalized <- rownames_to_column(as.data.frame(Mastertable_viral_rarefied_categories_euk_gen_normalized))
names(Mastertable_viral_rarefied_categories_euk_gen_normalized)[names(Mastertable_viral_rarefied_categories_euk_gen_normalized) == "rowname"] <- "sample"
rownames(Mastertable_viral_rarefied_categories_euk_gen_normalized) <- Mastertable_viral_rarefied_categories_euk_gen_normalized$sample

# (D) melt
# (D.1) families
Mastertable_viral_rarefied_categories_euk_fam_normalized <- melt(Mastertable_viral_rarefied_categories_euk_fam_normalized, id.vars = c("sample"))
colnames(Mastertable_viral_rarefied_categories_euk_fam_normalized) <- c("sample", "groups", "reads")

# (D.2) genera
Mastertable_viral_rarefied_categories_euk_gen_normalized <- melt(Mastertable_viral_rarefied_categories_euk_gen_normalized, id.vars = c("sample"))
colnames(Mastertable_viral_rarefied_categories_euk_gen_normalized) <- c("sample", "groups", "reads")

# (E) stacked barplot: relative abundances
# (E.1) families (two different color panels)
ggplot(Mastertable_viral_rarefied_categories_euk_fam_normalized, aes(x=sample, y= reads)) +
  geom_area(aes(y= reads, group = groups, fill = groups), position = 'stack', show.legend = T) +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_fill_manual(values = rev(c(inferno(length(unique(Mastertable_viral_rarefied_categories_euk_fam_normalized$groups))/2,1, begin = 0.2, end = 0.9), rev(plasma(length(unique(Mastertable_viral_rarefied_categories_euk_fam_normalized$groups))/2,1, begin = 0.3, end = 1))))) +
  theme(legend.title = element_blank(), legend.text = element_text(size = 9)) + 
  xlab("") +
  ylab("relative abundance") +
  theme(axis.title.y = element_text(size=9)) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  theme(axis.text.x = element_text(size = 5, angle = 45, vjust = 1 , hjust = 1, color = "black"), axis.title = element_text(size = 9)) +
  theme(axis.text.y = element_text(color = "black"))

ggplot(Mastertable_viral_rarefied_categories_euk_fam_normalized, aes(x=sample, y= reads)) +
  geom_area(aes(y= reads, group = groups, fill = groups), position = 'stack', show.legend = T) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_fill_manual(values = rev(c(viridis(length(unique(Mastertable_viral_rarefied_categories_euk_fam_normalized$groups))/2,1, begin = 0.2, end = 0.9), rev(inferno(length(unique(Mastertable_viral_rarefied_categories_euk_fam_normalized$groups))/2,1, begin = 0.3, end = 1))))) +
  theme(legend.title = element_blank(), legend.text = element_text(size = 9)) + 
  xlab("") +
  ylab("relative abundance") +
  theme(axis.title.y = element_text(size=9)) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  theme(axis.text = element_text(colour = "black")) +
  theme(axis.text.x = element_text(size = 5, angle = 45, vjust = 1 , hjust = 1, color = "black"), axis.title = element_text(size = 9))

# (E.2) genera
ggplot(Mastertable_viral_rarefied_categories_euk_gen_normalized, aes(x=sample, y= reads)) +
  geom_area(aes(y= reads, group = groups, fill = groups), position = 'stack', show.legend = T) +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_fill_manual(values = rev(c(plasma(length(unique(Mastertable_viral_rarefied_categories_euk_gen_normalized$groups))/2,1, begin = 0.2, end = 0.9), rev(viridis(length(unique(Mastertable_viral_rarefied_categories_euk_gen_normalized$groups))/2,1, begin = 0.3, end = 1))))) +
  theme(legend.title = element_blank(), legend.text = element_text(size = 9)) + 
  xlab("") +
  ylab("relative abundance") +
  theme(axis.title.y = element_text(size=9)) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  theme(axis.text.x = element_text(size = 5, angle = 45, vjust = 1 , hjust = 1, color = "black"), axis.title = element_text(size = 9)) +
  theme(axis.text.y = element_text(color = "black"))

ggplot(Mastertable_viral_rarefied_categories_euk_gen_normalized, aes(x=sample, y= reads)) +
  geom_area(aes(y= reads, group = groups, fill = groups), position = 'stack', show.legend = T) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_fill_manual(values = rev(c(viridis(length(unique(Mastertable_viral_rarefied_categories_euk_gen_normalized$groups))/2,1, begin = 0.2, end = 0.9), rev(inferno(length(unique(Mastertable_viral_rarefied_categories_euk_gen_normalized$groups))/2,1, begin = 0.3, end = 1))))) +
  theme(legend.title = element_blank(), legend.text = element_text(size = 9)) + 
  xlab("") +
  ylab("relative abundance") +
  theme(axis.title.y = element_text(size=9)) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  theme(axis.text = element_text(colour = "black")) +
  theme(axis.text.x = element_text(size = 5, angle = 45, vjust = 1 , hjust = 1, color = "black"), axis.title = element_text(size = 9))
####################################
# 4.3.2 Samples with eukaryotic virus presence
####################################
# (A) subset samples with a presence of euk. viruses
# (A.1) families
Mastertable_viral_rarefied_categories_euk_fam_normalized_vector <- Mastertable_viral_rarefied_categories_euk_fam_normalized %>% filter(reads > 0)
vector_euk <- unique(Mastertable_viral_rarefied_categories_euk_fam_normalized_vector$sample)
#view(vector_euk) # 38 samples with at least one euk. virus
Mastertable_viral_rarefied_categories_euk_fam_normalized <- Mastertable_viral_rarefied_categories_euk_fam_normalized[Mastertable_viral_rarefied_categories_euk_fam_normalized$sample %in% vector_euk,]

# (A.1) genera
Mastertable_viral_rarefied_categories_euk_gen_normalized <- Mastertable_viral_rarefied_categories_euk_gen_normalized[Mastertable_viral_rarefied_categories_euk_gen_normalized$sample %in% vector_euk,]

# (B) stacked barplot: relative abundances
# (B.1) families
ggplot(Mastertable_viral_rarefied_categories_euk_fam_normalized, aes(x=sample, y= reads)) +
  geom_area(aes(y= reads, group = groups, fill = groups), position = 'stack', show.legend = T) +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_fill_manual(values = rev(c(viridis(length(unique(Mastertable_viral_rarefied_categories_euk_fam_normalized$groups))/2,1, begin = 0.2, end = 0.9), rev(plasma(length(unique(Mastertable_viral_rarefied_categories_euk_fam_normalized$groups))/2,1, begin = 0.3, end = 1))))) +
  theme(legend.title = element_blank(), legend.text = element_text(size = 9)) + 
  xlab("") +
  ylab("relative abundance") +
  theme(axis.title.y = element_text(size=9)) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  theme(axis.text.x = element_text(size = 5, angle = 45, vjust = 1 , hjust = 1, color = "black"), axis.title = element_text(size = 9)) +
  theme(axis.text.y = element_text(color = "black"))

# (B.2) genera
ggplot(Mastertable_viral_rarefied_categories_euk_gen_normalized, aes(x=sample, y= reads)) +
  geom_area(aes(y= reads, group = groups, fill = groups), position = 'stack', show.legend = T) +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_fill_manual(values = rev(c(viridis(length(unique(Mastertable_viral_rarefied_categories_euk_gen_normalized$groups))/2,1, begin = 0.2, end = 0.9), rev(plasma(length(unique(Mastertable_viral_rarefied_categories_euk_gen_normalized$groups))/2,1, begin = 0.3, end = 1))))) +
  theme(legend.title = element_blank(), legend.text = element_text(size = 9)) + 
  xlab("") +
  ylab("relative abundance") +
  theme(axis.title.y = element_text(size=9)) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  theme(axis.text.x = element_text(size = 5, angle = 45, vjust = 1 , hjust = 1, color = "black"), axis.title = element_text(size = 9)) +
  theme(axis.text.y = element_text(color = "black"))

# (H) Go to Illustrator - finish & save figures for RA heatmap
####################################
# 4.4 Complex heatmap
####################################
# (1) ALL SAMPLES

Mastertable_viral_rarefied_euk_virus_host

# (A) subset eukaryotic viral families viruses and aggregate
# (A.1) families
Mastertable_global_rarefied_euk_fam_host <- Mastertable_viral_rarefied_euk_virus_host
Mastertable_global_rarefied_euk_fam_host <- aggregate(. ~Final_family, FUN = sum, data = Mastertable_viral_rarefied_euk_virus_host[,colnames(Mastertable_viral_rarefied_euk_virus_host) %in% Mastertable_viral_rarefied_categories_vector | colnames(Mastertable_viral_rarefied_euk_virus_host) == 'Final_family'])
rownames(Mastertable_global_rarefied_euk_fam_host) <- Mastertable_global_rarefied_euk_fam_host$Final_family
Mastertable_global_rarefied_euk_fam_host$Final_family <- NULL

# (A.2) genera
Mastertable_global_rarefied_euk_gen_host <- Mastertable_viral_rarefied_euk_virus_host
Mastertable_global_rarefied_euk_gen_host <- aggregate(. ~Final_genus, FUN = sum, data = Mastertable_global_rarefied_euk_gen_host[,colnames(Mastertable_global_rarefied_euk_gen_host) %in% Mastertable_viral_rarefied_categories_vector | colnames(Mastertable_global_rarefied_euk_gen_host) == 'Final_genus'])
rownames(Mastertable_global_rarefied_euk_gen_host) <- Mastertable_global_rarefied_euk_gen_host$Final_genus
Mastertable_global_rarefied_euk_gen_host$Final_genus <- NULL

# (B) convert reads to relative abundances + add sample ID
# (B.1) families
Mastertable_global_rarefied_euk_fam_host_normalized <- sweep(Mastertable_global_rarefied_euk_fam_host, 2, colSums(Mastertable_global_rarefied_euk_fam_host), '/') 
Mastertable_global_rarefied_euk_fam_host_normalized[is.na(Mastertable_global_rarefied_euk_fam_host_normalized)] <- 0
Mastertable_global_rarefied_euk_fam_host_normalized <- t(Mastertable_global_rarefied_euk_fam_host_normalized)

# (B.2) genera
Mastertable_global_rarefied_euk_gen_host_normalized <- sweep(Mastertable_global_rarefied_euk_gen_host, 2, colSums(Mastertable_global_rarefied_euk_gen_host), '/') 
Mastertable_global_rarefied_euk_gen_host_normalized[is.na(Mastertable_global_rarefied_euk_gen_host_normalized)] <- 0
Mastertable_global_rarefied_euk_gen_host_normalized <- t(Mastertable_global_rarefied_euk_gen_host_normalized)

# (C) melt & subset samples with presence of euk. viruses
# (C.1) families
Mastertable_global_rarefied_euk_fam_host_m <- melt(Mastertable_global_rarefied_euk_fam_host_normalized)
Mastertable_global_rarefied_euk_fam_host_m_subset <- Mastertable_global_rarefied_euk_fam_host_m[Mastertable_global_rarefied_euk_fam_host_m$Var1 %in% vector_euk,]

# (C.2) genera
Mastertable_global_rarefied_euk_gen_host_m <- melt(Mastertable_global_rarefied_euk_gen_host_normalized)
Mastertable_global_rarefied_euk_gen_host_m_subset <- Mastertable_global_rarefied_euk_gen_host_m[Mastertable_global_rarefied_euk_gen_host_m$Var1 %in% vector_euk,]

# (D) heatmap: relativer abundance
# (D.1) families
ggplot(Mastertable_global_rarefied_euk_fam_host_m_subset, aes(x= Var1, y= Var2)) + 
  geom_tile(aes(fill=value) , colour = "grey", show.legend = T) +
  scale_fill_gradient(low = 'white', high = '#3182bd', limits = c(0,1), breaks=c(0,0.25,0.5,0.75,1),labels=c("0%", "25%", "50%","75%" ,"100%")) +
  ggtitle("") +
  xlab("") +
  ylab("") +
  theme_minimal() +
  coord_equal() +
  theme(axis.text.y = element_text(color = "black", size = 7)) +
  scale_y_discrete(limits = unique(rev(Mastertable_global_rarefied_euk_fam_host_m_subset$Var2)), expand = c(0,0)) +  
  scale_x_discrete(expand = c(0,0)) +
  theme(panel.background = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.title = element_blank(), legend.text = element_text(size = 10), legend.key.size = unit(1,"line")) +
  theme(plot.title = element_text(hjust = -0.2 , vjust=2.12, size = 7, face ="bold")) +
  theme(axis.text.x = element_text(size = 5, angle = 45, vjust = 1 , hjust = 1, color = "black"), axis.title = element_text(size = 9)) 

# (D.2) genera
ggplot(Mastertable_global_rarefied_euk_gen_host_m_subset, aes(x= Var1, y= Var2)) + 
  geom_tile(aes(fill=value) , colour = "grey", show.legend = T) +
  scale_fill_gradient(low = 'white', high = '#3182bd', limits = c(0,1), breaks=c(0,0.25,0.5,0.75,1),labels=c("0%", "25%", "50%","75%" ,"100%")) +
  ggtitle("") +
  xlab("") +
  ylab("") +
  theme_minimal() +
  coord_equal() +
  theme(axis.text.y = element_text(color = "black", size = 7)) +
  scale_y_discrete(limits = unique(rev(Mastertable_global_rarefied_euk_gen_host_m_subset$Var2)), expand = c(0,0)) +  
  scale_x_discrete(expand = c(0,0)) +
  theme(panel.background = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.title = element_blank(), legend.text = element_text(size = 10), legend.key.size = unit(1,"line")) +
  theme(plot.title = element_text(hjust = -0.2 , vjust=2.12, size = 7, face ="bold")) +
  theme(axis.text.x = element_text(size = 5, angle = 45, vjust = 1 , hjust = 1, color = "black"), axis.title = element_text(size = 9))

# (2) only samples with eukaryotic viruses present
# (A) Create a vector for animal viruses, plant and fungal viruses, small circular viruses and unclassified viruses
# (A.1) families
unique(Mastertable_viral_rarefied_euk_virus_host$Host)
## FIX IF NEEDED
animal_viruses_fam_vector <-  Mastertable_viral_rarefied_euk_virus_host[Mastertable_viral_rarefied_euk_virus_host$Host == "animal viruses", c(1:vector_1_1)]
plant_and_fungal_viruses_fam_vector <- Mastertable_viral_rarefied_euk_virus_host$sample
small_circular_viruses_fam_vector <- Mastertable_viral_rarefied_euk_virus_host$sample
unclassified_viruses_fam_vector <- Mastertable_viral_rarefied_euk_virus_host$sample
View(animal_viruses_fam_vector)

# (A.2) genera
animal_viruses_gen_vector <- Mastertable_global_rarefied_euk_gen_host_animal_preabs$sample
plant_and_fungal_viruses_gen_vector <- Mastertable_global_rarefied_euk_gen_host_plant_preabs$sample
small_circular_viruses_gen_vector <- Mastertable_global_rarefied_euk_gen_host_scv_preabs$sample
unclassified_viruses_gen_vector <- Mastertable_global_rarefied_euk_gen_host_unc_preabs$sample

# (B) animal viruses
# (B.1) families
animal_viruses_fam_RA <- Mastertable_global_rarefied_euk_fam_host_m_subset[Mastertable_global_rarefied_euk_fam_host_m_subset$Var2 %in% animal_viruses_fam_vector,]

animal_viruses_fam_RA_plot <- ggplot(animal_viruses_fam_RA, aes(x= Var1, y= Var2)) + 
  geom_tile(aes(fill=value) , colour = "grey", show.legend = F) +
  scale_fill_gradient(low = 'white', high = '#3182bd', limits = c(0,1), breaks=c(0,0.25,0.5,0.75,1),labels=c("0%", "25%", "50%","75%" ,"100%")) +
  ggtitle("") +
  xlab("") +
  ylab("") +
  theme_minimal() +
  coord_equal() +
  theme(axis.text.y = element_text(color = "black", size = 7)) +
  scale_y_discrete(limits = unique(rev(animal_viruses_fam_RA$Var2)), expand = c(0,0)) +  
  scale_x_discrete(expand = c(0,0)) +
  theme(panel.background = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.title = element_blank(), legend.text = element_text(size = 10), legend.key.size = unit(1,"line")) +
  theme(plot.title = element_text(hjust = -0.2 , vjust=2.12, size = 7, face ="bold")) +
  theme(axis.text.x = element_text(size = 5, angle = 45, vjust = 1 , hjust = 1, color = "black"), axis.title = element_text(size = 9)) 

# (B.2) genera
animal_viruses_gen_RA <- Mastertable_global_rarefied_euk_gen_host_m_subset[Mastertable_global_rarefied_euk_gen_host_m_subset$Var2 %in% animal_viruses_gen_vector,]

animal_viruses_gen_RA_plot <- ggplot(animal_viruses_gen_RA, aes(x= Var1, y= Var2)) + 
  geom_tile(aes(fill=value) , colour = "grey", show.legend = F) +
  scale_fill_gradient(low = 'white', high = '#3182bd', limits = c(0,1), breaks=c(0,0.25,0.5,0.75,1),labels=c("0%", "25%", "50%","75%" ,"100%")) +
  ggtitle("") +
  xlab("") +
  ylab("") +
  theme_minimal() +
  coord_equal() +
  theme(axis.text.y = element_text(color = "black", size = 7)) +
  scale_y_discrete(limits = unique(rev(animal_viruses_gen_RA$Var2)), expand = c(0,0)) +  
  scale_x_discrete(expand = c(0,0)) +
  theme(panel.background = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.title = element_blank(), legend.text = element_text(size = 10), legend.key.size = unit(1,"line")) +
  theme(plot.title = element_text(hjust = -0.2 , vjust=2.12, size = 7, face ="bold")) +
  theme(axis.text.x = element_text(size = 5, angle = 45, vjust = 1 , hjust = 1, color = "black"), axis.title = element_text(size = 9)) 

# (C) plant and fungal viruses
# (C.1) families
plant_and_fungal_viruses_fam_RA <- Mastertable_global_rarefied_euk_fam_host_m_subset[Mastertable_global_rarefied_euk_fam_host_m_subset$Var2 %in% plant_and_fungal_viruses_fam_vector,]

plant_and_fungal_viruses_fam_RA_plot <- ggplot(plant_and_fungal_viruses_fam_RA, aes(x= Var1, y= Var2)) + 
  geom_tile(aes(fill=value) , colour = "grey", show.legend = F) +
  scale_fill_gradient(low = 'white', high = '#3182bd', limits = c(0,1), breaks=c(0,0.25,0.5,0.75,1),labels=c("0%", "25%", "50%","75%" ,"100%")) +
  ggtitle("") +
  xlab("") +
  ylab("") +
  theme_minimal() +
  coord_equal() +
  theme(axis.text.y = element_text(color = "black", size = 7)) +
  scale_y_discrete(limits = unique(rev(plant_and_fungal_viruses_fam_RA$Var2)), expand = c(0,0)) +  
  scale_x_discrete(expand = c(0,0)) +
  theme(panel.background = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.title = element_blank(), legend.text = element_text(size = 10), legend.key.size = unit(1,"line")) +
  theme(plot.title = element_text(hjust = -0.2 , vjust=2.12, size = 7, face ="bold")) +
  theme(axis.text.x = element_text(size = 5, angle = 45, vjust = 1 , hjust = 1, color = "black"), axis.title = element_text(size = 9)) 

# (C.2) genera
plant_and_fungal_viruses_gen_RA <- Mastertable_global_rarefied_euk_gen_host_m_subset[Mastertable_global_rarefied_euk_gen_host_m_subset$Var2 %in% plant_and_fungal_viruses_gen_vector,]

plant_and_fungal_viruses_gen_RA_plot <- ggplot(plant_and_fungal_viruses_gen_RA, aes(x= Var1, y= Var2)) + 
  geom_tile(aes(fill=value) , colour = "grey", show.legend = T) +
  scale_fill_gradient(low = 'white', high = '#3182bd', limits = c(0,1), breaks=c(0,0.25,0.5,0.75,1),labels=c("0%", "25%", "50%","75%" ,"100%")) +
  ggtitle("") +
  xlab("") +
  ylab("") +
  theme_minimal() +
  coord_equal() +
  theme(axis.text.y = element_text(color = "black", size = 7)) +
  scale_y_discrete(limits = unique(rev(plant_and_fungal_viruses_gen_RA$Var2)), expand = c(0,0)) +  
  scale_x_discrete(expand = c(0,0)) +
  theme(panel.background = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.title = element_blank(), legend.text = element_text(size = 10), legend.key.size = unit(1,"line")) +
  theme(plot.title = element_text(hjust = -0.2 , vjust=2.12, size = 7, face ="bold")) +
  theme(axis.text.x = element_text(size = 5, angle = 45, vjust = 1 , hjust = 1, color = "black"), axis.title = element_text(size = 9)) 

plant_and_fungal_viruses_gen_RA_plot

# (D) small circular viruses
# (D.1) families
small_circular_viruses_fam_RA <- Mastertable_global_rarefied_euk_fam_host_m_subset[Mastertable_global_rarefied_euk_fam_host_m_subset$Var2 %in% small_circular_viruses_fam_vector,]

small_circular_viruses_fam_RA_plot <- ggplot(small_circular_viruses_fam_RA, aes(x= Var1, y= Var2)) + 
  geom_tile(aes(fill=value) , colour = "grey", show.legend = F) +
  scale_fill_gradient(low = 'white', high = '#3182bd', limits = c(0,1), breaks=c(0,0.25,0.5,0.75,1),labels=c("0%", "25%", "50%","75%" ,"100%")) +
  ggtitle("") +
  xlab("") +
  ylab("") +
  theme_minimal() +
  coord_equal() +
  theme(axis.text.y = element_text(color = "black", size = 7)) +
  scale_y_discrete(limits = unique(rev(small_circular_viruses_fam_RA$Var2)), expand = c(0,0)) +  
  scale_x_discrete(expand = c(0,0)) +
  theme(panel.background = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.title = element_blank(), legend.text = element_text(size = 10), legend.key.size = unit(1,"line")) +
  theme(plot.title = element_text(hjust = -0.2 , vjust=2.12, size = 7, face ="bold")) +
  theme(axis.text.x = element_text(size = 5, angle = 45, vjust = 1 , hjust = 1, color = "black"), axis.title = element_text(size = 9)) 

# (D.2) genera
small_circular_viruses_gen_RA <- Mastertable_global_rarefied_euk_gen_host_m_subset[Mastertable_global_rarefied_euk_gen_host_m_subset$Var2 %in% small_circular_viruses_gen_vector,]

small_circular_viruses_gen_RA_plot <- ggplot(small_circular_viruses_gen_RA, aes(x= Var1, y= Var2)) + 
  geom_tile(aes(fill=value) , colour = "grey", show.legend = T) +
  scale_fill_gradient(low = 'white', high = '#3182bd', limits = c(0,1), breaks=c(0,0.25,0.5,0.75,1),labels=c("0%", "25%", "50%","75%" ,"100%")) +
  ggtitle("") +
  xlab("") +
  ylab("") +
  theme_minimal() +
  coord_equal() +
  theme(axis.text.y = element_text(color = "black", size = 7)) +
  scale_y_discrete(limits = unique(rev(small_circular_viruses_gen_RA$Var2)), expand = c(0,0)) +  
  scale_x_discrete(expand = c(0,0)) +
  theme(panel.background = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.title = element_blank(), legend.text = element_text(size = 10), legend.key.size = unit(1,"line")) +
  theme(plot.title = element_text(hjust = -0.2 , vjust=2.12, size = 7, face ="bold")) +
  theme(axis.text.x = element_text(size = 5, angle = 45, vjust = 1 , hjust = 1, color = "black"), axis.title = element_text(size = 9)) 

# (E) unclassified viruses
# (E.1) families
unclassified_viruses_fam_RA <- Mastertable_global_rarefied_euk_fam_host_m_subset[Mastertable_global_rarefied_euk_fam_host_m_subset$Var2 %in% unclassified_viruses_fam_vector,]

unclassified_viruses_fam_RA_plot <- ggplot(unclassified_viruses_fam_RA, aes(x= Var1, y= Var2)) + 
  geom_tile(aes(fill=value) , colour = "grey", show.legend = F) +
  scale_fill_gradient(low = 'white', high = '#3182bd', limits = c(0,1), breaks=c(0,0.25,0.5,0.75,1),labels=c("0%", "25%", "50%","75%" ,"100%")) +
  ggtitle("") +
  xlab("") +
  ylab("") +
  theme_minimal() +
  coord_equal() +
  theme(axis.text.y = element_text(color = "black", size = 7)) +
  scale_y_discrete(limits = unique(rev(unclassified_viruses_fam_RA$Var2)), expand = c(0,0)) +  
  scale_x_discrete(expand = c(0,0)) +
  theme(panel.background = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.title = element_blank(), legend.text = element_text(size = 10), legend.key.size = unit(1,"line")) +
  theme(plot.title = element_text(hjust = -0.2 , vjust=2.12, size = 7, face ="bold")) +
  theme(axis.text.x = element_text(size = 5, angle = 45, vjust = 1 , hjust = 1, color = "black"), axis.title = element_text(size = 9)) 

# (E.2) genera
unclassified_viruses_gen_RA <- Mastertable_global_rarefied_euk_gen_host_m_subset[Mastertable_global_rarefied_euk_gen_host_m_subset$Var2 %in% unclassified_viruses_gen_vector,]

unclassified_viruses_gen_RA_plot <- ggplot(unclassified_viruses_gen_RA, aes(x= Var1, y= Var2)) + 
  geom_tile(aes(fill=value) , colour = "grey", show.legend = T) +
  scale_fill_gradient(low = 'white', high = '#3182bd', limits = c(0,1), breaks=c(0,0.25,0.5,0.75,1),labels=c("0%", "25%", "50%","75%" ,"100%")) +
  ggtitle("") +
  xlab("") +
  ylab("") +
  theme_minimal() +
  coord_equal() +
  theme(axis.text.y = element_text(color = "black", size = 7)) +
  scale_y_discrete(limits = unique(rev(unclassified_viruses_gen_RA$Var2)), expand = c(0,0)) +  
  scale_x_discrete(expand = c(0,0)) +
  theme(panel.background = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.title = element_blank(), legend.text = element_text(size = 10), legend.key.size = unit(1,"line")) +
  theme(plot.title = element_text(hjust = -0.2 , vjust=2.12, size = 7, face ="bold")) +
  theme(axis.text.x = element_text(size = 5, angle = 45, vjust = 1 , hjust = 1, color = "black"), axis.title = element_text(size = 9)) 

# (F) save these plots
# (F.1) families
animal_viruses_fam_RA_plot
plant_and_fungal_viruses_fam_RA_plot
small_circular_viruses_fam_RA_plot
unclassified_viruses_fam_RA_plot

# (F.2) genera
animal_viruses_gen_RA_plot
plant_and_fungal_viruses_gen_RA_plot
small_circular_viruses_gen_RA_plot
unclassified_viruses_gen_RA_plot

# (G) One figure containing all heatmaps: we want the same format, which is done down below. This can be saved an imported in Illustrator to modify and finish figure.
# (G.1) families
animal_viruses_fam_RA_no_leg <- animal_viruses_fam_RA_plot + theme(legend.position='none',  axis.title.y=element_blank(),  axis.text.y=element_blank(), axis.text.x = element_blank())
plant_and_fungal_viruses_fam_RA_no_leg <- plant_and_fungal_viruses_fam_RA_plot + theme(legend.position='none', axis.title.y=element_blank(),  axis.text.y=element_blank(),axis.text.x = element_blank())
small_circular_viruses_fam_RA_no_leg <- small_circular_viruses_fam_RA_plot + theme(legend.position='none', axis.title.y=element_blank(), axis.text.y=element_blank(),axis.text.x = element_blank())
unclassified_viruses_fam_RA_plot_no_leg <- unclassified_viruses_fam_RA_plot + theme(legend.position='none', axis.title.y=element_blank(), axis.text.y=element_blank())

# (G.2) genera
animal_viruses_gen_RA_no_leg <- animal_viruses_gen_RA_plot + theme(legend.position='none',  axis.title.y=element_blank(),  axis.text.y=element_blank(),axis.text.x = element_blank())
plant_and_fungal_viruses_gen_RA_no_leg <- plant_and_fungal_viruses_gen_RA_plot + theme(legend.position='none', axis.title.y=element_blank(),  axis.text.y=element_blank(),axis.text.x = element_blank())
small_circular_viruses_gen_RA_no_leg <- small_circular_viruses_gen_RA_plot + theme(legend.position='none', axis.title.y=element_blank(), axis.text.y=element_blank(),axis.text.x = element_blank())
unclassified_viruses_gen_RA_plot_noleg <- unclassified_viruses_gen_RA_plot + theme(legend.position='none', axis.title.y=element_blank(), axis.text.y=element_blank())

# (G.X) keep figure with a legend (to import legend in Illustrator)
plant_and_fungal_viruses_gen_RA_plot

animal_viruses_fam_RA_no_leg
plant_and_fungal_viruses_fam_RA_no_leg
small_circular_viruses_fam_RA_no_leg
unclassified_viruses_fam_RA_plot_no_leg

animal_viruses_gen_RA_no_leg
plant_and_fungal_viruses_gen_RA_no_leg
small_circular_viruses_gen_RA_no_leg
unclassified_viruses_gen_RA_plot_noleg
####################################
# 5. Presence/absence
####################################
# Presence/absence data is qualitative date and not quantitative data.
# Therefore, we can work with the unrarefied table.
####################################
# 5.1 Heatmap
####################################
## Use unrarefied mastertable with hosts.
nrow(Mastertable_viral_unrarefied_euk_virus_host)
####################################
# 5.1.1 All samples
####################################
# (1) ALL SAMPLES

Mastertable_viral_unrarefied_euk_virus_host$Final_family

# (A) aggregate 
# (A.1) families
Mastertable_viral_rarefied_categories_euk_fam <- aggregate(. ~FinalAnnotationFamily, FUN = sum, data = Mastertable_global_rarefied_euk_host[,colnames(Mastertable_global_rarefied_euk_host) %in% Mastertable_viral_rarefied_categories_vector | colnames(Mastertable_global_rarefied_euk_host) == 'FinalAnnotationFamily'])
rownames(Mastertable_viral_rarefied_categories_euk_fam) <- Mastertable_viral_rarefied_categories_euk_fam$FinalAnnotationFamily
Mastertable_viral_rarefied_categories_euk_fam$FinalAnnotationFamily <- NULL

# (A.2) genera
Mastertable_viral_rarefied_categories_euk_gen <- aggregate(. ~FinalAnnotationGenus, FUN = sum, data = Mastertable_global_rarefied_euk_host[,colnames(Mastertable_global_rarefied_euk_host) %in% Mastertable_viral_rarefied_categories_vector | colnames(Mastertable_global_rarefied_euk_host) == 'FinalAnnotationGenus'])
rownames(Mastertable_viral_rarefied_categories_euk_gen) <- Mastertable_viral_rarefied_categories_euk_gen$FinalAnnotationGenus
Mastertable_viral_rarefied_categories_euk_gen$FinalAnnotationGenus <- NULL

# (B) convert reads to absence ('0) or presence ('1) + add sample ID
# (B.1) families
Mastertable_viral_rarefied_categories_euk_fam_preabs <- ifelse(Mastertable_viral_rarefied_categories_euk_fam[, colnames(Mastertable_viral_rarefied_categories_euk_fam) %in% Mastertable_viral_rarefied_categories_vector] < 1, 0, 1)
Mastertable_viral_rarefied_categories_euk_fam_preabs <- rownames_to_column(as.data.frame(Mastertable_viral_rarefied_categories_euk_fam_preabs))
names(Mastertable_viral_rarefied_categories_euk_fam_preabs)[names(Mastertable_viral_rarefied_categories_euk_fam_preabs) == "rowname"] <- "sample"
rownames(Mastertable_viral_rarefied_categories_euk_fam_preabs) <- Mastertable_viral_rarefied_categories_euk_fam_preabs$sample

# (B.2) genera
Mastertable_viral_rarefied_categories_euk_gen_preabs <- ifelse(Mastertable_viral_rarefied_categories_euk_gen[, colnames(Mastertable_viral_rarefied_categories_euk_gen) %in% Mastertable_viral_rarefied_categories_vector] < 1, 0, 1)
Mastertable_viral_rarefied_categories_euk_gen_preabs <- rownames_to_column(as.data.frame(Mastertable_viral_rarefied_categories_euk_gen_preabs))
names(Mastertable_viral_rarefied_categories_euk_gen_preabs)[names(Mastertable_viral_rarefied_categories_euk_gen_preabs) == "rowname"] <- "sample"
rownames(Mastertable_viral_rarefied_categories_euk_gen_preabs) <- Mastertable_viral_rarefied_categories_euk_gen_preabs$sample

# (C) melt
# (C.1) families
Mastertable_viral_rarefied_categories_euk_fam_preabs_m <- melt(Mastertable_viral_rarefied_categories_euk_fam_preabs, id.vars = ("sample"))
Mastertable_viral_rarefied_categories_euk_fam_preabs_m$value[Mastertable_viral_rarefied_categories_euk_fam_preabs_m$value == "0"] <- "absent" 
Mastertable_viral_rarefied_categories_euk_fam_preabs_m$value[Mastertable_viral_rarefied_categories_euk_fam_preabs_m$value == "1"] <- "present" 

# (C.2) genera
Mastertable_viral_rarefied_categories_euk_gen_preabs_m <- melt(Mastertable_viral_rarefied_categories_euk_gen_preabs, id.vars = ("sample"))
Mastertable_viral_rarefied_categories_euk_gen_preabs_m$value[Mastertable_viral_rarefied_categories_euk_gen_preabs_m$value == "0"] <- "absent" 
Mastertable_viral_rarefied_categories_euk_gen_preabs_m$value[Mastertable_viral_rarefied_categories_euk_gen_preabs_m$value == "1"] <- "present" 

# (D) heatmap: presence/absence
# (D.1) families
ggplot(Mastertable_viral_rarefied_categories_euk_fam_preabs_m, aes(x=variable, y=sample)) +
  geom_tile(aes(fill=factor(value)), colour = "grey", show.legend = T) +
  scale_fill_manual(breaks = c("absent", "present"),values = c("white", "#225ea8")) +
  ggtitle("") +
  xlab("") +
  ylab("") +
  theme_minimal() +
  coord_equal() +
  theme(axis.text.y = element_text(color = "black", size = 5)) +
  scale_y_discrete(limits = unique(rev(Mastertable_viral_rarefied_categories_euk_fam_preabs_m$sample)), expand = c(0,0)) + 
  scale_x_discrete(expand = c(0,0)) +
  theme(panel.background = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.title = element_blank(), legend.text = element_text(size = 10), legend.key.size = unit(1,"line")) +
  theme(plot.title = element_text(hjust = -0.2 , vjust=2.12, size = 7, face ="bold")) +
  theme(axis.text.x = element_text(size = 5, angle = 45, vjust = 1 , hjust = 1, color = "black"), axis.title = element_text(size = 9)) 

# (D.2) genera
ggplot(Mastertable_viral_rarefied_categories_euk_gen_preabs_m, aes(x=variable, y=sample)) +
  geom_tile(aes(fill=factor(value)), colour = "grey", show.legend = T) +
  scale_fill_manual(breaks = c("absent", "present"),values = c("white", "#225ea8")) +
  ggtitle("") +
  xlab("") +
  ylab("") +
  theme_minimal() +
  coord_equal() +
  theme(axis.text.y = element_text(color = "black", size = 5)) +
  scale_y_discrete(limits = unique(rev(Mastertable_viral_rarefied_categories_euk_gen_preabs_m$sample)), expand = c(0,0)) + 
  scale_x_discrete(expand = c(0,0)) +
  theme(panel.background = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.title = element_blank(), legend.text = element_text(size = 10), legend.key.size = unit(1,"line")) +
  theme(plot.title = element_text(hjust = -0.2 , vjust=2.12, size = 7, face ="bold")) +
  theme(axis.text.x = element_text(size = 5, angle = 45, vjust = 1 , hjust = 1, color = "black"), axis.title = element_text(size = 9)) 
####################################
# 5.1.2 Samples with eukaryotic viral presence
####################################
Mastertable_viral_unrarefied_euk_virus_host
vector_1_1 <- (which(names(Mastertable_viral_unrarefied_euk_virus_host)== "Virsorter")-1)
Mastertable_viral_rarefied_categories_vector <- colnames(Mastertable_viral_unrarefied_euk_virus_host[,1:vector_1_1])

# (A) aggregate 
# (A.1) families
Mastertable_viral_rarefied_categories_euk_fam <- aggregate(. ~Final_family, FUN = sum, data = Mastertable_viral_unrarefied_euk_virus_host[,colnames(Mastertable_viral_unrarefied_euk_virus_host) %in% Mastertable_viral_rarefied_categories_vector | colnames(Mastertable_viral_unrarefied_euk_virus_host) == 'Final_family'])
rownames(Mastertable_viral_rarefied_categories_euk_fam) <- Mastertable_viral_rarefied_categories_euk_fam$Final_family
Mastertable_viral_rarefied_categories_euk_fam$Final_family <- NULL
View(Mastertable_viral_rarefied_categories_euk_fam)

# (A.2) genera
Mastertable_viral_rarefied_categories_euk_gen <- aggregate(. ~Final_genus, FUN = sum, data = Mastertable_viral_unrarefied_euk_virus_host[,colnames(Mastertable_viral_unrarefied_euk_virus_host) %in% Mastertable_viral_rarefied_categories_vector | colnames(Mastertable_viral_unrarefied_euk_virus_host) == 'Final_genus'])
rownames(Mastertable_viral_rarefied_categories_euk_gen) <- Mastertable_viral_rarefied_categories_euk_gen$Final_genus
Mastertable_viral_rarefied_categories_euk_gen$Final_genus <- NULL
View(Mastertable_viral_rarefied_categories_euk_gen)

# (B) convert reads to absence ('0) or presence ('1) + add sample ID
# (B.1) families
Mastertable_viral_rarefied_categories_euk_fam_preabs <- ifelse(Mastertable_viral_rarefied_categories_euk_fam[, colnames(Mastertable_viral_rarefied_categories_euk_fam) %in% Mastertable_viral_rarefied_categories_vector] < 1, 0, 1)
Mastertable_viral_rarefied_categories_euk_fam_preabs <- rownames_to_column(as.data.frame(Mastertable_viral_rarefied_categories_euk_fam_preabs))
names(Mastertable_viral_rarefied_categories_euk_fam_preabs)[names(Mastertable_viral_rarefied_categories_euk_fam_preabs) == "rowname"] <- "sample"
rownames(Mastertable_viral_rarefied_categories_euk_fam_preabs) <- Mastertable_viral_rarefied_categories_euk_fam_preabs$sample
View(Mastertable_viral_rarefied_categories_euk_fam_preabs)

# (B.2) genera
Mastertable_viral_rarefied_categories_euk_gen_preabs <- ifelse(Mastertable_viral_rarefied_categories_euk_gen[, colnames(Mastertable_viral_rarefied_categories_euk_gen) %in% Mastertable_viral_rarefied_categories_vector] < 1, 0, 1)
Mastertable_viral_rarefied_categories_euk_gen_preabs <- rownames_to_column(as.data.frame(Mastertable_viral_rarefied_categories_euk_gen_preabs))
names(Mastertable_viral_rarefied_categories_euk_gen_preabs)[names(Mastertable_viral_rarefied_categories_euk_gen_preabs) == "rowname"] <- "sample"
rownames(Mastertable_viral_rarefied_categories_euk_gen_preabs) <- Mastertable_viral_rarefied_categories_euk_gen_preabs$sample
View(Mastertable_viral_rarefied_categories_euk_gen_preabs)

# (C) melt & subset samples with presence of euk. viruses
# (C.1) families
Mastertable_viral_rarefied_categories_euk_fam_preabs_m <- melt(Mastertable_viral_rarefied_categories_euk_fam_preabs, id.vars = ("sample"))
Mastertable_viral_rarefied_categories_euk_fam_normalized_vector <- Mastertable_viral_rarefied_categories_euk_fam_preabs_m %>% filter(value > 0)
vector_euk <- unique(Mastertable_viral_rarefied_categories_euk_fam_normalized_vector$variable)
Mastertable_viral_rarefied_categories_euk_fam_preabs_m_subset <- Mastertable_viral_rarefied_categories_euk_fam_preabs_m[Mastertable_viral_rarefied_categories_euk_fam_preabs_m$variable %in% vector_euk,]
Mastertable_viral_rarefied_categories_euk_fam_preabs_m_subset$value[Mastertable_viral_rarefied_categories_euk_fam_preabs_m_subset$value == "0"] <- "absent" 
Mastertable_viral_rarefied_categories_euk_fam_preabs_m_subset$value[Mastertable_viral_rarefied_categories_euk_fam_preabs_m_subset$value == "1"] <- "present" 
length(vector_euk) # 123 samples with presence of at least one eukaryotic virus
View(Mastertable_viral_rarefied_categories_euk_fam_preabs_m_subset)

# (C.2) genera
Mastertable_viral_rarefied_categories_euk_gen_preabs_m <- melt(Mastertable_viral_rarefied_categories_euk_gen_preabs, id.vars = ("sample"))
Mastertable_viral_rarefied_categories_euk_gen_preabs_m_subset <- Mastertable_viral_rarefied_categories_euk_gen_preabs_m[Mastertable_viral_rarefied_categories_euk_gen_preabs_m$variable %in% vector_euk,]
Mastertable_viral_rarefied_categories_euk_gen_preabs_m_subset$value[Mastertable_viral_rarefied_categories_euk_gen_preabs_m_subset$value == "0"] <- "absent" 
Mastertable_viral_rarefied_categories_euk_gen_preabs_m_subset$value[Mastertable_viral_rarefied_categories_euk_gen_preabs_m_subset$value == "1"] <- "present" 
View(Mastertable_viral_rarefied_categories_euk_gen_preabs_m_subset)

# (D) heatmap: presence/absence
# (D.1) families
ggplot(Mastertable_viral_rarefied_categories_euk_fam_preabs_m_subset, aes(x=variable, y=sample)) +
  geom_tile(aes(fill=factor(value)), colour = "grey", show.legend = T) +
  scale_fill_manual(breaks = c("absent", "present"),values = c("white", "#225ea8")) +
  ggtitle("") +
  xlab("") +
  ylab("") +
  theme_minimal() +
  coord_equal() +
  theme(axis.text.y = element_text(color = "black", size = 5)) +
  scale_y_discrete(limits = unique(rev(Mastertable_viral_rarefied_categories_euk_fam_preabs_m_subset$sample)), expand = c(0,0)) + 
  scale_x_discrete(expand = c(0,0)) +
  theme(panel.background = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.title = element_blank(), legend.text = element_text(size = 10), legend.key.size = unit(1,"line")) +
  theme(plot.title = element_text(hjust = -0.2 , vjust=2.12, size = 7, face ="bold")) +
  theme(axis.text.x = element_text(size = 5, angle = 45, vjust = 1 , hjust = 1, color = "black"), axis.title = element_text(size = 9)) 

# (D.2) genera
ggplot(Mastertable_viral_rarefied_categories_euk_gen_preabs_m_subset, aes(x=variable, y=sample)) +
  geom_tile(aes(fill=factor(value)), colour = "grey", show.legend = T) +
  scale_fill_manual(breaks = c("absent", "present"),values = c("white", "#225ea8")) +
  ggtitle("") +
  xlab("") +
  ylab("") +
  theme_minimal() +
  coord_equal() +
  theme(axis.text.y = element_text(color = "black", size = 5)) +
  scale_y_discrete(limits = unique(rev(Mastertable_viral_rarefied_categories_euk_gen_preabs_m$sample)), expand = c(0,0)) + 
  scale_x_discrete(expand = c(0,0)) +
  theme(panel.background = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.title = element_blank(), legend.text = element_text(size = 10), legend.key.size = unit(1,"line")) +
  theme(plot.title = element_text(hjust = -0.2 , vjust=2.12, size = 7, face ="bold")) +
  theme(axis.text.x = element_text(size = 5, angle = 45, vjust = 1 , hjust = 1, color = "black"), axis.title = element_text(size = 9)) 
####################################
# 5.2 Heatmap (group host annotation)
####################################
# (1) ANIMAL VIRUSES

# (A) subset animal viruses and aggregate
# (A.1) families
Mastertable_global_rarefied_euk_fam_host_animal <- Mastertable_global_rarefied_euk_host[Mastertable_global_rarefied_euk_host$FinalAnnotationHost == "animal viruses",]
Mastertable_global_rarefied_euk_fam_host_animal <- aggregate(. ~FinalAnnotationFamily, FUN = sum, data = Mastertable_global_rarefied_euk_fam_host_animal[,colnames(Mastertable_global_rarefied_euk_fam_host_animal) %in% Mastertable_viral_rarefied_categories_vector | colnames(Mastertable_global_rarefied_euk_fam_host_animal) == 'FinalAnnotationFamily'])
rownames(Mastertable_global_rarefied_euk_fam_host_animal) <- Mastertable_global_rarefied_euk_fam_host_animal$FinalAnnotationFamily
Mastertable_global_rarefied_euk_fam_host_animal$FinalAnnotationFamily <- NULL

# (A.2) genera
Mastertable_global_rarefied_euk_gen_host_animal <- Mastertable_global_rarefied_euk_host[Mastertable_global_rarefied_euk_host$FinalAnnotationHost == "animal viruses",]
Mastertable_global_rarefied_euk_gen_host_animal <- aggregate(. ~FinalAnnotationGenus, FUN = sum, data = Mastertable_global_rarefied_euk_gen_host_animal[,colnames(Mastertable_global_rarefied_euk_gen_host_animal) %in% Mastertable_viral_rarefied_categories_vector | colnames(Mastertable_global_rarefied_euk_gen_host_animal) == 'FinalAnnotationGenus'])
rownames(Mastertable_global_rarefied_euk_gen_host_animal) <- Mastertable_global_rarefied_euk_gen_host_animal$FinalAnnotationGenus
Mastertable_global_rarefied_euk_gen_host_animal$FinalAnnotationGenus <- NULL

# (B) convert reads to absence ('0) or presence ('1) + add sample ID
# (B.1) families
Mastertable_global_rarefied_euk_host_animal_preabs <- ifelse(Mastertable_global_rarefied_euk_host_animal[, colnames(Mastertable_global_rarefied_euk_host_animal) %in% Mastertable_viral_rarefied_categories_vector] < 1, 0, 1)
Mastertable_global_rarefied_euk_host_animal_preabs <- rownames_to_column(as.data.frame(Mastertable_global_rarefied_euk_host_animal_preabs))
names(Mastertable_global_rarefied_euk_host_animal_preabs)[names(Mastertable_global_rarefied_euk_host_animal_preabs) == "rowname"] <- "sample"
rownames(Mastertable_global_rarefied_euk_host_animal_preabs) <- Mastertable_global_rarefied_euk_host_animal_preabs$sample

# (B.2) genera
Mastertable_global_rarefied_euk_gen_host_animal_preabs <- ifelse(Mastertable_global_rarefied_euk_gen_host_animal[, colnames(Mastertable_global_rarefied_euk_gen_host_animal) %in% Mastertable_viral_rarefied_categories_vector] < 1, 0, 1)
Mastertable_global_rarefied_euk_gen_host_animal_preabs <- rownames_to_column(as.data.frame(Mastertable_global_rarefied_euk_gen_host_animal_preabs))
names(Mastertable_global_rarefied_euk_gen_host_animal_preabs)[names(Mastertable_global_rarefied_euk_gen_host_animal_preabs) == "rowname"] <- "sample"
rownames(Mastertable_global_rarefied_euk_gen_host_animal_preabs) <- Mastertable_global_rarefied_euk_gen_host_animal_preabs$sample

# (C) melt & subset samples with presence of euk. viruses
# (C.1) families
Mastertable_viral_rarefied_categories_euk_fam_preabs_m <- melt(Mastertable_global_rarefied_euk_host_animal_preabs, id.vars = ("sample"))
Mastertable_viral_rarefied_categories_euk_fam_preabs_m_subset <- Mastertable_viral_rarefied_categories_euk_fam_preabs_m[Mastertable_viral_rarefied_categories_euk_fam_preabs_m$variable %in% vector_euk,]
Mastertable_viral_rarefied_categories_euk_fam_preabs_m_subset$value[Mastertable_viral_rarefied_categories_euk_fam_preabs_m_subset$value == "0"] <- "absent" 
Mastertable_viral_rarefied_categories_euk_fam_preabs_m_subset$value[Mastertable_viral_rarefied_categories_euk_fam_preabs_m_subset$value == "1"] <- "present" 

# (C.2) genera
Mastertable_viral_rarefied_categories_euk_gen_preabs_m <- melt(Mastertable_global_rarefied_euk_gen_host_animal_preabs, id.vars = ("sample"))
Mastertable_viral_rarefied_categories_euk_gen_preabs_m_subset <- Mastertable_viral_rarefied_categories_euk_gen_preabs_m[Mastertable_viral_rarefied_categories_euk_gen_preabs_m$variable %in% vector_euk,]
Mastertable_viral_rarefied_categories_euk_gen_preabs_m_subset$value[Mastertable_viral_rarefied_categories_euk_gen_preabs_m_subset$value == "0"] <- "absent" 
Mastertable_viral_rarefied_categories_euk_gen_preabs_m_subset$value[Mastertable_viral_rarefied_categories_euk_gen_preabs_m_subset$value == "1"] <- "present" 

# (D) heatmap animal viruses: presence/absence
animal_viruses_family <- ggplot(Mastertable_viral_rarefied_categories_euk_fam_preabs_m_subset, aes(x=variable, y=sample)) +
  geom_tile(aes(fill=factor(value)), colour = "grey", show.legend = F) +
  scale_fill_manual(breaks = c("absent", "present"),values = c("white", "#225ea8")) +
  ggtitle("") +
  xlab("") +
  ylab("") +
  theme_minimal() +
  coord_equal() +
  theme(axis.text.y = element_text(color = "black", size = 7)) +
  scale_y_discrete(limits = unique(rev(Mastertable_viral_rarefied_categories_euk_fam_preabs_m_subset$sample)), expand = c(0,0)) + 
  scale_x_discrete(expand = c(0,0)) +
  theme(panel.background = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.title = element_blank(), legend.text = element_text(size = 10), legend.key.size = unit(1,"line")) +
  theme(plot.title = element_text(hjust = -0.2 , vjust=2.12, size = 7, face ="bold")) +
  theme(axis.text.x = element_text(size = 5, angle = 45, vjust = 1 , hjust = 1, color = "black"), axis.title = element_text(size = 9)) 

## Kwe's script
ggplot(Mastertable_viral_rarefied_categories_euk_fam_preabs_m_subset, aes(x=variable, y=sample)) +
  geom_tile(aes(fill=factor(value)), colour = "grey", show.legend = F) +
  scale_fill_manual(breaks = c("absent", "present"),values = c("white", "#225ea8")) +
  ggtitle("") +
  xlab("") +
  ylab("viral family") +
  theme_grey(base_size = 5) +
  coord_equal() + 
  scale_x_discrete(expand = c(0,0)) +
  theme(axis.ticks = element_blank()) +
  scale_y_discrete(limits = unique(rev(Mastertable_viral_rarefied_categories_euk_fam_preabs_m_subset$sample)), expand = c(0,0)) + 
  theme(axis.text.x = element_text(size = 5   +  0.8, angle = 45, vjust = 1 , hjust = 1, color = "grey50"), axis.title = element_text(size = 9)) +
  theme(axis.text.y = element_text(color = "grey50", size = 8), axis.title.y = element_text(size=9)) +
  theme(legend.title = element_blank()) 

# (D.2) genera
animal_viruses_genus <- ggplot(Mastertable_viral_rarefied_categories_euk_gen_preabs_m_subset, aes(x=variable, y=sample)) +
  geom_tile(aes(fill=factor(value)), colour = "grey", show.legend = F) +
  scale_fill_manual(breaks = c("absent", "present"),values = c("white", "#225ea8")) +
  ggtitle("") +
  xlab("") +
  ylab("") +
  theme_minimal() +
  coord_equal() +
  theme(axis.text.y = element_text(color = "black", size = 7)) +
  scale_y_discrete(limits = unique(rev(Mastertable_viral_rarefied_categories_euk_gen_preabs_m_subset$sample)), expand = c(0,0)) + 
  scale_x_discrete(expand = c(0,0)) +
  theme(panel.background = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.title = element_blank(), legend.text = element_text(size = 10), legend.key.size = unit(1,"line")) +
  theme(plot.title = element_text(hjust = -0.2 , vjust=2.12, size = 7, face ="bold")) +
  theme(axis.text.x = element_text(size = 5, angle = 45, vjust = 1 , hjust = 1, color = "black"), axis.title = element_text(size = 9)) 

## kwe's script
ggplot(Mastertable_viral_rarefied_categories_euk_gen_preabs_m_subset, aes(x=variable, y=sample)) +
  geom_tile(aes(fill=factor(value)), colour = "grey", show.legend = F) +
  scale_fill_manual(breaks = c("absent", "present"),values = c("white", "#225ea8")) +
  ggtitle("") +
  xlab("") +
  ylab("viral genus") +
  theme_grey(base_size = 5) +
  coord_equal() + 
  scale_x_discrete(expand = c(0,0)) +
  theme(axis.ticks = element_blank()) +
  scale_y_discrete(limits = unique(rev(Mastertable_viral_rarefied_categories_euk_gen_preabs_m_subset$sample)), expand = c(0,0)) + 
  theme(axis.text.x = element_text(size = 5   +  0.8, angle = 45, vjust = 1 , hjust = 1, color = "grey50"), axis.title = element_text(size = 9)) +
  theme(axis.text.y = element_text(color = "grey50", size = 8), axis.title.y = element_text(size=9)) +
  theme(legend.title = element_blank()) 

# (2) PLANT AND FUNGAL VIRUSES
# (A) subset animal viruses and aggregate
# (A.1) families
Mastertable_global_rarefied_euk_fam_host_plant <- Mastertable_global_rarefied_euk_host[Mastertable_global_rarefied_euk_host$FinalAnnotationHost == "plant & fungal viruses",]
Mastertable_global_rarefied_euk_fam_host_plant <- aggregate(. ~FinalAnnotationFamily, FUN = sum, data = Mastertable_global_rarefied_euk_fam_host_plant[,colnames(Mastertable_global_rarefied_euk_fam_host_plant) %in% Mastertable_viral_rarefied_categories_vector | colnames(Mastertable_global_rarefied_euk_fam_host_plant) == 'FinalAnnotationFamily'])
rownames(Mastertable_global_rarefied_euk_fam_host_plant) <- Mastertable_global_rarefied_euk_fam_host_plant$FinalAnnotationFamily
Mastertable_global_rarefied_euk_fam_host_plant$FinalAnnotationFamily <- NULL

# (A.2) genera
Mastertable_global_rarefied_euk_gen_host_plant <- Mastertable_global_rarefied_euk_host[Mastertable_global_rarefied_euk_host$FinalAnnotationHost == "plant & fungal viruses",]
Mastertable_global_rarefied_euk_gen_host_plant <- aggregate(. ~FinalAnnotationGenus, FUN = sum, data = Mastertable_global_rarefied_euk_gen_host_plant[,colnames(Mastertable_global_rarefied_euk_gen_host_plant) %in% Mastertable_viral_rarefied_categories_vector | colnames(Mastertable_global_rarefied_euk_gen_host_plant) == 'FinalAnnotationGenus'])
rownames(Mastertable_global_rarefied_euk_gen_host_plant) <- Mastertable_global_rarefied_euk_gen_host_plant$FinalAnnotationGenus
Mastertable_global_rarefied_euk_gen_host_plant$FinalAnnotationGenus <- NULL

# (B) convert reads to absence ('0) or presence ('1) + add sample ID
# (B.1) families
Mastertable_global_rarefied_euk_host_plant_preabs <- ifelse(Mastertable_global_rarefied_euk_fam_host_plant[, colnames(Mastertable_global_rarefied_euk_fam_host_plant) %in% Mastertable_viral_rarefied_categories_vector] < 1, 0, 1)
Mastertable_global_rarefied_euk_host_plant_preabs <- rownames_to_column(as.data.frame(Mastertable_global_rarefied_euk_host_plant_preabs))
names(Mastertable_global_rarefied_euk_host_plant_preabs)[names(Mastertable_global_rarefied_euk_host_plant_preabs) == "rowname"] <- "sample"
rownames(Mastertable_global_rarefied_euk_host_plant_preabs) <- Mastertable_global_rarefied_euk_host_plant_preabs$sample

# (B.2) genera
Mastertable_global_rarefied_euk_gen_host_plant_preabs <- ifelse(Mastertable_global_rarefied_euk_gen_host_plant[, colnames(Mastertable_global_rarefied_euk_gen_host_plant) %in% Mastertable_viral_rarefied_categories_vector] < 1, 0, 1)
Mastertable_global_rarefied_euk_gen_host_plant_preabs <- rownames_to_column(as.data.frame(Mastertable_global_rarefied_euk_gen_host_plant_preabs))
names(Mastertable_global_rarefied_euk_gen_host_plant_preabs)[names(Mastertable_global_rarefied_euk_gen_host_plant_preabs) == "rowname"] <- "sample"
rownames(Mastertable_global_rarefied_euk_gen_host_plant_preabs) <- Mastertable_global_rarefied_euk_gen_host_plant_preabs$sample

# (C) melt & subset samples with presence of euk. viruses
# (C.1) families
Mastertable_viral_rarefied_categories_euk_fam_preabs_m <- melt(Mastertable_global_rarefied_euk_host_plant_preabs, id.vars = ("sample"))
Mastertable_viral_rarefied_categories_euk_fam_preabs_m_subset <- Mastertable_viral_rarefied_categories_euk_fam_preabs_m[Mastertable_viral_rarefied_categories_euk_fam_preabs_m$variable %in% vector_euk,]
Mastertable_viral_rarefied_categories_euk_fam_preabs_m_subset$value[Mastertable_viral_rarefied_categories_euk_fam_preabs_m_subset$value == "0"] <- "absent" 
Mastertable_viral_rarefied_categories_euk_fam_preabs_m_subset$value[Mastertable_viral_rarefied_categories_euk_fam_preabs_m_subset$value == "1"] <- "present" 

# (C.2) genera
Mastertable_viral_rarefied_categories_euk_gen_preabs_m <- melt(Mastertable_global_rarefied_euk_gen_host_plant_preabs, id.vars = ("sample"))
Mastertable_viral_rarefied_categories_euk_gen_preabs_m_subset <- Mastertable_viral_rarefied_categories_euk_gen_preabs_m[Mastertable_viral_rarefied_categories_euk_gen_preabs_m$variable %in% vector_euk,]
Mastertable_viral_rarefied_categories_euk_gen_preabs_m_subset$value[Mastertable_viral_rarefied_categories_euk_gen_preabs_m_subset$value == "0"] <- "absent" 
Mastertable_viral_rarefied_categories_euk_gen_preabs_m_subset$value[Mastertable_viral_rarefied_categories_euk_gen_preabs_m_subset$value == "1"] <- "present" 

# (D) heatmap plant & fungal: presence/absence
# (D.1) families
plant_fungal_viruses_family <- ggplot(Mastertable_viral_rarefied_categories_euk_fam_preabs_m_subset, aes(x=variable, y=sample)) +
  geom_tile(aes(fill=factor(value)), colour = "grey", show.legend = F) +
  scale_fill_manual(breaks = c("absent", "present"),values = c("white", "#225ea8")) +
  ggtitle("") +
  xlab("") +
  ylab("") +
  theme_minimal() +
  coord_equal() +
  theme(axis.text.y = element_text(color = "black", size = 7)) +
  scale_y_discrete(limits = unique(rev(Mastertable_viral_rarefied_categories_euk_fam_preabs_m_subset$sample)), expand = c(0,0)) + 
  scale_x_discrete(expand = c(0,0)) +
  theme(panel.background = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.title = element_blank(), legend.text = element_text(size = 10), legend.key.size = unit(1,"line")) +
  theme(plot.title = element_text(hjust = -0.2 , vjust=2.12, size = 7, face ="bold")) +
  theme(axis.text.x = element_text(size = 5, angle = 45, vjust = 1 , hjust = 1, color = "black"), axis.title = element_text(size = 9)) 

# (D.2) genera
plant_fungal_viruses_genus <- ggplot(Mastertable_viral_rarefied_categories_euk_gen_preabs_m_subset, aes(x=variable, y=sample)) +
  geom_tile(aes(fill=factor(value)), colour = "grey", show.legend = F) +
  scale_fill_manual(breaks = c("absent", "present"),values = c("white", "#225ea8")) +
  ggtitle("") +
  xlab("") +
  ylab("") +
  theme_minimal() +
  coord_equal() +
  theme(axis.text.y = element_text(color = "black", size = 7)) +
  scale_y_discrete(limits = unique(rev(Mastertable_viral_rarefied_categories_euk_gen_preabs_m$sample)), expand = c(0,0)) + 
  scale_x_discrete(expand = c(0,0)) +
  theme(panel.background = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.title = element_blank(), legend.text = element_text(size = 10), legend.key.size = unit(1,"line")) +
  theme(plot.title = element_text(hjust = -0.2 , vjust=2.12, size = 7, face ="bold")) +
  theme(axis.text.x = element_text(size = 5, angle = 45, vjust = 1 , hjust = 1, color = "black"), axis.title = element_text(size = 9)) 

# (3) SMALL CIRCULAR VIRUSES
# (A) subset small circular viruses and aggregate
# (A.1) families
Mastertable_global_rarefied_euk_fam_host_scv <- Mastertable_global_rarefied_euk_host[Mastertable_global_rarefied_euk_host$FinalAnnotationHost == "small circular viruses",]
Mastertable_global_rarefied_euk_fam_host_scv <- aggregate(. ~FinalAnnotationFamily, FUN = sum, data = Mastertable_global_rarefied_euk_fam_host_scv[,colnames(Mastertable_global_rarefied_euk_fam_host_scv) %in% Mastertable_viral_rarefied_categories_vector | colnames(Mastertable_global_rarefied_euk_fam_host_scv) == 'FinalAnnotationFamily'])
rownames(Mastertable_global_rarefied_euk_fam_host_scv) <- Mastertable_global_rarefied_euk_fam_host_scv$FinalAnnotationFamily
Mastertable_global_rarefied_euk_fam_host_scv$FinalAnnotationFamily <- NULL

# (A.2) genera
Mastertable_global_rarefied_euk_gen_host_scv <- Mastertable_global_rarefied_euk_host[Mastertable_global_rarefied_euk_host$FinalAnnotationHost == "small circular viruses",]
Mastertable_global_rarefied_euk_gen_host_scv <- aggregate(. ~FinalAnnotationGenus, FUN = sum, data = Mastertable_global_rarefied_euk_gen_host_scv[,colnames(Mastertable_global_rarefied_euk_gen_host_scv) %in% Mastertable_viral_rarefied_categories_vector | colnames(Mastertable_global_rarefied_euk_gen_host_scv) == 'FinalAnnotationGenus'])
rownames(Mastertable_global_rarefied_euk_gen_host_scv) <- Mastertable_global_rarefied_euk_gen_host_scv$FinalAnnotationGenus
Mastertable_global_rarefied_euk_gen_host_scv$FinalAnnotationGenus <- NULL

# (B) convert reads to absence ('0) or presence ('1) + add sample ID
# (B.1) families
Mastertable_global_rarefied_euk_fam_host_scv_preabs <- ifelse(Mastertable_global_rarefied_euk_fam_host_scv[, colnames(Mastertable_global_rarefied_euk_fam_host_scv) %in% Mastertable_viral_rarefied_categories_vector] < 1, 0, 1)
Mastertable_global_rarefied_euk_fam_host_scv_preabs <- rownames_to_column(as.data.frame(Mastertable_global_rarefied_euk_fam_host_scv_preabs))
names(Mastertable_global_rarefied_euk_fam_host_scv_preabs)[names(Mastertable_global_rarefied_euk_fam_host_scv_preabs) == "rowname"] <- "sample"
rownames(Mastertable_global_rarefied_euk_fam_host_scv_preabs) <- Mastertable_global_rarefied_euk_fam_host_scv_preabs$sample

# (B.2) genera
Mastertable_global_rarefied_euk_gen_host_scv_preabs <- ifelse(Mastertable_global_rarefied_euk_gen_host_scv[, colnames(Mastertable_global_rarefied_euk_gen_host_scv) %in% Mastertable_viral_rarefied_categories_vector] < 1, 0, 1)
Mastertable_global_rarefied_euk_gen_host_scv_preabs <- rownames_to_column(as.data.frame(Mastertable_global_rarefied_euk_gen_host_scv_preabs))
names(Mastertable_global_rarefied_euk_gen_host_scv_preabs)[names(Mastertable_global_rarefied_euk_gen_host_scv_preabs) == "rowname"] <- "sample"
rownames(Mastertable_global_rarefied_euk_gen_host_scv_preabs) <- Mastertable_global_rarefied_euk_gen_host_scv_preabs$sample

# (C) melt & subset samples with presence of euk. viruses
# (C.1) families
Mastertable_viral_rarefied_categories_euk_fam_preabs_m <- melt(Mastertable_global_rarefied_euk_fam_host_scv_preabs, id.vars = ("sample"))
Mastertable_viral_rarefied_categories_euk_fam_preabs_m_subset <- Mastertable_viral_rarefied_categories_euk_fam_preabs_m[Mastertable_viral_rarefied_categories_euk_fam_preabs_m$variable %in% vector_euk,]
Mastertable_viral_rarefied_categories_euk_fam_preabs_m_subset$value[Mastertable_viral_rarefied_categories_euk_fam_preabs_m_subset$value == "0"] <- "absent" 
Mastertable_viral_rarefied_categories_euk_fam_preabs_m_subset$value[Mastertable_viral_rarefied_categories_euk_fam_preabs_m_subset$value == "1"] <- "present" 

# (C.2) genera
Mastertable_viral_rarefied_categories_euk_gen_preabs_m <- melt(Mastertable_global_rarefied_euk_gen_host_scv_preabs, id.vars = ("sample"))
Mastertable_viral_rarefied_categories_euk_gen_preabs_m_subset <- Mastertable_viral_rarefied_categories_euk_gen_preabs_m[Mastertable_viral_rarefied_categories_euk_gen_preabs_m$variable %in% vector_euk,]
Mastertable_viral_rarefied_categories_euk_gen_preabs_m_subset$value[Mastertable_viral_rarefied_categories_euk_gen_preabs_m_subset$value == "0"] <- "absent" 
Mastertable_viral_rarefied_categories_euk_gen_preabs_m_subset$value[Mastertable_viral_rarefied_categories_euk_gen_preabs_m_subset$value == "1"] <- "present" 

# (D) heatmap small circular viruses: presence/absence
# (D.1) families
small_circular_viruses_family <- ggplot(Mastertable_viral_rarefied_categories_euk_fam_preabs_m_subset, aes(x=variable, y=sample)) +
  geom_tile(aes(fill=factor(value)), colour = "grey", show.legend = F) +
  scale_fill_manual(breaks = c("absent", "present"),values = c("white", "#225ea8")) +
  ggtitle("") +
  xlab("") +
  ylab("") +
  theme_minimal() +
  coord_equal() +
  theme(axis.text.y = element_text(color = "black", size = 7)) +
  scale_y_discrete(limits = unique(rev(Mastertable_viral_rarefied_categories_euk_fam_preabs_m_subset$sample)), expand = c(0,0)) + 
  scale_x_discrete(expand = c(0,0)) +
  theme(panel.background = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.title = element_blank(), legend.text = element_text(size = 10), legend.key.size = unit(1,"line")) +
  theme(plot.title = element_text(hjust = -0.2 , vjust=2.12, size = 7, face ="bold")) +
  theme(axis.text.x = element_text(size = 5, angle = 45, vjust = 1 , hjust = 1, color = "black"), axis.title = element_text(size = 9)) 

# (D.2) genera
small_circular_viruses_genus <- ggplot(Mastertable_viral_rarefied_categories_euk_gen_preabs_m_subset, aes(x=variable, y=sample)) +
  geom_tile(aes(fill=factor(value)), colour = "grey", show.legend = F) +
  scale_fill_manual(breaks = c("absent", "present"),values = c("white", "#225ea8")) +
  ggtitle("") +
  xlab("") +
  ylab("") +
  theme_minimal() +
  coord_equal() +
  theme(axis.text.y = element_text(color = "black", size = 7)) +
  scale_y_discrete(limits = unique(rev(Mastertable_viral_rarefied_categories_euk_gen_preabs_m$sample)), expand = c(0,0)) + 
  scale_x_discrete(expand = c(0,0)) +
  theme(panel.background = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.title = element_blank(), legend.text = element_text(size = 10), legend.key.size = unit(1,"line")) +
  theme(plot.title = element_text(hjust = -0.2 , vjust=2.12, size = 7, face ="bold")) +
  theme(axis.text.x = element_text(size = 5, angle = 45, vjust = 1 , hjust = 1, color = "black"), axis.title = element_text(size = 9)) 

# (4) UNCLASSIFIEDVIRUSES
# (A) subset small circular viruses and aggregate
# (A.1) families
Mastertable_global_rarefied_euk_fam_host_unc <- Mastertable_global_rarefied_euk_host[Mastertable_global_rarefied_euk_host$FinalAnnotationHost == "unclassified viruses",]
Mastertable_global_rarefied_euk_fam_host_unc <- aggregate(. ~FinalAnnotationFamily, FUN = sum, data = Mastertable_global_rarefied_euk_fam_host_unc[,colnames(Mastertable_global_rarefied_euk_fam_host_unc) %in% Mastertable_viral_rarefied_categories_vector | colnames(Mastertable_global_rarefied_euk_fam_host_unc) == 'FinalAnnotationFamily'])
rownames(Mastertable_global_rarefied_euk_fam_host_unc) <- Mastertable_global_rarefied_euk_fam_host_unc$FinalAnnotationFamily
Mastertable_global_rarefied_euk_fam_host_unc$FinalAnnotationFamily <- NULL

# (A.2) genera
Mastertable_global_rarefied_euk_gen_host_unc <- Mastertable_global_rarefied_euk_host[Mastertable_global_rarefied_euk_host$FinalAnnotationHost == "unclassified viruses",]
Mastertable_global_rarefied_euk_gen_host_unc <- aggregate(. ~FinalAnnotationGenus, FUN = sum, data = Mastertable_global_rarefied_euk_gen_host_unc[,colnames(Mastertable_global_rarefied_euk_gen_host_unc) %in% Mastertable_viral_rarefied_categories_vector | colnames(Mastertable_global_rarefied_euk_gen_host_unc) == 'FinalAnnotationGenus'])
rownames(Mastertable_global_rarefied_euk_gen_host_unc) <- Mastertable_global_rarefied_euk_gen_host_unc$FinalAnnotationGenus
Mastertable_global_rarefied_euk_gen_host_unc$FinalAnnotationGenus <- NULL

# (B) convert reads to absence ('0) or presence ('1) + add sample ID
# (B.1) families
Mastertable_global_rarefied_euk_fam_host_unc_preabs <- ifelse(Mastertable_global_rarefied_euk_fam_host_unc[, colnames(Mastertable_global_rarefied_euk_fam_host_unc) %in% Mastertable_viral_rarefied_categories_vector] < 1, 0, 1)
Mastertable_global_rarefied_euk_fam_host_unc_preabs <- rownames_to_column(as.data.frame(Mastertable_global_rarefied_euk_fam_host_unc_preabs))
names(Mastertable_global_rarefied_euk_fam_host_unc_preabs)[names(Mastertable_global_rarefied_euk_fam_host_unc_preabs) == "rowname"] <- "sample"
rownames(Mastertable_global_rarefied_euk_fam_host_unc_preabs) <- Mastertable_global_rarefied_euk_fam_host_unc_preabs$sample

# (B.2) genera
Mastertable_global_rarefied_euk_gen_host_unc_preabs <- ifelse(Mastertable_global_rarefied_euk_gen_host_unc[, colnames(Mastertable_global_rarefied_euk_gen_host_unc) %in% Mastertable_viral_rarefied_categories_vector] < 1, 0, 1)
Mastertable_global_rarefied_euk_gen_host_unc_preabs <- rownames_to_column(as.data.frame(Mastertable_global_rarefied_euk_gen_host_unc_preabs))
names(Mastertable_global_rarefied_euk_gen_host_unc_preabs)[names(Mastertable_global_rarefied_euk_gen_host_unc_preabs) == "rowname"] <- "sample"
rownames(Mastertable_global_rarefied_euk_gen_host_unc_preabs) <- Mastertable_global_rarefied_euk_gen_host_unc_preabs$sample

# (C) melt & subset samples with presence of euk. viruses
# (C.1) families
Mastertable_viral_rarefied_categories_euk_fam_preabs_m <- melt(Mastertable_global_rarefied_euk_fam_host_unc_preabs, id.vars = ("sample"))
Mastertable_viral_rarefied_categories_euk_fam_preabs_m_subset <- Mastertable_viral_rarefied_categories_euk_fam_preabs_m[Mastertable_viral_rarefied_categories_euk_fam_preabs_m$variable %in% vector_euk,]
Mastertable_viral_rarefied_categories_euk_fam_preabs_m_subset$value[Mastertable_viral_rarefied_categories_euk_fam_preabs_m_subset$value == "0"] <- "absent" 
Mastertable_viral_rarefied_categories_euk_fam_preabs_m_subset$value[Mastertable_viral_rarefied_categories_euk_fam_preabs_m_subset$value == "1"] <- "present" 
# adapt long names so format doesn't change
Mastertable_viral_rarefied_categories_euk_fam_preabs_m_subset$sample[Mastertable_viral_rarefied_categories_euk_fam_preabs_m_subset$sample == "unclassified eukaryotic virus"] <- "unclassified virus"

# (C.2) genera
Mastertable_viral_rarefied_categories_euk_gen_preabs_m <- melt(Mastertable_global_rarefied_euk_gen_host_unc_preabs, id.vars = ("sample"))
Mastertable_viral_rarefied_categories_euk_gen_preabs_m_subset <- Mastertable_viral_rarefied_categories_euk_gen_preabs_m[Mastertable_viral_rarefied_categories_euk_gen_preabs_m$variable %in% vector_euk,]
Mastertable_viral_rarefied_categories_euk_gen_preabs_m_subset$value[Mastertable_viral_rarefied_categories_euk_gen_preabs_m_subset$value == "0"] <- "absent" 
Mastertable_viral_rarefied_categories_euk_gen_preabs_m_subset$value[Mastertable_viral_rarefied_categories_euk_gen_preabs_m_subset$value == "1"] <- "present" 
# adapt long names so format doesn't change
Mastertable_viral_rarefied_categories_euk_gen_preabs_m_subset$sample[Mastertable_viral_rarefied_categories_euk_gen_preabs_m_subset$sample == "Human lung-associated vientovirus LT"] <- "unclassified vientovirus" 

# (D) heatmap unclassified viruses: presence/absence
# (D.1) families
unclassified_viruses_family <- ggplot(Mastertable_viral_rarefied_categories_euk_fam_preabs_m_subset, aes(x=variable, y=sample)) +
  geom_tile(aes(fill=factor(value)), colour = "grey", show.legend = F) +
  scale_fill_manual(breaks = c("absent", "present"),values = c("white", "#225ea8")) +
  ggtitle("") +
  xlab("") +
  ylab("") +
  theme_minimal() +
  coord_equal() +
  theme(axis.text.y = element_text(color = "black", size = 7)) +
  scale_y_discrete(limits = unique(rev(Mastertable_viral_rarefied_categories_euk_fam_preabs_m_subset$sample)), expand = c(0,0)) + 
  scale_x_discrete(expand = c(0,0)) +
  theme(panel.background = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.title = element_blank(), legend.text = element_text(size = 10), legend.key.size = unit(1,"line")) +
  theme(plot.title = element_text(hjust = -0.2 , vjust=2.12, size = 7, face ="bold")) +
  theme(axis.text.x = element_text(size = 5, angle = 45, vjust = 1 , hjust = 1, color = "black"), axis.title = element_text(size = 9)) 

# (D.2) genera
unclassified_viruses_genus <- ggplot(Mastertable_viral_rarefied_categories_euk_gen_preabs_m_subset, aes(x=variable, y=sample)) +
  geom_tile(aes(fill=factor(value)), colour = "grey", show.legend = F) +
  scale_fill_manual(breaks = c("absent", "present"),values = c("white", "#225ea8")) +
  ggtitle("") +
  xlab("") +
  ylab("") +
  theme_minimal() +
  coord_equal() +
  theme(axis.text.y = element_text(color = "black", size = 7)) +
  scale_y_discrete(limits = unique(rev(Mastertable_viral_rarefied_categories_euk_gen_preabs_m_subset$sample)), expand = c(0,0)) + 
  scale_x_discrete(expand = c(0,0)) +
  theme(panel.background = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.title = element_blank(), legend.text = element_text(size = 10), legend.key.size = unit(1,"line")) +
  theme(plot.title = element_text(hjust = -0.2 , vjust=2.12, size = 7, face ="bold")) +
  theme(axis.text.x = element_text(size = 5, angle = 45, vjust = 1 , hjust = 1, color = "black"), axis.title = element_text(size = 9)) 

# (5) legend extraction
## SAVE plot only for the legend
ggplot(Mastertable_viral_rarefied_categories_euk_gen_preabs_m_subset, aes(x=variable, y=sample)) +
  geom_tile(aes(fill=factor(value)), colour = "grey", show.legend = T) +
  scale_fill_manual(breaks = c("absent", "present"),values = c("white", "#225ea8")) +
  ggtitle("") +
  xlab("") +
  ylab("") +
  theme_minimal() +
  coord_equal() +
  theme(axis.text.y = element_text(color = "black", size = 7)) +
  scale_y_discrete(limits = unique(rev(Mastertable_viral_rarefied_categories_euk_gen_preabs_m_subset$sample)), expand = c(0,0)) + 
  scale_x_discrete(expand = c(0,0)) +
  theme(panel.background = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.title = element_blank(), legend.text = element_text(size = 10), legend.key.size = unit(2,"line")) +
  theme(plot.title = element_text(hjust = -0.2 , vjust=2.12, size = 7, face ="bold")) +
  theme(axis.text.x = element_text(size = 5, angle = 45, vjust = 1 , hjust = 1, color = "black"), axis.title = element_text(size = 9))

# save these plots
animal_viruses_family
plant_fungal_viruses_family
small_circular_viruses_family
unclassified_viruses_family

animal_viruses_genus
plant_fungal_viruses_genus
small_circular_viruses_genus
unclassified_viruses_genus

dev.off()

# (7) One figure containing all heatmaps
# (7.1) We will plot the previous plots all together in one, it is not necessary that all plots have a legend and title so we will remove some
# families
animal_viruses_family_no_leg <- animal_viruses_family + theme(legend.position='none',  axis.title.y=element_blank(),  axis.text.y=element_blank(), axis.text.x = element_blank())
plant_fungal_viruses_family_no_leg <- plant_fungal_viruses_family + theme(legend.position='none', axis.title.y=element_blank(),  axis.text.y=element_blank(),axis.text.x = element_blank())
small_circular_viruses_family_no_leg <- small_circular_viruses_family + theme(legend.position='none', axis.title.y=element_blank(), axis.text.y=element_blank(),axis.text.x = element_blank())
unclassified_viruses_family_no_leg <- unclassified_viruses_family + theme(legend.position='none', axis.title.y=element_blank(), axis.text.y=element_blank())

# genera
animal_viruses_genus_no_leg <- animal_viruses_genus + theme(legend.position='none',  axis.title.y=element_blank(),  axis.text.y=element_blank(),axis.text.x = element_blank())
plant_fungal_viruses_genus_no_leg <- plant_fungal_viruses_genus + theme(legend.position='none', axis.title.y=element_blank(),  axis.text.y=element_blank(),axis.text.x = element_blank())
small_circular_viruses_genus_no_leg <- small_circular_viruses_genus + theme(legend.position='none', axis.title.y=element_blank(), axis.text.y=element_blank(),axis.text.x = element_blank())
unclassified_viruses_genus_no_leg <- unclassified_viruses_genus + theme(legend.position='none', axis.title.y=element_blank(), axis.text.y=element_blank())

animal_viruses_genus_no_leg
plant_fungal_viruses_genus_no_leg
small_circular_viruses_genus_no_leg
unclassified_viruses_genus_no_leg

### you can better save these here & add all together in Illustrator Adobe!

# (7.2) Define layout for the plots (5 rows, 1 column)
# families
layt_fam <-grid.layout(nrow=4,ncol=1,heights = c(1, 1.79, 1, 1)) 

# genera
layt_gen <-grid.layout(nrow=4,ncol=1,heights =c(2,3,1,1.2)) 

# (7.3) view the layout of plots
# families
grid.show.layout(layt_fam)

# genera
grid.show.layout(layt_gen)

# (7.4) Draw plots one by one in their positions
# families
grid.newpage()
pushViewport(viewport(layout=layt_fam))
print(animal_viruses_family_no_leg,vp=viewport(layout.pos.row=1,layout.pos.col=1))
print(plant_fungal_viruses_family_no_leg,vp=viewport(layout.pos.row=2,layout.pos.col=1))
print(small_circular_viruses_family_no_leg,vp=viewport(layout.pos.row=3,layout.pos.col=1))
print(unclassified_viruses_family_no_leg,vp=viewport(layout.pos.row=4,layout.pos.col=1))
dev.off()

# genera
grid.newpage()
pushViewport(viewport(layout=layt_gen))
print(animal_viruses_genus_no_leg,vp=viewport(layout.pos.row=1,layout.pos.col=1))
print(plant_fungal_viruses_genus_no_leg,vp=viewport(layout.pos.row=2,layout.pos.col=1))
print(small_circular_viruses_genus_no_leg,vp=viewport(layout.pos.row=3,layout.pos.col=1))
print(unclassified_viruses_genus_no_leg,vp=viewport(layout.pos.row=4,layout.pos.col=1))
dev.off()

# (7) finish figure in Illustrator 
#   ---> add titles there (animal viruses, plant & fungal viruses, small circular viruses, unclassified viruses)
#   ---> add y-axis names (families in italic)
#   ---> add legend there
#   ---> add long names that change format such as "human lung-associated vientovirus LT" (for now change to unclassified vientovirus)
#   ---> change "unclassified eukaryotic virus", to "unclassified virus"
####################################
# 5.3 Complex heatmap: hierarchical clustering based on euclidean/pearson
####################################
# Hierarchical clustering you have to evaluate what the best way to di it. Since Pearson thrives on linear connections, I'm not sure if this is the best option.
####################################
# 5.3.1 preparation of matrix
####################################
# info: https://jokergoo.github.io/ComplexHeatmap-reference/book/introduction.html
# https://bioconductor.org/packages/release/bioc/vignettes/ComplexHeatmap/inst/doc/most_probably_asked_questions.html
# check also Chenyen's thesis

# (1) Subset of samples containg eukaryotic viruses: 
#               --> single heatmap

Mastertable_viral_rarefied_euk_virus_host$Final_genus[Mastertable_viral_rarefied_euk_virus_host$Final_species == "Torque teno virus"] <- "Alphatorquevirus"
Mastertable_viral_rarefied_euk_virus_host$Final_genus[Mastertable_viral_rarefied_euk_virus_host$Final_species == "Chicken genomovirus mg2_75"] <- "unclassified Genomoviridae"

# (A) Create an abundance table: groups as rownames for y-axis & sample ID as colnames for x-axis
# (A.1) families
View(Mastertable_viral_rarefied_euk_virus_host)
table(Mastertable_viral_rarefied_euk_virus_host$Final_family)
table(Mastertable_viral_rarefied_euk_virus_host$Final_genus)

Mastertable_global_rarefied_fam_complex_heatm <- Mastertable_viral_rarefied_euk_virus_host
Mastertable_global_rarefied_fam_complex_heatm <- aggregate(. ~Final_family, FUN = sum, data = Mastertable_viral_rarefied_euk_virus_host[,colnames(Mastertable_viral_rarefied_euk_virus_host) %in% vector_euk | colnames(Mastertable_viral_rarefied_euk_virus_host) == 'Final_family'])
rownames(Mastertable_global_rarefied_fam_complex_heatm) <- Mastertable_global_rarefied_fam_complex_heatm$Final_family
Mastertable_global_rarefied_fam_complex_heatm$Final_family <- NULL
View(Mastertable_global_rarefied_fam_complex_heatm)

# (A.2) genera
Mastertable_global_rarefied_gen_complex_heatm <- Mastertable_viral_rarefied_euk_virus_host
Mastertable_global_rarefied_gen_complex_heatm <- aggregate(. ~Final_genus, FUN = sum, data = Mastertable_global_rarefied_gen_complex_heatm[,colnames(Mastertable_global_rarefied_gen_complex_heatm) %in% vector_euk | colnames(Mastertable_global_rarefied_gen_complex_heatm) == 'Final_genus'])
rownames(Mastertable_global_rarefied_gen_complex_heatm) <- Mastertable_global_rarefied_gen_complex_heatm$Final_genus
Mastertable_global_rarefied_gen_complex_heatm$Final_genus <- NULL
View(Mastertable_global_rarefied_gen_complex_heatm)

# (A.3) species
Mastertable_global_rarefied_species_complex_heatm <- Mastertable_viral_rarefied_euk_virus_host
Mastertable_global_rarefied_species_complex_heatm <- aggregate(. ~Final_species, FUN = sum, data = Mastertable_global_rarefied_species_complex_heatm[,colnames(Mastertable_global_rarefied_species_complex_heatm) %in% vector_euk | colnames(Mastertable_global_rarefied_species_complex_heatm) == 'Final_species'])
rownames(Mastertable_global_rarefied_species_complex_heatm) <- Mastertable_global_rarefied_species_complex_heatm$Final_species
Mastertable_global_rarefied_species_complex_heatm$Final_species <- NULL
View(Mastertable_global_rarefied_species_complex_heatm)

# (B) convert reads to relative abundances + add sample ID
# (B.1) families
Mastertable_global_rarefied_fam_complex_heatm <- sweep(Mastertable_global_rarefied_fam_complex_heatm, 2, colSums(Mastertable_global_rarefied_fam_complex_heatm), '/') 
Mastertable_global_rarefied_fam_complex_heatm[is.na(Mastertable_global_rarefied_fam_complex_heatm)] <- 0
matrix_family <- as.matrix(Mastertable_global_rarefied_fam_complex_heatm)
View(matrix_family)

# (B.2) genera
Mastertable_global_rarefied_gen_complex_heatm <- sweep(Mastertable_global_rarefied_gen_complex_heatm, 2, colSums(Mastertable_global_rarefied_gen_complex_heatm), '/') 
Mastertable_global_rarefied_gen_complex_heatm[is.na(Mastertable_global_rarefied_gen_complex_heatm)] <- 0
matrix_genus <- as.matrix(Mastertable_global_rarefied_gen_complex_heatm)
View(matrix_genus)

# (B.3) species
Mastertable_global_rarefied_species_complex_heatm <- sweep(Mastertable_global_rarefied_species_complex_heatm, 2, colSums(Mastertable_global_rarefied_species_complex_heatm), '/') 
Mastertable_global_rarefied_species_complex_heatm[is.na(Mastertable_global_rarefied_species_complex_heatm)] <- 0
matrix_species <- as.matrix(Mastertable_global_rarefied_species_complex_heatm)
View(matrix_species)

# (D) Create list for side annotations and create complex heatmap
####################################
# 5.3.2 families
####################################
## (i) Create a percentage of the presence of a specific eukaryotic viral family over all the samples - by this measure you can see which families is the most prevalent in the samples. This is obviously based on all the column & not only on this subset.
#           --> I did a few extra step because I thought I needed the order for the sample_ID's (x-axis) after dendrogram, but this is not needed

View(Mastertable_viral_rarefied_euk_virus_host)

Mastertable_global_rarefied_fam_complex_heatm <- Mastertable_viral_rarefied_euk_virus_host
used_for_heatmap_later <- aggregate(. ~Final_family, FUN = sum, data = Mastertable_global_rarefied_fam_complex_heatm[,colnames(Mastertable_global_rarefied_fam_complex_heatm) %in% vector_euk | colnames(Mastertable_global_rarefied_fam_complex_heatm) == 'Final_family'])
Mastertable_global_rarefied_fam_complex_heatm <- aggregate(. ~Final_family, FUN = sum, data = Mastertable_global_rarefied_fam_complex_heatm[,colnames(Mastertable_global_rarefied_fam_complex_heatm) %in% vector_euk | colnames(Mastertable_global_rarefied_fam_complex_heatm) == 'Final_family'])
rownames(Mastertable_global_rarefied_fam_complex_heatm) <- Mastertable_global_rarefied_fam_complex_heatm$Final_family
Mastertable_global_rarefied_fam_complex_heatm$Final_family <- NULL
rownames(used_for_heatmap_later) <- used_for_heatmap_later$Final_family
used_for_heatmap_later$Final_family <- NULL

vector_eukfamily_names_row <- rownames(Mastertable_global_rarefied_fam_complex_heatm)

percentage_families_row <- vector()  
for (l in 1: length(vector_eukfamily_names_row)){
  percentage_families_row[l] <- ((length(which(Mastertable_global_rarefied_fam_complex_heatm[l,] > 0)))/ncol(used_for_heatmap_later))*100
}

percentage_families_row

## (ii) Create a percentage of the number of different eukaryotic viral families over all families found within one sample
family_names_column <- colnames(matrix_family)
Mastertable_global_rarefied_fam_complex_heatm_dendo_reordered <- Mastertable_global_rarefied_fam_complex_heatm

richness_families_column <- vector()  
for (l in 1: length(family_names_column)){
  richness_families_column[l] <- (length(which(Mastertable_global_rarefied_fam_complex_heatm_dendo_reordered[,l] > 0)))
}

richness_families_column

## (iii) Add patients: UC & CD
View(matrix_family)
View(Mastertable_viral_unrarefied_euk_virus)

## (iv) Add host group
Mastertable_global_rarefied_fam_complex_heatm <- Mastertable_viral_rarefied_euk_virus_host
used_for_heatmap_later <- aggregate(. ~Final_family, FUN = sum, data = Mastertable_global_rarefied_fam_complex_heatm[,colnames(Mastertable_global_rarefied_fam_complex_heatm) %in% vector_euk | colnames(Mastertable_global_rarefied_fam_complex_heatm) == 'Final_family'])
Mastertable_global_rarefied_fam_complex_heatm <- aggregate(. ~Final_family, FUN = sum, data = Mastertable_global_rarefied_fam_complex_heatm[,colnames(Mastertable_global_rarefied_fam_complex_heatm) %in% vector_euk | colnames(Mastertable_global_rarefied_fam_complex_heatm) == 'Final_family'])
rownames(Mastertable_global_rarefied_fam_complex_heatm) <- Mastertable_global_rarefied_fam_complex_heatm$Final_family
Mastertable_global_rarefied_fam_complex_heatm$Final_family <- NULL
rownames(used_for_heatmap_later) <- used_for_heatmap_later$Final_family
used_for_heatmap_later$Final_family <- NULL

##  Make vector of host by subsetting instead of writing everything!
## Look for an easy way
## Also add disease annotation on top instead of richness/decide which one to go I say numbers
View(used_for_heatmap_later)
Mastertable_global_rarefied_fam_complex_heatm

## (iii) complex family heatmap: WITHOUT ROW DENDROGRAM (ALPHABETICALLY)
col_fun = colorRamp2(c(0,0.001,0.3,0.5,0.7,1), c("gray98","#ffdb00", "#ffa904","#ee7b06","#a12424","#400b0b"))
group_annotations <- c("animal virus","plant and fungal virus", "small circular virus", "plant and fungal virus", "plant and fungal virus","small circular virus","plant and fungal virus", "small circular virus", "small circular virus","plant and fungal virus","animal virus","animal virus","plant and fungal virus","plant and fungal virus","plant and fungal virus","unclassified virus","plant and fungal virus")
col = c("plant and fungal virus" = "lightgreen", "animal virus" = "orange", "small circular virus" = "#74a9cf", "unclassified virus" = "red")

View(matrix_family)
Heatmap_family_without_row_dend_ord <- Heatmap(matrix_family, name = "RA", col = col_fun,  heatmap_width = unit(20, "cm"), heatmap_height = unit(10, "cm"),
                                               column_title = "", column_title_gp = gpar(fontsize = 20, fontface = "bold"),
                                               cluster_rows = FALSE, cluster_columns = TRUE, clustering_distance_columns = "euclidean", column_dend_height = unit(1, "cm"),
                                               border = TRUE,
                                               row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 7),  column_names_rot = 45,
                                               top_annotation = columnAnnotation(richness = richness_families_column, gp = gpar(col = "black"), bar1 = anno_points(richness_families_column)),
                                               right_annotation = rowAnnotation(presence = percentage_families_row, groups = anno_simple(group_annotations, col = col, gp = gpar(col = "black")), bar2 = anno_points(percentage_families_row), gp = gpar(col = "black")))

Heatmap_family_without_row_dend_ord

View(matrix_family)
matrix_family

## 123 samples with at least one eukaryotic virus present (32.6% of all samples)
## The eukaryotic virome in these samples is dominated by a plant virus called Virgaviridae (~ 50%)

## Fix color for group annotations
## Put in illustrator on this figure in gray bar on right top the number of samples with at least one eukar. virus present

## Put later also in barplot the presence and absence of the three groups & total in (i) UC, CD, remission, non-remission, before & after treatment
## think about this?..


# Finish in illustrators
#       --> Increase annotation size: make sure round circiles fall in borders.
#       --> remove titles next to bar & re-order legend next to figures
#       --> Make sure everything is in Arial + Family level in Italic


View(matrix_family)

table(Mastertable_viral_unrarefied_euk_virus$Final_family)

## (iv) division of heatmap into group annotations
## decide row order first - group annotation based on unchanged 
rownames(matrix_family)
[1] "Alphaflexiviridae"  "Anelloviridae"      "Betaflexiviridae"  
[4] "Caliciviridae"      "Chrysoviridae"      "Circoviridae"      
[7] "Closteroviridae"    "Endornaviridae"     "Geminiviridae"     
[10] "Genomoviridae"      "Papillomaviridae"   "Parvoviridae"      
[13] "Picornaviridae"     "Secoviridae"        "Smacoviridae"      
[16] "unclassified CRESS" "Virgaviridae"      

## check col, group annotation, ord fam etc
col_fun = colorRamp2(c(0,0.001,0.3,0.5,0.7,1), c("gray98","#ffdb00", "#ffa904","#ee7b06","#a12424","#400b0b"))
order_fam <- c("Alphaflexiviridae","Betaflexiviridae","Chrysoviridae","Closteroviridae","Endornaviridae","Secoviridae","Virgaviridae","Caliciviridae","Papillomaviridae","Parvoviridae","Picornaviridae","Anelloviridae","Circoviridae","Genomoviridae","Geminiviridae","Smacoviridaee","unclassified CRESS")
group_annotations <- c("plant and fungal virus","small circular virus", "plant and fungal virus","animal virus","plant and fungal virus", "small circular virus","plant and fungal virus","plant and fungal virus","small circular virus","small circular virus","animal virus","animal virus","animal virus","plant and fungal virus", "small circular virus", "small circular virus","plant and fungal virus")
col_ = c("plant and fungal virus" = "lightgreen", "animal virus" = "orange", "small circular virus" = "#74a9cf")

#View(matrix_family)
Heatmap(matrix_family, name = "RA", col = col_fun,  heatmap_width = unit(20, "cm"), heatmap_height = unit(10, "cm"),
        column_title = "", column_title_gp = gpar(fontsize = 20, fontface = "bold"),
        cluster_rows = FALSE, cluster_columns = TRUE, clustering_distance_columns = "euclidean", column_dend_height = unit(1, "cm"),
        border = TRUE,
        row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 7),  column_names_rot = 45,
        top_annotation = columnAnnotation(richness = richness_families_column, gp = gpar(col = "black"), bar1 = anno_points(richness_families_column)),
        right_annotation = rowAnnotation(presence = percentage_families_row, groups = anno_simple(group_annotations, col = col_, gp = gpar(col = "black")), bar2 = anno_points(percentage_families_row), gp = gpar(col = "black")),
        row_split = group_annotations, row_title = NULL,row_gap = unit(1, "mm"))

ncol(matrix_family) # 124 samples have at least one euk. virus present

# same in illustrators
#       --> Increase annotation size: make sure round circiles fall in borders.
#       --> remove titles next to bar & re-order legend next to figures
#       --> Make sure everything is in Arial + Family level in Italic
#       ---> remove sample names

####################################
# 5.3.3 genera
####################################
Mastertable_global_rarefied_gen_complex_heatm_ <- Mastertable_viral_rarefied_euk_virus_host
used_for_heatmap_later <- aggregate(. ~Final_genus, FUN = sum, data = Mastertable_global_rarefied_gen_complex_heatm_[,colnames(Mastertable_global_rarefied_gen_complex_heatm_) %in% vector_euk | colnames(Mastertable_global_rarefied_gen_complex_heatm_) == 'Final_genus'])
Mastertable_global_rarefied_gen_complex_heatm_ <- aggregate(. ~Final_genus, FUN = sum, data = Mastertable_global_rarefied_gen_complex_heatm_[,colnames(Mastertable_global_rarefied_gen_complex_heatm_) %in% vector_euk | colnames(Mastertable_global_rarefied_gen_complex_heatm_) == 'Final_genus'])
rownames(Mastertable_global_rarefied_gen_complex_heatm_) <- Mastertable_global_rarefied_gen_complex_heatm_$Final_family
Mastertable_global_rarefied_gen_complex_heatm_$Final_genus <- NULL
rownames(used_for_heatmap_later) <- used_for_heatmap_later$Final_genus
used_for_heatmap_later$Final_genus <- NULL
genus_names_row <- rownames(used_for_heatmap_later)

percentage_genus_row <- vector()  
for (l in 1: length(genus_names_row)){
  percentage_genus_row[l] <- ((length(which(Mastertable_global_rarefied_gen_complex_heatm_[l,] > 0)))/ncol(used_for_heatmap_later))*100
}

View(Mastertable_global_rarefied_gen_complex_heatm_)
percentage_genus_row

## (ii) Create a percentage of the number of different eukaryotic viral families over all families found within one sample
genus_names_column <- colnames(matrix_genus)
Mastertable_global_rarefied_gen_complex_heatm_dendo_reordered <- Mastertable_global_rarefied_gen_complex_heatm_

richness_genera_column <- vector()  
for (l in 1: length(genus_names_column)){
  richness_genera_column[l] <- (length(which(Mastertable_global_rarefied_gen_complex_heatm_dendo_reordered[,l] > 0)))
}

## (iii) complex genera heatmap: WITHOUT ROW DENDROGRAM (ALPHABETICALLY)
col_fun = colorRamp2(c(0,0.0001,0.3,0.5,0.7,1), c("gray98","#ffdb00", "#ffa904","#ee7b06","#a12424","#400b0b"))
order_gen <- c("Alphapapillomavirus","Cosavirus","Enterovirus","Norovirus","Sapovirus","unclassified Astroviridae","unclassified Parvoviridae","Alphachrysovirus","Alphaendornavirus","Carlavirus","Crinivirus","Nepovirus","Potexvirus","Tobamovirus","unclassified Tombusviridae","Trichovirus","Alphatorquevirus","Begomovirus","Betatorquevirus","Circovirus","Gammatorquevirus","Gemycircularvirus","Gemykrogvirus","Gyrovirus","Huchismacovirus","unclassified Anelloviridae","unclassified Circoviridae","unclassified CRESS")
group_annotations <- c("plant and fungal virus","plant and fungal virus","animal virus","small circular virus","small circular virus","small circular virus","plant and fungal virus","small circular virus","animal virus","plant and fungal virus","animal virus","small circular virus","small circular virus","small circular virus","small circular virus","small circular virus","plant and fungal virus","animal virus","plant and fungal virus","animal virus","plant and fungal virus","plant and fungal virus","small circular virus","animal virus","small circular virus","small circular virus","animal virus","plant and fungal virus")
col = c("plant and fungal virus" = "lightgreen", "animal virus" = "orange", "small circular virus" = "#74a9cf")
Heatmap(matrix_genus, name = "RA", col = col_fun,  heatmap_width = unit(20, "cm"), heatmap_height = unit(12, "cm"),
        column_title = "", column_title_gp = gpar(fontsize = 20, fontface = "bold"),
        cluster_rows = FALSE, cluster_columns = TRUE, clustering_distance_columns = "euclidean", column_dend_height = unit(1, "cm"),
        border = TRUE,
        row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 7),  column_names_rot = 45,
        top_annotation = columnAnnotation(richness = richness_genera_column, gp = gpar(col = "black"), bar1 = anno_points(richness_genera_column)),
        right_annotation = rowAnnotation(presence = percentage_genus_row, groups = anno_simple(group_annotations, col = col, gp = gpar(col = "black")), bar2 = anno_points(percentage_genus_row), gp = gpar(col = "black")))

# Finish in illustrators
#       --> Increase annotation size: make sure round circiles fall in borders.
#       --> remove titles next to bar & re-order legend next to figures
#       --> Make sure everything is in Arial + Family level in Italic

rownames(matrix_genus)

## (iv) division of heatmap into group annotations
col_fun = colorRamp2(c(0,0.001,0.3,0.5,0.7,1), c("gray98","#ffdb00", "#ffa904","#ee7b06","#a12424","#400b0b"))
order_gen <- c("Alphapapillomavirus","Cosavirus","Enterovirus","Norovirus","Sapovirus","unclassified Parvoviridae","Alphachrysovirus","Alphaendornavirus","Carlavirus","Crinivirus","Nepovirus","Potexvirus","Tobamovirus","Trichovirus","Alphatorquevirus","Begomovirus","Betatorquevirus","Circovirus","Gammatorquevirus","Gemycircularvirus","Gemykrogvirus","Gyrovirus","Huchismacovirus","unclassified Circoviridae","unclassified CRESS", "unclassified Genomoviridae")
group_annotations <- c("plant and fungal virus","plant and fungal virus","animal virus","small circular virus","small circular virus","small circular virus","plant and fungal virus","small circular virus","animal virus","plant and fungal virus","animal virus","small circular virus","small circular virus","small circular virus","small circular virus","small circular virus","plant and fungal virus","animal virus","plant and fungal virus","animal virus","plant and fungal virus","plant and fungal virus","small circular virus","small circular virus","small circular virus","animal virus")
col = c("plant and fungal virus" = "lightgreen", "animal virus" = "orange", "small circular virus" = "#74a9cf", "unclassified virus" = "red")
group_annotations
order_gen

Heatmap(matrix_genus, name = "RA", col = col_fun,  heatmap_width = unit(20, "cm"), heatmap_height = unit(12, "cm"),
        column_title = "", column_title_gp = gpar(fontsize = 20, fontface = "bold"), 
        cluster_rows = FALSE, cluster_columns = TRUE, clustering_distance_columns = "euclidean", column_dend_height = unit(1, "cm"),
        border = TRUE,
        row_order = order_gen,
        row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 7),  column_names_rot = 45,
        top_annotation = columnAnnotation(richness = richness_genera_column, gp = gpar(col = "black"), bar1 = anno_points(richness_genera_column)),
        right_annotation = rowAnnotation(presence = percentage_genus_row, groups = anno_simple(group_annotations, col = col, gp = gpar(col = "black")), bar2 = anno_points(percentage_genus_row), gp = gpar(col = "black")),
        row_split = group_annotations, row_title = NULL,row_gap = unit(1, "mm"))
####################################
# 5.3.4 species
####################################
Mastertable_global_rarefied_gen_complex_heatm_ <- Mastertable_viral_rarefied_euk_virus_host
used_for_heatmap_later <- aggregate(. ~Final_species, FUN = sum, data = Mastertable_global_rarefied_gen_complex_heatm_[,colnames(Mastertable_global_rarefied_gen_complex_heatm_) %in% vector_euk | colnames(Mastertable_global_rarefied_gen_complex_heatm_) == 'Final_species'])
Mastertable_global_rarefied_gen_complex_heatm_ <- aggregate(. ~Final_species, FUN = sum, data = Mastertable_global_rarefied_gen_complex_heatm_[,colnames(Mastertable_global_rarefied_gen_complex_heatm_) %in% vector_euk | colnames(Mastertable_global_rarefied_gen_complex_heatm_) == 'Final_species'])
rownames(Mastertable_global_rarefied_gen_complex_heatm_) <- Mastertable_global_rarefied_gen_complex_heatm_$Final_species
Mastertable_global_rarefied_gen_complex_heatm_$Final_species <- NULL
rownames(used_for_heatmap_later) <- used_for_heatmap_later$Final_species
used_for_heatmap_later$Final_species <- NULL
genus_names_row <- rownames(used_for_heatmap_later)

percentage_genus_row <- vector()  
for (l in 1: length(genus_names_row)){
  percentage_genus_row[l] <- ((length(which(Mastertable_global_rarefied_gen_complex_heatm_[l,] > 0)))/ncol(used_for_heatmap_later))*100
}

View(Mastertable_viral_rarefied_euk_virus_host)
percentage_genus_row
View(percentage_genus_row)

rownames(Mastertable_global_rarefied_gen_complex_heatm_)
percentage_genus_row
## (ii) Create a percentage of the number of different eukaryotic viral families over all families found within one sample
genus_names_column <- colnames(matrix_genus)
Mastertable_global_rarefied_gen_complex_heatm_dendo_reordered <- Mastertable_global_rarefied_gen_complex_heatm_

richness_genera_column <- vector()  
for (l in 1: length(genus_names_column)){
  richness_genera_column[l] <- (length(which(Mastertable_global_rarefied_gen_complex_heatm_dendo_reordered[,l] > 0)))
}
####################################
# 6 Barplot: presence/absence euk. viruses
####################################
library("readxl")

setwd("/Users/daan/Desktop/Bioinformatics/pre-processing/Biological_study/Metadata/output_metadata")
getwd()
dir()

metadata <- read_excel("Final_Metadata_subselection_All.xlsx")
rownames(metadata) <- metadata$`Randomized Nr.`
vector_diagnosis <- (which(names(metadata) == "Diagnosis"))
vector_B <- (which(names(metadata) == "Randomized Nr."))
vector_R <- (which(names(metadata) == "Endoscopic_outcome_combined (NR/R)"))
vector_W <- (which(names(metadata) == "Week"))

metadata_diagnosis <- metadata[c(vector_diagnosis,vector_B, vector_R,vector_W)]
rownames(metadata_diagnosis) <- metadata_diagnosis$`Randomized Nr.`

Mastertable_global_rarefied_fam_complex_heatm <- Mastertable_viral_rarefied_euk_virus_host

Mastertable_global_rarefied_fam_complex_heatm <- aggregate(. ~Final_family, FUN = sum, data = Mastertable_global_rarefied_fam_complex_heatm[,colnames(Mastertable_global_rarefied_fam_complex_heatm) %in% vector_euk | colnames(Mastertable_global_rarefied_fam_complex_heatm) == 'Final_family'])
rownames(Mastertable_global_rarefied_fam_complex_heatm) <- Mastertable_global_rarefied_fam_complex_heatm$Final_family
Mastertable_global_rarefied_fam_complex_heatm$Final_family <- NULL

Mastertable_viral_unrarefied_euk_virus_host_t <- t(Mastertable_global_rarefied_fam_complex_heatm)
merged_euk_table <- merge(Mastertable_viral_unrarefied_euk_virus_host_t, metadata_diagnosis, by=0, all=F) 
rownames(merged_euk_table) <- merged_euk_table$`Randomized Nr.`
merged_euk_table$`Randomized Nr.` <- NULL
merged_euk_table$Row.names <- NULL

## Add host
merged_euk_table$`plant and fungal viruses` <- "0"
merged_euk_table$`small circular viruses` <- "0"
merged_euk_table$`animal viruses` <- "0"
str(merged_euk_table)
colnames(merged_euk_table)

merged_euk_table$`plant and fungal viruses`[merged_euk_table$Alphaflexiviridae > 0 | merged_euk_table$Betaflexiviridae > 0 | merged_euk_table$Chrysoviridae > 0 | merged_euk_table$Closteroviridae > 0 | merged_euk_table$Endornaviridae > 0 | merged_euk_table$Secoviridae > 0 | merged_euk_table$Virgaviridae > 0] <- 1
merged_euk_table$`animal viruses`[merged_euk_table$Caliciviridae > 0 | merged_euk_table$Papillomaviridae > 0 | merged_euk_table$Parvoviridae > 0 | merged_euk_table$Picornaviridae > 0] <- 1
merged_euk_table$`small circular viruses`[merged_euk_table$Anelloviridae > 0 | merged_euk_table$Circoviridae > 0 | merged_euk_table$`unclassified CRESS` > 0 | merged_euk_table$Geminiviridae > 0 | merged_euk_table$Genomoviridae > 0 | merged_euk_table$Smacoviridae > 0] <- 1

merged_euk_table$`Endoscopic_outcome_combined (NR/R)`[merged_euk_table$`Endoscopic_outcome_combined (NR/R)` == "1"] <- "remission"
merged_euk_table$`Endoscopic_outcome_combined (NR/R)`[merged_euk_table$`Endoscopic_outcome_combined (NR/R)` == "0"] <- "non-remission"
merged_euk_table$`Endoscopic_outcome_combined (NR/R)`[merged_euk_table$`Endoscopic_outcome_combined (NR/R)` == "NA"] <- "unknown"

merged_euk_table$`plant and fungal viruses`[merged_euk_table$`plant and fungal viruses` == "1"] <- "present"
merged_euk_table$`plant and fungal viruses`[merged_euk_table$`plant and fungal viruses` == "0"] <- "absent"

merged_euk_table$`small circular viruses`[merged_euk_table$`small circular viruses` == "1"] <- "present"
merged_euk_table$`small circular viruses`[merged_euk_table$`small circular viruses` == "0"] <- "absent"

merged_euk_table$`animal viruses`[merged_euk_table$`animal viruses` == "1"] <- "present"
merged_euk_table$`animal viruses`[merged_euk_table$`animal viruses` == "0"] <- "absent"
View(merged_euk_table)

esquisser(merged_euk_table)
(nrow(merged_euk_table[merged_euk_table$`plant and fungal viruses` == "present",])/nrow(merged_euk_table))*100
(nrow(merged_euk_table[merged_euk_table$`plant and fungal viruses` == "present",])/377)*100
nrow(merged_euk_table[merged_euk_table$`plant and fungal viruses` == "present",])

## (A) PLANT VIRUSES
# UC: plant
merged_euk_table_UC_plant_pres <- merged_euk_table %>%
 filter(Diagnosis %in% "UC") %>%
 filter(!(`Endoscopic_outcome_combined (NR/R)` %in% 
 "unknown")) %>%
 filter(`plant and fungal viruses` %in% "present") %>%
 ggplot() +
 aes(x = `plant and fungal viruses`, fill = `Endoscopic_outcome_combined (NR/R)`) +
 geom_bar(position = "dodge") +
 scale_fill_hue(direction = 1) +
 theme_bw() +
 theme(axis.text.x = element_blank(),axis.ticks.x = element_blank()) +
 facet_wrap(vars(Week),scales = "fixed")

merged_euk_table_UC_plant_pres

# CD :plant
merged_euk_table_CD_plant_pres <- merged_euk_table %>%
  filter(Diagnosis %in% "CD") %>%
  filter(!(Week %in% "w14")) %>%
  filter(!(`Endoscopic_outcome_combined (NR/R)` %in% 
             "unknown")) %>%
  filter(`plant and fungal viruses` %in% "present") %>%
  ggplot() +
  aes(x = `plant and fungal viruses`, fill = `Endoscopic_outcome_combined (NR/R)`) +
  geom_bar(position = "dodge") +
  scale_fill_hue(direction = 1) +
  theme_bw() +
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank()) +
  facet_wrap(vars(Week), scales = "fixed")

merged_euk_table_CD_plant_pres

## Merge both now in Illustrator: this is for unrarefied data
Mastertable_viral_unrarefied_euk_virus_tt <-t(Mastertable_viral_rarefied_euk_virus_host)
merged_euk_table_t <- merge(Mastertable_viral_unrarefied_euk_virus_tt, metadata_diagnosis, by=0, all=F) 
nrow(merged_euk_table_t) # 377 IBD samples

## UC
nrow(merged_euk_table_t[merged_euk_table_t$Diagnosis == "UC",]) # 95 UC samples
nrow(merged_euk_table_t[merged_euk_table_t$Diagnosis == "UC" & merged_euk_table_t$Week == "w0" & merged_euk_table_t$`Endoscopic_outcome_combined (NR/R)` == "1",]) # 22 w0
nrow(merged_euk_table_t[merged_euk_table_t$Diagnosis == "UC" & merged_euk_table_t$Week == "w0" & merged_euk_table_t$`Endoscopic_outcome_combined (NR/R)` == "0",]) # 22 w0
nrow(merged_euk_table[merged_euk_table$Diagnosis == "UC" & merged_euk_table$Week == "w0" & merged_euk_table$`Endoscopic_outcome_combined (NR/R)` == "remission" & merged_euk_table$`plant and fungal viruses` == "present",]) # 3 w0
nrow(merged_euk_table[merged_euk_table$Diagnosis == "UC" & merged_euk_table$Week == "w0" & merged_euk_table$`Endoscopic_outcome_combined (NR/R)` == "non-remission" & merged_euk_table$`plant and fungal viruses` == "present",]) # 7 w0
(3/22)*100 # 13.6%
(7/22)*100 # 31.8%

nrow(merged_euk_table_t[merged_euk_table_t$Diagnosis == "UC" & merged_euk_table_t$Week == "w14" & merged_euk_table_t$`Endoscopic_outcome_combined (NR/R)` == "1",]) # 31 w14
nrow(merged_euk_table_t[merged_euk_table_t$Diagnosis == "UC" & merged_euk_table_t$Week == "w14" & merged_euk_table_t$`Endoscopic_outcome_combined (NR/R)` == "0",]) # 20 w14
nrow(merged_euk_table[merged_euk_table$Diagnosis == "UC" & merged_euk_table$Week == "w14" & merged_euk_table$`Endoscopic_outcome_combined (NR/R)` == "remission" & merged_euk_table$`plant and fungal viruses` == "present",]) # 10 w0
nrow(merged_euk_table[merged_euk_table$Diagnosis == "UC" & merged_euk_table$Week == "w14" & merged_euk_table$`Endoscopic_outcome_combined (NR/R)` == "non-remission" & merged_euk_table$`plant and fungal viruses` == "present",]) # 7 w0
(10/31)*100 # 32.3%
(7/20)*100 # 35%

## 95 UC patients; 44 w0; 3/22 remission; Of these, 13.6% contain plant or fungal viruses
## 95 UC patients; 44 w0; 7/22 non-remission; Of these, 31.8% contain plant or fungal viruses
## 95 UC patients; 51 w14; 10/31 remission; Of these 32.3% contain plant or fungal viruses
## 95 UC patients; 51 w14; 7/20 non-remission; Of these 35% contain plant or fungal viruses

View(merged_euk_table_t)
## CD
nrow(merged_euk_table_t[merged_euk_table_t$Diagnosis == "CD",]) # 282 UC
nrow(merged_euk_table_t[merged_euk_table_t$Diagnosis == "CD" & merged_euk_table_t$Week == "w0" & merged_euk_table_t$`Endoscopic_outcome_combined (NR/R)` == "1",]) # 43 w0
nrow(merged_euk_table_t[merged_euk_table_t$Diagnosis == "CD" & merged_euk_table_t$Week == "w0" & merged_euk_table_t$`Endoscopic_outcome_combined (NR/R)` == "0",]) # 47 w0
nrow(merged_euk_table_t[merged_euk_table_t$Diagnosis == "CD" & merged_euk_table_t$Week == "w0" & merged_euk_table_t$`Endoscopic_outcome_combined (NR/R)` == "NA",]) # 23 w0

nrow(merged_euk_table[merged_euk_table$Diagnosis == "CD" & merged_euk_table$Week == "w0" & merged_euk_table$`Endoscopic_outcome_combined (NR/R)` == "remission" & merged_euk_table$`plant and fungal viruses` == "present",]) # 7 w0
nrow(merged_euk_table[merged_euk_table$Diagnosis == "CD" & merged_euk_table$Week == "w0" & merged_euk_table$`Endoscopic_outcome_combined (NR/R)` == "non-remission" & merged_euk_table$`plant and fungal viruses` == "present",]) # 10 w0
nrow(merged_euk_table[merged_euk_table$Diagnosis == "CD" & merged_euk_table$Week == "w0" & merged_euk_table$`Endoscopic_outcome_combined (NR/R)` == "unknown" & merged_euk_table$`plant and fungal viruses` == "present",]) # 5 w0
(7/43)*100 # 16.3%
(10/47)*100 # 21.3%
(5/23)*100 # 21.7

nrow(merged_euk_table_t[merged_euk_table_t$Diagnosis == "CD" & merged_euk_table_t$Week == "w24" & merged_euk_table_t$`Endoscopic_outcome_combined (NR/R)` == "1",]) # 46 w24
nrow(merged_euk_table_t[merged_euk_table_t$Diagnosis == "CD" & merged_euk_table_t$Week == "w24" & merged_euk_table_t$`Endoscopic_outcome_combined (NR/R)` == "0",]) # 48 w24
nrow(merged_euk_table_t[merged_euk_table_t$Diagnosis == "CD" & merged_euk_table_t$Week == "w24" & merged_euk_table_t$`Endoscopic_outcome_combined (NR/R)` == "NA",]) # 28 w24

nrow(merged_euk_table[merged_euk_table$Diagnosis == "CD" & merged_euk_table$Week == "w24" & merged_euk_table$`Endoscopic_outcome_combined (NR/R)` == "remission" & merged_euk_table$`plant and fungal viruses` == "present",]) # 11 w24
nrow(merged_euk_table[merged_euk_table$Diagnosis == "CD" & merged_euk_table$Week == "w24" & merged_euk_table$`Endoscopic_outcome_combined (NR/R)` == "non-remission" & merged_euk_table$`plant and fungal viruses` == "present",]) # 8 w24
nrow(merged_euk_table[merged_euk_table$Diagnosis == "CD" & merged_euk_table$Week == "w24" & merged_euk_table$`Endoscopic_outcome_combined (NR/R)` == "unknown" & merged_euk_table$`plant and fungal viruses` == "present",]) # 8/28 w24

(11/46)*100 # 23.9%
(8/48)*100 # 16.7%
(8/28)*100 # 28.6

View(merged_euk_table_t)
nrow(merged_euk_table_t[merged_euk_table_t$Diagnosis == "CD" & merged_euk_table_t$`Endoscopic_outcome_combined (NR/R)` == "NA",]) # 46 w24
nrow(merged_euk_table_t[merged_euk_table_t$Diagnosis == "CD"  & merged_euk_table_t$`Endoscopic_outcome_combined (NR/R)` == "NA",]) # 46 w24

## 282 CD
## 282 CD patients; 90 w0; 7/43 remission; Of these, 16.3% contain plant or fungal viruses
## 282 CD patients; 90 w0; 10/47 non-remission; Of these 21.3% contain plant or fungal viruses
## 282 CD patients; 94 w24; 11/46 remission; Of these 23.9% contain plant or fungal viruses
## 282 CD patients; 94 w24; 8/48 non-remission; Of these 16.7% contain plant or fungal viruses

Percent_euk <- as.data.frame(c("plant and fungal viruses","plant and fungal viruses","plant and fungal viruses","plant and fungal viruses","plant and fungal viruses","plant and fungal viruses","plant and fungal viruses","plant and fungal viruses"))
Percent_euk$host <- Percent_euk$`c("plant and fungal viruses", "plant and fungal viruses", "plant and fungal viruses", "plant and fungal viruses", "plant and fungal viruses", "plant and fungal viruses", "plant and fungal viruses", "plant and fungal viruses")`
Percent_euk$`c("plant and fungal viruses", "plant and fungal viruses", "plant and fungal viruses", "plant and fungal viruses", "plant and fungal viruses", "plant and fungal viruses", "plant and fungal viruses", "plant and fungal viruses")` <- NULL
Percent_euk$diagnosis <- c("UC", "UC", "UC", "UC", "CD", "CD", "CD","CD")
Percent_euk$remission <- c("remission","non-remission","remission","non-remission","remission","non-remission","remission","non-remission")
Percent_euk$Week <- c("w0","w0","w14","w14","w0","w0","w24","w24")
Percent_euk$Percentage <- as.numeric(c("10.0","26.9","28.9","26.9","15.5","20.0","23.9","15.1"))

View(Percent_euk)
esquisser(Percent_euk)
  
library(ggplot2)

ggplot(Percent_euk) +
 aes(x = Week, fill = remission, weight = Percentage) +
 geom_bar(position = "dodge") +
 scale_fill_hue(direction = 1) +
 theme_bw() +
 facet_wrap(vars(diagnosis))
####################################
# 7 Create contingency table for presence
####################################
## 95 UC patients; 44 w0; 3/22 remission; Of these, 13.6% contain plant or fungal viruses
## 95 UC patients; 44 w0; 7/22 non-remission; Of these, 31.8% contain plant or fungal viruses
## 95 UC patients; 51 w14; 10/31 remission; Of these 32.3% contain plant or fungal viruses
## 95 UC patients; 51 w14; 7/20 non-remission; Of these 35% contain plant or fungal viruses

## 282 CD
## 282 CD patients; 90 w0; 7/43 remission; Of these, 16.3% contain plant or fungal viruses
## 282 CD patients; 90 w0; 10/47 non-remission; Of these 21.3% contain plant or fungal viruses
## 282 CD patients; 94 w24; 11/46 remission; Of these 23.9% contain plant or fungal viruses
## 282 CD patients; 94 w24; 8/48 non-remission; Of these 16.7% contain plant or fungal viruses

merged_euk_table[is.na(merged_euk_table$`Endoscopic_outcome_combined (NR/R)`)] <- "unknown"
nrow(merged_euk_table[!merged_euk_table$`Endoscopic_outcome_combined (NR/R)` == "unknown", ]) # 98 are unknown
nrow(merged_euk_table[merged_euk_table$`Endoscopic_outcome_combined (NR/R)` == "remission" & merged_euk_table$Diagnosis == "UC" & !merged_euk_table$Week == "w14" & merged_euk_table$`plant and fungal viruses` == "present", ]) # 13
nrow(merged_euk_table[merged_euk_table$`Endoscopic_outcome_combined (NR/R)` == "non-remission" & merged_euk_table$Diagnosis == "UC"  & merged_euk_table$`plant and fungal viruses` == "present", ]) # 14
nrow(merged_euk_table[merged_euk_table$`Endoscopic_outcome_combined (NR/R)` == "remission" & merged_euk_table$Diagnosis == "CD" & !merged_euk_table$Week == "w14"  & merged_euk_table$`plant and fungal viruses` == "present", ]) # 18
nrow(merged_euk_table[merged_euk_table$`Endoscopic_outcome_combined (NR/R)` == "non-remission" & merged_euk_table$Diagnosis == "CD" & !merged_euk_table$Week == "w14"  & merged_euk_table$`plant and fungal viruses` == "present",]) # 18

# total: 63 IBD samples with presence of euk. virus in one of these 4 groups. 
# These groups have a total of 279 samples of UC, CD (not w14), in remission or non-remission?
# Therefore, there are still 98 samples that have no endoscopical outcome result, or are at w14 for CD.

# A) Create DF with pres/abs & groups
Contiguency_table <- as.data.frame(c("UC-w0-R", "UC-w0-NR", "UC-w14-R", "UC-w14-NR", "CD-w0-R", "CD-w0-NR", "CD-w24-R", "CD-w24-NR"))
Contiguency_table$groups <- Contiguency_table$`c("UC-w0-R", "UC-w0-NR", "UC-w14-R", "UC-w14-NR", "CD-w0-R", "CD-w0-NR", "CD-w24-R", "CD-w24-NR")`
Contiguency_table$`c("UC-w0-R", "UC-w0-NR", "UC-w14-R", "UC-w14-NR", "CD-w0-R", "CD-w0-NR", "CD-w24-R", "CD-w24-NR")` <- NULL
Contiguency_table$plant_presence <- c(3,7,10,7,7,10,11,8)
Contiguency_table$plant_absence <- c(19,15,21,13,36,37,35,40)
rownames(Contiguency_table) <- Contiguency_table$groups
Contiguency_table$groups <- NULL
View(Contiguency_table)
# 27 UC
# 36 CD
# 63 IBD

# B) Graphical display of contengency tables
# http://www.sthda.com/english/wiki/chi-square-test-of-independence-in-r
library("gplots")

# 1. convert the data as a table
dt <- as.table(as.matrix(Contiguency_table))

# 2. Graph
balloonplot(t(dt), main ="plant viral presence", xlab ="", ylab="",
            label = TRUE, show.margins = FALSE)

# C) mosaiccplot
install.packages("graphics")
library("graphics")
mosaicplot(dt, shade = TRUE, las=2,main = "")

# D) Chi-squared test
chisq <- chisq.test(Contiguency_table)
chisq
chisq$observed
round(chisq$expected,2)
round(chisq$residuals, 3)

# E) Visualize residuals
library(corrplot)
corrplot(chisq$residuals, is.cor = FALSE)

contrib <- 100*chisq$residuals^2/chisq$statistic
round(contrib, 3)
corrplot(contrib, is.cor = FALSE)

# F) P-value 
chisq$p.value
# printing the mean
chisq$estimate

## Subset groups: chi-squared
chisq <- chisq.test(Contiguency_table)
chisq

?chisq.test()

# UC
chisq <- chisq.test(Contiguency_table[1:4,])
chisq

# Baseline UC
chisq <- chisq.test(Contiguency_table[1:2,])
chisq

# w14 UC
chisq <- chisq.test(Contiguency_table[c(1,3),])
chisq
####################################
## PS: Great package for connections of contigs with Circlize in https://jokergoo.github.io/circlize_book/book/advanced-layout.html#arrange-multiple-plots
######### !!!!!!!!!!!!!!! #########



## Now, make a complex heatmap, both of phages and eukaryotic viruses in the same heatmap. or two? think how to proceed

# Monday: Complex heatmap for skin samples - end this & put in ppt on desktop - Discussion Tine Marcelis about pipeline (tell to make PubCrawler + Github account and add to our lab)
# Tuesday: Create template for powerpoint Research Seminar broadly.. (I) Descriptive analysis & new method of viral identification, annotation etc (II) RA,richness, diversity, (III) Host prediction, (IV) Predictive biomarker: is it possible?, (V) additionlaly try of clustering..
# Make Tuesday big lines, use perhaps wednesday thursday

# --> Put in powerpoint & make with illustrator some symbols (infect humans, infect diatoms, infect fungi, etc)
# --> Do ONLY for family level because most of them we cannot annotate on genus level anyway, likely to do with fragmentations.


# SAVE Global Environment #
######### !!!!!!!!!!!!!!! #########