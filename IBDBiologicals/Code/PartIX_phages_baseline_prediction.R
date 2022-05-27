####################################
# SCRIPT 11: PROKARYOTIC VIRUSES: qPCR baseline prediction
####################################
# ------------------> Load the "Input_qPCR_data" 
####################################
# 0. Packages: Install packages and load them wherever needed
####################################
#install.packages("readxl")
#install.packages("pscl")
#install.packages("lmtest")
#install.packages("pROC")
#install.packages("caret")
#install.packages("reshape2")
#install.packages("/Users/daan/Downloads/randomForest_4.6-1.tar.gz", repos = NULL, type = "source")
#install.packages("car")
library(randomForest)
library(readxl)
library(pscl)
library(lmtest)
library(pROC)
library(randomForest)
library(caret)
library(reshape2)
library(esquisse)
library(ggplot2)
library(gridExtra)
library(scales)
library(dplyr)
library(car)
####################################
# 1. Import qPCR data
####################################
# 1.1 Import DF
####################################
setwd("/Users/daan/Desktop/Bioinformatics/pre-processing/Biological_study/wetlab/qPCR")
getwd()
dir()

baseline_DF <- read_excel('Input_R_qPCR_model.xlsx')
baseline_DF <- as.data.frame(baseline_DF)
rownames(baseline_DF) <- baseline_DF$`Sample ID`

## Import NGS data
dir()
baseline_DF_NGS <- read_excel('Input_R_NGS_model.xlsx')
baseline_DF_NGS <- as.data.frame(baseline_DF_NGS)
rownames(baseline_DF_NGS) <- baseline_DF_NGS$`Sample ID`
####################################
# 1.2 Basic manipulations DF
####################################
baseline_DF$viral_copies_B26 <- as.numeric(baseline_DF$viral_copies_B26)
baseline_DF$viral_copies_B380 <- as.numeric(baseline_DF$viral_copies_B380)
baseline_DF$viral_copies_B85 <- as.numeric(baseline_DF$viral_copies_B85)
baseline_DF$viral_copies_B137 <- as.numeric(baseline_DF$viral_copies_B137)
baseline_DF$viral_copies_B261 <- as.numeric(baseline_DF$viral_copies_B261)

baseline_DF_NGS$NGS_B26 <- as.numeric(baseline_DF_NGS$NGS_B26)
baseline_DF_NGS$NGS_B380 <- as.numeric(baseline_DF_NGS$NGS_B380)
baseline_DF_NGS$NGS_B261 <- as.numeric(baseline_DF_NGS$NGS_B261)
baseline_DF_NGS$NGS_B85 <- as.numeric(baseline_DF_NGS$NGS_B85)
baseline_DF_NGS$NGS_B137 <- as.numeric(baseline_DF_NGS$NGS_B137)
####################################
# 2. Calculations data
####################################
baseline_DF_figures <- baseline_DF
baseline_DF_figures[is.na(baseline_DF_figures)] <- 0
# To visualize we convert the NA's to zero

# B26
baseline_DF_figures$viral_copies_B26 <- baseline_DF_figures$viral_copies_B26+1
baseline_DF_figures$`Sample ID`<- reorder(baseline_DF_figures$`Sample ID`, baseline_DF_figures$viral_copies_B26)
ggplot(baseline_DF_figures) +
 aes(x = `Sample ID`, y = viral_copies_B26) +
  geom_point(shape=21, fill="#499B78", colour="black", size=1) +
  theme_bw() +
 theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
 theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
  geom_hline(yintercept=100, linetype="dashed",color = "darkred", size=0.3) +
  scale_y_log10(limit = c(1e+00,1e+08), breaks = c(1e+00, 1e+01, 1e+02, 1e+02, 1e+03, 1e+04, 1e+05, 1e+06,1e+07,1e+08)) 

# B380
baseline_DF_figures$viral_copies_B380 <- baseline_DF_figures$viral_copies_B380+1
baseline_DF_figures$`Sample ID`<- reorder(baseline_DF_figures$`Sample ID`, baseline_DF_figures$viral_copies_B380)
ggplot(baseline_DF_figures) +
  aes(x = `Sample ID`, y = viral_copies_B380) +
  geom_point(shape=21, fill="#499B78", colour="black", size=1) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
  geom_hline(yintercept=100, linetype="dashed",color = "darkred", size=0.3) +
  scale_y_log10(limit = c(1e+00,1e+08), breaks = c(1e+00, 1e+01, 1e+02, 1e+02, 1e+03, 1e+04, 1e+05, 1e+06,1e+07,1e+08)) 

# B85
baseline_DF_figures$viral_copies_B85 <- baseline_DF_figures$viral_copies_B85+1
baseline_DF_figures$`Sample ID`<- reorder(baseline_DF_figures$`Sample ID`, baseline_DF_figures$viral_copies_B85)
ggplot(baseline_DF_figures) +
  aes(x = `Sample ID`, y = viral_copies_B85) +
  geom_point(shape=21, fill="#AD4B38", colour="black", size=1) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
  geom_hline(yintercept=100, linetype="dashed",color = "darkred", size=0.3) +
  scale_y_log10(limit = c(1e+00,1e+08), breaks = c(1e+00, 1e+01, 1e+02, 1e+02, 1e+03, 1e+04, 1e+05, 1e+06,1e+07,1e+08)) 

# B137
baseline_DF_figures$viral_copies_B137 <- baseline_DF_figures$viral_copies_B137+1
baseline_DF_figures$`Sample ID`<- reorder(baseline_DF_figures$`Sample ID`, baseline_DF_figures$viral_copies_B137)
ggplot(baseline_DF_figures) +
  aes(x = `Sample ID`, y = viral_copies_B137) +
  geom_point(shape=21, fill="#AD4B38", colour="black", size=1) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
  geom_hline(yintercept=100, linetype="dashed",color = "darkred", size=0.3) +
  scale_y_log10(limit = c(1e+00,1e+08), breaks = c(1e+00, 1e+01, 1e+02, 1e+02, 1e+03, 1e+04, 1e+05, 1e+06,1e+07,1e+08)) 

# B261
baseline_DF_figures$viral_copies_B261 <- baseline_DF_figures$viral_copies_B261+1
baseline_DF_figures$`Sample ID`<- reorder(baseline_DF_figures$`Sample ID`, baseline_DF_figures$viral_copies_B261)
ggplot(baseline_DF_figures) +
  aes(x = `Sample ID`, y = viral_copies_B261) +
  geom_point(shape=21, fill="#499B78", colour="black", size=1) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
  geom_hline(yintercept=100, linetype="dashed",color = "darkred", size=0.3) +
  scale_y_log10(limit = c(1e+00,1e+08), breaks = c(1e+00, 1e+01, 1e+02, 1e+02, 1e+03, 1e+04, 1e+05, 1e+06,1e+07,1e+08)) 

View(baseline_DF_figures)

## do later for all & add percentage NGS & qCPR +
####################################
# 2.1 Introduce cutoff
####################################
baseline_DF_500 <- baseline_DF
baseline_DF_500[,c(5:9)][baseline_DF_500[,c(5:9)] < 500] <- 0
table(baseline_DF_500$Diagnosis[baseline_DF_500$Diagnosis == "UC"])

# Percentage B26 presence
nrow(baseline_DF_500[baseline_DF_500$viral_copies_B26 > 0,])/nrow(baseline_DF_500)*100
nrow(baseline_DF_500[baseline_DF_500$Diagnosis == "UC" & baseline_DF_500$viral_copies_B26 > 0,])/nrow(baseline_DF_500[baseline_DF_500$Diagnosis == "UC",])*100 # 18,2%
nrow(baseline_DF_500[baseline_DF_500$Diagnosis == "CD" & baseline_DF_500$viral_copies_B26 > 0,])/nrow(baseline_DF_500[baseline_DF_500$Diagnosis == "CD",])*100 # 9,52%

# Percentage B380 presence
nrow(baseline_DF_500[baseline_DF_500$viral_copies_B380 > 0,])/nrow(baseline_DF_500)*100
nrow(baseline_DF_500[baseline_DF_500$Diagnosis == "UC" & baseline_DF_500$viral_copies_B380 > 0,])/nrow(baseline_DF_500[baseline_DF_500$Diagnosis == "UC",])*100 # 43,2%
nrow(baseline_DF_500[baseline_DF_500$Diagnosis == "CD" & baseline_DF_500$viral_copies_B380 > 0,])/nrow(baseline_DF_500[baseline_DF_500$Diagnosis == "CD",])*100 # 47,6%

# Percentage B85 presence
nrow(baseline_DF_500[baseline_DF_500$viral_copies_B85 > 0,])/nrow(baseline_DF_500)*100
nrow(baseline_DF_500[baseline_DF_500$Diagnosis == "UC" & baseline_DF_500$viral_copies_B85 > 0,])/nrow(baseline_DF_500[baseline_DF_500$Diagnosis == "UC",])*100 # 54,5%
nrow(baseline_DF_500[baseline_DF_500$Diagnosis == "CD" & baseline_DF_500$viral_copies_B85 > 0,])/nrow(baseline_DF_500[baseline_DF_500$Diagnosis == "CD",])*100 # 45,2%

# Percentage B137 presence
nrow(baseline_DF_500[baseline_DF_500$viral_copies_B137 > 0,])/nrow(baseline_DF_500)*100
nrow(baseline_DF_500[baseline_DF_500$Diagnosis == "UC" & baseline_DF_500$viral_copies_B137 > 0,])/nrow(baseline_DF_500[baseline_DF_500$Diagnosis == "UC",])*100 # 34,1%
nrow(baseline_DF_500[baseline_DF_500$Diagnosis == "CD" & baseline_DF_500$viral_copies_B137 > 0,])/nrow(baseline_DF_500[baseline_DF_500$Diagnosis == "CD",])*100 # 33,3%

# Percentage B261 presence
nrow(baseline_DF_500[baseline_DF_500$viral_copies_B261 > 0,])/nrow(baseline_DF_500)*100
nrow(baseline_DF_500[baseline_DF_500$Diagnosis == "UC" & baseline_DF_500$viral_copies_B261 > 0,])/nrow(baseline_DF_500[baseline_DF_500$Diagnosis == "UC",])*100 # 63,5%
nrow(baseline_DF_500[baseline_DF_500$Diagnosis == "CD" & baseline_DF_500$viral_copies_B261 > 0,])/nrow(baseline_DF_500[baseline_DF_500$Diagnosis == "CD",])*100 # 47,6%

#View(baseline_DF_500)
####################################
# 2.2 Proportion and log10
####################################
# 2.2.1 Proportion
####################################
baseline_DF_500_prop <- baseline_DF_500
baseline_DF_500_prop$proportion <- 1

# Proportion is calculated by summing up the absolute quantities of marker suposedly predicting remission divided by marker suposedly predicting non-remission
# Add proportion UC by formula (B26+B261+1)/(B85+1)  
baseline_DF_500_prop$proportion[baseline_DF_500_prop$Diagnosis == "UC"] <- ((baseline_DF_500_prop$viral_copies_B26[baseline_DF_500_prop$Diagnosis == "UC"]+baseline_DF_500_prop$viral_copies_B261[baseline_DF_500_prop$Diagnosis == "UC"]+1)/(baseline_DF_500_prop$viral_copies_B85[baseline_DF_500_prop$Diagnosis == "UC"]+1))

# Add proportion CD by formula (B26+B380+1)/(B137+1)  
baseline_DF_500_prop$proportion[baseline_DF_500_prop$Diagnosis == "CD"] <- ((baseline_DF_500_prop$viral_copies_B26[baseline_DF_500_prop$Diagnosis == "CD"]+baseline_DF_500_prop$viral_copies_B380[baseline_DF_500_prop$Diagnosis == "CD"]+1)/(baseline_DF_500_prop$viral_copies_B137[baseline_DF_500_prop$Diagnosis == "CD"]+1))

# Add proportion IBD by formula (B26+B261+1)/(B85+1)  
#baseline_DF_500_prop$proportion <- ((baseline_DF_500_prop$viral_copies_B26 + baseline_DF_500_prop$viral_copies_B380+1)/(baseline_DF_500_prop$viral_copies_B85+baseline_DF_500_prop$viral_copies_B137+1))
####################################
# 2.2.2 Log10
####################################
baseline_DF_500_prop$log_proportion <- log10(baseline_DF_500_prop$proportion)
####################################
# 2.2.3 Basic evaluation of formula and judgment TP/TN
####################################
# Convert NA's to zero's
baseline_DF_500_prop_ev <- baseline_DF_500_prop
baseline_DF_500_prop_ev$log_proportion[is.na(baseline_DF_500_prop_ev$log_proportion)] <- 0 # 128 samples
baseline_DF_500_prop_ev <- baseline_DF_500_prop_ev[!baseline_DF_500_prop_ev$log_proportion == 0,] # 82 samples
View(baseline_DF_500_prop_ev)

# Positive predictor (for remission) that is actually positive
table_postive <- as.data.frame(table(baseline_DF_500_prop_ev$`Endoscopic outcome`[baseline_DF_500_prop_ev$log_proportion > 0]))
sum(table_postive$Freq) # Positive predictions
(table_postive$Freq[table_postive$Var1 == "Yes"]/sum(table_postive$Freq))*100 # Actually positive
# 66.7% TP

# Negative predictor (for non-remission) that is actually negative
table_negative <- as.data.frame(table(baseline_DF_500_prop_ev$`Endoscopic outcome`[baseline_DF_500_prop_ev$log_proportion < 0]))
sum(table_negative$Freq) # Negative predictions
(table_negative$Freq[table_negative$Var1 == "No"]/sum(table_negative$Freq))*100 # Actually negative
# 88% TN

# No prediction able to make
table_neutral <- as.data.frame(table(baseline_DF_500_prop_ev$`Endoscopic outcome`[baseline_DF_500_prop_ev$log_proportion == 0]))
sum(table_neutral$Freq) # total non-predictions 
# 17 No
# 17 Yes

# total predictions made for remission, non-remission, neutral
sum(table_postive$Freq)/nrow(baseline_DF_500_prop)*100 # + predictions
sum(table_negative$Freq)/nrow(baseline_DF_500_prop)*100 # - predictions
sum(table_neutral$Freq)/nrow(baseline_DF_500_prop)*100 # neutral predictions
####################################
# 3. Logistic regression model to make predictions
####################################
# 3.1 Trim dataframe
####################################
# Remove redundant columns
baseline_DF_500_prop_model <- baseline_DF_500_prop_ev
baseline_DF_500_prop_model$`Sample ID` <- NULL
baseline_DF_500_prop_model$Week <- NULL
baseline_DF_500_prop_model$proportion <- NULL
baseline_DF_500_prop_model$Diagnosis <- NULL
baseline_DF_500_prop_model$viral_copies_B26 <- NULL
baseline_DF_500_prop_model$viral_copies_B380 <- NULL
baseline_DF_500_prop_model$viral_copies_B85 <- NULL
baseline_DF_500_prop_model$viral_copies_B137 <- NULL
baseline_DF_500_prop_model$viral_copies_B261 <- NULL

# Convert predicting class variables to binomial
table(baseline_DF_500_prop_model$`Endoscopic outcome`)
baseline_DF_500_prop_model$`Endoscopic outcome`[baseline_DF_500_prop_model$`Endoscopic outcome` == "Yes"] <- 1
baseline_DF_500_prop_model$`Endoscopic outcome`[baseline_DF_500_prop_model$`Endoscopic outcome` == "No"] <- 0
baseline_DF_500_prop_model$`Endoscopic outcome` <- as.numeric(baseline_DF_500_prop_model$`Endoscopic outcome`)
str(baseline_DF_500_prop_model)
table(baseline_DF_500_prop_model$`Endoscopic outcome`)

## Contigency table
baseline_DF_500_prop_model$prediction <- baseline_DF_500_prop_model$log_proportion
baseline_DF_500_prop_model$prediction[baseline_DF_500_prop_model$prediction > 1] <- "qPCR remission"
baseline_DF_500_prop_model$prediction[!baseline_DF_500_prop_model$prediction > 1] <- "qPCR non-remission"
baseline_DF_500_prop_model$`Endoscopic outcome`[baseline_DF_500_prop_model$`Endoscopic outcome` == 1] <- "endoscopic remission"
baseline_DF_500_prop_model$`Endoscopic outcome`[baseline_DF_500_prop_model$`Endoscopic outcome` == 0] <- "endoscopic non-remission"
baseline_DF_500_prop_model$log_proportion <- NULL
baseline_DF_500_prop_model <- baseline_DF_500_prop_model[c(2,1)]

baseline_DF_500_prop_model_table <- table(baseline_DF_500_prop_model)
baseline_DF_500_prop_model_tableDF <- as.data.frame(prop.table(baseline_DF_500_prop_model_table, margin =2))
baseline_DF_500_prop_model_tableDF$perc <- baseline_DF_500_prop_model_tableDF$Freq*100
chisq.test(baseline_DF_500_prop_model_table)
# Check with chisq if the two variables are independent.
# Meaning is endoscopic remission different significantly different in qPCR non-remission. group then remission.

View(baseline_DF_500_prop_model_table)
## or calcauclate true positives?  false positives?
# specificy, sensitivty?
# check: https://online.stat.psu.edu/stat507/lesson/10/10.3
# like this: https://www.jmp.com/en_us/statistics-knowledge-portal/exploratory-data-analysis/mosaic-plot.html

## put TP first, right FP, FP, TN 

# 66.6 qPCR NR = acutal NR; 66.6% qPCR R = actual R
mosaicplot(baseline_DF_500_prop_model_table, color = T)

str(baseline_DF_500_prop_model_tableDF)
esquisser(baseline_DF_500_prop_model_tableDF)

ggplot(baseline_DF_500_prop_model_tableDF) +
 aes(x = Endoscopic.outcome, y = perc, fill = prediction) +
 geom_col() +
 scale_fill_hue(direction = 1) +
 theme_classic()

View(baseline_DF_500_prop_model)
# Mosaic plot:
## 1) A mosaic plot is a special type of stacked bar chart. For two variables, the width of the columns is proportional to the number of observations in each level of the variable plotted on the horizontal axis. 
## 2) https://rpubs.com/tskam/mosaic_plot
####################################
# 3.2 Build logistic regression model
####################################
glm.fit <- glm(formula = `Endoscopic outcome` ~ log_proportion, family = binomial("logit"), data = baseline_DF_500_prop_model, method = "glm.fit", na.action(na.omit))
summary(glm.fit)
####################################
# 3.3 Characteristics model
####################################
# Analysis of variance for individual terms
Anova(glm.fit, type="II", test="Wald")
summary(glm.fit, test="Chisq")
?Anova

# Pseudo-R-squared
pR2(glm.fit)*100 # McFadden output*100 is best measure for R2
lrtest(glm.fit) # Significance of model

# Calculate relative risk to visualize.
RR <- exp(coef(glm.fit))[-1]
RR
exp(confint(glm.fit))[-1,]

## Simply show the significant ones out of the GLR model. Therefore in this case only show remission
glm.fit$coefficients

## RR figure
variable <- c("Phage markers")
RR_value <- c(1.312019)
lowerCI <- c(1.125002)
upperCI <- c(1.559694)
RR_df <- data.frame(variable, RR_value,lowerCI,upperCI)

ggplot(RR_df) +
  geom_point(aes(x = variable, y = RR_value), size = 4,  color = "black") +
  geom_errorbar(aes(x = variable, y = RR_value, ymin = lowerCI, ymax = upperCI), width = 0, color = "black", ) +
  geom_text(aes(x = variable, y = RR_value, label = round(RR_value, 2)), vjust = - 1.5, size = 2.5) +
  ylim(0.5,1.6) +
  geom_hline(yintercept=1, linetype='dotted', col = 'black') +
  coord_flip() +
  xlab("") +
  ylab("Community typ M Relative risk") +
  theme_bw()

## Phage markers 
# Create figure logistic regression
ggplot(baseline_DF_500_prop_model, aes(x=log_proportion, y=`Endoscopic outcome`)) + 
  stat_smooth(method="glm", se=TRUE, method.args = list(family=binomial("logit")), col="#b87a4d", fill="#b87a4d", lty=1) +
  ylab("Endoscopic remission prevalence") +
  xlab("log10(potential phage markers)") +
  theme_bw() +
  ylim(0,1)

# The larger the logarithm the more change of hosting a positive endoscopic outcome
# Thus, the more + phage/-phages the more likely for patients to lead to endoscopic response
####################################
# 3.4 ROC curve
####################################
# Convert NA's
baseline_DF_500_prop_model_ROC <- baseline_DF_500_prop_model
baseline_DF_500_prop_model_ROC[is.na(baseline_DF_500_prop_model_ROC)] <- 0 
nrow(baseline_DF_500_prop_model_ROC) # 128 samples
# baseline_DF_500_prop_model_ROC <- baseline_DF_500_prop_model_ROC[!baseline_DF_500_prop_model_ROC$log_proportion == 0,]
# nrow(baseline_DF_500_prop_model_ROC) # 82 samples
#View(baseline_DF_500_prop_model_ROC)

# Create vectors
endoscopic <- baseline_DF_500_prop_model_ROC$`Endoscopic outcome`
log_markers <- baseline_DF_500_prop_model_ROC$log_proportion
View(baseline_DF_500_prop_model_ROC)

glm.fit <- glm(formula = endoscopic ~ log_markers, family = binomial("logit"), method = "glm.fit")
summary(glm.fit)
plot(x = log_markers, y = endoscopic)
lines(log_markers, glm.fit$fitted.values)

# Figure
ggplot(baseline_DF_500_prop_model_ROC, aes(x=log_proportion, y=`Endoscopic outcome`)) + 
  stat_smooth(method="glm", se=TRUE, method.args = list(family=binomial("logit")), col="#b87a4d", fill="#b87a4d", lty=1) +
  ylab("Endoscopic remission prevalence") +
  xlab("log10(potential phage markers)") +
  theme_bw() +
  ylim(0,1) 

# ROC curve
par(pty = "s") # square
ROC_potential_biomarkers <- roc(endoscopic, glm.fit$fitted.values, plot = TRUE, legacy.axes = TRUE,
    percent = TRUE, xlab = "False positive rate",
    ylab = "True positive rate", col = "#377eb8", lwd = 2,
    print.auc = TRUE)

pROC_obj <- roc(endoscopic, glm.fit$fitted.values,
                smoothed = TRUE,
                # arguments for ci
                ci=TRUE, ci.alpha=0.9, stratified=FALSE,
                # arguments for plot
                plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
                print.auc=TRUE, show.thres=TRUE)

sens.ci <- ci.se(pROC_obj)
plot(sens.ci, type="shape", col="lightblue")

plot(sens.ci, type="bars")

ROC_potential_biomarkers
mean(ROC_potential_biomarkers$sensitivities)
mean(ROC_potential_biomarkers$specificities)

View(Mastertable_viral_rarefied)
View(Mastertable_viral_rarefied[Mastertable_viral_rarefied$names == "NODE_1_length_98031_cov_823_635982_B26",])

# B26 circular
# B380 circular
# B261 circular
  
# B85 not circular
# B137 not circular

# In general, an AUC of 0.5 suggests no discrimination (i.e., ability to diagnose patients with and without the disease or condition based on the test), 0.7 to 0.8 is considered acceptable, 0.8 to 0.9 is considered excellent, and more than 0.9 is considered outstanding.
####################################

######### !!!!!!!!!!!!!!! #########
# SAVE Global Environment #