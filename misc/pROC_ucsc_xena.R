library(tidyverse)
library(stringr)
library(UCSCXenaTools)
library(pROC)


# Codes of the samples
# https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/sample-type-codes

# Step 1: Load Xena Datasets
data("XenaData")
datasets <- XenaData

# Step 2: Check the datasets
datasets[grep("LUAD", datasets$XenaCohorts), ]  

# Step 3: Download Gene Expression Data
# HiSeqV2: Normalization only LUAD: RSEM
luad_expression <- XenaGenerate(subset = XenaData$XenaCohorts == "TCGA Lung Adenocarcinoma (LUAD)") %>%
  XenaFilter(filterDatasets = "TCGA.LUAD.sampleMap/HiSeqV2") %>%  # Change if using different expression dataset
  XenaQuery() %>%
  XenaDownload() %>%
  XenaPrepare()

# Step4: Check if DNAJC9 is present in the dataset
#gene_data
gene_data <- luad_expression$HiSeqV2.gz %>%
  filter(sample == "DNAJC9") %>%
  select(-sample) %>%
  pivot_longer(cols = everything(), names_to = "TCGA_ID", values_to = "Expression") %>%
  mutate(Expression = as.numeric(Expression),
         category = case_when(
           str_detect(TCGA_ID, "-01$") ~ "Tumor",  # Assign "Tumor" for samples ending in -01
           str_detect(TCGA_ID, "-11$") ~ "Normal", # Assign "Normal" for samples ending in -11
           TRUE ~ NA_character_))

plot.roc(gene_data$category, gene_data$Expression,          # data
         percent = TRUE,                    # show all values in percent
         partial.auc=c(100, 90), 
         partial.auc.correct=TRUE,          # define a partial AUC (pAUC)
         print.auc=TRUE,                    
         #display pAUC value on the plot with following options:
         print.auc.pattern = " Corrected pAUC \n (100-90%% SP):\n%.1f%%",
         print.auc.col = "#1c61b6",
         auc.polygon = TRUE, 
         auc.polygon.col = "#1c61b6",       # show pAUC as a polygon
         max.auc.polygon = TRUE, 
         max.auc.polygon.col = "#1c61b622", # also show the 100% polygon
         main = "Partial AUC (pAUC)")


# To remove the max AUC polygon as well:
plot.roc(gene_data$category, gene_data$Expression,
         percent = TRUE,
         partial.auc = c(100, 90),
         partial.auc.correct = TRUE,
         print.auc = TRUE,
         print.auc.pattern = "Corrected pAUC \n (100-90%% SP):\n%.1f%%",
         print.auc.col = "#1c61b6",
         auc.polygon = FALSE,  # Remove the pAUC polygon/rectangle
         #auc.polygon.col = "#1c61b6",  # No longer needed
         max.auc.polygon = FALSE,  # Remove the max AUC polygon
         #max.auc.polygon.col = "#1c61b622", # No longer needed
         main = "Partial AUC (pAUC)")
