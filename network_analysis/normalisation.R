# Normalizing the data for all subclusters together
# By Joy Otten

# Libraries:
library(cowplot)
library(tidyverse)
library(ggridges)
library(scales)
library(SingleCellExperiment)
library(Matrix)
library(knitr)
library(readxl)
library(writexl)
library(Seurat)
library(msigdbr)
library(clusterProfiler)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(enrichplot)
library(dplyr)
library(limma)
library(edgeR)
set.seed(1234)

# Data path
data_path <- "set your path"
output_path <- "set your path"
functions_path <- "set your path"
output_workenvironment <- "set your path"

# Load in data
total <- readRDS(paste0(output_workenvironment, "/total.RDS"))

regions <- SplitObject(total, split.by = "sub_regions")
samples <- SplitObject(total, split.by = "orig.ident")
sample_name <- c("35.1", "35.2", "35.3", "35.4", "61.1", "61.2", "61.3", "61.4",
                 "69.1", "69.2", "69.3", "69.4", "79.1", "79.2", "79.3", "79.4")
names(samples) <- sample_name

count_matrix <- DataFrame()
for(i in names(samples)){
  message(i)
  
  data <- samples[[i]]
  r <- SplitObject(data, split.by = "sub_regions")
  
  for(region in names(r)){
    name <- paste0(i, "_", region)
    reg <- r[[region]]
    all_count <- as.data.frame(reg@assays[["Spatial"]]@counts)
    count_matrix[name] <- as.data.frame(apply(all_count, 1, sum))
  }
}
count_matrix <- as.data.frame(count_matrix)

# Calculate the normalisation factors with Deseq2 and perform normalisation
dge2 <- DGEList(counts=count_matrix)
dge2 <- calcNormFactors(dge2, method = "TMM")
normalised_data <- cpm(dge2, log=TRUE, prior.count=0.5) 

saveRDS(normalised_data, paste0(output_workenvironment, "/visium_normalised_entire_brain.RDS"))

# Create a metadata file per region, subject and condition
# Formatting metadata
metadata <- data.frame(samples = colnames(logCPM_offset), stringsAsFactors = FALSE)
metadata$sample <- "NA"
metadata$region <- "NA"
metadata$slide_nr <- "NA"
metadata$group <- "NA"

# Extract sample and region information
for (i in seq_len(nrow(metadata))) {
  sample_info <- str_split(metadata$samples[i], "_")[[1]]
  metadata$sample[i] <- sample_info[1]
  metadata$region[i] <- sample_info[2]
}

# Assign slide numbers based on sample prefix
metadata$slide_nr <- case_when(
  str_starts(metadata$sample, "X35") ~ "1",
  str_starts(metadata$sample, "X61") ~ "2",
  str_starts(metadata$sample, "X69") ~ "3",
  str_starts(metadata$sample, "X79") ~ "4",
  TRUE ~ "NA"
)

# Assign group based on sample suffix
metadata$group <- case_when(
  str_ends(metadata$sample, "1") ~ "FC",
  str_ends(metadata$sample, "2") ~ "CTRL",
  str_ends(metadata$sample, "3") ~ "FC",
  str_ends(metadata$sample, "4") ~ "CTRL",
  TRUE ~ "NA"
)

# Calculate sequencing depth
metadata$seq_depth <- colSums(logCPM_offset)

# Save the formatted metadata
saveRDS(metadata, paste0(data_path, "/visium_metadata_normalised_entire_brain.RDS"))


