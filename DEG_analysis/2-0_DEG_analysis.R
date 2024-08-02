# Script to pseudo bulk the raw sequencinig count, we are using the non-normalised spatial counts to prevent
# for normalising the dataset twice. We pseudo-bulk the counts per cluster and perform a differential expression
# analysis with limma
# By: Joy otten & Shu Dan


# Libaries:
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(Matrix)
library(dplyr)
library(stringr)
library(tidyverse)
library(DESeq2)
library(RColorBrewer)
library(leiden)
library(igraph)
library(pheatmap)
library(clustree)
library(KEGGprofile)
library(clusterProfiler)
library(DOSE)
library(biomaRt)
library(KEGGREST)
library(msigdbr)
library(PCAtools)
library(plyr)
library(enrichR)
library(sva)
library(edgeR)
library(QCEWAS)
library(bacon)
library(RUVSeq)
set.seed(1234)

# Data path
data_path <- "set your path"
output_path <- "set your path"
functions_path <- "set your path"
output_workenvironment <- "set your path"

# Load in functions
source(paste0(functions_path, "/General_functions.R"))
ggPCA <- function(pca, metadata, variables){
  PCA_out <- as.data.frame(pca$x)
  percentage <- round(pca$sdev / sum(pca$sdev) * 100, 2)
  percentage <- paste0(colnames(PCA_out), " (", percentage, "%", ")")
  theme <- theme(panel.background = element_blank(), 
                 panel.border = element_rect(fill = NA), 
                 panel.grid.major = element_blank(), 
                 panel.grid.minor = element_blank(), 
                 strip.background = element_blank(), 
                 axis.text.x = element_text(colour = "black"), 
                 axis.text.y = element_text(colour = "black"),
                 axis.ticks = element_line(colour = "black"), 
                 plot.margin = unit(c(1, 1, 1, 1), "line"))
  for(i in variables){
    p <- ggplot(PCA_out, aes(x = PC1, y = PC2, color = metadata[, i])) + 
      geom_point() + theme + xlab(percentage[1]) + ylab(percentage[2]) + 
      labs(color = i)
    print(p)
  }
}
source(paste0(functions_path, "/customPCA.R"))

# Data obtained from the QC_analysis.R script
total <- readRDS(paste0(output_workenvironment, "/total.RDS"))
sample_name <- c("35.1", "35.2", "35.3", "35.4", "61.1", "61.2", "61.3", "61.4",
                 "69.1", "69.2", "69.3", "69.4", "79.1", "79.2", "79.3", "79.4")
# Split the data from Visium in the clusters that we observed
samples <- SplitObject(total, split.by = "orig.ident")
xsample_name <- paste0("X", sample_name, sep = "")
all_genes <- rownames(total)

# Load metadata
metadata <- read.csv("/PHShome/je637/Visium/data/metadata_Visium.csv", header = TRUE, sep = ",")
rownames(metadata) <- xsample_name
metadata$group <- relevel(as.factor(metadata$group), ref = "CTRL")
group <- relevel(as.factor(metadata$group), ref = "CTRL")
contrast <- c("group", "FC", "CTRL")

# Subclusters
# ------------
# we set the 
layers <- as.character(unique(total$sub_regions)) # we set the layers to the sub-clusters
test_methods <- "pseudo_bulk"

# pseudo bulk per cluster
bulk_count_all <- list()
for (i in test_methods){
  bulk_count <- list()
  message(i)
  
  for(j in sample_name){
    message(j)
    
    # strip sample name
    temp_name <- strsplit(j, '.', fixed = TRUE)
    name_i <- temp_name[[1]][1]
    panel_i <- temp_name[[1]][2]
    image_name_i <- paste("slide", name_i, "_panel_", panel_i, sep = "")
    samp_name_i <- paste("slide", name_i, "_panel", panel_i, sep = "")
    
    # subset samples
    message("subset samples")
    ICA_sample <- samples[[samp_name_i]]
    Idents(ICA_sample) <- ICA_sample$sub_regions
    sample_sc <- ICA_sample
    sample_sc$layer_assign <- sample_sc$sub_regions
    
    # Visualization of the layers
    
    # Pseudo bulk the counts in each layers
    message("pseudo bulking counts for each layer")
    sample_bulk_count <- pseudo_bulk_raw(layers, sample_sc, sample_name = j, verbose = FALSE)
    bulk_count <- append(bulk_count, sample_bulk_count)
  }
  
  bulk_count_all <- append(bulk_count_all, list(bulk_count))
}
names(bulk_count_all) <- test_methods
data <- as.data.frame(bulk_count_all[["pseudo_bulk"]])

# Checking the data for co-variables
#---------
## To check the how different the sequencing depths are per cluster, this is due to
## that each region contains a different number of spots and therefor it leads to more/less sequencing reads
seq_depth <- colSums(data)
barplot(seq_depth)

# Create metadata dataframe from column names
metadata <- as.data.frame(colnames(data))
colnames(metadata) <- "samples"

# Initialize metadata columns with default values
metadata$region <- "NA"
metadata$sample <- "NA"
metadata$slide_nr <- "NA"
metadata$group <- "NA"

# Populate 'region' and 'sample' columns
for (i in metadata$samples) {
  message(i)
  x <- str_split(i, "_")[[1]]  # Split sample name by '_'
  y <- which(metadata$samples == i)
  metadata[y, "region"] <- x[1]
  metadata[y, "sample"] <- x[2]
}

# Assign slide number based on the starting digits of 'sample'
metadata$slide_nr <- case_when(
  str_starts(metadata$sample, "35") ~ "1",
  str_starts(metadata$sample, "61") ~ "2",
  str_starts(metadata$sample, "69") ~ "3",
  str_starts(metadata$sample, "79") ~ "4",
  TRUE ~ metadata$slide_nr
)

# Assign group based on the ending digits of 'sample'
metadata$group <- case_when(
  str_ends(metadata$sample, "1") ~ "FC",
  str_ends(metadata$sample, "2") ~ "CTRL",
  str_ends(metadata$sample, "3") ~ "FC",
  str_ends(metadata$sample, "4") ~ "CTRL",
  TRUE ~ metadata$group
)

# Add sequencing depth column
metadata$seq_depth <- seq_depth

# PCA analysis to detect co-variables
PCA1 <- prcomp(t(data))
ggPCA(PCA1, metadata, colnames(metadata)[2:6])
log_counts <- log2(data +1)
custom.PCA(beta = log_counts, pd = metadata, adjust.p = TRUE)
# The most variation in the data is explained by the different regions that we detect 
# and the sequencing depth. Therefore we normalise the data per cluster. 
# Additionally at a pc of 17 we detect 0.16% variation in the data due to the spatial slides

# PCA on pseudo-bulk per region this is in order to check if you perform the DEG
# analysis seperately if there are additional co-variables significantly contributing to 
# the variation in the data.
layers <- str_replace(layers, "-", ".")
layers <- as.character(layers)
# Load metadata
metadata <- read.csv(paste0(data_path, "/metadata_Visium.csv"), header = TRUE, sep = ",")
rownames(metadata) <- xsample_name
metadata$group <- relevel(as.factor(metadata$group), ref = "CTRL")
group <- relevel(as.factor(metadata$group), ref = "CTRL")
contrast <- c("group", "FC", "CTRL")

pdf(paste0(image_path,"/limma_layers.pdf"))
for(c in layers){
  layer_count <- log_counts %>% dplyr::select(matches(paste("X",c, "_", sep = "")))
  custom.PCA(beta = layer_count, pd = metadata, adjust.p = TRUE, plot.title = c)
}
dev.off()
# ------ 

# Differential expression analysis with Limma
#-------
# replace layer and sample names for easier downstream analysis [numbers as sample name is not the best]
test_layers2p <- gsub(" ", ".", layers) 
test_layers2p <- gsub("/",".",test_layers2p)
xsample_name <- paste0("X", sample_name, sep = "")

# Load metadata
metadata <- read.csv(paste0(data_path, "/metadata_Visium.csv"), header = TRUE, sep = ",")
rownames(metadata) <- xsample_name
metadata$group <- relevel(as.factor(metadata$group), ref = "CTRL")
group <- relevel(as.factor(metadata$group), ref = "CTRL")
contrast <- c("group", "FC", "CTRL")

ICA_bulk_count <- as.data.frame(bulk_count_all[["pseudo_bulk"]])
layers <- str_replace(layers, "-", ".")
layers <- as.character(layers)
norm_counts <- list()
ICA_limma <- list()
for (c in layers){
  message(c)
  pdf(paste0(image_path, "/new_limma_", c, ".pdf"))
  layer_count <- ICA_bulk_count %>% dplyr::select(matches(paste("X",c, "_", sep = "")))
  colnames(layer_count) <- paste("X", sub(".*_", "", colnames(layer_count)), sep = "")
  rownames(layer_count) <- all_genes
  # if there's no spot detected in a layer for a sample, we exclude them from the metadata and group variables. 
  if (length(colnames(layer_count)) != length(rownames(metadata))){
    message("layer is not detected in all samples")
    metadata_layer <- metadata[match(colnames(layer_count), rownames(metadata)),]
    metadata_layer$group <-  relevel(as.factor(metadata_layer$group), ref = "CTRL")
  } else { 
    metadata_layer <- metadata
  }
  # Set variables for the limma analysis
  slide <- as.numeric(metadata_layer$slide)
  group <- metadata_layer$group
  design <- model.matrix(~group+slide) # We correct over here for the slide effect found earlier with the PCA analysis
  
  message("limma")
  # Normalisation is performed within the de_limma function with EdgeR
  limma_res <- de_limma(layer_count, group, design, coef = "groupFC",
                        weights = FALSE, contrast = NULL, plot = TRUE)
  # Bacon analysis to correct for inflation of p-values
  limma_output <- bacon_limma(limma_res$efit, limma_res$res, coef = "groupFC", plot = TRUE)
  limma_plots(efit = limma_res$efit, res = limma_output$res, d_norm = limma_res$d_norm, group = group)
  ICA_limma <- append(ICA_limma, list(limma_output))
  limma_norm <- as.data.frame(limma_res[["voom_count"]][["E"]])
  norm_counts <- append(norm_counts, list(limma_norm))
  
  # add title to the page for pdf
  mtext(paste("cluster", c, sep = ""), outer=TRUE,  cex=1, line = -0.5)
  dev.off()
}
dev.off()
names(ICA_limma) <- layers
names(norm_counts) <- layers

saveRDS(norm_counts, paste0(output_workenvironment, "/new_count_sub_regions_raw.RDS"))
saveRDS(ICA_limma, paste0(output_workenvironment, "/new_limmaresults_regions_raw.RDS"))


######
### Limma analysis with RUVseq variance testing before to detect outliers
#######

# Data obtained from the QC_analysis.R script
total <- readRDS(paste0(data_path, "/total.RDS"))

sample_name <- c("35.1", "35.2", "35.3", "35.4", "61.1", "61.2", "61.3", "61.4",
                 "69.1", "69.2", "69.3", "69.4", "79.1", "79.2", "79.3", "79.4")
# Split the data from Visium in the clusters that we observed
samples <- SplitObject(total, split.by = "orig.ident")
xsample_name <- paste0("X", sample_name, sep = "")
all_genes <- rownames(total)

# Load metadata
metadata <- read.csv("/PHShome/je637/Visium/data/metadata_Visium.csv", header = TRUE, sep = ",")
rownames(metadata) <- xsample_name
metadata$group <- relevel(as.factor(metadata$group), ref = "CTRL")
group <- relevel(as.factor(metadata$group), ref = "CTRL")
contrast <- c("group", "FC", "CTRL")

layers <- as.character(unique(total$sub_regions))
test_methods <- "pseudo_bulk"

# pseudo bulk per cluster
bulk_count_all <- list()
for (i in test_methods){
  bulk_count <- list()
  message(i)
  
  for(j in sample_name){
    message(j)
    
    # strip sample name
    temp_name <- strsplit(j, '.', fixed = TRUE)
    name_i <- temp_name[[1]][1]
    panel_i <- temp_name[[1]][2]
    image_name_i <- paste("slide", name_i, "_panel_", panel_i, sep = "")
    samp_name_i <- paste("slide", name_i, "_panel", panel_i, sep = "")
    
    # subset samples
    message("subset samples")
    ICA_sample <- samples[[samp_name_i]]
    Idents(ICA_sample) <- ICA_sample$sub_regions
    sample_sc <- ICA_sample
    sample_sc$layer_assign <- sample_sc$sub_regions
    
    # Visualization of the layers
    
    # Pseudo bulk the counts in each layers
    message("pseudo bulking counts for each layer")
    sample_bulk_count <- pseudo_bulk_raw(layers, sample_sc, sample_name = j, verbose = FALSE)
    bulk_count <- append(bulk_count, sample_bulk_count)
  }
  
  bulk_count_all <- append(bulk_count_all, list(bulk_count))
}
names(bulk_count_all) <- test_methods
data <- as.data.frame(bulk_count_all[["pseudo_bulk"]])
seq_depth <- colSums(data)
barplot(seq_depth)

# Load metadata
metadata <- read.csv(paste0(data_path, "/metadata_Visium.csv"), header = TRUE, sep = ",")
rownames(metadata) <- xsample_name
metadata$group <- relevel(as.factor(metadata$group), ref = "CTRL")
group <- relevel(as.factor(metadata$group), ref = "CTRL")
contrast <- c("group", "FC", "CTRL")

ICA_bulk_count <- as.data.frame(bulk_count_all[["pseudo_bulk"]])
layers <- str_replace(layers, "-", ".")
layers <- as.character(layers)
norm_counts <- list()
ICA_limma <- list()

for (c in layers){
  message(c)
  pdf(paste0(output_path, "/limma_", c, ".pdf"))
  layer_count <- ICA_bulk_count %>% dplyr::select(matches(paste("X",c, "_", sep = "")))
  colnames(layer_count) <- paste("X", sub(".*_", "", colnames(layer_count)), sep = "")
  rownames(layer_count) <- all_genes
  # if there's no spot detected in a layer for a sample, we exclude them from the metadata and group variables. 
  if (length(colnames(layer_count)) != length(rownames(metadata))){
    message("layer is not detected in all samples")
    metadata_layer <- metadata[match(colnames(layer_count), rownames(metadata)),]
    metadata_layer$group <-  relevel(as.factor(metadata_layer$group), ref = "CTRL")
  } else { 
    metadata_layer <- metadata
  }
  slide <- as.numeric(metadata_layer$slide)
  group <- metadata_layer$group
  design <- model.matrix(~group+slide)
  
  message("limma") # be aware that within the limma function the data will be normalised by EdgeR
  limma_res <- de_limma(layer_count, group, design, coef = "groupFC",
                        weights = FALSE, contrast = NULL, plot = TRUE, metadata_layer = metadata_layer)
  # Bacon is used to adjust for inflation of p-values
  limma_output <- bacon_limma(limma_res$efit, limma_res$res, coef = "groupFC", plot = TRUE)
  limma_plots(efit = limma_res$efit, res = limma_output$res, d_norm = limma_res$d_norm, group = group)
  ICA_limma <- append(ICA_limma, list(limma_output))
  limma_norm <- as.data.frame(limma_res[["voom_count"]][["E"]])
  set_limma <- newSeqExpressionSet(as.matrix(limma_norm), phenoData = metadata_layer)
  plotRLE(set_limma, outline = FALSE, ylim=c(-4,4))
  norm_counts <- append(norm_counts, list(limma_norm))
 
  # add title to the page for pdf
  mtext(paste("cluster", c, sep = ""), outer=TRUE,  cex=1, line = -0.5)
  dev.off()
}
dev.off()
# The RTSNE plots show that there is more variance in sample 69.1. This might indicate that sample 69.1 is an outlier.
# However, we decided to leave it in since if the sample would be fastly different there would have been additional clusters been present. 


