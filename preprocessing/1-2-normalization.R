# QC analysis and determining the number of PC's to use wherein the Hbb, Hba and Ttr genes are deleted on the broad regions
# Additionally the first clustering is performed.
# Joy Otten

# Libaries:
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(Matrix)
library(stringr)
library(tidyverse)
library(DESeq2)
library(RColorBrewer)
library(leiden)
library(igraph)
library(pheatmap)
library(clustree)
library(clusterProfiler)
library(DOSE)
library(biomaRt)
library(KEGGREST)
library(msigdbr)
library(plyr)
library(enrichR)
library(future)
library(CellChat)
set.seed(1234)

# data paths
data_path <- "set your path"
output_path <- "set your path"
functions_path <- "set your path"
output_workenvironment <- "set your path"

# Data
load(paste0(data_path, "/Visium-preprocessing.RData"))
sample_name <- c("35.1", "35.2", "35.3", "35.4", "61.1", "61.2", "61.3", "61.4",
                 "69.1", "69.2", "69.3", "69.4", "79.1", "79.2", "79.3", "79.4")
source(paste0(functions_path, "/DE_res_0*4_functions.R"))
source(paste0(functions_path, "/custom.PCA.R"))
# Spatial plots of certain genes.

# Perform SCTransform normalization on the Spatial assay and retain all genes
total <- SCTransform(total, assay = "Spatial", return.only.var.genes = FALSE, verbose = FALSE)
# Run PCA on the SCT assay
total <- RunPCA(total, assay = "SCT", verbose = FALSE)

# Plot the elbow plot for the first 30 principal components
ElbowPlot(total, ndims = 30)
# Find neighbors using the first 30 principal components from the PCA
total <- FindNeighbors(total, reduction = "pca", dims = 1:30)
# Calculate the percentage of variance explained by each PC
pct <- total[["pca"]]@stdev / sum(total[["pca"]]@stdev) * 100
# Calculate cumulative percentage of variance explained
cumu <- cumsum(pct)
# Identify the first PC where the cumulative percentage exceeds 90% 
# and the percentage of variation is less than 5%
co1 <- which(cumu > 90 & pct < 5)[1]
# Identify the first PC where the change in percentage variation exceeds 0.1%
co2 <- sort(which((pct[1:(length(pct) - 1)] - pct[2:length(pct)]) > 0.1), decreasing = TRUE)[1] + 1
# Select the minimum of the two calculated PCs
pcs <- min(co1, co2)
pcs

# This indicates that the dimensions we want to use is 18
total <- FindNeighbors(total, reduction = "pca", dims = 1:18)
VizDimLoadings(total, dims = 1:18, nfeatures = 10, reduction = "pca", combine = FALSE, col = 'black')
# according to the gene loading plots it seems that certain genes contribute to a high extent to the PC's, such as the Hbb genes.
# I want to check if these genes have a certain spatial pattern or not.

# Spatial plots of certain genes.
# Define the list of genes to plot
genes <- c("Ttr", "Hbb-bt", "Hbb-bs", "Hbb-y", "Hba-a1", "Hba-a2")
# Define custom color scales for specific genes
color_scales <- list(
  "Hbb-bt" = list(limits = c(0, 5), breaks = c(0, 2.5, 5), colours = c("grey60", "steelblue", "tomato")),
  "Hbb-bs" = list(limits = c(0, 10), breaks = c(0, 5, 10), colours = c("grey60", "steelblue", "tomato")),
  "Hbb-y"  = list(limits = c(0, 1), breaks = c(0, 0.5, 1), colours = c("grey60", "steelblue", "tomato")),
  "Hba-a1" = list(limits = c(0, 10), breaks = c(0, 5, 10), colours = c("grey60", "steelblue", "tomato")),
  "Hba-a2" = list(limits = c(0, 10), breaks = c(0, 5, 10), colours = c("grey60", "steelblue", "tomato")),
  "Ttr"    = list(limits = c(0, 10), breaks = c(0, 5, 10), colours = c("grey60", "steelblue", "tomato"))
)

# Iterate over each sample
for (j in sample_name) {
  message(j)
  
  # Extract sample and panel names
  temp_name <- strsplit(j, '.', fixed = TRUE)[[1]]
  name_i <- temp_name[1]
  panel_i <- temp_name[2]
  
  # Construct the image name
  image_name_i <- paste("slide", name_i, "_panel_", panel_i, sep = "")
  
  # Create a PDF to save plots for this sample
  pdf(file.path("/PHShome/je637/Visium/output/deleted_genes/", paste0(image_name_i, ".pdf")))
  
  # Iterate over each gene to plot
  for (i in genes) {
    # Generate spatial plot for the current gene
    p1 <- SpatialPlot(total, features = i, images = image_name_i)
    
    # Apply custom color scale if defined
    if (i %in% names(color_scales)) {
      p1 <- p1 + scale_fill_gradientn(
        limits = color_scales[[i]]$limits,
        breaks = color_scales[[i]]$breaks,
        colours = color_scales[[i]]$colours
      )
    }
    
    # Plot the result
    plot(p1)
  }
  
  # Close the PDF device
  dev.off()
}
# As you can see the Hbb genes do not show a spatial coherent pattern and therefore we deleted
# those genes in our data manually. Futhermore, for the Ttr gene you see that it is smeared out over the sample 
# which represents a technical artificat related to the slicing direction of the sample.

# Deletion of the Hbb-genes and Ttr
# Hbb-genes interfere with clustering which is non-specific therefore they are deleted. 
# Ttr gene semms to have a clustering effect which is non-specific therefore this is deleted (mainly weird signature in the cortical area)
counts <- GetAssayData(total, assay = "Spatial")
counts <- counts[-(which(rownames(counts) %in% c("Hbb-bs", "Hbb-bt", "Hbb-y", "Hba-a2", "Hba-a1", "Ttr"))),]
total <- subset(total, features = rownames(counts))
rm(counts)

# Normalization and clustering with PCA and the Leiden algorithm, we need to this again since we removed genes
total <- SCTransform(total, assay = "Spatial", return.only.var.genes = FALSE, verbose = FALSE)
total <- RunPCA(total, assay = "SCT", verbose = FALSE)

# Generate Elbow Plot for the first 30 principal components
ElbowPlot(total, ndims = 30)

# Find neighbors using the first 30 principal components from PCA
total <- FindNeighbors(total, reduction = "pca", dims = 1:30)
# Calculate the percentage of variance explained by each PC
pct <- total[["pca"]]@stdev / sum(total[["pca"]]@stdev) * 100
# Calculate the cumulative percentage of variance explained
cumu <- cumsum(pct)
# Find the first PC where cumulative variance exceeds 90% and individual variance is less than 5%
co1 <- which(cumu > 90 & pct < 5)[1]
# Find the last point where the change in percentage variance exceeds 0.1%
co2 <- which(diff(pct) > -0.1)[1] + 1
# Select the minimum of the two identified PCs
pcs <- min(co1, co2)
# Output the selected number of PCs
pcs

# Plot the variance explained by the PC's
# https://bioconductor.org/packages/release/bioc/vignettes/PCAtools/inst/doc/PCAtools.html
pca_data <- total[["pca"]]
screeplot(pca_data)

total <- FindNeighbors(total, reduction = "pca", dims = 1:18)
total <- FindClusters(total, algorithm = 4, verbose = FALSE, method = "igraph",
                      resolution = c(0, 0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 1, 1.2, 1.8))
saveRDS(total, paste0(data_path, "/total.RDS"))
        
# Umap plot visualisation of the clusters based on PCA analysis
# Define custom colors
my_cols <- c('3'='#F68282', '15'='#31C53F', '5'='#1FA195', '1'='#B95FBB', '13'='#D4D915',
             '14'='#28CECA', '9'='#ff9a36', '8'='#2FF18B', '11'='#aeadb3', '6'='#faf4cf',
             '2'='#CCB1F1', '12'='#25aff5', '7'='#A4DFF2', '4'='#4B4BF7', '16'='#AC8F14',
             '10'='#E6C122', '17'='#A6CEE3', '18'='#85B9D7', '19'='#64A4CC', '20'='#448FC0',
             '21'='#237AB5', '22'='#3D8DAB', '23'='#60A6A1', '24'='#84BF97', '25'='#A7D78C',
             '26'='#AF9C6F', '27'='#E42521', '0'='#5EB54C', '29'='#F26D6D', '30'='#FE9D35',
             '31'='#F765A3', '32'='#F68282', '33'='#B95FBB', '34'='#2FF18B', '35'='#448FC0',
             '36'='#B95FBB', '37'='#AC8F14', '38'='#D4D915', '39'='#25aff5', '40'='#CCB1F1',
             '41'='#25aff5', '42'='#A4DFF2', '43'='#4B4BF7', '44'='#AC8F14', '45'='#E6C122',
             '46'='#A6CEE3', '47'='#85B9D7', '48'='#64A4CC', '49'='#448FC0', '50'='#C0F6FD',
             '51'='#FDC0D8')

# Run UMAP on PCA-reduced data
total <- Seurat::RunUMAP(total, reduction = "pca", dims = 1:18)
# Set identities for clusters
Idents(total_new) <- total$SCT_snn_res.0.2
# Plot UMAP
DimPlot(total, reduction = "umap", label = FALSE, cols = my_cols)

# Split total object by the samples
samples <- SplitObject(total, split.by = "orig.ident")

# Data with dimensionality reduction of 18 identifies at a resolution of 0.2, 11 clusters
# -----------
# Assess the different resolution parameters with Clustree
# Reference: https://cran.r-project.org/web/packages/clustree/vignettes/clustree.html
# With clustree we can assess the clustering parameters and find out which clustering parameter
# is the best to use for further processing
p2 <- clustree(total_new, prefix = 'SCT_snn_res.')
plot(p2)
# In the clustree plot you can visualize that a resolution of 0.1 is the best since if you use a 
# bigger resolution you will get over clustering which will lead to overfitting and clusters that are 
# not present in the data. Although, resolution of 0.2 up until 0.4 are also accepted.

# ---------------
# Different resolutions
# ----------
color.use <- ggPalette(length(unique(clusters)))

# Resolution 0.1
pdf(paste0(output_path, "/resolution_0*1.pdf"))
for(j in sample_name){
  message(j)
  temp_name <- strsplit(j, '.', fixed = TRUE)
  name_i <- temp_name[[1]][1]
  panel_i <- temp_name[[1]][2]
  image_name_i <- paste("slide", name_i, "_panel_", panel_i, sep = "")
  samp_name_i <- paste("slide", name_i, "_panel", panel_i, sep = "")
  sample <- samples[[samp_name_i]]
  Idents(sample) <- sample$SCT_snn_res.0.1
  
  p1 <- SpatialDimPlot(sample, images = image_name_i, cols = color.use)
  plot(p1)
}
dev.off()

#Resolution 0.2
pdf(paste0(output_path, "/resolution_0*2.pdf"))
for(j in sample_name){
  message(j)
  temp_name <- strsplit(j, '.', fixed = TRUE)
  name_i <- temp_name[[1]][1]
  panel_i <- temp_name[[1]][2]
  image_name_i <- paste("slide", name_i, "_panel_", panel_i, sep = "")
  samp_name_i <- paste("slide", name_i, "_panel", panel_i, sep = "")
  sample <- samples[[samp_name_i]]
  Idents(sample) <- sample$SCT_snn_res.0.2
  
  p1 <- SpatialDimPlot(sample, images = image_name_i, cols = color.use)
  plot(p1)
}
dev.off()

#Resolution 0.3
pdf(paste0(output_path, "/resolution_0*3.pdf"))
for(j in sample_name){
  message(j)
  temp_name <- strsplit(j, '.', fixed = TRUE)
  name_i <- temp_name[[1]][1]
  panel_i <- temp_name[[1]][2]
  image_name_i <- paste("slide", name_i, "_panel_", panel_i, sep = "")
  samp_name_i <- paste("slide", name_i, "_panel", panel_i, sep = "")
  sample <- samples_new[[samp_name_i]]
  Idents(sample) <- sample$SCT_snn_res.0.3
  
  p1 <- SpatialDimPlot(sample, images = image_name_i, cols = my_cols)
  plot(p1)
}
dev.off()


#Resolution 0.4
pdf(paste0(output_path, "/resolution_0*4.pdf"))
for(j in sample_name){
  message(j)
  temp_name <- strsplit(j, '.', fixed = TRUE)
  name_i <- temp_name[[1]][1]
  panel_i <- temp_name[[1]][2]
  image_name_i <- paste("slide", name_i, "_panel_", panel_i, sep = "")
  samp_name_i <- paste("slide", name_i, "_panel", panel_i, sep = "")
  sample <- samples_new[[samp_name_i]]
  Idents(sample) <- sample$SCT_snn_res.0.4
  
  p1 <- SpatialDimPlot(sample, images = image_name_i, cols = my_cols)
  plot(p1)
}
dev.off()
# After inspection I would go with a resolution of 0.2 at dims 18.

# ---------
# General information about the clusters
Idents(total) <- total$SCT_snn_res.0.2

## Calculate the number of spots per cluster
n_spots <- FetchData(total,
                     vars = c("orig.ident", "SCT_snn_res.0.2")) %>% 
  dplyr::count(SCT_snn_res.0.2)

p1 <- n_spots %>% ggplot(aes(x = SCT_snn_res.0.2, y = n)) + geom_col() +
  ggtitle("Number of spots in clusters")+
  geom_text(aes(label = n), hjust = 0) +
  coord_flip() +
  xlab("Cluster ID") +
  ylab("Number of spots") + theme_classic()
pdf(paste0(output_path, "/number_of_spots_per_cluster.pdf"))
plot(p1)
dev.off()

## Calculate the number of genes per cluster over all samples
metadata <- total@meta.data

p1 <- metadata %>%
  ggplot(aes(SCT_snn_res.0.2, nFeature_Spatial)) + 
  geom_violin(aes(fill = SCT_snn_res.0.2)) + 
  ylab("number of genes") + xlab("clusters") +
  geom_vline(xintercept = 500, linetype = "dashed") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual(values=my_cols) + theme_classic()
pdf(paste0(output_path, "/number_of_genes_per_cluster.pdf"))
plot(p1)
dev.off()

## Calculate the number of genes per cluster per sample
clusters <- SplitObject(total, split.by = "SCT_snn_res.0.2")
pdf(paste0(output_path, "/number_of_genes_per_cluster_per_sample.pdf"))
for(i in 1:11){
  message(i)
  i <- as.character(i)
  cluster <- clusters[[i]]
  metadata <- cluster@meta.data
  p1 <- metadata %>%
    ggplot(aes(orig.ident, nFeature_Spatial)) + 
    geom_violin(aes(fill = orig.ident)) + 
    ylab("total number of genes") + xlab("samples") +
    geom_vline(xintercept = 500, linetype = "dashed") +
    theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="none") +
    ggtitle(paste0("cluster", i)) 
  plot(p1)
}
dev.off()

## Calculate the number of spots per cluster per sample
pdf(paste0(output_path, "/number_of_spots_per_cluster_per_sample.pdf"))
for(j in sample_name){
  message(j)
  temp_name <- strsplit(j, '.', fixed = TRUE)
  name_i <- temp_name[[1]][1]
  panel_i <- temp_name[[1]][2]
  image_name_i <- paste("slide", name_i, "_panel_", panel_i, sep = "")
  samp_name_i <- paste("slide", name_i, "_panel", panel_i, sep = "")
  sample <- samples[[samp_name_i]]
  n_spots <- FetchData(sample,
                       vars = c("orig.ident", "SCT_snn_res.0.2")) %>% 
    dplyr::count(SCT_snn_res.0.2)
  
  p1 <- n_spots %>% ggplot(aes(x = SCT_snn_res.0.2, y = n)) + geom_col() +
    ggtitle(paste0("Number of spots in clusters in sample", j))+
    geom_text(aes(label = n), hjust = 0) +
    coord_flip() +
    xlab("Cluster ID") +
    ylab("Number of spots") + theme_classic()
  plot(p1)
}
dev.off()

# --------------
# Differential expressed genes per region
# -----------
Idents(total) <- total@meta.data$SCT_snn_res.0.2
ST.markers_0.2 <- FindAllMarkers(total, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
write.xlsx2(ST.markers_0.2, "/PHShome/je637/Visium/output/cluster_specific_markers_broad_regions.xlsx")

ST.markers_0.2 %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top10

pdf(paste0(output_path, "/heatmap.pdf"))
DoHeatmap(total, features = top10$gene, group.colors = my_cols) + NoLegend()
dev.off()

pdf(paste0(output_path, "/heatmap_added_legend.pdf"))
DoHeatmap(total, features = top10$gene, group.colors = my_cols)
dev.off()

# Resolution 0.2
clusters = unique(total$SCT_snn_res.0.2)
DE_list_0.2 <- list()

for(i in clusters){
  x <- which(ST.markers_0.2$cluster == i)
  y <- ST.markers_0.2[x,]
  n <- which(y$p_val_adj <= 0.05)
  y <- y[n,]
  DE_list_0.2 <- append(DE_list_0.2, list(y))
}
names(DE_list_0.2) <- clusters

for(i in clusters){
  DE <- as.data.frame(DE_list_0.2[[i]])
  write.csv(DE, paste0(output_path, "/cluster_", i, ".csv"))
}

# -------------
# Find cluster specific genes
# ----------------------
top <- ST.markers_0.2 %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)
clusters = unique(total$SCT_snn_res.0.2)
pdf(paste0(output_path, "/dotplot_clusters_top30_res0*2.pdf"))
for(i in clusters){
  message(i)
  x <- which(top$cluster == i)
  y <- top[x,]
  p1 <- DotPlot(total, features = y$gene) + RotatedAxis()
  plot(p1)
}
dev.off()
# --------

# Check if there are any batch effects
options(future.globals.maxSize = 30000 * 1024^2) # more memory for each thread
plan(multisession, workers = 12) # 12 parallel

counts <- GetAssayData(total, assay = "Spatial")
test <- as.data.frame(counts)
log_counts <- log2(test + 1)
metadata <- total@meta.data
metadata <- metadata[,c(2:5)]
log_counts <- as.data.frame(log_counts)

PCA <- prcomp(t(counts))
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
pdf(paste0(output_path, "/batches_detected_raw_pseudo_bulked_spatial_data.pdf"))
ggPCA(PCA, metadata, colnames(metadata))
dev.off()

pdf(paste0(output_path, "/test_raw_data.pdf"))
custom.PCA(beta = log_counts, pd = metadata, adjust.p = TRUE)
dev.off()

# Testing batches on the SCT normalised data
counts <- GetAssayData(total, assay = "SCT")
test <- as.data.frame(counts)
log_counts <- log2(test + 1)
metadata <- total@meta.data
metadata <- metadata[,c(2:5)]
log_counts <- as.data.frame(log_counts)

PCA <- prcomp(t(counts))
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
pdf(paste0(output_path, "/batches_detected_SCT_spatial_data.pdf"))
ggPCA(PCA, metadata, colnames(metadata))
dev.off()

# save data
save.image(paste0(output_workenvironment, "/QC_analysis.RData")
saveRDS(total, paste0(output_workenvironment,"/total.RDS")
