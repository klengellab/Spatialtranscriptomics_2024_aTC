# Reclustering the first clustering
# Joy Otten

# Libraries:
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
library(clusterProfiler)
set.seed(1234)

# Data path
data_path <- "set your path"
output_path <- "set your path"
functions_path <- "set your path"
output_workenvironment <- "set your path"

# Data
total <- readRDS(paste0(output_workenvironment, "/total.RDS"))

# Splitting the data so that it will result in the regions detected earlier
regions_total <- list()
res <- c("0.2")
for(i in res){
  message(i)
  resolution <- paste0("SCT_snn_res.", i)
  Idents(total) <- resolution
  regions <- SplitObject(total, split.by = resolution)
  regions_total <- append(regions_total, list(regions))
}
names(regions_total) <- res

## QC analysis per cluster
## Elbow plot
res_region <- list()
clusters <- as.character(unique(Idents(total)))
res_temp <- list()
total_res <- regions_total[["0.2"]]
for(i in clusters){
  message(i)
  region <- total_res[[i]]
  cluster <- FindNeighbors(region, reduction = "pca", dims = 1:30)
  pct <- cluster[["pca"]]@stdev / sum(cluster[["pca"]]@stdev) * 100
  cumu <- cumsum(pct)
  co1 <- which(cumu > 90 & pct < 5)[1]
  message(co1)
  co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
  message(co2)
  res_temp <- append(res_temp, list(co2))
}
names(res_temp) <- clusters
# PCA 40 exhibits cumulative perfect greater than 90% for all regions. The last
# PCA where the change of % of variation is more than 0.1% is PC 17 for all regions.

# Resolution 0.2
# I did this on purpose per sample since this takes long
region <- regions_total[["0.2"]]
cluster_results_0.2 <- list()
for(i in clusters){
  message(i)
  number <- region[[i]]
  # treat each cluster separate 
  cluster <- FindNeighbors(number, reduction = "pca", dims = 1:17)
  cluster <- FindClusters(cluster, algorithm = 4, verbose = FALSE, method = "igraph",
                          resolution = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1))
  clustered <- Seurat::RunUMAP(cluster, reduction = "pca", dims = 1:17)
  cluster_results_0.2 <- append(cluster_results_0.2, list(clustered))
}
names(cluster_results_0.2) <- clusters

# Clustree
# Assess the different resolution parameters with Clustree
# Reference: https://cran.r-project.org/web/packages/clustree/vignettes/clustree.html
# With clustree we can assess the clustering parameters and find out which clustering parameter
# is the best to use for further processing
for(i in clusters){
  message(i)
  region <- cluster_results_0.2[[i]]
  p1 <- clustree(region, prefix = 'SCT_snn_res.')
  pdf(paste0(output_path, "/clustree_cluster_", i, ".pdf"))
  plot(p1)
  dev.off()
}
# Adding clusters together does not work since you will receive the same clusters again with the ones that you added together.

# Visualization of the clusters with different resolution parameters
sample_name <- c("35.1", "35.2", "35.3", "35.4", "61.1", "61.2", "61.3", "61.4",
                 "69.1", "69.2", "69.3", "69.4", "79.1", "79.2", "79.3", "79.4")
resolution <- as.character(c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8))
all_samples <- list()
output <- output_path
for(i in clusters){
  message(i)
  region <- cluster_results_0.2[[i]]
  samples <- SplitObject(region, split.by = "orig.ident")
  for(x in resolution){
    message(x)
    name <- paste0("SCT_snn_res.", x)
    message(name)
    pdf(paste0(output, "spatialplots_", i,"labels_resolution_", x, "_.pdf"))
    for(j in sample_name){
      message(j)
      temp_name <- strsplit(j, '.', fixed = TRUE)
      name_i <- temp_name[[1]][1]
      panel_i <- temp_name[[1]][2]
      image_name_i <- paste("slide", name_i, "_panel_", panel_i, sep = "")
      samp_name_i <- paste("slide", name_i, "_panel", panel_i, sep = "")
      sample <- samples[[samp_name_i]]
      
      Idents(sample) <- sample@meta.data[[name]]
      p1 <- SpatialDimPlot(sample, images = image_name_i, crop = FALSE, pt.size.factor = 1, label = TRUE)
      plot(p1)
    }
    dev.off()
    all_samples <- append(all_samples, list(samples))
    names(all_samples) <- x
  }
}

# Differential expression analysis per cluster we used the following resolutions for the first clustering to obtain the subclusters
# Cluster 1: Res.0.3
# Cluster 2: Res 0.3
# Cluster 3: Res 0.2
# Cluster 4: Res 0.2
# Cluster 5: Res 0.1
# Cluster 6: -
# Cluster 7: 0.2
# Cluster 8: 0.3
# Cluster 9: 0.2
# Cluster 10: -
# Cluster 11: -

# DE genes per region compared to all other regions
# Initialize the list to store DE results for each cluster
DE_clusters <- list()
clusters <- names(cluster_results_0.2)

# Loop through each cluster
for (i in clusters) {
  message(i)
  region <- cluster_results_0.2[[i]]
  
  # Determine the identity for each cluster and apply specific conditions
  DE_result <- NULL
  ident_resolution <- NULL
  subset_idents <- NULL
  
  if (i == "1") {
    ident_resolution <- "SCT_snn_res.0.3"
    subset_idents <- 5
  } else if (i == "2" || i == "8") {
    ident_resolution <- "SCT_snn_res.0.3"
    subset_idents <- ifelse(i == "8", 6, NULL)
  } else if (i == "3" || i == "4" || i == "9" || i == "7") {
    ident_resolution <- "SCT_snn_res.0.2"
  } else if (i == "5") {
    ident_resolution <- "SCT_snn_res.0.1"
  } 
  
  # Skip specific clusters (6, 10, and 11)
  if (i %in% c("6", "10", "11")) {
    next
  }
  
  # Set the identity and subset the region if needed
  Idents(region) <- region[[ident_resolution]]
  if (!is.null(subset_idents)) {
    region <- subset(region, idents = subset_idents, invert = TRUE)
  }
  
  # Find markers for the region
  DE_result <- FindAllMarkers(region, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
  
  # Store the result in the DE_clusters list
  DE_clusters[[as.character(i)]] <- list(DE_result)
  names(DE_clusters[[as.character(i)]]) <- ident_resolution
}
# Assign final names to DE_clusters based on the clusters used
names(DE_clusters) <- as.character(c(3, 4, 2, 8, 5, 1, 9, 7))

# Write data tables of DE genes per cluster vs all other clusters
output <- output_path
clusters <- c("7","9","3","1","2","4","5", "8")
for(i in clusters){
  message(i)
  cluster <- as.data.frame(DE_clusters[[i]])
  write.csv(cluster, paste0(output, "DE_reclusterd_cl_vs_cl_", i, ".csv"))
}

# Visualization of the DE genes per cluster per resolution
clusters <- as.character(c(1,2,3,4,5,7,8,9))
for(i in clusters){
  message(paste0("cluster_", i))
  DE <- DE_clusters[[i]]
  for(x in length(DE)){
    region <- DE[[x]]
    top_region <- region %>% group_by(cluster) %>% top_n(n = 25, wt = avg_log2FC)
    z <- unique(top_region$cluster)
    cluster <- cluster_results_0.2[[i]]
    temp <- strsplit(names(DE), split = "_")
    for(n in temp){
      resolution <- as.character(paste0("SCT_snn_res.", n[2]))
      message(resolution)
      Idents(cluster) <- cluster@meta.data[[resolution]]
      pdf(paste0(output, "cluster_", i, "_resolution_", resolution, ".pdf"))
      for(y in z){
        message(y)
        q <- which(top_region$cluster == y)
        w <- top_region[q,]
        plot(DotPlot(cluster, features = w$gene) + RotatedAxis())
      }
      dev.off()
    }
  }
}

# add the subclusters to the metadata file
# Initialize sub_regions column with 'NA'
total@meta.data$sub_regions <- 'NA'

# Define clusters to process
clusters <- as.character(1:11)

# Process each cluster
for (i in clusters) {
  region <- cluster_results_0.2[[i]]
  sub_clust_col <- ifelse(i %in% c("5"), "SCT_snn_res.0.1", 
                          ifelse(i %in% c("3", "4", "7", "9"), "SCT_snn_res.0.2", "SCT_snn_res.0.3"))
  
  # Get unique sub-clusters
  sub_clust <- unique(region[[sub_clust_col]])
  
  for (x in sub_clust) {
    message(paste0("cluster_", i, "_sub_clust_", x))
    # Identify cells in the sub-cluster
    clust <- which(region[[sub_clust_col]] == x)
    y <- region@meta.data[clust, ]
    n <- match(rownames(y), rownames(total@meta.data))
    # Assign sub-region labels
    total@meta.data[n, "sub_regions"] <- paste0(i, "-", x)
  }
}

# Assign specific labels for clusters 6, 10, and 11 based on broad regions
total@meta.data[total@meta.data$broad_regions == "6", "sub_regions"] <- "6"
total@meta.data[total@meta.data$broad_regions == "10", "sub_regions"] <- "10"
total@meta.data[total@meta.data$broad_regions == "11", "sub_regions"] <- "11"

# Visualisation of all subclusters on the samples
## get pallette
unique(total@meta.data$sub_regions)
colourCount = 35
getPalette = colorRampPalette(brewer.pal(12, "Set3"))
fill = getPalette(colourCount)
my_cols <- c('1-1'='#faf4cf','1-2'='#aeadb3', '1-3' = '#50d2fa','1-4'='#ff9a36','1-6'='#28CECA',
             '2-1'='#B95FBB','2-2'='#2FF18B','2-3' = '#ffbbee','3-1' = '#25aff5','3-2 '='#7F7F7F' ,'3-3'='#AB6B51',
             '4-1' = '#D4D915','4-2'='#1FA195','4-3' = '#f37735','5-1'='#6d5c91','5-2'='#E6C122','6'='#eec9d2',
             '7-1' = '#916e0f','7-2' = '#ef94ff','7-3' = '#c99789','7-4'='#c5f3fa','7-5'='#A4DFF2','8-1'='#def9ff',
             '8-2' = '#60A6A1','8-3' = '#448FC0','8-4' = '#F68282','8-5' = '#AF9C6F','9-1' = '#ff3377',
             '9-2' = '#00FFFF','9-3'='#CCB1F1','9-4' = '#A7D78C','10'='#AC8F14','11'='#0000FF')

Idents(total) <- total@meta.data$sub_regions

SpatialDimPlot(total, images = "slide35_panel_2", cols = my_cols)
total <- RunUMAP(total, reduction = "pca", dims = 1:17)
pdf(paste0(output_path, "/umap_reclustering.pdf")
DimPlot(total, reduction = "umap", cols = my_cols)
dev.off()

# Visualization of all the sub regions on the slides. 
output <- output_path
pdf(paste0(output, "/spatialplots_","subclustering",".pdf"))
for(j in sample_name){
  message(j)
  temp_name <- strsplit(j, '.', fixed = TRUE)
  name_i <- temp_name[[1]][1]
  panel_i <- temp_name[[1]][2]
  image_name_i <- paste("slide", name_i, "_panel_", panel_i, sep = "")
  samp_name_i <- paste("slide", name_i, "_panel", panel_i, sep = "")
  
  Idents(total) <- total@meta.data$sub_regions
  p1 <- SpatialDimPlot(total, images = image_name_i, crop = TRUE, cols = my_cols)
  plot(p1)
}
dev.off()

# Visualization of the subclusters but only within the broader region
Idents(total) <- total@meta.data$broad_regions
cluster <- as.character(c(1:11))
for(i in cluster){
  message(i)
  broad_region <- subset(total, idents = c(i))
  Idents(broad_region) <- broad_region@meta.data$sub_regions
  pdf(paste0(output, "/spatialplots_cluster_",i,".pdf"))
  for(j in sample_name){
    message(j)
    temp_name <- strsplit(j, '.', fixed = TRUE)
    name_i <- temp_name[[1]][1]
    panel_i <- temp_name[[1]][2]
    image_name_i <- paste("slide", name_i, "_panel_", panel_i, sep = "")
    samp_name_i <- paste("slide", name_i, "_panel", panel_i, sep = "")
    
    p1 <- SpatialDimPlot(broad_region, images = image_name_i, crop = FALSE, pt.size.factor = 1.2)
    plot(p1)
  }
  dev.off()
}

Idents(total) <- total@meta.data$sub_regions
pdf(paste0(output, "/umap_plot_subclusters",".pdf"))
DimPlot(total, reduction = "umap", label = TRUE, cols = my_cols)
dev.off()
Idents(total) <- total@meta.data$orig.ident
pdf(paste0(output, "/umap_plot_samples",".pdf"))
DimPlot(total, reduction = "umap", label = FALSE)
dev.off()
Idents(total) <- total@meta.data$percent.mt
pdf(paste0(output, "/umap_plot_mito",".pdf"))
DimPlot(total, reduction = "umap", label = FALSE)
dev.off()

# DE genes on all subclusters with one another
Idents(total) <- total$sub_regions
# removal of region 1-5 and 8-6 since they were only detected in a subset of the samples on purpose we didn't
# renormalise the data since that will lead also to entirely different clusters. The clusters that we found still carry
# meaningfull biological information however, we deleted those regions for downstream analysis
total <- subset(total, idents = c("1-5", "8-6"), invert = TRUE)
DE <- FindAllMarkers(total, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(DE, paste0(output_path, "/subcluster_specific_markers_clusters.csv"))

# Heatmap
DE %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top10
Idents(total) <- factor(x = Idents(total), levels=sort(levels(total)))
# I adjusted the order of clusters in excel/

pdf(paste0(output_path, "/heatmap.pdf"))
DoHeatmap(total, assay = "SCT", features = top10$gene, group.colors =  my_cols, size = 0) + NoLegend ()
dev.off()

pdf(paste0(output_path, "heatmap_legend.pdf")
DoHeatmap(total, assay = "SCT", features = top10$gene, group.colors =  my_cols, size = 0)
dev.off()

# QC matrixes of the subclustering
## Calculate the number of genes per cluster over all samples
metadata <- total@meta.data
p1 <- metadata %>%
  ggplot(aes(sub_regions, nFeature_Spatial)) + 
  geom_violin(aes(fill = sub_regions)) + 
  ylab("number of genes") + xlab("clusters") +
  geom_vline(xintercept = 500, linetype = "dashed") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual(values=my_cols) + theme_classic()

pdf(paste0(output_path, "/total_number_genes.pdf"))
plot(p1)
dev.off()

## Calculate the number of genes per cluster per sample
clusters <- SplitObject(total, split.by = "sub_regions")                                         
for(clust in names(clusters)){
  message(clust)
  cluster <- clusters[[clust]]
  meta <- cluster@meta.data
  count_results <- meta %>%
    group_by(orig.ident) %>%
    summarise(sum(nFeature_SCT))
  write.csv(count_results, paste0(output_path, "/number_genes_per_sample_cluster_", clust, ".csv"))
}

metadata <- total@meta.data                                                            
n_spots <- FetchData(total,
                     vars = c("orig.ident", "sub_regions")) %>% 
  dplyr::count(orig.ident, sub_regions)
write.csv(n_spots, paste0(output_path, "/total_spots_per_cluster_per_sample.csv"))           

save.image(paste0(workenvironment_path, "/reclustered_small_regions.RData")
saveRDS(total, paste0(workenvironment_path, "/total.RDS")

