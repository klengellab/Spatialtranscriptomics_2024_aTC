# Visium pre-processing

# libaries
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(Matrix)
library(tidyverse)
library(topGO)
library(knitr)
library(kableExtra)
library(magrittr)
set.seed(1234)

# Functions
set_metadata <- function(object, ident, metadata, slide_nr) {
  object@meta.data[["orig.ident"]] <- ident
  object@meta.data[["metadata"]] <- metadata
  object@meta.data[["slide_nr"]] <- slide_nr
  return(object)
}
# Define a function to update coordinates for a given slide panel
update_coordinates <- function(panel, count_data, total) {
  # Extract coordinates for the specified panel
  coordinates <- as.data.frame(total@images[[panel]]@coordinates)
  
  # Get all coordinates from count_data and filter for the current panel
  all_coordinates <- colnames(count_data)
  slide_panel <- grepl(panel, all_coordinates)
  selected_coordinates <- all_coordinates[slide_panel]
  
  # Convert selected coordinates to a data frame
  selected_df <- as.data.frame(selected_coordinates)
  
  # Match the selected coordinates with the row names in the coordinates data
  matched_indices <- match(selected_df$selected_coordinates, rownames(coordinates))
  updated_coordinates <- coordinates[matched_indices, ]
  
  # Update the coordinates in the total object
  total@images[[panel]]@coordinates <- updated_coordinates
  
  return(dim(updated_coordinates))
}
image_dir <- "/PHShome/je637/Visium/images/"

# Loading in the data
#-------
slide35_panel1 <- Load10X_Spatial(data.dir = "/data/klengellab/Visium/data/slide_35/200522_NS500668_0814_AHGHTHBGXF/S3__6-3-Slide-3-35-panel-1/outs",
                                  filename = "filtered_feature_bc_matrix.h5",
                                  assay = "Spatial",
                                  slice = "slide35_panel_1",
                                  filter.matrix = TRUE)

slide35_panel2 <- Load10X_Spatial(data.dir = "/data/klengellab/Visium/data/slide_35/200522_NS500668_0814_AHGHTHBGXF/S3__1-2-Slide-3-35-panel-2/outs",
                                  filename = "filtered_feature_bc_matrix.h5",
                                  assay = "Spatial",
                                  slice = "slide35_panel_2",
                                  filter.matrix = TRUE)

slide35_panel3 <- Load10X_Spatial(data.dir = "/data/klengellab/Visium/data/slide_35/200522_NS500668_0814_AHGHTHBGXF/S3__5-4-Slide-3-35-panel-3/outs",
                                  filename = "filtered_feature_bc_matrix.h5",
                                  assay = "Spatial",
                                  slice = "slide35_panel_3",
                                  filter.matrix = TRUE)

slide35_panel4 <- Load10X_Spatial(data.dir = "/data/klengellab/Visium/data/slide_35/200522_NS500668_0814_AHGHTHBGXF/S3__2-2-Slide-3-35-panel-4/outs",
                                  filename = "filtered_feature_bc_matrix.h5",
                                  assay = "Spatial",
                                  slice = "slide35_panel_4",
                                  filter.matrix = TRUE)

slide61_panel1 <- Load10X_Spatial(data.dir = "/data/klengellab/Visium/data/slide_61/200521_NS500668_0813_AHGHJLBGXF/S3__4-1-Slide-1-61-panel-1/outs",
                                  filename = "filtered_feature_bc_matrix.h5",
                                  assay = "Spatial",
                                  slice = "slide61_panel_1",
                                  filter.matrix = TRUE)

slide61_panel2 <- Load10X_Spatial(data.dir = "/data/klengellab/Visium/data/slide_61/200521_NS500668_0813_AHGHJLBGXF/S3__1-1-Slide-1-61-panel-2/outs",
                                  filename = "filtered_feature_bc_matrix.h5",
                                  assay = "Spatial",
                                  slice = "slide61_panel_2",
                                  filter.matrix = TRUE)

slide61_panel3 <- Load10X_Spatial(data.dir = "/data/klengellab/Visium/data/slide_61/200521_NS500668_0813_AHGHJLBGXF/S3__6-2-Slide-1-61-panel-3/outs",
                                  filename = "filtered_feature_bc_matrix.h5",
                                  assay = "Spatial",
                                  slice = "slide61_panel_3",
                                  filter.matrix = TRUE)

slide61_panel4 <- Load10X_Spatial(data.dir = "/data/klengellab/Visium/data/slide_61/200521_NS500668_0813_AHGHJLBGXF/S3__2-1-Slide-1-61-panel-4/outs",
                                  filename = "filtered_feature_bc_matrix.h5",
                                  assay = "Spatial",
                                  slice = "slide61_panel_4",
                                  filter.matrix = TRUE)

slide69_panel1 <- Load10X_Spatial(data.dir = "/data/klengellab/Visium/data/slide_69/200525_NS500668_0815_AHHCLVBGXF/S3__4-3-Slide-2-069-panel-1/outs",
                                  filename = "filtered_feature_bc_matrix.h5",
                                  assay = "Spatial",
                                  slice = "slide69_panel_1",
                                  filter.matrix = TRUE)

slide69_panel2 <- Load10X_Spatial(data.dir = "/data/klengellab/Visium/data/slide_69/200525_NS500668_0815_AHHCLVBGXF/S3__1-2-Slide-2-069-panel-2/outs",
                                  filename = "filtered_feature_bc_matrix.h5",
                                  assay = "Spatial",
                                  slice = "slide69_panel_2",
                                  filter.matrix = TRUE)

slide69_panel3 <- Load10X_Spatial(data.dir = "/data/klengellab/Visium/data/slide_69/200525_NS500668_0815_AHHCLVBGXF/S3__5-4-Slide-2-069-panel-3/outs",
                                  filename = "filtered_feature_bc_matrix.h5",
                                  assay = "Spatial",
                                  slice = "slide69_panel_3",
                                  filter.matrix = TRUE)

slide69_panel4 <- Load10X_Spatial(data.dir = "/data/klengellab/Visium/data/slide_69/200525_NS500668_0815_AHHCLVBGXF/S3__2-2-Slide-2-069-panel-4/outs",
                                  filename = "filtered_feature_bc_matrix.h5",
                                  assay = "Spatial",
                                  slice = "slide69_panel_4",
                                  filter.matrix = TRUE)

slide79_panel1 <- Load10X_Spatial(data.dir = "/data/klengellab/Visium/data/slide_79/200526_NS500668_0816_AHGJHCBGXF/S3__6-4-Slide-4-079-panel-1/outs",
                                  filename = "filtered_feature_bc_matrix.h5",
                                  assay = "Spatial",
                                  slice = "slide79_panel_1",
                                  filter.matrix = TRUE)

slide79_panel2 <- Load10X_Spatial(data.dir = "/data/klengellab/Visium/data/slide_79/200526_NS500668_0816_AHGJHCBGXF/S3__1-4-Slide-4-079-panel-2/outs",
                                  filename = "filtered_feature_bc_matrix.h5",
                                  assay = "Spatial",
                                  slice = "slide79_panel_2",
                                  filter.matrix = TRUE)

slide79_panel3 <- Load10X_Spatial(data.dir = "/data/klengellab/Visium/data/slide_79/200526_NS500668_0816_AHGJHCBGXF/S3__3-3-Slide-4-079-panel-3/outs",
                                  filename = "filtered_feature_bc_matrix.h5",
                                  assay = "Spatial",
                                  slice = "slide79_panel_3",
                                  filter.matrix = TRUE)

slide79_panel4 <- Load10X_Spatial(data.dir = "/data/klengellab/Visium/data/slide_79/200526_NS500668_0816_AHGJHCBGXF/S3__2-4-Slide-4-079-panel-4/outs",
                                  filename = "filtered_feature_bc_matrix.h5",
                                  assay = "Spatial",
                                  slice = "slide79_panel_4",
                                  filter.matrix = TRUE)
#-------
#metadata <- read.csv("/PHShome/je637/Visium/data/metadata_Visium.csv", header = TRUE, sep = ",")

slide35_panel1 <- set_metadata(slide35_panel1, "slide35_panel1", "FC", "1")
slide35_panel2 <- set_metadata(slide35_panel2, "slide35_panel2", "CTRL", "2")
slide35_panel3 <- set_metadata(slide35_panel3, "slide35_panel3", "FC", "3")
slide35_panel4 <- set_metadata(slide35_panel4, "slide35_panel4", "CTRL", "4")
slide61_panel1 <- set_metadata(slide61_panel1, "slide61_panel1", "FC", "1")
slide61_panel2 <- set_metadata(slide61_panel2, "slide61_panel2", "CTRL", "2")
slide61_panel3 <- set_metadata(slide61_panel3, "slide61_panel3", "FC", "3")
slide61_panel4 <- set_metadata(slide61_panel4, "slide61_panel4", "CTRL", "4")
slide69_panel1 <- set_metadata(slide69_panel1, "slide69_panel1", "FC", "1")
slide69_panel2 <- set_metadata(slide69_panel2, "slide69_panel2", "CTRL", "2")
slide69_panel3 <- set_metadata(slide69_panel3, "slide69_panel3", "FC", "3")
slide69_panel4 <- set_metadata(slide69_panel4, "slide69_panel4", "CTRL", "4")
slide79_panel1 <- set_metadata(slide79_panel1, "slide79_panel1", "FC", "1")
slide79_panel2 <- set_metadata(slide79_panel2, "slide79_panel2", "CTRL", "2")
slide79_panel3 <- set_metadata(slide79_panel3, "slide79_panel3", "FC", "3")
slide79_panel4 <- set_metadata(slide79_panel4, "slide79_panel4", "CTRL", "4")


# Merging all data together
## merging is based on the raw counts later on with SCTransform I will perform a normalization
total <- merge(x = slide35_panel1, y = c(slide35_panel2, slide35_panel3, slide35_panel4, slide61_panel1,
                                         slide61_panel2, slide61_panel3, slide61_panel4, slide69_panel1,
                                         slide69_panel2, slide69_panel3, slide69_panel4, slide79_panel1,
                                         slide79_panel2, slide79_panel3, slide79_panel4),
               add.cell.ids = c("35.1", "35.2", "35.3", "35.4", "61.1", "61.2", "61.3", "61.4", "69.1",
                                "69.2", "69.3", "69.4", "79.1", "79.2", "79.3", "79.4"))

# raw count data
count_data <- as.data.frame(total@assays[["Spatial"]]@counts) 

# How is the data looking before filtering
plot1 <- VlnPlot(total, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(total, images = "slide35_panel_2", features = "nCount_Spatial") + theme(legend.position = "right")
plot3 <- SpatialFeaturePlot(total, images = "slide35_panel_2", features = "nFeature_Spatial") + theme(legend.position = "right")

wrap_plots(plot1, plot2)

# filtering
## Each spot needs at least a gene count of 200
## Each gene needs at least a gene count of 10 and detected in more than 2 spots otherwise discarded
# Original dimensions of count_data: 31053 genes, 47062 spots
dim(count_data) # 31053 47062

# Filter out spots with fewer than 200 counts
keep <- apply(count_data, 2, function(x) any(x >= 200))
count_data <- count_data[, keep]

# Filter out genes with fewer than 10 total reads
keep <- rowSums(count_data)
keep1 <- keep >= 10
count_data <- count_data[keep1, ]

# Dimensions after filtering genes: 19427 genes, 46837 spots
dim(count_data) # 19427 46837

# Iterate through each gene to check for detection in spots
for (i in 1:nrow(count_data)) {
  message(i) # Print the current index
  gene <- count_data[i, ]
  detected_in_spots <- apply(gene, 2, function(x) any(x > 0))
  num_detected_spots <- sum(detected_in_spots)
  
  # If the gene is detected in 2 or fewer spots, mark it as NA
  if (num_detected_spots <= 2) {
    data[i, ] <- NA
  }
}

# Filter the spots with NA values out

count_data <- as.matrix(count_data)
count_data <- Matrix(count_data, sparse = TRUE)
total@assays[["Spatial"]]@counts <- count_data
total@assays[["Spatial"]]@data <- count_data

counts_per_spot <- Matrix::colSums(count_data)
counts_per_gene <- Matrix::rowSums(count_data)
genes_per_spot <- Matrix::colSums(count_data > 1)
cat("counts for non-zero genes: ", genes_per_spot[1:5])

hist(log10(counts_per_spot+1), main = 'counts per spot', col = 'wheat')
hist(log10(genes_per_spot+1), main='genes per spot', col='wheat')
plot(counts_per_spot, genes_per_spot, log='xy', col='wheat')
title('counts vs genes per spot')
plot(sort(genes_per_spot), xlab='cell', log='y', main='genes per spot (ordered)')

#methods for seurat object
utils::methods(class = 'Seurat')

# overview mt genes
total[["percent.mt"]] <- PercentageFeatureSet(object = total, pattern = "mt-")
VlnPlot(object = total, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), ncol = 3)

count_data <- as.data.frame(count_data)
rownames(total@assays[["Spatial"]]@meta.features) <- total@assays[["Spatial"]]@meta.features$`rownames(count_data)`
test <- total@assays[["Spatial"]]@meta.features
test <- test[,-c(1)]
total@assays[["Spatial"]]@meta.features <- test
x <- match(colnames(count_data),rownames(total@meta.data))
total@meta.data <- total@meta.data[x,]

# need to have the coordinates in the images correct with the spots in the filtered dataset
# since I filtered out some spots
#Update coordinates for each panel
panels <- c("slide35_panel_1", "slide35_panel_2", "slide35_panel_3", "slide35_panel_4",
            "slide61_panel_1", "slide61_panel_2", "slide61_panel_3", "slide61_panel_4",
            "slide69_panel_1", "slide69_panel_2", "slide69_panel_3", "slide69_panel_4",
            "slide79_panel_1", "slide79_panel_2", "slide79_panel_3", "slide79_panel_4")

# Iterate over all panels and update coordinates
for (panel in panels) {
  print(paste("Updating coordinates for", panel))
  print(update_coordinates(panel, count_data, total))
}
rm(list=setdiff(ls(), c("total")))
save.image("/data/humgen/klengellab/Visium/data/workenvironment/Visium-preprocessing.RData")

######
# only to visualise the different normalization methods. I did not use that for
# the saved file and further processing of the sample
######

# run normalization to store sctransform residuals for all genes
total <- SCTransform(total, assay = "Spatial", return.only.var.genes = FALSE, verbose = FALSE) %>% RunPCA(verbose = FALSE) %>% RunUMAP(dims = 1:30)
# run standard log normalization for comparison
total <- NormalizeData(total, verbose = FALSE, assay = "Spatial")

DimPlot(total, group.by = "slide_nr", label = FALSE)

# Computes the correlation of the log normalized data and sctransform residuals with
# the number of UMIs
total <- GroupCorrelation(total, group.assay = "Spatial", assay = "Spatial", slot = "data",
                          do.plot = FALSE)
total <- GroupCorrelation(total, group.assay = "Spatial", assay = "SCT", slot = "scale.data",
                          do.plot = FALSE)


p7 <- GroupCorrelationPlot(total, assay = "Spatial", cor = "nCount_Spatial_cor") +
  ggtitle("Log Normalization") +
  theme(plot.title = element_text(hjust = 0.5))
p8 = GroupCorrelationPlot(total, assay = "SCT", cor = "nCount_Spatial_cor") +
  ggtitle("SCTransform Normalization") +
  theme(plot.title = element_text(hjust = 0.5))

# Percentage mitochondria content
total[["percent.mt"]] <- PercentageFeatureSet(total, pattern = "mt-*")
plot4 <- VlnPlot(total, features = "percent.mt", pt.size = 0.1) + NoLegend()
plot5 <- SpatialFeaturePlot(total, images = "slide35_panel_1", features = "percent.mt") + theme(legend.position = "right")

# Check the mitochondrial genes in the spatial slides.
for(i in names(total@images)){
  message(i)
  p1 <- SpatialFeaturePlot(total, images = i, features = "percent.mt")
  p2 <- p1 + ggtitle(i)
  plot(p2)
}
