# Script to obtain the spatial coordinates for deconvolution analysis by CARD

# Libraries
library(Seurat)
library(Matrix)
set.seed(1234)

# Data path
data_path <- "set your path"
output_path <- "set your path"
functions_path <- "set your path"
output_workenvironment <- "set your path"


# Data
total <- readRDS("/data/klengellab/Visium/workenvironment/total_res.RDS")
counts <- GetAssayData(total, assay = "Spatial")

# To obtain the positions of the spots since all spatial coordinates are the same we can use one sample to extract all coordinates
slide35_panel1 <- read10xVisium(paste0(data_path, "/S3__6-3-Slide-3-35-panel-1/outs",
                                type = c("sparse"),
                                data = c("raw"))
pos <- SpatialExperiment::spatialCoords(slide35_panel1)
colnames(pos) <- c("x", "y")
all_pos <- pos

# Now I need to change the barcodes and remove the slide and panel specific code
sample_name <- c("35.1", "35.2", "35.3", "35.4", "61.1", "61.2", "61.3", "61.4", 
                 "69.1", "69.2" ,"69.3", "69.4", "79.1", "79.2", "79.3", "79.4")
pos_final <- DataFrame()
for(i in sample_name){
  message(i)
  temp_name <- strsplit(i, '.', fixed = TRUE)
  name_i <- temp_name[[1]][1]
  panel_i <- temp_name[[1]][2]
  name <- paste0("slide", name_i, "_panel", panel_i)
  x <- str_remove(rownames(counts[[name]]), paste0(i, "_"))
  y <- match(x, rownames(all_pos)) 
  pos <- as.data.frame(all_pos[y,])
  pos$slice <- as.character(paste0(i))
  rownames(pos) <- rownames(counts[[name]])
  pos_final <- rbind(pos_final, pos)
}
pos_final <- as.data.frame(pos_final)

saveRDS(pos_final, paste0(output_workenvironment, "/positions_visium.RDS"))

# Preparation of the count data from the spatial data
# Removal of MT genes ince it's influencing the cell type predictions
x <- rownames(counts) %>% startsWith("mt-") 
y <- which(x == TRUE)
counts <- counts[-y,]
total <- subset(total, features = rownames(counts))
rm(counts)
all_samples <- SplitObject(total, split.by = "orig.ident")

counts <- list()
for(i in names(all_samples)){
  message(i)
  count <- all_samples[[i]]@assays[["Spatial"]]@counts
  counts <- append(counts, list(t(count)))
}
names(counts) <- names(all_samples)
saveRDS(counts, paste0(output_workenvironment, "/counts_st.RDS"))

