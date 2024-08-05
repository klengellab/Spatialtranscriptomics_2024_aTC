# Script for the configuration of the dropviz data for the deconvolution of the spatial data
# ref: http://dropviz.org/ or article: Molecular Diversity and Specializations among the Cells of the Adult Mouse Brain Saunders et al. 2018
# By: Joy Otten

# Libraries:
library(Biobase)
library(limma)
library(marray)
library(convert)
library(biomaRt)
library(MuSiC)
library(ggplot2)
library(nnls)
library(xbioc)
library(reshape)
library(cowplot)
library(DropSeq.util)
set.seed(1234)

# Data path
data_path <- "set your path"
output_path <- "set your path"
functions_path <- "set your path"
output_workenvironment <- "set your path"

# Read in data
dge.path <- data_path
outcomes.path <- output_path
datasets <- c("F_GRCm38.81.P60Cortex_noRep5_FRONTALonly.raw.dge.txt.gz", "F_GRCm38.81.P60Cortex_noRep5_POSTERIORonly.raw.dge.txt.gz",
              "F_GRCm38.81.P60EntoPeduncular.raw.dge.txt.gz", "F_GRCm38.81.P60GlobusPallidus.raw.dge.txt.gz","F_GRCm38.81.P60Hippocampus.raw.dge.txt.gz",
              "F_GRCm38.81.P60Striatum.raw.dge.txt.gz", "F_GRCm38.81.P60Thalamus.raw.dge.txt.gz","H_1stRound_CrossTissue_Astrocytes_9-13-17.raw.dge.txt.gz", 
              "H_1stRound_CrossTissue_Endothelial_5-3-17.raw.dge.txt.gz","H_1stRound_CrossTissue_FibroblastLike_5-3-17.raw.dge.txt.gz", 
              "H_1stRound_CrossTissue_Microglia_Macrophage_5-3-17.raw.dge.txt.gz","H_1stRound_CrossTissue_Mural_5-3-17.raw.dge.txt.gz", 
              "H_1stRound_CrossTissue_Oligodendrocytes_5-3-17.raw.dge.txt.gz","H_1stRound_CrossTissue_Polydendrocytes_5-3-17.raw.dge.txt.gz"
)
cluster_outcomes <- c("F_GRCm38.81.P60Cortex_noRep5_FRONTALonly.cell_cluster_outcomes.RDS", "F_GRCm38.81.P60Cortex_noRep5_POSTERIORonly.cell_cluster_outcomes.RDS",
                      "F_GRCm38.81.P60EntoPeduncular.cell_cluster_outcomes.RDS", "F_GRCm38.81.P60GlobusPallidus.cell_cluster_outcomes.RDS",
                      "F_GRCm38.81.P60Hippocampus.cell_cluster_outcomes.RDS", "F_GRCm38.81.P60Striatum.cell_cluster_outcomes.RDS",
                      "F_GRCm38.81.P60Thalamus.cell_cluster_outcomes.RDS","H_1stRound_CrossTissue_Astrocytes_9-13-17.cell_cluster_outcomes.RDS", 
                      "H_1stRound_CrossTissue_Endothelial_5-3-17.cell_cluster_outcomes.RDS", "H_1stRound_CrossTissue_FibroblastLike_5-3-17.cell_cluster_outcomes.RDS", 
                      "H_1stRound_CrossTissue_Microglia_Macrophage_5-3-17.cell_cluster_outcomes.RDS","H_1stRound_CrossTissue_Mural_5-3-17.cell_cluster_outcomes.RDS",
                      "H_1stRound_CrossTissue_Oligodendrocytes_5-3-17.cell_cluster_outcomes.RDS","H_1stRound_CrossTissue_Polydendrocytes_5-3-17.cell_cluster_outcomes.RDS")

# Read in annotations
scRNA_anno = readRDS(paste0(data_path, "/annotation.BrainCellAtlas_Saunders_version_2018.04.01.RDS"))

# Put data into a list containing the snRNA-seq data
data <- list()
for(i in datasets){
  path <- paste0(dge.path, i)
  dge <- loadSparseDge(path)
  dge <- as.matrix(dge)
  data <- append(data, list(dge))
}
names(data) <- c("frontal_cortex", "posterior_cortex", "entopeduncular", "globus_pallidus", "hippocampus", "striatum",
                 "thalamus", "astrocytes", "endothelial", "fibroblast", "microglia",
                 "mural", "oligodendrocyte", "polydendrocyte")

outcomes <- list()
for(i in cluster_outcomes){
  path <- paste0(outcomes.path, i)
  outcome <- readRDS(path)
  outcomes <- append(outcomes, list(outcome))
}
names(outcomes) <- names(data)

data_expr <- list()
for(i in names(data)){
  d <- data[[i]]
  o <- outcomes[[i]]
  colnames(d) <- rownames(o)
  data_expr <- append(data_expr, list(d))
}
names(data_expr) <- names(data)
rm(data)

all_outcomes <- list()
for(i in names(outcomes)){
  # Filtering out all filter steps by selecting for NA if there are no reasons to delete out these genes
  o <- outcomes[[i]]
  x <- which(is.na(o$reason) == TRUE)
  d <- o[x,]
  all_outcomes <- append(all_outcomes, list(d))
}
names(all_outcomes) <- names(outcomes)
rm(outcomes)

# filter out all cells with doublets, outliers in expression data set
all_expr_filtered <- list()
for(i in names(data_expr)){
  expr <- data_expr[[i]]
  outcome <- all_outcomes[[i]]
  all(rownames(outcome) %in% colnames(expr))
  all(rownames(outcome) == colnames(expr))
  reorder <- match(rownames(outcome), colnames(expr))
  Expr_filtered <- expr[,reorder]
  all_expr_filtered <- append(all_expr_filtered, list(Expr_filtered))
}
names(all_expr_filtered) <- names(data_expr)

# Create metadata
all_metadata <- list()
for (i in names(all_outcomes)) {
  outcome <- all_outcomes[[i]]
  outcome$cell_type <- NA
  outcome$cell_labels <- rownames(outcome)
  
  # Define cell types based on cluster and subcluster assignments
  if (i == "frontal_cortex") {
    outcome$cell_type <- ifelse(outcome$cluster %in% c("3", "4", "5", "6", "7"), "GLUT", outcome$cell_type)
    outcome$cell_type <- ifelse(outcome$cluster %in% c("1", "2"), "GABA", outcome$cell_type)
    outcome$cell_type <- ifelse(outcome$cluster == "8", "ASTROCYTE", outcome$cell_type)
    outcome$cell_type <- ifelse(outcome$cluster == "9", "OLIGODENDROCYTE", outcome$cell_type)
    outcome$cell_type <- ifelse(outcome$cluster == "10", "POLYDENDROCYTE", outcome$cell_type)
    outcome$cell_type <- ifelse(outcome$cluster %in% c("12", "14"), "ENDOTHELIAL", outcome$cell_type)
    outcome$cell_type <- ifelse(outcome$cluster == "13", "MURAL", outcome$cell_type)
    
    outcome$cell_type <- ifelse(outcome$subcluster %in% c("11-1", "11-4"), "MICROGLIA", outcome$cell_type)
    outcome$cell_type <- ifelse(outcome$subcluster == "11-2", "GLUT", outcome$cell_type)
    outcome$cell_type <- ifelse(outcome$subcluster == "11-3", "MACROPHAGE", outcome$cell_type)
  }
  else if (i == "posterior_cortex") {
    outcome$cell_type <- ifelse(outcome$cluster %in% c("1", "2", "3", "6", "7"), "GLUT", outcome$cell_type)
    outcome$cell_type <- ifelse(outcome$cluster %in% c("4", "5"), "GABA", outcome$cell_type)
    outcome$cell_type <- ifelse(outcome$cluster == "8", "ASTROCYTE", outcome$cell_type)
    outcome$cell_type <- ifelse(outcome$cluster == "9", "OLIGODENDROCYTE", outcome$cell_type)
    outcome$cell_type <- ifelse(outcome$cluster == "10", "POLYDENDROCYTE", outcome$cell_type)
    outcome$cell_type <- ifelse(outcome$cluster %in% c("12", "14"), "ENDOTHELIAL", outcome$cell_type)
    outcome$cell_type <- ifelse(outcome$cluster == "13", "MURAL", outcome$cell_type)
    
    outcome$cell_type <- ifelse(outcome$subcluster %in% c("11-1", "11-3", "11-4"), "MICROGLIA", outcome$cell_type)
    outcome$cell_type <- ifelse(outcome$subcluster == "11-2", "MACROPHAGE", outcome$cell_type)
  }
  else if (i == "entopeduncular") {
    outcome$cell_type <- ifelse(outcome$cluster %in% c("1", "8"), "OLIGODENDROCYTE", outcome$cell_type)
    outcome$cell_type <- ifelse(outcome$cluster == "2", "ENDOTHELIAL", outcome$cell_type)
    outcome$cell_type <- ifelse(outcome$cluster %in% c("5", "6"), "POLYDENDROCYTE", outcome$cell_type)
    outcome$cell_type <- ifelse(outcome$cluster == "7", "ASTROCYTE", outcome$cell_type)
    
    outcome$cell_type <- ifelse(outcome$subcluster %in% c("3-1", "3-3", "3-4", "3-5", "3-6"), "MURAL", outcome$cell_type)
    outcome$cell_type <- ifelse(outcome$subcluster == "3-2", "ENDOTHELIAL", outcome$cell_type)
    outcome$cell_type <- ifelse(outcome$subcluster %in% c("4-1", "4-6"), "GLUT", outcome$cell_type)
    outcome$cell_type <- ifelse(outcome$subcluster %in% c("4-2", "4-3", "4-4", "4-5", "4-7", "4-8", "4-9", "4-10", "4-11", "4-12"), "GABA", outcome$cell_type)
  }
  else if (i == "globus_pallidus") {
    outcome$cell_type <- ifelse(outcome$cluster %in% c("1", "3"), "GABA", outcome$cell_type)
    outcome$cell_type <- ifelse(outcome$cluster == "5", "ASTROCYTE", outcome$cell_type)
    outcome$cell_type <- ifelse(outcome$cluster == "6", "EPENDYMAL", outcome$cell_type)
    outcome$cell_type <- ifelse(outcome$cluster %in% c("7", "9"), "ENDOTHELIAL", outcome$cell_type)
    outcome$cell_type <- ifelse(outcome$cluster == "8", "MURAL", outcome$cell_type)
    outcome$cell_type <- ifelse(outcome$cluster == "10", "OLIGODENDROCYTE", outcome$cell_type)
    
    outcome$cell_type <- ifelse(outcome$subcluster %in% c("4-1", "4-2", "4-3", "4-4", "4-5", "4-7"), "POLYDENDROCYTE", outcome$cell_type)
    outcome$cell_type <- ifelse(outcome$subcluster == "4-6", "MITOTIC", outcome$cell_type)
    outcome$cell_type <- ifelse(outcome$subcluster == "11-1", "MICROGLIA", outcome$cell_type)
    outcome$cell_type <- ifelse(outcome$subcluster == "11-2", "MACROPHAGE", outcome$cell_type)
    outcome$cell_type <- ifelse(outcome$subcluster %in% c("2-12", "2-11", "2-2", "2-8", "2-9"), "GLUT", outcome$cell_type)
    outcome$cell_type <- ifelse(outcome$subcluster %in% c("2-1", "2-3", "2-4", "2-5", "2-6", "2-7", "2-10", "2-13", "2-14", "2-15", "2-16", "2-17", "2-18", "2-19", "2-20", "2-21", "2-22", "2-23", "2-24", "2-25"), "GABA", outcome$cell_type)
  }
  else if (i == "hippocampus") {
    outcome$cell_type <- ifelse(outcome$cluster == "1", "GABA", outcome$cell_type)
    outcome$cell_type <- ifelse(outcome$cluster %in% c("2", "3", "4", "5", "6", "14"), "GLUT", outcome$cell_type)
    outcome$cell_type <- ifelse(outcome$cluster == "7", "ASTROCYTE", outcome$cell_type)
    outcome$cell_type <- ifelse(outcome$cluster == "8", "OLIGODENDROCYTE", outcome$cell_type)
    outcome$cell_type <- ifelse(outcome$cluster == "9", "POLYDENDROCYTE", outcome$cell_type)
    outcome$cell_type <- ifelse(outcome$cluster == "11", "EPENDYMAL", outcome$cell_type)
    outcome$cell_type <- ifelse(outcome$cluster == "12", "CHOROID_PLEXUS", outcome$cell_type)
    outcome$cell_type <- ifelse(outcome$cluster == "13", "NEUROGENESIS", outcome$cell_type)
    outcome$cell_type <- ifelse(outcome$cluster %in% c("15", "17"), "ENDOTHELIAL", outcome$cell_type)
    outcome$cell_type <- ifelse(outcome$cluster == "16", "MURAL", outcome$cell_type)
    
    outcome$cell_type <- ifelse(outcome$subcluster == "10-1", "MICROGLIA", outcome$cell_type)
    outcome$cell_type <- ifelse(outcome$subcluster == "10-2", "MACROPHAGE", outcome$cell_type)
  }
  else if (i == "striatum") {
    outcome$cell_type <- ifelse(outcome$cluster %in% c("10", "11", "12", "14", "15"), "GABA", outcome$cell_type)
    outcome$cell_type <- ifelse(outcome$cluster == "1", "EPENDYMAL", outcome$cell_type)
    outcome$cell_type <- ifelse(outcome$cluster == "2", "NEUROGENESIS", outcome$cell_type)
    outcome$cell_type <- ifelse(outcome$cluster == "3", "OLIGODENDROCYTE", outcome$cell_type)
    outcome$cell_type <- ifelse(outcome$cluster == "4", "ASTROCYTE", outcome$cell_type)
    outcome$cell_type <- ifelse(outcome$cluster == "5", "POLYDENDROCYTE", outcome$cell_type)
    outcome$cell_type <- ifelse(outcome$cluster %in% c("7", "8"), "ENDOTHELIAL", outcome$cell_type)
    outcome$cell_type <- ifelse(outcome$cluster == "9", "MURAL", outcome$cell_type)
    outcome$cell_type <- ifelse(outcome$cluster == "6", "GLUT", outcome$cell_type)
    outcome$cell_type <- ifelse(outcome$cluster == "13", "IMMUNE", outcome$cell_type)
  }
  
  all_metadata[[i]] <- outcome
}


# Assign cell types based on cluster and subcluster

# Assign cell types by cluster
for (number in unique(outcome$cluster)) {
  if (number %in% c("1", "2")) {
    outcome[outcome$cluster == number, 4] <- "GLUT"
  } else if (number == "4") {
    outcome[outcome$cluster == number, 4] <- "TH"
  } else if (number == "7") {
    outcome[outcome$cluster == number, 4] <- "ASTROCYTE"
  } else if (number == "8") {
    outcome[outcome$cluster == number, 4] <- "EPENDYMAL"
  } else if (number %in% c("5", "6")) {
    outcome[outcome$cluster == number, 4] <- "POLYDENDROCYTE"
  } else if (number %in% c("12", "14")) {
    outcome[outcome$cluster == number, 4] <- "ENDOTHELIAL"
  } else if (number %in% c("10", "11")) {
    outcome[outcome$cluster == number, 4] <- "OLIGODENDROCYTE"
  } else if (number == "13") {
    outcome[outcome$cluster == number, 4] <- "MURAL"
  }
}

# Assign cell types by subcluster
for (number in unique(outcome$subcluster)) {
  if (number == "9-2") {
    outcome[outcome$subcluster == number, 4] <- "MICROGLIA"
  } else if (number == "9-1") {
    outcome[outcome$subcluster == number, 4] <- "MACROPHAGE"
  } else if (number %in% c("3-1", "3-2", "3-3", "3-4", "3-5", "3-8", "3-13", "3-14", "3-16", "3-17", "3-18", "3-19")) {
    outcome[outcome$subcluster == number, 4] <- "GABA"
  } else if (number %in% c("3-6", "3-7", "3-9", "3-10", "3-11", "3-12", "3-15")) {
    outcome[outcome$subcluster == number, 4] <- "GLUT"
  }
}

# Matching the metadata
test <- rbindlist(all_metadata)
x <- unique(test$cell_labels)
y <- match(x, test$cell_labels)
test <- test[y,]
metadata <- test

write.csv(metadata, paste0(data_path, "metadata_dropviz_glut_gad.csv"))

# Make a sparse matrix
matrixes <- list()
for(i in names(all_expr_filtered)){
  d <- all_expr_filtered[[i]]
  data <- Matrix(d, sparse = TRUE)
  matrixes <- append(matrixes, list(data))
}
names(matrixes) <- names(all_expr_filtered)
rm(all_expr_filtered)

# Merge the sparse matrices
data <- RowMergeSparseMatrices(data, matrixes)

# Match the metadata and data
reorder <- match(metadata$cell_labels, colnames(data))
data <- data[,reorder]
all(colnames(data) %in% metadata$cell_labels)
all(metadata$cell_labels %in% colnames(data))

saveRDS(data, paste0(output_workenvironment, "/expr_data_dropviz.RDS"))

