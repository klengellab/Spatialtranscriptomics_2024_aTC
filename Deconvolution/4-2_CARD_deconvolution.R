# Deconvolution of the spatial transcriptomics data with the dropviz snRNA-seq data making use of the CARD deconvolution tool
# ref: https://github.com/YMa-lab/CARD
# By: Joy Otten

# Set environment variable
Sys.setenv('R_MAX_VSIZE'=100000000000)

# Libraries:
library(Biobase)
library(SummarizedExperiment)
library("SingleCellExperiment")
library(CARD)
library(stringr)
library(rdist)
library(rBeta2009)
library(gtools)
library(scatterpie)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
set.seed(1234)

# Data path
data_path <- "set your path"
output_path <- "set your path"
functions_path <- "set your path"
output_workenvironment <- "set your path"

# functions:
CARD_decon <- function (CARD_object) 
{
  ct.select = CARD_object@info_parameters$ct.select
  ct.varname = CARD_object@info_parameters$ct.varname
  sample.varname = NULL
  cat(paste0("## create reference matrix from scRNASeq...\n"))
  sc_eset = CARD_object@sc_eset
  Basis_ref = createscRef(sc_eset, ct.select, ct.varname, 
                          sample.varname)
  Basis = Basis_ref$basis
  Basis = Basis[, colnames(Basis) %in% ct.select]
  Basis = Basis[, match(ct.select, colnames(Basis))]
  spatial_count = CARD_object@spatial_countMat
  commonGene = intersect(rownames(spatial_count), rownames(Basis))
  commonGene = commonGene[!(commonGene %in% commonGene[grep("mt-", 
                                                            commonGene)])]
  cat(paste0("## Select Informative Genes! ...\n"))
  common = selectInfo(Basis, sc_eset, commonGene, ct.select, 
                      ct.varname)
  Xinput = spatial_count
  rm(spatial_count)
  B = Basis
  rm(Basis)
  Xinput = Xinput[order(rownames(Xinput)), ]
  B = B[order(rownames(B)), ]
  B = B[rownames(B) %in% common, ]
  Xinput = Xinput[rownames(Xinput) %in% common, ]
  Xinput = Xinput[rowSums(Xinput) > 0, ]
  Xinput = Xinput[, colSums(Xinput) > 0]
  colsumvec = colSums(Xinput)
  Xinput_norm = sweep(Xinput, 2, colsumvec, "/")
  B = B[rownames(B) %in% rownames(Xinput_norm), ]
  B = B[match(rownames(Xinput_norm), rownames(B)), ]
  spatial_location = CARD_object@spatial_location
  spatial_location = spatial_location[rownames(spatial_location) %in% 
                                        colnames(Xinput_norm), ]
  spatial_location = spatial_location[match(colnames(Xinput_norm), 
                                            rownames(spatial_location)), ]
  norm_cords = spatial_location[, c("x", "y")]
  norm_cords$x = norm_cords$x - min(norm_cords$x)
  norm_cords$y = norm_cords$y - min(norm_cords$y)
  scaleFactor = max(norm_cords$x, norm_cords$y)
  norm_cords$x = norm_cords$x/scaleFactor
  norm_cords$y = norm_cords$y/scaleFactor
  ED <- rdist(as.matrix(norm_cords))
  cat(paste0("## Deconvolution Starts! ...\n"))
  set.seed(20200107)
  Vint1 = as.matrix(rdirichlet(ncol(Xinput_norm), rep(10, 
                                                      ncol(B))))
  colnames(Vint1) = colnames(B)
  rownames(Vint1) = colnames(Xinput_norm)
  b = rep(0, length(ct.select))
  isigma = 0.1
  epsilon = 1e-04
  phi = c(0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99)
  kernel_mat <- exp(-ED^2/(2 * isigma^2))
  diag(kernel_mat) <- 0
  rm(ED)
  rm(Xinput)
  rm(norm_cords)
  gc()
  mean_X = mean(Xinput_norm)
  mean_B = mean(B)
  Xinput_norm = Xinput_norm * 0.1/mean_X
  B = B * 0.1/mean_B
  gc()
  ResList = list()
  Obj = c()
  for (iphi in 1:length(phi)) {
    res = CARDref(XinputIn = as.matrix(Xinput_norm), UIn = as.matrix(B), 
                  WIn = as.matrix(kernel_mat), phiIn = phi[iphi], max_iterIn = 1000, 
                  epsilonIn = epsilon, initV = as.matrix(Vint1), initb = rep(0, 
                                                                             ncol(B)), initSigma_e2 = 0.1, initLambda = rep(10, 
                                                                                                                            length(ct.select)))
    rownames(res$V) = colnames(Xinput_norm)
    colnames(res$V) = colnames(B)
    ResList[[iphi]] = res
    Obj = c(Obj, res$Obj)
  }
  Optimal = which(Obj == max(Obj))
  Optimal = Optimal[length(Optimal)]
  OptimalPhi = phi[Optimal]
  OptimalRes = ResList[[Optimal]]
  cat(paste0("## Deconvolution Finish! ...\n"))
  CARD_object@info_parameters$phi = OptimalPhi
  CARD_object@Proportion_CARD = sweep(OptimalRes$V, 1, rowSums(OptimalRes$V), 
                                      "/")
  CARD_object@algorithm_matrix = list(B = B * mean_B/0.1, 
                                      Xinput_norm = Xinput_norm * mean_X/0.1, Res = OptimalRes)
  CARD_object@spatial_location = spatial_location
  return(CARD_object)
}
CARD.visualize.pie.joy <- function (proportion, spatial_location, colors = NULL) 
{
  res_CARD = as.data.frame(proportion)
  res_CARD = res_CARD[, mixedsort(colnames(res_CARD))]
  location = as.data.frame(spatial_location)
  if (sum(rownames(res_CARD) == rownames(location)) != nrow(res_CARD)) {
    stop("The rownames of proportion data does not match with the rownames of spatial location data")
  }
  colorCandidate = c("#1e77b4", "#ff7d0b", "#ceaaa3", "#2c9f2c", 
                              "#babc22", "#d52828", "#9267bc", "#8b544c", "#e277c1", 
                              "#d42728", "#adc6e8", "#97df89", "#fe9795", "#4381bd", 
                              "#f2941f", "#5aa43a", "#cc4d2e", "#9f83c8", "#91675a", 
                              "#da8ec8", "#929292", "#c3c237", "#b4e0ea", "#bacceb", 
                              "#f7c685", "#dcf0d0", "#f4a99f", "#c8bad8", "#F56867", 
                              "#FEB915", "#C798EE", "#59BE86", "#7495D3", "#D1D1D1", 
                              "#6D1A9C", "#15821E", "#3A84E6", "#997273", "#787878", 
                              "#DB4C6C", "#9E7A7A", "#554236", "#AF5F3C", "#93796C", 
                              "#F9BD3F", "#DAB370", "#877F6C", "#268785")
                              if (is.null(colors)) {
                                if (ncol(res_CARD) > length(colorCandidate)) {
                                  colors = colorRampPalette(colorCandidate)(ncol(res_CARD))
                                }
                                else {
                                  colors = colorCandidate[sample(1:length(colorCandidate), 
                                                                 ncol(res_CARD))]
                                }
                              }
  else {
    colors = colors
  }
  data = cbind(res_CARD, location)
  ct.select = colnames(res_CARD)
  p = suppressMessages(ggplot() + geom_scatterpie(aes(x = x, 
                                                      y = y, r = 50), data = data, cols = ct.select, color = NA) + 
                         coord_fixed(ratio = 1) + scale_fill_manual(values = colors) + 
                         theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"), 
                               panel.background = element_blank(), plot.background = element_blank(), 
                               panel.border = element_rect(colour = "grey89", fill = NA, 
                                                           size = 0.5), axis.text = element_blank(), axis.ticks = element_blank(), 
                               axis.title = element_blank(), legend.title = element_text(size = 1, 
                                                                                         face = "bold"), legend.text = element_text(size = 8), 
                               legend.key = element_rect(colour = "transparent", 
                                                         fill = "white"), legend.key.size = unit(0.45, 
                                                                                                 "cm"), strip.text = element_text(size = 8, 
                                                                                                                                  face = "bold"), legend.position = "bottom") + 
                         guides(fill = guide_legend(title = "Cell Type")))
  return(p)
}

# Data
## special transcriptomics data
st_data <- readRDS(paste0(output_workenvironment, "/total.RDS"))
all_pos <- readRDS(paste0(output_workenvironment, "/positions_visium.RDS"))
counts <- readRDS(paste0(output_workenvironment,"/counts_st.RDS"))

## scRNA-seq data
scRNA <- readRDS(paste0(output_workenvironment, "/expr_data_dropviz.RDS"))
metadata <- read.csv(paste0(data_path, "/metadata_dropviz_glut_gad.csv"))
rownames(metadata) <- metadata$cell_labels # Change barcodes

# Remove NA values from metadata
unique(metadata$cell_type)
x <- which(is.na(metadata$cell_type) == "TRUE")
View(metadata[x,])

# Change the naming of the ASTROCYTES label
x <- which(metadata$cell_type == "ASTROCYTES")
metadata[x,5] <- "ASTROCYTE"
all(metadata$cell_labels %in% colnames(scRNA))
all(colnames(scRNA) %in% metadata$cell_labels)

# Reorder the snRNA-seq data and the metadata
reorder <- match(metadata$cell_labels, colnames(scRNA))
test <- scRNA[,reorder]
all(colnames(test) %in% metadata$cell_labels)
all(metadata$cell_labels %in% colnames(test))
scRNA <- test

# Deconvolution over all the samples
rownames(metadata) <- metadata$cell_labels
sample_name <- c("35.1", "35.2", "35.3", "35.4", "61.1", "61.2", "61.3", "61.4", 
                 "69.1", "69.2" ,"69.3", "69.4", "79.1", "79.2", "79.3", "79.4")
deconvolution <- list()
for(i in sample_name){
  message(i)
  temp_name <- strsplit(i, '.', fixed = TRUE)
  name_i <- temp_name[[1]][1]
  panel_i <- temp_name[[1]][2]
  name <- paste0("slide", name_i, "_panel", panel_i)
  st_expr <- counts[[name]]
  st_expr <- t(st_expr)
  pos <- all_pos[[name]]
  pos <- pos[,c(1:2)]
  rownames(pos) <- paste0(i, "_", rownames(pos))
  all(rownames(pos) %in% st_expr@Dimnames[[2]])
  x <- which(rownames(pos) %in% st_expr@Dimnames[[2]])
  pos <- as.data.frame(pos[x,])
  
  CARD_obj = createCARDObject(sc_count = scRNA,
                              sc_meta = metadata,
                              spatial_count = st_expr,
                              spatial_location = pos,
                              ct.varname = "cell_type",
                              sample.varname = NULL,
                              ct.select = c("GABA","GLUT", "OLIGODENDROCYTE","ASTROCYTE","ENDOTHELIAL","MICROGLIA", "POLYDENDROCYTE", "MURAL", "FIBROBLAST"),
                              minCountGene = 100,
                              minCountSpot = 5)
  CARD_obj_decon = CARD_decon(CARD_object = CARD_obj)
  deconvolution <- append(deconvolution, list(CARD_obj_decon))
}
names(deconvolution) <- sample_name

saveRDS(deconvolution, paste0(output_workenvironment, "/deconvolution_samples.RDS"))

# Visualisation of the deconvolution
colors = c("#eec9d2","#aeadb3","#B95FBB","#0096FF", "#D4D915","#31C53F", "#f37735", "#f5054f", "#00FFFF")
pdf(paste0(output_path, "/deconvolution_samples.pdf"))
for(i in names(deconvolution)){
  message(i)
  CARD_obj_decon <- deconvolution[[i]]
  p1 <- CARD.visualize.pie.joy(proportion = CARD_obj_decon@Proportion_CARD,spatial_location = CARD_obj_decon@spatial_location, colors = colors)
  plot(p1)
}
dev.off()


