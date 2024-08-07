# WGCNA directed on pseudo bulk count of the samples normalised and sequencing batch regressed out
# according to the Normalization_entire_dataset_V2.R script

# By Joy Otten

# Library
library(Seurat)
library(ggplot2)
library(dplyr)
library(WGCNA)
library(QCEWAS)
library(ggrepel)
library(corrplot)
library(EDASeq)
library(dplyr)
library(QCEWAS)
library(edgeR)
library(bacon)
library(biomaRt)
library(msigdbr)
library(clusterProfiler)
library(stringr)
library(parallel)
library(xlsx)
library(openxlsx)
set.seed(1234)

# Data path
data_path <- "set your path"
output_path <- "set your path"
functions_path <- "set your path"
output_workenvironment <- "set your path"
results_path <- output_path # need to set this for some of the WGCNA functions, might change this

# Functions
source(paste0(functions_path, "/WGCNA_functions.R"))
source(paste0(functions_path, "/customPCA.R"))
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

# Data 
normalized_data <- readRDS(paste0(output_workenvironment, "/visium_normalised_entire_brain.RDS"))
metadata <- readRDS(paste0(data_path, "/visium_metadata_normalised_entire_brain.RDS"))
total <- readRDS(paste0(data_path, "/total.RDS"))

# Running a PCA analysis to check if there are significant co-variables
custom.PCA(beta = normalised_data, pd = metadata plot.title = "PCA")
# Here you can see that the sequencing depth is still not entirely even between all samples and more specifically between the clusters.
# This is due to that the number of spots are very variable between the clusters contributing to the sequencing depth difference.
# Additionally, the spatial slide number and region/cluster is significantly co-variable. Therefor, we remove these effects with empiricalBayesLM for WGCNA. 
# If we don't remove the region/cluster component, we end up with WGCNA modules related to the clusters. 

# Removal of significant co-variables with empiricalBayesLM
data <- as.data.frame(normalized_data)
rownames(metadata) <- metadata$samples

covariates_df <- data.frame(
  seq_depth = metadata$seq_depth,
  slide_nr = metadata$slide_nr,
  region = metadata$region)

data.lm = empiricalBayesLM(data, removedCovariates = covariates_df, retainedCovariates = metadata$group)
batch_corrected_data <- data.lm[["adjustedData"]]
batch_corrected_data <- as.data.frame(batch_corrected_data)
saveRDS(batch_corrected_data, paste0(output_workenvironment, "/batch_corrected_data.RDS"))

# PCA analysis after bayesLM
custom.PCA(beta = t(batch_corrected_data), pd = metadata, plot.title = "PCA after bayesLM")
# Here we see indeed that the significant co-variables found earlier are removed.

# WGCNA

## paths, you need to set this in some of the WGCNA functions the output are images in your allocated folder
image_path <- output_path
results_path <- output_path

# prepare for WGCNA
metadata$group <- relevel(as.factor(metadata$group), ref = "CTRL")
group <- relevel(as.factor(metadata$group), ref = "CTRL")
contrast <- c("group", "FC", "CTRL")

# Direction of the network of WGCNA is a signed network
data <- as.data.frame(batch_corrected_data)
direction = "signed"
data <- data_preprocessing(data)

# clustering samples to detect outliers
data <- clust_outliers(data, filter == FALSE)

# Set-up of the metadata for WGCNA
traitData = metadata[,c(3:6)]
rownames(traitData) <- rownames(metadata)
traitData$group <- as.numeric(traitData$group)
traitData$region <- as.numeric(traitData$region)
traitData$slide_nr <- as.numeric(traitData$slide_nr)
traitData$seq_depth <- as.numeric(traitData$seq_depth)

# Cluster samples with traitdata
cluster_samples(data, traitData)

## power calculation
# check the power_calc.pdf image in the output path to set the power
power <- power_calc(data) # should be a power of 4
power <- 4

# Automatic module detection via dynamic tree cutting
modules <- WGCNA_module(Expr = data, power = power, direction = "signed", split = 2) 
table(modules$colors) # 28 modules with module 0 containing 14591 genes
moduleLabels <- modules$colors
moduleColors = labels2colors(moduleLabels)
MEs <- modules$MEs;
geneTree = modules$dendrograms[[1]]

# To identify which genes are in which modules located
num <- unique(modules$colors)
modulecolours_labels <- as.data.frame(moduleColors)
modulecolours_labels$label <- moduleLabels
rownames(modulecolours_labels) <- names(moduleLabels)

genes_in_module <- list()
for(color in modules_interest){
  message(color)
  x <- which(modulecolours_labels$moduleColors == color)
  tmp <- modulecolours_labels[x,]
  genes_module <- rownames(tmp)
  genes_in_module <- append(genes_in_module, list(genes_module))
}
names(genes_in_module) <- modules_interest

# To write in excel which genes are located in which modules
for(module in names(genes_in_module)){
  message(module)
  module_data <- genes_in_module[[module]]
  if(module == "darkgrey"){
    n <- module
    wb <- createWorkbook()
    addWorksheet(wb, n)
    writeData(wb, n, module_data)
    saveWorkbook(wb, file = paste0(output_path, "/genes_in_module.xlsx")), overwrite = TRUE)
  } else {
    n <- module
    wb <- loadWorkbook(paste0(output_path, "/genes_in_module.xlsx"))
    addWorksheet(wb, n)
    writeData(wb, n, module_data, startRow = 1, startCol = 1, colNames = TRUE)
    saveWorkbook(wb, file = paste0(output_path, "/genes_in_module.xlsx"), overwrite = TRUE)
  }
}

## plot the dendrogram and the module colors underneath
# Metadata variable
names <- colnames(traitData)
colors <- list()
groups <- list()
for(i in names){
  x <- which(colnames(traitData) == i)
  x <- as.data.frame(traitData[,x])
  names(x) = i
  # next use this trait to define a gene significance variable
  x <- as.numeric(stats::cor(data, x, use = "p", method = c("pearson")))
  groups <- append(groups, list(x))
  color <- as.character(numbers2colors(x, signed = T))
  colors <- append(colors, list(color))
}
names(colors) <- names
names(groups) <- names

# plot the dendrogram and the module colors underneath 
## change each colors variable for every experiment. # also want to have the regions as an extra variable
blocknumber = 1
datColors = data.frame(moduleColors, colors[["region"]], colors[["slide_nr"]], colors[["group"]], colors[["seq_depth"]])[modules$blockGenes[[1]],]
pdf(paste0(output_path, "cluster_dendrogram_signed.pdf"))
plotDendroAndColors(modules$dendrograms[[1]], colors = datColors,
                    groupLabels = c("Module colors", "region","slide nr","group","seq depth"),
                    dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
dev.off()

## Calculation eigengenes
ME1 <- eigengene_network(Expr = data, module = modules, metadata = traitData)

## Relate modules to phenotype data
MEs0 <- module_trait(data, traitData, 87) 

# Gene relationship to traint and important modules: Gene significance and module membership
## calculate module membership values (aka. module eigengene based connectivity kME)
datKME = signedKME(data, MEs0, corOptions = "use = 'p', method = 'pearson'")

# groups of interest/metadata
group = colnames(traitData)
intramodular_analyses(kme = datKME, group = group, variables = groups)

# To identify which genes are in which modules located
num <- c(0:27)
modulecolours_labels <- as.data.frame(moduleColors)
modulecolours_labels$label <- moduleLabels
rownames(modulecolours_labels) <- names(moduleLabels)
# check for each experiment if the order of columns to numbers is the same so for example for grey the number is 0 and for turquoise it's 1 etc.
col <- as.list(c("grey", "turquoise", "blue", "brown", "yellow", "green", "red", "black", "pink",
                       "magenta", "purple", "greenyellow", "tan", "salmon", "cyan", "midnightblue", 
                       "lightcyan", "grey60", "lightgreen", "lightyellow", "royalblue", "darkred",
                       "darkgreen", "darkturquoise", "darkgrey", "orange", "darkorange", "white"))
                       
names(col) <- num

# Hub genes
hubgenes(Expr = data, moduleColors = moduleColors, power = power, direction = direction, MEs1 = ME1, output = output_path)
# The hubgenes script is adapted so that there is a correlation an p-value for all genes within the module. The one with the highest correlation value within the module is the hubgene

# intersect the hubgenes list with the genes in the module.
hubgenes_list <- read.csv(paste0(output_path, "list_hubgenes.csv"))
hubgenes <- list()
modules_interest <- c("darkgrey", "lightgreen")
for(module in modules_interest){
  message(module)
  intersect <- genes_in_module[[module]]
  x <- which(grepl(paste0("^ME", module), colnames(hubgenes_list)))
  hub_list <- hubgenes_list[,x]
  rownames(hub_list) <- hubgenes_list$X
  x <- which(rownames(hub_list) %in% intersect)
  hubs <- hub_list[x,]
  hubgenes <- append(hubgenes, list(hubs))
}
names(hubgenes) <- modules_interest

for(module in names(hubgenes)){
  data <- hubgenes[[module]]
  data$genes <- rownames(data)
  if(module == "darkgrey"){
    n <- module
    wb <- createWorkbook()
    addWorksheet(wb, n)
    writeData(wb, n, data)
    saveWorkbook(wb, file = paste0(output_path, "/hubgenes_WGCNA.xlsx"), overwrite = TRUE)
  }
  else {
    n <- module
    wb <- loadWorkbook(paste0(output_path, "/hubgenes_WGCNA.xlsx"))
    addWorksheet(wb, n)
    writeData(wb, n, data, startRow = 1, startCol = 1, colNames = TRUE)
    saveWorkbook(wb, file = paste0(output_path, "/hubgenes_WGCNA.xlsx"), overwrite = TRUE)
  }
}

# GO, Kegg, reactome and TF enrichment analysis
# I adapted the normal enrichment analysis function to include the transcription factor enrichment analysis

modules_interest <- c("darkgrey", "lightgreen", "white", "royalblue", "orange")
Label = as.data.frame(moduleLabels)
# Running the GSEA analysis
enriched <- enrichment_analysis(moduleLabels = moduleLabels, moduleColors = moduleColors, Expr = data, Labels1 = Label, modules_interest = modules_interest)

# Plotting enrichment analysis results
for(module in names(enriched)){
  message(module)
  module_data <- enriched[[module]]
  plot_data <- DataFrame()
  for(d in names(module_data)){
    if(d == "GO"){
      data <- module_data[[d]]
      if(nrow(data) == 0){
        message("no enriched GO terms")
      } else {
        x <- data[c(1:5), ]
        plot_data <- rbind(plot_data, x)
        plot_data <- as.data.frame(plot_data)
      } 
    }
    if(d == "KEGG"){
      data <- module_data[[d]]
      if(nrow(data) == 0){
        message("no enriched KEGG terms")
      } else {
        x <- data[c(1:5), ]
        plot_data <- rbind(plot_data, x)
        plot_data <- as.data.frame(plot_data)
      } 
    }
    if(d == "Reactome"){
      data <- module_data[[d]]
      if(nrow(data) == 0){
        message("no enriched Reactome terms")
      } else {
        x <- data[c(1:5), ]
        plot_data <- rbind(plot_data, x)
        plot_data <- as.data.frame(plot_data)
      } 
    }
    if(d == "TF"){
      data <- module_data[[d]]
      if(nrow(data) == 0){
        message("no enriched TF")
      } else {
        x <- data[c(1:5),]
        plot_data <- rbind(plot_data, x)
        plot_data <- as.data.frame(plot_data)
      }
    }
  }
  if(nrow(plot_data) >= 1){
    x <- which(is.na(plot_data) == TRUE)
    if(length(x) == 0){
      p <- ggplot(data=plot_data, aes(x = rownames(plot_data), y = Count, fill = p.adjust)) +
        geom_bar(stat = "identity") + theme_classic(base_size = 6) + coord_flip() + 
        labs(x = "", y = "Count") + guides(fill=guide_legend(title = "p.value"))
      pdf(paste0(output_path, "enrichment_module_", module, ".pdf"))
      plot(p)
      dev.off()
    } else {
      plot_data <- plot_data[-x,]
      p <- ggplot(data=plot_data,
                  aes(x = rownames(plot_data), y = Count, fill = p.adjust)) +
        geom_bar(stat = "identity") + theme_classic(base_size = 6) + coord_flip() + 
        labs(x = "", y = "Count") + guides(fill=guide_legend(title = "p.value"))
      pdf(paste0(output_path, "enrichment_module_", module, ".pdf"))
      plot(p)
      dev.off()
    }
  }
}
# Save the enrichment analysis results in excel files.
for(module in names(enriched)){
  message(module)
  module_data <- enriched[[module]]
  for(database in names(module_data)){
    message(database)
    data <- module_data[[database]]
    if(database == "GO"){
      n <- database
      wb <- createWorkbook()
      addWorksheet(wb, n)
      writeData(wb, n, data)
      saveWorkbook(wb, file = paste0(output_path, "/Enrichment_", module, ".xlsx"), overwrite = TRUE)
    }
    else {
      n <- database
      wb <- loadWorkbook(paste0(output_path, "/Enrichment_", module, ".xlsx"))
      addWorksheet(wb, n)
      writeData(wb, n, data, startRow = 1, startCol = 1, colNames = TRUE)
      saveWorkbook(wb, file = paste0(output_path, "/Enrichment_", module, ".xlsx"), overwrite = TRUE)
    }
  }
}

# Save workenvironment
save.image(paste0(workenvironment_path, "/all_regions_WGCNA.RData"))

