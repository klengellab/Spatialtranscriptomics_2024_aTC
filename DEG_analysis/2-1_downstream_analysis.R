# Downstream analysis of the DEGs analysis such as enrichment analysis and the correlation analysis
# By: Joy Otten


# Libraries:
library(Seurat)
library(ggplot2)
library(patchwork)
library(Matrix)
library(dplyr)
library(stringr)
library(tidyverse)
library(DESeq2)
library(leiden)
library(igraph)
library(pheatmap)
library(gridExtra)
library(grid)
library(tximport)
library(qqman)
library(QCEWAS)
library(org.Mm.eg.db)
library(car)
library(ensembldb)
library(sva)
library(bacon)
library(ggfortify)
library(ggpmisc)
library(readxl)
library(edgeR)
library(gridBase)
library(VennDiagram)
library(enrichR)
library(S4Vectors)
library(msigdbr)
library(clusterProfiler)
library(RRHO2)
library(RColorBrewer)
library(pheatmap)
library(scales)
library(reshape2)
library(circlize)
library(ggrepel)
library(openxlsx)
library(stringr)
set.seed(1234)

# Data path
data_path <- "set your path"
output_path <- "set your path"
functions_path <- "set your path"
output_workenvironment <- "set your path"

# Removal of two clusters that are not detected in all samples for downstream analysis
total <- readRDS(paste0(data_path, "/total.RDS"))
Idents(total) <- total@meta.data$sub_regions
total <- subset(total, idents = c("1-5", "8-6"), invert = TRUE) # deletion of 2 clusters not detected in all samples
total <- saveRDS(total, paste0(data_path, "/total.RDS"))

# Set up the environment
source(paste0(functions_path, "General_functions.R"))

# Data
sample_name <- c("35.1", "35.2", "35.3", "35.4", "61.1", "61.2", "61.3",
                 "61.4", "69.1", "69.2", "69.3", "69.4", "79.1", "79.2",
                 "79.3", "79.4")
test_methods <- c("ICA")
all_genes <- rownames(total)
layers <- unique(total$sub_regions) 

# Split data samples
samples <- SplitObject(total, split.by = "orig.ident")

# pseudo bulk counts, normalisation and differential expression analaysis performed in the 2-0_DEG_analysis.R script
# Read in the normalised counts and the limma results
ICA_limma <- readRDS(paste0(workenvironment_path, "/new_limmaresults_regions_raw.RDS"))
norm_counts <- readRDS(paste0(workenvironment_path,"/new_count_sub_regions_raw.RDS"))

# obtain DE_genes
DE_results <- list()
layers <- str_replace(layers, "-", ".")
for(i in layers){
  x <- which(ICA_limma[[i]][["res"]]$baconFDR <= 0.05)
  y <- ICA_limma[[i]][["res"]][x,]
  DE_results <- append(DE_results, list(y))
}
names(DE_results) <- layers

# Write output DE analysis 
for(i in layers){
  message(i)
  if(i == layers[1]){
    result <- as.data.frame(DE_results[[i]])
    write.xlsx2(result, paste0(output_path, "/DEGs_sub_regions.xlsx"), sheetName = i,
                append=FALSE, showNA = FALSE)
  } else {
    result <- as.data.frame(DE_results[[i]])
    write.xlsx2(result, paste0(output_path, "/DEGs_sub_regions.xlsx"), sheetName = i,
                append=TRUE, showNA = FALSE)
  }
}

all_length <- list()
for(i in names(DE_results)){
  x <- length(rownames(DE_results[[i]]))
  all_length <- append(all_length, x)
}
all_length <- as.numeric(all_length)
sum(all_length) # Total number of differential expressed genes across all brain regions


# Calculate ratio of DEGs wherein we devide the number of spots by the number of DEGs
# We have for each region 16 samples and therefore we divide the total number of spots per cluster by
# the number of DEGs.
# Initialize an empty list to store ratios
all_ratio <- list()

# Get unique sub-region layers from metadata
layers <- unique(total@meta.data$sub_regions)

# Loop through each layer
for (i in layers) {
  message(i)  # Print the current layer
  y <- str_replace(i, "-", ".")  # Replace "-" with "." in the layer name
  
  # Get the number of spots and their indices for the current layer
  number_spots <- length(which(total@meta.data$sub_regions == i))
  x <- which(total@meta.data$sub_regions == i)
  
  # Get metadata for the spots in the current layer
  cluster <- total@meta.data[x,]
  samples <- length(unique(cluster$orig.ident))  # Count unique samples
  
  # Calculate the number of spots per sample
  num <- number_spots / samples
  
  # Get the number of DEGs for the current layer
  DEGs <- nrow(DE_results[[y]])
  
  # Calculate the ratio of DEGs per spot and scale it
  r <- DEGs / num
  ratio <- r * 10  # You can modify the scaling factor as needed
  print(ratio)
  
  # Append the ratio to the list
  all_ratio <- append(all_ratio, ratio)
}
# Set names for the ratios based on layer names
names(all_ratio) <- layers
# Save the results to CSV and RDS files
write.csv(all_ratio, paste0(output_path,"/ratio_DEGs_subregions.csv"))
saveRDS(DE_results, paste0(workenvironment_path, "/DEGs_subregions.RDS"))

# Volcano plots for the DE genes
# reference https://bioconductor.org/packages/release/bioc/vignettes/EnhancedVolcano/inst/doc/EnhancedVolcano.html
# Define output path for the PDF
output <- output_path
pdf(paste0(output_path, "/volcanoplots_DEGs_small_regions.pdf"))

# Loop through each layer and generate volcano plots
for (i in layers) {
  message(i)  # Print the current layer
  
  # Replace "-" with "." in the layer name
  a <- str_replace(i, "-", ".")
  
  # Extract the results dataframe for the current layer
  test <- as.data.frame(ICA_limma[[a]][["res"]])
  
  # Initialize differential expression status
  test$dif <- "NO"
  
  # Mark differentially expressed genes (DEGs) with nominal and adjusted significance
  test$dif[test$logFC > 0 & test$baconP < 0.05] <- "up_nominal"
  test$dif[test$logFC < 0 & test$baconP < 0.05] <- "down_nominal"
  test$dif[test$logFC > 0 & test$baconFDR < 0.05] <- "UP"
  test$dif[test$logFC < 0 & test$baconFDR < 0.05] <- "DOWN"
  
  # Define colors for each category
  mycolors <- c("#0272bd", "#41a1e1", "#AFADB3", "#d47d7d", "#B22222")
  names(mycolors) <- c("DOWN", "down_nominal", "NO", "up_nominal", "UP")
  
  # Set rownames as symbol column for labeling
  test$symbol <- rownames(test)
  
  # Initialize label column for significant DEGs
  test$delabel <- NA
  test$delabel[!test$dif %in% c("NO", "up_nominal", "down_nominal")] <- test$symbol[!test$dif %in% c("NO", "up_nominal", "down_nominal")]
  
  # Plot volcano plot
  p <- ggplot(data = test, aes(x = logFC, y = -log10(baconFDR), col = dif, label = delabel)) +
    geom_point() +
    theme_classic()
  
  # Customize plot with colors and labels
  p2 <- p + scale_color_manual(values = mycolors) + geom_text_repel(aes(label = delabel), max.overlaps = Inf)
  p3 <- p2 + labs(title = paste0("Cluster ", i)) + xlab("Log2FC") + ylab("-log10 FDR")
  
  # Print the plot
  plot(p3)
}
# Close the PDF device
dev.off()

# File output of the DEG results in an excel file
for(a in layers){
  message(a)
  i <- str_replace(a, "-", ".")
  if(a == layers[1]){
    result <- as.data.frame(DE_results[[i]])
    write.xlsx2(result, paste0(output_path, "/DEGs_subclusters.xlsx"), sheetName = i,
                append=FALSE, showNA = FALSE)
  } else {
    result <- as.data.frame(DE_results[[i]])
    write.xlsx2(result, paste0(output_path, "/DEGs_subclusters.xlsx"), sheetName = i,
                append=TRUE, showNA = FALSE)
  }
}

# Enrichment analysis with GSEA on the nominal significant genes

# Obtain the nominal significant genes
output <- output_path

DE_results_nom <- list()
for(i in layers){
  x <- which(ICA_limma[[i]][["res"]]$baconP <= 0.05)
  y <- ICA_limma[[i]][["res"]][x,]
  DE_results_nom <- append(DE_results_nom, list(y))
}
names(DE_results_nom) <- layers

all_genes <- rownames(total) # set the universe
all_results <- list()
for(i in layers){
  message(i)
  DE_region <- DE_results_nom[[i]]
  
  if(nrow(DE_region) == 0){
    message(paste0("no differential expressed genes detected in cluster ", i))
    results <- DataFrame()
  } else {
    DE <- rownames(DE_region)
    results <- enrichment_analysis(genes = DE, universe = all_genes)
  }
  all_results <- append(all_results, list(results))
}
names(all_results) <- layers

# Z score for bubble plot
library(GOplot)
for(i in layers){
  message(i)
  region <- all_results[[i]]
  term <- DE_results_nom[[i]]
  term$ID <- rownames(term)
  term <- term[,c(1:6,9)]
  for(database in names(region)){
    message(database)
    if(database == "GO"){
      data <- region[[database]]
      l <- length(rownames(data))
      l <- l >= 1
      if(l == TRUE){
        colnames(data) <- c("term", "ID", "GeneRatio" ,"BgRatio","pvalue", "adj_pval", "qvalue", "Genes", "Count")
        data$Genes <- gsub("/", ",", data$Genes, fixed = TRUE)
        data$Category <- "GO"
        output <- circle_dat(data, term)
        write.csv(output, paste0(output_path, "/bubble_plot_region_", i, ".csv"))
      } else {
        message("no terms")
      }
    }
  }
}

# Write the nominal enriched terms to excel sheets
for (c in layers) {
  data <- all_results[[c]]
  tmp_data <- list()
  tmp_name <- list()
  
  for (i in names(data)) {
    if (nrow(data[[i]]) == 0) {
      message("empty")
    } else {
      d <- data[[i]]
      tmp_name <- append(tmp_name, i)
      tmp_data <- append(tmp_data, list(d))
    }
  }
  
  tmp_name <- as.character(tmp_name)
  names(tmp_data) <- tmp_name
  
  if (length(tmp_data) == 0) {
    print("No enrichment")
  } else {
    for (databases in 1:length(tmp_data)) {
      name <- names(tmp_data)
      message(databases)
      
      if (databases == 1) {
        n <- name[[databases]]
        wb <- createWorkbook()
        addWorksheet(wb, n)
        writeData(wb, n, tmp_data[[n]])
        saveWorkbook(wb, file = paste0(output_path, "/Enrichment_analysis_region_nominal_", c, ".xlsx"), overwrite = TRUE)
      } else {
        n <- name[[databases]]
        wb <- loadWorkbook(paste0(output_path, "/Enrichment_analysis_region_nominal_", c, ".xlsx"))
        addWorksheet(wb, n)
        writeData(wb, n, tmp_data[[n]], startRow = 1, startCol = 1, colNames = TRUE)
        saveWorkbook(wb, file = paste0(output_path, "/Enrichment_analysis_region_nominal_", c, ".xlsx"), overwrite = TRUE)
      }
    }
  }
}
save.image(paste0(workenvironment_path, "/DE_reclustered.RData"))


#### Correlation analysis on freezing with linear model
# ----------
# Define the list of FC samples and load metadata
fc_samples <- c("X35.1", "X35.3", "X61.1", "X61.3", "X69.1", "X69.3", "X79.1", "X79.3")
metadata <- read.csv(paste0(data_path, "/metadata_Visium.csv"))
output <- output_path

# Filter metadata for the FC group
metadata_fc <- metadata[metadata$group == "FC", ]
rownames(metadata_fc) <- paste0("X", metadata_fc$Samples)

# Clean up layer names
layers <- str_replace(layers, "-", ".")

# Initialize a dataframe to store all gene correlations
all_gene_correlations <- data.frame()

# Loop through each layer
for (i in layers) {
  message(i)
  
  # Extract DEGs and normalized counts for the current layer
  DE <- DE_results[[i]]
  genes <- rownames(DE)
  ICA_bulk_count_sub <- as.data.frame(norm_counts[[i]])
  x <- which(colnames(ICA_bulk_count_sub) %in% fc_samples)
  ICA_bulk_count_FC <- ICA_bulk_count_sub[, x]
  
  # Loop through each gene in the DEGs
  pdf(paste0(output, "cluster_boxplots_mean_", i,".pdf"))
  for (gene in genes) {
    gene_count <- as.data.frame(ICA_bulk_count_FC[gene, ])
    
    # Filter metadata for the current gene's samples
    x <- which(rownames(metadata_fc) %in% colnames(gene_count))
    metadata_s <- metadata_fc[x, ]
    gene_count <- as.data.frame(cbind(log(t(gene_count)), metadata_s$CSUS5))
    gene_count$sample <- rownames(gene_count)
    colnames(gene_count) <- c("count", "freezing", "s", "slide")
    gene_count <- as.data.frame(gene_count %>% arrange(freezing))
    
    # Define group labels and plot boxplot
    group <- factor(rep(c("FC", "CTRL"), 8), levels = c("FC", "CTRL"))
    names(group) <- paste0("X", metadata_s$Samples)
    x <- which(names(group) %in% colnames(ICA_bulk_count_sub))
    group_t <- group[x]
    p2 <- gene_boxplot(ICA_bulk_count_sub, i, gene, group_t)
    
    # Pearson correlation test
    pearson <- cor.test(gene_count$count, gene_count$freezing, method = "pearson")
    pearson_cor <- pearson$estimate
    pearson_pval <- pearson$p.value
    padj_pearson <- p.adjust(pearson$p.value, method = "BH", n = length(genes))
    
    # Create the ggplot for Pearson correlation
    p1 <- ggplot(gene_count, aes(x = freezing, y = count)) +
      geom_point() +
      geom_smooth(method = 'lm', se = FALSE, color = 'grey') +
      geom_text(aes(label = s), nudge_x = 0.25) +
      labs(x = "% Freezing according to CSUS5", 
           y = "Log2(normalized count)", 
           title = paste(gene, i, sep = " "),
           subtitle = paste("Pearson correlation:", round(pearson_cor, 3), 
                            "p.adjust:", round(padj_pearson, 5))) +
      theme_classic()
    
    # If Pearson correlation p-value is significant, record the results
    if (pearson_pval < 0.05) {
      # Combine Pearson correlation results
      data_correlations <- data.frame(
        gene = gene,
        layer = i,
        pearson_cor = pearson_cor,
        pearson_pval = pearson_pval,
        pearson_padj = padj_pearson
      )
      
      # Append to the main dataframe
      all_gene_correlations <- rbind(all_gene_correlations, data_correlations)
      plot(wrap_plots(p1, p2))
    }
    dev.off()
  }
}
# Convert all_gene_correlations to a data frame and save
all_gene_correlations <- as.data.frame(all_gene_correlations)
# Write correlations to a csv file
write.csv(all_gene_correlations, paste0(output_path, "/correlations.csv"))





# Want to plot the average %freezing along gene expression
fc_samples <- c("X35.1", "X35.3", "X61.1", "X61.3", "X69.1", "X69.3", "X79.1", "X79.3")
metadata <- read_csv("/PHShome/je637/Visium/data/metadata_Visium.csv")
x <- which(metadata$group == "FC")
metadata_fc <- metadata[x,]
fc_data <- metadata_fc[,c(6:10)]
test <- apply(fc_data, 1, mean)
metadata_fc$freezing_mean <- test
regions_of_int <- c("9.2", "8.4", "8.5", "9.1", "9.3")
for(i in regions_of_int){
  message(i)
  DE <- DE_results[[i]]
  genes <- rownames(DE)
  ICA_bulk_count_sub <- as.data.frame(norm_counts[[i]])
  x <- which(colnames(ICA_bulk_count_sub) %in% fc_samples)
  ICA_bulk_count_FC <- ICA_bulk_count_sub[,x]
  
  pdf(paste0(output, "cluster_boxplots_mean_", i,".pdf"))
  for(gene in genes){
    gene_count <- ICA_bulk_count_FC[gene,]
    
    gene_count <- cbind(log(t(gene_count)), metadata_fc$freezing_mean)
    colnames(gene_count) <- c("counts","freezing")
    # p <- boxplot(layer_count[gene,:]~group)
    gene_count <- as.data.frame(gene_count)
    gene_count$slide <- metadata_fc$slide_nr
    gene_count <- gene_count %>% arrange(g)
    
    test.lm <- lm(gene_count$count ~ gene_count$freezing + gene_count$slide)
    sum.lm <- summary(test.lm)
    p.val <- test.lm$coefficients[2]
    p1 <- ggplot(gene_count, aes(x=g, y = c)) +
      geom_point() + geom_smooth(method = 'lm') +
      geom_text(label = rownames(gene_count), nudge_x = 0.25) +
      labs(x = "mean %Freezing", y = "Log2(normalized count)", title = paste(gene,i,sep = " "),
           subtitle = paste("R-squared ", round(sum.lm$r.squared, 3), "", "p.value ", round(p.val, 3))) +
      theme_classic()
    
    
    p2 <- gene_boxplot(ICA_bulk_count_sub, i, gene, group)
    
    plot(wrap_plots(p1, p2))
  }
  dev.off()
}
save.image(paste0(output_path, "/downstream_analysis_DEGs.RData"))

