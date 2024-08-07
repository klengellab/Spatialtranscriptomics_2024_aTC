# Script that describes the functions for WGCNA with default settings

# By: Joy Otten


data_preprocessing <- function(Expr, verbose = FALSE){
  # Expr data is dataframe columns are the samples and rows the genes
  # This function filters on that the median variance is higher than 0,
  # if not this is filtered out. 
  
  # checks for genes and samples with too many missing values. 
  print(dim(Expr))
  gsg = goodSamplesGenes(Expr, verbose = 3)
  print(gsg$allOK)
  
  # removes genes with too many missing values
  if (!gsg$allOK){
    # Optionally, print the gene and sample names that were removed:
    if (sum(!gsg$goodGenes)>0)
      printFlush(paste("Removing genes:", paste(names(Expr)[!gsg$goodGenes], collapse = ", ")));
    if (sum(!gsg$goodSamples)>0)
      printFlush(paste("Removing samples:", paste(rownames(Expr)[!gsg$goodSamples], collapse = ", ")));
    # Remove the offending genes and samples from the data:
    Expr = Expr[gsg$goodSamples, gsg$goodGenes]
    dim(Expr)
  }
  return(Expr)
}

clust_outliers <- function(Expr, filter, verbose = FALSE){
  # filter is an object wherein you choose if you want to filter out the outliers
  # from the rest of the data or that you want to take them along.
  # plots the outliers in the samples
  sampleTree <- hclust(dist(Expr), method = "average")
  sizeGrWindow(12,9)
  par(cex = 0.6);
  par(mar = c(0,4,2,0))
  pdf(paste0(image_path, "sample_outliers.pdf"))
  plot(sampleTree, main = "Sample clustering to detect outliers", sub="",
       xlab = "", cex.lab = 1.5,
       cex.axis = 1.5, cex.main = 2)
  abline(h = 15, col = "red")
  dev.off()
  clust = cutreeStatic(sampleTree, cutHeight = 15, minSize = 10)
  print(table(clust))
  if(filter == TRUE){
    keepSamples = (clust==1)
    Expr <- Expr[keepSamples, ]
  }
  nGenes = ncol(Expr)
  nSamples = nrow(Expr)
  
  return(Expr)
}

power_calc <- function(Expr){
  # Expr: Dataframe with samples on the rows and genes on the columns
  # Calculates the power needed for WGCNA check always the plots
  # it could be that the power estimate gives another value back than 
  # that you think according to the plots.
  powers = c(c(1:10), seq(from = 12, to = 20, by =2))
  # Call the network topology analysis function
  sft = pickSoftThreshold(Expr, powerVector = powers, verbose = 5, cor =)
  #Plot the results
  sizeGrWindow(9, 5)
  par(mfrow = c(1,2));
  cex1 = 0.9
  pdf(paste0(image_path, "power_calc.pdf"))
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",
       ylab="Scale Free Topology Model Fit,signed R^2", type="n",
       main = paste("Scale independence"));
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red");
  # this line corresponds to using an R^2 cut-off of h
  abline(h=0.90,col="red")
  dev.off()
  pdf(paste0(image_path, "mean_connectivity.pdf"))
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab ="Soft Threshold (power)", ylab = "Mean connectivity", type = "n",
       main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
  dev.off()
  
  print(sft$fitIndices)
  power <- sft$powerEstimate
  
  return(power)
}

WGCNA_module <- function(Expr, power, direction, split){
  # Plotting the modules with pearson correlation
  # Expr is dataframe with the samples on the rows and genes on the columns
  # Power is the power you obtained after the power_calc function
  # Direction is the direction of network you want either directed or undirected in character
  # Split is a numerical value of 0, 2, 4
  cor <- WGCNA::cor
  enableWGCNAThreads(8)
  modules = blockwiseModules(Expr, power = power, cortype = "pearson", 
                             networkType = direction,
                             TOMType = direction, minModuleSize = 30,
                             maxBlockSize = 30000, deepSplit = split,
                             reassignThreshold = 0, mergeCutHeight = 0.25,
                             numericLabels = TRUE, pamRespectsDendro = FALSE,
                             verbose = 3)
  disableWGCNAThreads()
  cor <- stats::cor
  
  # plotting the figures
  plot_modules <- function(module){
    mergedColors = labels2colors(module$colors)
    pdf(paste0(image_path, "modules.pdf"))
    plotDendroAndColors(module$dendrogram[[1]],
                        mergedColors[module$blockGenes[[1]]],
                        "Module colors",
                        dendroLabels = FALSE, hang = 0.03,
                        addGuide = TRUE, guideHang = 0.05)
    dev.off()
  }
  
  return(modules)
}

cluster_samples <- function(Expr, datTraits, verbose = FALSE){
  # Expr is dataframe with rows containing the samples and columns the genes
  # DatTraits is a dataframe containing the metadata only numerical values
  sampleTree2 = hclust(dist(Expr), method = "average")
  traitColors = numbers2colors(datTraits, signed = FALSE);
  # Plot the sample dendrogram and the colors underneath.
  pdf(paste0(image_path, "sampledendro.pdf"))
  plotDendroAndColors(sampleTree2, traitColors,
                      groupLabels = names(datTraits),
                      main = "Sample dendrogram and trait heatmap")
  dev.off()
}

eigengene_network <- function(Expr, module, metadata){
  # Expr is the dataframe containing samples on the rows and genes on the columns
  # metadata: is a dataframe containing the phenotype/metadata with only numerical values
  # module: is the output from the WGCNA_module function
  moduleLabels <- module$colors
  moduleColors = labels2colors(moduleLabels)
  # plot eigen gene network
  ME1 <- moduleEigengenes(Expr, moduleColors)$eigengenes
  MET = orderMEs(cbind(ME1, metadata))
  pdf(paste0(image_path, "eigengene_network.pdf"))
  plotEigengeneNetworks(MET, "", marDendro = c(0.5,2,0.5,2), marHeatmap = c(4,6,2,6),
                        cex.lab = 0.4, xLabelsAngle=90)
  dev.off()
  
  return(ME1)
}

module_trait <- function(Expr, metadata, num_comparisons){
  # Expr is the dataframe with samples on the rows and genes on the columns
  # Metadata is the traitData dataframe with the phenotype containing only numerical values
  # num_comparisons, the number of comparisons that you eventually take to perform the BH correction
  # Visualization of modules with the phenotype data
  nGenes = ncol(Expr)
  nSamples = nrow(Expr)
  MEs0 <- moduleEigengenes(Expr, moduleColors)$eigengenes
  MEs0 <- orderMEs(MEs0)
  modTraitCor = stats::cor(MEs0, metadata, use = "p", method = "pearson")
  modTraitP = corPvalueStudent(modTraitCor, nSamples)
  modTraitP <- p.adjust(modTraitP, method = "BH", n = num_comparisons)
  textMatrix = paste(signif(modTraitCor, 2), "\n(",
                     signif(modTraitP, 1), ")", se = "")
  dim(textMatrix) = dim(modTraitCor)
  pdf(paste0(image_path, "module_trait_relationship_raw_signed_corrected_p.pdf"), width = 12, height = 10)
  par(mar = c(4,8,2,2))
  labeledHeatmap(Matrix = modTraitCor, xLabels = names(metadata), 
                 yLabels = names(MEs0), ySymbols = names(MEs0),
                 colorLabels = F, colors = blueWhiteRed(50),
                 textMatrix = textMatrix, setStdMargins = F, cex.text = 0.3, cex.lab = 0.3,
                 zlim = c(-1,1), main = paste("Module-trait relationships"))
  dev.off()
  
  return(MEs0)
}

hubgenes <- function(Expr, moduleColors, power, direction, MEs1, ensembl = FALSE, output){
  # This is the code from WGCNA adapted so that it put out the entire list of how well a
  # gene correlates to the other genes within the module. The gene with the highest correlation
  # is the hubgene
  # output is a path were the results will be saved
  # ensembl: If there are ensembl gene id's
  # Expr: Dataframe containing samples on rows and genes on columns
  # power: is the power calculated from the power_calc function
  # direction: Is the direction of the network either signed or unsigned
  # MEs1: Is the output from the eigengene_network function
  hub_genes <- as.data.frame(chooseTopHubInEachModule(Expr, moduleColors, omitColors = "grey", 
                                                      power = power, type = direction))
  names(hub_genes) <- "genes"
  
  list_hubgenes <- t(as.data.frame(mclapply(colnames(Expr), function(x){
    y <- numeric()
    for(i in colnames(MEs1)){
      y <- append(y, unlist(cor.test(Expr[, x], MEs1[, i], method = "pearson")[4:3]))
    }
    return(y)
  })))
  rownames(list_hubgenes) <- colnames(Expr)
  colnames(list_hubgenes) <- paste0(rep(colnames(MEs1), each = 2), rep(c(".r", ".p"), ncol(MEs1)))
  
  if(ensembl == TRUE){
    # convert Ensembl ID to NCBI gene ID
    ensembl <- useMart("ensembl")
    ensembl = useDataset("mmusculus_gene_ensembl", mart = ensembl)
    attributes = listAttributes(ensembl)
    
    x <- getBM(attributes=c("ensembl_gene_id", "external_gene_name"),
               filters = "ensembl_gene_id",
               values = hub_genes,
               mart = ensembl)
    y <- match(x$ensembl_gene_id, hub_genes$genes)
    x <- x[y,]
    hub_genes$gene_names <- x$external_gene_name
    
    print(hub_genes)
    write.csv(hub_genes, paste0(output, "hubgenes.csv"))
    
    x <- getBM(attributes=c("ensembl_gene_id", "external_gene_name"),
               filters = "ensembl_gene_id",
               values = rownames(list_hubgenes),
               mart = ensembl)
    
    y <- match(x$ensembl_gene_id, rownames(list_hubgenes))
    message("It could be that in your list of hubgenes are less genes compared to
          the expr data. This is due to that those ensembl names don't have a gene
          name")
    list_hubgenes <- list_hubgenes[y,]
    rownames(list_hubgenes) <- x$external_gene_name
  }
  
  write.csv(list_hubgenes, paste0(output, "list_hubgenes.csv"))
}

enrichment_analysis <- function(moduleLabels, moduleColors, Expr, Labels1, output, universe, modules_interest){
  # Performing enrichment analysis with GO, Reactome, KEGG and Transcription factor databases
  # Expr: Dataframe containing the samples on the rows and genes on columns
  # output is the path on where you want to have the output
  # Universe: is a character string of genes
  # Labels1: is the output from modules$colors in dataframe format
  # modules_interest: character string with the names of the modules that you are interested in to look further
  # moduleColors: is the output from labels2colors
  
  enrichment <- list()
  # keytype: character to identify the gene id's such as ENSEMBL, external_gene_name etc
  Labels1$moduleColors <- moduleColors
  Labels1$moduleColors = moduleColors
  Labels1$genes <- colnames(Expr)
  select_genes <- function(x){
    require(dplyr)
    y <- dplyr::filter(Labels1, moduleColors == x)
    genes <- y$genes
    return(as.data.frame(genes))
  }
  
  # Getting the genes per module
  Colours <- as.list(unique(Labels1$moduleColors))
  names(Colours) <- unique(Labels1$moduleColors)
  list_modules <- lapply(Colours, select_genes)
  for(i in modules_interest){
    enriched <- list()
    message("GO analysis")
    GO_datasets <- msigdbr(species = "Mus musculus", category = "C5")
    # Select for only the BP, CC and MF category in GO analysis
    filtered_GO_datasets <- GO_datasets[GO_datasets$gs_subcat %in% c("GO:CC", "GO:BP", "GO:MF"), ]
    
    x <- Colours[[i]]
    genes <- list_modules[[x]][["genes"]]
    enriched_GO <- as.data.frame(enricher(genes, universe = colnames(Expr), pvalueCutoff = 0.1,
                                          qvalueCutoff = 0.1,pAdjustMethod = "BH",TERM2GENE = filtered_GO_datasets[, c("gs_name", "gene_symbol")],
                                          TERM2NAME = filtered_GO_datasets[, c("gs_name", "gs_exact_source")]
    ))
    enriched <- append(enriched, list(enriched_GO))
    message("KEGG analysis")
    KEGG_datasets <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:KEGG")
    enriched_KEGG <- as.data.frame(enricher(genes, universe = colnames(Expr), pvalueCutoff = 0.1, 
                                            qvalueCutoff = 0.1, pAdjustMethod = "BH", TERM2GENE =
                                              KEGG_datasets[,c("gs_name", "gene_symbol")]))
    enriched_KEGG$Description <- str_sub(enriched_KEGG$ID,6) 
    enriched_KEGG$Description <- str_replace_all(enriched_KEGG$Description, "_", " ")
    
    enriched <- append(enriched, list(enriched_KEGG))
    # Reactome analysis
    message("Reactome analysis")
    reactome_datasets <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:REACTOME")
    enriched_reactome <- as.data.frame(enricher(genes, universe = colnames(Expr), pvalueCutoff = 0.1, 
                                                qvalueCutoff = 0.1, pAdjustMethod = "BH", TERM2GENE =
                                                  reactome_datasets[,c("gs_name", "gene_symbol")]))
    enriched_reactome$Description <- str_sub(enriched_reactome$Description, 10)
    enriched_reactome$ID <- str_sub(enriched_reactome$ID, start = 1, end = 8)
    
    enriched <- append(enriched, list(enriched_reactome))
    
    # Transcription factor enrichment analysis
    message("TF analysis")
    TF_datasets <- msigdbr(species = "Mus musculus", category = "C3")
    x <- which(TF_datasets$gs_subcat %in% c("TFT:TFT_Legacy", "TFT:GTRD"))
    TF_datasets <- TF_datasets[x,]
    enriched_tf <- as.data.frame(enricher(genes, universe = colnames(Expr), pvalueCutoff = 0.05, 
                                          qvalueCutoff = 0.1, pAdjustMethod = "BH", TERM2GENE =
                                            TF_datasets[,c("gs_name", "gene_symbol")]))
    enriched_tf$Description <- str_sub(enriched_tf$Description, 10)
    enriched_tf$ID <- str_sub(enriched_tf$ID, start = 1, end = 8)
    
    enriched <- append(enriched, list(enriched_tf))
    names(enriched) <- c("GO", "KEGG", "Reactome", "TF")
    enrichment <- append(enrichment, list(enriched))
  }
  names(enrichment) <- as.character(modules_interest)
  return(enrichment)
}
