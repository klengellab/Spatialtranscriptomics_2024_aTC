# Transcriptional changes across brain regions according to the method developed by Michael Gandal
# Reference: https://github.com/dhglab/Broad-transcriptomic-dysregulation-across-the-cerebral-cortex-in-ASD

# Libraries:
library(cowplot)
library(tidyverse)
library(ggridges)
library(scales)
library(SingleCellExperiment)
library(Matrix)
library(knitr)
library(readxl)
library(writexl)
library(Seurat)
library(msigdbr)
library(clusterProfiler)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(enrichplot)
library(dplyr)
library(limma)
library(edgeR)
library(VennDiagram)

# Data path
data_path <- "set your path"
output_path <- "set your path"
functions_path <- "set your path"
output_workenvironment <- "set your path"


### Now, make a histogram for every comparison of the number of times a 'true' gene was in a scrambled permutation
load(paste0(output_workenvironment, "/TrueGenes_in_PermutedGenes_All_Compiled.RData")

pdf(file=paste0(output_path, "attenuation/Permutation-Count-Density.pdf"),width=10)
par(mfrow=c(1,2))
for(comp in names(true_diff_genes_fc_all)){
  if(length(true_diff_genes_count_fc_all[[comp]]) > 0){
    plot(density(as.numeric(true_diff_genes_count_fc_all[[comp]]),bw=100), main=paste(comp,"; FC",sep=""))
    abline(v=500,col="red")
  }else{
    plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
    text(x = 0.5, y = 0.5, paste(comp,"; \n No True Diff Genes; \n FC",sep=""), 
         cex = 1, col = "black")
  }
  if(length(true_diff_genes_count_ctrl_all[[comp]]) > 0){
    plot(density(as.numeric(true_diff_genes_count_ctrl_all[[comp]]),bw=100),main=paste(comp,"; CTRL",sep=""))
    abline(v=500,col="red")
  }else{
    plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
    text(x = 0.5, y = 0.5, paste(comp,"; \n No True Diff Genes; \n CTRL",sep=""), 
         cex = 1, col = "black")
  }
}
dev.off()

### Choose filter of 5% (highly conservative)
### Filer is for a gene to be kept as a 'true' cp gene, from perm test; needs to be present in less than 500 out of the 10,000 perms
true_diff_genes_final_fc=list()
true_diff_genes_final_ctrl=list()

for(comp in names(true_diff_genes_fc_all)){
  idx = which(true_diff_genes_count_fc_all[[comp]] < 500)
  true_diff_genes_final_fc[[comp]] = true_diff_genes_fc_all[[comp]][idx]
  idx = which(true_diff_genes_count_ctrl_all[[comp]] < 500)
  true_diff_genes_final_ctrl[[comp]] = true_diff_genes_ctrl_all[[comp]][idx]
}
save(true_diff_genes_final_fc,true_diff_genes_final_ctrl,file=paste0(output_workenvironment, "/All_ARI_Genes.RData"))

# Separate the lists for attenuated and over-patterned
overpatterned <- c("4.2_v_9.4", "4.2_v_7.5", "3.2_v_1.1", "1.1_v_3.3", "4.2_v_1.1",
                   "2.2_v_3.3", "3.1_v_1.1", "8.3_v_3.3")
attenuated <- c("10_v_8.5", "5.1_v_10", "6_v_10", "6_v_7.4", "10_v_7.3", "10_v_7.1")
total <- readRDS(paste0(output_workenvironment, "/total.RDS"))
total_genes <- rownames(total)
rm(total)

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

# Perform enrichment analysis per attenuated or over patterned region containing two brain regions
enrichment_attenuated <- list()
enrichment_overpatterned <- list()
for(i in names(true_diff_genes_final_ctrl)){
  message(i)
  if(i %in% overpatterned){
    data <- true_diff_genes_final_fc[[i]]
    tmp <- enrichment_analysis(data, total_genes)
    enrichment_overpatterned <- append(enrichment_overpatterned, list(tmp))
  }
  if(i %in% attenuated){
    data <- true_diff_genes_final_ctrl[[i]]
    tmp <- enrichment_analysis(data, total_genes)
    enrichment_attenuated <- append(enrichment_attenuated, list(tmp))
  }
}
names(enrichment_attenuated) <- attenuated
names(enrichment_overpatterned) <- overpatterned

for(i in names(enrichment_attenuated)){
  data <- enrichment_attenuated[[i]]
  for(database in names(data)){
    if(database == "GO"){
      wb = createWorkbook()
      n = database
      addWorksheet(wb, n)
      writeData(wb, n, data[[database]])
      saveWorkbook(wb, file = paste0(output_path, "/attenuated_regions_together_", i, ".xlsx"), overwrite = TRUE)
    }
    else {
      n = database
      addWorksheet(wb, n)
      writeData(wb, n, data[[database]])
      saveWorkbook(wb, file = paste0(output_path, "/attenuated_regions_together_", i, ".xlsx"), overwrite = TRUE)
    }
  }
}

for(i in names(enrichment_overpatterned)){
  data <- enrichment_overpatterned[[i]]
  for(database in names(data)){
    if(database == "GO"){
      wb = createWorkbook()
      n = database
      addWorksheet(wb, n)
      writeData(wb, n, data[[database]])
      saveWorkbook(wb, file = paste0(output_path, "/overpatterned_regions_together_", i, ".xlsx"), overwrite = TRUE)
    }
    else {
      n = database
      addWorksheet(wb, n)
      writeData(wb, n, data[[database]])
      saveWorkbook(wb, file = paste0(output_path, "/overpatterned_regions_together_", i, ".xlsx"), overwrite = TRUE)
    }
  }
}

load(paste0(output_workenvironment, "/True_TRI_DGE_FC_v_Control.RData")) # DGE results from TRI analysis with the true samples

### from this point on - just look at the regional comparisons significantly attenuated in FC
sig_comps = attenuated
names(sig_comps) = rep("CTRL", length(sig_comps))
sig_comps = gsub("-", "_v_", sig_comps)

sig_comp_list_ctrl <- true_diff_genes_final_ctrl[which(names(true_diff_genes_final_ctrl) %in% sig_comps[which(names(sig_comps)=="CTRL")])] 

### remove any genes from the controls that are present in FC
sig_comp_list_ctrl_unique = list()

for(comp in names(sig_comp_list_ctrl)){
  print(comp)
  rm = true_diff_genes_final_fc[[comp]][which(true_diff_genes_final_fc[[comp]] %in% sig_comp_list_ctrl[[comp]])]
  print(length(rm))
  if(length(rm)==0){
    sig_comp_list_ctrl_unique[[comp]]=sig_comp_list_ctrl[[comp]]
  }else{
    sig_comp_list_ctrl_unique[[comp]]=sig_comp_list_ctrl[[comp]][-match(rm,sig_comp_list_ctrl[[comp]])]
  }
}
sig_comp_list_ctrl = sig_comp_list_ctrl_unique

names(sig_comp_list_ctrl) <- paste(names(sig_comp_list_ctrl),"_CTRL-gt-FC",sep="")
sig_comp_list = sig_comp_list_ctrl


### determine which region each comparison is expressed in
datExpr.reg <- readRDS(paste0(output_workenvironment, "/visium_normalised_entire_brain.RDS"))
datMeta <- readRDS(paste0(data_path, "/visium_metadata_normalised_entire_brain.RDS"))
load(paste0(output_workenvironment, "/Permutation_Index.RData"))

# Here we are going to sort the genes based on higher/lower expression in region 1 or region 2
cp_genes = list()

for(comp in names(sig_comp_list)){
  
  fc_att_genes = sig_comp_list[[comp]]
  
  idx_name = gsub("_CTRL-gt-FC", "", comp)
  idx_name = gsub("_v_", "-",idx_name)
  
  reg1 = strsplit(comp, split="_")[[1]][1]
  reg2 = strsplit(comp,split="_")[[1]][3]
  
  idx_name = paste0(reg1, " - ", reg2)
  fc_att_expr = datExpr.reg[which(rownames(datExpr.reg) %in% fc_att_genes),idx_comp_all[[idx_name]]]
  datMeta_fc_att = datMeta[idx_comp_all[[idx_name]],]
  datMeta_fc_att$Dx_Reg = paste(datMeta_fc_att$group,datMeta_fc_att$region,sep="_")
  
  ## sort the genes which are more highly expressed in reg1 or reg2
  
  count=rep(NA,length(fc_att_genes))
  names(count) <- rownames(fc_att_expr)
  for(gene in rownames(fc_att_expr)){
    avg_1=median(fc_att_expr[gene,which(datMeta_fc_att$Dx_Reg==paste("CTRL",reg1,sep="_"))])
    avg_2=median(fc_att_expr[gene,which(datMeta_fc_att$Dx_Reg==paste("CTRL",reg2,sep="_"))])
    if(avg_1 > avg_2){
      count[gene]=reg1
    }else{
      count[gene]=reg2
    }
  }
  
  reg1_genes = names(count)[which(count==reg1)]
  reg2_genes = names(count)[which(count==reg2)]
  cp_genes[[comp]][[reg1]]=reg1_genes
  cp_genes[[comp]][[reg2]]=reg2_genes
  
}
save(cp_genes,file=paste0(output_workenvironment, "/ARI_genes_from_RegComparisons_Attenuated_inFC_RegSorted.RData")

# Additional enrichment analysis on the genes assigned to either region 1 or region 2 for the attenuated network
## enrichment analysis
enrichment_ARI_genes <- list()
for(i in names(cp_genes)){
  message(i)
  data <- cp_genes[[i]]
  enrichment <- list()
  for(region in names(data)){
    reg <- data[[region]]
    enriched <- enrichment_analysis(reg, total_genes)
    enrichment <- append(enrichment, list(enriched))
  }
  names(enrichment) <- names(data)
  enrichment_ARI_genes <- append(enrichment_ARI_genes, list(enrichment))
}
names(enrichment_ARI_genes) <- names(cp_genes)
enrichment_ARI_genes_attenuated <- enrichment_ARI_genes
rm(enrichment_ARI_genes)

for(i in names(enrichment_ARI_genes_attenuated)){
  data <- enrichment_ARI_genes_attenuated[[i]]
  for(region in names(data)){
    reg <- data[[region]]
    for(database in names(reg)){
      if(region == names(data)[1]){
        if(database == "GO"){
          wb = createWorkbook()
          n = paste0(region, "-", database)
          addWorksheet(wb, n)
          writeData(wb, n, reg[[database]])
          saveWorkbook(wb, file = paste0("output_path","enrichment_region_", i, ".xlsx"), overwrite = TRUE)
        }
        else {
          n = paste0(region, "-", database)
          addWorksheet(wb, n)
          writeData(wb, n, reg[[database]])
          saveWorkbook(wb, file = paste0("output_path", "enrichment_region_", i, ".xlsx"), overwrite = TRUE)
        }
      }
      else {
        n = paste0(region, "-", database)
        addWorksheet(wb, n)
        writeData(wb, n, reg[[database]])
        saveWorkbook(wb, file = paste0("output_path", "enrichment_region_", i, ".xlsx"), overwrite = TRUE)
      }
    }
  }
}

### from this point on - just look at the regional comparisons significantly overpatterned in FC
sig_comps = overpatterned
names(sig_comps) = rep("FC", length(sig_comps))
sig_comps = gsub("-", "_v_", sig_comps)

sig_comp_list_fc <- true_diff_genes_final_fc[which(names(true_diff_genes_final_fc) %in% sig_comps[which(names(sig_comps)=="FC")])] 

### remove any genes from the FC that are present in ctrl
sig_comp_list_fc_unique = list()

for(comp in names(sig_comp_list_fc)){
  print(comp)
  rm = true_diff_genes_final_ctrl[[comp]][which(true_diff_genes_final_ctrl[[comp]] %in% sig_comp_list_fc[[comp]])]
  print(length(rm))
  if(length(rm)==0){
    sig_comp_list_fc_unique[[comp]]=sig_comp_list_fc[[comp]]
  }else{
    sig_comp_list_fc_unique[[comp]]=sig_comp_list_fc[[comp]][-match(rm,sig_comp_list_fc[[comp]])]
  }
}
sig_comp_list_fc = sig_comp_list_fc_unique

names(sig_comp_list_fc) <- paste(names(sig_comp_list_fc),"_FC-gt-CTRL",sep="")
sig_comp_list = sig_comp_list_fc

cp_genes_fc = list()
for(comp in names(sig_comp_list)){
  
  fc_att_genes = sig_comp_list[[comp]]
  
  idx_name = gsub("_FC-gt-CTRL", "", comp)
  idx_name = gsub("_v_", "-",idx_name)
  
  reg1 = strsplit(comp, split="_")[[1]][1]
  reg2 = strsplit(comp,split="_")[[1]][3]
  
  idx_name = paste0(reg1, " - ", reg2)
  fc_att_expr = datExpr.reg[which(rownames(datExpr.reg) %in% fc_att_genes),idx_comp_all[[idx_name]]]
  datMeta_fc_att = datMeta[idx_comp_all[[idx_name]],]
  datMeta_fc_att$Dx_Reg = paste(datMeta_fc_att$group,datMeta_fc_att$region,sep="_")
  
  ## sort the genes which are more highly expressed in reg1 or reg2
  
  count=rep(NA,length(fc_att_genes))
  names(count) <- rownames(fc_att_expr)
  for(gene in rownames(fc_att_expr)){
    avg_1=median(fc_att_expr[gene,which(datMeta_fc_att$Dx_Reg==paste("FC",reg1,sep="_"))])
    avg_2=median(fc_att_expr[gene,which(datMeta_fc_att$Dx_Reg==paste("FC",reg2,sep="_"))])
    if(avg_1 > avg_2){
      count[gene]=reg1
    }else{
      count[gene]=reg2
    }
  }
  
  reg1_genes = names(count)[which(count==reg1)]
  reg2_genes = names(count)[which(count==reg2)]
  cp_genes_fc[[comp]][[reg1]]=reg1_genes
  cp_genes_fc[[comp]][[reg2]]=reg2_genes
  
}

save(cp_genes_fc,file=paste0(output_workenvironment, "/ARI_genes_from_RegComparisons_Overpatterned_inFC_RegSorted.RData"))

# Enrichment analysis
enrichment_ARI_genes <- list()
for(i in names(cp_genes_fc)){
  message(i)
  data <- cp_genes_fc[[i]]
  enrichment <- list()
  for(region in names(data)){
    reg <- data[[region]]
    enriched <- enrichment_analysis(reg, total_genes)
    enrichment <- append(enrichment, list(enriched))
  }
  names(enrichment) <- names(data)
  enrichment_ARI_genes <- append(enrichment_ARI_genes, list(enrichment))
}
names(enrichment_ARI_genes) <- names(cp_genes_fc)
enrichment_ARI_genes_overpatterned <- enrichment_ARI_genes
rm(enrichment_ARI_genes)

for(i in names(enrichment_ARI_genes_overpatterned)){
  data <- enrichment_ARI_genes_overpatterned[[i]]
  for(region in names(data)){
    reg <- data[[region]]
    for(database in names(reg)){
      if(region == names(data)[1]){
        if(database == "GO"){
          wb = createWorkbook()
          n = paste0(region, "-", database)
          addWorksheet(wb, n)
          writeData(wb, n, reg[[database]])
          saveWorkbook(wb, file = paste0(output_workenvironment, "enrichment_overpatterned_region_", i, ".xlsx"), overwrite = TRUE)
        }
        else {
          n = paste0(region, "-", database)
          addWorksheet(wb, n)
          writeData(wb, n, reg[[database]])
          saveWorkbook(wb, file = paste0(output_workenvironment, "enrichment_overpatterned_region_", i, ".xlsx"), overwrite = TRUE)
        }
      }
      else {
        n = paste0(region, "-", database)
        addWorksheet(wb, n)
        writeData(wb, n, reg[[database]])
        saveWorkbook(wb, file = paste0(output_workenvironment, "enrichment_overpatterned_region_", i, ".xlsx"), overwrite = TRUE)
      }
    }
  }
}

