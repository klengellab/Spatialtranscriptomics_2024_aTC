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

# Read in Data
data <- readRDS(paste0(output_workenvironment, "/visium_normalised_entire_brain.RDS"))
metadata <- readRDS(paste0(data_path, "/visium_metadata_normalised_entire_brain.RDS"))

# Contrast matrix
cont_mat=t(combn(paste(unique(metadata$region),sep=""), 2))
colnames(cont_mat)=c("Region1","Region2")

idx_comp_all=list()
for(i in c(1:nrow(cont_mat))){
  idx=which(metadata$region==cont_mat[i,1] | metadata$region==cont_mat[i,2])
  metadata_comp=metadata[idx,]
  idx_comp=list()
  for(dx in c("FC","CTRL")){
    idx2=which(metadata_comp$group==dx)
    metadata_comp_dx=metadata_comp[idx2,]
    idx_keep=which(duplicated(metadata_comp_dx$sample)==TRUE)
    subs_keep=metadata_comp_dx$sample[idx_keep]
    idx_comp[[dx]]=which(metadata$sample %in% subs_keep & (metadata$region==cont_mat[i,1] | metadata$region==cont_mat[i,2]) )
    tmp = metadata[idx_comp[[dx]],]
    tmp=tmp[order(tmp$sample),]
    idx_comp[[dx]]=idx_comp[[dx]]
  }
  idx_comp_all[[paste(cont_mat[i,1],cont_mat[i,2],sep=" - ")]]=c(idx_comp[["CTRL"]],idx_comp[["FC"]])
}

n_vec_fc <- n_vec_ctrl <- rep(NA,528)
names(n_vec_fc) <- names(n_vec_ctrl) <- names(idx_comp_all)

for(comp in names(idx_comp_all)){
  datMeta_comp=metadata[idx_comp_all[[comp]],]
  datMeta_unq=datMeta_comp[!duplicated(datMeta_comp$group),]
  fc_subs=length(datMeta_unq$group[which(datMeta_unq$Diagnosis=="FC")])
  ctrl_subs=length(datMeta_unq$group[which(datMeta_unq$Diagnosis=="CTRL")])
  n_vec_fc[comp]=fc_subs
  n_vec_ctrl[comp]=ctrl_subs
}
save(idx_comp_all,n_vec_fc,n_vec_ctrl,file=paste0(output_workenvironment, "/Permutation_Index.RData"))

p_list=list()
metadata$group <- as.factor(metadata$group)
for(comp in names(idx_comp_all)){
  print(comp)
  for(dx in c("FC","CTRL")){
    print(dx)
    idx_dx=which(metadata$group[idx_comp_all[[comp]]]==dx)
    
    comp_meta_mat=metadata[idx_comp_all[[comp]][idx_dx],]
    comp_expr_mat=logCPM_offset[,idx_comp_all[[comp]][idx_dx]]
    names=unique(comp_meta_mat$sample)
    
    reg1_name=strsplit(comp,split=" - ")[[1]][1]
    reg2_name=strsplit(comp,split=" - ")[[1]][2]
    
    reg1=comp_expr_mat[,which(comp_meta_mat$region==reg1_name)]
    reg2=comp_expr_mat[,which(comp_meta_mat$region==reg2_name)]
    colnames(reg1)=comp_meta_mat$sample[which(comp_meta_mat$region==reg1_name)]
    colnames(reg2)=comp_meta_mat$sample[which(comp_meta_mat$region==reg2_name)]
    
    reg1=reg1[,match(names,colnames(reg1))]
    reg2=reg2[,match(names,colnames(reg2))]
    
    ttable=data.frame("V"=rep(NA,nrow(comp_expr_mat)),"PVal"=rep(NA,nrow(comp_expr_mat)))
    rownames(ttable)=rownames(comp_expr_mat)
    
    for(i in c(1:nrow(comp_expr_mat))){
      if (i%%1000 == 0) {print(i)}
      tmp=wilcox.test(reg1[i,],reg2[i,],paired=TRUE)
      ttable[i,1]=tmp$statistic
      ttable[i,2]=tmp$p.value
    }
    cp_list[[comp]][[dx]]$data=ttable
    cp_list[[comp]][[dx]]$n=ncol(reg1)
  }
}
saveRDS(cp_list, paste0(output_workenvironment, "/03_01_A_01_True_TRI_DGE_CompleteList.RDS"))

### Compile results into a matrix
cp_mat <- data.frame(matrix(NA,nrow=length(names(cp_list)),ncol=5))
rownames(cp_mat) <- names(cp_list)
colnames(cp_mat) <- c("Comparison","Num_DGE_FC","Num_DGE_CTRL","Num_Samps_FC","Num_Samps_CTRL")
cp_mat$Comparison <- names(cp_list)
for(comp in cp_mat$Comparison){
  tmp=cp_list[[which(names(cp_list)==comp)]]$FC$data
  tmp$FDR=p.adjust(tmp$PVal,method="fdr")
  cp_mat$Num_DGE_FC[which(cp_mat$Comparison==comp)]=length(which(tmp$FDR < 0.05))
  tmp=cp_list[[which(names(cp_list)==comp)]]$CTRL$data
  tmp$FDR=p.adjust(tmp$PVal,method="fdr")
  cp_mat$Num_DGE_CTRL[which(cp_mat$Comparison==comp)]=length(which(tmp$FDR < 0.05))
  cp_mat$Num_Samps_FC[which(cp_mat$Comparison==comp)]=cp_list[[which(names(cp_list)==comp)]]$FC$n
  cp_mat$Num_Samps_CTRL[which(cp_mat$Comparison==comp)]=cp_list[[which(names(cp_list)==comp)]]$CTRL$n
}
saveRDS(cp_mat, paste0(output_workenvironment, "/True_TRI_DGE_SummarizedMatrix.RDS"))

### Get the difference in the absolute number of DE genes between regions for aTC and CTRL
diff_list=list()
### take difference of num dge CTRL - num DGE aTC. Positive diff = CTRL greater, Neg diff = aTC greater.
for(comp in names(idx_comp_all)){
  diff_list[[comp]]$Difference=cp_mat[comp,3]-cp_mat[comp,2]
  diff_list[[comp]]$FC_N=cp_mat[comp,4]
  diff_list[[comp]]$CTRL_N=cp_mat[comp,5]
}
saveRDS(diff_list, paste0(output_workenvironment, "/True_TRI_DGE_FC_v_Control.RData"))

# Permutation analysis
### It is recommended to run individual permutations in parallel on a high-performance computing cluster.
## In total we ran 10.000 permutations
## Make sure that your output workenvironment path is different than the general output workenvironment path since you later on
## Need to coutn the number of files stored there
for(iter in c(1:10000)){
  set.seed(iter)
  
  datMeta_model_perm=list()
  
  for(comp in names(idx_comp_all)){
    datMeta_tmp=metadata[idx_comp_all[[comp]],]
    datMeta_tmp_unq=metadata[!duplicated(metadata$sample),]
    dx_permute=c(rep("CTRL",length(which(datMeta_tmp_unq$group=="CTRL"))),
                 rep("FC",length(which(datMeta_tmp_unq$group =="FC"))))
    dx_permute=sample(dx_permute)
    datMeta_tmp_unq$group=dx_permute
    datMeta_tmp$group=datMeta_tmp_unq$group[match(datMeta_tmp$sample,datMeta_tmp_unq$sample)]
    datMeta_model_perm[[comp]]=datMeta_tmp
  }
  
  datExpr.reg_perm=logCPM_offset
  colnames(datExpr.reg_perm)=paste(metadata$sample,metadata$region,sep="_")
  
  cp_list=list()
  
  for(comp in names(idx_comp_all)){
    print(comp)
    for(dx in c("FC","CTRL")){
      print(dx)
      idx_dx=which(datMeta_model_perm[[comp]]$group==dx)
      comp_meta_mat=datMeta_model_perm[[comp]][idx_dx,]
      names_expr=paste(comp_meta_mat$sample,comp_meta_mat$region,sep="_")
      
      comp_expr_mat=datExpr.reg_perm[,match(names_expr,colnames(datExpr.reg_perm))]
      
      names=unique(comp_meta_mat$sample)
      
      reg1_name=strsplit(comp,split=" - ")[[1]][1]
      reg2_name=strsplit(comp,split=" - ")[[1]][2]
      
      reg1=comp_expr_mat[,which(comp_meta_mat$region==reg1_name)]
      reg2=comp_expr_mat[,which(comp_meta_mat$region==reg2_name)]
      colnames(reg1)=comp_meta_mat$sample[which(comp_meta_mat$region==reg1_name)]
      colnames(reg2)=comp_meta_mat$sample[which(comp_meta_mat$region==reg2_name)]
      
      reg1=reg1[,match(names,colnames(reg1))]
      reg2=reg2[,match(names,colnames(reg2))]
      
      ttable=data.frame("V"=rep(NA,nrow(comp_expr_mat)),"PVal"=rep(NA,nrow(comp_expr_mat)))
      rownames(ttable)=rownames(comp_expr_mat)
      
      for(i in c(1:nrow(comp_expr_mat))){
        if (i%%1000 == 0) {print(i)}
        tmp=wilcox.test(reg1[i,],reg2[i,],paired=TRUE)
        ttable[i,1]=tmp$statistic
        ttable[i,2]=tmp$p.value
      }
      
      cp_list[[comp]][[dx]]$data=ttable
      cp_list[[comp]][[dx]]$n=ncol(reg1)
    }
  }
  
  save(cp_list,file=paste("your output workenvironment path", "/TRI_list_permutation",iter,".RData", sep=""))
  
}
# Check if you have indeed 10.000 files

### Compile permutations
"data_user/03_01_A_01_permutation/"

cp_mat_list=list()
for(i in c(1:length(list.files("your output workenvironment path")))){
  if (i%%50 == 0) {print(paste(i,"/",length(list.files("your output workenvironment path")),sep=""))}
  load(paste("your output workenvironment path", "/TRI_list_permutation",i,".RData",sep=""))
  cp_mat_list[[i]] <- data.frame(matrix(NA,nrow=length(names(cp_list)),ncol=5))
  rownames(cp_mat_list[[i]]) <- names(cp_list)
  colnames(cp_mat_list[[i]]) <- c("Comparison","Num_DGE_FC","Num_DGE_CTRL","Num_Samps_FC","Num_Samps_CTRL")
  cp_mat_list[[i]]$Comparison <- names(cp_list)
  
  for(comp in cp_mat_list[[i]]$Comparison){
    tmp=cp_list[[which(names(cp_list)==comp)]]$FC$data
    tmp$FDR=p.adjust(tmp$PVal,method="fdr")
    cp_mat_list[[i]]$Num_DGE_FC[which(cp_mat_list[[i]]$Comparison==comp)]=length(which(tmp$FDR < 0.05))
    tmp=cp_list[[which(names(cp_list)==comp)]]$CTRL$data
    tmp$FDR=p.adjust(tmp$PVal,method="fdr")
    cp_mat_list[[i]]$Num_DGE_CTRL[which(cp_mat_list[[i]]$Comparison==comp)]=length(which(tmp$FDR < 0.05))
    cp_mat_list[[i]]$Num_Samps_FC[which(cp_mat_list[[i]]$Comparison==comp)]=cp_list[[which(names(cp_list)==comp)]]$FC$n
    cp_mat_list[[i]]$Num_Samps_CTRL[which(cp_mat_list[[i]]$Comparison==comp)]=cp_list[[which(names(cp_list)==comp)]]$CTRL$n
  }
}
saveRDS(cp_mat_list, paste0(output_workenvironment, "cp_matrix_list.RDS"))

## Now, use permutations to create a null distribution to test if 
### the # of DE genes between aTC and CTRLS is significant for each regional pair.

perm_diff_list=list()

for(comp in names(idx_comp_all)){
  perm_diff_list[[comp]]$data=data.frame(matrix(NA,nrow=10000,ncol=4))
  colnames(perm_diff_list[[comp]]$data)=c("Iteration","Difference","FC_Num_DGE","CTRL_Num_DGE")
  rownames(perm_diff_list[[comp]]$data)=seq(1:10000)
}

for(i in c(1:10000)){
  if (i%%50 == 0) {print(paste(i,"/10000",sep=""))}
  cp_mat_tmp=cp_mat_list[[i]]
  for(comp in names(idx_comp_all)){
    perm_diff_list[[comp]]$data$Iteration[i]=i
    perm_diff_list[[comp]]$data$Difference[i]=cp_mat_tmp[comp,3]-cp_mat_tmp[comp,2]
    perm_diff_list[[comp]]$data$FC_Num_DGE[i]=cp_mat_tmp[comp,2]
    perm_diff_list[[comp]]$data$CTRL_Num_DGE[i]=cp_mat_tmp[comp,3]
  }
}
saveRDS(perm_diff_list, paste0(output_workenvironment, "/perm_diff_list.RDS"))



ttable=data.frame(matrix(NA,nrow=length(cp_mat_list[[1]]$Comparison),ncol=6))
colnames(ttable) = c("Comparison","CTRL_FC_NumDGE_MeanDiff","Perm_PVal","FDR_adj_Perm_PVal","FC_N","CTRL_N")
rownames(ttable) <- ttable$Comparison <- cp_mat_list[[1]]$Comparison

pdf(paste0(output_path, "/03_01_A_01_TRI_Permutation.pdf"))
for( comp in names(idx_comp_all)){
  cent=round(scale(perm_diff_list[[comp]]$data$Difference,scale=FALSE),0)
  mn=round(mean(perm_diff_list[[comp]]$data$Difference),0)
  diff=diff_list[[comp]]$Difference-mn
  ttable$CTRL_FC_NumDGE_MeanDiff[which(rownames(ttable)==comp)]=diff
  
  more=length(which(abs(cent) >= abs(diff)))
  p = (more)/10001
  ttable$Perm_PVal[which(rownames(ttable)==comp)]=signif(p,5)
  ttable$FC_N[which(rownames(ttable)==comp)]=diff_list[[comp]]$FC_N
  ttable$CTRL_N[which(rownames(ttable)==comp)]=diff_list[[comp]]$CTRL_N
  
  plot.cent=data.frame(cent)
  
  print(ggplot(plot.cent,aes(x=cent))+
          geom_density()+
          ggtitle(paste(comp,"; Real_Diff=",diff,"; p=",signif(p,3),"; N_CTRL=",diff_list[[comp]]$CTRL_N,
                        "; N_FC=",diff_list[[comp]]$FC_N,sep=""))+
          xlab("CTRL_NDiff-FC_NDiff (Mean Centered)")+
          theme(plot.title = element_text(hjust = 0.5))+
          geom_vline(xintercept=diff,col="red"))
}

dev.off()
ttable$FDR_adj_Perm_PVal=p.adjust(ttable$Perm_PVal,method="fdr")
write.csv(ttable,file=paste0(output_path, "/Permutation_RegComp-TTable.csv"))


### Obtain genes from the Wilcoxon test (with true samples) that are missing in 95% of permutations
### It is recommended to run the iterations below in parallel on a high-performance computing cluster.

### First, extract the permuted genes for each comparison
idx_list=list(c(1:1000),c(1001:2000),c(2001:3000),c(3001:4000),c(4001:5000),
              c(5001:6000),c(6001:7000),c(7001:8000),c(8001:9000),c(9001:10000))

for(iter_idx in c(1:10)){
  
  # This is that they filter on FDR significant genes
  idx = idx_list[[as.numeric(iter_idx)]]
  
  perm_diff_genes_ctrl=list()
  perm_diff_genes_fc=list()
  
  for(i in c(1:length(idx))){
    if (i%%50 == 0) {print(paste(i,"/1000",sep=""))}
    load(paste0("your output workenvironment path", "/TRI_list_permutation",as.character(idx[i]),".RData",sep=""))
    names(cp_list) <- gsub(" - ","_v_",names(cp_list))
    for(comp in names(cp_list)){
      if(i==1){
        perm_diff_genes_ctl[[comp]] <- list()
        perm_diff_genes_fc[[comp]] <- list()
      }
      tmp=cp_list[[which(names(cp_list)==comp)]]$FC$data
      tmp$FDR=p.adjust(tmp$PVal,method="fdr")
      perm_diff_genes_fc[[comp]][[i]]=rownames(tmp)[which(tmp$FDR < 0.05)]
      tmp=cp_list[[which(names(cp_list)==comp)]]$CTRL$data
      tmp$FDR=p.adjust(tmp$PVal,method="fdr")
      perm_diff_genes_ctrl[[comp]][[i]]=rownames(tmp)[which(tmp$FDR < 0.05)]
    }
  }
  
  save(perm_diff_genes_ctrl,perm_diff_genes_fc,
       file=paste("your output workenvironment path", "/TRI_permutation_genes_Perm",iter_idx,".RData",sep=""))
  
}

### now compile based on comparison
perm_diff_genes_fc_comp <- list()
perm_diff_genes_ctrl_comp <- list()
for(j in c(1:10)){
  print(j)
  load(paste("your output workenvironment path", "/TRI_permutation_genes_Perm",j,".RData",sep=""))
  for(comp in names(cp_list)){
    print(comp)
    if(j==1){
      perm_diff_genes_fc_comp[[comp]] <- vector(mode="list",length=10000)
      perm_diff_genes_ctrl_comp[[comp]] <- vector(mode="list",length=10000)
    }
    perm_diff_genes_fc_comp[[comp]][idx_list[[j]]] = perm_diff_genes_fc[[comp]]
    perm_diff_genes_ctrl_comp[[comp]][idx_list[[j]]] = perm_diff_genes_ctrl[[comp]]
  }
}

for(comp in names(perm_diff_genes_fc_comp)){
  perm_diff_genes_fc_comp_single = perm_diff_genes_fc_comp[[comp]]
  perm_diff_genes_ctrl_comp_single = perm_diff_genes_ctrl_comp[[comp]]
  save(perm_diff_genes_fc_comp_single,perm_diff_genes_ctrl_comp_single,
       file=paste("your output workenvironment path","/Compiled_TRI_permutation_genes_",comp,".RData",sep=""))
}

### Now, count the number of true diff genes in the permutations 

for(iter_idx in c(1:10)){
  for(dx in c("CTRL","FC")){
    
    true_diff_genes_ctrl=list()
    true_diff_genes_fc=list()
    
    true_diff_genes_count_ctrl=list()
    true_diff_genes_count_fc=list()
    
    load(paste0(output_workenvironment, "/03_01_A_01_True_TRI_DGE_CompleteList.RData"))

    names(cp_list) <- gsub(" - ","_v_",names(cp_list))
    
    comp = names(cp_list)[as.numeric(iter_idx)]
    dx = as.character(dx)  
    
    "Getting True CP Genes"
    
    tmp=cp_list[[which(names(cp_list)==comp)]]$FC$data
    tmp$FDR=p.adjust(tmp$PVal,method="fdr")
    true_diff_genes_fc[[comp]]=rownames(tmp)[which(tmp$FDR < 0.05)]
    tmp=cp_list[[which(names(cp_list)==comp)]]$CTRL$data
    tmp$FDR=p.adjust(tmp$PVal,method="fdr")
    true_diff_genes_ctrl[[comp]]=rownames(tmp)[which(tmp$FDR < 0.05)]
    
    "Counting Number of 'True Gene' Occurences in Permutations"
    
    true_diff_genes_count_fc[[comp]]=rep(NA,length(true_diff_genes_fc[[comp]]))
    names(true_diff_genes_count_fc[[comp]]) <- true_diff_genes_fc[[comp]]
    true_diff_genes_count_ctrl[[comp]]=rep(NA,length(true_diff_genes_ctrl[[comp]]))
    names(true_diff_genes_count_ctrl[[comp]]) <- true_diff_genes_ctrl[[comp]]
    
    load(paste("your output workenvironment path","/Compiled_TRI_permutation_genes_",comp,".RData",sep=""))
    
    if(dx == "FC"){
      
      print("Counting for FC")
      
      ind=1
      for(gene in true_diff_genes_fc[[comp]]){
        count=0
        if (ind%%50 == 0) {print(paste(ind,"/",length(true_diff_genes_fc[[comp]]),sep=""))}
        for(i in c(1:length(perm_diff_genes_fc_comp_single))){
          iter=length(grep(gene,perm_diff_genes_fc_comp_single[[i]]))
          count=count+iter
        }
        true_diff_genes_count_fc[[comp]][gene]=count
        ind=ind+1
      }
      
      save(true_diff_genes_count_fc,true_diff_genes_fc,
           file=paste("your output workenvironment", "/TrueGenes_in_PermutedGenes_Compiled_",as.character(comp),"_FC.RData",sep=""))
      
      
    }else if(dx == "CTRL"){
      
      print("Counting for CTRL")
      
      ind=1
      for(gene in true_diff_genes_ctrl[[comp]]){
        count=0
        if (ind%%50 == 0) {print(paste(ind,"/",length(true_diff_genes_ctrl[[comp]]),sep=""))}
        for(i in c(1:length(perm_diff_genes_ctrl_comp_single))){
          iter=length(grep(gene,perm_diff_genes_ctrl_comp_single[[i]]))
          count=count+iter
        }
        true_diff_genes_count_ctrl[[comp]][gene]=count
        ind=ind+1
      }
      
      save(true_diff_genes_count_ctrl,true_diff_genes_ctrl,
           file=paste("your output workenvironment", "/TrueGenes_in_PermutedGenes_Compiled_",as.character(comp),"_CTRL.RData",sep=""))
    } 
  }
}

### Combine CTRL and FC data
true_diff_genes_ctrl_all <- vector(mode="list",length=55)
true_diff_genes_fc_all <- vector(mode="list",length=55)

true_diff_genes_count_ctrl_all <- vector(mode="list",length=55)
true_diff_genes_count_fc_all <- vector(mode="list",length=55)

load(paste0(output_workenvironment, "/03_01_A_01_True_TRI_DGE_CompleteList.RDS"))

names(cp_list) <- gsub(" - ","_v_",names(cp_list))

names(true_diff_genes_ctrl_all) <- names(true_diff_genes_fc_all) <- names(cp_list)
names(true_diff_genes_count_ctrl_all) <- names(true_diff_genes_count_fc_all) <- names(cp_list)

for(comp in names(cp_list)){
  print(comp)
  for(dx in c("FC","CTRL")){
    print(dx)
    if(dx == "FC"){
      
      load(paste("your output workenvironment", "/TrueGenes_in_PermutedGenes_Compiled_",comp,"_FC.RData",sep=""))
     
      true_diff_genes_fc_all[[comp]] <- true_diff_genes_fc[[comp]]        
      true_diff_genes_count_fc_all[[comp]] <- true_diff_genes_count_ctrl[[comp]]
      
    }else if(dx == "CTRL"){
      
      load(paste("your output workenvironment", "/TrueGenes_in_PermutedGenes_Compiled_",comp,"_CTRL.RData",sep=""))
      
      true_diff_genes_ctrl_all[[comp]] <- true_diff_genes_ctrl[[comp]]        
      true_diff_genes_count_ctrl_all[[comp]] <- true_diff_genes_count_ctrl[[comp]]
      
    }
  }
}
save(true_diff_genes_count_fc_all,true_diff_genes_count_ctrl_all,
     true_diff_genes_fc_all,true_diff_genes_ctrl_all,
     file=paste0(output_workenvironment, "/TrueGenes_in_PermutedGenes_All_Compiled.RData"))
