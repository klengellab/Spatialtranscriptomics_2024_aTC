**Spatial transcriptomics reveals modulation of transcriptional networks across 
brain regions after auditory threat conditioning**

**Abstract**

**Data availability:**

* Raw data: The stRNAseq fastq files can be accessed through GEO accession GSE243140
  
**Guide through the data analysis scripts:**

pipeline stages:
1. preprocessing
2. DEG_analysis
3. Deconvolution
4. network_analysis
   
The scripts are ordered by number. To replicate our analysis the start is with preprocessing of the data that entails data loading, cleaning, normalisation and clustering.
Accordingly, we performed differential expression analysis to obtain the differences between ctrl and aTC per cluster, which can be found in the DEG_analysis folder. There 
you start with pseudo-bulk of the data, differential gene expression (DEG) analysis carried out by Limma and in script 2-1 you have a closer look at the output of the DEG analysis
and perform gene set enrichment analysis (GSEA). The Deconvolution and network_analysis can be performed simultaneously but please start with the [number]-0 script and follow accordingly.
We also have a folder entailing general functions that you need to load in to get through the scripts.


