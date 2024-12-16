#!/usr/local/bin/Rscript
# title: "RBM33 is a unique m(6)A RNA-binding protein that regulates ALKBH5 demethylase activity and substrate selectivity"
# author: "Ruihan Luo"
# date: "April 19th,2024"
# rm(list=ls())
library(DESeq2)
library(edgeR)
library(data.table)
library(stringr)
library(tidyverse) 
library(GEOquery)
library(readxl)
Sys.setenv("VROOM_CONNECTION_SIZE" = 999999999)
getOption('timeout')
options(timeout=100000000)
path <- c('/home/rluo4/R/x86_64-pc-linux-gnu-library/4.3', '/opt/R/4.3.1/lib64/R/library', '/home/rluo4/R/x86_64-conda-linux-gnu-library/4.3', '/data2/rluo4/bin/miniconda3/lib/R/library')
.libPaths(path)
raw_directory <- '/data2/rluo4/EpiTrans/RMDatasets/Raw'
data_directory <- '/data2/rluo4/EpiTrans/RMDatasets/GEO'
setwd(data_directory)
plot_dir <- '/data2/rluo4/RPMfunc/Output/GEO'
# 1) load the rdata from RMdeg.R
load(paste0('/data2/rluo4/EpiTrans/RMDatasets/allDat.RData'))
load(paste0('/data2/rluo4/EpiTrans/RMDatasets/RMP_allDat.RData'))
load(paste0('/data2/rluo4/EpiTrans/RMDatasets/RMdeg_allDat.RData'))
RMP_update <- read_excel('/data2/rluo4/EpiTrans/DataCollection/RMP_update.xlsx')
asp_sc <- fread(file="/data2/rluo4/EpiTrans/DataCollection/asp_RMP_epitrans.link.txt", header = F)
# remove.packages("MatrixGenerics", lib="/data2/rluo4/lorihan/R/site-library")
# BiocManager::install("MatrixGenerics")
# remove.packages("Rgraphviz", lib="/data2/rluo4/lorihan/R/site-library")
# BiocManager::install("Rgraphviz")
# BiocManager::install("DiffBind")
# 2) load the functions
extract_names <- function(df) {
  colnames(df)
}
# Function to transform the strings
transform_path <- function(x) {
  x <- gsub("^X\\.", "/", x)
  x <- gsub("\\.(?=.*_)", "/", x, perl = TRUE)
  return(x)
}
# # Replace 'X.' with '/'
# modified_string <- gsub("^X\\.", "/", original_string)
# # Replace the last '_' with '.'
# modified_string <- gsub("_(?=.*\\.bam)", ".", modified_string, perl = TRUE)
# Define a function to extract gene_name from the attribute string
extract_gene_name <- function(attribute_string) {
  # Use a regular expression to match the gene_name
  gene_name <- str_extract(attribute_string, "gene_name [^;]+")
  # Remove the "gene_name " prefix to get the actual gene name
  gene_name <- str_replace(gene_name, "gene_name ", "")
  gene_name <- str_replace_all(gene_name, '"', '')
  return(gene_name)
}
process_celine <- function(c, resM, resP, cellName, data_directory, db, genome_version, i, RMPdeg_datasets, n) {
  print(c)
  cl <- paste0(c, '_vs_WT')
  index = grepl('DESeq2', cellName)
  
  if (length(cellName) >= 2 & length(unique(index)) > 1) {
    a <- rownames(resM[[paste0(cl, '_DESeq2')]])
    b <- rownames(resM[[paste0(cl, '_limma')]])
    peaks <- union(a, b)
  } else {
    b <- rownames(resM[[paste0(cl, '_limma')]])
    peaks <- b
  }
  
  merge_peak <- resP[[cl]]
  if(i=='GSE122961'){
    merge_bed <- merge_peak[, c('chr', 'start', 'end', 'name', 'strand')]
    # merge_bed$Stat <- merge_peak$score
    colnames(merge_bed) <- c('chr', 'start', 'end', 'Stat', 'strand')
  }  else if (i=='GSE38957' & c =='Hela_NSUN2PD'){
    merge_bed <- merge_peak[, c('V2', 'V3', 'V4', 'V1', 'V5')]
    # merge_bed$Stat <- merge_peak$score
    colnames(merge_bed) <- c('chr', 'start', 'end', 'Stat', 'strand')
  } else {
    merge_bed <- merge_peak[-1, match(c('chr', 'start', 'end', 'Stat', 'strand'), merge_peak[1,])]
    colnames(merge_bed) <- c('chr', 'start', 'end', 'Stat', 'strand')
    merge_bed$Stat <- merge_peak$V1[-1]
  }
  i <- RMPdeg_datasets[n, 'GSE']
  
  output_path <- file.path(data_directory, i, "metaPlotR", paste0(c, "_merge_peak.bed"))
  write.table(merge_bed, output_path, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
  
  HG <- genome_version$HG[genome_version$GSE == i]
  # cmd <- paste0("bedtools intersect -wa -wb -a ", output_path, " -b ", db, HG, ".refGene.gtf.gz > ", 
  #               file.path(data_directory, i, "metaPlotR", paste0(c, "_merge_peak.anno.bed")))
  cmd <- paste0("bedtools intersect -wa -wb -a ", output_path, " -b ", db, "gencode.v19.annotation.gtf.gz > ", 
                file.path(data_directory, i, "metaPlotR", paste0(c, "_merge_peak.anno.bed")))
  if(HG=='hg38'){
    cmd <- paste0("bedtools intersect -wa -wb -a ", output_path, " -b ", db, "gencode.v36.annotation.gtf.gz > ", 
                  file.path(data_directory, i, "metaPlotR", paste0(c, "_merge_peak.anno.bed")))
  }
  merge_bed_anno_path <- file.path(data_directory, i, "metaPlotR", paste0(c, "_merge_peak.anno.bed"))
  # if( !file.exists(merge_bed_anno_path) ){
  system(cmd)
  # }
  merge_bed_anno <- fread(merge_bed_anno_path, sep = '\t')
  # merge_bed_anno$V14 <- gsub('gene_id ', '', str_split(merge_bed_anno$V14, ';', simplify = TRUE)[, 1]) #for annotation of hg19.refGene.gtf.gz
  # merge_bed_anno$V14 <- gsub(' gene_name ', '', str_split(merge_bed_anno$V14, ';', simplify = TRUE)[, 4]) 
  # merge_bed_anno$V14 <- sapply(merge_bed_anno$V14, extract_gene_name)
  merge_bed_anno$V14 <- extract_gene_name(merge_bed_anno$V14)
  # Print the results
  merge_bed_anno <- merge_bed_anno[!duplicated(merge_bed_anno$V4), c(1:4, 14)]
  colnames(merge_bed_anno) <- c('chr', 'start', 'end', 'peak_id', 'geneSymbol')
  print(length(unique(merge_bed_anno$geneSymbol)))
  df_peak <- merge_bed_anno[merge_bed_anno$peak_id %in% peaks, ]
  df_gene <- unique(df_peak$geneSymbol)
  print(paste0('diff methylated genes in ', cl, ': ', length(df_gene)))
  return(df_peak)
}
process_samples <- function(n, RMPdeg_datasets, data_directory, RMdem_list, RMdep_list, genome_version) {
  # Extract the current dataset information
  print(RMPdeg_datasets[n, ])
  i <- RMPdeg_datasets[n, 'GSE']
  # Fetch RNA modification and deposition lists for the current dataset
  resM <- RMdem_list[[i]]
  resP <- RMdep_list[[i]]
  # if( is.null(resM) ){
  #   print(paste0(i, ": no RNA modification data !"))
  #   next;
  # } This works for for loops 
  # Check if RNA modification data is available
  if (is.null(resM)|is.null(resP)) {
    print(paste0(i, ": no RNA modification/Peak data!"))
    return(NULL)  # Exit the function early if no data
  }
  setwd(file.path(data_directory, i))
  # Extract unique cell lines from cell names
  cellName <- names(resM)
  celine <- unique(str_split(cellName, '_vs', simplify = TRUE)[, 1])
  print(celine)
  
  db <- '/data2/rluo4/hg38/'
  resMeR_list <- vector("list", length(celine))
  names(resMeR_list) <- celine
  
  for (c in celine) {
    resMeR_list[[c]] <- process_celine(c, resM, resP, cellName, data_directory, db, genome_version, i, RMPdeg_datasets, n)
  }
  return(resMeR_list)
}

############
#gene info#
############
library(org.Hs.eg.db)
eg2symbol=toTable(org.Hs.egSYMBOL)
eg2ensembl=toTable(org.Hs.egENSEMBL)
length(unique(eg2ensembl$ensembl_id))
eg2name=toTable(org.Hs.egGENENAME)
eg2alias=toTable(org.Hs.egALIAS2EG)
eg2ncbi=toTable(org.Hs.egREFSEQ)
eg2alis_list=lapply(split(eg2alias,eg2alias$gene_id),function(x){paste0(x[,2],collapse = ";")})
GeneList=mappedLkeys(org.Hs.egSYMBOL)
if( GeneList[1] %in% eg2symbol$symbol ){
  symbols=GeneList
  geneIds=eg2symbol[match(symbols,eg2symbol$symbol),'gene_id']
}else{
  geneIds=GeneList
  symbols=eg2symbol[match(geneIds,eg2symbol$gene_id),'symbol']
}
geneNames=eg2name[match(geneIds,eg2name$gene_id),'gene_name']
geneEnsembl=eg2ensembl[match(geneIds,eg2ensembl$gene_id),'ensembl_id']
head(eg2ncbi)
geneRefSeq=eg2ncbi[match(geneIds,eg2ncbi$gene_id),'accession']
geneAlias=sapply(geneIds,function(x){ifelse(is.null(eg2alis_list[[x]]),"no_alias",eg2alis_list[[x]])})
createLink <- function(base,val) {
  sprintf('<a href="%s" class="btn btn-link" target="_blank" >%s</a>',base,val)
}
gene_info = data.frame( symbols=symbols,
                        geneIds=geneIds,#createLink(paste0("http://www.ncbi.nlm.nih.gov/gene/",geneIds),geneIds),
                        geneNames=geneNames,
                        geneEnsembls=geneEnsembl,
                        geneAlias=geneAlias,
                        geneRefSeqs=geneRefSeq,
                        stringsAsFactors = F
)
eg2CHR <- toTable(org.Hs.egMAP)
eg2TYPE <- toTable(org.Hs.egGENETYPE)
eg2GO <- toTable(org.Hs.egGO2ALLEGS)
eg2UniProt <- toTable(org.Hs.egUNIPROT)
gene_info$CHR <- eg2CHR$cytogenetic_location[match(gene_info$geneIds,eg2CHR$gene_id)]
# gene_info$geneName <- eg2name$gene_name[match(gene_info$geneID,eg2name$gene_id)]
gene_info$geneType <- eg2TYPE$gene_type[match(gene_info$geneIds,eg2TYPE$gene_id)]
gene_info$GO <- eg2GO$go_id[match(gene_info$geneIds,eg2GO$gene_id)]
gene_info$UniProtAcc <- eg2UniProt$uniprot_id[match(gene_info$geneIds,eg2UniProt$gene_id)]
# Patt <- read.table('/data2/rluo4/All/Output/organ13-deg.txt',sep = '\t',fill = TRUE,header = TRUE) #PCTanno database DEGs
############
#RMDatasets#
############
RMDatasets <- read_excel('/data2/rluo4/EpiTrans/DataCollection/Regulatory-enzymes-RNA-modifications.xlsx', sheet = 1)
unique(RMDatasets$Modifcation)
RMDatasets$GSE[duplicated(RMDatasets$GSE)]
table(RMDatasets$GSE %in% names(allDat))
setdiff(RMDatasets$GSE, names(allDat)) #PRJNA498900
table(RMDatasets$Layout)
table(RMDatasets$Modifcation)
table(RMDatasets$Treatment.Condition)
length(unique(RMDatasets$GSE))
# RMDatasets$GSE[duplicated(RMDatasets$GSE)]
RMDatasets$SRP[duplicated(RMDatasets$SRP)]
unique(RMDatasets$PMID)#168
unique(RMDatasets$GSE)#177
RMP <- read.csv(file = '/data2/rluo4/EpiTrans/DataCollection/RMP.csv')
# Identify the rows with specific characters in the "Gene" column
matching_rows <- grepl(paste(RMP$RMP_update, collapse = "|"), RMDatasets$Treatment.Condition)
# Print the dataframe with only the matching rows
# View(RMDatasets[matching_rows, ])
# View(RMDatasets[! matching_rows, ])
matching_rows1 <- grepl(paste(RMP$RMP_update, collapse = "|"), RMDatasets$Sample)
# Print the dataframe with only the matching rows
# View(RMDatasets[! (matching_rows|matching_rows1), ])
RMPdeg_datasets <- RMDatasets[ matching_rows|matching_rows1, ] #120 datasets with RMPs
unique(RMPdeg_datasets$GSE)
table(RMPdeg_datasets$Layout)
Layout <- gsub('count|FPKM|RPKM|result','matrix',RMPdeg_datasets$Layout)
table(Layout)
grep("GSE|and|PAIRED matrix", Layout, value = TRUE)
index = grep("GSE|and|PAIRED matrix", Layout)#, value = TRUE)
length(index) #77
# View(RMP_datasets[-index,])
setwd(data_directory)
GEO_dir <- unique(RMPdeg_datasets$GSE)
# GEO_dir <- gsub('./', '',list.dirs('.'))
# GEO_dir <- RMP_datasets$GSE[RMP_datasets$GSE %in% GEO_dir]
# RMP_datasets are GEO datasets with peak callings
RMP_GEO <- vector("list", length(GEO_dir))
names(RMP_GEO) <- GEO_dir
for (i in GEO_dir) {
  RMP_GEO[[i]] <- list.files(i)
}
RMPdeg_datasets$peak <- apply(RMPdeg_datasets, 1, function(x){
  i <- x['GSE']
  level3 <- grep("eak|graph|bed|bw|wig", RMP_GEO[[i]], value = TRUE)
  level3 <- ifelse(length(level3)>=1, 'peakcalling','nopeak')
  return(level3)
})
table(RMPdeg_datasets$peak)
RMPdeg_datasets <- as.data.frame(RMPdeg_datasets)
RMPeak_datasets <- RMPdeg_datasets[RMPdeg_datasets$peak=='peakcalling',]
Layout <- gsub('count|FPKM|RPKM|result','matrix',RMPeak_datasets$Layout)
table(Layout)
grep("GSE|and|PAIRED matrix", Layout, value = TRUE)
index = grep("GSE|and|PAIRED matrix", Layout)#, value = TRUE)
length(index) #47
RMPeak_datasets <- RMPeak_datasets[index,]

############
## RMdeg ##
############
RMdeg_list <- vector("list", 120)
names(RMdeg_list) <- RMPdeg_datasets$GSE
RMdem_list <- vector("list", 120)
names(RMdem_list) <- RMPdeg_datasets$GSE
RMdep_list <- vector("list", 120)
names(RMdep_list) <- RMPdeg_datasets$GSE

contrast_input_list <- vector("list", 120)
names(contrast_input_list) <- RMPdeg_datasets$GSE
contrast_IP_list <- vector("list", 120)
names(contrast_IP_list) <- RMPdeg_datasets$GSE

condition_input_list <- vector("list", 120)
names(condition_input_list) <- RMPdeg_datasets$GSE
condition_IP_list <- vector("list", 120)
names(condition_IP_list) <- RMPdeg_datasets$GSE

extract_names <- function(df) {
  colnames(df)
}

############
n = 1
############
#GSE102113# with RIP and Input
############
print(RMPdeg_datasets[n,])
i = 'GSE102113'
cl = 'Hela'
RMp = 'NAT10'
RM = 'ac4C'
condition = c(rep('NAT10KO', 2), rep('WT', 2))
# View(allDat[[i]])
ENA <- as.data.frame(RMP_allDat[grepl('GSE102113', RMP_allDat$study_alias ) & ( !grepl('IgG|RNA-Seq|Ribo-seq|BRICseq', RMP_allDat$sample_title))
                                ,c("study_alias", "sample_alias", "sample_title")])
ENA <- ENA[order(ENA$sample_alias),]
ENA
skip = ENA
paste0(skip$sample_alias, collapse = "|")

skip <- RMP_allDat[grepl('GSE102113', RMP_allDat$study_alias ) & ( !grepl('IgG|RNA-Seq|Ribo-seq|BRICseq', RMP_allDat$sample_title))
                   &  grepl('Input', RMP_allDat$sample_title)
                   , "sample_alias"]
paste0(skip$sample_alias, collapse = " ")
skip <- RMP_allDat[grepl('GSE102113', RMP_allDat$study_alias ) & ( !grepl('IgG|RNA-Seq|Ribo-seq|BRICseq', RMP_allDat$sample_title))
                   &  grepl('RIP', RMP_allDat$sample_title)
                   , "sample_alias"]
paste0(skip$sample_alias, collapse = " ")

files <- list.files(file.path(data_directory,i))
files
setwd(data_directory)
setwd(i)
# ###################################
# #         RNA-seq analysis        #
# ###################################
# df_in <- read.csv(paste0(data_directory,i, "/GSE102113_RNA-Seq_variance_stabilized_data.csv"),comment="#")
# # rownames(df_in) <- df_in$X
# DEG <- read_excel(file.path(data_directory,i '/DEG-',i,'.xlsx'),sheet = 2)
# summary(10^(-DEM$log10q))
# colnames(DEG)[1] <- 'geneSymbol'
# table(DEM$geneSymbol %in% DEP$geneSymbol)
# table(is.na(DEG$padj))
# DEG$padj <- as.numeric(DEG$padj)
# # DEG <- DEG[DEG$pvalue < 0.05,]
# table(DEG$padj < 0.05)
# length(which(DEG$padj < 0.05 & abs(DEG$log2FoldChange) > log2(1.5))) #how many are significant?
# DEG$label <- with(DEG, ifelse(DEG$baseMean>10 & abs(DEG$log2FoldChange)>log2(1.5) & DEG$padj<0.05,"sig","nosig"))
# DEG$change <- with(DEG, ifelse(DEG$padj < 0.05 & abs(DEG$log2FoldChange) > log2(1.5),
#                                ifelse(DEG$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
# table(DEG$change)
# DEG = as.data.frame(DEG[DEG$change %in% c('UP', 'DOWN'), ])
# resInP = DEG
# ###################################
# #         RIP-seq analysis        #
# ###################################
# peaks <- read.table(file.path(data_directory,i,'GSE102113_acRIP_peaks.csv'),sep = ',', header = T)
# # peak score (MACS)
# # Fold change (MACS)
# # -log10 p-value (MACS) 
# # -log10 q-value (FDR adjustment, MACS) 
# summary(peaks$log10p)
# summary(10^(-peaks$log10p))
# summary(10^(-peaks$log10q))
# DEP <- read_excel(file.path(data_directory,i, 'DEP-',i,'.xlsx'),sheet = 2)
# dim(DEP) #pooled peaks
# # [1] 4250   31
# DEM <- read_excel(file.path(data_directory,i, 'DEM-',i,'.xlsx'),sheet = 2)#  ac4C(+) Targets
# dim(DEM)
# # [1] 2135   24
# resMeR <- DEM
# 
# setdiff(unique(DEP$geneSymbol), DEM$geneSymbol)
# length(unique(DEP$geneSymbol)) #2139
# table(DEM$peak_id %in% DEP$peak_id)
# df_peak <- DEP[, c('KO.t','KO.c','WT.t','WT.c')]
# #判断数据是否需要转换
# ex <- df_peak#df_in[,-1]
# qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
# LogC <- (qx[5] > 100) ||
#   (qx[6]-qx[1] > 50 && qx[2] > 0) ||
#   (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
# if (LogC) { ex[which(ex <= 0)] <- NaN
# print("log2 transform unfinished")}else{print("log2 transform not needed")}
# 
# df_peak <- voom(df_peak)
# df_peak <- as.data.frame(df_peak$E)
# df_peak$peak_id <- DEP$peak_id
# Norm.RIP <- left_join(DEP[, c('peak_id','chr', 'start', 'end','geneSymbol')], df_peak, by="peak_id")
# colnames(df_peak)
# colnames(Norm.RIP) <- c('peak_id','chr', 'start', 'end','geneSymbol', paste0(c("RIP",'Input'), rep(1:2, each = 2)))
# Norm.RIP$Ave.input <- log2(rowMeans(2^Norm.RIP[, grepl('Input', colnames(Norm.RIP))]))
# # RIP1/Input1: KO; RIP2/Input2: WT
# # Adjusting RIP values for gene expression
# # All RIP values are corrected by dividing by the ratio (Input.sample/Input.average) (note: all values are in log2 format so this appears as "-Input#+Ave.Input")
# Norm.RIP <- mutate(Norm.RIP,
#                    Ratio1 = RIP1-Input1+Ave.input,
#                    Ratio2 = RIP2-Input2+Ave.input,
# )
# # Formatting the "Ratio" table (i.e. the adjusted RIP peaks used for the final analysis)
# # No filtering is performed as both input and RIP data are already filtered for low counts
# Ratio <- cbind(Norm.RIP[,c('chr', 'start', 'end','geneSymbol')], 2^(Norm.RIP[,grepl('Ratio', colnames(Norm.RIP))]) )
# Ratio$RefID <- paste("Ref", 1:nrow(Ratio)) # used to bring back peak information later
# rownames(Ratio) <- Ratio$RefID
# ex <- Ratio[, grepl('Ratio', colnames(Ratio))]
# qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
# LogC <- (qx[5] > 100) ||
#   (qx[6]-qx[1] > 50 && qx[2] > 0) ||
#   (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
# if (LogC) { ex[which(ex <= 0)] <- NaN
# print("log2 transform unfinished")}else{print("log2 transform not needed")}
# # 2. Differential gene analysis with edgeR
# # Without replicates, you cannot estimate which genes are differentially expressed using EdgeR or DESeq2. You can only calculate fold changes based on normalized read counts (preferentially CPM normalized by TMM method included in any EdgeR analysis) and apply a stringent fold-change cut-off to determine which genes are more or less expressed depending on the condition.
# y <- DGEList(ex)
# # Calculate fold change
# # norm.factor <- calcNormFactors(y, method = c('TMM'))
# # lib.size <- colSums(ex)
# # final.factor <- norm.factor*lib.size
# # perMillion.factor <- (final.factor/1000000)*-1
# y <- calcNormFactors(y, method = c('TMM'))
# y$samples$lib.size <- colSums(y$counts)
# logcpm <- cpm(y, prior.count=2, log=TRUE)
# logcpm[1:5,]
# cpm <- cpm(y, prior.count=2, log=FALSE)
# cpm[1:5, ]
# # Calculate fold changes between conditions (e.g., condition A vs. condition B)
# # fold_change <- logcpm[,'Ratio1']/logcpm[,'Ratio2']
# fold_change <- cpm[,'Ratio1']/cpm[,'Ratio2']
# Ratio$logFC <- log2(fold_change)
# Ratio$FoldChange <- fold_change
# All <- Ratio
# All$change <- with(All, ifelse(abs(All$logFC) > log2(1.5),
#                  ifelse(All$logFC > log2(1.5),'UP','DOWN'),'NOT'))
# table(All$change)
# target_res <- All[All$change %in% c('UP', 'DOWN'), ]
# length(unique(target_res$geneSymbol))
# target_res <- target_res[!duplicated(target_res$geneSymbol),]
# table(target_res$geneSymbol %in% resInP$geneSymbol)
# table(target_res$geneSymbol %in% resMeR$geneSymbol)
# target_res <- target_res[target_res$geneSymbol %in% resInP$geneSymbol, ]
# target_res <- target_res[target_res$geneSymbol %in% resMeR$geneSymbol, ]

###################################
#  Cell number-normalized (CNN)   #
#      MeRIP INPUT analysis       #
###################################
c <- unique(condition[!condition %in% 'WT'])
countGTF <- read.table( paste0("metaPlotR/Input_counts_", c, ".txt"),header=T)
dim(countGTF)
colnames(countGTF)
rownames(countGTF) <- countGTF$Geneid
countGTF <- countGTF[, !colnames(countGTF) %in% c('Geneid','Chr','Start','End','Strand','Length')]
dim(countGTF) #28265     4
## Sample name
## Sample information
sample=data.frame(
  sample=colnames(countGTF),
  condition=condition)
## SummarizedExperiment input
dds<-DESeqDataSetFromMatrix(countData=countGTF,colData=sample,design=~condition)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 10
table(keep)
dds <- dds[keep,]
## Note on factor levels
dds$condition <- factor(dds$condition, levels = unique(condition))
dds$condition <- relevel(dds$condition, ref = "WT")
# Differential expression analysis-------------------------------------------
dds <- DESeq(dds)

res <-  results(dds, contrast=c("condition",unique(condition)))
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
table(res$change)
table(is.na(res$change))
# resOrdered <- as.data.frame(res[order(res$padj), ])
resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
resInP$geneSymbol <- rownames(resInP)
name = res@elementMetadata$description[2]
name = paste0(cl, '_',str_split(name, "condition ", simplify = T)[,2])
RMdeg_list[[i]] <- list( resInP)
names(RMdeg_list[[i]] ) <- paste0(gsub(' ','_',(name)))
countGTF_list <- list(colnames(countGTF)); names(countGTF_list) <- c
condition_list <- list(condition); names(condition_list) <- c
contrast_input_list[[i]] <- countGTF_list
condition_input_list[[i]] <- condition_list

###################################
#  Cell number-normalized (CNN)   #
#   MeRIP RIP analysis - Part1    #
###################################
## Sample information
merge_peak <- read.delim(paste0("metaPlotR/", c, "_newPeakFile.txt"),header=F)
countMatrix <- read.table(paste0("metaPlotR/ac4C_counts_", c, ".txt"),header=T)
table(countMatrix$Geneid == merge_peak[-1,1])
table(duplicated(merge_peak$V1[-1]))
dim(countMatrix)
colnames(countMatrix)
rownames(countMatrix) <- countMatrix$Geneid
countMatrix <- countMatrix[, !colnames(countMatrix) %in% c('Geneid','Chr','Start','End','Strand','Length')]
# 1. 做好分组因子即可
groupList<-factor(condition)
colData <- data.frame(row.names=colnames(countMatrix), groupList=groupList)
# 2. Exploring Method 1: DESeq2-------------------------------------------
# 构建 DESeqDataSet 对象
dds <- DESeqDataSetFromMatrix(countData = countMatrix,
                              colData = colData,
                              design = ~ groupList)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 1
table(keep)
dds <- dds[keep,]
# 接着计算每个样本的 library size 用于后续的标准化。library size 等于落在该样本所有peaks上的reads的总数，即 counts 矩阵中每列的加和。如果bFullLibrarySize设为TRUE，则会使用library的总reads数（根据原BAM/BED文件统计获得）。然后使用estimateDispersions估计统计分布，需要将参数fitType设为’local’
dds<-estimateSizeFactors(dds)
dds<-estimateDispersions(dds,fitType="local")
dds <- nbinomWaldTest(dds)
res <- results(dds)
## 提取你想要的差异分析结果，我们这里是treated组对untreated组进行比较
res <-  results(dds, contrast=c("groupList", unique(condition)))
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
res_DESeq2 = as.data.frame(res)
table(res_DESeq2$change)#Not only
res_DESeq2 = res_DESeq2[res_DESeq2$change %in% c('UP', 'DOWN'), ]
# 2. Exploring Method 2: limma-------------------------------------------
library(edgeR)
d0 <- DGEList(counts = countMatrix)
dim(d0)
# Filtering low count peaks
cutoff <- 1           
drop <- which(apply(cpm(d0), 1, max) < cutoff)
print(length(drop))
if (length(drop) >0) {
  d <- d0[-drop,] 
} else{ d <- d0}
dim(d) #number of peaks left: 8501      
d <- calcNormFactors(d, method = c('TMM'))
colnames(countMatrix)
genotype <- groupList
table(genotype)
mm <- model.matrix(~0 + genotype)#
v <- voom(d, mm, plot=F, normalize="quantile")#voom from limma:Transform count data to log2-counts per million (logtpm)
fit <- lmFit(v, mm)
head(coef(fit))
contr <- makeContrasts(paste0(colnames(coef(fit)), collapse  = ' - '), levels = colnames(coef(fit))) #Change if other comparison
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)
length(which(top.table$adj.P.Val < 0.05 & abs(top.table$logFC) > log2(1.5))) #how many are significant?
res_limma <- as.data.frame(top.table)
res_limma$change <- with(res_limma, ifelse(res_limma$adj.P.Val < 0.05 & abs(res_limma$logFC) > log2(1.5),
                                           ifelse(res_limma$logFC > log2(1.5),'UP','DOWN'),'NOT'))
table(res_limma$change)#NOT only
res_limma = res_limma[res_limma$change %in% c('UP', 'DOWN'), ]

RMdem_list[[i]] <- list(res_DESeq2, res_limma)
name = res@elementMetadata$description[2]
name = paste0(cl, '_',str_split(name, "groupList ", simplify = T)[,2])
name= paste0(gsub(' ','_',(name)))
names(RMdem_list[[i]] ) <- paste0(name, '_', c('DESeq2', 'limma'))

RMdep_list[[i]] <- list( merge_peak )
names(RMdep_list[[i]] ) <- name

countMatrix_list <- list(colnames(countMatrix)); names(countMatrix_list) <- c
condition_list <- list(condition); names(condition_list) <- c
contrast_IP_list[[i]] <- countMatrix_list
condition_IP_list[[i]] <- condition_list

############
n = 2
############
#GSE90963#  with RIP and Input
############
print(RMPdeg_datasets[n,])
i = RMPdeg_datasets[n,'GSE']
cl = 'Hela'
RMp = 'DKC1'
RM = 'Psi'
condition=c(rep("DKC1KD",2),rep("WT",2))
# View(allDat[[i]])
ENA <- as.data.frame(RMP_allDat[grepl('PRJNA355164', RMP_allDat$study_alias ) & grepl('DKC1|Ctrl', RMP_allDat$sample_alias)
                                ,c("study_alias", "sample_alias", "sample_title")])
ENA[order(ENA$sample_alias),]

skip = ENA
paste0(skip$sample_alias, collapse = "|")

skip <- RMP_allDat[grepl('PRJNA355164', RMP_allDat$study_alias ) & grepl('DKC1|Ctrl', RMP_allDat$sample_alias)
                   &  grepl('NBS', RMP_allDat$sample_alias)
                   , "sample_alias"]
paste0(skip$sample_alias, collapse = " ")
skip <- RMP_allDat[grepl('PRJNA355164', RMP_allDat$study_alias ) & grepl('DKC1|Ctrl', RMP_allDat$sample_alias)
                   &  (! grepl('NBS', RMP_allDat$sample_alias))
                   , "sample_alias"]
paste0(skip$sample_alias, collapse = " ")

# Pseudouridine positions were called if the following criteria were met: bisulfite (BS) ≥5 deletions, BS fraction deletion ≥0.01, and BS coverage ≥10 reads. Adjacent positions were merged into individual deletion groups, which were further pruned by removing positions with fraction deletions less than half the maximum observed faction deletion in the group
DEM <- read_excel(file.path(data_directory,i, '/GSE90963_Table_S8-PseudoU_candidate_site_validation.xlsx'),sheet = 2)#
DEM <- read_excel(file.path(data_directory,i, '/GSE90963_Table_S5-PseudoU_candidate_sites.xlsx'),sheet = 2)# 
DEM <- read_excel(file.path(data_directory,i, '/pnas.1817334116.sd09.xlsx'),sheet = 2)# 
# Comparison of the DKC1-siRNA with control-siRNA datasets revealed significant reduction (>25% reduction, FDR < 0.01) of the deletion signature levels in 227 sites; most reside within rRNAs, although 18 sites were observed within mRNAs 
table(DEM$NBSPval...42<0.05)
target_res <- DEM[DEM$FDR < 0.05,]
setwd(data_directory)
setwd(i)
###################################
#  Cell number-normalized (CNN)   #
#      MeRIP INPUT analysis       #
###################################
condition=c(rep("DKC1KD",1),rep("WT",1))
c <- unique(condition[!condition %in% 'WT'])
condition <- unique(condition)
countGTF <- read.table( paste0("metaPlotR/Input_counts_", c, ".txt"),header=T)
dim(countGTF)
colnames(countGTF)
rownames(countGTF) <- countGTF$Geneid
countGTF <- countGTF[, !colnames(countGTF) %in% c('Geneid','Chr','Start','End','Strand','Length')]
dim(countGTF) # 28265     4
# 2. Differential gene analysis with edgeR
# Without replicates, you cannot estimate which genes are differentially expressed using EdgeR or DESeq2. You can only calculate fold changes based on normalized read counts (preferentially CPM normalized by TMM method included in any EdgeR analysis) and apply a stringent fold-change cut-off to determine which genes are more or less expressed depending on the condition.
y <- DGEList(countGTF)
# Calculate fold change
y <- calcNormFactors(y, method = c('TMM'))
y$samples$lib.size <- colSums(y$counts)
logcpm <- cpm(y, prior.count=2, log=TRUE)
logcpm[1:5,]
cpm <- cpm(y, prior.count=2, log=FALSE)
cpm[1:5, ]
fold_change <- cpm[,1]/cpm[,2]
res <- countGTF
res$logFC <- log2(fold_change)
res$FoldChange <- fold_change
res$change <- with(res, ifelse(abs(res$logFC) > log2(1.5),
                               ifelse(res$logFC > log2(1.5),'UP','DOWN'),'NOT'))
table(res$change)
resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
resInP$geneSymbol <- rownames(resInP)
table(resInP$change)
table(resInP$FoldChange %in% c(0, Inf))
table(resInP$logFC %in% c(-Inf, Inf))
# resInP <- resInP[! resInP$FoldChange %in% c(0, Inf), ]
RMdeg_list[[i]] <- list( resInP)
name = "Hela_DKC1KD_vs_WT"
RMdeg_list[[i]] <- list( resInP)
names(RMdeg_list[[i]] ) <- paste0(gsub(' ','_',(name)))

countGTF_list <- list(colnames(countGTF)); names(countGTF_list) <- c
condition_list <- list(condition); names(condition_list) <- c
contrast_input_list[[i]] <- countGTF_list
condition_input_list[[i]] <- condition_list

###################################
#  Cell number-normalized (CNN)   #
#   MeRIP RIP analysis - Part1    #
###################################
## Sample information
# # 1) WT
# macs2 callpeak -t $outdir/PRJNA355164/Ctrl_PolyA_BS_IP.bam -c $outdir/PRJNA355164/Ctrl_Ribo_NBS_input.bam --nomodel --extsize 200 -f BAM -q 0.01 -n WT_rep1 --outdir $outfolder
# macs2 callpeak -t $outdir/PRJNA355164/Ctrl_Ribo_BS_IP.bam -c $outdir/PRJNA355164/Ctrl_Ribo_NBS_input.bam --nomodel --extsize 200 -f BAM -q 0.01 -n WT_rep2 --outdir $outfolder
# # 2) DKC1KD
# macs2 callpeak -t $outdir/PRJNA355164/DKC1_PolyA_BS_IP.bam -c $outdir/PRJNA355164/DKC1_Ribo_NBS_input.bam --nomodel --extsize 200 -f BAM -q 0.01 -n DKC1KD_rep1 --outdir $outfolder
# macs2 callpeak -t $outdir/PRJNA355164/DKC1_Ribo_BS_IP.bam -c $outdir/PRJNA355164/DKC1_Ribo_NBS_input.bam --nomodel --extsize 200 -f BAM -q 0.01 -n DKC1KD_rep2 --outdir $outfolder
condition=c(rep("DKC1KD",1),rep("WT",1))
c <- unique(condition[!condition %in% 'WT'])
merge_peak <- read.delim(paste0("metaPlotR/", c, "_newPeakFile.txt"),header=F)
countMatrix <- read.table(paste0("metaPlotR/Psi_counts_", c, ".txt"),header=T)
table(countMatrix$Geneid == merge_peak[-1,1])
table(duplicated(merge_peak$V1[-1]))
dim(countMatrix)
colnames(countMatrix)
rownames(countMatrix) <- countMatrix$Geneid
countMatrix <- countMatrix[, !colnames(countMatrix) %in% c('Geneid','Chr','Start','End','Strand','Length')]
# countMatrix <- countMatrix[, ! grepl('Ribo', colnames(countMatrix))]
countMatrix <- countMatrix[, ! grepl('PolyA', colnames(countMatrix))]
# 2. Differential gene analysis with edgeR
# Without replicates, you cannot estimate which genes are differentially expressed using EdgeR or DESeq2. You can only calculate fold changes based on normalized read counts (preferentially CPM normalized by TMM method included in any EdgeR analysis) and apply a stringent fold-change cut-off to determine which genes are more or less expressed depending on the condition.
y <- DGEList(countMatrix)
# Calculate fold change
y <- calcNormFactors(y, method = c('TMM'))
y$samples$lib.size <- colSums(y$counts)
logcpm <- cpm(y, prior.count=2, log=TRUE)
logcpm[1:5,]
cpm <- cpm(y, prior.count=2, log=FALSE)
cpm[1:5, ]
# fold_change <- logcpm[,'Ratio1']/logcpm[,'Ratio2']
fold_change <- cpm[,1]/cpm[,2]
resMeR <- countMatrix
resMeR$logFC <- log2(fold_change)
resMeR$FoldChange <- fold_change
resMeR$change <- with(resMeR, ifelse(abs(resMeR$logFC) > log2(1.5),
                                     ifelse(resMeR$logFC > log2(1.5),'UP','DOWN'),'NOT'))
table(resMeR$change)
res_limma <- resMeR[resMeR$change %in% c('UP','DOWN'),]
RMdem_list[[i]] <- list(res_limma)
name = paste0(c, ' vs WT')
name = paste0(gsub(' ','_', name ))
names(RMdem_list[[i]]) <- paste0(cl, '_', name, '_limma')
RMdep_list[[i]] <- list( merge_peak )
names(RMdep_list[[i]] ) <-  paste0(cl, '_', name)

countMatrix_list <- list(colnames(countMatrix)); names(countMatrix_list) <- c
condition_list <- list(condition); names(condition_list) <- c
contrast_IP_list[[i]] <- countMatrix_list
condition_IP_list[[i]] <- condition_list

###################################
#  Cell number-normalized (CNN)   #
#   MeRIP RIP analysis - Part1    #
###################################
## Sample information
# condition=c(rep("DKC1KD",2),rep("WT",2))
# c <- unique(condition[!condition %in% 'WT'])
# merge_peak <- read.delim(paste0("metaPlotR/", c, "_newPeakFile.txt"),header=F)
# countMatrix <- read.table(paste0("metaPlotR/Psi_counts_", c, ".txt"),header=T)
# table(countMatrix$Geneid == merge_peak[-1,1])
# table(duplicated(merge_peak$V1[-1]))
# dim(countMatrix)
# colnames(countMatrix)
# rownames(countMatrix) <- countMatrix$Geneid
# countMatrix <- countMatrix[, !colnames(countMatrix) %in% c('Geneid','Chr','Start','End','Strand','Length')]
# # 1. 做好分组因子即可
# groupList<-factor(condition)
# colData <- data.frame(row.names=colnames(countMatrix), groupList=groupList)
# # 2. Exploring Method 1: DESeq2-------------------------------------------
# # 构建 DESeqDataSet 对象
# dds <- DESeqDataSetFromMatrix(countData = countMatrix,
#                               colData = colData,
#                               design = ~ groupList)
# ## Pre-filtering
# keep <- rowSums(counts(dds)) >= 1
# table(keep)
# dds <- dds[keep,]
# # 接着计算每个样本的 library size 用于后续的标准化。library size 等于落在该样本所有peaks上的reads的总数，即 counts 矩阵中每列的加和。如果bFullLibrarySize设为TRUE，则会使用library的总reads数（根据原BAM/BED文件统计获得）。然后使用estimateDispersions估计统计分布，需要将参数fitType设为’local’
# dds<-estimateSizeFactors(dds)
# dds<-estimateDispersions(dds,fitType="local")
# dds <- nbinomWaldTest(dds)
# res <- results(dds)
# ## 提取你想要的差异分析结果，我们这里是treated组对untreated组进行比较
# res <-  results(dds, contrast=c("groupList", unique(condition)))
# res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
# length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
# res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
#                                ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
# res_DESeq2 = as.data.frame(res)
# table(res_DESeq2$change)#Not only
# res_DESeq2 = res_DESeq2[res_DESeq2$change %in% c('UP', 'DOWN'), ]
# # 2. Exploring Method 2: limma-------------------------------------------
# library(edgeR)
# d0 <- DGEList(counts = countMatrix)
# dim(d0)
# # Filtering low count peaks
# cutoff <- 1           
# drop <- which(apply(cpm(d0), 1, max) < cutoff)
# print(length(drop))
# if (length(drop) >0) {
#   d <- d0[-drop,] 
# } else{ d <- d0}
# dim(d) #number of peaks left: 8501      
# d <- calcNormFactors(d, method = c('TMM'))
# colnames(countMatrix)
# genotype <- groupList
# table(genotype)
# mm <- model.matrix(~0 + genotype)#
# v <- voom(d, mm, plot=F, normalize="quantile")#voom from limma:Transform count data to log2-counts per million (logtpm)
# fit <- lmFit(v, mm)
# head(coef(fit))
# contr <- makeContrasts(paste0(colnames(coef(fit)), collapse  = ' - '), levels = colnames(coef(fit))) #Change if other comparison
# tmp <- contrasts.fit(fit, contr)
# tmp <- eBayes(tmp)
# top.table <- topTable(tmp, sort.by = "P", n = Inf)
# head(top.table, 20)
# length(which(top.table$adj.P.Val < 0.05 & abs(top.table$logFC) > log2(1.5))) #how many are significant?
# res_limma <- as.data.frame(top.table)
# res_limma$change <- with(res_limma, ifelse(res_limma$adj.P.Val < 0.05 & abs(res_limma$logFC) > log2(1.5),
#                                            ifelse(res_limma$logFC > log2(1.5),'UP','DOWN'),'NOT'))
# table(res_limma$change)#NOT only
# res_limma = res_limma[res_limma$change %in% c('UP', 'DOWN'), ]
# RMdem_list[[i]] <- list(res_DESeq2, res_limma)
# name = res@elementMetadata$description[2]
# name = paste0(cl, '_',str_split(name, "groupList ", simplify = T)[,2])
# name= paste0(gsub(' ','_',(name)))
# names(RMdem_list[[i]] ) <- paste0(name, '_', c('DESeq2', 'limma'))
# RMdep_list[[i]] <- list( merge_peak )
# names(RMdep_list[[i]] ) <- name
# 
# countMatrix_list <- list(colnames(countMatrix)); names(countMatrix_list) <- c
# condition_list <- list(condition); names(condition_list) <- c
# contrast_IP_list[[i]] <- countMatrix_list
# condition_IP_list[[i]] <- condition_list


############
n = 3
############
#GSE97419# with RIP and Input
############
print(RMPdeg_datasets[n,])
i = RMPdeg_datasets[n,'GSE']
cl = 'HEK293T'
RMp = 'TRMT61A'
RM = 'm1A'
condition=c(rep("TRMT61AOE",2),rep("WT",2))

# View(allDat[[i]])
files <- list.files(file.path(data_directory,i))
files
setwd(data_directory)
setwd(i)
###################################
#  Cell number-normalized (CNN)   #
#      MeRIP INPUT analysis       #
###################################
c <- unique(condition[!condition %in% 'WT'])
count_in <- read.table("metaPlotR/Input_counts.txt",header=T)
countGTF <- count_in[, !colnames(count_in) %in% c('Geneid','Chr','Start','End','Strand','Length')]
rownames(countGTF) <- count_in$Geneid
dim(countGTF) #28265     4
## Sample name
## Sample information
sample=data.frame(
  sample=colnames(countGTF),
  condition=condition)
## SummarizedExperiment input
dds<-DESeqDataSetFromMatrix(countData=countGTF,colData=sample,design=~condition)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 10
table(keep)
dds <- dds[keep,]
## Note on factor levels
dds$condition <- factor(dds$condition, levels = unique(condition))
dds$condition <- relevel(dds$condition, ref = "WT")
# Differential expression analysis-------------------------------------------
dds <- DESeq(dds)
# res<-results(dds,name = "condition_TRMT61AOE_vs_WT")
res <-  results(dds, contrast=c("condition",unique(condition)))
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
table(res$change)
table(is.na(res$change))
# resOrdered <- as.data.frame(res[order(res$padj), ])
resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
resInP$geneSymbol <- rownames(resInP)
name = res@elementMetadata$description[2]
name = paste0(cl, '_',str_split(name, "condition ", simplify = T)[,2])
RMdeg_list[[i]] <- list( resInP)
names(RMdeg_list[[i]] ) <- paste0(gsub(' ','_',(name)))

countGTF_list <- list(colnames(countGTF)); names(countGTF_list) <- c
condition_list <- list(condition); names(condition_list) <- c
contrast_input_list[[i]] <- countGTF_list
condition_input_list[[i]] <- condition_list
###################################
#  Cell number-normalized (CNN)   #
#   MeRIP RIP analysis - Part1    #
###################################
# count_table<-read.delim("metaPlotR/counts_multicov.txt",header=F)
# countMatrix<-as.matrix(count_table[,-c(1:6)])
# rownames(countMatrix)<-count_table$V4
## Sample information
merge_peak <- read.delim("metaPlotR/newPeakFile.txt",header=F)
countMatrix <- read.table("metaPlotR/m1A_counts.txt",header=T)
table(countMatrix$Geneid == merge_peak$V1[-1])
table(duplicated(merge_peak$V1[-1]))
dim(countMatrix)
colnames(countMatrix)
rownames(countMatrix) <- countMatrix$Geneid
countMatrix <- countMatrix[, !colnames(countMatrix) %in% c('Geneid','Chr','Start','End','Strand','Length')]
# 1. 做好分组因子即可
groupList<-factor(condition)
colData <- data.frame(row.names=colnames(countMatrix), groupList=groupList)
# 2. Exploring Method 1: DESeq2-------------------------------------------
# 构建 DESeqDataSet 对象
dds <- DESeqDataSetFromMatrix(countData = countMatrix,
                              colData = colData,
                              design = ~ groupList)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 1
# keep <- rowSums(counts(dds)>1) >= 1
table(keep)
dds <- dds[keep,]
# 接着计算每个样本的 library size 用于后续的标准化。library size 等于落在该样本所有peaks上的reads的总数，即 counts 矩阵中每列的加和。如果bFullLibrarySize设为TRUE，则会使用library的总reads数（根据原BAM/BED文件统计获得）。然后使用estimateDispersions估计统计分布，需要将参数fitType设为’local’
dds<-estimateSizeFactors(dds)
dds<-estimateDispersions(dds,fitType="local")
dds <- nbinomWaldTest(dds)
res <- results(dds)
## 提取你想要的差异分析结果，我们这里是treated组对untreated组进行比较
res <-  results(dds, contrast=c("groupList", unique(condition)))
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
# resOrdered <- res[order(res$padj), ]
res_DESeq2 = as.data.frame(res)
table(res_DESeq2$change)
res_DESeq2 = res_DESeq2[res_DESeq2$change %in% c('UP', 'DOWN'), ]
# 2. Exploring Method 2: limma-------------------------------------------
library(edgeR)
d0 <- DGEList(counts = countMatrix)
dim(d0)
# Filtering low count peaks
cutoff <- 1           
drop <- which(apply(cpm(d0), 1, max) < cutoff)
print(length(drop))
if (length(drop) >0) {
  d <- d0[-drop,] 
} else{ d <- d0}
dim(d) #number of peaks left: 8501      
d <- calcNormFactors(d, method = c('TMM'))
# d$samples$lib.size <- colSums(d$counts)
# N <- colSums(countMatrix)
# nf <- calcNormFactors(ercc_expr, lib.size=N) #calculating normalizing factors (nf) based on ERCCs
colnames(countMatrix)
genotype <- groupList
table(genotype)
mm <- model.matrix(~0 + genotype)#
v <- voom(d, mm, plot=F, normalize="quantile")#voom from limma:Transform count data to log2-counts per million (logtpm)
fit <- lmFit(v, mm)
head(coef(fit))
contr <- makeContrasts(paste0(colnames(coef(fit)), collapse  = ' - '), levels = colnames(coef(fit))) #Change if other comparison
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)
length(which(top.table$adj.P.Val < 0.05 & abs(top.table$logFC) > log2(1.5))) #how many are significant?
res_limma <- as.data.frame(top.table)
res_limma$change <- with(res_limma, ifelse(res_limma$adj.P.Val < 0.05 & abs(res_limma$logFC) > log2(1.5),
                                           ifelse(res_limma$logFC > log2(1.5),'UP','DOWN'),'NOT'))
table(res_limma$change)
res_limma = res_limma[res_limma$change %in% c('UP', 'DOWN'), ]

RMdem_list[[i]] <- list(res_DESeq2, res_limma)
name = res@elementMetadata$description[2]
name = paste0(cl, '_',str_split(name, "groupList ", simplify = T)[,2])
name= paste0(gsub(' ','_',(name)))
names(RMdem_list[[i]] ) <- paste0(name, '_', c('DESeq2', 'limma'))
RMdep_list[[i]] <- list( merge_peak )
names(RMdep_list[[i]] ) <- name

countMatrix_list <- list(colnames(countMatrix)); names(countMatrix_list) <- c
condition_list <- list(condition); names(condition_list) <- c
contrast_IP_list[[i]] <- countMatrix_list
condition_IP_list[[i]] <- condition_list

# ###################################
# #  Cell number-normalized (CNN)   #
# #   MeRIP RIP analysis - Part2    #
# ###################################
# library(tidyverse)
# library(edgeR)
# # 1. Formatting tables as needed for downstream analysis
# d0 <- DGEList(counts = countGTF)
# # Filtering low count genes
# cutoff <- 1           
# drop <- which(apply(cpm(d0), 1, max) < cutoff)
# length(drop)
# d <- d0[-drop,] 
# d <- calcNormFactors(d)
# N <- colSums(countGTF)
# colnames(countGTF)
# genotype <- condition %>% as.factor()
# mm <- model.matrix(~0 + genotype)
# v <- voom(d, mm, plot=T, normalize="quantile")#voom from limma:Transform count data to log2-counts per million (logtpm)
# Input <- as.data.frame(v$E)
# Input$geneSymbol <- rownames(Input)
# 
# d0 <- DGEList(countMatrix)
# N <- colSums(countMatrix)  
# cutoff <- 1
# drop <- which(apply(cpm(d0), 1, max) < cutoff)
# length(drop)
# if(length(drop)>0){
#   d <- d0[-drop,] 
# }else{
#   d <- d0 
# }
# d <- calcNormFactors(d)
# genotype <- condition %>% as.factor()##c(rep(c("WT", "KO"), 6))
# mm <- model.matrix(~0 + genotype)#
# v <- voom(d, mm, plot=T, normalize="quantile")#voom from limma:Transform count data to log2-counts per million (logtpm)
# 
# RIP <- data.frame(v$E)
# RIP$peak_id <- rownames(RIP)
# RIP <- inner_join(RIP, merge_bed.anno, by=c("peak_id"))
# 
# Norm.RIP <- inner_join(RIP, Input, by="geneSymbol")
# Norm.RIP$Ave.input <- log2(rowMeans(2^Norm.RIP[,colnames(Norm.RIP) %in% colnames(countGTF)]))
# colnames(Norm.RIP) <- c(paste0("RIP", 1:length(condition)), 'peak_id','chr', 'start', 'end','geneSymbol',
#                         paste0("Input", 1:length(condition)), "Ave.input")
# # Adjusting RIP values for gene expression
# # All RIP values are corrected by dividing by the ratio (Input.sample/Input.average) (note: all values are in log2 format so this appears as "-Input#+Ave.Input")
# Norm.RIP <- mutate(Norm.RIP, 
#                    Ratio1 = RIP1-Input1+Ave.input,
#                    Ratio2 = RIP2-Input2+Ave.input,
#                    Ratio3 = RIP3-Input3+Ave.input,
#                    Ratio4 = RIP4-Input4+Ave.input,
# )
# # Formatting the "Ratio" table (i.e. the adjusted RIP peaks used for the final analysis)
# # No filtering is performed as both input and RIP data are already filtered for low counts
# Ratio <- cbind(Norm.RIP[,c('chr', 'start', 'end','geneSymbol')], 2^(Norm.RIP[,grepl('Ratio', colnames(Norm.RIP))]) )
# Ratio$RefID <- paste("Ref", 1:nrow(Ratio)) # used to bring back peak information later
# rownames(Ratio) <- Ratio$RefID
# # 2. Differential gene analysis with edgeR
# # Building the model for limma-voom
# # As the data is already fully normalized, further re-scaling is prevented by fixing the library size in the voom function
# d0 <- DGEList(Ratio[, grepl('Ratio', colnames(Ratio))])
# d <- calcNormFactors(d0)
# N <- colSums(Ratio[,grepl('Ratio', colnames(Ratio))])
# genotype <- condition %>% as.factor()
# mm <- model.matrix(~0 + genotype)
# y <- voom(d, mm, plot = T, lib.size = rep(10^6,length(condition)))  #lib.size fixed as 10^6 for all samples
# # v <- voom(d, mm, plot=T, normalize="quantile")#voom from limma:Transform count data to log2-counts per million (logtpm)
# fit <- lmFit(y, mm)
# head(coef(fit))
# # Differential analysis for inputs (Paused vs FBS)
# contr <- makeContrasts(genotypeTRMT61AOE-genotypeWT, levels = colnames(coef(fit))) #Change if other comparison
# tmp <- contrasts.fit(fit, contr)
# tmp <- eBayes(tmp)
# top.table <- topTable(tmp, sort.by = "P", n = Inf)
# head(top.table, 20)
# length(which(top.table$adj.P.Val < 0.05 & abs(top.table$logFC) > log2(1.5))) #how many are significant?
# # Bring back peak information and export
# top.table$RefID <- rownames(top.table)
# top.table <- inner_join(top.table, Ratio, by=c("RefID"))
# # top.table <- top.table[,-c(7,18)]
# # write.table(top.table, file = paste0(export_folder,"/m6A.CNN.toptable.txt"))
# All <- inner_join(Norm.RIP, top.table[,!grepl('Ratio', colnames(top.table))], by=c('chr', 'start', 'end','geneSymbol'))
# # write.table(All, file = paste0(export_folder,"/m6A.CNN.all_results.txt"))
# All$change <- with(All, ifelse(All$adj.P.Val < 0.05 & abs(All$logFC) > log2(1.5),
#                                ifelse(All$logFC > log2(1.5),'UP','DOWN'),'NOT'))
# table(All$change)
# target_res <- All[All$change %in% c('UP', 'DOWN'), ]
# length(unique(target_res$geneSymbol))
# target_res <- target_res[!duplicated(target_res$geneSymbol),]
# table(target_res$geneSymbol %in% resMeR$geneSymbol)
# table(target_res$geneSymbol %in% resInP$geneSymbol)
# target_res <- target_res[target_res$geneSymbol %in% intersect(resMeR$geneSymbol, resInP$geneSymbol), ]

############
n = 4
############
#GSE133671# with RIP and RNA
############
print(RMPdeg_datasets[n,])
i = RMPdeg_datasets[n,'GSE']
cl = 'T24'
RMp = 'NSUN2'
RM = 'm5C'
condition = c(rep('NSUN2KD', 2), rep('WT', 2))
ENA <- as.data.frame(RMP_allDat[grepl('GSE133621|GSE133623', RMP_allDat$study_alias )
                                ,c("study_alias", "sample_alias", "sample_title")])
ENA[order(ENA$sample_alias),]
# View(allDat[[i]])
files <- list.files(file.path(data_directory,i))
files
setwd(data_directory)
setwd(i)
###################################
#      RNA/Ribo-seq analysis      #
###################################
condition = c(rep('NSUN2KD', 2), rep('WT', 2))
df_in <- read.table(file.path(data_directory,i, '/GSE133621_reads-count-cell-line.txt.gz'), header = T)
## Raw counts
#判断数据是否需要转换
ex <- df_in[,-1]
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
# Calculate LogC
LogC <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0) || (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
# Check if LogC is TRUE
if (LogC) {
  # Replace values in 'ex' that are less than or equal to 0 with NaN
  ex[ex <= 0] <- NaN
  print("log2 transform unfinished")
} else {
  print("log2 transform not needed")
}
##Gene count matrix
countGTF=df_in[,2:ncol(df_in)]
## Gene name
library(data.table)
id_mapped <- read.table("/data2/rluo4/hg38/gencode.v44lift37.annotation.gene.probeMap",header = T)
id_mapped <- read.table("/data2/rluo4/hg38/gencode.v36.annotation.gene.probeMap",header = T)
dim(id_mapped) #  60675     6
table(str_split(df_in$GeneID, '[.]', simplify = T)[,1] %in% gene_info$geneEnsembls)
table(df_in$GeneID %in% id_mapped$id)
table(str_split(df_in$GeneID, '[.]', simplify = T)[,1] %in% str_split(id_mapped$id, '[.]', simplify = T)[,1])
countGTF$id <- str_split(df_in$GeneID, '[.]', simplify = T)[,1]
id_mapped$id <- str_split(id_mapped$id, '[.]', simplify = T)[,1]
# id_mapped <- unique(id_mapped) #%>% dplyr::select(id,gene) # 去重
id_mapped <- id_mapped[which(id_mapped$id %in% countGTF$id),]
dim(id_mapped) # 51343     6 for v36 # 54667     6 for v44 
table(duplicated(id_mapped$gene))
countGTF <- left_join(id_mapped,countGTF,by="id")
countGTF <- aggregate(x = countGTF[,(ncol(id_mapped)+1):ncol(countGTF)],   #此时exprSet的第三列开始是表达矩阵内容
                      by = list(symbol = countGTF$gene),   #按照相同symbol分组，在组内计算
                      FUN = max) %>%   #原文中是计算最大值（max），也可以计算平均值（mean）或者中位数（median）
  column_to_rownames(var = 'symbol')
dim(countGTF) # 50460     4 for v36 # 45344     4 for v44
## Sample name
colnames(countGTF)
countGTF <- countGTF[, c('siNSUN2_rep1','siNSUN2_rep2', 'siControl.rep1','siControl.rep2')]
## Sample information
sample=data.frame(
  sample=colnames(countGTF),
  condition=condition)
## SummarizedExperiment input
dds<-DESeqDataSetFromMatrix(countData=countGTF,colData=sample,design=~condition)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 10
table(keep)
dds <- dds[keep,]
## Note on factor levels
dds$condition <- factor(dds$condition, levels = unique(condition))
dds$condition <- relevel(dds$condition, ref = "WT")
# Differential expression analysis-------------------------------------------
dds <- DESeq(dds)
# res<-results(dds,name = "condition_TRMT61AOE_vs_WT")
res <-  results(dds, contrast=c("condition",unique(condition)))
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
table(res$change)
table(is.na(res$change))
# resOrdered <- as.data.frame(res[order(res$padj), ])
resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
resInP$geneSymbol <- rownames(resInP)
name = res@elementMetadata$description[2]
name = paste0(cl, '_',str_split(name, "condition ", simplify = T)[,2])
RMdeg_list[[i]] <- list( resInP)
names(RMdeg_list[[i]] ) <- paste0(gsub(' ','_',(name)))
c <- unique(condition[!condition %in% 'WT'])
countGTF_list <- list(colnames(countGTF)); names(countGTF_list) <- c
condition_list <- list(condition); names(condition_list) <- c
contrast_input_list[[i]] <- countGTF_list
condition_input_list[[i]] <- condition_list

###################################
#  Cell number-normalized (CNN)   #
#   MeRIP RIP analysis - Part1    #
###################################
## Sample information
condition = c(rep('NSUN2KD', 1), rep('WT', 1))
merge_peak <- read.delim("metaPlotR/newPeakFile.txt",header=F)
countMatrix <- read.table("metaPlotR/m5C_counts.txt",header=T)
table(countMatrix$Geneid == merge_peak$V1[-1])
table(duplicated(merge_peak$V1[-1]))
dim(countMatrix)
colnames(countMatrix)
rownames(countMatrix) <- countMatrix$Geneid
countMatrix <- countMatrix[, !colnames(countMatrix) %in% c('Geneid','Chr','Start','End','Strand','Length')]
# 1. 做好分组因子即可
groupList<-factor(unique(condition))
ex <- countMatrix
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
# Calculate LogC
LogC <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0) || (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
# Check if LogC is TRUE
if (LogC) {
  # Replace values in 'ex' that are less than or equal to 0 with NaN
  ex[ex <= 0] <- NaN
  print("log2 transform unfinished")
} else {
  print("log2 transform not needed")
}
# 2. Differential gene analysis with edgeR
# Without replicates, you cannot estimate which genes are differentially expressed using EdgeR or DESeq2. You can only calculate fold changes based on normalized read counts (preferentially CPM normalized by TMM method included in any EdgeR analysis) and apply a stringent fold-change cut-off to determine which genes are more or less expressed depending on the condition.
y <- DGEList(countMatrix)
# Calculate fold change
y <- calcNormFactors(y, method = c('TMM'))
y$samples$lib.size <- colSums(y$counts)
logcpm <- cpm(y, prior.count=2, log=TRUE)
logcpm[1:5,]
cpm <- cpm(y, prior.count=2, log=FALSE)
cpm[1:5, ]
# fold_change <- logcpm[,'Ratio1']/logcpm[,'Ratio2']
fold_change <- cpm[,1]/cpm[,2]
resMeR <- countMatrix
resMeR$logFC <- log2(fold_change)
resMeR$FoldChange <- fold_change
resMeR$change <- with(resMeR, ifelse(abs(resMeR$logFC) > log2(1.5),
                                     ifelse(resMeR$logFC > log2(1.5),'UP','DOWN'),'NOT'))
table(resMeR$change)
res_limma <- resMeR[resMeR$change %in% c('UP','DOWN'),]
RMdem_list[[i]] <- list(res_limma)
name = paste0(c, ' vs WT')
name = paste0(gsub(' ','_', name ))
names(RMdem_list[[i]]) <- paste0(cl, '_', name, '_limma')
RMdep_list[[i]] <- list( merge_peak )
names(RMdep_list[[i]] ) <- paste0(cl, '_', name)

c <- unique(condition[!condition %in% 'WT'])
countMatrix_list <- list(colnames(countMatrix)); names(countMatrix_list) <- c
condition_list <- list(condition); names(condition_list) <- c
contrast_IP_list[[i]] <- countMatrix_list
condition_IP_list[[i]] <- condition_list

# ###################################
# #  Cell number-normalized (CNN)   #
# #   MeRIP RIP analysis - Part2    #
# ###################################
# library(tidyverse)
# library(edgeR)
# Input <- countGTF[, c('siNSUN2_rep1','siControl.rep1')]
# Input$geneSymbol <- rownames(Input)
# RIP <- left_join(df_peak, Input, by="geneSymbol")
# RIP <- na.omit(RIP)
# peak_id <- RIP$peak_id
# RIP <- RIP[, c(1:length(groupList), (ncol(RIP)-(length(groupList)-1)):ncol(RIP) )]
# RIP <- voom(RIP)
# RIP <- as.data.frame(RIP$E)
# RIP$peak_id <- peak_id
# Norm.RIP <- left_join(df_peak[, c('peak_id','chr', 'start', 'end','geneSymbol')], RIP, by="peak_id")
# Norm.RIP <- na.omit(Norm.RIP)
# colnames(Norm.RIP)
# colnames(Norm.RIP) <- c('peak_id','chr', 'start', 'end','geneSymbol', paste0(rep(c("RIP",'Input'), each = 2), rep(1:2, 2)))
#                         #paste0(c("RIP",'Input'), rep(1:2, each = 2)))
# Norm.RIP$Ave.input <- log2(rowMeans(2^Norm.RIP[, grepl('Input', colnames(Norm.RIP))]))
# # Adjusting RIP values for gene expression
# # All RIP values are corrected by dividing by the ratio (Input.sample/Input.average) (note: all values are in log2 format so this appears as "-Input#+Ave.Input")
# Norm.RIP <- mutate(Norm.RIP,
#                    Ratio1 = RIP1-Input1+Ave.input,
#                    Ratio2 = RIP2-Input2+Ave.input,
# )
# # Formatting the "Ratio" table (i.e. the adjusted RIP peaks used for the final analysis)
# # No filtering is performed as both input and RIP data are already filtered for low counts
# Ratio <- cbind(Norm.RIP[,c('chr', 'start', 'end','geneSymbol')], 2^(Norm.RIP[,grepl('Ratio', colnames(Norm.RIP))]) )
# Ratio$RefID <- paste("Ref", 1:nrow(Ratio)) # used to bring back peak information later
# rownames(Ratio) <- Ratio$RefID
# ex <- Ratio[, grepl('Ratio', colnames(Ratio))]
# qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
# LogC <- (qx[5] > 100) ||
#   (qx[6]-qx[1] > 50 && qx[2] > 0) ||
#   (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
# if (LogC) { ex[which(ex <= 0)] <- NaN
# print("log2 transform unfinished")}else{print("log2 transform not needed")}
# y <- DGEList(ex)
# y <- calcNormFactors(y, method = c('TMM'))
# y$samples$lib.size <- colSums(y$counts)
# logcpm <- cpm(y, prior.count=2, log=TRUE)
# logcpm[1:5,]
# cpm <- cpm(y, prior.count=2, log=FALSE)
# cpm[1:5, ]
# # Calculate fold changes between conditions (e.g., condition A vs. condition B)
# # fold_change <- logcpm[,'Ratio1']/logcpm[,'Ratio2']
# fold_change <- cpm[,'Ratio1']/cpm[,'Ratio2']
# Ratio$logFC <- log2(fold_change)
# Ratio$FoldChange <- fold_change
# All <- Ratio
# All$change <- with(All, ifelse(abs(All$logFC) > log2(1.5),
#                                ifelse(All$logFC > log2(1.5),'UP','DOWN'),'NOT'))
# table(All$change)
# target_res <- All[All$change %in% c('UP', 'DOWN'), ]
# length(unique(target_res$geneSymbol))# 140
# target_res <- target_res[!duplicated(target_res$geneSymbol),]
# table(target_res$geneSymbol %in% resMeR$geneSymbol)
# # FALSE  TRUE 
# # 31   109 
# table(target_res$geneSymbol %in% resInP$geneSymbol)
# # FALSE  TRUE 
# # 116    24 
# target_res <- target_res[target_res$geneSymbol %in% resInP$geneSymbol, ]
# target_res <- target_res[target_res$geneSymbol %in% resMeR$geneSymbol, ]
# # 17 -->15
############
n = 5
############
#GSE122254# with RIP without Input 
############
print(RMPdeg_datasets[n,])
i = RMPdeg_datasets[n,'GSE']
cl = 'Hela' #(use data from GSE93750/GSE93749 as Input)
RMp = 'NSUN2'
RM = 'm5C'
condition=c(rep("NSUN2KO",2),rep("WT",2))
# View(allDat[[i]])
ENA <- as.data.frame(RMP_allDat[grepl('GSE122254', RMP_allDat$study_alias ) & grepl('HeLa', RMP_allDat$sample_title)
                                ,c("study_alias", "sample_alias", "sample_title")])
ENA[order(ENA$sample_alias),]
files <- list.files(file.path(data_directory,i))
files
setwd(data_directory)
setwd(i)
library(readxl)
x <-  read_excel(file.path(data_directory,i, 'GSE122254_human_sites.xls'),sheet = 2)
###################################
#  Cell number-normalized (CNN)   #
#      MeRIP INPUT analysis       #
###################################
count_in <- read.table("metaPlotR/Input_counts_NSUN2KO.txt",header=T)
countGTF <- count_in[, !colnames(count_in) %in% c('Geneid','Chr','Start','End','Strand','Length')]
rownames(countGTF) <- count_in$Geneid
dim(countGTF) #56752     5
## Sample name
## Sample information
sample=data.frame(
  sample=colnames(countGTF),
  condition=condition)
## SummarizedExperiment input
dds<-DESeqDataSetFromMatrix(countData=countGTF,colData=sample,design=~condition)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 10
table(keep)
dds <- dds[keep,]
## Note on factor levels
dds$condition <- factor(dds$condition, levels = unique(condition))
dds$condition <- relevel(dds$condition, ref = "WT")
# Differential expression analysis-------------------------------------------
dds <- DESeq(dds)

res <-  results(dds, contrast=c("condition",unique(condition)))
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
table(res$change)
table(is.na(res$change))
# resOrdered <- as.data.frame(res[order(res$padj), ])
resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
resInP$geneSymbol <- rownames(resInP)
name = res@elementMetadata$description[2]
name = paste0(cl, '_',str_split(name, "condition ", simplify = T)[,2])
RMdeg_list[[i]] <- list( resInP)
names(RMdeg_list[[i]] ) <- paste0(gsub(' ','_',(name)))

c <- unique(condition[!condition %in% 'WT'])
countGTF_list <- list(colnames(countGTF)); names(countGTF_list) <- c
condition_list <- list(condition); names(condition_list) <- c
contrast_input_list[[i]] <- countGTF_list
condition_input_list[[i]] <- condition_list
###################################
#  Cell number-normalized (CNN)   #
#   MeRIP RIP analysis - Part1    #
###################################
## Sample information
merge_peak <- read.delim("metaPlotR/NSUN2KO_newPeakFile.txt",header=F)
countMatrix <- read.table("metaPlotR/m5C_counts_NSUN2KO.txt",header=T)
table(countMatrix$Geneid == merge_peak$V1[-1])
table(duplicated(merge_peak$V1[-1]))
dim(countMatrix)
colnames(countMatrix)
rownames(countMatrix) <- countMatrix$Geneid
countMatrix <- countMatrix[, !colnames(countMatrix) %in% c('Geneid','Chr','Start','End','Strand','Length')]
# 1. 做好分组因子即可
groupList<-factor(condition)
colData <- data.frame(row.names=colnames(countMatrix), groupList=groupList)
# 2. Exploring Method 1: DESeq2-------------------------------------------
# 构建 DESeqDataSet 对象
dds <- DESeqDataSetFromMatrix(countData = countMatrix,
                              colData = colData,
                              design = ~ groupList)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 1
table(keep)
dds <- dds[keep,]
# 接着计算每个样本的 library size 用于后续的标准化。library size 等于落在该样本所有peaks上的reads的总数，即 counts 矩阵中每列的加和。如果bFullLibrarySize设为TRUE，则会使用library的总reads数（根据原BAM/BED文件统计获得）。然后使用estimateDispersions估计统计分布，需要将参数fitType设为’local’
dds<-estimateSizeFactors(dds)
dds<-estimateDispersions(dds,fitType="local")
dds <- nbinomWaldTest(dds)
res <- results(dds)
## 提取你想要的差异分析结果，我们这里是treated组对untreated组进行比较
res <-  results(dds, contrast=c("groupList", unique(condition)))
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
res_DESeq2 = as.data.frame(res)
table(res_DESeq2$change)
res_DESeq2 = res_DESeq2[res_DESeq2$change %in% c('UP', 'DOWN'), ]
# 2. Exploring Method 2: limma-------------------------------------------
library(edgeR)
d0 <- DGEList(counts = countMatrix)
dim(d0)
# Filtering low count peaks
cutoff <- 1           
drop <- which(apply(cpm(d0), 1, max) < cutoff)
print(length(drop))
if (length(drop) >0) {
  d <- d0[-drop,] 
} else{ d <- d0}
dim(d) #number of peaks left: 8501      
d <- calcNormFactors(d, method = c('TMM'))
colnames(countMatrix)
genotype <- groupList
table(genotype)
mm <- model.matrix(~0 + genotype)#
v <- voom(d, mm, plot=F, normalize="quantile")#voom from limma:Transform count data to log2-counts per million (logtpm)
fit <- lmFit(v, mm)
head(coef(fit))
contr <- makeContrasts(paste0(colnames(coef(fit)), collapse  = ' - '), levels = colnames(coef(fit))) #Change if other comparison
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)
length(which(top.table$adj.P.Val < 0.05 & abs(top.table$logFC) > log2(1.5))) #how many are significant?
res_limma <- as.data.frame(top.table)
res_limma$change <- with(res_limma, ifelse(res_limma$adj.P.Val < 0.05 & abs(res_limma$logFC) > log2(1.5),
                                           ifelse(res_limma$logFC > log2(1.5),'UP','DOWN'),'NOT'))
table(res_limma$change)
res_limma = res_limma[res_limma$change %in% c('UP', 'DOWN'), ]

RMdem_list[[i]] <- list(res_DESeq2, res_limma)
name = res@elementMetadata$description[2]
name = paste0(cl, '_',str_split(name, "groupList ", simplify = T)[,2])
name= paste0(gsub(' ','_',(name)))
names(RMdem_list[[i]] ) <- paste0(name, '_', c('DESeq2', 'limma'))
RMdep_list[[i]] <- list( merge_peak )
names(RMdep_list[[i]] ) <- name

c <- unique(condition[!condition %in% 'WT'])
countMatrix_list <- list(colnames(countMatrix)); names(countMatrix_list) <- c
condition_list <- list(condition); names(condition_list) <- c
contrast_IP_list[[i]] <- countMatrix_list
condition_IP_list[[i]] <- condition_list
# ###################################
# #  Cell number-normalized (CNN)   #
# #   MeRIP RIP analysis - Part2    #
# ###################################
# library(tidyverse)
# library(edgeR)
# # 1. Formatting tables as needed for downstream analysis
# d0 <- DGEList(counts = countGTF)
# # Filtering low count genes
# cutoff <- 1           
# drop <- which(apply(cpm(d0), 1, max) < cutoff)
# length(drop)
# d <- d0[-drop,] 
# d <- calcNormFactors(d)
# N <- colSums(countGTF)
# colnames(countGTF)
# genotype <- condition %>% as.factor()
# mm <- model.matrix(~0 + genotype)
# v <- voom(d, mm, plot=T, normalize="quantile")#voom from limma:Transform count data to log2-counts per million (logtpm)
# Input <- as.data.frame(v$E)
# Input$geneSymbol <- rownames(Input)
# 
# d0 <- DGEList(countMatrix)
# N <- colSums(countMatrix)  
# cutoff <- 1
# drop <- which(apply(cpm(d0), 1, max) < cutoff)
# length(drop)
# if(length(drop)>0){
#   d <- d0[-drop,] 
# }else{
#   d <- d0 
# }
# d <- calcNormFactors(d)
# genotype <- condition %>% as.factor()##c(rep(c("WT", "KO"), 6))
# mm <- model.matrix(~0 + genotype)#
# v <- voom(d, mm, plot=T, normalize="quantile")#voom from limma:Transform count data to log2-counts per million (logtpm)
# 
# RIP <- data.frame(v$E)
# RIP$peak_id <- rownames(RIP)
# RIP <- inner_join(RIP, merge_bed.anno, by=c("peak_id"))
# 
# Norm.RIP <- inner_join(RIP, Input, by="geneSymbol")
# Norm.RIP$Ave.input <- log2(rowMeans(2^Norm.RIP[,colnames(Norm.RIP) %in% colnames(countGTF)]))
# colnames(Norm.RIP) <- c(paste0("RIP", 1:length(condition)), 'peak_id','chr', 'start', 'end','geneSymbol',
#                         paste0("Input", 1:length(condition)), "Ave.input")
# # Adjusting RIP values for gene expression
# # All RIP values are corrected by dividing by the ratio (Input.sample/Input.average) (note: all values are in log2 format so this appears as "-Input#+Ave.Input")
# Norm.RIP <- mutate(Norm.RIP, 
#                    Ratio1 = RIP1-Input1+Ave.input,
#                    Ratio2 = RIP2-Input2+Ave.input,
#                    Ratio3 = RIP3-Input3+Ave.input,
#                    Ratio4 = RIP4-Input4+Ave.input,
# )
# # Formatting the "Ratio" table (i.e. the adjusted RIP peaks used for the final analysis)
# # No filtering is performed as both input and RIP data are already filtered for low counts
# Ratio <- cbind(Norm.RIP[,c('chr', 'start', 'end','geneSymbol')], 2^(Norm.RIP[,grepl('Ratio', colnames(Norm.RIP))]) )
# Ratio$RefID <- paste("Ref", 1:nrow(Ratio)) # used to bring back peak information later
# rownames(Ratio) <- Ratio$RefID
# # 2. Differential gene analysis with edgeR
# # Building the model for limma-voom
# # As the data is already fully normalized, further re-scaling is prevented by fixing the library size in the voom function
# d0 <- DGEList(Ratio[, grepl('Ratio', colnames(Ratio))])
# d <- calcNormFactors(d0)
# N <- colSums(Ratio[,grepl('Ratio', colnames(Ratio))])
# genotype <- condition %>% as.factor()
# mm <- model.matrix(~0 + genotype)
# y <- voom(d, mm, plot = T, lib.size = rep(10^6,length(condition)))  #lib.size fixed as 10^6 for all samples
# # v <- voom(d, mm, plot=T, normalize="quantile")#voom from limma:Transform count data to log2-counts per million (logtpm)
# fit <- lmFit(y, mm)
# head(coef(fit))
# # Differential analysis for inputs (Paused vs FBS)
# contr <- makeContrasts(paste0(colnames(coef(fit)), collapse  = ' - '), levels = colnames(coef(fit))) #Change if other comparison
# tmp <- contrasts.fit(fit, contr)
# tmp <- eBayes(tmp)
# top.table <- topTable(tmp, sort.by = "P", n = Inf)
# head(top.table, 20)
# length(which(top.table$adj.P.Val < 0.05 & abs(top.table$logFC) > log2(1.5))) #how many are significant?
# # Bring back peak information and export
# top.table$RefID <- rownames(top.table)
# top.table <- inner_join(top.table, Ratio, by=c("RefID"))
# # top.table <- top.table[,-c(7,18)]
# # write.table(top.table, file = paste0(export_folder,"/m6A.CNN.toptable.txt"))
# All <- inner_join(Norm.RIP, top.table[,!grepl('Ratio', colnames(top.table))], by=c('chr', 'start', 'end','geneSymbol'))
# # write.table(All, file = paste0(export_folder,"/m6A.CNN.all_results.txt"))
# All$change <- with(All, ifelse(All$adj.P.Val < 0.05 & abs(All$logFC) > log2(1.5),
#                                ifelse(All$logFC > log2(1.5),'UP','DOWN'),'NOT'))
# table(All$change)
# target_res <- All[All$change %in% c('UP', 'DOWN'), ]
# length(unique(target_res$geneSymbol))# 37
# target_res <- target_res[!duplicated(target_res$geneSymbol),]
# table(target_res$geneSymbol %in% resMeR$geneSymbol)
# # FALSE  TRUE 
# # 20   17 
# table(target_res$geneSymbol %in% resInP$geneSymbol)
# # target_res <- target_res[target_res$geneSymbol %in% resInP$geneSymbol, ]
# # target_res <- target_res[target_res$geneSymbol %in% resMeR$geneSymbol, ]

############
n = 6
############
#GSE148764## with RIP and Input
############
print(RMPdeg_datasets[n,])
i = RMPdeg_datasets[n,'GSE']
cl = 'HEK293T'
RMp = 'NSUN6'
RM = 'm5C'
condition = c(rep('NSUN6KO', 2), rep('WT', 2))
# View(allDat[[i]])
RMP_allDat[RMP_allDat$study_alias %in% c('GSE155203') | RMP_allDat$sample_alias %in% c('GSM4478618', 'GSM4478619', 'GSM4478620', 'GSM4478621')
           ,c("study_alias", "sample_alias", "sample_title")]
files <- list.files(file.path(data_directory,i))
files
setwd(data_directory)
setwd(i)
###################################
#  Cell number-normalized (CNN)   #
#      MeRIP INPUT analysis       #
###################################
count_in <- read.table("metaPlotR/Input_counts.txt",header=T)
countGTF <- count_in[, !colnames(count_in) %in% c('Geneid','Chr','Start','End','Strand','Length')]
rownames(countGTF) <- count_in$Geneid
dim(countGTF) #56752     5
## Sample name
## Sample information
sample=data.frame(
  sample=colnames(countGTF),
  condition=condition)
## SummarizedExperiment input
dds<-DESeqDataSetFromMatrix(countData=countGTF,colData=sample,design=~condition)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 10
table(keep)
dds <- dds[keep,]
## Note on factor levels
dds$condition <- factor(dds$condition, levels = unique(condition))
dds$condition <- relevel(dds$condition, ref = "WT")
# Differential expression analysis-------------------------------------------
dds <- DESeq(dds)

res <-  results(dds, contrast=c("condition",unique(condition)))
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
table(res$change)
table(is.na(res$change))
# resOrdered <- as.data.frame(res[order(res$padj), ])
resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
resInP$geneSymbol <- rownames(resInP)
name = res@elementMetadata$description[2]
name = paste0(cl, '_',str_split(name, "condition ", simplify = T)[,2])
RMdeg_list[[i]] <- list( resInP)
names(RMdeg_list[[i]] ) <- paste0(gsub(' ','_',(name)))

c <- unique(condition[!condition %in% 'WT'])
countGTF_list <- list(colnames(countGTF)); names(countGTF_list) <- c
condition_list <- list(condition); names(condition_list) <- c
contrast_input_list[[i]] <- countGTF_list
condition_input_list[[i]] <- condition_list
###################################
#  Cell number-normalized (CNN)   #
#   MeRIP RIP analysis - Part1    #
###################################
## Sample information
merge_peak <- read.delim("metaPlotR/newPeakFile.txt",header=F)
countMatrix <- read.table("metaPlotR/m5C_counts.txt",header=T)
table(countMatrix$Geneid == merge_peak$V1[-1])
table(duplicated(merge_peak$V1[-1]))
dim(countMatrix)
colnames(countMatrix)
rownames(countMatrix) <- countMatrix$Geneid
countMatrix <- countMatrix[, !colnames(countMatrix) %in% c('Geneid','Chr','Start','End','Strand','Length')]
# 1. 做好分组因子即可
groupList<-factor(condition)
colData <- data.frame(row.names=colnames(countMatrix), groupList=groupList)
# 2. Exploring Method 1: DESeq2-------------------------------------------
# 构建 DESeqDataSet 对象
dds <- DESeqDataSetFromMatrix(countData = countMatrix,
                              colData = colData,
                              design = ~ groupList)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 1
table(keep)
dds <- dds[keep,]
# 接着计算每个样本的 library size 用于后续的标准化。library size 等于落在该样本所有peaks上的reads的总数，即 counts 矩阵中每列的加和。如果bFullLibrarySize设为TRUE，则会使用library的总reads数（根据原BAM/BED文件统计获得）。然后使用estimateDispersions估计统计分布，需要将参数fitType设为’local’
dds<-estimateSizeFactors(dds)
dds<-estimateDispersions(dds,fitType="local")
dds <- nbinomWaldTest(dds)
res <- results(dds)
## 提取你想要的差异分析结果，我们这里是treated组对untreated组进行比较
res <-  results(dds, contrast=c("groupList", unique(condition)))
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
res_DESeq2 = as.data.frame(res)
table(res_DESeq2$change)
res_DESeq2 = res_DESeq2[res_DESeq2$change %in% c('UP', 'DOWN'), ]
# 2. Exploring Method 2: limma-------------------------------------------
library(edgeR)
d0 <- DGEList(counts = countMatrix)
dim(d0)
# Filtering low count peaks
cutoff <- 1           
drop <- which(apply(cpm(d0), 1, max) < cutoff)
print(length(drop))
if (length(drop) >0) {
  d <- d0[-drop,] 
} else{ d <- d0}
dim(d) #number of peaks left: 8501      
d <- calcNormFactors(d, method = c('TMM'))
colnames(countMatrix)
genotype <- groupList
table(genotype)
mm <- model.matrix(~0 + genotype)#
v <- voom(d, mm, plot=F, normalize="quantile")#voom from limma:Transform count data to log2-counts per million (logtpm)
fit <- lmFit(v, mm)
head(coef(fit))
contr <- makeContrasts(paste0(colnames(coef(fit)), collapse  = ' - '), levels = colnames(coef(fit))) #Change if other comparison
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)
length(which(top.table$adj.P.Val < 0.05 & abs(top.table$logFC) > log2(1.5))) #how many are significant?
res_limma <- as.data.frame(top.table)
res_limma$change <- with(res_limma, ifelse(res_limma$adj.P.Val < 0.05 & abs(res_limma$logFC) > log2(1.5),
                                           ifelse(res_limma$logFC > log2(1.5),'UP','DOWN'),'NOT'))
table(res_limma$change)#NOT only
res_limma = res_limma[res_limma$change %in% c('UP', 'DOWN'), ]

RMdem_list[[i]] <- list(res_DESeq2, res_limma)
name = res@elementMetadata$description[2]
name = paste0(cl, '_',str_split(name, "groupList ", simplify = T)[,2])
name= paste0(gsub(' ','_',(name)))
names(RMdem_list[[i]] ) <- paste0(name, '_', c('DESeq2', 'limma'))
RMdep_list[[i]] <- list( merge_peak )
names(RMdep_list[[i]] ) <- name

c <- unique(condition[!condition %in% 'WT'])
countMatrix_list <- list(colnames(countMatrix)); names(countMatrix_list) <- c
condition_list <- list(condition); names(condition_list) <- c
contrast_IP_list[[i]] <- countMatrix_list
condition_IP_list[[i]] <- condition_list
# ###################################
# #  Cell number-normalized (CNN)   #
# #   MeRIP RIP analysis - Part2    #
# ###################################
# library(tidyverse)
# library(edgeR)
# # 1. Formatting tables as needed for downstream analysis
# d0 <- DGEList(counts = countGTF)
# # Filtering low count genes
# cutoff <- 1           
# drop <- which(apply(cpm(d0), 1, max) < cutoff)
# length(drop)
# d <- d0[-drop,] 
# d <- calcNormFactors(d)
# N <- colSums(countGTF)
# colnames(countGTF)
# genotype <- condition %>% as.factor()
# mm <- model.matrix(~0 + genotype)
# v <- voom(d, mm, plot=T, normalize="quantile")#voom from limma:Transform count data to log2-counts per million (logtpm)
# Input <- as.data.frame(v$E)
# Input$geneSymbol <- rownames(Input)
# 
# d0 <- DGEList(countMatrix)
# N <- colSums(countMatrix)  
# cutoff <- 1
# drop <- which(apply(cpm(d0), 1, max) < cutoff)
# length(drop)
# if(length(drop)>0){
#   d <- d0[-drop,] 
# }else{
#   d <- d0 
# }
# d <- calcNormFactors(d)
# genotype <- condition %>% as.factor()##c(rep(c("WT", "KO"), 6))
# mm <- model.matrix(~0 + genotype)#
# v <- voom(d, mm, plot=T, normalize="quantile")#voom from limma:Transform count data to log2-counts per million (logtpm)
# 
# RIP <- data.frame(v$E)
# RIP$peak_id <- rownames(RIP)
# RIP <- inner_join(RIP, merge_bed.anno, by=c("peak_id"))
# 
# Norm.RIP <- inner_join(RIP, Input, by="geneSymbol")
# Norm.RIP$Ave.input <- log2(rowMeans(2^Norm.RIP[,colnames(Norm.RIP) %in% colnames(countGTF)]))
# colnames(Norm.RIP) <- c(paste0("RIP", 1:length(condition)), 'peak_id','chr', 'start', 'end','geneSymbol',
#                         paste0("Input", 1:length(condition)), "Ave.input")
# # Adjusting RIP values for gene expression
# # All RIP values are corrected by dividing by the ratio (Input.sample/Input.average) (note: all values are in log2 format so this appears as "-Input#+Ave.Input")
# Norm.RIP <- mutate(Norm.RIP, 
#                    Ratio1 = RIP1-Input1+Ave.input,
#                    Ratio2 = RIP2-Input2+Ave.input,
#                    Ratio3 = RIP3-Input3+Ave.input,
#                    Ratio4 = RIP4-Input4+Ave.input,
# )
# # Formatting the "Ratio" table (i.e. the adjusted RIP peaks used for the final analysis)
# # No filtering is performed as both input and RIP data are already filtered for low counts
# Ratio <- cbind(Norm.RIP[,c('chr', 'start', 'end','geneSymbol')], 2^(Norm.RIP[,grepl('Ratio', colnames(Norm.RIP))]) )
# Ratio$RefID <- paste("Ref", 1:nrow(Ratio)) # used to bring back peak information later
# rownames(Ratio) <- Ratio$RefID
# # 2. Differential gene analysis with edgeR
# # Building the model for limma-voom
# # As the data is already fully normalized, further re-scaling is prevented by fixing the library size in the voom function
# d0 <- DGEList(Ratio[, grepl('Ratio', colnames(Ratio))])
# d <- calcNormFactors(d0)
# N <- colSums(Ratio[,grepl('Ratio', colnames(Ratio))])
# genotype <- condition %>% as.factor()
# mm <- model.matrix(~0 + genotype)
# y <- voom(d, mm, plot = T, lib.size = rep(10^6,length(condition)))  #lib.size fixed as 10^6 for all samples
# # v <- voom(d, mm, plot=T, normalize="quantile")#voom from limma:Transform count data to log2-counts per million (logtpm)
# fit <- lmFit(y, mm)
# head(coef(fit))
# # Differential analysis for inputs (Paused vs FBS)
# contr <- makeContrasts(paste0(colnames(coef(fit)), collapse  = ' - '), levels = colnames(coef(fit))) #Change if other comparison
# tmp <- contrasts.fit(fit, contr)
# tmp <- eBayes(tmp)
# top.table <- topTable(tmp, sort.by = "P", n = Inf)
# head(top.table, 20)
# length(which(top.table$adj.P.Val < 0.05 & abs(top.table$logFC) > log2(1.5))) #how many are significant?
# # Bring back peak information and export
# top.table$RefID <- rownames(top.table)
# top.table <- inner_join(top.table, Ratio, by=c("RefID"))
# # top.table <- top.table[,-c(7,18)]
# # write.table(top.table, file = paste0(export_folder,"/m6A.CNN.toptable.txt"))
# All <- inner_join(Norm.RIP, top.table[,!grepl('Ratio', colnames(top.table))], by=c('chr', 'start', 'end','geneSymbol'))
# # write.table(All, file = paste0(export_folder,"/m6A.CNN.all_results.txt"))
# All$change <- with(All, ifelse(All$adj.P.Val < 0.05 & abs(All$logFC) > log2(1.5),
#                                ifelse(All$logFC > log2(1.5),'UP','DOWN'),'NOT'))
# table(All$change)
# target_res <- All[All$change %in% c('UP', 'DOWN'), ]
# length(unique(target_res$geneSymbol))# 37
# target_res <- target_res[!duplicated(target_res$geneSymbol),]
# table(target_res$geneSymbol %in% resMeR$geneSymbol)
# # FALSE  TRUE 
# # 20   17 
# table(target_res$geneSymbol %in% resInP$geneSymbol)
# # FALSE  TRUE 
# # 29    8 
# target_res <- target_res[target_res$geneSymbol %in% resInP$geneSymbol, ]
# target_res <- target_res[target_res$geneSymbol %in% resMeR$geneSymbol, ]

############
n = 7
############
#GSE93749# with RIP and Input
############
print(RMPdeg_datasets[n,])
i = RMPdeg_datasets[n,'GSE']
cl = 'Hela'
RMp = c('NSUN2', 'ALYREFKD')
RM = 'm5C'
condition_list <- list(NSUN2KD = c(rep('NSUN2KD', 2), rep('WT', 2)),
                       ALYREFPD = c(rep('ALYREFPD', 2), rep('WT', 2))
)
# View(allDat[[i]])
ENA <- as.data.frame(RMP_allDat[grep('GSE937', RMP_allDat$study_alias ) 
                                ,c("study_alias", "sample_alias", "sample_title")])
ENA[order(ENA$sample_alias),]
files <- list.files(file.path(data_directory,i))
files
setwd(data_directory)
setwd(i)
# ALYREF-RIP
DEP <-  read_excel(file.path(data_directory,i, 'GSE93749_m5C_sites_in_ALYREF-RIP.xls'),sheet = 1)
DEG <-  read_excel(file.path(data_directory,i, 'GSE93750_all-samples-reads-count.xls'),sheet = 1)
countfiles <- list.files('./metaPlotR/')
countfiles <- countfiles[grep('.txt', countfiles)]
countfiles <- countfiles[! grepl('summary', countfiles)] #mouse embryonic fibroblasts (MEFs)
countfiles
###################################
#  Cell number-normalized (CNN)   #
#      MeRIP INPUT analysis       #
###################################
resInP_list <- vector("list", 2)
names(resInP_list) <- names(condition_list)
countGTF_list <- vector("list", 2)
names(countGTF_list) <- names(condition_list)
for(c in names(condition_list)){
  countGTF <- read.table( paste0("metaPlotR/Input_counts_", c, ".txt"),header=T)
  dim(countGTF)
  colnames(countGTF)
  rownames(countGTF) <- countGTF$Geneid
  countGTF <- countGTF[, !colnames(countGTF) %in% c('Geneid','Chr','Start','End','Strand','Length')]
  dim(countGTF) # 28265     4
  condition <- condition_list[[c]]
  ## Sample information
  sample=data.frame(
    sample=colnames(countGTF),
    condition=condition)
  ## SummarizedExperiment input
  dds<-DESeqDataSetFromMatrix(countData=countGTF,colData=sample,design=~condition)
  ## Pre-filtering
  keep <- rowSums(counts(dds)) >= 10
  table(keep)
  dds <- dds[keep,]
  ## Note on factor levels
  dds$condition <- factor(dds$condition, levels = unique(condition))
  dds$condition <- relevel(dds$condition, ref = "WT")
  # Differential expression analysis-------------------------------------------
  dds <- DESeq(dds)
  res <-  results(dds, contrast=c("condition",unique(condition)))
  res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
  length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
  res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                                 ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
  table(res$change)
  table(is.na(res$change))
  # resOrdered <- as.data.frame(res[order(res$padj), ])
  resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
  resInP$geneSymbol <- rownames(resInP)
  table(resInP$change)
  resInP_list[[c]] <- resInP
  countGTF_list[[c]] <- countGTF
}
name = paste0(names(resInP_list), ' vs WT')
name = paste0(cl, '_', name)
RMdeg_list[[i]] <- resInP_list 
names(RMdeg_list[[i]] ) <- paste0(gsub(' ','_',(name)))

contrast_input_list[[i]] <- lapply(countGTF_list, extract_names)
condition_input_list[[i]] <- condition_list

###################################
#  Cell number-normalized (CNN)   #
#   MeRIP RIP analysis - Part1    #
###################################
resMeR_list <- vector("list", length(condition_list))
names(resMeR_list) <- names(condition_list)
countMatrix_list <- vector("list", length(condition_list))
names(countMatrix_list) <- names(condition_list)
resPeak_list <- vector("list",length(condition_list))
names(resPeak_list) <- names(condition_list)
for(c in names(condition_list)){
  merge_peak <- read.delim(paste0("metaPlotR/", c, "_newPeakFile.txt"),header=F)
  countMatrix <- read.table(paste0("metaPlotR/m5C_counts_", c, ".txt"),header=T)
  table(countMatrix$Geneid == merge_peak[-1,1])
  table(duplicated(merge_peak[-1,1]))
  dim(countMatrix)
  colnames(countMatrix)
  rownames(countMatrix) <- countMatrix$Geneid
  countMatrix <- countMatrix[, !colnames(countMatrix) %in% c('Geneid','Chr','Start','End','Strand','Length')]
  # 1. 做好分组因子即可
  condition <- condition_list[[c]]
  groupList<-factor(condition)
  colData <- data.frame(row.names=colnames(countMatrix), groupList=groupList)
  # 2. Exploring Method 1: DESeq2-------------------------------------------
  # 构建 DESeqDataSet 对象
  dds <- DESeqDataSetFromMatrix(countData = countMatrix,
                                colData = colData,
                                design = ~ groupList)
  ## Pre-filtering
  keep <- rowSums(counts(dds)) >= 1
  table(keep)
  dds <- dds[keep,]
  # 接着计算每个样本的 library size 用于后续的标准化。library size 等于落在该样本所有peaks上的reads的总数，即 counts 矩阵中每列的加和。如果bFullLibrarySize设为TRUE，则会使用library的总reads数（根据原BAM/BED文件统计获得）。然后使用estimateDispersions估计统计分布，需要将参数fitType设为’local’
  dds<-estimateSizeFactors(dds)
  dds<-estimateDispersions(dds,fitType="local")
  dds <- nbinomWaldTest(dds)
  res <- results(dds)
  ## 提取你想要的差异分析结果，我们这里是treated组对untreated组进行比较
  res <-  results(dds, contrast=c("groupList", unique(condition)))
  res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
  length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
  res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                                 ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
  res_DESeq2 = as.data.frame(res)
  table(res_DESeq2$change)
  res_DESeq2 = res_DESeq2[res_DESeq2$change %in% c('UP', 'DOWN'), ]
  # 2. Exploring Method 2: limma-------------------------------------------
  library(edgeR)
  d0 <- DGEList(counts = countMatrix)
  dim(d0)
  # Filtering low count peaks
  cutoff <- 1           
  drop <- which(apply(cpm(d0), 1, max) < cutoff)
  print(length(drop))
  if (length(drop) >0) {
    d <- d0[-drop,] 
  } else{ d <- d0}
  dim(d) #number of peaks left: 8501      
  d <- calcNormFactors(d, method = c('TMM'))
  colnames(countMatrix)
  genotype <- groupList
  table(genotype)
  mm <- model.matrix(~0 + genotype)#
  v <- voom(d, mm, plot=F, normalize="quantile")#voom from limma:Transform count data to log2-counts per million (logtpm)
  fit <- lmFit(v, mm)
  head(coef(fit))
  contr <- makeContrasts(paste0(colnames(coef(fit)), collapse  = ' - '), levels = colnames(coef(fit))) #Change if other comparison
  tmp <- contrasts.fit(fit, contr)
  tmp <- eBayes(tmp)
  top.table <- topTable(tmp, sort.by = "P", n = Inf)
  head(top.table, 20)
  length(which(top.table$adj.P.Val < 0.05 & abs(top.table$logFC) > log2(1.5))) #how many are significant?
  res_limma <- as.data.frame(top.table)
  res_limma$change <- with(res_limma, ifelse(res_limma$adj.P.Val < 0.05 & abs(res_limma$logFC) > log2(1.5),
                                             ifelse(res_limma$logFC > log2(1.5),'UP','DOWN'),'NOT'))
  table(res_limma$change)
  res_limma = res_limma[res_limma$change %in% c('UP', 'DOWN'), ]
  resMeR_list[[c]] <- list(res_DESeq2, res_limma)
  countMatrix_list[[c]] <- countMatrix
  resPeak_list[[c]] <- merge_peak
  name = res@elementMetadata$description[2]
  name = paste0(cl, '_',str_split(name, "groupList ", simplify = T)[,2])
  name = paste0(gsub(' ','_',(name)))
  names(resMeR_list[[c]] ) <- paste0(name, '_', c('DESeq2', 'limma'))
}
# Convert the list of lists to a single list
# L_flat <- unlist(resMeR_list, recursive = FALSE)
# # Convert the flattened list to a list of 4 elements
# L_final <- as.list(L_flat)
RMdem_list[[i]] <- as.list(unlist(resMeR_list, recursive = FALSE))
name <- str_split(names(RMdem_list[[i]] ), "[.]", simplify = T)[,2]
names(RMdem_list[[i]] ) <- name

name = paste0(names(resPeak_list), ' vs WT')
name = paste0(cl, '_', name)
RMdep_list[[i]] <- resPeak_list 
names(RMdep_list[[i]] ) <- paste0(gsub(' ','_',(name)))

contrast_IP_list[[i]] <- lapply(countMatrix_list, extract_names)
condition_IP_list[[i]] <- condition_list

# ###################################
# #  Cell number-normalized (CNN)   #
# #   MeRIP RIP analysis - Part2    #
# ###################################
# target_res <- vector("list", 2)
# names(target_res) <- names(condition_list)
# for(c in names(condition_list)){  
#   library(tidyverse)
#   library(edgeR)
#   countGTF <- countGTF_list[[c]]
#   # 1. Formatting tables as needed for downstream analysis
#   d0 <- DGEList(counts = countGTF)
#   # Filtering low count genes
#   cutoff <- 1           
#   drop <- which(apply(cpm(d0), 1, max) < cutoff)
#   length(drop)
#   d <- d0[-drop,] 
#   d <- calcNormFactors(d)
#   N <- colSums(countGTF)
#   colnames(countGTF)
#   genotype <- condition %>% as.factor()
#   mm <- model.matrix(~0 + genotype)
#   v <- voom(d, mm, plot=T, normalize="quantile")#voom from limma:Transform count data to log2-counts per million (logtpm)
#   Input <- as.data.frame(v$E)
#   Input$geneSymbol <- rownames(Input)
#   
#   countMatrix <- countMatrix_list[[c]]
#   d0 <- DGEList(countMatrix)
#   N <- colSums(countMatrix)  
#   cutoff <- 1
#   drop <- which(apply(cpm(d0), 1, max) < cutoff)
#   length(drop)
#   if(length(drop)>0){
#     d <- d0[-drop,] 
#   }else{
#     d <- d0 
#   }
#   d <- calcNormFactors(d)
#   genotype <- condition %>% as.factor()##c(rep(c("WT", "KO"), 6))
#   mm <- model.matrix(~0 + genotype)#
#   v <- voom(d, mm, plot=T, normalize="quantile")#voom from limma:Transform count data to log2-counts per million (logtpm)
#   
#   RIP <- data.frame(v$E)
#   RIP$peak_id <- rownames(RIP)
#   RIP <- inner_join(RIP, merge_bed.anno, by=c("peak_id"))
#   
#   Norm.RIP <- inner_join(RIP, Input, by="geneSymbol")
#   Norm.RIP$Ave.input <- log2(rowMeans(2^Norm.RIP[,colnames(Norm.RIP) %in% colnames(countGTF)]))
#   colnames(Norm.RIP) <- c(paste0("RIP", 1:length(condition)), 'peak_id','chr', 'start', 'end','geneSymbol',
#                           paste0("Input", 1:length(condition)), "Ave.input")
#   # Adjusting RIP values for gene expression
#   # All RIP values are corrected by dividing by the ratio (Input.sample/Input.average) (note: all values are in log2 format so this appears as "-Input#+Ave.Input")
#   Norm.RIP <- mutate(Norm.RIP, 
#                      Ratio1 = RIP1-Input1+Ave.input,
#                      Ratio2 = RIP2-Input2+Ave.input,
#                      Ratio3 = RIP3-Input3+Ave.input,
#                      Ratio4 = RIP4-Input4+Ave.input,
#   )
#   # Formatting the "Ratio" table (i.e. the adjusted RIP peaks used for the final analysis)
#   # No filtering is performed as both input and RIP data are already filtered for low counts
#   Ratio <- cbind(Norm.RIP[,c('chr', 'start', 'end','geneSymbol')], 2^(Norm.RIP[,grepl('Ratio', colnames(Norm.RIP))]) )
#   Ratio$RefID <- paste("Ref", 1:nrow(Ratio)) # used to bring back peak information later
#   rownames(Ratio) <- Ratio$RefID
#   # 2. Differential gene analysis with edgeR
#   # Building the model for limma-voom
#   # As the data is already fully normalized, further re-scaling is prevented by fixing the library size in the voom function
#   d0 <- DGEList(Ratio[, grepl('Ratio', colnames(Ratio))])
#   d <- calcNormFactors(d0)
#   N <- colSums(Ratio[,grepl('Ratio', colnames(Ratio))])
#   genotype <- condition %>% as.factor()
#   mm <- model.matrix(~0 + genotype)
#   y <- voom(d, mm, plot = T, lib.size = rep(10^6,length(condition)))  #lib.size fixed as 10^6 for all samples
#   # v <- voom(d, mm, plot=T, normalize="quantile")#voom from limma:Transform count data to log2-counts per million (logtpm)
#   fit <- lmFit(y, mm)
#   head(coef(fit))
#   # Differential analysis for inputs (Paused vs FBS)
#   contr <- makeContrasts(paste0(colnames(coef(fit)), collapse  = ' - '), levels = colnames(coef(fit))) #Change if other comparison
#   tmp <- contrasts.fit(fit, contr)
#   tmp <- eBayes(tmp)
#   top.table <- topTable(tmp, sort.by = "P", n = Inf)
#   head(top.table, 20)
#   length(which(top.table$adj.P.Val < 0.05 & abs(top.table$logFC) > log2(1.5))) #how many are significant?
#   # Bring back peak information and export
#   top.table$RefID <- rownames(top.table)
#   top.table <- inner_join(top.table, Ratio, by=c("RefID"))
#   # top.table <- top.table[,-c(7,18)]
#   # write.table(top.table, file = paste0(export_folder,"/m6A.CNN.toptable.txt"))
#   All <- inner_join(Norm.RIP, top.table[,!grepl('Ratio', colnames(top.table))], by=c('chr', 'start', 'end','geneSymbol'))
#   # write.table(All, file = paste0(export_folder,"/m6A.CNN.all_results.txt"))
#   All$change <- with(All, ifelse(All$adj.P.Val < 0.05 & abs(All$logFC) > log2(1.5),
#                                  ifelse(All$logFC > log2(1.5),'UP','DOWN'),'NOT'))
#   table(All$change)
#   All <- All[All$change %in% c('UP', 'DOWN'), ]
#   print(length(unique(All$geneSymbol)))# 67
#   table(All$geneSymbol %in% resMeR$geneSymbol)
#   table(All$geneSymbol %in% resInP$geneSymbol)
#   target_res[[c]] <- All
#   }
# target_res <- target_res[target_res$geneSymbol %in% resInP$geneSymbol, ]
# target_res <- target_res[target_res$geneSymbol %in% resMeR$geneSymbol, ]
# ###################################################################
# # remove duplicated genes (according to the baseMean of DESeq2 results)
# ###################################################################
# for(c in names(condition_list)){
#   resMeR <- resMeR_list[[c]]
#   print(dim(resMeR))
#   length(unique(resMeR$geneSymbol))
#   # tmp=by(gene_summary,gene_summary$Symbol,function(x) rownames(x)[which.min(x$Pvalue_Adj)])
#   tmp=by(resMeR,resMeR$geneSymbol,function(x) rownames(x)[which.max(x$baseMean)])
#   probes = as.character(tmp)
#   length(probes) # 7785
#   target_res[[c]] <- resMeR[rownames(resMeR) %in% probes,]
#   index = target_res[[c]]$geneSymbol %in% resInP_list[[c]]$geneSymbol
#   print( table(index) )
#   target_res[[c]] <- target_res[[c]][index,]
# }

############
n = 8
############
#GSE74771# with RIP without Input
############
print(RMPdeg_datasets[n,])
i = RMPdeg_datasets[n,'GSE']
cl = 'Hela' # (use data from GSE93750/GSE93749 as Input)
RMp = 'NSUN2'
RM = 'm5C'
condition = c(rep('NSUN2KO', 2), rep('WT', 2))
ENA <- as.data.frame(RMP_allDat[grepl('GSE74771', RMP_allDat$study_alias ) & grepl('HeLa', RMP_allDat$sample_title)
                                ,c("study_alias", "sample_alias", "sample_title")])
ENA[order(ENA$sample_alias),]
# View(allDat[[i]])
files <- list.files(file.path(data_directory,i))
files
setwd(data_directory)
setwd(i)
###################################
#  Cell number-normalized (CNN)   #
#      MeRIP INPUT analysis       #
###################################
count_in <- read.table("metaPlotR/Input_counts_NSUN2KO.txt",header=T)
countGTF <- count_in[, !colnames(count_in) %in% c('Geneid','Chr','Start','End','Strand','Length')]
rownames(countGTF) <- count_in$Geneid
dim(countGTF) #28265     4
## Sample name
## Sample information
sample=data.frame(
  sample=colnames(countGTF),
  condition=condition)
## SummarizedExperiment input
dds<-DESeqDataSetFromMatrix(countData=countGTF,colData=sample,design=~condition)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 10
table(keep)
dds <- dds[keep,]
## Note on factor levels
dds$condition <- factor(dds$condition, levels = unique(condition))
dds$condition <- relevel(dds$condition, ref = "WT")
# Differential expression analysis-------------------------------------------
dds <- DESeq(dds)

res <-  results(dds, contrast=c("condition",unique(condition)))
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
table(res$change)
table(is.na(res$change))
# resOrdered <- as.data.frame(res[order(res$padj), ])
resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
resInP$geneSymbol <- rownames(resInP)
name = res@elementMetadata$description[2]
name = paste0(cl, '_',str_split(name, "condition ", simplify = T)[,2])
RMdeg_list[[i]] <- list( resInP)
names(RMdeg_list[[i]] ) <- paste0(gsub(' ','_',(name)))
c <- unique(condition[!condition %in% 'WT'])
countGTF_list <- list(colnames(countGTF)); names(countGTF_list) <- c
condition_list <- list(condition); names(condition_list) <- c
contrast_input_list[[i]] <- countGTF_list
condition_input_list[[i]] <- condition_list
###################################
#  Cell number-normalized (CNN)   #
#   MeRIP RIP analysis - Part1    #
###################################
## Sample information
merge_peak <- read.delim("metaPlotR/NSUN2KO_newPeakFile.txt",header=F)
countMatrix <- read.table("metaPlotR/m5C_counts_NSUN2KO.txt",header=T)
table(countMatrix$Geneid == merge_peak$V1[-1])
table(duplicated(merge_peak$V1[-1]))
dim(countMatrix)
colnames(countMatrix)
rownames(countMatrix) <- countMatrix$Geneid
countMatrix <- countMatrix[, !colnames(countMatrix) %in% c('Geneid','Chr','Start','End','Strand','Length')]
# 1. 做好分组因子即可
groupList<-factor(condition)
colData <- data.frame(row.names=colnames(countMatrix), groupList=groupList)
# 2. Exploring Method 1: DESeq2-------------------------------------------
# 构建 DESeqDataSet 对象
dds <- DESeqDataSetFromMatrix(countData = countMatrix,
                              colData = colData,
                              design = ~ groupList)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 1
table(keep)
dds <- dds[keep,]
# 接着计算每个样本的 library size 用于后续的标准化。library size 等于落在该样本所有peaks上的reads的总数，即 counts 矩阵中每列的加和。如果bFullLibrarySize设为TRUE，则会使用library的总reads数（根据原BAM/BED文件统计获得）。然后使用estimateDispersions估计统计分布，需要将参数fitType设为’local’
dds<-estimateSizeFactors(dds)
dds<-estimateDispersions(dds,fitType="local")
dds <- nbinomWaldTest(dds)
res <- results(dds)
## 提取你想要的差异分析结果，我们这里是treated组对untreated组进行比较
res <-  results(dds, contrast=c("groupList", unique(condition)))
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
res_DESeq2 = as.data.frame(res)
table(res_DESeq2$change)#Not only
res_DESeq2 = res_DESeq2[res_DESeq2$change %in% c('UP', 'DOWN'), ]
# 2. Exploring Method 2: limma-------------------------------------------
library(edgeR)
d0 <- DGEList(counts = countMatrix)
dim(d0)
# Filtering low count peaks
cutoff <- 1           
drop <- which(apply(cpm(d0), 1, max) < cutoff)
print(length(drop))
if (length(drop) >0) {
  d <- d0[-drop,] 
} else{ d <- d0}
dim(d) #number of peaks left: 8501      
d <- calcNormFactors(d, method = c('TMM'))
colnames(countMatrix)
genotype <- groupList
table(genotype)
mm <- model.matrix(~0 + genotype)#
v <- voom(d, mm, plot=F, normalize="quantile")#voom from limma:Transform count data to log2-counts per million (logtpm)
fit <- lmFit(v, mm)
head(coef(fit))
contr <- makeContrasts(paste0(colnames(coef(fit)), collapse  = ' - '), levels = colnames(coef(fit))) #Change if other comparison
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)
length(which(top.table$adj.P.Val < 0.05 & abs(top.table$logFC) > log2(1.5))) #how many are significant?
res_limma <- as.data.frame(top.table)
res_limma$change <- with(res_limma, ifelse(res_limma$adj.P.Val < 0.05 & abs(res_limma$logFC) > log2(1.5),
                                           ifelse(res_limma$logFC > log2(1.5),'UP','DOWN'),'NOT'))
table(res_limma$change)#NOT only
res_limma = res_limma[res_limma$change %in% c('UP', 'DOWN'), ]

RMdem_list[[i]] <- list(res_DESeq2, res_limma)
name = res@elementMetadata$description[2]
name = paste0(cl, '_',str_split(name, "groupList ", simplify = T)[,2])
name= paste0(gsub(' ','_',(name)))
names(RMdem_list[[i]] ) <- paste0(name, '_', c('DESeq2', 'limma'))
RMdep_list[[i]] <- list( merge_peak )
names(RMdep_list[[i]] ) <- name

c <- unique(condition[!condition %in% 'WT'])
countMatrix_list <- list(colnames(countMatrix)); names(countMatrix_list) <- c
condition_list <- list(condition); names(condition_list) <- c
contrast_IP_list[[i]] <- countMatrix_list
condition_IP_list[[i]] <- condition_list
# ###################################
# #  Cell number-normalized (CNN)   #
# #   MeRIP RIP analysis - Part2    #
# ###################################
# library(tidyverse)
# library(edgeR)
# # 1. Formatting tables as needed for downstream analysis
# d0 <- DGEList(counts = countGTF)
# # Filtering low count genes
# cutoff <- 1           
# drop <- which(apply(cpm(d0), 1, max) < cutoff)
# length(drop)
# d <- d0[-drop,] 
# d <- calcNormFactors(d)
# N <- colSums(countGTF)
# colnames(countGTF)
# genotype <- condition %>% as.factor()
# mm <- model.matrix(~0 + genotype)
# v <- voom(d, mm, plot=T, normalize="quantile")#voom from limma:Transform count data to log2-counts per million (logtpm)
# Input <- as.data.frame(v$E)
# Input$geneSymbol <- rownames(Input)
# 
# d0 <- DGEList(countMatrix)
# N <- colSums(countMatrix)  
# cutoff <- 1
# drop <- which(apply(cpm(d0), 1, max) < cutoff)
# length(drop)
# if(length(drop)>0){
#   d <- d0[-drop,] 
# }else{
#   d <- d0 
# }
# d <- calcNormFactors(d)
# genotype <- condition %>% as.factor()##c(rep(c("WT", "KO"), 6))
# mm <- model.matrix(~0 + genotype)#
# v <- voom(d, mm, plot=T, normalize="quantile")#voom from limma:Transform count data to log2-counts per million (logtpm)
# 
# RIP <- data.frame(v$E)
# RIP$peak_id <- rownames(RIP)
# RIP <- inner_join(RIP, merge_bed.anno, by=c("peak_id"))
# 
# Norm.RIP <- inner_join(RIP, Input, by="geneSymbol")
# Norm.RIP$Ave.input <- log2(rowMeans(2^Norm.RIP[,colnames(Norm.RIP) %in% colnames(countGTF)]))
# colnames(Norm.RIP) <- c(paste0("RIP", 1:length(condition)), 'peak_id','chr', 'start', 'end','geneSymbol',
#                         paste0("Input", 1:length(condition)), "Ave.input")
# # Adjusting RIP values for gene expression
# # All RIP values are corrected by dividing by the ratio (Input.sample/Input.average) (note: all values are in log2 format so this appears as "-Input#+Ave.Input")
# Norm.RIP <- mutate(Norm.RIP, 
#                    Ratio1 = RIP1-Input1+Ave.input,
#                    Ratio2 = RIP2-Input2+Ave.input,
#                    Ratio3 = RIP3-Input3+Ave.input,
#                    Ratio4 = RIP4-Input4+Ave.input,
# )
# # Formatting the "Ratio" table (i.e. the adjusted RIP peaks used for the final analysis)
# # No filtering is performed as both input and RIP data are already filtered for low counts
# Ratio <- cbind(Norm.RIP[,c('chr', 'start', 'end','geneSymbol')], 2^(Norm.RIP[,grepl('Ratio', colnames(Norm.RIP))]) )
# Ratio$RefID <- paste("Ref", 1:nrow(Ratio)) # used to bring back peak information later
# rownames(Ratio) <- Ratio$RefID
# # 2. Differential gene analysis with edgeR
# # Building the model for limma-voom
# # As the data is already fully normalized, further re-scaling is prevented by fixing the library size in the voom function
# d0 <- DGEList(Ratio[, grepl('Ratio', colnames(Ratio))])
# d <- calcNormFactors(d0)
# N <- colSums(Ratio[,grepl('Ratio', colnames(Ratio))])
# genotype <- condition %>% as.factor()
# mm <- model.matrix(~0 + genotype)
# y <- voom(d, mm, plot = T, lib.size = rep(10^6,length(condition)))  #lib.size fixed as 10^6 for all samples
# # v <- voom(d, mm, plot=T, normalize="quantile")#voom from limma:Transform count data to log2-counts per million (logtpm)
# fit <- lmFit(y, mm)
# head(coef(fit))
# # Differential analysis for inputs (Paused vs FBS)
# contr <- makeContrasts(paste0(colnames(coef(fit)), collapse  = ' - '), levels = colnames(coef(fit))) #Change if other comparison
# tmp <- contrasts.fit(fit, contr)
# tmp <- eBayes(tmp)
# top.table <- topTable(tmp, sort.by = "P", n = Inf)
# head(top.table, 20)
# length(which(top.table$adj.P.Val < 0.05 & abs(top.table$logFC) > log2(1.5))) #how many are significant?
# # Bring back peak information and export
# top.table$RefID <- rownames(top.table)
# top.table <- inner_join(top.table, Ratio, by=c("RefID"))
# # top.table <- top.table[,-c(7,18)]
# # write.table(top.table, file = paste0(export_folder,"/m6A.CNN.toptable.txt"))
# All <- inner_join(Norm.RIP, top.table[,!grepl('Ratio', colnames(top.table))], by=c('chr', 'start', 'end','geneSymbol'))
# # write.table(All, file = paste0(export_folder,"/m6A.CNN.all_results.txt"))
# All$change <- with(All, ifelse(All$adj.P.Val < 0.05 & abs(All$logFC) > log2(1.5),
#                                ifelse(All$logFC > log2(1.5),'UP','DOWN'),'NOT'))
# table(All$change)
# target_res <- All[All$change %in% c('UP', 'DOWN'), ]
# length(unique(target_res$geneSymbol))# 37
# target_res <- target_res[!duplicated(target_res$geneSymbol),]
# table(target_res$geneSymbol %in% resMeR$geneSymbol)
# # FALSE  TRUE 
# # 20   17 
# table(target_res$geneSymbol %in% resInP$geneSymbol)
# # target_res <- target_res[target_res$geneSymbol %in% resInP$geneSymbol, ]
# # target_res <- target_res[target_res$geneSymbol %in% resMeR$geneSymbol, ]

############
n = 9
############
#GSE125046# with RIP without Input
############
print(RMPdeg_datasets[n,]) # hg38
i = RMPdeg_datasets[n,'GSE']
cl = c('H9','HEK293') # human embryonic fibroblasts/stem cells (H9) and HEK293 cells
RMp = 'NSUN6'
RM = 'm5C'
condition_list <- list(H9.KO.1 = c(rep('NSUN6KO', 8), rep('WT', 8)), #for H9: Crispr /Cas9 mediated (HDR) deletion of EXON 2 of NSUN6, clone 1
                       H9.KO.2 = c(rep('NSUN6KO', 8), rep('WT', 8)), #for H9: Crispr /Cas9 mediated (HDR) deletion of EXON 2 of NSUN6, clone 2
                       H9.KO.3 = c(rep('NSUN6KO', 4), rep('WT', 8)), #for H9: Crispr /Cas9 mediated (NHEJ) InDel of exon9 of NSUN6 (deletion of 7nt on allele1 and 50 nt on allele 2 causes framshift and premature stop)
                       H9.OEX = c(rep('NSUN6OE', 8), rep('WT', 8)), #for H9: WT H9 carrying insertions of a piggybac vector coding for a flagged version of NSUN6, puromycin selection, clonal
                       
                       HEK.KO.1 = c(rep('NSUN6KO', 4), rep('WT', 4)), #for HEK: Crispr /Cas9 mediated (NHEJ) InDel of exon2 of NSUN6 (deletion causes frameshift and premature stop), clone1
                       HEK.KO.2 = c(rep('NSUN6KO', 4), rep('WT', 4)) #for HEK: Crispr /Cas9 mediated (NHEJ) InDel of exon2 of NSUN6 (deletion causes frameshift and premature stop), clone2
)
RMP_allDat[grep('GSE125046', RMP_allDat$study_alias ) 
           ,c("study_alias", "sample_alias", "sample_title")]
ENA <- as.data.frame(RMP_allDat[grepl('GSE125046', RMP_allDat$study_alias ) & (! grepl('rescue|EV|crispr', RMP_allDat$sample_title))
                                ,c("study_alias", "sample_alias", "sample_title")])
ENA[order(ENA$sample_alias),]

RMP_allDat[grepl('GSE125046', RMP_allDat$study_alias ) & (! grepl('rescue|EV|crispr', RMP_allDat$sample_title))
           ,c("study_alias", "sample_alias", "sample_title")]
skip <- RMP_allDat[grepl('GSE125046', RMP_allDat$study_alias ) & grepl('rescue|EV|crispr', RMP_allDat$sample_title)
                   , "sample_alias"]
paste0(skip$sample_alias, collapse = "|")

# View(allDat[[i]]) # batch 1-Illumina HiSeq 2500
# setwd(data_directory)
# gse <- getGEO('GSE125046',destdir=".", AnnotGPL=F, getGPL=F)
# pdata <- pData(gse[[2]]) # batch 2,3,4 - Illumina HiSeq 4000
# allDat[[i]] <- rbind(allDat[[i]], pdata)
# The bismark coverage reports are tab delimited and contain methylation level at each covered cytosine in the following format: <chromosome> <start position> <end position> <methylation percentage> <count methylated> <count unmethylated>
files <- list.files(file.path(data_directory,i))
files 
getwd()
setwd(data_directory)
setwd(i)
bedfiles <- list.files('./metaPlotR/')
bedfiles <- bedfiles[grep('.bed', bedfiles)]
# bedfiles <- bedfiles[grep('KO.1|H9.WT', bedfiles)] # clone 1-Crispr /Cas9 mediated (HDR) deletion of EXON 2 of NSUN6
# bedfiles <- bedfiles[grep('KO.2|H9.WT', bedfiles)] # clone 2-Crispr /Cas9 mediated (HDR) deletion of EXON 2 of NSUN6
# bedfiles <- bedfiles[grep('KO.3|H9.WT', bedfiles)] # Crispr /Cas9 mediated (NHEJ) InDel of exon9 of NSUN6 (deletion of 7nt on allele1 and 50 nt on allele 2 causes framshift and premature stop)
# bedfiles <- bedfiles[grep('HEK.KO1|HEK.WT', bedfiles)]
bedfiles <- bedfiles[grep('HEK.KO2|HEK.WT', bedfiles)]
paste0(bedfiles, collapse = ' ')
x <- 'GSM3561740_H9.KO.1.1.bam GSM3561741_H9.KO.1.2.bam GSM3561742_H9.KO.1.3.bam GSM3561743_H9.KO.1.4.bam GSM3561744_H9.KO.1o.1.bam GSM3561745_H9.KO.1o.2.bam GSM3561746_H9.KO.1o.3.bam GSM3561747_H9.KO.1o.4.bam GSM3561768_H9.WT.1.1.bam GSM3561769_H9.WT.1.2.bam GSM3561770_H9.WT.1.3.bam GSM3561771_H9.WT.1.4.bam GSM3561772_H9.WTo.1.bam GSM3561773_H9.WTo.2.bam GSM3561774_H9.WTo.3.bam GSM3561775_H9.WTo.4.bam'
x <- as.vector(str_split(x, ' ', simplify = T))
x <- str_split(x, '_H9', simplify = T)[,1]
x <- paste0('$outdir/GSE125046/', x, '_IP.bam')
paste(x, collapse = ' ')
ENA <- ENA[order(ENA$sample_alias),]
y <- ENA[grepl('H9.KO.2|H9.WT',ENA$sample_title),]
y <- ENA[grepl('H9.KO.3|H9.WT',ENA$sample_title),]
y <- ENA[grepl('H9.OEX|H9.WT',ENA$sample_title),]
y <- ENA[grepl('HEK.KO.1|HEK.WT',ENA$sample_title),]
y <- ENA[grepl('HEK.KO.2|HEK.WT',ENA$sample_title),]
y
x <- paste0('$outdir/GSE125046/', y$sample_alias, '_IP.bam')
paste(x, collapse = ' ')
# rm *.bam*; rm *.bw
countfiles <- list.files('./metaPlotR/')
countfiles <- countfiles[grep('counts', countfiles)]
countfiles <- countfiles[! grepl('summary', countfiles)] 
countfiles
###################################
#         RIP-seq analysis        #
###################################
resMeR_list <- vector("list", length(condition_list))
names(resMeR_list) <- names(condition_list)
countMatrix_list <- vector("list", length(condition_list))
names(countMatrix_list) <- names(condition_list)
resPeak_list <- vector("list",length(condition_list))
names(resPeak_list) <- names(condition_list)
for(c in names(condition_list)){
  # merge_peak <- read.delim(paste0("metaPlotR/", c, "_newPeakFile.txt"),header=F)
  # countMatrix <- read.table(paste0("metaPlotR/m5C_counts_", c, ".txt"),header=T)
  merge_peak <- read.delim(paste0("metaPlotR/", c, "_filteredPeak.txt"),header=F)
  countMatrix <- read.table(paste0("metaPlotR/", c, "_filteredCount.txt"),header=T)
  table(countMatrix$Geneid == merge_peak[-1,1])
  table(duplicated(merge_peak[-1,1]))
  dim(countMatrix)
  colnames(countMatrix)
  rownames(countMatrix) <- countMatrix$Geneid
  countMatrix <- countMatrix[, !colnames(countMatrix) %in% c('Geneid','Chr','Start','End','Strand','Length')]
  # 1. 做好分组因子即可
  condition <- condition_list[[c]]
  groupList <- factor(condition)
  colData <- data.frame(row.names=colnames(countMatrix), groupList=groupList)
  # 2. Exploring Method 1: DESeq2-------------------------------------------
  # 构建 DESeqDataSet 对象
  dds <- DESeqDataSetFromMatrix(countData = countMatrix,
                                colData = colData,
                                design = ~ groupList)
  ## Pre-filtering
  keep <- rowSums(counts(dds)) >= 1
  table(keep)
  dds <- dds[keep,]
  # 接着计算每个样本的 library size 用于后续的标准化。library size 等于落在该样本所有peaks上的reads的总数，即 counts 矩阵中每列的加和。如果bFullLibrarySize设为TRUE，则会使用library的总reads数（根据原BAM/BED文件统计获得）。然后使用estimateDispersions估计统计分布，需要将参数fitType设为’local’
  dds<-estimateSizeFactors(dds)
  dds<-estimateDispersions(dds,fitType="local")
  dds <- nbinomWaldTest(dds)
  res <- results(dds)
  ## 提取你想要的差异分析结果，我们这里是treated组对untreated组进行比较
  res <-  results(dds, contrast=c("groupList", unique(condition)))
  res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
  length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
  res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                                 ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
  res_DESeq2 = as.data.frame(res)
  table(res_DESeq2$change)
  res_DESeq2 = res_DESeq2[res_DESeq2$change %in% c('UP', 'DOWN'), ]
  # 2. Exploring Method 2: limma-------------------------------------------
  library(edgeR)
  d0 <- DGEList(counts = countMatrix)
  dim(d0)
  # Filtering low count peaks
  cutoff <- 1           
  drop <- which(apply(cpm(d0), 1, max) < cutoff)
  print(length(drop))
  if (length(drop) >0) {
    d <- d0[-drop,] 
  } else{ d <- d0}
  dim(d) #number of peaks left: 8501      
  d <- calcNormFactors(d, method = c('TMM'))
  colnames(countMatrix)
  genotype <- groupList
  table(genotype)
  mm <- model.matrix(~0 + genotype)#
  v <- voom(d, mm, plot=F, normalize="quantile")#voom from limma:Transform count data to log2-counts per million (logtpm)
  fit <- lmFit(v, mm)
  head(coef(fit))
  contr <- makeContrasts(paste0(colnames(coef(fit)), collapse  = ' - '), levels = colnames(coef(fit))) #Change if other comparison
  tmp <- contrasts.fit(fit, contr)
  tmp <- eBayes(tmp)
  top.table <- topTable(tmp, sort.by = "P", n = Inf)
  head(top.table, 20)
  length(which(top.table$adj.P.Val < 0.05 & abs(top.table$logFC) > log2(1.5))) #how many are significant?
  res_limma <- as.data.frame(top.table)
  res_limma$change <- with(res_limma, ifelse(res_limma$adj.P.Val < 0.05 & abs(res_limma$logFC) > log2(1.5),
                                             ifelse(res_limma$logFC > log2(1.5),'UP','DOWN'),'NOT'))
  table(res_limma$change)
  res_limma = res_limma[res_limma$change %in% c('UP', 'DOWN'), ]
  resMeR_list[[c]] <- list(res_DESeq2, res_limma)
  countMatrix_list[[c]] <- countMatrix
  resPeak_list[[c]] <- merge_peak
  name = res@elementMetadata$description[2]
  # name = paste0(cl, '_',str_split(name, "groupList ", simplify = T)[,2])
  name = str_split(name, "groupList ", simplify = T)[,2]
  name = paste0(gsub(' ','_',(name)))
  names(resMeR_list[[c]] ) <- paste0(name, '_', c('DESeq2', 'limma'))
}
RMdem_list[[i]] <- as.list(unlist(resMeR_list, recursive = FALSE))
name <- gsub('[.]', '_', names(RMdem_list[[i]] ))
name <- gsub('H9_KO_', 'H9.', name)
name <- gsub('H9_OEX_', 'H9_', name)
name <- gsub('HEK_KO_', 'HEK293.', name)
names(RMdem_list[[i]] ) <- name
name
RMdep_list[[i]] <- resPeak_list 
name <- unique(gsub('(_limma|_DESeq2)', '', name))
names(RMdep_list[[i]] ) <- name
contrast_IP_list[[i]] <- lapply(countMatrix_list, extract_names)
names(contrast_IP_list[[i]]) <- gsub('_vs_WT', '',names(RMdep_list[[i]] ) )
condition_IP_list[[i]] <- condition_list
names(condition_IP_list[[i]]) <- gsub('_vs_WT', '',names(RMdep_list[[i]] ) )
############
n = 10
############
#GSE109183# with miCLIPs without Ctrl 
############
print(RMPdeg_datasets[n,])#
i = RMPdeg_datasets[n,'GSE']
cl = c('HEK293')#'HAP1', (use data from GSE66012 as WT)
RMp = 'TRMT2A'
RM = 'm5U'
condition = c(rep('TRMT2APD', 2), rep('WT', 2))
ENA <- as.data.frame(RMP_allDat[grepl('GSE109183', RMP_allDat$study_alias ) #& (! grepl('rescue|EV|crispr', RMP_allDat$sample_title))
                                ,c("study_alias", "sample_alias", "sample_title")])
ENA[order(ENA$sample_alias),]
skip <- RMP_allDat[grepl('GSE109183', RMP_allDat$study_alias ) & grepl('miCLIP', RMP_allDat$sample_title)
                   , "sample_alias"]
paste0(skip$sample_alias, collapse = "|")
# View(allDat[[i]])
files <- list.files(file.path(data_directory,i))
files # The processed bedgraph files contains the start position of unique genome mapping reads. Each position represent accumulation of TRMT2A binding sites.
setwd(data_directory)
setwd(i)
countfiles <- list.files('./metaPlotR/')
countfiles <- countfiles[grep('counts', countfiles)]
countfiles <- countfiles[! grepl('summary', countfiles)] 
countfiles
###################################
#         RIP-seq analysis        #
###################################
merge_peak <- read.delim("metaPlotR/newPeakFile.txt",header=F)
countMatrix <- read.table("metaPlotR/m5U_counts.txt",header=T)
table(countMatrix$Geneid == merge_peak$V1[-1])
table(duplicated(merge_peak$V1[-1]))
dim(countMatrix)
colnames(countMatrix)
rownames(countMatrix) <- countMatrix$Geneid
countMatrix <- countMatrix[, !colnames(countMatrix) %in% c('Geneid','Chr','Start','End','Strand','Length')]
# 1. 做好分组因子即可
groupList<-factor(condition)
colData <- data.frame(row.names=colnames(countMatrix), groupList=groupList)
# 2. Exploring Method 1: DESeq2-------------------------------------------
# 构建 DESeqDataSet 对象
dds <- DESeqDataSetFromMatrix(countData = countMatrix,
                              colData = colData,
                              design = ~ groupList)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 1
table(keep)
dds <- dds[keep,]
# 接着计算每个样本的 library size 用于后续的标准化。library size 等于落在该样本所有peaks上的reads的总数，即 counts 矩阵中每列的加和。如果bFullLibrarySize设为TRUE，则会使用library的总reads数（根据原BAM/BED文件统计获得）。然后使用estimateDispersions估计统计分布，需要将参数fitType设为’local’
dds<-estimateSizeFactors(dds)
dds<-estimateDispersions(dds,fitType="local")
dds <- nbinomWaldTest(dds)
res <- results(dds)
## 提取你想要的差异分析结果，我们这里是treated组对untreated组进行比较
res <-  results(dds, contrast=c("groupList", unique(condition)))
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
res_DESeq2 = as.data.frame(res)
table(res_DESeq2$change)
res_DESeq2 = res_DESeq2[res_DESeq2$change %in% c('UP', 'DOWN'), ]
# 2. Exploring Method 2: limma-------------------------------------------
library(edgeR)
d0 <- DGEList(counts = countMatrix)
dim(d0)
# Filtering low count peaks
cutoff <- 1           
drop <- which(apply(cpm(d0), 1, max) < cutoff)
print(length(drop))
if (length(drop) >0) {
  d <- d0[-drop,] 
} else{ d <- d0}
dim(d) #number of peaks left: 8501      
d <- calcNormFactors(d, method = c('TMM'))
colnames(countMatrix)
genotype <- groupList
table(genotype)
mm <- model.matrix(~0 + genotype)#
v <- voom(d, mm, plot=F, normalize="quantile")#voom from limma:Transform count data to log2-counts per million (logtpm)
fit <- lmFit(v, mm)
head(coef(fit))
contr <- makeContrasts(paste0(colnames(coef(fit)), collapse  = ' - '), levels = colnames(coef(fit))) #Change if other comparison
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)
length(which(top.table$adj.P.Val < 0.05 & abs(top.table$logFC) > log2(1.5))) #how many are significant?
res_limma <- as.data.frame(top.table)
res_limma$change <- with(res_limma, ifelse(res_limma$adj.P.Val < 0.05 & abs(res_limma$logFC) > log2(1.5),
                                           ifelse(res_limma$logFC > log2(1.5),'UP','DOWN'),'NOT'))
table(res_limma$change)#NOT only
res_limma = res_limma[res_limma$change %in% c('UP', 'DOWN'), ]

RMdem_list[[i]] <- list(res_DESeq2, res_limma)
name = res@elementMetadata$description[2]
name = paste0(cl, '_',str_split(name, "groupList ", simplify = T)[,2])
name= paste0(gsub(' ','_',(name)))
names(RMdem_list[[i]] ) <- paste0(name, '_', c('DESeq2', 'limma'))
RMdep_list[[i]] <- list( merge_peak )
names(RMdep_list[[i]] ) <- name
c <- unique(condition[!condition %in% 'WT'])
countMatrix_list <- list(colnames(countMatrix)); names(countMatrix_list) <- c
condition_list <- list(condition); names(condition_list) <- c
contrast_IP_list[[i]] <- countMatrix_list
condition_IP_list[[i]] <- condition_list

############
n = 11
############
#GSE122961# with RIP, without Input 
############
print(RMPdeg_datasets[n,])
i = RMPdeg_datasets[n,'GSE']
cl = c('hESCs', 'HEK293T')
RMp = c('FTO', 'ALKBH5')
RM = 'm6A'
condition_list <- list(FTOKO = c(rep('FTOKO', 3), rep('WT', 3)), #for hESCs
                       ALKBH5OE = c(rep('ALKBH5OE', 3), rep('WT', 3)), #for HEK293T
                       FTOOE = c(rep('FTOOE', 3), rep('WT', 2))#for HEK293T
)
RMP_allDat[grep('GSE129842', RMP_allDat$study_alias ) 
           ,c("study_alias", "sample_alias", "sample_title")]
# View(allDat[[i]])
files <- list.files(file.path(data_directory,i))
files
setwd(data_directory)
setwd(i)
###################################
#         RIP-seq analysis        #
###################################
## Sample information
merge_peak <- read.delim("GSE129842_Table_S6.txt.gz",header=T)
countfiles <- list.files('./metaPlotR/')
countfiles <- countfiles[grep('.txt', countfiles)]
countfiles <- countfiles[! grepl('summary', countfiles)] #mouse embryonic fibroblasts (MEFs)
countfiles
resMeR_list <- vector("list", 3)
names(resMeR_list) <- names(condition_list)
countMatrix_list <- vector("list", length(condition_list))
names(countMatrix_list) <- names(condition_list)
resPeak_list <- vector("list",length(condition_list))
names(resPeak_list) <- names(condition_list)
for(c in names(condition_list)){
  if(c=='FTOKO'){
    countMatrix <- read.table( paste0("metaPlotR/m6A_counts_", c, "_hESCs.txt"),header=T)
  } else{
    countMatrix <- read.table( paste0("metaPlotR/m6A_counts_", c, "_HEK293T.txt"),header=T)
  }
  if(c=='FTOOE'){
    FTOWT <- read.table( paste0("metaPlotR/m6A_counts_FTOWT_HEK293T.txt"),header=T)
    countMatrix <- left_join(countMatrix, FTOWT[, !colnames(FTOWT) %in% c('Chr','Start','End','Strand','Length')], by = 'Geneid')
  }
  table(countMatrix$Geneid == merge_peak$name)
  table(duplicated(merge_peak$name))
  dim(countMatrix)
  colnames(countMatrix)
  rownames(countMatrix) <- countMatrix$Geneid
  countMatrix <- countMatrix[, !colnames(countMatrix) %in% c('Geneid','Chr','Start','End','Strand','Length')]
  # 1. 做好分组因子即可
  condition <- condition_list[[c]]
  groupList <- factor(condition)
  colData <- data.frame(row.names=colnames(countMatrix), groupList=groupList)
  # 2. Exploring Method 1: DESeq2-------------------------------------------
  # 构建 DESeqDataSet 对象
  dds <- DESeqDataSetFromMatrix(countData = countMatrix,
                                colData = colData,
                                design = ~ groupList)
  ## Pre-filtering
  keep <- rowSums(counts(dds)) >= 1
  table(keep)
  dds <- dds[keep,]
  # 接着计算每个样本的 library size 用于后续的标准化。library size 等于落在该样本所有peaks上的reads的总数，即 counts 矩阵中每列的加和。如果bFullLibrarySize设为TRUE，则会使用library的总reads数（根据原BAM/BED文件统计获得）。然后使用estimateDispersions估计统计分布，需要将参数fitType设为’local’
  dds<-estimateSizeFactors(dds)
  dds<-estimateDispersions(dds,fitType="local")
  dds <- nbinomWaldTest(dds)
  res <- results(dds)
  ## 提取你想要的差异分析结果，我们这里是treated组对untreated组进行比较
  res <-  results(dds, contrast=c("groupList", unique(condition)))
  res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
  length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
  res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                                 ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
  res_DESeq2 = as.data.frame(res)
  table(res_DESeq2$change)
  res_DESeq2 = res_DESeq2[res_DESeq2$change %in% c('UP', 'DOWN'), ]
  # 2. Exploring Method 2: limma-------------------------------------------
  library(edgeR)
  d0 <- DGEList(counts = countMatrix)
  dim(d0)
  # Filtering low count peaks
  cutoff <- 1           
  drop <- which(apply(cpm(d0), 1, max) < cutoff)
  print(length(drop))
  if (length(drop) >0) {
    d <- d0[-drop,] 
  } else{ d <- d0}
  dim(d) #number of peaks left: 8501      
  d <- calcNormFactors(d, method = c('TMM'))
  colnames(countMatrix)
  genotype <- groupList
  table(genotype)
  mm <- model.matrix(~0 + genotype)#
  v <- voom(d, mm, plot=F, normalize="quantile")#voom from limma:Transform count data to log2-counts per million (logtpm)
  fit <- lmFit(v, mm)
  head(coef(fit))
  contr <- makeContrasts(paste0(colnames(coef(fit)), collapse  = ' - '), levels = colnames(coef(fit))) #Change if other comparison
  tmp <- contrasts.fit(fit, contr)
  tmp <- eBayes(tmp)
  top.table <- topTable(tmp, sort.by = "P", n = Inf)
  head(top.table, 20)
  length(which(top.table$adj.P.Val < 0.05 & abs(top.table$logFC) > log2(1.5))) #how many are significant?
  res_limma <- as.data.frame(top.table)
  res_limma$change <- with(res_limma, ifelse(res_limma$adj.P.Val < 0.05 & abs(res_limma$logFC) > log2(1.5),
                                             ifelse(res_limma$logFC > log2(1.5),'UP','DOWN'),'NOT'))
  table(res_limma$change)
  res_limma = res_limma[res_limma$change %in% c('UP', 'DOWN'), ]
  resMeR_list[[c]] <- list(res_DESeq2, res_limma)
  countMatrix_list[[c]] <- countMatrix
  resPeak_list[[c]] <- merge_peak
  name = res@elementMetadata$description[2]
  name = str_split(name, "groupList ", simplify = T)[,2]
  name = paste0(gsub(' ','_',(name)))
  names(resMeR_list[[c]] ) <- paste0(name, '_', c('DESeq2', 'limma'))
}
RMdem_list[[i]] <- as.list(unlist(resMeR_list, recursive = FALSE))
name <- str_split(names(RMdem_list[[i]] ), "[.]", simplify = T)[,2]
print(cl)
name = paste0('HEK293T', '_', name)
name = gsub('HEK293T_FTOKO', 'hESCs_FTOKO', name)
names(RMdem_list[[i]] ) <- name
RMdep_list[[i]] <- resPeak_list 
name <- unique(gsub('(_limma|_DESeq2)', '', name))
names(RMdep_list[[i]] ) <- name

contrast_IP_list[[i]] <- lapply(countMatrix_list, extract_names)
names(condition_list) <-  gsub('_vs_WT', '', names(RMdep_list[[i]]))
names(contrast_IP_list[[i]]) <- names(condition_list) 
condition_IP_list[[i]] <- condition_list
# ###################################################################
# # remove duplicated genes (according to the baseMean of DESeq2 results)
# ###################################################################
# target_res <- vector("list", 3)
# names(target_res) <- names(condition_list)
# for(c in names(condition_list)){
#   resMeR <- resMeR_list[[c]]
#   length(unique(resMeR$geneSymbol))
# # tmp=by(gene_summary,gene_summary$Symbol,function(x) rownames(x)[which.min(x$Pvalue_Adj)])
#   tmp=by(resMeR,resMeR$geneSymbol,function(x) rownames(x)[which.max(x$baseMean)])
#   probes = as.character(tmp)
#   length(probes) # 7785
#   target_res[[c]] <- resMeR[rownames(resMeR) %in% probes,]
# }

############
n = 12
############
#GSE122948# with miCLIP without Input
############
print(RMPdeg_datasets[n,]) # hg38
i = RMPdeg_datasets[n,'GSE']
cl = 'HEK293T'
RMp = 'PCIF1'
RM = 'm6Am'
condition = c(rep('PCIF1KO', 2), rep('WT', 1))
RMP_allDat[grep('GSE122948', RMP_allDat$study_alias ) 
           ,c("study_alias", "sample_alias", "sample_title")]
RMP_allDat[grepl('GSM3489018|GSM3489019|GSM3489020', RMP_allDat$sample_alias ) 
           ,c("study_alias", "sample_alias", "sample_title")]

# View(allDat[[i]])
files <- list.files(file.path(data_directory,i))
files
getwd()
setwd(data_directory)
setwd(i)
bedfiles <- list.files('./metaPlotR/')
# bedfiles <- bedfiles[grep('.bedgraph', bedfiles)]
# bedfiles <- bedfiles[! grepl('annot|MEF', bedfiles)]
bedfiles <- bedfiles[grep('.bw', bedfiles)]
bedfiles <- bedfiles[! grepl('bed.gz|annot|MEF', bedfiles)] #mouse embryonic fibroblasts (MEFs)
bedfiles
###################################
#      RNA/Ribo-seq analysis      #
###################################
df_in <- read.table(file.path(data_directory,i, 'GSE122948_riboseq_input.txt.gz'), header = T)
## Raw counts
# ## Gene-level variance-stabilized counts from DESeq2 for RNA-Seq samples are normalized expression values that have been transformed to stabilize the variance across samples. These counts are obtained through a variance stabilizing transformation (VST), which is a common step in RNA-Seq data analysis pipelines.
setwd(data_directory)
setwd(i)
#判断数据是否需要转换
ex <- df_in[,-1]
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
# Calculate LogC
LogC <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0) || (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
# Log2 transform if conditions are met
if (LogC) {
  # Convert elements less than or equal to 0 to NaN
  ex[ex <= 0] <- NaN
  print("log2 transform unfinished")
} else {
  print("log2 transform not needed")
}

# Exploring Method 1: DESeq2-------------------------------------------
##Gene count matrix
count=df_in[,2:ncol(df_in)]
condition = c(rep('PCIF1KO', 3),rep('WT', 2))
## Gene name
library(data.table)
# id_mapped <- read.table("/data2/rluo4/hg38/gencode.v22.annotation.gene.probeMap",header = T)
# id_mapped$id <- str_split(id_mapped$id, '[.]', simplify = T)[,1]
id_mapped <- read.table("/data2/rluo4/hg38/gencode.v36.annotation.gene.probeMap",header = T)
dim(id_mapped) #  60675     6
table(str_split(df_in$gene_name, '[.]', simplify = T)[,1] %in% gene_info$geneEnsembls)
table(df_in$gene_name %in% id_mapped$id)
table(str_split(df_in$gene_name, '[.]', simplify = T)[,1] %in% str_split(id_mapped$id, '[.]', simplify = T)[,1])

count$id <- str_split(df_in$gene_name, '[.]', simplify = T)[,1]
# id_mapped <- unique(id_mapped) #%>% dplyr::select(id,gene) # 去重
id_mapped <- id_mapped[which(id_mapped$id %in% count$id),]
dim(id_mapped) #57930     6
table(duplicated(id_mapped$gene))
countGTF <- left_join(id_mapped,count,by="id")
countGTF <- aggregate(x = countGTF[,(ncol(id_mapped)+1):ncol(countGTF)],   #此时exprSet的第三列开始是表达矩阵内容
                      by = list(symbol = countGTF$gene),   #按照相同symbol分组，在组内计算
                      FUN = max) %>%   #原文中是计算最大值（max），也可以计算平均值（mean）或者中位数（median）
  column_to_rownames(var = 'symbol')
dim(countGTF) #56752     5
## Sample name
colnames(countGTF)<-gsub('riboseq_input_','',colnames(countGTF))
## Sample information
sample=data.frame(
  sample=colnames(countGTF),
  condition=c(rep("WT",2),rep("PCIF1KO",3)))
# dds <- DESeqDataSetFromMatrix(...)
# dds <- estimateSizeFactors(dds)
# # normalized expression
# normalized_counts <- counts(dds, normalized = TRUE)
# head(normalized_counts)
## SummarizedExperiment input
dds<-DESeqDataSetFromMatrix(countData=countGTF,colData=sample,design=~condition)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 10
table(keep)
dds <- dds[keep,]
## Note on factor levels
dds$condition <- factor(dds$condition, levels = c("PCIF1KO","WT"))
dds$condition <- relevel(dds$condition, ref = "WT")
# Differential expression analysis-------------------------------------------
dds <- DESeq(dds)
# res<-results(dds,name = "condition_PCIF1KO_vs_WT")
res <-  results(dds, contrast=c("condition",unique(condition))) # c("condition", rev(levels(genotype)))
# rev(levels(Group))
# res$label<-with(res,ifelse(res$baseMean>10&abs(res$log2FoldChange)>0.58&res$padj<0.05,"sig","nosig"))
# ress=as.data.frame(res)
# res1=subset(ress,ress$label=="sig")
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
table(res$change)
resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
resInP$geneSymbol <- rownames(resInP)
name = res@elementMetadata$description[2]
name = paste0(cl, '_',str_split(name, "condition ", simplify = T)[,2])
RMdeg_list[[i]] <- list( resInP)
names(RMdeg_list[[i]] ) <- paste0(gsub(' ','_',(name)))

c <- unique(condition[!condition %in% 'WT'])
countGTF_list <- list(colnames(countGTF)); names(countGTF_list) <- c
condition_list <- list(condition); names(condition_list) <- c
contrast_input_list[[i]] <- countGTF_list
condition_input_list[[i]] <- condition_list
# # Exploring Method 2: edgeR/limma-------------------------------------------
# library(edgeR)
# # 1. Formatting tables as needed for downstream analysis
# d0 <- DGEList(counts = countGTF)
# # Filtering low count genes
# cutoff <- 1
# drop <- which(apply(cpm(d0), 1, max) < cutoff)
# length(drop)
# d <- d0[-drop,]
# d <- calcNormFactors(d)
# N <- colSums(countGTF)
# # nf <- calcNormFactors(ercc_expr, lib.size=N) #calculating normalizing factors (nf) based on ERCCs
# dim(d) #number of genes left: 17454
# colnames(countGTF)
# genotype <- c(rep(c('WT'),2), rep(c('PCIF1KO'),3)) %>% as.factor()##c(rep(c("WT", "KO"), 6))
# table(genotype)
# # group <- interaction(ttmt, genotype)
# # 2. Differential gene analysis with edgeR
# # Building the model for limma-voom
# # keep <- rowSums(cpm(d0)>1) >= 2#至少在两个样本里cpm大于1
# # d0$samples$lib.size <- colSums(d0$counts)
# # y <- d0[keep, , keep.lib.sizes=FALSE]
# # y <- calcNormFactors(y)
# # #y$samples
# # logcpm <- cpm(y, prior.count=2, log=TRUE)
# # logcpm[1:5,1:5]
# # dim(logcpm)
# mm <- model.matrix(~0 + genotype)#
# # colnames(mm)=levels(factor(genotype))
# # rownames(mm)=colnames(countGTF)
# # y <- voom(d, mm, lib.size = N * nf, plot = T)  #lib.size adjusted by nf for CNN purposes
# v <- voom(d, mm, plot=T, normalize="quantile")#voom from limma:Transform count data to log2-counts per million (logtpm)
# #v <- voom(expr,design,plot=FALSE, normalize="quantile")#:Error in RStudioGD() Shadow graphics device error: r error 4 (R code execution error)
# fit <- lmFit(v, mm)
# head(coef(fit))
# contr <- makeContrasts(genotypePCIF1KO - genotypeWT, levels = colnames(coef(fit))) #Change if other comparison
# tmp <- contrasts.fit(fit, contr)
# tmp <- eBayes(tmp)
# top.table <- topTable(tmp, sort.by = "P", n = Inf)
# head(top.table, 20)
# length(which(top.table$adj.P.Val < 0.05 & abs(top.table$logFC) > log2(1.5))) #how many are significant?
# # KOtop.table <- top.table
# # KOtop.table$GeneID <- row.names(KOtop.table)
# # logCPM <- data.frame(v$E) #saving normalized counts for all samples
# # top.table <- full_join(WTtop.table, KOtop.table, by="GeneID") #merging the 2 top tables
# res_limma <- as.data.frame(top.table)
# res_limma$change <- with(res_limma, ifelse(res_limma$adj.P.Val < 0.05 & abs(res_limma$logFC) > log2(1.5),
#                                            ifelse(res_limma$logFC > log2(1.5),'UP','DOWN'),'NOT'))
# table(res_limma$change)
# res_1 = as.data.frame(res_limma[res_limma$change %in% c('UP', 'DOWN'), ])
# table(res_1$change)
# ##or
# y <- d0[-drop,]
# # design<-model.matrix(~genotype)
# y <- calcNormFactors(y) # # 估计dispersion
# # # edgeR内部进行了以下三步调用，有兴趣的同学可以查阅文档，看一看它们分别在做什么事情 # y <- estimateGLMCommonDisp(y,design = design) # y <- estimateGLMTrendedDisp(y,design = design) # y <- estimateGLMTagwiseDisp(y,design = design)
# # # 拟合广义线性模型
# design <- model.matrix(~0+genotype)
# rownames(design) <- colnames(countGTF)
# # colnames(design) <- levels(genotype)
# y <- estimateDisp(y,design = design)
# fit <- glmFit(y, design = design)
# result <- glmLRT(fit, contrast=c(1, -1))#c(-1,1)) -1-normal 1-tumor
# # 似然比检验 # coef = 2指的是对design矩阵的第二列（即是否照光）对应的系数进行检验
# # # result <- glmLRT(fit, coef=c(1:ncol(fit$design))) #glmLRT(fit,coef=2)
# # GLM testing for differential expression
# # Just as we used a GLM to fit the trend line above, we can also use this in finding the tags that are interesting by using a likelihood ratio test.
# # design.mat
# # ##   MM WM WW
# # ## 1  1  0  0
# # ## 2  1  0  0
# # ## 3  0  1  0
# # ## 4  0  1  0
# # ## 5  0  0  1
# # ## 6  0  0  1
# # ## attr(,"assign")
# # ## [1] 1 1 1
# # ## attr(,"contrasts")
# # ## attr(,"contrasts")$`d$samples$group`
# # ## [1] "contr.treatment"
# # fit <- glmFit(d2, design.mat)
# # # compare (group 1 - group 2) to 0:
# # # this is equivalent to comparing group 1 to group 2
# # lrt12 <- glmLRT(fit, contrast=c(1,-1,0))
# # lrt13 <- glmLRT(fit, contrast=c(1,0,-1))
# # lrt23 <- glmLRT(fit, contrast=c(0,1,-1))
# # topTags(lrt12, n=10)
# ## Coefficient:  1*MM -1*WM
# top<- topTags(result,n=Inf,adjust.method="BH", sort.by="logFC")
# topTags(result, n=10)
# degs<-top@.Data[[1]]
# summary(de <- decideTestsDGE(result))
# res_edgeR <- as.data.frame(top)
# res_edgeR$change <- with(res_edgeR, ifelse(res_edgeR$FDR < 0.05 & abs(res_edgeR$logFC) > log2(1.5),
#                                            ifelse(res_edgeR$logFC > log2(1.5),'UP','DOWN'),'NOT'))
# table(res_edgeR$change)
# res_2 = as.data.frame(res_edgeR[res_edgeR$change %in% c('UP', 'DOWN'), ])
# table(res_2$change)

###################################
#         RIP-seq analysis        #
###################################
## Sample information
condition = c(rep('PCIF1KO', 2), rep('WT', 1))
merge_peak <- read.delim("metaPlotR/newPeakFile.txt",header=F)
countMatrix <- read.table("metaPlotR/m6Am_counts.txt",header=T)
table(countMatrix$Geneid == merge_peak$V1[-1])
table(duplicated(merge_peak$V1[-1]))
dim(countMatrix)
colnames(countMatrix)
rownames(countMatrix) <- countMatrix$Geneid
countMatrix <- countMatrix[, !colnames(countMatrix) %in% c('Geneid','Chr','Start','End','Strand','Length')]
# 1. 做好分组因子即可
groupList<-factor(condition)
colData <- data.frame(row.names=colnames(countMatrix), groupList=groupList)
# 2. Exploring Method 1: DESeq2-------------------------------------------
# 构建 DESeqDataSet 对象
dds <- DESeqDataSetFromMatrix(countData = countMatrix,
                              colData = colData,
                              design = ~ groupList)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 1
table(keep)
dds <- dds[keep,]
# 接着计算每个样本的 library size 用于后续的标准化。library size 等于落在该样本所有peaks上的reads的总数，即 counts 矩阵中每列的加和。如果bFullLibrarySize设为TRUE，则会使用library的总reads数（根据原BAM/BED文件统计获得）。然后使用estimateDispersions估计统计分布，需要将参数fitType设为’local’
dds<-estimateSizeFactors(dds)
dds<-estimateDispersions(dds,fitType="local")
dds <- nbinomWaldTest(dds)
res <- results(dds)
## 提取你想要的差异分析结果，我们这里是treated组对untreated组进行比较
res <-  results(dds, contrast=c("groupList", unique(condition)))
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
res_DESeq2 = as.data.frame(res)
table(res_DESeq2$change)
res_DESeq2 = res_DESeq2[res_DESeq2$change %in% c('UP', 'DOWN'), ]
# 2. Exploring Method 2: limma-------------------------------------------
library(edgeR)
d0 <- DGEList(counts = countMatrix)
dim(d0)
# Filtering low count peaks
cutoff <- 1           
drop <- which(apply(cpm(d0), 1, max) < cutoff)
print(length(drop))
if (length(drop) >0) {
  d <- d0[-drop,] 
} else{ d <- d0}
dim(d) #number of peaks left: 8501      
d <- calcNormFactors(d, method = c('TMM'))
colnames(countMatrix)
genotype <- groupList
table(genotype)
mm <- model.matrix(~0 + genotype)#
v <- voom(d, mm, plot=F, normalize="quantile")#voom from limma:Transform count data to log2-counts per million (logtpm)
fit <- lmFit(v, mm)
head(coef(fit))
contr <- makeContrasts(paste0(colnames(coef(fit)), collapse  = ' - '), levels = colnames(coef(fit))) #Change if other comparison
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)
length(which(top.table$adj.P.Val < 0.05 & abs(top.table$logFC) > log2(1.5))) #how many are significant?
res_limma <- as.data.frame(top.table)
res_limma$change <- with(res_limma, ifelse(res_limma$adj.P.Val < 0.05 & abs(res_limma$logFC) > log2(1.5),
                                           ifelse(res_limma$logFC > log2(1.5),'UP','DOWN'),'NOT'))
table(res_limma$change)#NOT only
res_limma = res_limma[res_limma$change %in% c('UP', 'DOWN'), ]

RMdem_list[[i]] <- list(res_DESeq2, res_limma)
name = res@elementMetadata$description[2]
name = paste0(cl, '_',str_split(name, "groupList ", simplify = T)[,2])
name= paste0(gsub(' ','_',(name)))
names(RMdem_list[[i]] ) <- paste0(name, '_', c('DESeq2', 'limma'))
RMdep_list[[i]] <- list( merge_peak )
names(RMdep_list[[i]] ) <- name

c <- unique(condition[!condition %in% 'WT'])
countMatrix_list <- list(colnames(countMatrix)); names(countMatrix_list) <- c
condition_list <- list(condition); names(condition_list) <- c
contrast_IP_list[[i]] <- countMatrix_list
condition_IP_list[[i]] <- condition_list

############
n = 13
############
#GSE121942# with RIP and Input
############
print(RMPdeg_datasets[n,])# hg38 (paper used: hg19; I used: hg38 on server 57)
i = RMPdeg_datasets[n,'GSE']
cl = 'HepG2'
RMp = 'METTL3'
RM = 'm6A'
ENA <- as.data.frame(RMP_allDat[grepl('GSE110320', RMP_allDat$study_alias ) 
                                ,c("study_alias", "sample_alias", "sample_title")])
ENA[order(ENA$sample_alias),]
skip <- RMP_allDat[grepl('GSE110320', RMP_allDat$study_alias ) & grepl('input', RMP_allDat$sample_title)
                   , "sample_alias"]
skip <- RMP_allDat[grepl('GSE110320', RMP_allDat$study_alias ) & grepl('m6a-', RMP_allDat$sample_title)
                   , "sample_alias"]
paste0(skip$sample_alias, collapse = " ")
RMP_allDat[grepl('GSE110320', RMP_allDat$study_alias ) & grepl('m6a-', RMP_allDat$sample_title)
           ,c("study_alias", "sample_alias", "sample_title")]
# View(allDat[[i]])
condition_list <- list(METTL3KD = c(rep('METTL3KD', 3), rep('WT', 3)),
                       METTL14KD = c(rep('METTL14KD', 3), rep('WT', 3)),
                       SETD2KD = c(rep('SETD2KD', 3), rep('WT', 3)),
                       WTAPKD = c(rep('WTAPKD', 3), rep('WT', 3)),
                       SETD2_2KD = c(rep('SETD2KD', 3), rep('WT', 3)))

# View(allDat[[i]])
files <- list.files(file.path(data_directory,i))
files
getwd()
setwd(data_directory)
setwd(i)
ENA <- ENA[order(ENA$sample_alias),]
y <- ENA[grepl('shCtrl input',ENA$sample_title),]
y <- ENA[grepl('shCtrl m6a',ENA$sample_title),]
y <- ENA[grepl('shSETD2#2) input',ENA$sample_title),]
y <- ENA[grepl('shSETD2#2) m6a',ENA$sample_title),]

y <- ENA[grepl('shMETTL3 input',ENA$sample_title),]
y <- ENA[grepl('shMETTL14 input',ENA$sample_title),]
y <- ENA[grepl('shSETD2#1 input',ENA$sample_title),]
y <- ENA[grepl('shWTAP input',ENA$sample_title),]
y <- ENA[grepl('shSETD2#2 input',ENA$sample_title),]
y
x <- paste0('$outdir/GSE110320/', y$sample_alias, '_input.bam')
paste(x, collapse = ' ')
countfiles <- list.files('./metaPlotR/')
countfiles <- countfiles[grep('counts', countfiles)]
countfiles <- countfiles[! grepl('summary', countfiles)] 
countfiles
# ###################################
# #      RNA/Ribo-seq analysis      #
# ###################################
# # 1) RNA-seq 
# # GSE121949	RNA stability profiling in SETD2 and METTL14 knockdown HepG2 cells
# # GSE121952	Active ribosome profiling in SETD2 or METTL14 knockdown HepG2 cells
# df_in <- read.table(file.path(data_directory,i, 'GSE121949_RNAseq_Raw_and_RPKM.txt.gz'), header = T)
# # df_in <- df_in[, ! grepl('RAW', colnames(df_in))]
# df_in <- df_in[, ! grepl('RPKM', colnames(df_in))]
# ###################################
# #   more than 1 condition here   #
# ###################################
# resInP_list <- list()
# ##Gene count matrix
# countGTF=df_in[,c(9:ncol(df_in))]
# ## Gene name
# table(duplicated(df_in$symbol))
# dim(countGTF) #56752     5
# ## Sample name
# colnames(countGTF)#<-gsub('riboseq_input_','',colnames(countGTF))
# ## Sample information
# sample_meta=data.frame(
#   sample=colnames(countGTF),
#   condition=c(rep("WT",8),rep("METTL14KD",8), rep("SETD2KD",8)))
# #####condition_SETD2KD_vs_WT#####
# count <- countGTF[, !grepl('MTL14', colnames(countGTF))]
# sample <- sample_meta[sample_meta$sample %in% colnames(count), ]
# ## SummarizedExperiment input
# condition = c(rep('SETD2KD', 8), rep('WT', 8))
# dds<-DESeqDataSetFromMatrix(countData=count, colData=sample,design=~condition)
# ## Pre-filtering
# keep <- rowSums(counts(dds)) >= 10
# table(keep)
# dds <- dds[keep,]
# ## Note on factor levels
# dds$condition <- factor(dds$condition, levels = c("SETD2KD","WT"))
# dds$condition <- relevel(dds$condition, ref = "WT")
# dds <- DESeq(dds)
# res<-results(dds,name = "condition_SETD2KD_vs_WT")
# res <-  results(dds, contrast=c("condition",unique(condition))) # c("condition", rev(levels(genotype)))
# length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
# res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
# res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
#                                ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
# table(res$change)
# res$geneSymbol <- df_in$symbol[keep]
# resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
# table(duplicated(resInP$geneSymbol))
# # tmp=by(resInP,resInP$geneSymbol,function(x) rownames(x)[which.max(x$baseMean)])
# # probes = as.character(tmp)
# # length(probes) # 7785
# # resInP <- resInP[rownames(resInP) %in% probes,]
# resInP_list[['SETD2KD_vs_WT']] <- resInP
# #####condition_METTL14KD_vs_WT#####
# count <- countGTF[, !grepl('SETD2', colnames(countGTF))]
# sample <- sample_meta[sample_meta$sample %in% colnames(count), ]
# ## SummarizedExperiment input
# condition = c(rep('METTL14KD', 8), rep('WT', 8))
# dds<-DESeqDataSetFromMatrix(countData=count, colData=sample,design=~condition)
# ## Pre-filtering
# keep <- rowSums(counts(dds)) >= 10
# table(keep)
# dds <- dds[keep,]
# ## Note on factor levels
# dds$condition <- factor(dds$condition, levels = c("METTL14KD","WT"))
# dds$condition <- relevel(dds$condition, ref = "WT")
# dds <- DESeq(dds)
# res<-results(dds,name = "condition_METTL14KD_vs_WT")
# res <-  results(dds, contrast=c("condition",unique(condition))) # c("condition", rev(levels(genotype)))
# length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
# res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
# res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
#                                ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
# table(res$change)
# res$geneSymbol <- df_in$symbol[keep]
# resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
# table(duplicated(resInP$geneSymbol))
# # tmp=by(resInP,resInP$geneSymbol,function(x) rownames(x)[which.max(x$baseMean)])
# # probes = as.character(tmp)
# # length(probes) # 7785
# # resInP <- resInP[rownames(resInP) %in% probes,]
# resInP_list[['METTL14KD_vs_WT']] <- resInP
# 
# RMdeg_list[[i]] <- resInP_list
# name <- paste0(cl,'_', names(resInP_list))
# names(RMdeg_list[[i]] ) <- name

###################################
# 2) Ribo-seq (with input and ribo, not knowing how to analyze it)
###################################
df_in <- read.delim(file.path(data_directory,i, 'GSE121952_HepG2_ribosome-profiling_featureCount_results.txt.gz'), header = T)
colnames(df_in)
###################################
#         RIP-seq analysis        #
###################################
peaks <- read.table(file.path(data_directory,i,'GSE110320_HepG2_shMettl3_m6A-seq_rep1-3_exomePeak_peak.txt.gz'),header = T)
# peaks <- read.table(file.path(data_directory,i,'metaPlotR/m6Am_counts_METTL3KD.txt'),header = T)
# peak score (MACS)
# Fold change (MACS)
# -log10 p-value (MACS) 
# -log10 q-value (FDR adjustment, MACS) 
summary(peaks$lg.fdr)
summary(10^(peaks$lg.fdr))
###################################
#  Cell number-normalized (CNN)   #
#      MeRIP INPUT analysis       #
###################################
resInP_list <- vector("list", length(condition_list))
names(resInP_list) <- names(condition_list)
countGTF_list <- vector("list", length(condition_list))
names(countGTF_list) <- names(condition_list)
for(c in names(condition_list)){
  countGTF <- read.table( paste0("metaPlotR/Input_counts_", c, ".txt"),header=T)
  dim(countGTF)
  colnames(countGTF)
  rownames(countGTF) <- countGTF$Geneid
  countGTF <- countGTF[, !colnames(countGTF) %in% c('Geneid','Chr','Start','End','Strand','Length')]
  dim(countGTF) # 28265     4
  condition <- condition_list[[c]]
  ## Sample information
  sample=data.frame(
    sample=colnames(countGTF),
    condition=condition)
  ## SummarizedExperiment input
  dds<-DESeqDataSetFromMatrix(countData=countGTF,colData=sample,design=~condition)
  ## Pre-filtering
  keep <- rowSums(counts(dds)) >= 10
  table(keep)
  dds <- dds[keep,]
  ## Note on factor levels
  dds$condition <- factor(dds$condition, levels = unique(condition))
  dds$condition <- relevel(dds$condition, ref = "WT")
  # Differential expression analysis-------------------------------------------
  dds <- DESeq(dds)
  res <-  results(dds, contrast=c("condition",unique(condition)))
  res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
  length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
  res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                                 ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
  table(res$change)
  table(is.na(res$change))
  # resOrdered <- as.data.frame(res[order(res$padj), ])
  resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
  resInP$geneSymbol <- rownames(resInP)
  table(resInP$change)
  resInP_list[[c]] <- resInP
  countGTF_list[[c]] <- countGTF
}
name = paste0(names(resInP_list), ' vs WT')
name = paste0(cl, '_', name)
name <- gsub('HepG2_SETD2_2KD vs WT', 'HepG2.2_SETD2KD vs WT', name)
RMdeg_list[[i]] <- resInP_list 
names(RMdeg_list[[i]] ) <- paste0(gsub(' ','_',(name)))

contrast_input_list[[i]] <- lapply(countGTF_list, extract_names)
condition_input_list[[i]] <- condition_list
names(contrast_input_list[[i]]) <- gsub('_vs_WT', '',names(RMdeg_list[[i]] ) )
names(condition_input_list[[i]]) <- gsub('_vs_WT', '',names(RMdeg_list[[i]] ) )

###################################
#  Cell number-normalized (CNN)   #
#   MeRIP RIP analysis - Part1    #
###################################
resMeR_list <- vector("list", length(condition_list))
names(resMeR_list) <- names(condition_list)
countMatrix_list <- vector("list", length(condition_list))
names(countMatrix_list) <- names(condition_list)
resPeak_list <- vector("list",length(condition_list))
names(resPeak_list) <- names(condition_list)

for(c in names(condition_list)){
  merge_peak <- read.delim(paste0("metaPlotR/", c, "_newPeakFile.txt"),header=F)
  countMatrix <- read.table(paste0("metaPlotR/m6A_counts_", c, ".txt"),header=T)
  table(countMatrix$Geneid == merge_peak[-1,1])
  table(duplicated(merge_peak[-1,1]))
  dim(countMatrix)
  colnames(countMatrix)
  rownames(countMatrix) <- countMatrix$Geneid
  countMatrix <- countMatrix[, !colnames(countMatrix) %in% c('Geneid','Chr','Start','End','Strand','Length')]
  # 1. 做好分组因子即可
  condition <- condition_list[[c]]
  groupList<-factor(condition)
  colData <- data.frame(row.names=colnames(countMatrix), groupList=groupList)
  # 2. Exploring Method 1: DESeq2-------------------------------------------
  # 构建 DESeqDataSet 对象
  dds <- DESeqDataSetFromMatrix(countData = countMatrix,
                                colData = colData,
                                design = ~ groupList)
  ## Pre-filtering
  keep <- rowSums(counts(dds)) >= 1
  table(keep)
  dds <- dds[keep,]
  # 接着计算每个样本的 library size 用于后续的标准化。library size 等于落在该样本所有peaks上的reads的总数，即 counts 矩阵中每列的加和。如果bFullLibrarySize设为TRUE，则会使用library的总reads数（根据原BAM/BED文件统计获得）。然后使用estimateDispersions估计统计分布，需要将参数fitType设为’local’
  dds<-estimateSizeFactors(dds)
  dds<-estimateDispersions(dds,fitType="local")
  dds <- nbinomWaldTest(dds)
  res <- results(dds)
  ## 提取你想要的差异分析结果，我们这里是treated组对untreated组进行比较
  res <-  results(dds, contrast=c("groupList", unique(condition)))
  res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
  length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
  res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                                 ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
  res_DESeq2 = as.data.frame(res)
  table(res_DESeq2$change)
  res_DESeq2 = res_DESeq2[res_DESeq2$change %in% c('UP', 'DOWN'), ]
  # 2. Exploring Method 2: limma-------------------------------------------
  library(edgeR)
  d0 <- DGEList(counts = countMatrix)
  dim(d0)
  # Filtering low count peaks
  cutoff <- 1           
  drop <- which(apply(cpm(d0), 1, max) < cutoff)
  print(length(drop))
  if (length(drop) >0) {
    d <- d0[-drop,] 
  } else{ d <- d0}
  dim(d) #number of peaks left: 8501      
  d <- calcNormFactors(d, method = c('TMM'))
  colnames(countMatrix)
  genotype <- groupList
  table(genotype)
  mm <- model.matrix(~0 + genotype)#
  v <- voom(d, mm, plot=F, normalize="quantile")#voom from limma:Transform count data to log2-counts per million (logtpm)
  fit <- lmFit(v, mm)
  head(coef(fit))
  contr <- makeContrasts(paste0(colnames(coef(fit)), collapse  = ' - '), levels = colnames(coef(fit))) #Change if other comparison
  tmp <- contrasts.fit(fit, contr)
  tmp <- eBayes(tmp)
  top.table <- topTable(tmp, sort.by = "P", n = Inf)
  head(top.table, 20)
  length(which(top.table$adj.P.Val < 0.05 & abs(top.table$logFC) > log2(1.5))) #how many are significant?
  res_limma <- as.data.frame(top.table)
  res_limma$change <- with(res_limma, ifelse(res_limma$adj.P.Val < 0.05 & abs(res_limma$logFC) > log2(1.5),
                                             ifelse(res_limma$logFC > log2(1.5),'UP','DOWN'),'NOT'))
  table(res_limma$change)
  res_limma = res_limma[res_limma$change %in% c('UP', 'DOWN'), ]
  resMeR_list[[c]] <- list(res_DESeq2, res_limma)
  countMatrix_list[[c]] <- countMatrix
  resPeak_list[[c]] <- merge_peak
  name = res@elementMetadata$description[2]
  name = paste0(cl, '_',str_split(name, "groupList ", simplify = T)[,2])
  name = paste0(gsub(' ','_',(name)))
  names(resMeR_list[[c]] ) <- paste0(name, '_', c('DESeq2', 'limma'))
}
names(resMeR_list)
RMdem_list[[i]] <- as.list(unlist(resMeR_list, recursive = FALSE))
name <- gsub('SETD2_2KD.HepG2', 'SETD2KD.HepG2_2', names(RMdem_list[[i]] ))
name <- str_split(name, "[.]", simplify = T)[,2]
name <- gsub('HepG2_2', 'HepG2.2', name)
names(RMdem_list[[i]] ) <- name
RMdep_list[[i]] <- resPeak_list 
name <- unique(gsub('(_limma|_DESeq2)', '', name))
names(RMdep_list[[i]] ) <- name
contrast_IP_list[[i]] <- lapply(countMatrix_list, extract_names)
condition_IP_list[[i]] <- condition_list
names(contrast_IP_list[[i]]) <- gsub('_vs_WT', '',names(RMdep_list[[i]] ) )
names(condition_IP_list[[i]]) <- gsub('_vs_WT', '',names(RMdep_list[[i]] ) )

# ###################################
# #  Cell number-normalized (CNN)   #
# #   MeRIP RIP analysis - Part2    #
# ###################################
# target_res <- vector("list", 2)
# names(target_res) <- names(condition_list)
# for(c in names(condition_list)){  
#   library(tidyverse)
#   library(edgeR)
#   countGTF <- countGTF_list[[c]]
#   # 1. Formatting tables as needed for downstream analysis
#   d0 <- DGEList(counts = countGTF)
#   # Filtering low count genes
#   cutoff <- 1           
#   drop <- which(apply(cpm(d0), 1, max) < cutoff)
#   length(drop)
#   d <- d0[-drop,] 
#   d <- calcNormFactors(d)
#   N <- colSums(countGTF)
#   colnames(countGTF)
#   genotype <- condition %>% as.factor()
#   mm <- model.matrix(~0 + genotype)
#   v <- voom(d, mm, plot=T, normalize="quantile")#voom from limma:Transform count data to log2-counts per million (logtpm)
#   Input <- as.data.frame(v$E)
#   Input$geneSymbol <- rownames(Input)
#   
#   countMatrix <- countMatrix_list[[c]]
#   d0 <- DGEList(countMatrix)
#   N <- colSums(countMatrix)  
#   cutoff <- 1
#   drop <- which(apply(cpm(d0), 1, max) < cutoff)
#   length(drop)
#   if(length(drop)>0){
#     d <- d0[-drop,] 
#   }else{
#     d <- d0 
#   }
#   d <- calcNormFactors(d)
#   genotype <- condition %>% as.factor()##c(rep(c("WT", "KO"), 6))
#   mm <- model.matrix(~0 + genotype)#
#   v <- voom(d, mm, plot=T, normalize="quantile")#voom from limma:Transform count data to log2-counts per million (logtpm)
#   
#   RIP <- data.frame(v$E)
#   RIP$peak_id <- rownames(RIP)
#   RIP <- inner_join(RIP, merge_bed.anno, by=c("peak_id"))
#   
#   Norm.RIP <- inner_join(RIP, Input, by="geneSymbol")
#   Norm.RIP$Ave.input <- log2(rowMeans(2^Norm.RIP[,colnames(Norm.RIP) %in% colnames(countGTF)]))
#   colnames(Norm.RIP) <- c(paste0("RIP", 1:length(condition)), 'peak_id','chr', 'start', 'end','geneSymbol',
#                           paste0("Input", 1:length(condition)), "Ave.input")
#   # Adjusting RIP values for gene expression
#   # All RIP values are corrected by dividing by the ratio (Input.sample/Input.average) (note: all values are in log2 format so this appears as "-Input#+Ave.Input")
#   Norm.RIP <- mutate(Norm.RIP, 
#                      Ratio1 = RIP1-Input1+Ave.input,
#                      Ratio2 = RIP2-Input2+Ave.input,
#                      Ratio3 = RIP3-Input3+Ave.input,
#                      Ratio4 = RIP4-Input4+Ave.input,
#   )
#   # Formatting the "Ratio" table (i.e. the adjusted RIP peaks used for the final analysis)
#   # No filtering is performed as both input and RIP data are already filtered for low counts
#   Ratio <- cbind(Norm.RIP[,c('chr', 'start', 'end','geneSymbol')], 2^(Norm.RIP[,grepl('Ratio', colnames(Norm.RIP))]) )
#   Ratio$RefID <- paste("Ref", 1:nrow(Ratio)) # used to bring back peak information later
#   rownames(Ratio) <- Ratio$RefID
#   # 2. Differential gene analysis with edgeR
#   # Building the model for limma-voom
#   # As the data is already fully normalized, further re-scaling is prevented by fixing the library size in the voom function
#   d0 <- DGEList(Ratio[, grepl('Ratio', colnames(Ratio))])
#   d <- calcNormFactors(d0)
#   N <- colSums(Ratio[,grepl('Ratio', colnames(Ratio))])
#   genotype <- condition %>% as.factor()
#   mm <- model.matrix(~0 + genotype)
#   y <- voom(d, mm, plot = T, lib.size = rep(10^6,length(condition)))  #lib.size fixed as 10^6 for all samples
#   # v <- voom(d, mm, plot=T, normalize="quantile")#voom from limma:Transform count data to log2-counts per million (logtpm)
#   fit <- lmFit(y, mm)
#   head(coef(fit))
#   # Differential analysis for inputs (Paused vs FBS)
#   contr <- makeContrasts(paste0(colnames(coef(fit)), collapse  = ' - '), levels = colnames(coef(fit))) #Change if other comparison
#   tmp <- contrasts.fit(fit, contr)
#   tmp <- eBayes(tmp)
#   top.table <- topTable(tmp, sort.by = "P", n = Inf)
#   head(top.table, 20)
#   length(which(top.table$adj.P.Val < 0.05 & abs(top.table$logFC) > log2(1.5))) #how many are significant?
#   # Bring back peak information and export
#   top.table$RefID <- rownames(top.table)
#   top.table <- inner_join(top.table, Ratio, by=c("RefID"))
#   # top.table <- top.table[,-c(7,18)]
#   # write.table(top.table, file = paste0(export_folder,"/m6A.CNN.toptable.txt"))
#   All <- inner_join(Norm.RIP, top.table[,!grepl('Ratio', colnames(top.table))], by=c('chr', 'start', 'end','geneSymbol'))
#   # write.table(All, file = paste0(export_folder,"/m6A.CNN.all_results.txt"))
#   All$change <- with(All, ifelse(All$adj.P.Val < 0.05 & abs(All$logFC) > log2(1.5),
#                                  ifelse(All$logFC > log2(1.5),'UP','DOWN'),'NOT'))
#   table(All$change)
#   All <- All[All$change %in% c('UP', 'DOWN'), ]
#   print(length(unique(All$geneSymbol)))# 67
#   table(All$geneSymbol %in% resMeR$geneSymbol)
#   table(All$geneSymbol %in% resInP$geneSymbol)
#   target_res[[c]] <- All
# }
############
n = 14
############
#GSE98623# with miCLIP and RNA-seq, without case RIP 
############
print(RMPdeg_datasets[n,])
i = RMPdeg_datasets[n,'GSE']
cl = 'MOLM13'
RMp = 'METTL3'
RM = 'm6A'
condition = c(rep('METTL3KD', 3), rep('WT', 2))
# Two independent shRNAs showed efficient knockdown of METTL3 (shRNA #9 and #12) (Fig. 1a) and resulted in a reduction of global m6A levels
RMP_allDat[grep('GSE98623', RMP_allDat$study_alias ) 
           ,c("study_alias", "sample_alias", "sample_title")]
# View(allDat[[i]])
files <- list.files(file.path(data_directory,i))
files
getwd()
setwd(data_directory)
setwd(i)
###################################
#      RNA/Ribo-seq analysis      #
###################################
# 1) RNA-seq
df_in <- read_excel(file.path(data_directory,i, 'GSE98623_rnaseq.counts.xls'), sheet = 1)
## Raw counts
#判断数据是否需要转换
ex <- df_in[,-(1:3)]
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
# Calculate LogC
LogC <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0) || (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
# Log2 transform if conditions are met
if (LogC) {
  # Convert elements less than or equal to 0 to NaN
  ex[ex <= 0] <- NaN
  print("log2 transform unfinished")
} else {
  print("log2 transform not needed")
}
# [1] "log2 transform unfinished" for RNA
###################################
#   more than 1 condition here   #
###################################
resInP_list <- list()
##Gene count matrix
countGTF=df_in[,c(4:ncol(df_in))]
## Gene name
table(duplicated(df_in[1]))
table(duplicated(df_in$external_gene_name))
dim(countGTF) #1867    9
## Sample name
colnames(countGTF)
## Sample information
sample_meta=data.frame(
  sample=colnames(countGTF),
  condition=c(rep("WT",3),rep("METTL3_9KD",3), rep("METTL3_12KD",3)))
#####condition_METTL3_9KD_vs_WT#####
count <- as.matrix(countGTF[, !grepl('m3kd12', colnames(countGTF))])
table(is.na(count))
sample <- sample_meta[sample_meta$sample %in% colnames(count), ]
## SummarizedExperiment input
condition = c(rep('METTL3_9KD', 3), rep('WT', 3))
dds<-DESeqDataSetFromMatrix(countData=count, colData=sample,design=~condition)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 10
table(keep)
dds <- dds[keep,]
## Note on factor levels
dds$condition <- factor(dds$condition, levels = c("METTL3_9KD","WT"))
dds$condition <- relevel(dds$condition, ref = "WT")
dds <- DESeq(dds)
res<-results(dds,name = "condition_METTL3_9KD_vs_WT")
res <-  results(dds, contrast=c("condition",unique(condition))) # c("condition", rev(levels(genotype)))
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
table(res$change)
res$geneSymbol <-df_in$external_gene_name
resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
table(duplicated(resInP$geneSymbol))
# tmp=by(resInP,resInP$geneSymbol,function(x) rownames(x)[which.max(x$baseMean)])
# probes = as.character(tmp)
# length(probes) # 7785
# resInP <- resInP[rownames(resInP) %in% probes,]
resInP_list[['METTL3_9KD_vs_WT']] <- resInP
#####condition_METTL3_12KD_vs_WT#####
count <- as.matrix(countGTF[, !grepl('m3kd9', colnames(countGTF))])
table(is.na(count))
sample <- sample_meta[sample_meta$sample %in% colnames(count), ]
## SummarizedExperiment input
condition = c(rep('METTL3_12KD', 3), rep('WT', 3))
dds<-DESeqDataSetFromMatrix(countData=count, colData=sample,design=~condition)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 10
table(keep)
dds <- dds[keep,]
## Note on factor levels
dds$condition <- factor(dds$condition, levels = c("METTL3_12KD","WT"))
dds$condition <- relevel(dds$condition, ref = "WT")
dds <- DESeq(dds)
res<-results(dds,name = "condition_METTL3_12KD_vs_WT")
res <-  results(dds, contrast=c("condition",unique(condition))) # c("condition", rev(levels(genotype)))
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
table(res$change)
res$geneSymbol <-df_in$external_gene_name
resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
table(duplicated(resInP$geneSymbol))
# tmp=by(resInP,resInP$geneSymbol,function(x) rownames(x)[which.max(x$baseMean)])
# probes = as.character(tmp)
# length(probes) # 7785
# resInP <- resInP[rownames(resInP) %in% probes,]
resInP_list[['METTL3_12KD_vs_WT']] <- resInP
RMdeg_list[[i]] <- resInP_list
name <- paste0(cl,'_', names(resInP_list))
name <- c('MOLM13.9_METTL3KD_vs_WT', 'MOLM13.12_METTL3KD_vs_WT')
names(RMdeg_list[[i]] ) <- name
names(condition_list) <- gsub('_vs_WT', '', name)
contrast_input_list[[i]] <- NULL
condition_input_list[[i]] <- condition_list

# ###################################
# # 2) Ribo-seq (with input and ribo, not knowing how to analyze it)
# df_in <- read_excel(file.path(data_directory,i, 'GSE98623_rpf.counts.xls'), sheet = 1)

############
n = 15
############
#GSE128699#  with miCLIP and RIP, without RNA-seq
############
print(RMPdeg_datasets[n,])# hg38
i = RMPdeg_datasets[n,'GSE']
cl = 'HCT116'
RMp = c('METTL5','ZCCHC4')
RM = 'm6A'
condition_list <- list(
  METTL5KO = c(rep('METTL5KO', 2), rep('WT', 2)),
  ZCCHC4KO = c(rep('ZCCHC4KO', 2), rep('WT', 2))
)
# View(allDat[[i]])
ENA <- as.data.frame(RMP_allDat[grepl('GSE128699', RMP_allDat$study_alias ) 
                                ,c("study_alias", "sample_alias", "sample_title")])
ENA[order(ENA$sample_alias),]
files <- list.files(file.path(data_directory,i))
files
getwd()
setwd(data_directory)
setwd(i)
ENA <- ENA[order(ENA$sample_alias),]
y <- ENA[grepl('METTL5|WT',ENA$sample_title),]
y <- ENA[grepl('ZCCHC4|WT',ENA$sample_title),]
y
x <- paste0('$outdir/GSE128699/', y$sample_alias, '_IP.bam')
paste(x, collapse = ' ')
countfiles <- list.files('./metaPlotR/')
countfiles <- countfiles[grep('counts', countfiles)]
countfiles <- countfiles[! grepl('summary', countfiles)] 
countfiles
###################################
#         RIP-seq analysis        #
###################################
resMeR_list <- vector("list", length(condition_list))
names(resMeR_list) <- names(condition_list)
countMatrix_list <- vector("list", length(condition_list))
names(countMatrix_list) <- names(condition_list)
resPeak_list <- vector("list",length(condition_list))
names(resPeak_list) <- names(condition_list)

for(c in names(condition_list)){
  merge_peak <- read.delim(paste0("metaPlotR/", c, "_newPeakFile.txt"),header=F)
  countMatrix <- read.table(paste0("metaPlotR/m6A_counts_", c, ".txt"),header=T)
  table(countMatrix$Geneid == merge_peak[-1,1])
  table(duplicated(merge_peak[-1,1]))
  dim(countMatrix)
  colnames(countMatrix)
  rownames(countMatrix) <- countMatrix$Geneid
  countMatrix <- countMatrix[, !colnames(countMatrix) %in% c('Geneid','Chr','Start','End','Strand','Length')]
  # 1. 做好分组因子即可
  condition <- condition_list[[c]]
  groupList <- factor(condition)
  colData <- data.frame(row.names=colnames(countMatrix), groupList=groupList)
  # 2. Exploring Method 1: DESeq2-------------------------------------------
  # 构建 DESeqDataSet 对象
  dds <- DESeqDataSetFromMatrix(countData = countMatrix,
                                colData = colData,
                                design = ~ groupList)
  ## Pre-filtering
  keep <- rowSums(counts(dds)) >= 1
  table(keep)
  dds <- dds[keep,]
  # 接着计算每个样本的 library size 用于后续的标准化。library size 等于落在该样本所有peaks上的reads的总数，即 counts 矩阵中每列的加和。如果bFullLibrarySize设为TRUE，则会使用library的总reads数（根据原BAM/BED文件统计获得）。然后使用estimateDispersions估计统计分布，需要将参数fitType设为’local’
  dds<-estimateSizeFactors(dds)
  dds<-estimateDispersions(dds,fitType="local")
  dds <- nbinomWaldTest(dds)
  res <- results(dds)
  ## 提取你想要的差异分析结果，我们这里是treated组对untreated组进行比较
  res <-  results(dds, contrast=c("groupList", unique(condition)))
  res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
  length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
  res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                                 ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
  res_DESeq2 = as.data.frame(res)
  table(res_DESeq2$change)
  res_DESeq2 = res_DESeq2[res_DESeq2$change %in% c('UP', 'DOWN'), ]
  # 2. Exploring Method 2: limma-------------------------------------------
  library(edgeR)
  d0 <- DGEList(counts = countMatrix)
  dim(d0)
  # Filtering low count peaks
  cutoff <- 1           
  drop <- which(apply(cpm(d0), 1, max) < cutoff)
  print(length(drop))
  if (length(drop) >0) {
    d <- d0[-drop,] 
  } else{ d <- d0}
  dim(d) #number of peaks left: 8501      
  d <- calcNormFactors(d, method = c('TMM'))
  colnames(countMatrix)
  genotype <- groupList
  table(genotype)
  mm <- model.matrix(~0 + genotype)#
  v <- voom(d, mm, plot=F, normalize="quantile")#voom from limma:Transform count data to log2-counts per million (logtpm)
  fit <- lmFit(v, mm)
  head(coef(fit))
  contr <- makeContrasts(paste0(colnames(coef(fit)), collapse  = ' - '), levels = colnames(coef(fit))) #Change if other comparison
  tmp <- contrasts.fit(fit, contr)
  tmp <- eBayes(tmp)
  top.table <- topTable(tmp, sort.by = "P", n = Inf)
  head(top.table, 20)
  length(which(top.table$adj.P.Val < 0.05 & abs(top.table$logFC) > log2(1.5))) #how many are significant?
  res_limma <- as.data.frame(top.table)
  res_limma$change <- with(res_limma, ifelse(res_limma$adj.P.Val < 0.05 & abs(res_limma$logFC) > log2(1.5),
                                             ifelse(res_limma$logFC > log2(1.5),'UP','DOWN'),'NOT'))
  table(res_limma$change)
  res_limma = res_limma[res_limma$change %in% c('UP', 'DOWN'), ]
  resMeR_list[[c]] <- list(res_DESeq2, res_limma)
  countMatrix_list[[c]] <- countMatrix
  resPeak_list[[c]] <- merge_peak
  name = res@elementMetadata$description[2]
  name = paste0(cl, '_',str_split(name, "groupList ", simplify = T)[,2])
  name = paste0(gsub(' ','_',(name)))
  names(resMeR_list[[c]] ) <- paste0(name, '_', c('DESeq2', 'limma'))
}
RMdem_list[[i]] <- as.list(unlist(resMeR_list, recursive = FALSE))
name <- str_split(names(RMdem_list[[i]] ), "[.]", simplify = T)[,2]
names(RMdem_list[[i]] ) <- name
name
RMdep_list[[i]] <- resPeak_list 
name <- unique(gsub('(_limma|_DESeq2)', '', name))
names(RMdep_list[[i]] ) <- name
contrast_IP_list[[i]] <- lapply(countMatrix_list, extract_names)
condition_IP_list[[i]] <- condition_list
############
n = 16
############
#GSE78040# with miCLIP and RNA-seq, without case RIP 
############
print(RMPdeg_datasets[n,])
i = RMPdeg_datasets[n,'GSE']
cl = 'HEK293T'
RMp = c('FTO','ALKBH5','DCP2')
RM = 'm6Am'
condition_list = list(
  ALKBH5KD = c(rep('ALKBH5KD', 2), rep('WT', 2)),
  DCP2KD = c(rep('DCP2KD', 2), rep('WT', 2)),
  FTOKD = c(rep('FTOKD', 2), rep('WT', 2)))
# View(allDat[[i]])
files <- list.files(file.path(data_directory,i))
files
getwd()
setwd(data_directory)
setwd(i)
###################################
#      RNA/Ribo-seq analysis      #
###################################
# 1) RNA-seq
# Filter files that end with 'gene.counts.results.txt.gz'
count_files <- files[grepl('gene.counts.results.txt.gz', files)]
###################################
# solutions <- lapply(count_files, function(x){
#   fread(file.path(data_directory,i,'',x)) # read.table(paste0(data_directory,i,x),header = T)
#   })
# names(solutions) <- gsub('_gene.counts.results.txt.gz','',count_files)
# df_empty <- data.frame(gene_id = solutions[[1]]$gene_id)
# # Process count data and merge into df_empty
# countGTF <- lapply(names(solutions), function(y) {
#   # Extract required columns
#   df_new <- solutions[[y]][, c(1, 5)]
#   # Set column names same as file names
#   colnames(df_new) <- y
#   # Merge with df_empty by gene_id
#   df_empty <<- merge(df_empty, df_new, by = 'gene_id', all.x = TRUE)
#   return(df_empty)  
# })
###################################
# Read and process each count file
count <- lapply(count_files, function(x) {
  # Read data using fread for faster reading
  data <- fread(file.path(data_directory, i, x))
  # Extract required columns
  df_new <- data[, c(1, 5)]
  # Set column names same as file names
  colnames(df_new) <- c('geneSymbol',gsub('_gene.counts.results.txt.gz', '', x))
  return(df_new)
})
# Merge dataframes in countGTF list by 'gene_id'
df_in <- Reduce(function(x, y) merge(x, y, by = 'geneSymbol', all = TRUE), count)
colSums(df_in[,-1])
table(colSums(df_in[,-1])=="1e+06")
###################################
#   more than 1 condition here   #
###################################
resInP_list <- list()
##Gene count matrix
countGTF=as.matrix(ceiling(df_in[, -1]))
## Gene name
table(duplicated(df_in$geneSymbol))
dim(countGTF) #25343    12
## Sample name
colnames(countGTF)
## Sample information
sample_meta=data.frame(
  sample=colnames(countGTF),
  condition=c(rep("WT_ALKBH5",2),rep("ALKBH5KD",2), rep("WT_DCP2",2),rep("DCP2KD",2), rep("WT_FTO",2),rep("FTOKD",2)))
#####condition_ALKBH5KD_vs_WT#####
count <- as.matrix(countGTF[, !grepl('_d|_f', colnames(countGTF))])
table(is.na(count))
sample <- sample_meta[sample_meta$sample %in% colnames(count), ]
## SummarizedExperiment input
condition = c(rep('ALKBH5KD', 2), rep('WT_ALKBH5', 2))
dds<-DESeqDataSetFromMatrix(countData=count, colData=sample,design=~condition)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 10
table(keep)
dds <- dds[keep,]
## Note on factor levels
dds$condition <- factor(dds$condition, levels = c("ALKBH5KD","WT_ALKBH5"))
dds$condition <- relevel(dds$condition, ref = "WT_ALKBH5")
dds <- DESeq(dds)
res<-results(dds,name = "condition_ALKBH5KD_vs_WT_ALKBH5")
res <-  results(dds, contrast=c("condition",unique(condition))) # c("condition", rev(levels(genotype)))
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
table(res$change)
res$geneSymbol <-df_in$geneSymbol[keep]
resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
table(duplicated(resInP$geneSymbol))
resInP_list[['ALKBH5KD_vs_WT']] <- resInP
#####condition_DCP2KD_vs_WT#####
count <- as.matrix(countGTF[, !grepl('_a|_f', colnames(countGTF))])
table(is.na(count))
sample <- sample_meta[sample_meta$sample %in% colnames(count), ]
## SummarizedExperiment input
condition = c(rep('DCP2KD', 2), rep('WT_DCP2', 2))
dds<-DESeqDataSetFromMatrix(countData=count, colData=sample,design=~condition)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 10
table(keep)
dds <- dds[keep,]
## Note on factor levels
dds$condition <- factor(dds$condition, levels = c("DCP2KD","WT_DCP2"))
dds$condition <- relevel(dds$condition, ref = "WT_DCP2")
dds <- DESeq(dds)
res<-results(dds,name = "condition_DCP2KD_vs_WT_DCP2")
res <-  results(dds, contrast=c("condition",unique(condition))) # c("condition", rev(levels(genotype)))
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
table(res$change)
res$geneSymbol <-df_in$geneSymbol[keep]
resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
table(duplicated(resInP$geneSymbol))
resInP_list[['DCP2KD_vs_WT']] <- resInP

#####condition_FTOKD_vs_WT#####
count <- as.matrix(countGTF[, !grepl('_a|_d', colnames(countGTF))])
table(is.na(count))
sample <- sample_meta[sample_meta$sample %in% colnames(count), ]
## SummarizedExperiment input
condition = c(rep('FTOKD', 2), rep('WT_FTO', 2))
dds<-DESeqDataSetFromMatrix(countData=count, colData=sample,design=~condition)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 10
table(keep)
dds <- dds[keep,]
## Note on factor levels
dds$condition <- factor(dds$condition, levels = c("FTOKD","WT_FTO"))
dds$condition <- relevel(dds$condition, ref = "WT_FTO")
dds <- DESeq(dds)
res<-results(dds,name = "condition_FTOKD_vs_WT_FTO")
res <-  results(dds, contrast=c("condition",unique(condition))) # c("condition", rev(levels(genotype)))
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
table(res$change)
res$geneSymbol <-df_in$geneSymbol[keep]
resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
table(duplicated(resInP$geneSymbol))
resInP_list[['FTOKD_vs_WT']] <- resInP

RMdeg_list[[i]] <- resInP_list
name <- paste0(cl,'_', names(resInP_list))
names(RMdeg_list[[i]] ) <- name

# names(condition_list) <- gsub('_vs_WT', '', names(RMdeg_list[[i]] ))
contrast_input_list[[i]] <- NULL
condition_input_list[[i]] <- condition_list

# ###################################
# #         RIP-seq analysis        #
# ###################################
# # m6Am_regions.bed is a .bed format file and contains the coordinates of m6Am-containing regions in HEK293T cells. hek293t_rnaseq.bed is a .bed format file that contains the RNA-Seq for normalization of m6Am site relative stoichiometry. The *counts* files contain raw STAR aligner counts.
# peaks <- read.table(file.path(data_directory,i,'GSM2065370_m6am_regions.bed.gz'),header = F)
# peaks <- read.table(file.path(data_directory,i,'GSM2340158_hek293t_rnaseq.bed.gz'),header = F)

############
n = 17
############
#GSE112276# with RIP and Input
############
print(RMPdeg_datasets[n,])
i = RMPdeg_datasets[n,'GSE']
cl = c('HeLa', 'HepG2')
RMp = 'METTL1'
RM = 'm7G'
condition_list = list(Hela_METTL1KD = c(rep('METTL1KD', 2), rep('WT', 2)),
                      HepG2_METTL1KD = c(rep('METTL1KD', 2), rep('WT', 2)))
RMP_allDat[grep('GSE112276', RMP_allDat$study_alias ) 
           ,c("study_alias", "sample_alias", "sample_title")]
ENA <- as.data.frame(RMP_allDat[grepl('GSE112276', RMP_allDat$study_alias ) 
                                ,c("study_alias", "sample_alias", "sample_title")])
ENA[order(ENA$sample_alias),]
ENA <- ENA[order(ENA$sample_alias),]
ENA[grepl('GSE112276', ENA$study_alias ) & grepl('Control|METTL1', ENA$sample_title) &  grepl('dataset 4', ENA$sample_title)
    ,c("study_alias", "sample_alias", "sample_title")]
skip <- RMP_allDat[grepl('GSE112276', RMP_allDat$study_alias ) & grepl('Control|METTL1', RMP_allDat$sample_title) &  grepl('dataset 4', RMP_allDat$sample_title)
                   , "sample_alias"]
paste0(skip$sample_alias, collapse = "|")

skip <- RMP_allDat[grepl('GSE112276', RMP_allDat$study_alias ) & grepl('Control-IP|METTL1-IP', RMP_allDat$sample_title) &  grepl('dataset 4', RMP_allDat$sample_title)
                   , "sample_alias"]
paste0(skip$sample_alias, collapse = " ")

skip <- RMP_allDat[grepl('GSE112276', RMP_allDat$study_alias ) & grepl('Control-input|METTL1-input', RMP_allDat$sample_title) &  grepl('dataset 4', RMP_allDat$sample_title)
                   , "sample_alias"]
paste0(skip$sample_alias, collapse = " ")

# View(allDat[[i]])
y <- ENA[grepl('HeLa_siControl-IP',ENA$sample_title),]
y <- ENA[grepl('HeLa_siMETTL1-IP',ENA$sample_title),]
y <- ENA[grepl('HepG2_shControl-IP',ENA$sample_title),]
y <- ENA[grepl('HepG2_shMETTL1-IP',ENA$sample_title),]
y
x <- paste0('$outdir/GSE112276/', y$sample_alias, '_IP.bam')
paste(x, collapse = ' ')

files <- list.files(file.path(data_directory,i))
files
getwd()
setwd(data_directory)
setwd(i)
countfiles <- list.files('./metaPlotR/')
countfiles <- countfiles[grep('.txt', countfiles)]
countfiles <- countfiles[! grepl('summary', countfiles)] #mouse embryonic fibroblasts (MEFs)
countfiles
###################################
#  Cell number-normalized (CNN)   #
#      MeRIP INPUT analysis       #
###################################
resInP_list <- vector("list", length(condition_list))
names(resInP_list) <- names(condition_list)
countGTF_list <- vector("list", length(condition_list))
names(countGTF_list) <- names(condition_list)
for(c in names(condition_list)){
  countGTF <- read.table( paste0("metaPlotR/Input_counts_", c, ".txt"),header=T)
  dim(countGTF)
  colnames(countGTF)
  rownames(countGTF) <- countGTF$Geneid
  countGTF <- countGTF[, !colnames(countGTF) %in% c('Geneid','Chr','Start','End','Strand','Length')]
  dim(countGTF) # 28265     4
  condition <- condition_list[[c]]
  ## Sample information
  sample=data.frame(
    sample=colnames(countGTF),
    condition=condition)
  ## SummarizedExperiment input
  dds<-DESeqDataSetFromMatrix(countData=countGTF,colData=sample,design=~condition)
  ## Pre-filtering
  keep <- rowSums(counts(dds)) >= 10
  table(keep)
  dds <- dds[keep,]
  ## Note on factor levels
  dds$condition <- factor(dds$condition, levels = unique(condition))
  dds$condition <- relevel(dds$condition, ref = "WT")
  # Differential expression analysis-------------------------------------------
  dds <- DESeq(dds)
  res <-  results(dds, contrast=c("condition",unique(condition)))
  res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
  length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
  res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                                 ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
  table(res$change)
  table(is.na(res$change))
  # resOrdered <- as.data.frame(res[order(res$padj), ])
  resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
  resInP$geneSymbol <- rownames(resInP)
  table(resInP$change)
  resInP_list[[c]] <- resInP
  countGTF_list[[c]] <- countGTF
}
name = paste0(names(resInP_list), ' vs WT')
# name = paste0(cl, '_', name)
RMdeg_list[[i]] <- resInP_list 
names(RMdeg_list[[i]] ) <- paste0(gsub(' ','_',(name)))
contrast_input_list[[i]] <- lapply(countGTF_list, extract_names)
condition_input_list[[i]] <- condition_list

###################################
#  Cell number-normalized (CNN)   #
#   MeRIP RIP analysis - Part1    #
###################################
resMeR_list <- vector("list", length(condition_list))
names(resMeR_list) <- names(condition_list)
countMatrix_list <- vector("list", length(condition_list))
names(countMatrix_list) <- names(condition_list)
resPeak_list <- vector("list",length(condition_list))
names(resPeak_list) <- names(condition_list)

for(c in names(condition_list)){
  merge_peak <- read.delim(paste0("metaPlotR/", c, "_newPeakFile.txt"),header=F)
  countMatrix <- read.table(paste0("metaPlotR/m7G_counts_", c, ".txt"),header=T)
  table(countMatrix$Geneid == merge_peak[-1,1])
  table(duplicated(merge_peak[-1,1]))
  dim(countMatrix)
  colnames(countMatrix)
  rownames(countMatrix) <- countMatrix$Geneid
  countMatrix <- countMatrix[, !colnames(countMatrix) %in% c('Geneid','Chr','Start','End','Strand','Length')]
  # 1. 做好分组因子即可
  condition <- condition_list[[c]]
  groupList<-factor(condition)
  colData <- data.frame(row.names=colnames(countMatrix), groupList=groupList)
  # 2. Exploring Method 1: DESeq2-------------------------------------------
  # 构建 DESeqDataSet 对象
  dds <- DESeqDataSetFromMatrix(countData = countMatrix,
                                colData = colData,
                                design = ~ groupList)
  ## Pre-filtering
  keep <- rowSums(counts(dds)) >= 1
  table(keep)
  dds <- dds[keep,]
  # 接着计算每个样本的 library size 用于后续的标准化。library size 等于落在该样本所有peaks上的reads的总数，即 counts 矩阵中每列的加和。如果bFullLibrarySize设为TRUE，则会使用library的总reads数（根据原BAM/BED文件统计获得）。然后使用estimateDispersions估计统计分布，需要将参数fitType设为’local’
  dds<-estimateSizeFactors(dds)
  dds<-estimateDispersions(dds,fitType="local")
  dds <- nbinomWaldTest(dds)
  res <- results(dds)
  ## 提取你想要的差异分析结果，我们这里是treated组对untreated组进行比较
  res <-  results(dds, contrast=c("groupList", unique(condition)))
  res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
  length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
  res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                                 ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
  res_DESeq2 = as.data.frame(res)
  table(res_DESeq2$change)
  res_DESeq2 = res_DESeq2[res_DESeq2$change %in% c('UP', 'DOWN'), ]
  # 2. Exploring Method 2: limma-------------------------------------------
  library(edgeR)
  d0 <- DGEList(counts = countMatrix)
  dim(d0)
  # Filtering low count peaks
  cutoff <- 1           
  drop <- which(apply(cpm(d0), 1, max) < cutoff)
  print(length(drop))
  if (length(drop) >0) {
    d <- d0[-drop,] 
  } else{ d <- d0}
  dim(d) #number of peaks left: 8501      
  d <- calcNormFactors(d, method = c('TMM'))
  colnames(countMatrix)
  genotype <- groupList
  table(genotype)
  mm <- model.matrix(~0 + genotype)#
  v <- voom(d, mm, plot=F, normalize="quantile")#voom from limma:Transform count data to log2-counts per million (logtpm)
  fit <- lmFit(v, mm)
  head(coef(fit))
  contr <- makeContrasts(paste0(colnames(coef(fit)), collapse  = ' - '), levels = colnames(coef(fit))) #Change if other comparison
  tmp <- contrasts.fit(fit, contr)
  tmp <- eBayes(tmp)
  top.table <- topTable(tmp, sort.by = "P", n = Inf)
  head(top.table, 20)
  length(which(top.table$adj.P.Val < 0.05 & abs(top.table$logFC) > log2(1.5))) #how many are significant?
  res_limma <- as.data.frame(top.table)
  res_limma$change <- with(res_limma, ifelse(res_limma$adj.P.Val < 0.05 & abs(res_limma$logFC) > log2(1.5),
                                             ifelse(res_limma$logFC > log2(1.5),'UP','DOWN'),'NOT'))
  table(res_limma$change)
  res_limma = res_limma[res_limma$change %in% c('UP', 'DOWN'), ]
  resMeR_list[[c]] <- list(res_DESeq2, res_limma)
  countMatrix_list[[c]] <- countMatrix
  resPeak_list[[c]] <- merge_peak
  name = res@elementMetadata$description[2]
  # name = paste0(cl, '_',str_split(name, "groupList ", simplify = T)[,2])
  name = str_split(name, "groupList ", simplify = T)[,2]
  name = paste0(gsub(' ','_',(name)))
  names(resMeR_list[[c]] ) <- paste0(name, '_', c('DESeq2', 'limma'))
}
RMdem_list[[i]] <- as.list(unlist(resMeR_list, recursive = FALSE))
name <- gsub("_METTL1KD.", '_', names(RMdem_list[[i]] ))#str_split(names(RMdem_list[[i]] ), "[.]", simplify = T)[,2]
names(RMdem_list[[i]] ) <- name
RMdep_list[[i]] <- resPeak_list 
name <- unique(gsub('(_limma|_DESeq2)', '', name))
names(RMdep_list[[i]] ) <- name
contrast_IP_list[[i]] <- lapply(countMatrix_list, extract_names)
condition_IP_list[[i]] <- condition_list

############
n = 18
############
#GSE60047# with RIP and Input
############
print(RMPdeg_datasets[n,])
i = RMPdeg_datasets[n,'GSE']
cl = 'HEK293'
RMp = 'DKC1'
RM = 'Psi'
condition = c(rep('DKC1KD', 2), rep('WT', 2))
RMP_allDat[grep('GSE60047', RMP_allDat$study_alias ) 
           ,c("study_alias", "sample_alias", "sample_title")]
ENA <- as.data.frame(RMP_allDat[grepl('GSE60047', RMP_allDat$study_alias ) 
                                ,c("study_alias", "sample_alias", "sample_title")])
ENA[order(ENA$sample_alias),]
ENA <- ENA[order(ENA$sample_alias),]
ENA[grepl('GSE60047', ENA$study_alias ) & grepl('Lib72PUdkc1|Lib72PUcontrol', ENA$sample_title)
    ,c("study_alias", "sample_alias", "sample_title")]

RMP_allDat[grepl('GSE60047', RMP_allDat$study_alias ) & grepl('Lib72PUdkc1|Lib72PUcontrol', RMP_allDat$sample_title)
           ,c("study_alias", "sample_alias", "sample_title")]
skip <- RMP_allDat[grepl('GSE60047', RMP_allDat$study_alias ) & grepl('Lib72PUdkc1|Lib72PUcontrol', RMP_allDat$sample_title)
                   , "sample_alias"]
paste0(skip$sample_alias, collapse = "|")

skip <- RMP_allDat[grepl('GSE60047', RMP_allDat$study_alias ) & grepl('Lib72PUdkc1|Lib72PUcontrol', RMP_allDat$sample_title) & 
                     grepl('_input', RMP_allDat$sample_title)
                   , "sample_alias"]
paste0(skip$sample_alias, collapse = " ")
skip <- RMP_allDat[grepl('GSE60047', RMP_allDat$study_alias ) & grepl('Lib72PUdkc1|Lib72PUcontrol', RMP_allDat$sample_title) & 
                     grepl('CMC_treatment', RMP_allDat$sample_title)
                   , "sample_alias"]
paste0(skip$sample_alias, collapse = " ")

# View(allDat[[i]])
files <- list.files(file.path(data_directory,i))
files
getwd()
setwd(data_directory)
setwd(i)
c <- unique(condition[!condition %in% 'WT'])
###################################
#  Cell number-normalized (CNN)   #
#      MeRIP INPUT analysis       #
###################################
countGTF <- read.table( paste0("metaPlotR/Input_counts_", c, ".txt"),header=T)
rownames(countGTF) <- countGTF$Geneid
countGTF <- countGTF[, !colnames(countGTF) %in% c('Geneid','Chr','Start','End','Strand','Length')]
dim(countGTF) #28265     4
## Sample name
## Sample information
sample=data.frame(
  sample=colnames(countGTF),
  condition=condition)
## SummarizedExperiment input
dds<-DESeqDataSetFromMatrix(countData=countGTF,colData=sample,design=~condition)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 10
table(keep)
dds <- dds[keep,]
## Note on factor levels
dds$condition <- factor(dds$condition, levels = unique(condition))
dds$condition <- relevel(dds$condition, ref = "WT")
# Differential expression analysis-------------------------------------------
dds <- DESeq(dds)

res <-  results(dds, contrast=c("condition",unique(condition)))
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
table(res$change)
table(is.na(res$change))
# resOrdered <- as.data.frame(res[order(res$padj), ])
resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
resInP$geneSymbol <- rownames(resInP)
name = res@elementMetadata$description[2]
name = paste0(cl, '_',str_split(name, "condition ", simplify = T)[,2])
RMdeg_list[[i]] <- list( resInP)
names(RMdeg_list[[i]] ) <- paste0(gsub(' ','_',(name)))
c <- unique(condition[!condition %in% 'WT'])
countGTF_list <- list(colnames(countGTF)); names(countGTF_list) <- c
condition_list <- list(condition); names(condition_list) <- c
contrast_input_list[[i]] <- countGTF_list
condition_input_list[[i]] <- condition_list

###################################
#  Cell number-normalized (CNN)   #
#   MeRIP RIP analysis - Part1    #
###################################
## Sample information
merge_peak <- read.delim(paste0("metaPlotR/", c, "_newPeakFile.txt"),header=F)
countMatrix <- read.table(paste0("metaPlotR/Psi_counts_", c, ".txt"),header=T)
table(countMatrix$Geneid == merge_peak[-1,1])
table(duplicated(merge_peak$V1[-1]))
dim(countMatrix)
colnames(countMatrix)
rownames(countMatrix) <- countMatrix$Geneid
countMatrix <- countMatrix[, !colnames(countMatrix) %in% c('Geneid','Chr','Start','End','Strand','Length')]
# 1. 做好分组因子即可
groupList<-factor(condition)
colData <- data.frame(row.names=colnames(countMatrix), groupList=groupList)
# 2. Exploring Method 1: DESeq2-------------------------------------------
# 构建 DESeqDataSet 对象
dds <- DESeqDataSetFromMatrix(countData = countMatrix,
                              colData = colData,
                              design = ~ groupList)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 1
table(keep)
dds <- dds[keep,]
# 接着计算每个样本的 library size 用于后续的标准化。library size 等于落在该样本所有peaks上的reads的总数，即 counts 矩阵中每列的加和。如果bFullLibrarySize设为TRUE，则会使用library的总reads数（根据原BAM/BED文件统计获得）。然后使用estimateDispersions估计统计分布，需要将参数fitType设为’local’
dds<-estimateSizeFactors(dds)
dds<-estimateDispersions(dds,fitType="local")
dds <- nbinomWaldTest(dds)
res <- results(dds)
## 提取你想要的差异分析结果，我们这里是treated组对untreated组进行比较
res <-  results(dds, contrast=c("groupList", unique(condition)))
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
res_DESeq2 = as.data.frame(res)
table(res_DESeq2$change)#Not only
res_DESeq2 = res_DESeq2[res_DESeq2$change %in% c('UP', 'DOWN'), ]
# 2. Exploring Method 2: limma-------------------------------------------
library(edgeR)
d0 <- DGEList(counts = countMatrix)
dim(d0)
# Filtering low count peaks
cutoff <- 1           
drop <- which(apply(cpm(d0), 1, max) < cutoff)
print(length(drop))
if (length(drop) >0) {
  d <- d0[-drop,] 
} else{ d <- d0}
dim(d) #number of peaks left: 8501      
d <- calcNormFactors(d, method = c('TMM'))
colnames(countMatrix)
genotype <- groupList
table(genotype)
mm <- model.matrix(~0 + genotype)#
v <- voom(d, mm, plot=F, normalize="quantile")#voom from limma:Transform count data to log2-counts per million (logtpm)
fit <- lmFit(v, mm)
head(coef(fit))
contr <- makeContrasts(paste0(colnames(coef(fit)), collapse  = ' - '), levels = colnames(coef(fit))) #Change if other comparison
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)
length(which(top.table$adj.P.Val < 0.05 & abs(top.table$logFC) > log2(1.5))) #how many are significant?
res_limma <- as.data.frame(top.table)
res_limma$change <- with(res_limma, ifelse(res_limma$adj.P.Val < 0.05 & abs(res_limma$logFC) > log2(1.5),
                                           ifelse(res_limma$logFC > log2(1.5),'UP','DOWN'),'NOT'))
table(res_limma$change)#NOT only
res_limma = res_limma[res_limma$change %in% c('UP', 'DOWN'), ]

RMdem_list[[i]] <- list(res_DESeq2, res_limma)
name = res@elementMetadata$description[2]
name = paste0(cl, '_',str_split(name, "groupList ", simplify = T)[,2])
name= paste0(gsub(' ','_',(name)))
names(RMdem_list[[i]] ) <- paste0(name, '_', c('DESeq2', 'limma'))
RMdep_list[[i]] <- list( merge_peak )
names(RMdep_list[[i]] ) <- name

countMatrix_list <- list(colnames(countMatrix)); names(countMatrix_list) <- c
condition_list <- list(condition); names(condition_list) <- c
contrast_IP_list[[i]] <- countMatrix_list
condition_IP_list[[i]] <- condition_list

############
n = 19
############
#GSE102493# with RIP and Input
############
print(RMPdeg_datasets[n,])#a characterization of the full m6A methyltransferase complex in HeLa cells identifying METTL3/METTL14/WTAP/VIRMA/HAKAI/ZC3H13 as the key components, and we show that VIRMA mediates preferential mRNA methylation in 3`UTR and near stop codon. Biochemical studies reveal that VIRMA recruits the catalytic core components METTL3/METTL14/WTAP to guide region-selective methylations. Around 60% of VIRMA mRNA immunoprecipitation targets manifest strong m6A enrichment in 3`UTR. Depletions of VIRMA and METTL3 induce 3`UTR lengthening of several hundred mRNAs with over 50% targets in common. 
i = RMPdeg_datasets[n,'GSE']
cl = 'HeLa'
RMp = c('ZC3H13','HAKAI','VIRMA', 'METTL3', 'METTL14')
RM = 'm6A'
condition = c(rep('CPSF5', 3), rep('WT', 2))
RMP_allDat[grep('GSE102493', RMP_allDat$study_alias ) 
           ,c("study_alias", "sample_alias", "sample_title")]

RMP_allDat[grepl('GSE102493', RMP_allDat$study_alias )
           ,c("study_alias", "sample_alias", "sample_title")]

skip <- RMP_allDat[grepl('GSE102493', RMP_allDat$study_alias ) &  grepl('_Input', RMP_allDat$sample_title)
                   , "sample_alias"]
paste0(skip$sample_alias, collapse = " ")
skip <- RMP_allDat[grepl('GSE102493', RMP_allDat$study_alias ) &  grepl('_IP', RMP_allDat$sample_title)
                   , "sample_alias"]
paste0(skip$sample_alias, collapse = " ")
ENA <- as.data.frame(RMP_allDat[grepl('GSE102493', RMP_allDat$study_alias ) 
                                ,c("study_alias", "sample_alias", "sample_title")])
ENA <- ENA[order(ENA$sample_alias),]
condition_list <- list(VIRMAPD = c(rep('VIRMAPD', 2), rep('WT', 4)),
                       METTL3MU = c(rep('METTL3MU', 2), rep('WT', 4)),
                       METTL14MU = c(rep('METTL14MU', 2), rep('WT', 4)),
                       ZC3H13KD = c(rep('ZC3H13KD', 2), rep('WT', 8)),
                       HAKAIKD = c(rep('HAKAIKD', 2), rep('WT', 8)),
                       VIRMAKD = c(rep('VIRMAKD', 2), rep('WT', 8)),
                       VIRMA_1MU = c(rep('VIRMAMU', 2), rep('WT', 2)), #single end
                       VIRMA_2MU = c(rep('VIRMAMU', 2), rep('WT', 4)) #paied end
)

# View(allDat[[i]])
files <- list.files(file.path(data_directory,i))
files
getwd()
setwd(data_directory)
setwd(i)
ENA <- ENA[order(ENA$sample_alias),]
# countfiles <- list.files('./metaPlotR/')
metaPlotR <- file.path(raw_directory, "Out/macs2", i)
countfiles <- list.files(metaPlotR)
countfiles <- countfiles[grep('counts', countfiles)]
countfiles <- countfiles[! grepl('summary', countfiles)] 
countfiles
###################################
#  Cell number-normalized (CNN)   #
#      MeRIP INPUT analysis       #
###################################
resInP_list <- vector("list", length(condition_list))
names(resInP_list) <- names(condition_list)
countGTF_list <- vector("list", length(condition_list))
names(countGTF_list) <- names(condition_list)
for(c in names(condition_list)){
  countGTF <- read.table( paste0(metaPlotR,"/Input_counts_", c, ".txt"),header=T)
  dim(countGTF)
  colnames(countGTF)
  rownames(countGTF) <- countGTF$Geneid
  countGTF <- countGTF[, !colnames(countGTF) %in% c('Geneid','Chr','Start','End','Strand','Length')]
  dim(countGTF) # 28265     4
  condition <- condition_list[[c]]
  ## Sample information
  sample=data.frame(
    sample=colnames(countGTF),
    condition=condition)
  ## SummarizedExperiment input
  dds<-DESeqDataSetFromMatrix(countData=countGTF,colData=sample,design=~condition)
  ## Pre-filtering
  keep <- rowSums(counts(dds)) >= 10
  table(keep)
  dds <- dds[keep,]
  ## Note on factor levels
  dds$condition <- factor(dds$condition, levels = unique(condition))
  dds$condition <- relevel(dds$condition, ref = "WT")
  # Differential expression analysis-------------------------------------------
  dds <- DESeq(dds)
  res <-  results(dds, contrast=c("condition",unique(condition)))
  res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
  length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
  res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                                 ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
  table(res$change)
  table(is.na(res$change))
  # resOrdered <- as.data.frame(res[order(res$padj), ])
  resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
  resInP$geneSymbol <- rownames(resInP)
  table(resInP$change)
  resInP_list[[c]] <- resInP
  countGTF_list[[c]] <- countGTF
  
}
name = paste0(names(resInP_list), ' vs WT')
RMdeg_list[[i]] <- resInP_list 
name <- paste0(gsub(' ','_',(name)))
names(RMdeg_list[[i]] ) <-  paste0(cl, '_', name)

contrast_input_list[[i]] <- lapply(countGTF_list, extract_names)
condition_input_list[[i]] <- condition_list

###################################
#  Cell number-normalized (CNN)   #
#   MeRIP RIP analysis - Part1    #
###################################
resMeR_list <- vector("list", length(condition_list))
names(resMeR_list) <- names(condition_list)
countMatrix_list <- vector("list", length(condition_list))
names(countMatrix_list) <- names(condition_list)
resPeak_list <- vector("list",length(condition_list))
names(resPeak_list) <- names(condition_list)
for(c in names(condition_list)){
  merge_peak <- read.delim(paste0(metaPlotR, '/', c, "_newPeakFile.txt"),header=F)
  countMatrix <- read.table(paste0(metaPlotR, "/m6A_counts_", c, ".txt"),header=T)
  table(countMatrix$Geneid == merge_peak[-1,1])
  table(duplicated(merge_peak[-1,1]))
  dim(countMatrix)
  colnames(countMatrix)
  rownames(countMatrix) <- countMatrix$Geneid
  countMatrix <- countMatrix[, !colnames(countMatrix) %in% c('Geneid','Chr','Start','End','Strand','Length')]
  # 1. 做好分组因子即可
  condition <- condition_list[[c]]
  groupList<-factor(condition)
  colData <- data.frame(row.names=colnames(countMatrix), groupList=groupList)
  # 2. Exploring Method 1: DESeq2-------------------------------------------
  # 构建 DESeqDataSet 对象
  dds <- DESeqDataSetFromMatrix(countData = countMatrix,
                                colData = colData,
                                design = ~ groupList)
  ## Pre-filtering
  keep <- rowSums(counts(dds)) >= 1
  table(keep)
  dds <- dds[keep,]
  # 接着计算每个样本的 library size 用于后续的标准化。library size 等于落在该样本所有peaks上的reads的总数，即 counts 矩阵中每列的加和。如果bFullLibrarySize设为TRUE，则会使用library的总reads数（根据原BAM/BED文件统计获得）。然后使用estimateDispersions估计统计分布，需要将参数fitType设为’local’
  dds<-estimateSizeFactors(dds)
  dds<-estimateDispersions(dds,fitType="local")
  dds <- nbinomWaldTest(dds)
  res <- results(dds)
  ## 提取你想要的差异分析结果，我们这里是treated组对untreated组进行比较
  res <-  results(dds, contrast=c("groupList", unique(condition)))
  res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
  length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
  res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                                 ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
  res_DESeq2 = as.data.frame(res)
  table(res_DESeq2$change)
  res_DESeq2 = res_DESeq2[res_DESeq2$change %in% c('UP', 'DOWN'), ]
  # 2. Exploring Method 2: limma-------------------------------------------
  library(edgeR)
  d0 <- DGEList(counts = countMatrix)
  dim(d0)
  # Filtering low count peaks
  cutoff <- 1           
  drop <- which(apply(cpm(d0), 1, max) < cutoff)
  print(length(drop))
  if (length(drop) >0) {
    d <- d0[-drop,] 
  } else{ d <- d0}
  dim(d) #number of peaks left: 8501      
  d <- calcNormFactors(d, method = c('TMM'))
  colnames(countMatrix)
  genotype <- groupList
  table(genotype)
  mm <- model.matrix(~0 + genotype)#
  v <- voom(d, mm, plot=F, normalize="quantile")#voom from limma:Transform count data to log2-counts per million (logtpm)
  fit <- lmFit(v, mm)
  head(coef(fit))
  contr <- makeContrasts(paste0(colnames(coef(fit)), collapse  = ' - '), levels = colnames(coef(fit))) #Change if other comparison
  tmp <- contrasts.fit(fit, contr)
  tmp <- eBayes(tmp)
  top.table <- topTable(tmp, sort.by = "P", n = Inf)
  head(top.table, 20)
  length(which(top.table$adj.P.Val < 0.05 & abs(top.table$logFC) > log2(1.5))) #how many are significant?
  res_limma <- as.data.frame(top.table)
  res_limma$change <- with(res_limma, ifelse(res_limma$adj.P.Val < 0.05 & abs(res_limma$logFC) > log2(1.5),
                                             ifelse(res_limma$logFC > log2(1.5),'UP','DOWN'),'NOT'))
  table(res_limma$change)
  res_limma = res_limma[res_limma$change %in% c('UP', 'DOWN'), ]
  resMeR_list[[c]] <- list(res_DESeq2, res_limma)
  countMatrix_list[[c]] <- countMatrix
  resPeak_list[[c]] <- merge_peak
  name = paste0(c, ' vs WT')
  name = paste0(gsub(' ','_', name ))
  names(resMeR_list[[c]] ) <- paste0(name, '_', c('DESeq2', 'limma'))
}
names(resMeR_list) <- NULL
RMdem_list[[i]] <- as.list(unlist(resMeR_list, recursive = FALSE))
names(RMdem_list[[i]] )
names(RMdem_list[[i]]) <- paste0(cl, '_', names(RMdem_list[[i]] ))
RMdep_list[[i]] <- resPeak_list 
name <- unique(gsub('(_limma|_DESeq2)', '', names(RMdem_list[[i]] )))
names(RMdep_list[[i]] ) <- name# paste0(cl, '_', name)

contrast_IP_list[[i]] <- lapply(countMatrix_list, extract_names)
condition_IP_list[[i]] <- condition_list


############
n = 20
############
#GSE90914# with RIP and Input 
############
print(RMPdeg_datasets[n,])# hg38
i = RMPdeg_datasets[n,'GSE']
cl = 'HEK293A-TOA'
RMp = c('METTL16')
RM = 'm6A'
condition = c(rep('METTL16KD', 3), rep('WT', 3))
RMP_allDat[grep('GSE90914', RMP_allDat$study_alias ) 
           ,c("study_alias", "sample_alias", "sample_title")]

skip <- RMP_allDat[grepl('GSE90914', RMP_allDat$study_alias ) &  grepl('input', RMP_allDat$sample_title)
                   , "sample_alias"]
paste0(skip$sample_alias, collapse = " ")
skip <- RMP_allDat[grepl('GSE90914', RMP_allDat$study_alias ) &  grepl('IP', RMP_allDat$sample_title)
                   , "sample_alias"]
paste0(skip$sample_alias, collapse = " ")
ENA <- as.data.frame(RMP_allDat[grepl('GSE90914', RMP_allDat$study_alias ) 
                                ,c("study_alias", "sample_alias", "sample_title")])
ENA[order(ENA$sample_alias),]
ENA <- ENA[order(ENA$sample_alias),]
# y <- ENA[grepl('siMETT16 m6A',ENA$sample_title),]
# y <- ENA[grepl('siCtrl m6A',ENA$sample_title),]
y <- ENA[grepl('siMETTL16 input_',ENA$sample_title),]
y <- ENA[grepl('siCtrl input',ENA$sample_title),]
y
# x <- paste0('$outdir/GSE90914/', y$sample_alias, '_IP.bam')
x <- paste0('$outdir/GSE90914/', y$sample_alias, '_input.bam')
paste(x, collapse = ' ')
# View(allDat[[i]])
files <- list.files(file.path(data_directory,i))
files
getwd()
setwd(data_directory)
setwd(i)
c <- unique(condition[!condition %in% 'WT'])
###################################
#  Cell number-normalized (CNN)   #
#      MeRIP INPUT analysis       #
###################################
countGTF <- read.table( paste0("metaPlotR/Input_counts_", c, ".txt"),header=T)
rownames(countGTF) <- countGTF$Geneid
countGTF <- countGTF[, !colnames(countGTF) %in% c('Geneid','Chr','Start','End','Strand','Length')]
dim(countGTF) #28265     4
## Sample name
## Sample information
sample=data.frame(
  sample=colnames(countGTF),
  condition=condition)
## SummarizedExperiment input
dds<-DESeqDataSetFromMatrix(countData=countGTF,colData=sample,design=~condition)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 10
table(keep)
dds <- dds[keep,]
## Note on factor levels
dds$condition <- factor(dds$condition, levels = unique(condition))
dds$condition <- relevel(dds$condition, ref = "WT")
# Differential expression analysis-------------------------------------------
dds <- DESeq(dds)

res <-  results(dds, contrast=c("condition",unique(condition)))
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
table(res$change)
table(is.na(res$change))
# resOrdered <- as.data.frame(res[order(res$padj), ])
resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
resInP$geneSymbol <- rownames(resInP)
name = res@elementMetadata$description[2]
name = paste0(cl, '_',str_split(name, "condition ", simplify = T)[,2])
RMdeg_list[[i]] <- list( resInP)
names(RMdeg_list[[i]] ) <- paste0(gsub(' ','_',(name)))
c <- unique(condition[!condition %in% 'WT'])
countGTF_list <- list(colnames(countGTF)); names(countGTF_list) <- c
condition_list <- list(condition); names(condition_list) <- c
contrast_input_list[[i]] <- countGTF_list
condition_input_list[[i]] <- condition_list

###################################
#  Cell number-normalized (CNN)   #
#   MeRIP RIP analysis - Part1    #
###################################
## Sample information
merge_peak <- read.delim(paste0("metaPlotR/", c, "_newPeakFile.txt"),header=F)
countMatrix <- read.table(paste0("metaPlotR/m6A_counts_", c, ".txt"),header=T)
table(countMatrix$Geneid == merge_peak[-1,1])
table(duplicated(merge_peak$V1[-1]))
dim(countMatrix)
colnames(countMatrix)
rownames(countMatrix) <- countMatrix$Geneid
countMatrix <- countMatrix[, !colnames(countMatrix) %in% c('Geneid','Chr','Start','End','Strand','Length')]
# 1. 做好分组因子即可
groupList<-factor(condition)
colData <- data.frame(row.names=colnames(countMatrix), groupList=groupList)
# 2. Exploring Method 1: DESeq2-------------------------------------------
# 构建 DESeqDataSet 对象
dds <- DESeqDataSetFromMatrix(countData = countMatrix,
                              colData = colData,
                              design = ~ groupList)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 1
table(keep)
dds <- dds[keep,]
# 接着计算每个样本的 library size 用于后续的标准化。library size 等于落在该样本所有peaks上的reads的总数，即 counts 矩阵中每列的加和。如果bFullLibrarySize设为TRUE，则会使用library的总reads数（根据原BAM/BED文件统计获得）。然后使用estimateDispersions估计统计分布，需要将参数fitType设为’local’
dds<-estimateSizeFactors(dds)
dds<-estimateDispersions(dds,fitType="local")
dds <- nbinomWaldTest(dds)
res <- results(dds)
## 提取你想要的差异分析结果，我们这里是treated组对untreated组进行比较
res <-  results(dds, contrast=c("groupList", unique(condition)))
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
res_DESeq2 = as.data.frame(res)
table(res_DESeq2$change)#Not only
res_DESeq2 = res_DESeq2[res_DESeq2$change %in% c('UP', 'DOWN'), ]
# 2. Exploring Method 2: limma-------------------------------------------
library(edgeR)
d0 <- DGEList(counts = countMatrix)
dim(d0)
# Filtering low count peaks
cutoff <- 1           
drop <- which(apply(cpm(d0), 1, max) < cutoff)
print(length(drop))
if (length(drop) >0) {
  d <- d0[-drop,] 
} else{ d <- d0}
dim(d) #number of peaks left: 8501      
d <- calcNormFactors(d, method = c('TMM'))
colnames(countMatrix)
genotype <- groupList
table(genotype)
mm <- model.matrix(~0 + genotype)#
v <- voom(d, mm, plot=F, normalize="quantile")#voom from limma:Transform count data to log2-counts per million (logtpm)
fit <- lmFit(v, mm)
head(coef(fit))
contr <- makeContrasts(paste0(colnames(coef(fit)), collapse  = ' - '), levels = colnames(coef(fit))) #Change if other comparison
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)
length(which(top.table$adj.P.Val < 0.05 & abs(top.table$logFC) > log2(1.5))) #how many are significant?
res_limma <- as.data.frame(top.table)
res_limma$change <- with(res_limma, ifelse(res_limma$adj.P.Val < 0.05 & abs(res_limma$logFC) > log2(1.5),
                                           ifelse(res_limma$logFC > log2(1.5),'UP','DOWN'),'NOT'))
table(res_limma$change)#NOT only
res_limma = res_limma[res_limma$change %in% c('UP', 'DOWN'), ]

RMdem_list[[i]] <- list(res_DESeq2, res_limma)
name = res@elementMetadata$description[2]
name = paste0(cl, '_',str_split(name, "groupList ", simplify = T)[,2])
name= paste0(gsub(' ','_',(name)))
names(RMdem_list[[i]] ) <- paste0(name, '_', c('DESeq2', 'limma'))
RMdep_list[[i]] <- list( merge_peak )
names(RMdep_list[[i]] ) <- name
countMatrix_list <- list(colnames(countMatrix)); names(countMatrix_list) <- c
condition_list <- list(condition); names(condition_list) <- c
contrast_IP_list[[i]] <- countMatrix_list
condition_IP_list[[i]] <- condition_list
############
n = 21
############
#GSE93911# with RIP and Input 
############
print(RMPdeg_datasets[n,]) # hg38
i = RMPdeg_datasets[n,'GSE']
cl = 'HEC-1-A'
RMp = c('METTL14')
RM = 'm6A'
condition_list = list(METTL14MU = c(rep('METTL14MU', 1), rep('WT', 1)),
                      METTL3KD = c(rep('METTL3KD', 1), rep('WT', 1)))

RMP_allDat[grep('GSE93911', RMP_allDat$study_alias ) 
           ,c("study_alias", "sample_alias", "sample_title")]
RMP_allDat[grepl('GSE93911', RMP_allDat$study_alias ) & grepl('HEC-1-A', RMP_allDat$sample_title)
           ,c("study_alias", "sample_alias", "sample_title")]
skip <- RMP_allDat[grepl('GSE93911', RMP_allDat$study_alias ) & grepl('HEC-1-A', RMP_allDat$sample_title)
                   , "sample_alias"]
paste0(skip$sample_alias, collapse = "|")

skip <- RMP_allDat[grepl('GSE93911', RMP_allDat$study_alias ) & grepl('HEC-1-A', RMP_allDat$sample_title)
                   &  grepl('input', RMP_allDat$sample_title)
                   , "sample_alias"]
paste0(skip$sample_alias, collapse = " ")
skip <- RMP_allDat[grepl('GSE93911', RMP_allDat$study_alias ) & grepl('HEC-1-A', RMP_allDat$sample_title)
                   &  grepl('IP', RMP_allDat$sample_title)
                   , "sample_alias"]
paste0(skip$sample_alias, collapse = " ")

ENA <- as.data.frame(RMP_allDat[grepl('GSE93911', RMP_allDat$study_alias ) & grepl('HEC-1-A', RMP_allDat$sample_title)
                                ,c("study_alias", "sample_alias", "sample_title")])
ENA[order(ENA$sample_alias),]
ENA <- ENA[order(ENA$sample_alias),]
# View(allDat[[i]])
files <- list.files(file.path(data_directory,i))
files
getwd()
setwd(data_directory)
setwd(i)
countfiles <- list.files('./metaPlotR/')
countfiles <- countfiles[grep('.txt', countfiles)]
countfiles <- countfiles[! grepl('summary', countfiles)] #mouse embryonic fibroblasts (MEFs)
countfiles
###################################
#  Cell number-normalized (CNN)   #
#      MeRIP INPUT analysis       #
###################################
resInP_list <- vector("list", length(condition_list))
names(resInP_list) <- names(condition_list)
countGTF_list <- vector("list", length(condition_list))
names(countGTF_list) <- names(condition_list)
for(c in names(condition_list)){
  countGTF <- read.table( paste0("metaPlotR/Input_counts_", c, ".txt"),header=T)
  dim(countGTF)
  colnames(countGTF)
  rownames(countGTF) <- countGTF$Geneid
  countGTF <- countGTF[, !colnames(countGTF) %in% c('Geneid','Chr','Start','End','Strand','Length')]
  dim(countGTF) # 28265     4
  condition <- condition_list[[c]]
  # 2. Differential gene analysis with edgeR
  # Without replicates, you cannot estimate which genes are differentially expressed using EdgeR or DESeq2. You can only calculate fold changes based on normalized read counts (preferentially CPM normalized by TMM method included in any EdgeR analysis) and apply a stringent fold-change cut-off to determine which genes are more or less expressed depending on the condition.
  y <- DGEList(countGTF)
  # Calculate fold change
  y <- calcNormFactors(y, method = c('TMM'))
  y$samples$lib.size <- colSums(y$counts)
  logcpm <- cpm(y, prior.count=2, log=TRUE)
  logcpm[1:5,]
  cpm <- cpm(y, prior.count=2, log=FALSE)
  cpm[1:5, ]
  fold_change <- cpm[,1]/cpm[,2]
  res <- countGTF
  res$logFC <- log2(fold_change)
  res$FoldChange <- fold_change
  res$change <- with(res, ifelse(abs(res$logFC) > log2(1.5),
                                 ifelse(res$logFC > log2(1.5),'UP','DOWN'),'NOT'))
  table(res$change)
  resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
  resInP$geneSymbol <- rownames(resInP)
  table(resInP$change)
  resInP_list[[c]] <- resInP
  countGTF_list[[c]] <- countGTF
}
name = paste0(names(resInP_list), ' vs WT')
name = paste0(cl, '_', name)
RMdeg_list[[i]] <- resInP_list 
names(RMdeg_list[[i]] ) <- paste0(gsub(' ','_',(name)))
contrast_input_list[[i]] <- lapply(countGTF_list, extract_names)
condition_input_list[[i]] <- condition_list

###################################
#  Cell number-normalized (CNN)   #
#   MeRIP RIP analysis - Part1    #
###################################
resMeR_list <- vector("list", length(condition_list))
names(resMeR_list) <- names(condition_list)
countMatrix_list <- vector("list", length(condition_list))
names(countMatrix_list) <- names(condition_list)
resPeak_list <- vector("list",length(condition_list))
names(resPeak_list) <- names(condition_list)
for(c in names(condition_list)){
  merge_peak <- read.delim(paste0("metaPlotR/", c, "_newPeakFile.txt"),header=F)
  countMatrix <- read.table(paste0("metaPlotR/m6A_counts_", c, ".txt"),header=T)
  table(countMatrix$Geneid == merge_peak[-1,1])
  table(duplicated(merge_peak[-1,1]))
  dim(countMatrix)
  colnames(countMatrix)
  rownames(countMatrix) <- countMatrix$Geneid
  countMatrix <- countMatrix[, !colnames(countMatrix) %in% c('Geneid','Chr','Start','End','Strand','Length')]
  # 2. Differential gene analysis with edgeR
  # Without replicates, you cannot estimate which genes are differentially expressed using EdgeR or DESeq2. You can only calculate fold changes based on normalized read counts (preferentially CPM normalized by TMM method included in any EdgeR analysis) and apply a stringent fold-change cut-off to determine which genes are more or less expressed depending on the condition.
  y <- DGEList(countMatrix)
  # Calculate fold change
  y <- calcNormFactors(y, method = c('TMM'))
  y$samples$lib.size <- colSums(y$counts)
  logcpm <- cpm(y, prior.count=2, log=TRUE)
  logcpm[1:5,]
  cpm <- cpm(y, prior.count=2, log=FALSE)
  cpm[1:5, ]
  # fold_change <- logcpm[,'Ratio1']/logcpm[,'Ratio2']
  fold_change <- cpm[,1]/cpm[,2]
  resMeR <- countMatrix
  resMeR$logFC <- log2(fold_change)
  resMeR$FoldChange <- fold_change
  resMeR$change <- with(resMeR, ifelse(abs(resMeR$logFC) > log2(1.5),
                                       ifelse(resMeR$logFC > log2(1.5),'UP','DOWN'),'NOT'))
  table(resMeR$change)
  resMeR = as.data.frame(resMeR[resMeR$change %in% c('UP', 'DOWN'), ])
  resMeR_list[[c]] <- resMeR
  countMatrix_list[[c]] <- countMatrix
  resPeak_list[[c]] <- merge_peak
}
RMdem_list[[i]] <- resMeR_list
name = paste0(names(resMeR_list), ' vs WT')#res@elementMetadata$description[2]
name = paste0(cl, '_', name)
name = paste0(gsub(' ','_', name ))
names(RMdem_list[[i]] ) <- paste0(name, '_', c('limma'))
RMdep_list[[i]] <- resPeak_list 
name <- unique(gsub('(_limma|_DESeq2)', '', name))
names(RMdep_list[[i]] ) <- name
contrast_IP_list[[i]] <- lapply(countMatrix_list, extract_names)
condition_IP_list[[i]] <- condition_list

############
n = 22
############
#GSE52600# without case RIP with WT
############
print(RMPdeg_datasets[n,]) #hg38
i = RMPdeg_datasets[n,'GSE']
cl = 'hESCs' #m6A-seq data was used as WT for GSE210867
RMp = c('METTL3')
RM = 'm6A'
# condition = c(rep('METTL3KD', 3), rep('WT', 3))
RMP_allDat[grep('GSE52600', RMP_allDat$study_alias ) 
           ,c("study_alias", "sample_alias", "sample_title")]
# View(allDat[[i]])
files <- list.files(file.path(data_directory,i))
files
getwd()
setwd(data_directory)
setwd(i)


# ############
# n = 23
# ############
# #GSE29714# without case RIP and Input 
# ############
# print(RMPdeg_datasets[n,]) #skip?
# i = RMPdeg_datasets[n,'GSE']
# cl = 'HEK293T'
# RMp = c('FTO')
# RM = 'm6A'
# # condition = c(rep('FTO', 3), rep('WT', 3))
# 
# RMP_allDat[grep('GSE29714', RMP_allDat$study_alias ) 
#            ,c("study_alias", "sample_alias", "sample_title")]
# 
# # View(allDat[[i]])
# files <- list.files(file.path(data_directory,i))
# files
# getwd()
# setwd(data_directory)
# setwd(i)

############
n = 24
############
#GSE97408# with RIP and Input 
############
print(RMPdeg_datasets[n,]) # hg38
i = RMPdeg_datasets[n,'GSE']
cl = c('MM6', 'NB4')# AML
RMp = c('METTL14')
RM = 'm6A'
condition_list = list( NB4_METTL14KD = c(rep('METTL14KD', 2), rep('WT', 3)),
                       MM6_METTL14KD = c(rep('METTL14KD', 2), rep('WT', 3))
) # RNA
RMP_allDat[grep('GSE97408|GSE97443', RMP_allDat$study_alias ) 
           ,c("study_alias", "sample_alias", "sample_title")]
skip <- RMP_allDat[grepl('GSE97408|GSE97443', RMP_allDat$study_alias ) &  grepl('input|RNA-Seq', RMP_allDat$sample_title)
                   , "sample_alias"]
paste0(skip$sample_alias, collapse = " ")
skip <- RMP_allDat[grepl('GSE97408|GSE97443', RMP_allDat$study_alias ) &  grepl('IP', RMP_allDat$sample_title)
                   , "sample_alias"]
paste0(skip$sample_alias, collapse = " ")

ENA <- as.data.frame(RMP_allDat[grepl('GSE97408|GSE97443', RMP_allDat$study_alias )
                                ,c("study_alias", "sample_alias", "sample_title")])
ENA[order(ENA$sample_alias),]
ENA <- ENA[order(ENA$sample_alias),]
# View(allDat[[i]])
files <- list.files(file.path(data_directory,i))
files
getwd()
setwd(data_directory)
setwd(i)
y <- ENA[grepl('NB4',ENA$sample_title) & grepl('IP',ENA$sample_title),]
y <- ENA[grepl('MM6',ENA$sample_title) & grepl('IP',ENA$sample_title),]
y <- ENA[grepl('NB4',ENA$sample_title) & grepl('RNA|input',ENA$sample_title),]
y <- ENA[grepl('MM6',ENA$sample_title) & grepl('RNA|input',ENA$sample_title),]
y
# x <- paste0('$outdir/GSE97408/', y$sample_alias, '_IP.bam')
x <- paste0('$outdir/', y$study_alias, '/', y$sample_alias, '_input.bam')
paste(x, collapse = ' ')

countfiles <- list.files('./metaPlotR/')
countfiles <- countfiles[grep('counts', countfiles)]
countfiles <- countfiles[! grepl('summary', countfiles)] 
countfiles
###################################
#  Cell number-normalized (CNN)   #
#      MeRIP INPUT analysis       #
###################################
resInP_list <- vector("list", length(condition_list))
names(resInP_list) <- names(condition_list)
countGTF_list <- vector("list", length(condition_list))
names(countGTF_list) <- names(condition_list)
for(c in names(condition_list)){
  countGTF <- read.table( paste0("metaPlotR/Input_counts_", c, ".txt"),header=T)
  dim(countGTF)
  colnames(countGTF)
  rownames(countGTF) <- countGTF$Geneid
  countGTF <- countGTF[, !colnames(countGTF) %in% c('Geneid','Chr','Start','End','Strand','Length')]
  dim(countGTF) # 28265     4
  condition <- condition_list[[c]]
  ## Sample information
  sample=data.frame(
    sample=colnames(countGTF),
    condition=condition)
  ## SummarizedExperiment input
  dds<-DESeqDataSetFromMatrix(countData=countGTF,colData=sample,design=~condition)
  ## Pre-filtering
  keep <- rowSums(counts(dds)) >= 10
  table(keep)
  dds <- dds[keep,]
  ## Note on factor levels
  dds$condition <- factor(dds$condition, levels = unique(condition))
  dds$condition <- relevel(dds$condition, ref = "WT")
  # Differential expression analysis-------------------------------------------
  dds <- DESeq(dds)
  res <-  results(dds, contrast=c("condition",unique(condition)))
  res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
  length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
  res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                                 ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
  table(res$change)
  table(is.na(res$change))
  # resOrdered <- as.data.frame(res[order(res$padj), ])
  resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
  resInP$geneSymbol <- rownames(resInP)
  table(resInP$change)
  resInP_list[[c]] <- resInP
  countGTF_list[[c]] <- countGTF
}
name = paste0(names(resInP_list), ' vs WT')
# name = paste0(cl, '_', name)
RMdeg_list[[i]] <- resInP_list 
names(RMdeg_list[[i]] ) <- paste0(gsub(' ','_',(name)))
contrast_input_list[[i]] <- lapply(countGTF_list, extract_names)
condition_input_list[[i]] <- condition_list

###################################
#  Cell number-normalized (CNN)   #
#   MeRIP RIP analysis - Part1    #
###################################
condition_list = list( NB4_METTL14KD = c(rep('METTL14KD', 2), rep('WT', 1)),
                       MM6_METTL14KD = c(rep('METTL14KD', 2), rep('WT', 1))
) # RIP 
resMeR_list <- vector("list", length(condition_list))
names(resMeR_list) <- names(condition_list)
countMatrix_list <- vector("list", length(condition_list))
names(countMatrix_list) <- names(condition_list)
resPeak_list <- vector("list",length(condition_list))
names(resPeak_list) <- names(condition_list)
for(c in names(condition_list)){
  merge_peak <- read.delim(paste0("metaPlotR/", c, "_newPeakFile.txt"),header=F)
  countMatrix <- read.table(paste0("metaPlotR/m6A_counts_", c, ".txt"),header=T)
  table(countMatrix$Geneid == merge_peak[-1,1])
  table(duplicated(merge_peak[-1,1]))
  dim(countMatrix)
  colnames(countMatrix)
  rownames(countMatrix) <- countMatrix$Geneid
  countMatrix <- countMatrix[, !colnames(countMatrix) %in% c('Geneid','Chr','Start','End','Strand','Length')]
  # 1. 做好分组因子即可
  condition <- condition_list[[c]]
  groupList<-factor(condition)
  colData <- data.frame(row.names=colnames(countMatrix), groupList=groupList)
  # 2. Exploring Method 1: DESeq2-------------------------------------------
  # 构建 DESeqDataSet 对象
  dds <- DESeqDataSetFromMatrix(countData = countMatrix,
                                colData = colData,
                                design = ~ groupList)
  ## Pre-filtering
  keep <- rowSums(counts(dds)) >= 1
  table(keep)
  dds <- dds[keep,]
  # 接着计算每个样本的 library size 用于后续的标准化。library size 等于落在该样本所有peaks上的reads的总数，即 counts 矩阵中每列的加和。如果bFullLibrarySize设为TRUE，则会使用library的总reads数（根据原BAM/BED文件统计获得）。然后使用estimateDispersions估计统计分布，需要将参数fitType设为’local’
  dds<-estimateSizeFactors(dds)
  dds<-estimateDispersions(dds,fitType="local")
  dds <- nbinomWaldTest(dds)
  res <- results(dds)
  ## 提取你想要的差异分析结果，我们这里是treated组对untreated组进行比较
  res <-  results(dds, contrast=c("groupList", unique(condition)))
  res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
  length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
  res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                                 ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
  res_DESeq2 = as.data.frame(res)
  table(res_DESeq2$change)
  res_DESeq2 = res_DESeq2[res_DESeq2$change %in% c('UP', 'DOWN'), ]
  # 2. Exploring Method 2: limma-------------------------------------------
  library(edgeR)
  d0 <- DGEList(counts = countMatrix)
  dim(d0)
  # Filtering low count peaks
  cutoff <- 1           
  drop <- which(apply(cpm(d0), 1, max) < cutoff)
  print(length(drop))
  if (length(drop) >0) {
    d <- d0[-drop,] 
  } else{ d <- d0}
  dim(d) #number of peaks left: 8501      
  d <- calcNormFactors(d, method = c('TMM'))
  colnames(countMatrix)
  genotype <- groupList
  table(genotype)
  mm <- model.matrix(~0 + genotype)#
  v <- voom(d, mm, plot=F, normalize="quantile")#voom from limma:Transform count data to log2-counts per million (logtpm)
  fit <- lmFit(v, mm)
  head(coef(fit))
  contr <- makeContrasts(paste0(colnames(coef(fit)), collapse  = ' - '), levels = colnames(coef(fit))) #Change if other comparison
  tmp <- contrasts.fit(fit, contr)
  tmp <- eBayes(tmp)
  top.table <- topTable(tmp, sort.by = "P", n = Inf)
  head(top.table, 20)
  length(which(top.table$adj.P.Val < 0.05 & abs(top.table$logFC) > log2(1.5))) #how many are significant?
  res_limma <- as.data.frame(top.table)
  res_limma$change <- with(res_limma, ifelse(res_limma$adj.P.Val < 0.05 & abs(res_limma$logFC) > log2(1.5),
                                             ifelse(res_limma$logFC > log2(1.5),'UP','DOWN'),'NOT'))
  table(res_limma$change)
  res_limma = res_limma[res_limma$change %in% c('UP', 'DOWN'), ]
  resMeR_list[[c]] <- list(res_DESeq2, res_limma)
  countMatrix_list[[c]] <- countMatrix
  resPeak_list[[c]] <- merge_peak
  name = res@elementMetadata$description[2]
  # name = paste0(cl, '_',str_split(name, "groupList ", simplify = T)[,2])
  name = str_split(name, "groupList ", simplify = T)[,2]
  name = paste0(gsub(' ','_',(name)))
  names(resMeR_list[[c]] ) <- paste0(name, '_', c('DESeq2', 'limma'))
}
RMdem_list[[i]] <- as.list(unlist(resMeR_list, recursive = FALSE))
name <- gsub("_METTL14KD.", '_', names(RMdem_list[[i]] ))
names(RMdem_list[[i]] ) <- name
RMdep_list[[i]] <- resPeak_list 
name <- unique(gsub('(_limma|_DESeq2)', '', name))
names(RMdep_list[[i]] ) <- name
contrast_IP_list[[i]] <- lapply(countMatrix_list, extract_names)
condition_IP_list[[i]] <- condition_list

############
n = 25
############
#GSE120229# with RIP and Input 
############
print(RMPdeg_datasets[n,]) 
i = RMPdeg_datasets[n,'GSE']
cl = c('HEK293T')
RMp = c('PCIF1')
RM = 'm6Am'
condition = c(rep('PCIF1KD', 2), rep('WT', 2))
RMP_allDat[grep('GSE120229', RMP_allDat$study_alias ) 
           ,c("study_alias", "sample_alias", "sample_title")]
skip <- RMP_allDat[grepl('GSE120229', RMP_allDat$study_alias ) &  grepl('input', RMP_allDat$sample_title)
                   , "sample_alias"]
paste0(skip$sample_alias, collapse = " ")
skip <- RMP_allDat[grepl('GSE120229', RMP_allDat$study_alias ) &  grepl('IP', RMP_allDat$sample_title)
                   , "sample_alias"]
paste0(skip$sample_alias, collapse = " ")
ENA <- as.data.frame(RMP_allDat[grepl('GSE120229', RMP_allDat$study_alias )
                                ,c("study_alias", "sample_alias", "sample_title")])
ENA[order(ENA$sample_alias),]
ENA <- ENA[order(ENA$sample_alias),]

# View(allDat[[i]])
files <- list.files(file.path(data_directory,i))
files
getwd()
setwd(data_directory)
setwd(i)
c <- unique(condition[!condition %in% 'WT'])
###################################
#  Cell number-normalized (CNN)   #
#      MeRIP INPUT analysis       #
###################################
countGTF <- read.table( paste0("metaPlotR/Input_counts_", c, ".txt"),header=T)
rownames(countGTF) <- countGTF$Geneid
countGTF <- countGTF[, !colnames(countGTF) %in% c('Geneid','Chr','Start','End','Strand','Length')]
dim(countGTF) #28265     4
## Sample information
sample=data.frame(
  sample=colnames(countGTF),
  condition=condition)
## SummarizedExperiment input
dds<-DESeqDataSetFromMatrix(countData=countGTF,colData=sample,design=~condition)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 10
table(keep)
dds <- dds[keep,]
## Note on factor levels
dds$condition <- factor(dds$condition, levels = unique(condition))
dds$condition <- relevel(dds$condition, ref = "WT")
# Differential expression analysis-------------------------------------------
dds <- DESeq(dds)

res <-  results(dds, contrast=c("condition",unique(condition)))
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
table(res$change)
table(is.na(res$change))
# resOrdered <- as.data.frame(res[order(res$padj), ])
resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
resInP$geneSymbol <- rownames(resInP)
name = res@elementMetadata$description[2]
name = paste0(cl, '_',str_split(name, "condition ", simplify = T)[,2])
RMdeg_list[[i]] <- list( resInP)
names(RMdeg_list[[i]] ) <- paste0(gsub(' ','_',(name)))
c <- unique(condition[!condition %in% 'WT'])
countGTF_list <- list(colnames(countGTF)); names(countGTF_list) <- c
condition_list <- list(condition); names(condition_list) <- c
contrast_input_list[[i]] <- countGTF_list
condition_input_list[[i]] <- condition_list

###################################
#  Cell number-normalized (CNN)   #
#   MeRIP RIP analysis - Part1    #
###################################
## Sample information
merge_peak <- read.delim(paste0("metaPlotR/", c, "_newPeakFile.txt"),header=F)
countMatrix <- read.table(paste0("metaPlotR/m6Am_counts_", c, ".txt"),header=T)
table(countMatrix$Geneid == merge_peak[-1,1])
table(duplicated(merge_peak$V1[-1]))
dim(countMatrix)
colnames(countMatrix)
rownames(countMatrix) <- countMatrix$Geneid
countMatrix <- countMatrix[, !colnames(countMatrix) %in% c('Geneid','Chr','Start','End','Strand','Length')]
# 1. 做好分组因子即可
groupList<-factor(condition)
colData <- data.frame(row.names=colnames(countMatrix), groupList=groupList)
# 2. Exploring Method 1: DESeq2-------------------------------------------
# 构建 DESeqDataSet 对象
dds <- DESeqDataSetFromMatrix(countData = countMatrix,
                              colData = colData,
                              design = ~ groupList)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 1
table(keep)
dds <- dds[keep,]
# 接着计算每个样本的 library size 用于后续的标准化。library size 等于落在该样本所有peaks上的reads的总数，即 counts 矩阵中每列的加和。如果bFullLibrarySize设为TRUE，则会使用library的总reads数（根据原BAM/BED文件统计获得）。然后使用estimateDispersions估计统计分布，需要将参数fitType设为’local’
dds<-estimateSizeFactors(dds)
dds<-estimateDispersions(dds,fitType="local")
dds <- nbinomWaldTest(dds)
res <- results(dds)
## 提取你想要的差异分析结果，我们这里是treated组对untreated组进行比较
res <-  results(dds, contrast=c("groupList", unique(condition)))
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
res_DESeq2 = as.data.frame(res)
table(res_DESeq2$change)#Not only
res_DESeq2 = res_DESeq2[res_DESeq2$change %in% c('UP', 'DOWN'), ]
# 2. Exploring Method 2: limma-------------------------------------------
library(edgeR)
d0 <- DGEList(counts = countMatrix)
dim(d0)
# Filtering low count peaks
cutoff <- 1           
drop <- which(apply(cpm(d0), 1, max) < cutoff)
print(length(drop))
if (length(drop) >0) {
  d <- d0[-drop,] 
} else{ d <- d0}
dim(d) #number of peaks left: 8501      
d <- calcNormFactors(d, method = c('TMM'))
colnames(countMatrix)
genotype <- groupList
table(genotype)
mm <- model.matrix(~0 + genotype)#
v <- voom(d, mm, plot=F, normalize="quantile")#voom from limma:Transform count data to log2-counts per million (logtpm)
fit <- lmFit(v, mm)
head(coef(fit))
contr <- makeContrasts(paste0(colnames(coef(fit)), collapse  = ' - '), levels = colnames(coef(fit))) #Change if other comparison
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)
length(which(top.table$adj.P.Val < 0.05 & abs(top.table$logFC) > log2(1.5))) #how many are significant?
res_limma <- as.data.frame(top.table)
res_limma$change <- with(res_limma, ifelse(res_limma$adj.P.Val < 0.05 & abs(res_limma$logFC) > log2(1.5),
                                           ifelse(res_limma$logFC > log2(1.5),'UP','DOWN'),'NOT'))
table(res_limma$change)#NOT only
res_limma = res_limma[res_limma$change %in% c('UP', 'DOWN'), ]

RMdem_list[[i]] <- list(res_DESeq2, res_limma)
name = res@elementMetadata$description[2]
name = paste0(cl, '_',str_split(name, "groupList ", simplify = T)[,2])
name= paste0(gsub(' ','_',(name)))
names(RMdem_list[[i]] ) <- paste0(name, '_', c('DESeq2', 'limma'))
RMdep_list[[i]] <- list( merge_peak )
names(RMdep_list[[i]] ) <- name
countMatrix_list <- list(colnames(countMatrix)); names(countMatrix_list) <- c
condition_list <- list(condition); names(condition_list) <- c
contrast_IP_list[[i]] <- countMatrix_list
condition_IP_list[[i]] <- condition_list
# ############
# n = 26
# ############
# #GSE83438# without case RIP and Input 
# ############
# print(RMPdeg_datasets[n,]) #skip?
# i = RMPdeg_datasets[n,'GSE']
# cl = c('Huh7')
# RMp = c('METTL3','METTL14','FTO')
# RM = 'm6A'
# # condition = c(rep('METTL3+14', 3), rep('WT', 3))
# # condition = c(rep('FTO', 3), rep('WT', 3))
# RMP_allDat[grep('GSE83438', RMP_allDat$study_alias ) 
#            ,c("study_alias", "sample_alias", "sample_title")]
# RMP_allDat[grepl('GSE83438', RMP_allDat$study_alias ) & grepl('_MeRIP', RMP_allDat$sample_title)
#            ,c("study_alias", "sample_alias", "sample_title")]
# skip <- RMP_allDat[grepl('GSE83438', RMP_allDat$study_alias ) &  grepl('Input', RMP_allDat$sample_title)
#                    , "sample_alias"]
# paste0(skip$sample_alias, collapse = " ")
# skip <- RMP_allDat[grepl('GSE83438', RMP_allDat$study_alias ) &  grepl('IP', RMP_allDat$sample_title)
#                    , "sample_alias"]
# paste0(skip$sample_alias, collapse = " ")
# 
# # View(allDat[[i]])
# # setwd(data_directory)
# # gse <- getGEO('GSE83438',destdir=".", AnnotGPL=F, getGPL=F)
# # pdata <- pData(gse[[2]])
# # allDat[[i]] <- rbind(allDat[[i]], pdata)
# files <- list.files(file.path(data_directory,i))
# files
# getwd()
# setwd(data_directory)
# setwd(i)
# 
# 
# ############
# n = 27
# ############
# #GSE70299# without case RIP and Input 
# ############
# print(RMPdeg_datasets[n,]) #skip?
# i = RMPdeg_datasets[n,'GSE']
# cl = c('Lymphoblastoid')
# RMp = c('METTL3','METTL14','WTAP')
# RM = 'm6A'
# # condition = c(rep('METTL3+14', 3), rep('WT', 3))
# RMP_allDat[grep('GSE70299', RMP_allDat$study_alias ) 
#            ,c("study_alias", "sample_alias", "sample_title")]
# 
# # View(allDat[[i]])
# # H:human C/R: Other species
# files <- list.files(file.path(data_directory,i))
# files
# getwd()
# setwd(data_directory)
# setwd(i)

############
n = 28
############
#GSE55572# with RIP and Input 
############
print(RMPdeg_datasets[n,]) #
i = RMPdeg_datasets[n,'GSE']
cl = c('HEK293','A549')
RMp = c('METTL3','METTL14','WTAP','KIAA1429')
RM = 'm6A'
condition_list = list(
  HEK293_WTAPKD = c(rep('WTAPKD', 1), rep('WT', 1)),
  HEK293_METTL3KD = c(rep('METTL3KD', 1), rep('WT', 1)),
  
  A549_WTAPKD = c(rep('WTAPKD', 3), rep('WT', 4)),
  A549_METTL14KD = c(rep('METTL14KD', 5), rep('WT', 4)),
  A549_METTL3KD = c(rep('METTL3KD', 3), rep('WT', 4)),
  A549_KIAA1429KD = c(rep('KIAA1429KD', 1), rep('WT', 4)),
  A549_METTL3_14KD = c(rep('METTL3_14KD', 1), rep('WT', 4))
)
RMP_allDat[grep('GSE55572', RMP_allDat$study_alias ) 
           ,c("study_alias", "sample_alias", "sample_title")]
RMP_allDat[grepl('GSE55572', RMP_allDat$study_alias ) & (!grepl('human_h|human_OKMS', RMP_allDat$sample_title ))
           ,c("study_alias", "sample_alias", "sample_title")]


# View(allDat[[i]])
# pdata <- pData(gse[[2]])
# allDat[[i]] <- pdata
ENA <- as.data.frame(RMP_allDat[grepl('GSE55572', RMP_allDat$study_alias ) & (!grepl('human_h|human_OKMS', RMP_allDat$sample_title ))
                                ,c("study_alias", "sample_alias", "sample_title")])
ENA[order(ENA$sample_alias),]
ENA <- ENA[order(ENA$sample_alias),]
# View(allDat[[i]])
files <- list.files(file.path(data_directory,i))
files
getwd()
setwd(data_directory)
setwd(i)
y <- ENA[grepl('NB4',ENA$sample_title) & grepl('IP',ENA$sample_title),]
y <- ENA[grepl('MM6',ENA$sample_title) & grepl('IP',ENA$sample_title),]
y <- ENA[grepl('NB4',ENA$sample_title) & grepl('RNA|input',ENA$sample_title),]
y <- ENA[grepl('MM6',ENA$sample_title) & grepl('RNA|input',ENA$sample_title),]
y
# x <- paste0('$outdir/GSE97408/', y$sample_alias, '_IP.bam')
x <- paste0('$outdir/', y$study_alias, '/', y$sample_alias, '_input.bam')
paste(x, collapse = ' ')

countfiles <- list.files('./metaPlotR/')
countfiles <- countfiles[grep('counts', countfiles)]
countfiles <- countfiles[! grepl('summary', countfiles)] 
countfiles
###################################
#  Cell number-normalized (CNN)   #
#      MeRIP INPUT analysis       #
###################################
resInP_list <- vector("list", length(condition_list))
names(resInP_list) <- names(condition_list)
countGTF_list <- vector("list", length(condition_list))
names(countGTF_list) <- names(condition_list)
for(c in names(condition_list)){
  countGTF <- read.table( paste0("metaPlotR/Input_counts_", c, ".txt"),header=T)
  dim(countGTF)
  colnames(countGTF)
  rownames(countGTF) <- countGTF$Geneid
  countGTF <- countGTF[, !colnames(countGTF) %in% c('Geneid','Chr','Start','End','Strand','Length')]
  dim(countGTF) # 28265     4
  condition <- condition_list[[c]]
  if(c %in% c("HEK293_WTAPKD","HEK293_METTL3KD")){
    # 1. Differential gene analysis with edgeR
    # Without replicates, you cannot estimate which genes are differentially expressed using EdgeR or DESeq2. You can only calculate fold changes based on normalized read counts (preferentially CPM normalized by TMM method included in any EdgeR analysis) and apply a stringent fold-change cut-off to determine which genes are more or less expressed depending on the condition.
    y <- DGEList(countGTF)
    # Calculate fold change
    y <- calcNormFactors(y, method = c('TMM'))
    y$samples$lib.size <- colSums(y$counts)
    logcpm <- cpm(y, prior.count=2, log=TRUE)
    logcpm[1:5,]
    cpm <- cpm(y, prior.count=2, log=FALSE)
    cpm[1:5, ]
    fold_change <- cpm[,1]/cpm[,2]
    res <- countGTF
    res$logFC <- log2(fold_change)
    res$FoldChange <- fold_change
    res$change <- with(res, ifelse(abs(res$logFC) > log2(1.5),
                                   ifelse(res$logFC > log2(1.5),'UP','DOWN'),'NOT'))
    table(res$change)
    resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
    resInP$geneSymbol <- rownames(resInP)
    table(resInP$change)
    resInP_list[[c]] <- resInP
    countGTF_list[[c]] <- countGTF
    
  } else{
    
    ## Sample information
    sample=data.frame(
      sample=colnames(countGTF),
      condition=condition)
    ## SummarizedExperiment input
    dds<-DESeqDataSetFromMatrix(countData=countGTF,colData=sample,design=~condition)
    ## Pre-filtering
    keep <- rowSums(counts(dds)) >= 10
    table(keep)
    dds <- dds[keep,]
    ## Note on factor levels
    dds$condition <- factor(dds$condition, levels = unique(condition))
    dds$condition <- relevel(dds$condition, ref = "WT")
    # Differential expression analysis-------------------------------------------
    dds <- DESeq(dds)
    res <-  results(dds, contrast=c("condition",unique(condition)))
    res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
    length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
    res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                                   ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
    table(res$change)
    table(is.na(res$change))
    # resOrdered <- as.data.frame(res[order(res$padj), ])
    resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
    resInP$geneSymbol <- rownames(resInP)
    table(resInP$change)
    resInP_list[[c]] <- resInP
    countGTF_list[[c]] <- countGTF
  }
  
}
name = paste0(names(resInP_list), ' vs WT')
RMdeg_list[[i]] <- resInP_list 
names(RMdeg_list[[i]] ) <- paste0(gsub(' ','_',(name)))
contrast_input_list[[i]] <- lapply(countGTF_list, extract_names)
condition_input_list[[i]] <- condition_list

###################################
#  Cell number-normalized (CNN)   #
#   MeRIP RIP analysis - Part1    #
###################################
resMeR_list <- vector("list", length(condition_list))
names(resMeR_list) <- names(condition_list)
countMatrix_list <- vector("list", length(condition_list))
names(countMatrix_list) <- names(condition_list)
resPeak_list <- vector("list",length(condition_list))
names(resPeak_list) <- names(condition_list)
for(c in names(condition_list)){
  merge_peak <- read.delim(paste0("metaPlotR/", c, "_newPeakFile.txt"),header=F)
  countMatrix <- read.table(paste0("metaPlotR/m6A_counts_", c, ".txt"),header=T)
  table(countMatrix$Geneid == merge_peak[-1,1])
  table(duplicated(merge_peak[-1,1]))
  dim(countMatrix)
  colnames(countMatrix)
  rownames(countMatrix) <- countMatrix$Geneid
  countMatrix <- countMatrix[, !colnames(countMatrix) %in% c('Geneid','Chr','Start','End','Strand','Length')]
  # 1. 做好分组因子即可
  condition <- condition_list[[c]]
  groupList<-factor(condition)
  colData <- data.frame(row.names=colnames(countMatrix), groupList=groupList)
  if(c %in% c("HEK293_WTAPKD","HEK293_METTL3KD")){
    # 2. Differential gene analysis with edgeR
    # Without replicates, you cannot estimate which genes are differentially expressed using EdgeR or DESeq2. You can only calculate fold changes based on normalized read counts (preferentially CPM normalized by TMM method included in any EdgeR analysis) and apply a stringent fold-change cut-off to determine which genes are more or less expressed depending on the condition.
    y <- DGEList(countMatrix)
    # Calculate fold change
    y <- calcNormFactors(y, method = c('TMM'))
    y$samples$lib.size <- colSums(y$counts)
    logcpm <- cpm(y, prior.count=2, log=TRUE)
    logcpm[1:5,]
    cpm <- cpm(y, prior.count=2, log=FALSE)
    cpm[1:5, ]
    # fold_change <- logcpm[,'Ratio1']/logcpm[,'Ratio2']
    fold_change <- cpm[,1]/cpm[,2]
    resMeR <- countMatrix
    resMeR$logFC <- log2(fold_change)
    resMeR$FoldChange <- fold_change
    resMeR$change <- with(resMeR, ifelse(abs(resMeR$logFC) > log2(1.5),
                                         ifelse(resMeR$logFC > log2(1.5),'UP','DOWN'),'NOT'))
    table(resMeR$change)
    resMeR = as.data.frame(resMeR[resMeR$change %in% c('UP', 'DOWN'), ])
    resMeR_list[[c]] <- list(resMeR)
    countMatrix_list[[c]] <- countMatrix
    resPeak_list[[c]] <- merge_peak
    name = paste0(c, ' vs WT')
    name = paste0(gsub(' ','_', name ))
    names(resMeR_list[[c]]) <- paste0(name, '_', c('limma'))
    
  } else{
    # 2. Exploring Method 1: DESeq2-------------------------------------------
    # 构建 DESeqDataSet 对象
    dds <- DESeqDataSetFromMatrix(countData = countMatrix,
                                  colData = colData,
                                  design = ~ groupList)
    ## Pre-filtering
    keep <- rowSums(counts(dds)) >= 1
    table(keep)
    dds <- dds[keep,]
    # 接着计算每个样本的 library size 用于后续的标准化。library size 等于落在该样本所有peaks上的reads的总数，即 counts 矩阵中每列的加和。如果bFullLibrarySize设为TRUE，则会使用library的总reads数（根据原BAM/BED文件统计获得）。然后使用estimateDispersions估计统计分布，需要将参数fitType设为’local’
    dds<-estimateSizeFactors(dds)
    dds<-estimateDispersions(dds,fitType="local")
    dds <- nbinomWaldTest(dds)
    res <- results(dds)
    ## 提取你想要的差异分析结果，我们这里是treated组对untreated组进行比较
    res <-  results(dds, contrast=c("groupList", unique(condition)))
    res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
    length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
    res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                                   ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
    res_DESeq2 = as.data.frame(res)
    table(res_DESeq2$change)
    res_DESeq2 = res_DESeq2[res_DESeq2$change %in% c('UP', 'DOWN'), ]
    # 2. Exploring Method 2: limma-------------------------------------------
    library(edgeR)
    d0 <- DGEList(counts = countMatrix)
    dim(d0)
    # Filtering low count peaks
    cutoff <- 1           
    drop <- which(apply(cpm(d0), 1, max) < cutoff)
    print(length(drop))
    if (length(drop) >0) {
      d <- d0[-drop,] 
    } else{ d <- d0}
    dim(d) #number of peaks left: 8501      
    d <- calcNormFactors(d, method = c('TMM'))
    colnames(countMatrix)
    genotype <- groupList
    table(genotype)
    mm <- model.matrix(~0 + genotype)#
    v <- voom(d, mm, plot=F, normalize="quantile")#voom from limma:Transform count data to log2-counts per million (logtpm)
    fit <- lmFit(v, mm)
    head(coef(fit))
    contr <- makeContrasts(paste0(colnames(coef(fit)), collapse  = ' - '), levels = colnames(coef(fit))) #Change if other comparison
    tmp <- contrasts.fit(fit, contr)
    tmp <- eBayes(tmp)
    top.table <- topTable(tmp, sort.by = "P", n = Inf)
    head(top.table, 20)
    length(which(top.table$adj.P.Val < 0.05 & abs(top.table$logFC) > log2(1.5))) #how many are significant?
    res_limma <- as.data.frame(top.table)
    res_limma$change <- with(res_limma, ifelse(res_limma$adj.P.Val < 0.05 & abs(res_limma$logFC) > log2(1.5),
                                               ifelse(res_limma$logFC > log2(1.5),'UP','DOWN'),'NOT'))
    table(res_limma$change)
    res_limma = res_limma[res_limma$change %in% c('UP', 'DOWN'), ]
    resMeR_list[[c]] <- list(res_DESeq2, res_limma)
    countMatrix_list[[c]] <- countMatrix
    resPeak_list[[c]] <- merge_peak
    name = paste0(c, ' vs WT')
    name = paste0(gsub(' ','_', name ))
    names(resMeR_list[[c]] ) <- paste0(name, '_', c('DESeq2', 'limma'))
    
  }
}
names(resMeR_list) <- NULL
RMdem_list[[i]] <- as.list(unlist(resMeR_list, recursive = FALSE))
names(RMdem_list[[i]] )
RMdep_list[[i]] <- resPeak_list 
name <- unique(gsub('(_limma|_DESeq2)', '', names(RMdem_list[[i]] )))
names(RMdep_list[[i]] ) <- name
contrast_IP_list[[i]] <- lapply(countMatrix_list, extract_names)
condition_IP_list[[i]] <- condition_list

############
n = 29
############
#GSE37002# without case RIP, with case RNA 
############
print(RMPdeg_datasets[n,]) #
i = RMPdeg_datasets[n,'GSE']
cl = c('HepG2')
RMp = c('METTL3')
RM = 'm6A'
condition = c(rep('METTL3KD', 2), rep('WT', 2))
RMP_allDat[grep('GSE37002', RMP_allDat$study_alias ) 
           ,c("study_alias", "sample_alias", "sample_title")]
# View(allDat[[i]])
files <- list.files(file.path(data_directory,i))
files
getwd()
setwd(data_directory)
setwd(i)
countfiles <- files[grep('txt', files)]
countfiles <- countfiles[! grepl('summary', countfiles)] 
countfiles
###################################
#      RNA/Ribo-seq analysis      #
###################################
df_in <- read.table(file.path(data_directory,i, 'GSE37001_DEgenes.txt.gz'), header = T)
# df_in <- read.csv(file.path(data_directory,i, 'GSE37002_Supplementary_data_HepG2_treatments.txt.gz'),
#                   sep = '\t',header = T)
# df_in <- read.table(file.path(data_directory,i, 'GSE37001_DEexons.txt.gz'), header = T)
## Gene name
library(data.table)
id_mapped <- read.table("/data2/rluo4/hg38/gencode.v44lift37.annotation.gene.probeMap",header = T)
# id_mapped <- read.table("/data2/rluo4/hg38/gencode.v36.annotation.gene.probeMap",header = T)
dim(id_mapped) #  60675     6
table(str_split(df_in$id, '[.]', simplify = T)[,1] %in% gene_info$geneEnsembls)
table(df_in$id %in% id_mapped$id)
table(str_split(df_in$id, '[.]', simplify = T)[,1] %in% str_split(id_mapped$id, '[.]', simplify = T)[,1])
df_in$id <- str_split(df_in$id, '[.]', simplify = T)[,1]
id_mapped$id <- str_split(id_mapped$id, '[.]', simplify = T)[,1]
id_mapped <- id_mapped[which(id_mapped$id %in% df_in$id),]
dim(id_mapped) # 51343     6 for v36 # 54667     6 for v44 
table(duplicated(id_mapped$gene))
res <- left_join(id_mapped,df_in,by="id")
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
table(duplicated(res$gene))
tmp=by(res,res$gene,function(x) rownames(x)[which.max(x$baseMean)])
probes = as.character(tmp)
length(probes) # 7785
resInP = as.data.frame(res[rownames(res) %in% probes, ])
resInP$geneSymbol <- resInP$gene
name = 'METTL3KD vs WT'
name = paste0(cl, '_', name)
RMdeg_list[[i]] <- list( resInP)
names(RMdeg_list[[i]] ) <- paste0(gsub(' ','_',(name)))


############
n = 30
############
#GSE46705# with RIP and Input 
############
print(RMPdeg_datasets[n,]) # hg38 (paper used: hg18; I used: hg38 on server 57)
i = RMPdeg_datasets[n,'GSE']
cl = c('Hela')
RMp = c('METTL3','METTL14','WTAP')
RM = 'm6A'
condition = c(rep('METTL3OE', 2), rep('WT', 2))
condition = c(rep('METTL3OE', 2), rep('WT', 2))
condition = c(rep('WTAPOE', 2), rep('WT', 2))
condition = c(rep('METTL14OE', 2), rep('WT', 2))
condition_list = list(
  METTL14KD = c(rep('METTL14KD', 2), rep('WT', 4)),
  METTL3KD = c(rep('METTL3KD', 2), rep('WT', 4)),
  WTAPKD = c(rep('WTAPKD', 2), rep('WT', 4))
)
RMP_allDat[grep('GSE46705', RMP_allDat$study_alias ) 
           ,c("study_alias", "sample_alias", "sample_title")]
RMP_allDat[grepl('GSE46705', RMP_allDat$study_alias ) & (!grepl('_parClip', RMP_allDat$sample_title))
           ,c("study_alias", "sample_alias", "sample_title")]
skip = RMP_allDat[grepl('GSE46705', RMP_allDat$study_alias ) & (!grepl('_parClip', RMP_allDat$sample_title))
                  ,c("study_alias", "sample_alias", "sample_title")]
paste0(skip$sample_alias, collapse = "|")

skip <- RMP_allDat[grepl('GSE46705', RMP_allDat$study_alias ) & (!grepl('_parClip', RMP_allDat$sample_title))
                   &  grepl('Input', RMP_allDat$sample_title)
                   , "sample_alias"]
paste0(skip$sample_alias, collapse = " ")
skip <- RMP_allDat[grepl('GSE46705', RMP_allDat$study_alias ) &  grepl('IP', RMP_allDat$sample_title)
                   , "sample_alias"]
paste0(skip$sample_alias, collapse = " ")

# View(allDat[[i]])

files <- list.files(file.path(data_directory,i))
files
getwd()
setwd(data_directory)
setwd(i)
ENA <- as.data.frame(RMP_allDat[grepl('GSE46705', RMP_allDat$study_alias ) & (!grepl('_parClip', RMP_allDat$sample_title))
                                ,c("study_alias", "sample_alias", "sample_title")])
ENA[order(ENA$sample_alias),]
ENA <- ENA[order(ENA$sample_alias),]
# View(allDat[[i]])
files <- list.files(file.path(data_directory,i))
files
getwd()
setwd(data_directory)
setwd(i)
y <- ENA[grepl('M14|C',ENA$sample_title) & grepl('IP',ENA$sample_title),]
y <- ENA[grepl('M3|C',ENA$sample_title) & grepl('IP',ENA$sample_title),]
y <- ENA[grepl('WTAP|C',ENA$sample_title) & grepl('IP',ENA$sample_title),]

y <- ENA[grepl('M14|C',ENA$sample_title) & grepl('Input',ENA$sample_title),]
y <- ENA[grepl('M3|C',ENA$sample_title) & grepl('Input',ENA$sample_title),]
y <- ENA[grepl('WTAP|C',ENA$sample_title) & grepl('Input',ENA$sample_title),]
y
# x <- paste0('$outdir/', y$study_alias, '/', y$sample_alias, '_IP.bam')
x <- paste0('$outdir/', y$study_alias, '/', y$sample_alias, '_input.bam')
paste(x, collapse = ' ')

countfiles <- list.files('./metaPlotR/')
countfiles <- countfiles[grep('counts', countfiles)]
countfiles <- countfiles[! grepl('summary', countfiles)] 
countfiles
###################################
#  Cell number-normalized (CNN)   #
#      MeRIP INPUT analysis       #
###################################
resInP_list <- vector("list", length(condition_list))
names(resInP_list) <- names(condition_list)
countGTF_list <- vector("list", length(condition_list))
names(countGTF_list) <- names(condition_list)
for(c in names(condition_list)){
  countGTF <- read.table( paste0("metaPlotR/Input_counts_", c, ".txt"),header=T)
  dim(countGTF)
  colnames(countGTF)
  rownames(countGTF) <- countGTF$Geneid
  countGTF <- countGTF[, !colnames(countGTF) %in% c('Geneid','Chr','Start','End','Strand','Length')]
  dim(countGTF) # 28265     4
  condition <- condition_list[[c]]
  ## Sample information
  sample=data.frame(
    sample=colnames(countGTF),
    condition=condition)
  ## SummarizedExperiment input
  dds<-DESeqDataSetFromMatrix(countData=countGTF,colData=sample,design=~condition)
  ## Pre-filtering
  keep <- rowSums(counts(dds)) >= 10
  table(keep)
  dds <- dds[keep,]
  ## Note on factor levels
  dds$condition <- factor(dds$condition, levels = unique(condition))
  dds$condition <- relevel(dds$condition, ref = "WT")
  # Differential expression analysis-------------------------------------------
  dds <- DESeq(dds)
  res <-  results(dds, contrast=c("condition",unique(condition)))
  res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
  length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
  res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                                 ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
  table(res$change)
  table(is.na(res$change))
  # resOrdered <- as.data.frame(res[order(res$padj), ])
  resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
  resInP$geneSymbol <- rownames(resInP)
  table(resInP$change)
  resInP_list[[c]] <- resInP
  countGTF_list[[c]] <- countGTF 
}

name = paste0(cl, '_',names(resInP_list), ' vs WT')
RMdeg_list[[i]] <- resInP_list 
names(RMdeg_list[[i]] ) <- paste0(gsub(' ','_',(name)))
contrast_input_list[[i]] <- lapply(countGTF_list, extract_names)
condition_input_list[[i]] <- condition_list

###################################
#  Cell number-normalized (CNN)   #
#   MeRIP RIP analysis - Part1    #
###################################
resMeR_list <- vector("list", length(condition_list))
names(resMeR_list) <- names(condition_list)
countMatrix_list <- vector("list", length(condition_list))
names(countMatrix_list) <- names(condition_list)
resPeak_list <- vector("list",length(condition_list))
names(resPeak_list) <- names(condition_list)
for(c in names(condition_list)){
  merge_peak <- read.delim(paste0("metaPlotR/", c, "_newPeakFile.txt"),header=F)
  countMatrix <- read.table(paste0("metaPlotR/m6A_counts_", c, ".txt"),header=T)
  table(countMatrix$Geneid == merge_peak[-1,1])
  table(duplicated(merge_peak[-1,1]))
  dim(countMatrix)
  colnames(countMatrix)
  rownames(countMatrix) <- countMatrix$Geneid
  countMatrix <- countMatrix[, !colnames(countMatrix) %in% c('Geneid','Chr','Start','End','Strand','Length')]
  # 1. 做好分组因子即可
  condition <- condition_list[[c]]
  groupList<-factor(condition)
  colData <- data.frame(row.names=colnames(countMatrix), groupList=groupList)
  
  # 2. Exploring Method 1: DESeq2-------------------------------------------
  # 构建 DESeqDataSet 对象
  dds <- DESeqDataSetFromMatrix(countData = countMatrix,
                                colData = colData,
                                design = ~ groupList)
  ## Pre-filtering
  keep <- rowSums(counts(dds)) >= 1
  table(keep)
  dds <- dds[keep,]
  # 接着计算每个样本的 library size 用于后续的标准化。library size 等于落在该样本所有peaks上的reads的总数，即 counts 矩阵中每列的加和。如果bFullLibrarySize设为TRUE，则会使用library的总reads数（根据原BAM/BED文件统计获得）。然后使用estimateDispersions估计统计分布，需要将参数fitType设为’local’
  dds<-estimateSizeFactors(dds)
  dds<-estimateDispersions(dds,fitType="local")
  dds <- nbinomWaldTest(dds)
  res <- results(dds)
  ## 提取你想要的差异分析结果，我们这里是treated组对untreated组进行比较
  res <-  results(dds, contrast=c("groupList", unique(condition)))
  res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
  length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
  res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                                 ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
  res_DESeq2 = as.data.frame(res)
  table(res_DESeq2$change)
  res_DESeq2 = res_DESeq2[res_DESeq2$change %in% c('UP', 'DOWN'), ]
  # 2. Exploring Method 2: limma-------------------------------------------
  library(edgeR)
  d0 <- DGEList(counts = countMatrix)
  dim(d0)
  # Filtering low count peaks
  cutoff <- 1           
  drop <- which(apply(cpm(d0), 1, max) < cutoff)
  print(length(drop))
  if (length(drop) >0) {
    d <- d0[-drop,] 
  } else{ d <- d0}
  dim(d) #number of peaks left: 8501      
  d <- calcNormFactors(d, method = c('TMM'))
  colnames(countMatrix)
  genotype <- groupList
  table(genotype)
  mm <- model.matrix(~0 + genotype)#
  v <- voom(d, mm, plot=F, normalize="quantile")#voom from limma:Transform count data to log2-counts per million (logtpm)
  fit <- lmFit(v, mm)
  head(coef(fit))
  contr <- makeContrasts(paste0(colnames(coef(fit)), collapse  = ' - '), levels = colnames(coef(fit))) #Change if other comparison
  tmp <- contrasts.fit(fit, contr)
  tmp <- eBayes(tmp)
  top.table <- topTable(tmp, sort.by = "P", n = Inf)
  head(top.table, 20)
  length(which(top.table$adj.P.Val < 0.05 & abs(top.table$logFC) > log2(1.5))) #how many are significant?
  res_limma <- as.data.frame(top.table)
  res_limma$change <- with(res_limma, ifelse(res_limma$adj.P.Val < 0.05 & abs(res_limma$logFC) > log2(1.5),
                                             ifelse(res_limma$logFC > log2(1.5),'UP','DOWN'),'NOT'))
  table(res_limma$change)
  res_limma = res_limma[res_limma$change %in% c('UP', 'DOWN'), ]
  resMeR_list[[c]] <- list(res_DESeq2, res_limma)
  countMatrix_list[[c]] <- countMatrix
  resPeak_list[[c]] <- merge_peak
  name = paste0(c, ' vs WT')
  name = paste0(gsub(' ','_', name ))
  names(resMeR_list[[c]] ) <- paste0(name, '_', c('DESeq2', 'limma'))
}
names(resMeR_list) <- NULL
RMdem_list[[i]] <- as.list(unlist(resMeR_list, recursive = FALSE))
names(RMdem_list[[i]] ) <- paste0(cl, '_', names(RMdem_list[[i]]))
RMdep_list[[i]] <- resPeak_list 
name <- unique(gsub('(_limma|_DESeq2)', '', names(RMdem_list[[i]] )))
names(RMdep_list[[i]] ) <- name
contrast_IP_list[[i]] <- lapply(countMatrix_list, extract_names)
condition_IP_list[[i]] <- condition_list

# ############
# n = 31
# ############
# #GSE112795# without case RIP and Input 
# ############
# print(RMPdeg_datasets[n,]) #skip?
# i = RMPdeg_datasets[n,'GSE']
# cl = c('Hela')
# RMp = c('METTL3')
# RM = 'm6A'
# # condition = c(rep('METTL3KO', 2), rep('WT', 2))
# 
# RMP_allDat[grep('GSE112795', RMP_allDat$study_alias ) 
#            ,c("study_alias", "sample_alias", "sample_title")]
# 
# # View(allDat[[i]])
# 
# files <- list.files(file.path(data_directory,i))
# files
# getwd()
# setwd(data_directory)
# setwd(i)


# ############
# n = 32
# ############
# #GSE85724# without case RIP and Input 
# ############
# print(RMPdeg_datasets[n,]) #skip?
# i = RMPdeg_datasets[n,'GSE']
# cl = c('Jurkat', 'HEK293T')
# RMp = c('YTHDF1', 'YTHDF2', 'YTHDF3')
# RM = 'm6A'
# # condition = c(rep('METTL3OE', 2), rep('WT', 2))
# RMP_allDat[grep('GSE85724', RMP_allDat$study_alias ) 
#            ,c("study_alias", "sample_alias", "sample_title")]
# # View(allDat[[i]])
# 
# files <- list.files(file.path(data_directory,i))
# files
# getwd()
# setwd(data_directory)
# setwd(i)

############
n = 33
############
#GSE87190# with RIP and Input 
############
print(RMPdeg_datasets[n,]) #
i = RMPdeg_datasets[n,'GSE']
cl = c('NOMO-1', 'MA9.3RAS', 'MA9.3ITD')
RMp = c('FTO')
RM = 'm6A'
# Changes of m6A peaks on MYC transcripts in PBS- or R-2HG-treated sensitive cells (MA9.3ITD) upon FTO knockdown, or resistant cells (MA9.3RAS) upon FTO overexpression.
RMP_allDat[grep('GSE87190', RMP_allDat$study_alias ) 
           ,c("study_alias", "sample_alias", "sample_title")]
# View(allDat[[i]])
skip <- RMP_allDat[grepl('GSE87190', RMP_allDat$study_alias )  &  grepl('input', RMP_allDat$sample_title)
                   , "sample_alias"]
paste0(skip$sample_alias, collapse = " ")
skip <- RMP_allDat[grepl('GSE87190', RMP_allDat$study_alias ) &  grepl('m6A-seq', RMP_allDat$sample_title)
                   , "sample_alias"]
paste0(skip$sample_alias, collapse = " ")
ENA <- as.data.frame(RMP_allDat[grepl('GSE87190', RMP_allDat$study_alias ) & (!grepl('NOMO-1', RMP_allDat$sample_title))
                                ,c("study_alias", "sample_alias", "sample_title")])
ENA[order(ENA$sample_alias),]
ENA <- ENA[order(ENA$sample_alias),]

files <- list.files(file.path(data_directory,i))
files
getwd()
setwd(data_directory)
setwd(i)
extract_last_element <- function(x) {
  split_string <- strsplit(x, " ")
  last_element <- sapply(split_string, function(y) tail(y, n = 1))
  return(last_element)
}
# for (n in 1:nrow(allDat[[i]])) {
#   gsm = allDat[[i]]$geo_accession[n]
#   title = allDat[[i]]$title[n]
#   title <-gsub('m6A-seq', 'IP', extract_last_element(title))
#   
#   cmd = paste0("mv -i ", "metaPlotR/", gsm, "*.bam ", "metaPlotR/", gsm, "_", title, ".bam")
#   print(cmd)
#   system(cmd)
#   cmd = paste0("mv -i ", "metaPlotR/", gsm, "*.bam.bai ", "metaPlotR/", gsm, "_", title, ".bam.bai")
#   print(cmd)
#   system(cmd)
#   
# }
# nrows <- nrow(allDat[[i]])
# for (n in seq(1, nrows, 2)) {
#   gsm = allDat[[i]]$geo_accession[n]
#   gsm_IP = allDat[[i]]$geo_accession[n+1]
#   title = allDat[[i]]$title[n]
#   title <-gsub('m6A-seq', 'IP', extract_last_element(title))
#   
#   cmd = paste0("mv -i ", "metaPlotR/", gsm, "_input.bam ", "metaPlotR/", gsm_IP, "_input.bam")
#   print(cmd)
#   system(cmd)
#   cmd = paste0("mv -i ", "metaPlotR/", gsm, "_input.bam.bai ", "metaPlotR/", gsm_IP, "_input.bam.bai")
#   print(cmd)
#   system(cmd)
#   
# }
countfiles <- list.files('./metaPlotR/')
countfiles <- countfiles[grep('counts', countfiles)]
countfiles <- countfiles[! grepl('summary', countfiles)] 
countfiles
condition_list = list(MA9.3ITD_FTOKD = c(rep('FTOKD', 1), rep('WT', 1)),
                      MA9.3RAS_FTOOE = c(rep('FTOOE', 1), rep('WT', 1)))

###################################
#  Cell number-normalized (CNN)   #
#      MeRIP INPUT analysis       #
###################################
resInP_list <- vector("list", length(condition_list))
names(resInP_list) <- names(condition_list)
countGTF_list <- vector("list", length(condition_list))
names(countGTF_list) <- names(condition_list)
for(c in names(condition_list)){
  countGTF <- read.table( paste0("metaPlotR/Input_counts_", c, ".txt"),header=T)
  dim(countGTF)
  colnames(countGTF)
  rownames(countGTF) <- countGTF$Geneid
  countGTF <- countGTF[, !colnames(countGTF) %in% c('Geneid','Chr','Start','End','Strand','Length')]
  dim(countGTF) # 28265     4
  countGTF <- countGTF[, c(1,3)] # exclude the R-2HG (R-2-hydroxyglutarate)-treated samples
  condition <- condition_list[[c]]
  # if(c %in% c("HEK293_WTAPKD","HEK293_METTL3KD")){
  # 1. Differential gene analysis with edgeR
  # Without replicates, you cannot estimate which genes are differentially expressed using EdgeR or DESeq2. You can only calculate fold changes based on normalized read counts (preferentially CPM normalized by TMM method included in any EdgeR analysis) and apply a stringent fold-change cut-off to determine which genes are more or less expressed depending on the condition.
  y <- DGEList(countGTF)
  # Calculate fold change
  y <- calcNormFactors(y, method = c('TMM'))
  y$samples$lib.size <- colSums(y$counts)
  logcpm <- cpm(y, prior.count=2, log=TRUE)
  logcpm[1:5,]
  cpm <- cpm(y, prior.count=2, log=FALSE)
  cpm[1:5, ]
  fold_change <- cpm[,1]/cpm[,2]
  res <- countGTF
  res$logFC <- log2(fold_change)
  res$FoldChange <- fold_change
  res$change <- with(res, ifelse(abs(res$logFC) > log2(1.5),
                                 ifelse(res$logFC > log2(1.5),'UP','DOWN'),'NOT'))
  table(res$change)
  resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
  resInP$geneSymbol <- rownames(resInP)
  table(resInP$change)
  resInP_list[[c]] <- resInP
  countGTF_list[[c]] <- countGTF
}
name = paste0(names(resInP_list), ' vs WT')
RMdeg_list[[i]] <- resInP_list 
names(RMdeg_list[[i]] ) <- paste0(gsub(' ','_',(name)))
contrast_input_list[[i]] <- lapply(countGTF_list, extract_names)
condition_input_list[[i]] <- condition_list

###################################
#  Cell number-normalized (CNN)   #
#   MeRIP RIP analysis - Part1    #
###################################
resMeR_list <- vector("list", length(condition_list))
names(resMeR_list) <- names(condition_list)
countMatrix_list <- vector("list", length(condition_list))
names(countMatrix_list) <- names(condition_list)
resPeak_list <- vector("list",length(condition_list))
names(resPeak_list) <- names(condition_list)
for(c in names(condition_list)){
  merge_peak <- read.delim(paste0("metaPlotR/", c, "_newPeakFile.txt"),header=F)
  countMatrix <- read.table(paste0("metaPlotR/m6A_counts_", c, ".txt"),header=T)
  table(countMatrix$Geneid == merge_peak[-1,1])
  table(duplicated(merge_peak[-1,1]))
  dim(countMatrix)
  colnames(countMatrix)
  rownames(countMatrix) <- countMatrix$Geneid
  countMatrix <- countMatrix[, !colnames(countMatrix) %in% c('Geneid','Chr','Start','End','Strand','Length')]
  countMatrix <- countMatrix[, c(1,3)] # exclude the R-2HG (R-2-hydroxyglutarate)-treated samples
  # 1. 做好分组因子即可
  condition <- condition_list[[c]]
  groupList<-factor(condition)
  colData <- data.frame(row.names=colnames(countMatrix), groupList=groupList)
  # if(c %in% c("HEK293_WTAPKD","HEK293_METTL3KD")){
  # 2. Differential gene analysis with edgeR
  # Without replicates, you cannot estimate which genes are differentially expressed using EdgeR or DESeq2. You can only calculate fold changes based on normalized read counts (preferentially CPM normalized by TMM method included in any EdgeR analysis) and apply a stringent fold-change cut-off to determine which genes are more or less expressed depending on the condition.
  y <- DGEList(countMatrix)
  # Calculate fold change
  y <- calcNormFactors(y, method = c('TMM'))
  y$samples$lib.size <- colSums(y$counts)
  logcpm <- cpm(y, prior.count=2, log=TRUE)
  logcpm[1:5,]
  cpm <- cpm(y, prior.count=2, log=FALSE)
  cpm[1:5, ]
  # fold_change <- logcpm[,'Ratio1']/logcpm[,'Ratio2']
  fold_change <- cpm[,1]/cpm[,2]
  resMeR <- countMatrix
  resMeR$logFC <- log2(fold_change)
  resMeR$FoldChange <- fold_change
  resMeR$change <- with(resMeR, ifelse(abs(resMeR$logFC) > log2(1.5),
                                       ifelse(resMeR$logFC > log2(1.5),'UP','DOWN'),'NOT'))
  table(resMeR$change)
  resMeR = as.data.frame(resMeR[resMeR$change %in% c('UP', 'DOWN'), ])
  resMeR_list[[c]] <- list(resMeR)
  countMatrix_list[[c]] <- countMatrix
  resPeak_list[[c]] <- merge_peak
  name = paste0(c, ' vs WT')
  name = paste0(gsub(' ','_', name ))
  names(resMeR_list[[c]]) <- paste0(name, '_', c('limma'))
}
names(resMeR_list) <- NULL
RMdem_list[[i]] <- as.list(unlist(resMeR_list, recursive = FALSE))
names(RMdem_list[[i]] )
RMdep_list[[i]] <- resPeak_list 
name <- unique(gsub('(_limma|_DESeq2)', '', names(RMdem_list[[i]] )))
names(RMdep_list[[i]] ) <- name
contrast_IP_list[[i]] <- lapply(countMatrix_list, extract_names)
condition_IP_list[[i]] <- condition_list

############
n = 34
############
#GSE87515# with RIP and Input 
############
print(RMPdeg_datasets[n,]) # hg38
i = RMPdeg_datasets[n,'GSE']
cl = c('GSC11')#Glioblastoma
RMp = c('ALKBH5')
RM = 'm6A'
condition = c(rep('ALKBH5KD', 2), rep('WT', 2))
RMP_allDat[grep('GSE87515', RMP_allDat$study_alias ) 
           ,c("study_alias", "sample_alias", "sample_title")]
# View(allDat[[i]])
skip <- RMP_allDat[grepl('GSE87515', RMP_allDat$study_alias )  &  grepl('input', RMP_allDat$sample_title)
                   , "sample_alias"]
paste0(skip$sample_alias, collapse = " ")
skip <- RMP_allDat[grepl('GSE87515', RMP_allDat$study_alias ) &  grepl('IP', RMP_allDat$sample_title)
                   , "sample_alias"]
paste0(skip$sample_alias, collapse = " ")
ENA <- as.data.frame(RMP_allDat[grepl('GSE87515', RMP_allDat$study_alias )
                                ,c("study_alias", "sample_alias", "sample_title")])
ENA[order(ENA$sample_alias),]
ENA <- ENA[order(ENA$sample_alias),]
y <- ENA[ grepl('IP',ENA$sample_title),]
y <- ENA[ grepl('input',ENA$sample_title),]
y
x <- paste0('$outdir/', y$study_alias, '/', y$sample_alias, '_IP.bam')
x <- paste0('$outdir/', y$study_alias, '/', y$sample_alias, '_input.bam')
paste(x, collapse = ' ')

files <- list.files(file.path(data_directory,i))
files
getwd()
setwd(data_directory)
setwd(i)
c <- unique(condition[! condition %in% 'WT'])
###################################
#  Cell number-normalized (CNN)   #
#      MeRIP INPUT analysis       #
###################################
countGTF <- read.table( paste0("metaPlotR/Input_counts_", c, ".txt"),header=T)
rownames(countGTF) <- countGTF$Geneid
countGTF <- countGTF[, !colnames(countGTF) %in% c('Geneid','Chr','Start','End','Strand','Length')]
dim(countGTF) #28265     4
## Sample name
## Sample information
sample=data.frame(
  sample=colnames(countGTF),
  condition=condition)
## SummarizedExperiment input
dds<-DESeqDataSetFromMatrix(countData=countGTF,colData=sample,design=~condition)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 10
table(keep)
dds <- dds[keep,]
## Note on factor levels
dds$condition <- factor(dds$condition, levels = unique(condition))
dds$condition <- relevel(dds$condition, ref = "WT")
# Differential expression analysis-------------------------------------------
dds <- DESeq(dds)

res <-  results(dds, contrast=c("condition",unique(condition)))
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
table(res$change)
table(is.na(res$change))
# resOrdered <- as.data.frame(res[order(res$padj), ])
resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
resInP$geneSymbol <- rownames(resInP)
name = res@elementMetadata$description[2]
name = paste0(cl, '_',str_split(name, "condition ", simplify = T)[,2])
RMdeg_list[[i]] <- list( resInP)
names(RMdeg_list[[i]] ) <- paste0(gsub(' ','_',(name)))
c <- unique(condition[!condition %in% 'WT'])
countGTF_list <- list(colnames(countGTF)); names(countGTF_list) <- c
condition_list <- list(condition); names(condition_list) <- c
contrast_input_list[[i]] <- countGTF_list
condition_input_list[[i]] <- condition_list

###################################
#  Cell number-normalized (CNN)   #
#   MeRIP RIP analysis - Part1    #
###################################
## Sample information
merge_peak <- read.delim(paste0("metaPlotR/", c, "_newPeakFile.txt"),header=F)
countMatrix <- read.table(paste0("metaPlotR/m6A_counts_", c, ".txt"),header=T)
table(countMatrix$Geneid == merge_peak[-1,1])
table(duplicated(merge_peak$V1[-1]))
dim(countMatrix)
colnames(countMatrix)
rownames(countMatrix) <- countMatrix$Geneid
countMatrix <- countMatrix[, !colnames(countMatrix) %in% c('Geneid','Chr','Start','End','Strand','Length')]
# 1. 做好分组因子即可
groupList<-factor(condition)
colData <- data.frame(row.names=colnames(countMatrix), groupList=groupList)
# 2. Exploring Method 1: DESeq2-------------------------------------------
# 构建 DESeqDataSet 对象
dds <- DESeqDataSetFromMatrix(countData = countMatrix,
                              colData = colData,
                              design = ~ groupList)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 1
table(keep)
dds <- dds[keep,]
# 接着计算每个样本的 library size 用于后续的标准化。library size 等于落在该样本所有peaks上的reads的总数，即 counts 矩阵中每列的加和。如果bFullLibrarySize设为TRUE，则会使用library的总reads数（根据原BAM/BED文件统计获得）。然后使用estimateDispersions估计统计分布，需要将参数fitType设为’local’
dds<-estimateSizeFactors(dds)
dds<-estimateDispersions(dds,fitType="local")
dds <- nbinomWaldTest(dds)
res <- results(dds)
## 提取你想要的差异分析结果，我们这里是treated组对untreated组进行比较
res <-  results(dds, contrast=c("groupList", unique(condition)))
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
res_DESeq2 = as.data.frame(res)
table(res_DESeq2$change)#
res_DESeq2 = res_DESeq2[res_DESeq2$change %in% c('UP', 'DOWN'), ]
# 2. Exploring Method 2: limma-------------------------------------------
library(edgeR)
d0 <- DGEList(counts = countMatrix)
dim(d0)
# Filtering low count peaks
cutoff <- 1           
drop <- which(apply(cpm(d0), 1, max) < cutoff)
print(length(drop))
if (length(drop) >0) {
  d <- d0[-drop,] 
} else{ d <- d0}
dim(d) #number of peaks left: 8501      
d <- calcNormFactors(d, method = c('TMM'))
colnames(countMatrix)
genotype <- groupList
table(genotype)
mm <- model.matrix(~0 + genotype)#
v <- voom(d, mm, plot=F, normalize="quantile")#voom from limma:Transform count data to log2-counts per million (logtpm)
fit <- lmFit(v, mm)
head(coef(fit))
contr <- makeContrasts(paste0(colnames(coef(fit)), collapse  = ' - '), levels = colnames(coef(fit))) #Change if other comparison
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)
length(which(top.table$adj.P.Val < 0.05 & abs(top.table$logFC) > log2(1.5))) #how many are significant?
res_limma <- as.data.frame(top.table)
res_limma$change <- with(res_limma, ifelse(res_limma$adj.P.Val < 0.05 & abs(res_limma$logFC) > log2(1.5),
                                           ifelse(res_limma$logFC > log2(1.5),'UP','DOWN'),'NOT'))
table(res_limma$change)#
res_limma = res_limma[res_limma$change %in% c('UP', 'DOWN'), ]

RMdem_list[[i]] <- list(res_DESeq2, res_limma)
name = res@elementMetadata$description[2]
name = paste0(cl, '_',str_split(name, "groupList ", simplify = T)[,2])
name= paste0(gsub(' ','_',(name)))
names(RMdem_list[[i]] ) <- paste0(name, '_', c('DESeq2', 'limma'))
RMdep_list[[i]] <- list( merge_peak )
names(RMdep_list[[i]] ) <- name
countMatrix_list <- list(colnames(countMatrix)); names(countMatrix_list) <- c
condition_list <- list(condition); names(condition_list) <- c
contrast_IP_list[[i]] <- countMatrix_list
condition_IP_list[[i]] <- condition_list
############
n = 35
############
#GSE94808# with RIP and Input 
############
print(RMPdeg_datasets[n,]) # hg38
i = RMPdeg_datasets[n,'GSE']
cl = c('PBT003')#,'PBT707', 'PBT726') #Glioblastoma
RMp = c('METTL14','METTL3')
RM = 'm6A'
condition_list <- list(
  METTL3KD = c(rep('METTL3KD', 1), rep('WT', 1)),
  METTL14KD = c(rep('METTL14KD', 2), rep('WT', 1))
  
)
RMP_allDat[grep('GSE94808', RMP_allDat$study_alias ) 
           ,c("study_alias", "sample_alias", "sample_title")]
# View(allDat[[i]])
skip <- RMP_allDat[grepl('GSE94808', RMP_allDat$study_alias )  &  grepl('Input', RMP_allDat$sample_title)
                   , "sample_alias"]
paste0(skip$sample_alias, collapse = " ")
skip <- RMP_allDat[grepl('GSE94808', RMP_allDat$study_alias ) &  grepl('IP', RMP_allDat$sample_title)
                   , "sample_alias"]
paste0(skip$sample_alias, collapse = " ")

files <- list.files(file.path(data_directory,i))
files
getwd()
setwd(data_directory)
setwd(i)
ENA <- as.data.frame(RMP_allDat[grepl('GSE94808', RMP_allDat$study_alias )
                                ,c("study_alias", "sample_alias", "sample_title")])
ENA[order(ENA$sample_alias),]
ENA <- ENA[order(ENA$sample_alias),]
y <- ENA[ grepl('METTL14|control',ENA$sample_title) & grepl('IP',ENA$sample_title),]
y <- ENA[ grepl('METTL3|control',ENA$sample_title) & grepl('IP',ENA$sample_title),]
y <- ENA[ grepl('METTL14|control',ENA$sample_title) & grepl('Input',ENA$sample_title),]
y <- ENA[ grepl('METTL3|control',ENA$sample_title) & grepl('Input',ENA$sample_title),]
y
x <- paste0('$outdir/', y$study_alias, '/', y$sample_alias, '_IP.bam')
x <- paste0('$outdir/', y$study_alias, '/', y$sample_alias, '_input.bam')
paste(x, collapse = ' ')

countfiles <- list.files('./metaPlotR/')
countfiles <- countfiles[grep('counts', countfiles)]
countfiles <- countfiles[! grepl('summary', countfiles)] 
countfiles
###################################
#  Cell number-normalized (CNN)   #
#      MeRIP INPUT analysis       #
###################################
resInP_list <- vector("list", length(condition_list))
names(resInP_list) <- names(condition_list)
countGTF_list <- vector("list", length(condition_list))
names(countGTF_list) <- names(condition_list)
for(c in names(condition_list)){
  countGTF <- read.table( paste0("metaPlotR/Input_counts_", c, ".txt"),header=T)
  dim(countGTF)
  colnames(countGTF)
  rownames(countGTF) <- countGTF$Geneid
  countGTF <- countGTF[, !colnames(countGTF) %in% c('Geneid','Chr','Start','End','Strand','Length')]
  dim(countGTF) # 28265     4
  condition <- condition_list[[c]]
  if(c %in% c("METTL3KD")){
    # 1. Differential gene analysis with edgeR
    # Without replicates, you cannot estimate which genes are differentially expressed using EdgeR or DESeq2. You can only calculate fold changes based on normalized read counts (preferentially CPM normalized by TMM method included in any EdgeR analysis) and apply a stringent fold-change cut-off to determine which genes are more or less expressed depending on the condition.
    y <- DGEList(countGTF)
    # Calculate fold change
    y <- calcNormFactors(y, method = c('TMM'))
    y$samples$lib.size <- colSums(y$counts)
    logcpm <- cpm(y, prior.count=2, log=TRUE)
    logcpm[1:5,]
    cpm <- cpm(y, prior.count=2, log=FALSE)
    cpm[1:5, ]
    fold_change <- cpm[,1]/cpm[,2]
    res <- countGTF
    res$logFC <- log2(fold_change)
    res$FoldChange <- fold_change
    res$change <- with(res, ifelse(abs(res$logFC) > log2(1.5),
                                   ifelse(res$logFC > log2(1.5),'UP','DOWN'),'NOT'))
    table(res$change)
    resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
    resInP$geneSymbol <- rownames(resInP)
    table(resInP$change)
    resInP_list[[c]] <- resInP
    countGTF_list[[c]] <- countGTF
    
  } else{
    ## Sample information
    sample=data.frame(
      sample=colnames(countGTF),
      condition=condition)
    ## SummarizedExperiment input
    dds<-DESeqDataSetFromMatrix(countData=countGTF,colData=sample,design=~condition)
    ## Pre-filtering
    keep <- rowSums(counts(dds)) >= 10
    table(keep)
    dds <- dds[keep,]
    ## Note on factor levels
    dds$condition <- factor(dds$condition, levels = unique(condition))
    dds$condition <- relevel(dds$condition, ref = "WT")
    # Differential expression analysis-------------------------------------------
    dds <- DESeq(dds)
    res <-  results(dds, contrast=c("condition",unique(condition)))
    res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
    length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
    res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                                   ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
    table(res$change)
    table(is.na(res$change))
    # resOrdered <- as.data.frame(res[order(res$padj), ])
    resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
    resInP$geneSymbol <- rownames(resInP)
    table(resInP$change)
    resInP_list[[c]] <- resInP
    countGTF_list[[c]] <- countGTF
  }
  
}
name = paste0(cl, '_',names(resInP_list), ' vs WT')
RMdeg_list[[i]] <- resInP_list 
names(RMdeg_list[[i]] ) <- paste0(gsub(' ','_',(name)))
contrast_input_list[[i]] <- lapply(countGTF_list, extract_names)
condition_input_list[[i]] <- condition_list

###################################
#  Cell number-normalized (CNN)   #
#   MeRIP RIP analysis - Part1    #
###################################
resMeR_list <- vector("list", length(condition_list))
names(resMeR_list) <- names(condition_list)
countMatrix_list <- vector("list", length(condition_list))
names(countMatrix_list) <- names(condition_list)
resPeak_list <- vector("list",length(condition_list))
names(resPeak_list) <- names(condition_list)
for(c in names(condition_list)){
  merge_peak <- read.delim(paste0("metaPlotR/", c, "_newPeakFile.txt"),header=F)
  countMatrix <- read.table(paste0("metaPlotR/m6A_counts_", c, ".txt"),header=T)
  table(countMatrix$Geneid == merge_peak[-1,1])
  table(duplicated(merge_peak[-1,1]))
  dim(countMatrix)
  colnames(countMatrix)
  rownames(countMatrix) <- countMatrix$Geneid
  countMatrix <- countMatrix[, !colnames(countMatrix) %in% c('Geneid','Chr','Start','End','Strand','Length')]
  # 1. 做好分组因子即可
  condition <- condition_list[[c]]
  groupList<-factor(condition)
  colData <- data.frame(row.names=colnames(countMatrix), groupList=groupList)
  if(c %in% c("METTL3KD")){
    # 2. Differential gene analysis with edgeR
    # Without replicates, you cannot estimate which genes are differentially expressed using EdgeR or DESeq2. You can only calculate fold changes based on normalized read counts (preferentially CPM normalized by TMM method included in any EdgeR analysis) and apply a stringent fold-change cut-off to determine which genes are more or less expressed depending on the condition.
    y <- DGEList(countMatrix)
    # Calculate fold change
    y <- calcNormFactors(y, method = c('TMM'))
    y$samples$lib.size <- colSums(y$counts)
    logcpm <- cpm(y, prior.count=2, log=TRUE)
    logcpm[1:5,]
    cpm <- cpm(y, prior.count=2, log=FALSE)
    cpm[1:5, ]
    # fold_change <- logcpm[,'Ratio1']/logcpm[,'Ratio2']
    fold_change <- cpm[,1]/cpm[,2]
    resMeR <- countMatrix
    resMeR$logFC <- log2(fold_change)
    resMeR$FoldChange <- fold_change
    resMeR$change <- with(resMeR, ifelse(abs(resMeR$logFC) > log2(1.5),
                                         ifelse(resMeR$logFC > log2(1.5),'UP','DOWN'),'NOT'))
    table(resMeR$change)
    resMeR = as.data.frame(resMeR[resMeR$change %in% c('UP', 'DOWN'), ])
    resMeR_list[[c]] <- list(resMeR)
    countMatrix_list[[c]] <- countMatrix
    resPeak_list[[c]] <- merge_peak
    name = paste0(c, ' vs WT')
    name = paste0(gsub(' ','_', name ))
    names(resMeR_list[[c]]) <- paste0(name, '_', c('limma'))
    
  } else{
    # 2. Exploring Method 1: DESeq2-------------------------------------------
    # 构建 DESeqDataSet 对象
    dds <- DESeqDataSetFromMatrix(countData = countMatrix,
                                  colData = colData,
                                  design = ~ groupList)
    ## Pre-filtering
    keep <- rowSums(counts(dds)) >= 1
    table(keep)
    dds <- dds[keep,]
    # 接着计算每个样本的 library size 用于后续的标准化。library size 等于落在该样本所有peaks上的reads的总数，即 counts 矩阵中每列的加和。如果bFullLibrarySize设为TRUE，则会使用library的总reads数（根据原BAM/BED文件统计获得）。然后使用estimateDispersions估计统计分布，需要将参数fitType设为’local’
    dds<-estimateSizeFactors(dds)
    dds<-estimateDispersions(dds,fitType="local")
    dds <- nbinomWaldTest(dds)
    res <- results(dds)
    ## 提取你想要的差异分析结果，我们这里是treated组对untreated组进行比较
    res <-  results(dds, contrast=c("groupList", unique(condition)))
    res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
    length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
    res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                                   ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
    res_DESeq2 = as.data.frame(res)
    table(res_DESeq2$change)
    res_DESeq2 = res_DESeq2[res_DESeq2$change %in% c('UP', 'DOWN'), ]
    # 2. Exploring Method 2: limma-------------------------------------------
    library(edgeR)
    d0 <- DGEList(counts = countMatrix)
    dim(d0)
    # Filtering low count peaks
    cutoff <- 1           
    drop <- which(apply(cpm(d0), 1, max) < cutoff)
    print(length(drop))
    if (length(drop) >0) {
      d <- d0[-drop,] 
    } else{ d <- d0}
    dim(d) #number of peaks left: 8501      
    d <- calcNormFactors(d, method = c('TMM'))
    colnames(countMatrix)
    genotype <- groupList
    table(genotype)
    mm <- model.matrix(~0 + genotype)#
    v <- voom(d, mm, plot=F, normalize="quantile")#voom from limma:Transform count data to log2-counts per million (logtpm)
    fit <- lmFit(v, mm)
    head(coef(fit))
    contr <- makeContrasts(paste0(colnames(coef(fit)), collapse  = ' - '), levels = colnames(coef(fit))) #Change if other comparison
    tmp <- contrasts.fit(fit, contr)
    tmp <- eBayes(tmp)
    top.table <- topTable(tmp, sort.by = "P", n = Inf)
    head(top.table, 20)
    length(which(top.table$adj.P.Val < 0.05 & abs(top.table$logFC) > log2(1.5))) #how many are significant?
    res_limma <- as.data.frame(top.table)
    res_limma$change <- with(res_limma, ifelse(res_limma$adj.P.Val < 0.05 & abs(res_limma$logFC) > log2(1.5),
                                               ifelse(res_limma$logFC > log2(1.5),'UP','DOWN'),'NOT'))
    table(res_limma$change)
    res_limma = res_limma[res_limma$change %in% c('UP', 'DOWN'), ]
    resMeR_list[[c]] <- list(res_DESeq2, res_limma)
    countMatrix_list[[c]] <- countMatrix
    resPeak_list[[c]] <- merge_peak
    name = paste0(c, ' vs WT')
    name = paste0(gsub(' ','_', name ))
    names(resMeR_list[[c]] ) <- paste0(name, '_', c('DESeq2', 'limma'))
    
  }
}
names(resMeR_list) <- NULL
RMdem_list[[i]] <- as.list(unlist(resMeR_list, recursive = FALSE))
names(RMdem_list[[i]] ) <- paste0(cl, '_', names(RMdem_list[[i]] ))
RMdep_list[[i]] <- resPeak_list 
name <- unique(gsub('(_limma|_DESeq2)', '', names(RMdem_list[[i]] )))
names(RMdep_list[[i]] ) <- name
contrast_IP_list[[i]] <- lapply(countMatrix_list, extract_names)
condition_IP_list[[i]] <- condition_list

############
n = 36
############
#GSE76414# with RIP and Input 
############
print(RMPdeg_datasets[n,]) # hg38
i = RMPdeg_datasets[n,'GSE']
cl = c('MONO-MAC-6')#,'MA9.3ITD')
RMp = c('FTO')
RM = 'm6A'
condition = c(rep('FTOOE', 2), rep('WT', 2))#Sequencing of messenger RNAs with N6-methyladenosine modifications in acute myeloid leukemia (AML) with and without forced expression of FTO
# condition = c(rep('FTOKD', 1), rep('WT', 1)) # RNA-seq
RMP_allDat[grep('GSE76414', RMP_allDat$study_alias ) 
           ,c("study_alias", "sample_alias", "sample_title")]
# View(allDat[[i]])
skip <- RMP_allDat[grepl('GSE76414', RMP_allDat$study_alias )  &  grepl('input', RMP_allDat$sample_title)
                   , "sample_alias"]
paste0(skip$sample_alias, collapse = " ")
skip <- RMP_allDat[grepl('GSE76414', RMP_allDat$study_alias ) &  grepl('IP', RMP_allDat$sample_title)
                   , "sample_alias"]
paste0(skip$sample_alias, collapse = " ")
ENA <- as.data.frame(RMP_allDat[grepl('GSE76414', RMP_allDat$study_alias )
                                ,c("study_alias", "sample_alias", "sample_title")])
ENA[order(ENA$sample_alias),]
ENA <- ENA[order(ENA$sample_alias),]
files <- list.files(file.path(data_directory,i))
files
getwd()
setwd(data_directory)
setwd(i)
c <- unique(condition[! condition %in% 'WT'])
###################################
#  Cell number-normalized (CNN)   #
#      MeRIP INPUT analysis       #
###################################
countGTF <- read.table( paste0("metaPlotR/Input_counts_", c, ".txt"),header=T)
rownames(countGTF) <- countGTF$Geneid
countGTF <- countGTF[, !colnames(countGTF) %in% c('Geneid','Chr','Start','End','Strand','Length')]
dim(countGTF) #28265     4
## Sample name
## Sample information
sample=data.frame(
  sample=colnames(countGTF),
  condition=condition)
## SummarizedExperiment input
dds<-DESeqDataSetFromMatrix(countData=countGTF,colData=sample,design=~condition)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 10
table(keep)
dds <- dds[keep,]
## Note on factor levels
dds$condition <- factor(dds$condition, levels = unique(condition))
dds$condition <- relevel(dds$condition, ref = "WT")
# Differential expression analysis-------------------------------------------
dds <- DESeq(dds)

res <-  results(dds, contrast=c("condition",unique(condition)))
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
table(res$change)
table(is.na(res$change))
# resOrdered <- as.data.frame(res[order(res$padj), ])
resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
resInP$geneSymbol <- rownames(resInP)
name = res@elementMetadata$description[2]
name = paste0(cl, '_',str_split(name, "condition ", simplify = T)[,2])
RMdeg_list[[i]] <- list( resInP)
names(RMdeg_list[[i]] ) <- paste0(gsub(' ','_',(name)))
c <- unique(condition[!condition %in% 'WT'])
countGTF_list <- list(colnames(countGTF)); names(countGTF_list) <- c
condition_list <- list(condition); names(condition_list) <- c
contrast_input_list[[i]] <- countGTF_list
condition_input_list[[i]] <- condition_list

###################################
#  Cell number-normalized (CNN)   #
#   MeRIP RIP analysis - Part1    #
###################################
## Sample information
merge_peak <- read.delim(paste0("metaPlotR/", c, "_newPeakFile.txt"),header=F)
countMatrix <- read.table(paste0("metaPlotR/m6A_counts_", c, ".txt"),header=T)
table(countMatrix$Geneid == merge_peak[-1,1])
table(duplicated(merge_peak$V1[-1]))
dim(countMatrix)
colnames(countMatrix)
rownames(countMatrix) <- countMatrix$Geneid
countMatrix <- countMatrix[, !colnames(countMatrix) %in% c('Geneid','Chr','Start','End','Strand','Length')]
# 1. 做好分组因子即可
groupList<-factor(condition)
colData <- data.frame(row.names=colnames(countMatrix), groupList=groupList)
# 2. Exploring Method 1: DESeq2-------------------------------------------
# 构建 DESeqDataSet 对象
dds <- DESeqDataSetFromMatrix(countData = countMatrix,
                              colData = colData,
                              design = ~ groupList)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 1
table(keep)
dds <- dds[keep,]
# 接着计算每个样本的 library size 用于后续的标准化。library size 等于落在该样本所有peaks上的reads的总数，即 counts 矩阵中每列的加和。如果bFullLibrarySize设为TRUE，则会使用library的总reads数（根据原BAM/BED文件统计获得）。然后使用estimateDispersions估计统计分布，需要将参数fitType设为’local’
dds<-estimateSizeFactors(dds)
dds<-estimateDispersions(dds,fitType="local")
dds <- nbinomWaldTest(dds)
res <- results(dds)
## 提取你想要的差异分析结果，我们这里是treated组对untreated组进行比较
res <-  results(dds, contrast=c("groupList", unique(condition)))
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
res_DESeq2 = as.data.frame(res)
table(res_DESeq2$change)#
res_DESeq2 = res_DESeq2[res_DESeq2$change %in% c('UP', 'DOWN'), ]
# 2. Exploring Method 2: limma-------------------------------------------
library(edgeR)
d0 <- DGEList(counts = countMatrix)
dim(d0)
# Filtering low count peaks
cutoff <- 1           
drop <- which(apply(cpm(d0), 1, max) < cutoff)
print(length(drop))
if (length(drop) >0) {
  d <- d0[-drop,] 
} else{ d <- d0}
dim(d) #number of peaks left: 8501      
d <- calcNormFactors(d, method = c('TMM'))
colnames(countMatrix)
genotype <- groupList
table(genotype)
mm <- model.matrix(~0 + genotype)#
v <- voom(d, mm, plot=F, normalize="quantile")#voom from limma:Transform count data to log2-counts per million (logtpm)
fit <- lmFit(v, mm)
head(coef(fit))
contr <- makeContrasts(paste0(colnames(coef(fit)), collapse  = ' - '), levels = colnames(coef(fit))) #Change if other comparison
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)
length(which(top.table$adj.P.Val < 0.05 & abs(top.table$logFC) > log2(1.5))) #how many are significant?
res_limma <- as.data.frame(top.table)
res_limma$change <- with(res_limma, ifelse(res_limma$adj.P.Val < 0.05 & abs(res_limma$logFC) > log2(1.5),
                                           ifelse(res_limma$logFC > log2(1.5),'UP','DOWN'),'NOT'))
table(res_limma$change)#
res_limma = res_limma[res_limma$change %in% c('UP', 'DOWN'), ]

RMdem_list[[i]] <- list(res_DESeq2, res_limma)
name = res@elementMetadata$description[2]
name = paste0(cl, '_',str_split(name, "groupList ", simplify = T)[,2])
name= paste0(gsub(' ','_',(name)))
names(RMdem_list[[i]] ) <- paste0(name, '_', c('DESeq2', 'limma'))
RMdep_list[[i]] <- list( merge_peak )
names(RMdep_list[[i]] ) <- name
countMatrix_list <- list(colnames(countMatrix)); names(countMatrix_list) <- c
condition_list <- list(condition); names(condition_list) <- c
contrast_IP_list[[i]] <- countMatrix_list
condition_IP_list[[i]] <- condition_list

contrast_input_list[['GSE76414']]$FTOOE <- gsub('XC.', 'XC-', contrast_input_list[['GSE76414']]$FTOOE )
contrast_IP_list[['GSE76414']]$FTOOE <- gsub('XC.', 'XC-', contrast_IP_list[['GSE76414']]$FTOOE )

############
n = 37
############
#GSE76367# without case RIP, with case RNA 
############
print(RMPdeg_datasets[n,]) #
i = RMPdeg_datasets[n,'GSE']
cl = c('A549','H1299')
RMp = c('METTL3')
RM = 'm6A'
condition = c(rep('METTL3KD', 2), rep('WT', 2))
RMP_allDat[grep('GSE76367', RMP_allDat$study_alias ) 
           ,c("study_alias", "sample_alias", "sample_title")]
# View(allDat[[i]])
files <- list.files(file.path(data_directory,i))
files
getwd()
setwd(data_directory)
setwd(i)
###################################
#      RNA/Ribo-seq analysis      #
###################################
# 1) RNA-seq 
df_in <- list(A549 =  read.table(file.path(data_directory,i, 'GSE76367_A549_RNAseq.txt.gz'), header = T),
              H1299 =  read.table(file.path(data_directory,i, 'GSE76367_H1299_RNAseq.txt.gz'), header = T))
# a matrix table, which contained four columns including the gene's name, reads numbers in different samples
###################################
## Raw counts
#判断数据是否需要转换
ex <- df_in$H1299[,-1]
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
# Calculate LogC
LogC <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0) || (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)

# Log2 transform if conditions are met
if (LogC) {
  # Convert elements less than or equal to 0 to NaN
  ex[ex <= 0] <- NaN
  print("log2 transform unfinished")
} else {
  print("log2 transform not needed")
}
###################################
#   more than 1 cell line here   #
###################################
resInP_list <- list()
#####A549#####
##Gene count matrix
countGTF=ceiling(df_in$A549[,-1])
## Gene name
table(duplicated(df_in$A549$Gene))
dim(countGTF) #21316     3
## Sample name
colnames(countGTF)#<-gsub('riboseq_input_','',colnames(countGTF))
## Sample information
sample_meta=data.frame(
  sample=colnames(countGTF),
  condition=c(rep("WT",1),rep("METTL3KD", 2)))
## SummarizedExperiment input
condition = sample_meta$condition
dds<-DESeqDataSetFromMatrix(countData=countGTF, colData=sample_meta, design=~condition)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 10
table(keep)
dds <- dds[keep,]
## Note on factor levels
dds$condition <- factor(dds$condition, levels = c("METTL3KD","WT"))
dds$condition <- relevel(dds$condition, ref = "WT")
dds <- DESeq(dds)
res<-results(dds,name = "condition_METTL3KD_vs_WT")
# res <-  results(dds, contrast=c("condition",unique(condition))) # c("condition", rev(levels(genotype)))
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
# res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
#                                ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
res$change <- with(res, ifelse( abs(res$log2FoldChange) > log2(1.5),
                                ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
table(res$change)
res$geneSymbol <- df_in$A549$Gene[keep]
resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
table(duplicated(resInP$geneSymbol))
resInP_list[['A549_METTL3KD_vs_WT']] <- resInP

#####H1299#####
##Gene count matrix
countGTF=ceiling(df_in$H1299[,-1])
## Gene name
table(duplicated(df_in$A549$Gene))
dim(countGTF) #21316     3
## Sample name
colnames(countGTF)#<-gsub('riboseq_input_','',colnames(countGTF))
## Sample information
sample_meta=data.frame(
  sample=colnames(countGTF),
  condition=c(rep("WT",1),rep("METTL3KD", 2)))
## SummarizedExperiment input
condition = sample_meta$condition
dds<-DESeqDataSetFromMatrix(countData=countGTF, colData=sample_meta, design=~condition)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 10
table(keep)
dds <- dds[keep,]
## Note on factor levels
dds$condition <- factor(dds$condition, levels = c("METTL3KD","WT"))
dds$condition <- relevel(dds$condition, ref = "WT")
dds <- DESeq(dds)
res<-results(dds,name = "condition_METTL3KD_vs_WT")
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
# res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
#                                ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
# There's no fix because there's nothing to fix. I'm a little surprised you got adjusted p-values that low even. There are no meaningful p-values for your experiment because computing them would be impossible. No replicates means no meaningful p-values.
res$change <- with(res, ifelse( abs(res$log2FoldChange) > log2(1.5),
                                ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
table(res$change)
res$geneSymbol <- df_in$H1299$gene[keep]
resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
table(duplicated(resInP$geneSymbol))
resInP_list[['H1299_METTL3KD_vs_WT']] <- resInP
RMdeg_list[[i]] <- resInP_list
names(RMdeg_list[[i]] )
############
n = 38
############
#GSE107954# without case RIP, with case RNA 
############
print(RMPdeg_datasets[n,]) #
i = RMPdeg_datasets[n,'GSE']
cl = c('hUCB')#	human UCB CD34+ HSPCs
RMp = c('YTHDF2')
RM = 'm6A'
condition = c(rep('YTHDF2KD', 3), rep('WT', 3))
RMP_allDat[grep('GSE107954', RMP_allDat$study_alias ) 
           ,c("study_alias", "sample_alias", "sample_title")]
# View(allDat[[i]])
files <- list.files(file.path(data_directory,i))
files
getwd()
setwd(data_directory)
setwd(i)
###################################
#      RNA/Ribo-seq analysis      #
###################################
# 1) RNA-seq 
df_in <- read.table(file.path(data_directory,i, 'GSE107956_norm_counts.txt.gz'), header = T)
# a matrix table, which contained four columns including the gene's name, reads numbers in different samples
###################################
## Raw counts
#判断数据是否需要转换
ex <- df_in[,-1]
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
# Calculate LogC
LogC <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0) || (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)

# Log2 transform if conditions are met
if (LogC) {
  # Convert elements less than or equal to 0 to NaN
  ex[ex <= 0] <- NaN
  print("log2 transform unfinished")
} else {
  print("log2 transform not needed")
}
# [1] "log2 transform unfinished" for RAW
countGTF=as.data.frame(ceiling(df_in[, -1]))
id_mapped <- read.table("/data2/rluo4/hg38/gencode.v36.annotation.gene.probeMap",header = T)
dim(id_mapped) #  60675     6
countGTF$id <- str_split(df_in$ensid, '[.]', simplify = T)[,1]
id_mapped <- id_mapped[which(id_mapped$id %in% countGTF$id),]
dim(id_mapped) #57930     6
table(duplicated(id_mapped$gene))
countGTF <- left_join(id_mapped,countGTF,by="id")
countGTF <- aggregate(x = countGTF[,(ncol(id_mapped)+1):ncol(countGTF)],   #此时exprSet的第三列开始是表达矩阵内容
                      by = list(symbol = countGTF$gene),   #按照相同symbol分组，在组内计算
                      FUN = max) %>%   #原文中是计算最大值（max），也可以计算平均值（mean）或者中位数（median）
  column_to_rownames(var = 'symbol')
dim(countGTF) #56752     5
##Gene count matrix
## Gene name
table(duplicated(df_in$ensid))
dim(countGTF) #57992     6
## Sample name
colnames(countGTF)#<-gsub('riboseq_input_','',colnames(countGTF))
## Sample information
sample_meta=data.frame(
  sample=colnames(countGTF),
  condition=c(rep("WT",3),rep("YTHDF2KD", 3)))
## SummarizedExperiment input
condition = sample_meta$condition
dds<-DESeqDataSetFromMatrix(countData=countGTF, colData=sample_meta, design=~condition)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 10
table(keep)
dds <- dds[keep,]
## Note on factor levels
dds$condition <- factor(dds$condition, levels = c("YTHDF2KD","WT"))
dds$condition <- relevel(dds$condition, ref = "WT")
dds <- DESeq(dds)
res<-results(dds,name = "condition_YTHDF2KD_vs_WT")
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
table(res$change)
res$geneSymbol <- rownames(countGTF)[keep]
resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
table(duplicated(resInP$geneSymbol))
c <- unique(condition[!condition %in% 'WT'])
name = res@elementMetadata$description[2]
name = paste0(cl, '_',str_split(name, "condition ", simplify = T)[,2])
RMdeg_list[[i]] <- list( resInP)
names(RMdeg_list[[i]] ) <- paste0(gsub(' ','_',(name)))
condition_list <- list(condition); names(condition_list) <- c
contrast_input_list[[i]] <- NULL
condition_input_list[[i]] <- condition_list

############
n = 39
############
#GSE133132# with RIP and Input 
############
print(RMPdeg_datasets[n,]) #
i = RMPdeg_datasets[n,'GSE']
cl = c('BGC823')
RMp = c('METTL3')
RM = 'm6A'
condition = c(rep('METTL3KD', 1), rep('WT', 1))
RMP_allDat[grep('GSE133132', RMP_allDat$study_alias ) 
           ,c("study_alias", "sample_alias", "sample_title")]
# View(allDat[[i]])
skip <- RMP_allDat[grepl('GSE133132', RMP_allDat$study_alias )  &  grepl('Input', RMP_allDat$sample_title)
                   , "sample_alias"]
paste0(skip$sample_alias, collapse = " ")
skip <- RMP_allDat[grepl('GSE133132', RMP_allDat$study_alias ) &  grepl('IP', RMP_allDat$sample_title)
                   , "sample_alias"]
paste0(skip$sample_alias, collapse = " ")

files <- list.files(file.path(data_directory,i))
files
getwd()
setwd(data_directory)
setwd(i)
ENA <- as.data.frame(RMP_allDat[grepl('GSE133132', RMP_allDat$study_alias )
                                ,c("study_alias", "sample_alias", "sample_title")])
ENA[order(ENA$sample_alias),]
ENA <- ENA[order(ENA$sample_alias),]
setwd(data_directory)
setwd(i)
###################################
#  Cell number-normalized (CNN)   #
#      MeRIP INPUT analysis       #
###################################
c <- unique(condition[!condition %in% 'WT'])
countGTF <- read.table( paste0("metaPlotR/Input_counts_", c, ".txt"),header=T)
dim(countGTF)
colnames(countGTF)
rownames(countGTF) <- countGTF$Geneid
countGTF <- countGTF[, !colnames(countGTF) %in% c('Geneid','Chr','Start','End','Strand','Length')]
dim(countGTF) # 28265     4
# 2. Differential gene analysis with edgeR
# Without replicates, you cannot estimate which genes are differentially expressed using EdgeR or DESeq2. You can only calculate fold changes based on normalized read counts (preferentially CPM normalized by TMM method included in any EdgeR analysis) and apply a stringent fold-change cut-off to determine which genes are more or less expressed depending on the condition.
y <- DGEList(countGTF)
# Calculate fold change
y <- calcNormFactors(y, method = c('TMM'))
y$samples$lib.size <- colSums(y$counts)
logcpm <- cpm(y, prior.count=2, log=TRUE)
logcpm[1:5,]
cpm <- cpm(y, prior.count=2, log=FALSE)
cpm[1:5, ]
fold_change <- cpm[,1]/cpm[,2]
res <- countGTF
res$logFC <- log2(fold_change)
res$FoldChange <- fold_change
res$change <- with(res, ifelse(abs(res$logFC) > log2(1.5),
                               ifelse(res$logFC > log2(1.5),'UP','DOWN'),'NOT'))
table(res$change)
resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
resInP$geneSymbol <- rownames(resInP)
table(resInP$change)
table(resInP$FoldChange %in% c(0, Inf))
# resInP <- resInP[! resInP$FoldChange %in% c(0, Inf), ]
name = paste0(c, ' vs WT')
name = paste0(cl, '_', name)
RMdeg_list[[i]] <- list( resInP)
names(RMdeg_list[[i]] ) <- paste0(gsub(' ','_',(name)))

countGTF_list <- list(colnames(countGTF)); names(countGTF_list) <- c
condition_list <- list(condition); names(condition_list) <- c
contrast_input_list[[i]] <- countGTF_list
condition_input_list[[i]] <- condition_list

###################################
#  Cell number-normalized (CNN)   #
#   MeRIP RIP analysis - Part1    #
###################################
## Sample information
merge_peak <- read.delim(paste0("metaPlotR/", c, "_newPeakFile.txt"),header=F)
countMatrix <- read.table(paste0("metaPlotR/m6A_counts_", c, ".txt"),header=T)
table(countMatrix$Geneid == merge_peak[-1,1])
table(duplicated(merge_peak$V1[-1]))
dim(countMatrix)
colnames(countMatrix)
rownames(countMatrix) <- countMatrix$Geneid
countMatrix <- countMatrix[, !colnames(countMatrix) %in% c('Geneid','Chr','Start','End','Strand','Length')]
# 2. Differential gene analysis with edgeR
# Without replicates, you cannot estimate which genes are differentially expressed using EdgeR or DESeq2. You can only calculate fold changes based on normalized read counts (preferentially CPM normalized by TMM method included in any EdgeR analysis) and apply a stringent fold-change cut-off to determine which genes are more or less expressed depending on the condition.
y <- DGEList(countMatrix)
# Calculate fold change
y <- calcNormFactors(y, method = c('TMM'))
y$samples$lib.size <- colSums(y$counts)
logcpm <- cpm(y, prior.count=2, log=TRUE)
logcpm[1:5,]
cpm <- cpm(y, prior.count=2, log=FALSE)
cpm[1:5, ]
# fold_change <- logcpm[,'Ratio1']/logcpm[,'Ratio2']
fold_change <- cpm[,1]/cpm[,2]
resMeR <- countMatrix
resMeR$logFC <- log2(fold_change)
resMeR$FoldChange <- fold_change
resMeR$change <- with(resMeR, ifelse(abs(resMeR$logFC) > log2(1.5),
                                     ifelse(resMeR$logFC > log2(1.5),'UP','DOWN'),'NOT'))
table(resMeR$change)
res_limma <- resMeR[resMeR$change %in% c('UP','DOWN'),]
RMdem_list[[i]] <- list(res_limma)
name = paste0(c, ' vs WT')
name = paste0(gsub(' ','_', name ))
names(RMdem_list[[i]]) <- paste0(cl, '_', name, '_limma')
RMdep_list[[i]] <- list( merge_peak )
names(RMdep_list[[i]] ) <-  paste0(cl, '_', name)

countMatrix_list <- list(colnames(countMatrix)); names(countMatrix_list) <- c
condition_list <- list(condition); names(condition_list) <- c
contrast_IP_list[[i]] <- countMatrix_list
condition_IP_list[[i]] <- condition_list

############
n = 40
############
#GSE38957# with Aza-IP and Input
############
print(RMPdeg_datasets[n,]) #
i = RMPdeg_datasets[n,'GSE']
cl = c('Hela')
RMp = c('DNMT2','NSUN2')
RM = 'm5C'
# Proteini	TRNA aspartic acid methyltransferase 1	
# Gene namei	TRDMT1 (DNMT2, RNMT1)
condition_list = list(NSUN2PD = c(rep('NSUN2PD', 2), rep('WT', 1)),
                      DNMT2PD = c(rep('DNMT2PD', 1), rep('WT', 1)))

RMP_allDat[grep('GSE38957', RMP_allDat$study_alias ) 
           ,c("study_alias", "sample_alias", "sample_title")]
# View(allDat[[i]])
ENA <- as.data.frame(RMP_allDat[grepl('GSE38957', RMP_allDat$study_alias )
                                ,c("study_alias", "sample_alias", "sample_title")])
ENA[order(ENA$sample_alias),]

skip <- RMP_allDat[grepl('GSE38957', RMP_allDat$study_alias ) &  grepl('Input', RMP_allDat$sample_title)
                   , "sample_alias"]
paste0(skip$sample_alias, collapse = " ")
skip <- RMP_allDat[grepl('GSE38957', RMP_allDat$study_alias ) &  (!grepl('Input', RMP_allDat$sample_title))
                   , "sample_alias"]
paste0(skip$sample_alias, collapse = " ")

files <- list.files(file.path(data_directory,i))
files
getwd()
setwd(data_directory)
setwd(i)
for(c in names(condition_list)){
  if(c=='NSUN2PD'){
    DEP <- read_excel(file.path(data_directory,i, 'GSE38957_03-NSUN2_Aza-IP_All_Genes.xls'),sheet = 2)
    # DEP_sig <- read_excel(file.path(data_directory,i, 'GSE38957_03-NSUN2_Aza-IP_All_Genes.xls'),sheet = 3)
  } else{
    DEP <- read_excel(file.path(data_directory,i, 'GSE38957_02-DNMT2_Aza-IP_All_Genes.xls'),sheet = 1)
  }
  
  DEP <- DEP[!is.na(DEP$qValFDR), ]
  q_values <- 10^(-DEP$qValFDR / 10)
  DEP$qval <- q_values
  summary(q_values)
  table(DEP$qval < 0.01)
  DEP_FDR <- as.data.frame(DEP[DEP$qval < 0.01, ])
  
  DEP_FDR[,1] <- paste0(DEP_FDR[,1],'-peak-', rownames(DEP_FDR))
  #[-10Log10(q-values) - Window level bionomial p-values converted into q-value FDRs using John Storey's R package. (See the USeq package descriptions) ]
  # write.table(DEP_FDR[, c(1,2,4,5,3,6)], paste0("metaPlotR/", c, "_newPeakFile.txt"),sep="\t",col.names=F,row.names=F,quote=F)
  # DEP <- read_excel(file.path(data_directory,i, 'GSE38957_01-DNMT2_Aza-IP_tRNAs_and_signature_analysis.xls'),sheet = 2)
  
}

countfiles <- list.files('./metaPlotR/')
countfiles <- countfiles[grep('counts', countfiles)]
countfiles <- countfiles[! grepl('summary', countfiles)] 
countfiles
###################################
#         RIP-seq analysis        #
###################################
resMeR_list <- vector("list", length(condition_list))
names(resMeR_list) <- names(condition_list)
countMatrix_list <- vector("list", length(condition_list))
names(countMatrix_list) <- names(condition_list)
resPeak_list <- vector("list",length(condition_list))
names(resPeak_list) <- names(condition_list)

for(c in names(condition_list)){
  merge_peak <- read.delim(paste0("metaPlotR/", c, "_newPeakFile.txt"),header=F)
  countMatrix <- read.table(paste0("metaPlotR/m5C_counts_", c, ".txt"),header=T,sep = '\t')
  table(countMatrix$Geneid == merge_peak[-1,1])
  table(duplicated(merge_peak[-1,1]))
  dim(countMatrix)
  colnames(countMatrix)
  rownames(countMatrix) <- countMatrix$Geneid
  countMatrix <- countMatrix[, !colnames(countMatrix) %in% c('Geneid','Chr','Start','End','Strand','Length')]
  # 1. 做好分组因子即可
  condition <- condition_list[[c]]
  groupList<-factor(condition)
  colData <- data.frame(row.names=colnames(countMatrix), groupList=groupList)
  if(c %in% c("DNMT2PD")){
    # 2. Differential gene analysis with edgeR
    # Without replicates, you cannot estimate which genes are differentially expressed using EdgeR or DESeq2. You can only calculate fold changes based on normalized read counts (preferentially CPM normalized by TMM method included in any EdgeR analysis) and apply a stringent fold-change cut-off to determine which genes are more or less expressed depending on the condition.
    y <- DGEList(countMatrix)
    # Calculate fold change
    y <- calcNormFactors(y, method = c('TMM'))
    y$samples$lib.size <- colSums(y$counts)
    logcpm <- cpm(y, prior.count=2, log=TRUE)
    logcpm[1:5,]
    cpm <- cpm(y, prior.count=2, log=FALSE)
    cpm[1:5, ]
    # fold_change <- logcpm[,'Ratio1']/logcpm[,'Ratio2']
    fold_change <- cpm[,1]/cpm[,2]
    resMeR <- countMatrix
    resMeR$logFC <- log2(fold_change)
    resMeR$FoldChange <- fold_change
    resMeR$change <- with(resMeR, ifelse(abs(resMeR$logFC) > log2(1.5),
                                         ifelse(resMeR$logFC > log2(1.5),'UP','DOWN'),'NOT'))
    table(resMeR$change)
    resMeR = as.data.frame(resMeR[resMeR$change %in% c('UP', 'DOWN'), ])
    resMeR_list[[c]] <- list(resMeR)
    countMatrix_list[[c]] <- countMatrix
    resPeak_list[[c]] <- merge_peak
    name = paste0(c, ' vs WT')
    name = paste0(gsub(' ','_', name ))
    names(resMeR_list[[c]]) <- paste0(name, '_', c('limma'))
    
  } else{
    # 2. Exploring Method 1: DESeq2-------------------------------------------
    # 构建 DESeqDataSet 对象
    dds <- DESeqDataSetFromMatrix(countData = countMatrix,
                                  colData = colData,
                                  design = ~ groupList)
    ## Pre-filtering
    keep <- rowSums(counts(dds)) >= 1
    table(keep)
    dds <- dds[keep,]
    # 接着计算每个样本的 library size 用于后续的标准化。library size 等于落在该样本所有peaks上的reads的总数，即 counts 矩阵中每列的加和。如果bFullLibrarySize设为TRUE，则会使用library的总reads数（根据原BAM/BED文件统计获得）。然后使用estimateDispersions估计统计分布，需要将参数fitType设为’local’
    dds<-estimateSizeFactors(dds)
    dds<-estimateDispersions(dds,fitType="local")
    dds <- nbinomWaldTest(dds)
    res <- results(dds)
    ## 提取你想要的差异分析结果，我们这里是treated组对untreated组进行比较
    res <-  results(dds, contrast=c("groupList", unique(condition)))
    res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
    length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
    res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                                   ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
    res_DESeq2 = as.data.frame(res)
    table(res_DESeq2$change)
    res_DESeq2 = res_DESeq2[res_DESeq2$change %in% c('UP', 'DOWN'), ]
    # 2. Exploring Method 2: limma-------------------------------------------
    library(edgeR)
    d0 <- DGEList(counts = countMatrix)
    dim(d0)
    # Filtering low count peaks
    cutoff <- 1           
    drop <- which(apply(cpm(d0), 1, max) < cutoff)
    print(length(drop))
    if (length(drop) >0) {
      d <- d0[-drop,] 
    } else{ d <- d0}
    dim(d) #number of peaks left: 8501      
    d <- calcNormFactors(d, method = c('TMM'))
    colnames(countMatrix)
    genotype <- groupList
    table(genotype)
    mm <- model.matrix(~0 + genotype)#
    v <- voom(d, mm, plot=F, normalize="quantile")#voom from limma:Transform count data to log2-counts per million (logtpm)
    fit <- lmFit(v, mm)
    head(coef(fit))
    contr <- makeContrasts(paste0(colnames(coef(fit)), collapse  = ' - '), levels = colnames(coef(fit))) #Change if other comparison
    tmp <- contrasts.fit(fit, contr)
    tmp <- eBayes(tmp)
    top.table <- topTable(tmp, sort.by = "P", n = Inf)
    head(top.table, 20)
    length(which(top.table$adj.P.Val < 0.05 & abs(top.table$logFC) > log2(1.5))) #how many are significant?
    res_limma <- as.data.frame(top.table)
    res_limma$change <- with(res_limma, ifelse(res_limma$adj.P.Val < 0.05 & abs(res_limma$logFC) > log2(1.5),
                                               ifelse(res_limma$logFC > log2(1.5),'UP','DOWN'),'NOT'))
    table(res_limma$change)
    res_limma = res_limma[res_limma$change %in% c('UP', 'DOWN'), ]
    resMeR_list[[c]] <- list(res_DESeq2, res_limma)
    countMatrix_list[[c]] <- countMatrix
    resPeak_list[[c]] <- merge_peak
    name = paste0(c, ' vs WT')
    name = paste0(gsub(' ','_', name ))
    names(resMeR_list[[c]] ) <- paste0(name, '_', c('DESeq2', 'limma'))
    
  }
}
names(resMeR_list) <- NULL
RMdem_list[[i]] <- as.list(unlist(resMeR_list, recursive = FALSE))
names(RMdem_list[[i]] ) <- paste0(cl, '_',names(RMdem_list[[i]]) )
RMdep_list[[i]] <- resPeak_list 
name <- unique(gsub('(_limma|_DESeq2)', '', names(RMdem_list[[i]] )))
names(RMdep_list[[i]] ) <- name
contrast_IP_list[[i]] <- lapply(countMatrix_list, extract_names)
condition_IP_list[[i]] <- condition_list

############
n = 41
############
#GSE112181# with RIP and Input 
############
print(RMPdeg_datasets[n,]) # hg38
i = RMPdeg_datasets[n,'GSE']
cl = c('A549')#,'Caco-2')
RMp = c('METTL1')
RM = 'm5C'
condition = c(rep('METTL1KD', 3), rep('WT', 5))
RMP_allDat[grep('GSE112181', RMP_allDat$study_alias ) 
           ,c("study_alias", "sample_alias", "sample_title")]
RMP_allDat[grep('GSE120455', RMP_allDat$study_alias ) 
           ,c("study_alias", "sample_alias", "sample_title")]
skip = RMP_allDat[grepl('GSE112181|GSE120455', RMP_allDat$study_alias ) & (grepl('BoRed-Seq', RMP_allDat$sample_title))
                  ,c("study_alias", "sample_alias", "sample_title")]
paste0(skip$sample_alias, collapse = "|")
skip <- RMP_allDat[grepl('GSE112181|GSE120455', RMP_allDat$study_alias )  & (! grepl('BoRed-Seq', RMP_allDat$sample_title)) 
                   &  grepl('Input|METTL1 KD Rep', RMP_allDat$sample_title)
                   , "sample_alias"]
paste0(skip$sample_alias, collapse = " ")
skip <- RMP_allDat[grepl('GSE112181|GSE120455', RMP_allDat$study_alias )  & (! grepl('BoRed-Seq', RMP_allDat$sample_title)) 
                   &  grepl('IP Rep', RMP_allDat$sample_title)
                   , "sample_alias"]
paste0(skip$sample_alias, collapse = " ")
ENA <- as.data.frame(RMP_allDat[grepl('GSE112181|GSE120455', RMP_allDat$study_alias ) 
                                ,c("study_alias", "sample_alias", "sample_title")])
ENA <- ENA[order(ENA$sample_alias),]
ENA

# View(allDat[[i]])
# setwd(data_directory)
# gse <- getGEO('GSE112181',destdir=".", AnnotGPL=F, getGPL=F)
# pdata <- pData(gse[[2]])
files <- list.files(file.path(data_directory,i))
files
getwd()
setwd(data_directory)
setwd(i)
c <- unique(condition[!condition %in% 'WT'])
###################################
#  Cell number-normalized (CNN)   #
#      MeRIP INPUT analysis       #
###################################
countGTF <- read.table( paste0("metaPlotR/Input_counts_", c, ".txt"),header=T)
rownames(countGTF) <- countGTF$Geneid
countGTF <- countGTF[, !colnames(countGTF) %in% c('Geneid','Chr','Start','End','Strand','Length')]
dim(countGTF) #28265     4
## Sample name
## Sample information
sample=data.frame(
  sample=colnames(countGTF),
  condition=condition)
## SummarizedExperiment input
dds<-DESeqDataSetFromMatrix(countData=countGTF,colData=sample,design=~condition)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 10
table(keep)
dds <- dds[keep,]
## Note on factor levels
dds$condition <- factor(dds$condition, levels = unique(condition))
dds$condition <- relevel(dds$condition, ref = "WT")
# Differential expression analysis-------------------------------------------
dds <- DESeq(dds)

res <-  results(dds, contrast=c("condition",unique(condition)))
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
table(res$change)
table(is.na(res$change))
# resOrdered <- as.data.frame(res[order(res$padj), ])
resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
resInP$geneSymbol <- rownames(resInP)
name = res@elementMetadata$description[2]
name = paste0(cl, '_',str_split(name, "condition ", simplify = T)[,2])
RMdeg_list[[i]] <- list( resInP)
names(RMdeg_list[[i]] ) <- paste0(gsub(' ','_',(name)))
countGTF_list <- list(colnames(countGTF)); names(countGTF_list) <- c
condition_list <- list(condition); names(condition_list) <- c
contrast_input_list[[i]] <- countGTF_list
condition_input_list[[i]] <- condition_list

###################################
#  Cell number-normalized (CNN)   #
#   MeRIP RIP analysis - Part1    #
###################################
## Sample information
merge_peak <- read.delim(paste0("metaPlotR/", c, "_newPeakFile.txt"),header=F)
countMatrix <- read.table(paste0("metaPlotR/m7G_counts_", c, ".txt"),header=T)
table(countMatrix$Geneid == merge_peak[-1,1])
table(duplicated(merge_peak$V1[-1]))
dim(countMatrix)
colnames(countMatrix)
rownames(countMatrix) <- countMatrix$Geneid
countMatrix <- countMatrix[, !colnames(countMatrix) %in% c('Geneid','Chr','Start','End','Strand','Length')]
# 1. 做好分组因子即可
groupList<-factor(condition)
colData <- data.frame(row.names=colnames(countMatrix), groupList=groupList)
# 2. Exploring Method 1: DESeq2-------------------------------------------
# 构建 DESeqDataSet 对象
dds <- DESeqDataSetFromMatrix(countData = countMatrix,
                              colData = colData,
                              design = ~ groupList)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 1
table(keep)
dds <- dds[keep,]
# 接着计算每个样本的 library size 用于后续的标准化。library size 等于落在该样本所有peaks上的reads的总数，即 counts 矩阵中每列的加和。如果bFullLibrarySize设为TRUE，则会使用library的总reads数（根据原BAM/BED文件统计获得）。然后使用estimateDispersions估计统计分布，需要将参数fitType设为’local’
dds<-estimateSizeFactors(dds)
dds<-estimateDispersions(dds,fitType="local")
dds <- nbinomWaldTest(dds)
res <- results(dds)
## 提取你想要的差异分析结果，我们这里是treated组对untreated组进行比较
res <-  results(dds, contrast=c("groupList", unique(condition)))
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
res_DESeq2 = as.data.frame(res)
table(res_DESeq2$change)#Not only
res_DESeq2 = res_DESeq2[res_DESeq2$change %in% c('UP', 'DOWN'), ]
# 2. Exploring Method 2: limma-------------------------------------------
library(edgeR)
d0 <- DGEList(counts = countMatrix)
dim(d0)
# Filtering low count peaks
cutoff <- 1           
drop <- which(apply(cpm(d0), 1, max) < cutoff)
print(length(drop))
if (length(drop) >0) {
  d <- d0[-drop,] 
} else{ d <- d0}
dim(d) #number of peaks left: 8501      
d <- calcNormFactors(d, method = c('TMM'))
colnames(countMatrix)
genotype <- groupList
table(genotype)
mm <- model.matrix(~0 + genotype)#
v <- voom(d, mm, plot=F, normalize="quantile")#voom from limma:Transform count data to log2-counts per million (logtpm)
fit <- lmFit(v, mm)
head(coef(fit))
contr <- makeContrasts(paste0(colnames(coef(fit)), collapse  = ' - '), levels = colnames(coef(fit))) #Change if other comparison
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)
length(which(top.table$adj.P.Val < 0.05 & abs(top.table$logFC) > log2(1.5))) #how many are significant?
res_limma <- as.data.frame(top.table)
res_limma$change <- with(res_limma, ifelse(res_limma$adj.P.Val < 0.05 & abs(res_limma$logFC) > log2(1.5),
                                           ifelse(res_limma$logFC > log2(1.5),'UP','DOWN'),'NOT'))
table(res_limma$change)#NOT only
res_limma = res_limma[res_limma$change %in% c('UP', 'DOWN'), ]

RMdem_list[[i]] <- list(res_DESeq2, res_limma)
name = res@elementMetadata$description[2]
name = paste0(cl, '_',str_split(name, "groupList ", simplify = T)[,2])
name= paste0(gsub(' ','_',(name)))
names(RMdem_list[[i]] ) <- paste0(name, '_', c('DESeq2', 'limma'))
RMdep_list[[i]] <- list( merge_peak )
names(RMdep_list[[i]] ) <- name
countMatrix_list <- list(colnames(countMatrix)); names(countMatrix_list) <- c
condition_list <- list(condition); names(condition_list) <- c
contrast_IP_list[[i]] <- countMatrix_list
condition_IP_list[[i]] <- condition_list
############
n = 42
############
#GSE94613# with RIP and Input 
############
print(RMPdeg_datasets[n,])
i = RMPdeg_datasets[n,'GSE']
cl = 'MOLM13'
RMp = 'METTL3'
RM = 'm6A'
condition = c(rep('METTL3KD', 4), rep('WT', 2))
RMP_allDat[grep('GSE94613', RMP_allDat$study_alias ) 
           ,c("study_alias", "sample_alias", "sample_title")]
RMP_allDat[grepl('GSE94613', RMP_allDat$study_alias ) & grepl('IP KD|Input KD|IP Ctrl|Input Ctrl', RMP_allDat$sample_title)
           , c("study_alias", "sample_alias", "sample_title")]
skip <- RMP_allDat[grepl('GSE94613', RMP_allDat$study_alias ) & grepl('IP KD|Input KD|IP Ctrl|Input Ctrl', RMP_allDat$sample_title) & (! grepl('Input Ctrl Day', RMP_allDat$sample_title))
                   , "sample_alias"]
paste0(skip$sample_alias, collapse = "|")
skip <- RMP_allDat[grepl('GSE94613', RMP_allDat$study_alias ) & grepl('IP KD|IP Ctrl', RMP_allDat$sample_title) & (! grepl('Input Ctrl Day', RMP_allDat$sample_title))
                   , "sample_alias"]
paste0(skip$sample_alias, collapse = " ")
skip <- RMP_allDat[grepl('GSE94613', RMP_allDat$study_alias ) & grepl('Input KD|Input Ctrl', RMP_allDat$sample_title) & (! grepl('Input Ctrl Day', RMP_allDat$sample_title))
                   , "sample_alias"]
paste0(skip$sample_alias, collapse = " ")
# View(allDat[[i]])
files <- list.files(file.path(data_directory,i))
files
# m6A RNA IP was run in duplicate for two METTL3 shRNA inducible knockdown cell lines (KD1 and KD2) and a scramble shRNA control at day 8 following doxycyclin induction
getwd()
setwd(data_directory)
setwd(i)
ENA <- as.data.frame(RMP_allDat[grepl('GSE94613', RMP_allDat$study_alias ) & grepl('IP KD|Input KD|IP Ctrl|Input Ctrl', RMP_allDat$sample_title)
                                , c("study_alias", "sample_alias", "sample_title")])
ENA <- ENA[order(ENA$sample_alias),]
ENA
c <- unique(condition[!condition %in% 'WT'])
###################################
#  Cell number-normalized (CNN)   #
#      MeRIP INPUT analysis       #
###################################
countGTF <- read.table( paste0("metaPlotR/Input_counts_", c, ".txt"),header=T)
rownames(countGTF) <- countGTF$Geneid
countGTF <- countGTF[, !colnames(countGTF) %in% c('Geneid','Chr','Start','End','Strand','Length')]
dim(countGTF) #28265     4
## Sample name
## Sample information
sample=data.frame(
  sample=colnames(countGTF),
  condition=condition)
## SummarizedExperiment input
dds<-DESeqDataSetFromMatrix(countData=countGTF,colData=sample,design=~condition)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 10
table(keep)
dds <- dds[keep,]
## Note on factor levels
dds$condition <- factor(dds$condition, levels = unique(condition))
dds$condition <- relevel(dds$condition, ref = "WT")
# Differential expression analysis-------------------------------------------
dds <- DESeq(dds)

res <-  results(dds, contrast=c("condition",unique(condition)))
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
table(res$change)
table(is.na(res$change))
# resOrdered <- as.data.frame(res[order(res$padj), ])
resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
resInP$geneSymbol <- rownames(resInP)
name = res@elementMetadata$description[2]
name = paste0(cl, '_',str_split(name, "condition ", simplify = T)[,2])
RMdeg_list[[i]] <- list( resInP)
names(RMdeg_list[[i]] ) <- paste0(gsub(' ','_',(name)))
countGTF_list <- list(colnames(countGTF)); names(countGTF_list) <- c
condition_list <- list(condition); names(condition_list) <- c
contrast_input_list[[i]] <- countGTF_list
condition_input_list[[i]] <- condition_list

###################################
#  Cell number-normalized (CNN)   #
#   MeRIP RIP analysis - Part1    #
###################################
## Sample information
merge_peak <- read.delim(paste0("metaPlotR/", c, "_newPeakFile.txt"),header=F)
countMatrix <- read.table(paste0("metaPlotR/m6A_counts_", c, ".txt"),header=T)
table(countMatrix$Geneid == merge_peak[-1,1])
table(duplicated(merge_peak$V1[-1]))
dim(countMatrix)
colnames(countMatrix)
rownames(countMatrix) <- countMatrix$Geneid
countMatrix <- countMatrix[, !colnames(countMatrix) %in% c('Geneid','Chr','Start','End','Strand','Length')]
# 1. 做好分组因子即可
groupList<-factor(condition)
colData <- data.frame(row.names=colnames(countMatrix), groupList=groupList)
# 2. Exploring Method 1: DESeq2-------------------------------------------
# 构建 DESeqDataSet 对象
dds <- DESeqDataSetFromMatrix(countData = countMatrix,
                              colData = colData,
                              design = ~ groupList)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 1
table(keep)
dds <- dds[keep,]
# 接着计算每个样本的 library size 用于后续的标准化。library size 等于落在该样本所有peaks上的reads的总数，即 counts 矩阵中每列的加和。如果bFullLibrarySize设为TRUE，则会使用library的总reads数（根据原BAM/BED文件统计获得）。然后使用estimateDispersions估计统计分布，需要将参数fitType设为’local’
dds<-estimateSizeFactors(dds)
dds<-estimateDispersions(dds,fitType="local")
dds <- nbinomWaldTest(dds)
res <- results(dds)
## 提取你想要的差异分析结果，我们这里是treated组对untreated组进行比较
res <-  results(dds, contrast=c("groupList", unique(condition)))
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
res_DESeq2 = as.data.frame(res)
table(res_DESeq2$change)#Not only
res_DESeq2 = res_DESeq2[res_DESeq2$change %in% c('UP', 'DOWN'), ]
# 2. Exploring Method 2: limma-------------------------------------------
library(edgeR)
d0 <- DGEList(counts = countMatrix)
dim(d0)
# Filtering low count peaks
cutoff <- 1           
drop <- which(apply(cpm(d0), 1, max) < cutoff)
print(length(drop))
if (length(drop) >0) {
  d <- d0[-drop,] 
} else{ d <- d0}
dim(d) #number of peaks left: 8501      
d <- calcNormFactors(d, method = c('TMM'))
colnames(countMatrix)
genotype <- groupList
table(genotype)
mm <- model.matrix(~0 + genotype)#
v <- voom(d, mm, plot=F, normalize="quantile")#voom from limma:Transform count data to log2-counts per million (logtpm)
fit <- lmFit(v, mm)
head(coef(fit))
contr <- makeContrasts(paste0(colnames(coef(fit)), collapse  = ' - '), levels = colnames(coef(fit))) #Change if other comparison
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)
length(which(top.table$adj.P.Val < 0.05 & abs(top.table$logFC) > log2(1.5))) #how many are significant?
res_limma <- as.data.frame(top.table)
res_limma$change <- with(res_limma, ifelse(res_limma$adj.P.Val < 0.05 & abs(res_limma$logFC) > log2(1.5),
                                           ifelse(res_limma$logFC > log2(1.5),'UP','DOWN'),'NOT'))
table(res_limma$change)#NOT only
res_limma = res_limma[res_limma$change %in% c('UP', 'DOWN'), ]

RMdem_list[[i]] <- list(res_DESeq2, res_limma)
name = res@elementMetadata$description[2]
name = paste0(cl, '_',str_split(name, "groupList ", simplify = T)[,2])
name= paste0(gsub(' ','_',(name)))
names(RMdem_list[[i]] ) <- paste0(name, '_', c('DESeq2', 'limma'))
RMdep_list[[i]] <- list( merge_peak )
names(RMdep_list[[i]] ) <- name
countMatrix_list <- list(colnames(countMatrix)); names(countMatrix_list) <- c
condition_list <- list(condition); names(condition_list) <- c
contrast_IP_list[[i]] <- countMatrix_list
condition_IP_list[[i]] <- condition_list

############
n = 43
############
#GSE90684# with RIP and Input 
############
print(RMPdeg_datasets[n,]) 
i = RMPdeg_datasets[n,'GSE']
cl = c('HepG2')#,'HEK293T')
RMp = c('METTL3')
RM = 'm6A'
# condition = c(rep('IGF2BP1', 2), rep('WT', 2))
# condition = c(rep('IGF2BP2', 2), rep('WT', 2))
# condition = c(rep('IGF2BP3', 2), rep('WT', 2))
condition = c(rep('METTL14KD', 2), rep('WT', 2))
RMP_allDat[grep('GSE90684', RMP_allDat$study_alias ) 
           ,c("study_alias", "sample_alias", "sample_title")]
RMP_allDat[grep('GSE90642', RMP_allDat$study_alias ) 
           ,c("study_alias", "sample_alias", "sample_title")]
skip <- RMP_allDat[grepl('GSE90642', RMP_allDat$study_alias )  &  grepl('input', RMP_allDat$sample_title)
                   , "sample_alias"]
paste0(skip$sample_alias, collapse = " ")
skip <- RMP_allDat[grepl('GSE90642', RMP_allDat$study_alias ) &  grepl('IP', RMP_allDat$sample_title)
                   , "sample_alias"]
paste0(skip$sample_alias, collapse = " ")

# View(allDat[[i]])
files <- list.files(file.path(data_directory,i))
files
getwd()
setwd(data_directory)
setwd(i)
ENA <- as.data.frame(RMP_allDat[grepl('GSE90642', RMP_allDat$study_alias ) 
                                , c("study_alias", "sample_alias", "sample_title")])
ENA <- ENA[order(ENA$sample_alias),]
ENA
c <- unique(condition[!condition %in% 'WT'])
###################################
#  Cell number-normalized (CNN)   #
#      MeRIP INPUT analysis       #
###################################
countGTF <- read.table( paste0("metaPlotR/Input_counts_", c, ".txt"),header=T)
rownames(countGTF) <- countGTF$Geneid
countGTF <- countGTF[, !colnames(countGTF) %in% c('Geneid','Chr','Start','End','Strand','Length')]
dim(countGTF) #28265     4
## Sample name
## Sample information
sample=data.frame(
  sample=colnames(countGTF),
  condition=condition)
## SummarizedExperiment input
dds<-DESeqDataSetFromMatrix(countData=countGTF,colData=sample,design=~condition)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 10
table(keep)
dds <- dds[keep,]
## Note on factor levels
dds$condition <- factor(dds$condition, levels = unique(condition))
dds$condition <- relevel(dds$condition, ref = "WT")
# Differential expression analysis-------------------------------------------
dds <- DESeq(dds)

res <-  results(dds, contrast=c("condition",unique(condition)))
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
table(res$change)
table(is.na(res$change))
# resOrdered <- as.data.frame(res[order(res$padj), ])
resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
resInP$geneSymbol <- rownames(resInP)
name = res@elementMetadata$description[2]
name = paste0(cl, '_',str_split(name, "condition ", simplify = T)[,2])
RMdeg_list[[i]] <- list( resInP)
names(RMdeg_list[[i]] ) <- paste0(gsub(' ','_',(name)))
c <- unique(condition[!condition %in% 'WT'])
countGTF_list <- list(colnames(countGTF)); names(countGTF_list) <- c
condition_list <- list(condition); names(condition_list) <- c
contrast_input_list[[i]] <- countGTF_list
condition_input_list[[i]] <- condition_list

###################################
#  Cell number-normalized (CNN)   #
#   MeRIP RIP analysis - Part1    #
###################################
## Sample information
merge_peak <- read.delim(paste0("metaPlotR/", c, "_newPeakFile.txt"),header=F)
countMatrix <- read.table(paste0("metaPlotR/m6A_counts_", c, ".txt"),header=T)
table(countMatrix$Geneid == merge_peak[-1,1])
table(duplicated(merge_peak$V1[-1]))
dim(countMatrix)
colnames(countMatrix)
rownames(countMatrix) <- countMatrix$Geneid
countMatrix <- countMatrix[, !colnames(countMatrix) %in% c('Geneid','Chr','Start','End','Strand','Length')]
# 1. 做好分组因子即可
groupList<-factor(condition)
colData <- data.frame(row.names=colnames(countMatrix), groupList=groupList)
# 2. Exploring Method 1: DESeq2-------------------------------------------
# 构建 DESeqDataSet 对象
dds <- DESeqDataSetFromMatrix(countData = countMatrix,
                              colData = colData,
                              design = ~ groupList)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 1
table(keep)
dds <- dds[keep,]
# 接着计算每个样本的 library size 用于后续的标准化。library size 等于落在该样本所有peaks上的reads的总数，即 counts 矩阵中每列的加和。如果bFullLibrarySize设为TRUE，则会使用library的总reads数（根据原BAM/BED文件统计获得）。然后使用estimateDispersions估计统计分布，需要将参数fitType设为’local’
dds<-estimateSizeFactors(dds)
dds<-estimateDispersions(dds,fitType="local")
dds <- nbinomWaldTest(dds)
res <- results(dds)
## 提取你想要的差异分析结果，我们这里是treated组对untreated组进行比较
res <-  results(dds, contrast=c("groupList", unique(condition)))
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
res_DESeq2 = as.data.frame(res)
table(res_DESeq2$change)#Not only
res_DESeq2 = res_DESeq2[res_DESeq2$change %in% c('UP', 'DOWN'), ]
# 2. Exploring Method 2: limma-------------------------------------------
library(edgeR)
d0 <- DGEList(counts = countMatrix)
dim(d0)
# Filtering low count peaks
cutoff <- 1           
drop <- which(apply(cpm(d0), 1, max) < cutoff)
print(length(drop))
if (length(drop) >0) {
  d <- d0[-drop,] 
} else{ d <- d0}
dim(d) #number of peaks left: 8501      
d <- calcNormFactors(d, method = c('TMM'))
colnames(countMatrix)
genotype <- groupList
table(genotype)
mm <- model.matrix(~0 + genotype)#
v <- voom(d, mm, plot=F, normalize="quantile")#voom from limma:Transform count data to log2-counts per million (logtpm)
fit <- lmFit(v, mm)
head(coef(fit))
contr <- makeContrasts(paste0(colnames(coef(fit)), collapse  = ' - '), levels = colnames(coef(fit))) #Change if other comparison
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)
length(which(top.table$adj.P.Val < 0.05 & abs(top.table$logFC) > log2(1.5))) #how many are significant?
res_limma <- as.data.frame(top.table)
res_limma$change <- with(res_limma, ifelse(res_limma$adj.P.Val < 0.05 & abs(res_limma$logFC) > log2(1.5),
                                           ifelse(res_limma$logFC > log2(1.5),'UP','DOWN'),'NOT'))
table(res_limma$change)#NOT only
res_limma = res_limma[res_limma$change %in% c('UP', 'DOWN'), ]

RMdem_list[[i]] <- list(res_DESeq2, res_limma)
name = res@elementMetadata$description[2]
name = paste0(cl, '_',str_split(name, "groupList ", simplify = T)[,2])
name= paste0(gsub(' ','_',(name)))
names(RMdem_list[[i]] ) <- paste0(name, '_', c('DESeq2', 'limma'))
RMdep_list[[i]] <- list( merge_peak )
names(RMdep_list[[i]] ) <- name
countMatrix_list <- list(colnames(countMatrix)); names(countMatrix_list) <- c
condition_list <- list(condition); names(condition_list) <- c
contrast_IP_list[[i]] <- countMatrix_list
condition_IP_list[[i]] <- condition_list

############
n = 44
############
#GSE102336# with RIP and Input 
############
print(RMPdeg_datasets[n,]) # hg38
i = RMPdeg_datasets[n,'GSE']
cl = c('HepG2')#,'Hela')
RMp = c('ZCCHC4')
RM = 'm6A'
condition = c(rep('ZCCHC4KO', 2), rep('WT', 2))
RMP_allDat[grep('GSE102336', RMP_allDat$study_alias ) 
           ,c("study_alias", "sample_alias", "sample_title")]
skip = RMP_allDat[grepl('GSE102336', RMP_allDat$study_alias ) & (!grepl('polysome|PARCLIP', RMP_allDat$sample_title))
                  ,c("study_alias", "sample_alias", "sample_title")]
paste0(skip$sample_alias, collapse = "|")

skip <- RMP_allDat[grepl('GSE102336', RMP_allDat$study_alias )  & (!grepl('polysome|PARCLIP', RMP_allDat$sample_title))
                   &  grepl('Input', RMP_allDat$sample_title)
                   , "sample_alias"]
paste0(skip$sample_alias, collapse = " ")
skip <- RMP_allDat[grepl('GSE102336', RMP_allDat$study_alias ) & (!grepl('polysome|PARCLIP', RMP_allDat$sample_title))
                   &  grepl('IP', RMP_allDat$sample_title)
                   , "sample_alias"]
paste0(skip$sample_alias, collapse = " ")
ENA <- as.data.frame(RMP_allDat[grepl('GSE102336', RMP_allDat$study_alias ) & (!grepl('polysome|PARCLIP', RMP_allDat$sample_title))
                                , c("study_alias", "sample_alias", "sample_title")])
ENA <- ENA[order(ENA$sample_alias),]
ENA

# View(allDat[[i]])
files <- list.files(file.path(data_directory,i))
files
getwd()
setwd(data_directory)
setwd(i)
c <- unique(condition[!condition %in% 'WT'])
###################################
#  Cell number-normalized (CNN)   #
#      MeRIP INPUT analysis       #
###################################
countGTF <- read.table( paste0("metaPlotR/Input_counts_", c, ".txt"),header=T)
rownames(countGTF) <- countGTF$Geneid
countGTF <- countGTF[, !colnames(countGTF) %in% c('Geneid','Chr','Start','End','Strand','Length')]
dim(countGTF) #28265     4
## Sample name
## Sample information
sample=data.frame(
  sample=colnames(countGTF),
  condition=condition)
## SummarizedExperiment input
dds<-DESeqDataSetFromMatrix(countData=countGTF,colData=sample,design=~condition)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 10
table(keep)
dds <- dds[keep,]
## Note on factor levels
dds$condition <- factor(dds$condition, levels = unique(condition))
dds$condition <- relevel(dds$condition, ref = "WT")
# Differential expression analysis-------------------------------------------
dds <- DESeq(dds)

res <-  results(dds, contrast=c("condition",unique(condition)))
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
table(res$change)
table(is.na(res$change))
# resOrdered <- as.data.frame(res[order(res$padj), ])
resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
resInP$geneSymbol <- rownames(resInP)
name = res@elementMetadata$description[2]
name = paste0(cl, '_',str_split(name, "condition ", simplify = T)[,2])
RMdeg_list[[i]] <- list( resInP)
names(RMdeg_list[[i]] ) <- paste0(gsub(' ','_',(name)))
countGTF_list <- list(colnames(countGTF)); names(countGTF_list) <- c
condition_list <- list(condition); names(condition_list) <- c
contrast_input_list[[i]] <- countGTF_list
condition_input_list[[i]] <- condition_list

###################################
#  Cell number-normalized (CNN)   #
#   MeRIP RIP analysis - Part1    #
###################################
## Sample information
merge_peak <- read.delim(paste0("metaPlotR/", c, "_newPeakFile.txt"),header=F)
countMatrix <- read.table(paste0("metaPlotR/m6A_counts_", c, ".txt"),header=T)
table(countMatrix$Geneid == merge_peak[-1,1])
table(duplicated(merge_peak$V1[-1]))
dim(countMatrix)
colnames(countMatrix)
rownames(countMatrix) <- countMatrix$Geneid
countMatrix <- countMatrix[, !colnames(countMatrix) %in% c('Geneid','Chr','Start','End','Strand','Length')]
# 1. 做好分组因子即可
groupList<-factor(condition)
colData <- data.frame(row.names=colnames(countMatrix), groupList=groupList)
# 2. Exploring Method 1: DESeq2-------------------------------------------
# 构建 DESeqDataSet 对象
dds <- DESeqDataSetFromMatrix(countData = countMatrix,
                              colData = colData,
                              design = ~ groupList)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 1
table(keep)
dds <- dds[keep,]
# 接着计算每个样本的 library size 用于后续的标准化。library size 等于落在该样本所有peaks上的reads的总数，即 counts 矩阵中每列的加和。如果bFullLibrarySize设为TRUE，则会使用library的总reads数（根据原BAM/BED文件统计获得）。然后使用estimateDispersions估计统计分布，需要将参数fitType设为’local’
dds<-estimateSizeFactors(dds)
dds<-estimateDispersions(dds,fitType="local")
dds <- nbinomWaldTest(dds)
res <- results(dds)
## 提取你想要的差异分析结果，我们这里是treated组对untreated组进行比较
res <-  results(dds, contrast=c("groupList", unique(condition)))
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
res_DESeq2 = as.data.frame(res)
table(res_DESeq2$change)#Not only
res_DESeq2 = res_DESeq2[res_DESeq2$change %in% c('UP', 'DOWN'), ]
# 2. Exploring Method 2: limma-------------------------------------------
library(edgeR)
d0 <- DGEList(counts = countMatrix)
dim(d0)
# Filtering low count peaks
cutoff <- 1           
drop <- which(apply(cpm(d0), 1, max) < cutoff)
print(length(drop))
if (length(drop) >0) {
  d <- d0[-drop,] 
} else{ d <- d0}
dim(d) #number of peaks left: 8501      
d <- calcNormFactors(d, method = c('TMM'))
colnames(countMatrix)
genotype <- groupList
table(genotype)
mm <- model.matrix(~0 + genotype)#
v <- voom(d, mm, plot=F, normalize="quantile")#voom from limma:Transform count data to log2-counts per million (logtpm)
fit <- lmFit(v, mm)
head(coef(fit))
contr <- makeContrasts(paste0(colnames(coef(fit)), collapse  = ' - '), levels = colnames(coef(fit))) #Change if other comparison
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)
length(which(top.table$adj.P.Val < 0.05 & abs(top.table$logFC) > log2(1.5))) #how many are significant?
res_limma <- as.data.frame(top.table)
res_limma$change <- with(res_limma, ifelse(res_limma$adj.P.Val < 0.05 & abs(res_limma$logFC) > log2(1.5),
                                           ifelse(res_limma$logFC > log2(1.5),'UP','DOWN'),'NOT'))
table(res_limma$change)#NOT only
res_limma = res_limma[res_limma$change %in% c('UP', 'DOWN'), ]

RMdem_list[[i]] <- list(res_DESeq2, res_limma)
name = res@elementMetadata$description[2]
name = paste0(cl, '_',str_split(name, "groupList ", simplify = T)[,2])
name= paste0(gsub(' ','_',(name)))
names(RMdem_list[[i]] ) <- paste0(name, '_', c('DESeq2', 'limma'))
RMdep_list[[i]] <- list( merge_peak )
names(RMdep_list[[i]] ) <- name
countMatrix_list <- list(colnames(countMatrix)); names(countMatrix_list) <- c
condition_list <- list(condition); names(condition_list) <- c
contrast_IP_list[[i]] <- countMatrix_list
condition_IP_list[[i]] <- condition_list

############
n = 45
############
#GSE103497# with RIP and Input 
############
print(RMPdeg_datasets[n,]) # hg38
i = RMPdeg_datasets[n,'GSE']
cl = c('NB4', 'MONOMAC6')
cl = 'NB4'
RMp = c('FTO')
RM = 'm6A'
condition = c(rep('FTOKD', 2), rep('WT', 2)) #RNA-seq
condition = c(rep('FTOIB', 3), rep('WT', 3)) #RNA-seq:FB23-1;FB23-2
condition = c(rep('FTOIB', 1), rep('WT', 1)) #m6A-seq:FB23-2
condition_list = list(FTOKD =  c(rep('FTOKD', 2), rep('WT', 2)), 
                      FTO.1IB = c(rep('FTO.1IB', 3), rep('WT', 3)),
                      FTO.2IB = c(rep('FTO.2IB', 3), rep('WT', 3))
)
RMP_allDat[grep('GSE103496', RMP_allDat$study_alias ) 
           ,c("study_alias", "sample_alias", "sample_title")]
# View(allDat[[i]])
# setwd(data_directory)
# gse <- getGEO('GSE103497',destdir=".", AnnotGPL=F, getGPL=F)
# pdata <- pData(gse[[2]])
# RNA-seq and m6A-seq of AML cells with FTO knockdown or inhibition
ENA <- as.data.frame(RMP_allDat[grepl('GSE10349', RMP_allDat$study_alias ) 
                                , c("study_alias", "sample_alias", "sample_title")])
ENA <- ENA[order(ENA$sample_alias),]
ENA
skip = ENA[grepl('FB23-2|input|m6A-seq', ENA$sample_title)
           ,c("study_alias", "sample_alias", "sample_title")]
skip = ENA[grepl('FB23-2', ENA$sample_title)
           ,c("study_alias", "sample_alias", "sample_title")]
paste0(skip$sample_alias, collapse = "|")

skip <- ENA[ grepl('FB23-2|input', ENA$sample_title)  
             , c("study_alias", "sample_alias", "sample_title")]
paste0(skip$sample_alias, collapse = " ")
skip <- ENA[ grepl('m6A-seq', ENA$sample_title)  
             , c("study_alias", "sample_alias", "sample_title")]
paste0(skip$sample_alias, collapse = " ")

files <- list.files(file.path(data_directory,i))
files
getwd()
setwd(data_directory)
setwd(i)
###################################
#      RNA/Ribo-seq analysis      #
###################################
# 1) RNA-seq
count_files <- files[grepl('.genes.results.txt.gz', files)]
###################################
# Read and process each count file
count <- lapply(count_files, function(x) {
  # Read data using fread for faster reading
  data <- fread(file.path(data_directory, i, x))
  # Extract required columns
  df_new <- data[, c(1, 5)]# expected_count
  # Set column names same as file names
  colnames(df_new) <- c('geneSymbol',gsub('.genes.results.txt.gz', '', x))
  return(df_new)
})
# Merge dataframes in countGTF list by 'gene_id'
df_in <- Reduce(function(x, y) merge(x, y, by = 'geneSymbol', all = TRUE), count)
colSums(df_in[,-1])
table(colSums(df_in[,-1])=="1e+06")
###################################
#   more than 1 condition here   #
###################################
resInP_list <- list()
##Gene count matrix
countGTF=as.matrix(ceiling(df_in[, -1]))
## Gene name
table(duplicated(df_in$geneSymbol))
dim(countGTF) #25343    12
## Sample name
colnames(countGTF)
## Sample information
sample_meta=data.frame(
  sample=colnames(countGTF),
  condition=c(rep("WT",2),rep("FTOKD",2), rep("WT",3),rep("FTO.1IB",3), rep("FTO.2IB",3)))
#####condition_FTOKD_vs_WT#####
count <- as.matrix(countGTF[, !grepl('_WT|_23', colnames(countGTF))])
table(is.na(count))
sample <- sample_meta[sample_meta$sample %in% colnames(count), ]
## SummarizedExperiment input
condition = c(rep('FTOKD', 2), rep('WT', 2))
dds<-DESeqDataSetFromMatrix(countData=count, colData=sample,design=~condition)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 10
table(keep)
dds <- dds[keep,]
## Note on factor levels
dds$condition <- factor(dds$condition, levels = c("FTOKD","WT"))
dds$condition <- relevel(dds$condition, ref = "WT")
dds <- DESeq(dds)
res<-results(dds,name = "condition_FTOKD_vs_WT")
res <-  results(dds, contrast=c("condition",unique(condition))) # c("condition", rev(levels(genotype)))
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
table(res$change)
res$geneSymbol <-df_in$geneSymbol[keep]
resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
table(duplicated(resInP$geneSymbol))
resInP_list[['FTOKD_vs_WT']] <- resInP
#####condition_FTOIB.1_vs_WT#####
count <- as.matrix(countGTF[, !grepl('_sh|_232', colnames(countGTF))])
table(is.na(count))
sample <- sample_meta[sample_meta$sample %in% colnames(count), ]
## SummarizedExperiment input
condition = c(rep('FTO.1IB', 3), rep('WT', 3))
dds<-DESeqDataSetFromMatrix(countData=count, colData=sample,design=~condition)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 10
table(keep)
dds <- dds[keep,]
## Note on factor levels
dds$condition <- factor(dds$condition, levels = c("FTO.1IB","WT"))
dds$condition <- relevel(dds$condition, ref = "WT")
dds <- DESeq(dds)
res<-results(dds,name = "condition_FTO.1IB_vs_WT")
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
table(res$change)
res$geneSymbol <-df_in$geneSymbol[keep]
resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
table(duplicated(resInP$geneSymbol))
resInP_list[['FTO.1IB_vs_WT']] <- resInP
#####condition_FTOIB.2_vs_WT#####
count <- as.matrix(countGTF[, !grepl('_sh|_231', colnames(countGTF))])
table(is.na(count))
sample <- sample_meta[sample_meta$sample %in% colnames(count), ]
## SummarizedExperiment input
condition = c(rep('FTO.2IB', 3), rep('WT', 3))
dds<-DESeqDataSetFromMatrix(countData=count, colData=sample,design=~condition)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 10
table(keep)
dds <- dds[keep,]
## Note on factor levels
dds$condition <- factor(dds$condition, levels = c("FTO.2IB","WT"))
dds$condition <- relevel(dds$condition, ref = "WT")
dds <- DESeq(dds)
res<-results(dds,name = "condition_FTO.2IB_vs_WT")
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
table(res$change)
res$geneSymbol <-df_in$geneSymbol[keep]
resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
table(duplicated(resInP$geneSymbol))
resInP_list[['FTO.2IB_vs_WT']] <- resInP

RMdeg_list[[i]] <- resInP_list
name <- paste0(cl[1],'_', names(resInP_list))
names(RMdeg_list[[i]] ) <- name
contrast_input_list[[i]] <- NULL#lapply(countGTF_list, extract_names)
condition_input_list[[i]] <- condition_list
###################################
#         RIP-seq analysis        #
###################################
## Sample information
condition = c(rep('FTO.2IB', 1), rep('WT', 1)) #m6A-seq:FB23-2
c <- unique(condition[!condition %in% 'WT'])
merge_peak <- read.delim("metaPlotR/FTOIB_newPeakFile.txt",header=F)
countMatrix <- read.table("metaPlotR/m6A_counts_FTOIB.txt",header=T)
table(countMatrix$Geneid == merge_peak$V1[-1])
table(duplicated(merge_peak$V1[-1]))
dim(countMatrix)
colnames(countMatrix)
rownames(countMatrix) <- countMatrix$Geneid
countMatrix <- countMatrix[, !colnames(countMatrix) %in% c('Geneid','Chr','Start','End','Strand','Length')]
# 1. 做好分组因子即可
groupList<-factor(unique(condition))
# 2. Differential gene analysis with edgeR
# Without replicates, you cannot estimate which genes are differentially expressed using EdgeR or DESeq2. You can only calculate fold changes based on normalized read counts (preferentially CPM normalized by TMM method included in any EdgeR analysis) and apply a stringent fold-change cut-off to determine which genes are more or less expressed depending on the condition.
y <- DGEList(countMatrix)
# Calculate fold change
y <- calcNormFactors(y, method = c('TMM'))
y$samples$lib.size <- colSums(y$counts)
logcpm <- cpm(y, prior.count=2, log=TRUE)
logcpm[1:5,]
cpm <- cpm(y, prior.count=2, log=FALSE)
cpm[1:5, ]
# fold_change <- logcpm[,'Ratio1']/logcpm[,'Ratio2']
fold_change <- cpm[,1]/cpm[,2]
resMeR <- countMatrix
resMeR$logFC <- log2(fold_change)
resMeR$FoldChange <- fold_change
resMeR$change <- with(resMeR, ifelse(abs(resMeR$logFC) > log2(1.5),
                                     ifelse(resMeR$logFC > log2(1.5),'UP','DOWN'),'NOT'))
table(resMeR$change)
res_limma <- resMeR[resMeR$change %in% c('UP','DOWN'),]
RMdem_list[[i]] <- list(res_limma)
name = paste0(c, ' vs WT')
name = paste0(gsub(' ','_', name ))
names(RMdem_list[[i]]) <- paste0(cl, '_', name, '_limma')
RMdep_list[[i]] <- list( merge_peak )
names(RMdep_list[[i]] ) <-  paste0(cl, '_', name)
countMatrix_list <- list(colnames(countMatrix)); names(countMatrix_list) <- c
condition_list <- list(condition); names(condition_list) <- c
contrast_IP_list[[i]] <- countMatrix_list
condition_IP_list[[i]] <- list(FTO.2IB=condition)
contrast_IP_list[['GSE103497']]$FTO.2IB <- gsub('X.data2.rluo4.EpiTrans.RMDatasets.Raw.Out.macs2.GSE103495.', '', contrast_IP_list[['GSE103497']]$FTO.2IB )
contrast_input_list[[i]]  <- contrast_IP_list[[i]] 
contrast_input_list[[i]]$FTO.2IB <- gsub('_IP.bam', '_input.bam', contrast_input_list[[i]]$FTO.2IB )
############
n = 46
############
#GSE122803# with RIP and Input 
############
print(RMPdeg_datasets[n,]) # hg38
i = RMPdeg_datasets[n,'GSE']
cl = c('MEL624')
RMp = c('PCIF1')
RM = 'm6Am'
condition_list = list(MEL624.1_PCIF1KO = c(rep('PCIF1KO', 3), rep('WT', 3)), 
                      MEL624.2_PCIF1KO = c(rep('PCIF1KO', 3), rep('WT', 3)))

RMP_allDat[grep('GSE12280', RMP_allDat$study_alias ) 
           ,c("study_alias", "sample_alias", "sample_title")]
# View(allDat[[i]])
ENA <- as.data.frame(RMP_allDat[grepl('GSE122801', RMP_allDat$study_alias ) 
                                , c("study_alias", "sample_alias", "sample_title")])
ENA <- ENA[order(ENA$sample_alias),]
ENA
# paste0(ENA$sample_alias, collapse = "|")

skip <- ENA[ grepl('input', ENA$sample_title)  
             , c("study_alias", "sample_alias", "sample_title")]
paste0(skip$sample_alias, collapse = " ")
skip <- ENA[ grepl('m6A', ENA$sample_title)  
             , c("study_alias", "sample_alias", "sample_title")]
paste0(skip$sample_alias, collapse = " ")

files <- list.files(file.path(data_directory,i))
files
getwd()
setwd(data_directory)
setwd(i)
y <- ENA[grepl('PCIF1 KO1|control',ENA$sample_title) & grepl('IP',ENA$sample_title),]
y <- ENA[grepl('PCIF1 KO2|control',ENA$sample_title) & grepl('IP',ENA$sample_title),]
y <- ENA[grepl('PCIF1 KO1|control',ENA$sample_title) & grepl('input',ENA$sample_title),]
y <- ENA[grepl('PCIF1 KO2|control',ENA$sample_title) & grepl('input',ENA$sample_title),]
y
x <- paste0('$outdir/', y$study_alias, '/', y$sample_alias, '_IP.bam')
x <- paste0('$outdir/', y$study_alias, '/', y$sample_alias, '_input.bam')
paste(x, collapse = ' ')
countfiles <- list.files('./metaPlotR/')
countfiles <- countfiles[grep('counts', countfiles)]
countfiles <- countfiles[! grepl('summary', countfiles)] 
countfiles

###################################
#  Cell number-normalized (CNN)   #
#      MeRIP INPUT analysis       #
###################################
resInP_list <- vector("list", length(condition_list))
names(resInP_list) <- names(condition_list)
countGTF_list <- vector("list", length(condition_list))
names(countGTF_list) <- names(condition_list)
for(c in names(condition_list)){
  countGTF <- read.table( paste0("metaPlotR/Input_counts_", c, ".txt"),header=T)
  dim(countGTF)
  colnames(countGTF)
  rownames(countGTF) <- countGTF$Geneid
  countGTF <- countGTF[, !colnames(countGTF) %in% c('Geneid','Chr','Start','End','Strand','Length')]
  dim(countGTF) # 28265     4
  condition <- condition_list[[c]]
  ## Sample information
  sample=data.frame(
    sample=colnames(countGTF),
    condition=condition)
  ## SummarizedExperiment input
  dds<-DESeqDataSetFromMatrix(countData=countGTF,colData=sample,design=~condition)
  ## Pre-filtering
  keep <- rowSums(counts(dds)) >= 10
  table(keep)
  dds <- dds[keep,]
  ## Note on factor levels
  dds$condition <- factor(dds$condition, levels = unique(condition))
  dds$condition <- relevel(dds$condition, ref = "WT")
  # Differential expression analysis-------------------------------------------
  dds <- DESeq(dds)
  res <-  results(dds, contrast=c("condition",unique(condition)))
  res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
  length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
  res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                                 ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
  table(res$change)
  table(is.na(res$change))
  # resOrdered <- as.data.frame(res[order(res$padj), ])
  resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
  resInP$geneSymbol <- rownames(resInP)
  table(resInP$change)
  resInP_list[[c]] <- resInP
  countGTF_list[[c]] <- countGTF 
}
name = paste0(names(resInP_list), ' vs WT')
RMdeg_list[[i]] <- resInP_list 
names(RMdeg_list[[i]] ) <- paste0(gsub(' ','_',(name)))
contrast_input_list[[i]] <- lapply(countGTF_list, extract_names)
condition_input_list[[i]] <- condition_list

###################################
#  Cell number-normalized (CNN)   #
#   MeRIP RIP analysis - Part1    #
###################################
resMeR_list <- vector("list", length(condition_list))
names(resMeR_list) <- names(condition_list)
countMatrix_list <- vector("list", length(condition_list))
names(countMatrix_list) <- names(condition_list)
resPeak_list <- vector("list",length(condition_list))
names(resPeak_list) <- names(condition_list)
for(c in names(condition_list)){
  merge_peak <- read.delim(paste0("metaPlotR/", c, "_newPeakFile.txt"),header=F)
  countMatrix <- read.table(paste0("metaPlotR/m6A_counts_", c, ".txt"),header=T)
  table(countMatrix$Geneid == merge_peak[-1,1])
  table(duplicated(merge_peak[-1,1]))
  dim(countMatrix)
  colnames(countMatrix)
  rownames(countMatrix) <- countMatrix$Geneid
  countMatrix <- countMatrix[, !colnames(countMatrix) %in% c('Geneid','Chr','Start','End','Strand','Length')]
  # 1. 做好分组因子即可
  condition <- condition_list[[c]]
  groupList<-factor(condition)
  colData <- data.frame(row.names=colnames(countMatrix), groupList=groupList)
  
  # 2. Exploring Method 1: DESeq2-------------------------------------------
  # 构建 DESeqDataSet 对象
  dds <- DESeqDataSetFromMatrix(countData = countMatrix,
                                colData = colData,
                                design = ~ groupList)
  ## Pre-filtering
  keep <- rowSums(counts(dds)) >= 1
  table(keep)
  dds <- dds[keep,]
  # 接着计算每个样本的 library size 用于后续的标准化。library size 等于落在该样本所有peaks上的reads的总数，即 counts 矩阵中每列的加和。如果bFullLibrarySize设为TRUE，则会使用library的总reads数（根据原BAM/BED文件统计获得）。然后使用estimateDispersions估计统计分布，需要将参数fitType设为’local’
  dds<-estimateSizeFactors(dds)
  dds<-estimateDispersions(dds,fitType="local")
  dds <- nbinomWaldTest(dds)
  res <- results(dds)
  ## 提取你想要的差异分析结果，我们这里是treated组对untreated组进行比较
  res <-  results(dds, contrast=c("groupList", unique(condition)))
  res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
  length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
  res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                                 ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
  res_DESeq2 = as.data.frame(res)
  table(res_DESeq2$change)
  res_DESeq2 = res_DESeq2[res_DESeq2$change %in% c('UP', 'DOWN'), ]
  # 2. Exploring Method 2: limma-------------------------------------------
  library(edgeR)
  d0 <- DGEList(counts = countMatrix)
  dim(d0)
  # Filtering low count peaks
  cutoff <- 1           
  drop <- which(apply(cpm(d0), 1, max) < cutoff)
  print(length(drop))
  if (length(drop) >0) {
    d <- d0[-drop,] 
  } else{ d <- d0}
  dim(d) #number of peaks left: 8501      
  d <- calcNormFactors(d, method = c('TMM'))
  colnames(countMatrix)
  genotype <- groupList
  table(genotype)
  mm <- model.matrix(~0 + genotype)#
  v <- voom(d, mm, plot=F, normalize="quantile")#voom from limma:Transform count data to log2-counts per million (logtpm)
  fit <- lmFit(v, mm)
  head(coef(fit))
  contr <- makeContrasts(paste0(colnames(coef(fit)), collapse  = ' - '), levels = colnames(coef(fit))) #Change if other comparison
  tmp <- contrasts.fit(fit, contr)
  tmp <- eBayes(tmp)
  top.table <- topTable(tmp, sort.by = "P", n = Inf)
  head(top.table, 20)
  length(which(top.table$adj.P.Val < 0.05 & abs(top.table$logFC) > log2(1.5))) #how many are significant?
  res_limma <- as.data.frame(top.table)
  res_limma$change <- with(res_limma, ifelse(res_limma$adj.P.Val < 0.05 & abs(res_limma$logFC) > log2(1.5),
                                             ifelse(res_limma$logFC > log2(1.5),'UP','DOWN'),'NOT'))
  table(res_limma$change)
  res_limma = res_limma[res_limma$change %in% c('UP', 'DOWN'), ]
  resMeR_list[[c]] <- list(res_DESeq2, res_limma)
  countMatrix_list[[c]] <- countMatrix
  resPeak_list[[c]] <- merge_peak
  name = paste0(c, ' vs WT')
  name = paste0(gsub(' ','_', name ))
  names(resMeR_list[[c]] ) <- paste0(name, '_', c('DESeq2', 'limma'))
}
names(resMeR_list) <- NULL
RMdem_list[[i]] <- as.list(unlist(resMeR_list, recursive = FALSE))
# names(RMdem_list[[i]] ) <- paste0(cl, '_', names(RMdem_list[[i]]))
RMdep_list[[i]] <- resPeak_list 
name <- unique(gsub('(_limma|_DESeq2)', '', names(RMdem_list[[i]] )))
names(RMdep_list[[i]] ) <- name
contrast_IP_list[[i]] <- lapply(countMatrix_list, extract_names)
condition_IP_list[[i]] <- condition_list

############
n = 47
############
#GSE120024# with RIP and Input
############
print(RMPdeg_datasets[n,]) # hg38
i = RMPdeg_datasets[n,'GSE']
cl = c('EndoC-bH1')#'T2D islets', human pancreatic beta cell line 
RMp = c('METTL3', 'METTL14')
RM = 'm6A'
condition_list <- list(METTL3KD = c(rep('METTL3KD', 3), rep('WT', 3)), 
                       METTL14KD = c(rep('METTL14KD', 3), rep('WT', 3)) 
)
RMP_allDat[grep('GSE120024', RMP_allDat$study_alias ) 
           ,c("study_alias", "sample_alias", "sample_title")]
RMP_allDat[grep('GSE132306', RMP_allDat$study_alias ) 
           ,c("study_alias", "sample_alias", "sample_title")]
skip <- RMP_allDat[grepl('GSE132306', RMP_allDat$study_alias )  &  grepl('input', RMP_allDat$sample_title)
                   , "sample_alias"]
paste0(skip$sample_alias, collapse = " ")
skip <- RMP_allDat[grepl('GSE132306', RMP_allDat$study_alias ) &  grepl('IP', RMP_allDat$sample_title)
                   , "sample_alias"]
paste0(skip$sample_alias, collapse = " ")
ENA <- as.data.frame(RMP_allDat[grepl('GSE132306', RMP_allDat$study_alias ) 
                                , c("study_alias", "sample_alias", "sample_title")])
ENA <- ENA[order(ENA$sample_alias),]
ENA
# paste0(ENA$sample_alias, collapse = "|")
# View(allDat[[i]])

files <- list.files(file.path(data_directory,i))
files
getwd()
setwd(data_directory)
setwd(i)
y <- ENA[grepl('METTL3|Control',ENA$sample_title) & grepl('IP',ENA$sample_title),]
y <- ENA[grepl('METTL14|Control',ENA$sample_title) & grepl('IP',ENA$sample_title),]
y <- ENA[grepl('METTL3|Control',ENA$sample_title) & grepl('input',ENA$sample_title),]
y <- ENA[grepl('METTL14|Control',ENA$sample_title) & grepl('input',ENA$sample_title),]
y
x <- paste0('$outdir/', y$study_alias, '/', y$sample_alias, '_IP.bam')
x <- paste0('$outdir/', y$study_alias, '/', y$sample_alias, '_input.bam')
paste(x, collapse = ' ')
countfiles <- list.files('./metaPlotR/')
countfiles <- countfiles[grep('counts', countfiles)]
countfiles <- countfiles[! grepl('summary', countfiles)] 
countfiles
###################################
#  Cell number-normalized (CNN)   #
#      MeRIP INPUT analysis       #
###################################
resInP_list <- vector("list", length(condition_list))
names(resInP_list) <- names(condition_list)
countGTF_list <- vector("list", length(condition_list))
names(countGTF_list) <- names(condition_list)
for(c in names(condition_list)){
  countGTF <- read.table( paste0("metaPlotR/Input_counts_", c, ".txt"),header=T)
  dim(countGTF)
  colnames(countGTF)
  rownames(countGTF) <- countGTF$Geneid
  countGTF <- countGTF[, !colnames(countGTF) %in% c('Geneid','Chr','Start','End','Strand','Length')]
  dim(countGTF) # 28265     4
  condition <- condition_list[[c]]
  ## Sample information
  sample=data.frame(
    sample=colnames(countGTF),
    condition=condition)
  ## SummarizedExperiment input
  dds<-DESeqDataSetFromMatrix(countData=countGTF,colData=sample,design=~condition)
  ## Pre-filtering
  keep <- rowSums(counts(dds)) >= 10
  table(keep)
  dds <- dds[keep,]
  ## Note on factor levels
  dds$condition <- factor(dds$condition, levels = unique(condition))
  dds$condition <- relevel(dds$condition, ref = "WT")
  # Differential expression analysis-------------------------------------------
  dds <- DESeq(dds)
  res <-  results(dds, contrast=c("condition",unique(condition)))
  res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
  length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
  res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                                 ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
  table(res$change)
  table(is.na(res$change))
  # resOrdered <- as.data.frame(res[order(res$padj), ])
  resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
  resInP$geneSymbol <- rownames(resInP)
  table(resInP$change)
  resInP_list[[c]] <- resInP
  countGTF_list[[c]] <- countGTF 
}

name = paste0(cl, '_',names(resInP_list), ' vs WT')
RMdeg_list[[i]] <- resInP_list 
names(RMdeg_list[[i]] ) <- paste0(gsub(' ','_',(name)))
contrast_input_list[[i]] <- lapply(countGTF_list, extract_names)
condition_input_list[[i]] <- condition_list

###################################
#  Cell number-normalized (CNN)   #
#   MeRIP RIP analysis - Part1    #
###################################
resMeR_list <- vector("list", length(condition_list))
names(resMeR_list) <- names(condition_list)
countMatrix_list <- vector("list", length(condition_list))
names(countMatrix_list) <- names(condition_list)
resPeak_list <- vector("list",length(condition_list))
names(resPeak_list) <- names(condition_list)
for(c in names(condition_list)){
  merge_peak <- read.delim(paste0("metaPlotR/", c, "_newPeakFile.txt"),header=F)
  countMatrix <- read.table(paste0("metaPlotR/m6A_counts_", c, ".txt"),header=T)
  table(countMatrix$Geneid == merge_peak[-1,1])
  table(duplicated(merge_peak[-1,1]))
  dim(countMatrix)
  colnames(countMatrix)
  rownames(countMatrix) <- countMatrix$Geneid
  countMatrix <- countMatrix[, !colnames(countMatrix) %in% c('Geneid','Chr','Start','End','Strand','Length')]
  # 1. 做好分组因子即可
  condition <- condition_list[[c]]
  groupList<-factor(condition)
  colData <- data.frame(row.names=colnames(countMatrix), groupList=groupList)
  
  # 2. Exploring Method 1: DESeq2-------------------------------------------
  # 构建 DESeqDataSet 对象
  dds <- DESeqDataSetFromMatrix(countData = countMatrix,
                                colData = colData,
                                design = ~ groupList)
  ## Pre-filtering
  keep <- rowSums(counts(dds)) >= 1
  table(keep)
  dds <- dds[keep,]
  # 接着计算每个样本的 library size 用于后续的标准化。library size 等于落在该样本所有peaks上的reads的总数，即 counts 矩阵中每列的加和。如果bFullLibrarySize设为TRUE，则会使用library的总reads数（根据原BAM/BED文件统计获得）。然后使用estimateDispersions估计统计分布，需要将参数fitType设为’local’
  dds<-estimateSizeFactors(dds)
  dds<-estimateDispersions(dds,fitType="local")
  dds <- nbinomWaldTest(dds)
  res <- results(dds)
  ## 提取你想要的差异分析结果，我们这里是treated组对untreated组进行比较
  res <-  results(dds, contrast=c("groupList", unique(condition)))
  res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
  length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
  res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                                 ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
  res_DESeq2 = as.data.frame(res)
  table(res_DESeq2$change)
  res_DESeq2 = res_DESeq2[res_DESeq2$change %in% c('UP', 'DOWN'), ]
  # 2. Exploring Method 2: limma-------------------------------------------
  library(edgeR)
  d0 <- DGEList(counts = countMatrix)
  dim(d0)
  # Filtering low count peaks
  cutoff <- 1           
  drop <- which(apply(cpm(d0), 1, max) < cutoff)
  print(length(drop))
  if (length(drop) >0) {
    d <- d0[-drop,] 
  } else{ d <- d0}
  dim(d) #number of peaks left: 8501      
  d <- calcNormFactors(d, method = c('TMM'))
  colnames(countMatrix)
  genotype <- groupList
  table(genotype)
  mm <- model.matrix(~0 + genotype)#
  v <- voom(d, mm, plot=F, normalize="quantile")#voom from limma:Transform count data to log2-counts per million (logtpm)
  fit <- lmFit(v, mm)
  head(coef(fit))
  contr <- makeContrasts(paste0(colnames(coef(fit)), collapse  = ' - '), levels = colnames(coef(fit))) #Change if other comparison
  tmp <- contrasts.fit(fit, contr)
  tmp <- eBayes(tmp)
  top.table <- topTable(tmp, sort.by = "P", n = Inf)
  head(top.table, 20)
  length(which(top.table$adj.P.Val < 0.05 & abs(top.table$logFC) > log2(1.5))) #how many are significant?
  res_limma <- as.data.frame(top.table)
  res_limma$change <- with(res_limma, ifelse(res_limma$adj.P.Val < 0.05 & abs(res_limma$logFC) > log2(1.5),
                                             ifelse(res_limma$logFC > log2(1.5),'UP','DOWN'),'NOT'))
  table(res_limma$change)
  res_limma = res_limma[res_limma$change %in% c('UP', 'DOWN'), ]
  resMeR_list[[c]] <- list(res_DESeq2, res_limma)
  countMatrix_list[[c]] <- countMatrix
  resPeak_list[[c]] <- merge_peak
  name = paste0(c, ' vs WT')
  name = paste0(gsub(' ','_', name ))
  names(resMeR_list[[c]] ) <- paste0(name, '_', c('DESeq2', 'limma'))
}
names(resMeR_list) <- NULL
RMdem_list[[i]] <- as.list(unlist(resMeR_list, recursive = FALSE))
names(RMdem_list[[i]] ) <- paste0(cl, '_', names(RMdem_list[[i]]))
RMdep_list[[i]] <- resPeak_list 
name <- unique(gsub('(_limma|_DESeq2)', '', names(RMdem_list[[i]] )))
names(RMdep_list[[i]] ) <- name
contrast_IP_list[[i]] <- lapply(countMatrix_list, extract_names)
condition_IP_list[[i]] <- condition_list

############
n = 48
############
#GSE132306# same with n = 47
############

############
n = 49
############
#GSE106122# without case RIP, with case RNA
############
print(RMPdeg_datasets[n,]) 
i = RMPdeg_datasets[n,'GSE']
cl = c('HEL')#Erythroleukemia
RMp = c('METTL3', 'METTL14', 'FTO')
RM = 'm6A'
condition_list = list( METTL3KO = c(rep('METTL3KO', 4), rep('WT', 3)), 
                       WTAPKO =  c(rep('WTAPKO', 4), rep('WT', 3)) # 	CRISPR KO of WTAP with sg1 replicate 1
)
RMP_allDat[grep('GSE106122', RMP_allDat$study_alias ) 
           ,c("study_alias", "sample_alias", "sample_title")]
# View(allDat[[i]])
# setwd(data_directory)
# gse <- getGEO('GSE112181',destdir=".", AnnotGPL=F, getGPL=F)
# pdata <- pData(gse[[2]])
files <- list.files(file.path(data_directory,i))
files
getwd()
setwd(data_directory)
setwd(i)
###################################
#         RNA-seq analysis        #
###################################
files <- files[grepl('cpm_data_NTC|cpm_Sample_METTL3|cpm_Sample_WTAP', files)]
df_in <- read.delim(file.path(data_directory, i, 'GSM3661970_cpm_Sample_WTAP_sg1_rep1.txt.gz'))
df_in <- read.delim(file.path(data_directory, i, 'GSM3661966_cpm_Sample_METTL3_sg1_rep1.txt.gz'))
# keep =  rowSums(CPM > 1) >= 2
# v   = log2(data[keep,] + 0.1)
# fit  = lmFit(v, design)
# cont.matrix = makeContrasts(contrasts = Disease-Control, levels = design)
# fit2 = contrasts.fit(fit, cont.matrix)
# fit2  = eBayes(fit2, trend = TRUE)
###################################
# Read and process each count file
count <- lapply(files, function(x) {
  data <- read.delim(file.path(data_directory, i, x))
  df_new <- data[, grepl('gene|VALUE|cpm_value', colnames(data))]
  # Set column names same as file names
  colnames(df_new) <- c('geneSymbol',gsub('.txt.gz', '', x))
  df_new <- df_new[df_new[,1] !='-', ]
  df_new[,1] <- str_split(df_new[,1], ',', simplify = T)[,1]
  return(df_new)
})
names(count) <- gsub('.txt.gz', '', files)
# Merge dataframes in countGTF list by 'gene_id'
df_in <- Reduce(function(x, y) merge(x, y, by = 'geneSymbol', all = TRUE), count)
dim(df_in) #
## Sample name
colnames(df_in)
colSums(df_in[,-1])
table(colSums(df_in[,-1])=="1e+06")

## CPM ######
#判断数据是否需要转换
ex <- df_in[,-1]
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
# Calculate LogC
LogC <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0) || (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)

# Log2 transform if conditions are met
if (LogC) {
  # Convert elements less than or equal to 0 to NaN
  ex[ex <= 0] <- NaN
  print("log2 transform unfinished")
} else {
  print("log2 transform not needed")
}
colnames(countGTF)
library(tidyverse)
# Function to convert CPM to raw counts
cpm_to_counts <- function(cpm_matrix) {
  # Get the sample (column) names
  sample_names <- colnames(cpm_matrix)[-1] # Assuming the first column is gene_id
  # Calculate the total read counts for each sample
  total_read_counts <- colSums(cpm_matrix[, -1]) / 1e6 # Sum of CPM values for each sample column and divide by 1e6
  # Initialize a matrix to store raw counts
  raw_counts_matrix <- cpm_matrix
  # Convert CPM to raw counts
  for (sample in sample_names) {
    raw_counts_matrix[[sample]] <- cpm_matrix[[sample]] * total_read_counts[sample]
  }
  # Round the values to the nearest integer since read counts are integers
  raw_counts_matrix <- raw_counts_matrix %>% mutate(across(all_of(sample_names), round))
  return(raw_counts_matrix)
}
# Apply the function
count <- df_in
count <- na.omit(count)
countReads <- cpm_to_counts(count)
# Save the raw counts matrix to a file# write_tsv(raw_counts_matrix, "raw_counts_matrix.tsv")
###################################
#   more than 1 condition here   #
###################################
resInP_list <- list()
##Gene count matrix
countGTF=countReads[,-1]
## Gene name
table(duplicated(countReads$geneSymbol))
dim(countGTF) #56752     5
## Sample name
colnames(countGTF)#<-gsub('riboseq_input_','',colnames(countGTF))
## Sample information
sample_meta=data.frame(
  sample=colnames(countGTF),
  condition=c(rep("WT",3),rep("METTL3KO",4), rep("WTAPKO",4)))
#####condition_METTL3KO_vs_WT#####
count <- countGTF[, !grepl('WTAP', colnames(countGTF))]
count <- na.omit(count)
sample <- sample_meta[sample_meta$sample %in% colnames(count), ]
## SummarizedExperiment input
condition = c(rep('WT', 3), rep('METTL3KO', 4))
dds<-DESeqDataSetFromMatrix(countData=count, colData=sample,design=~condition)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 10
table(keep)
dds <- dds[keep,]
## Note on factor levels
dds$condition <- factor(dds$condition, levels = c("METTL3KO","WT"))
dds$condition <- relevel(dds$condition, ref = "WT")
dds <- DESeq(dds)
res<-results(dds,name = "condition_METTL3KO_vs_WT")
res <-  results(dds, contrast=c("condition",unique(condition))) # c("condition", rev(levels(genotype)))
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
table(res$change)
res$geneSymbol <- countReads$geneSymbol[keep]
resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
table(duplicated(resInP$geneSymbol))
resInP_list[['METTL3KO_vs_WT']] <- resInP
#####condition_WTAPKO_vs_WT#####
count <- countGTF[, !grepl('METTL3', colnames(countGTF))]
sample <- sample_meta[sample_meta$sample %in% colnames(count), ]
## SummarizedExperiment input
condition = c(rep('WTAPKO', 4), rep('WT', 3))
dds<-DESeqDataSetFromMatrix(countData=count, colData=sample,design=~condition)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 10
table(keep)
dds <- dds[keep,]
## Note on factor levels
dds$condition <- factor(dds$condition, levels = c("WTAPKO","WT"))
dds$condition <- relevel(dds$condition, ref = "WT")
dds <- DESeq(dds)
res<-results(dds,name = "condition_WTAPKO_vs_WT")
res <-  results(dds, contrast=c("condition",unique(condition))) # c("condition", rev(levels(genotype)))
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
table(res$change)
res$geneSymbol <- countReads$geneSymbol[keep]
resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
table(duplicated(resInP$geneSymbol))
resInP_list[['WTAPKO_vs_WT']] <- resInP
RMdeg_list[[i]] <- resInP_list
name <- paste0(cl,'_', names(resInP_list))
names(RMdeg_list[[i]] ) <- name
contrast_input_list[[i]] <- NULL
condition_input_list[[i]] <- condition_list

############
n = 50
############
#GSE117299#  without case RIP, with case RNA
############
print(RMPdeg_datasets[n,]) 
i = RMPdeg_datasets[n,'GSE']
cl = c('Hela')#, 'A549', 'H1299')#
RMp = c('METTL3')#, 'YTHDF1')
RM = 'm6A'
condition = c(rep('METTL3KD', 2), rep('WT', 1)) 
# condition = c(rep('FTOKO', 2), rep('WT', 3)) 
RMP_allDat[grep('GSE117299', RMP_allDat$study_alias ) 
           ,c("study_alias", "sample_alias", "sample_title")]
# View(allDat[[i]])
RMP_allDat[grepl('GSE117299', RMP_allDat$study_alias ) & grepl('lungcancer', RMP_allDat$sample_title)
           ,c("study_alias", "sample_alias", "sample_title")]
skip = RMP_allDat[grepl('GSE117299', RMP_allDat$study_alias ) & grepl('lungcancer', RMP_allDat$sample_title)
                  ,c("study_alias", "sample_alias", "sample_title")]
paste0(skip$sample_alias, collapse = "|")

skip <- RMP_allDat[grepl('GSE117299', RMP_allDat$study_alias ) & grepl('lungcancer', RMP_allDat$sample_title)
                   &  grepl('Input', RMP_allDat$sample_title)
                   , "sample_alias"]
paste0(skip$sample_alias, collapse = " ")
skip <- RMP_allDat[grepl('GSE117299', RMP_allDat$study_alias ) & grepl('lungcancer', RMP_allDat$sample_title)
                   &  grepl('IP', RMP_allDat$sample_title)
                   , "sample_alias"]
paste0(skip$sample_alias, collapse = " ")

files <- list.files(file.path(data_directory,i))
files
getwd()
setwd(data_directory)
setwd(i)
ENA <- as.data.frame(RMP_allDat[grepl('GSE117299', RMP_allDat$study_alias ) 
                                , c("study_alias", "sample_alias", "sample_title")])
ENA <- ENA[order(ENA$sample_alias),]
ENA
peaks <- read_excel(file.path(data_directory,i,'GSE117299_m6A_gene_list.xlsx'), sheet = 1)
peaks <- read_excel(file.path(data_directory,i,'GSE117299_Polysome_profiling.xlsx'), sheet = 1)
peaks <- read_excel(file.path(data_directory,i,'GSE117299_RNA-seq_Polysome_profiling_20160428.xlsx'), sheet = 1)
###################################
#      RNA/Ribo-seq analysis      #
###################################
df_in <- read_excel(file.path(data_directory,i,'GSE117299_RNA-seq_Polysome_profiling_20160428.xlsx'), sheet = 1)
## Raw counts
# ## Gene-level variance-stabilized counts from DESeq2 for RNA-Seq samples are normalized expression values that have been transformed to stabilize the variance across samples. These counts are obtained through a variance stabilizing transformation (VST), which is a common step in RNA-Seq data analysis pipelines.
setwd(data_directory)
setwd(i)
##Gene count matrix
count=as.data.frame(df_in[-(1:2), c(6, 9, 12)])  %>% mutate(across(everything(), as.numeric))
## Sample name
df_in <- df_in[-(1:2), 1] %>% as.data.frame()
df_in <- cbind(df_in, count)
colnames(df_in)<- c('geneSymbol','ShGFP_Input','ShML3.2_Input', 'ShML3.3_Input')
## Gene name
table(duplicated(df_in$geneSymbol))
dim(df_in) #56752     5
# x <- aggregate(x = df_in,   #此时exprSet的第三列开始是表达矩阵内容
#                       by = list(geneSymbol = countGTF$geneSymbol),   #按照相同symbol分组，在组内计算
#                       FUN = max) %>%   #原文中是计算最大值（max），也可以计算平均值（mean）或者中位数（median）
#   column_to_rownames(var = 'geneSymbol')
dim(df_in) # 50460     4 for v36 # 45344     4 for v44
countGTF <- df_in[,-1]
# countGTF <- unique(df_in) %>% dplyr::select(geneSymbol) # 去重
# countGTF <- countGTF[!duplicated(df_in$geneSymbol), ]
## Sample information
sample=data.frame(
  sample=colnames(countGTF),
  condition=c(rep("WT",1),rep("METTL3KD",2)))
## SummarizedExperiment input
dds<-DESeqDataSetFromMatrix(countData=countGTF,colData=sample,design=~condition)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 10
table(keep)
dds <- dds[keep,]
## Note on factor levels
dds$condition <- factor(dds$condition, levels = c("METTL3KD","WT"))
dds$condition <- relevel(dds$condition, ref = "WT")
# Differential expression analysis-------------------------------------------
dds <- DESeq(dds)
# res <- results(dds,name = "condition_METTL3KD_vs_WT")
res <-  results(dds, contrast=c("condition",unique(condition))) 
# rev(levels(Group))
# res$label<-with(res,ifelse(res$baseMean>10&abs(res$log2FoldChange)>0.58&res$padj<0.05,"sig","nosig"))
# ress=as.data.frame(res)
# res1=subset(ress,ress$label=="sig")
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
table(res$change)
res$geneSymbol <- df_in$geneSymbol[keep]
resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
tmp=by(resInP,resInP$geneSymbol,function(x) rownames(x)[which.max(x$baseMean)])
probes = as.character(tmp)
length(probes) # 7785
resInP <- resInP[rownames(resInP) %in% probes,]
name = res@elementMetadata$description[2]
name = paste0(cl, '_',str_split(name, "condition ", simplify = T)[,2])
RMdeg_list[[i]] <- list( resInP)
names(RMdeg_list[[i]] ) <- paste0(gsub(' ','_',(name)))

############
n = 51
############
#GSE128443# without case RIP and RNA
############
print(RMPdeg_datasets[n,])# hg38 skip?
i = RMPdeg_datasets[n,'GSE']
cl = c('Hela')# m6A-seq data was used as WT for GSE162357
RMp = c('FTO','SFPQ') #only SFPQ, not FTO,
RM = 'm6A'
# condition = c(rep('FTOOE', 2), rep('WT', 3)) 
condition = c(rep('SFPQOE', 2), rep('WT', 3)) 
# major RNA binding protein-SFPQ as a direct interaction partner of FTO. Our study showed that FTO and SFPQ were located in close proximity throughout the transcriptome and overexpression of SFPQ led to the demethylation of adjacent N6-methyladenosine on RN
RMP_allDat[grep('GSE128443', RMP_allDat$study_alias ) 
           ,c("study_alias", "sample_alias", "sample_title")]
ENA <- as.data.frame(RMP_allDat[grepl('GSE128443', RMP_allDat$study_alias ) & grepl('vector', RMP_allDat$sample_title)
                                , c("study_alias", "sample_alias", "sample_title")])
ENA <- ENA[order(ENA$sample_alias),]
ENA# View(allDat[[i]])
paste0(ENA$sample_alias, collapse = "|")
files <- list.files(file.path(data_directory,i))
files
getwd()
setwd(data_directory)
setwd(i)

# ############
# n = 52
# ############
# #GSE99017# without case RIP and RNA
# ############
# print(RMPdeg_datasets[n,])# skip?
# i = RMPdeg_datasets[n,'GSE']
# cl = c('Hela')#
# RMp = c('METTL14') 
# RM = 'm6A'
# # condition = c(rep('METTL14KD', 2), rep('WT', 3)) 
# RMP_allDat[grep('GSE99017', RMP_allDat$study_alias ) 
#            ,c("study_alias", "sample_alias", "sample_title")]
# # View(allDat[[i]])
# files <- list.files(file.path(data_directory,i))
# files
# getwd()
# setwd(data_directory)
# setwd(i)

############
n = 53
############
#GSE129716# with RIP and RNA
############
print(RMPdeg_datasets[n,])# hg38
i = RMPdeg_datasets[n,'GSE']
cl = c('HCT116')#
RMp = c('YTHDF3') 
RM = 'm6A'
condition = c(rep('YTHDF3KD', 2), rep('WT', 2)) 
RMP_allDat[grep('GSE129716', RMP_allDat$study_alias ) 
           ,c("study_alias", "sample_alias", "sample_title")]
# View(allDat[[i]])
skip <- RMP_allDat[grepl('GSE129716', RMP_allDat$study_alias )  &  grepl('input', RMP_allDat$sample_title)
                   , "sample_alias"]
paste0(skip$sample_alias, collapse = " ")
skip <- RMP_allDat[grepl('GSE129716', RMP_allDat$study_alias )  &  grepl('IP', RMP_allDat$sample_title)
                   , "sample_alias"]
paste0(skip$sample_alias, collapse = " ")

files <- list.files(file.path(data_directory,i))
files
getwd()
setwd(data_directory)
setwd(i)
ENA <- as.data.frame(RMP_allDat[grepl('GSE129716', RMP_allDat$study_alias ) 
                                , c("study_alias", "sample_alias", "sample_title")])
ENA <- ENA[order(ENA$sample_alias),]
ENA
c <- unique(condition[!condition %in% 'WT'])
###################################
#  Cell number-normalized (CNN)   #
#      MeRIP INPUT analysis       #
###################################
countGTF <- read.table( paste0("metaPlotR/Input_counts_", c, ".txt"),header=T)
rownames(countGTF) <- countGTF$Geneid
countGTF <- countGTF[, !colnames(countGTF) %in% c('Geneid','Chr','Start','End','Strand','Length')]
dim(countGTF) #28265     4
## Sample name
## Sample information
sample=data.frame(
  sample=colnames(countGTF),
  condition=condition)
## SummarizedExperiment input
dds<-DESeqDataSetFromMatrix(countData=countGTF,colData=sample,design=~condition)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 10
table(keep)
dds <- dds[keep,]
## Note on factor levels
dds$condition <- factor(dds$condition, levels = unique(condition))
dds$condition <- relevel(dds$condition, ref = "WT")
# Differential expression analysis-------------------------------------------
dds <- DESeq(dds)

res <-  results(dds, contrast=c("condition",unique(condition)))
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
table(res$change)
table(is.na(res$change))
# resOrdered <- as.data.frame(res[order(res$padj), ])
resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
resInP$geneSymbol <- rownames(resInP)
name = res@elementMetadata$description[2]
name = paste0(cl, '_',str_split(name, "condition ", simplify = T)[,2])
RMdeg_list[[i]] <- list( resInP)
names(RMdeg_list[[i]] ) <- paste0(gsub(' ','_',(name)))
countGTF_list <- list(colnames(countGTF)); names(countGTF_list) <- c
condition_list <- list(condition); names(condition_list) <- c
contrast_input_list[[i]] <- countGTF_list
condition_input_list[[i]] <- condition_list

###################################
#  Cell number-normalized (CNN)   #
#   MeRIP RIP analysis - Part1    #
###################################
## Sample information
merge_peak <- read.delim(paste0("metaPlotR/", c, "_newPeakFile.txt"),header=F)
countMatrix <- read.table(paste0("metaPlotR/m6A_counts_", c, ".txt"),header=T)
table(countMatrix$Geneid == merge_peak[-1,1])
table(duplicated(merge_peak$V1[-1]))
dim(countMatrix)
colnames(countMatrix)
rownames(countMatrix) <- countMatrix$Geneid
countMatrix <- countMatrix[, !colnames(countMatrix) %in% c('Geneid','Chr','Start','End','Strand','Length')]
# 1. 做好分组因子即可
groupList<-factor(condition)
colData <- data.frame(row.names=colnames(countMatrix), groupList=groupList)
# 2. Exploring Method 1: DESeq2-------------------------------------------
# 构建 DESeqDataSet 对象
dds <- DESeqDataSetFromMatrix(countData = countMatrix,
                              colData = colData,
                              design = ~ groupList)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 1
table(keep)
dds <- dds[keep,]
# 接着计算每个样本的 library size 用于后续的标准化。library size 等于落在该样本所有peaks上的reads的总数，即 counts 矩阵中每列的加和。如果bFullLibrarySize设为TRUE，则会使用library的总reads数（根据原BAM/BED文件统计获得）。然后使用estimateDispersions估计统计分布，需要将参数fitType设为’local’
dds<-estimateSizeFactors(dds)
dds<-estimateDispersions(dds,fitType="local")
dds <- nbinomWaldTest(dds)
res <- results(dds)
## 提取你想要的差异分析结果，我们这里是treated组对untreated组进行比较
res <-  results(dds, contrast=c("groupList", unique(condition)))
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
res_DESeq2 = as.data.frame(res)
table(res_DESeq2$change)#Not only
res_DESeq2 = res_DESeq2[res_DESeq2$change %in% c('UP', 'DOWN'), ]
# 2. Exploring Method 2: limma-------------------------------------------
library(edgeR)
d0 <- DGEList(counts = countMatrix)
dim(d0)
# Filtering low count peaks
cutoff <- 1           
drop <- which(apply(cpm(d0), 1, max) < cutoff)
print(length(drop))
if (length(drop) >0) {
  d <- d0[-drop,] 
} else{ d <- d0}
dim(d) #number of peaks left: 8501      
d <- calcNormFactors(d, method = c('TMM'))
colnames(countMatrix)
genotype <- groupList
table(genotype)
mm <- model.matrix(~0 + genotype)#
v <- voom(d, mm, plot=F, normalize="quantile")#voom from limma:Transform count data to log2-counts per million (logtpm)
fit <- lmFit(v, mm)
head(coef(fit))
contr <- makeContrasts(paste0(colnames(coef(fit)), collapse  = ' - '), levels = colnames(coef(fit))) #Change if other comparison
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)
length(which(top.table$adj.P.Val < 0.05 & abs(top.table$logFC) > log2(1.5))) #how many are significant?
res_limma <- as.data.frame(top.table)
res_limma$change <- with(res_limma, ifelse(res_limma$adj.P.Val < 0.05 & abs(res_limma$logFC) > log2(1.5),
                                           ifelse(res_limma$logFC > log2(1.5),'UP','DOWN'),'NOT'))
table(res_limma$change)#NOT only
res_limma = res_limma[res_limma$change %in% c('UP', 'DOWN'), ]

RMdem_list[[i]] <- list(res_DESeq2, res_limma)
name = res@elementMetadata$description[2]
name = paste0(cl, '_',str_split(name, "groupList ", simplify = T)[,2])
name= paste0(gsub(' ','_',(name)))
names(RMdem_list[[i]] ) <- paste0(name, '_', c('DESeq2', 'limma'))
RMdep_list[[i]] <- list( merge_peak )
names(RMdep_list[[i]] ) <- name
countMatrix_list <- list(colnames(countMatrix)); names(countMatrix_list) <- c
condition_list <- list(condition); names(condition_list) <- c
contrast_IP_list[[i]] <- countMatrix_list
condition_IP_list[[i]] <- condition_list

# ############
# n = 54
# ############
# #GSE119963# without case RIP and RNA
# ############
# print(RMPdeg_datasets[n,])# skip?
# i = RMPdeg_datasets[n,'GSE']
# cl = c('PEO1')#
# RMp = c('FTO') 
# RM = 'm6A'
# # condition = c(rep('FTOKD', 2), rep('WT', 2)) 
# RMP_allDat[grep('GSE119963', RMP_allDat$study_alias ) 
#            ,c("study_alias", "sample_alias", "sample_title")]
# # View(allDat[[i]])
# files <- list.files(file.path(data_directory,i))
# files
# getwd()
# setwd(data_directory)
# setwd(i)
# ############
# n = 55
# ############
# #GSE112970# without case RIP and RNA
# ############
# print(RMPdeg_datasets[n,])# skip?
# i = RMPdeg_datasets[n,'GSE']
# cl = c('SV-HUC-1 and RWPE-1 and 16HBE')#
# RMp = c('METTL3', 'YTHDF1') 
# RM = 'm6A'
# # condition = c(rep('YTHDF1KD', 2), rep('WT', 2)) 
# RMP_allDat[grep('GSE112970', RMP_allDat$study_alias ) 
#            ,c("study_alias", "sample_alias", "sample_title")]
# # View(allDat[[i]])
# files <- list.files(file.path(data_directory,i))
# files
# getwd()
# setwd(data_directory)
# setwd(i)

############
n = 56
############
#GSE120860# without case RIP, with case RNA
############
print(RMPdeg_datasets[n,])#
i = RMPdeg_datasets[n,'GSE']
cl = c('SMMC7721')#Hepatocellular carcinoma
RMp = c('YTHDF1') 
RM = 'm6A'
condition_list <- list( SMMC7721_YTHDF1KD = c(rep('YTHDF1KD', 1), rep('WT', 1)) ,
                        SMMC7721.Hx_YTHDF1OE = c(rep('YTHDF1OE', 1), rep('WT', 1)) ,
                        SMMC7721.Nx_YTHDF1OE = c(rep('YTHDF1OE', 1), rep('WT', 1))
)
RMP_allDat[grep('GSE120860', RMP_allDat$study_alias ) 
           ,c("study_alias", "sample_alias", "sample_title")]
# View(allDat[[i]])
files <- list.files(file.path(data_directory,i))
files
getwd()
setwd(data_directory)
setwd(i)
###################################
#         RNA-seq analysis        #
###################################
# df_new <- read.table(paste0(data_directory,i, "/GSM3407085_SMMC7721-shYTHDF2.txt.gz"),comment="#", header = T)
# head(df_new)
# df_new <- df_new[df_new$gene_name !='-', ]
# df_new$gene_name <- str_split(df_new$gene_name, ',', simplify = T)[,1]
# df_new <- df_new[df_new$gene_name !='-', ]
# df_new <- aggregate(x = df_new,   #此时exprSet的第三列开始是表达矩阵内容
#                     by = list(geneSymbol = df_new$gene_name),   #按照相同symbol分组，在组内计算
#                     FUN = max) %>%   #原文中是计算最大值（max），也可以计算平均值（mean）或者中位数（median）
#   column_to_rownames(var = 'geneSymbol')
# df_in <- read.delim(paste0(data_directory,i, "/GSM3405798_SMMC7721-OE+Nx.txt.gz"),comment="#")#, header = T)
# head(df_in)
files <- files[grepl('.txt.gz', files)]
###################################
# Read and process each count file
count <- lapply(files, function(x) {
  data <- read.delim(file.path(data_directory, i, x))
  df_new <- data[, grepl('gene_name|fpkm|count', colnames(data))]
  # Set column names same as file names
  colnames(df_new) <- c('geneSymbol',gsub('.txt.gz', '', x))
  df_new <- df_new[df_new[,1] !='-', ]
  df_new[,1] <- str_split(df_new[,1], ',', simplify = T)[,1]
  df_new <- df_new[df_new[,1] !='-', ]
  df_new <- aggregate(x = df_new,   #此时exprSet的第三列开始是表达矩阵内容
                      by = list(geneSymbol = df_new$geneSymbol),   #按照相同symbol分组，在组内计算
                      FUN = max) %>%   #原文中是计算最大值（max），也可以计算平均值（mean）或者中位数（median）
    column_to_rownames(var = 'geneSymbol')
  
  return(df_new)
})
names(count) <- gsub('.txt.gz', '', files)
# Remove any NULL elements from the list
# count <- Filter(Negate(is.null), count)
# # Check if there are any non-empty data frames to merge
# if (length(count) > 0) {
#   countGTF <- Reduce(function(x, y) merge(x, y, by = 'geneSymbol', all = TRUE), count[5:6])
# } else {
#   message("No valid data frames to merge.")
# }
# fpkm_to_counts = function(bg=NULL, mat=NULL, tlengths=NULL, mean_rps=100e6, 
#                           threshold=0){
#   if(is.null(mat)){
#     tmeas = as.matrix(ballgown::texpr(bg, 'FPKM'))
#     tlengths = sapply(width(ballgown::structure(bg)$trans), sum)
#   }else{
#     tmeas = mat
#     stopifnot(!is.null(tlengths))
#   }
#   index1 = which(rowMeans(tmeas) >= threshold)
#   tlengths = tlengths[index1]
#   counts = tlengths*tmeas[index1,]/1000
#   counts = round(counts*mean_rps/1e6)
#   return(counts)    
# }
fpkm_to_counts <- function(mat, tlengths, mean_rps = 100e6, threshold = 0) {
  # Ensure that the inputs are provided
  stopifnot(!is.null(mat), !is.null(tlengths))
  # Ensure the matrix and transcript lengths are compatible
  stopifnot(nrow(mat) == length(tlengths))
  # Convert to matrix if not already
  tmeas <- as.matrix(mat)
  # Filter based on the threshold
  index1 <- which(rowMeans(tmeas) >= threshold)
  # Filter the lengths and the measurements
  tlengths <- tlengths[index1]
  tmeas <- tmeas[index1, , drop = FALSE]
  # Calculate the counts
  counts <- tlengths * tmeas / 1000
  counts <- round(counts * mean_rps / 1e6)
  return(counts)
}
# library(polyester)
# help(fpkm_to_counts)
## Sample name
colnames(countGTF)
###################################
#   more than 1 condition here   #
###################################
resInP_list <- list()
for (c in names(condition_list)) {
  # c <- 'SMMC7721_Hx_YTHDF1OE'
  condition <- condition_list[[c]]
  ##Gene count matrix
  if(c == 'SMMC7721_YTHDF1KD'){
    df_in <- Reduce(function(x, y) merge(x, y, by = 'geneSymbol', all = TRUE), count[5:6])
    # library(ballgown)
    # data(bg)
    # countmat = fpkm_to_counts(bg, mean_rps=400000)
    # countGTF <- fpkm_to_counts(df_in[, -1], tlengths = df_new$length)
    countGTF <- ceiling(df_in[, -1])
  } else{
    df_in <- Reduce(function(x, y) merge(x, y, by = 'geneSymbol', all = TRUE), count[1:4])
    index <- str_split(c, '[.]', simplify = T)[,2]
    index <- str_split(index, '[_]', simplify = T)[,1]
    countGTF <- df_in[, grepl(index, colnames(df_in))]
  }
  # 2. Differential gene analysis with edgeR
  # Without replicates, you cannot estimate which genes are differentially expressed using EdgeR or DESeq2. You can only calculate fold changes based on normalized read counts (preferentially CPM normalized by TMM method included in any EdgeR analysis) and apply a stringent fold-change cut-off to determine which genes are more or less expressed depending on the condition.
  y <- DGEList(countGTF)
  # Calculate fold change
  y <- calcNormFactors(y, method = c('TMM'))
  y$samples$lib.size <- colSums(y$counts)
  logcpm <- cpm(y, prior.count=2, log=TRUE)
  logcpm[1:5,]
  cpm <- cpm(y, prior.count=2, log=FALSE)
  cpm[1:5, ]
  fold_change <- cpm[,1]/cpm[,2]
  res <- countGTF
  res$logFC <- log2(fold_change)
  res$FoldChange <- fold_change
  res$change <- with(res, ifelse(abs(res$logFC) > log2(1.5),
                                 ifelse(res$logFC > log2(1.5),'UP','DOWN'),'NOT'))
  table(res$change)
  res$geneSymbol <- df_in[,1]
  resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
  table(resInP$change)
  table(duplicated(resInP$geneSymbol))
  resInP_list[[c]] <- resInP
}
name = paste0(names(resInP_list), ' vs WT')
RMdeg_list[[i]] <- resInP_list 
name <- paste0(gsub(' ','_',(name)))
names(RMdeg_list[[i]] ) <- name

contrast_input_list[[i]] <- NULL
condition_input_list[[i]] <- condition_list


############
n = 57
############
#GSE149510# with RIP and Input
############
print(RMPdeg_datasets[n,])#
i = RMPdeg_datasets[n,'GSE']
cl = c('HCCLM3')#Liver cancer
RMp = c('ALKBH5') 
RM = 'm6A'
condition = c(rep('ALKBH5OE', 2), rep('WT', 2)) 
RMP_allDat[grep('GSE149510', RMP_allDat$study_alias ) 
           ,c("study_alias", "sample_alias", "sample_title")]
# View(allDat[[i]])
files <- list.files(file.path(data_directory,i))
files
getwd()
setwd(data_directory)
setwd(i)
skip <- RMP_allDat[grepl('GSE149510', RMP_allDat$study_alias )  &  grepl('Input', RMP_allDat$sample_title)
                   , "sample_alias"]
paste0(skip$sample_alias, collapse = " ")
skip <- RMP_allDat[grepl('GSE149510', RMP_allDat$study_alias )  &  grepl('IP', RMP_allDat$sample_title)
                   , "sample_alias"]
paste0(skip$sample_alias, collapse = " ")
ENA <- as.data.frame(RMP_allDat[grepl('GSE149510', RMP_allDat$study_alias ) 
                                , c("study_alias", "sample_alias", "sample_title")])
ENA <- ENA[order(ENA$sample_alias),]
ENA

c <- unique(condition[!condition %in% 'WT'])
###################################
#  Cell number-normalized (CNN)   #
#      MeRIP INPUT analysis       #
###################################
countGTF <- read.table( paste0("metaPlotR/Input_counts_", c, ".txt"),header=T)
rownames(countGTF) <- countGTF$Geneid
countGTF <- countGTF[, !colnames(countGTF) %in% c('Geneid','Chr','Start','End','Strand','Length')]
dim(countGTF) #28265     4
## Sample name
## Sample information
sample=data.frame(
  sample=colnames(countGTF),
  condition=condition)
## SummarizedExperiment input
dds<-DESeqDataSetFromMatrix(countData=countGTF,colData=sample,design=~condition)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 10
table(keep)
dds <- dds[keep,]
## Note on factor levels
dds$condition <- factor(dds$condition, levels = unique(condition))
dds$condition <- relevel(dds$condition, ref = "WT")
# Differential expression analysis-------------------------------------------
dds <- DESeq(dds)

res <-  results(dds, contrast=c("condition",unique(condition)))
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
table(res$change)
table(is.na(res$change))
# resOrdered <- as.data.frame(res[order(res$padj), ])
resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
resInP$geneSymbol <- rownames(resInP)
name = res@elementMetadata$description[2]
name = paste0(cl, '_',str_split(name, "condition ", simplify = T)[,2])
RMdeg_list[[i]] <- list( resInP)
names(RMdeg_list[[i]] ) <- paste0(gsub(' ','_',(name)))
c <- unique(condition[!condition %in% 'WT'])
countGTF_list <- list(colnames(countGTF)); names(countGTF_list) <- c
condition_list <- list(condition); names(condition_list) <- c
contrast_input_list[[i]] <- countGTF_list
condition_input_list[[i]] <- condition_list

###################################
#  Cell number-normalized (CNN)   #
#   MeRIP RIP analysis - Part1    #
###################################
## Sample information
merge_peak <- read.delim(paste0("metaPlotR/", c, "_newPeakFile.txt"),header=F)
countMatrix <- read.table(paste0("metaPlotR/m6A_counts_", c, ".txt"),header=T)
table(countMatrix$Geneid == merge_peak[-1,1])
table(duplicated(merge_peak$V1[-1]))
dim(countMatrix)
colnames(countMatrix)
rownames(countMatrix) <- countMatrix$Geneid
countMatrix <- countMatrix[, !colnames(countMatrix) %in% c('Geneid','Chr','Start','End','Strand','Length')]
# 1. 做好分组因子即可
groupList<-factor(condition)
colData <- data.frame(row.names=colnames(countMatrix), groupList=groupList)
# 2. Exploring Method 1: DESeq2-------------------------------------------
# 构建 DESeqDataSet 对象
dds <- DESeqDataSetFromMatrix(countData = countMatrix,
                              colData = colData,
                              design = ~ groupList)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 1
table(keep)
dds <- dds[keep,]
# 接着计算每个样本的 library size 用于后续的标准化。library size 等于落在该样本所有peaks上的reads的总数，即 counts 矩阵中每列的加和。如果bFullLibrarySize设为TRUE，则会使用library的总reads数（根据原BAM/BED文件统计获得）。然后使用estimateDispersions估计统计分布，需要将参数fitType设为’local’
dds<-estimateSizeFactors(dds)
dds<-estimateDispersions(dds,fitType="local")
dds <- nbinomWaldTest(dds)
res <- results(dds)
## 提取你想要的差异分析结果，我们这里是treated组对untreated组进行比较
res <-  results(dds, contrast=c("groupList", unique(condition)))
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
res_DESeq2 = as.data.frame(res)
table(res_DESeq2$change)#Not only
res_DESeq2 = res_DESeq2[res_DESeq2$change %in% c('UP', 'DOWN'), ]
# 2. Exploring Method 2: limma-------------------------------------------
library(edgeR)
d0 <- DGEList(counts = countMatrix)
dim(d0)
# Filtering low count peaks
cutoff <- 1
drop <- which(apply(cpm(d0), 1, max) < cutoff)
print(length(drop))
if (length(drop) >0) {
  d <- d0[-drop,]
} else{ d <- d0}
dim(d) #number of peaks left: 8501
d <- calcNormFactors(d, method = c('TMM'))
colnames(countMatrix)
genotype <- groupList
table(genotype)
mm <- model.matrix(~0 + genotype)#
v <- voom(d, mm, plot=F, normalize="quantile")#voom from limma:Transform count data to log2-counts per million (logtpm)
fit <- lmFit(v, mm)
head(coef(fit))
contr <- makeContrasts(paste0(colnames(coef(fit)), collapse  = ' - '), levels = colnames(coef(fit))) #Change if other comparison
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)
length(which(top.table$adj.P.Val < 0.05 & abs(top.table$logFC) > log2(1.5))) #how many are significant?
res_limma <- as.data.frame(top.table)
res_limma$change <- with(res_limma, ifelse(res_limma$adj.P.Val < 0.05 & abs(res_limma$logFC) > log2(1.5),
                                           ifelse(res_limma$logFC > log2(1.5),'UP','DOWN'),'NOT'))
table(res_limma$change)#NOT only
res_limma = res_limma[res_limma$change %in% c('UP', 'DOWN'), ]

RMdem_list[[i]] <- list(res_DESeq2, res_limma)
name = res@elementMetadata$description[2]
name = paste0(cl, '_',str_split(name, "groupList ", simplify = T)[,2])
name= paste0(gsub(' ','_',(name)))
names(RMdem_list[[i]] ) <- paste0(name, '_', c('DESeq2', 'limma'))
RMdep_list[[i]] <- list( merge_peak )
names(RMdep_list[[i]] ) <- name
countMatrix_list <- list(colnames(countMatrix)); names(countMatrix_list) <- c
condition_list <- list(condition); names(condition_list) <- c
contrast_IP_list[[i]] <- countMatrix_list
condition_IP_list[[i]] <- condition_list

############
n = 58
############
#GSE128575# with RIP and Input
############
print(RMPdeg_datasets[n,])# hg38
i = RMPdeg_datasets[n,'GSE']
cl = 'THP1'#c('MOLM-13 and THP1 and MV4:11')#AML
RMp = c('ALKBH5') 
RM = 'm6A'
condition_list = list(  
  ALKBH5KD = c(rep('ALKBH5KD', 2), rep('WT', 2)),
  ALKBH5RS = c(rep('ALKBH5RS', 2), rep('WT', 2)) 
)
RMP_allDat[grep('GSE128520', RMP_allDat$study_alias ) 
           ,c("study_alias", "sample_alias", "sample_title")]
ENA <- as.data.frame(RMP_allDat[grepl('GSE128520|GSE128574', RMP_allDat$study_alias ) 
                                , c("study_alias", "sample_alias", "sample_title")])
ENA <- ENA[order(ENA$sample_alias),]
ENA

# View(allDat[[i]])
skip <- ENA[ grepl('RNA-seq', ENA$sample_title) ,]
paste0(skip$sample_alias, collapse = " ")
skip <- ENA[ grepl('m6A', ENA$sample_title) ,]
paste0(skip$sample_alias, collapse = " ")

y <- ENA[ grepl('A5-res|shCtrl',ENA$sample_title) & grepl('RNA',ENA$sample_title),]
y <- ENA[ grepl('shA5|shCtrl',ENA$sample_title) & grepl('RNA',ENA$sample_title),]
y <- ENA[ grepl('A5-res|shCtrl',ENA$sample_title) & grepl('m6A',ENA$sample_title),]
y <- ENA[ grepl('shA5|shCtrl',ENA$sample_title) & grepl('m6A',ENA$sample_title),]
y
x <- paste0('$outdir/', y$study_alias, '/', y$sample_alias, '_IP.bam')
x <- paste0('$outdir/', y$study_alias, '/', y$sample_alias, '_input.bam')
paste(x, collapse = ' ')

files <- list.files(file.path(data_directory,i))
files
getwd()
setwd(data_directory)
setwd(i)
c <- unique(condition[!condition %in% 'WT'])

countfiles <- list.files('./metaPlotR/')
countfiles <- countfiles[grep('counts', countfiles)]
countfiles <- countfiles[! grepl('summary', countfiles)] 
countfiles
###################################
#  Cell number-normalized (CNN)   #
#      MeRIP INPUT analysis       #
###################################
resInP_list <- vector("list", length(condition_list))
names(resInP_list) <- names(condition_list)
countGTF_list <- vector("list", length(condition_list))
names(countGTF_list) <- names(condition_list)
for(c in names(condition_list)){
  countGTF <- read.table( paste0("metaPlotR/Input_counts_", c, ".txt"),header=T)
  dim(countGTF)
  colnames(countGTF)
  rownames(countGTF) <- countGTF$Geneid
  countGTF <- countGTF[, !colnames(countGTF) %in% c('Geneid','Chr','Start','End','Strand','Length')]
  dim(countGTF) # 28265     4
  condition <- condition_list[[c]]
  ## Sample information
  sample=data.frame(
    sample=colnames(countGTF),
    condition=condition)
  ## SummarizedExperiment input
  dds<-DESeqDataSetFromMatrix(countData=countGTF,colData=sample,design=~condition)
  ## Pre-filtering
  keep <- rowSums(counts(dds)) >= 10
  table(keep)
  dds <- dds[keep,]
  ## Note on factor levels
  dds$condition <- factor(dds$condition, levels = unique(condition))
  dds$condition <- relevel(dds$condition, ref = "WT")
  # Differential expression analysis-------------------------------------------
  dds <- DESeq(dds)
  res <-  results(dds, contrast=c("condition",unique(condition)))
  res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
  length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
  res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                                 ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
  table(res$change)
  table(is.na(res$change))
  # resOrdered <- as.data.frame(res[order(res$padj), ])
  resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
  resInP$geneSymbol <- rownames(resInP)
  table(resInP$change)
  resInP_list[[c]] <- resInP
  countGTF_list[[c]] <- countGTF 
}
name = paste0(cl, '_',names(resInP_list), ' vs WT')
RMdeg_list[[i]] <- resInP_list 
names(RMdeg_list[[i]] ) <- paste0(gsub(' ','_',(name)))
contrast_input_list[[i]] <- lapply(countGTF_list, extract_names)
condition_input_list[[i]] <- condition_list

###################################
#  Cell number-normalized (CNN)   #
#   MeRIP RIP analysis - Part1    #
###################################
resMeR_list <- vector("list", length(condition_list))
names(resMeR_list) <- names(condition_list)
countMatrix_list <- vector("list", length(condition_list))
names(countMatrix_list) <- names(condition_list)
resPeak_list <- vector("list",length(condition_list))
names(resPeak_list) <- names(condition_list)
for(c in names(condition_list)){
  merge_peak <- read.delim(paste0("metaPlotR/", c, "_newPeakFile.txt"),header=F)
  countMatrix <- read.table(paste0("metaPlotR/m6A_counts_", c, ".txt"),header=T)
  table(countMatrix$Geneid == merge_peak[-1,1])
  table(duplicated(merge_peak[-1,1]))
  dim(countMatrix)
  colnames(countMatrix)
  rownames(countMatrix) <- countMatrix$Geneid
  countMatrix <- countMatrix[, !colnames(countMatrix) %in% c('Geneid','Chr','Start','End','Strand','Length')]
  # 1. 做好分组因子即可
  condition <- condition_list[[c]]
  groupList<-factor(condition)
  colData <- data.frame(row.names=colnames(countMatrix), groupList=groupList)
  
  # 2. Exploring Method 1: DESeq2-------------------------------------------
  # 构建 DESeqDataSet 对象
  dds <- DESeqDataSetFromMatrix(countData = countMatrix,
                                colData = colData,
                                design = ~ groupList)
  ## Pre-filtering
  keep <- rowSums(counts(dds)) >= 1
  table(keep)
  dds <- dds[keep,]
  # 接着计算每个样本的 library size 用于后续的标准化。library size 等于落在该样本所有peaks上的reads的总数，即 counts 矩阵中每列的加和。如果bFullLibrarySize设为TRUE，则会使用library的总reads数（根据原BAM/BED文件统计获得）。然后使用estimateDispersions估计统计分布，需要将参数fitType设为’local’
  dds<-estimateSizeFactors(dds)
  dds<-estimateDispersions(dds,fitType="local")
  dds <- nbinomWaldTest(dds)
  res <- results(dds)
  ## 提取你想要的差异分析结果，我们这里是treated组对untreated组进行比较
  res <-  results(dds, contrast=c("groupList", unique(condition)))
  res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
  length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
  res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                                 ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
  res_DESeq2 = as.data.frame(res)
  table(res_DESeq2$change)
  res_DESeq2 = res_DESeq2[res_DESeq2$change %in% c('UP', 'DOWN'), ]
  # 2. Exploring Method 2: limma-------------------------------------------
  library(edgeR)
  d0 <- DGEList(counts = countMatrix)
  dim(d0)
  # Filtering low count peaks
  cutoff <- 1           
  drop <- which(apply(cpm(d0), 1, max) < cutoff)
  print(length(drop))
  if (length(drop) >0) {
    d <- d0[-drop,] 
  } else{ d <- d0}
  dim(d) #number of peaks left: 8501      
  d <- calcNormFactors(d, method = c('TMM'))
  colnames(countMatrix)
  genotype <- groupList
  table(genotype)
  mm <- model.matrix(~0 + genotype)#
  v <- voom(d, mm, plot=F, normalize="quantile")#voom from limma:Transform count data to log2-counts per million (logtpm)
  fit <- lmFit(v, mm)
  head(coef(fit))
  contr <- makeContrasts(paste0(colnames(coef(fit)), collapse  = ' - '), levels = colnames(coef(fit))) #Change if other comparison
  tmp <- contrasts.fit(fit, contr)
  tmp <- eBayes(tmp)
  top.table <- topTable(tmp, sort.by = "P", n = Inf)
  head(top.table, 20)
  length(which(top.table$adj.P.Val < 0.05 & abs(top.table$logFC) > log2(1.5))) #how many are significant?
  res_limma <- as.data.frame(top.table)
  res_limma$change <- with(res_limma, ifelse(res_limma$adj.P.Val < 0.05 & abs(res_limma$logFC) > log2(1.5),
                                             ifelse(res_limma$logFC > log2(1.5),'UP','DOWN'),'NOT'))
  table(res_limma$change)
  res_limma = res_limma[res_limma$change %in% c('UP', 'DOWN'), ]
  resMeR_list[[c]] <- list(res_DESeq2, res_limma)
  countMatrix_list[[c]] <- countMatrix
  resPeak_list[[c]] <- merge_peak
  name = paste0(c, ' vs WT')
  name = paste0(gsub(' ','_', name ))
  names(resMeR_list[[c]] ) <- paste0(name, '_', c('DESeq2', 'limma'))
}
names(resMeR_list) <- NULL
RMdem_list[[i]] <- as.list(unlist(resMeR_list, recursive = FALSE))
names(RMdem_list[[i]] ) <- paste0(cl, '_', names(RMdem_list[[i]]))
RMdep_list[[i]] <- resPeak_list 
name <- unique(gsub('(_limma|_DESeq2)', '', names(RMdem_list[[i]] )))
names(RMdep_list[[i]] ) <- name
contrast_IP_list[[i]] <- lapply(countMatrix_list, extract_names)
condition_IP_list[[i]] <- condition_list

############
n = 59
############
#GSE130172# with RIP and Input
############
print(RMPdeg_datasets[n,])# hg38
i = RMPdeg_datasets[n,'GSE']
cl = c('MDA-MB-231')#Breast Cancer Cell
RMp = c('YTHDF3') 
RM = 'm6A'
condition = c(rep('YTHDF3PD', 2), rep('WT', 2)) 
RMP_allDat[grep('GSE130172', RMP_allDat$study_alias ) 
           ,c("study_alias", "sample_alias", "sample_title")]
RMP_allDat[grep('GSE130171', RMP_allDat$study_alias ) 
           ,c("study_alias", "sample_alias", "sample_title")]
skip <- RMP_allDat[grepl('GSE130171', RMP_allDat$study_alias )  &  grepl('Input', RMP_allDat$sample_title)
                   , "sample_alias"]
paste0(skip$sample_alias, collapse = " ")
skip <- RMP_allDat[grepl('GSE130171', RMP_allDat$study_alias )  &  grepl('IP', RMP_allDat$sample_title)
                   , "sample_alias"]
paste0(skip$sample_alias, collapse = " ")

# View(allDat[[i]])
files <- list.files(file.path(data_directory,i))
files
getwd()
setwd(data_directory)
setwd(i)
###################################
#       RNA/RIP-seq analysis      #
###################################
condition = c(rep('YTHDF3PD', 2), rep('WT', 4)) 
DEG <- read_excel(file.path(data_directory,i, 'GSE130171_RIPseq_GeneExp.xlsx'),sheet = 1)
df_in <- read.table(file.path(data_directory,i, 'metaPlotR/Input_counts_YTHDF3PD.txt'),header = T)
## Raw counts
#判断数据是否需要转换
ex <- DEG[,-(1:3)]
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
# Check if LogC is TRUE
if (LogC) {
  # Replace values in 'ex' that are less than or equal to 0 with NaN
  ex[ex <= 0] <- NaN
  print("log2 transform unfinished")
} else {
  print("log2 transform not needed")
}
##Gene count matrix
countGTF=DEG[,c('RIPSZ','RIPZ', 'inputSZ', 'inputZ')]
countGTF <- df_in[, !colnames(df_in) %in% c('Geneid','Chr','Start','End','Strand','Length')]
rownames(countGTF) <- df_in$Geneid
dim(countGTF) #28278     4
colSums(countGTF)
table(colSums(countGTF)=="1e+06")
## Gene name
library(data.table)
table(duplicated(df_in$Hgnc_symbol))
table(duplicated(df_in$Geneid))
dim(countGTF) # 
## Sample name
colnames(countGTF)
## Sample information
sample=data.frame(
  sample=colnames(countGTF),
  condition=condition)
## SummarizedExperiment input
dds<-DESeqDataSetFromMatrix(countData=countGTF,colData=sample,design=~condition)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 10
table(keep)
dds <- dds[keep,]
## Note on factor levels
dds$condition <- factor(dds$condition, levels = unique(condition))
dds$condition <- relevel(dds$condition, ref = "WT")
# Differential expression analysis-------------------------------------------
dds <- DESeq(dds)
# res<-results(dds,name = "condition_TRMT61AOE_vs_WT")
res <-  results(dds, contrast=c("condition",unique(condition)))
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
table(res$change)
table(is.na(res$change))
# resOrdered <- as.data.frame(res[order(res$padj), ])
resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
resInP$geneSymbol <- rownames(resInP)
name = res@elementMetadata$description[2]
name = paste0(cl, '_',str_split(name, "condition ", simplify = T)[,2])
RMdeg_list[[i]] <- list( resInP)
names(RMdeg_list[[i]] ) <- paste0(gsub(' ','_',(name)))
c <- unique(condition[!condition %in% 'WT'])
countGTF_list <- list(colnames(countGTF)); names(countGTF_list) <- c
condition_list <- list(condition); names(condition_list) <- c
contrast_input_list[[i]] <- countGTF_list
condition_input_list[[i]] <- condition_list

###################################
#  Cell number-normalized (CNN)   #
#   MeRIP RIP analysis - Part1    #
###################################
condition = c(rep('YTHDF3PD', 2), rep('WT', 2)) 
## Sample information
merge_peak <- read.delim(paste0("metaPlotR/", c, "_newPeakFile.txt"),header=F)
countMatrix <- read.table(paste0("metaPlotR/m6A_counts_", c, ".txt"),header=T)
table(countMatrix$Geneid == merge_peak[-1,1])
table(duplicated(merge_peak$V1[-1]))
dim(countMatrix)
colnames(countMatrix)
rownames(countMatrix) <- countMatrix$Geneid
countMatrix <- countMatrix[, !colnames(countMatrix) %in% c('Geneid','Chr','Start','End','Strand','Length')]
# 1. 做好分组因子即可
groupList<-factor(condition)
colData <- data.frame(row.names=colnames(countMatrix), groupList=groupList)
# 2. Exploring Method 1: DESeq2-------------------------------------------
# 构建 DESeqDataSet 对象
dds <- DESeqDataSetFromMatrix(countData = countMatrix,
                              colData = colData,
                              design = ~ groupList)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 1
table(keep)
dds <- dds[keep,]
# 接着计算每个样本的 library size 用于后续的标准化。library size 等于落在该样本所有peaks上的reads的总数，即 counts 矩阵中每列的加和。如果bFullLibrarySize设为TRUE，则会使用library的总reads数（根据原BAM/BED文件统计获得）。然后使用estimateDispersions估计统计分布，需要将参数fitType设为’local’
dds<-estimateSizeFactors(dds)
dds<-estimateDispersions(dds,fitType="local")
dds <- nbinomWaldTest(dds)
res <- results(dds)
## 提取你想要的差异分析结果，我们这里是treated组对untreated组进行比较
res <-  results(dds, contrast=c("groupList", unique(condition)))
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
res_DESeq2 = as.data.frame(res)
table(res_DESeq2$change)#Not only
res_DESeq2 = res_DESeq2[res_DESeq2$change %in% c('UP', 'DOWN'), ]
# 2. Exploring Method 2: limma-------------------------------------------
library(edgeR)
d0 <- DGEList(counts = countMatrix)
dim(d0)
# Filtering low count peaks
cutoff <- 1
drop <- which(apply(cpm(d0), 1, max) < cutoff)
print(length(drop))
if (length(drop) >0) {
  d <- d0[-drop,]
} else{ d <- d0}
dim(d) #number of peaks left: 8501
d <- calcNormFactors(d, method = c('TMM'))
colnames(countMatrix)
genotype <- groupList
table(genotype)
mm <- model.matrix(~0 + genotype)#
v <- voom(d, mm, plot=F, normalize="quantile")#voom from limma:Transform count data to log2-counts per million (logtpm)
fit <- lmFit(v, mm)
head(coef(fit))
contr <- makeContrasts(paste0(colnames(coef(fit)), collapse  = ' - '), levels = colnames(coef(fit))) #Change if other comparison
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)
length(which(top.table$adj.P.Val < 0.05 & abs(top.table$logFC) > log2(1.5))) #how many are significant?
res_limma <- as.data.frame(top.table)
res_limma$change <- with(res_limma, ifelse(res_limma$adj.P.Val < 0.05 & abs(res_limma$logFC) > log2(1.5),
                                           ifelse(res_limma$logFC > log2(1.5),'UP','DOWN'),'NOT'))
table(res_limma$change)#NOT only
res_limma = res_limma[res_limma$change %in% c('UP', 'DOWN'), ]

RMdem_list[[i]] <- list(res_DESeq2, res_limma)
name = res@elementMetadata$description[2]
name = paste0(cl, '_',str_split(name, "groupList ", simplify = T)[,2])
name= paste0(gsub(' ','_',(name)))
names(RMdem_list[[i]] ) <- paste0(name, '_', c('DESeq2', 'limma'))
RMdep_list[[i]] <- list( merge_peak )
names(RMdep_list[[i]] ) <- name
countMatrix_list <- list(colnames(countMatrix)); names(countMatrix_list) <- c
condition_list <- list(condition); names(condition_list) <- c
contrast_IP_list[[i]] <- countMatrix_list
condition_IP_list[[i]] <- condition_list

# ############
# n = 60
# ############
# #GSE122744# without RIP and Input
# ############
# print(RMPdeg_datasets[n,])# skip?
# i = RMPdeg_datasets[n,'GSE']
# cl = c('HEK293T')#
# RMp = c('METTL3') 
# RM = 'm6A'
# # condition = c(rep('METTL3KD', 1), rep('WT', 1)) 
# RMP_allDat[grep('GSE122744', RMP_allDat$study_alias ) 
#            ,c("study_alias", "sample_alias", "sample_title")]
# # View(allDat[[i]])
# files <- list.files(file.path(data_directory,i))
# files
# getwd()
# setwd(data_directory)
# setwd(i)
#
# ############
# n = 61
# ############
# #PRJNA498900# without RIP and Input
# ############
# print(RMPdeg_datasets[n,])# skip?
# i = RMPdeg_datasets[n,'GSE']
# cl = c('BCa5637 and UM-UC-3 and SV-HUC-1')#
# RMp = c('METTL3') 
# RM = 'm6A'
# condition = c(rep('METTL3KD', 1), rep('WT', 1)) 
# RMP_allDat[grep('PRJNA498900', RMP_allDat$study_alias ) 
#            ,c("study_alias", "sample_alias", "sample_title")]
# # View(allDat[[i]])
# files <- list.files(file.path(data_directory,i))
# files
# getwd()
# setwd(data_directory)
# setwd(i)

############
n = 62
############
#GSE141994# with RIP and Input
############
print(RMPdeg_datasets[n,])# 
i = RMPdeg_datasets[n,'GSE']
cl = c('IMR90')#Embryonic lung fibroblasts
RMp = c('METTL3','METTL14') 
RM = 'm6A'
condition_list = list(METTL14KD = c(rep('METTL14KD', 1), rep('WT', 1)), 
                      METTL3KD = c(rep('METTL3KD', 1), rep('WT', 1)) 
)
RMP_allDat[grep('GSE141994', RMP_allDat$study_alias ) 
           ,c("study_alias", "sample_alias", "sample_title")]
# View(allDat[[i]])
files <- list.files(file.path(data_directory,i))
files
getwd()
setwd(data_directory)
setwd(i)
countfiles <- list.files('./metaPlotR/')
countfiles <- countfiles[grep('counts', countfiles)]
countfiles <- countfiles[! grepl('summary', countfiles)] 
countfiles
###################################
#  Cell number-normalized (CNN)   #
#      MeRIP INPUT analysis       #
###################################
resInP_list <- vector("list", length(condition_list))
names(resInP_list) <- names(condition_list)
countGTF_list <- vector("list", length(condition_list))
names(countGTF_list) <- names(condition_list)
for(c in names(condition_list)){
  countGTF <- read.table( paste0("metaPlotR/Input_counts_", c, ".txt"),header=T)
  dim(countGTF)
  colnames(countGTF)
  rownames(countGTF) <- countGTF$Geneid
  countGTF <- countGTF[, !colnames(countGTF) %in% c('Geneid','Chr','Start','End','Strand','Length')]
  dim(countGTF) # 28265     4
  condition <- condition_list[[c]]
  # 1. Differential gene analysis with edgeR
  # Without replicates, you cannot estimate which genes are differentially expressed using EdgeR or DESeq2. You can only calculate fold changes based on normalized read counts (preferentially CPM normalized by TMM method included in any EdgeR analysis) and apply a stringent fold-change cut-off to determine which genes are more or less expressed depending on the condition.
  y <- DGEList(countGTF)
  # Calculate fold change
  y <- calcNormFactors(y, method = c('TMM'))
  y$samples$lib.size <- colSums(y$counts)
  logcpm <- cpm(y, prior.count=2, log=TRUE)
  logcpm[1:5,]
  cpm <- cpm(y, prior.count=2, log=FALSE)
  cpm[1:5, ]
  fold_change <- cpm[,1]/cpm[,2]
  res <- countGTF
  res$logFC <- log2(fold_change)
  res$FoldChange <- fold_change
  res$change <- with(res, ifelse(abs(res$logFC) > log2(1.5),
                                 ifelse(res$logFC > log2(1.5),'UP','DOWN'),'NOT'))
  table(res$change)
  resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
  resInP$geneSymbol <- rownames(resInP)
  table(resInP$change)
  resInP_list[[c]] <- resInP
  countGTF_list[[c]] <- countGTF
  
}
name = paste0(cl, '_' ,names(resInP_list), ' vs WT')
RMdeg_list[[i]] <- resInP_list 
names(RMdeg_list[[i]] ) <- paste0(gsub(' ','_',(name)))
contrast_input_list[[i]] <- lapply(countGTF_list, extract_names)
condition_input_list[[i]] <- condition_list
###################################
#  Cell number-normalized (CNN)   #
#   MeRIP RIP analysis - Part1    #
###################################
resMeR_list <- vector("list", length(condition_list))
names(resMeR_list) <- names(condition_list)
countMatrix_list <- vector("list", length(condition_list))
names(countMatrix_list) <- names(condition_list)
resPeak_list <- vector("list",length(condition_list))
names(resPeak_list) <- names(condition_list)
for(c in names(condition_list)){
  merge_peak <- read.delim(paste0("metaPlotR/", c, "_newPeakFile.txt"),header=F)
  countMatrix <- read.table(paste0("metaPlotR/m6A_counts_", c, ".txt"),header=T)
  table(countMatrix$Geneid == merge_peak[-1,1])
  table(duplicated(merge_peak[-1,1]))
  dim(countMatrix)
  colnames(countMatrix)
  rownames(countMatrix) <- countMatrix$Geneid
  countMatrix <- countMatrix[, !colnames(countMatrix) %in% c('Geneid','Chr','Start','End','Strand','Length')]
  # 1. 做好分组因子即可
  condition <- condition_list[[c]]
  groupList<-factor(condition)
  colData <- data.frame(row.names=colnames(countMatrix), groupList=groupList)
  # if(c %in% c("HEK293_WTAPKD","HEK293_METTL3KD")){
  # 2. Differential gene analysis with edgeR
  # Without replicates, you cannot estimate which genes are differentially expressed using EdgeR or DESeq2. You can only calculate fold changes based on normalized read counts (preferentially CPM normalized by TMM method included in any EdgeR analysis) and apply a stringent fold-change cut-off to determine which genes are more or less expressed depending on the condition.
  y <- DGEList(countMatrix)
  # Calculate fold change
  y <- calcNormFactors(y, method = c('TMM'))
  y$samples$lib.size <- colSums(y$counts)
  logcpm <- cpm(y, prior.count=2, log=TRUE)
  logcpm[1:5,]
  cpm <- cpm(y, prior.count=2, log=FALSE)
  cpm[1:5, ]
  # fold_change <- logcpm[,'Ratio1']/logcpm[,'Ratio2']
  fold_change <- cpm[,1]/cpm[,2]
  resMeR <- countMatrix
  resMeR$logFC <- log2(fold_change)
  resMeR$FoldChange <- fold_change
  resMeR$change <- with(resMeR, ifelse(abs(resMeR$logFC) > log2(1.5),
                                       ifelse(resMeR$logFC > log2(1.5),'UP','DOWN'),'NOT'))
  table(resMeR$change)
  resMeR = as.data.frame(resMeR[resMeR$change %in% c('UP', 'DOWN'), ])
  resMeR_list[[c]] <- list(resMeR)
  countMatrix_list[[c]] <- countMatrix
  resPeak_list[[c]] <- merge_peak
  name = paste0(cl,'_', c, ' vs WT')
  name = paste0(gsub(' ','_', name ))
  names(resMeR_list[[c]]) <- paste0(name, '_', c('limma'))
}
names(resMeR_list) <- NULL
RMdem_list[[i]] <- as.list(unlist(resMeR_list, recursive = FALSE))
names(RMdem_list[[i]] )
RMdep_list[[i]] <- resPeak_list 
name <- unique(gsub('(_limma|_DESeq2)', '', names(RMdem_list[[i]] )))
names(RMdep_list[[i]] ) <- name
contrast_IP_list[[i]] <- lapply(countMatrix_list, extract_names)
condition_IP_list[[i]] <- condition_list


############
n = 63
############
#GSE128580# with RIP and Input
############
print(RMPdeg_datasets[n,])#
i = RMPdeg_datasets[n,'GSE']
cl = c('PC-3','SK-BR-3')#Breast cancer
RMp = c('FTO')
RM = 'm6A'
condition_list = list(FTO.1KD = c(rep('FTOKD', 1), rep('WT', 1)), #PC-3
                      FTO.2KD = c(rep('FTOKD', 2), rep('WT', 2))) #SK-BR-2
RMP_allDat[grep('GSE128580', RMP_allDat$study_alias )
           ,c("study_alias", "sample_alias", "sample_title")]
# View(allDat[[i]])
RMP_allDat[grepl('GSE128580', RMP_allDat$study_alias ) & (!grepl('biopsy', RMP_allDat$sample_title))
           ,c("study_alias", "sample_alias", "sample_title")]

skip = RMP_allDat[grepl('GSE128580', RMP_allDat$study_alias ) & grepl('biopsy', RMP_allDat$sample_title)
                  ,c("study_alias", "sample_alias", "sample_title")]
paste0(skip$sample_alias, collapse = "|")

skip <- RMP_allDat[grepl('GSE128580', RMP_allDat$study_alias ) & (!grepl('biopsy', RMP_allDat$sample_title))
                   &  grepl('m6A-RIP Input', RMP_allDat$sample_title)
                   , "sample_alias"]
paste0(skip$sample_alias, collapse = " ")
skip <- RMP_allDat[grepl('GSE128580', RMP_allDat$study_alias ) & (!grepl('biopsy', RMP_allDat$sample_title))
                   &  grepl('m6A-RIP IP', RMP_allDat$sample_title)
                   , "sample_alias"]
paste0(skip$sample_alias, collapse = " ")
ENA <- as.data.frame(RMP_allDat[grepl('GSE128580', RMP_allDat$study_alias ) & (!grepl('biopsy', RMP_allDat$sample_title))
                                , c("study_alias", "sample_alias", "sample_title")])
ENA <- ENA[order(ENA$sample_alias),]
ENA
files <- list.files(file.path(data_directory,i))
files
getwd()
setwd(data_directory)
setwd(i)

countfiles <- list.files('./metaPlotR/')
countfiles <- countfiles[grep('counts', countfiles)]
countfiles <- countfiles[! grepl('summary', countfiles)] 
countfiles
###################################
#  Cell number-normalized (CNN)   #
#      MeRIP INPUT analysis       #
###################################
resInP_list <- vector("list", length(condition_list))
names(resInP_list) <- names(condition_list)
countGTF_list <- vector("list", length(condition_list))
names(countGTF_list) <- names(condition_list)
for(c in names(condition_list)){
  countGTF <- read.table( paste0("metaPlotR/Input_counts_", c, ".txt"),header=T)
  dim(countGTF)
  colnames(countGTF)
  rownames(countGTF) <- countGTF$Geneid
  countGTF <- countGTF[, !colnames(countGTF) %in% c('Geneid','Chr','Start','End','Strand','Length')]
  dim(countGTF) # 28265     4
  condition <- condition_list[[c]]
  if(c %in% c("FTO.1KD")){
    # 1. Differential gene analysis with edgeR
    # Without replicates, you cannot estimate which genes are differentially expressed using EdgeR or DESeq2. You can only calculate fold changes based on normalized read counts (preferentially CPM normalized by TMM method included in any EdgeR analysis) and apply a stringent fold-change cut-off to determine which genes are more or less expressed depending on the condition.
    y <- DGEList(countGTF)
    # Calculate fold change
    y <- calcNormFactors(y, method = c('TMM'))
    y$samples$lib.size <- colSums(y$counts)
    logcpm <- cpm(y, prior.count=2, log=TRUE)
    logcpm[1:5,]
    cpm <- cpm(y, prior.count=2, log=FALSE)
    cpm[1:5, ]
    fold_change <- cpm[,1]/cpm[,2]
    res <- countGTF
    res$logFC <- log2(fold_change)
    res$FoldChange <- fold_change
    res$change <- with(res, ifelse(abs(res$logFC) > log2(1.5),
                                   ifelse(res$logFC > log2(1.5),'UP','DOWN'),'NOT'))
    table(res$change)
    resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
    resInP$geneSymbol <- rownames(resInP)
    table(resInP$change)
    resInP_list[[c]] <- resInP
    countGTF_list[[c]] <- countGTF
    
  } else{
    
    ## Sample information
    sample=data.frame(
      sample=colnames(countGTF),
      condition=condition)
    ## SummarizedExperiment input
    dds<-DESeqDataSetFromMatrix(countData=countGTF,colData=sample,design=~condition)
    ## Pre-filtering
    keep <- rowSums(counts(dds)) >= 10
    table(keep)
    dds <- dds[keep,]
    ## Note on factor levels
    dds$condition <- factor(dds$condition, levels = unique(condition))
    dds$condition <- relevel(dds$condition, ref = "WT")
    # Differential expression analysis-------------------------------------------
    dds <- DESeq(dds)
    res <-  results(dds, contrast=c("condition",unique(condition)))
    res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
    length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
    res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                                   ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
    table(res$change)
    table(is.na(res$change))
    # resOrdered <- as.data.frame(res[order(res$padj), ])
    resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
    resInP$geneSymbol <- rownames(resInP)
    table(resInP$change)
    resInP_list[[c]] <- resInP
    countGTF_list[[c]] <- countGTF
  }
  
}
name = paste0(names(resInP_list), ' vs WT')
RMdeg_list[[i]] <- resInP_list 
name <-  paste0(cl,'_', gsub(".1KD|.2KD","KD", name))
names(RMdeg_list[[i]] ) <- gsub(' ','_', name) 
contrast_input_list[[i]] <- lapply(countGTF_list, extract_names)
name = names(countGTF_list)
names(contrast_input_list[[i]]) <- paste0(cl,'_', gsub(".1KD|.2KD","KD", name))
condition_input_list[[i]] <- condition_list
names(condition_input_list[[i]] ) <- names(contrast_input_list[[i]]) 
###################################
#  Cell number-normalized (CNN)   #
#   MeRIP RIP analysis - Part1    #
###################################
resMeR_list <- vector("list", length(condition_list))
names(resMeR_list) <- names(condition_list)
countMatrix_list <- vector("list", length(condition_list))
names(countMatrix_list) <- names(condition_list)
resPeak_list <- vector("list",length(condition_list))
names(resPeak_list) <- names(condition_list)
for(c in names(condition_list)){
  merge_peak <- read.delim(paste0("metaPlotR/", c, "_newPeakFile.txt"),header=F)
  countMatrix <- read.table(paste0("metaPlotR/m6A_counts_", c, ".txt"),header=T)
  table(countMatrix$Geneid == merge_peak[-1,1])
  table(duplicated(merge_peak[-1,1]))
  dim(countMatrix)
  colnames(countMatrix)
  rownames(countMatrix) <- countMatrix$Geneid
  countMatrix <- countMatrix[, !colnames(countMatrix) %in% c('Geneid','Chr','Start','End','Strand','Length')]
  # 1. 做好分组因子即可
  condition <- condition_list[[c]]
  groupList<-factor(condition)
  colData <- data.frame(row.names=colnames(countMatrix), groupList=groupList)
  if(c %in% c("FTO.1KD")){
    # 2. Differential gene analysis with edgeR
    # Without replicates, you cannot estimate which genes are differentially expressed using EdgeR or DESeq2. You can only calculate fold changes based on normalized read counts (preferentially CPM normalized by TMM method included in any EdgeR analysis) and apply a stringent fold-change cut-off to determine which genes are more or less expressed depending on the condition.
    y <- DGEList(countMatrix)
    # Calculate fold change
    y <- calcNormFactors(y, method = c('TMM'))
    y$samples$lib.size <- colSums(y$counts)
    logcpm <- cpm(y, prior.count=2, log=TRUE)
    logcpm[1:5,]
    cpm <- cpm(y, prior.count=2, log=FALSE)
    cpm[1:5, ]
    # fold_change <- logcpm[,'Ratio1']/logcpm[,'Ratio2']
    fold_change <- cpm[,1]/cpm[,2]
    resMeR <- countMatrix
    resMeR$logFC <- log2(fold_change)
    resMeR$FoldChange <- fold_change
    resMeR$change <- with(resMeR, ifelse(abs(resMeR$logFC) > log2(1.5),
                                         ifelse(resMeR$logFC > log2(1.5),'UP','DOWN'),'NOT'))
    table(resMeR$change)
    resMeR = as.data.frame(resMeR[resMeR$change %in% c('UP', 'DOWN'), ])
    resMeR_list[[c]] <- list(resMeR)
    countMatrix_list[[c]] <- countMatrix
    resPeak_list[[c]] <- merge_peak
    name = paste0(c, ' vs WT')
    name = paste0(gsub(' ','_', name ))
    names(resMeR_list[[c]]) <- paste0(name, '_', c('limma'))
    
  } else{
    # 2. Exploring Method 1: DESeq2-------------------------------------------
    # 构建 DESeqDataSet 对象
    dds <- DESeqDataSetFromMatrix(countData = countMatrix,
                                  colData = colData,
                                  design = ~ groupList)
    ## Pre-filtering
    keep <- rowSums(counts(dds)) >= 1
    table(keep)
    dds <- dds[keep,]
    # 接着计算每个样本的 library size 用于后续的标准化。library size 等于落在该样本所有peaks上的reads的总数，即 counts 矩阵中每列的加和。如果bFullLibrarySize设为TRUE，则会使用library的总reads数（根据原BAM/BED文件统计获得）。然后使用estimateDispersions估计统计分布，需要将参数fitType设为’local’
    dds<-estimateSizeFactors(dds)
    dds<-estimateDispersions(dds,fitType="local")
    dds <- nbinomWaldTest(dds)
    res <- results(dds)
    ## 提取你想要的差异分析结果，我们这里是treated组对untreated组进行比较
    res <-  results(dds, contrast=c("groupList", unique(condition)))
    res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
    length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
    res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                                   ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
    res_DESeq2 = as.data.frame(res)
    table(res_DESeq2$change)
    res_DESeq2 = res_DESeq2[res_DESeq2$change %in% c('UP', 'DOWN'), ]
    # 2. Exploring Method 2: limma-------------------------------------------
    library(edgeR)
    d0 <- DGEList(counts = countMatrix)
    dim(d0)
    # Filtering low count peaks
    cutoff <- 1           
    drop <- which(apply(cpm(d0), 1, max) < cutoff)
    print(length(drop))
    if (length(drop) >0) {
      d <- d0[-drop,] 
    } else{ d <- d0}
    dim(d) #number of peaks left: 8501      
    d <- calcNormFactors(d, method = c('TMM'))
    colnames(countMatrix)
    genotype <- groupList
    table(genotype)
    mm <- model.matrix(~0 + genotype)#
    v <- voom(d, mm, plot=F, normalize="quantile")#voom from limma:Transform count data to log2-counts per million (logtpm)
    fit <- lmFit(v, mm)
    head(coef(fit))
    contr <- makeContrasts(paste0(colnames(coef(fit)), collapse  = ' - '), levels = colnames(coef(fit))) #Change if other comparison
    tmp <- contrasts.fit(fit, contr)
    tmp <- eBayes(tmp)
    top.table <- topTable(tmp, sort.by = "P", n = Inf)
    head(top.table, 20)
    length(which(top.table$adj.P.Val < 0.05 & abs(top.table$logFC) > log2(1.5))) #how many are significant?
    res_limma <- as.data.frame(top.table)
    res_limma$change <- with(res_limma, ifelse(res_limma$adj.P.Val < 0.05 & abs(res_limma$logFC) > log2(1.5),
                                               ifelse(res_limma$logFC > log2(1.5),'UP','DOWN'),'NOT'))
    table(res_limma$change)
    res_limma = res_limma[res_limma$change %in% c('UP', 'DOWN'), ]
    resMeR_list[[c]] <- list(res_DESeq2, res_limma)
    countMatrix_list[[c]] <- countMatrix
    resPeak_list[[c]] <- merge_peak
    name = paste0(c, ' vs WT')
    name = paste0(gsub(' ','_', name ))
    names(resMeR_list[[c]] ) <- paste0(name, '_', c('DESeq2', 'limma'))
    
  }
}
names(resMeR_list) <- NULL
RMdem_list[[i]] <- as.list(unlist(resMeR_list, recursive = FALSE))
names(RMdem_list[[i]] )
name <-  paste0( c( "PC-3","SK-BR-3","SK-BR-3"),'_', 
                 gsub(".1KD|.2KD","KD", names(RMdem_list[[i]])) )
names(RMdem_list[[i]] ) <- name
RMdep_list[[i]] <- resPeak_list 
name <- unique(gsub('(_limma|_DESeq2)', '', names(RMdem_list[[i]] )))
names(RMdep_list[[i]] ) <- name
contrast_IP_list[[i]] <- lapply(countMatrix_list, extract_names)
names(contrast_IP_list[[i]]) <- gsub('_vs_WT', '',names(RMdep_list[[i]] ) )
condition_IP_list[[i]] <- condition_list
names(condition_IP_list[[i]]) <- gsub('_vs_WT', '',names(RMdep_list[[i]] ) )

############
n = 64
############
#GSE166972# with RIP and Input
############
print(RMPdeg_datasets[n,])#
i = RMPdeg_datasets[n,'GSE']
cl = c('AGS')#Gastric cancer
RMp = c('YTHDF1')
RM = 'm6A'
condition = c(rep('YTHDF1KD', 3), rep('WT', 3))
RMP_allDat[grep('GSE166972', RMP_allDat$study_alias )
           ,c("study_alias", "sample_alias", "sample_title")]
# View(allDat[[i]])
RMP_allDat[grepl('GSE166972', RMP_allDat$study_alias ) & (!grepl('RNA-seq', RMP_allDat$sample_title))
           ,c("study_alias", "sample_alias", "sample_title")]

skip <- RMP_allDat[grepl('GSE166972', RMP_allDat$study_alias ) & (!grepl('RNA-seq', RMP_allDat$sample_title))
                   &  grepl('INPUT', RMP_allDat$sample_title)
                   , "sample_alias"]
paste0(skip$sample_alias, collapse = " ")
skip <- RMP_allDat[grepl('GSE166972', RMP_allDat$study_alias ) & (!grepl('RNA-seq', RMP_allDat$sample_title))
                   &  (!grepl('INPUT', RMP_allDat$sample_title))
                   , "sample_alias"]
paste0(skip$sample_alias, collapse = " ")
ENA <- as.data.frame(RMP_allDat[grepl('GSE166972', RMP_allDat$study_alias ) & (!grepl('RNA-seq', RMP_allDat$sample_title))
                                , c("study_alias", "sample_alias", "sample_title")])
ENA <- ENA[order(ENA$sample_alias),]
ENA
files <- list.files(file.path(data_directory,i))
files
getwd()
setwd(data_directory)
setwd(i)

c <- unique(condition[! condition %in% 'WT'])
###################################
#  Cell number-normalized (CNN)   #
#      MeRIP INPUT analysis       #
###################################
countGTF <- read.table( paste0("metaPlotR/Input_counts_", c, ".txt"),header=T)
rownames(countGTF) <- countGTF$Geneid
countGTF <- countGTF[, !colnames(countGTF) %in% c('Geneid','Chr','Start','End','Strand','Length')]
dim(countGTF) #28265     4
## Sample name
## Sample information
sample=data.frame(
  sample=colnames(countGTF),
  condition=condition)
## SummarizedExperiment input
dds<-DESeqDataSetFromMatrix(countData=countGTF,colData=sample,design=~condition)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 10
table(keep)
dds <- dds[keep,]
## Note on factor levels
dds$condition <- factor(dds$condition, levels = unique(condition))
dds$condition <- relevel(dds$condition, ref = "WT")
# Differential expression analysis-------------------------------------------
dds <- DESeq(dds)

res <-  results(dds, contrast=c("condition",unique(condition)))
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
table(res$change)
table(is.na(res$change))
# resOrdered <- as.data.frame(res[order(res$padj), ])
resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
resInP$geneSymbol <- rownames(resInP)
name = res@elementMetadata$description[2]
name = paste0(cl, '_',str_split(name, "condition ", simplify = T)[,2])
RMdeg_list[[i]] <- list( resInP)
names(RMdeg_list[[i]] ) <- paste0(gsub(' ','_',(name)))
c <- unique(condition[!condition %in% 'WT'])
countGTF_list <- list(colnames(countGTF)); names(countGTF_list) <- c
condition_list <- list(condition); names(condition_list) <- c
contrast_input_list[[i]] <- countGTF_list
condition_input_list[[i]] <- condition_list

###################################
#  Cell number-normalized (CNN)   #
#   MeRIP RIP analysis - Part1    #
###################################
## Sample information
merge_peak <- read.delim(paste0("metaPlotR/", c, "_newPeakFile.txt"),header=F)
countMatrix <- read.table(paste0("metaPlotR/m6A_counts_", c, ".txt"),header=T)
table(countMatrix$Geneid == merge_peak[-1,1])
table(duplicated(merge_peak$V1[-1]))
dim(countMatrix)
colnames(countMatrix)
rownames(countMatrix) <- countMatrix$Geneid
countMatrix <- countMatrix[, !colnames(countMatrix) %in% c('Geneid','Chr','Start','End','Strand','Length')]
# 1. 做好分组因子即可
groupList<-factor(condition)
colData <- data.frame(row.names=colnames(countMatrix), groupList=groupList)
# 2. Exploring Method 1: DESeq2-------------------------------------------
# 构建 DESeqDataSet 对象
dds <- DESeqDataSetFromMatrix(countData = countMatrix,
                              colData = colData,
                              design = ~ groupList)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 1
table(keep)
dds <- dds[keep,]
# 接着计算每个样本的 library size 用于后续的标准化。library size 等于落在该样本所有peaks上的reads的总数，即 counts 矩阵中每列的加和。如果bFullLibrarySize设为TRUE，则会使用library的总reads数（根据原BAM/BED文件统计获得）。然后使用estimateDispersions估计统计分布，需要将参数fitType设为’local’
dds<-estimateSizeFactors(dds)
dds<-estimateDispersions(dds,fitType="local")
dds <- nbinomWaldTest(dds)
res <- results(dds)
## 提取你想要的差异分析结果，我们这里是treated组对untreated组进行比较
res <-  results(dds, contrast=c("groupList", unique(condition)))
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
res_DESeq2 = as.data.frame(res)
table(res_DESeq2$change)#
res_DESeq2 = res_DESeq2[res_DESeq2$change %in% c('UP', 'DOWN'), ]
# 2. Exploring Method 2: limma-------------------------------------------
library(edgeR)
d0 <- DGEList(counts = countMatrix)
dim(d0)
# Filtering low count peaks
cutoff <- 1           
drop <- which(apply(cpm(d0), 1, max) < cutoff)
print(length(drop))
if (length(drop) >0) {
  d <- d0[-drop,] 
} else{ d <- d0}
dim(d) #number of peaks left: 8501      
d <- calcNormFactors(d, method = c('TMM'))
colnames(countMatrix)
genotype <- groupList
table(genotype)
mm <- model.matrix(~0 + genotype)#
v <- voom(d, mm, plot=F, normalize="quantile")#voom from limma:Transform count data to log2-counts per million (logtpm)
fit <- lmFit(v, mm)
head(coef(fit))
contr <- makeContrasts(paste0(colnames(coef(fit)), collapse  = ' - '), levels = colnames(coef(fit))) #Change if other comparison
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)
length(which(top.table$adj.P.Val < 0.05 & abs(top.table$logFC) > log2(1.5))) #how many are significant?
res_limma <- as.data.frame(top.table)
res_limma$change <- with(res_limma, ifelse(res_limma$adj.P.Val < 0.05 & abs(res_limma$logFC) > log2(1.5),
                                           ifelse(res_limma$logFC > log2(1.5),'UP','DOWN'),'NOT'))
table(res_limma$change)#
res_limma = res_limma[res_limma$change %in% c('UP', 'DOWN'), ]

RMdem_list[[i]] <- list(res_DESeq2, res_limma)
name = res@elementMetadata$description[2]
name = paste0(cl, '_',str_split(name, "groupList ", simplify = T)[,2])
name= paste0(gsub(' ','_',(name)))
names(RMdem_list[[i]] ) <- paste0(name, '_', c('DESeq2', 'limma'))
RMdep_list[[i]] <- list( merge_peak )
names(RMdep_list[[i]] ) <- name
countMatrix_list <- list(colnames(countMatrix)); names(countMatrix_list) <- c
condition_list <- list(condition); names(condition_list) <- c
contrast_IP_list[[i]] <- countMatrix_list
condition_IP_list[[i]] <- condition_list

############
n = 65
############
#GSE144984# with RIP and Input
############
print(RMPdeg_datasets[n,])# hg38
i = RMPdeg_datasets[n,'GSE']
cl = c('NOMO1','MOLM13')#AML
cl = 'NOMO1' #only RNA-seq for MOLM13
RMp = c('ALKBH5')
RM = 'm6A'
condition_list = list(
  ALKBH5KD = c(rep('ALKBH5KD', 1), rep('WT', 1)), #NOMO1 cells
  ALKBH5OE = c(rep('ALKBH5OE', 2), rep('WT', 1)) )#AML cell line NOMO-1
RMP_allDat[grep('GSE144963', RMP_allDat$study_alias )
           ,c("study_alias", "sample_alias", "sample_title")]
# View(allDat[[i]])
RMP_allDat[grep('GSE134763|GSE144963', RMP_allDat$study_alias )
           ,c("study_alias", "sample_alias", "sample_title")]
# RMP_allDat[grepl('GSE144968', RMP_allDat$study_alias )
#              ,c("study_alias", "sample_alias", "sample_title")]
RMP_allDat[grepl('GSE144963', RMP_allDat$study_alias )
           ,c("study_alias", "sample_alias", "sample_title")]

skip <- RMP_allDat[grepl('GSE134763|GSE144963', RMP_allDat$study_alias )
                   &  grepl('input|Input', RMP_allDat$sample_title)
                   , "sample_alias"]
paste0(skip$sample_alias, collapse = " ")
skip <- RMP_allDat[grepl('GSE134763|GSE144963', RMP_allDat$study_alias )
                   &  grepl('IP', RMP_allDat$sample_title)
                   , "sample_alias"]
paste0(skip$sample_alias, collapse = " ")
ENA <- as.data.frame(RMP_allDat[grep('GSE134763|GSE144963', RMP_allDat$study_alias )
                                , c("study_alias", "sample_alias", "sample_title")])
ENA <- ENA[order(ENA$sample_alias),]
ENA

files <- list.files(file.path(data_directory,i))
files
getwd()
setwd(data_directory)
setwd(i)

###################################
#  Cell number-normalized (CNN)   #
#      MeRIP INPUT analysis       #
###################################
resInP_list <- vector("list", length(condition_list))
names(resInP_list) <- names(condition_list)
countGTF_list <- vector("list", length(condition_list))
names(countGTF_list) <- names(condition_list)
for(c in names(condition_list)){
  countGTF <- read.table( paste0("metaPlotR/Input_counts_", c, ".txt"),header=T)
  dim(countGTF)
  colnames(countGTF)
  rownames(countGTF) <- countGTF$Geneid
  countGTF <- countGTF[, !colnames(countGTF) %in% c('Geneid','Chr','Start','End','Strand','Length')]
  dim(countGTF) # 28265     4
  condition <- condition_list[[c]]
  if(c %in% c("ALKBH5KD")){
    # 1. Differential gene analysis with edgeR
    # Without replicates, you cannot estimate which genes are differentially expressed using EdgeR or DESeq2. You can only calculate fold changes based on normalized read counts (preferentially CPM normalized by TMM method included in any EdgeR analysis) and apply a stringent fold-change cut-off to determine which genes are more or less expressed depending on the condition.
    y <- DGEList(countGTF)
    # Calculate fold change
    y <- calcNormFactors(y, method = c('TMM'))
    y$samples$lib.size <- colSums(y$counts)
    logcpm <- cpm(y, prior.count=2, log=TRUE)
    logcpm[1:5,]
    cpm <- cpm(y, prior.count=2, log=FALSE)
    cpm[1:5, ]
    fold_change <- cpm[,1]/cpm[,2]
    res <- countGTF
    res$logFC <- log2(fold_change)
    res$FoldChange <- fold_change
    res$change <- with(res, ifelse(abs(res$logFC) > log2(1.5),
                                   ifelse(res$logFC > log2(1.5),'UP','DOWN'),'NOT'))
    table(res$change)
    resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
    resInP$geneSymbol <- rownames(resInP)
    table(resInP$change)
    resInP_list[[c]] <- resInP
    countGTF_list[[c]] <- countGTF
    
  } else{
    ## Sample information
    sample=data.frame(
      sample=colnames(countGTF),
      condition=condition)
    ## SummarizedExperiment input
    dds<-DESeqDataSetFromMatrix(countData=countGTF,colData=sample,design=~condition)
    ## Pre-filtering
    keep <- rowSums(counts(dds)) >= 10
    table(keep)
    dds <- dds[keep,]
    ## Note on factor levels
    dds$condition <- factor(dds$condition, levels = unique(condition))
    dds$condition <- relevel(dds$condition, ref = "WT")
    # Differential expression analysis-------------------------------------------
    dds <- DESeq(dds)
    res <-  results(dds, contrast=c("condition",unique(condition)))
    res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
    length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
    res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                                   ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
    table(res$change)
    table(is.na(res$change))
    # resOrdered <- as.data.frame(res[order(res$padj), ])
    resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
    resInP$geneSymbol <- rownames(resInP)
    table(resInP$change)
    resInP_list[[c]] <- resInP
    countGTF_list[[c]] <- countGTF
  }
  
}
name = paste0(cl, '_',names(resInP_list), ' vs WT')
RMdeg_list[[i]] <- resInP_list 
names(RMdeg_list[[i]] ) <- paste0(gsub(' ','_',(name)))
contrast_input_list[[i]] <- lapply(countGTF_list, extract_names)
condition_input_list[[i]] <- condition_list

###################################
#  Cell number-normalized (CNN)   #
#   MeRIP RIP analysis - Part1    #
###################################
resMeR_list <- vector("list", length(condition_list))
names(resMeR_list) <- names(condition_list)
countMatrix_list <- vector("list", length(condition_list))
names(countMatrix_list) <- names(condition_list)
resPeak_list <- vector("list",length(condition_list))
names(resPeak_list) <- names(condition_list)
for(c in names(condition_list)){
  merge_peak <- read.delim(paste0("metaPlotR/", c, "_newPeakFile.txt"),header=F)
  countMatrix <- read.table(paste0("metaPlotR/m6A_counts_", c, ".txt"),header=T)
  table(countMatrix$Geneid == merge_peak[-1,1])
  table(duplicated(merge_peak[-1,1]))
  dim(countMatrix)
  colnames(countMatrix)
  rownames(countMatrix) <- countMatrix$Geneid
  countMatrix <- countMatrix[, !colnames(countMatrix) %in% c('Geneid','Chr','Start','End','Strand','Length')]
  # 1. 做好分组因子即可
  condition <- condition_list[[c]]
  groupList<-factor(condition)
  colData <- data.frame(row.names=colnames(countMatrix), groupList=groupList)
  if(c %in% c("ALKBH5KD")){
    # 2. Differential gene analysis with edgeR
    # Without replicates, you cannot estimate which genes are differentially expressed using EdgeR or DESeq2. You can only calculate fold changes based on normalized read counts (preferentially CPM normalized by TMM method included in any EdgeR analysis) and apply a stringent fold-change cut-off to determine which genes are more or less expressed depending on the condition.
    y <- DGEList(countMatrix)
    # Calculate fold change
    y <- calcNormFactors(y, method = c('TMM'))
    y$samples$lib.size <- colSums(y$counts)
    logcpm <- cpm(y, prior.count=2, log=TRUE)
    logcpm[1:5,]
    cpm <- cpm(y, prior.count=2, log=FALSE)
    cpm[1:5, ]
    # fold_change <- logcpm[,'Ratio1']/logcpm[,'Ratio2']
    fold_change <- cpm[,1]/cpm[,2]
    resMeR <- countMatrix
    resMeR$logFC <- log2(fold_change)
    resMeR$FoldChange <- fold_change
    resMeR$change <- with(resMeR, ifelse(abs(resMeR$logFC) > log2(1.5),
                                         ifelse(resMeR$logFC > log2(1.5),'UP','DOWN'),'NOT'))
    table(resMeR$change)
    resMeR = as.data.frame(resMeR[resMeR$change %in% c('UP', 'DOWN'), ])
    resMeR_list[[c]] <- list(resMeR)
    countMatrix_list[[c]] <- countMatrix
    resPeak_list[[c]] <- merge_peak
    name = paste0(c, ' vs WT')
    name = paste0(gsub(' ','_', name ))
    names(resMeR_list[[c]]) <- paste0(name, '_', c('limma'))
    
  } else{
    # 2. Exploring Method 1: DESeq2-------------------------------------------
    # 构建 DESeqDataSet 对象
    dds <- DESeqDataSetFromMatrix(countData = countMatrix,
                                  colData = colData,
                                  design = ~ groupList)
    ## Pre-filtering
    keep <- rowSums(counts(dds)) >= 1
    table(keep)
    dds <- dds[keep,]
    # 接着计算每个样本的 library size 用于后续的标准化。library size 等于落在该样本所有peaks上的reads的总数，即 counts 矩阵中每列的加和。如果bFullLibrarySize设为TRUE，则会使用library的总reads数（根据原BAM/BED文件统计获得）。然后使用estimateDispersions估计统计分布，需要将参数fitType设为’local’
    dds<-estimateSizeFactors(dds)
    dds<-estimateDispersions(dds,fitType="local")
    dds <- nbinomWaldTest(dds)
    res <- results(dds)
    ## 提取你想要的差异分析结果，我们这里是treated组对untreated组进行比较
    res <-  results(dds, contrast=c("groupList", unique(condition)))
    res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
    length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
    res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                                   ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
    res_DESeq2 = as.data.frame(res)
    table(res_DESeq2$change)
    res_DESeq2 = res_DESeq2[res_DESeq2$change %in% c('UP', 'DOWN'), ]
    # 2. Exploring Method 2: limma-------------------------------------------
    library(edgeR)
    d0 <- DGEList(counts = countMatrix)
    dim(d0)
    # Filtering low count peaks
    cutoff <- 1           
    drop <- which(apply(cpm(d0), 1, max) < cutoff)
    print(length(drop))
    if (length(drop) >0) {
      d <- d0[-drop,] 
    } else{ d <- d0}
    dim(d) #number of peaks left: 8501      
    d <- calcNormFactors(d, method = c('TMM'))
    colnames(countMatrix)
    genotype <- groupList
    table(genotype)
    mm <- model.matrix(~0 + genotype)#
    v <- voom(d, mm, plot=F, normalize="quantile")#voom from limma:Transform count data to log2-counts per million (logtpm)
    fit <- lmFit(v, mm)
    head(coef(fit))
    contr <- makeContrasts(paste0(colnames(coef(fit)), collapse  = ' - '), levels = colnames(coef(fit))) #Change if other comparison
    tmp <- contrasts.fit(fit, contr)
    tmp <- eBayes(tmp)
    top.table <- topTable(tmp, sort.by = "P", n = Inf)
    head(top.table, 20)
    length(which(top.table$adj.P.Val < 0.05 & abs(top.table$logFC) > log2(1.5))) #how many are significant?
    res_limma <- as.data.frame(top.table)
    res_limma$change <- with(res_limma, ifelse(res_limma$adj.P.Val < 0.05 & abs(res_limma$logFC) > log2(1.5),
                                               ifelse(res_limma$logFC > log2(1.5),'UP','DOWN'),'NOT'))
    table(res_limma$change)
    res_limma = res_limma[res_limma$change %in% c('UP', 'DOWN'), ]
    resMeR_list[[c]] <- list(res_DESeq2, res_limma)
    countMatrix_list[[c]] <- countMatrix
    resPeak_list[[c]] <- merge_peak
    name = paste0(c, ' vs WT')
    name = paste0(gsub(' ','_', name ))
    names(resMeR_list[[c]] ) <- paste0(name, '_', c('DESeq2', 'limma'))
    
  }
}
names(resMeR_list) <- NULL
RMdem_list[[i]] <- as.list(unlist(resMeR_list, recursive = FALSE))
names(RMdem_list[[i]] ) <- paste0(cl, '_', names(RMdem_list[[i]] ))
RMdep_list[[i]] <- resPeak_list 
name <- unique(gsub('(_limma|_DESeq2)', '', names(RMdem_list[[i]] )))
names(RMdep_list[[i]] ) <- name
contrast_IP_list[[i]] <- lapply(countMatrix_list, extract_names)
condition_IP_list[[i]] <- condition_list

############
n = 66
############
#GSE137675# without case RIP, with Input and RNA
############
print(RMPdeg_datasets[n,])#
i = RMPdeg_datasets[n,'GSE']
cl = c('CM2005.1 and CRMM1 and OCM1 and OCM1a and OM431 and PIG1')#Ocular melanoma cell
cl = c('OCM1', 'PIG1')
RMp = c('ALKBH5', 'METTL3')
RM = 'm6A'
condition = c(rep('ALKBH5KD', 1), rep('WT', 1))
condition = c(rep('METTL3KD', 1), rep('WT', 1))
RMP_allDat[grep('GSE137675', RMP_allDat$study_alias )
           ,c("study_alias", "sample_alias", "sample_title")]
# View(allDat[[i]])
RMP_allDat[grepl('GSE137675', RMP_allDat$study_alias )
           ,c("study_alias", "sample_alias", "sample_title")]
files <- list.files(file.path(data_directory,i))
files
getwd()
setwd(data_directory)
setwd(i)
###################################
#   more than 1 condition here   #
###################################
resInP_list <- list()
##Gene count matrix
# DEG <- read_xlsx('GSE137675_OCM1.sgALKBH5.and.OCM1.sgNC.RNA-seq.processed.data.xlsx', sheet = 1)#FPKM
# countGTF <- as.data.frame(DEG[-(1:8),c(2,6,7)])
# colnames(countGTF) <- DEG[8,c(2,6,7)]
# table(duplicated(countGTF$Gene_Name))
# countGTF <- countGTF[!duplicated(countGTF$Gene_Name), ]
# rownames(countGTF) <- countGTF$Gene_Name
# countGTF <-  fpkm_to_counts(countGTF[, -1], mean_rps=400000) 
DEG <- read_xls('GSE137675_PIG1.shMETTL3.and.PIG1.shNC.RNA-seq.processed.data.xls', sheet = 1)
countGTF <- as.data.frame(DEG[,c(2,4,7)])
table(duplicated(countGTF$`Gene Name`))
countGTF <- countGTF[!duplicated(countGTF$`Gene Name`), ]
rownames(countGTF) <- countGTF$`Gene Name`
countGTF <- countGTF[, -1]
# 2. Differential gene analysis with edgeR
# Without replicates, you cannot estimate which genes are differentially expressed using EdgeR or DESeq2. You can only calculate fold changes based on normalized read counts (preferentially CPM normalized by TMM method included in any EdgeR analysis) and apply a stringent fold-change cut-off to determine which genes are more or less expressed depending on the condition.
y <- DGEList(countGTF)
# Calculate fold change
y <- calcNormFactors(y, method = c('TMM'))
y$samples$lib.size <- colSums(y$counts)
logcpm <- cpm(y, prior.count=2, log=TRUE)
logcpm[1:5,]
cpm <- cpm(y, prior.count=2, log=FALSE)
cpm[1:5, ]
fold_change <- cpm[,1]/cpm[,2]
res <- countGTF
res$logFC <- log2(fold_change)
res$FoldChange <- fold_change
res$change <- with(res, ifelse(abs(res$logFC) > log2(1.5),
                               ifelse(res$logFC > log2(1.5),'UP','DOWN'),'NOT'))
table(res$change)
resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
resInP$geneSymbol <- rownames(resInP)
table(resInP$change)
table(resInP$FoldChange %in% c(0, Inf))
resInP_list[['PIG1_METTL3KD_vs_WT']] <- resInP

RMdeg_list[[i]] <- resInP_list
c <- unique(condition[!condition %in% 'WT'])
contrast_input_list[[i]] <- NULL#lapply(countGTF_list, extract_names)
condition_list <- list(condition); names(condition_list) <- c
condition_input_list[[i]] <- condition_list

############
n = 67
############
#GSE137752# without case RIP, with case RNA
############
print(RMPdeg_datasets[n,])#
i = RMPdeg_datasets[n,'GSE']
cl = c('Hela and K562 and MCF7 and HMC3')#CML
cl = c('K562', 'hNPC')
RMp = c('SAFB2')
RM = 'm6A'
condition_list = list(hNPC_SAFB2KD = c(rep('SAFB2KD', 2), rep('WT', 2)),
                      K562_SAFB2KD = c(rep('SAFB2KD', 2), rep('WT', 1))
)
RMP_allDat[grep('GSE137752', RMP_allDat$study_alias )
           ,c("study_alias", "sample_alias", "sample_title")]
# View(allDat[[i]])
files <- list.files(file.path(data_directory,i))
files
getwd()
setwd(data_directory)
setwd(i)
###################################
#         RNA-seq analysis        #
###################################
resInP_list <- vector("list", length(condition_list))
names(resInP_list) <- names(condition_list)
countGTF_list <- vector("list", length(condition_list))
names(countGTF_list) <- names(condition_list)
for(c in names(condition_list)){
  countGTF <- read.table( paste0("metaPlotR/Input_counts_", c, ".txt"),header=T)
  dim(countGTF)
  colnames(countGTF)
  rownames(countGTF) <- countGTF$Geneid
  countGTF <- countGTF[, !colnames(countGTF) %in% c('Geneid','Chr','Start','End','Strand','Length')]
  dim(countGTF) # 28265     4
  condition <- condition_list[[c]]
  ## Sample information
  sample=data.frame(
    sample=colnames(countGTF),
    condition=condition)
  ## SummarizedExperiment input
  dds<-DESeqDataSetFromMatrix(countData=countGTF,colData=sample,design=~condition)
  ## Pre-filtering
  keep <- rowSums(counts(dds)) >= 10
  table(keep)
  dds <- dds[keep,]
  ## Note on factor levels
  dds$condition <- factor(dds$condition, levels = unique(condition))
  dds$condition <- relevel(dds$condition, ref = "WT")
  # Differential expression analysis-------------------------------------------
  dds <- DESeq(dds)
  res <-  results(dds, contrast=c("condition",unique(condition)))
  res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
  length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
  res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                                 ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
  table(res$change)
  table(is.na(res$change))
  # resOrdered <- as.data.frame(res[order(res$padj), ])
  resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
  resInP$geneSymbol <- rownames(resInP)
  table(resInP$change)
  resInP_list[[c]] <- resInP
  countGTF_list[[c]] <- countGTF 
}
name = paste0(names(resInP_list), ' vs WT')
RMdeg_list[[i]] <- resInP_list 
names(RMdeg_list[[i]] ) <- paste0(gsub(' ','_',(name)))
contrast_input_list[[i]] <- lapply(countGTF_list, extract_names)
condition_input_list[[i]] <- condition_list

############
n = 68
############
#GSE129469# with RIP and Input
############
print(RMPdeg_datasets[n,])
i = RMPdeg_datasets[n,'GSE']
cl = c('SU-DHL-8')#DLBCL:Diffuse large B-cell lymphoma
RMp = c('WTAP')
RM = 'm6A'
condition = c(rep('WTAPKD', 1), rep('WT', 1))
RMP_allDat[grep('GSE129469', RMP_allDat$study_alias )
           ,c("study_alias", "sample_alias", "sample_title")]
# View(allDat[[i]])
ENA <- as.data.frame(RMP_allDat[grep('GSE129469', RMP_allDat$study_alias )
                                , c("study_alias", "sample_alias", "sample_title")])
ENA <- ENA[order(ENA$sample_alias),]
ENA

files <- list.files(file.path(data_directory,i))
files
getwd()
setwd(data_directory)
setwd(i)
###################################
#  Cell number-normalized (CNN)   #
#      MeRIP INPUT analysis       #
###################################
c <- unique(condition[!condition %in% 'WT'])
countGTF <- read.table( paste0("metaPlotR/Input_counts_", c, ".txt"),header=T)
dim(countGTF)
colnames(countGTF)
rownames(countGTF) <- countGTF$Geneid
countGTF <- countGTF[, !colnames(countGTF) %in% c('Geneid','Chr','Start','End','Strand','Length')]
dim(countGTF) # 28265     4
# 2. Differential gene analysis with edgeR
# Without replicates, you cannot estimate which genes are differentially expressed using EdgeR or DESeq2. You can only calculate fold changes based on normalized read counts (preferentially CPM normalized by TMM method included in any EdgeR analysis) and apply a stringent fold-change cut-off to determine which genes are more or less expressed depending on the condition.
y <- DGEList(countGTF)
# Calculate fold change
y <- calcNormFactors(y, method = c('TMM'))
y$samples$lib.size <- colSums(y$counts)
logcpm <- cpm(y, prior.count=2, log=TRUE)
logcpm[1:5,]
cpm <- cpm(y, prior.count=2, log=FALSE)
cpm[1:5, ]
fold_change <- cpm[,1]/cpm[,2]
res <- countGTF
res$logFC <- log2(fold_change)
res$FoldChange <- fold_change
res$change <- with(res, ifelse(abs(res$logFC) > log2(1.5),
                               ifelse(res$logFC > log2(1.5),'UP','DOWN'),'NOT'))
table(res$change)
resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
resInP$geneSymbol <- rownames(resInP)
table(resInP$change)
table(resInP$FoldChange %in% c(0, Inf))
# resInP <- resInP[! resInP$FoldChange %in% c(0, Inf), ]
name = paste0(c, ' vs WT')
name = paste0(cl, '_', name)
RMdeg_list[[i]] <- list( resInP)
names(RMdeg_list[[i]] ) <- paste0(gsub(' ','_',(name)))

countGTF_list <- list(colnames(countGTF)); names(countGTF_list) <- c
condition_list <- list(condition); names(condition_list) <- c
contrast_input_list[[i]] <- countGTF_list
condition_input_list[[i]] <- condition_list

###################################
#  Cell number-normalized (CNN)   #
#   MeRIP RIP analysis - Part1    #
###################################
## Sample information
merge_peak <- read.delim(paste0("metaPlotR/", c, "_newPeakFile.txt"),header=F)
countMatrix <- read.table(paste0("metaPlotR/m6A_counts_", c, ".txt"),header=T)
table(countMatrix$Geneid == merge_peak[-1,1])
table(duplicated(merge_peak$V1[-1]))
dim(countMatrix)
colnames(countMatrix)
rownames(countMatrix) <- countMatrix$Geneid
countMatrix <- countMatrix[, !colnames(countMatrix) %in% c('Geneid','Chr','Start','End','Strand','Length')]
# 2. Differential gene analysis with edgeR
# Without replicates, you cannot estimate which genes are differentially expressed using EdgeR or DESeq2. You can only calculate fold changes based on normalized read counts (preferentially CPM normalized by TMM method included in any EdgeR analysis) and apply a stringent fold-change cut-off to determine which genes are more or less expressed depending on the condition.
y <- DGEList(countMatrix)
# Calculate fold change
y <- calcNormFactors(y, method = c('TMM'))
y$samples$lib.size <- colSums(y$counts)
logcpm <- cpm(y, prior.count=2, log=TRUE)
logcpm[1:5,]
cpm <- cpm(y, prior.count=2, log=FALSE)
cpm[1:5, ]
# fold_change <- logcpm[,'Ratio1']/logcpm[,'Ratio2']
fold_change <- cpm[,1]/cpm[,2]
resMeR <- countMatrix
resMeR$logFC <- log2(fold_change)
resMeR$FoldChange <- fold_change
resMeR$change <- with(resMeR, ifelse(abs(resMeR$logFC) > log2(1.5),
                                     ifelse(resMeR$logFC > log2(1.5),'UP','DOWN'),'NOT'))
table(resMeR$change)
res_limma <- resMeR[resMeR$change %in% c('UP','DOWN'),]
RMdem_list[[i]] <- list(res_limma)
name = paste0(c, ' vs WT')
name = paste0(gsub(' ','_', name ))
names(RMdem_list[[i]]) <- paste0(cl, '_', name, '_limma')
RMdep_list[[i]] <- list( merge_peak )
names(RMdep_list[[i]] ) <-  paste0(cl, '_', name)

countMatrix_list <- list(colnames(countMatrix)); names(countMatrix_list) <- c
condition_list <- list(condition); names(condition_list) <- c
contrast_IP_list[[i]] <- countMatrix_list
condition_IP_list[[i]] <- condition_list

############
n = 69
############
#GSE171472# with RIP and Input
############
print(RMPdeg_datasets[n,])#
i = RMPdeg_datasets[n,'GSE']
cl = c('HEK293T and H322 and H358')#CML
cl = 'H322'
RMp = c('FTO')
RM = 'm6A'
condition = c(rep('FTOKD', 1), rep('WT', 1))
RMP_allDat[grep('GSE171472', RMP_allDat$study_alias )
           ,c("study_alias", "sample_alias", "sample_title")]
# View(allDat[[i]])
ENA <- as.data.frame(RMP_allDat[grep('GSE171472', RMP_allDat$study_alias )
                                , c("study_alias", "sample_alias", "sample_title")])
ENA <- ENA[order(ENA$sample_alias),]
ENA
files <- list.files(file.path(data_directory,i))
files
getwd()
setwd(data_directory)
setwd(i)

###################################
#  Cell number-normalized (CNN)   #
#      MeRIP INPUT analysis       #
###################################
c <- unique(condition[!condition %in% 'WT'])
countGTF <- read.table( paste0("metaPlotR/Input_counts_", c, ".txt"),header=T)
dim(countGTF)
colnames(countGTF)
rownames(countGTF) <- countGTF$Geneid
countGTF <- countGTF[, !colnames(countGTF) %in% c('Geneid','Chr','Start','End','Strand','Length')]
dim(countGTF) # 28265     4
# 2. Differential gene analysis with edgeR
# Without replicates, you cannot estimate which genes are differentially expressed using EdgeR or DESeq2. You can only calculate fold changes based on normalized read counts (preferentially CPM normalized by TMM method included in any EdgeR analysis) and apply a stringent fold-change cut-off to determine which genes are more or less expressed depending on the condition.
y <- DGEList(countGTF)
# Calculate fold change
y <- calcNormFactors(y, method = c('TMM'))
y$samples$lib.size <- colSums(y$counts)
logcpm <- cpm(y, prior.count=2, log=TRUE)
logcpm[1:5,]
cpm <- cpm(y, prior.count=2, log=FALSE)
cpm[1:5, ]
fold_change <- cpm[,1]/cpm[,2]
res <- countGTF
res$logFC <- log2(fold_change)
res$FoldChange <- fold_change
res$change <- with(res, ifelse(abs(res$logFC) > log2(1.5),
                               ifelse(res$logFC > log2(1.5),'UP','DOWN'),'NOT'))
table(res$change)
resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
resInP$geneSymbol <- rownames(resInP)
table(resInP$change)
table(resInP$FoldChange %in% c(0, Inf))
# resInP <- resInP[! resInP$FoldChange %in% c(0, Inf), ]
name = paste0(c, ' vs WT')
name = paste0(cl, '_', name)
RMdeg_list[[i]] <- list( resInP)
names(RMdeg_list[[i]] ) <- paste0(gsub(' ','_',(name)))

countGTF_list <- list(colnames(countGTF)); names(countGTF_list) <- c
condition_list <- list(condition); names(condition_list) <- c
contrast_input_list[[i]] <- countGTF_list
condition_input_list[[i]] <- condition_list

###################################
#  Cell number-normalized (CNN)   #
#   MeRIP RIP analysis - Part1    #
###################################
## Sample information
merge_peak <- read.delim(paste0("metaPlotR/", c, "_newPeakFile.txt"),header=F)
countMatrix <- read.table(paste0("metaPlotR/m6A_counts_", c, ".txt"),header=T)
table(countMatrix$Geneid == merge_peak[-1,1])
table(duplicated(merge_peak$V1[-1]))
dim(countMatrix)
colnames(countMatrix)
rownames(countMatrix) <- countMatrix$Geneid
countMatrix <- countMatrix[, !colnames(countMatrix) %in% c('Geneid','Chr','Start','End','Strand','Length')]
# 2. Differential gene analysis with edgeR
# Without replicates, you cannot estimate which genes are differentially expressed using EdgeR or DESeq2. You can only calculate fold changes based on normalized read counts (preferentially CPM normalized by TMM method included in any EdgeR analysis) and apply a stringent fold-change cut-off to determine which genes are more or less expressed depending on the condition.
y <- DGEList(countMatrix)
# Calculate fold change
y <- calcNormFactors(y, method = c('TMM'))
y$samples$lib.size <- colSums(y$counts)
logcpm <- cpm(y, prior.count=2, log=TRUE)
logcpm[1:5,]
cpm <- cpm(y, prior.count=2, log=FALSE)
cpm[1:5, ]
# fold_change <- logcpm[,'Ratio1']/logcpm[,'Ratio2']
fold_change <- cpm[,1]/cpm[,2]
resMeR <- countMatrix
resMeR$logFC <- log2(fold_change)
resMeR$FoldChange <- fold_change
resMeR$change <- with(resMeR, ifelse(abs(resMeR$logFC) > log2(1.5),
                                     ifelse(resMeR$logFC > log2(1.5),'UP','DOWN'),'NOT'))
table(resMeR$change)
res_limma <- resMeR[resMeR$change %in% c('UP','DOWN'),]
RMdem_list[[i]] <- list(res_limma)
name = paste0(c, ' vs WT')
name = paste0(gsub(' ','_', name ))
names(RMdem_list[[i]]) <- paste0(cl, '_', name, '_limma')
RMdep_list[[i]] <- list( merge_peak )
names(RMdep_list[[i]] ) <-  paste0(cl, '_', name)

countMatrix_list <- list(colnames(countMatrix)); names(countMatrix_list) <- c
condition_list <- list(condition); names(condition_list) <- c
contrast_IP_list[[i]] <- countMatrix_list
condition_IP_list[[i]] <- condition_list

# ############
# n = 70
# ############
# #GSE168778# without case RIP and Input
# ############
# print(RMPdeg_datasets[n,]) # skip?
# i = RMPdeg_datasets[n,'GSE']
# cl = c('HEK293T and H322 and H358')#CML
# RMp = c('FTO')
# RM = 'm6A'
# condition = c(rep('FTOKD', 1), rep('WT', 1))
# RMP_allDat[grep('GSE168778', RMP_allDat$study_alias )
#            ,c("study_alias", "sample_alias", "sample_title")]
# # View(allDat[[i]])
# files <- list.files(file.path(data_directory,i))
# files
# getwd()
# setwd(data_directory)
# setwd(i)
# 
# 
############
n = 71
############
#GSE154555# with RIP and Input
############
print(RMPdeg_datasets[n,])#
i = RMPdeg_datasets[n,'GSE']
cl = c('KYSE180 and KYSE450')#Human ESCC
cl = 'K180'
RMp = c('METTL3')
RM = 'm6A'
condition = c(rep('METTL3KD', 1), rep('WT', 1))
RMP_allDat[grep('GSE154555', RMP_allDat$study_alias )
           ,c("study_alias", "sample_alias", "sample_title")]
# View(allDat[[i]])
skip <- RMP_allDat[grepl('GSE154555', RMP_allDat$study_alias )
                   &  grepl('Input', RMP_allDat$sample_title)
                   , "sample_alias"]
paste0(skip$sample_alias, collapse = " ")
skip <- RMP_allDat[grepl('GSE154555', RMP_allDat$study_alias )
                   &  grepl('IP', RMP_allDat$sample_title)
                   , "sample_alias"]
paste0(skip$sample_alias, collapse = " ")
ENA <- as.data.frame(RMP_allDat[grep('GSE154555', RMP_allDat$study_alias )
                                , c("study_alias", "sample_alias", "sample_title")])
ENA <- ENA[order(ENA$sample_alias),]
ENA
files <- list.files(file.path(data_directory,i))
files
getwd()
setwd(data_directory)
setwd(i)
# DEP <- read_excel(file.path(data_directory,i, 'GSE154555_Methylated_RNA_sites.mRNA.xlsx'),sheet = 2)
###################################
#  Cell number-normalized (CNN)   #
#      MeRIP INPUT analysis       #
###################################
c <- unique(condition[!condition %in% 'WT'])
countGTF <- read.table( paste0("metaPlotR/Input_counts_", c, ".txt"),header=T)
dim(countGTF)
colnames(countGTF)
rownames(countGTF) <- countGTF$Geneid
countGTF <- countGTF[, !colnames(countGTF) %in% c('Geneid','Chr','Start','End','Strand','Length')]
dim(countGTF) # 28265     4
# 2. Differential gene analysis with edgeR
# Without replicates, you cannot estimate which genes are differentially expressed using EdgeR or DESeq2. You can only calculate fold changes based on normalized read counts (preferentially CPM normalized by TMM method included in any EdgeR analysis) and apply a stringent fold-change cut-off to determine which genes are more or less expressed depending on the condition.
y <- DGEList(countGTF)
# Calculate fold change
y <- calcNormFactors(y, method = c('TMM'))
y$samples$lib.size <- colSums(y$counts)
logcpm <- cpm(y, prior.count=2, log=TRUE)
logcpm[1:5,]
cpm <- cpm(y, prior.count=2, log=FALSE)
cpm[1:5, ]
fold_change <- cpm[,1]/cpm[,2]
res <- countGTF
res$logFC <- log2(fold_change)
res$FoldChange <- fold_change
res$change <- with(res, ifelse(abs(res$logFC) > log2(1.5),
                               ifelse(res$logFC > log2(1.5),'UP','DOWN'),'NOT'))
table(res$change)
resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
resInP$geneSymbol <- rownames(resInP)
table(resInP$change)
table(resInP$FoldChange %in% c(0, Inf))
# resInP <- resInP[! resInP$FoldChange %in% c(0, Inf), ]
name = paste0(c, ' vs WT')
name = paste0(cl, '_', name)
RMdeg_list[[i]] <- list( resInP)
names(RMdeg_list[[i]] ) <- paste0(gsub(' ','_',(name)))

countGTF_list <- list(colnames(countGTF)); names(countGTF_list) <- c
condition_list <- list(condition); names(condition_list) <- c
contrast_input_list[[i]] <- countGTF_list
condition_input_list[[i]] <- condition_list

###################################
#  Cell number-normalized (CNN)   #
#   MeRIP RIP analysis - Part1    #
###################################
## Sample information
merge_peak <- read.delim(paste0("metaPlotR/", c, "_newPeakFile.txt"),header=F)
countMatrix <- read.table(paste0("metaPlotR/m6A_counts_", c, ".txt"),header=T)
table(countMatrix$Geneid == merge_peak[-1,1])
table(duplicated(merge_peak$V1[-1]))
dim(countMatrix)
colnames(countMatrix)
rownames(countMatrix) <- countMatrix$Geneid
countMatrix <- countMatrix[, !colnames(countMatrix) %in% c('Geneid','Chr','Start','End','Strand','Length')]
# 2. Differential gene analysis with edgeR
# Without replicates, you cannot estimate which genes are differentially expressed using EdgeR or DESeq2. You can only calculate fold changes based on normalized read counts (preferentially CPM normalized by TMM method included in any EdgeR analysis) and apply a stringent fold-change cut-off to determine which genes are more or less expressed depending on the condition.
y <- DGEList(countMatrix)
# Calculate fold change
y <- calcNormFactors(y, method = c('TMM'))
y$samples$lib.size <- colSums(y$counts)
logcpm <- cpm(y, prior.count=2, log=TRUE)
logcpm[1:5,]
cpm <- cpm(y, prior.count=2, log=FALSE)
cpm[1:5, ]
# fold_change <- logcpm[,'Ratio1']/logcpm[,'Ratio2']
fold_change <- cpm[,1]/cpm[,2]
resMeR <- countMatrix
resMeR$logFC <- log2(fold_change)
resMeR$FoldChange <- fold_change
resMeR$change <- with(resMeR, ifelse(abs(resMeR$logFC) > log2(1.5),
                                     ifelse(resMeR$logFC > log2(1.5),'UP','DOWN'),'NOT'))
table(resMeR$change)
res_limma <- resMeR[resMeR$change %in% c('UP','DOWN'),]
RMdem_list[[i]] <- list(res_limma)
name = paste0(c, ' vs WT')
name = paste0(gsub(' ','_', name ))
names(RMdem_list[[i]]) <- paste0(cl, '_', name, '_limma')
RMdep_list[[i]] <- list( merge_peak )
names(RMdep_list[[i]] ) <-  paste0(cl, '_', name)

countMatrix_list <- list(colnames(countMatrix)); names(countMatrix_list) <- c
condition_list <- list(condition); names(condition_list) <- c
contrast_IP_list[[i]] <- countMatrix_list
condition_IP_list[[i]] <- condition_list


# ############
# n = 72
# ############
# #GSE158742# without case RIP and Input
# ############
# print(RMPdeg_datasets[n,])# skip?
# i = RMPdeg_datasets[n,'GSE']
# cl = c('GSC387')#Glioblastoma
# RMp = c('YTHDF2')
# RM = 'm6A'
# # setwd(data_directory)
# # gse <- getGEO('GSE158020',destdir=".", AnnotGPL=F, getGPL=F)
# # setwd('GSE158742')
# # pdata <- pData(gse[[1]])
# # condition = c(rep('METTL3KD', 1), rep('WT', 1))
# RMP_allDat[grep('GSE158020', RMP_allDat$study_alias )
#            ,c("study_alias", "sample_alias", "sample_title")]
# # View(allDat[[i]])
# paste0(skip$sample_alias, collapse = " ")
# files <- list.files(file.path(data_directory,i))
# files
# getwd()
# setwd(data_directory)
# setwd(i)
# 
# 
############
n = 73
############
#GSE155413# without case RIP, with case RNA
############
print(RMPdeg_datasets[n,])#
i = RMPdeg_datasets[n,'GSE']
cl = c('Huh7')#Liver cancer
RMp = c('METTL3', 'METTL14')
RM = 'm6A'
condition_list <- list(Huh7.IFN_METTL3_14KD = c(rep('METTL3_14KD', 3), rep('WT', 3)),
                         Huh7_METTL3_14KD = c(rep('METTL3_14KD', 3), rep('WT', 3))
)
RMP_allDat[grep('GSE155413', RMP_allDat$study_alias )
           ,c("study_alias", "sample_alias", "sample_title")]
# View(allDat[[i]])
paste0(skip$sample_alias, collapse = " ")
files <- list.files(file.path(data_directory,i))
files
getwd()
setwd(data_directory)
setwd(i)
# 1) RNA-seq
count_files <- files[grepl('_countMatrix.txt.gz', files)]
###################################
# Read and process each count file
count <- lapply(count_files[!grepl('RIBO', count_files)], function(x) {
  # Read data using fread for faster reading
  data <- fread(file.path(data_directory, i, x))
  # Extract required columns
  df_new <- data[, c(1, 2:7)]
  # Set column names same as file names
  colnames(df_new)[1] <- c('geneSymbol')
  return(df_new)
})
# Merge dataframes in countGTF list by 'gene_id'
df_in <- Reduce(function(x, y) merge(x, y, by = 'geneSymbol', all = TRUE), count)
colSums(df_in[,-1])
colnames(df_in) <- str_split(colnames(df_in), '[.]', simplify = T)[,1]
df_in <- as.data.frame(df_in)
index = match(unique(colnames(df_in)), colnames(df_in))
df_in <- df_in[, index]
table(colSums(df_in[,-1])=="1e+06")
## Gene name
library(data.table)
table(duplicated(df_in$geneSymbol))
countGTF=df_in[, -1]
rownames(countGTF) <- df_in$geneSymbol
###################################
#   more than 1 condition here   #
###################################
resInP_list <- list()
##Gene count matrix
# countGTF=as.matrix(ceiling(df_in[, -1]))
## Gene name
table(duplicated(df_in$geneSymbol))
dim(countGTF) #25343    12
## Sample name
colnames(countGTF)
## Sample information
sample_meta=data.frame(
  sample=colnames(countGTF),
  condition=c(rep("Huh7.IFN_WT",3), rep("Huh7.IFN_METTL3_14KD",3), 
              rep("Huh7_WT",3), rep("Huh7_METTL3_14KD",3)))
#####condition_Huh7.IFN_METTL314_vs_WT#####
count <- as.matrix(countGTF[, grepl('.IFN', sample_meta$condition)])
table(is.na(count))
sample <- sample_meta[sample_meta$sample %in% colnames(count), ]
## SummarizedExperiment input
condition = c(rep('Huh7.IFN_METTL3_14KD', 3), rep('Huh7.IFN_WT', 3))
dds<-DESeqDataSetFromMatrix(countData=count, colData=sample,design=~condition)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 10
table(keep)
dds <- dds[keep,]
## Note on factor levels
dds$condition <- factor(dds$condition, levels = c("Huh7.IFN_METTL3_14KD","Huh7.IFN_WT"))
dds$condition <- relevel(dds$condition, ref = "Huh7.IFN_WT")
dds <- DESeq(dds)
res <-  results(dds, contrast=c("condition",unique(condition))) # c("condition", rev(levels(genotype)))
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
table(res$change)
res$geneSymbol <- rownames(countGTF)[keep]
resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
table(duplicated(resInP$geneSymbol))
resInP_list[['Huh7.IFN_METTL3_14KD_vs_WT']] <- resInP
#####condition_K562_shNSUN2_vs_WT#####
count <- as.matrix(countGTF[, grepl('Huh7_', sample_meta$condition)])
table(is.na(count))
sample <- sample_meta[sample_meta$sample %in% colnames(count), ]
## SummarizedExperiment input
condition = c(rep('Huh7_METTL3_14KD', 3), rep('Huh7_WT', 3))
dds<-DESeqDataSetFromMatrix(countData=count, colData=sample,design=~condition)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 10
table(keep)
dds <- dds[keep,]
## Note on factor levels
dds$condition <- factor(dds$condition, levels = c("Huh7_METTL3_14KD","Huh7_WT"))
dds$condition <- relevel(dds$condition, ref = "Huh7_WT")
dds <- DESeq(dds)
res <-  results(dds, contrast=c("condition",unique(condition))) # c("condition", rev(levels(genotype)))
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
table(res$change)
res$geneSymbol <- rownames(countGTF)[keep]
resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
table(duplicated(resInP$geneSymbol))
resInP_list[['Huh7_METTL3_14KD_vs_WT']] <- resInP


RMdeg_list[[i]] <- resInP_list
contrast_input_list[[i]] <- NULL#lapply(countGTF_list, extract_names)
condition_input_list[[i]] <- condition_list


############
n = 74
############
#GSE148962# with RIP and Input
############
print(RMPdeg_datasets[n,])#
i = RMPdeg_datasets[n,'GSE']
cl = c('HTR8.Svneo')#Embryo placenta
RMp = c('METTL3')
RM = 'm6A'
condition = c(rep('METTL3KD', 2), rep('WT', 2))
RMP_allDat[grep('GSE148962', RMP_allDat$study_alias )
           ,c("study_alias", "sample_alias", "sample_title")]
# View(allDat[[i]])
RMP_allDat[grepl('GSE148962', RMP_allDat$study_alias ) & (! grepl('RNC', RMP_allDat$sample_title))
           ,c("study_alias", "sample_alias", "sample_title")]
skip = RMP_allDat[grepl('GSE148962', RMP_allDat$study_alias ) & (! grepl('RNC', RMP_allDat$sample_title))
                  ,c("study_alias", "sample_alias", "sample_title")]
paste0(skip$sample_alias, collapse = "|")

skip <- RMP_allDat[grepl('GSE148962', RMP_allDat$study_alias ) & (! grepl('RNC', RMP_allDat$sample_title))
                   &  grepl('Input', RMP_allDat$sample_title)
                   , "sample_alias"]
paste0(skip$sample_alias, collapse = " ")
skip <- RMP_allDat[grepl('GSE148962', RMP_allDat$study_alias )  & (! grepl('RNC', RMP_allDat$sample_title))
                   &  grepl('IP', RMP_allDat$sample_title)
                   , "sample_alias"]
paste0(skip$sample_alias, collapse = " ")
ENA <- as.data.frame(RMP_allDat[grep('GSE148962', RMP_allDat$study_alias )
                                , c("study_alias", "sample_alias", "sample_title")])
ENA <- ENA[order(ENA$sample_alias),]
ENA
files <- list.files(file.path(data_directory,i))
files
getwd()
setwd(data_directory)
setwd(i)

c <- unique(condition[! condition %in% 'WT'])
###################################
#  Cell number-normalized (CNN)   #
#      MeRIP INPUT analysis       #
###################################
countGTF <- read.table( paste0("metaPlotR/Input_counts_", c, ".txt"),header=T)
rownames(countGTF) <- countGTF$Geneid
countGTF <- countGTF[, !colnames(countGTF) %in% c('Geneid','Chr','Start','End','Strand','Length')]
dim(countGTF) #28265     4
## Sample name
## Sample information
sample=data.frame(
  sample=colnames(countGTF),
  condition=condition)
## SummarizedExperiment input
dds<-DESeqDataSetFromMatrix(countData=countGTF,colData=sample,design=~condition)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 10
table(keep)
dds <- dds[keep,]
## Note on factor levels
dds$condition <- factor(dds$condition, levels = unique(condition))
dds$condition <- relevel(dds$condition, ref = "WT")
# Differential expression analysis-------------------------------------------
dds <- DESeq(dds)

res <-  results(dds, contrast=c("condition",unique(condition)))
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
table(res$change)
table(is.na(res$change))
# resOrdered <- as.data.frame(res[order(res$padj), ])
resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
resInP$geneSymbol <- rownames(resInP)
name = res@elementMetadata$description[2]
name = paste0(cl, '_',str_split(name, "condition ", simplify = T)[,2])
RMdeg_list[[i]] <- list( resInP)
names(RMdeg_list[[i]] ) <- paste0(gsub(' ','_',(name)))
c <- unique(condition[!condition %in% 'WT'])
countGTF_list <- list(colnames(countGTF)); names(countGTF_list) <- c
condition_list <- list(condition); names(condition_list) <- c
contrast_input_list[[i]] <- countGTF_list
condition_input_list[[i]] <- condition_list

###################################
#  Cell number-normalized (CNN)   #
#   MeRIP RIP analysis - Part1    #
###################################
## Sample information
merge_peak <- read.delim(paste0("metaPlotR/", c, "_newPeakFile.txt"),header=F)
countMatrix <- read.table(paste0("metaPlotR/m6A_counts_", c, ".txt"),header=T)
table(countMatrix$Geneid == merge_peak[-1,1])
table(duplicated(merge_peak$V1[-1]))
dim(countMatrix)
colnames(countMatrix)
rownames(countMatrix) <- countMatrix$Geneid
countMatrix <- countMatrix[, !colnames(countMatrix) %in% c('Geneid','Chr','Start','End','Strand','Length')]
# 1. 做好分组因子即可
groupList<-factor(condition)
colData <- data.frame(row.names=colnames(countMatrix), groupList=groupList)
# 2. Exploring Method 1: DESeq2-------------------------------------------
# 构建 DESeqDataSet 对象
dds <- DESeqDataSetFromMatrix(countData = countMatrix,
                              colData = colData,
                              design = ~ groupList)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 1
table(keep)
dds <- dds[keep,]
# 接着计算每个样本的 library size 用于后续的标准化。library size 等于落在该样本所有peaks上的reads的总数，即 counts 矩阵中每列的加和。如果bFullLibrarySize设为TRUE，则会使用library的总reads数（根据原BAM/BED文件统计获得）。然后使用estimateDispersions估计统计分布，需要将参数fitType设为’local’
dds<-estimateSizeFactors(dds)
dds<-estimateDispersions(dds,fitType="local")
dds <- nbinomWaldTest(dds)
res <- results(dds)
## 提取你想要的差异分析结果，我们这里是treated组对untreated组进行比较
res <-  results(dds, contrast=c("groupList", unique(condition)))
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
res_DESeq2 = as.data.frame(res)
table(res_DESeq2$change)#
res_DESeq2 = res_DESeq2[res_DESeq2$change %in% c('UP', 'DOWN'), ]
# 2. Exploring Method 2: limma-------------------------------------------
library(edgeR)
d0 <- DGEList(counts = countMatrix)
dim(d0)
# Filtering low count peaks
cutoff <- 1           
drop <- which(apply(cpm(d0), 1, max) < cutoff)
print(length(drop))
if (length(drop) >0) {
  d <- d0[-drop,] 
} else{ d <- d0}
dim(d) #number of peaks left: 8501      
d <- calcNormFactors(d, method = c('TMM'))
colnames(countMatrix)
genotype <- groupList
table(genotype)
mm <- model.matrix(~0 + genotype)#
v <- voom(d, mm, plot=F, normalize="quantile")#voom from limma:Transform count data to log2-counts per million (logtpm)
fit <- lmFit(v, mm)
head(coef(fit))
contr <- makeContrasts(paste0(colnames(coef(fit)), collapse  = ' - '), levels = colnames(coef(fit))) #Change if other comparison
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)
length(which(top.table$adj.P.Val < 0.05 & abs(top.table$logFC) > log2(1.5))) #how many are significant?
res_limma <- as.data.frame(top.table)
res_limma$change <- with(res_limma, ifelse(res_limma$adj.P.Val < 0.05 & abs(res_limma$logFC) > log2(1.5),
                                           ifelse(res_limma$logFC > log2(1.5),'UP','DOWN'),'NOT'))
table(res_limma$change)#
res_limma = res_limma[res_limma$change %in% c('UP', 'DOWN'), ]

RMdem_list[[i]] <- list(res_DESeq2, res_limma)
name = res@elementMetadata$description[2]
name = paste0(cl, '_',str_split(name, "groupList ", simplify = T)[,2])
name= paste0(gsub(' ','_',(name)))
names(RMdem_list[[i]] ) <- paste0(name, '_', c('DESeq2', 'limma'))
RMdep_list[[i]] <- list( merge_peak )
names(RMdep_list[[i]] ) <- name
countMatrix_list <- list(colnames(countMatrix)); names(countMatrix_list) <- c
condition_list <- list(condition); names(condition_list) <- c
contrast_IP_list[[i]] <- countMatrix_list
condition_IP_list[[i]] <- condition_list

############
n = 75
############
#GSE145924# with RIP and Input
############
print(RMPdeg_datasets[n,])#
i = RMPdeg_datasets[n,'GSE']
cl = c('NHEK and HaCaT and HaCaT-Ras and A431')#Skin cancer
cl = 'HaCaT'
RMp = c('METTL14')
RM = 'm6A'
condition = c(rep('METTL14KD', 2), rep('WT', 2))
RMP_allDat[grep('GSE145924', RMP_allDat$study_alias )
           ,c("study_alias", "sample_alias", "sample_title")]
# View(allDat[[i]])
files <- list.files(file.path(data_directory,i))
files
getwd()
setwd(data_directory)
setwd(i)
ENA <- as.data.frame(RMP_allDat[grep('GSE145924', RMP_allDat$study_alias )
                                , c("study_alias", "sample_alias", "sample_title")])
ENA <- ENA[order(ENA$sample_alias),]
ENA

c <- unique(condition[! condition %in% 'WT'])
###################################
#  Cell number-normalized (CNN)   #
#      MeRIP INPUT analysis       #
###################################
countGTF <- read.table( paste0("metaPlotR/Input_counts_", c, ".txt"),header=T)
rownames(countGTF) <- countGTF$Geneid
countGTF <- countGTF[, !colnames(countGTF) %in% c('Geneid','Chr','Start','End','Strand','Length')]
dim(countGTF) #28265     4
## Sample name
## Sample information
sample=data.frame(
  sample=colnames(countGTF),
  condition=condition)
## SummarizedExperiment input
dds<-DESeqDataSetFromMatrix(countData=countGTF,colData=sample,design=~condition)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 10
table(keep)
dds <- dds[keep,]
## Note on factor levels
dds$condition <- factor(dds$condition, levels = unique(condition))
dds$condition <- relevel(dds$condition, ref = "WT")
# Differential expression analysis-------------------------------------------
dds <- DESeq(dds)

res <-  results(dds, contrast=c("condition",unique(condition)))
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
table(res$change)
table(is.na(res$change))
# resOrdered <- as.data.frame(res[order(res$padj), ])
resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
resInP$geneSymbol <- rownames(resInP)
name = res@elementMetadata$description[2]
name = paste0(cl, '_',str_split(name, "condition ", simplify = T)[,2])
RMdeg_list[[i]] <- list( resInP)
names(RMdeg_list[[i]] ) <- paste0(gsub(' ','_',(name)))
c <- unique(condition[!condition %in% 'WT'])
countGTF_list <- list(colnames(countGTF)); names(countGTF_list) <- c
condition_list <- list(condition); names(condition_list) <- c
contrast_input_list[[i]] <- countGTF_list
condition_input_list[[i]] <- condition_list

###################################
#  Cell number-normalized (CNN)   #
#   MeRIP RIP analysis - Part1    #
###################################
## Sample information
merge_peak <- read.delim(paste0("metaPlotR/", c, "_newPeakFile.txt"),header=F)
countMatrix <- read.table(paste0("metaPlotR/m6A_counts_", c, ".txt"),header=T)
table(countMatrix$Geneid == merge_peak[-1,1])
table(duplicated(merge_peak$V1[-1]))
dim(countMatrix)
colnames(countMatrix)
rownames(countMatrix) <- countMatrix$Geneid
countMatrix <- countMatrix[, !colnames(countMatrix) %in% c('Geneid','Chr','Start','End','Strand','Length')]
# 1. 做好分组因子即可
groupList<-factor(condition)
colData <- data.frame(row.names=colnames(countMatrix), groupList=groupList)
# 2. Exploring Method 1: DESeq2-------------------------------------------
# 构建 DESeqDataSet 对象
dds <- DESeqDataSetFromMatrix(countData = countMatrix,
                              colData = colData,
                              design = ~ groupList)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 1
table(keep)
dds <- dds[keep,]
# 接着计算每个样本的 library size 用于后续的标准化。library size 等于落在该样本所有peaks上的reads的总数，即 counts 矩阵中每列的加和。如果bFullLibrarySize设为TRUE，则会使用library的总reads数（根据原BAM/BED文件统计获得）。然后使用estimateDispersions估计统计分布，需要将参数fitType设为’local’
dds<-estimateSizeFactors(dds)
dds<-estimateDispersions(dds,fitType="local")
dds <- nbinomWaldTest(dds)
res <- results(dds)
## 提取你想要的差异分析结果，我们这里是treated组对untreated组进行比较
res <-  results(dds, contrast=c("groupList", unique(condition)))
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
res_DESeq2 = as.data.frame(res)
table(res_DESeq2$change)#
res_DESeq2 = res_DESeq2[res_DESeq2$change %in% c('UP', 'DOWN'), ]
# 2. Exploring Method 2: limma-------------------------------------------
library(edgeR)
d0 <- DGEList(counts = countMatrix)
dim(d0)
# Filtering low count peaks
cutoff <- 1           
drop <- which(apply(cpm(d0), 1, max) < cutoff)
print(length(drop))
if (length(drop) >0) {
  d <- d0[-drop,] 
} else{ d <- d0}
dim(d) #number of peaks left: 8501      
d <- calcNormFactors(d, method = c('TMM'))
colnames(countMatrix)
genotype <- groupList
table(genotype)
mm <- model.matrix(~0 + genotype)#
v <- voom(d, mm, plot=F, normalize="quantile")#voom from limma:Transform count data to log2-counts per million (logtpm)
fit <- lmFit(v, mm)
head(coef(fit))
contr <- makeContrasts(paste0(colnames(coef(fit)), collapse  = ' - '), levels = colnames(coef(fit))) #Change if other comparison
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)
length(which(top.table$adj.P.Val < 0.05 & abs(top.table$logFC) > log2(1.5))) #how many are significant?
res_limma <- as.data.frame(top.table)
res_limma$change <- with(res_limma, ifelse(res_limma$adj.P.Val < 0.05 & abs(res_limma$logFC) > log2(1.5),
                                           ifelse(res_limma$logFC > log2(1.5),'UP','DOWN'),'NOT'))
table(res_limma$change)#
res_limma = res_limma[res_limma$change %in% c('UP', 'DOWN'), ]

RMdem_list[[i]] <- list(res_DESeq2, res_limma)
name = res@elementMetadata$description[2]
name = paste0(cl, '_',str_split(name, "groupList ", simplify = T)[,2])
name= paste0(gsub(' ','_',(name)))
names(RMdem_list[[i]] ) <- paste0(name, '_', c('DESeq2', 'limma'))
RMdep_list[[i]] <- list( merge_peak )
names(RMdep_list[[i]] ) <- name
countMatrix_list <- list(colnames(countMatrix)); names(countMatrix_list) <- c
condition_list <- list(condition); names(condition_list) <- c
contrast_IP_list[[i]] <- countMatrix_list
condition_IP_list[[i]] <- condition_list

contrast_input_list[['GSE145924']]$METTL14KD <- gsub('M14.input.', 'M14-input-', contrast_input_list[['GSE145924']]$METTL14KD )
contrast_IP_list[['GSE145924']]$METTL14KD <- gsub('M14.ip.', 'M14-ip-', contrast_IP_list[['GSE145924']]$METTL14KD )

############
n = 76
############
#GSE147891# without case RIP, with RNA
############
print(RMPdeg_datasets[n,])#
i = RMPdeg_datasets[n,'GSE']
cl = c('LNCaP')#, 'RWPE-1')#Prostate cancer
RMp = c('METTL3')
RM = 'm6A'
condition_list = list(LNCaP_METTL3KD = c(rep('METTL3KD', 2), rep('WT', 2)),
                      LNCaP.shM3_METTL3MU = c(rep('METTL3MU', 2), rep('WT', 2))
                      )
RMP_allDat[grep('GSE147891', RMP_allDat$study_alias )
           ,c("study_alias", "sample_alias", "sample_title")]
# View(allDat[[i]])
files <- list.files(file.path(data_directory,i))
files
getwd()
setwd(data_directory)
setwd(i)
df_in <- read.table('GSE147884_Counts.txt.gz')
## Gene name
library(data.table)
countGTF=df_in[, -1]
###################################
#   more than 1 condition here   #
###################################
resInP_list <- list()
##Gene count matrix
# countGTF=as.matrix(ceiling(df_in[, -1]))
## Gene name
table(duplicated(df_in$geneSymbol))
dim(countGTF) #25343    12
## Sample name
colnames(countGTF)
## Sample information
sample_meta=data.frame(
  sample=colnames(countGTF),
  condition=c(rep("LNCaP_WT",2), rep("LNCaP.shM3_METTL3MU",2),
              rep("LNCaP.shM3_WT",2), rep("LNCaP_METTL3KD",2)))
#####condition_LNCaP.shM3_METTL3MU_vs_WT#####
count <- as.matrix(countGTF[, grepl('.shM3', sample_meta$condition)])
table(is.na(count))
sample <- sample_meta[sample_meta$sample %in% colnames(count), ]
## SummarizedExperiment input
condition = c(rep('LNCaP.shM3_METTL3MU', 2), rep('LNCaP.shM3_WT', 2))
dds<-DESeqDataSetFromMatrix(countData=count, colData=sample,design=~condition)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 10
table(keep)
dds <- dds[keep,]
## Note on factor levels
dds$condition <- factor(dds$condition, levels = c("LNCaP.shM3_METTL3MU","LNCaP.shM3_WT"))
dds$condition <- relevel(dds$condition, ref = "LNCaP.shM3_WT")
dds <- DESeq(dds)
res <-  results(dds, contrast=c("condition",unique(condition))) # c("condition", rev(levels(genotype)))
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
table(res$change)
res$geneSymbol <- rownames(countGTF)[keep]
resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
table(duplicated(resInP$geneSymbol))
resInP_list[['LNCaP.shM3_METTL3MU_vs_WT']] <- resInP
#####condition_LNCaP_METTL3KD_vs_WT#####
count <- as.matrix(countGTF[, !grepl('.shM3', sample_meta$condition)])
table(is.na(count))
sample <- sample_meta[sample_meta$sample %in% colnames(count), ]
## SummarizedExperiment input
condition = c(rep('LNCaP_METTL3KD', 2), rep('LNCaP_WT', 2))
dds<-DESeqDataSetFromMatrix(countData=count, colData=sample,design=~condition)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 10
table(keep)
dds <- dds[keep,]
## Note on factor levels
dds$condition <- factor(dds$condition, levels = c("LNCaP_METTL3KD","LNCaP_WT"))
dds$condition <- relevel(dds$condition, ref = "LNCaP_WT")
dds <- DESeq(dds)
res <-  results(dds, contrast=c("condition",unique(condition))) # c("condition", rev(levels(genotype)))
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
table(res$change)
res$geneSymbol <- rownames(countGTF)[keep]
resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
table(duplicated(resInP$geneSymbol))
resInP_list[['LNCaP_METTL3KD_vs_WT']] <- resInP
RMdeg_list[[i]] <- resInP_list
contrast_input_list[[i]] <- NULL#lapply(countGTF_list, extract_names)
condition_input_list[[i]] <- condition_list


############
n = 77
############
#GSE119026# without case RIP, with RNA
############
print(RMPdeg_datasets[n,])##undone! for SND1PD_vs_m6A
i = RMPdeg_datasets[n,'GSE']
cl = c('BCBL1-Rta')#Kaposi's sarcoma
RMp = c('SND1')
RM = 'm6A'
condition_list = list( BCBL1.la_SND1KO = c(rep('SND1KO', 2), rep('WT', 2)),
                       BCBL1.ly_SND1KO = c(rep('SND1KO', 2), rep('WT', 2)),
                       BCBL1.TREX.la_SND1KO = c(rep('SND1KO', 2), rep('WT', 2)),
                       BCBL1.TREX.ly_SND1KO = c(rep('SND1KO', 2), rep('WT', 2))
)
RMP_allDat[grep('GSE119026', RMP_allDat$study_alias )
           ,c("study_alias", "sample_alias", "sample_title")]
# View(allDat[[i]])
paste0(skip$sample_alias, collapse = "|")
files <- list.files(file.path(data_directory,i))
files
getwd()
setwd(data_directory)
setwd(i)
## Gene name
library(data.table)
countGTF=df_in[, -1]
###################################
#   more than 1 condition here   #
###################################
resInP_list <- list()
##Gene count matrix
df_in1 <- read_excel('GSE119026_RAW_READ_COUNTS_SND1_KO_BCBL.xlsx',sheet = 1)
df_in2 <- read_excel('GSE119026_RAW_READ_COUNTS_SND1_KO_BCBL.xlsx',sheet = 2)
colnames(df_in1);colnames(df_in2)
df_in <- left_join(df_in1[,-(2:6)], df_in2[,-c(1,3:7)], by = 'GeneID')
countGTF=as.data.frame(ceiling(df_in[, -1]))
id_mapped <- read.table("/data2/rluo4/hg38/gencode.v36.annotation.gene.probeMap",header = T)
dim(id_mapped) #  60675     6
table(str_split(df_in$GeneID, '[.]', simplify = T)[,1] %in% str_split(id_mapped$id, '[.]', simplify = T)[,1])
countGTF$id <- str_split(df_in$GeneID, '[.]', simplify = T)[,1]
id_mapped <- id_mapped[which(id_mapped$id %in% countGTF$id),]
dim(id_mapped) #57930     6
table(duplicated(id_mapped$gene))
countGTF <- left_join(id_mapped,countGTF,by="id")
countGTF <- aggregate(x = countGTF[,(ncol(id_mapped)+1):ncol(countGTF)],   #此时exprSet的第三列开始是表达矩阵内容
                      by = list(symbol = countGTF$gene),   #按照相同symbol分组，在组内计算
                      FUN = max) %>%   #原文中是计算最大值（max），也可以计算平均值（mean）或者中位数（median）
  column_to_rownames(var = 'symbol')
dim(countGTF) #56752     5
## Sample name
colnames(countGTF)
## Sample information
sample_meta=data.frame(
  sample=colnames(countGTF),
  condition=c("BCBL1.la_WT","BCBL1.la_SND1KO",  "BCBL1.ly_WT","BCBL1.ly_SND1KO",
              "BCBL1.la_SND1KO", "BCBL1.la_WT", "BCBL1.ly_SND1KO", "BCBL1.ly_WT"))
#####condition_BCBL1.la_SND1KO_vs_WT#####
count <- as.matrix(countGTF[, grepl('.la', sample_meta$condition)])
table(is.na(count))
sample <- sample_meta[sample_meta$sample %in% colnames(count), ]
## SummarizedExperiment input
condition = c(rep('BCBL1.la_SND1KO', 2), rep('BCBL1.la_WT', 2))
dds<-DESeqDataSetFromMatrix(countData=count, colData=sample,design=~condition)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 10
table(keep)
dds <- dds[keep,]
## Note on factor levels
dds$condition <- factor(dds$condition, levels = c("BCBL1.la_SND1KO","BCBL1.la_WT"))
dds$condition <- relevel(dds$condition, ref = "BCBL1.la_WT")
dds <- DESeq(dds)
res <-  results(dds, contrast=c("condition",unique(condition))) # c("condition", rev(levels(genotype)))
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
table(res$change)
res$geneSymbol <- rownames(countGTF)[keep]
resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
table(duplicated(resInP$geneSymbol))
resInP_list[['BCBL1.la_SND1KO_vs_WT']] <- resInP
#####condition_BCBL1.ly_SND1KO_vs_WT#####
count <- as.matrix(countGTF[, grepl('.ly', sample_meta$condition)])
table(is.na(count))
sample <- sample_meta[sample_meta$sample %in% colnames(count), ]
## SummarizedExperiment input
condition = c(rep('BCBL1.ly_SND1KO', 2), rep('BCBL1.ly_WT', 2))
dds<-DESeqDataSetFromMatrix(countData=count, colData=sample,design=~condition)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 10
table(keep)
dds <- dds[keep,]
## Note on factor levels
dds$condition <- factor(dds$condition, levels = c("BCBL1.ly_SND1KO","BCBL1.ly_WT"))
dds$condition <- relevel(dds$condition, ref = "BCBL1.ly_WT")
dds <- DESeq(dds)
res <-  results(dds, contrast=c("condition",unique(condition))) # c("condition", rev(levels(genotype)))
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
table(res$change)
res$geneSymbol <- rownames(countGTF)[keep]
resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
table(duplicated(resInP$geneSymbol))
resInP_list[['BCBL1.ly_SND1KO_vs_WT']] <- resInP

##Gene count matrix
df_in <- fread('GSE119026_RAW_READ_COUNTS_SND1_KO_TREX.txt.gz')
dim(df_in)
df_in <- df_in[df_in$ENS !='']
colnames(df_in)
countGTF=as.data.frame(ceiling(df_in[, -(1:7)]))
id_mapped <- read.table("/data2/rluo4/hg38/gencode.v36.annotation.gene.probeMap",header = T)
dim(id_mapped) #  60675     6
table(str_split(df_in$GeneID, '[.]', simplify = T)[,1] %in% str_split(id_mapped$id, '[.]', simplify = T)[,1])
countGTF$id <- str_split(df_in$GeneID, '[.]', simplify = T)[,1]
id_mapped <- id_mapped[which(id_mapped$id %in% countGTF$id),]
dim(id_mapped) #57930     6
table(duplicated(id_mapped$gene))
countGTF <- left_join(id_mapped,countGTF,by="id")
countGTF <- aggregate(x = countGTF[,(ncol(id_mapped)+1):ncol(countGTF)],   #此时exprSet的第三列开始是表达矩阵内容
                      by = list(symbol = countGTF$gene),   #按照相同symbol分组，在组内计算
                      FUN = max) %>%   #原文中是计算最大值（max），也可以计算平均值（mean）或者中位数（median）
  column_to_rownames(var = 'symbol')
dim(countGTF) #56752     5
## Sample name
colnames(countGTF)
## Sample information
sample_meta=data.frame(
  sample=colnames(countGTF),
  condition=c(rep("BCBL1.TREX.la_SND1KO",2), rep( "BCBL1.TREX.la_WT", 2),
              rep("BCBL1.TREX.ly_SND1KO",2), rep( "BCBL1.TREX.ly_WT", 2)))
#####condition_BCBL1.TREX.la_SND1KO_vs_WT#####
count <- as.matrix(countGTF[, grepl('.la', sample_meta$condition)])
table(is.na(count))
sample <- sample_meta[sample_meta$sample %in% colnames(count), ]
## SummarizedExperiment input
condition = c(rep('BCBL1.TREX.la_SND1KO', 2), rep('BCBL1.TREX.la_WT', 2))
dds<-DESeqDataSetFromMatrix(countData=count, colData=sample,design=~condition)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 10
table(keep)
dds <- dds[keep,]
## Note on factor levels
dds$condition <- factor(dds$condition, levels = c("BCBL1.TREX.la_SND1KO","BCBL1.TREX.la_WT"))
dds$condition <- relevel(dds$condition, ref = "BCBL1.TREX.la_WT")
dds <- DESeq(dds)
res <-  results(dds, contrast=c("condition",unique(condition))) # c("condition", rev(levels(genotype)))
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
table(res$change)
res$geneSymbol <- rownames(countGTF)[keep]
resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
table(duplicated(resInP$geneSymbol))
resInP_list[['BCBL1.TREX.la_SND1KO_vs_WT']] <- resInP
#####condition_BCBL1.TREX.ly_SND1KO_vs_WT#####
count <- as.matrix(countGTF[, grepl('.ly', sample_meta$condition)])
table(is.na(count))
sample <- sample_meta[sample_meta$sample %in% colnames(count), ]
## SummarizedExperiment input
condition = c(rep('BCBL1.TREX.ly_SND1KO', 2), rep('BCBL1.TREX.ly_WT', 2))
dds<-DESeqDataSetFromMatrix(countData=count, colData=sample,design=~condition)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 10
table(keep)
dds <- dds[keep,]
## Note on factor levels
dds$condition <- factor(dds$condition, levels = c("BCBL1.TREX.ly_SND1KO","BCBL1.TREX.ly_WT"))
dds$condition <- relevel(dds$condition, ref = "BCBL1.TREX.ly_WT")
dds <- DESeq(dds)
res <-  results(dds, contrast=c("condition",unique(condition))) # c("condition", rev(levels(genotype)))
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
table(res$change)
res$geneSymbol <- rownames(countGTF)[keep]
resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
table(duplicated(resInP$geneSymbol))
resInP_list[['BCBL1.TREX.ly_SND1KO_vs_WT']] <- resInP

RMdeg_list[[i]] <- resInP_list
contrast_input_list[[i]] <- NULL#lapply(countGTF_list, extract_names)
condition_input_list[[i]] <- condition_list

############
n = 78
############
#GSE114019# with RIP and Input
############
print(RMPdeg_datasets[n,])#
i = RMPdeg_datasets[n,'GSE']
cl = c('Fibroblasts and THP1')#Foreskin and blood #post infection (hpi)
cl = 'Fibroblasts'
RMp = c('METTL3')
RM = 'm6A'
condition = c(rep('METTL3KD', 3), rep('WT', 3))#METTL3-depleted sgRNA
RMP_allDat[grep('GSE114019', RMP_allDat$study_alias )
           ,c("study_alias", "sample_alias", "sample_title")]
# View(allDat[[i]])
RMP_allDat[grepl('GSE114019', RMP_allDat$study_alias ) & (!grepl('hpi|uv', RMP_allDat$sample_title))
           ,c("study_alias", "sample_alias", "sample_title")]
skip = RMP_allDat[grepl('GSE114019', RMP_allDat$study_alias ) & (!grepl('hpi|uv', RMP_allDat$sample_title))
                  ,c("study_alias", "sample_alias", "sample_title")]
paste0(skip$sample_alias, collapse = "|")

skip <- RMP_allDat[grepl('GSE114019', RMP_allDat$study_alias ) & (!grepl('hpi|uv', RMP_allDat$sample_title))
                   &  grepl('input', RMP_allDat$sample_title)
                   , "sample_alias"]
paste0(skip$sample_alias, collapse = " ")
skip <- RMP_allDat[grepl('GSE114019', RMP_allDat$study_alias ) & (!grepl('hpi|uv', RMP_allDat$sample_title))
                   &  grepl('IP', RMP_allDat$sample_title)
                   , "sample_alias"]
paste0(skip$sample_alias, collapse = " ")
ENA <- as.data.frame(RMP_allDat[grepl('GSE114019', RMP_allDat$study_alias ) & (!grepl('hpi|uv', RMP_allDat$sample_title))
                                , c("study_alias", "sample_alias", "sample_title")])
ENA <- ENA[order(ENA$sample_alias),]
ENA
files <- list.files(file.path(data_directory,i))
files
getwd()
setwd(data_directory)
setwd(i)
c <- unique(condition[! condition %in% 'WT'])
###################################
#  Cell number-normalized (CNN)   #
#      MeRIP INPUT analysis       #
###################################
countGTF <- read.table( paste0("metaPlotR/Input_counts_", c, ".txt"),header=T)
rownames(countGTF) <- countGTF$Geneid
countGTF <- countGTF[, !colnames(countGTF) %in% c('Geneid','Chr','Start','End','Strand','Length')]
dim(countGTF) #28265     4
## Sample name
## Sample information
sample=data.frame(
  sample=colnames(countGTF),
  condition=condition)
## SummarizedExperiment input
dds<-DESeqDataSetFromMatrix(countData=countGTF,colData=sample,design=~condition)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 10
table(keep)
dds <- dds[keep,]
## Note on factor levels
dds$condition <- factor(dds$condition, levels = unique(condition))
dds$condition <- relevel(dds$condition, ref = "WT")
# Differential expression analysis-------------------------------------------
dds <- DESeq(dds)

res <-  results(dds, contrast=c("condition",unique(condition)))
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
table(res$change)
table(is.na(res$change))
# resOrdered <- as.data.frame(res[order(res$padj), ])
resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
resInP$geneSymbol <- rownames(resInP)
name = res@elementMetadata$description[2]
name = paste0(cl, '_',str_split(name, "condition ", simplify = T)[,2])
RMdeg_list[[i]] <- list( resInP)
names(RMdeg_list[[i]] ) <- paste0(gsub(' ','_',(name)))
c <- unique(condition[!condition %in% 'WT'])
countGTF_list <- list(colnames(countGTF)); names(countGTF_list) <- c
condition_list <- list(condition); names(condition_list) <- c
contrast_input_list[[i]] <- countGTF_list
condition_input_list[[i]] <- condition_list

###################################
#  Cell number-normalized (CNN)   #
#   MeRIP RIP analysis - Part1    #
###################################
## Sample information
merge_peak <- read.delim(paste0("metaPlotR/", c, "_newPeakFile.txt"),header=F)
countMatrix <- read.table(paste0("metaPlotR/m6A_counts_", c, ".txt"),header=T)
table(countMatrix$Geneid == merge_peak[-1,1])
table(duplicated(merge_peak$V1[-1]))
dim(countMatrix)
colnames(countMatrix)
rownames(countMatrix) <- countMatrix$Geneid
countMatrix <- countMatrix[, !colnames(countMatrix) %in% c('Geneid','Chr','Start','End','Strand','Length')]
# 1. 做好分组因子即可
groupList<-factor(condition)
colData <- data.frame(row.names=colnames(countMatrix), groupList=groupList)
# 2. Exploring Method 1: DESeq2-------------------------------------------
# 构建 DESeqDataSet 对象
dds <- DESeqDataSetFromMatrix(countData = countMatrix,
                              colData = colData,
                              design = ~ groupList)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 1
table(keep)
dds <- dds[keep,]
# 接着计算每个样本的 library size 用于后续的标准化。library size 等于落在该样本所有peaks上的reads的总数，即 counts 矩阵中每列的加和。如果bFullLibrarySize设为TRUE，则会使用library的总reads数（根据原BAM/BED文件统计获得）。然后使用estimateDispersions估计统计分布，需要将参数fitType设为’local’
dds<-estimateSizeFactors(dds)
dds<-estimateDispersions(dds,fitType="local")
dds <- nbinomWaldTest(dds)
res <- results(dds)
## 提取你想要的差异分析结果，我们这里是treated组对untreated组进行比较
res <-  results(dds, contrast=c("groupList", unique(condition)))
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
res_DESeq2 = as.data.frame(res)
table(res_DESeq2$change)#
res_DESeq2 = res_DESeq2[res_DESeq2$change %in% c('UP', 'DOWN'), ]
# 2. Exploring Method 2: limma-------------------------------------------
library(edgeR)
d0 <- DGEList(counts = countMatrix)
dim(d0)
# Filtering low count peaks
cutoff <- 1           
drop <- which(apply(cpm(d0), 1, max) < cutoff)
print(length(drop))
if (length(drop) >0) {
  d <- d0[-drop,] 
} else{ d <- d0}
dim(d) #number of peaks left: 8501      
d <- calcNormFactors(d, method = c('TMM'))
colnames(countMatrix)
genotype <- groupList
table(genotype)
mm <- model.matrix(~0 + genotype)#
v <- voom(d, mm, plot=F, normalize="quantile")#voom from limma:Transform count data to log2-counts per million (logtpm)
fit <- lmFit(v, mm)
head(coef(fit))
contr <- makeContrasts(paste0(colnames(coef(fit)), collapse  = ' - '), levels = colnames(coef(fit))) #Change if other comparison
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)
length(which(top.table$adj.P.Val < 0.05 & abs(top.table$logFC) > log2(1.5))) #how many are significant?
res_limma <- as.data.frame(top.table)
res_limma$change <- with(res_limma, ifelse(res_limma$adj.P.Val < 0.05 & abs(res_limma$logFC) > log2(1.5),
                                           ifelse(res_limma$logFC > log2(1.5),'UP','DOWN'),'NOT'))
table(res_limma$change)#
res_limma = res_limma[res_limma$change %in% c('UP', 'DOWN'), ]

RMdem_list[[i]] <- list(res_DESeq2, res_limma)
name = res@elementMetadata$description[2]
name = paste0(cl, '_',str_split(name, "groupList ", simplify = T)[,2])
name= paste0(gsub(' ','_',(name)))
names(RMdem_list[[i]] ) <- paste0(name, '_', c('DESeq2', 'limma'))
RMdep_list[[i]] <- list( merge_peak )
names(RMdep_list[[i]] ) <- name
countMatrix_list <- list(colnames(countMatrix)); names(countMatrix_list) <- c
condition_list <- list(condition); names(condition_list) <- c
contrast_IP_list[[i]] <- countMatrix_list
condition_IP_list[[i]] <- condition_list
# 
############
n = 79
############
#GSE167075# with RIP and Input
############
print(RMPdeg_datasets[n,])#hg38
i = RMPdeg_datasets[n,'GSE']
cl = c('Caco-2;SARS-CoV-2')#Human colorectal adenocarcinoma Caco-2 cells;SARS-CoV-2 virus
cl = 'Caco-2'
RMp = c('METTL3')
RM = 'm6A'
condition = c(rep('METTL3KD', 2), rep('WT', 2))#
RMP_allDat[grep('GSE167075', RMP_allDat$study_alias )
           ,c("study_alias", "sample_alias", "sample_title")]
# View(allDat[[i]])
ENA <- as.data.frame(RMP_allDat[grepl('GSE167075', RMP_allDat$study_alias ) & (!grepl('shNTC|shM3', RMP_allDat$sample_title))
                                , c("study_alias", "sample_alias", "sample_title")])
ENA <- ENA[order(ENA$sample_alias),]
ENA
paste0(ENA$sample_alias, collapse = "|")
skip <- ENA[ grepl('Input', ENA$sample_title)
             , ]
paste0(skip$sample_alias, collapse = " ")
skip <- ENA[ grepl('IP', ENA$sample_title)
             , ]
paste0(skip$sample_alias, collapse = " ")

files <- list.files(file.path(data_directory,i))
files
getwd()
setwd(data_directory)
setwd(i)
DEG <- read_excel(file.path(data_directory,i, 'GSE167075_shM3_vs_NTC_human_.xlsx'),sheet = 1)
c <- unique(condition[!condition %in% 'WT'])
###################################
#  Cell number-normalized (CNN)   #
#      MeRIP INPUT analysis       #
###################################
countGTF <- read.table( paste0("metaPlotR/Input_counts_", c, ".txt"),header=T)
rownames(countGTF) <- countGTF$Geneid
countGTF <- countGTF[, !colnames(countGTF) %in% c('Geneid','Chr','Start','End','Strand','Length')]
dim(countGTF) #28265     4
## Sample name
## Sample information
sample=data.frame(
  sample=colnames(countGTF),
  condition=condition)
## SummarizedExperiment input
dds<-DESeqDataSetFromMatrix(countData=countGTF,colData=sample,design=~condition)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 10
table(keep)
dds <- dds[keep,]
## Note on factor levels
dds$condition <- factor(dds$condition, levels = unique(condition))
dds$condition <- relevel(dds$condition, ref = "WT")
# Differential expression analysis-------------------------------------------
dds <- DESeq(dds)

res <-  results(dds, contrast=c("condition",unique(condition)))
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
table(res$change)
table(is.na(res$change))
# resOrdered <- as.data.frame(res[order(res$padj), ])
resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
resInP$geneSymbol <- rownames(resInP)
name = res@elementMetadata$description[2]
name = paste0(cl, '_',str_split(name, "condition ", simplify = T)[,2])
RMdeg_list[[i]] <- list( resInP)
names(RMdeg_list[[i]] ) <- paste0(gsub(' ','_',(name)))
countGTF_list <- list(colnames(countGTF)); names(countGTF_list) <- c
condition_list <- list(condition); names(condition_list) <- c
contrast_input_list[[i]] <- countGTF_list
condition_input_list[[i]] <- condition_list

###################################
#  Cell number-normalized (CNN)   #
#   MeRIP RIP analysis - Part1    #
###################################
## Sample information
merge_peak <- read.delim(paste0("metaPlotR/", c, "_newPeakFile.txt"),header=F)
countMatrix <- read.table(paste0("metaPlotR/m6A_counts_", c, ".txt"),header=T)
table(countMatrix$Geneid == merge_peak[-1,1])
table(duplicated(merge_peak$V1[-1]))
dim(countMatrix)
colnames(countMatrix)
rownames(countMatrix) <- countMatrix$Geneid
countMatrix <- countMatrix[, !colnames(countMatrix) %in% c('Geneid','Chr','Start','End','Strand','Length')]
# 1. 做好分组因子即可
groupList<-factor(condition)
colData <- data.frame(row.names=colnames(countMatrix), groupList=groupList)
# 2. Exploring Method 1: DESeq2-------------------------------------------
# 构建 DESeqDataSet 对象
dds <- DESeqDataSetFromMatrix(countData = countMatrix,
                              colData = colData,
                              design = ~ groupList)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 1
table(keep)
dds <- dds[keep,]
# 接着计算每个样本的 library size 用于后续的标准化。library size 等于落在该样本所有peaks上的reads的总数，即 counts 矩阵中每列的加和。如果bFullLibrarySize设为TRUE，则会使用library的总reads数（根据原BAM/BED文件统计获得）。然后使用estimateDispersions估计统计分布，需要将参数fitType设为’local’
dds<-estimateSizeFactors(dds)
dds<-estimateDispersions(dds,fitType="local")
dds <- nbinomWaldTest(dds)
res <- results(dds)
## 提取你想要的差异分析结果，我们这里是treated组对untreated组进行比较
res <-  results(dds, contrast=c("groupList", unique(condition)))
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
res_DESeq2 = as.data.frame(res)
table(res_DESeq2$change)#Not only
res_DESeq2 = res_DESeq2[res_DESeq2$change %in% c('UP', 'DOWN'), ]
# 2. Exploring Method 2: limma-------------------------------------------
library(edgeR)
d0 <- DGEList(counts = countMatrix)
dim(d0)
# Filtering low count peaks
cutoff <- 1           
drop <- which(apply(cpm(d0), 1, max) < cutoff)
print(length(drop))
if (length(drop) >0) {
  d <- d0[-drop,] 
} else{ d <- d0}
dim(d) #number of peaks left: 8501      
d <- calcNormFactors(d, method = c('TMM'))
colnames(countMatrix)
genotype <- groupList
table(genotype)
mm <- model.matrix(~0 + genotype)#
v <- voom(d, mm, plot=F, normalize="quantile")#voom from limma:Transform count data to log2-counts per million (logtpm)
fit <- lmFit(v, mm)
head(coef(fit))
contr <- makeContrasts(paste0(colnames(coef(fit)), collapse  = ' - '), levels = colnames(coef(fit))) #Change if other comparison
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)
length(which(top.table$adj.P.Val < 0.05 & abs(top.table$logFC) > log2(1.5))) #how many are significant?
res_limma <- as.data.frame(top.table)
res_limma$change <- with(res_limma, ifelse(res_limma$adj.P.Val < 0.05 & abs(res_limma$logFC) > log2(1.5),
                                           ifelse(res_limma$logFC > log2(1.5),'UP','DOWN'),'NOT'))
table(res_limma$change)#NOT only
res_limma = res_limma[res_limma$change %in% c('UP', 'DOWN'), ]

RMdem_list[[i]] <- list(res_DESeq2, res_limma)
name = res@elementMetadata$description[2]
name = paste0(cl, '_',str_split(name, "groupList ", simplify = T)[,2])
name= paste0(gsub(' ','_',(name)))
names(RMdem_list[[i]] ) <- paste0(name, '_', c('DESeq2', 'limma'))
RMdep_list[[i]] <- list( merge_peak )
names(RMdep_list[[i]] ) <- name
countMatrix_list <- list(colnames(countMatrix)); names(countMatrix_list) <- c
condition_list <- list(condition); names(condition_list) <- c
contrast_IP_list[[i]] <- countMatrix_list
condition_IP_list[[i]] <- condition_list

############
n = 80
############
#GSE124509# with RIP and Input
############
print(RMPdeg_datasets[n,])#hg38
i = RMPdeg_datasets[n,'GSE']
cl = c('HEK293T')#
RMp = c('METTL3','METTL16','FTO','ALKBH5','PCIF1')
RM = 'm6A'
condition_list = list(
  METTL3KO = c(rep('METTL3KO', 3), rep('WT', 9)),
  METTL16KD = c(rep('METTL16KD', 3), rep('WT', 9)),
  FTOKO = c(rep('FTOKO', 3), rep('WT', 9)),
  ALKBH5KO = c(rep('ALKBH5KO', 3), rep('WT', 9)),
  PCIF1KO = c(rep('PCIF1KO', 3), rep('WT', 9)),
  FTOOE = c(rep('FTOOE', 3), rep('WT', 9))
)
RMP_allDat[grep('GSE124509', RMP_allDat$study_alias )
           ,c("study_alias", "sample_alias", "sample_title")]
# View(allDat[[i]])
ENA <- as.data.frame(RMP_allDat[grepl('GSE124509', RMP_allDat$study_alias ) & (!grepl(' Mettl3-KO ', RMP_allDat$sample_title))
                                , c("study_alias", "sample_alias", "sample_title")])
ENA <- ENA[order(ENA$sample_alias),]
ENA

files <- list.files(file.path(data_directory,i))
files
getwd()
setwd(data_directory)
setwd(i)
###################################
#  Cell number-normalized (CNN)   #
#      MeRIP INPUT analysis       #
###################################
resInP_list <- vector("list", length(condition_list))
names(resInP_list) <- names(condition_list)
countGTF_list <- vector("list", length(condition_list))
names(countGTF_list) <- names(condition_list)
for(c in names(condition_list)){
  countGTF <- read.table( paste0("metaPlotR/Input_counts_", c, ".txt"),header=T)
  dim(countGTF)
  colnames(countGTF)
  rownames(countGTF) <- countGTF$Geneid
  countGTF <- countGTF[, !colnames(countGTF) %in% c('Geneid','Chr','Start','End','Strand','Length')]
  dim(countGTF) # 28265     4
  condition <- condition_list[[c]]
  ## Sample information
  sample=data.frame(
    sample=colnames(countGTF),
    condition=condition)
  ## SummarizedExperiment input
  dds<-DESeqDataSetFromMatrix(countData=countGTF,colData=sample,design=~condition)
  ## Pre-filtering
  keep <- rowSums(counts(dds)) >= 10
  table(keep)
  dds <- dds[keep,]
  ## Note on factor levels
  dds$condition <- factor(dds$condition, levels = unique(condition))
  dds$condition <- relevel(dds$condition, ref = "WT")
  # Differential expression analysis-------------------------------------------
  dds <- DESeq(dds)
  res <-  results(dds, contrast=c("condition",unique(condition)))
  res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
  length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
  res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                                 ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
  table(res$change)
  table(is.na(res$change))
  # resOrdered <- as.data.frame(res[order(res$padj), ])
  resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
  resInP$geneSymbol <- rownames(resInP)
  table(resInP$change)
  resInP_list[[c]] <- resInP
  countGTF_list[[c]] <- countGTF 
}
name = paste0(cl, '_',names(resInP_list), ' vs WT')
RMdeg_list[[i]] <- resInP_list 
names(RMdeg_list[[i]] ) <- paste0(gsub(' ','_',(name)))
contrast_input_list[[i]] <- lapply(countGTF_list, extract_names)
condition_input_list[[i]] <- condition_list

###################################
#  Cell number-normalized (CNN)   #
#   MeRIP RIP analysis - Part1    #
###################################
resMeR_list <- vector("list", length(condition_list))
names(resMeR_list) <- names(condition_list)
countMatrix_list <- vector("list", length(condition_list))
names(countMatrix_list) <- names(condition_list)
resPeak_list <- vector("list",length(condition_list))
names(resPeak_list) <- names(condition_list)
for(c in names(condition_list)){
  merge_peak <- read.delim(paste0("metaPlotR/", c, "_newPeakFile.txt"),header=F)
  countMatrix <- read.table(paste0("metaPlotR/m6A_counts_", c, ".txt"),header=T)
  table(countMatrix$Geneid == merge_peak[-1,1])
  table(duplicated(merge_peak[-1,1]))
  dim(countMatrix)
  colnames(countMatrix)
  rownames(countMatrix) <- countMatrix$Geneid
  countMatrix <- countMatrix[, !colnames(countMatrix) %in% c('Geneid','Chr','Start','End','Strand','Length')]
  # 1. 做好分组因子即可
  condition <- condition_list[[c]]
  groupList<-factor(condition)
  colData <- data.frame(row.names=colnames(countMatrix), groupList=groupList)
  
  # 2. Exploring Method 1: DESeq2-------------------------------------------
  # 构建 DESeqDataSet 对象
  dds <- DESeqDataSetFromMatrix(countData = countMatrix,
                                colData = colData,
                                design = ~ groupList)
  ## Pre-filtering
  keep <- rowSums(counts(dds)) >= 1
  table(keep)
  dds <- dds[keep,]
  # 接着计算每个样本的 library size 用于后续的标准化。library size 等于落在该样本所有peaks上的reads的总数，即 counts 矩阵中每列的加和。如果bFullLibrarySize设为TRUE，则会使用library的总reads数（根据原BAM/BED文件统计获得）。然后使用estimateDispersions估计统计分布，需要将参数fitType设为’local’
  dds<-estimateSizeFactors(dds)
  dds<-estimateDispersions(dds,fitType="local")
  dds <- nbinomWaldTest(dds)
  res <- results(dds)
  ## 提取你想要的差异分析结果，我们这里是treated组对untreated组进行比较
  res <-  results(dds, contrast=c("groupList", unique(condition)))
  res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
  length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
  res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                                 ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
  res_DESeq2 = as.data.frame(res)
  table(res_DESeq2$change)
  res_DESeq2 = res_DESeq2[res_DESeq2$change %in% c('UP', 'DOWN'), ]
  # 2. Exploring Method 2: limma-------------------------------------------
  library(edgeR)
  d0 <- DGEList(counts = countMatrix)
  dim(d0)
  # Filtering low count peaks
  cutoff <- 1           
  drop <- which(apply(cpm(d0), 1, max) < cutoff)
  print(length(drop))
  if (length(drop) >0) {
    d <- d0[-drop,] 
  } else{ d <- d0}
  dim(d) #number of peaks left: 8501      
  d <- calcNormFactors(d, method = c('TMM'))
  colnames(countMatrix)
  genotype <- groupList
  table(genotype)
  mm <- model.matrix(~0 + genotype)#
  v <- voom(d, mm, plot=F, normalize="quantile")#voom from limma:Transform count data to log2-counts per million (logtpm)
  fit <- lmFit(v, mm)
  head(coef(fit))
  contr <- makeContrasts(paste0(colnames(coef(fit)), collapse  = ' - '), levels = colnames(coef(fit))) #Change if other comparison
  tmp <- contrasts.fit(fit, contr)
  tmp <- eBayes(tmp)
  top.table <- topTable(tmp, sort.by = "P", n = Inf)
  head(top.table, 20)
  length(which(top.table$adj.P.Val < 0.05 & abs(top.table$logFC) > log2(1.5))) #how many are significant?
  res_limma <- as.data.frame(top.table)
  res_limma$change <- with(res_limma, ifelse(res_limma$adj.P.Val < 0.05 & abs(res_limma$logFC) > log2(1.5),
                                             ifelse(res_limma$logFC > log2(1.5),'UP','DOWN'),'NOT'))
  table(res_limma$change)
  res_limma = res_limma[res_limma$change %in% c('UP', 'DOWN'), ]
  resMeR_list[[c]] <- list(res_DESeq2, res_limma)
  countMatrix_list[[c]] <- countMatrix
  resPeak_list[[c]] <- merge_peak
  name = paste0(c, ' vs WT')
  name = paste0(gsub(' ','_', name ))
  names(resMeR_list[[c]] ) <- paste0(name, '_', c('DESeq2', 'limma'))
}
names(resMeR_list) <- NULL
RMdem_list[[i]] <- as.list(unlist(resMeR_list, recursive = FALSE))
names(RMdem_list[[i]] ) <- paste0(cl, '_', names(RMdem_list[[i]]))
RMdep_list[[i]] <- resPeak_list 
name <- unique(gsub('(_limma|_DESeq2)', '', names(RMdem_list[[i]] )))
names(RMdep_list[[i]] ) <- name
contrast_IP_list[[i]] <- lapply(countMatrix_list, extract_names)
condition_IP_list[[i]] <- condition_list

############
n = 81
############
#GSE161304# without case RIP, with RNA
############
print(RMPdeg_datasets[n,])#
i = RMPdeg_datasets[n,'GSE']
cl = c('LNCaP and RWPE')#
cl = 'LNCaP'
RMp = c('METTL3')
RM = 'm6A'
condition_list = list(LNCaP.DOX.0_METTL3KD = c(rep('METTL3KD', 6), rep('WT', 3)),#
                      LNCaP.DOX.200_METTL3KD = c(rep('METTL3KD', 6), rep('WT', 3)),#
                      LNCaP.DMSO_METTL3KD = c(rep('METTL3KD', 6), rep('WT', 3)),#
                      LNCaP.ENZ_METTL3KD = c(rep('METTL3KD', 6), rep('WT', 2))#
)

RMP_allDat[grep('GSE161304', RMP_allDat$study_alias )
           ,c("study_alias", "sample_alias", "sample_title")]
# View(allDat[[i]])
files <- list.files(file.path(data_directory,i))
files
getwd()
setwd(data_directory)
setwd(i)
###################################
#   more than 1 condition here   #
###################################
resInP_list <- list()
## Gene name
df_in <- fread(file.path(data_directory,i, 'GSE161301_rna-seq_STAR_18sample_normCounts.txt.gz'))
library(data.table)
countGTF=as.data.frame(df_in[, -1])
##Gene count matrix
# countGTF=as.matrix(ceiling(df_in[, -1]))
## Gene name
table(duplicated(df_in$geneSymbol))
dim(countGTF) #25343    12
## Sample name
colnames(countGTF)
## Sample information
sample_meta=data.frame(
  sample=colnames(countGTF),
  condition=c(rep("LNCaP.DOX.0_WT",3), rep("LNCaP.DOX.200_WT",3),
              rep("LNCaP.DOX.0_METTL3KD",3), rep("LNCaP.DOX.200_METTL3KD",3),
               rep("LNCaP.DOX.0_METTL3KD",3), rep("LNCaP.DOX.200_METTL3KD",3))
)
#####condition_LNCaP.DOX.200_METTL3KD_vs_WT#####
count <- as.matrix(ceiling(countGTF[, grepl('LNCaP.DOX.200_', sample_meta$condition)]))
table(is.na(count))
sample <- sample_meta[sample_meta$sample %in% colnames(count), ]
## SummarizedExperiment input
condition = c(rep('LNCaP.DOX.200_METTL3KD', 6), rep('LNCaP.DOX.200_WT', 3))
dds<-DESeqDataSetFromMatrix(countData=count, colData=sample,design=~condition)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 10
table(keep)
dds <- dds[keep,]
## Note on factor levels
dds$condition <- factor(dds$condition, levels = c("LNCaP.DOX.200_METTL3KD","LNCaP.DOX.200_WT"))
dds$condition <- relevel(dds$condition, ref = "LNCaP.DOX.200_WT")
dds <- DESeq(dds)
res <-  results(dds, contrast=c("condition",unique(condition))) # c("condition", rev(levels(genotype)))
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
table(res$change)
res$geneSymbol <- df_in$V1[keep]
resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
table(duplicated(resInP$geneSymbol))
resInP_list[['LNCaP.DOX.200_METTL3KD_vs_WT']] <- resInP

df_in <- fread(file.path(data_directory,i, 'GSE161302_allSample_Log2NormCounts.txt.gz'))
library(data.table)
countGTF=as.matrix(2^(df_in[, -c(1,19)]))
rownames(countGTF) <- df_in$gene_name
# count=2^count-1
# countGTF[is.infinite(countGTF)] <- NA
# [, !grepl('V1|gene_name', colnames(df_in))])
##Gene count matrix
## Gene name
## Sample name
colnames(countGTF)
## Sample information
sample_meta=data.frame(
  sample=colnames(countGTF),
  condition=c(rep("LNCaP.DMSO_WT",3), rep("LNCaP.ENZ_WT",3),
              rep("LNCaP.DMSO_METTL3KD",3), rep("LNCaP.ENZ_METTL3KD",3),
              rep("LNCaP.DMSO_METTL3KD",3), rep("LNCaP.ENZ_METTL3KD",2))
)
countGTF <- na.omit(ceiling(countGTF))
#####condition_LNCaP.DMSO_METTL3KD_vs_WT#####
count <- as.matrix(countGTF[, grepl('LNCaP.DMSO_', sample_meta$condition)])
table(is.na(count))
sample <- sample_meta[sample_meta$sample %in% colnames(count), ]
## SummarizedExperiment input
condition = c(rep('LNCaP.DMSO_METTL3KD', 6), rep('LNCaP.DMSO_WT', 3))
dds<-DESeqDataSetFromMatrix(countData=count, colData=sample,design=~condition)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 10
table(keep)
dds <- dds[keep,]
## Note on factor levels
dds$condition <- factor(dds$condition, levels = c("LNCaP.DMSO_METTL3KD","LNCaP.DMSO_WT"))
dds$condition <- relevel(dds$condition, ref = "LNCaP.DMSO_WT")
dds <- DESeq(dds)
res <-  results(dds, contrast=c("condition",unique(condition))) # c("condition", rev(levels(genotype)))
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
table(res$change)
res$geneSymbol <- df_in$V1[keep]
resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
table(duplicated(resInP$geneSymbol))
resInP_list[['LNCaP.DMSO_METTL3KD_vs_WT']] <- resInP

#####condition_LNCaP.ENZ_METTL3KD_vs_WT#####
count <- as.matrix(countGTF[, grepl('LNCaP.ENZ_', sample_meta$condition)])
table(is.na(count))
sample <- sample_meta[sample_meta$sample %in% colnames(count), ]
## SummarizedExperiment input
condition = c(rep('LNCaP.ENZ_METTL3KD', 5), rep('LNCaP.ENZ_WT', 3))
dds<-DESeqDataSetFromMatrix(countData=count, colData=sample,design=~condition)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 10
table(keep)
dds <- dds[keep,]
## Note on factor levels
dds$condition <- factor(dds$condition, levels = c("LNCaP.ENZ_METTL3KD","LNCaP.ENZ_WT"))
dds$condition <- relevel(dds$condition, ref = "LNCaP.ENZ_WT")
dds <- DESeq(dds)
res <-  results(dds, contrast=c("condition",unique(condition))) # c("condition", rev(levels(genotype)))
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
table(res$change)
res$geneSymbol <- df_in$V1[keep]
resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
table(duplicated(resInP$geneSymbol))
resInP_list[['LNCaP.ENZ_METTL3KD_vs_WT']] <- resInP

RMdeg_list[[i]] <- resInP_list
contrast_input_list[[i]] <- NULL#lapply(countGTF_list, extract_names)
condition_input_list[[i]] <- condition_list


############
n = 82
############
#GSE134103# with RIP and Input
############
print(RMPdeg_datasets[n,])#
i = RMPdeg_datasets[n,'GSE']
cl = c('Hela')#
RMp = c('METTL3')
RM = 'm6A'
condition = c(rep('METTL3KD', 1), rep('WT', 1))#

RMP_allDat[grep('GSE134103', RMP_allDat$study_alias )
           ,c("study_alias", "sample_alias", "sample_title")]
# View(allDat[[i]])
ENA <- as.data.frame(RMP_allDat[grepl('GSE134103', RMP_allDat$study_alias ) & ( !grepl('_miclip', RMP_allDat$sample_title) ) 
                                , c("study_alias", "sample_alias", "sample_title")])
ENA <- ENA[order(ENA$sample_alias),]
ENA
paste0(ENA$sample_alias, collapse = "|")

files <- list.files(file.path(data_directory,i))
files
getwd()
setwd(data_directory)
setwd(i)
###################################
#  Cell number-normalized (CNN)   #
#      MeRIP INPUT analysis       #
###################################
c <- unique(condition[!condition %in% 'WT'])
countGTF <- read.table( paste0("metaPlotR/Input_counts_", c, ".txt"),header=T)
dim(countGTF)
colnames(countGTF)
rownames(countGTF) <- countGTF$Geneid
countGTF <- countGTF[, !colnames(countGTF) %in% c('Geneid','Chr','Start','End','Strand','Length')]
dim(countGTF) # 28265     4
# 2. Differential gene analysis with edgeR
# Without replicates, you cannot estimate which genes are differentially expressed using EdgeR or DESeq2. You can only calculate fold changes based on normalized read counts (preferentially CPM normalized by TMM method included in any EdgeR analysis) and apply a stringent fold-change cut-off to determine which genes are more or less expressed depending on the condition.
y <- DGEList(countGTF)
# Calculate fold change
y <- calcNormFactors(y, method = c('TMM'))
y$samples$lib.size <- colSums(y$counts)
logcpm <- cpm(y, prior.count=2, log=TRUE)
logcpm[1:5,]
cpm <- cpm(y, prior.count=2, log=FALSE)
cpm[1:5, ]
fold_change <- cpm[,1]/cpm[,2]
res <- countGTF
res$logFC <- log2(fold_change)
res$FoldChange <- fold_change
res$change <- with(res, ifelse(abs(res$logFC) > log2(1.5),
                               ifelse(res$logFC > log2(1.5),'UP','DOWN'),'NOT'))
table(res$change)
resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
resInP$geneSymbol <- rownames(resInP)
table(resInP$change)
table(resInP$FoldChange %in% c(0, Inf))
# resInP <- resInP[! resInP$FoldChange %in% c(0, Inf), ]
name = paste0(c, ' vs WT')
name = paste0(cl, '_', name)
RMdeg_list[[i]] <- list( resInP)
names(RMdeg_list[[i]] ) <- paste0(gsub(' ','_',(name)))

countGTF_list <- list(colnames(countGTF)); names(countGTF_list) <- c
condition_list <- list(condition); names(condition_list) <- c
contrast_input_list[[i]] <- countGTF_list
condition_input_list[[i]] <- condition_list

###################################
#  Cell number-normalized (CNN)   #
#   MeRIP RIP analysis - Part1    #
###################################
## Sample information
merge_peak <- read.delim(paste0("metaPlotR/", c, "_newPeakFile.txt"),header=F)
countMatrix <- read.table(paste0("metaPlotR/m6A_counts_", c, ".txt"),header=T)
table(countMatrix$Geneid == merge_peak[-1,1])
table(duplicated(merge_peak$V1[-1]))
dim(countMatrix)
colnames(countMatrix)
rownames(countMatrix) <- countMatrix$Geneid
countMatrix <- countMatrix[, !colnames(countMatrix) %in% c('Geneid','Chr','Start','End','Strand','Length')]
# 2. Differential gene analysis with edgeR
# Without replicates, you cannot estimate which genes are differentially expressed using EdgeR or DESeq2. You can only calculate fold changes based on normalized read counts (preferentially CPM normalized by TMM method included in any EdgeR analysis) and apply a stringent fold-change cut-off to determine which genes are more or less expressed depending on the condition.
y <- DGEList(countMatrix)
# Calculate fold change
y <- calcNormFactors(y, method = c('TMM'))
y$samples$lib.size <- colSums(y$counts)
logcpm <- cpm(y, prior.count=2, log=TRUE)
logcpm[1:5,]
cpm <- cpm(y, prior.count=2, log=FALSE)
cpm[1:5, ]
# fold_change <- logcpm[,'Ratio1']/logcpm[,'Ratio2']
fold_change <- cpm[,1]/cpm[,2]
resMeR <- countMatrix
resMeR$logFC <- log2(fold_change)
resMeR$FoldChange <- fold_change
resMeR$change <- with(resMeR, ifelse(abs(resMeR$logFC) > log2(1.5),
                                     ifelse(resMeR$logFC > log2(1.5),'UP','DOWN'),'NOT'))
table(resMeR$change)
res_limma <- resMeR[resMeR$change %in% c('UP','DOWN'),]
RMdem_list[[i]] <- list(res_limma)
name = paste0(c, ' vs WT')
name = paste0(gsub(' ','_', name ))
names(RMdem_list[[i]]) <- paste0(cl, '_', name, '_limma')
RMdep_list[[i]] <- list( merge_peak )
names(RMdep_list[[i]] ) <-  paste0(cl, '_', name)

countMatrix_list <- list(colnames(countMatrix)); names(countMatrix_list) <- c
condition_list <- list(condition); names(condition_list) <- c
contrast_IP_list[[i]] <- countMatrix_list
condition_IP_list[[i]] <- condition_list


# ############
# n = 83
# ############
# #GSE163491# without case RIP, with RNA
# ############
# print(RMPdeg_datasets[n,])# skip?
# i = RMPdeg_datasets[n,'GSE']
# cl = c('HEK293T')#
# RMp = c('METTL3')
# RM = 'm6A'
# condition = c(rep('METTL3KO', 2), rep('WT', 2))#
# 
# RMP_allDat[grep('GSE163491', RMP_allDat$study_alias )
#            ,c("study_alias", "sample_alias", "sample_title")]
# # View(allDat[[i]])
# 
# files <- list.files(file.path(data_directory,i))
# files
# getwd()
# setwd(data_directory)
# setwd(i)
# 

############
n = 84
############
#GSE162357# with RIP and Input
############
print(RMPdeg_datasets[n,])# hg38
i = RMPdeg_datasets[n,'GSE']
cl = c('HeLa and HEK293 and HepG2 and HSPC')#
cl = 'Hela' # use m6A-seq data from GSE128443 as WT
RMp = c('METTL3')
RM = 'm6A'
condition_list = list( FTOKD = c(rep('FTOKD', 2), rep('WT', 3)),
                       FTOOE = c(rep('FTOOE', 2), rep('WT', 3)))
# #only the intersection of the two comparisons, FTO- vs Input, and FTO- vs FTO+, were kept for further analysis
RMP_allDat[grep('GSE162356', RMP_allDat$study_alias )
           ,c("study_alias", "sample_alias", "sample_title")]
# View(allDat[[i]])
RMP_allDat[grepl('GSE162356', RMP_allDat$study_alias ) & grepl('HeLa', RMP_allDat$sample_title)
           ,c("study_alias", "sample_alias", "sample_title")]
skip = RMP_allDat[grepl('GSE162356', RMP_allDat$study_alias ) & grepl('HeLa', RMP_allDat$sample_title)
                  ,c("study_alias", "sample_alias", "sample_title")]
paste0(skip$sample_alias, collapse = "|")

skip <- RMP_allDat[grepl('GSE162356', RMP_allDat$study_alias ) & grepl('HeLa', RMP_allDat$sample_title)
                   & (! grepl('RIP_polyA_FTO', RMP_allDat$sample_title))
                   , "sample_alias"]
paste0(skip$sample_alias, collapse = " ")
skip <- RMP_allDat[grepl('GSE162356', RMP_allDat$study_alias ) & grepl('HeLa', RMP_allDat$sample_title)
                   &  grepl('RIP_polyA_FTO', RMP_allDat$sample_title)
                   , "sample_alias"]
paste0(skip$sample_alias, collapse = " ")
ENA <- as.data.frame(RMP_allDat[grepl('GSE162356', RMP_allDat$study_alias ) & grepl('HeLa', RMP_allDat$sample_title)
                                , c("study_alias", "sample_alias", "sample_title")])
ENA <- ENA[order(ENA$sample_alias),]
ENA
paste0(ENA$sample_alias, collapse = "|")

files <- list.files(file.path(data_directory,i))
files
getwd()
setwd(data_directory)
setwd(i)

###################################
#  Cell number-normalized (CNN)   #
#      MeRIP INPUT analysis       #
###################################
resInP_list <- vector("list", length(condition_list))
names(resInP_list) <- names(condition_list)
countGTF_list <- vector("list", length(condition_list))
names(countGTF_list) <- names(condition_list)
for(c in names(condition_list)){
  countGTF <- read.table( paste0("metaPlotR/Input_counts_", c, ".txt"),header=T)
  dim(countGTF)
  colnames(countGTF)
  rownames(countGTF) <- countGTF$Geneid
  countGTF <- countGTF[, !colnames(countGTF) %in% c('Geneid','Chr','Start','End','Strand','Length')]
  dim(countGTF) # 28265     4
  condition <- condition_list[[c]]
  ## Sample information
  sample=data.frame(
    sample=colnames(countGTF),
    condition=condition)
  ## SummarizedExperiment input
  dds<-DESeqDataSetFromMatrix(countData=countGTF,colData=sample,design=~condition)
  ## Pre-filtering
  keep <- rowSums(counts(dds)) >= 10
  table(keep)
  dds <- dds[keep,]
  ## Note on factor levels
  dds$condition <- factor(dds$condition, levels = unique(condition))
  dds$condition <- relevel(dds$condition, ref = "WT")
  # Differential expression analysis-------------------------------------------
  dds <- DESeq(dds)
  res <-  results(dds, contrast=c("condition",unique(condition)))
  res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
  length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
  res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                                 ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
  table(res$change)
  table(is.na(res$change))
  # resOrdered <- as.data.frame(res[order(res$padj), ])
  resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
  resInP$geneSymbol <- rownames(resInP)
  table(resInP$change)
  resInP_list[[c]] <- resInP
  countGTF_list[[c]] <- countGTF 
}
name = paste0(cl, '_',names(resInP_list), ' vs WT')
RMdeg_list[[i]] <- resInP_list 
names(RMdeg_list[[i]] ) <- paste0(gsub(' ','_',(name)))
contrast_input_list[[i]] <- lapply(countGTF_list, extract_names)
condition_input_list[[i]] <- condition_list

###################################
#  Cell number-normalized (CNN)   #
#   MeRIP RIP analysis - Part1    #
###################################
resMeR_list <- vector("list", length(condition_list))
names(resMeR_list) <- names(condition_list)
countMatrix_list <- vector("list", length(condition_list))
names(countMatrix_list) <- names(condition_list)
resPeak_list <- vector("list",length(condition_list))
names(resPeak_list) <- names(condition_list)
for(c in names(condition_list)){
  merge_peak <- read.delim(paste0("metaPlotR/", c, "_newPeakFile.txt"),header=F)
  countMatrix <- read.table(paste0("metaPlotR/m6A_counts_", c, ".txt"),header=T)
  table(countMatrix$Geneid == merge_peak[-1,1])
  table(duplicated(merge_peak[-1,1]))
  dim(countMatrix)
  colnames(countMatrix)
  rownames(countMatrix) <- countMatrix$Geneid
  countMatrix <- countMatrix[, !colnames(countMatrix) %in% c('Geneid','Chr','Start','End','Strand','Length')]
  # 1. 做好分组因子即可
  condition <- condition_list[[c]]
  groupList<-factor(condition)
  colData <- data.frame(row.names=colnames(countMatrix), groupList=groupList)
  
  # 2. Exploring Method 1: DESeq2-------------------------------------------
  # 构建 DESeqDataSet 对象
  dds <- DESeqDataSetFromMatrix(countData = countMatrix,
                                colData = colData,
                                design = ~ groupList)
  ## Pre-filtering
  keep <- rowSums(counts(dds)) >= 1
  table(keep)
  dds <- dds[keep,]
  # 接着计算每个样本的 library size 用于后续的标准化。library size 等于落在该样本所有peaks上的reads的总数，即 counts 矩阵中每列的加和。如果bFullLibrarySize设为TRUE，则会使用library的总reads数（根据原BAM/BED文件统计获得）。然后使用estimateDispersions估计统计分布，需要将参数fitType设为’local’
  dds<-estimateSizeFactors(dds)
  dds<-estimateDispersions(dds,fitType="local")
  dds <- nbinomWaldTest(dds)
  res <- results(dds)
  ## 提取你想要的差异分析结果，我们这里是treated组对untreated组进行比较
  res <-  results(dds, contrast=c("groupList", unique(condition)))
  res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
  length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
  res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                                 ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
  res_DESeq2 = as.data.frame(res)
  table(res_DESeq2$change)
  res_DESeq2 = res_DESeq2[res_DESeq2$change %in% c('UP', 'DOWN'), ]
  # 2. Exploring Method 2: limma-------------------------------------------
  library(edgeR)
  d0 <- DGEList(counts = countMatrix)
  dim(d0)
  # Filtering low count peaks
  cutoff <- 1           
  drop <- which(apply(cpm(d0), 1, max) < cutoff)
  print(length(drop))
  if (length(drop) >0) {
    d <- d0[-drop,] 
  } else{ d <- d0}
  dim(d) #number of peaks left: 8501      
  d <- calcNormFactors(d, method = c('TMM'))
  colnames(countMatrix)
  genotype <- groupList
  table(genotype)
  mm <- model.matrix(~0 + genotype)#
  v <- voom(d, mm, plot=F, normalize="quantile")#voom from limma:Transform count data to log2-counts per million (logtpm)
  fit <- lmFit(v, mm)
  head(coef(fit))
  contr <- makeContrasts(paste0(colnames(coef(fit)), collapse  = ' - '), levels = colnames(coef(fit))) #Change if other comparison
  tmp <- contrasts.fit(fit, contr)
  tmp <- eBayes(tmp)
  top.table <- topTable(tmp, sort.by = "P", n = Inf)
  head(top.table, 20)
  length(which(top.table$adj.P.Val < 0.05 & abs(top.table$logFC) > log2(1.5))) #how many are significant?
  res_limma <- as.data.frame(top.table)
  res_limma$change <- with(res_limma, ifelse(res_limma$adj.P.Val < 0.05 & abs(res_limma$logFC) > log2(1.5),
                                             ifelse(res_limma$logFC > log2(1.5),'UP','DOWN'),'NOT'))
  table(res_limma$change)
  res_limma = res_limma[res_limma$change %in% c('UP', 'DOWN'), ]
  resMeR_list[[c]] <- list(res_DESeq2, res_limma)
  countMatrix_list[[c]] <- countMatrix
  resPeak_list[[c]] <- merge_peak
  name = paste0(c, ' vs WT')
  name = paste0(gsub(' ','_', name ))
  names(resMeR_list[[c]] ) <- paste0(name, '_', c('DESeq2', 'limma'))
}
names(resMeR_list) <- NULL
RMdem_list[[i]] <- as.list(unlist(resMeR_list, recursive = FALSE))
names(RMdem_list[[i]] ) <- paste0(cl, '_', names(RMdem_list[[i]]))
RMdep_list[[i]] <- resPeak_list 
name <- unique(gsub('(_limma|_DESeq2)', '', names(RMdem_list[[i]] )))
names(RMdep_list[[i]] ) <- name
contrast_IP_list[[i]] <- lapply(countMatrix_list, extract_names)
condition_IP_list[[i]] <- condition_list

############
n = 85
############
#GSE66012# with miCLIP and RIP
############
print(RMPdeg_datasets[n,])#
i = RMPdeg_datasets[n,'GSE']
cl = c('HEK293')#
RMp = c('NSUN3')
RM = 'm5C'
condition_list <- list(NSUN3MU = c(rep('NSUN3MU', 3), rep('WT', 2)),
                       NSUN6MU = c(rep('NSUN6MU', 3), rep('WT', 2)))
RMP_allDat[grep('GSE66011', RMP_allDat$study_alias )
           ,c("study_alias", "sample_alias", "sample_title")]
# View(allDat[[i]])
ENA <- as.data.frame(RMP_allDat[grepl('GSE66011', RMP_allDat$study_alias )
                                ,c("study_alias", "sample_alias", "sample_title")])
ENA <- ENA[order(ENA$sample_alias),]
skip <- RMP_allDat[grepl('GSE66011', RMP_allDat$study_alias ) & grepl('miCLIP', RMP_allDat$sample_title)
                   , "sample_alias"]
paste0(skip$sample_alias, collapse = "|")
# View(allDat[[i]])
y <- ENA[grepl('miCLIP_NSUN3|BS-seq',ENA$sample_title),]
y <- ENA[grepl('miCLIP_NSUN6|BS-seq',ENA$sample_title),]
y
x <- paste0('$outdir/GSE66011/', y$sample_alias, '_IP.bam')
paste(x, collapse = ' ')

files <- list.files(file.path(data_directory,i))
files
getwd()
setwd(data_directory)
setwd(i)
countfiles <- list.files('./metaPlotR/')
countfiles <- countfiles[grep('counts', countfiles)]
countfiles <- countfiles[! grepl('summary', countfiles)]
countfiles
###################################
#         RIP-seq analysis        #
###################################
resMeR_list <- vector("list", length(condition_list))
names(resMeR_list) <- names(condition_list)
countMatrix_list <- vector("list", length(condition_list))
names(countMatrix_list) <- names(condition_list)
resPeak_list <- vector("list",length(condition_list))
names(resPeak_list) <- names(condition_list)
for(c in names(condition_list)){
  merge_peak <- read.delim(paste0("metaPlotR/", c, "_newPeakFile.txt"),header=F)
  countMatrix <- read.table(paste0("metaPlotR/m5C_counts_", c, ".txt"),header=T)
  table(countMatrix$Geneid == merge_peak[-1,1])
  table(duplicated(merge_peak[-1,1]))
  dim(countMatrix)
  colnames(countMatrix)
  rownames(countMatrix) <- countMatrix$Geneid
  countMatrix <- countMatrix[, !colnames(countMatrix) %in% c('Geneid','Chr','Start','End','Strand','Length')]
  # 1. 做好分组因子即可
  condition <- condition_list[[c]]
  groupList <- factor(condition)
  colData <- data.frame(row.names=colnames(countMatrix), groupList=groupList)
  # 2. Exploring Method 1: DESeq2-------------------------------------------
  # 构建 DESeqDataSet 对象
  dds <- DESeqDataSetFromMatrix(countData = countMatrix,
                                colData = colData,
                                design = ~ groupList)
  ## Pre-filtering
  keep <- rowSums(counts(dds)) >= 1
  table(keep)
  dds <- dds[keep,]
  if(c=='NSUN6MU'){
    cts <- counts(dds)
    # geoMeans <- apply(cts, 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))
    geoMeans <- apply(cts, 1, function(row) if (all(row == 0)) 0 else exp(sum(log(row[row != 0]))/length(row)))
    dds <- estimateSizeFactors(dds, geoMeans=geoMeans)
  } else{
    dds<-estimateSizeFactors(dds)
  }
  # 接着计算每个样本的 library size 用于后续的标准化。library size 等于落在该样本所有peaks上的reads的总数，即 counts 矩阵中每列的加和。如果bFullLibrarySize设为TRUE，则会使用library的总reads数（根据原BAM/BED文件统计获得）。然后使用estimateDispersions估计统计分布，需要将参数fitType设为’local’
  dds<-estimateDispersions(dds,fitType="local")
  dds <- nbinomWaldTest(dds)
  res <- results(dds)
  ## 提取你想要的差异分析结果，我们这里是treated组对untreated组进行比较
  res <-  results(dds, contrast=c("groupList", unique(condition)))
  res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
  length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
  res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                                 ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
  res_DESeq2 = as.data.frame(res)
  table(res_DESeq2$change)
  res_DESeq2 = res_DESeq2[res_DESeq2$change %in% c('UP', 'DOWN'), ]
  # 2. Exploring Method 2: limma-------------------------------------------
  library(edgeR)
  d0 <- DGEList(counts = countMatrix)
  dim(d0)
  # Filtering low count peaks
  cutoff <- 1
  drop <- which(apply(cpm(d0), 1, max) < cutoff)
  print(length(drop))
  if (length(drop) >0) {
    d <- d0[-drop,]
  } else{ d <- d0}
  dim(d) #number of peaks left: 8501
  d <- calcNormFactors(d, method = c('TMM'))
  colnames(countMatrix)
  genotype <- groupList
  table(genotype)
  mm <- model.matrix(~0 + genotype)#
  v <- voom(d, mm, plot=F, normalize="quantile")#voom from limma:Transform count data to log2-counts per million (logtpm)
  fit <- lmFit(v, mm)
  head(coef(fit))
  contr <- makeContrasts(paste0(colnames(coef(fit)), collapse  = ' - '), levels = colnames(coef(fit))) #Change if other comparison
  tmp <- contrasts.fit(fit, contr)
  tmp <- eBayes(tmp)
  top.table <- topTable(tmp, sort.by = "P", n = Inf)
  head(top.table, 20)
  length(which(top.table$adj.P.Val < 0.05 & abs(top.table$logFC) > log2(1.5))) #how many are significant?
  res_limma <- as.data.frame(top.table)
  res_limma$change <- with(res_limma, ifelse(res_limma$adj.P.Val < 0.05 & abs(res_limma$logFC) > log2(1.5),
                                             ifelse(res_limma$logFC > log2(1.5),'UP','DOWN'),'NOT'))
  table(res_limma$change)
  res_limma = res_limma[res_limma$change %in% c('UP', 'DOWN'), ]
  
  resMeR_list[[c]] <- list(res_DESeq2, res_limma)
  countMatrix_list[[c]] <- countMatrix
  resPeak_list[[c]] <- merge_peak
  name = res@elementMetadata$description[2]
  name = paste0(cl, '_',str_split(name, "groupList ", simplify = T)[,2])
  name = paste0(gsub(' ','_',(name)))
  names(resMeR_list[[c]] ) <- paste0(name, '_', c('DESeq2', 'limma'))
}
# RMdem_list[[i]] <- as.list(unlist(resMeR_list, recursive = FALSE))
# name <- str_split(names(RMdem_list[[i]] ), "[.]", simplify = T)[,2]
# names(RMdem_list[[i]] ) <- name
names(resMeR_list) <- NULL
RMdem_list[[i]] <- as.list(unlist(resMeR_list, recursive = FALSE))
names(RMdem_list[[i]] ) #<- paste0(cl, '_', names(RMdem_list[[i]]))
RMdep_list[[i]] <- resPeak_list 
name <- unique(gsub('(_limma|_DESeq2)', '', names(RMdem_list[[i]] )))
names(RMdep_list[[i]] ) <- name
contrast_IP_list[[i]] <- lapply(countMatrix_list, extract_names)
condition_IP_list[[i]] <- condition_list

############
n = 86
############
#GSE122413# with RIP without Input
############
print(RMPdeg_datasets[n,])# 
i = RMPdeg_datasets[n,'GSE']
cl = c('HEK293')# use RNA-seq data from GSE44384 as Input (despite of NSUN2-KD)
RMp = c('NSUN2')
RM = 'm5C'
condition = c(rep('NSUN2KO', 1), rep('WT', 1))#
RMP_allDat[grep('GSE122413|GSE44384', RMP_allDat$study_alias )
           ,c("study_alias", "sample_alias", "sample_title")]
# View(allDat[[i]])
files <- list.files(file.path(data_directory,i))
files
getwd()
setwd(data_directory)
setwd(i)
DEP <- read_excel(file.path(data_directory,i, 'GSE122413_Methylated_RNA_sites.xlsx'),sheet = 2)
###################################
#  Cell number-normalized (CNN)   #
#      MeRIP INPUT analysis       #
###################################
c <- unique(condition[!condition %in% 'WT'])
countGTF <- read.table( paste0("metaPlotR/Input_counts_", c, ".txt"),header=T)
dim(countGTF)
colnames(countGTF)
rownames(countGTF) <- countGTF$Geneid
countGTF <- countGTF[, !colnames(countGTF) %in% c('Geneid','Chr','Start','End','Strand','Length')]
dim(countGTF) # 28265     4
# 2. Differential gene analysis with edgeR
# Without replicates, you cannot estimate which genes are differentially expressed using EdgeR or DESeq2. You can only calculate fold changes based on normalized read counts (preferentially CPM normalized by TMM method included in any EdgeR analysis) and apply a stringent fold-change cut-off to determine which genes are more or less expressed depending on the condition.
y <- DGEList(countGTF)
# Calculate fold change
y <- calcNormFactors(y, method = c('TMM'))
y$samples$lib.size <- colSums(y$counts)
logcpm <- cpm(y, prior.count=2, log=TRUE)
logcpm[1:5,]
cpm <- cpm(y, prior.count=2, log=FALSE)
cpm[1:5, ]
fold_change <- cpm[,1]/cpm[,2]
res <- countGTF
res$logFC <- log2(fold_change)
res$FoldChange <- fold_change
res$change <- with(res, ifelse(abs(res$logFC) > log2(1.5),
                               ifelse(res$logFC > log2(1.5),'UP','DOWN'),'NOT'))
table(res$change)
resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
resInP$geneSymbol <- rownames(resInP)
table(resInP$change)
table(resInP$FoldChange %in% c(0, Inf))
# resInP <- resInP[! resInP$FoldChange %in% c(0, Inf), ]
name = paste0(c, ' vs WT')
name = paste0(cl, '_', name)
RMdeg_list[[i]] <- list( resInP)
names(RMdeg_list[[i]] ) <- paste0(gsub(' ','_',(name)))

countGTF_list <- list(colnames(countGTF)); names(countGTF_list) <- c
condition_list <- list(condition); names(condition_list) <- c
contrast_input_list[[i]] <- countGTF_list
condition_input_list[[i]] <- condition_list

###################################
#  Cell number-normalized (CNN)   #
#   MeRIP RIP analysis - Part1    #
###################################
## Sample information
merge_peak <- read.delim(paste0("metaPlotR/", c, "_newPeakFile.txt"),header=F)
countMatrix <- read.table(paste0("metaPlotR/m5C_counts_", c, ".txt"),header=T)
table(countMatrix$Geneid == merge_peak[-1,1])
table(duplicated(merge_peak$V1[-1]))
dim(countMatrix)
colnames(countMatrix)
rownames(countMatrix) <- countMatrix$Geneid
countMatrix <- countMatrix[, !colnames(countMatrix) %in% c('Geneid','Chr','Start','End','Strand','Length')]
# 2. Differential gene analysis with edgeR
# Without replicates, you cannot estimate which genes are differentially expressed using EdgeR or DESeq2. You can only calculate fold changes based on normalized read counts (preferentially CPM normalized by TMM method included in any EdgeR analysis) and apply a stringent fold-change cut-off to determine which genes are more or less expressed depending on the condition.
y <- DGEList(countMatrix)
# Calculate fold change
y <- calcNormFactors(y, method = c('TMM'))
y$samples$lib.size <- colSums(y$counts)
logcpm <- cpm(y, prior.count=2, log=TRUE)
logcpm[1:5,]
cpm <- cpm(y, prior.count=2, log=FALSE)
cpm[1:5, ]
# fold_change <- logcpm[,'Ratio1']/logcpm[,'Ratio2']
fold_change <- cpm[,1]/cpm[,2]
resMeR <- countMatrix
resMeR$logFC <- log2(fold_change)
resMeR$FoldChange <- fold_change
resMeR$change <- with(resMeR, ifelse(abs(resMeR$logFC) > log2(1.5),
                                     ifelse(resMeR$logFC > log2(1.5),'UP','DOWN'),'NOT'))
table(resMeR$change)
res_limma <- resMeR[resMeR$change %in% c('UP','DOWN'),]
RMdem_list[[i]] <- list(res_limma)
name = paste0(c, ' vs WT')
name = paste0(gsub(' ','_', name ))
names(RMdem_list[[i]]) <- paste0(cl, '_', name, '_limma')
RMdep_list[[i]] <- list( merge_peak )
names(RMdep_list[[i]] ) <-  paste0(cl, '_', name)

countMatrix_list <- list(colnames(countMatrix)); names(countMatrix_list) <- c
condition_list <- list(condition); names(condition_list) <- c
contrast_IP_list[[i]] <- countMatrix_list
condition_IP_list[[i]] <- condition_list


############
n = 87
############
#GSE169589# without RIP, with RNA/Ribo
############
print(RMPdeg_datasets[n,])#
i = RMPdeg_datasets[n,'GSE']
cl = '5-8F'#c('Nasopharyngeal carcinoma')#
RMp = c('METTL1')
RM = 'm7G'
condition = c(rep('METTL1KO', 1), rep('WT', 1))#METTL1 was knocked out by sgRNA.
RMP_allDat[grep('GSE169589', RMP_allDat$study_alias )
           ,c("study_alias", "sample_alias", "sample_title")]
# View(allDat[[i]])
# setwd(data_directory)
# gse <- getGEO('GSE169589',destdir=".", AnnotGPL=F, getGPL=F)
# setwd('GSE169589')
# pdata <- pData(gse[[3]])
files <- list.files(file.path(data_directory,i))
files
getwd()
setwd(data_directory)
setwd(i)
# DEG <- read.table(file.path(data_directory,i, 'GSE169589_TRAC-seq_expression.txt.gz'))
# # TRAC-seq represents the ratio and cleavage score for every tRNA in each sample.
df_in <- fread(file.path(data_directory,i, 'GSE169589_Ribo-seq_reads.txt.gz'))
# DEG <- fread(file.path(data_directory,i, 'GSE169589_Polyribosome_mRNA_expression.txt.gz'))
## Gene name
library(data.table)
countGTF=as.data.frame(df_in[, c(8, 3, 6)])
countGTF <- aggregate(x = countGTF[,(2:ncol(countGTF))],   #此时exprSet的第三列开始是表达矩阵内容
                      by = list(symbol = countGTF$Symbol),   #按照相同symbol分组，在组内计算
                      FUN = max) %>%   #原文中是计算最大值（max），也可以计算平均值（mean）或者中位数（median）
  column_to_rownames(var = 'symbol')
dim(countGTF) # 51242    12
c <- unique(condition[!condition %in% 'WT'])
condition <- unique(condition)
# 2. Differential gene analysis with edgeR
# Without replicates, you cannot estimate which genes are differentially expressed using EdgeR or DESeq2. You can only calculate fold changes based on normalized read counts (preferentially CPM normalized by TMM method included in any EdgeR analysis) and apply a stringent fold-change cut-off to determine which genes are more or less expressed depending on the condition.
y <- DGEList(countGTF)
# Calculate fold change
y <- calcNormFactors(y, method = c('TMM'))
y$samples$lib.size <- colSums(y$counts)
logcpm <- cpm(y, prior.count=2, log=TRUE)
logcpm[1:5,]
cpm <- cpm(y, prior.count=2, log=FALSE)
cpm[1:5, ]
fold_change <- cpm[,1]/cpm[,2]
res <- countGTF
res$logFC <- log2(fold_change)
res$FoldChange <- fold_change
res$change <- with(res, ifelse(abs(res$logFC) > log2(1.5),
                               ifelse(res$logFC > log2(1.5),'UP','DOWN'),'NOT'))
table(res$change)
resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
resInP$geneSymbol <- rownames(resInP)
table(resInP$change)
table(resInP$FoldChange %in% c(0, Inf))
table(resInP$logFC %in% c(-Inf, Inf))
# resInP <- resInP[! resInP$FoldChange %in% c(0, Inf), ]
name = '5-8F_METTL1KO_vs_WT'
RMdeg_list[[i]] <- list( resInP)
names(RMdeg_list[[i]] ) <- name
condition_list <- list(condition); names(condition_list) <- c
contrast_input_list[[i]] <- NULL
condition_input_list[[i]] <- condition_list

############
n = 88
############
#GSE215095# without RIP, with RNA
############
print(RMPdeg_datasets[n,])#
i = RMPdeg_datasets[n,'GSE']
cl = '92.1' # 'Ocular melanoma cell'
RMp = c('METTL14')
RM = 'm6A'
condition = c(rep('METTL14OE', 3), rep('WT', 3))#
RMP_allDat[grep('GSE215095', RMP_allDat$study_alias )
           ,c("study_alias", "sample_alias", "sample_title")]
# View(allDat[[i]])
files <- list.files(file.path(data_directory,i))
files
getwd()
setwd(data_directory)
setwd(i)
DEG <- read.table(file.path(data_directory,i, 'GSE215095_All.HTSeq.counts.txt.gz'), header = T)
DEG1<- read.table(file.path(data_directory,i, 'GSE215095_All.HTSeq.counts1.txt.gz'), header = T)
countGTF <- left_join(DEG1,  DEG,   by = 'AccID')
countGTF <- aggregate(x = countGTF[,2:ncol(countGTF)],   #此时exprSet的第三列开始是表达矩阵内容
                      by = list(symbol = countGTF$AccID),   #按照相同symbol分组，在组内计算
                      FUN = max) %>%   #原文中是计算最大值（max），也可以计算平均值（mean）或者中位数（median）
  column_to_rownames(var = 'symbol')
dim(countGTF) # 51242    12
geneSymbol <- rownames(countGTF)
## Gene name
dim(countGTF) #25343    12
## Sample name
colnames(countGTF)
## Sample information
sample=data.frame(
  sample=colnames(countGTF),
  condition=c( 'WT', 'METTL14OE', rep(c('WT','METTL14OE'), each = 2)))
## SummarizedExperiment input
dds<-DESeqDataSetFromMatrix(countData=countGTF,colData=sample,design=~condition)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 10
table(keep)
dds <- dds[keep,]
## Note on factor levels
dds$condition <- factor(dds$condition, levels = unique(condition))
dds$condition <- relevel(dds$condition, ref = "WT")
# Differential expression analysis-------------------------------------------
dds <- DESeq(dds)
res <-  results(dds, contrast=c("condition",unique(condition)))
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
table(res$change)
table(is.na(res$change))
# resOrdered <- as.data.frame(res[order(res$padj), ])
resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
resInP$geneSymbol <- rownames(resInP)
c <- unique(condition[!condition %in% 'WT'])
name = res@elementMetadata$description[2]
name = paste0(cl, '_',str_split(name, "condition ", simplify = T)[,2])
RMdeg_list[[i]] <- list( resInP)
names(RMdeg_list[[i]] ) <- paste0(gsub(' ','_',(name)))
condition_list <- list(condition); names(condition_list) <- c
contrast_input_list[[i]] <- NULL
condition_input_list[[i]] <- condition_list

############
n = 89
############
#GSE210867# with RIP and RNA
############
print(RMPdeg_datasets[n,])#hg38
i = RMPdeg_datasets[n,'GSE']
cl = c('hESCs')# use m6A-seq data from GSE52600 as WT (H1A/B; Endoderm)
RMp = c('YTHDC2')
RM = 'm6A'
condition = c(rep('YTHDC2PD', 2), rep('WT', 2))#
RMP_allDat[grep('GSE210867', RMP_allDat$study_alias )
           ,c("study_alias", "sample_alias", "sample_title")]
# View(allDat[[i]])
files <- list.files(file.path(data_directory,i))
files
getwd()
setwd(data_directory)
setwd(i)
###################################
#       RNA/RIP-seq analysis      #
###################################
condition = c(rep('YTHDC2PD', 2), rep('WT', 8))#
df_in <- read.table(file.path(data_directory,i, 'metaPlotR/Input_counts_YTHDC2PD.txt'),header = T)
##Gene count matrix
countGTF <- df_in[, !colnames(df_in) %in% c('Geneid','Chr','Start','End','Strand','Length')]
rownames(countGTF) <- df_in$Geneid
dim(countGTF) #28278     4
colSums(countGTF)
table(colSums(countGTF)=="1e+06")
## Gene name
library(data.table)
table(duplicated(df_in$Hgnc_symbol))
table(duplicated(df_in$Geneid))
dim(countGTF) # 
## Sample name
colnames(countGTF)
## Sample information
sample=data.frame(
  sample=colnames(countGTF),
  condition=condition)
## SummarizedExperiment input
dds<-DESeqDataSetFromMatrix(countData=countGTF,colData=sample,design=~condition)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 10
table(keep)
dds <- dds[keep,]
## Note on factor levels
dds$condition <- factor(dds$condition, levels = unique(condition))
dds$condition <- relevel(dds$condition, ref = "WT")
# Differential expression analysis-------------------------------------------
dds <- DESeq(dds)
# res<-results(dds,name = "condition_TRMT61AOE_vs_WT")
res <-  results(dds, contrast=c("condition",unique(condition)))
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
table(res$change)
table(is.na(res$change))
# resOrdered <- as.data.frame(res[order(res$padj), ])
resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
resInP$geneSymbol <- rownames(resInP)
name = res@elementMetadata$description[2]
name = paste0(cl, '_',str_split(name, "condition ", simplify = T)[,2])
RMdeg_list[[i]] <- list( resInP)
names(RMdeg_list[[i]] ) <- paste0(gsub(' ','_',(name)))
c <- unique(condition[!condition %in% 'WT'])
countGTF_list <- list(colnames(countGTF)); names(countGTF_list) <- c
condition_list <- list(condition); names(condition_list) <- c
contrast_input_list[[i]] <- countGTF_list
condition_input_list[[i]] <- condition_list

###################################
#  Cell number-normalized (CNN)   #
#   MeRIP RIP analysis - Part1    #
###################################
condition = c(rep('YTHDC2PD', 2), rep('WT', 6))#
## Sample information
merge_peak <- read.delim(paste0("metaPlotR/", c, "_newPeakFile.txt"),header=F)
countMatrix <- read.table(paste0("metaPlotR/m6A_counts_", c, ".txt"),header=T)
table(countMatrix$Geneid == merge_peak[-1,1])
table(duplicated(merge_peak$V1[-1]))
dim(countMatrix)
colnames(countMatrix)
rownames(countMatrix) <- countMatrix$Geneid
countMatrix <- countMatrix[, !colnames(countMatrix) %in% c('Geneid','Chr','Start','End','Strand','Length')]
# 1. 做好分组因子即可
groupList<-factor(condition)
colData <- data.frame(row.names=colnames(countMatrix), groupList=groupList)
# 2. Exploring Method 1: DESeq2-------------------------------------------
# 构建 DESeqDataSet 对象
dds <- DESeqDataSetFromMatrix(countData = countMatrix,
                              colData = colData,
                              design = ~ groupList)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 1
table(keep)
dds <- dds[keep,]
# 接着计算每个样本的 library size 用于后续的标准化。library size 等于落在该样本所有peaks上的reads的总数，即 counts 矩阵中每列的加和。如果bFullLibrarySize设为TRUE，则会使用library的总reads数（根据原BAM/BED文件统计获得）。然后使用estimateDispersions估计统计分布，需要将参数fitType设为’local’
dds<-estimateSizeFactors(dds)
dds<-estimateDispersions(dds,fitType="local")
dds <- nbinomWaldTest(dds)
res <- results(dds)
## 提取你想要的差异分析结果，我们这里是treated组对untreated组进行比较
res <-  results(dds, contrast=c("groupList", unique(condition)))
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
res_DESeq2 = as.data.frame(res)
table(res_DESeq2$change)#Not only
res_DESeq2 = res_DESeq2[res_DESeq2$change %in% c('UP', 'DOWN'), ]
# 2. Exploring Method 2: limma-------------------------------------------
library(edgeR)
d0 <- DGEList(counts = countMatrix)
dim(d0)
# Filtering low count peaks
cutoff <- 1
drop <- which(apply(cpm(d0), 1, max) < cutoff)
print(length(drop))
if (length(drop) >0) {
  d <- d0[-drop,]
} else{ d <- d0}
dim(d) #number of peaks left: 8501
d <- calcNormFactors(d, method = c('TMM'))
colnames(countMatrix)
genotype <- groupList
table(genotype)
mm <- model.matrix(~0 + genotype)#
v <- voom(d, mm, plot=F, normalize="quantile")#voom from limma:Transform count data to log2-counts per million (logtpm)
fit <- lmFit(v, mm)
head(coef(fit))
contr <- makeContrasts(paste0(colnames(coef(fit)), collapse  = ' - '), levels = colnames(coef(fit))) #Change if other comparison
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)
length(which(top.table$adj.P.Val < 0.05 & abs(top.table$logFC) > log2(1.5))) #how many are significant?
res_limma <- as.data.frame(top.table)
res_limma$change <- with(res_limma, ifelse(res_limma$adj.P.Val < 0.05 & abs(res_limma$logFC) > log2(1.5),
                                           ifelse(res_limma$logFC > log2(1.5),'UP','DOWN'),'NOT'))
table(res_limma$change)#NOT only
res_limma = res_limma[res_limma$change %in% c('UP', 'DOWN'), ]

RMdem_list[[i]] <- list(res_DESeq2, res_limma)
name = res@elementMetadata$description[2]
name = paste0(cl, '_',str_split(name, "groupList ", simplify = T)[,2])
name= paste0(gsub(' ','_',(name)))
names(RMdem_list[[i]] ) <- paste0(name, '_', c('DESeq2', 'limma'))
RMdep_list[[i]] <- list( merge_peak )
names(RMdep_list[[i]] ) <- name
countMatrix_list <- list(colnames(countMatrix)); names(countMatrix_list) <- c
condition_list <- list(condition); names(condition_list) <- c
contrast_IP_list[[i]] <- countMatrix_list
condition_IP_list[[i]] <- condition_list

contrast_input_list[['GSE210867']]$YTHDC2PD <- gsub('X.data2.rluo4.EpiTrans.RMDatasets.GEO..GSE52600.metaPlotR.', '', contrast_input_list[['GSE210867']]$YTHDC2PD )
contrast_IP_list[['GSE210867']]$YTHDC2PD <- gsub('X.data2.rluo4.EpiTrans.RMDatasets.GEO..GSE52600.metaPlotR.', '', contrast_IP_list[['GSE210867']]$YTHDC2PD )

############
n = 90
############
#GSE202815# with RIP and Input
############
print(RMPdeg_datasets[n,])# hg38
i = RMPdeg_datasets[n,'GSE']
cl = c('Flp-In TRex 293 and HEK293')#
cl = 'HEK293'
RMp = c('ALKBH1')
RM = 'm5C'
condition = c(rep('ALKBH1KO', 2), rep('WT', 2))#
RMP_allDat[grep('GSE202814', RMP_allDat$study_alias )
           ,c("study_alias", "sample_alias", "sample_title")]
# View(allDat[[i]])
skip <- RMP_allDat[grepl('GSE202814', RMP_allDat$study_alias )
                   &  grepl('input', RMP_allDat$sample_title)
                   , "sample_alias"]
paste0(skip$sample_alias, collapse = " ")
skip <- RMP_allDat[grepl('GSE202814', RMP_allDat$study_alias )
                   &  grepl('ip', RMP_allDat$sample_title)
                   , "sample_alias"]
paste0(skip$sample_alias, collapse = " ")
ENA <- as.data.frame(RMP_allDat[grep('GSE202814', RMP_allDat$study_alias )
                                , c("study_alias", "sample_alias", "sample_title")])
ENA <- ENA[order(ENA$sample_alias),]
ENA
files <- list.files(file.path(data_directory,i))
files
getwd()
setwd(data_directory)
setwd(i)
c <- unique(condition[!condition %in% 'WT'])
###################################
#  Cell number-normalized (CNN)   #
#      MeRIP INPUT analysis       #
###################################
countGTF <- read.table( paste0("metaPlotR/Input_counts_", c, ".txt"),header=T)
rownames(countGTF) <- countGTF$Geneid
countGTF <- countGTF[, !colnames(countGTF) %in% c('Geneid','Chr','Start','End','Strand','Length')]
dim(countGTF) #28265     4
## Sample name
## Sample information
sample=data.frame(
  sample=colnames(countGTF),
  condition=condition)
## SummarizedExperiment input
dds<-DESeqDataSetFromMatrix(countData=countGTF,colData=sample,design=~condition)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 10
table(keep)
dds <- dds[keep,]
## Note on factor levels
dds$condition <- factor(dds$condition, levels = unique(condition))
dds$condition <- relevel(dds$condition, ref = "WT")
# Differential expression analysis-------------------------------------------
dds <- DESeq(dds)

res <-  results(dds, contrast=c("condition",unique(condition)))
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
table(res$change)
table(is.na(res$change))
# resOrdered <- as.data.frame(res[order(res$padj), ])
resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
resInP$geneSymbol <- rownames(resInP)
name = res@elementMetadata$description[2]
name = paste0(cl, '_',str_split(name, "condition ", simplify = T)[,2])
RMdeg_list[[i]] <- list( resInP)
names(RMdeg_list[[i]] ) <- paste0(gsub(' ','_',(name)))
countGTF_list <- list(colnames(countGTF)); names(countGTF_list) <- c
condition_list <- list(condition); names(condition_list) <- c
contrast_input_list[[i]] <- countGTF_list
condition_input_list[[i]] <- condition_list

###################################
#  Cell number-normalized (CNN)   #
#   MeRIP RIP analysis - Part1    #
###################################
## Sample information
merge_peak <- read.delim(paste0("metaPlotR/", c, "_newPeakFile.txt"),header=F)
countMatrix <- read.table(paste0("metaPlotR/m5C_counts_", c, ".txt"),header=T)
table(countMatrix$Geneid == merge_peak[-1,1])
table(duplicated(merge_peak$V1[-1]))
dim(countMatrix)
colnames(countMatrix)
rownames(countMatrix) <- countMatrix$Geneid
countMatrix <- countMatrix[, !colnames(countMatrix) %in% c('Geneid','Chr','Start','End','Strand','Length')]
# 1. 做好分组因子即可
groupList<-factor(condition)
colData <- data.frame(row.names=colnames(countMatrix), groupList=groupList)
# 2. Exploring Method 1: DESeq2-------------------------------------------
# 构建 DESeqDataSet 对象
dds <- DESeqDataSetFromMatrix(countData = countMatrix,
                              colData = colData,
                              design = ~ groupList)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 1
table(keep)
dds <- dds[keep,]
# 接着计算每个样本的 library size 用于后续的标准化。library size 等于落在该样本所有peaks上的reads的总数，即 counts 矩阵中每列的加和。如果bFullLibrarySize设为TRUE，则会使用library的总reads数（根据原BAM/BED文件统计获得）。然后使用estimateDispersions估计统计分布，需要将参数fitType设为’local’
dds<-estimateSizeFactors(dds)
dds<-estimateDispersions(dds,fitType="local")
dds <- nbinomWaldTest(dds)
res <- results(dds)
## 提取你想要的差异分析结果，我们这里是treated组对untreated组进行比较
res <-  results(dds, contrast=c("groupList", unique(condition)))
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
res_DESeq2 = as.data.frame(res)
table(res_DESeq2$change)#Not only
res_DESeq2 = res_DESeq2[res_DESeq2$change %in% c('UP', 'DOWN'), ]
# 2. Exploring Method 2: limma-------------------------------------------
library(edgeR)
d0 <- DGEList(counts = countMatrix)
dim(d0)
# Filtering low count peaks
cutoff <- 1           
drop <- which(apply(cpm(d0), 1, max) < cutoff)
print(length(drop))
if (length(drop) >0) {
  d <- d0[-drop,] 
} else{ d <- d0}
dim(d) #number of peaks left: 8501      
d <- calcNormFactors(d, method = c('TMM'))
colnames(countMatrix)
genotype <- groupList
table(genotype)
mm <- model.matrix(~0 + genotype)#
v <- voom(d, mm, plot=F, normalize="quantile")#voom from limma:Transform count data to log2-counts per million (logtpm)
fit <- lmFit(v, mm)
head(coef(fit))
contr <- makeContrasts(paste0(colnames(coef(fit)), collapse  = ' - '), levels = colnames(coef(fit))) #Change if other comparison
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)
length(which(top.table$adj.P.Val < 0.05 & abs(top.table$logFC) > log2(1.5))) #how many are significant?
res_limma <- as.data.frame(top.table)
res_limma$change <- with(res_limma, ifelse(res_limma$adj.P.Val < 0.05 & abs(res_limma$logFC) > log2(1.5),
                                           ifelse(res_limma$logFC > log2(1.5),'UP','DOWN'),'NOT'))
table(res_limma$change)#NOT only
res_limma = res_limma[res_limma$change %in% c('UP', 'DOWN'), ]

RMdem_list[[i]] <- list(res_DESeq2, res_limma)
name = res@elementMetadata$description[2]
name = paste0(cl, '_',str_split(name, "groupList ", simplify = T)[,2])
name= paste0(gsub(' ','_',(name)))
names(RMdem_list[[i]] ) <- paste0(name, '_', c('DESeq2', 'limma'))
RMdep_list[[i]] <- list( merge_peak )
names(RMdep_list[[i]] ) <- name
countMatrix_list <- list(colnames(countMatrix)); names(countMatrix_list) <- c
condition_list <- list(condition); names(condition_list) <- c
contrast_IP_list[[i]] <- countMatrix_list
condition_IP_list[[i]] <- condition_list

# 
# ############
# n = 91
# ############
# #GSE240879# without case RIP and Input
# ############
# print(RMPdeg_datasets[n,])# skip?
# i = RMPdeg_datasets[n,'GSE']
# cl = c('H460')#
# RMp = c('FTO')
# RM = 'm6A'
# # condition = c(rep('FTOKD', 2), rep('WT', 2))#
# RMP_allDat[grep('GSE240879', RMP_allDat$study_alias )
#            ,c("study_alias", "sample_alias", "sample_title")]
# # View(allDat[[i]])
# files <- list.files(file.path(data_directory,i))
# files
# getwd()
# setwd(data_directory)
# setwd(i)

############
n = 92
############
#GSE40132# without RIP, with RNA
############
print(RMPdeg_datasets[n,])#
i = RMPdeg_datasets[n,'GSE']
cl = c('Hela')#
RMp = c('ALKBH5')
RM = 'm6A'
condition = c(rep('ALKBH5KD', 1), rep('WT', 1))#
RMP_allDat[grep('GSE40132', RMP_allDat$study_alias )
           ,c("study_alias", "sample_alias", "sample_title")]
# View(allDat[[i]])
files <- list.files(file.path(data_directory,i))
files
getwd()
setwd(data_directory)
setwd(i)
# library(GenomicFeatures)
# hg19.ens <- makeTxDbFromUCSC(genome="hg19", tablename="ensGene")
# exonic <- exonsBy(hg19.ens, by="gene")
# red.exonic <- reduce(exonic)
# exon.lengths <- sum(width(red.exonic))
###################################
#         RNA-seq analysis        #
###################################
# Reads Per Kilobase of exon per Megabase of library size (RPKM) were calculated using cuffdiff
# 1) RNA-seq
count_files <- files[grepl('.gene.txt.gz', files)]
###################################
# Read and process each count file
count <- lapply(count_files[grepl('HeLa', count_files)], function(x) {
  # Read data using fread for faster reading
  data <- fread(file.path(data_directory, i, x))
  # Extract required columns
  df_new <- data[, c(2, 6)]
  df_new <- df_new[df_new$Symbol !='',]
  # Set column names same as file names
  colnames(df_new) <- c('geneSymbol',gsub('.gene.txt.gz', '', x))
  return(df_new)
})
# Merge dataframes in countGTF list by 'gene_id'
df_in <- Reduce(function(x, y) merge(x, y, by = 'geneSymbol', all = TRUE), count)
colSums(df_in[,-1])
table(colSums(df_in[,-1])=="1e+06")
# 
# rpkm_data <- read.table(paste0(data_directory,"/GSE120860/GSM3407085_SMMC7721-shYTHDF2.txt.gz"),comment="#", header = T)
# rpkm_data <- rpkm_data[, c('gene_name', 'length')]
# colnames(rpkm_data)[1] <- 'geneSymbol'
# table(rpkm_data$geneSymbol %in% df_in$geneSymbol)
# rpkm_data <- rpkm_data[rpkm_data$geneSymbol %in% df_in$geneSymbol, ]
# # df_new <- read.table(paste0(data_directory,"/GSE120860/GSM3405798_SMMC7721-OE+Nx.txt.gz"),comment="#", header = T, sep = '\t')# Total mapped reads in the sample (in millions)
# total_mapped_reads_million <- 25  # For example, 25 million reads (This is a reasonable typical value)
# # Add a column for Gene Length in Kilobases
# rpkm_data$Gene_Length_KB <- rpkm_data$length / 1000
# rpkm_data <- left_join(rpkm_data, df_in, by = 'geneSymbol')
# # Calculate Counts from RPKM
# rpkm_data$GSM986109_HeLa.Conc_counts <- (rpkm_data$GSM986109_HeLa.Conc * total_mapped_reads_million * rpkm_data$Gene_Length_KB) / 10^6
# rpkm_data$GSM986109_HeLa.ALKBH5_counts <- (rpkm_data$GSM986110_HeLa.ALKBH5 * total_mapped_reads_million * rpkm_data$Gene_Length_KB) / 10^6
# # 10^3标准化了基因长度的影响，10^6标准化了测序深度的影响。FPKM方法与RPKM类似，主要针对双末端RNA-seq实验的转录本定量。在双末端RNA-seq实验中，有左右两个对应的read来自相同的DNA片段。在进行双末端read进行比对时，来自同一DNA片段的高质量的一对或单个read可以定位到参考序列上。为避免混淆或多次计数，统计一对或单个read比对上的参考序列片段（Fragment），来计算FPKM，计算方法同RPKM。
# # RPKM与FPKM的区别：RPKM值适用于单末端RNA-seq实验数据，FPKM适用于双末端RNA-seq测序数据。
## Gene name
library(data.table)
table(duplicated(id_mapped$gene))
countGTF <- df_in#rpkm_data[, grepl('counts|gene', colnames(rpkm_data))]
countGTF <- aggregate(x = countGTF[,2:ncol(countGTF)],   #此时exprSet的第三列开始是表达矩阵内容
                      by = list(symbol = countGTF$geneSymbol),   #按照相同symbol分组，在组内计算
                      FUN = max) %>%   #原文中是计算最大值（max），也可以计算平均值（mean）或者中位数（median）
  column_to_rownames(var = 'symbol')
dim(countGTF) # 51242    12
# countGTF <- as.matrix(ceiling(countGTF))
# ##Gene count matrix
# y <- DGEList(countGTF)
# # Calculate fold change
# y <- calcNormFactors(y, method = c('TMM'))
# y$samples$lib.size <- colSums(y$counts)
logcpm <- log2(countGTF+1)#cpm(y, prior.count=2, log=TRUE)
logcpm[1:5,]
# cpm <- cpm(y, prior.count=2, log=FALSE)
# cpm[1:5, ]
fold_change <- countGTF[,1]/countGTF[,2]#cpm[,1]/cpm[,2]
res <- countGTF
res$logFC <- log2(fold_change)
res$FoldChange <- fold_change
res$change <- with(res, ifelse(abs(res$logFC) > log2(1.5),
                               ifelse(res$logFC > log2(1.5),'UP','DOWN'),'NOT'))
table(res$change)
resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
resInP$geneSymbol <- rownames(resInP)
table(resInP$change)
table(resInP$FoldChange %in% c(0, Inf))
table(resInP$logFC %in% c(-Inf, Inf))
# resInP <- resInP[! resInP$FoldChange %in% c(0, Inf), ]
name = "Hela_ALKBH5KD_vs_WT"
RMdeg_list[[i]] <- list( resInP)
names(RMdeg_list[[i]] ) <- name
condition_list <- list(condition); names(condition_list) <- 'ALKBH5KD'
contrast_input_list[[i]] <- NULL#countGTF_list
condition_input_list[[i]] <- condition_list

############
n = 93
############
#GSE49339# with RIP and Input
############
print(RMPdeg_datasets[n,])#
i = RMPdeg_datasets[n,'GSE']
cl = c('Hela')#use m6A-seq data from GSE144620 as WT
RMp = c('YTHDF2')
RM = 'm6A'
# condition = c(rep('METTL3KD', 2), rep('WT', 2))#
condition = c(rep('YTHDF2PD', 2), rep('WT', 2))#

RMP_allDat[grepl('GSE49339', RMP_allDat$study_alias ) & grepl('YTHDF2-RIP-', RMP_allDat$sample_title )
           ,c("study_alias", "sample_alias", "sample_title")]
# View(allDat[[i]])
ENA <- as.data.frame(RMP_allDat[grepl('GSE49339', RMP_allDat$study_alias ) & grepl('YTHDF2-RIP-', RMP_allDat$sample_title )
                                ,c("study_alias", "sample_alias", "sample_title")])
ENA[order(ENA$sample_alias),]
skip <- ENA[order(ENA$sample_alias),]
paste0(skip$sample_alias, collapse = "|")

skip <- RMP_allDat[grepl('GSE49339', RMP_allDat$study_alias ) & grepl('YTHDF2-RIP-', RMP_allDat$sample_title )
                   & grepl('-input', RMP_allDat$sample_title)
                   , "sample_alias"]
paste0(skip$sample_alias, collapse = " ")
skip <- RMP_allDat[grepl('GSE49339', RMP_allDat$study_alias ) & grepl('YTHDF2-RIP-', RMP_allDat$sample_title )
                   &  grepl('-IP', RMP_allDat$sample_title)
                   , "sample_alias"]
paste0(skip$sample_alias, collapse = " ")
files <- list.files(file.path(data_directory,i))
files
getwd()
setwd(data_directory)
setwd(i)
###################################
#       RNA/RIP-seq analysis      #
###################################
condition = c(rep('YTHDF2PD', 2), rep('WT', 4))#
c <- unique(condition[!condition %in% 'WT'])
countGTF <- read.table( paste0("metaPlotR/Input_counts_", c, ".txt"),header=T)
rownames(countGTF) <- countGTF$Geneid
countGTF <- countGTF[, !colnames(countGTF) %in% c('Geneid','Chr','Start','End','Strand','Length')]
dim(countGTF) #28265     4
## Sample name
## Sample information
sample=data.frame(
  sample=colnames(countGTF),
  condition=condition)
## SummarizedExperiment input
dds<-DESeqDataSetFromMatrix(countData=countGTF,colData=sample,design=~condition)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 10
table(keep)
dds <- dds[keep,]
## Note on factor levels
dds$condition <- factor(dds$condition, levels = unique(condition))
dds$condition <- relevel(dds$condition, ref = "WT")
# Differential expression analysis-------------------------------------------
dds <- DESeq(dds)

res <-  results(dds, contrast=c("condition",unique(condition)))
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
table(res$change)
table(is.na(res$change))
# resOrdered <- as.data.frame(res[order(res$padj), ])
resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
resInP$geneSymbol <- rownames(resInP)
name = res@elementMetadata$description[2]
name = paste0(cl, '_',str_split(name, "condition ", simplify = T)[,2])
RMdeg_list[[i]] <- list( resInP)
names(RMdeg_list[[i]] ) <- paste0(gsub(' ','_',(name)))
countGTF_list <- list(colnames(countGTF)); names(countGTF_list) <- c
condition_list <- list(condition); names(condition_list) <- c
contrast_input_list[[i]] <- countGTF_list
condition_input_list[[i]] <- condition_list

###################################
#  Cell number-normalized (CNN)   #
#   MeRIP RIP analysis - Part1    #
###################################
condition = c(rep('YTHDF2PD', 2), rep('WT', 2))#
## Sample information
merge_peak <- read.delim(paste0("metaPlotR/", c, "_newPeakFile.txt"),header=F)
countMatrix <- read.table(paste0("metaPlotR/m6A_counts_", c, ".txt"),header=T)
table(countMatrix$Geneid == merge_peak[-1,1])
table(duplicated(merge_peak$V1[-1]))
dim(countMatrix)
colnames(countMatrix)
rownames(countMatrix) <- countMatrix$Geneid
countMatrix <- countMatrix[, !colnames(countMatrix) %in% c('Geneid','Chr','Start','End','Strand','Length')]
# 1. 做好分组因子即可
groupList<-factor(condition)
colData <- data.frame(row.names=colnames(countMatrix), groupList=groupList)
# 2. Exploring Method 1: DESeq2-------------------------------------------
# 构建 DESeqDataSet 对象
dds <- DESeqDataSetFromMatrix(countData = countMatrix,
                              colData = colData,
                              design = ~ groupList)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 1
table(keep)
dds <- dds[keep,]
# 接着计算每个样本的 library size 用于后续的标准化。library size 等于落在该样本所有peaks上的reads的总数，即 counts 矩阵中每列的加和。如果bFullLibrarySize设为TRUE，则会使用library的总reads数（根据原BAM/BED文件统计获得）。然后使用estimateDispersions估计统计分布，需要将参数fitType设为’local’
dds<-estimateSizeFactors(dds)
dds<-estimateDispersions(dds,fitType="local")
dds <- nbinomWaldTest(dds)
res <- results(dds)
## 提取你想要的差异分析结果，我们这里是treated组对untreated组进行比较
res <-  results(dds, contrast=c("groupList", unique(condition)))
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
res_DESeq2 = as.data.frame(res)
table(res_DESeq2$change)#Not only
res_DESeq2 = res_DESeq2[res_DESeq2$change %in% c('UP', 'DOWN'), ]
# 2. Exploring Method 2: limma-------------------------------------------
library(edgeR)
d0 <- DGEList(counts = countMatrix)
dim(d0)
# Filtering low count peaks
cutoff <- 1           
drop <- which(apply(cpm(d0), 1, max) < cutoff)
print(length(drop))
if (length(drop) >0) {
  d <- d0[-drop,] 
} else{ d <- d0}
dim(d) #number of peaks left: 8501      
d <- calcNormFactors(d, method = c('TMM'))
colnames(countMatrix)
genotype <- groupList
table(genotype)
mm <- model.matrix(~0 + genotype)#
v <- voom(d, mm, plot=F, normalize="quantile")#voom from limma:Transform count data to log2-counts per million (logtpm)
fit <- lmFit(v, mm)
head(coef(fit))
contr <- makeContrasts(paste0(colnames(coef(fit)), collapse  = ' - '), levels = colnames(coef(fit))) #Change if other comparison
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)
length(which(top.table$adj.P.Val < 0.05 & abs(top.table$logFC) > log2(1.5))) #how many are significant?
res_limma <- as.data.frame(top.table)
res_limma$change <- with(res_limma, ifelse(res_limma$adj.P.Val < 0.05 & abs(res_limma$logFC) > log2(1.5),
                                           ifelse(res_limma$logFC > log2(1.5),'UP','DOWN'),'NOT'))
table(res_limma$change)#NOT only
res_limma = res_limma[res_limma$change %in% c('UP', 'DOWN'), ]

RMdem_list[[i]] <- list(res_DESeq2, res_limma)
name = res@elementMetadata$description[2]
name = paste0(cl, '_',str_split(name, "groupList ", simplify = T)[,2])
name= paste0(gsub(' ','_',(name)))
names(RMdem_list[[i]] ) <- paste0(name, '_', c('DESeq2', 'limma'))
RMdep_list[[i]] <- list( merge_peak )
names(RMdep_list[[i]] ) <- name
countMatrix_list <- list(colnames(countMatrix)); names(countMatrix_list) <- c
condition_list <- list(condition); names(condition_list) <- c
contrast_IP_list[[i]] <- countMatrix_list
condition_IP_list[[i]] <- condition_list

############
n = 94
############
#GSE63591# without RIP, with RNA
############
print(RMPdeg_datasets[n,])#
i = RMPdeg_datasets[n,'GSE']
cl = c('Hela')#
RMp = c('METTL3')
RM = 'm6A'
condition = c(rep('METTL3KD', 1), rep('WT', 1))#
RMP_allDat[grep('GSE63591', RMP_allDat$study_alias )
           ,c("study_alias", "sample_alias", "sample_title")]
# View(allDat[[i]])
files <- list.files(file.path(data_directory,i))
files
getwd()
setwd(data_directory)
setwd(i)
df_in <- read_excel('GSE63591_D-METTL3-Ribosome_profiling.xlsx')
colN <- as.data.frame(df_in[1,c(1,5,6,7)])
colN$rpkm<- 'geneSymbol'; colN$...7 <- 'logFC'
countGTF <- as.data.frame(df_in[-1,])[, c(1,5,6,7)]
colnames(countGTF) <- colN[1,]
fold_change <- as.numeric(countGTF[,3])/as.numeric(countGTF[,2])
res <- countGTF
res$logFC <- as.numeric(res$logFC)
res$FoldChange <- fold_change
res$change <- with(res, ifelse(abs(res$logFC) > log2(1.5),
                               ifelse(res$logFC > log2(1.5),'UP','DOWN'),'NOT'))
table(res$change)
resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
table(resInP$change)
table(resInP$FoldChange %in% c(0, Inf))
table(resInP$logFC %in% c(-Inf, Inf))
# resInP <- resInP[! resInP$FoldChange %in% c(0, Inf), ]
name = "Hela_METTL3KD_vs_WT"
RMdeg_list[[i]] <- list( resInP)
names(RMdeg_list[[i]] ) <- name
condition_list <- list(condition); names(condition_list) <- 'METTL3KD'
contrast_input_list[[i]] <- NULL#countGTF_list
condition_input_list[[i]] <- condition_list

############
n = 95
############
#GSE86214# with RIP and Input
############
print(RMPdeg_datasets[n,])# hg38
i = RMPdeg_datasets[n,'GSE']
cl = c('Hela')#(use m6A-seq data from GSE46705 as WT)
RMp = c('YTHDF3')
RM = 'm6A'
condition_list = list(
  Hela_YTHDF3PD = c(rep('YTHDF3PD', 2), rep('WT', 6)),#
  Hela.YTHDF3KD_YTHDF1PD = c(rep('YTHDF3KD_YTHDF1PD', 2), rep('WT', 1)),#WT:YTHDF1PD
  Hela.YTHDF3KD_YTHDF2PD = c(rep('YTHDF3KD_YTHDF2PD', 2), rep('WT', 1))#WT:YTHDF2PD
)
RMP_allDat[grep('GSE86214', RMP_allDat$study_alias )
           ,c("study_alias", "sample_alias", "sample_title")]
# View(allDat[[i]])
# setwd(data_directory)
# gse <- getGEO('GSE86214',destdir=".", AnnotGPL=F, getGPL=F)
# setwd('GSE86214')
# pdata <- pData(gse[[2]])
RMP_allDat[grepl('GSE86214', RMP_allDat$study_alias ) & grepl('HeLa', RMP_allDat$sample_title)
           ,c("study_alias", "sample_alias", "sample_title")]
skip = RMP_allDat[grepl('GSE86214', RMP_allDat$study_alias ) & grepl('RIP-input|RIP-IP', RMP_allDat$sample_title)
                  ,c("study_alias", "sample_alias", "sample_title")]
paste0(skip$sample_alias, collapse = "|")

skip <- RMP_allDat[grepl('GSE86214', RMP_allDat$study_alias ) &  grepl('RIP-input|RIP-IP', RMP_allDat$sample_title)
                   &  grepl('RIP-input', RMP_allDat$sample_title)
                   , "sample_alias"]
paste0(skip$sample_alias, collapse = " ")
skip <- RMP_allDat[grepl('GSE86214', RMP_allDat$study_alias ) & grepl('RIP-input|RIP-IP', RMP_allDat$sample_title)
                   &  grepl('RIP-IP', RMP_allDat$sample_title)
                   , "sample_alias"]
paste0(skip$sample_alias, collapse = " ")
ENA <- as.data.frame(RMP_allDat[grepl('GSE86214', RMP_allDat$study_alias ) & grepl('RIP-input|RIP-IP', RMP_allDat$sample_title)
                                ,c("study_alias", "sample_alias", "sample_title")])
ENA[order(ENA$sample_alias),]

files <- list.files(file.path(data_directory,i))
files
getwd()
setwd(data_directory)
setwd(i)

###################################
#  Cell number-normalized (CNN)   #
#      MeRIP INPUT analysis       #
###################################
resInP_list <- vector("list", length(condition_list))
names(resInP_list) <- names(condition_list)
countGTF_list <- vector("list", length(condition_list))
names(countGTF_list) <- names(condition_list)
for(c in names(condition_list)){
  d = gsub('Hela.|Hela_', '', c)
  countGTF <- read.table( paste0("metaPlotR/Input_counts_", d, ".txt"),header=T)
  dim(countGTF)
  colnames(countGTF)
  rownames(countGTF) <- countGTF$Geneid
  countGTF <- countGTF[, !colnames(countGTF) %in% c('Geneid','Chr','Start','End','Strand','Length')]
  dim(countGTF) # 28265     4
  condition <- condition_list[[c]]
  ## Sample information
  sample=data.frame(
    sample=colnames(countGTF),
    condition=condition)
  ## SummarizedExperiment input
  dds<-DESeqDataSetFromMatrix(countData=countGTF,colData=sample,design=~condition)
  ## Pre-filtering
  keep <- rowSums(counts(dds)) >= 10
  table(keep)
  dds <- dds[keep,]
  ## Note on factor levels
  dds$condition <- factor(dds$condition, levels = unique(condition))
  dds$condition <- relevel(dds$condition, ref = "WT")
  # Differential expression analysis-------------------------------------------
  dds <- DESeq(dds)
  res <-  results(dds, contrast=c("condition",unique(condition)))
  res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
  length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
  res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                                 ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
  table(res$change)
  table(is.na(res$change))
  # resOrdered <- as.data.frame(res[order(res$padj), ])
  resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
  resInP$geneSymbol <- rownames(resInP)
  table(resInP$change)
  resInP_list[[c]] <- resInP
  countGTF_list[[c]] <- countGTF
}
# name = paste0(cl, '_',names(resInP_list), ' vs WT')
name = paste0(names(resInP_list), ' vs WT')
RMdeg_list[[i]] <- resInP_list 
names(RMdeg_list[[i]] ) <- paste0(gsub(' ','_',(name)))
contrast_input_list[[i]] <- lapply(countGTF_list, extract_names)
condition_input_list[[i]] <- condition_list

###################################
#  Cell number-normalized (CNN)   #
#   MeRIP RIP analysis - Part1    #
###################################
condition_list = list(
  Hela_YTHDF3PD = c(rep('YTHDF3PD', 2), rep('WT', 4)),#
  Hela.YTHDF3KD_YTHDF1PD = c(rep('YTHDF3KD_YTHDF1PD', 2), rep('WT', 1)),#WT:YTHDF1PD
  Hela.YTHDF3KD_YTHDF2PD = c(rep('YTHDF3KD_YTHDF2PD', 2), rep('WT', 1))#WT:YTHDF2PD
)
resMeR_list <- vector("list", length(condition_list))
names(resMeR_list) <- names(condition_list)
countMatrix_list <- vector("list", length(condition_list))
names(countMatrix_list) <- names(condition_list)
resPeak_list <- vector("list",length(condition_list))
names(resPeak_list) <- names(condition_list)
for(c in names(condition_list)){
  d = gsub('Hela.|Hela_', '', c)
  merge_peak <- read.delim(paste0("metaPlotR/", d, "_newPeakFile.txt"),header=F)
  countMatrix <- read.table(paste0("metaPlotR/m6A_counts_", d, ".txt"),header=T)
  table(countMatrix$Geneid == merge_peak[-1,1])
  table(duplicated(merge_peak[-1,1]))
  dim(countMatrix)
  colnames(countMatrix)
  rownames(countMatrix) <- countMatrix$Geneid
  countMatrix <- countMatrix[, !colnames(countMatrix) %in% c('Geneid','Chr','Start','End','Strand','Length')]
  # 1. 做好分组因子即可
  condition <- condition_list[[c]]
  groupList<-factor(condition)
  colData <- data.frame(row.names=colnames(countMatrix), groupList=groupList)
  # 2. Exploring Method 1: DESeq2-------------------------------------------
  # 构建 DESeqDataSet 对象
  dds <- DESeqDataSetFromMatrix(countData = countMatrix,
                                colData = colData,
                                design = ~ groupList)
  ## Pre-filtering
  keep <- rowSums(counts(dds)) >= 1
  table(keep)
  dds <- dds[keep,]
  # 接着计算每个样本的 library size 用于后续的标准化。library size 等于落在该样本所有peaks上的reads的总数，即 counts 矩阵中每列的加和。如果bFullLibrarySize设为TRUE，则会使用library的总reads数（根据原BAM/BED文件统计获得）。然后使用estimateDispersions估计统计分布，需要将参数fitType设为’local’
  dds<-estimateSizeFactors(dds)
  dds<-estimateDispersions(dds,fitType="local")
  dds <- nbinomWaldTest(dds)
  res <- results(dds)
  ## 提取你想要的差异分析结果，我们这里是treated组对untreated组进行比较
  res <-  results(dds, contrast=c("groupList", unique(condition)))
  res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
  length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
  res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                                 ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
  res_DESeq2 = as.data.frame(res)
  table(res_DESeq2$change)
  res_DESeq2 = res_DESeq2[res_DESeq2$change %in% c('UP', 'DOWN'), ]
  # 2. Exploring Method 2: limma-------------------------------------------
  library(edgeR)
  d0 <- DGEList(counts = countMatrix)
  dim(d0)
  # Filtering low count peaks
  cutoff <- 1           
  drop <- which(apply(cpm(d0), 1, max) < cutoff)
  print(length(drop))
  if (length(drop) >0) {
    d <- d0[-drop,] 
  } else{ d <- d0}
  dim(d) #number of peaks left: 8501      
  d <- calcNormFactors(d, method = c('TMM'))
  colnames(countMatrix)
  genotype <- groupList
  table(genotype)
  mm <- model.matrix(~0 + genotype)#
  v <- voom(d, mm, plot=F, normalize="quantile")#voom from limma:Transform count data to log2-counts per million (logtpm)
  fit <- lmFit(v, mm)
  head(coef(fit))
  contr <- makeContrasts(paste0(colnames(coef(fit)), collapse  = ' - '), levels = colnames(coef(fit))) #Change if other comparison
  tmp <- contrasts.fit(fit, contr)
  tmp <- eBayes(tmp)
  top.table <- topTable(tmp, sort.by = "P", n = Inf)
  head(top.table, 20)
  length(which(top.table$adj.P.Val < 0.05 & abs(top.table$logFC) > log2(1.5))) #how many are significant?
  res_limma <- as.data.frame(top.table)
  res_limma$change <- with(res_limma, ifelse(res_limma$adj.P.Val < 0.05 & abs(res_limma$logFC) > log2(1.5),
                                             ifelse(res_limma$logFC > log2(1.5),'UP','DOWN'),'NOT'))
  table(res_limma$change)
  res_limma = res_limma[res_limma$change %in% c('UP', 'DOWN'), ]
  resMeR_list[[c]] <- list(res_DESeq2, res_limma)
  countMatrix_list[[c]] <- countMatrix
  resPeak_list[[c]] <- merge_peak
  name = paste0(c, ' vs WT')
  name = paste0(gsub(' ','_', name ))
  names(resMeR_list[[c]] ) <- paste0(name, '_', c('DESeq2', 'limma'))
}
names(resMeR_list) <- NULL
RMdem_list[[i]] <- as.list(unlist(resMeR_list, recursive = FALSE))
# names(RMdem_list[[i]] ) <- paste0(cl, '_', names(RMdem_list[[i]] ))
RMdep_list[[i]] <- resPeak_list 
name <- unique(gsub('(_limma|_DESeq2)', '', names(RMdem_list[[i]] )))
names(RMdep_list[[i]] ) <- name
contrast_IP_list[[i]] <- lapply(countMatrix_list, extract_names)
condition_IP_list[[i]] <- condition_list

############
n = 96
############
#GSE79577# without case RIP, with RNA
############
print(RMPdeg_datasets[n,])#
i = RMPdeg_datasets[n,'GSE']
cl = c('HEK293')#
RMp = c('FTO')
RM = 'm6A'
condition = c(rep('FTOKO', 3), rep('WT', 3))#
RMP_allDat[grep('GSE79577', RMP_allDat$study_alias )
           ,c("study_alias", "sample_alias", "sample_title")]
# View(allDat[[i]])
files <- list.files(file.path(data_directory,i))
files
getwd()
setwd(data_directory)
setwd(i)
df_in <- read.csv('GSE79577_DESeq2_FTO_KO.csv.gz', sep = ',')
countGTF <- df_in[, 2:7]
id_mapped <- read.table("/data2/rluo4/hg38/gencode.v36.annotation.gene.probeMap",header = T)
countGTF$id <- str_split(df_in$X, '[.]', simplify = T)[,1]
id_mapped$id <- str_split(id_mapped$id, '[.]', simplify = T)[,1]
id_mapped <- id_mapped[which(id_mapped$id %in% countGTF$id ),]
dim(id_mapped) #57930     6
table(duplicated(id_mapped$gene))
countGTF <- left_join(id_mapped[,1:2],countGTF,by="id")
tmp=by(countGTF,countGTF$gene,function(x) rownames(x)[which.max(x$baseMean)])
probes = as.character(tmp)
length(probes) # 7785
res <- countGTF[rownames(countGTF) %in% probes,-1]
colnames(res)[1] <- 'geneSymbol'
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
table(res$change)
resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
table(duplicated(resInP$geneSymbol))
name = "HEK293_FTOKO_vs_WT"
RMdeg_list[[i]] <- list( resInP)
names(RMdeg_list[[i]] ) <- name
condition_list <- list(condition); names(condition_list) <- 'FTOKO'
contrast_input_list[[i]] <- NULL#countGTF_list
condition_input_list[[i]] <- condition_list

############
n = 97
############
#GSE142386# with RIP and Input
############
print(RMPdeg_datasets[n,])#
i = RMPdeg_datasets[n,'GSE']
cl = c('HUVECs')#Human Umbilical Vein Endothelial Cells
RMp = c('METTL3', "WTAP")
RM = 'm6A'
condition_list = list(METTL3KD = c(rep('METTL3KD', 2), rep('WT', 2)),
                      WTAPKD = c(rep('WTAPKD', 2), rep('WT', 2)))
RMP_allDat[grep('GSE142386', RMP_allDat$study_alias )
           ,c("study_alias", "sample_alias", "sample_title")]
# View(allDat[[i]])
skip <- RMP_allDat[grepl('GSE142386', RMP_allDat$study_alias )
                   &  grepl('Input', RMP_allDat$sample_title)
                   , "sample_alias"]
paste0(skip$sample_alias, collapse = " ")
skip <- RMP_allDat[grepl('GSE142386', RMP_allDat$study_alias )
                   &  grepl('IP', RMP_allDat$sample_title)
                   , "sample_alias"]
paste0(skip$sample_alias, collapse = " ")
ENA <- as.data.frame(RMP_allDat[grepl('GSE142386', RMP_allDat$study_alias ) 
                                ,c("study_alias", "sample_alias", "sample_title")])
ENA[order(ENA$sample_alias),]

files <- list.files(file.path(data_directory,i))
files
getwd()
setwd(data_directory)
setwd(i)

###################################
#  Cell number-normalized (CNN)   #
#      MeRIP INPUT analysis       #
###################################
resInP_list <- vector("list", length(condition_list))
names(resInP_list) <- names(condition_list)
countGTF_list <- vector("list", length(condition_list))
names(countGTF_list) <- names(condition_list)
for(c in names(condition_list)){
  countGTF <- read.table( paste0("metaPlotR/Input_counts_", c, ".txt"),header=T)
  dim(countGTF)
  colnames(countGTF)
  rownames(countGTF) <- countGTF$Geneid
  countGTF <- countGTF[, !colnames(countGTF) %in% c('Geneid','Chr','Start','End','Strand','Length')]
  dim(countGTF) # 28265     4
  condition <- condition_list[[c]]
  ## Sample information
  sample=data.frame(
    sample=colnames(countGTF),
    condition=condition)
  ## SummarizedExperiment input
  dds<-DESeqDataSetFromMatrix(countData=countGTF,colData=sample,design=~condition)
  ## Pre-filtering
  keep <- rowSums(counts(dds)) >= 10
  table(keep)
  dds <- dds[keep,]
  ## Note on factor levels
  dds$condition <- factor(dds$condition, levels = unique(condition))
  dds$condition <- relevel(dds$condition, ref = "WT")
  # Differential expression analysis-------------------------------------------
  dds <- DESeq(dds)
  res <-  results(dds, contrast=c("condition",unique(condition)))
  res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
  length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
  res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                                 ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
  table(res$change)
  table(is.na(res$change))
  # resOrdered <- as.data.frame(res[order(res$padj), ])
  resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
  resInP$geneSymbol <- rownames(resInP)
  table(resInP$change)
  resInP_list[[c]] <- resInP
  countGTF_list[[c]] <- countGTF
}
name = paste0(cl, '_',names(resInP_list), ' vs WT')
RMdeg_list[[i]] <- resInP_list 
names(RMdeg_list[[i]] ) <- paste0(gsub(' ','_',(name)))
contrast_input_list[[i]] <- lapply(countGTF_list, extract_names)
condition_input_list[[i]] <- condition_list

###################################
#  Cell number-normalized (CNN)   #
#   MeRIP RIP analysis - Part1    #
###################################
resMeR_list <- vector("list", length(condition_list))
names(resMeR_list) <- names(condition_list)
countMatrix_list <- vector("list", length(condition_list))
names(countMatrix_list) <- names(condition_list)
resPeak_list <- vector("list",length(condition_list))
names(resPeak_list) <- names(condition_list)
for(c in names(condition_list)){
  merge_peak <- read.delim(paste0("metaPlotR/", c, "_newPeakFile.txt"),header=F)
  countMatrix <- read.table(paste0("metaPlotR/m6A_counts_", c, ".txt"),header=T)
  table(countMatrix$Geneid == merge_peak[-1,1])
  table(duplicated(merge_peak[-1,1]))
  dim(countMatrix)
  colnames(countMatrix)
  rownames(countMatrix) <- countMatrix$Geneid
  countMatrix <- countMatrix[, !colnames(countMatrix) %in% c('Geneid','Chr','Start','End','Strand','Length')]
  # 1. 做好分组因子即可
  condition <- condition_list[[c]]
  groupList<-factor(condition)
  colData <- data.frame(row.names=colnames(countMatrix), groupList=groupList)
  # 2. Exploring Method 1: DESeq2-------------------------------------------
  # 构建 DESeqDataSet 对象
  dds <- DESeqDataSetFromMatrix(countData = countMatrix,
                                colData = colData,
                                design = ~ groupList)
  ## Pre-filtering
  keep <- rowSums(counts(dds)) >= 1
  table(keep)
  dds <- dds[keep,]
  # 接着计算每个样本的 library size 用于后续的标准化。library size 等于落在该样本所有peaks上的reads的总数，即 counts 矩阵中每列的加和。如果bFullLibrarySize设为TRUE，则会使用library的总reads数（根据原BAM/BED文件统计获得）。然后使用estimateDispersions估计统计分布，需要将参数fitType设为’local’
  dds<-estimateSizeFactors(dds)
  dds<-estimateDispersions(dds,fitType="local")
  dds <- nbinomWaldTest(dds)
  res <- results(dds)
  ## 提取你想要的差异分析结果，我们这里是treated组对untreated组进行比较
  res <-  results(dds, contrast=c("groupList", unique(condition)))
  res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
  length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
  res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                                 ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
  res_DESeq2 = as.data.frame(res)
  table(res_DESeq2$change)
  res_DESeq2 = res_DESeq2[res_DESeq2$change %in% c('UP', 'DOWN'), ]
  # 2. Exploring Method 2: limma-------------------------------------------
  library(edgeR)
  d0 <- DGEList(counts = countMatrix)
  dim(d0)
  # Filtering low count peaks
  cutoff <- 1           
  drop <- which(apply(cpm(d0), 1, max) < cutoff)
  print(length(drop))
  if (length(drop) >0) {
    d <- d0[-drop,] 
  } else{ d <- d0}
  dim(d) #number of peaks left: 8501      
  d <- calcNormFactors(d, method = c('TMM'))
  colnames(countMatrix)
  genotype <- groupList
  table(genotype)
  mm <- model.matrix(~0 + genotype)#
  v <- voom(d, mm, plot=F, normalize="quantile")#voom from limma:Transform count data to log2-counts per million (logtpm)
  fit <- lmFit(v, mm)
  head(coef(fit))
  contr <- makeContrasts(paste0(colnames(coef(fit)), collapse  = ' - '), levels = colnames(coef(fit))) #Change if other comparison
  tmp <- contrasts.fit(fit, contr)
  tmp <- eBayes(tmp)
  top.table <- topTable(tmp, sort.by = "P", n = Inf)
  head(top.table, 20)
  length(which(top.table$adj.P.Val < 0.05 & abs(top.table$logFC) > log2(1.5))) #how many are significant?
  res_limma <- as.data.frame(top.table)
  res_limma$change <- with(res_limma, ifelse(res_limma$adj.P.Val < 0.05 & abs(res_limma$logFC) > log2(1.5),
                                             ifelse(res_limma$logFC > log2(1.5),'UP','DOWN'),'NOT'))
  table(res_limma$change)
  res_limma = res_limma[res_limma$change %in% c('UP', 'DOWN'), ]
  resMeR_list[[c]] <- list(res_DESeq2, res_limma)
  countMatrix_list[[c]] <- countMatrix
  resPeak_list[[c]] <- merge_peak
  name = paste0(c, ' vs WT')
  name = paste0(gsub(' ','_', name ))
  names(resMeR_list[[c]] ) <- paste0(name, '_', c('DESeq2', 'limma'))
}
names(resMeR_list) <- NULL
RMdem_list[[i]] <- as.list(unlist(resMeR_list, recursive = FALSE))
names(RMdem_list[[i]] ) <- paste0(cl, '_', names(RMdem_list[[i]] ))
RMdep_list[[i]] <- resPeak_list 
name <- unique(gsub('(_limma|_DESeq2)', '', names(RMdem_list[[i]] )))
names(RMdep_list[[i]] ) <- name
contrast_IP_list[[i]] <- lapply(countMatrix_list, extract_names)
condition_IP_list[[i]] <- condition_list

############
n = 98
############
#GSE133517# without case RIP, with RNA
############
print(RMPdeg_datasets[n,])#
i = RMPdeg_datasets[n,'GSE']
cl = c('Hela')#
RMp = c('FTO')
RM = 'm6A'
condition = c(rep('FTOKD', 2), rep('WT', 2))#
RMP_allDat[grep('GSE133517', RMP_allDat$study_alias )
           ,c("study_alias", "sample_alias", "sample_title")]
# View(allDat[[i]])
files <- list.files(file.path(data_directory,i))
files
getwd()
setwd(data_directory)
setwd(i)
###################################
#      RNA/Ribo-seq analysis      #
###################################
# 1) RNA-seq
count_files <- files[grepl('.genes.results.txt.gz', files)]
###################################
# Read and process each count file
count <- lapply(files, function(x) {
  data <- read.delim(file.path(data_directory, i, x))
  df_new <- data[, grepl('gene_id|count', colnames(data))]
  # Set column names same as file names
  colnames(df_new) <- c('geneSymbol',gsub('.genes.results.txt.gz', '', x))
  df_new <- df_new[df_new[,1] !='-', ]
  df_new[,1] <- str_split(df_new[,1], ',', simplify = T)[,1]
  df_new <- df_new[df_new[,1] !='-', ]
  df_new <- aggregate(x = df_new,   #此时exprSet的第三列开始是表达矩阵内容
                      by = list(geneSymbol = df_new$geneSymbol),   #按照相同symbol分组，在组内计算
                      FUN = max) %>%   #原文中是计算最大值（max），也可以计算平均值（mean）或者中位数（median）
    column_to_rownames(var = 'geneSymbol')
  
  return(df_new)
})
# Merge dataframes in countGTF list by 'gene_id'
df_in <- Reduce(function(x, y) merge(x, y, by = 'geneSymbol', all = TRUE), count)
colSums(df_in[,-1])
table(colSums(df_in[,-1])=="1e+06")
## Gene name
library(data.table)
countGTF= ceiling(df_in[,2:ncol(df_in)])
## Gene name
table(duplicated(df_in$geneSymbol))
dim(countGTF) #25343    12
rownames(countGTF) <- df_in$geneSymbol
## Sample name
colnames(countGTF)
## Sample information
sample=data.frame(
  sample=colnames(countGTF),
  condition = c(rep('WT', 2), rep('FTOKD', 2)))
## SummarizedExperiment input
dds<-DESeqDataSetFromMatrix(countData=countGTF,colData=sample,design=~condition)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 10
table(keep)
dds <- dds[keep,]
## Note on factor levels
dds$condition <- factor(dds$condition, levels = unique(condition))
dds$condition <- relevel(dds$condition, ref = "WT")
# Differential expression analysis-------------------------------------------
dds <- DESeq(dds)
res <-  results(dds, contrast=c("condition",unique(condition)))
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
table(res$change)
table(is.na(res$change))
# resOrdered <- as.data.frame(res[order(res$padj), ])
resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
resInP$geneSymbol <- rownames(resInP)
name = res@elementMetadata$description[2]
name = paste0(cl, '_',str_split(name, "condition ", simplify = T)[,2])
RMdeg_list[[i]] <- list( resInP)
names(RMdeg_list[[i]] ) <- paste0(gsub(' ','_',(name)))
c <- unique(condition[!condition %in% 'WT'])
condition_list <- list(condition); names(condition_list) <- c
contrast_input_list[[i]] <- NULL
condition_input_list[[i]] <- condition_list


############
n = 99
############
#GSE144620# with RIP and Input
############
print(RMPdeg_datasets[n,])#
i = RMPdeg_datasets[n,'GSE']
cl = c('Hela')#
RMp = c('ALKBH5')
RM = 'm6A'
condition_list = list( ALKBH5KO = c(rep('ALKBH5KO', 2), rep('WT', 2)),
                       ALKBH5MU = c(rep('ALKBH5MU', 2), rep('WT', 2)))#
RMP_allDat[grep('GSE144620', RMP_allDat$study_alias )
           ,c("study_alias", "sample_alias", "sample_title")]
# View(allDat[[i]])
skip <- RMP_allDat[grepl('GSE144620', RMP_allDat$study_alias )
                   &  grepl('Input', RMP_allDat$sample_title)
                   , "sample_alias"]
paste0(skip$sample_alias, collapse = " ")
skip <- RMP_allDat[grepl('GSE144620', RMP_allDat$study_alias )
                   &  grepl('IP', RMP_allDat$sample_title)
                   , "sample_alias"]
paste0(skip$sample_alias, collapse = " ")
files <- list.files(file.path(data_directory,i))
files
getwd()
setwd(data_directory)
setwd(i)
ENA <- as.data.frame(RMP_allDat[grepl('GSE144620', RMP_allDat$study_alias )
                                ,c("study_alias", "sample_alias", "sample_title")])
ENA[order(ENA$sample_alias),]

files <- list.files(file.path(data_directory,i))
files
getwd()
setwd(data_directory)
setwd(i)

###################################
#  Cell number-normalized (CNN)   #
#      MeRIP INPUT analysis       #
###################################
resInP_list <- vector("list", length(condition_list))
names(resInP_list) <- names(condition_list)
countGTF_list <- vector("list", length(condition_list))
names(countGTF_list) <- names(condition_list)
for(c in names(condition_list)){
  countGTF <- read.table( paste0("metaPlotR/Input_counts_", c, ".txt"),header=T)
  dim(countGTF)
  colnames(countGTF)
  rownames(countGTF) <- countGTF$Geneid
  countGTF <- countGTF[, !colnames(countGTF) %in% c('Geneid','Chr','Start','End','Strand','Length')]
  dim(countGTF) # 28265     4
  condition <- condition_list[[c]]
  ## Sample information
  sample=data.frame(
    sample=colnames(countGTF),
    condition=condition)
  ## SummarizedExperiment input
  dds<-DESeqDataSetFromMatrix(countData=countGTF,colData=sample,design=~condition)
  ## Pre-filtering
  keep <- rowSums(counts(dds)) >= 10
  table(keep)
  dds <- dds[keep,]
  ## Note on factor levels
  dds$condition <- factor(dds$condition, levels = unique(condition))
  dds$condition <- relevel(dds$condition, ref = "WT")
  # Differential expression analysis-------------------------------------------
  dds <- DESeq(dds)
  res <-  results(dds, contrast=c("condition",unique(condition)))
  res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
  length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
  res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                                 ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
  table(res$change)
  table(is.na(res$change))
  # resOrdered <- as.data.frame(res[order(res$padj), ])
  resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
  resInP$geneSymbol <- rownames(resInP)
  table(resInP$change)
  resInP_list[[c]] <- resInP
  countGTF_list[[c]] <- countGTF
}
name = paste0(cl, '_',names(resInP_list), ' vs WT')
RMdeg_list[[i]] <- resInP_list 
names(RMdeg_list[[i]] ) <- paste0(gsub(' ','_',(name)))
contrast_input_list[[i]] <- lapply(countGTF_list, extract_names)
condition_input_list[[i]] <- condition_list

###################################
#  Cell number-normalized (CNN)   #
#   MeRIP RIP analysis - Part1    #
###################################
resMeR_list <- vector("list", length(condition_list))
names(resMeR_list) <- names(condition_list)
countMatrix_list <- vector("list", length(condition_list))
names(countMatrix_list) <- names(condition_list)
resPeak_list <- vector("list",length(condition_list))
names(resPeak_list) <- names(condition_list)
for(c in names(condition_list)){
  merge_peak <- read.delim(paste0("metaPlotR/", c, "_newPeakFile.txt"),header=F)
  countMatrix <- read.table(paste0("metaPlotR/m6A_counts_", c, ".txt"),header=T)
  table(countMatrix$Geneid == merge_peak[-1,1])
  table(duplicated(merge_peak[-1,1]))
  dim(countMatrix)
  colnames(countMatrix)
  rownames(countMatrix) <- countMatrix$Geneid
  countMatrix <- countMatrix[, !colnames(countMatrix) %in% c('Geneid','Chr','Start','End','Strand','Length')]
  # 1. 做好分组因子即可
  condition <- condition_list[[c]]
  groupList<-factor(condition)
  colData <- data.frame(row.names=colnames(countMatrix), groupList=groupList)
  # 2. Exploring Method 1: DESeq2-------------------------------------------
  # 构建 DESeqDataSet 对象
  dds <- DESeqDataSetFromMatrix(countData = countMatrix,
                                colData = colData,
                                design = ~ groupList)
  ## Pre-filtering
  keep <- rowSums(counts(dds)) >= 1
  table(keep)
  dds <- dds[keep,]
  # 接着计算每个样本的 library size 用于后续的标准化。library size 等于落在该样本所有peaks上的reads的总数，即 counts 矩阵中每列的加和。如果bFullLibrarySize设为TRUE，则会使用library的总reads数（根据原BAM/BED文件统计获得）。然后使用estimateDispersions估计统计分布，需要将参数fitType设为’local’
  dds<-estimateSizeFactors(dds)
  dds<-estimateDispersions(dds,fitType="local")
  dds <- nbinomWaldTest(dds)
  res <- results(dds)
  ## 提取你想要的差异分析结果，我们这里是treated组对untreated组进行比较
  res <-  results(dds, contrast=c("groupList", unique(condition)))
  res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
  length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
  res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                                 ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
  res_DESeq2 = as.data.frame(res)
  table(res_DESeq2$change)
  res_DESeq2 = res_DESeq2[res_DESeq2$change %in% c('UP', 'DOWN'), ]
  # 2. Exploring Method 2: limma-------------------------------------------
  library(edgeR)
  d0 <- DGEList(counts = countMatrix)
  dim(d0)
  # Filtering low count peaks
  cutoff <- 1           
  drop <- which(apply(cpm(d0), 1, max) < cutoff)
  print(length(drop))
  if (length(drop) >0) {
    d <- d0[-drop,] 
  } else{ d <- d0}
  dim(d) #number of peaks left: 8501      
  d <- calcNormFactors(d, method = c('TMM'))
  colnames(countMatrix)
  genotype <- groupList
  table(genotype)
  mm <- model.matrix(~0 + genotype)#
  v <- voom(d, mm, plot=F, normalize="quantile")#voom from limma:Transform count data to log2-counts per million (logtpm)
  fit <- lmFit(v, mm)
  head(coef(fit))
  contr <- makeContrasts(paste0(colnames(coef(fit)), collapse  = ' - '), levels = colnames(coef(fit))) #Change if other comparison
  tmp <- contrasts.fit(fit, contr)
  tmp <- eBayes(tmp)
  top.table <- topTable(tmp, sort.by = "P", n = Inf)
  head(top.table, 20)
  length(which(top.table$adj.P.Val < 0.05 & abs(top.table$logFC) > log2(1.5))) #how many are significant?
  res_limma <- as.data.frame(top.table)
  res_limma$change <- with(res_limma, ifelse(res_limma$adj.P.Val < 0.05 & abs(res_limma$logFC) > log2(1.5),
                                             ifelse(res_limma$logFC > log2(1.5),'UP','DOWN'),'NOT'))
  table(res_limma$change)
  res_limma = res_limma[res_limma$change %in% c('UP', 'DOWN'), ]
  resMeR_list[[c]] <- list(res_DESeq2, res_limma)
  countMatrix_list[[c]] <- countMatrix
  resPeak_list[[c]] <- merge_peak
  name = paste0(c, ' vs WT')
  name = paste0(gsub(' ','_', name ))
  names(resMeR_list[[c]] ) <- paste0(name, '_', c('DESeq2', 'limma'))
}
names(resMeR_list) <- NULL
RMdem_list[[i]] <- as.list(unlist(resMeR_list, recursive = FALSE))
names(RMdem_list[[i]] ) <- paste0(cl, '_', names(RMdem_list[[i]] ))
RMdep_list[[i]] <- resPeak_list 
name <- unique(gsub('(_limma|_DESeq2)', '', names(RMdem_list[[i]] )))
names(RMdep_list[[i]] ) <- name
contrast_IP_list[[i]] <- lapply(countMatrix_list, extract_names)
condition_IP_list[[i]] <- condition_list

############
n = 100
############
#GSE195637# with RIP and Input
############
print(RMPdeg_datasets[n,])# paper used: hg38; I used: hg19 on server 57
i = RMPdeg_datasets[n,'GSE']
cl = c('Hela')#
RMp = c('ALKBH3')
RM = c('m1A', 'm6A')
condition_list = list(m1A_ALKBH3KO = c(rep('ALKBH3KO', 2), rep('WT', 2)),
                      m6A_ALKBH3KO = c(rep('ALKBH3KO', 2), rep('WT', 2)))#
RMP_allDat[grep('GSE195637', RMP_allDat$study_alias )
           ,c("study_alias", "sample_alias", "sample_title")]
# View(allDat[[i]])
skip <- RMP_allDat[grepl('GSE195637', RMP_allDat$study_alias )
                   &  grepl('input', RMP_allDat$sample_title)
                   , "sample_alias"]
paste0(skip$sample_alias, collapse = " ")
skip <- RMP_allDat[grepl('GSE195637', RMP_allDat$study_alias )
                   &  grepl('-IP', RMP_allDat$sample_title)
                   , "sample_alias"]
paste0(skip$sample_alias, collapse = " ")
files <- list.files(file.path(data_directory,i))
files
getwd()
setwd(data_directory)
setwd(i)
ENA <- as.data.frame(RMP_allDat[grepl('GSE195637', RMP_allDat$study_alias )
                                ,c("study_alias", "sample_alias", "sample_title")])
ENA[order(ENA$sample_alias),]

###################################
#  Cell number-normalized (CNN)   #
#      MeRIP INPUT analysis       #
###################################
condition = c(rep('ALKBH3KO', 2), rep('WT', 2))#
c <- unique(condition[!condition %in% 'WT'])
countGTF <- read.table( paste0("metaPlotR/Input_counts_", c, ".txt"),header=T)
rownames(countGTF) <- countGTF$Geneid
countGTF <- countGTF[, !colnames(countGTF) %in% c('Geneid','Chr','Start','End','Strand','Length')]
dim(countGTF) #28265     4
## Sample name
## Sample information
sample=data.frame(
  sample=colnames(countGTF),
  condition=condition)
## SummarizedExperiment input
dds<-DESeqDataSetFromMatrix(countData=countGTF,colData=sample,design=~condition)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 10
table(keep)
dds <- dds[keep,]
## Note on factor levels
dds$condition <- factor(dds$condition, levels = unique(condition))
dds$condition <- relevel(dds$condition, ref = "WT")
# Differential expression analysis-------------------------------------------
dds <- DESeq(dds)

res <-  results(dds, contrast=c("condition",unique(condition)))
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
table(res$change)
table(is.na(res$change))
# resOrdered <- as.data.frame(res[order(res$padj), ])
resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
resInP$geneSymbol <- rownames(resInP)
resInP_list <- list()
name = res@elementMetadata$description[2]
name = paste0(cl, '.m1A_',str_split(name, "condition ", simplify = T)[,2])
resInP_list[[name]] <- resInP
name = res@elementMetadata$description[2]
name = paste0(cl, '.m6A_',str_split(name, "condition ", simplify = T)[,2])
resInP_list[[name]] <- resInP
names(resInP_list ) <- gsub(' ','_', names(resInP_list))

RMdeg_list[[i]] <- resInP_list
names(RMdeg_list[[i]] )
countGTF_list <- list(m1A_ALKBH3KO = colnames(countGTF), m6A_ALKBH3KO = colnames(countGTF)); 
names(countGTF_list) <- paste0(cl, '.', names(countGTF_list))
# condition_list <- list(condition); names(condition_list) <- c
names(condition_list) <- paste0(cl, '.', names(condition_list))
contrast_input_list[[i]] <- countGTF_list
condition_input_list[[i]] <- condition_list

###################################
#  Cell number-normalized (CNN)   #
#   MeRIP RIP analysis - Part1    #
###################################
condition_list = list(m1A_ALKBH3KO = c(rep('ALKBH3KO', 2), rep('WT', 2)),
                      m6A_ALKBH3KO = c(rep('ALKBH3KO', 2), rep('WT', 2)))#
resMeR_list <- vector("list", length(condition_list))
names(resMeR_list) <- names(condition_list)
countMatrix_list <- vector("list", length(condition_list))
names(countMatrix_list) <- names(condition_list)
resPeak_list <- vector("list",length(condition_list))
names(resPeak_list) <- names(condition_list)
for(c in names(condition_list)){
  merge_peak <- read.delim(paste0("metaPlotR/", c, "_newPeakFile.txt"),header=F)
  countMatrix <- read.table(paste0("metaPlotR/counts_", c, ".txt"),header=T)
  table(countMatrix$Geneid == merge_peak[-1,1])
  table(duplicated(merge_peak[-1,1]))
  dim(countMatrix)
  colnames(countMatrix)
  rownames(countMatrix) <- countMatrix$Geneid
  countMatrix <- countMatrix[, !colnames(countMatrix) %in% c('Geneid','Chr','Start','End','Strand','Length')]
  # 1. 做好分组因子即可
  condition <- condition_list[[c]]
  groupList<-factor(condition)
  colData <- data.frame(row.names=colnames(countMatrix), groupList=groupList)
  # 2. Exploring Method 1: DESeq2-------------------------------------------
  # 构建 DESeqDataSet 对象
  dds <- DESeqDataSetFromMatrix(countData = countMatrix,
                                colData = colData,
                                design = ~ groupList)
  ## Pre-filtering
  keep <- rowSums(counts(dds)) >= 1
  table(keep)
  dds <- dds[keep,]
  # 接着计算每个样本的 library size 用于后续的标准化。library size 等于落在该样本所有peaks上的reads的总数，即 counts 矩阵中每列的加和。如果bFullLibrarySize设为TRUE，则会使用library的总reads数（根据原BAM/BED文件统计获得）。然后使用estimateDispersions估计统计分布，需要将参数fitType设为’local’
  dds<-estimateSizeFactors(dds)
  dds<-estimateDispersions(dds,fitType="local")
  dds <- nbinomWaldTest(dds)
  res <- results(dds)
  ## 提取你想要的差异分析结果，我们这里是treated组对untreated组进行比较
  res <-  results(dds, contrast=c("groupList", unique(condition)))
  res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
  length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
  res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                                 ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
  res_DESeq2 = as.data.frame(res)
  table(res_DESeq2$change)
  res_DESeq2 = res_DESeq2[res_DESeq2$change %in% c('UP', 'DOWN'), ]
  # 2. Exploring Method 2: limma-------------------------------------------
  library(edgeR)
  d0 <- DGEList(counts = countMatrix)
  dim(d0)
  # Filtering low count peaks
  cutoff <- 1           
  drop <- which(apply(cpm(d0), 1, max) < cutoff)
  print(length(drop))
  if (length(drop) >0) {
    d <- d0[-drop,] 
  } else{ d <- d0}
  dim(d) #number of peaks left: 8501      
  d <- calcNormFactors(d, method = c('TMM'))
  colnames(countMatrix)
  genotype <- groupList
  table(genotype)
  mm <- model.matrix(~0 + genotype)#
  v <- voom(d, mm, plot=F, normalize="quantile")#voom from limma:Transform count data to log2-counts per million (logtpm)
  fit <- lmFit(v, mm)
  head(coef(fit))
  contr <- makeContrasts(paste0(colnames(coef(fit)), collapse  = ' - '), levels = colnames(coef(fit))) #Change if other comparison
  tmp <- contrasts.fit(fit, contr)
  tmp <- eBayes(tmp)
  top.table <- topTable(tmp, sort.by = "P", n = Inf)
  head(top.table, 20)
  length(which(top.table$adj.P.Val < 0.05 & abs(top.table$logFC) > log2(1.5))) #how many are significant?
  res_limma <- as.data.frame(top.table)
  res_limma$change <- with(res_limma, ifelse(res_limma$adj.P.Val < 0.05 & abs(res_limma$logFC) > log2(1.5),
                                             ifelse(res_limma$logFC > log2(1.5),'UP','DOWN'),'NOT'))
  table(res_limma$change)
  res_limma = res_limma[res_limma$change %in% c('UP', 'DOWN'), ]
  resMeR_list[[c]] <- list(res_DESeq2, res_limma)
  countMatrix_list[[c]] <- countMatrix
  resPeak_list[[c]] <- merge_peak
  name = paste0(c, ' vs WT')
  name = paste0(gsub(' ','_', name ))
  names(resMeR_list[[c]] ) <- paste0(name, '_', c('DESeq2', 'limma'))
}
names(resMeR_list) <- NULL
RMdem_list[[i]] <- as.list(unlist(resMeR_list, recursive = FALSE))
names(RMdem_list[[i]] ) <- paste0(cl, '.', names(RMdem_list[[i]] ))
RMdep_list[[i]] <- resPeak_list 
name <- unique(gsub('(_limma|_DESeq2)', '', names(RMdem_list[[i]] )))
names(RMdep_list[[i]] ) <- name
names(countMatrix_list) <- paste0(cl, '.', names(countMatrix_list))
contrast_IP_list[[i]] <- lapply(countMatrix_list, extract_names)
names(condition_list) <- paste0(cl, '.', names(condition_list))
condition_IP_list[[i]] <- condition_list

############
n = 101
############
#GSE163310# without case RIP, with RNA
############
print(RMPdeg_datasets[n,])#
i = RMPdeg_datasets[n,'GSE']
cl = c('MGC-803')#Gastric cancer
RMp = c('METTL3')
RM = 'm6A'
condition_list = list(METTL3KD =  c(rep('METTL3KD', 1), rep('WT', 1)),#
                      #METTL3WT =  c(rep('METTL3WT', 1), rep('OE', 1)),
                      METTL3MU =  c(rep('METTL3MU', 1), rep('WT', 1))
                     )

RMP_allDat[grep('GSE163310', RMP_allDat$study_alias )
           ,c("study_alias", "sample_alias", "sample_title")]
# View(allDat[[i]])
files <- list.files(file.path(data_directory,i))
files
getwd()
setwd(data_directory)
setwd(i)
###################################
#         RNA-seq analysis        #
###################################
# Reads Per Kilobase of exon per Megabase of library size (RPKM) were calculated using cuffdiff
# 1) RNA-seq
count_files <- files[grepl('_input_RNA_seq.txt.gz', files)]
###################################
# Read and process each count file
count <- lapply(count_files[grepl('shMETTL3|METTL3_mut', count_files)], function(x) {
  # Read data using fread for faster reading
  data <- fread(file.path(data_directory, i, x))
  # Extract required columns
  df_new <- data[, 1:3]
  df_new <- df_new[!duplicated(df_new[,1])]
  # Set column names same as file names
  colnames(df_new)[1] <- c('id')
  return(df_new)
})
# Merge dataframes in countGTF list by 'gene_id'
df_in <- Reduce(function(x, y) merge(x, y, by = 'id', all = TRUE), count)
colSums(df_in[,-1])
table(colSums(df_in[,-1])=="1e+06")
## Gene name
library(data.table)
countGTF=df_in
# id_mapped <- read.table("/data2/rluo4/hg38/gencode.v44lift37.annotation.gene.probeMap",header = T)
id_mapped <- read.table("/data2/rluo4/hg38/gencode.v36.annotation.gene.probeMap",header = T)
id_mapped$id <- str_split(id_mapped$id, '[.]', simplify = T)[,1]
id_mapped <- id_mapped[which(id_mapped$id %in% countGTF$id ),]
dim(id_mapped) #57930     6
table(duplicated(id_mapped$gene))
countGTF <- left_join(id_mapped,countGTF,by="id")
countGTF <- aggregate(x = countGTF[,(ncol(id_mapped)+1):ncol(countGTF)],   #此时exprSet的第三列开始是表达矩阵内容
                      by = list(symbol = countGTF$gene),   #按照相同symbol分组，在组内计算
                      FUN = max) %>%   #原文中是计算最大值（max），也可以计算平均值（mean）或者中位数（median）
  column_to_rownames(var = 'symbol')
dim(countGTF) # 51242    12
###################################
#   more than 1 condition here   #
###################################
resInP_list <- list()
##Gene count matrix
## Sample name
colnames(countGTF)
#####condition_MGC-803_METTL3MU_vs_WT#####
count <- as.matrix(countGTF[, grepl('METTL3_mut|OE_Ctrl', sample_meta$condition)])
fold_change <- countGTF[,1]/countGTF[,2]#cpm[,1]/cpm[,2]
res <- countGTF
res$logFC <- log2(fold_change)
res$FoldChange <- fold_change
res$change <- with(res, ifelse(abs(res$logFC) > log2(1.5),
                               ifelse(res$logFC > log2(1.5),'UP','DOWN'),'NOT'))
table(res$change)
resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
resInP$geneSymbol <- rownames(resInP)
table(resInP$change)
table(resInP$FoldChange %in% c(0, Inf))
table(resInP$logFC %in% c(-Inf, Inf))
# resInP <- resInP[! resInP$FoldChange %in% c(0, Inf), ]
resInP_list[['MGC-803_METTL3KD_vs_WT']] <- resInP
#####condition_MGC-803_METTL3KD_vs_WT#####
count <- as.matrix(countGTF[, grepl('shMETTL3|KD_Ctrl', sample_meta$condition)])
fold_change <- countGTF[,1]/countGTF[,2]#cpm[,1]/cpm[,2]
res <- countGTF
res$logFC <- log2(fold_change)
res$FoldChange <- fold_change
res$change <- with(res, ifelse(abs(res$logFC) > log2(1.5),
                               ifelse(res$logFC > log2(1.5),'UP','DOWN'),'NOT'))
table(res$change)
resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
resInP$geneSymbol <- rownames(resInP)
table(resInP$change)
table(resInP$FoldChange %in% c(0, Inf))
table(resInP$logFC %in% c(-Inf, Inf))
# resInP <- resInP[! resInP$FoldChange %in% c(0, Inf), ]
resInP_list[['MGC-803_METTL3KD_vs_WT']] <- resInP

RMdeg_list[[i]] <- resInP_list
names(RMdeg_list[[i]] )
contrast_input_list[[i]] <- NULL#countGTF_list
condition_input_list[[i]] <- condition_list



############
n = 102
############
#GSE226129# with RIP and Input
############
print(RMPdeg_datasets[n,])#
i = RMPdeg_datasets[n,'GSE']
cl = c('HCT15')#Gastric cancer
RMp = c('NSUN2' )
RM = 'm5C'
condition = c(rep('NSUN2KD', 2), rep('WT', 2))#
RMP_allDat[grep('GSE226129', RMP_allDat$study_alias )
           ,c("study_alias", "sample_alias", "sample_title")]
# View(allDat[[i]])
skip <- RMP_allDat[grepl('GSE226129', RMP_allDat$study_alias )
                   &  grepl('Input', RMP_allDat$sample_title)
                   , "sample_alias"]
paste0(skip$sample_alias, collapse = " ")
skip <- RMP_allDat[grepl('GSE226129', RMP_allDat$study_alias )
                   &  grepl('IP', RMP_allDat$sample_title)
                   , "sample_alias"]
paste0(skip$sample_alias, collapse = " ")
ENA <- as.data.frame(RMP_allDat[grepl('GSE226129', RMP_allDat$study_alias )
                                ,c("study_alias", "sample_alias", "sample_title")])
ENA[order(ENA$sample_alias),]

files <- list.files(file.path(data_directory,i))
files
getwd()
setwd(data_directory)
setwd(i)
c <- unique(condition[!condition %in% 'WT'])
###################################
#  Cell number-normalized (CNN)   #
#      MeRIP INPUT analysis       #
###################################
countGTF <- read.table( paste0("metaPlotR/Input_counts_", c, ".txt"),header=T)
rownames(countGTF) <- countGTF$Geneid
countGTF <- countGTF[, !colnames(countGTF) %in% c('Geneid','Chr','Start','End','Strand','Length')]
dim(countGTF) #28265     4
## Sample name
## Sample information
sample=data.frame(
  sample=colnames(countGTF),
  condition=condition)
## SummarizedExperiment input
dds<-DESeqDataSetFromMatrix(countData=countGTF,colData=sample,design=~condition)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 10
table(keep)
dds <- dds[keep,]
## Note on factor levels
dds$condition <- factor(dds$condition, levels = unique(condition))
dds$condition <- relevel(dds$condition, ref = "WT")
# Differential expression analysis-------------------------------------------
dds <- DESeq(dds)

res <-  results(dds, contrast=c("condition",unique(condition)))
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
table(res$change)
table(is.na(res$change))
# resOrdered <- as.data.frame(res[order(res$padj), ])
resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
resInP$geneSymbol <- rownames(resInP)
name = res@elementMetadata$description[2]
name = paste0(cl, '_',str_split(name, "condition ", simplify = T)[,2])
RMdeg_list[[i]] <- list( resInP)
names(RMdeg_list[[i]] ) <- paste0(gsub(' ','_',(name)))
countGTF_list <- list(colnames(countGTF)); names(countGTF_list) <- c
condition_list <- list(condition); names(condition_list) <- c
contrast_input_list[[i]] <- countGTF_list
condition_input_list[[i]] <- condition_list

###################################
#  Cell number-normalized (CNN)   #
#   MeRIP RIP analysis - Part1    #
###################################
## Sample information
merge_peak <- read.delim(paste0("metaPlotR/", c, "_newPeakFile.txt"),header=F)
countMatrix <- read.table(paste0("metaPlotR/m5C_counts_", c, ".txt"),header=T)
table(countMatrix$Geneid == merge_peak[-1,1])
table(duplicated(merge_peak$V1[-1]))
dim(countMatrix)
colnames(countMatrix)
rownames(countMatrix) <- countMatrix$Geneid
countMatrix <- countMatrix[, !colnames(countMatrix) %in% c('Geneid','Chr','Start','End','Strand','Length')]
# 1. 做好分组因子即可
groupList<-factor(condition)
colData <- data.frame(row.names=colnames(countMatrix), groupList=groupList)
# 2. Exploring Method 1: DESeq2-------------------------------------------
# 构建 DESeqDataSet 对象
dds <- DESeqDataSetFromMatrix(countData = countMatrix,
                              colData = colData,
                              design = ~ groupList)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 1
table(keep)
dds <- dds[keep,]
# 接着计算每个样本的 library size 用于后续的标准化。library size 等于落在该样本所有peaks上的reads的总数，即 counts 矩阵中每列的加和。如果bFullLibrarySize设为TRUE，则会使用library的总reads数（根据原BAM/BED文件统计获得）。然后使用estimateDispersions估计统计分布，需要将参数fitType设为’local’
dds<-estimateSizeFactors(dds)
dds<-estimateDispersions(dds,fitType="local")
dds <- nbinomWaldTest(dds)
res <- results(dds)
## 提取你想要的差异分析结果，我们这里是treated组对untreated组进行比较
res <-  results(dds, contrast=c("groupList", unique(condition)))
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
res_DESeq2 = as.data.frame(res)
table(res_DESeq2$change)#Not only
res_DESeq2 = res_DESeq2[res_DESeq2$change %in% c('UP', 'DOWN'), ]
# 2. Exploring Method 2: limma-------------------------------------------
library(edgeR)
d0 <- DGEList(counts = countMatrix)
dim(d0)
# Filtering low count peaks
cutoff <- 1           
drop <- which(apply(cpm(d0), 1, max) < cutoff)
print(length(drop))
if (length(drop) >0) {
  d <- d0[-drop,] 
} else{ d <- d0}
dim(d) #number of peaks left: 8501      
d <- calcNormFactors(d, method = c('TMM'))
colnames(countMatrix)
genotype <- groupList
table(genotype)
mm <- model.matrix(~0 + genotype)#
v <- voom(d, mm, plot=F, normalize="quantile")#voom from limma:Transform count data to log2-counts per million (logtpm)
fit <- lmFit(v, mm)
head(coef(fit))
contr <- makeContrasts(paste0(colnames(coef(fit)), collapse  = ' - '), levels = colnames(coef(fit))) #Change if other comparison
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)
length(which(top.table$adj.P.Val < 0.05 & abs(top.table$logFC) > log2(1.5))) #how many are significant?
res_limma <- as.data.frame(top.table)
res_limma$change <- with(res_limma, ifelse(res_limma$adj.P.Val < 0.05 & abs(res_limma$logFC) > log2(1.5),
                                           ifelse(res_limma$logFC > log2(1.5),'UP','DOWN'),'NOT'))
table(res_limma$change)#NOT only
res_limma = res_limma[res_limma$change %in% c('UP', 'DOWN'), ]

RMdem_list[[i]] <- list(res_DESeq2, res_limma)
name = res@elementMetadata$description[2]
name = paste0(cl, '_',str_split(name, "groupList ", simplify = T)[,2])
name= paste0(gsub(' ','_',(name)))
names(RMdem_list[[i]] ) <- paste0(name, '_', c('DESeq2', 'limma'))
RMdep_list[[i]] <- list( merge_peak )
names(RMdep_list[[i]] ) <- name
countMatrix_list <- list(colnames(countMatrix)); names(countMatrix_list) <- c
condition_list <- list(condition); names(condition_list) <- c
contrast_IP_list[[i]] <- countMatrix_list
condition_IP_list[[i]] <- condition_list


############
n = 103
############
#GSE223731# with RIP and Input
############
print(RMPdeg_datasets[n,])# hg38
i = RMPdeg_datasets[n,'GSE']
cl = c('NALM6 and Mino and JeKo-1 and HG3 and MEC1 (ACC497, DSMZ) and MOLM13')#
cl = 'HG3'
RMp = c('METTL3')
RM = 'm6A'
condition = c(rep('METTL3KO', 2), rep('WT', 2))#
RMP_allDat[grep('GSE223728', RMP_allDat$study_alias )
           ,c("study_alias", "sample_alias", "sample_title")]
# View(allDat[[i]])
RMP_allDat[grep('GSE223729', RMP_allDat$study_alias )
           ,c("study_alias", "sample_alias", "sample_title")]
skip <- RMP_allDat[grepl('GSE223729', RMP_allDat$study_alias )
                   &  grepl('Input', RMP_allDat$sample_title)
                   , "sample_alias"]
paste0(skip$sample_alias, collapse = " ")
skip <- RMP_allDat[grepl('GSE223729', RMP_allDat$study_alias )
                   &  grepl('MeRIP', RMP_allDat$sample_title)
                   , "sample_alias"]
paste0(skip$sample_alias, collapse = " ")
ENA <- as.data.frame(RMP_allDat[grepl('GSE223729', RMP_allDat$study_alias )
                                ,c("study_alias", "sample_alias", "sample_title")])
ENA[order(ENA$sample_alias),]

files <- list.files(file.path(data_directory,i))
files
getwd()
setwd(data_directory)
setwd(i)
c <- unique(condition[!condition %in% 'WT'])
###################################
#  Cell number-normalized (CNN)   #
#      MeRIP INPUT analysis       #
###################################
countGTF <- read.table( paste0("metaPlotR/Input_counts_", c, ".txt"),header=T)
rownames(countGTF) <- countGTF$Geneid
countGTF <- countGTF[, !colnames(countGTF) %in% c('Geneid','Chr','Start','End','Strand','Length')]
dim(countGTF) #28265     4
## Sample name
## Sample information
sample=data.frame(
  sample=colnames(countGTF),
  condition=condition)
## SummarizedExperiment input
dds<-DESeqDataSetFromMatrix(countData=countGTF,colData=sample,design=~condition)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 10
table(keep)
dds <- dds[keep,]
## Note on factor levels
dds$condition <- factor(dds$condition, levels = unique(condition))
dds$condition <- relevel(dds$condition, ref = "WT")
# Differential expression analysis-------------------------------------------
dds <- DESeq(dds)

res <-  results(dds, contrast=c("condition",unique(condition)))
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
table(res$change)
table(is.na(res$change))
# resOrdered <- as.data.frame(res[order(res$padj), ])
resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
resInP$geneSymbol <- rownames(resInP)
name = res@elementMetadata$description[2]
name = paste0(cl, '_',str_split(name, "condition ", simplify = T)[,2])
RMdeg_list[[i]] <- list( resInP)
names(RMdeg_list[[i]] ) <- paste0(gsub(' ','_',(name)))
countGTF_list <- list(colnames(countGTF)); names(countGTF_list) <- c
condition_list <- list(condition); names(condition_list) <- c
contrast_input_list[[i]] <- countGTF_list
condition_input_list[[i]] <- condition_list

###################################
#  Cell number-normalized (CNN)   #
#   MeRIP RIP analysis - Part1    #
###################################
## Sample information
merge_peak <- read.delim(paste0("metaPlotR/", c, "_newPeakFile.txt"),header=F)
countMatrix <- read.table(paste0("metaPlotR/m6A_counts_", c, ".txt"),header=T)
table(countMatrix$Geneid == merge_peak[-1,1])
table(duplicated(merge_peak$V1[-1]))
dim(countMatrix)
colnames(countMatrix)
rownames(countMatrix) <- countMatrix$Geneid
countMatrix <- countMatrix[, !colnames(countMatrix) %in% c('Geneid','Chr','Start','End','Strand','Length')]
# 1. 做好分组因子即可
groupList<-factor(condition)
colData <- data.frame(row.names=colnames(countMatrix), groupList=groupList)
# 2. Exploring Method 1: DESeq2-------------------------------------------
# 构建 DESeqDataSet 对象
dds <- DESeqDataSetFromMatrix(countData = countMatrix,
                              colData = colData,
                              design = ~ groupList)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 1
table(keep)
dds <- dds[keep,]
# 接着计算每个样本的 library size 用于后续的标准化。library size 等于落在该样本所有peaks上的reads的总数，即 counts 矩阵中每列的加和。如果bFullLibrarySize设为TRUE，则会使用library的总reads数（根据原BAM/BED文件统计获得）。然后使用estimateDispersions估计统计分布，需要将参数fitType设为’local’
dds<-estimateSizeFactors(dds)
dds<-estimateDispersions(dds,fitType="local")
dds <- nbinomWaldTest(dds)
res <- results(dds)
## 提取你想要的差异分析结果，我们这里是treated组对untreated组进行比较
res <-  results(dds, contrast=c("groupList", unique(condition)))
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
res_DESeq2 = as.data.frame(res)
table(res_DESeq2$change)#Not only
res_DESeq2 = res_DESeq2[res_DESeq2$change %in% c('UP', 'DOWN'), ]
# 2. Exploring Method 2: limma-------------------------------------------
library(edgeR)
d0 <- DGEList(counts = countMatrix)
dim(d0)
# Filtering low count peaks
cutoff <- 1           
drop <- which(apply(cpm(d0), 1, max) < cutoff)
print(length(drop))
if (length(drop) >0) {
  d <- d0[-drop,] 
} else{ d <- d0}
dim(d) #number of peaks left: 8501      
d <- calcNormFactors(d, method = c('TMM'))
colnames(countMatrix)
genotype <- groupList
table(genotype)
mm <- model.matrix(~0 + genotype)#
v <- voom(d, mm, plot=F, normalize="quantile")#voom from limma:Transform count data to log2-counts per million (logtpm)
fit <- lmFit(v, mm)
head(coef(fit))
contr <- makeContrasts(paste0(colnames(coef(fit)), collapse  = ' - '), levels = colnames(coef(fit))) #Change if other comparison
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)
length(which(top.table$adj.P.Val < 0.05 & abs(top.table$logFC) > log2(1.5))) #how many are significant?
res_limma <- as.data.frame(top.table)
res_limma$change <- with(res_limma, ifelse(res_limma$adj.P.Val < 0.05 & abs(res_limma$logFC) > log2(1.5),
                                           ifelse(res_limma$logFC > log2(1.5),'UP','DOWN'),'NOT'))
table(res_limma$change)#NOT only
res_limma = res_limma[res_limma$change %in% c('UP', 'DOWN'), ]

RMdem_list[[i]] <- list(res_DESeq2, res_limma)
name = res@elementMetadata$description[2]
name = paste0(cl, '_',str_split(name, "groupList ", simplify = T)[,2])
name= paste0(gsub(' ','_',(name)))
names(RMdem_list[[i]] ) <- paste0(name, '_', c('DESeq2', 'limma'))
RMdep_list[[i]] <- list( merge_peak )
names(RMdep_list[[i]] ) <- name
countMatrix_list <- list(colnames(countMatrix)); names(countMatrix_list) <- c
condition_list <- list(condition); names(condition_list) <- c
contrast_IP_list[[i]] <- countMatrix_list
condition_IP_list[[i]] <- condition_list


############
n = 104
############
#GSE242276# with RIP and Input
############
print(RMPdeg_datasets[n,])# hg38
#Head and Neck Squamous Carcinoma
i = RMPdeg_datasets[n,'GSE']
cl = c('UM-SCC-1' )#and UM-SCC-2 and SCC090 and JHU022 and 93-VU-147T and SqCC/Y1 and 1483')#
RMp = c('ALKBH5') 
RM = 'm6A'
condition_list <- list(
  'UM-SCC-1.1_ALKBH5KD' = c(rep('ALKBH5KD', 3), rep('WT', 3)),
  'UM-SCC-1.2_ALKBH5KD' = c(rep('ALKBH5KD', 3), rep('WT', 3))#
)
condition = c(rep('ALKBH5KD', 6), rep('WT', 3))
RMP_allDat[grep('GSE242276', RMP_allDat$study_alias ) 
           ,c("study_alias", "sample_alias", "sample_title")]
# View(allDat[[i]])
RMP_allDat[grepl('GSE242276', RMP_allDat$study_alias ) & grepl('ALKBH5', RMP_allDat$sample_title)
           ,c("study_alias", "sample_alias", "sample_title")]
skip = RMP_allDat[grepl('GSE242276', RMP_allDat$study_alias ) & grepl('ALKBH5', RMP_allDat$sample_title)
                  ,c("study_alias", "sample_alias", "sample_title")]
paste0(skip$sample_alias, collapse = "|")

skip <- RMP_allDat[grepl('GSE242276', RMP_allDat$study_alias ) &  grepl('ALKBH5', RMP_allDat$sample_title)
                   &  grepl('RNA-seq', RMP_allDat$sample_title)
                   , "sample_alias"]
paste0(skip$sample_alias, collapse = " ")
skip <- RMP_allDat[grepl('GSE242276', RMP_allDat$study_alias ) & grepl('ALKBH5', RMP_allDat$sample_title)
                   &  grepl('m6A-seq', RMP_allDat$sample_title)
                   , "sample_alias"]
paste0(skip$sample_alias, collapse = " ")
files <- list.files(file.path(data_directory,i))
files
getwd()
setwd(data_directory)
setwd(i)
ENA <- as.data.frame(RMP_allDat[grepl('GSE242276', RMP_allDat$study_alias ) & ( !grepl('RBM33', RMP_allDat$sample_title))
                                ,c("study_alias", "sample_alias", "sample_title")])
ENA <- ENA[order(ENA$sample_alias),]
ENA

countfiles <- list.files('./metaPlotR/')
countfiles <- countfiles[grep('.txt', countfiles)]
countfiles <- countfiles[! grepl('summary', countfiles)] #mouse embryonic fibroblasts (MEFs)
countfiles
# Method1: 
c <- unique(condition[!condition %in% 'WT'])
###################################
#  Cell number-normalized (CNN)   #
#      MeRIP INPUT analysis       #
###################################
countGTF <- read.table( paste0("metaPlotR/Input_counts_", c, ".txt"),header=T)
rownames(countGTF) <- countGTF$Geneid
countGTF <- countGTF[, !colnames(countGTF) %in% c('Geneid','Chr','Start','End','Strand','Length')]
dim(countGTF) #28265     4
## Sample name
## Sample information
sample=data.frame(
  sample=colnames(countGTF),
  condition=condition)
## SummarizedExperiment input
dds<-DESeqDataSetFromMatrix(countData=countGTF,colData=sample,design=~condition)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 10
table(keep)
dds <- dds[keep,]
## Note on factor levels
dds$condition <- factor(dds$condition, levels = unique(condition))
dds$condition <- relevel(dds$condition, ref = "WT")
# Differential expression analysis-------------------------------------------
dds <- DESeq(dds)

res <-  results(dds, contrast=c("condition",unique(condition)))
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
table(res$change)
table(is.na(res$change))
# resOrdered <- as.data.frame(res[order(res$padj), ])
resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
resInP$geneSymbol <- rownames(resInP)
name = res@elementMetadata$description[2]
name = paste0(cl, '_',str_split(name, "condition ", simplify = T)[,2])
RMdeg_list[[i]] <- list( resInP)
names(RMdeg_list[[i]] ) <- paste0(gsub(' ','_',(name)))
countGTF_list <- list(colnames(countGTF)); names(countGTF_list) <- c
condition_list <- list(condition); names(condition_list) <- c
contrast_input_list[[i]] <- countGTF_list
condition_input_list[[i]] <- condition_list

###################################
#  Cell number-normalized (CNN)   #
#   MeRIP RIP analysis - Part1    #
###################################
## Sample information
merge_peak <- read.delim(paste0("metaPlotR/", c, "_newPeakFile.txt"),header=F)
countMatrix <- read.table(paste0("metaPlotR/m6A_counts_", c, ".txt"),header=T)
table(countMatrix$Geneid == merge_peak[-1,1])
table(duplicated(merge_peak$V1[-1]))
dim(countMatrix)
colnames(countMatrix)
rownames(countMatrix) <- countMatrix$Geneid
countMatrix <- countMatrix[, !colnames(countMatrix) %in% c('Geneid','Chr','Start','End','Strand','Length')]
# 1. 做好分组因子即可
groupList<-factor(condition)
colData <- data.frame(row.names=colnames(countMatrix), groupList=groupList)
# 2. Exploring Method 1: DESeq2-------------------------------------------
# 构建 DESeqDataSet 对象
dds <- DESeqDataSetFromMatrix(countData = countMatrix,
                              colData = colData,
                              design = ~ groupList)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 1
table(keep)
dds <- dds[keep,]
# 接着计算每个样本的 library size 用于后续的标准化。library size 等于落在该样本所有peaks上的reads的总数，即 counts 矩阵中每列的加和。如果bFullLibrarySize设为TRUE，则会使用library的总reads数（根据原BAM/BED文件统计获得）。然后使用estimateDispersions估计统计分布，需要将参数fitType设为’local’
dds<-estimateSizeFactors(dds)
dds<-estimateDispersions(dds,fitType="local")
dds <- nbinomWaldTest(dds)
res <- results(dds)
## 提取你想要的差异分析结果，我们这里是treated组对untreated组进行比较
res <-  results(dds, contrast=c("groupList", unique(condition)))
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
res_DESeq2 = as.data.frame(res)
table(res_DESeq2$change)#Not only
res_DESeq2 = res_DESeq2[res_DESeq2$change %in% c('UP', 'DOWN'), ]
# 2. Exploring Method 2: limma-------------------------------------------
library(edgeR)
d0 <- DGEList(counts = countMatrix)
dim(d0)
# Filtering low count peaks
cutoff <- 1           
drop <- which(apply(cpm(d0), 1, max) < cutoff)
print(length(drop))
if (length(drop) >0) {
  d <- d0[-drop,] 
} else{ d <- d0}
dim(d) #number of peaks left: 8501      
d <- calcNormFactors(d, method = c('TMM'))
colnames(countMatrix)
genotype <- groupList
table(genotype)
mm <- model.matrix(~0 + genotype)#
v <- voom(d, mm, plot=F, normalize="quantile")#voom from limma:Transform count data to log2-counts per million (logtpm)
fit <- lmFit(v, mm)
head(coef(fit))
contr <- makeContrasts(paste0(colnames(coef(fit)), collapse  = ' - '), levels = colnames(coef(fit))) #Change if other comparison
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)
length(which(top.table$adj.P.Val < 0.05 & abs(top.table$logFC) > log2(1.5))) #how many are significant?
res_limma <- as.data.frame(top.table)
res_limma$change <- with(res_limma, ifelse(res_limma$adj.P.Val < 0.05 & abs(res_limma$logFC) > log2(1.5),
                                           ifelse(res_limma$logFC > log2(1.5),'UP','DOWN'),'NOT'))
table(res_limma$change)#NOT only
res_limma = res_limma[res_limma$change %in% c('UP', 'DOWN'), ]

RMdem_list[[i]] <- list(res_DESeq2, res_limma)
name = res@elementMetadata$description[2]
name = paste0(cl, '_',str_split(name, "groupList ", simplify = T)[,2])
name= paste0(gsub(' ','_',(name)))
names(RMdem_list[[i]] ) <- paste0(name, '_', c('DESeq2', 'limma'))
RMdep_list[[i]] <- list( merge_peak )
names(RMdep_list[[i]] ) <- name
countMatrix_list <- list(colnames(countMatrix)); names(countMatrix_list) <- c
condition_list <- list(condition); names(condition_list) <- c
contrast_IP_list[[i]] <- countMatrix_list
condition_IP_list[[i]] <- condition_list
# # Method2: 
# ###################################
# #  Cell number-normalized (CNN)   #
# #      MeRIP INPUT analysis       #
# ###################################
# resInP_list <- vector("list", length(condition_list))
# names(resInP_list) <- names(condition_list)
# countGTF_list <- vector("list", length(condition_list))
# names(countGTF_list) <- names(condition_list)
# for(c in names(condition_list)){
#   countGTF <- read.table( paste0("metaPlotR/Input_counts_", c, ".txt"),header=T)
#   dim(countGTF)
#   colnames(countGTF)
#   rownames(countGTF) <- countGTF$Geneid
#   countGTF <- countGTF[, !colnames(countGTF) %in% c('Geneid','Chr','Start','End','Strand','Length')]
#   dim(countGTF) # 28265     4
#   condition <- condition_list[[c]]
#   ## Sample information
#   sample=data.frame(
#     sample=colnames(countGTF),
#     condition=condition)
#   ## SummarizedExperiment input
#   dds<-DESeqDataSetFromMatrix(countData=countGTF,colData=sample,design=~condition)
#   ## Pre-filtering
#   keep <- rowSums(counts(dds)) >= 10
#   table(keep)
#   dds <- dds[keep,]
#   ## Note on factor levels
#   dds$condition <- factor(dds$condition, levels = unique(condition))
#   dds$condition <- relevel(dds$condition, ref = "WT")
#   # Differential expression analysis-------------------------------------------
#   dds <- DESeq(dds)
#   res <-  results(dds, contrast=c("condition",unique(condition)))
#   res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
#   length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
#   res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
#                                  ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
#   table(res$change)
#   table(is.na(res$change))
#   # resOrdered <- as.data.frame(res[order(res$padj), ])
#   resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
#   resInP$geneSymbol <- rownames(resInP)
#   table(resInP$change)
#   resInP_list[[c]] <- resInP
#   countGTF_list[[c]] <- countGTF
# }
# name = paste0(names(resInP_list), ' vs WT')
# # name = paste0(cl, '_', name)
# RMdeg_list[[i]] <- resInP_list 
# names(RMdeg_list[[i]] ) <- paste0(gsub(' ','_',(name)))
# contrast_input_list[[i]] <- lapply(countGTF_list, extract_names)
# condition_input_list[[i]] <- condition_list
# 
# ###################################
# #  Cell number-normalized (CNN)   #
# #   MeRIP RIP analysis - Part1    #
# ###################################
# resMeR_list <- vector("list", length(condition_list))
# names(resMeR_list) <- names(condition_list)
# countMatrix_list <- vector("list", length(condition_list))
# names(countMatrix_list) <- names(condition_list)
# resPeak_list <- vector("list",length(condition_list))
# names(resPeak_list) <- names(condition_list)
# 
# for(c in names(condition_list)){
#   merge_peak <- read.delim(paste0("metaPlotR/", c, "_newPeakFile.txt"),header=F)
#   countMatrix <- read.table(paste0("metaPlotR/m6A_counts_", c, ".txt"),header=T)
#   table(countMatrix$Geneid == merge_peak[-1,1])
#   table(duplicated(merge_peak[-1,1]))
#   dim(countMatrix)
#   colnames(countMatrix)
#   rownames(countMatrix) <- countMatrix$Geneid
#   countMatrix <- countMatrix[, !colnames(countMatrix) %in% c('Geneid','Chr','Start','End','Strand','Length')]
#   # 1. 做好分组因子即可
#   condition <- condition_list[[c]]
#   groupList<-factor(condition)
#   colData <- data.frame(row.names=colnames(countMatrix), groupList=groupList)
#   # 2. Exploring Method 1: DESeq2-------------------------------------------
#   # 构建 DESeqDataSet 对象
#   dds <- DESeqDataSetFromMatrix(countData = countMatrix,
#                                 colData = colData,
#                                 design = ~ groupList)
#   ## Pre-filtering
#   keep <- rowSums(counts(dds)) >= 1
#   table(keep)
#   dds <- dds[keep,]
#   # 接着计算每个样本的 library size 用于后续的标准化。library size 等于落在该样本所有peaks上的reads的总数，即 counts 矩阵中每列的加和。如果bFullLibrarySize设为TRUE，则会使用library的总reads数（根据原BAM/BED文件统计获得）。然后使用estimateDispersions估计统计分布，需要将参数fitType设为’local’
#   dds<-estimateSizeFactors(dds)
#   dds<-estimateDispersions(dds,fitType="local")
#   dds <- nbinomWaldTest(dds)
#   res <- results(dds)
#   ## 提取你想要的差异分析结果，我们这里是treated组对untreated组进行比较
#   res <-  results(dds, contrast=c("groupList", unique(condition)))
#   res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
#   length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
#   res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
#                                  ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
#   res_DESeq2 = as.data.frame(res)
#   table(res_DESeq2$change)
#   res_DESeq2 = res_DESeq2[res_DESeq2$change %in% c('UP', 'DOWN'), ]
#   # 2. Exploring Method 2: limma-------------------------------------------
#   library(edgeR)
#   d0 <- DGEList(counts = countMatrix)
#   dim(d0)
#   # Filtering low count peaks
#   cutoff <- 1           
#   drop <- which(apply(cpm(d0), 1, max) < cutoff)
#   print(length(drop))
#   if (length(drop) >0) {
#     d <- d0[-drop,] 
#   } else{ d <- d0}
#   dim(d) #number of peaks left: 8501      
#   d <- calcNormFactors(d, method = c('TMM'))
#   colnames(countMatrix)
#   genotype <- groupList
#   table(genotype)
#   mm <- model.matrix(~0 + genotype)#
#   v <- voom(d, mm, plot=F, normalize="quantile")#voom from limma:Transform count data to log2-counts per million (logtpm)
#   fit <- lmFit(v, mm)
#   head(coef(fit))
#   contr <- makeContrasts(paste0(colnames(coef(fit)), collapse  = ' - '), levels = colnames(coef(fit))) #Change if other comparison
#   tmp <- contrasts.fit(fit, contr)
#   tmp <- eBayes(tmp)
#   top.table <- topTable(tmp, sort.by = "P", n = Inf)
#   head(top.table, 20)
#   length(which(top.table$adj.P.Val < 0.05 & abs(top.table$logFC) > log2(1.5))) #how many are significant?
#   res_limma <- as.data.frame(top.table)
#   res_limma$change <- with(res_limma, ifelse(res_limma$adj.P.Val < 0.05 & abs(res_limma$logFC) > log2(1.5),
#                                              ifelse(res_limma$logFC > log2(1.5),'UP','DOWN'),'NOT'))
#   table(res_limma$change)
#   res_limma = res_limma[res_limma$change %in% c('UP', 'DOWN'), ]
#   resMeR_list[[c]] <- list(res_DESeq2, res_limma)
#   countMatrix_list[[c]] <- countMatrix
#   resPeak_list[[c]] <- merge_peak
#   name = res@elementMetadata$description[2]
#   # name = paste0(cl, '_',str_split(name, "groupList ", simplify = T)[,2])
#   name = str_split(name, "groupList ", simplify = T)[,2]
#   name = paste0(gsub(' ','_',(name)))
#   names(resMeR_list[[c]] ) <- paste0(name, '_', c('DESeq2', 'limma'))
# }
# RMdem_list[[i]] <- as.list(unlist(resMeR_list, recursive = FALSE))
# name <- gsub("_ALKBH5KD.", '_', names(RMdem_list[[i]] ))#str_split(names(RMdem_list[[i]] ), "[.]", simplify = T)[,2]
# names(RMdem_list[[i]] ) <- name
# RMdep_list[[i]] <- resPeak_list 
# name <- unique(gsub('(_limma|_DESeq2)', '', name))
# names(RMdep_list[[i]] ) <- name
# contrast_IP_list[[i]] <- lapply(countMatrix_list, extract_names)
# condition_IP_list[[i]] <- condition_list
############
n = 105
############
#GSE252752# without case RIP, with RNA
############
print(RMPdeg_datasets[n,])# undone !
i = RMPdeg_datasets[n,'GSE']
cl = c('MHCC97H/Huh7/Hep3B')#Liver cancer
RMp = c('YTHDF1')
RM = 'm6A'
condition = c(rep('YTHDF1KD', 2), rep('WT', 2))#
RMP_allDat[grep('GSE252752', RMP_allDat$study_alias )
           ,c("study_alias", "sample_alias", "sample_title")]
# View(allDat[[i]])
files <- list.files(file.path(data_directory,i))
files
getwd()
setwd(data_directory)
setwd(i)

df_in <- read_excel('GSE252752_processed_data_file.xlsx', sheet = 2)
df_in <- read_excel('GSE252752_processed_data_file.xlsx', sheet = 3)
countGTF <- as.data.frame(df_in[, 12:15])
table(duplicated(df_in$Symbol))
tmp=by(countGTF, countGTF$Symbol,function(x) rownames(x)[which.min(x$FDR)])
probes = as.character(tmp)
length(probes) # 7785
res <- countGTF[rownames(countGTF) %in% probes,]
colnames(res) <- c('log2FoldChange', 'pvalue', 'padj', 'geneSymbol')
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$label <- with(res, ifelse(abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
table(res$change)
resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
table(duplicated(resInP$geneSymbol))
name = "MHCC97H_YTHDF1KD_vs_WT"
RMdeg_list[[i]] <- list( resInP)
names(RMdeg_list[[i]] ) <- name
condition_list <- list(condition); names(condition_list) <- 'YTHDF1KD'
contrast_input_list[[i]] <- NULL#countGTF_list
condition_input_list[[i]] <- condition_list

############
n = 106
############
#GSE254232# without case RIP, with RNA
############
print(RMPdeg_datasets[n,])#
i = RMPdeg_datasets[n,'GSE']
cl = c('EC109')#Esophageal squamous cell carcinoma
RMp = c('METTL3', 'METTL14')
RM = 'm6A'
condition_list = list(EC109_METTL3KD =  c(rep('METTL3KD', 2), rep('WT', 2)),
                      EC109.hs_METTL3KD =  c(rep('METTL3KD', 2), rep('WT', 2)),
                      EC109_METTL14KD =  c(rep('METTL14KD', 2), rep('WT', 2)))#
RMP_allDat[grep('GSE254232', RMP_allDat$study_alias )
           ,c("study_alias", "sample_alias", "sample_title")]
# View(allDat[[i]])
files <- list.files(file.path(data_directory,i))
files
getwd()
setwd(data_directory)
setwd(i)
df_in <- fread('GSE254232_gene_count.xls.gz')
###################################
#   more than 1 condition here   #
###################################
resInP_list <- list()
##Gene count matrix
countGTF=as.data.frame(df_in[, 2:12])
## Gene name
table(duplicated(df_in$gene_name))
countGTF <- aggregate(x = countGTF[,-ncol(countGTF)],   #此时exprSet的第三列开始是表达矩阵内容
                      by = list(symbol = countGTF$gene_name),   #按照相同symbol分组，在组内计算
                      FUN = max) %>%   #原文中是计算最大值（max），也可以计算平均值（mean）或者中位数（median）
  column_to_rownames(var = 'symbol')
dim(countGTF) # 51242    12
## Sample name
colnames(countGTF)
## Sample information
sample_meta=data.frame(
  sample=colnames(countGTF),
  condition=c(rep("EC109_WT",2),rep("EC109.hs_WT",2), rep("EC109_METTL3KD",2), rep("EC109_METTL14KD",2), rep("EC109.hs_METTL3KD",2)))
#####condition_EC109_METTL3KD_vs_WT#####
count <- as.matrix(countGTF[, grepl('EC109_WT|EC109_METTL3KD', sample_meta$condition)])
table(is.na(count))
sample <- sample_meta[sample_meta$sample %in% colnames(count), ]
## SummarizedExperiment input
condition = c(rep('EC109_METTL3KD', 2), rep('EC109_WT', 2))
dds<-DESeqDataSetFromMatrix(countData=count, colData=sample,design=~condition)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 10
table(keep)
dds <- dds[keep,]
## Note on factor levels
dds$condition <- factor(dds$condition, levels = c("EC109_METTL3KD","EC109_WT"))
dds$condition <- relevel(dds$condition, ref = "EC109_WT")
dds <- DESeq(dds)
res <-  results(dds, contrast=c("condition",unique(condition))) # c("condition", rev(levels(genotype)))
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
table(res$change)
res$geneSymbol <- rownames(res)
resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
table(duplicated(resInP$geneSymbol))
resInP_list[['EC109_METTL3KD_vs_WT']] <- resInP
#####condition_EC109.heatshock_METTL3KD_vs_WT#####
count <- as.matrix(countGTF[, grepl('EC109.hs_WT|EC109.hs_METTL3KD', sample_meta$condition)])
table(is.na(count))
sample <- sample_meta[sample_meta$sample %in% colnames(count), ]
## SummarizedExperiment input
condition = c(rep('EC109.hs_METTL3KD', 2), rep('EC109.hs_WT', 2))
dds<-DESeqDataSetFromMatrix(countData=count, colData=sample,design=~condition)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 10
table(keep)
dds <- dds[keep,]
## Note on factor levels
dds$condition <- factor(dds$condition, levels = c("EC109.hs_METTL3KD","EC109.hs_WT"))
dds$condition <- relevel(dds$condition, ref = "EC109.hs_WT")
dds <- DESeq(dds)
res <-  results(dds, contrast=c("condition",unique(condition))) # c("condition", rev(levels(genotype)))
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
table(res$change)
res$geneSymbol <- rownames(res)
resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
table(duplicated(resInP$geneSymbol))
resInP_list[['EC109.hs_METTL3KD_vs_WT']] <- resInP
#####condition_EC109_METTL14KD_vs_WT#####
count <- as.matrix(countGTF[, grepl('EC109_WT|EC109_METTL14KD', sample_meta$condition)])
table(is.na(count))
sample <- sample_meta[sample_meta$sample %in% colnames(count), ]
## SummarizedExperiment input
condition = c(rep('EC109_METTL14KD', 2), rep('EC109_WT', 2))
dds<-DESeqDataSetFromMatrix(countData=count, colData=sample,design=~condition)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 10
table(keep)
dds <- dds[keep,]
## Note on factor levels
dds$condition <- factor(dds$condition, levels = c("EC109_METTL14KD","EC109_WT"))
dds$condition <- relevel(dds$condition, ref = "EC109_WT")
dds <- DESeq(dds)
res <-  results(dds, contrast=c("condition",unique(condition))) # c("condition", rev(levels(genotype)))
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
table(res$change)
res$geneSymbol <- rownames(res)
resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
table(duplicated(resInP$geneSymbol))
resInP_list[['EC109_METTL14KD_vs_WT']] <- resInP
RMdeg_list[[i]] <- resInP_list
names(RMdeg_list[[i]] ) <- names(resInP_list)#name
contrast_input_list[[i]] <- NULL
condition_input_list[[i]] <- condition_list

############
n = 107
############
#GSE171497# with RIP and Input
############
print(RMPdeg_datasets[n,])# undone ! see http://localhost:8888/edit/Shell/UTH_RPMfunc/RPMdb-raw.sh
i = RMPdeg_datasets[n,'GSE']
cl = c('U87')#Glioblastoma
RMp = c('ALKBH5')
RM = 'm6A'
condition = c(rep('ALKBH5KD', 1), rep('WT', 1))#
RMP_allDat[grep('GSE171497', RMP_allDat$study_alias )
           ,c("study_alias", "sample_alias", "sample_title")]
# View(allDat[[i]])
files <- list.files(file.path(data_directory,i))
files
getwd()
setwd(data_directory)
setwd(i)
df_in <- read_excel(file.path(data_directory,i, 'GSE171227_ALKBH5-RNA-seq-FPKM.xls'),sheet = 2)
countGTF <- as.data.frame(df_in[, 2:4])
table(duplicated(df_in$gene_name))
countGTF <- aggregate(x = countGTF[,1:2],   #此时exprSet的第三列开始是表达矩阵内容
                      by = list(symbol = countGTF$gene_name),   #按照相同symbol分组，在组内计算
                      FUN = max) %>%   #原文中是计算最大值（max），也可以计算平均值（mean）或者中位数（median）
  column_to_rownames(var = 'symbol')
dim(countGTF) # 51242    12
res <- countGTF
fold_change <- countGTF[,1]/countGTF[,2]#cpm[,1]/cpm[,2]
res <- countGTF
res$logFC <- log2(fold_change)
res$FoldChange <- fold_change
res$change <- with(res, ifelse(abs(res$logFC) > log2(1.5),
                               ifelse(res$logFC > log2(1.5),'UP','DOWN'),'NOT'))
table(res$change)
res$geneSymbol <- rownames(res)
resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
table(resInP$change)
table(resInP$FoldChange %in% c(0, Inf))
table(resInP$logFC %in% c(-Inf, Inf))
# resInP <- resInP[! resInP$FoldChange %in% c(0, Inf), ]
name = "U87_ALKBH5KD_vs_WT"
RMdeg_list[[i]] <- list( resInP)
names(RMdeg_list[[i]] ) <- name
condition_list <- list(condition); names(condition_list) <- 'ALKBH5KD'
contrast_input_list[[i]] <- NULL#countGTF_list
condition_input_list[[i]] <- condition_list
############
n = 108
############
#GSE240674# with RIP and Input
############
print(RMPdeg_datasets[n,])#
i = RMPdeg_datasets[n,'GSE']
cl = c('THP1')#AML
RMp = c('METTL3')
RM = 'm6A'
condition_list = list(THP1.N_METTL3KO = c(rep('METTL3KO', 2), rep('WT', 2)),
                      THP1.NL_METTL3KO = c(rep('METTL3KO', 2), rep('WT', 2)),
                      THP1.TL_METTL3KO = c(rep('METTL3KO', 2), rep('WT', 2))
)
RMP_allDat[grep('GSE240674', RMP_allDat$study_alias )
           ,c("study_alias", "sample_alias", "sample_title")]
# View(allDat[[i]])
RMP_allDat[grepl('GSE240674', RMP_allDat$study_alias ) &  (!grepl('G9a', RMP_allDat$sample_title))
           ,c("study_alias", "sample_alias", "sample_title")]
skip = RMP_allDat[grepl('GSE240674', RMP_allDat$study_alias ) & grepl('G9a', RMP_allDat$sample_title)
                  ,c("study_alias", "sample_alias", "sample_title")]
paste0(skip$sample_alias, collapse = "|")

skip <- RMP_allDat[grepl('GSE240674', RMP_allDat$study_alias ) & (!grepl('G9a', RMP_allDat$sample_title))
                   &  grepl('Input', RMP_allDat$sample_title)
                   , "sample_alias"]
paste0(skip$sample_alias, collapse = " ")
skip <- RMP_allDat[grepl('GSE240674', RMP_allDat$study_alias ) & (!grepl('G9a', RMP_allDat$sample_title))
                   &  grepl('m6A-IP', RMP_allDat$sample_title)
                   , "sample_alias"]
paste0(skip$sample_alias, collapse = " ")
files <- list.files(file.path(data_directory,i))
files
getwd()
setwd(data_directory)
setwd(i)
ENA <- as.data.frame(RMP_allDat[grepl('GSE240674', RMP_allDat$study_alias ) & (!grepl('G9a', RMP_allDat$sample_title))
                                ,c("study_alias", "sample_alias", "sample_title")])
ENA <- ENA[order(ENA$sample_alias),]
ENA

countfiles <- list.files('./metaPlotR/')
countfiles <- countfiles[grep('.txt', countfiles)]
countfiles <- countfiles[! grepl('summary', countfiles)] #mouse embryonic fibroblasts (MEFs)
countfiles
###################################
#  Cell number-normalized (CNN)   #
#      MeRIP INPUT analysis       #
###################################
resInP_list <- vector("list", length(condition_list))
names(resInP_list) <- names(condition_list)
countGTF_list <- vector("list", length(condition_list))
names(countGTF_list) <- names(condition_list)
for(c in names(condition_list)){
  d = str_split(c, '[.]', simplify = T)[,2]
  countGTF <- read.table( paste0("metaPlotR/Input_counts_", d, ".txt"),header=T)
  dim(countGTF)
  colnames(countGTF)
  rownames(countGTF) <- countGTF$Geneid
  countGTF <- countGTF[, !colnames(countGTF) %in% c('Geneid','Chr','Start','End','Strand','Length')]
  dim(countGTF) # 28265     4
  condition <- condition_list[[c]]
  ## Sample information
  sample=data.frame(
    sample=colnames(countGTF),
    condition=condition)
  ## SummarizedExperiment input
  dds<-DESeqDataSetFromMatrix(countData=countGTF,colData=sample,design=~condition)
  ## Pre-filtering
  keep <- rowSums(counts(dds)) >= 10
  table(keep)
  dds <- dds[keep,]
  ## Note on factor levels
  dds$condition <- factor(dds$condition, levels = unique(condition))
  dds$condition <- relevel(dds$condition, ref = "WT")
  # Differential expression analysis-------------------------------------------
  dds <- DESeq(dds)
  res <-  results(dds, contrast=c("condition",unique(condition)))
  res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
  length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
  res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                                 ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
  table(res$change)
  table(is.na(res$change))
  # resOrdered <- as.data.frame(res[order(res$padj), ])
  resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
  resInP$geneSymbol <- rownames(resInP)
  table(resInP$change)
  resInP_list[[c]] <- resInP
  countGTF_list[[c]] <- countGTF
}
name = paste0(names(resInP_list), ' vs WT')
# name = paste0(cl, '-', name)
RMdeg_list[[i]] <- resInP_list
names(RMdeg_list[[i]] ) <- paste0(gsub(' ','_',(name)))
contrast_input_list[[i]] <- lapply(countGTF_list, extract_names)
condition_input_list[[i]] <- condition_list

###################################
#  Cell number-normalized (CNN)   #
#   MeRIP RIP analysis - Part1    #
###################################
resMeR_list <- vector("list", length(condition_list))
names(resMeR_list) <- names(condition_list)
countMatrix_list <- vector("list", length(condition_list))
names(countMatrix_list) <- names(condition_list)
resPeak_list <- vector("list",length(condition_list))
names(resPeak_list) <- names(condition_list)
for(c in names(condition_list)){
  d= str_split(c, '[.]', simplify = T)[,2]
  merge_peak <- read.delim(paste0("metaPlotR/", d, "_newPeakFile.txt"),header=F)
  countMatrix <- read.table(paste0("metaPlotR/m6A_counts_", d, ".txt"),header=T)
  table(countMatrix$Geneid == merge_peak[-1,1])
  table(duplicated(merge_peak[-1,1]))
  dim(countMatrix)
  colnames(countMatrix)
  rownames(countMatrix) <- countMatrix$Geneid
  countMatrix <- countMatrix[, !colnames(countMatrix) %in% c('Geneid','Chr','Start','End','Strand','Length')]
  # 1. 做好分组因子即可
  condition <- condition_list[[c]]
  groupList<-factor(condition)
  colData <- data.frame(row.names=colnames(countMatrix), groupList=groupList)
  # 2. Exploring Method 1: DESeq2-------------------------------------------
  # 构建 DESeqDataSet 对象
  dds <- DESeqDataSetFromMatrix(countData = countMatrix,
                                colData = colData,
                                design = ~ groupList)
  ## Pre-filtering
  keep <- rowSums(counts(dds)) >= 1
  table(keep)
  dds <- dds[keep,]
  # 接着计算每个样本的 library size 用于后续的标准化。library size 等于落在该样本所有peaks上的reads的总数，即 counts 矩阵中每列的加和。如果bFullLibrarySize设为TRUE，则会使用library的总reads数（根据原BAM/BED文件统计获得）。然后使用estimateDispersions估计统计分布，需要将参数fitType设为’local’
  dds<-estimateSizeFactors(dds)
  dds<-estimateDispersions(dds,fitType="local")
  dds <- nbinomWaldTest(dds)
  res <- results(dds)
  ## 提取你想要的差异分析结果，我们这里是treated组对untreated组进行比较
  res <-  results(dds, contrast=c("groupList", unique(condition)))
  res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
  length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
  res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                                 ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
  res_DESeq2 = as.data.frame(res)
  table(res_DESeq2$change)
  res_DESeq2 = res_DESeq2[res_DESeq2$change %in% c('UP', 'DOWN'), ]
  # 2. Exploring Method 2: limma-------------------------------------------
  library(edgeR)
  d0 <- DGEList(counts = countMatrix)
  dim(d0)
  # Filtering low count peaks
  cutoff <- 1
  drop <- which(apply(cpm(d0), 1, max) < cutoff)
  print(length(drop))
  if (length(drop) >0) {
    d <- d0[-drop,]
  } else{ d <- d0}
  dim(d) #number of peaks left: 8501
  d <- calcNormFactors(d, method = c('TMM'))
  colnames(countMatrix)
  genotype <- groupList
  table(genotype)
  mm <- model.matrix(~0 + genotype)#
  v <- voom(d, mm, plot=F, normalize="quantile")#voom from limma:Transform count data to log2-counts per million (logtpm)
  fit <- lmFit(v, mm)
  head(coef(fit))
  contr <- makeContrasts(paste0(colnames(coef(fit)), collapse  = ' - '), levels = colnames(coef(fit))) #Change if other comparison
  tmp <- contrasts.fit(fit, contr)
  tmp <- eBayes(tmp)
  top.table <- topTable(tmp, sort.by = "P", n = Inf)
  head(top.table, 20)
  length(which(top.table$adj.P.Val < 0.05 & abs(top.table$logFC) > log2(1.5))) #how many are significant?
  res_limma <- as.data.frame(top.table)
  res_limma$change <- with(res_limma, ifelse(res_limma$adj.P.Val < 0.05 & abs(res_limma$logFC) > log2(1.5),
                                             ifelse(res_limma$logFC > log2(1.5),'UP','DOWN'),'NOT'))
  table(res_limma$change)
  res_limma = res_limma[res_limma$change %in% c('UP', 'DOWN'), ]
  resMeR_list[[c]] <- list(res_DESeq2, res_limma)
  countMatrix_list[[c]] <- countMatrix
  resPeak_list[[c]] <- merge_peak
  name = res@elementMetadata$description[2]
  # name = paste0(cl, '_',str_split(name, "groupList ", simplify = T)[,2])
  name = str_split(name, "groupList ", simplify = T)[,2]
  name = paste0(gsub(' ','_',(name)))
  names(resMeR_list[[c]] ) <- paste0(name, '_', c('DESeq2', 'limma'))
}
RMdem_list[[i]] <- as.list(unlist(resMeR_list, recursive = FALSE))
name <- gsub("_METTL3KO.", '_', names(RMdem_list[[i]] ))#str_split(names(RMdem_list[[i]] ), "[.]", simplify = T)[,2]
names(RMdem_list[[i]] ) <- name
RMdep_list[[i]] <- resPeak_list
name <- unique(gsub('(_limma|_DESeq2)', '', name))
names(RMdep_list[[i]] ) <- name
contrast_IP_list[[i]] <- lapply(countMatrix_list, extract_names)
condition_IP_list[[i]] <- condition_list

############
n = 109
############
#GSE174492# without case RIP, with RNA
############
print(RMPdeg_datasets[n,])#
i = RMPdeg_datasets[n,'GSE']
cl = c('MHCC97H and SNU449')#Liver cancer
RMp = c('METTL1')
RM = 'm7G'
# condition = c(rep('METTL1KD', 2), rep('WT', 2))#
condition = c(rep('METTL1OE',1), rep('WT', 1))#
RMP_allDat[grep('GSE174492', RMP_allDat$study_alias )
           ,c("study_alias", "sample_alias", "sample_title")]
# View(allDat[[i]])
files <- list.files(file.path(data_directory,i))
files
getwd()
setwd(data_directory)
setwd(i)
df_in <- read.table(file.path(data_directory,i, 'GSE174492_RPKM-TRAC-Seq-SNU449.txt'), header = T)
# DEG <- fread(file.path(data_directory,i, 'GSE174492_RPKM_TRAC-Seq_MHCC97H.txt'))# file is corrupt ?
# df_in <- read.table(file.path(data_directory,i, 'GSE174492_FPKM_polyribosome_associated_mRNA-seq.txt'), header = T)
countGTF <- as.data.frame(df_in[, c(2, 4)])
table(duplicated(df_in$tRNA.name))
res <- countGTF
fold_change <- countGTF[,1]/countGTF[,2]#cpm[,1]/cpm[,2]
res <- countGTF
res$logFC <- log2(fold_change)
res$FoldChange <- fold_change
res$change <- with(res, ifelse(abs(res$logFC) > log2(1.5),
                               ifelse(res$logFC > log2(1.5),'UP','DOWN'),'NOT'))
table(res$change)
res$geneSymbol <- df_in$tRNA.name
resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
table(resInP$change)
table(resInP$FoldChange %in% c(0, Inf))
table(resInP$logFC %in% c(-Inf, Inf))
# resInP <- resInP[! resInP$FoldChange %in% c(0, Inf), ]
name = "SNU449_METTL1OE_vs_WT"
RMdeg_list[[i]] <- list( resInP)
names(RMdeg_list[[i]] ) <- name
condition_list <- list(condition); names(condition_list) <- 'METTL1OE'
contrast_input_list[[i]] <- NULL#countGTF_list
condition_input_list[[i]] <- condition_list

countMatrix <- as.data.frame(df_in[, c(3, 5)])
table(duplicated(df_in$tRNA.name))
resMeR <- countMatrix
fold_change <- countMatrix[,1]/countMatrix[,2]#cpm[,1]/cpm[,2]
resMeR <- countMatrix
resMeR$logFC <- log2(fold_change)
resMeR$FoldChange <- fold_change
resMeR$change <- with(resMeR, ifelse(abs(resMeR$logFC) > log2(1.5),
                               ifelse(resMeR$logFC > log2(1.5),'UP','DOWN'),'NOT'))
table(resMeR$change)
resMeR$geneSymbol <- df_in$tRNA.name
res_limma <- resMeR[resMeR$change %in% c('UP','DOWN'),]
RMdem_list[[i]] <- list(res_limma)
name = "SNU449_METTL1OE_vs_WT"
RMdem_list[[i]] <- list( resMeR)
names(RMdem_list[[i]] ) <- name
RMdep_list[[i]] <- NULL #list( res_limma ) no peak analysis
# names(RMdep_list[[i]] ) <- name
condition_list <- list(condition); names(condition_list) <- 'METTL1OE'
contrast_IP_list[[i]] <- NULL
condition_IP_list[[i]] <- condition_list

############
n = 110
############
#GSE71096# without case RIP, with RNA
############
print(RMPdeg_datasets[n,])#
i = RMPdeg_datasets[n,'GSE']
cl = c('Hela')#
RMp = c('YTHDC1','SRSF10')
RM = 'm6A'
# condition = c(rep('SRSF10OE', 2), rep('WT', 2))#
condition = c(rep('YTHDC1OE', 2), rep('WT', 3))#
RMP_allDat[grep('GSE71096', RMP_allDat$study_alias )
           ,c("study_alias", "sample_alias", "sample_title")]
# View(allDat[[i]])
files <- list.files(file.path(data_directory,i))
files
getwd()
setwd(data_directory)
setwd(i)
df_in <- read.table(file.path(data_directory,i, 'GSE71095_reads-count.txt.gz'), header = T)
library(data.table)
colnames(df_in)
countGTF=df_in[, grepl('CTRL|YTHDC1', colnames(df_in))]
# id_mapped <- read.table("/data2/rluo4/hg38/gencode.v44lift37.annotation.gene.probeMap",header = T)
id_mapped <- read.table("/data2/rluo4/hg38/gencode.v36.annotation.gene.probeMap",header = T)
countGTF$id <- str_split(rownames(df_in), '[.]', simplify = T)[,1]
# id_mapped <- unique(id_mapped) #%>% dplyr::select(id,gene) # 去重
id_mapped$id <- str_split(id_mapped$id, '[.]', simplify = T)[,1]
id_mapped <- id_mapped[which(id_mapped$id %in% countGTF$id ),]
dim(id_mapped) #57930     6
table(duplicated(id_mapped$gene))
countGTF <- left_join(id_mapped,countGTF,by="id")
countGTF <- aggregate(x = countGTF[,(ncol(id_mapped)+1):ncol(countGTF)],   #此时exprSet的第三列开始是表达矩阵内容
                      by = list(symbol = countGTF$gene),   #按照相同symbol分组，在组内计算
                      FUN = max) %>%   #原文中是计算最大值（max），也可以计算平均值（mean）或者中位数（median）
  column_to_rownames(var = 'symbol')
dim(countGTF) # 51242    12
geneSymbol <- rownames(countGTF)
## Gene name
dim(countGTF) #25343    12
## Sample name
colnames(countGTF)
## Sample information
sample=data.frame(
  sample=colnames(countGTF),
  condition=c('WT', 'YTHDC1OE', 'WT', 'WT', 'YTHDC1OE'))
## SummarizedExperiment input
dds<-DESeqDataSetFromMatrix(countData=countGTF,colData=sample,design=~condition)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 10
table(keep)
dds <- dds[keep,]
## Note on factor levels
dds$condition <- factor(dds$condition, levels = unique(condition))
dds$condition <- relevel(dds$condition, ref = "WT")
# Differential expression analysis-------------------------------------------
dds <- DESeq(dds)
res <-  results(dds, contrast=c("condition",unique(condition)))
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
table(res$change)
table(is.na(res$change))
# resOrdered <- as.data.frame(res[order(res$padj), ])
resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
resInP$geneSymbol <- rownames(resInP)
c <- unique(condition[!condition %in% 'WT'])
name = res@elementMetadata$description[2]
name = paste0(cl, '_',str_split(name, "condition ", simplify = T)[,2])
RMdeg_list[[i]] <- list( resInP)
names(RMdeg_list[[i]] ) <- paste0(gsub(' ','_',(name)))
condition_list <- list(condition); names(condition_list) <- c
contrast_input_list[[i]] <- NULL
condition_input_list[[i]] <- condition_list

############
n = 111
############
#GSE207643# with RIP and Input
############
print(RMPdeg_datasets[n,])#
i = RMPdeg_datasets[n,'GSE']
cl = c('K562', 'Hela')#CMML
RMp = c('NSUN2')
RM = 'm5C'
condition_list = list(Hela_NSUN2KD = c(rep('NSUN2KD', 2), rep('WT', 2)), 
                      K562_NSUN2KD = c(rep('NSUN2KD', 2), rep('WT', 2))#
)
RMP_allDat[grep('GSE207643', RMP_allDat$study_alias )
           ,c("study_alias", "sample_alias", "sample_title")]
# View(allDat[[i]])
ENA <- as.data.frame(RMP_allDat[grepl('GSE207643', RMP_allDat$study_alias ) & (!grepl('PAR-CLIP|Biopsy', RMP_allDat$sample_title))
                                ,c("study_alias", "sample_alias", "sample_title")])
ENA <- ENA[order(ENA$sample_alias),]
ENA
files <- list.files(file.path(data_directory,i))
files
getwd()
setwd(data_directory)
setwd(i)
###################################
#      RNA/Ribo-seq analysis      #
###################################
# 1) RNA-seq
count_files <- files[grepl('_htseqcount.tsv.gz', files)]
###################################
# Read and process each count file
count <- lapply(count_files, function(x) {
  # Read data using fread for faster reading
  data <- fread(file.path(data_directory, i, x))
  # Extract required columns
  df_new <- data[, c(1, 2)]
  # Set column names same as file names
  colnames(df_new) <- c('geneSymbol',gsub('_htseqcount.tsv.gz', '', x))
  return(df_new)
})
# Merge dataframes in countGTF list by 'gene_id'
df_in <- Reduce(function(x, y) merge(x, y, by = 'geneSymbol', all = TRUE), count)
colSums(df_in[,-1])
table(colSums(df_in[,-1])=="1e+06")
## Gene name
library(data.table)
countGTF=df_in[,2:ncol(df_in)]
# id_mapped <- read.table("/data2/rluo4/hg38/gencode.v44lift37.annotation.gene.probeMap",header = T)
id_mapped <- read.table("/data2/rluo4/hg38/gencode.v36.annotation.gene.probeMap",header = T)
countGTF$id <- str_split(df_in$geneSymbol, '[.]', simplify = T)[,1]
# id_mapped <- unique(id_mapped) #%>% dplyr::select(id,gene) # 去重
id_mapped$id <- str_split(id_mapped$id, '[.]', simplify = T)[,1]
id_mapped <- id_mapped[which(id_mapped$id %in% countGTF$id ),]
dim(id_mapped) #57930     6
table(duplicated(id_mapped$gene))
countGTF <- left_join(id_mapped,countGTF,by="id")
countGTF <- aggregate(x = countGTF[,(ncol(id_mapped)+1):ncol(countGTF)],   #此时exprSet的第三列开始是表达矩阵内容
                      by = list(symbol = countGTF$gene),   #按照相同symbol分组，在组内计算
                      FUN = max) %>%   #原文中是计算最大值（max），也可以计算平均值（mean）或者中位数（median）
  column_to_rownames(var = 'symbol')
dim(countGTF) # 51242    12
geneSymbol <- rownames(countGTF)
# countGTF_1=df_in[!df_in$geneSymbol %in% id_mapped$id, 2:ncol(df_in)]
# geneSymbol_1 <- df_in[!df_in$geneSymbol %in% id_mapped$id, ]
# table(geneSymbol_1$geneSymbol %in% geneSymbol)
# intersect(geneSymbol_1$geneSymbol, geneSymbol)
###################################
#   more than 1 condition here   #
###################################
resInP_list <- list()
##Gene count matrix
# countGTF=as.matrix(ceiling(df_in[, -1]))
## Gene name
table(duplicated(df_in$geneSymbol))
dim(countGTF) #25343    12
## Sample name
colnames(countGTF)
## Sample information
sample_meta=data.frame(
  sample=colnames(countGTF),
  condition=c(rep("Hela_WT",2),rep("Hela_NSUN2KD",2), rep("Hela_SRSF2KD",2),rep("K562_WT",2), rep("K562_NSUN2KD",2),rep("K562_SRSF2MU",2)))
#####condition_HeLa_siNSUN2_vs_WT#####
count <- as.matrix(countGTF[, grepl('Hela_WT|Hela_NSUN2KD', sample_meta$condition)])
table(is.na(count))
sample <- sample_meta[sample_meta$sample %in% colnames(count), ]
## SummarizedExperiment input
condition = c(rep('Hela_NSUN2KD', 2), rep('Hela_WT', 2))
dds<-DESeqDataSetFromMatrix(countData=count, colData=sample,design=~condition)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 10
table(keep)
dds <- dds[keep,]
## Note on factor levels
dds$condition <- factor(dds$condition, levels = c("Hela_NSUN2KD","Hela_WT"))
dds$condition <- relevel(dds$condition, ref = "Hela_WT")
dds <- DESeq(dds)
res <-  results(dds, contrast=c("condition",unique(condition))) # c("condition", rev(levels(genotype)))
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
table(res$change)
res$geneSymbol <- geneSymbol[keep]
resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
table(duplicated(resInP$geneSymbol))
resInP_list[['Hela_NSUN2KD_vs_WT']] <- resInP
#####condition_K562_shNSUN2_vs_WT#####
count <- as.matrix(countGTF[, grepl('K562_WT|K562_NSUN2KD', sample_meta$condition)])
table(is.na(count))
sample <- sample_meta[sample_meta$sample %in% colnames(count), ]
## SummarizedExperiment input
condition = c(rep('K562_NSUN2KD', 2), rep('K562_WT', 2))
dds<-DESeqDataSetFromMatrix(countData=count, colData=sample,design=~condition)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 10
table(keep)
dds <- dds[keep,]
## Note on factor levels
dds$condition <- factor(dds$condition, levels = c("K562_NSUN2KD","K562_WT"))
dds$condition <- relevel(dds$condition, ref = "K562_WT")
dds <- DESeq(dds)
res <-  results(dds, contrast=c("condition",unique(condition))) # c("condition", rev(levels(genotype)))
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
table(res$change)
res$geneSymbol <- geneSymbol[keep]
resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
table(duplicated(resInP$geneSymbol))
resInP_list[['K562_NSUN2KD_vs_WT']] <- resInP

###################################
#  Cell number-normalized (CNN)   #
#      MeRIP INPUT analysis       #
###################################
cl = 'Hela'
condition = c(rep('NSUN2KD', 2), rep('WT', 2))#
c <- unique(condition[! condition %in% 'WT'])
countGTF <- read.table( paste0("metaPlotR/Input_counts_", c, ".txt"),header=T)
rownames(countGTF) <- countGTF$Geneid
countGTF <- countGTF[, !colnames(countGTF) %in% c('Geneid','Chr','Start','End','Strand','Length')]
dim(countGTF) #28265     4
## Sample name
## Sample information
sample=data.frame(
  sample=colnames(countGTF),
  condition=condition)
## SummarizedExperiment input
dds<-DESeqDataSetFromMatrix(countData=countGTF,colData=sample,design=~condition)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 10
table(keep)
dds <- dds[keep,]
## Note on factor levels
dds$condition <- factor(dds$condition, levels = unique(condition))
dds$condition <- relevel(dds$condition, ref = "WT")
# Differential expression analysis-------------------------------------------
dds <- DESeq(dds)
res <-  results(dds, contrast=c("condition",unique(condition)))
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
table(res$change)
table(is.na(res$change))
# resOrdered <- as.data.frame(res[order(res$padj), ])
resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
resInP$geneSymbol <- rownames(resInP)
countGTF_list <- list(colnames(countGTF)); names(countGTF_list) <- c
resInP_list[['Hela_NSUN2KD_vs_WT']] <- resInP
RMdeg_list[[i]] <- resInP_list
# name <- paste0(cl[1],'_', names(resInP_list))
names(RMdeg_list[[i]] ) <- names(resInP_list)#name
contrast_input_list[[i]] <- countGTF_list
names(contrast_input_list[[i]]) <- names(condition_list)[1]
condition_input_list[[i]] <- condition_list

###################################
#  Cell number-normalized (CNN)   #
#   MeRIP RIP analysis - Part1    #
###################################
## Sample information
merge_peak <- read.delim(paste0("metaPlotR/", c, "_newPeakFile.txt"),header=F)
countMatrix <- read.table(paste0("metaPlotR/m5C_counts_", c, ".txt"),header=T)
table(countMatrix$Geneid == merge_peak[-1,1])
table(duplicated(merge_peak$V1[-1]))
dim(countMatrix)
colnames(countMatrix)
rownames(countMatrix) <- countMatrix$Geneid
countMatrix <- countMatrix[, !colnames(countMatrix) %in% c('Geneid','Chr','Start','End','Strand','Length')]
# 1. 做好分组因子即可
groupList<-factor(condition)
colData <- data.frame(row.names=colnames(countMatrix), groupList=groupList)
# 2. Exploring Method 1: DESeq2-------------------------------------------
# 构建 DESeqDataSet 对象
dds <- DESeqDataSetFromMatrix(countData = countMatrix,
                              colData = colData,
                              design = ~ groupList)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 1
table(keep)
dds <- dds[keep,]
# 接着计算每个样本的 library size 用于后续的标准化。library size 等于落在该样本所有peaks上的reads的总数，即 counts 矩阵中每列的加和。如果bFullLibrarySize设为TRUE，则会使用library的总reads数（根据原BAM/BED文件统计获得）。然后使用estimateDispersions估计统计分布，需要将参数fitType设为’local’
dds<-estimateSizeFactors(dds)
dds<-estimateDispersions(dds,fitType="local")
dds <- nbinomWaldTest(dds)
res <- results(dds)
## 提取你想要的差异分析结果，我们这里是treated组对untreated组进行比较
res <-  results(dds, contrast=c("groupList", unique(condition)))
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
res_DESeq2 = as.data.frame(res)
table(res_DESeq2$change)#
res_DESeq2 = res_DESeq2[res_DESeq2$change %in% c('UP', 'DOWN'), ]
# 2. Exploring Method 2: limma-------------------------------------------
library(edgeR)
d0 <- DGEList(counts = countMatrix)
dim(d0)
# Filtering low count peaks
cutoff <- 1           
drop <- which(apply(cpm(d0), 1, max) < cutoff)
print(length(drop))
if (length(drop) >0) {
  d <- d0[-drop,] 
} else{ d <- d0}
dim(d) #number of peaks left: 8501      
d <- calcNormFactors(d, method = c('TMM'))
colnames(countMatrix)
genotype <- groupList
table(genotype)
mm <- model.matrix(~0 + genotype)#
v <- voom(d, mm, plot=F, normalize="quantile")#voom from limma:Transform count data to log2-counts per million (logtpm)
fit <- lmFit(v, mm)
head(coef(fit))
contr <- makeContrasts(paste0(colnames(coef(fit)), collapse  = ' - '), levels = colnames(coef(fit))) #Change if other comparison
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)
length(which(top.table$adj.P.Val < 0.05 & abs(top.table$logFC) > log2(1.5))) #how many are significant?
res_limma <- as.data.frame(top.table)
res_limma$change <- with(res_limma, ifelse(res_limma$adj.P.Val < 0.05 & abs(res_limma$logFC) > log2(1.5),
                                           ifelse(res_limma$logFC > log2(1.5),'UP','DOWN'),'NOT'))
table(res_limma$change)#
res_limma = res_limma[res_limma$change %in% c('UP', 'DOWN'), ]

RMdem_list[[i]] <- list(res_DESeq2, res_limma)
name = res@elementMetadata$description[2]
name = paste0(cl, '_',str_split(name, "groupList ", simplify = T)[,2])
name= paste0(gsub(' ','_',(name)))
names(RMdem_list[[i]] ) <- paste0(name, '_', c('DESeq2', 'limma'))
RMdep_list[[i]] <- list( merge_peak )
names(RMdep_list[[i]] ) <- name
countMatrix_list <- list(colnames(countMatrix)); names(countMatrix_list) <-  paste0(cl, '_', c)
condition_list <- list(condition); names(condition_list) <-  paste0(cl, '_', c)
contrast_IP_list[[i]] <- countMatrix_list
condition_IP_list[[i]] <- condition_list
# Define a function to transform the filename
transform_filename <- function(filename) {
  # Use regular expressions to replace patterns
  # Replace .siNSUN2. or .siCTRL. with -siNSUN2- or -siCTRL-
  filename <- gsub("\\.si([A-Za-z0-9]+)\\.", "-si\\1-", filename)
  # Replace remaining dots with hyphens
  filename <- gsub("\\.", "-", filename)
  return(filename)
}
contrast_IP_list[['GSE207643']]$Hela_NSUN2KD <- gsub('-bam', '.bam', transform_filename(contrast_IP_list[['GSE207643']]$Hela_NSUN2KD))
contrast_input_list[['GSE207643']]$Hela_NSUN2KD <- gsub('-bam', '.bam', transform_filename(contrast_input_list[['GSE207643']]$Hela_NSUN2KD))

############
n = 112
############
#GSE174374# with RIP and Input
############
print(RMPdeg_datasets[n,])#
i = RMPdeg_datasets[n,'GSE']
cl = c('A549')#
RMp = c('NSUN2') 
RM = 'm5C'

condition_list = list(A549.RSV_NSUN2KO = c(rep('NSUN2KO', 2), rep('WT', 2)),#
                      A549.VSV_NSUN2KO = c(rep('NSUN2KO', 2), rep('WT', 2)))
RMP_allDat[grep('GSE174374', RMP_allDat$study_alias ) 
           ,c("study_alias", "sample_alias", "sample_title")]
# View(allDat[[i]])
RMP_allDat[grepl('GSE174374', RMP_allDat$study_alias ) &  grepl('m5C_', RMP_allDat$sample_title)
           ,c("study_alias", "sample_alias", "sample_title")]
skip = RMP_allDat[grepl('GSE174374', RMP_allDat$study_alias ) & grepl('m5C_', RMP_allDat$sample_title)
                  ,c("study_alias", "sample_alias", "sample_title")]
paste0(skip$sample_alias, collapse = "|")

skip <- RMP_allDat[grepl('GSE174374', RMP_allDat$study_alias ) & grepl('m5C_', RMP_allDat$sample_title)
                   &  grepl('Input', RMP_allDat$sample_title)
                   , "sample_alias"]
paste0(skip$sample_alias, collapse = " ")
skip <- RMP_allDat[grepl('GSE174374', RMP_allDat$study_alias ) & grepl('m5C_', RMP_allDat$sample_title)
                   &  grepl('BS', RMP_allDat$sample_title)
                   , "sample_alias"]
paste0(skip$sample_alias, collapse = " ")
ENA <- as.data.frame(RMP_allDat[grepl('GSE174374', RMP_allDat$study_alias ) &  grepl('m5C_', RMP_allDat$sample_title)
                                ,c("study_alias", "sample_alias", "sample_title")])
ENA <- ENA[order(ENA$sample_alias),]
ENA
files <- list.files(file.path(data_directory,i))
files
getwd()
setwd(data_directory)
setwd(i)
###################################
#  Cell number-normalized (CNN)   #
#      MeRIP INPUT analysis       #
###################################
resInP_list <- vector("list", length(condition_list))
names(resInP_list) <- names(condition_list)
countGTF_list <- vector("list", length(condition_list))
names(countGTF_list) <- names(condition_list)
for(c in names(condition_list)){
  countGTF <- read.table( paste0("metaPlotR/Input_counts_", c, ".txt"),header=T)
  dim(countGTF)
  colnames(countGTF)
  rownames(countGTF) <- countGTF$Geneid
  countGTF <- countGTF[, !colnames(countGTF) %in% c('Geneid','Chr','Start','End','Strand','Length')]
  dim(countGTF) # 28265     4
  condition <- condition_list[[c]]
  ## Sample information
  sample=data.frame(
    sample=colnames(countGTF),
    condition=condition)
  ## SummarizedExperiment input
  dds<-DESeqDataSetFromMatrix(countData=countGTF,colData=sample,design=~condition)
  ## Pre-filtering
  keep <- rowSums(counts(dds)) >= 10
  table(keep)
  dds <- dds[keep,]
  ## Note on factor levels
  dds$condition <- factor(dds$condition, levels = unique(condition))
  dds$condition <- relevel(dds$condition, ref = "WT")
  # Differential expression analysis-------------------------------------------
  dds <- DESeq(dds)
  res <-  results(dds, contrast=c("condition",unique(condition)))
  res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
  length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
  res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                                 ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
  table(res$change)
  table(is.na(res$change))
  # resOrdered <- as.data.frame(res[order(res$padj), ])
  resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
  resInP$geneSymbol <- rownames(resInP)
  table(resInP$change)
  resInP_list[[c]] <- resInP
  countGTF_list[[c]] <- countGTF
}
name = paste0(names(resInP_list), ' vs WT')
# name = paste0(cl, '-', name)
RMdeg_list[[i]] <- resInP_list
names(RMdeg_list[[i]] ) <- paste0(gsub(' ','_',(name)))
contrast_input_list[[i]] <- lapply(countGTF_list, extract_names)
condition_input_list[[i]] <- condition_list

###################################
#  Cell number-normalized (CNN)   #
#   MeRIP RIP analysis - Part1    #
###################################
resMeR_list <- vector("list", length(condition_list))
names(resMeR_list) <- names(condition_list)
countMatrix_list <- vector("list", length(condition_list))
names(countMatrix_list) <- names(condition_list)
resPeak_list <- vector("list",length(condition_list))
names(resPeak_list) <- names(condition_list)
for(c in names(condition_list)){
  merge_peak <- read.delim(paste0("metaPlotR/", c, "_newPeakFile.txt"),header=F)
  countMatrix <- read.table(paste0("metaPlotR/m5C_counts_", c, ".txt"),header=T)
  table(countMatrix$Geneid == merge_peak[-1,1])
  table(duplicated(merge_peak[-1,1]))
  dim(countMatrix)
  colnames(countMatrix)
  rownames(countMatrix) <- countMatrix$Geneid
  countMatrix <- countMatrix[, !colnames(countMatrix) %in% c('Geneid','Chr','Start','End','Strand','Length')]
  # 1. 做好分组因子即可
  condition <- condition_list[[c]]
  groupList<-factor(condition)
  colData <- data.frame(row.names=colnames(countMatrix), groupList=groupList)
  # 2. Exploring Method 1: DESeq2-------------------------------------------
  # 构建 DESeqDataSet 对象
  dds <- DESeqDataSetFromMatrix(countData = countMatrix,
                                colData = colData,
                                design = ~ groupList)
  ## Pre-filtering
  keep <- rowSums(counts(dds)) >= 1
  table(keep)
  dds <- dds[keep,]
  # 接着计算每个样本的 library size 用于后续的标准化。library size 等于落在该样本所有peaks上的reads的总数，即 counts 矩阵中每列的加和。如果bFullLibrarySize设为TRUE，则会使用library的总reads数（根据原BAM/BED文件统计获得）。然后使用estimateDispersions估计统计分布，需要将参数fitType设为’local’
  dds<-estimateSizeFactors(dds)
  dds<-estimateDispersions(dds,fitType="mean")
  dds <- nbinomWaldTest(dds)
  res <- results(dds)
  ## 提取你想要的差异分析结果，我们这里是treated组对untreated组进行比较
  res <-  results(dds, contrast=c("groupList", unique(condition)))
  res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
  length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
  res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                                 ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
  res_DESeq2 = as.data.frame(res)
  table(res_DESeq2$change)
  res_DESeq2 = res_DESeq2[res_DESeq2$change %in% c('UP', 'DOWN'), ]
  # 2. Exploring Method 2: limma-------------------------------------------
  library(edgeR)
  d0 <- DGEList(counts = countMatrix)
  dim(d0)
  # Filtering low count peaks
  cutoff <- 1
  drop <- which(apply(cpm(d0), 1, max) < cutoff)
  print(length(drop))
  if (length(drop) >0) {
    d <- d0[-drop,]
  } else{ d <- d0}
  dim(d) #number of peaks left: 8501
  d <- calcNormFactors(d, method = c('TMM'))
  colnames(countMatrix)
  genotype <- groupList
  table(genotype)
  mm <- model.matrix(~0 + genotype)#
  v <- voom(d, mm, plot=F, normalize="quantile")#voom from limma:Transform count data to log2-counts per million (logtpm)
  fit <- lmFit(v, mm)
  head(coef(fit))
  contr <- makeContrasts(paste0(colnames(coef(fit)), collapse  = ' - '), levels = colnames(coef(fit))) #Change if other comparison
  tmp <- contrasts.fit(fit, contr)
  tmp <- eBayes(tmp)
  top.table <- topTable(tmp, sort.by = "P", n = Inf)
  head(top.table, 20)
  length(which(top.table$adj.P.Val < 0.05 & abs(top.table$logFC) > log2(1.5))) #how many are significant?
  res_limma <- as.data.frame(top.table)
  res_limma$change <- with(res_limma, ifelse(res_limma$adj.P.Val < 0.05 & abs(res_limma$logFC) > log2(1.5),
                                             ifelse(res_limma$logFC > log2(1.5),'UP','DOWN'),'NOT'))
  table(res_limma$change)
  res_limma = res_limma[res_limma$change %in% c('UP', 'DOWN'), ]
  resMeR_list[[c]] <- list(res_DESeq2, res_limma)
  countMatrix_list[[c]] <- countMatrix
  resPeak_list[[c]] <- merge_peak
  name = res@elementMetadata$description[2]
  # name = paste0(cl, '_',str_split(name, "groupList ", simplify = T)[,2])
  name = str_split(name, "groupList ", simplify = T)[,2]
  name = paste0(gsub(' ','_',(name)))
  names(resMeR_list[[c]] ) <- paste0(name, '_', c('DESeq2', 'limma'))
}
RMdem_list[[i]] <- as.list(unlist(resMeR_list, recursive = FALSE))
name <- gsub("_NSUN2KO.", '_', names(RMdem_list[[i]] ))#str_split(names(RMdem_list[[i]] ), "[.]", simplify = T)[,2]
names(RMdem_list[[i]] ) <- name
RMdep_list[[i]] <- resPeak_list
name <- unique(gsub('(_limma|_DESeq2)', '', name))
names(RMdep_list[[i]] ) <- name
contrast_IP_list[[i]] <- lapply(countMatrix_list, extract_names)
condition_IP_list[[i]] <- condition_list

############
n = 113
############
#GSE253795# with RIP and Input
############
print(RMPdeg_datasets[n,])#
i = RMPdeg_datasets[n,'GSE']
cl = c('A375')#
RMp = c('ALKBH5') 
RM = 'm6A'
condition = c(rep('ALKBH5KD', 3), rep('WT', 3))#
RMP_allDat[grep('GSE253795', RMP_allDat$study_alias ) 
           ,c("study_alias", "sample_alias", "sample_title")]
# View(allDat[[i]])

skip <- RMP_allDat[grepl('GSE253795', RMP_allDat$study_alias ) &  grepl('input', RMP_allDat$sample_title)
                   , "sample_alias"]
paste0(skip$sample_alias, collapse = " ")
skip <- RMP_allDat[grepl('GSE253795', RMP_allDat$study_alias ) &  grepl('IP', RMP_allDat$sample_title)
                   , "sample_alias"]
paste0(skip$sample_alias, collapse = " ")
ENA <- as.data.frame(RMP_allDat[grepl('GSE253795', RMP_allDat$study_alias )
                                ,c("study_alias", "sample_alias", "sample_title")])
ENA <- ENA[order(ENA$sample_alias),]
ENA
files <- list.files(file.path(data_directory,i))
files
getwd()
setwd(data_directory)
setwd(i)
c <- unique(condition[!condition %in% 'WT'])
###################################
#  Cell number-normalized (CNN)   #
#      MeRIP INPUT analysis       #
###################################
countGTF <- read.table( paste0("metaPlotR/Input_counts_", c, ".txt"),header=T)
rownames(countGTF) <- countGTF$Geneid
countGTF <- countGTF[, !colnames(countGTF) %in% c('Geneid','Chr','Start','End','Strand','Length')]
dim(countGTF) #28265     4
## Sample name
## Sample information
sample=data.frame(
  sample=colnames(countGTF),
  condition=condition)
## SummarizedExperiment input
dds<-DESeqDataSetFromMatrix(countData=countGTF,colData=sample,design=~condition)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 10
table(keep)
dds <- dds[keep,]
## Note on factor levels
dds$condition <- factor(dds$condition, levels = unique(condition))
dds$condition <- relevel(dds$condition, ref = "WT")
# Differential expression analysis-------------------------------------------
dds <- DESeq(dds)

res <-  results(dds, contrast=c("condition",unique(condition)))
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
table(res$change)
table(is.na(res$change))
# resOrdered <- as.data.frame(res[order(res$padj), ])
resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
resInP$geneSymbol <- rownames(resInP)
name = res@elementMetadata$description[2]
name = paste0(cl, '_',str_split(name, "condition ", simplify = T)[,2])
RMdeg_list[[i]] <- list( resInP)
names(RMdeg_list[[i]] ) <- paste0(gsub(' ','_',(name)))
countGTF_list <- list(colnames(countGTF)); names(countGTF_list) <- c
condition_list <- list(condition); names(condition_list) <- c
contrast_input_list[[i]] <- countGTF_list
condition_input_list[[i]] <- condition_list

###################################
#  Cell number-normalized (CNN)   #
#   MeRIP RIP analysis - Part1    #
###################################
## Sample information
merge_peak <- read.delim(paste0("metaPlotR/", c, "_newPeakFile.txt"),header=F)
countMatrix <- read.table(paste0("metaPlotR/m6A_counts_", c, ".txt"),header=T)
table(countMatrix$Geneid == merge_peak[-1,1])
table(duplicated(merge_peak$V1[-1]))
dim(countMatrix)
colnames(countMatrix)
rownames(countMatrix) <- countMatrix$Geneid
countMatrix <- countMatrix[, !colnames(countMatrix) %in% c('Geneid','Chr','Start','End','Strand','Length')]
# 1. 做好分组因子即可
groupList<-factor(condition)
colData <- data.frame(row.names=colnames(countMatrix), groupList=groupList)
# 2. Exploring Method 1: DESeq2-------------------------------------------
# 构建 DESeqDataSet 对象
dds <- DESeqDataSetFromMatrix(countData = countMatrix,
                              colData = colData,
                              design = ~ groupList)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 1
table(keep)
dds <- dds[keep,]
# 接着计算每个样本的 library size 用于后续的标准化。library size 等于落在该样本所有peaks上的reads的总数，即 counts 矩阵中每列的加和。如果bFullLibrarySize设为TRUE，则会使用library的总reads数（根据原BAM/BED文件统计获得）。然后使用estimateDispersions估计统计分布，需要将参数fitType设为’local’
dds<-estimateSizeFactors(dds)
dds<-estimateDispersions(dds,fitType="local")
dds <- nbinomWaldTest(dds)
res <- results(dds)
## 提取你想要的差异分析结果，我们这里是treated组对untreated组进行比较
res <-  results(dds, contrast=c("groupList", unique(condition)))
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
res_DESeq2 = as.data.frame(res)
table(res_DESeq2$change)#Not only
res_DESeq2 = res_DESeq2[res_DESeq2$change %in% c('UP', 'DOWN'), ]
# 2. Exploring Method 2: limma-------------------------------------------
library(edgeR)
d0 <- DGEList(counts = countMatrix)
dim(d0)
# Filtering low count peaks
cutoff <- 1
drop <- which(apply(cpm(d0), 1, max) < cutoff)
print(length(drop))
if (length(drop) >0) {
  d <- d0[-drop,]
} else{ d <- d0}
dim(d) #number of peaks left: 8501
d <- calcNormFactors(d, method = c('TMM'))
colnames(countMatrix)
genotype <- groupList
table(genotype)
mm <- model.matrix(~0 + genotype)#
v <- voom(d, mm, plot=F, normalize="quantile")#voom from limma:Transform count data to log2-counts per million (logtpm)
fit <- lmFit(v, mm)
head(coef(fit))
contr <- makeContrasts(paste0(colnames(coef(fit)), collapse  = ' - '), levels = colnames(coef(fit))) #Change if other comparison
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)
length(which(top.table$adj.P.Val < 0.05 & abs(top.table$logFC) > log2(1.5))) #how many are significant?
res_limma <- as.data.frame(top.table)
res_limma$change <- with(res_limma, ifelse(res_limma$adj.P.Val < 0.05 & abs(res_limma$logFC) > log2(1.5),
                                           ifelse(res_limma$logFC > log2(1.5),'UP','DOWN'),'NOT'))
table(res_limma$change)#NOT only
res_limma = res_limma[res_limma$change %in% c('UP', 'DOWN'), ]

RMdem_list[[i]] <- list(res_DESeq2, res_limma)
name = res@elementMetadata$description[2]
name = paste0(cl, '_',str_split(name, "groupList ", simplify = T)[,2])
name= paste0(gsub(' ','_',(name)))
names(RMdem_list[[i]] ) <- paste0(name, '_', c('DESeq2', 'limma'))
RMdep_list[[i]] <- list( merge_peak )
names(RMdep_list[[i]] ) <- name
countMatrix_list <- list(colnames(countMatrix)); names(countMatrix_list) <- c
condition_list <- list(condition); names(condition_list) <- c
contrast_IP_list[[i]] <- countMatrix_list
condition_IP_list[[i]] <- condition_list

############
n = 114
############
#GSE211442# with RIP and Input
############
print(RMPdeg_datasets[n,])# hg38
i = RMPdeg_datasets[n,'GSE']
cl = c('HaCAT')#Human keratinocytes
RMp = c('ALKBH5') 
RM = 'm6A'
condition = c(rep('ALKBH5KD', 2), rep('WT', 2))#
RMP_allDat[grep('GSE211076', RMP_allDat$study_alias ) 
           ,c("study_alias", "sample_alias", "sample_title")]
# View(allDat[[i]])
RMP_allDat[grepl('GSE211442', RMP_allDat$study_alias ) & grepl('siControl|siALKBH5', RMP_allDat$sample_title)
           ,c("study_alias", "sample_alias", "sample_title")]

skip = RMP_allDat[grepl('GSE211442', RMP_allDat$study_alias ) & grepl('siControl|siALKBH5', RMP_allDat$sample_title)
                  ,c("study_alias", "sample_alias", "sample_title")]
paste0(skip$sample_alias, collapse = "|")

skip <- RMP_allDat[grepl('GSE211442', RMP_allDat$study_alias ) & grepl('siControl|siALKBH5', RMP_allDat$sample_title)
                   &  grepl('Input', RMP_allDat$sample_title)
                   , "sample_alias"]
paste0(skip$sample_alias, collapse = " ")
skip <- RMP_allDat[grepl('GSE211442', RMP_allDat$study_alias ) &  grepl('siControl|siALKBH5', RMP_allDat$sample_title)
                   &  grepl('IP', RMP_allDat$sample_title)
                   , "sample_alias"]
paste0(skip$sample_alias, collapse = " ")
ENA <- as.data.frame(RMP_allDat[grepl('GSE211442', RMP_allDat$study_alias ) & grepl('siControl|siALKBH5', RMP_allDat$sample_title)
                                ,c("study_alias", "sample_alias", "sample_title")])
ENA <- ENA[order(ENA$sample_alias),]
ENA
files <- list.files(file.path(data_directory,i))
files
getwd()
setwd(data_directory)
setwd(i)
c <- unique(condition[!condition %in% 'WT'])
###################################
#  Cell number-normalized (CNN)   #
#      MeRIP INPUT analysis       #
###################################
countGTF <- read.table( paste0("metaPlotR/Input_counts_", c, ".txt"),header=T)
rownames(countGTF) <- countGTF$Geneid
countGTF <- countGTF[, !colnames(countGTF) %in% c('Geneid','Chr','Start','End','Strand','Length')]
dim(countGTF) #28265     4
## Sample name
## Sample information
sample=data.frame(
  sample=colnames(countGTF),
  condition=condition)
## SummarizedExperiment input
dds<-DESeqDataSetFromMatrix(countData=countGTF,colData=sample,design=~condition)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 10
table(keep)
dds <- dds[keep,]
## Note on factor levels
dds$condition <- factor(dds$condition, levels = unique(condition))
dds$condition <- relevel(dds$condition, ref = "WT")
# Differential expression analysis-------------------------------------------
dds <- DESeq(dds)

res <-  results(dds, contrast=c("condition",unique(condition)))
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
table(res$change)
table(is.na(res$change))
# resOrdered <- as.data.frame(res[order(res$padj), ])
resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
resInP$geneSymbol <- rownames(resInP)
name = res@elementMetadata$description[2]
name = paste0(cl, '_',str_split(name, "condition ", simplify = T)[,2])
RMdeg_list[[i]] <- list( resInP)
names(RMdeg_list[[i]] ) <- paste0(gsub(' ','_',(name)))
countGTF_list <- list(colnames(countGTF)); names(countGTF_list) <- c
condition_list <- list(condition); names(condition_list) <- c
contrast_input_list[[i]] <- countGTF_list
condition_input_list[[i]] <- condition_list

###################################
#  Cell number-normalized (CNN)   #
#   MeRIP RIP analysis - Part1    #
###################################
## Sample information
merge_peak <- read.delim(paste0("metaPlotR/", c, "_newPeakFile.txt"),header=F)
countMatrix <- read.table(paste0("metaPlotR/m6A_counts_", c, ".txt"),header=T)
table(countMatrix$Geneid == merge_peak[-1,1])
table(duplicated(merge_peak$V1[-1]))
dim(countMatrix)
colnames(countMatrix)
rownames(countMatrix) <- countMatrix$Geneid
countMatrix <- countMatrix[, !colnames(countMatrix) %in% c('Geneid','Chr','Start','End','Strand','Length')]
# 1. 做好分组因子即可
groupList<-factor(condition)
colData <- data.frame(row.names=colnames(countMatrix), groupList=groupList)
# 2. Exploring Method 1: DESeq2-------------------------------------------
# 构建 DESeqDataSet 对象
dds <- DESeqDataSetFromMatrix(countData = countMatrix,
                              colData = colData,
                              design = ~ groupList)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 1
table(keep)
dds <- dds[keep,]
# 接着计算每个样本的 library size 用于后续的标准化。library size 等于落在该样本所有peaks上的reads的总数，即 counts 矩阵中每列的加和。如果bFullLibrarySize设为TRUE，则会使用library的总reads数（根据原BAM/BED文件统计获得）。然后使用estimateDispersions估计统计分布，需要将参数fitType设为’local’
dds<-estimateSizeFactors(dds)
dds<-estimateDispersions(dds,fitType="local")
dds <- nbinomWaldTest(dds)
res <- results(dds)
## 提取你想要的差异分析结果，我们这里是treated组对untreated组进行比较
res <-  results(dds, contrast=c("groupList", unique(condition)))
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
res_DESeq2 = as.data.frame(res)
table(res_DESeq2$change)#Not only
res_DESeq2 = res_DESeq2[res_DESeq2$change %in% c('UP', 'DOWN'), ]
# 2. Exploring Method 2: limma-------------------------------------------
library(edgeR)
d0 <- DGEList(counts = countMatrix)
dim(d0)
# Filtering low count peaks
cutoff <- 1
drop <- which(apply(cpm(d0), 1, max) < cutoff)
print(length(drop))
if (length(drop) >0) {
  d <- d0[-drop,]
} else{ d <- d0}
dim(d) #number of peaks left: 8501
d <- calcNormFactors(d, method = c('TMM'))
colnames(countMatrix)
genotype <- groupList
table(genotype)
mm <- model.matrix(~0 + genotype)#
v <- voom(d, mm, plot=F, normalize="quantile")#voom from limma:Transform count data to log2-counts per million (logtpm)
fit <- lmFit(v, mm)
head(coef(fit))
contr <- makeContrasts(paste0(colnames(coef(fit)), collapse  = ' - '), levels = colnames(coef(fit))) #Change if other comparison
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)
length(which(top.table$adj.P.Val < 0.05 & abs(top.table$logFC) > log2(1.5))) #how many are significant?
res_limma <- as.data.frame(top.table)
res_limma$change <- with(res_limma, ifelse(res_limma$adj.P.Val < 0.05 & abs(res_limma$logFC) > log2(1.5),
                                           ifelse(res_limma$logFC > log2(1.5),'UP','DOWN'),'NOT'))
table(res_limma$change)#NOT only
res_limma = res_limma[res_limma$change %in% c('UP', 'DOWN'), ]

RMdem_list[[i]] <- list(res_DESeq2, res_limma)
name = res@elementMetadata$description[2]
name = paste0(cl, '_',str_split(name, "groupList ", simplify = T)[,2])
name= paste0(gsub(' ','_',(name)))
names(RMdem_list[[i]] ) <- paste0(name, '_', c('DESeq2', 'limma'))
RMdep_list[[i]] <- list( merge_peak )
names(RMdep_list[[i]] ) <- name
countMatrix_list <- list(colnames(countMatrix)); names(countMatrix_list) <- c
condition_list <- list(condition); names(condition_list) <- c
contrast_IP_list[[i]] <- countMatrix_list
condition_IP_list[[i]] <- condition_list

############
n = 115
############
#GSE44386# without case RIP, with RNA
############
print(RMPdeg_datasets[n,])#
i = RMPdeg_datasets[n,'GSE']
cl = c('HEK293')
RMp = c('NSUN2')
RM = 'm5C'
condition = c(rep('NSUN2KD', 4), rep('WT', 4))#
RMP_allDat[grep('GSE44384', RMP_allDat$study_alias )
           ,c("study_alias", "sample_alias", "sample_title")]
# View(allDat[[i]])
files <- list.files(file.path(data_directory,i))
files
getwd()
setwd(data_directory)
setwd(i)
###################################
#         RNA-seq analysis        #
###################################
# Reads Per Kilobase of exon per Megabase of library size (RPKM) were calculated using cuffdiff
# 1) RNA-seq
###################################
# Read and process each count file
count <- lapply(files[grepl('Nsun_|Scr', files)], function(x) {
  # Read data using fread for faster reading
  data <- fread(file.path(data_directory, i, x))
  # Extract required columns
  df_new <- data[, 1:2]
  df_new <- df_new[!duplicated(df_new[,1])]
  # Set column names same as file names
  colnames(df_new) <- c('id',gsub('.txt.gz', '', x))
  return(df_new)
})
# Merge dataframes in countGTF list by 'gene_id'
df_in <- Reduce(function(x, y) merge(x, y, by = 'id', all = TRUE), count)
colSums(df_in[,-1])
table(colSums(df_in[,-1])=="1e+06")
## Gene name
library(data.table)
countGTF=df_in
# id_mapped <- read.table("/data2/rluo4/hg38/gencode.v44lift37.annotation.gene.probeMap",header = T)
id_mapped <- read.table("/data2/rluo4/hg38/gencode.v36.annotation.gene.probeMap",header = T)
id_mapped$id <- str_split(id_mapped$id, '[.]', simplify = T)[,1]
id_mapped <- id_mapped[which(id_mapped$id %in% countGTF$id ),]
dim(id_mapped) #57930     6
table(duplicated(id_mapped$gene))
countGTF <- left_join(id_mapped,countGTF,by="id")
countGTF <- aggregate(x = countGTF[,(ncol(id_mapped)+1):ncol(countGTF)],   #此时exprSet的第三列开始是表达矩阵内容
                      by = list(symbol = countGTF$gene),   #按照相同symbol分组，在组内计算
                      FUN = max) %>%   #原文中是计算最大值（max），也可以计算平均值（mean）或者中位数（median）
  column_to_rownames(var = 'symbol')
dim(countGTF) # 51242    12
geneSymbol <- rownames(countGTF)

## Sample name
colnames(countGTF)
## Sample information
sample=data.frame(
  sample=colnames(countGTF),
  condition=c(rep('WT', 4), rep('NSUN2KD', 4)))
## SummarizedExperiment input
dds<-DESeqDataSetFromMatrix(countData=countGTF,colData=sample,design=~condition)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 10
table(keep)
dds <- dds[keep,]
## Note on factor levels
dds$condition <- factor(dds$condition, levels = unique(condition))
dds$condition <- relevel(dds$condition, ref = "WT")
# Differential expression analysis-------------------------------------------
dds <- DESeq(dds)
res <-  results(dds, contrast=c("condition",unique(condition)))
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
table(res$change)
table(is.na(res$change))
# resOrdered <- as.data.frame(res[order(res$padj), ])
resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
resInP$geneSymbol <- rownames(resInP)
c <- unique(condition[!condition %in% 'WT'])

name = res@elementMetadata$description[2]
name = paste0(cl, '_',str_split(name, "condition ", simplify = T)[,2])
RMdeg_list[[i]] <- list( resInP)
names(RMdeg_list[[i]] ) <- paste0(gsub(' ','_',(name)))
condition_list <- list(condition); names(condition_list) <- c
contrast_input_list[[i]] <- NULL
condition_input_list[[i]] <- condition_list

############
n = 116
############
#GSE73941# with RIR and Input
############
print(RMPdeg_datasets[n,])#
i = RMPdeg_datasets[n,'GSE']
cl = c('HEK293T')
RMp = c('ALKBH3')
RM = 'm1A'
condition = c(rep('ALKBH3KO', 2), rep('WT', 2))#
RMP_allDat[grep('GSE73941', RMP_allDat$study_alias )
           ,c("study_alias", "sample_alias", "sample_title")]
# View(allDat[[i]])
RMP_allDat[grepl('GSE73941', RMP_allDat$study_alias ) & grepl('WT_|KO_', RMP_allDat$sample_title)
           ,c("study_alias", "sample_alias", "sample_title")]

skip = RMP_allDat[grepl('GSE73941', RMP_allDat$study_alias ) & grepl('WT_|KO_', RMP_allDat$sample_title)
                  ,c("study_alias", "sample_alias", "sample_title")]
paste0(skip$sample_alias, collapse = "|")

skip <- RMP_allDat[grepl('GSE73941', RMP_allDat$study_alias ) & grepl('WT_|KO_', RMP_allDat$sample_title)
                   &  grepl('Input', RMP_allDat$sample_title)
                   , "sample_alias"]
paste0(skip$sample_alias, collapse = " ")
skip <- RMP_allDat[grepl('GSE73941', RMP_allDat$study_alias ) &  grepl('WT_|KO_', RMP_allDat$sample_title)
                   &  (! grepl('Input', RMP_allDat$sample_title))
                   , "sample_alias"]
paste0(skip$sample_alias, collapse = " ")
files <- list.files(file.path(data_directory,i))
files
getwd()
setwd(data_directory)
setwd(i)
ENA <- as.data.frame(RMP_allDat[grepl('GSE73941', RMP_allDat$study_alias ) & grepl('WT_|KO_', RMP_allDat$sample_title)
                                ,c("study_alias", "sample_alias", "sample_title")])
ENA <- ENA[order(ENA$sample_alias),]
ENA
c <- unique(condition[!condition %in% 'WT'])
###################################
#  Cell number-normalized (CNN)   #
#      MeRIP INPUT analysis       #
###################################
countGTF <- read.table( paste0("metaPlotR/Input_counts_", c, ".txt"),header=T)
rownames(countGTF) <- countGTF$Geneid
countGTF <- countGTF[, !colnames(countGTF) %in% c('Geneid','Chr','Start','End','Strand','Length')]
dim(countGTF) #28265     4
## Sample name
## Sample information
sample=data.frame(
  sample=colnames(countGTF),
  condition=condition)
## SummarizedExperiment input
dds<-DESeqDataSetFromMatrix(countData=countGTF,colData=sample,design=~condition)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 10
table(keep)
dds <- dds[keep,]
## Note on factor levels
dds$condition <- factor(dds$condition, levels = unique(condition))
dds$condition <- relevel(dds$condition, ref = "WT")
# Differential expression analysis-------------------------------------------
dds <- DESeq(dds)

res <-  results(dds, contrast=c("condition",unique(condition)))
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
table(res$change)
table(is.na(res$change))
# resOrdered <- as.data.frame(res[order(res$padj), ])
resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
resInP$geneSymbol <- rownames(resInP)
name = res@elementMetadata$description[2]
name = paste0(cl, '_',str_split(name, "condition ", simplify = T)[,2])
RMdeg_list[[i]] <- list( resInP)
names(RMdeg_list[[i]] ) <- paste0(gsub(' ','_',(name)))
countGTF_list <- list(colnames(countGTF)); names(countGTF_list) <- c
condition_list <- list(condition); names(condition_list) <- c
contrast_input_list[[i]] <- countGTF_list
condition_input_list[[i]] <- condition_list

###################################
#  Cell number-normalized (CNN)   #
#   MeRIP RIP analysis - Part1    #
###################################
## Sample information
merge_peak <- read.delim(paste0("metaPlotR/", c, "_newPeakFile.txt"),header=F)
countMatrix <- read.table(paste0("metaPlotR/m1A_counts_", c, ".txt"),header=T)
table(countMatrix$Geneid == merge_peak[-1,1])
table(duplicated(merge_peak$V1[-1]))
dim(countMatrix)
colnames(countMatrix)
rownames(countMatrix) <- countMatrix$Geneid
countMatrix <- countMatrix[, !colnames(countMatrix) %in% c('Geneid','Chr','Start','End','Strand','Length')]
# 1. 做好分组因子即可
groupList<-factor(condition)
colData <- data.frame(row.names=colnames(countMatrix), groupList=groupList)
# 2. Exploring Method 1: DESeq2-------------------------------------------
# 构建 DESeqDataSet 对象
dds <- DESeqDataSetFromMatrix(countData = countMatrix,
                              colData = colData,
                              design = ~ groupList)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 1
table(keep)
dds <- dds[keep,]
# 接着计算每个样本的 library size 用于后续的标准化。library size 等于落在该样本所有peaks上的reads的总数，即 counts 矩阵中每列的加和。如果bFullLibrarySize设为TRUE，则会使用library的总reads数（根据原BAM/BED文件统计获得）。然后使用estimateDispersions估计统计分布，需要将参数fitType设为’local’
dds<-estimateSizeFactors(dds)
dds<-estimateDispersions(dds,fitType="local")
dds <- nbinomWaldTest(dds)
res <- results(dds)
## 提取你想要的差异分析结果，我们这里是treated组对untreated组进行比较
res <-  results(dds, contrast=c("groupList", unique(condition)))
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
res_DESeq2 = as.data.frame(res)
table(res_DESeq2$change)#Not only
res_DESeq2 = res_DESeq2[res_DESeq2$change %in% c('UP', 'DOWN'), ]
# 2. Exploring Method 2: limma-------------------------------------------
library(edgeR)
d0 <- DGEList(counts = countMatrix)
dim(d0)
# Filtering low count peaks
cutoff <- 1
drop <- which(apply(cpm(d0), 1, max) < cutoff)
print(length(drop))
if (length(drop) >0) {
  d <- d0[-drop,]
} else{ d <- d0}
dim(d) #number of peaks left: 8501
d <- calcNormFactors(d, method = c('TMM'))
colnames(countMatrix)
genotype <- groupList
table(genotype)
mm <- model.matrix(~0 + genotype)#
v <- voom(d, mm, plot=F, normalize="quantile")#voom from limma:Transform count data to log2-counts per million (logtpm)
fit <- lmFit(v, mm)
head(coef(fit))
contr <- makeContrasts(paste0(colnames(coef(fit)), collapse  = ' - '), levels = colnames(coef(fit))) #Change if other comparison
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)
length(which(top.table$adj.P.Val < 0.05 & abs(top.table$logFC) > log2(1.5))) #how many are significant?
res_limma <- as.data.frame(top.table)
res_limma$change <- with(res_limma, ifelse(res_limma$adj.P.Val < 0.05 & abs(res_limma$logFC) > log2(1.5),
                                           ifelse(res_limma$logFC > log2(1.5),'UP','DOWN'),'NOT'))
table(res_limma$change)#NOT only
res_limma = res_limma[res_limma$change %in% c('UP', 'DOWN'), ]

RMdem_list[[i]] <- list(res_DESeq2, res_limma)
name = res@elementMetadata$description[2]
name = paste0(cl, '_',str_split(name, "groupList ", simplify = T)[,2])
name= paste0(gsub(' ','_',(name)))
names(RMdem_list[[i]] ) <- paste0(name, '_', c('DESeq2', 'limma'))
RMdep_list[[i]] <- list( merge_peak )
names(RMdep_list[[i]] ) <- name
countMatrix_list <- list(colnames(countMatrix)); names(countMatrix_list) <- c
condition_list <- list(condition); names(condition_list) <- c
contrast_IP_list[[i]] <- countMatrix_list
condition_IP_list[[i]] <- condition_list

############
n = 117
############
#GSE198643# with RIP and Input
############
print(RMPdeg_datasets[n,])#
i = RMPdeg_datasets[n,'GSE']
cl = c('HK2')#Normal kidney
RMp = c('METTL3')
RM = 'm6A'
condition = c(rep('METTL3KD', 2), rep('WT', 2))#
RMP_allDat[grep('GSE198643', RMP_allDat$study_alias )
           ,c("study_alias", "sample_alias", "sample_title")]
# View(allDat[[i]])
files <- list.files(file.path(data_directory,i))
files
getwd()
setwd(data_directory)
setwd(i)
c <- unique(condition[!condition %in% 'WT'])
###################################
#  Cell number-normalized (CNN)   #
#      MeRIP INPUT analysis       #
###################################
countGTF <- read.table( paste0("metaPlotR/Input_counts_", c, ".txt"),header=T)
rownames(countGTF) <- countGTF$Geneid
countGTF <- countGTF[, !colnames(countGTF) %in% c('Geneid','Chr','Start','End','Strand','Length')]
dim(countGTF) #28265     4
## Sample name
## Sample information
sample=data.frame(
  sample=colnames(countGTF),
  condition=condition)
## SummarizedExperiment input
dds<-DESeqDataSetFromMatrix(countData=countGTF,colData=sample,design=~condition)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 10
table(keep)
dds <- dds[keep,]
## Note on factor levels
dds$condition <- factor(dds$condition, levels = unique(condition))
dds$condition <- relevel(dds$condition, ref = "WT")
# Differential expression analysis-------------------------------------------
dds <- DESeq(dds)

res <-  results(dds, contrast=c("condition",unique(condition)))
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
table(res$change)
table(is.na(res$change))
# resOrdered <- as.data.frame(res[order(res$padj), ])
resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
resInP$geneSymbol <- rownames(resInP)
name = res@elementMetadata$description[2]
name = paste0(cl, '_',str_split(name, "condition ", simplify = T)[,2])
RMdeg_list[[i]] <- list( resInP)
names(RMdeg_list[[i]] ) <- paste0(gsub(' ','_',(name)))
countGTF_list <- list(colnames(countGTF)); names(countGTF_list) <- c
condition_list <- list(condition); names(condition_list) <- c
contrast_input_list[[i]] <- countGTF_list
condition_input_list[[i]] <- condition_list

###################################
#  Cell number-normalized (CNN)   #
#   MeRIP RIP analysis - Part1    #
###################################
## Sample information
merge_peak <- read.delim(paste0("metaPlotR/", c, "_newPeakFile.txt"),header=F)
countMatrix <- read.table(paste0("metaPlotR/m6A_counts_", c, ".txt"),header=T)
table(countMatrix$Geneid == merge_peak[-1,1])
table(duplicated(merge_peak$V1[-1]))
dim(countMatrix)
colnames(countMatrix)
rownames(countMatrix) <- countMatrix$Geneid
countMatrix <- countMatrix[, !colnames(countMatrix) %in% c('Geneid','Chr','Start','End','Strand','Length')]
# 1. 做好分组因子即可
groupList<-factor(condition)
colData <- data.frame(row.names=colnames(countMatrix), groupList=groupList)
# 2. Exploring Method 1: DESeq2-------------------------------------------
# 构建 DESeqDataSet 对象
dds <- DESeqDataSetFromMatrix(countData = countMatrix,
                              colData = colData,
                              design = ~ groupList)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 1
table(keep)
dds <- dds[keep,]
# 接着计算每个样本的 library size 用于后续的标准化。library size 等于落在该样本所有peaks上的reads的总数，即 counts 矩阵中每列的加和。如果bFullLibrarySize设为TRUE，则会使用library的总reads数（根据原BAM/BED文件统计获得）。然后使用estimateDispersions估计统计分布，需要将参数fitType设为’local’
dds<-estimateSizeFactors(dds)
dds<-estimateDispersions(dds,fitType="local")
dds <- nbinomWaldTest(dds)
res <- results(dds)
## 提取你想要的差异分析结果，我们这里是treated组对untreated组进行比较
res <-  results(dds, contrast=c("groupList", unique(condition)))
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
res_DESeq2 = as.data.frame(res)
table(res_DESeq2$change)#Not only
res_DESeq2 = res_DESeq2[res_DESeq2$change %in% c('UP', 'DOWN'), ]
# 2. Exploring Method 2: limma-------------------------------------------
library(edgeR)
d0 <- DGEList(counts = countMatrix)
dim(d0)
# Filtering low count peaks
cutoff <- 1
drop <- which(apply(cpm(d0), 1, max) < cutoff)
print(length(drop))
if (length(drop) >0) {
  d <- d0[-drop,]
} else{ d <- d0}
dim(d) #number of peaks left: 8501
d <- calcNormFactors(d, method = c('TMM'))
colnames(countMatrix)
genotype <- groupList
table(genotype)
mm <- model.matrix(~0 + genotype)#
v <- voom(d, mm, plot=F, normalize="quantile")#voom from limma:Transform count data to log2-counts per million (logtpm)
fit <- lmFit(v, mm)
head(coef(fit))
contr <- makeContrasts(paste0(colnames(coef(fit)), collapse  = ' - '), levels = colnames(coef(fit))) #Change if other comparison
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)
length(which(top.table$adj.P.Val < 0.05 & abs(top.table$logFC) > log2(1.5))) #how many are significant?
res_limma <- as.data.frame(top.table)
res_limma$change <- with(res_limma, ifelse(res_limma$adj.P.Val < 0.05 & abs(res_limma$logFC) > log2(1.5),
                                           ifelse(res_limma$logFC > log2(1.5),'UP','DOWN'),'NOT'))
table(res_limma$change)#NOT only
res_limma = res_limma[res_limma$change %in% c('UP', 'DOWN'), ]

RMdem_list[[i]] <- list(res_DESeq2, res_limma)
name = res@elementMetadata$description[2]
name = paste0(cl, '_',str_split(name, "groupList ", simplify = T)[,2])
name= paste0(gsub(' ','_',(name)))
names(RMdem_list[[i]] ) <- paste0(name, '_', c('DESeq2', 'limma'))
RMdep_list[[i]] <- list( merge_peak )
names(RMdep_list[[i]] ) <- name
countMatrix_list <- list(colnames(countMatrix)); names(countMatrix_list) <- c
condition_list <- list(condition); names(condition_list) <- c
contrast_IP_list[[i]] <- countMatrix_list
condition_IP_list[[i]] <- condition_list

contrast_input_list[['GSE198643']]$METTL3KD <- gsub('M.KD', 'M-KD', contrast_input_list[['GSE198643']]$METTL3KD )
contrast_IP_list[['GSE198643']]$METTL3KD <- gsub('M.KD', 'M-KD', contrast_IP_list[['GSE198643']]$METTL3KD )
############
n = 118
############
#GSE78509# without RIP, with RNA
############
print(RMPdeg_datasets[n,])#
i = RMPdeg_datasets[n,'GSE']
cl = c('hESCs')
RMp = c('IGF2BP1', 'IGF2BP2', 'IGF2BP3')
RM = 'm6A'
condition = c(rep('IGF2BP1KD', 3), rep('WT', 2))#
RMP_allDat[grep('GSE78508', RMP_allDat$study_alias )
           ,c("study_alias", "sample_alias", "sample_title")]
# View(allDat[[i]])
# gse <- getGEO('GSE78509',destdir=".", AnnotGPL=F, getGPL=F)
# pdata <- pData(gse[[2]])
RMP_allDat[grep('GSE78507', RMP_allDat$study_alias )
           ,c("study_alias", "sample_alias", "sample_title")]
files <- list.files(file.path(data_directory,i))
files
getwd()
setwd(data_directory)
setwd(i)
###################################
# Read and process each count file
count <- lapply(files[grepl('_esc_', files)], function(x) {
  # Read data using fread for faster reading
  data <- fread(file.path(data_directory, i, x))
  # Extract required columns
  # Set column names same as file names
  colnames(data)[1] <- c('geneSymbol')
  return(data)
})
# Merge dataframes in countGTF list by 'gene_id'
df_in <- Reduce(function(x, y) merge(x, y, by = 'geneSymbol', all = TRUE), count)
colSums(df_in[,-1])
library(data.table)
colnames(df_in)
countGTF=df_in[, -1]
# id_mapped <- read.table("/data2/rluo4/hg38/gencode.v44lift37.annotation.gene.probeMap",header = T)
id_mapped <- read.table("/data2/rluo4/hg38/gencode.v36.annotation.gene.probeMap",header = T)
countGTF$id <- str_split(df_in$geneSymbol, '[.]', simplify = T)[,1]
# id_mapped <- unique(id_mapped) #%>% dplyr::select(id,gene) # 去重
id_mapped$id <- str_split(id_mapped$id, '[.]', simplify = T)[,1]
id_mapped <- id_mapped[which(id_mapped$id %in% countGTF$id ),]
dim(id_mapped) #57930     6
table(duplicated(id_mapped$gene))
countGTF <- left_join(id_mapped,countGTF,by="id")
countGTF <- aggregate(x = countGTF[,(ncol(id_mapped)+1):ncol(countGTF)],   #此时exprSet的第三列开始是表达矩阵内容
                      by = list(symbol = countGTF$gene),   #按照相同symbol分组，在组内计算
                      FUN = max) %>%   #原文中是计算最大值（max），也可以计算平均值（mean）或者中位数（median）
  column_to_rownames(var = 'symbol')
dim(countGTF) # 51242    12
geneSymbol <- rownames(countGTF)
## Gene name
dim(countGTF) #25343    12
## Sample name
colnames(countGTF)
## Sample information
sample=data.frame(
  sample=colnames(countGTF),
  condition=c( rep('WT', 2),rep('IGF2BP1KD', 3)))
## SummarizedExperiment input
dds<-DESeqDataSetFromMatrix(countData=countGTF,colData=sample,design=~condition)
## Pre-filtering
keep <- rowSums(counts(dds)) >= 10
table(keep)
dds <- dds[keep,]
## Note on factor levels
dds$condition <- factor(dds$condition, levels = unique(condition))
dds$condition <- relevel(dds$condition, ref = "WT")
# Differential expression analysis-------------------------------------------
dds <- DESeq(dds)
res <-  results(dds, contrast=c("condition",unique(condition)))
res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                               ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
table(res$change)
table(is.na(res$change))
# resOrdered <- as.data.frame(res[order(res$padj), ])
resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
resInP$geneSymbol <- rownames(resInP)
c <- unique(condition[!condition %in% 'WT'])
name = res@elementMetadata$description[2]
name = paste0(cl, '_',str_split(name, "condition ", simplify = T)[,2])
RMdeg_list[[i]] <- list( resInP)
names(RMdeg_list[[i]] ) <- paste0(gsub(' ','_',(name)))
condition_list <- list(condition); names(condition_list) <- c
contrast_input_list[[i]] <- NULL
condition_input_list[[i]] <- condition_list

############
n = 119
############
#GSE149989# without RIP with Input
############
print(RMPdeg_datasets[n,])#
i = RMPdeg_datasets[n,'GSE']
cl = c('LNZ308') #Giloma
RMp = c('METTL1')
RM = 'm7G'
condition = c(rep('METTL1KD', 2), rep('WT', 2))#
RMP_allDat[grep('GSE78509', RMP_allDat$study_alias )
           ,c("study_alias", "sample_alias", "sample_title")]
# View(allDat[[i]])
files <- list.files(file.path(data_directory,i))
files
getwd()
setwd(data_directory)
setwd(i)
df_in <- read_excel('GSE149989_LNZ308_METTL1_KD_Trac-seq.xlsx')
countGTF <- as.data.frame(df_in[, -ncol(df_in)])
table(duplicated(df_in$tRNAs))
res <- countGTF
colnames(res)[8:10] <- c('FoldChange','logFC', 'pvalue')
table(countGTF$`p value`<0.05)
res$change <- with(res, ifelse(abs(res$logFC) > log2(1.5),
                               ifelse(res$logFC > log2(1.5),'UP','DOWN'),'NOT'))
table(res$change)
resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
resInP$geneSymbol <- resInP$tRNAs
table(resInP$change)
table(resInP$FoldChange %in% c(0, Inf))
table(resInP$logFC %in% c(-Inf, Inf))
resInP_list <- list()
resInP_list[['LNZ308_METTL1KD_vs_WT']] <- resInP
RMdeg_list[[i]] <- resInP_list 
condition_list <- list(condition); names(condition_list) <- 'METTL1KD'
contrast_input_list[[i]] <- NULL#countGTF_list
condition_input_list[[i]] <- condition_list

############
n = 120
############
#GSE210563# without RIP and Input
############
print(RMPdeg_datasets[n,])#
i = RMPdeg_datasets[n,'GSE']
cl = c('Hela and HEK293T')
cl = 'Hela'
RMp = c('METTL3', 'METTL14', 'WTAP', 'YTHDF2')
RM = 'm6A'
condition_list <- list(METTL3KD = c(rep('METTL3KD', 2), rep('WT', 4)),
                      METTL14KD = c(rep('METTL14KD', 2), rep('WT', 4)),
                      WTAPKD = c(rep('WTAPKD', 2), rep('WT', 4)),
                      YTHDF2KD = c(rep('YTHDF2KD', 2), rep('WT', 4)))
RMP_allDat[grep('GSE210563', RMP_allDat$study_alias )
           ,c("study_alias", "sample_alias", "sample_title")]
# View(allDat[[i]])
# gse <- getGEO('GSE210563',destdir=".", AnnotGPL=F, getGPL=F)
# pdata <- pData(gse[[3]])
files <- list.files(file.path(data_directory,i))
files
getwd()
setwd(data_directory)
setwd(i)
count_files <- files[grepl('_expression.csv.gz', files)]
count_files <- count_files[grepl('-si', count_files)]
###################################
# Read and process each count file
count <- lapply(count_files, function(x) {
  data <- read.delim(file.path(data_directory, i, x))
  df_new <- data[, c('Gene','Counts')]
  # Set column names same as file names
  colnames(df_new) <- c('geneSymbol',gsub('_expression.csv.gz', '', x))
  df_new <- df_new[df_new[,1] !='-', ]
  df_new[,1] <- str_split(df_new[,1], ',', simplify = T)[,1]
  df_new <- df_new[df_new[,1] !='-', ]
  df_new <- aggregate(x = df_new,   #此时exprSet的第三列开始是表达矩阵内容
                      by = list(geneSymbol = df_new$geneSymbol),   #按照相同symbol分组，在组内计算
                      FUN = max) %>%   #原文中是计算最大值（max），也可以计算平均值（mean）或者中位数（median）
    column_to_rownames(var = 'geneSymbol')
  
  return(df_new)
})
names(count) <- gsub('_expression.csv.gz', '', count_files)
# Merge dataframes in countGTF list by 'gene_id'
df_in <- Reduce(function(x, y) merge(x, y, by = 'geneSymbol', all = TRUE), count)
###################################
#   more than 1 condition here   #
###################################
resInP_list <- list()
colnames(df_in)
##Gene count matrix
countGTF=as.matrix(df_in[, -1])
## Gene name
table(duplicated(df_in$geneSymbol))
dim(countGTF) #25343    12
## Sample name
colnames(countGTF)
## Sample information
sample_meta=data.frame(
  sample=colnames(countGTF),
  condition=c(rep("METTL3KD",2),rep("METTL14KD",2), rep("WTAPKD",2),rep("WT",2), rep("YTHDF2KD",2),rep("WT",2)))
###################################
resInP_list <- vector("list", length(condition_list))
names(resInP_list) <- names(condition_list)
# countGTF_list <- vector("list", length(condition_list))
# names(countGTF_list) <- names(condition_list)
for(c in names(condition_list)){
  condition <- condition_list[[c]]
  ## SummarizedExperiment input
  count <- as.matrix(countGTF[, grepl(paste0(substr(c, 1, nchar(c)-2), '|siNC'), colnames(countGTF))])
  colnames(count)
  table(is.na(count))
  sample <- sample_meta[sample_meta$sample %in% colnames(count), ]
  dds<-DESeqDataSetFromMatrix(countData=count, colData=sample,design=~condition)
  ## Pre-filtering
  keep <- rowSums(counts(dds)) >= 10
  table(keep)
  dds <- dds[keep,]
  ## Note on factor levels
  dds$condition <- factor(dds$condition, levels = unique(condition))
  dds$condition <- relevel(dds$condition, ref = "WT")
  # Differential expression analysis-------------------------------------------
  dds <- DESeq(dds)
  res <-  results(dds, contrast=c("condition",unique(condition)))
  res$label <- with(res, ifelse(res$baseMean>10 & abs(res$log2FoldChange)>log2(1.5) & res$padj<0.05,"sig","nosig"))
  length(which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))) #how many are significant?
  res$change <- with(res, ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5),
                                 ifelse(res$log2FoldChange > log2(1.5),'UP','DOWN'),'NOT'))
  res$geneSymbol <- df_in$geneSymbol[keep]
  table(res$change)
  table(is.na(res$change))
  # resOrdered <- as.data.frame(res[order(res$padj), ])
  resInP = as.data.frame(res[res$change %in% c('UP', 'DOWN'), ])
  table(resInP$change)
  resInP_list[[c]] <- resInP
}
name = paste0(cl, '_',names(resInP_list), ' vs WT')
RMdeg_list[[i]] <- resInP_list 
names(RMdeg_list[[i]] ) <- paste0(gsub(' ','_',(name)))
contrast_input_list[[i]] <- NULL
condition_input_list[[i]] <- condition_list

save(RMPdeg_datasets, gene_info,
     RMdeg_list, RMdem_list, RMdep_list, 
     contrast_input_list, contrast_IP_list, condition_input_list, condition_IP_list,
     file = paste0('/data2/rluo4/EpiTrans/RMDatasets/RMdeg_allDat.RData'))


# 3) the details about RMPdeg_datasets
genome_version <- data.frame(GSE = RMPdeg_datasets$GSE, HG = 'hg19', skip = '')
hg38_list <- c(9,12,13,15,20,21,22,24,30,34,35,36,41,44,45,46,47,53,58,59,65,84,90,95,103,104,108,112,114)
skip_list <- c(23, 26, 27, 31, 32, 52, 54, 55, 60, 61, 70, 72, 83, 91 )
genome_version$HG[hg38_list] <- 'hg38'
genome_version$skip[skip_list] <- 'NoData'
table(genome_version$skip[skip_list])
genome_version$macs2_dir <- genome_version$GSE
genome_version$macs2_dir[genome_version$GSE=='GSE90963'] = 'PRJNA355164'
genome_version$macs2_dir[genome_version$GSE=='GSE121942'] = 'GSE110320'
genome_version$macs2_dir[genome_version$GSE=='GSE103497'] = 'GSE103495'
genome_version$macs2_dir[genome_version$GSE=='GSE122803'] = 'GSE122801'
genome_version$macs2_dir[genome_version$GSE=='GSE120024'] = 'GSE132306'
genome_version$macs2_dir[genome_version$GSE=='GSE128575'] = 'GSE128520'
genome_version$macs2_dir[genome_version$GSE=='GSE202815'] = 'GSE202814'
genome_version$macs2_dir[genome_version$GSE=='GSE223731'] = 'GSE223729'
genome_version$macs2_dir[genome_version$GSE=='GSE90684'] = 'GSE90642'
# load the target annotation results #
# load(file.path(data_directory,'../target_results.RData'))
target_results <-  vector("list", 120)
names(target_results) <- RMPdeg_datasets$GSE
for (n in 1:nrow(RMPdeg_datasets)) {#n = 11, 36, 40, 60
  print(n)
  print(RMPdeg_datasets[n, 'GSE'])
  i <- RMPdeg_datasets[n, 'GSE']
  # outdir <- file.path(plot_dir, i)
  # if(! dir.exists( outdir )) {
  #   dir.create(outdir)
  # }
  result <- process_samples(n, RMPdeg_datasets, data_directory, RMdem_list, RMdep_list, genome_version)
  target_results[[i]] <- result
}
save(target_results, file = paste0('/data2/rluo4/EpiTrans/RMDatasets/target_results.RData'))
