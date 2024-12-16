#!/usr/local/bin/Rscript
# title: "RMzyme, a comprehensive platform of regulations of RNA modifying enzymes in human based on multiomics data collections"
# author: "Ruihan Luo"
# date: "April 19th,2024"
# rm(list=ls())
library(dplyr)
library(tidyr)
library(stringr)
library(data.table)
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

##################################
# 1) load the Rdata from RMdeg.R #
##################################
load(paste0('/data2/rluo4/EpiTrans/RMDatasets/allDat.RData')) # pdata from RMDatasets.R
load(paste0('/data2/rluo4/EpiTrans/RMDatasets/RMP_allDat.RData')) # asp metadata from RMDatasets_UTH36.R
load(paste0('/data2/rluo4/EpiTrans/RMDatasets/RMdeg_allDat.RData')) # from RMdeg.R
RMP_update <- read_excel('/data2/rluo4/EpiTrans/DataCollection/RMP_update.xlsx') # from RM.R
asp_sc <- fread(file="/data2/rluo4/EpiTrans/DataCollection/asp_RMP_epitrans.link.txt", header = F)

#########################
# 2) load the functions #
#########################
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
#####################################
#      nohup Rscript RMdeg.R &      #
#####################################
# target_results <-  vector("list", 120) # final length is 73
# names(target_results) <- RMPdeg_datasets$GSE
# for (n in 104:nrow(RMPdeg_datasets)) {#n = 11, 36, 40, 60
#   print(n)
#   print(RMPdeg_datasets[n, 'GSE'])
#   i <- RMPdeg_datasets[n, 'GSE']
#   
#   outdir <- file.path(plot_dir, i)
#   if(! dir.exists( outdir )) {
#     dir.create(outdir)
#   }
#   result <- process_samples(n, RMPdeg_datasets, data_directory, RMdem_list, RMdep_list, genome_version)
#   target_results[[i]] <- result
# }
# save(target_results, file = paste0('/data2/rluo4/EpiTrans/RMDatasets/target_results.RData'))
#####################################
# load the target annotation results #
load(file.path(data_directory,'../target_results.RData'))
########################################
# 3) the details about RMPdeg_datasets #
########################################
library(stringr)
library(dplyr)
genome_version <- data.frame(GSE = RMPdeg_datasets$GSE, HG = 'hg19', skip = '')
hg38_list <- c(9,12,13,15,20,21,22,24,30,34,35,36,41,44,45,46,47,53,58,59,65,84,90,95,103,104,108,112,114)
skip_list <- c(23, 26, 27, 31, 32, 52, 54, 55, 60, 61, 70, 72, 83, 91 )
genome_version$HG[hg38_list] <- 'hg38'
genome_version$skip[skip_list] <- 'NoData'
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
# ideogram_hg19 <- geom_ideogram(genome = 'hg19', plot.space = 0, highlight.centromere = TRUE)
# ideogram_hg38 <- geom_ideogram(genome = 'hg38', plot.space = 0, highlight.centromere = TRUE)
# save(ideogram_hg19, ideogram_hg38, file = paste0('/data2/rluo4/EpiTrans/RMDatasets/getIdeogram.RData') )
# load(paste0('/data2/rluo4/EpiTrans/RMDatasets/getIdeogram.RData') )
final_datasets <- RMPdeg_datasets[-skip_list,]
# write.xlsx(x = final_datasets, file = paste0('/data2/rluo4/EpiTrans/DataCollection/final_datasets.xlsx'),   sheetName = 'RNA-modification-Datasets'#,append = TRUE
#  )
library(writexl)
# write_xlsx(x = final_datasets, path = paste0('/data2/rluo4/EpiTrans/DataCollection/final_datasets.xlsx'))

######################################################
# 4) The visualization of differential binding peaks #
######################################################
###################################
#  4.1 MeRIP RIP visualization    #
#    nohup Rscript PeakVis.R &    #
###################################
library(GenomicFeatures)
library(Guitar)
library(ggcoverage)
# hg38.bed <- read.delim("/data2/rluo4/lorihan/hg38/hg38.txt",sep ="",header = FALSE)
# rownames(hg38.bed)<-hg38.bed$V4
# hg19.bed <- read.delim("/data2/rluo4/lorihan/hg38/hg19.txt",sep ="",header = FALSE)
# rownames(hg19.bed) <- hg19.bed$V4
# To add gene annotation, the gtf file should contain gene_type and gene_name attributes in column 9; to add transcript annotation, the gtf file should contain transcript_name attribute in column 9.
library(Gviz)
library(CAGEfightR)
library(rtracklayer)
library("ggcoverage")
library("ggpattern")
library("BSgenome.Hsapiens.UCSC.hg19")
library("BSgenome.Hsapiens.UCSC.hg38")
BSv19 = BSgenome.Hsapiens.UCSC.hg19
BSv38 = BSgenome.Hsapiens.UCSC.hg38

original <- c('GSE97419', 'GSE55572', 'GSE87190', 'GSE76414',
              'GSE103497', 'GSE141994', 'GSE145924', 'GSE124509',
              'GSE210867',  'GSE207643', 'GSE171497', 'GSE198643')
bw.rank <- match(original, RMPdeg_datasets$GSE)
final_datasets <- read_excel(paste0('/data2/rluo4/EpiTrans/DataCollection/final_datasets.xlsx'))
final_datasets$Technique[final_datasets$GSE=='GSE122413'] = 'RNA-Bisulfite-Seq'
final_datasets$Tissue <- (str_split(final_datasets$Cell.Line.Tissue, ';', simplify = TRUE)[, 2])
final_datasets$Tissue <- gsub('MA9.3ITD','Leukemia', final_datasets$Tissue )
final_datasets$Modifcation <- (str_split(final_datasets$Modifcation, ';', simplify = TRUE)[, 1])
filtered_deg <- Filter(Negate(is.null), RMdeg_list)
filtered_dem <- Filter(Negate(is.null), RMdem_list)
filtered_dep <- Filter(Negate(is.null), RMdep_list)

analyzed_gse <- unique(union(names(filtered_dem ), names(filtered_deg)))
setdiff(RMPdeg_datasets$GSE,  analyzed_gse)
setdiff(names(filtered_dem ), names(filtered_dep ))#[1] "GSE174492"

###############################################
# 4.6 Volcano visualization for each condition #
###############################################
library(EnhancedVolcano) #BiocManager::install("EnhancedVolcano")
# Loop over each dataset row to process the results and save them to PNG files
library(openxlsx)
plot_dir <- '/data2/rluo4/RPMfunc/Output/GEO_top' # rm -r */*/Guitar
for (i in seq_along(filtered_deg)) {
  GSE <- names(filtered_deg)[i]
  deg_results <- filtered_deg[[GSE]]
  
  for (contrast in names(deg_results)) {
    # contrast = names(deg_results)[1]
    deg_result <- deg_results[[contrast]]
    j = gsub('_vs_WT', '', contrast)
    outdir <- file.path(plot_dir, GSE, j)
    if(! dir.exists( outdir )) {
      dir.create(outdir)
    }
    print(outdir)
    print(j)
    # res <- deg_result
    print(colnames(deg_result))
    
    if( 'logFC' %in% colnames(deg_result)){
      colnames(deg_result) <- gsub('logFC', 'log2FoldChange', colnames(deg_result))
    }
    target_columns <- c("baseMean", "log2FoldChange", "lfcSE", "stat", "pval", "padj")#, "label", "change", "geneSymbol")
    actual_columns <- colnames(deg_result)
    # Find the partial matches between the target column names and actual column names in your second file
    matched_columns <- sapply(target_columns, function(target) {
      # Use grepl to find columns that partially match the target column names
      match <- actual_columns[grepl(target, actual_columns, ignore.case = TRUE)]
      if (length(match) == 0) return(NA)  # Return NA if no match is found
      return(match[1])  # Return the first matched column
    })
    # Create a named list for renaming matched columns
    matched_columns <- matched_columns[!is.na(matched_columns)]  # Remove NAs if any
    # Extract the matched columns from your dataframe
    res <- deg_result %>%
      dplyr::select(all_of(c(matched_columns, "change", "geneSymbol")))
    print(colnames(res))
    # If necessary, fill in missing columns with NA, such as 'lfcSE' and 'stat' (if not matched)
    if (!"baseMean" %in% colnames(res)) {
      res$baseMean <- NA
    }
    if (!"lfcSE" %in% colnames(res)) {
      res$lfcSE <- NA
    }
    if (!"stat" %in% colnames(res)) {
      res$stat <- NA
    }
    # if (!"pval" %in% colnames(res)) {
    #   res$pval <- 0.01
    # }
    if (!"padj" %in% colnames(res)) {
      res$padj <- NA
    }
    
    if ( nrow(res) > 0 && "pval" %in% colnames(res)) {
      target_columns <- c("baseMean", "log2FoldChange", "lfcSE", "stat", "pval", "padj", "change", "geneSymbol")
      res <- res[,match(target_columns, colnames(res))]
      # Rename the extracted columns to match the target names
      names(res) <- gsub('pval', 'pvalue', names(res))
      # View the final dataframe
      head(res)
      print(table(res$pvalue==0))
      pvalues <- sort(unique(res$pvalue))
      log2FCs <- sort(unique(res$log2FoldChange))
      res$pvalue[res$pvalue==0] <- pvalues[2]/2
      res$log2FoldChange[res$log2FoldChange==-Inf] <- log2FCs[2]*2
      res$log2FoldChange[res$log2FoldChange==Inf] <- log2FCs[length(log2FCs)-1]*2
      
      index <- abs(res$log2FoldChange) > log2(2) & res$pvalue < 0.05
      index1 <- order(abs(res$log2FoldChange), decreasing = T)
      index2 <- order(abs(res$pvalue), decreasing = F)
      key1_genes <- res$geneSymbol[index1][index & index1][1:15]
      key2_genes <- res$geneSymbol[index2][index & index2][1:15]
      # key_genes <- res$geneSymbol[index1][index & index1]
      # genes_to_label <-c("ANXA1", 'FPR1', key_genes[1:50])  # Replace with actual gene names
      genes_to_label <- unique(c(key1_genes, key2_genes)) 
      print(genes_to_label)
      if(GSE %in% c("GSE144984", "GSE94613")){
        genes_to_label <- unique(c("ANXA1", 'FPR1', 'LINC00324', key1_genes, key2_genes))  # Replace with actual gene names
      }
      hei <- 7.5; wid <- 9
      output_files  =  file.path(outdir, 'show_DEG_volcano.png')
      print( paste0("drawing volcano plot: ", output_files) )
      # Show a volcano plot and only label the key genes
      png(file = output_files , width = wid, height = hei, res = 500,units = 'in')
      p <- EnhancedVolcano(res,
                      lab = res$geneSymbol,
                      x = 'log2FoldChange',
                      y = 'pvalue',
                      selectLab = genes_to_label,
                      # pCutoff = 10e-14,
                      # FCcutoff = 1.0,
                      xlab = bquote(~Log[2]~ 'fold change'),
                      # ylab = "",
                      axisLabSize = 15,
                      pCutoff = 0.01,
                      FCcutoff = 1,
                      pointSize = 1.0,
                      labSize = 5,
                      colAlpha = 1,
                      legendPosition = 'top',
                      legendLabSize = 13, # legendInset = 0.1,          # Adjust the distance of the legend from the plot (increase/decrease this value)
                      legendIconSize = 5.0,
                      title = "",#paste0('Volcano plot for DEGs of ', gsub('_', ' ', contrast)),
                      titleLabSize = 13,    # Adjust the title size here
                      # Remove the subtitle and footnote
                      subtitle = NULL,     # Remove subtitle
                      caption = NULL,      # Remove footnote
                      #添加Connectors
                      drawConnectors = TRUE,
                      widthConnectors = 0.75)
      print(p)
      dev.off()
      
    }
  }
}
