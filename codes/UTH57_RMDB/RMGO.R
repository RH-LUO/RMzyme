#!/usr/local/bin/Rscript
# title: "RBM33 is a unique m(6)A RNA-binding protein that regulates ALKBH5 demethylase activity and substrate selectivity"
# author: "Ruihan Luo"
# date: "April 19th,2024"
# rm(list=ls())
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
# Function to extract the relevant columns and combine them into a summary table
extract_gene_summary <- function(filtered_deg_list, dataset_summary) {
  # Initialize an empty list to store results
  summary_list <- list()
  # Iterate over each dataset in the list
  for (dataset_name in names(filtered_deg_list)) {
    dataset <- filtered_deg_list[[dataset_name]]
    # Iterate over each condition in the dataset
    for (condition_name in names(dataset)) {
      condition_data <- dataset[[condition_name]]
      if(nrow(condition_data)==0){
        print(paste0(dataset_name, '-', condition_name, ": no DEG data in RMdeg_list !"))
        next;
      }
      # Check for the columns "log2FoldChange" or "logFC"
      if ("log2FoldChange" %in% colnames(condition_data)) {
        logFC_col <- "log2FoldChange"
      } else if ("logFC" %in% colnames(condition_data)) {
        logFC_col <- "logFC"
      } else {
        # Skip this condition if neither column is present
        next
      }
      # Extract the required columns
      extracted_data <- condition_data[, c("geneSymbol", "change", logFC_col)]
      colnames(extracted_data) <- c('GeneSymbol', 'Regulation', 'log2FoldChange')
      # Add additional columns for dataset and condition name
      extracted_data$Condition <- condition_name
      extracted_data$Dataset <- dataset_name
      Enzyme <-  gsub('_vs_WT', '', condition_name)
      remove_first_element <- function(char_vector) {
        char_vector_parts <- strsplit(char_vector, "_")
        char_vector_parts <- lapply(char_vector_parts, function(parts) {
          if (length(parts) > 1) {
            return(paste(parts[-1], collapse = "_"))
          } else {
            return(parts[1])
          }
        })
        return(unlist(char_vector_parts))
      }
      extracted_data$Enzyme <-  remove_first_element(Enzyme) #(str_split(Enzyme, '_', simplify = TRUE)[, 2])
      extracted_data$Treatment <-  substr(extracted_data$Enzyme, nchar(extracted_data$Enzyme)-1, nchar(extracted_data$Enzyme))
      extracted_data$Treatment <-  gsub('KD', 'knock-down',  gsub('MU', 'mutant', gsub('RS', 'rescue',  gsub('IB', 'inhibition',gsub('PD', 'pull-down',
                                                                                                                                     gsub('KO', 'knock-out', gsub('OE', 'over-expression', extracted_data$Treatment)
                                                                                                                                     ))))))
      extracted_data$Enzyme <- substr(extracted_data$Enzyme, 1, nchar(extracted_data$Enzyme)-2)
      extracted_data$Treatment <- paste0(extracted_data$Enzyme, ' ', extracted_data$Treatment)
      extracted_data$Modification_type <- dataset_summary$Modifcation[match(dataset_name, dataset_summary$GSE)]
      extracted_data$Cell_line <- (str_split(condition_name, '_', simplify = TRUE)[, 1])
      extracted_data$Tissue_disease_type <- dataset_summary$Tissue[match(dataset_name, dataset_summary$GSE)]
      extracted_data$Technique <- dataset_summary$Technique[match(dataset_name, dataset_summary$GSE)]
      # Append the extracted data to the list
      summary_list[[paste(dataset_name, condition_name, sep = "_")]] <- extracted_data
    }
  }
  # Combine all extracted data into a single data frame
  summary_df <- do.call(rbind, summary_list)
  # Return the summary data frame
  return(summary_df)
}
extract_peak_summary <- function(filtered_dem_list, dataset_summary) {
  # Initialize an empty list to store results
  summary_list <- list()
  # Iterate over each dataset in the list
  for (dataset_name in names(filtered_dem_list)) {
    target_res <- target_results[[dataset_name]]
    length(target_res)
    if( length(target_res) ==0 ){
      print(paste0(dataset_name, ": no data in target_results !"))
      next;
    }
    resM <- filtered_dem_list[[dataset_name]]
    resP <- filtered_dep[[dataset_name]]
    # table(target_res)
    cl <- unique(str_split( names(target_res), '_', simplify = TRUE)[, 1])
    # if (is.null(resM) | is.null(resP)) {
    #   print(paste0(i, ": no RNA Peak data!"))
    #   # condition_name = names(resM)
    #   # condition_data <- resM[[condition_name]]
    #   # if ("log2FoldChange" %in% colnames(condition_data)) {
    #   #   logFC_col <- "log2FoldChange"
    #   # } else if ("logFC" %in% colnames(condition_data)) {
    #   #   logFC_col <- "logFC"
    #   # } else {
    #   #   # Skip this condition if neither column is present
    #   #   next
    #   # }
    #   # extracted_data <- condition_data[, c("geneSymbol", "change", logFC_col)]
    #   # colnames(extracted_data) <- c('GeneSymbol', 'Regulation', 'log2FoldChange')
    #   return(NULL)  # Exit the function early if no data
    # }
    for (contrast in names(target_res)) {
      # contrast <- names(target_res)[1]
      c <- unique(str_split(contrast, '_', simplify = TRUE)[, 2])
      if( length(cl) >1 | dataset_name=='GSE207643'){
        c <- contrast
      }
      condition_name <- paste0(contrast, '_vs_WT')
      print(condition_name)
      target <- target_res[[contrast]]
      DEM <- resM[[paste0(condition_name, '_DESeq2')]]
      index = grepl('DESeq2', names(resM))
      if (length(unique(index)) > 1 & !is.null(DEM) ) {
        DEM <- resM[[paste0(condition_name, '_DESeq2')]]
      } else{
        DEM <- resM[[paste0(condition_name, '_limma')]]
      }
      if (length(unique(index)) > 1 & !is.null(DEM) & nrow(DEM) == 0) {
        DEM <- resM[[paste0(condition_name, '_limma')]]
      }
      DEM$peak_id <- rownames(DEM)
      DEM <- left_join(target, DEM, by = 'peak_id')
      condition_data <- as.data.frame(DEM[!is.na(DEM$change),])
      if(nrow(condition_data)==0){
        print(paste0(dataset_name, '-', condition_name, ": no dem data in RMdem_list !"))
        next;
      }
      # Check for the columns "log2FoldChange" or "logFC"
      if ("log2FoldChange" %in% colnames(condition_data)) {
        logFC_col <- "log2FoldChange"
      } else if ("logFC" %in% colnames(condition_data)) {
        logFC_col <- "logFC"
      } else {
        # Skip this condition if neither column is present
        next
      }
      # Extract the required columns
      cols <-  c('geneSymbol', 'chr', 'start', 'end', 'peak_id', logFC_col, "change")
      extracted_data <- condition_data[, match(cols, colnames(condition_data))]
      colnames(extracted_data) <- c('GeneSymbol','Chr', 'Start', 'End', 'Peak_id',  'log2FoldChange', 'Regulation')
      # Add additional columns for dataset and condition name
      extracted_data$Condition <- condition_name
      extracted_data$Dataset <- dataset_name
      Enzyme <-  gsub('_vs_WT', '', condition_name)
      remove_first_element <- function(char_vector) {
        char_vector_parts <- strsplit(char_vector, "_")
        char_vector_parts <- lapply(char_vector_parts, function(parts) {
          if (length(parts) > 1) {
            return(paste(parts[-1], collapse = "_"))
          } else {
            return(parts[1])
          }
        })
        return(unlist(char_vector_parts))
      }
      extracted_data$Enzyme <-  remove_first_element(Enzyme) #(str_split(Enzyme, '_', simplify = TRUE)[, 2])
      extracted_data$Treatment <-  substr(extracted_data$Enzyme, nchar(extracted_data$Enzyme)-1, nchar(extracted_data$Enzyme))
      extracted_data$Treatment <-  gsub('KD', 'knock-down',  gsub('MU', 'mutant', gsub('RS', 'rescue', gsub('IB', 'inhibition',gsub('PD', 'pull-down',
                                                                                                                                    gsub('KO', 'knock-out', gsub('OE', 'over-expression', extracted_data$Treatment)
                                                                                                                                    ))))))
      extracted_data$Enzyme <- substr(extracted_data$Enzyme, 1, nchar(extracted_data$Enzyme)-2)
      extracted_data$Treatment <- paste0(extracted_data$Enzyme, ' ', extracted_data$Treatment)
      extracted_data$Modification_type <- dataset_summary$Modifcation[match(dataset_name, dataset_summary$GSE)]
      extracted_data$Cell_line <- (str_split(condition_name, '_', simplify = TRUE)[, 1])
      extracted_data$Tissue_disease_type <- dataset_summary$Tissue[match(dataset_name, dataset_summary$GSE)]
      extracted_data$Technique <- dataset_summary$Technique[match(dataset_name, dataset_summary$GSE)]
      # Append the extracted data to the list
      summary_list[[paste(dataset_name, condition_name, sep = "_")]] <- extracted_data
    }
  }
  # Combine all extracted data into a single data frame
  summary_df <- do.call(rbind, summary_list)
  # Return the summary data frame
  return(summary_df)
}

#######################################
# 4.2 GO analysis for each condition  #
#      nohup Rscript RMGO.R &         #
#######################################
# Apply the function to the filtered_deg list
DEG_summary <- extract_gene_summary(filtered_deg_list = filtered_deg, dataset_summary = final_datasets)
# Save the summary table to a CSV file if needed
DEG_summary <- DEG_summary[!DEG_summary$GeneSymbol %in% c('-','--'),]
DEG_summary$Cell_line <- gsub('HeLa', 'Hela', DEG_summary$Cell_line)
DEG_summary$Cell_line <- gsub('HaCAT', 'HaCaT', DEG_summary$Cell_line)
index <- DEG_summary$Condition == 'Hela.m6A_ALKBH3KO_vs_WT'
DEG_summary$Modification_type[index] <- 'm6A'
write.table(DEG_summary, "/data2/rluo4/RPMfunc/Output/summary/DEG_summary.txt", row.names = FALSE, sep = '\t', quote = F)
# Apply the function to the filtered_dem list
DEM_summary <- extract_peak_summary(filtered_dem_list = filtered_dem, dataset_summary = final_datasets)
DEM_summary$Cell_line <- gsub('HeLa', 'Hela', DEM_summary$Cell_line)
DEM_summary$Cell_line <- gsub('HaCAT', 'HaCaT', DEM_summary$Cell_line)
index <- DEM_summary$Condition == 'Hela.m6A_ALKBH3KO_vs_WT'
DEM_summary$Modification_type[index] <- 'm6A'
write.table(DEM_summary, "/data2/rluo4/RPMfunc/Output/summary/DEM_summary.txt", row.names = FALSE, sep = '\t', quote = F)
table(DEG_summary$GeneSymbol %in% gene_info$symbols)
geneSymbol = unique(union(DEG_summary$GeneSymbol, DEM_summary$GeneSymbol))
setdiff(unique(RMP_update$RMP), geneSymbol)
geneSymbol =  union( geneSymbol, unique(RMP_update$RMP))
gene_summary <- data.frame(Index = 1:length(geneSymbol), 
                           geneSymbol = geneSymbol)
colnames(gene_info)[1] <- 'geneSymbol'
gene_summary <- left_join(gene_summary, gene_info, by = 'geneSymbol')
write.table(gene_summary,'/data2/rluo4/RPMfunc/Output/summary/gene_summary.txt',sep='\t',quote = F,row.names = F)
library(writexl)
write_xlsx(final_datasets[, -ncol(final_datasets)],'/data2/rluo4/RPMfunc/Output/summary/dataset_summary.xlsx')#,sep='\t',quote = F,row.names = F)
###################################
# GO analysis for each condition  #
###################################
set.seed(211)
library(clusterProfiler)
library(msigdbr)
genesets <- msigdbr_collections()
genesets[1:20,]
msigdbr_species() #列出有的物种
genesets <- msigdbr(species = "Homo sapiens", category = "C5")
genesets <- genesets[! genesets$gs_subcat %in% 'HPO',]
genesets <- subset(genesets, select = c('gs_name','gene_symbol')) %>% as.data.frame()
go_genesets <- genesets
# go_genesets <- split(genesets$gene_symbol, genesets$gs_name)
genesets <- msigdbr(species = "Homo sapiens", category = "C2")
genesets <- subset(genesets, select = c('gs_name','gene_symbol')) %>% as.data.frame()
kegg_genesets <-  genesets[grep("KEGG",genesets$gs_name),]
# ###############################################################################
# genesets <- msigdbr(species = "Homo sapiens", category = "C8")
# genesets <- subset(genesets, select = c('gs_name','gene_symbol')) %>% as.data.frame()
# cell_genesets <- genesets
# # cell_genesets  <-  split(genesets$gene_symbol, genesets$gs_name)
# genesets <- msigdbr(species = 'Homo sapiens', category = 'H')
# genesets <- subset(genesets, select = c('gs_name','gene_symbol')) %>% as.data.frame()
# H_genesets <- genesets
# SenMayo = read.table("/data2/rluo4/All/SenMayo.txt",fill=T,skip = F,header = T,sep = '\t')
# FRIDMAN_SENESCENCE_UP <- read.csv('/data2/rluo4/All/FRIDMAN_SENESCENCE_UP.txt',header = F)
# EpiSen <- read.csv('/data2/rluo4/All/EpiSen.txt',header = F)
# SenMayo$Term <- rep('SenMayo',nrow(SenMayo))
# colnames(SenMayo)[1] <- 'V1'
# FRIDMAN_SENESCENCE_UP$Term <- rep('FRIDMAN_SENESCENCE_UP',nrow(FRIDMAN_SENESCENCE_UP))
# EpiSen$Term <- rep('EpiSen',nrow(EpiSen))
# 
# Senescense <- rbind(EpiSen[,c(2,1)],FRIDMAN_SENESCENCE_UP[,c(2,1)],SenMayo[,c(4,1)])
# colnames(Senescense) <- colnames(H_genesets)
# 
# Metaplasia <- read.csv('/data2/rluo4/All/Metaplasia.txt',header = F)
# Metaplasia$Term <- rep('Metaplasia',nrow(Metaplasia))
# Metaplasia <- Metaplasia[,c(2,1)]
# colnames(Metaplasia) <- colnames(H_genesets)
# 
# CancerG0Arrest.set <- read.csv('/home/lorihan/SCI-supp/CancerG0Arrest/QuiescenceBiomarkers.csv')
# reactome.gmt <- read.gmt('/data2/rluo4/NSCLC/c2.cp.reactome.v7.0.symbols.gmt.txt')
# cellcycle.gmt <- reactome.gmt[grep("G1",reactome.gmt$term),]
# 
# head(CancerG0Arrest.set)
# CancerG0Arrest.set$Term <- rep('CancerG0Arrest',nrow(CancerG0Arrest.set))
# G0Arrest.gmt <- CancerG0Arrest.set[,c(4,1)]
# colnames(G0Arrest.gmt) <- c("term","gene")
# cellcycle.gmt <- G0Arrest.gmt#rbind(G0Arrest.gmt,cellcycle.gmt)
# cellcycle.gmt$term <- as.factor(cellcycle.gmt$term)
# unique(cellcycle.gmt$term)
# #cellcycle.gmt <-  split(cellcycle.gmt$gene, cellcycle.gmt$term)
# colnames(cellcycle.gmt) <- colnames(H_genesets)
# colnames(reactome.gmt) <- colnames(H_genesets)
# my_genesets <- rbind(go_genesets, Senescense, Metaplasia, H_genesets, cellcycle.gmt, reactome.gmt)
# py_genesets <- rbind(kegg_genesets, Senescense, Metaplasia, H_genesets, cellcycle.gmt, reactome.gmt)
# write.csv(py_genesets,'/data2/rluo4/All/my_genesets.csv')
# ###############################################################################
my_genesets <- read.csv('/data2/rluo4/All/my_genesets.csv')
my_genesets <- my_genesets[-grep("KEGG", my_genesets$gs_name),-1]
my_genesets <- rbind(go_genesets, my_genesets)
# # my_genesets <- rbind(Senescense, Metaplasia, H_genesets, cellcycle.gmt, reactome.gmt)
# colnames(my_genesets) <- c("term","gene")
# my_genesets$term <- as.factor(my_genesets$term)
# unique(my_genesets$term)
# my.gmt <- split(my_genesets$gene, my_genesets$term)
#############################################################################################################################
for (n in 1:nrow(RMPdeg_datasets)) {
  print(n)
  print(RMPdeg_datasets[n, 'GSE'])
  i <- RMPdeg_datasets[n, 'GSE']
  print(i)
  if( genome_version$skip[genome_version$GSE==i] =='NoData' ){
    print(paste0(i, ": no data finally, skip  !"))
    next;
  }
  outdir <- file.path(plot_dir, i)
  print(outdir)
  if(! dir.exists( outdir )) {
    dir.create(outdir)
  }
  ###################################
  DEG_subset <- DEG_summary %>% filter(Dataset ==i)
  print(unique(DEG_subset$Condition))
  DEM_subset <- DEM_summary %>% filter(Dataset ==i)
  print(unique(DEG_subset$Condition))
  index1 <- DEG_subset %>% nrow()
  index2 <- DEM_subset %>% nrow()
  if( index1 ==0 & index2 ==0 ){
    print(paste0(i, ": no data in DEG/DEM_summary !"))
    next;
  }
  condition = union(DEG_subset$Condition, DEM_subset$Condition)

  for (contrast in condition) {
    # contrast <- condition[1]
    DEG_res <- DEG_subset %>% filter(Condition ==contrast)
    DEM_res <- DEM_subset %>% filter(Condition ==contrast)
    index1 <- DEG_res %>% nrow()
    index2 <- DEM_res %>% nrow()
    start_time <- Sys.time()
    # print(paste("start GO analysis:", contrast, 'from', i, 'at', start_time))
    dir_name <- gsub('_vs_WT', '', contrast)
    outdir <- file.path(plot_dir, i, dir_name)
    if(! dir.exists( outdir )) {
      dir.create(outdir)
    }
    print(outdir)
    inter = intersect(DEG_res$GeneSymbol, DEM_res$GeneSymbol)
    print(paste0('intersection of DEG and DEM: ',length(unique(inter))))
    if( index1 !=0 ){
      write.table(DEG_res, file = paste0(outdir, '/show_DEG.txt'), sep = '\t', quote = F)
     }else{
      print(paste0(i, ": no data in show_DEG.txt !"))
    }
    if( index2 !=0 ){
      write.table(DEM_res, file = paste0(outdir, '/show_DEM.txt'), sep = '\t', quote = F)
   }else{
      print(paste0(i, ": no data in show_DEM.txt !"))
    }
  }
}

# Load required libraries
library(AnnotationDbi)
library(org.Hs.eg.db)
library(clusterProfiler)
library(dplyr)
library(ggplot2)
# library(openxlsx)  # For writing Excel files
library(Hmisc)

# Function to perform GO enrichment analysis for each dataset
perform_GO_enrichment <- function(dataset_row, DEG_summary, DEM_summary, plot_dir, genome_version, my_genesets, cutoff = 0.05) {
  GSE <- dataset_row[,'GSE']
  print(paste("Processing:", GSE))
  # # Check if the dataset is marked to be skipped
  # if (genome_version$skip[genome_version$GSE == GSE] == 'NoData') {
  #   print(paste0(GSE, ": no data available, skipping!"))
  #   return(NULL)
  # }
  # 
  # Create output directory if it doesn't exist
  outdir <- file.path(plot_dir, GSE)
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }
  print(outdir)
  
  # Filter DEG and DEM data for the current dataset
  DEG_subset <- DEG_summary %>% filter(Dataset == GSE)
  DEM_subset <- DEM_summary %>% filter(Dataset == GSE)
  
  if (nrow(DEG_subset) == 0 && nrow(DEM_subset) == 0) {
    print(paste0(GSE, ": no data in DEG/DEM summary!"))
    return(NULL)
  }
  
  conditions <- union(DEG_subset$Condition, DEM_subset$Condition)
  modules_DEG <- list()
  PID_DEG <- list()
  modules_DEM <- list()
  PID_DEM <- list()
  
  # Process each condition within the dataset
  for (condition in conditions) {
    DEG_res <- DEG_subset %>% filter(Condition == condition)
    DEM_res <- DEM_subset %>% filter(Condition == condition)
    print(paste("Start GO analysis for condition:", condition, "from", GSE, "at", Sys.time()))
    
    dir_name <- gsub('_vs_WT', '', condition)
    condition_outdir <- file.path(outdir, dir_name)
    if (!dir.exists(condition_outdir)) {
      dir.create(condition_outdir)
    }
    print(condition_outdir)
    
    inter_genes <- intersect(DEG_res$GeneSymbol, DEM_res$GeneSymbol)
    print(paste0('Intersection of DEG and DEM: ', length(unique(inter_genes))))
    
    # Perform enrichment analysis for DEG
    if (nrow(DEG_res) != 0) {
      modules_DEG[[dir_name]] <- unique(DEG_res$GeneSymbol)
      set.seed(123)
      PID_DEG[[dir_name]] <- enricher(gene = modules_DEG[[dir_name]], TERM2GENE = my_genesets, pvalueCutoff = 1)
    } else {
      print(paste0(GSE, ": no data in DEG for condition ", condition))
    }
    
    # Perform enrichment analysis for DEM
    if (nrow(DEM_res) != 0) {
      modules_DEM[[dir_name]] <- unique(DEM_res$GeneSymbol)
      set.seed(123)
      PID_DEM[[dir_name]] <- enricher(gene = modules_DEM[[dir_name]], TERM2GENE = my_genesets, pvalueCutoff = 1)
    } else {
      print(paste0(GSE, ": no data in DEM for condition ", condition))
    }
  }
  
  list(PID_DEG = PID_DEG, PID_DEM = PID_DEM)
}

# Initialize result lists
PID_DEG_res <- list()
PID_DEM_res <- list()
# Loop over each dataset row to perform the GO enrichment analysis
for (n in 1:nrow(RMPdeg_datasets)) {
  print(n)
  print(RMPdeg_datasets[n, 'GSE'])
  i <- RMPdeg_datasets[n, 'GSE']
  print(i)
  if( genome_version$skip[genome_version$GSE==i] =='NoData' ){
    print(paste0(i, ": no data finally, skip  !"))
    next;
  }
  result <- perform_GO_enrichment(RMPdeg_datasets[n, ], DEG_summary, DEM_summary, plot_dir, genome_version, my_genesets)
  if (!is.null(result)) {
    GSE <- RMPdeg_datasets[n, 'GSE']
    PID_DEG_res[[GSE]] <- result$PID_DEG
    PID_DEM_res[[GSE]] <- result$PID_DEM
  }
}
save(DEG_summary, DEM_summary, PID_DEG_res, PID_DEM_res, file = paste0('/data2/rluo4/EpiTrans/RMDatasets/GOanalysis.RData'))
# load(paste0('/data2/rluo4/EpiTrans/RMDatasets/GOanalysis.RData'))
# Process the results and save them to Excel files
cutoff <- 0.05
shorten_names <- function(x, n_word=4, n_char=40){
  if (length(strsplit(x, " ")[[1]]) > n_word || (nchar(x) > 40))
  {
    if (nchar(x) > 40) x <- substr(x, 1, 40)
    x <- paste(paste(strsplit(x, " ")[[1]][1:min(length(strsplit(x," ")[[1]]), n_word)],
                     collapse=" "), "...", sep="")
    return(x)
  }
  else
  {
    return(x)
  }
}

library(openxlsx)
for (i in seq_along(PID_DEG_res)) {
  GSE <- names(PID_DEG_res)[i]
  deg_results <- PID_DEG_res[[GSE]]
  
  for (j in names(deg_results)) {
    # j = names(deg_results)[1]
    fullname <- paste0(j ,'_vs_WT')
    outdir <- file.path(plot_dir, GSE, j)
    if(! dir.exists( outdir )) {
      dir.create(outdir)
    }
    deg_result <- deg_results[[j]]
    if (!is.null(deg_result) && nrow(deg_result@result) > 0) {
      # sheet_name <- paste("MP", i, j, sep = '_')
      GO_final <- deg_result@result[deg_result@result$pvalue < cutoff, -c(1)]
      write.table(GO_final, file = file.path(outdir, 'show_DEG_GO.txt'), sep = '\t', quote = F, row.names = F)
      # write.xlsx(GO_final,
      #            file = file.path(outdir, 'show_DEG_GO.xlsx'),
      #            # sheetName = sheet_name,
      #            append = TRUE)
      go_dat=GO_final[order(GO_final$GeneRatio, decreasing = T),]
      go_dat = go_dat[1:10,]
      go_dat$Description = sapply(go_dat$Description ,shorten_names)
      g_ego<-ggplot(data=go_dat)+
        geom_bar(aes(x=reorder(Description,-log10(pvalue)),y=-log10(pvalue), fill=-log10(pvalue)), stat='identity') +
        coord_flip() +
        scale_fill_gradient(expression(-log["10"](pvalue)),low="#7FFFD4", high = "#FF6A6A") +
        labs(x="",y="Enrichment score(-log10(pvalue))",title=paste0("Top enriched  terms for ", fullname )) +#xlab("")+ ylab("Enrichment score(-log10(pvalue))") +
        scale_y_continuous(expand=c(0, 0))+
        theme_bw()+theme(panel.grid=element_blank(),axis.ticks=element_line(color='black'),
                         #axis.text=element_text(colour="black",size=8),
                         panel.background = element_rect(fill='transparent'),
                         #panel.border=element_rect(fill='transparent'),
                         axis.line=element_line(color='black'),
                         axis.title=element_text(color='black',size=16),
                         plot.title = element_text(size=16,hjust = 0.5))+
        theme(axis.text.x=element_text(color="black",size=rel(1.6)),
              axis.text.y=element_text(color="black", size=rel(1.6)),
              axis.title.x = element_text(color="black", size=rel(1)),
              # legend.text=element_text(color="black",size=rel(1.0)),
              # legend.title = element_text(color="black",size=rel(1.1))
              legend.text = element_text(size = 16),
              legend.title=element_text(size = 16))
      
      hei <- nrow(go_dat)*0.6; wid <- nrow(go_dat)*1.2
      file =  file.path(outdir, 'show_DEG_GO.png')
      png(file, width = wid, height = hei, res = 500,units = 'in')
      print(g_ego)
      dev.off()
      
    }
  }
}
for (i in seq_along(PID_DEM_res)) {
  GSE <- names(PID_DEM_res)[i]
  DEM_results <- PID_DEM_res[[GSE]]
  
  for (j in names(DEM_results)) {
    # j = names(DEM_results)[1]
    fullname <- paste0(j ,'_vs_WT')
    outdir <- file.path(plot_dir, GSE, j)
    if(! dir.exists( outdir )) {
      dir.create(outdir)
    }
    DEM_result <- DEM_results[[j]]
    if (!is.null(DEM_result) && nrow(DEM_result@result) > 0) {
      GO_final <- DEM_result@result[DEM_result@result$pvalue < cutoff, -c(1)]
      write.table(GO_final, file = file.path(outdir, 'show_DEM_GO.txt'), sep = '\t', quote = F, row.names = F)
      # write.xlsx(GO_final,
      #            file = file.path(outdir, 'show_DEM_GO.xlsx'),
      #            # sheetName = sheet_name,
      #            append = TRUE)
      go_dat=GO_final[order(GO_final$GeneRatio, decreasing = T),]
      go_dat = go_dat[1:10,]
      go_dat$Description = sapply(go_dat$Description ,shorten_names)
      g_ego<-ggplot(data=go_dat)+
        geom_bar(aes(x=reorder(Description,-log10(pvalue)),y=-log10(pvalue), fill=-log10(pvalue)), stat='identity') +
        coord_flip() +
        scale_fill_gradient(expression(-log["10"](pvalue)),low="#20B2AA", high = "#CD5C5C") +
        labs(x="",y="Enrichment score(-log10(pvalue))",title=paste0("Top enriched terms for ", fullname )) +#xlab("")+ ylab("Enrichment score(-log10(pvalue))") +
        scale_y_continuous(expand=c(0, 0))+
        theme_bw()+theme(panel.grid=element_blank(),axis.ticks=element_line(color='black'),
                         #axis.text=element_text(colour="black",size=8),
                         panel.background = element_rect(fill='transparent'),
                         #panel.border=element_rect(fill='transparent'),
                         axis.line=element_line(color='black'),
                         axis.title=element_text(color='black',size=16),
                         plot.title = element_text(size=16,hjust = 0.5))+
        theme(axis.text.x=element_text(color="black",size=rel(1.6)),
              axis.text.y=element_text(color="black", size=rel(1.6)),
              axis.title.x = element_text(color="black", size=rel(1)),
              # legend.text=element_text(color="black",size=rel(1.0)),
              # legend.title = element_text(color="black",size=rel(1.1))
              legend.text = element_text(size = 16),
              legend.title=element_text(size = 16))
      
      hei <- nrow(go_dat)*0.6; wid <- nrow(go_dat)*1.2
      file =  file.path(outdir, 'show_DEM_GO.png')
      png(file, width = wid, height = hei, res = 500,units = 'in')
      print(g_ego)
      dev.off()
    }
  }
}

