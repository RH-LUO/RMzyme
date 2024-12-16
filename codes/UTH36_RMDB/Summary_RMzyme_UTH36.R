#!/usr/local/bin/Rscript
# title: "RMzyme, a comprehensive platform of regulations of RNA modifying enzymes in human based on multiomics data collections"
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
# path <- c('/home/rluo4/R/x86_64-pc-linux-gnu-library/4.3', '/opt/R/4.3.1/lib64/R/library', '/home/rluo4/R/x86_64-conda-linux-gnu-library/4.3', '/data/rluo4/bin/miniconda3/lib/R/library')
# .libPaths(path)
raw_directory <- '/data/rluo4/EpiTrans/RMDatasets/Raw'
data_directory <- '/data/rluo4/EpiTrans/RMDatasets/GEO'
setwd(data_directory)
plot_dir <- '/data/rluo4/RPMfunc/Output/GEO'


##################################
# 1) load the Rdata from RMdeg.R #
##################################
load(paste0('/data/rluo4/EpiTrans/RMDatasets/allDat.RData')) # pdata from RMDatasets.R
load(paste0('/data/rluo4/EpiTrans/RMDatasets/RMP_allDat.RData')) # asp metadata from RMDatasets_UTH36.R
load(paste0('/data/rluo4/EpiTrans/RMDatasets/RMdeg_allDat.RData')) # from RMdeg.R
RMP_update <- read_excel('/data/rluo4/EpiTrans/DataCollection/RMP_update.xlsx') # from RM.R
asp_sc <- fread(file="/data/rluo4/EpiTrans/DataCollection/asp_RMP_epitrans.link.txt", header = F)

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
  
  db <- '/data/rluo4/hg38/'
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
# save(target_results, file = paste0('/data/rluo4/EpiTrans/RMDatasets/target_results.RData'))
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
# save(ideogram_hg19, ideogram_hg38, file = paste0('/data/rluo4/EpiTrans/RMDatasets/getIdeogram.RData') )
# load(paste0('/data/rluo4/EpiTrans/RMDatasets/getIdeogram.RData') )
final_datasets <- RMPdeg_datasets[-skip_list,]
# write.xlsx(x = final_datasets, file = paste0('/data/rluo4/EpiTrans/DataCollection/final_datasets.xlsx'),   sheetName = 'RNA-modification-Datasets'#,append = TRUE
#  )
library(writexl)
# write_xlsx(x = final_datasets, path = paste0('/data/rluo4/EpiTrans/DataCollection/final_datasets.xlsx'))

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
# hg38.bed <- read.delim("/data/rluo4/lorihan/hg38/hg38.txt",sep ="",header = FALSE)
# rownames(hg38.bed)<-hg38.bed$V4
# hg19.bed <- read.delim("/data/rluo4/lorihan/hg38/hg19.txt",sep ="",header = FALSE)
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
final_datasets <- read_excel(paste0('/data/rluo4/EpiTrans/DataCollection/final_datasets.xlsx'))
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
library(writexl)
# write_xlsx(final_datasets[, -ncol(final_datasets)],'/data/rluo4/RPMfunc/Output/summary/dataset_summary.xlsx')#,sep='\t',quote = F,row.names = F)
#######################################
# 4.2 GO analysis for each condition  #
#      nohup Rscript RMGO.R &         #
#######################################
# Process the results and save them to Excel files
# save(PID_DEG_res, PID_DEM_res, file = paste0('/data/rluo4/EpiTrans/RMDatasets/GOanalysis.RData'))
load(paste0('/data/rluo4/EpiTrans/RMDatasets/GOanalysis.RData'))
DEM_statistic <- data.frame(table(DEM_summary$Dataset))
# table(DEM_summary$Dataset, DEM_summary$Condition)
# Process the results and save them to Excel files
table(DEG_summary$GeneSymbol %in% gene_info$symbols)
geneSymbol = unique(union(DEG_summary$GeneSymbol, DEM_summary$GeneSymbol))
setdiff(unique(RMP_update$RMP), geneSymbol)
geneSymbol =  union( geneSymbol, unique(RMP_update$RMP))
gene_summary <- data.frame(Index = 1:length(geneSymbol), 
                           geneSymbol = geneSymbol)
colnames(gene_info)[1] <- 'geneSymbol'
gene_summary <- left_join(gene_summary, gene_info, by = 'geneSymbol')
uniprot <- read.csv('/data/rluo4/All/uniprot-hs.tsv',fill = T,header = T,sep = '\t')
table(gene_summary$UniProtAcc %in% uniprot$Entry)
table(gene_summary$geneSymbol %in% uniprot$Gene.Names)
uniprot$geneSymbol <- str_split(uniprot$Gene.Names, ' ', simplify = T)[,1]
#######################################
# # gene_loc$localization <- uniprot$Subcellular.location..CC.[match(gene_loc$UniProtAcc,uniprot$Entry)]
# table(gene_summary$geneSymbol %in% uniprot$geneSymbol)
# gene_loc <- left_join(gene_summary, uniprot, by = 'geneSymbol')
# gene_loc$uniprotLoc <- gene_loc$Subcellular.location..CC. #uniprot$Subcellular.location..CC.[match(gene_loc$geneSymbol,uniprot$Gene.Names)]
# table(gene_loc$uniprotLoc=='')#133382
# gene_loc$uniprotLoc <- gsub('SUBCELLULAR LOCATION: ','',gene_loc$uniprotLoc)
# gene_loc$Nucleus <- ifelse(grepl('ucleus',gene_loc$uniprotLoc),1,0)
# gene_loc$Secreted <- ifelse(grepl('ecreted',gene_loc$uniprotLoc),1,0)
# gene_loc$Cytoplasm <- ifelse(grepl('ytoplasm',gene_loc$uniprotLoc),1,0)
# gene_loc$Surface <- ifelse(grepl('ell membrane',gene_loc$uniprotLoc),1,0)
# multi.loc <- apply(gene_loc[, c('Nucleus','Secreted','Cytoplasm','Surface')],1,function(gene_loc){
#   s <- sum(gene_loc)
#   return(s)
# })
# dim(gene_loc);length(multi.loc)
# gene_loc$multi.loc <- multi.loc  
# table(gene_loc$multi.loc)
# table(gene_loc$multi.loc %in% c(0,1))
# gene_loc <- gene_loc[gene_loc$multi.loc %in% c(0,1),]
# gene_loc$uniprotLoc[is.na(gene_loc$uniprotLoc)] = ''
# table(gene_loc$uniprotLoc=='')#3638
# 
# length(strsplit(gene_loc$uniprotLoc[100], "[;]")[[1]])
# gene_loc$proteinLoc[gene_loc$Nucleus==1] <- 'Nucleus'
# gene_loc$proteinLoc[gene_loc$Secreted==1] <- 'Secreted'
# gene_loc$proteinLoc[gene_loc$Cytoplasm==1] <- 'Cytoplasm'
# gene_loc$proteinLoc[gene_loc$Surface==1] <- 'Surface'
# # table(gene_loc$proteinLoc)
# gene_loc$proteinLoc <- str_split(gene_loc$proteinLoc,'[;]',simplify = T)[,1]
# gene_loc$proteinLoc <- str_split(gene_loc$proteinLoc,'[.]',simplify = T)[,1]
# gene_loc$proteinLoc <- str_split(gene_loc$proteinLoc,'[{]',simplify = T)[,1]
# gene_loc$proteinLoc <- str_split(gene_loc$proteinLoc,'[,]',simplify = T)[,1]
# table(gene_loc$proteinLoc)
# table(gene_loc$proteinLoc[gene_loc$proteinLoc %ni% c('Nucleus','Secreted','Cytoplasm','Surface')])
# gene_loc$proteinLoc[gene_loc$proteinLoc %ni% c('Nucleus','Secreted','Cytoplasm','Surface')] <- 'Others'
# table(gene_loc$proteinLoc)
# colnames(gene_loc)
####################################
# 4.3 Top 1000 Peak visualization  #
#  nohup Rscript PeakVis-top.R &   #
####################################
####################################
# 4.4 Top 200 Peak visualization   #
#  nohup Rscript Copu-PeakVis.R &  #
####################################
# Example: n = 65/GSE144984 and n = 42 / GSE94613 
# DEG_res <- DEG_subset %>% filter(Condition == contrast) %>% arrange(desc(abs(log2FoldChange))) #%>%  mutate(row_id = row_number()) 
# DEM_res <- DEM_subset %>% filter(Condition == contrast) %>% arrange(desc(abs(log2FoldChange))) #%>%  mutate(row_id = row_number()) 
# DEG_res <- head(DEG_res, 200)#1000)
# DEM_res <- head(DEM_res, 200)#1000)
###############################################
# 4.5 Guitar visualization for each condition #
###############################################
# Load GTF
# gtf_v36 <- rtracklayer::import.gff(con = '/data/rluo4/hg38/gencode.v36.annotation.gtf.gz', format = "gtf")
# gtf_v19 <- rtracklayer::import.gff(con = '/data/rluo4/hg38/gencode.v19.annotation.gtf.gz', format = "gtf")
load( paste0('/data/rluo4/EpiTrans/RMDatasets/txdb_Guitar.RData') )
txdb_hg19 <- makeTxDbFromGFF(file = '/data/rluo4/hg38/gencode.v19.annotation.gtf.gz',
                             format="gtf",  dataSource="Ensembl", organism="Homo sapiens")
txdb_hg38 <- makeTxDbFromGFF(file = '/data/rluo4/hg38/gencode.v36.annotation.gtf.gz',
                             format="gtf",  dataSource="Ensembl", organism="Homo sapiens")
# save(gtf_v19, gtf_v36, file = paste0('/data/rluo4/EpiTrans/RMDatasets/txdb_Guitar.RData') )
options(stringsAsFactors = F)
# load package
library(Guitar)
package.version("Guitar")
# [1] "2.10.0"
getwd()
setwd(data_directory)
# i='GSE122948'
# setwd(i)
summit_bed_files <- data.frame()
for (n in 1:nrow(RMPdeg_datasets)) {
  # for (n in bw.rank) {
  print(n)
  print(RMPdeg_datasets[n, 'GSE'])
  i <- RMPdeg_datasets[n, 'GSE']
  files <- list.files(file.path(data_directory, i, 'metaPlotR'))
  # files <- files[grep('summits.bed', files)]
  files <- files[grep('summits.bed', files)]
  if(length(files)==0){
    print('no summit.bed !')
    next;
  }
  summit_bed <- data.frame(rank = n, GSE = i, summit.bed = files)
  summit_bed_files <<- rbind(summit_bed_files, summit_bed)
}
RMP_pattern <- paste(c(RMP_update$RMP,'WT'), collapse = "|")
# Extract the matching words
matches <- str_extract(summit_bed_files$summit.bed, RMP_pattern)
# Add the matches as a new column in the dataframe
summit_bed_files$RMP <- matches
summit_bed_files$RMP[grepl('TGIRT',summit_bed_files$summit.bed)] <- 'TRMT61A'
summit_bed_files$RMP[grepl('KIAA1429',summit_bed_files$summit.bed)] <- 'KIAA1429' #VIRMA=KIAA1429;
summit_bed_files$RMP[grepl('HAKAI',summit_bed_files$summit.bed)] <- 'HAKAI' #HAKAI=CBLL1
summit_bed_files$RMP[grepl('SETD2',summit_bed_files$summit.bed)] <- 'SETD2' # not RMP
summit_bed_files$RMP[grepl('DNMT2',summit_bed_files$summit.bed)] <- 'DNMT2' # not RMP
summit_bed_files$RMP[grepl('ontrol|GFP|m6A_rep|METTL14WT',summit_bed_files$summit.bed)] <- 'WT'
summit_bed_files[is.na(summit_bed_files$RMP),]
# summit_bed_files <- summit_bed_files[!grepl('GSM2324292|GSM2324294|GSM2324302|GSM2324298|GSM2324310|GSM2324306',summit_bed_files$summit.bed),]
summit_bed_files <- summit_bed_files[!grepl('GSM2324292|GSM2324294',summit_bed_files$summit.bed),] # remove the files of NOMO-1 cell line
# GSM2324291	NOMO-1 PBS input
# GSM2324292	NOMO-1 PBS m6A-seq
# GSM2324293	NOMO-1 R-2HG input
# GSM2324294	NOMO-1 R-2HG m6A-seq
# summit_bed_files$RMP[is.na(summit_bed_files$RMP)] <- c('FTO','WT','FTO','WT')
summit_bed_files$RMP[summit_bed_files$rank==33] <- c('WT','WT','FTO','FTO', 'WT','WT','FTO','FTO')
flatten_condition_IP_list <- function(condition_list) {
  unlist(unique(unlist(condition_list)))
}
conditions <- flatten_condition_IP_list(condition_IP_list)
conditions <- conditions[conditions!='WT']
conditions <- c(conditions, 'WT')
conditions_pattern <- paste(conditions, collapse = "|")
matches <- str_extract(summit_bed_files$summit.bed, conditions_pattern)
# Add the matches as a new column in the dataframe
summit_bed_files$condition <- matches
summit_bed_files$condition[grepl('TGIRT',summit_bed_files$summit.bed)] <- 'TRMT61AOE'
summit_bed_files$condition[grepl('siNSUN2',summit_bed_files$summit.bed)] <- 'NSUN2KD'
summit_bed_files$condition[grepl('siMETTL3|shMETTL3',summit_bed_files$summit.bed)] <- 'METTL3KD'
summit_bed_files$condition[grepl('siMETTL14|shMETTL14',summit_bed_files$summit.bed)] <- 'METTL14KD'
summit_bed_files$condition[grepl('siMETTL3andsiMETTL14',summit_bed_files$summit.bed)] <- 'METTL3_14KD'

summit_bed_files$condition[grepl('siWTAP|shWTAP',summit_bed_files$summit.bed)] <- 'WTAPKD'
summit_bed_files$condition[grepl('FTOIB',summit_bed_files$summit.bed)] <- 'FTO.2IB' 
summit_bed_files$condition[grepl('siKIAA1429',summit_bed_files$summit.bed)] <- 'KIAA1429KD'#VIRMA=KIAA1429;
summit_bed_files$condition[grepl('SETD2_2KD',summit_bed_files$summit.bed)] <- 'SETD2_2KD' # not RMP condition
summit_bed_files$condition[grepl('DNMT2',summit_bed_files$summit.bed)] <- 'DNMT2PD' # not RMP condition
summit_bed_files$condition[grepl('ontrol|GFP|m6A_rep|METTL14WT',summit_bed_files$summit.bed)] <- 'WT'
summit_bed_files$condition[grepl('YTHDF3KD_YTHDF',summit_bed_files$summit.bed)] <- c('YTHDF3KD_YTHDF1PD', 'YTHDF3KD_YTHDF1PD', 'YTHDF3KD_YTHDF2PD', 'YTHDF3KD_YTHDF2PD')

summit_bed_files[is.na(summit_bed_files$condition),]
summit_bed_files$condition[summit_bed_files$rank==33] <- c('WT','WT','FTOKD','FTOKD','WT','WT','FTOOE','FTOOE')
# summit_bed_files$condition[grepl('METTL3KO_N_',summit_bed_files$summit.bed)] <- 'METTL3KO_N' 
# summit_bed_files$condition[grepl('METTL3KO_NL_',summit_bed_files$summit.bed)] <- 'METTL3KO_NL' 
# summit_bed_files$condition[grepl('METTL3KO_TL_',summit_bed_files$summit.bed)] <- 'METTL3KO_TL' 
summit_bed_files <- summit_bed_files[!grepl('GSE132306',summit_bed_files$GSE),] # rm GSE132306(n = 48) whose result is the same with GSE120024(n = 47)
# summit_bed_files <- summit_bed_files[!grepl('GSE130172',summit_bed_files$GSE),]
summit_bed_files <- summit_bed_files[!grepl('GSE171497',summit_bed_files$GSE),] # rm GSE171497(n = 107) whose result is not good for now
setdiff(names(target_results), unique(summit_bed_files$GSE))#[1] "GSE125046" "GSE109183" "GSE122961" "GSE122948" "GSE128699"
# with RIP/miCLIP without Input, so no summit.bed
summit_bed_files$path = file.path(data_directory, summit_bed_files$GSE, 'metaPlotR', summit_bed_files$summit.bed)
write.table(summit_bed_files, quote = F, row.names = F, sep = '\t', col.names = F,
            file = paste0('/data/rluo4/EpiTrans/RMDatasets/summit_bed_files.txt'))
# http://localhost:8888/edit/Shell/UTH_RPMfunc/RPMdb-bw.sh
# original <- c('GSE97419', 'GSE55572', 'GSE87190', 'GSE103497', 'GSE141994', 'GSE145924')
# 'GSE97419', 'GSE55572', 'GSE87190', 'GSE76414',
# 'GSE103497', 'GSE141994', 'GSE145924', 'GSE124509', 
# 'GSE210867',  'GSE207643', 'GSE171497', 'GSE198643'
# paste0(original, collapse =  "' | '")
# 'GSE97419' | 'GSE55572' | 'GSE87190' | 'GSE103497' | 'GSE141994'| 'GSE145924' | 'GSE207643' | 'GSE198643'
summit_bed_files$newpath <- file.path('/data/rluo4/RPMfunc/GEO_metaPlotR',summit_bed_files$GSE,
                                      summit_bed_files$summit.bed)
# see Guitar.R and RM.R
# for (n in 1:nrow(RMPdeg_datasets)) {# 1, 7, 13, 19, 25, 28, 43, 45, 48, c(95, 108)){#
###########################################
# 5) Summary of RNA moditication results  #
###########################################
# peak_summary <- do.call(rbind, lapply(filtered_dep, function(x) x[[1]]))
# Function to combine nested lists of data frames
# combine_filtered_dep <- function(filtered_dep) {
#   # Loop over the outer list (GSE entries)
#   do.call(rbind, lapply(filtered_dep, function(gse_list) {
#     # For each GSE entry, loop over the nested list and extract the data frames
#     do.call(rbind, lapply(gse_list, function(df) {
#       # You can add a column to keep track of the comparison name if needed
#       df$comparison <- names(gse_list)[1]  # Or use the relevant identifier
#       return(df)
#     }))
#   }))
# }
combine_filtered_dep <- function(filtered_dep_list, dataset_summary) {
  # Initialize an empty list to store results
  summary_list <- list()
  # Iterate over each dataset in the list
  for (dataset_name in names(filtered_dep_list)) {
    target_res <- target_results[[dataset_name]]
    length(target_res)
    if( length(target_res) ==0 ){
      print(paste0(dataset_name, ": no data in target_results !"))
      next;
    }
    # resM <- filtered_dem_list[[dataset_name]]
    resP <- filtered_dep[[dataset_name]]
    # table(target_res)
    cl <- unique(str_split( names(target_res), '_', simplify = TRUE)[, 1])
    for (contrast in names(target_res)) {
      # contrast <- names(target_res)[1]
      c <- unique(str_split(contrast, '_', simplify = TRUE)[, 2])
      if( length(cl) >1 | dataset_name=='GSE207643'){
        c <- contrast
      }
      condition_name <- paste0(contrast, '_vs_WT')
      print(condition_name)
      target <- target_res[[contrast]]
      # DEM <- resM[[paste0(condition_name, '_DESeq2')]]
      DEP <- resP[[condition_name]]
      
      extracted_data <- DEP[-1,1:6]
      if(dataset_name=='GSE122961'){
        extracted_data <- DEP[, c(4,1:3,5:6)]
      }
      colnames(extracted_data) <- c('Peak_id', 'chr', 'start', 'end', 'strand', 'Stat')
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
# Apply the function to combine all dataframes
peak_summary <- combine_filtered_dep(filtered_dep_list = filtered_dep, dataset_summary = final_datasets)
rownames(peak_summary) <- 1:nrow(peak_summary)
unique(RMP_update$RNA_modifications)
length(unique(DEM_summary$Condition))
length(unique(peak_summary$Condition))
setdiff(unique(peak_summary$Condition), unique(DEM_summary$Condition))
cell_lines <- str_split(peak_summary$Cell_line, '_', simplify = T)[,1]
cell_lines <- str_split(cell_lines, '[.]', simplify = T)[,1]
cell_lines <- tolower(cell_lines)
sort(unique(cell_lines)) # 44 --> 42
dim(peak_summary) # [1] 2972866      14
identifiers <- paste(peak_summary$Peak_id, peak_summary$Dataset, peak_summary$Condition)
length(unique(identifiers)) #[1] 506055 (DEM_summary) --> [1] 2972866 (peak_summary)
# peak_summary$identifiers <- identifiers
# index <- duplicated(identifiers)
# dups <- peak_summary[index, ]
# unique(dups$Dataset)
table(peak_summary$Cell_line)
peak_summary$Cell_line <- gsub('HeLa', 'Hela', peak_summary$Cell_line)
peak_summary$Cell_line <- gsub('HaCAT', 'HaCaT', peak_summary$Cell_line)
index <- peak_summary$Condition == 'Hela.m6A_ALKBH3KO_vs_WT'
peak_summary$Modification_type[index] <- 'm6A'
table(peak_summary$Modification_type)
CL_RMT <- peak_summary
CL_RMT$Cell_line<- str_split(peak_summary$Cell_line, '_', simplify = T)[,1]
CL_RMT <- as.data.frame(table(CL_RMT$Cell_line, CL_RMT$Modification_type))
CL_RMT <- CL_RMT[CL_RMT$Freq!=0,]
ggplot(melt(CL_RMT), aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  labs(x = "Variable 1", y = "Variable 2", fill = "Count")
# List of dataframes
df_list <- list(DEG_summary,DEM_summary)
# Using Reduce and merge for full outer join
bulk_result <- Reduce(function(x, y) merge(x, y, all = TRUE), df_list)
length(unique(bulk_result$GeneSymbol)) # 57274
unique(bulk_result$Modification_type)
# [1] "m6A"  "m7G"  "m5C"  "m1A"  "m6Am" "Psi"  "ac4C"
# [8] "m5U" 
unique(bulk_result$Enzyme) # 39 - 5 = 34
Enzymes <- unique(str_split(bulk_result$Enzyme, '_1|_2', simplify = T)[,1])
Enzymes <- Enzymes[ ! Enzymes %in% c('FTO.1', 'FTO.2')] # 35
length(unique(bulk_result$GeneSymbol))
unique(bulk_result$Technique) 
table(bulk_result$Treatment)
write.table(bulk_result, "/data/rluo4/RPMfunc/Output/summary/bulk_result.txt", row.names = FALSE, sep = '\t', quote = F)
length(unique(DEG_summary$GeneSymbol))
length(unique(DEM_summary$GeneSymbol))
length(unique(DEM_summary$Peak_id))
inter_genes <- (intersect(unique(DEG_summary$GeneSymbol),unique(DEM_summary$GeneSymbol)))
table(inter_genes %in% c('ANXA1'))
# target = intersect(unique(DEG$geneSymbol), unique(DEM$geneSymbol))
# http://localhost:8585/view/RPMfunc/Output/GEO/GSE94613/MOLM13_METTL3KD/PeakVis/FPR1_peak.png
# http://localhost:8585/view/RPMfunc/Output/GEO/GSE144984/NOMO1_ALKBH5OE/PeakVis/ANXA1_peak.png
# http://localhost:8585/view/RPMfunc/Output/GEO/GSE122803/MEL624.1_PCIF1KO/PeakVis/NECTIN2_peak.png
# http://localhost:8585/view/RPMfunc/Output/GEO/GSE122803/MEL624.1_PCIF1KO/PeakVis/TGFB1_peak.png
# n %in% c(12, 19, 56) DEG_res need to be modified
bulk_summary <- bulk_result[, c(4:9)]
bulk_summary <- bulk_summary[!duplicated(paste(bulk_summary$Condition, bulk_summary$Dataset)),]
bulk_summary$Folder <- file.path('/data/rluo4/RPMfunc/Output/GEO', bulk_summary$Dataset,gsub('_vs_WT', '',bulk_summary$Condition))
write.table(bulk_summary, "/data/rluo4/RPMfunc/Output/summary/bulk_summary.txt", row.names = FALSE, sep = '\t', quote = F)
# write.table(CCI_summary_All[,c(7,1:6)], '/home/lorihan/lrh/All/Output/Index_LR_All.txt',#paste0(TissueType,'_Index_CellChat.txt'),
#             sep='\t',quote=F,row.name=F,col.name=T)
CCI_summary_all <- read.table('/data/rluo4/summary/Index_LR_All.txt', fill = T, header = T)
# intersect(target, CCI_summary_all$Receptor)
inter <- intersect(DEG_summary$GeneSymbol, DEM_summary$GeneSymbol)
# inter <- bulk_result[bulk_result$GeneSymbol %in% inter, ]
inter <- DEM_summary[DEM_summary$GeneSymbol %in% inter, ]
inter <- inter[order(inter$Dataset),]
# # View(inter[inter$GeneSymbol %in% c('FN1','SDC1'),])
# # View(inter[inter$GeneSymbol %in% c('LGALS9','CD44','CD45'),])
# # View(inter[inter$GeneSymbol %in% c('TGFB1','TGFBR1'),])
# View(inter[inter$GeneSymbol %in% c('SAA1','COL6A3','ACTA2','FN1','SDC1','LGALS9','CD44','CD45','NECTIN2','NECTIN3','TIGIT','FPR1', 'ANXA1','TGFB1','TGFBR1'),])
# View(DEG_summary[DEG_summary$GeneSymbol %in% c('SAA1','COL6A3','ACTA2','FN1','SDC1','LGALS9','CD44','CD45','NECTIN2','NECTIN3','TIGIT','FPR1', 'ANXA1','TGFB1','TGFBR1'),])
# 
# View(bulk_result[grepl('NECTIN2',bulk_result$GeneSymbol),])
# View(bulk_result[grepl('TIGIT',bulk_result$GeneSymbol),])
# 
# View(bulk_result[grepl('FPR',bulk_result$GeneSymbol),])
# View(bulk_result[bulk_result$GeneSymbol == 'ANXA1',])
# View(bulk_result[bulk_result$GeneSymbol == 'TIGIT',])
# 
# # View(DEM_summary[DEM_summary$Enzyme=='VIR',])
# # View(DEM_summary[DEM_summary$Enzyme=='METT',])
##################################################
# 6) PTM of RNA modification protein regulators  #
##################################################
# PTM_info <- read_excel('/data/rluo4/EpiTrans/DataCollection/RBP-PTM-2024.5.31.xlsx', sheet = 8)
PTM_info <- read_excel('/data/rluo4/RPMfunc/Output/summary/ALL_PTM0629.xlsx', sheet = 1)
PTM_info$Disease <- gsub(' MB', ' Medulloblastoma', PTM_info$Disease)
table(PTM_info$Disease == ' ')
PTM_info$Disease[is.na(PTM_info$Disease)] <- 'unknown'
PTM_info$`Upstream regularory proteins`[is.na(PTM_info$`Upstream regularory proteins`)] <- 'unknown'
PTM_info$`Interacting Proteins`[is.na(PTM_info$`Interacting Proteins`)] <- 'unknown'
table(PTM_info$`Data Sources`)
PTM_info$Link_db <- '-'
PTM_info$Link_db[PTM_info$`Data Sources`=='Cell'] <- 'https://pdc.cancer.gov/pdc/cptac-pancancer'
PTM_info$Link_db[PTM_info$`Data Sources`=='dbPTM'] <- 'https://awi.cuhk.edu.cn/dbPTM/'
PTM_info$Link_db[PTM_info$`Data Sources`=='PTM-RBP-ATLAS'] <- 'http://ptm-rbp-atlas.igb.uci.edu/'
PTM_info$Link_PMID <- paste0('https://pubmed.ncbi.nlm.nih.gov/', PTM_info$PMID, '/' )
write_xlsx(x = PTM_info, path = paste0('/data/rluo4/RPMfunc/Output/summary/ALL_PTM0703.xlsx'))

sites <- paste(PTM_info$Regulators, PTM_info$PTMs, PTM_info$`Identified residues`)
Freq <- data.frame(table(PTM_info$PTMs, sites))
Freq_database <- Freq[Freq$Freq !=0,]
unique(PTM_info$Regulators)
unique(PTM_info$PTMs)
PTM_PMID <- read_excel('/data/rluo4/RPMfunc/Output/summary/PTM-Table23.xlsx', sheet = 1)
unique(PTM_PMID$PTMs)
PTM_PMID$PTMs <- gsub('SUMOylation', 'Sumoylation', PTM_PMID$PTMs)
index = PTM_PMID$`Identified residues`=="–"
PTM_PMID$`Identified residues`[index] <- paste('no.',rownames(PTM_PMID)[index],'-')
PTMs <- unique(c(PTM_info$PTMs,PTM_PMID$PTMs))
PTMs
sites <- paste(PTM_PMID$Regulators, PTM_PMID$PTMs, PTM_PMID$`Identified residues`)
Freq <- data.frame(table(PTM_PMID$PTMs, sites))
Freq_PMID <- Freq[Freq$Freq !=0,]
proteomics_stats <- read_excel('/data/rluo4/EpiTrans/DataCollection/diff_corgene_stat.xlsx',sheet = 1)
sum(proteomics_stats$diff_gene)

proteomics <- read_excel('/data/rluo4/EpiTrans/DataCollection/Proteomics-summary.xlsx', sheet = 1)
proteomics$omics <- rep('proteomics', nrow(proteomics))
phosphoproteomics  <- read_excel('/data/rluo4/EpiTrans/DataCollection/Proteomics-summary.xlsx', sheet = 2)
phosphoproteomics$omics <- rep('phosphoproteomics', nrow(phosphoproteomics))
acetylomics   <- read_excel('/data/rluo4/EpiTrans/DataCollection/Proteomics-summary.xlsx', sheet = 3)
acetylomics$omics <- rep('acetylomics', nrow(acetylomics))
bulkRNA   <- read_excel('/data/rluo4/EpiTrans/DataCollection/Proteomics-summary.xlsx', sheet = 4)
bulkRNA$omics <- rep('transcriptomics', nrow(bulkRNA))
bulkDNA   <- read_excel('/data/rluo4/EpiTrans/DataCollection/Proteomics-summary.xlsx', sheet = 5)
bulkDNA$Normal <- rep(0, nrow(bulkDNA))
bulkDNA$omics <- rep('genomics', nrow(bulkDNA))
glycosylomics  <- read_excel('/data/rluo4/EpiTrans/DataCollection/Proteomics-summary.xlsx', sheet = 6)
glycosylomics$omics <- rep('glycosylomics', nrow(glycosylomics))
ubiquitylomics  <- read_excel('/data/rluo4/EpiTrans/DataCollection/Proteomics-summary.xlsx', sheet = 7)
ubiquitylomics$omics <- rep('ubiquitylomics', nrow(ubiquitylomics))

df_list <- list(proteomics,phosphoproteomics,acetylomics,bulkRNA,bulkDNA, glycosylomics, ubiquitylomics)
# Using Reduce and merge for full outer join
Proteogenomic_datasets <- Reduce(function(x, y) merge(x, y, all = TRUE), df_list)
length(unique(Proteogenomic_datasets$Cancertype)) #11
sum(Proteogenomic_datasets$Tumor) #4642 -->5838
cancercellcli <- read_xlsx("/data/rluo4/EpiTrans/PTM/CPTAC/inputdata/cancer_cell_clinical.xlsx",sheet = 2)
load("/data/rluo4/EpiTrans/PTM/CPTAC/inputdata/acetylist.rdata")
acety_pro <- colnames(acetylist$BRCA)
acety_pro.uniq <- str_split(unique(acety_pro), ' ', simplify = T)[,1]
length(unique(acety_pro.uniq)) # 2241
# We identified 9446 cancer mutations that coincide with the position of categorized PTM modified sites in 1727 RBPs.
# We also characterized mutations surrounding PTM sites, as the local structure of proteins is known to be critical for the recognition of deposition enzymes to control PTM stoichiometry(62). 
# From this analysis we identified 100 002 cancer mutations that fall within 10 residues of PTM modified sites in 2285 RBPs
##################################################
# 7) #        Single-cell NMF analysis           #
##################################################
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(ggpubr)
library(ggsignif)
# library(ggstatsplot)
library(ggplot2)
library(ggalluvial)
library(RColorBrewer)
library(ggthemes)
library(maftools)
# out_dir <-  '/data/rluo4/RPMfunc/Output/scRNA/PCTanno_SNV'
data_path = '/data/rluo4/All/Output'
# load('/data/rluo4/All/Output/pheno_all.rds')#Cervix_Epi
# cell_summary <- read.csv('/data/rluo4/All/Output/cell_summary.txt',sep = "\t",header = T)
tissue_summary <- read.table('/data/rluo4/RPMfunc/PCTanno_pdata/tissue_summary.txt')
# disco_original <- readRDS(paste0('/data/rluo4/RPMfunc/disco_pdata/disco_all_meta.rds'))
disco <- readRDS(paste0('/data/rluo4/RPMfunc/disco_pdata/disco_adjustmeta_0827.rds'))#disco_adjustmeta_0720.rds'))
disco_summary <- disco[, c(1,2,4,7)]
disco_summary <- disco_summary[! duplicated(paste(disco_summary$patient_id, disco_summary$tissue, 
                                                  disco_summary$project_id)),]
colnames(disco_summary)
PCTanno <- tissue_summary[, c(10,2:4)]
colnames(PCTanno) <- colnames(disco_summary)
intersect(disco_summary$project_id, PCTanno$project_id)
# [1] "GSE134355" "GSE161529" "GSE148673" "GSE164898" "GSE195861"
# [6] "GSE164241" "GSE134520" "GSE136103" "GSE112271" "GSE208653"
# [11] "GSE172357" "GSE167297" "GSE154778" "GSE137829" "GSE172577"
# [16] "GSE144236" "GSE130973" "GSE193304"
sc_summary <- rbind(disco_summary[!disco_summary$project_id %in% PCTanno$project_id, ], PCTanno)
length(unique(sc_summary$project_id)) #396 --> 379
sc_summary$omics <- 'scRNA'
####################################################
# 7.1. Organization of scRNA datasets: disco + PCTanno
# Skip below since it's time-consuming:
####################################################
# # 3.1 Adjust sub_type in sc_majortype of diso & PCTanno
# disco <- readRDS(paste0('/data/rluo4/RPMfunc/disco_pdata/disco_adjustmeta_0827.rds'))#0809.rds'))
# disco_majortype <- readRDS('/data/rluo4/RPMfunc/Output/scRNA/disco_all_adjustmeta_0827.rds')#0809.rds') # made by DiscoDatasets.R:
# disco_pheno <- as.data.frame(table(disco_majortype$sample_type, disco_majortype$disease))
# disco_pheno <- disco_pheno[disco_pheno$Freq!=0,]
# # sample_type inconsistent with disease: GSE168453 (Critical/Moderate -- Healthy)
# # sample.index <- c(grep("ormal", disco_pheno$Var1),grep("ealthy", disco_pheno$Var1),
# #                   grep("control", disco_pheno$Var1))
# # View(disco_pheno[sample.index,])
# # disco changes of disease column
# sample.index <- grepl("ormal|ealthy|control", disco_majortype$sample_type)
# table(disco_majortype$disease[sample.index])
# disco_majortype$disease[sample.index] <- 'Healthy'
# covid.index <- disco_majortype$project_id == 'GSE168453'
# table(disco_majortype$disease[covid.index])
# table(disco_majortype$disease[covid.index], disco_majortype$sample_type[covid.index])
# disco_majortype$disease[covid.index] <- disco_majortype$sample_type[covid.index]# critical, moderate and severe are COVID negative with NIH clinical symptomes
# GSE168453_meta <- read.csv('/data/rluo4/EpiTrans/DataCollection/GSE168453_covars.csv')
# # GSE168453$NIH_clinical <- GSE168453_meta$NIH_clinical[match(GSE168453$patient_id, GSE168453_meta$donor)]
# # GSE168453$COVID_status <- GSE168453_meta$COVID_status[match(GSE168453$patient_id, GSE168453_meta$donor)]
# # GSE168453$pulmonary_infection <- GSE168453_meta$pulmonary_infection[match(GSE168453$patient_id, GSE168453_meta$donor)]
# # GSE168453_healty <- read.table('/data/rluo4/EpiTrans/DataCollection/COVID_healthy.txt')
# # GSE168453_healty$COVID_status <- GSE168453_meta$COVID_status[match(GSE168453_healty$V1, GSE168453_meta$donor)]
# # GSE168453_healty$NIH_clinical <- GSE168453_meta$NIH_clinical[match(GSE168453_healty$V1, GSE168453_meta$donor)]
# # GSE168453_healty$sample_type <- GSE168453$sample_type[match(GSE168453_healty$V1, GSE168453$patient_id)]
# # # write.table(GSE168453_disease$COVID_disease,file = '/data/rluo4/EpiTrans/DataCollection/COVID_disease.txt',
# # #             col.names = F, row.names = F, quote = F)
# # GSE168453_disease <- read.table('/data/rluo4/EpiTrans/DataCollection/COVID_disease.txt')
# # GSE168453_disease$COVID_status <- GSE168453_meta$COVID_status[match(GSE168453_disease$V1, GSE168453_meta$donor)]
# # GSE168453_disease$NIH_clinical <- GSE168453_meta$NIH_clinical[match(GSE168453_disease$V1, GSE168453_meta$donor)]
# # GSE168453_disease$sample_type <- GSE168453$sample_type[match(GSE168453_disease$V1, GSE168453$patient_id)]
# CT <- as.data.frame(table(disco_majortype$ct, disco_majortype$sub_type))
# CT <- CT[CT$Freq!=0,]
# # View(CT[duplicated(CT$Var1),])
# index1 <- grepl('eutrophil',disco_majortype$ct)
# disco_majortype$sub_type[index1] <- 'Neutrophils'
# index2 <- grepl('EMT',disco_majortype$ct)
# disco_majortype$sub_type[index2] <- 'Tissue_Specific_Cells'
# index3 <- grepl('NKT',disco_majortype$ct)
# disco_majortype$sub_type[index3] <- 'Tcells'
# index4 <- grepl(' EC',disco_majortype$ct) # from BM_MAST_meta_p1
# disco_majortype$sub_type[index4] <- 'Endothelial_Cells'
# index5 <- grepl('Granulocyte|Myelocyte',disco_majortype$ct)
# disco_majortype$sub_type[index5] <- 'OtherImmunecells'
# # Granulocyte <- disco_majortype[disco_majortype$ct=='Granulocyte',]
# # table(Granulocyte$project_id, Granulocyte$sub_type)#GSE185381--BM
# # table(Granulocyte$project_id, Granulocyte$tissue)#GSE185381--BM
# # HSC <- disco_majortype[disco_majortype$ct=='HSC',]
# # table(HSC$project_id, HSC$sub_type)#GSE181989--BM
# # table(HSC$project_id, HSC$tissue)#GSE181989--BM
# # Placental mesenchymal stem/stromal cells (PMSCs) are fibroblast-like, self-renewing, and multipotent cells that reside in the chorionic plate, umbilical cord, and amnion of the placenta
# index6 <- grepl('HSC',disco_majortype$ct) & disco_majortype$tissue %in% c('blood','bone marrow', 'Buffy coat','PBMC','thymus')
# disco_majortype$sub_type[index6] <- 'Stem_Cells'
# index7 <- grepl('myeloid',disco_majortype$ct)
# disco_majortype$sub_type[index7] <- 'Myeloids'
# index8 <- grepl(' T/NK',disco_majortype$ct)
# disco_majortype$sub_type[index8] <- 'NKcells'
# index9 <- grepl('Mast',disco_majortype$ct)
# disco_majortype$sub_type[index9] <- 'MastCell'
# index10 <- grepl('Oligodendrocyte',disco_majortype$ct)
# disco_majortype$sub_type[index10] <- 'Neuro_Cells'
# index11 <- grepl('epithelia',disco_majortype$ct)
# disco_majortype$sub_type[index11] <- 'Tissue_Specific_Cells'
# # saveRDS(disco_majortype,  file=paste0('/data/rluo4/RPMfunc/Output/scRNA/disco_majortype_0828.rds'))
# 
# PCT_majortype <- readRDS('/data/rluo4/RPMfunc/Output/scRNA/PCT_all_meta.rds') # made by Haixia
# table(unique(PCT_majortype$sample_name) %in% tissue_summary$Cohort)
# PCT_majortype$project_id <- tissue_summary$Dataset[match(PCT_majortype$sample_name, tissue_summary$Cohort)]
# celltype_summary <- read_excel('/data/rluo4/All/Output/table_xlsx/TableS1.xlsx')
# table(PCT_majortype$Cell_Type %in% celltype_summary$Minor.Celltype)
# setdiff(PCT_majortype$Cell_Type, celltype_summary$Minor.Celltype)#  "MSC.ADIPO" "MSC.MVA"   "MSC.SEC"   "MSC.PVA"
# index <- match(PCT_majortype$Cell_Type, celltype_summary$Minor.Celltype)
# PCT_majortype$FullName <- celltype_summary$FullName[index]
# table(is.na(PCT_majortype$FullName))
# unique(PCT_majortype$Cell_Type[is.na(PCT_majortype$FullName)])
# PCT_majortype$ct[ is.na(PCT_majortype$FullName) ] <- 'Mesenchymal stem cells'
# PCT_majortype$ct[! is.na(PCT_majortype$FullName) ] <- PCT_majortype$FullName[! is.na(PCT_majortype$FullName) ]
# unique(PCT_majortype$ct)
# unique(PCT_majortype$sub_type)
# unique(disco_majortype$sub_type)
# index1 <- PCT_majortype$Major_type=='Epithelial'
# PCT_majortype$sub_type[index1] <- 'Tissue_Specific_Cells'
# index2 <- PCT_majortype$Cell_Type=='MDSC'
# PCT_majortype$sub_type[index2] <- 'Myeloids'
# index3 <- PCT_majortype$Cell_Type=='STM'
# PCT_majortype$sub_type[index3] <- 'Stem_Cells'
# index4 <- grepl('Mast',PCT_majortype$ct)
# PCT_majortype$sub_type[index4] <- 'MastCell'
# index5 <- grepl('endritic',PCT_majortype$ct)
# PCT_majortype$sub_type[index5] <- 'OtherImmunecells'
# index6 <- grepl('uppfer',PCT_majortype$ct)
# PCT_majortype$sub_type[index6] <- 'Macrophages'
# index7 <- grepl('Megakaryocytes|Erythrocytes',PCT_majortype$ct)
# PCT_majortype$sub_type[index7] <- 'Blood_Cells'
# CT <- as.data.frame(table(PCT_majortype$ct, PCT_majortype$sub_type))
# CT <- CT[CT$Freq!=0,]
# # View(CT[duplicated(CT$Var1),])
# PCT_majortype$Major_type <- gsub('Epithelial','Tissue_Specific_Cells',
#                                  gsub('Stromal','Mesenchymal_Cells',
#                                       gsub('Immune','Immune_Cells',PCT_majortype$Major_type)))
# table(PCT_majortype$project_id %in% tissue_summary$Dataset)
# PCT_majortype$platform <- tissue_summary$Platform[match(PCT_majortype$project_id, tissue_summary$Dataset)]
# # saveRDS(PCT_majortype[, -10], file = '/data/rluo4/RPMfunc/Output/scRNA/PCT_all_adjustmeta_0715.rds')
# # write_rds() will save a big file.
# #Langerhans cell histiocytosis (LCH) cells are a type of dendritic cell that normally helps the body fight infection. Sometimes mutations (changes) develop in genes that control how dendritic cells function. These include mutations of the BRAF, MAP2K1, RAS, and ARAF genes. These mutations may cause too many LCH cells to grow and build up in certain parts of the body, where they can damage tissue or form lesions.
# colnames(disco_majortype)
# # disco_majortype <- disco_majortype[, c(1:4,5,7:10)]
# sc_majortype <- PCT_majortype %>% dplyr::select(barcode, project_id, patient_id = orig.ident, ct ,
#                                          tissue = Organ, platform, sample_type = Tissue, disease = Tissue, Major_type, sub_type)
# # sc_majortype$project_id <- PCTanno$project_id[match(sc_majortype$patient_id, PCTanno$patient_id)]
# sc_majortype <- rbind(sc_majortype, disco_majortype[! disco_majortype$project_id %in% sc_majortype$project_id,])
# unique(sc_majortype$sub_type) # 17
# unique(sc_majortype$project_id) # 340 --> 312 --> 296
# # sc_majortype$sub_type[sc_majortype$sub_type=='Fibroblast'] = 'Fibroblasts'
# # sc_majortype$sub_type[sc_majortype$sub_type=='Macrophage'] = 'Macrophages'
# table(sc_majortype$sub_type)
# table(sc_majortype$tissue)
# unique(sc_majortype$tissue) #79
# CT <- as.data.frame(table(sc_majortype$ct, sc_majortype$sub_type))
# CT <- CT[CT$Freq!=0,]
# # View(CT[duplicated(CT$Var1),])
# # 3.2 Adjust tissue and disease in sc_majortype of diso & PCTanno
# sc_majortype$tissue <-  gsub('Acute meyloid leukemia|Acute leukemia|Acute monocytic leukemia|Acute myeloid leukemia|Chronic lymphocytic leukemia|Chronic myeloid leukemia|Leukemia|Myelogenous leukemia', 'bone marrow',
#                               gsub('Breast and Brain metastases|Breast|BRCA', 'breast',
#                               gsub('Cervical', 'cervix',
#                               gsub('Cervical and Liver', 'cervix and liver',
#                               gsub('CCRCC|Human embryonic kidney|Normal kidney', 'kidney',
#                               gsub('Colorectal adenocarcinoma and SARS-CoV-2|COAD|Colon|Colorectum', 'colon',
#                               gsub('Diffuse large B-cell lymphoma|Erythroleukemia', 'blood',
#                               gsub('PDAC', 'pancreas',
#                               gsub('OV|HGSC', 'ovary',
#                               gsub('Myelogenous leukemia and Human embryonic microglia', 'bone marrow and brain',
#                               gsub('Embryo placenta', 'placenta',
#                               gsub('Endometrium|UCEC|Endometrial', 'endometrium',
#                               gsub('GBM|Glioblastoma|Giloma|Human embryonic microglia', 'brain',
#                               gsub('HGSC|OV', 'ovary',
#                               gsub('LSCC|LUAD|Embryonic lung fibroblasts', 'lung',
#                               gsub('HNSCC|Oral cavity|Head and Neck Squamous Carcinoma', 'head and neck',
#                               gsub('Skin cutaneous melanoma|Skin|Human keratinocytes|Foreskin|Melanoma', 'skin',
#                               gsub("Kaposi's sarcoma", 'intestine',
#                               gsub("Human Umbilical Vein Endothelial Cell|Human umblilical cord blood", 'blood',
#                               gsub('Nasopharyngeal carcinoma','nasopharynx',
#                               gsub('Uveal melanoma|Ocular melanoma', 'eye',
#                               gsub('T2D islets', 'pancreatic islet',
#                               gsub('Gastric', 'stomach',
#                               gsub('Esophageal squamous cell carcinoma', 'esophagus',
#                               gsub('Human embryonic stem cell|Human embryonic stem cells', 'endoderm',
#                                    sc_majortype$tissue,
#                                    )))))))))))))))))))))))))
# sc_majortype$tissue <- tolower(sc_majortype$tissue)
# table(sc_majortype$tissue) # 72 kinds of tissues in scRNA datasets
# sc_majortype$tissue <- gsub('head and neck', "mouth",sc_majortype$tissue)# all oral cavity --> head and neck --> mouth, sc_majortype$tissue <- gsub("oral cavity", 'mouth',sc_majortype$tissue)
# sc_majortype$tissue <- gsub("head", 'head and neck',sc_majortype$tissue) # sc_majortype$tissue <- gsub("teeth", 'gingiva',sc_majortype$tissue)
# 
# library(dplyr)
# # Create a named vector with abbreviations and their corresponding full names (PCTanno)
# abbreviations <- c(
#   ADJ = "Adjacent tumor/disease section",
#   IDC = "Invasive ductal carcinoma",
#   DCIS = "Ductal carcinoma in situ",
#   `BRCA1-mut` = "Precancerous lesion from BRCA1 mutation carriers",
#   CC = "Cervix cancer",
#   N_HPV = "HPV-infected normal cervix",
#   HSIL_HPV = "HPV-infected high-grade squamous intraepithelial lesions",
#   FAP = "Familial adenomatous polyposis",
#   CRC = "Colorectal cancer",
#   AD = "Adenomas",
#   SER = "Sessile serrated lesions",
#   `MSI-H` = "Microsatellite-high colorectal cancer",
#   MSS = "Microsatellite stable colorectal cancer",
#   AEH = "Atypical endometrial hyperplasia",
#   EEC = "Endometrioid cancer",
#   HGIN = "High-grade intraepithelial neoplasias",
#   LGIN = "High-grade intraepithelial neoplasias",
#   ESCC = "Esophageal squamous cell carcinoma",
#   HCC = "Hepatocellular carcinoma",
#   NAFLD = "Non-alcoholic fatty liver disease",
#   AAH = "Atypical adenomatous hyperplasia",
#   AIS = "Adenocarcinoma in situ",
#   MIAC = "Minimally invasive adenocarcinoma",
#   IAC = "Invasive lung adenocarcinoma",
#   LP = "Leukoplakia",
#   NEOLP = "Non-erosive oral lichen planus",
#   EOLP = "Erosive oral lichen planus",
#   OSCC = "Oral squamous cell carcinoma",
#   PDAC = 'Pancreatic ductal adenocarcinoma',
#   PanIN = 'Pancreatic Intraepithelial Neoplasia',
#   BPH = "Benign prostatic hyperplasia",
#   PRAD = "Prostate adenocarcinoma",
#   AK = "Actinic keratosis",
#   SCCIS = "Squamous cell carcinoma in situ",
#   cSCC = "Cutaneous squamous cell carcinoma",
#   CAG = "Chronic atrophic gastritis",
#   `CAG with IM` = "Chronic atrophic gastritis with intestinal metaplasia",
#   CSG = "Chronic superficial gastritis",
#   GC = "Gastric cancer",
#   SIM = "Severe intestinal metaplasia",
#   WIM = "Wild intestinal metaplasia",
#   PTC = "Papillary thyroid cancer",
#   ATC = "Anaplastic thyroid cancer",
#   HT = "Hashimoto's thyroiditis"
# )
# # Replace the abbreviations with full names
# sc_majortype$disease <- recode(sc_majortype$disease, !!!abbreviations)
# sc_majortype$disease <- gsub("AdjNorm", 'Adjacent tumor/disease section',
#   gsub('AML', 'Acute myeloid leukemia', gsub("’s", "'s", sc_majortype$disease)) )
# sc_majortype$disease[sc_majortype$disease == "Infection"] <- 'COVID-19'
# # sc_majortype$disease[sc_majortype$disease == "severe"] <- 'Severe COVID-19'
# # sc_majortype$disease[sc_majortype$disease == "critical"] <- 'Critical COVID-19'
# # sc_majortype$disease[sc_majortype$disease == "moderate"] <- 'Moderate COVID-19'
# sc_majortype$disease[sc_majortype$project_id=='GSE145927'] <- 'acute antibody-mediated rejection' # kidney transplant
# sc_majortype$disease <- tolower(sc_majortype$disease)
# # Check the result
# head(sc_majortype$disease)
# sort(unique(sc_majortype$disease)) # ~ 154 kinds --> 142 -->173
# # saveRDS(sc_majortype, file = '/data/rluo4/RPMfunc/Output/scRNA/sc_majortype_0725.rds')#sc_majortype.rds')
# # saveRDS(sc_majortype, file = '/data/rluo4/RPMfunc/Output/scRNA/sc_majortype_0808.rds')# BM_MAST_p1
# # saveRDS(sc_majortype, file = '/data/rluo4/RPMfunc/Output/scRNA/sc_majortype_0809.rds')# BM_MAST_p1+p2+p3(GSE181294)
# table(sc_majortype$project_id[is.na(sc_majortype$patient_id)])
# # GSE181989 <- subset(sc_majortype, project_id == 'GSE181989')
# # GSM <- str_split(rownames(GSE181989), '--', simplify = T)[,2]
# index <- is.na(sc_majortype$patient_id)
# sc_majortype$patient_id[index] <- str_split(rownames(sc_majortype[index, ]), '--', simplify = T)[,2]
# table(sc_majortype$project_id[is.na(sc_majortype$patient_id)])
# saveRDS(sc_majortype, file = '/data/rluo4/RPMfunc/Output/scRNA/sc_majortype_0828.rds')#3.rds')# BM_MAST_p1+p2+p3(GSE181294) + last2 + wrong3 + last1

##################################################
# 7.2. Organization of merge nmf metadata
##################################################
# sc_majortype <- readRDS(file = '/data/rluo4/RPMfunc/Output/scRNA/sc_majortype_0725.rds') # 340 --> 312 datasets
sc_majortype <- readRDS(file = '/data/rluo4/RPMfunc/Output/scRNA/sc_majortype_0828.rds') # 340 --> 312 datasets
merge_nmf_meta <- readRDS(file = '/data/rluo4/EpiTrans/DataCollection/merge_nmf_meta_0925.rds')
unique(merge_nmf_meta$project_id) # 260 --> 223 -->263 -->268-->270
unique(merge_nmf_meta$tissue)#70 --> 64
unique(merge_nmf_meta$disease)#171
# table(is.na(merge_nmf_meta$patient_id))
length(unique(merge_nmf_meta$patient_id)) #3578 -->3586
unique(sc_majortype$project_id) #327
setdiff(unique(merge_nmf_meta$project_id), unique(sc_majortype$project_id))
# [1] "E-MTAB-8495" "GSE159929"   "GSE126030"   "GSE121080"  
# [5] "GSE150430"
# merge_nmf_meta$project_id <- gsub('https://ngdc.cncb.ac.cn/bioproject/browse/PRJCA001063',
#                                   'PRJCA001063', merge_nmf_meta$project_id)
diff_DT <- setdiff(unique(merge_nmf_meta$project_id), unique(sc_majortype$project_id))
# "PRJNA001063"
table(unique(merge_nmf_meta$patient_id) %in% unique(sc_majortype$patient_id))
setdiff(unique(merge_nmf_meta$patient_id), unique(sc_majortype$patient_id))
diff_pt <- setdiff(unique(merge_nmf_meta$patient_id), unique(sc_majortype$patient_id))#52
##################################################
# change 1:
##################################################
table(merge_nmf_meta$disease[merge_nmf_meta$patient_id %in% diff_pt])  
# aplastic anemia         healthy     oral cancer 
# 29088          197466            8686 
index <- merge_nmf_meta$disease %in% c('aplastic anemia', 'oral cancer')
# View(merge_nmf_meta[index, ])
unique(merge_nmf_meta$project_id[index]) # [1] "GSE181989" "GSE172577"
merge_nmf_meta$barcode[index] <- rownames(merge_nmf_meta)[index]#str_split(rownames(merge_nmf_meta)[index], '[...]', simplify = T)[,1]
merge_nmf_meta$patient_id[index] <- str_split(merge_nmf_meta$barcode[index], '--', simplify = T)[,2]
# table(merge_nmf_meta$tissue[merge_nmf_meta$project_id %in% diff_pt]) # healthy 
# setdiff(unique(merge_nmf_meta$disease), unique(sc_majortype$disease))
unique(sc_majortype$patient_id[sc_majortype$project_id=='GSE181989'])

##################################################
# change 2:
##################################################
index = merge_nmf_meta$project_id %in% diff_DT
unique(merge_nmf_meta$project_id[index])
patient_id <- str_split(merge_nmf_meta$barcode[index], '--', simplify = T)[,2]
table(patient_id == '')
table(patient_id != merge_nmf_meta$patient_id[index])
index1 <- patient_id == ''
table(merge_nmf_meta$project_id[index][index1])
# E-MTAB-8495   GSE126030   GSE150430
# 13041        9865       28910
diff_dt <- unique(merge_nmf_meta$project_id[index][index1])
table(sc_majortype$project_id %in% diff_dt)#c('E-MTAB-8495', 'GSE126030', 'GSE150430'))
# View(merge_nmf_meta[index,][index1,])
table(merge_nmf_meta$patient_id[index][index1])
table(merge_nmf_meta$orig.ident[index][index1])
patient_id[index1] <- merge_nmf_meta$patient_id[index][index1]
index2 <-  patient_id != merge_nmf_meta$patient_id[index]
# View(merge_nmf_meta[index,][index2,])
merge_nmf_meta$patient_id[index][index2] <- patient_id[index2]
barcode <- str_split(merge_nmf_meta$barcode[index], '--', simplify = T)[,1]
table(merge_nmf_meta$patient_id[index] == patient_id)
merge_nmf_meta$barcode[index] <- paste0(barcode, '--', patient_id)

##################################################
# change 3:
##################################################
diff_pt <- setdiff(unique(merge_nmf_meta$patient_id), unique(sc_majortype$patient_id))#52 --> 64
index <- merge_nmf_meta$patient_id %in% diff_pt
# View(merge_nmf_meta[index, ])
patient_id <- str_split(merge_nmf_meta$barcode[index], '--', simplify = T)[,2]
table(patient_id == '')
table(patient_id != merge_nmf_meta$patient_id[index])
index1 <- patient_id != merge_nmf_meta$patient_id[index]
table(merge_nmf_meta$project_id[index][index1])
# GSE156625 GSE161529 
# 1278     13571 
table(merge_nmf_meta$patient_id[index][index1])
table(patient_id)
diff_dt <- unique(merge_nmf_meta$project_id[index][index1])
table(sc_majortype$project_id %in% diff_dt)
unique(sc_majortype$patient_id[sc_majortype$project_id %in% diff_dt])
unique(merge_nmf_meta$patient_id[merge_nmf_meta$project_id %in% diff_dt])
merge_nmf_meta$patient_id[index][index1] <- patient_id[index1]
# View(merge_nmf_meta[index,][index1,])
table(merge_nmf_meta$patient_id[index][index1])
table(merge_nmf_meta$orig.ident[index][index1])
# View(merge_nmf_meta[index,])
# unique(merge_nmf_meta[index, 'project_id'])
# unique(merge_nmf_meta[index, 'patient_id'])
diff_pt <- setdiff(unique(merge_nmf_meta$patient_id), unique(sc_majortype$patient_id))#47 -->22
table(sc_majortype$project_id %in% diff_dt)
diff_dt <- unique(merge_nmf_meta$project_id[merge_nmf_meta$patient_id %in% diff_pt])#10
# [1] "E-MTAB-8495" "GSE164898"   "GSE161529"   "GSE159929"   "GSE136103"  
# [6] "GSE156625"   "GSE126030"   "GSE121080"   "GSE150430"   "GSE130973"  
# [11] "GSE172577"   "GSE181989"
diff_DT
setdiff(diff_dt, diff_DT)
# [1] "GSE164898" "GSE161529" "GSE136103" "GSE156625"
# [5] "GSE130973" "GSE172577" "GSE181989"
table(merge_nmf_meta$barcode %in% sc_majortype$barcode)
# View(sc_samples[sc_samples$project_id %in% diff_dt, ])
table(merge_nmf_meta$disease[merge_nmf_meta$project_id %in% diff_dt]) # healthy 
#          aplastic anemia                cirrhotic                  healthy 
#                    29088                    12891                   275018 
# hepatocellular carcinoma              oral cancer 
#                    47534                     8686 
table(unique(merge_nmf_meta$project_id) %in% unique(sc_majortype$project_id))
diff_majortype <- subset(merge_nmf_meta, project_id %in% diff_DT) %>% dplyr::select(colnames(sc_majortype))
# rownames(diff_majortype) <- str_split(rownames(diff_majortype), '[...]', simplify = T)[, 1]
table(diff_majortype$tissue, diff_majortype$disease)
sc_majortype <- rbind(sc_majortype, diff_majortype)

table(unique(sc_majortype$project_id) %in% unique(merge_nmf_meta$project_id))
exclude_DT <- setdiff(unique(sc_majortype$project_id), unique(merge_nmf_meta$project_id)) #69 --> 62
index <- sc_majortype$project_id %in% merge_nmf_meta$project_id
unique(sc_majortype$disease[!index ])
exclude_CT <- as.data.frame(table(sc_majortype$disease[!index ], sc_majortype$project_id[!index ]))
exclude_CT <- exclude_CT[exclude_CT$Freq !=0,]
unique(exclude_CT$Var2) # 62
RMzyme_pheno <- as.data.frame(table(sc_majortype$tissue, sc_majortype$disease))
RMzyme_pheno <- RMzyme_pheno[RMzyme_pheno$Freq!=0,]
# sc_majortype <- sc_majortype[! sc_majortype$project_id %in% exclude_DT, ]
sort(unique(sc_majortype$disease))
sc_majortype$disease[sc_majortype$disease == "severe"] <- 'severe covid-19'
sc_majortype$disease[sc_majortype$disease == "critical"] <- 'critical covid-19'
sc_majortype$disease[sc_majortype$disease == "moderate"] <- 'moderate covid-19'
sc_majortype$disease <- gsub('esophagus squamous cell carcinoma', 'esophageal squamous cell carcinoma',
                             gsub('lung cancer', 'lung carcinoma', 
                                  gsub('inflamed', 'inflammatory bowel disease',
                                       gsub('head and neck squamous carcinoma', 'head and neck squamous cell carcinoma', sc_majortype$disease))))
sc_majortype$tissue[sc_majortype$project_id %in% c('GSE197677', 'GSE222078')] <- 'esophagus'

sc_summary <- sc_majortype[! duplicated(paste(sc_majortype$patient_id, sc_majortype$tissue, 
                                              sc_majortype$project_id)),]
sort(unique(sc_summary$disease)) #171
length(unique(sc_summary$project_id)) #396 --> 340 --> 327--> 332
table(unique(sc_summary$project_id) %in% exclude_DT)
sc_summary <- sc_summary[! sc_summary$project_id %in% exclude_DT,]

sample.index <- c(grep("ormal", sc_summary$sample_type),grep("ealthy", sc_summary$sample_type),
                  grep("ontrol", sc_summary$sample_type))
sc_Diseased_Samples <- sc_summary[ -sample.index, c(2,3,5:8)] # 219 --> 222 disease datasets
unique(sc_Diseased_Samples$project_id) # 222
sc_summary$omics <- 'scRNA'
sc_patient_id <- unique(paste(sc_summary$project_id, sc_summary$patient_id)) #4463 (sc_majortype) --> 4048

##################################################
# change 4:
##################################################
# View(sc_summary[sc_summary$project_id %in% diff_dt, ])
table(is.na(merge_nmf_meta$barcode))
table(grepl('--', merge_nmf_meta$barcode))
index <- ! grepl('--', merge_nmf_meta$barcode)
index1 <- ! grepl('--', rownames(merge_nmf_meta))
merge_nmf_meta$barcode[index & (! index1)] <- rownames(merge_nmf_meta)[index & (! index1)]
index <- ! grepl('--', merge_nmf_meta$barcode)
patient_id <- merge_nmf_meta$patient_id[index&index1]
extract_last_element <- function(x) {
  split_string <- strsplit(x, "_")
  last_element <- sapply(split_string, function(y) tail(y, n = 1))
  return(last_element)
}
barcode.index <- extract_last_element(rownames(merge_nmf_meta)[index & index1])
index2 <- ! grepl('-1...', barcode.index)
remove_first_element <- function(char_vector) {
  # Split the string using "." as the delimiter (properly escaped as "\\.")
  char_vector_parts <- strsplit(char_vector, "\\.")
  # Remove the first element from each part, then join the remaining parts back
  char_vector_parts <- lapply(char_vector_parts, function(parts) {
    if (length(parts) > 1) {
      return(paste(parts[-1], collapse = "."))
    } else {
      return(parts[1])  # If only one part, return it unchanged
    }
  })
  return(unlist(char_vector_parts))
}
barcode.index[index2] <- remove_first_element(barcode.index[index2])
index3 <- grepl("Via.", barcode.index)
barcode.index[index3] <- remove_first_element(barcode.index[index3])
# barcode.index[grepl('GSE172357', merge_nmf_meta$project_id[index & index1])]
# barcode.index <- str_split(barcode.index, '[...]', simplify = T)[, 1]
# barcode.index[grepl('E-MTAB-11948', merge_nmf_meta$project_id)]
barcode.index <- paste0(barcode.index, '--', patient_id)
merge_nmf_meta$barcode[index & index1] <- barcode.index
# duplicated 1: GSE188222 GGGATGAGTTATCACG-1--GSM5672599 ? 

##################################################
# change 5:
##################################################
sc_samples <- merge_nmf_meta[! duplicated(paste(merge_nmf_meta$patient_id, merge_nmf_meta$tissue, 
                                                merge_nmf_meta$project_id)),]
patterns <- c(
  "^([ACGT]+)-\\d.*",                # Starts with [ACGT]+ followed by a dash and digits
  "^([ACGT]+)_\\d.*",                # Starts with [ACGT]+ followed by underscore and digits
  "^([ACGT]+)-[A-Za-z0-9]+$",        # Starts with [ACGT]+ followed by a dash and alphanumeric characters
  "^([ACGT]+)$",                     # Matches only [ACGT]+ sequence
  "^([ACGT]+)-[A-Za-z0-9]+$",        # Starts with [ACGT]+ followed by alphanumeric characters and dash
  "^([ACGT]+)-.*",                   # General pattern starting with [ACGT]+ and a dash
  "^([ACGT]+)_.*",                   # Starts with [ACGT]+ and an underscore
  "^([ACGT]+)\\.\\d+$"               # Starts with [ACGT]+ and a dot followed by digits
)
new_vector <- sc_samples$barcode
# Check if any pattern matches for each element in the new_vector
logical_vector <- sapply(new_vector, function(x) {
  any(sapply(patterns, function(pattern) grepl(pattern, x)))
})
# View(sc_samples[!logical_vector,])
diff_bc <- sc_samples[!logical_vector, 'patient_id']
index <- merge_nmf_meta$patient_id %in% unique(diff_bc)
# View(merge_nmf_meta[index, ])
extract_clean_barcode <- function(barcodes) {
  library(stringr)
  patterns <- c(
    ".*_([ACGT]+)-\\d.*",
    ".*_([ACGT]+)\\.\\d.*",              # Pattern for "GSM5320155_scRNA_ctrl4_COH_Patients_CA1140_Ctrl_ACCAGTATCCTCTAGC.1"
    "^([ACGT]+)-\\d.*",                # Pattern for "AAATCAGAAGTGATGC-0"
    ".*_([ACGT]+)_\\d.*",                # General pattern with _ and digits
    ".*_([ACGT]+)-[A-Za-z0-9]+$",        # Pattern ending with alphanumeric
    ".*_([ACGT]+)$",                     # General pattern ending with barcode
    ".*([ACGT]+)-[A-Za-z0-9]+$",         # Pattern ending with alphanumeric
    ".*-([ACGT]+)$",                     # Pattern for "P1T-E-AAACCTGCACTCTGTC"
    "([ACGT]+).*",                       # New pattern for cases like "AAAACCTCCAATGACCC-6034-YX-3-Cancer-2021"
    # "^([ACGT]+)-.*$",                    # Pattern for barcodes at the start followed by hyphen
    "^([ACGT]+)-[A-Za-z0-9]+$",           # Pattern for "TTTGTGTCTTATCTGT-MSS_A_2-2020"
    "^[A-Za-z0-9]*_([ACGT]+)-.*$",        # Pattern for mixed content ending with barcode
    # ".*_([ACGT]+)-\\d_\\d$",              # Pattern for "CA_HPV_1_AAGCATCTCTGAGAAA-1_1"
    ".*\\.([ACGT]+)-\\d+$"                # Pattern for "GSM5252126_BPH283PrGF_Via.AAACCTGGTCCGAATT-1"
  )
  
  extract_barcode <- function(barcode) {
    for (pattern in patterns) {
      match <- str_match(barcode, pattern)
      if (!is.na(match[2])) {
        return(match[2])
      }
    }
    return(NA)
  }
  
  clean_barcodes <- sapply(barcodes, extract_barcode)
  return(clean_barcodes)
}
barcode.index <-  extract_clean_barcode(rownames(merge_nmf_meta)[index])
table(is.na(barcode.index))
table(merge_nmf_meta$project_id[index][is.na(barcode.index)])# GSE149512 
index1 <- is.na(barcode.index)
barcode.index[index1] <- rownames(merge_nmf_meta)[index][index1]
merge_nmf_meta$barcode[index] <-  paste0(barcode.index, '--', merge_nmf_meta$patient_id[index])

patient_id <- merge_nmf_meta$patient_id
barcode.index <- rownames(merge_nmf_meta)
# barcode.index[grepl('PRJCA001063', merge_nmf_meta$project_id)]
index4 <- grepl('PRJCA001063', merge_nmf_meta$project_id)
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
barcode.index[index4] <- remove_first_element(barcode.index[index4])
# barcode.index[grepl('GSE158291', merge_nmf_meta$project_id)]
index5 <- grepl('GSE158291', merge_nmf_meta$project_id)
rownames(merge_nmf_meta)[index5]
barcode.index[index5] <- remove_first_element(rownames(merge_nmf_meta)[index5])
# barcode.index[grepl('GSE146115', merge_nmf_meta$project_id)]
index6 <- grepl('GSE146115', merge_nmf_meta$project_id)
rownames(merge_nmf_meta)[index6]
barcode.index[index6] <- remove_first_element(rownames(merge_nmf_meta)[index6])
index7 <- grepl('GSE203612', merge_nmf_meta$project_id) & ( ! grepl("_Vis", merge_nmf_meta$patient_id) ) # GSE149512
rownames(merge_nmf_meta)[index7]
barcode.index[index7] <- str_split(rownames(merge_nmf_meta)[index7], '.gz_', simplify = T)[,2]
index8 <- grepl('HTA11', merge_nmf_meta$project_id)
rownames(merge_nmf_meta)[index8]
barcode.index[index8] <- str_split(rownames(merge_nmf_meta)[index8], '[-]', simplify = T)[,1]

index <- index4 | index5 | index6 | index7 | index8
barcode.index[index] <- paste0(barcode.index[index], '--', patient_id[index])
merge_nmf_meta$barcode[index] <- barcode.index[index]

##################################################
# # Check for matches using str_detect
##################################################
# patterns <- c(
#   ".*_([ACGT]+)-\\d.*",
#   ".*_([ACGT]+)\\.\\d.*", 
#   "^([ACGT]+)-\\d.*",
#   ".*_([ACGT]+)_\\d.*",
#   ".*_([ACGT]+)-[A-Za-z0-9]+$",
#   ".*_([ACGT]+)$",
#   ".*([ACGT]+)-[A-Za-z0-9]+$",
#   ".*-([ACGT]+)$",
#   "([ACGT]+).*",
#   "^([ACGT]+)-[A-Za-z0-9]+$",
#   "^[A-Za-z0-9]*_([ACGT]+)-.*$",
#   ".*\\.([ACGT]+)-\\d+$"
# )
# new_vector <- sc_samples$barcode
# # Check if any pattern matches for each element in the new_vector
# logical_vector <- sapply(new_vector, function(x) {
#   any(sapply(patterns, function(pattern) grepl(pattern, x)))
# })
# View(sc_samples[!logical_vector,])
# diff_bc <- sc_samples[!logical_vector, 'patient_id']
# index <- merge_nmf_meta$patient_id %in% unique(diff_bc)
# View(merge_nmf_meta[index, ])
# barcode.index <- extract_last_element(rownames(merge_nmf_meta)[index])
# barcode.index <- str_split(barcode.index, '[...]', simplify = T)[, 1]
# barcode.index <- str_split(rownames(merge_nmf_meta)[index], "_", simplify = T)[,1]
# barcode.index <- str_split(barcode.index, "-Pt", simplify = T)[,1]
# merge_nmf_meta$barcode[index] <-  paste0(barcode.index, '--', merge_nmf_meta$patient_id[index])
# index <- ! merge_nmf_meta$barcode %in% sc_majortype$barcode
# df_barcode <- merge_nmf_meta[ index, ]
# df_DT <- unique(df_barcode$project_id) # 34 -> 29
# length(unique(df_barcode$patient_id)) # 408

##################################################
# change 6:
##################################################
scMajortype_ct <- as.data.frame(table(sc_majortype$sub_type, sc_majortype$ct))
scMajortype_ct <- scMajortype_ct[scMajortype_ct$Freq!=0,]
index <- is.na(merge_nmf_meta$ct)
dts <- unique(merge_nmf_meta$project_id[index]) # [1] "GSE185991" "GSE185381"
table(unique(sc_majortype$project_id) %in% dts)
ct_dts <- subset(sc_majortype, project_id %in% dts)
table(merge_nmf_meta$barcode[index] %in% ct_dts$barcode)
index1 <- ! merge_nmf_meta$barcode[index] %in% ct_dts$barcode
dts <- unique(merge_nmf_meta[index,'project_id'][index1]) # [1] "GSE185991"
unique(merge_nmf_meta[index,'tissue'][index1])
unique(merge_nmf_meta[index,'patient_id'][index1])
# # View(sc_summary[sc_summary$project_id %in% unique(merge_nmf_meta[index,'project_id'][index1]),])
# ct_dts <- subset(sc_majortype, project_id %in% dts)
# unique(ct_dts$tissue) #[1] "bone marrow"
# # View(ct_dts[ct_dts$project_id %in% unique(merge_nmf_meta[index,'project_id'][index1]),])
# # View( merge_nmf_meta[index,][index1,] )
patient_id <- merge_nmf_meta$patient_id[index][index1]
barcode.index <- str_split(merge_nmf_meta$barcode[index][index1], '--', simplify = T)[,1]
barcode.index <- str_split(barcode.index, '[...]', simplify = T)[, 1]
barcode.index <- paste0(barcode.index, '--', patient_id)
table(barcode.index %in% ct_dts$barcode)
barcode.index <- c(barcode.index, merge_nmf_meta$barcode[index][!index1])
table( merge_nmf_meta$barcode[index][index1] %in% ct_dts$barcode)

index2 <- match(barcode.index, ct_dts$barcode)
merge_nmf_meta$ct[index] <- ct_dts$ct[index2]
table(is.na(merge_nmf_meta$ct))
merge_nmf_meta$Major_type[index] <- ct_dts$Major_type[index2]
merge_nmf_meta$sub_type[index] <- ct_dts$sub_type[index2]
RMzyme_ct <- as.data.frame(table(merge_nmf_meta$sub_type, merge_nmf_meta$ct))
RMzyme_ct <- RMzyme_ct[RMzyme_ct$Freq!=0,]
index <- merge_nmf_meta$sub_type=='MastCells'
merge_nmf_meta$sub_type[index] <- 'MastCell'
##################################################
# # changes on old version of merge_nmf_rds
##################################################
# [1] "GSE208653"    "SCP1950"      "GSE168652"    "E-MTAB-11948"
# [5] "E-MTAB-12305" "S-BSST1035"   "GSE185991"   
# remove_last_element <- function(char_vector) {
#   char_vector_parts <- strsplit(char_vector, "_")
#   char_vector_parts <- lapply(char_vector_parts, function(parts) {
#     if (length(parts) > 1) {
#       return(paste(parts[-length(parts)], collapse = "_"))
#     } else {
#       return(parts[1])  # If only one part exists, return it
#     }
#   })
#   return(unlist(char_vector_parts))
# }
# index2 <- grepl('-1_', ct_dts$barcode)
# ct_dts$barcode[index2] <- remove_last_element(ct_dts$barcode[index2])
# extract_last_element <- function(x) {
#   split_string <- strsplit(x, "_")
#   last_element <- sapply(split_string, function(y) tail(y, n = 1))
#   return(last_element)
# }
# index2 <- grepl('_', ct_dts$barcode)
# ct_dts$barcode[index2] <- extract_last_element(ct_dts$barcode[index2])
##################################################
# # View(merge_nmf_meta[merge_nmf_meta$sub_type=='Tcells' & merge_nmf_meta$ct=='Mast cells', ])
# # View(merge_nmf_meta[merge_nmf_meta$sub_type=='Tissue_Specific_Cells' & merge_nmf_meta$ct=='Stem-like cells', ])
# index <- merge_nmf_meta$sub_type=='Tissue_Specific_Cells' & merge_nmf_meta$ct=='Stem-like cells'
# table(merge_nmf_meta$seurat_clusters[index] );table(merge_nmf_meta$project_id[index] )
# # GSE174748 GSE183676 GSE212046 
# # 1135      5820       255
# # merge_nmf_meta$sub_type[index] <- 'Stem_Cells'
# 
# index <- merge_nmf_meta$sub_type=='Tcells' & merge_nmf_meta$ct=='Mast cells'
# table(merge_nmf_meta$seurat_clusters[index] );table(merge_nmf_meta$project_id[index] )
# # GSE183676     HTA12 
# # 1354     11283
# # merge_nmf_meta$sub_type[index] <- 'MastCell'
##################################################
# change 7:
##################################################
table(is.na(merge_nmf_meta$sub_type))
table(is.na(merge_nmf_meta$Major_type))
table(merge_nmf_meta$sub_type)
##################################################
index1 <- grepl('Bcells|Macrophages|Neutrophils|NKcells|OtherImmunecells|Tcells|MastCell|Myeloids',merge_nmf_meta$sub_type)
merge_nmf_meta$Major_type[index1] <- 'Immune_Cells'
index2 <- grepl('OtherMesenchymals|Fibroblasts',merge_nmf_meta$sub_type)
merge_nmf_meta$Major_type[index2] <- 'Mesenchymal_Cells'
index3 <- grepl('Endothelial_Cells',merge_nmf_meta$sub_type)
merge_nmf_meta$Major_type[index3] <- 'Endothelial_Cells'
index4 <- grepl('Tissue_Specific_Cells',merge_nmf_meta$sub_type)
merge_nmf_meta$Major_type[index4] <- 'Tissue_Specific_Cells'
index5 <- grepl('Stem_Cells',merge_nmf_meta$sub_type)
merge_nmf_meta$Major_type[index5] <- 'Stem_Cells'
index6 <- grepl('Blood_Cells',merge_nmf_meta$sub_type)
merge_nmf_meta$Major_type[index6] <- 'Blood_Cells'
table(merge_nmf_meta$Major_type)
table(is.na(merge_nmf_meta$Major_type))
merge_nmf_meta$Major_type[merge_nmf_meta$Major_type=='Epithelial'] <- 'Tissue_Specific_Cells'
merge_nmf_meta$Major_type[merge_nmf_meta$Major_type=='Stromal'] <- 'Mesenchymal_Cells'
merge_nmf_meta$Major_type[merge_nmf_meta$Major_type=='Immune'] <- 'Immune_Cells'
index7 <- merge_nmf_meta$sub_type=='Stem_Cells'
table(merge_nmf_meta$Major_type[index7])
table(merge_nmf_meta$ct[index7 & merge_nmf_meta$Major_type=='Tissue_Specific_Cells'])
# # View(merge_nmf_meta[merge_nmf_meta$ct=='Metaplasia stem cells', ])
# merge_nmf_meta$Major_type[index7] <- gsub('Tissue_Specific_Cells','Stem_Cells', merge_nmf_meta$Major_type[index7])
# index <- is.na(merge_nmf_meta$Major_type)
# table(merge_nmf_meta$sub_type[index])

##################################################
# change 8:
##################################################
# merge_nmf_meta <- readRDS(file = '/data/rluo4/RPMfunc/Output/scRNA/merge_nmf_meta_0915.rds')
setdiff(unique(merge_nmf_meta$disease), unique(sc_majortype$disease))
setdiff(unique(merge_nmf_meta$tissue), unique(sc_majortype$tissue))
sort(unique(merge_nmf_meta$disease))
sort(unique(sc_majortype$disease))
setdiff(unique(sc_majortype$disease), unique(merge_nmf_meta$disease))
setdiff(unique(sc_majortype$tissue), unique(merge_nmf_meta$tissue))
setdiff(unique(sc_summary$disease), unique(merge_nmf_meta$disease))
setdiff(unique(sc_summary$tissue), unique(merge_nmf_meta$tissue))
setdiff(unique(sc_summary$patient_id), unique(merge_nmf_meta$patient_id))#
merge_nmf_meta$patient_id <- str_split(merge_nmf_meta$patient_id, '[...]', simplify = T)[, 1]
merge_patient <- unique(paste(merge_nmf_meta$project_id, merge_nmf_meta$patient_id)) #4463 (sc_majortype) --> 4048
table(merge_patient %in% sc_patient_id)
diff_pt <- setdiff(merge_patient, sc_patient_id)#setdiff(unique(merge_nmf_meta$patient_id),unique(sc_summary$patient_id))
index <- paste(merge_nmf_meta$project_id, merge_nmf_meta$patient_id) %in% diff_pt
diff_dt <- merge_nmf_meta[index,]
unique(diff_dt$project_id)
# [1] "GSE164898" "GSE171555" "GSE130973"
# [4] "GSE167297" "GSE172577" "GSE144236"
# [7] "GSE137829" "GSE136103" "GSE156625"
# [10] "GSE112271"
# View(sc_summary[sc_summary$project_id %in% unique(diff_dt$project_id),])
diff_dt$patient_id <- str_split(diff_dt$barcode, '--', simplify = T)[,2]
table(unique(diff_dt$project_id) %in% tissue_summary$Dataset)
unique(diff_dt$disease)
merge_nmf_meta[index,'patient_id'] <- diff_dt$patient_id
merge_patient <- paste(merge_nmf_meta$project_id, merge_nmf_meta$patient_id)
sc_patient_id <- paste(sc_summary$project_id, sc_summary$patient_id)
index1 <- merge_patient %in% sc_patient_id
index2 <- match(merge_patient[index1], sc_patient_id)
merge_nmf_meta$disease[index1] <- sc_summary$disease[index2]
merge_nmf_meta$sample_type[index1] <- sc_summary$sample_type[index2]
merge_nmf_meta$tissue[index1] <- sc_summary$tissue[index2]
# View(merge_nmf_meta[ ! index1, ])
setdiff(unique(sc_summary$disease), unique(merge_nmf_meta$disease)) # [1] "hashimoto's thyroiditis"
setdiff(unique(sc_summary$sample_type), unique(merge_nmf_meta$sample_type)) # [1] "HT"
# test <- subset(merge_nmf_meta, project_id == 'GSE136103')
# unique(test$patient_id)
# test_major <- subset(sc_majortype, project_id == 'GSE136103')
# unique(test_major$patient_id)
# setdiff(unique(test$patient_id), unique(test_major$patient_id))
setdiff( unique(merge_nmf_meta$disease), unique(sc_summary$disease))
merge_nmf_meta$disease <- gsub('oral cancer', 'oral squamous cell carcinoma', merge_nmf_meta$disease)
# unique(merge_nmf_meta$patient_id)
# saveRDS(merge_nmf_meta, file = '/data/rluo4/RPMfunc/Output/scRNA/merge_nmf_meta_0925.rds')

sc_metadata <- merge_nmf_meta %>% dplyr::select( c(colnames(sc_majortype), 'seurat_clusters') )
unique(sc_metadata$project_id) #270
unique(sc_metadata$patient_id) #3596 --> 2600
unique(sc_metadata$platform)
index1 <- is.na(sc_metadata$platform)
# table(sc_metadata$patient_id %in% c( disco$patient_id, tissue_summary$Replicates) )
index2 <- match(sc_metadata$patient_id[index1], disco$patient_id)
table(is.na(index2))
sc_metadata$platform[index1] <- disco$platform[index2]
index1 <- is.na(sc_metadata$platform)
unique(sc_metadata$project_id[index1]) #[1] "CRA002497"   "GSE155673"   "covid_atlas" "GSE155224"
sc_metadata$platform[sc_metadata$project_id=='CRA002497'] = "10x5'v2"
sc_metadata$platform[sc_metadata$project_id=='GSE155673'] = "10x3'v3"
sc_metadata$platform[sc_metadata$project_id=='covid_atlas'] = "10x5'v2"
sc_metadata$platform[sc_metadata$project_id=='GSE155224'] = "10x5'"#"10x3'v3" #only pbmc? where are the blood tissues?
sc_metadata$platform[sc_metadata$platform=='10x3'] = "10x3'"
unique(sc_metadata$platform)
# [1] "10x3'v2"       "10x3'v3"       "10x3'"        
# [4] "10x5'"         "10x5'v2"       "10X Genomics" 
# [7] "10x"           "10x3'v1"       "inDrop V2"    
# [10] "BD Rhapsody"   "Microwell-seq" "Seq-well"     
# [13] "MARS-seq"
unique(sc_metadata$project_id[sc_metadata$platform=='10x'] )
sc_metadata$platform[sc_metadata$project_id=='GSE149689'] = "10x3'v3"
sc_metadata$platform[sc_metadata$project_id=='GSE168453'] = "10x5'"
sc_metadata$platform[sc_metadata$project_id=='GSE166992'] = "10x5'"
sc_metadata$platform[sc_metadata$project_id=='GSE171555'] = "10x5'v1"
sc_metadata$platform[sc_metadata$project_id=='GSE192391'] = "10x5'"
sc_metadata$platform[sc_metadata$project_id=='GSE154567'] = "10x3'v3"
sc_metadata$platform[sc_metadata$project_id=='GSE208337'] = "10x3'v3"
sc_metadata$platform[sc_metadata$project_id=='GSE150861'] = "10x3'v2"
unique(sc_metadata$project_id[sc_metadata$platform=='10X Genomics'] )
sc_metadata$platform[sc_metadata$platform=='10X Genomics'] = "10x"
NMF_clusters <- unique(sc_metadata$seurat_clusters)
NMF_RMPs <- str_split(NMF_clusters, "[+]", simplify = T)[,1]
sort(unique(NMF_RMPs))
# index <- grepl("NA+", sc_metadata$seurat_clusters)
# View(sc_metadata[index,])
# NMF_clusters <- gsub('^NA\\+$', 'NoneRMP+', NMF_clusters)
sc_metadata$seurat_clusters <- gsub('NA+', 'NoneRMP+', sc_metadata$seurat_clusters, fixed = TRUE)
sc_metadata$seurat_clusters <- gsub('MastCell', 'MastCells',  sc_metadata$seurat_clusters, fixed = TRUE)
sc_metadata$sub_type <- gsub('MastCell', 'MastCells',  sc_metadata$sub_type, fixed = TRUE)
table(is.na(sc_metadata$seurat_clusters))
# saveRDS(sc_metadata, file = '/data/rluo4/RPMfunc/Output/scRNA/sc_metadata_0925.rds')

sc_samples <- sc_metadata[! duplicated(paste(sc_metadata$patient_id, sc_metadata$tissue, 
                                             sc_metadata$project_id)),]
# save(disco, tissue_summary, sc_summary, sc_samples, file = '/data/rluo4/RPMfunc/Output/scRNA/sc_samplemeta_0925.RData')
sort(unique(sc_samples$tissue))
sort(unique(sc_samples$disease))

merge_patient <- paste(sc_samples$project_id, sc_samples$patient_id)
index1 <- merge_patient %in% sc_patient_id
# View(sc_samples[ ! index1, ])
unique(sc_samples[ ! index1, 'project_id'])
# [1] "GSE164898" "GSE171555" "GSE130973" "GSE167297"
# [5] "GSE172577" "GSE144236" "GSE137829" "GSE136103"
# [9] "GSE156625"
setdiff(unique(sc_metadata$project_id), unique(sc_majortype$project_id))
setdiff(unique(sc_metadata$patient_id), unique(sc_majortype$patient_id))
unique(sc_metadata$project_id[! sc_metadata$patient_id %in% unique(sc_majortype$patient_id)])
# [1] "GSE164898" "GSE171555" "GSE130973" "GSE167297" "GSE172577" "GSE136103"
# [7] "GSE156625"
# table(unique(sc_majortype$project_id) %in% sc_summary$project_id)
diff_dt <- setdiff(unique(sc_samples$project_id),unique(sc_Diseased_Samples$project_id))
table(merge_nmf_meta$disease[merge_nmf_meta$project_id %in% diff_dt]) # healthy 579268
table(merge_nmf_meta$tissue[merge_nmf_meta$project_id %in% diff_dt])  
setdiff(unique(sc_Diseased_Samples$project_id), unique(merge_nmf_meta$project_id))
# View(merge_nmf_meta[merge_nmf_meta$project_id %in% diff_dt,])
# unique(sc_majortype[! grepl('adjacent|ealthy', sc_majortype$disease),]$project_id)

##################################################
# 8) #        Statistics and Overview            #
##################################################
######################################################################
# 8.1. Summary of RNA modifications across single-cell transcriptome #
######################################################################
library(Seurat)
library(dplyr)
library(stringr)
RMP_update <- readxl::read_excel('/data/rluo4/EpiTrans/DataCollection/RMP_update.xlsx')

sc_allmeta <- readRDS(file = '/data/rluo4/RPMfunc/Output/scRNA/sc_allmeta_1122.rds')#sc_allmeta_0925.rds') 
RMzyme_alltissues <- sort(unique(sc_allmeta$tissue));
RMzyme_alldatasets <- sort(unique(sc_allmeta$project_id));
load(file = '/data/rluo4/RPMfunc/Output/scRNA/sc_allsample_1122.RData', verbose = T)
table(sc_allsample$patient_id %in% sc_summary$patient_id)
diff_pt <- setdiff(sc_allsample$patient_id, sc_summary$patient_id)
# View(sc_allsample[sc_allsample$patient_id %in% diff_pt,])
sc_allmeta$omics <- 'scRNA'
sc_allsample$omics <- 'scRNA'
RMzyme_pheno <- as.data.frame(table(sc_allmeta$tissue, sc_allmeta$disease))
RMzyme_pheno <- RMzyme_pheno[RMzyme_pheno$Freq!=0,]
RMzyme_ct <- as.data.frame(table(sc_allmeta$sub_type, sc_allmeta$ct))
RMzyme_ct <- RMzyme_ct[RMzyme_ct$Freq!=0,]
NMF_clusters <- unique(sc_allmeta$seurat_clusters)
length(NMF_clusters)#[1] 7321 --> 7304
NMF_RMPs <- str_split(NMF_clusters, "[+]", simplify = T)[,1]
sort(unique(NMF_RMPs))#178=177+1 NoneRMP
setdiff(RMP_update$RMP, NMF_RMPs) #"RET1" "RET2"
NMF_clusters <- sc_allmeta%>% dplyr::select(ddt_tis, ddt_id, split_ids, project_id, patient_id,  tissue, disease,
                                            barcode,   Major_type, sub_type, ct, NMF_clusters = seurat_clusters)
NMF_clusters$RMP <- str_split(NMF_clusters$NMF_clusters, "[+]", simplify = T)[,1]
NMF_clusters$disease <- gsub('adjacent tumor/disease section', 'healthy', NMF_clusters$disease)
remove_first_element <- function(char_vector) {
  char_vector_parts <- strsplit(char_vector, "[+]")
  char_vector_parts <- lapply(char_vector_parts, function(parts) {
    if (length(parts) > 1) {
      return(paste(parts[-1], collapse = "[+]"))
    } else {
      return(parts[1])
    }
  })
  return(unlist(char_vector_parts))
}
# NMF_clusters$clusters <-  remove_first_element(NMF_clusters$NMF_clusters);NMF_clusters <- NMF_clusters %>% dplyr::select(-clusters)
# Optimized function using stringi
library(stringi)
library(data.table)
remove_first_element_fast <- function(char_vector) {
  # Split based on "+" without treating it as a special character
  char_vector_parts <- stri_split_fixed(char_vector, "+", omit_empty = TRUE)
  
  # Apply vectorized operation to return all but the first part
  remaining_parts <- sapply(char_vector_parts, function(parts) {
    if (length(parts) > 1) {
      return(stri_join(parts[-1], collapse = "+"))  # Join all but the first with "+"
    } else {
      return(parts[1])  # If only one element, return it as-is
    }
  })
  
  return(remaining_parts)
}
# Use data.table for fast operations with large data
NMF_clusters <- data.table(NMF_clusters)
NMF_clusters[, clusters := remove_first_element_fast(NMF_clusters)]
table(is.na(NMF_clusters$clusters))
NMF_clusters <- as.data.frame(NMF_clusters)
# View(NMF_clusters[NMF_clusters$sub_type=='Bcells' & NMF_clusters$clusters=='Tcells' ,])
NMF_clusters$cluster_type <- ifelse(NMF_clusters$clusters %in% NMF_clusters$sub_type, 'Major', 'Minor')
index = NMF_clusters$cluster_type=='Major'
testMa = (NMF_clusters[index & NMF_clusters$sub_type != NMF_clusters$clusters, ])
unique(testMa$project_id) #[1] "GSE185991"
index = NMF_clusters$cluster_type=='Minor'
testMi = (NMF_clusters[index & gsub(' |/','_',NMF_clusters$ct) != NMF_clusters$clusters, ])
unique(testMi$project_id)
testTis = (NMF_clusters[NMF_clusters$ddt_tis!=NMF_clusters$tissue, ])
error_tis <-  as.data.frame(table(testTis$ddt_tis, testTis$tissue))
error_tis <- error_tis[error_tis$Freq>0, ]
unique(testTis$project_id[grepl('esophagus', testTis$tissue) | grepl('esophagus', testTis$ddt_tis)])
unique(testTis$project_id[grepl('cartilage', testTis$tissue) | grepl('cartilage', testTis$ddt_tis)])
unique(testTis$ddt_id[grepl('blood', testTis$tissue) & grepl('bone marrow', testTis$ddt_tis)])
error_idx <- grepl('blood', testTis$tissue) & grepl('lung', testTis$ddt_tis)
unique(testTis$patient_id[error_idx])
# View(table(testTis$ddt_tis, testTis$tissue))
# View(test[!duplicated(test$NMF_clusters), ])
NMF_clusters$disease_type <- NMF_clusters$disease
NMF_clusters$disease_type[grepl('cancer|carcinoma|leukemia|lymphoma|blastoma|melanoma|myeloma', NMF_clusters$disease_type)] ='cancer'
NMF_clusters$disease_type[! grepl('cancer|carcinoma|leukemia|lymphoma|blastoma|melanoma|myeloma|healthy', NMF_clusters$disease_type)] ='non-cancer'
NMF_pheno <- as.data.frame(table(NMF_clusters$disease_type, NMF_clusters$disease))
NMF_pheno <- NMF_pheno[NMF_pheno$Freq!=0,]

unique(NMF_clusters$clusters[NMF_clusters$cluster_type == 'Minor'])
unique(NMF_clusters$clusters[NMF_clusters$cluster_type == 'Major'])
table(unique(NMF_clusters$RMP) %in% RMP_update$RMP)
index <- match(NMF_clusters$RMP, RMP_update$RMP)
# NMF_clusters <- left_join(NMF_clusters, RMP_update, by = 'RMP' )
NMF_clusters$RNA_modifcations <- RMP_update$RNA_modifications[index]
NMF_clusters$RMP_type <- RMP_update$RMP_type[index]
length(unique(NMF_clusters$NMF_clusters)) #7135 --> 7329 --> 7210
colnames(NMF_clusters)
tableS8 <- NMF_clusters[, c("ddt_tis", "ddt_id", 'project_id', #'patient_id', #'disease',  'disease_type',
                            'NMF_clusters', 'clusters', 'sub_type', 'RMP', 'RNA_modifcations', 'RMP_type')]
colnames(tableS8) <-  c("disease_tissue", "disease_dataset_id", 'project_id',
                        'NMF_clusters', 'clusters', 'sub_type', 'RMP', 'RNA_modifcations', 'RMP_type')
identifiers <- paste(tableS8$tissue, tableS8$project_id, tableS8$NMF_clusters)
table(duplicated(identifiers))
tableS8 <- tableS8[! duplicated(identifiers), ]
write.table(tableS8, '/data/rluo4/RPMfunc/Output/summary/NMF_clusters.txt', sep = '\t', row.names = F, quote = F)

# allres_path <- "/data/rluo4/EpiTrans/DataCollection/RMzyme_scRNA"
# allres_files <- list.files(path = allres_path, pattern = "_CCI_pairLR.use.rds",#"aggregateNet_cellchat.rds$", 
#                            full.names = TRUE, recursive = T)
# strsplit(allres_files, '/')[[1]]
# extract_element <- function(x) {
#   split_string <- strsplit(x, "/")
#   len_element <- sapply(split_string, function(y) length(y))
#   return(len_element)
# }
# allres_files <- data.frame( file = allres_files, sep_len = extract_element(allres_files))
# allres_files <- allres_files %>% dplyr::mutate(
#   tissue = str_split(file, '/', simplify = T)[,7],
#   ddt_id = str_split(file, '/', simplify = T)[,8],
#   dataset = ifelse(sep_len==10, str_split(file, '/', simplify = T)[,8], str_split(file, '/', simplify = T)[,9]),
#   sub_type = ifelse(sep_len==10, str_split(file, '/', simplify = T)[,9], str_split(file, '/', simplify = T)[,10])
# )
# large_split <- unique(allres_files$ddt_id[allres_files$sep_len>10])

################################################################################
# 8.1. Overview of RMP annotated NMF clusters   #
################################################################################
disease_freq <- as.data.frame(table(NMF_clusters$disease_type))
colnames(disease_freq)[2] <- 'disease_freq'
RMzyme_RMP <- as.data.frame(table(NMF_clusters$disease_type, NMF_clusters$RMP))
RMzyme_RMP$RNAmodType <- RMP_update$RNA_modifications[match(RMzyme_RMP$Var2, RMP_update$RMP)]
RMzyme_RMP <- RMzyme_RMP[RMzyme_RMP$Freq!=0,]
RMzyme_RMP <- left_join(RMzyme_RMP, disease_freq, by = 'Var1')
RMzyme_RMP$frequency <- RMzyme_RMP$Freq/RMzyme_RMP$disease_freq#
RMzyme_RMP$frequency <- round(RMzyme_RMP$Freq/RMzyme_RMP$disease_freq,4)*100
# View(RMzyme_RMP[order(RMzyme_RMP$frequency,decreasing = T),])
RMzyme_RM <- as.data.frame(table(NMF_clusters$disease_type, NMF_clusters$RNA_modifcations))
RMzyme_RM <- RMzyme_RM[RMzyme_RM$Freq!=0,]
RMzyme_RM <- left_join(RMzyme_RM, disease_freq, by = 'Var1')
RMzyme_RM$frequency <- RMzyme_RM$Freq/RMzyme_RM$disease_freq#
RMzyme_RM$frequency <- round(RMzyme_RM$Freq/RMzyme_RM$disease_freq,4)*100
# View(RMzyme_RM[order(RMzyme_RM$frequency,decreasing = T),])

cell_freq <- as.data.frame(table(NMF_clusters$sub_type))
colnames(cell_freq)[2] <- 'cell_freq'
RMzyme_RMP <- as.data.frame(table(NMF_clusters$sub_type, NMF_clusters$RMP))
RMzyme_RMP$RNAmodType <- RMP_update$RNA_modifications[match(RMzyme_RMP$Var2, RMP_update$RMP)]
RMzyme_RMP <- RMzyme_RMP[RMzyme_RMP$Freq!=0,]
RMzyme_RMP <- left_join(RMzyme_RMP, cell_freq, by = 'Var1')
RMzyme_RMP$frequency <- RMzyme_RMP$Freq/RMzyme_RMP$cell_freq#
RMzyme_RMP$frequency <- round(RMzyme_RMP$Freq/RMzyme_RMP$cell_freq,4)*100
# View(RMzyme_RMP[order(RMzyme_RMP$frequency,decreasing = T),])

RMzyme_RM <- as.data.frame(table(NMF_clusters$sub_type, NMF_clusters$RNA_modifcations))
RMzyme_RM <- RMzyme_RM[RMzyme_RM$Freq!=0,]
RMzyme_RM <- left_join(RMzyme_RM, cell_freq, by = 'Var1')
RMzyme_RM$frequency <- RMzyme_RM$Freq/RMzyme_RM$cell_freq#
RMzyme_RM$frequency <- round(RMzyme_RM$Freq/RMzyme_RM$cell_freq,4)*100
# RMzyme_RM$Var1 <- gsub('MastCell', 'MastCells', gsub('_', '', RMzyme_RM$Var1))
RMzyme_RM$Var1 <- as.factor(RMzyme_RM$Var1)
colnames(RMzyme_RM)[1:2] <- c('Celltype', 'RNAmodType')
# View(RMzyme_RM[order(RMzyme_RM$frequency,decreasing = T),])

library(echarts4r)
library(htmlwidgets)
library(webshot2)  # Only if you want to convert HTML to PNG
# Create a unique list of cell types
unique_celltypes <- unique(RMzyme_RM$Celltype)
# Create a folder to save images (optional)
dir.create("/data2/rluo4/RPMfunc/Output/scRNA/figures/")
dir.create("/data2/rluo4/RPMfunc/Output/scRNA/figures/celltype_pies")
####################################
#            Figure 3c             #
####################################
# Loop through each Celltype and save the charts as HTML
for (celltype in unique_celltypes) {
  # Create the pie chart for each celltype
  chart <- RMzyme_RM %>%
    filter(Celltype == celltype) %>%
    e_charts(RNAmodType) %>%
    e_pie(
      Freq,
      radius = c("50%", "70%"),
      itemStyle = list(
        borderRadius = 20,
        borderColor = '#fff',
        borderWidth = 2 # 新版本该设置下画出来的是圆角环形图
      )
    ) %>%
    # e_title(text = paste("Pie chart for", celltype)) %>%
    e_show_loading()
  
  # Save the chart as an HTML file
  html_file <- paste0("/data2/rluo4/RPMfunc/Output/scRNA/figures/celltype_pies/", celltype, "_chart.html")
  saveWidget(chart, file = html_file)
  
  # Convert the HTML file to a PNG image using webshot2
  # png_file <- paste0("/data2/rluo4/RPMfunc/Output/scRNA/figures/celltype_pies/", celltype, "_chart.png")
  # Now run your webshot2 code
  # webshot2::webshot(html_file, file = png_file, vwidth = 1500, vheight = 1500, delay = 3)
  
}
# Sys.setenv(CHROMOTE_CHROME = "/path/to/chrome/or/chromium")
# # webshot::install_phantomjs(force = TRUE);system("phantomjs --version")
# index = RMzyme_RM$frequency >10
# View(RMzyme_RM[index,][order(RMzyme_RM$frequency[index],decreasing = T),])
# # View(RMzyme_RM[RMzyme_RM$frequency >10,])
# # View(as.data.frame(table(NMF_clusters$RNA_modifcations)))
color_define=c("#f89e81","#99a9cc","#dd9bc5","#84c7b3","#acd485","#1B9E77","#f6d573","#66A61E","skyblue")
custom_colors = c(RColorBrewer::brewer.pal(n = 12, name = 'Paired'),color_define)
rna_mod_colors <- c("#f89e81","#99a9cc","#dd9bc5","#84c7b3","#acd485","#1B9E77","#f6d573","skyblue",
                    '#A6CEE3','#1F78B4','#B2DF8A','#33A02C','#FB9A99','#E31A1C','#FDBF6F','#FF7F00','#CAB2D6','#6A3D9A','#FFFF99','#B15928' #RColorBrewer::brewer.pal(n = 12, name = 'Paired')
)
####################################
#            Figure 3d             #
####################################
library(circlize)
df.cir <- RMzyme_RMP[order(RMzyme_RMP$frequency, decreasing = T), ]
colnames(df.cir)[1:2] <- c('sectors','RMP')
df.cir$sectors <- gsub('_','', df.cir$sectors)
df.cir$RNAmodType[is.na(df.cir$RNAmodType)] <- 'None'
df.cir$RNAmodType <- as.factor(df.cir$RNAmodType)
df.cir$RMP <- as.factor(df.cir$RMP)
set.seed(123)
n = nrow(df.cir) #nrow(NMF_clusters)#
df.cir$x <- rnorm(n);  df.cir$y = runif(n)
# df.cir$value <- ifelse(df.cir$RNAmodType=='None', -1, 1)
# df.cir = data.frame(sectors = sample(unique(RMzyme_RM$sectors), n, replace = TRUE),
#                 x = rnorm(n), y = runif(n))
output_file <- "/data2/rluo4/RPMfunc/Output/scRNA/figures/circos_plot.png" 
# Set up the PNG device and define dimensions (e.g., 1000x1000 pixels)
png(output_file, width = 7, height = 7, res = 500,units = 'in')
circos.par("track.height" = 0.1)
# Initialize the sectors (first track)
# First Track: Sectors (e.g., cell types)
circos.initialize(df.cir$sectors, x = df.cir$x)
# circos.initialize(factors = df.cir$sectors, x = df.cir$frequency)
set.seed(123)
# First Track: Sector Background and Labels
n_sectors <- length(unique(df.cir$sectors))
if (length(custom_colors) < n_sectors) {
  custom_colors <- rep(custom_colors, length.out = n_sectors)
}
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  chr = CELL_META$sector.index
  xlim = CELL_META$xlim
  ylim = CELL_META$ylim
  # Get the color for the current sector
  sector_color <- custom_colors[which(unique(df.cir$sectors) == chr)]
  circos.rect(xlim[1], 0, xlim[2], 1, col = sector_color)
  circos.text(mean(xlim), mean(ylim), chr, cex = 0.6, col = "black",
              facing = "inside", niceFacing = TRUE)
}, track.height = 0.15, bg.border = NA)
# Second Track: Display Top 3 RMPs by frequency, arranged side by side (skip "NoneRMP")
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  sector_name = CELL_META$sector.index
  xlim = CELL_META$xlim
  # Get all rows for the current sector, excluding "NoneRMP"
  sector_data <- df.cir[df.cir$sectors == sector_name & df.cir$RMP != "NoneRMP", ]
  # Sort by frequency and take the top 3 RMPs (excluding "NoneRMP")
  top_rmp <- head(sector_data[order(-sector_data$frequency), "RMP"], 3)
  # Divide the x-axis range of each sector into three parts (for side-by-side placement)
  x_positions <- seq(xlim[1], xlim[2], length.out = 4)  # 3 sections between the 4 boundaries
  # Display each RMP in its own x segment, with vertical text orientation
  for (i in seq_along(top_rmp)) {
    circos.text(mean(x_positions[i:(i+1)]), 0.5, top_rmp[i], cex = 0.5, 
                facing = "reverse.clockwise", niceFacing = TRUE)
  }
}, bg.border = NA)
# Third Track: Frequency as a bar plot
bgcol = rep(c("#EFEFEF", "#CCCCCC"), 8)
circos.trackHist(df.cir$sectors, df.cir$frequency, bin.size = 0.2, bg.col = bgcol, col = NA)

# Fourth Track: Count distribution per sector using dots or boxplot (illustrative)
# Aggregate counts of RMPs per cell type
# rmp_counts <- aggregate(RMP ~ sectors, data = df.cir, FUN = function(x) length(unique(x)))
# rmp_counts$counts <- rmp_counts$RMP
# rmp_counts$RMP <- NULL
# Initialize the fourth track
rmp_counts <- df.cir %>%
  group_by(sectors, RNAmodType) %>%
  summarise(counts = n()) %>%
  ungroup()
# Create a color palette for RNAmodTypes
# rna_mod_colors <- setNames(rand_color(length(unique(rmp_counts$RNAmodType))), unique(rmp_counts$RNAmodType))
rna_mod_colors <- setNames(rna_mod_colors,levels(rmp_counts$RNAmodType))
# Plot dots for each RNAmodType without labels
circos.trackPlotRegion(ylim = c(0, max(rmp_counts$counts) + 5), track.height = 0.1, 
                       panel.fun = function(x, y) {
                         sector_name = CELL_META$sector.index
                         xlim = CELL_META$xlim
                         # Filter the counts for the current sector
                         sector_counts <- rmp_counts[rmp_counts$sectors == sector_name, ]
                         # Plot dots for each RNAmodType in this sector
                         if (nrow(sector_counts) > 0) {
                           for (i in 1:nrow(sector_counts)) {
                             # x is distributed evenly within the sector for different RNAmodTypes
                             x_pos = xlim[1] + (xlim[2] - xlim[1]) * (i / (nrow(sector_counts) + 1))
                             # Plot dot for the RNAmodType
                             circos.points(x = x_pos, 
                                           y = sector_counts$counts[i],  # y is based on the count of RMPs
                                           pch = 16, 
                                           cex = sector_counts$counts[i] / max(rmp_counts$counts) * 2,  # Size reflects the counts
                                           col = rna_mod_colors[sector_counts$RNAmodType[i]])  # Color per RNAmodType
                           }
                         }
                       })

dev.off()

output_file <- "/data2/rluo4/RPMfunc/Output/scRNA/figures/circos_legend.png" 
# Set up the PNG device and define dimensions (e.g., 1000x1000 pixels)
png(output_file, width = 5, height = 3, res = 500,units = 'in')
circos.par("track.height" = 0.1)
circos.initialize(df.cir$sectors, x = df.cir$x)
set.seed(123)
# Define legend details
sizes <- c(1, 2, 3.5)  # Sizes of the points to be represented in the legend
size_labels <- c(1, 11, 40) #c("Small", "Medium", "Large")  # Labels for size categories
# Plot the color legend
legend("center",                            # Position legend in the top (adjust as needed)
       legend = names(rna_mod_colors),   # RNAmodType labels
       col = rna_mod_colors,             # Corresponding colors
       pch = 16,                         # Use solid circles for the legend
       pt.cex = rep(sizes, each = length(rna_mod_colors)), # Set size
       ncol = 3,                         # Make the legend horizontal
       cex = 0.9,                        # Text size for the legend
       bty = "n",                        # No border around the legend
       xpd = TRUE,                       # Allow drawing outside the plot region
       inset = c(-0.1, 0.1),             # Move the legend further away from the plot
       y.intersp = 1,
       x.intersp = 1.5)                    # Increase space between symbol and labels
# Plot the size legend with more space between circles and labels
legend("top",                         # Position legend in the bottom (adjust as needed)
       legend = size_labels,             # Labels for size categories
       pch = 16,                         # Use solid circles for the legend
       pt.cex = sizes,                   # Set sizes for the legend
       ncol = 3,    
       cex = 1.2,                        # Text size for the legend
       xpd = TRUE,                       # Allow drawing outside the plot region
       bty = "n",                        # No border around the legend
       y.intersp = 1.5,
       x.intersp = 1)                    # Increase space between symbol and labels
dev.off()

# legend("bottom",                            # Position legend below the plot
#        legend = names(rna_mod_colors),      # RNAmodType labels
#        col = rna_mod_colors,                # Corresponding colors
#        pch = 16,                            # Use solid circles for the legend
#        xpd = TRUE,                          # Allow drawing outside the plot region
#        inset = c(-0.1, 0.1),                # Move the legend farther below (adjust -0.2 for more distance)
#        ncol = length(rna_mod_colors),       # Make the legend horizontal by setting ncol equal to the number of items
#        cex = 0.8)    
# circos.link("NKcells", c(-0.5, 0.5), "OtherImmunecells", c(-0.5,0.5),  border = "#f6d573",
#             col = "#f89e81", h = 0.2)
# circos.link("Neutrophils", c(-0.5, 0.5), "OtherImmunecells", c(-0.5,0.5), border = "#f6d573",
#             col = "#f89e81", h = 0.2)
# circos.link("MastCells", c(-0.5, 0.5), "OtherImmunecells", c(-0.5,0.5),  border = "#f6d573",
#             col = "#f89e81", h = 0.2)
# circos.link("Myeloids", c(-0.5, 0.5), "OtherImmunecells", c(-0.5,0.5), border =  "#f6d573",
#             col = "#f89e81", h = 0.2)
# circos.link("Fibroblasts", c(0, 0), "OtherMesenchymals", c(0,0.1), border = "#acd485",
#             col = '#33A02C', h = 0.1)
# # circos.link("Tcells", 0, "Bcells", 0, h = 0.4, lwd = 2, lty = 2)



################################################################################
# 8.2. Summary of RNA modifications nmf proportion across healthy vs disease   #
################################################################################
# Load required packages
library(dplyr)
library(ggplot2)
library(ggpubr)
NMF_colors <- c("#f89e81","#99a9cc","#dd9bc5","#84c7b3","#acd485","#1B9E77","#f6d573","orange1",
                '#CAB2D6','#6A3D9A','#FFFF99','#B15928', '#FB9A99','#E31A1C','#FDBF6F','#FF7F00',
                '#A6CEE3','#1F78B4','#33A02C', "#66A61E","skyblue","salmon",'#B2DF8A'
)
outdir = '/data2/rluo4/RPMfunc/Output/scRNA/figures/NMFcluster_bars'
dir.create(outdir)
# Function to calculate proportion, generate plot, and perform statistical tests
generate_plots_per_cluster <- function(ct_nmf, plot_dir) {
  # Get unique clusters
  unique_clusters <- unique(ct_nmf$clusters)
  # Iterate through each cluster
  for (cluster in unique_clusters) {
    # Subset data for the current cluster
    # cluster = 'Stem_Cells'
    # Subset data for the current cluster
    subset_data <- ct_nmf %>%
      filter(clusters == cluster) 
    print(table(subset_data$disease))
    # Create the contingency table for actual frequencies
    contingency_table <- table(subset_data$disease, subset_data$NMF_clusters)
    # Perform the chi-square test on the contingency table
    chisq_test <- chisq.test(contingency_table)
    p_value <- chisq_test$p.value
    print(p_value)
    # Prepare data for plotting proportions
    plot_data <- subset_data %>%
      group_by(disease, NMF_clusters) %>%
      summarise(count = n(), .groups = "drop") %>%
      group_by(disease) %>%
      mutate(proportion = count / sum(count))  # Calculate proportions for visualization
    # Generate the plot with position "fill" to normalize proportions
    plot <- ggplot(plot_data, aes(x = disease, y = proportion, fill = NMF_clusters)) +
      geom_bar(stat = "identity", position = "fill", width = 0.7) +  # Fill bars to represent 100% proportions
      scale_y_continuous(labels = scales::percent_format()) +        # Display y-axis as percentages
      labs(title = paste("Cell Proportion Comparison in: ", cluster), y = "Proportion (%)", x = "",#"Disease", 
           fill = "NMF Clusters") +
      scale_fill_manual(values = NMF_colors) +  # Use custom color palette
      # theme_minimal() +
      stat_compare_means(method = "chisq.test", label = "p.format") +
      annotate("text", x = 1.5, y = 1.05, 
               label = paste("p-value:", format(p_value, digits = 3)),
               size = 5, color = "black") +
      guides(fill = guide_legend(title = "NMF Clusters")) +  # Improve legend readability
      theme(axis.ticks = element_line(color = "black"),
            # panel.border = element_rect(color = "black"),
            # panel.grid.minor = element_blank(),
            # strip.background = element_blank(),
            legend.position = "right", #'none',
            legend.key.width = unit(0.3, "cm"),
            # legend.box = "vertical",
            # legend.direction = "vertical",
            # plot.title = element_blank(),
            # strip.placement = "outside",
            strip.text.x = element_text(size = 14),
            axis.text.x = element_text(angle = 40, hjust = 1),
            axis.text.y = element_text(face = "italic"),
            axis.text = element_text(color = "black", size = 12),
            plot.margin = unit(c(0,1,0,0), "cm"),
            panel.grid.major = element_line(size=0.25, colour="grey80", linetype = "dashed")) 
    # labs(size = "P value (-log10)") +
    # guides(size = guide_legend(override.aes = list(pch = 16), title.position = "top", title.hjust = 0.5)) 
    # Print the plot
    print(plot)
    # Save the plot as an image (optional)
    if( nchar(cluster) > 10){
      wid = 9
    }else{
      wid = 7
    }
    ggsave(filename = file.path(plot_dir, paste0(cluster, "_proportion.png")),
           plot = plot, width = wid, height = 6)
  }
}

# Apply the function to your dataset
for (tissue in unique(NMF_clusters$ddt_tis)) {
  print(tissue)
  tis_dir <- file.path(outdir, tissue)
  if(! dir.exists( tis_dir )) {
    dir.create(tis_dir)
  }
  tis_nmf <- NMF_clusters[ NMF_clusters$ddt_tis %in% tissue,]
  
  for (ddt_id in unique(tis_nmf$ddt_id)) {
    ct_nmf <- tis_nmf[ tis_nmf$ddt_id == ddt_id, ]
    print(table(ct_nmf$disease))
    table(ct_nmf$clusters)
    table(ct_nmf$split_ids)
    plot_dir <- file.path(tis_dir, ddt_id)
    
    if(! dir.exists( plot_dir )) {
      dir.create(plot_dir)
    }
    
    if( NA %in% ct_nmf$split_ids ) {
      generate_plots_per_cluster(ct_nmf, plot_dir)
    } else{
      for (split_id in unique(ct_nmf$split_ids)) {
        ct_nmf_split <- ct_nmf[ ct_nmf$split_ids == split_id, ]
        print(table(ct_nmf_split$disease))
        table(ct_nmf_split$clusters)
        plot_dir_split <- file.path(plot_dir, split_id)
        if(! dir.exists( plot_dir_split )) {
          dir.create(plot_dir_split)
        }
        generate_plots_per_cluster(ct_nmf = ct_nmf_split, plot_dir = plot_dir_split)
      }
    }
  }
}

View(DEM_summary[DEM_summary$GeneSymbol=='FOXO1',])
View(DEM_summary[DEM_summary$GeneSymbol=='ANXA1',])
DT <- as.data.frame(table(sc_allmeta$disease, sc_allmeta$tissue))
DT <- DT[DT$Freq!=0,]
index = grepl('ALKBH5', NMF_clusters$RMP) & grepl('leukemia', NMF_clusters$disease)
# index = grepl('ALKBH5', GSE185381_allmeta$RMP) & grepl('leukemia', GSE185381_allmeta$disease)
index = grepl('WTAP', NMF_clusters$RMP) & grepl('lung', NMF_clusters$disease)
index = grepl('ZCCHC4', NMF_clusters$RMP) & grepl('liver', NMF_clusters$disease)
# index = grepl('METTL14', NMF_clusters$RMP) & grepl('endometrial carcinoma', NMF_clusters$disease)
# index = grepl('METTL3', NMF_clusters$RMP) & grepl('acute myeloid leukemia', NMF_clusters$disease)
table(NMF_clusters$ddt_id[index])
table(NMF_clusters$ddt_id[index], NMF_clusters$sub_type[index])
# Example for scRNA-seq data analysis:
GSE185381_allmeta <- subset(sc_allmeta, ddt_id == 'GSE185381')
table(GSE185381_allmeta$split_ids)
GSE185381_allmeta <- subset(GSE185381_allmeta, split_ids == 'sc_adult_AML')
sc_adult_AML <- unique(GSE185381_allmeta$patient_id)

tissue = 'bone marrow'
ddt_id = 'GSE185381'
# ddt_id = 'GSE139369'
split_id = 'sc_adult_AML'
indir = '/data2/rluo4/EpiTrans/DataCollection/RMzyme_scRNA'
# indir = '/data2/rluo4/EpiTrans/DataCollection/GSE185381/GSE185381_sc_stem_cell/sc_adult_AML/Stem_Cells'
cluster = 'Stem_Cells'
# cluster = 'Tcells'
nmffile = file.path(indir, tissue, ddt_id, split_id,  cluster, paste0(cluster, '_gse.nmf.rds'))
nmf <- readRDS(nmffile)

DEGfile = file.path(indir, tissue, ddt_id, split_id,  cluster, paste0(cluster, '_gse.nmf_healthy_disease_DEGs.rds'))
# DEGfile = file.path(indir, 'Stem_Cells_gse.nmf.rds')#'Stem_Cells_gse.nmf_healthy_disease_DEGs.rds')
DEGs <- readRDS(DEGfile)
DEGfile = file.path(indir, tissue, ddt_id, split_id,  cluster, paste0(cluster, '_gse.nmf_RMP_celltype_DEGs.rds'))
Stem_Cells <- readRDS(DEGfile)
table(Stem_Cells$cluster)
ALKBH5_DEGs <- Stem_Cells %>% filter(cluster == 'ALKBH5+Stem_Cells')
bulk_result <- fread( "/data/rluo4/RPMfunc/Output/summary/bulk_result.txt")
NOMO1 <- subset(bulk_result, Dataset == 'GSE144984')
table(NOMO1$Condition)
NOMO1_markers <- intersect(NOMO1$GeneSymbol, ALKBH5_DEGs$gene)
length(NOMO1_markers)/nrow(ALKBH5_DEGs)
# NOMO1 <- subset(NOMO1, Condition == 'NOMO1_ALKBH5OE_vs_WT')
NOMO1 <- subset(NOMO1, Condition == 'NOMO1_ALKBH5KD_vs_WT')
markers <- intersect(NOMO1$GeneSymbol, ALKBH5_DEGs$gene)
length(markers)/nrow(ALKBH5_DEGs)
NOMO1_DEG <- NOMO1[ is.na(NOMO1$Peak_id), ]
NOMO1_DEM <- NOMO1[! is.na(NOMO1$Peak_id) & NOMO1$Regulation=='DOWN', ]
NOMO1_DEM <- NOMO1[! is.na(NOMO1$Peak_id) & NOMO1$Regulation=='UP', ]
length(unique(NOMO1$GeneSymbol))
markers <- intersect(NOMO1_DEG$GeneSymbol, ALKBH5_DEGs$gene)
markers <- intersect(NOMO1_DEM$GeneSymbol, markers)
markers
# "EMP3", "SZRD1", "UBE2G2" --> 66
# View(DEGs[DEGs$gene %in% markers,])
intersect(markers, DEGs$gene)

KEGGfile = file.path(indir, tissue, ddt_id,split_id,  cluster, paste0(cluster, '_KEGG.rds'))
# library(KEGG.db)
KEGG <- readRDS(KEGGfile)
# all_path = as.data.frame(KEGG.db::KEGGPATHID2NAME)
# KEGG@compareClusterResult$Description = all_path$path_name[match(KEGG@compareClusterResult$ID,
#                                                                         all_path$path_id)]
dotplot(KEGG) + 
  theme(axis.text.x = element_text(
    angle = 45,
    vjust = 1, hjust = 1
  ))
# Enhanced dotplot
# Extract data for plotting
# Extract the top 10 pathways for each cluster based on the lowest p-value
plot_data <- KEGG@compareClusterResult %>%
  # mutate( pvalue = round(pvalue, 3)) %>%
  # mutate(logPval = sapply(pvalue, function(x) -log10(x))) %>%
  mutate(GeneRatio = sapply(GeneRatio, function(x) eval(parse(text = x)))) %>%
  group_by(Cluster) %>%
  arrange(pvalue) %>%
  slice_head(n = 6) %>%  # Take top 10 pathways per cluster
  ungroup()
summary(plot_data$pvalue)
# Prepare ggplot
KEGG_p <- ggplot(plot_data, aes(x = Cluster, y = reorder(Description,  GeneRatio))) +
  geom_point(aes(size = GeneRatio, color = pvalue)) +  # Size by GeneRatio, color by p-value
  # scale_color_gradient(high = "#84c7b3", low = "salmon", name = "P-value") +  # Color gradient based on p-value
  scale_color_gradientn(
    colors = c( '#CAB2D6', "#99a9cc", "#1B9E77"),  # Low, middle, and high colors
    values = scales::rescale(c(0, mean(plot_data$pvalue), max(plot_data$pvalue))),  # Rescale within the small p-value range
    name = "P Value") +
  scale_size(range = c(2, 6), name = "Gene Ratio") +  # Adjust size range of dots for better visibility
  labs(title = '',#"KEGG Pathway Enrichment Analysis",
       y = "", x = "") +  # Customize labels
  theme_minimal(base_size = 14) +  # Clean and professional theme
  theme(axis.ticks = element_line(color = "black"),
        # panel.border = element_rect(color = "black"),
        # panel.grid.minor = element_blank(),
        # strip.background = element_blank(),
        legend.position = "right", #'none',
        legend.key.width = unit(0.3, "cm"),
        # legend.box = "vertical",
        # legend.direction = "vertical",
        # plot.title = element_blank(),
        # strip.placement = "outside",
        strip.text.x = element_text(size = 14),
        axis.text.x = element_text(angle = 30, hjust = 1),
        axis.text.y = element_text(face = "italic"),
        axis.text = element_text(color = "black", size = 12),
        plot.margin = unit(c(0,1,0,0), "cm"),
        panel.grid.major = element_line(size=0.25, colour="grey80", linetype = "dashed")) 

tis_dir <- file.path(outdir, tissue)
plot_dir <- file.path(tis_dir, ddt_id)
plot_dir_split = file.path(plot_dir, split_id)
ggsave(filename = file.path(plot_dir_split, paste0(cluster, "_KEGG.png")),
       plot = KEGG_p,  width = 11, height = 8)

CCIfile = file.path(indir, tissue, ddt_id,split_id,  cluster, paste0(cluster, '_CCI_pairLR.use.rds'))
CCI <- readRDS(CCIfile)
intersect(NOMO1_markers, CCI$Receptor)
intersect(NOMO1_markers, CCI$Ligand)
CCIfile = file.path(indir, tissue, ddt_id,split_id,  cluster, paste0(cluster, '_aggregateNet_cellchat.rds'))
cellchat <- readRDS(CCIfile)
print(table(cellchat@idents))
signaling.name <- cellchat@netP$pathways
if (is.null(signaling.name)) {
  signaling.name <- signaling
}
print(signaling.name)
data.cci=NULL
data.cci <- subsetCommunication(cellchat)
data.cci <- data.cci[, c(1:7, 9)]#data.cci[,c(3,4,1,2)]
LR_genes <- strsplit(as.character(data.cci$interaction_name), '_')
LR_genes <- sort(unique(unlist(LR_genes)))
inter_genes <- intersect(LR_genes, ALKBH5_DEGs$gene)
intersect(inter_genes, DEGs$gene)
idx <- str_detect(as.character(data.cci$interaction_name), str_c(inter_genes, collapse = "|"))
View(data.cci[idx,])
inter_genes <- intersect(LR_genes, DEGs$gene)


# Define the base directory
indir <- '/data2/rluo4/EpiTrans/DataCollection/RMzyme_scRNA'
# Recursively list all files in the base directory
all_files <- list.files(indir, pattern = "_gse.nmf_.*_DEGs\\.rds$", full.names = TRUE, recursive = TRUE)
# Separate into two lists based on pattern
healthy_disease_DEGs <- grep("_gse.nmf_healthy_disease_DEGs\\.rds$", all_files, value = TRUE)
RMP_celltype_DEGs <- grep("_gse.nmf_RMP_celltype_DEGs\\.rds$", all_files, value = TRUE)
# Load the files if needed
healthy_disease_list <- lapply(healthy_disease_DEGs, readRDS)
names(healthy_disease_list) <- gsub('/data2/rluo4/EpiTrans/DataCollection/RMzyme_scRNA', '', healthy_disease_DEGs)
RMP_celltype_list <- lapply(RMP_celltype_DEGs, readRDS)
names(RMP_celltype_list) <- gsub('/data2/rluo4/EpiTrans/DataCollection/RMzyme_scRNA', '', RMP_celltype_DEGs)

# Check the results
print(healthy_disease_DEGs)
print(RMP_celltype_DEGs)



library(dittoSeq)
DT <- as.data.frame(table(NMF_clusters$disease, NMF_clusters$tissue))
DT <- DT[DT$Freq!=0,]
ct_nmf <- NMF_clusters %>% filter( ddt_tis == 'bone marrow', ddt_id == "GSE185991", #"GSE185381",#'E-MTAB-9139',
                                   disease == c('acute myeloid leukemia', 'healthy')#, split_ids == "DX", #'sc_large_tissue_split' #'cite_GSE185381_adult_AML' 
)
table(ct_nmf$clusters)
table(ct_nmf$split_ids)
ct <- 'Bcells'
# ct_nmf <- NMF_clusters %>% filter( tissue == 'liver',  disease_dataset_id == 'GSE156625', #clusters == ct,
#                                    disease %in% c('hepatocellular carcinoma', 'healthy') )
table(ct_nmf$disease_dataset_id, ct_nmf$disease)
table(ct_nmf$project_id, ct_nmf$disease)
table(ct_nmf$clusters, ct_nmf$disease)

proportion_data <- ct_nmf %>%
  group_by(disease, clusters, NMF_clusters) %>%
  summarise(cell_count = n()) %>%
  group_by(disease, clusters) %>%
  mutate(proportion = cell_count / sum(cell_count))

# Create a bar plot to visualize proportions across disease states for each cluster
ggplot(proportion_data, aes(x = disease, y = proportion, fill = NMF_clusters)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ clusters, scales = "free") +
  theme_minimal() +
  labs(title = "Cell Proportion by NMF Clusters across Disease States",
       x = "Disease State", y = "Proportion of Cells") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))




# organize the diseased data of disco named as disease_sc_summary
extract_first_element <- function(x) {
  split_string <- strsplit(x, "[--]")
  first_element <- sapply(split_string, function(y) head(y, n = 1))
  return(first_element)
}
sc_disease_metadata <- sc_metadata[! grepl('adjacent|ealthy', sc_metadata$disease),] # 219 disease datasets
sc_disease_metadata <- sc_disease_metadata[ ! sc_disease_metadata$project_id %in% tissue_summary$Dataset, ] # 169 --> 172 disease datasets
first_elements <- extract_first_element(sc_disease_metadata$barcode)
CB = str_split(first_elements,'[-]',simplify = T)[,1]
sc_disease_metadata$CB <- CB #str_split(sc_disease_metadata$barcode, '-', simplify = T)[,1]
all_characters <- unique(unlist(strsplit(sc_disease_metadata$ct, "")))
potential_separators_and_whitespace <- all_characters[grep("[[:punct:][:space:]]", all_characters)]
# Define the replacement character
replacement <- "_"
# Replace each potential separator with the replacement character
modified_vector <- sc_disease_metadata$ct
for (sep in potential_separators_and_whitespace) {
  modified_vector <- gsub(sep, replacement, modified_vector, fixed = TRUE)
}
# Replace consecutive occurrences of the replacement character with a single underscore
sc_disease_metadata$Cell_Type <- gsub("_+", replacement, modified_vector)

##################################################################
# 8.2. Summary of RNA modifications across multo-omics datasets  #
##################################################################
bulk <- bulk_result[, c(5, 4, 10, 9)]
colnames(bulk) <- colnames(disco_summary)
bulk <- bulk[!duplicated(paste(bulk$disease, bulk$tissue , bulk$project_id, bulk$patient_id)),]
bulk <- bulk[order(bulk$project_id),]
multissue <- bulk[grepl('and',bulk$tissue),]
bulk[grepl('and',bulk$tissue),'tissue'] <- c('Liver cancer', 'Cervical cancer', rep('Human embryonic kidney',2), rep('Human embryonic stem cell', 3),
                                             'Human embryonic kidney', rep('Human embryonic stem cell', 2), 'Breast cancer','Ocular melanoma',
                                             'Human embryonic microglia','Myelogenous leukemia','Colorectal adenocarcinoma and SARS-CoV-2',
                                             'Head and Neck Squamous Carcinoma',rep('Human embryonic kidney',2), rep('Lung cancer',5))
# bulk[grepl('and',bulk$tissue),'tissue'] <-  c('Cervical cancer','Liver cancer', 'Human embryonic kidney', rep('Human embryonic stem cell', 3),
#   'Human embryonic kidney', rep('Human embryonic stem cell', 2), 'Breast cancer','Ocular melanoma',
#   'Human embryonic microglia','Myelogenous leukemia','Colorectal adenocarcinoma and SARS-CoV-2',
#   'Head and Neck Squamous Carcinoma','Human embryonic kidney','Lung cancer')
cell_lines <- str_split(bulk$disease, '_', simplify = T)[,1]
cell_lines <- str_split(cell_lines, '[.]', simplify = T)[,1]
sort(unique(cell_lines))
bulk$omics <- 'bulk'
# CPTAC <- Proteogenomic_datasets[, c(2,3,1,4)]
# CPTAC$Normal <- 'CPTAC'
# colnames(CPTAC) <- colnames(disco_summary)
# CPTAC$omics <- 'CPTAC'
CPTAC <- Proteogenomic_datasets[, c(2,3,1,4)]
CPTAC$Normal <- 'CPTAC'
colnames(CPTAC) <- colnames(disco_summary)
CPTAC$omics <- CPTAC$disease #'CPTAC'
CPTAC$disease <- CPTAC$tissue
##################################################
# # df_list <- list(disco_summary,PCTanno,bulk,CPTAC)
# # df_list <- list(sc_summary,bulk,CPTAC)
# df_list <- list(sc_samples,bulk,CPTAC)
# # Using Reduce and merge for full outer join
# multiomics_summary <- Reduce(function(x, y) merge(x, y, all = TRUE), df_list)
# multiomics_summary$tissue <- gsub(' cancer', '', multiomics_summary$tissue)
# tissue <- sort(unique(multiomics_summary$tissue))
# multiomics_summary$tissue <-  gsub('Acute meyloid leukemia|Acute leukemia|Acute monocytic leukemia|Acute myeloid leukemia|Chronic lymphocytic leukemia|Chronic myeloid leukemia|Leukemia|Myelogenous leukemia', 'bone marrow',
#                                    gsub('Breast and Brain metastases|Breast|BRCA', 'breast',
#                                         gsub('Cervical', 'cervix',
#                                              gsub('Cervical and Liver', 'cervix and liver',
#                                                   gsub('CCRCC|Human embryonic kidney|Normal kidney', 'kidney',     
#                                                        gsub('Colorectal adenocarcinoma and SARS-CoV-2|COAD|Colon|Colorectum', 'colon',
#                                                             gsub('Diffuse large B-cell lymphoma|Erythroleukemia', 'blood',
#                                                                  gsub('PDAC', 'pancreas',
#                                                                       gsub('OV|HGSC', 'ovary',
#                                                                            gsub('Myelogenous leukemia and Human embryonic microglia', 'bone marrow and brain',     
#                                                                                 gsub('Embryo placenta', 'placenta',
#                                                                                      gsub('Endometrium|UCEC|Endometrial', 'endometrium',     
#                                                                                           gsub('GBM|Glioblastoma|Giloma|Human embryonic microglia', 'brain',  
#                                                                                                gsub('HGSC|OV', 'ovary',
#                                                                                                     gsub('LSCC|LUAD|Embryonic lung fibroblasts', 'lung',
#                                                                                                          gsub('HNSCC|Oral cavity|Head and Neck Squamous Carcinoma', 'head and neck', 
#                                                                                                               gsub('Skin cutaneous melanoma|Skin|Human keratinocytes|Foreskin|Melanoma', 'skin',  
#                                                                                                                    gsub("Kaposi's sarcoma", 'intestine',  
#                                                                                                                         gsub("Human Umbilical Vein Endothelial Cell|Human umblilical cord blood", 'blood',  
#                                                                                                                              gsub('Nasopharyngeal carcinoma','nasopharynx',        
#                                                                                                                                   gsub('Uveal melanoma|Ocular melanoma', 'eye',
#                                                                                                                                        gsub('T2D islets', 'pancreatic islet', 
#                                                                                                                                             gsub('Gastric', 'stomach',
#                                                                                                                                                  gsub('Esophageal squamous cell carcinoma', 'esophagus',
#                                                                                                                                                       gsub('Human embryonic stem cell|Human embryonic stem cells', 'endoderm',
#                                                                                                                                                            multiomics_summary$tissue,         
#                                                                                                                                                       )))))))))))))))))))))))))
# multiomics_summary$tissue <- tolower(multiomics_summary$tissue)
# tissue <- sort(unique(multiomics_summary$tissue))
# tolower(unique( disco$tissue))
# setdiff(tolower(unique( disco$tissue)), tissue)
# tissue
# length(unique(multiomics_summary$project_id))
# # [1] 499 --> 468
# unique(multiomics_summary$tissue[multiomics_summary$omics=='scRNA'])
##################################################
CPTAC$disease <-  gsub('Acute meyloid leukemia|Acute meyloid leukemia|Acute monocytic leukemia|Acute myeloid leukemia|Chronic lymphocytic leukemia|Chronic myeloid leukemia|Leukemia|Myelogenous leukemia', 'bone marrow',
                       gsub('Breast and Brain metastases|Breast|BRCA', 'breast',
                            gsub('Cervical', 'cervix',
                                 gsub('Cervical and Liver', 'cervix and liver',
                                      gsub('CCRCC|Human embryonic kidney|Normal kidney', 'kidney',     
                                           gsub('Colorectal adenocarcinoma and SARS-CoV-2|COAD|Colon|Colorectum', 'colon',
                                                gsub('Diffuse large B-cell lymphoma|Erythroleukemia', 'blood',
                                                     gsub('PDAC', 'pancreas',
                                                          gsub('OV|HGSC', 'ovary',
                                                               gsub('Myelogenous leukemia and Human embryonic microglia', 'bone marrow and brain',     
                                                                    gsub('Embryo placenta', 'placenta',
                                                                         gsub('Endometrium|UCEC|Endometrial', 'endometrium',     
                                                                              gsub('GBM|Glioblastoma|Giloma|Human embryonic microglia', 'brain',  
                                                                                   gsub('HGSC|OV', 'ovary',
                                                                                        gsub('LSCC|LUAD|Embryonic lung fibroblasts', 'lung',
                                                                                             gsub('HNSCC|Oral cavity|Head and Neck Squamous Carcinoma', 'head and neck', 
                                                                                                  gsub('Skin cutaneous melanoma|Skin|Human keratinocytes|Foreskin|Melanoma', 'skin',  
                                                                                                       gsub("Kaposi's sarcoma", 'intestine',  
                                                                                                            gsub("Human Umbilical Vein Endothelial Cell|Human umblilical cord blood", 'blood',  
                                                                                                                 gsub('Nasopharyngeal carcinoma','nasopharynx',        
                                                                                                                      gsub('Uveal melanoma|Ocular melanoma', 'eye',
                                                                                                                           gsub('T2D islets', 'pancreatic islet', 
                                                                                                                                gsub('Gastric', 'stomach',
                                                                                                                                     gsub('Esophageal squamous cell carcinoma', 'esophagus',
                                                                                                                                          gsub('Human embryonic stem cell|Human embryonic stem cells', 'endoderm',
                                                                                                                                               CPTAC$tissue,         
                                                                                                                                          )))))))))))))))))))))))))
write.table(CPTAC, "/data/rluo4/RPMfunc/Output/summary/cptac_summary.txt", row.names = FALSE, sep = '\t', quote = F)

library(RColorBrewer)  
library(reshape2)
library(ggplot2)
library(ggprism)
# 1. Statistics for RMPs of different RNA modifications: /data/rluo4/RPMfunc/WebSite/RNA_modification_type.png 
dfm <- as.data.frame(table(RMP_update$RNA_modifications))#tissue_summary$Tissue))
colnames(dfm) <- c('RM_type', 'Freq')#c('Tissue','Freq')
dfm <- dfm[order(dfm$Freq),]
RColorBrewer::brewer.pal(n = 12, name = 'Paired')
color_cluster=c("#1B9E77","#66A61E","skyblue","#99a9cc","#f89e81","#acd485","#dd9bc5","#f6d573","#84c7b3")
color_cluster=c(RColorBrewer::brewer.pal(n = 12, name = 'Paired'),color_cluster)
color_cluster <- color_cluster[1:length(unique(dfm$RM_type))]
names(color_cluster)=unique(dfm$RM_type)#c("B cells","Epithelial cells","Myeloid cells","Neurons","Stromal cells","T cells")
# color_cluster <- color_cluster[match(dfm$RM_type, names(color_cluster))]
ggplot(dfm, aes(x = Freq, y = RM_type))+#, group = group)) +
  geom_col(aes(fill=RM_type)) + ylab('') + xlab('No. of RNA Modifying proteins') +
  # theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.25)) +
  geom_text(aes(label = Freq), position = position_stack(vjust = .5), size = 6) + # labels inside the bar segments
  scale_fill_manual(values = color_cluster)+
  # scale_y_reverse("Datasets by Primary site",expand = c(0,0),position = 'left',limits=c(250,0))+#,position = "right",limits=c(5000,0))+  #yÖáÔÚÓÒ
  theme(axis.text=element_text(face = "bold", color="black")) + #labs(title = "Predicted neopeptides in TCGA-SKCM")+ 
  # coord_flip() +
  theme_prism(#palette = "winter_bright",
    base_fontface = "plain", # ×ÖÌåÑùÊ½£¬¿ÉÑ¡ bold, plain, italic
    base_family = "serif", # ×ÖÌå¸ñÊ½£¬¿ÉÑ¡ serif, sans, mono, ArialµÈ
    base_size = 17,  # Í¼ÐÎµÄ×ÖÌå´óÐ¡
    base_line_size = 0.8, # ×ø±êÖáµÄ´ÖÏ¸
    axis_text_angle = 0) + theme(legend.position='none')
ggsave(file = paste0('/data/rluo4/RPMfunc/WebSite/RNA_modification_type.png'),height = 6, width = 8, dpi = 500)

# 2. Overview of modified RMPs that contain multiple PTMs of the same type:/data/rluo4/RPMfunc/WebSite/PTMs.png
Freq_PTM <- rbind(Freq_database, Freq_PMID)
colnames(Freq_PTM)[1] <- 'PTMs'
dfm <- Freq_PTM %>% aggregate(Freq~PTMs, FUN="sum")
RColorBrewer::brewer.pal(n = 12, name = 'Paired')
colors <- paletteer::palettes_d
colors <- colors$wesanderson$BottleRocket1
color_cluster=c("orange3", "#1B9E77", "salmon", "#1B9E77","#66A61E","skyblue","#99a9cc","#f89e81","#acd485","#dd9bc5","#f6d573","#84c7b3")
color_cluster=c(colors,RColorBrewer::brewer.pal(n = 12, name = 'Paired'), color_cluster)
color_cluster <- color_cluster[1:length(unique(dfm$PTMs))]
names(color_cluster)=unique(dfm$PTMs)
ggplot(dfm,aes(x = Freq, y = PTMs))+
  geom_col(aes(fill=PTMs)) + ylab('') + xlab('No. of PTM sites') +
  # theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.25)) +
  geom_text(aes(label = Freq), position = position_stack(vjust = 1.05), size = 4) + # labels inside the bar segments
  # scale_colour_gradientn(colors =  'salmon') + # scale_colour_continuous(dfm$Freq, type = getOption("ggplot2.continuous.colour")) +
  scale_fill_manual(values = color_cluster)+ #scale_fill_brewer(palette = "Blues")+
  # scale_y_reverse("Datasets by Primary site",expand = c(0,0),position = 'left',limits=c(250,0))+#,position = "right",limits=c(5000,0))+  #yÖáÔÚÓÒ
  theme(axis.text=element_text(face = "bold", color="black")) + #labs(title = "Predicted neopeptides in TCGA-SKCM")+ 
  coord_flip() +
  theme_prism(#palette = "winter_bright",
    base_fontface = "plain", # ×ÖÌåÑùÊ½£¬¿ÉÑ¡ bold, plain, italic
    base_family = "serif", # ×ÖÌå¸ñÊ½£¬¿ÉÑ¡ serif, sans, mono, ArialµÈ
    base_size = 16,  # Í¼ÐÎµÄ×ÖÌå´óÐ¡
    base_line_size = 0.8, # ×ø±êÖáµÄ´ÖÏ¸
    axis_text_angle = 45) + theme(legend.position='none')
ggsave('/data/rluo4/RPMfunc/WebSite/PTMs.png', height = 5, width = 11, dpi = 500)

# 3. Overview of major cell types for scRNA-seq datasets in RMzyme.
DT <- as.data.frame(table(sc_metadata$disease, sc_metadata$tissue))
DT <- DT[DT$Freq!=0,]
library(ArchR)
bar.df <- sc_metadata
title <- paste0('Major cell type fraction (%)')
set.seed(123)
colors <- paletteer::palettes_d
colors <- colors$wesanderson$BottleRocket1
color_cluster=c("orange3", "#1B9E77", "salmon", "#f6d573","#99a9cc","skyblue","#66A61E","#f89e81","#acd485","#dd9bc5")
color_cluster=c(RColorBrewer::brewer.pal(n = 12, name = 'Paired'), color_cluster)
color_cluster <- color_cluster[1:length(unique(bar.df$sub_type))]
names(color_cluster)=unique(bar.df$sub_type)
# color_cluster = paletteDiscrete(values = unique(bar.df$Stage))
p <- ggplot(bar.df,aes(x=tissue))+
  geom_bar(aes(fill=sub_type),position = "fill",width = .8)+
  scale_x_discrete("")+
  scale_y_continuous(title,expand = c(0,0),labels = scales::label_percent(),position = "right")+
  scale_fill_manual("",values = color_cluster)+
  theme_ArchR(baseSize = 10) + #ggtitle('') +
  # theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())+
  theme(
    plot.title    = element_text(color = 'black', size   = 16, hjust = 0.5),
    plot.subtitle = element_text(color = 'black', size   = 16,hjust = 0.5),
    plot.caption  = element_text(color = 'black', size   = 16,face = 'italic', hjust = 1),
    axis.text.x   = element_text(color = 'black', size = 12, angle = 0),
    axis.text.y   = element_text(color = 'black', size = 11, angle = 0),
    axis.title.x  = element_text(color = 'black', size = 12, angle = 0),
    axis.title.y  = element_text(color = 'black', size = 12, angle = 90),
    legend.title  = element_text(color = 'black', size  = 10),
    legend.text   = element_text(color = 'black', size   = 10),
    axis.line.y = element_line(color = 'black', linetype = 'solid'), # y轴线特征
    axis.line.x = element_line (color = 'black',linetype = 'solid'), # x轴线特征
    # panel.border = element_rect(linetype = 'solid', size = 1.2,fill = NA) # 图四周框起来
  )+   coord_flip() #让条形图横过来
# hei <- ceiling(length(unique(bar.df$Tissue))*1.5)
f = paste0('/data/rluo4/RPMfunc/WebSite/scRNA_subtype.png')
p
ggsave(plot=p,width = 10,height =12, filename=f, dpi = 500, device = "png")

# 4. Overview of multiomics datasets by tissue type collected in RMzyme:/data/rluo4/RPMfunc/WebSite/Omics.png
bulk <- bulk_result[, c(5, 4, 10, 9)]
colnames(bulk) <- colnames(disco_summary)
bulk <- bulk[!duplicated(paste(bulk$disease, bulk$tissue, bulk$project_id, bulk$patient_id)),]
bulk <- bulk[order(bulk$project_id),]
multissue <- bulk[grepl('and',bulk$tissue),]
bulk[grepl('and',bulk$tissue),'tissue'] <- c('Liver cancer', 'Cervical cancer', rep('Human embryonic kidney',2), rep('Human embryonic stem cell', 3),
                                             'Human embryonic kidney', rep('Human embryonic stem cell', 2), 'Breast cancer','Ocular melanoma',
                                             'Human embryonic microglia','Myelogenous leukemia','Colorectal adenocarcinoma and SARS-CoV-2',
                                             'Head and Neck Squamous Carcinoma',rep('Human embryonic kidney',2), rep('Lung cancer',5))
cell_lines <- str_split(bulk$disease, '_', simplify = T)[,1]
cell_lines <- str_split(cell_lines, '[.]', simplify = T)[,1]
sort(unique(cell_lines))
bulk.celine_RNA <- bulk[bulk$project_id %in% DEG_summary$Dataset,]
bulk.celine_RNA$omics <- 'transcriptomics'
bulk.celine_epitrans <- bulk[bulk$project_id %in% DEM_summary$Dataset,]
bulk.celine_epitrans$omics <- 'epitranscriptomics'
bulk <- rbind(bulk.celine_epitrans,bulk.celine_RNA)

CPTAC <- Proteogenomic_datasets[, c(2,3,1,4)]
CPTAC$Normal <- 'CPTAC'
colnames(CPTAC) <- colnames(disco_summary)
CPTAC$omics <- CPTAC$disease #'CPTAC'
CPTAC$disease <- CPTAC$tissue
# df_list <- list(disco_summary,PCTanno,bulk,CPTAC)
df_list <- list(sc_samples,bulk,CPTAC) #df_list <- list(sc_summary,bulk,CPTAC)
# Using Reduce and merge for full outer join
multi_omics <- Reduce(function(x, y) merge(x, y, all = TRUE), df_list)
multi_omics$tissue <- gsub(' cancer', '', multi_omics$tissue)
tissue <- sort(unique(multi_omics$tissue))
multi_omics$tissue <- gsub('Acute meyloid leukemia|Acute leukemia|Acute monocytic leukemia|Acute myeloid leukemia|Chronic lymphocytic leukemia|Chronic myeloid leukemia|Leukemia|Myelogenous leukemia', 'bone marrow',
                           gsub('Breast and Brain metastases|Breast|BRCA', 'breast',
                                gsub('Cervical', 'cervix',
                                     gsub('Cervical and Liver', 'cervix and liver',
                                          gsub('CCRCC|Human embryonic kidney|Normal kidney', 'kidney',     
                                               gsub('Colorectal adenocarcinoma and SARS-CoV-2|COAD|Colon|Colorectum', 'colon',
                                                    gsub('Diffuse large B-cell lymphoma|Erythroleukemia', 'blood',
                                                         gsub('PDAC', 'pancreas',
                                                              gsub('OV|HGSC', 'ovary',
                                                                   gsub('Myelogenous leukemia and Human embryonic microglia', 'bone marrow and brain',     
                                                                        gsub('Embryo placenta', 'placenta',
                                                                             gsub('Endometrium|UCEC|Endometrial', 'endometrium',     
                                                                                  gsub('GBM|Glioblastoma|Giloma|Human embryonic microglia', 'brain',  
                                                                                       gsub('HGSC|OV', 'ovary',
                                                                                            gsub('LSCC|LUAD|Embryonic lung fibroblasts', 'lung',
                                                                                                 gsub('HNSCC|Oral cavity|Head and Neck Squamous Carcinoma', 'head and neck', 
                                                                                                      gsub('Skin cutaneous melanoma|Skin|Human keratinocytes|Foreskin|Melanoma', 'skin',  
                                                                                                           gsub("Kaposi's sarcoma", 'intestine',  
                                                                                                                gsub("Human Umbilical Vein Endothelial Cell|Human umblilical cord blood", 'blood',  
                                                                                                                     gsub('Nasopharyngeal carcinoma','nasopharynx',        
                                                                                                                          gsub('Uveal melanoma|Ocular melanoma', 'eye',
                                                                                                                               gsub('T2D islets', 'pancreatic islet', 
                                                                                                                                    gsub('Gastric', 'stomach',
                                                                                                                                         gsub('Esophageal squamous cell carcinoma', 'esophagus',
                                                                                                                                              gsub('Human embryonic stem cell|Human embryonic stem cells', 'endoderm',
                                                                                                                                                   multi_omics$tissue,         
                                                                                                                                              )))))))))))))))))))))))))
multi_omics$tissue <- tolower(multi_omics$tissue)
multi_omics$omics <- gsub('scRNA', 'single-cell transcriptomics', multi_omics$omics)
write.table(multi_omics[,1:5], "/data/rluo4/RPMfunc/Output/summary/multiomics_summary.txt", row.names = FALSE, sep = '\t', quote = F)
dim(multi_omics)
unique(multi_omics$tissue)#72-->74-->72 kinds of tissues
bar.df <- multi_omics[!duplicated(paste(multi_omics$project_id, multi_omics$tissue, multi_omics$omics)),]#disease_summary
level <- unique(multi_omics$tissue)#unique(bar.df$Tissue)
bar.df <- mutate(bar.df,name=factor(bar.df$tissue, levels=level))
text.df <- as.data.frame(table(bar.df$name))
table(multi_omics$omics)
color_cluster=c("#99a9cc","#acd485","#f89e81","#1B9E77","#84c7b3","#dd9bc5","skyblue","#A42820","#f6d573" )
names(color_cluster)=unique(multi_omics$omics)#c('Human','Mus musculus')#c("B cells","Epithelial cells","Myeloid cells","Neurons","Stromal cells","T cells")
dfm = data.frame(table(bar.df$tissue,bar.df$omics))
dfm <- dfm[dfm$Freq!=0,]
colnames(dfm)[1:2] <- c('Tissue','Omics')#bar.df <- bar.df[order(bar.df$Organ,decreasing = T),]
library(ggprism)
ggplot(bar.df,aes(x=tissue))+
  geom_bar(aes(fill=omics),position = "stack",width = .6)+
  scale_x_discrete("",position = "top")+  #x轴在上
  scale_fill_manual(values = color_cluster)+ #ggtitle('Omics data')+
  scale_y_reverse("No. of Omics Datasets by Tissue",expand = c(0,0),position = 'left',limits=c(46,0))+#,position = "right",limits=c(5000,0))+  #y轴在右
  theme(axis.text=element_text(face = "bold", color="black")) + #labs(title = "Predicted neopeptides in TCGA-SKCM")+ 
  coord_flip() +
  theme_prism(#palette = "winter_bright",
    base_fontface = "plain", # 字体样式，可选 bold, plain, italic
    base_family = "serif", # 字体格式，可选 serif, sans, mono, Arial等
    base_size = 17,  # 图形的字体大小
    base_line_size = 0.8, # 坐标轴的粗细
    axis_text_angle = 0) + theme(legend.position='left')# 可选值有 0，45，90，270
# ggsave('/data/rluo4/RPMfunc/WebSite/Omics.png', height = 12, width = 10, dpi = 500)
ggplot(bar.df,aes(x=tissue))+
  geom_bar(aes(fill=omics),position = "stack",width = .7)+
  scale_x_discrete("",position = "bottom")+  ylab('No. of Omics Datasets by Tissue') + #x轴在上
  scale_fill_manual(values = color_cluster)+ #ggtitle('Omics data')+
  # scale_y_reverse("No. of Omics Datasets by Tissue",expand = c(0,0),position = 'left',limits=c(46,0))+#,position = "right",limits=c(5000,0))+  #y轴在右
  theme(axis.text=element_text(face = "bold", color="black")) + #labs(title = "Predicted neopeptides in TCGA-SKCM")+ 
  # coord_flip() +
  theme_prism(#palette = "winter_bright",
    base_fontface = "plain", # 字体样式，可选 bold, plain, italic
    base_family = "serif", # 字体格式，可选 serif, sans, mono, Arial等
    base_size = 14,  # 图形的字体大小
    base_line_size = 0.8, # 坐标轴的粗细
    axis_text_angle = 45) + theme(legend.position='top')# 可选值有 0，45，90，270
ggsave('/data/rluo4/RPMfunc/WebSite/Omics.png', height = 6, width = 15, dpi = 500)

datasets_info = final_datasets[, c(1:4, 8)]
colnames(datasets_info)[1] = 'Dataset'
datasets_info <- left_join(bulk_summary, datasets_info, by = 'Dataset')
index = paste( datasets_info$Dataset, datasets_info$Condition); length(unique(index))
multiomics_index = paste( multiomics_summary$project_id, multiomics_summary$patient_id )
length(unique(multiomics_index))
setdiff(unique(index), unique(multiomics_index))
table(index %in% multiomics_index)
datasets_info$Tissue <- multiomics_summary$tissue[match(index, multiomics_index)]
datasets_info$HyperLink <- paste0('https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=', datasets_info$Dataset)
write.table(datasets_info[, -c(1,7)], "/data/rluo4/RPMfunc/Output/summary/datasets_info.txt", row.names = FALSE, sep = '\t', quote = F)


#####################################################################################################
# 9) #   the RNA modification site is located in RBP binding/miRNA target/SNV/SNP site of Human     #
#####################################################################################################
################################################################################
# 9.1 A tab separated table of various of RNA modification sites related to RBP
################################################################################
# RBP_files <- list.files('/data/rluo4/EpiTrans/DataCollection', pattern = '.bed')
# RBP <- fread('/data/rluo4/EpiTrans/DataCollection/human.hg38.modrbp.m5C.reader.bed', sep = '\t')
# colN <- c("rbpSiteID", "rbpName", "narrowPeak", "broadPeak", "rbpRand", 
#           "clipExpNum", "dataSetIdList", "RNAmodID", "RNAmodLoc", 
#           "RNAmodType", "geneID", "transcriptID", "geneName", 
#           "geneType", "Region", "Conservation", "rbpType")
# # Read and process each count file
# RBP <- lapply(RBP_files, function(x) {
#   # Read data using fread for faster reading
#   data <- fread(file.path('/data/rluo4/EpiTrans/DataCollection', x))
#   # Extract required columns
#   # Set column names same as file names
#   colnames(data) <- colN# Display the first few rows with column names
#   return(data)
# })
# df_in <- Reduce(function(x, y) merge(x, y, by = colN, all = TRUE), RBP)
# # Define the output file
# output_file="merged_hg38_modrbp.bed"
# # Define the header line with column names
# header="rbpSiteID\trbpName\tnarrowPeak\tbroadPeak\trbpRand\tclipExpNum\tdataSetIdList\tRNAmodID\tRNAmodLoc\tRNAmodType\tgeneID\ttranscriptID\tgeneName\tgeneType\tRegion\tConservation\trbpType"
# # Create the output file and add the header
# echo -e $header > $output_file
# # Loop through each BED file and append its content to the output file
# for file in human.hg38.modrbp.*.bed; do
# echo "Processing $file..."
# # Append the content of the file, skipping the header if exists
# tail -n +1 $file >> $output_file
# done
# echo "All files have been merged into $output_file"
################################################################################
RBP_site <- fread('/data/rluo4/EpiTrans/DataCollection/merged_hg38_modrbp.bed')
table(RBP_site$RNAmodType)
RBP_Name <- unique(RBP_site$rbpName)
table(RBP_Name %in% RMP_update$RMP)
setdiff( RMP_update$RMP,  RBP_Name)
setdiff(bulk_summary$Enzyme,  RBP_Name)
intersect(Enzymes,  RBP_Name)
# [1] "YTHDF1"  "METTL3"  "ALKBH5"  "NSUN2"   "FTO"    
# [6] "IGF2BP1" "YTHDF3"  "METTL14" "YTHDF2"  "WTAP"   
# [11] "DKC1"    "YTHDC2"  "ALYREF"  "SAFB2"   "SND1"   
# [16] "YTHDC1" 
RBP_site$geneSymbol <- str_split(RBP_site$geneName, ',', simplify = T)[,1]
RBP_site$RNAmodType <- gsub('Y', 'Psi', RBP_site$RNAmodType)
# write.table(RBP_site[, c(1,2,3,8,9,10,17,18)], "/data/rluo4/RPMfunc/Output/summary/RBP_site.txt", row.names = FALSE, sep = '\t', quote = F)
length(unique(intersect(bulk_result$GeneSymbol, RBP_site$geneSymbol))) # 21494
unique(bulk_result$Enzyme)
unique(bulk_result$Modification_type) # 8 kinds
table(unique(bulk_result$Modification_type) %in% unique(RBP_site$RNAmodType) )
setdiff(unique(bulk_result$Modification_type), unique(RBP_site$RNAmodType) )
# RBP_site_less <- RBP_site[RBP_site$rbpName %in% Enzymes,]
# View(RBP_site[RBP_site$geneSymbol %in% c('ANXA1','FPR1'),])
# The default is hg38 and WTAP of m6A writer.
# ♥ The RBP Peak is the intersection peak that derived from different datasets.
# Parse narrowPeak column to extract chr, start, end
RBP_site[, `:=`(
  Chr = sub("^(chr[0-9XYM]+):.*", "\\1", narrowPeak),
  Start = as.numeric(sub("^chr[0-9XYM]+:([0-9]+)-[0-9]+:.*", "\\1", narrowPeak)),
  End = as.numeric(sub("^chr[0-9XYM]+:[0-9]+-([0-9]+):.*", "\\1", narrowPeak))
)]

# Save DEM_summary as a .bed file
DEM_summary.bed <- DEM_summary[, c("Chr", "Start", "End", "Peak_id", 'GeneSymbol','Dataset', 'Enzyme', 'Modification_type')]
DEM_summary.bed$hg_version <- genome_version$HG[match(DEM_summary.bed$Dataset, genome_version$GSE)] #"/data/rluo4/EpiTrans/DataCollection/RBP/DEM_summary.bed"
# for (i in unique(DEM_summary.bed$Dataset)) {
#     bed.path <- file.path('/data/rluo4/EpiTrans/DataCollection/RBP', paste0(i, '_DEM.bed'))
#     GSE.bed <- DEM_summary.bed[DEM_summary.bed$Dataset == i, 1:4]
#     GSE.bed$Peak_id <- gsub(' ', '-', GSE.bed$Peak_id)
#     write.table(GSE.bed, bed.path, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
# }
RBP_dir='/data/rluo4/EpiTrans/DataCollection/RBP'
for (n in 1:nrow(RMPdeg_datasets)) {
  print(n)
  print(RMPdeg_datasets[n, 'GSE'])
  i <- RMPdeg_datasets[n, 'GSE']
  print(i)
  if( genome_version$skip[genome_version$GSE==i] =='NoData' ){
    print(paste0(i, ": no data finally, skip  !"))
    next;
  }
  outdir <- file.path(RBP_dir, i)
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
    DEG_sub <- DEG_subset %>% filter(Condition == contrast) %>% arrange(desc(abs(log2FoldChange))) #%>%  mutate(row_id = row_number()) 
    DEM_sub <- DEM_subset %>% filter(Condition == contrast) %>% arrange(desc(abs(log2FoldChange))) #%>%  mutate(row_id = row_number()) 
    DEG_res <- head(DEG_sub, 200)#1000)
    DEM_res <- head(DEM_sub, 200)#1000)
    index1 <- DEG_res %>% nrow()
    index2 <- DEM_res %>% nrow()
    start_time <- Sys.time()
    if( i == 'GSE144984' & contrast=='NOMO1_ALKBH5OE_vs_WT'){
      ANXA1 = DEM_sub[DEM_sub$GeneSymbol=='ANXA1',]; DEM_res[133,] <- ANXA1 #AC211486.3
    }
    if( i == 'GSE94613' & contrast=='MOLM13_METTL3KD_vs_WT'){
      FPR1 = DEM_sub[DEM_sub$GeneSymbol=='FPR1',]; DEM_res[176,] <- FPR1 #AC133749.1
    }
    # print(paste("start GO analysis:", contrast, 'from', i, 'at', start_time))
    dir_name <- gsub('_vs_WT', '', contrast)
    outdir <- file.path(RBP_dir, i, dir_name)
    if(! dir.exists( outdir )) {
      dir.create(outdir)
    }
    print(outdir)
    inter = intersect(DEG_res$GeneSymbol, DEM_res$GeneSymbol)
    print(paste0('intersection of DEG and DEM: ',length(unique(inter))))
    
    if( index1 !=0 ){
      rownames(DEG_res) <- 1:index1
      # write.table(DEG_res, file = paste0(gsub('GEO/','GEO_top/',outdir), '/show_top_DEG.txt'), sep = '\t', quote = F, row.names = F)
    }else{
      print(paste0(i, ": no data in show_DEG.txt !"))
    }
    all_peak <- unique(DEM_subset$GeneSymbol)
    print(paste0('no. of DEM peaks: ',length(all_peak)))
    target_genes <- unique(DEM_res$GeneSymbol)
    print(paste0('no. of DEM targets: ',length(target_genes)))
    if( index2 !=0 ){
      rownames(DEM_res) <- 1:index2
      DEM_res <- DEM_res[, c("Chr", "Start", "End", "Peak_id", 'GeneSymbol','Dataset', 'Enzyme', 'Modification_type')]
      DEM_res$hg_version <- genome_version$HG[match(DEM_res$Dataset, genome_version$GSE)] #"/data/rluo4/EpiTrans/DataCollection/RBP/DEM_summary.bed"
      bed.path <- file.path(outdir, paste0('DEM.bed'))
      GSE.bed <- DEM_res[, 1:4]#DEM_summary.bed[DEM_summary.bed$Dataset == i, 1:4]
      GSE.bed$Peak_id <- gsub(' ', '-', GSE.bed$Peak_id)
      write.table(GSE.bed, bed.path, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
    }
  }
}
hg38_datasets <- unique(DEM_summary.bed$Dataset[DEM_summary.bed$hg_version=='hg38'])
paste0(hg38_datasets, collapse =  "' | '")
# Save RBP_site as a .bed file
RBP_site.bed <- RBP_site[, c("Chr", "Start", "End", "rbpSiteID", 'rbpName', 'rbpType', 'geneSymbol', 'RNAmodType')]
# write.table(RBP_site.bed, "/data/rluo4/EpiTrans/DataCollection/RBP_site.bed", 
#             sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
# cat RBP_site.bed | sortBed -i >RBP_site.sorted.bed
# Error: malformed BED entry at line 207. Start was greater than end. Exiting.
# more RBP_site.bed|sed -n '200, 207p'
# chr12   57166998        57167030        RBP_site_355205 LRP1    m6A
# chr5    74640673        74640716        RBP_site_236550 ENC1    m6A
# chrX    150763362       150763400       RBP_site_614845 MTMR1   m6A
# chr16   30765237        30765326        RBP_site_79288  RNF40   m6A
# chr11   61795976        61796029        RBP_site_490940 FADS2   m6A
# chr12   46204331        46204416        RBP_site_475354 SLC38A1 m6A
# chr22   30894012        30894153        RBP_site_181176 OSBP2   m6A
# chrM:15078-15128:+      NA      NA      RBP_site_357539 MT-CYB  m6A
dsets <- list.files(RBP_dir)
j <- NULL; k <- NULL
for (i in dsets) {
  print(i)
  file = list.files(file.path(RBP_dir, i), pattern = 'rbp.bed')
  j <<- c(j, file)
  # if(length(file) !=0 ){
  k <<- c(k, rep(i, length(file)))
  # }
}
bed_files <- data.frame(GSE=k, bed=j)
unique(bed_files$GSE)
colN <- c(colnames(DEM_summary.bed)[1:4], "chr", "start", "end", "rbpSiteID", 'rbpName', 'rbpType', 'geneSymbol', 'RNAmodType')
DEM_merge.bed <- apply(bed_files, 1, function(x) {
  # Read data using fread for faster reading
  print(x)
  data <- fread(file.path(RBP_dir, x['GSE'], x['bed']))
  print(head(data))
  print(dim(data))
  colnames(data) <- colN
  return(data)
})
names(DEM_merge.bed) <- paste0(bed_files$GSE, '_', gsub('_rbp.bed', '', bed_files$bed))

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
plot_dir <- '/data/rluo4/RPMfunc/Output/GEO_top/'
# for (i in names(DEM_merge.bed)) {
changes <- c('GSE94613_MOLM13_METTL3KD', 'GSE144984_NOMO1_ALKBH5KD', 'GSE144984_NOMO1_ALKBH5OE')
for (i in changes) {
  GSE <- str_split(i, '_', simplify = T)[,1]
  inter.bed <- DEM_merge.bed[[i]]
  GSE.bed <- DEM_summary.bed[DEM_summary.bed$Dataset == GSE, ]
  GSE.bed$Peak_id <- gsub(' ', '-', GSE.bed$Peak_id)
  print(length(GSE.bed$Peak_id))
  print(dim(GSE.bed))
  condition = str_split(rownames(GSE.bed), '_vs', simplify = T)[,1]
  index = condition %in% i
  table(index)
  GSE.bed <- GSE.bed[index,]
  inter.bed <- left_join( inter.bed, GSE.bed[, 4:9], by = 'Peak_id' )
  
  index = inter.bed$geneSymbol == inter.bed$GeneSymbol
  print(table(index))
  inter.bed <- inter.bed[index,]
  index = inter.bed$RNAmodType == inter.bed$Modification_type
  print(table(index))
  inter.bed <- inter.bed[index,]
  # RBP_site_805332 vs. RBP_site_572197
  rmdup <- paste(inter.bed$Peak_id, inter.bed$rbpSiteID, inter.bed$rbpName, inter.bed$geneSymbol, inter.bed$geneSymbol )
  index = !duplicated(rmdup)
  DEM_merge.bed[[i]] <-  inter.bed[index,]
  GSM <- unique( str_split(rownames(GSE.bed), '_vs', simplify = T)[,1])
  GSM <- remove_first_element(GSM)
  bed.path <- file.path(plot_dir, GSE, GSM,  paste0( 'show_top_rbp.txt'))
  bed <- DEM_merge.bed[[i]][, ]#c(5:12, 4, 14, 15)
  colnames(bed)
  print(dim(bed))
  print(unique(bed$RNAmodType))
  bed$rbpNarrowPeak <- paste0(bed$chr,":", bed$start, "-", bed$end)
  bed$RNAmodSite <- paste0(bed$Chr,":", bed$Start, "-", bed$End)
  # RMBase <- RBP_site[RBP_site$rbpSiteID %in% bed$rbpSiteID & RBP_site$RNAmodType %in% bed$RNAmodType, ]
  # print(dim(RMBase))
  # bed <- left_join(bed, RMBase, by = 'rbpSiteID')
  if(nrow(DEM_merge.bed[[i]]) !=0){
    write.table(bed[, c("rbpNarrowPeak","rbpSiteID","rbpName","rbpType","geneSymbol","RNAmodType","Dataset","Enzyme","RNAmodSite")], #c(18, 8:12, 14:15, 19)],
                bed.path, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
  }
}
# RBP_DEM.bed <- Reduce(function(x, y) merge(x, y, by = colN, all = TRUE), DEM_merge.bed)

############################################################################################
# 9.2 # A tab separated table of various of RNA modification sites related to miRNA target #
############################################################################################
# (the RNA modification site is located in miRNA target) among Homo sapien. The named format of file is "genome.mirTarget.tar.gz".
# The compressed file contains two files, that are mod-mirTar file and miRNA-target align file respectively. The modID (like "m6A_site_185574") of the column named "RNAmodList" in the mod-mirTar file corresponds the "modID" (4th column) of the RNA Modifications file. 
miR_files <- list.files('/data/rluo4/EpiTrans/DataCollection', pattern = 'human.miR')
miR <- fread('/data/rluo4/EpiTrans/DataCollection/human.miRmRNA.related.allRNAmod.list.txt', sep = '\t')
colN <- c("mirTarRand", "miRNAID", "mirRNAName", "mirTarChr", "mirTarStart",
          "mirTarEnd", "mirTarStrand", "predictIdNum", "PITAID",
          "PicTarID", "RNA22ID", "TargetScanID", "miRandaID",
          "miRmapID", "microTID", "RNAmodSiteNum", "RNAmodList", "RNAmodLoc")
# Read and process each count file
miR <- lapply(miR_files, function(x) {
  # Read data using fread for faster reading
  data <- fread(file.path('/data/rluo4/EpiTrans/DataCollection', x))
  data <- data[,1:18]
  # Extract required columns
  # Set column names same as file names
  colnames(data) <- colN# Display the first few rows with column names
  return(data)
})
df_in <- Reduce(function(x, y) merge(x, y, by = colN, all = TRUE), miR)
# write.table(df_in, '/data/rluo4/EpiTrans/DataCollection/human.miR.allRNAmod.list.txt', sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
# # Define the output file
# output_file="merged_hg38_miR.bed"
# # Define the header line with column names
# header="bindID\tmiRNAseq\talign\ttargetSeq"
# # Create the output file and add the header
# echo -e $header > $output_file
# # Loop through each BED file and append its content to the output file
# for file in hg38_miR*align_seq.txt; do
# echo "Processing $file..."
# # Append the content of the file, skipping the header if exists
# tail -n +2 $file >> $output_file
# done
# echo "All files have been merged into $output_file"
# Save miR_site as a .bed file
miR_merge.bed <- fread(file.path(miR_dir, '../merged_hg38_miR.bed'))
miR_merge.bed <- as.data.frame(miR_merge.bed)
colnames(miR_merge.bed)[1] <- 'mirTarRand'
miR_site.bed <- df_in[, c("mirTarChr", "mirTarStart", "mirTarEnd", "mirTarStrand", 'mirTarRand', "miRNAID", "mirRNAName")]
write.table(miR_site.bed, "/data/rluo4/EpiTrans/DataCollection/miR_site.bed",
            sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
# cp -rp /data/rluo4/EpiTrans/DataCollection/RBP /data/rluo4/EpiTrans/DataCollection/miRNA
# nohup bash ../miR-bed.sh >>./miR-bed.sh.o &
# find . -type f -name "*.bed" -size 0 -exec rm -f {} \;
miR_dir <- '/data/rluo4/EpiTrans/DataCollection/miRNA'
dsets <- list.files(miR_dir)
j <- NULL; k <- NULL
for (i in dsets) {
  print(i)
  file = list.files(file.path(miR_dir, i), pattern = 'miR.bed')
  j <<- c(j, file)
  # if(length(file) !=0 ){
  k <<- c(k, rep(i, length(file)))
  # }
}
bed_files <- data.frame(GSE=k, bed=j)
unique(bed_files$GSE)
colN <- c(colnames(DEM_summary.bed)[1:4], "mirTarChr", "mirTarStart", "mirTarEnd", "mirTarStrand", 'mirTarRand', "miRNAID", "mirRNAName")
DEM_merge.bed <- apply(bed_files, 1, function(x) {
  # Read data using fread for faster reading
  print(x)
  data <- fread(file.path(miR_dir, x['GSE'], x['bed']))
  print(head(data))
  print(dim(data))
  colnames(data) <- colN
  return(data)
})
names(DEM_merge.bed) <- paste0(bed_files$GSE, '_', gsub('_miR.bed', '', bed_files$bed))
plot_dir <- '/data/rluo4/RPMfunc/Output/GEO_top/'
# for (i in names(DEM_merge.bed)) {
changes <- c('GSE94613_MOLM13_METTL3KD')
for (i in changes) {
  GSE <- str_split(i, '_', simplify = T)[,1]
  inter.bed <- DEM_merge.bed[[i]]
  GSE.bed <- DEM_summary.bed[DEM_summary.bed$Dataset == GSE, ]
  GSE.bed$Peak_id <- gsub(' ', '-', GSE.bed$Peak_id)
  print(length(GSE.bed$Peak_id))
  print(dim(GSE.bed))
  condition = str_split(rownames(GSE.bed), '_vs', simplify = T)[,1]
  index = condition %in% i
  table(index)
  GSE.bed <- GSE.bed[index,]
  bed <- left_join( inter.bed, GSE.bed[, 4:9], by = 'Peak_id' )
  GSM <- unique( str_split(rownames(GSE.bed), '_vs', simplify = T)[,1])
  GSM <- remove_first_element(GSM)
  bed.path <- file.path(plot_dir, GSE, GSM,  paste0( 'show_top_miR.txt'))
  # bed <- DEM_merge.bed[[i]]
  colnames(bed)
  print(dim(bed))
  # print(unique(bed$mirRNAName))
  bed$miRTarSite <- paste0(bed$mirTarChr,":", bed$mirTarStart, "-", bed$mirTarEnd, ":", bed$mirTarStrand)
  bed$RNAmodSite <- paste0(bed$Chr,":", bed$Start, "-", bed$End)
  bed <- left_join( bed, miR_merge.bed, by = 'mirTarRand' )
  if(nrow(bed) !=0){
    write.table(bed[, c('miRTarSite', 'miRNAID', 'mirTarRand', 'mirRNAName',"GeneSymbol","Modification_type","Dataset","Enzyme","RNAmodSite")], bed.path, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
  }
}
# miR_DEM.bed <- Reduce(function(x, y) merge(x, y, by = colN, all = TRUE), DEM_merge.bed)
# i = 'GSE144984_NOMO1_ALKBH5KD'
i = 'GSE94613_MOLM13_METTL3KD'
DEG_subset <- DEG_summary[str_split(rownames(DEG_summary), '_vs', simplify = T)[,1] == i,]
DEM_subset <- DEM_summary[str_split(rownames(DEM_summary), '_vs', simplify = T)[,1] == i,]
inter_genes <- intersect(DEG_subset$GeneSymbol, DEM_subset$GeneSymbol)
intersect(bed$GeneSymbol, DEM_subset$GeneSymbol) # "C22orf24" [10] "TINF2"    "FZD9"
# "CPEB4"     "TMEM194A" 
# [4] "LINC00324" "RPL37"     "TBK1"  
intersect(inter_genes, bed$GeneSymbol)    
# [1] "LINC00324" "IL7R"      "SEMA6B"    "RPL37" 
targets <- c("CPEB4", "TMEM194A", "LINC00324", 'TBK1' )
# View(bed[bed$GeneSymbol %in% targets,])
########################################################################################
# 9.3 # A tab separated table of various of RNA modification sites related to SNV/SNP  #
########################################################################################
# 13: shiftPos -> The offset of SNV loci relative to the modification site. The range is from -10 to 10, the positive number means the SNV site appears in the downstream of RNA modification site and the negative means upstream.
SNV <- fread('/data/rluo4/EpiTrans/DataCollection/hg38.modSNV.tar.gz', sep = '\t')
SNV$V1[1] <- 'SNV_site_1'
SNV_site.bed <- SNV[, c(2:4,7,1,5:6,8:10)]
# write.table(SNV_site.bed, "/data/rluo4/EpiTrans/DataCollection/SNV_site.bed",
#             sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
SNV_dir <- '/data/rluo4/EpiTrans/DataCollection/SNV'
dsets <- list.files(SNV_dir)
j <- NULL; k <- NULL
for (i in dsets) {
  print(i)
  file = list.files(file.path(SNV_dir, i), pattern = 'SNV.bed')
  j <<- c(j, file)
  # if(length(file) !=0 ){
  k <<- c(k, rep(i, length(file)))
  # }
}
bed_files <- data.frame(GSE=k, bed=j)
unique(bed_files$GSE)
colN <- c(colnames(DEM_summary.bed)[1:4], "SNVChr", "SNVStart", "SNVEnd", "SNVStrand", "SNVsiteID", "SNVID", "SNVScore", "geneName", "Cancer", "Pubmed")
DEM_merge.bed <- apply(bed_files, 1, function(x) {
  # Read data using fread for faster reading
  print(x)
  data <- fread(file.path(SNV_dir, x['GSE'], x['bed']))
  print(head(data))
  print(dim(data))
  colnames(data) <- colN
  return(data)
})
names(DEM_merge.bed) <- paste0(bed_files$GSE, '_', gsub('_SNV.bed', '', bed_files$bed))
plot_dir <- '/data/rluo4/RPMfunc/Output/GEO_top/'
# for (i in names(DEM_merge.bed)) {
changes <- c('GSE94613_MOLM13_METTL3KD', 'GSE144984_NOMO1_ALKBH5KD', 'GSE144984_NOMO1_ALKBH5OE')
for (i in changes) {
  GSE <- str_split(i, '_', simplify = T)[,1]
  inter.bed <- DEM_merge.bed[[i]]
  GSE.bed <- DEM_summary.bed[DEM_summary.bed$Dataset == GSE, ]
  GSE.bed$Peak_id <- gsub(' ', '-', GSE.bed$Peak_id)
  print(length(GSE.bed$Peak_id))
  print(dim(GSE.bed))
  condition = str_split(rownames(GSE.bed), '_vs', simplify = T)[,1]
  index = condition %in% i
  table(index)
  GSE.bed <- GSE.bed[index,]
  bed <- left_join( inter.bed, GSE.bed[, 4:9], by = 'Peak_id' )
  GSM <- unique( str_split(rownames(GSE.bed), '_vs', simplify = T)[,1])
  GSM <- remove_first_element(GSM)
  bed.path <- file.path(plot_dir, GSE, GSM,  paste0( 'show_top_SNV.txt'))
  # bed <- DEM_merge.bed[[i]]
  colnames(bed)
  print(dim(bed))
  # print(unique(bed$SNVRNAName))
  bed$SNVTarSite <- paste0(bed$SNVChr,":", bed$SNVStart, "-", bed$SNVEnd, ":", bed$SNVStrand)
  bed$RNAmodSite <- paste0(bed$Chr,":", bed$Start, "-", bed$End)
  # bed <- left_join( bed, SNV_merge.bed, by = 'SNVTarRand' )
  if(nrow(bed) !=0){
    write.table(bed[, c('SNVTarSite', "SNVsiteID", "SNVID", "SNVScore", "geneName", "Cancer", "Pubmed","GeneSymbol","Modification_type","Dataset","Enzyme","RNAmodSite")], bed.path, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
  }
}
# i = 'GSE144984_NOMO1_ALKBH5KD'
# i = 'GSE94613_MOLM13_METTL3KD'
DEG_subset <- DEG_summary[str_split(rownames(DEG_summary), '_vs', simplify = T)[,1] == i,]
DEM_subset <- DEM_summary[str_split(rownames(DEM_summary), '_vs', simplify = T)[,1] == i,]
inter_genes <- intersect(DEG_subset$GeneSymbol, DEM_subset$GeneSymbol)
intersect(bed$GeneSymbol, DEM_subset$GeneSymbol)  
intersect(inter_genes, bed$GeneSymbol)    
# [1] "LINC00324" "IL7R"      "SEMA6B"    "RPL37" 
targets <- intersect(bed$GeneSymbol, DEM_subset$GeneSymbol[DEM_subset$Regulation=='DOWN'])#c("CPEB4", "TMEM194A", "LINC00324", 'TBK1' )
# View(bed[bed$GeneSymbol %in% targets,])
unique(sc_meta_SNV$RMP)
table(sc_meta_SNV$RMP %in% DEM_summary$Enzyme)

test <- sc_meta_SNV %>% filter(RMP=='ALKBH5')
maf_sub <- maf %>% filter(maf$Tumor_Sample_Barcode )
# test <- (maf[maf$Tumor_Sample_Barcode=='GSM5613768', c(1:15, 130:140)])

intersect(test$Hugo_Symbol, DEM_subset$GeneSymbol)

rs <-  paste0(RM.disease$seqnames,":", RM.disease$SNP_ChromStart, "-", RM.disease$SNP_ChromEnd, ":", RM.disease$SNP_Strand, ":", RM.disease$ref, "-",RM.disease$alt,":",RM.disease$type,":", RM.disease$Gene,":", RM.disease$refseq, ":", RM.disease$altseq)#, RM.disease$MD_Source )
RM.disease$identifier <- rs
# length(unique(rs))
# SNP_site.bed <- na.omit(RM.disease[, 1:4])#[,c(1:8, 26, 30)]
# SNP_site.bed$SNP_ChromEnd <- SNP_site.bed$SNP_ChromEnd +1
# Remove rows where SNP_ChromStart has NA values
SNP_site.bed <- RM.disease[,c(1:8, 22:23, 26, 30)]#'type', 'Gene', 'refseq', 'altseq')]
SNP_site.bed <- SNP_site.bed[complete.cases(SNP_site.bed$SNP_ChromStart), ]
# write.table(SNP_site.bed, "/data/rluo4/EpiTrans/DataCollection/SNP_site.bed",
#             sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
SNP_dir <- '/data/rluo4/EpiTrans/DataCollection/SNP'
dsets <- list.files(SNP_dir)
j <- NULL; k <- NULL
for (i in dsets) {
  print(i)
  file = list.files(file.path(SNP_dir, i), pattern = 'SNP.bed')
  j <<- c(j, file)
  # if(length(file) !=0 ){
  k <<- c(k, rep(i, length(file)))
  # }
}
bed_files <- data.frame(GSE=k, bed=j)
unique(bed_files$GSE)
colN <- c(colnames(DEM_summary.bed)[1:4], "SNPChr", "SNPStart", "SNPEnd", "SNPStrand", colnames(SNP_site.bed)[5:12])
DEM_merge.bed <- apply(bed_files, 1, function(x) {
  # Read data using fread for faster reading
  print(x)
  data <- fread(file.path(SNP_dir, x['GSE'], x['bed']))
  print(head(data))
  print(dim(data))
  colnames(data) <- colN
  return(data)
})
names(DEM_merge.bed) <- paste0(bed_files$GSE, '_', gsub('_SNP.bed', '', bed_files$bed))
for (i in names(DEM_merge.bed)) {
  GSE <- str_split(i, '_', simplify = T)[,1]
  GSM <- remove_first_element(i)
  bed.path <- file.path(plot_dir, GSE, GSM,  paste0( 'show_top_SNP.txt'))
  
  inter.bed <- DEM_merge.bed[[i]]
  inter.bed$identifier <-  paste0(inter.bed$SNPChr,":", inter.bed$SNPStart, "-", inter.bed$SNPEnd, ":", inter.bed$SNPStrand, ":", inter.bed$ref, "-",inter.bed$alt,":",inter.bed$type,":", inter.bed$Gene,":", inter.bed$refseq, ":", inter.bed$altseq)
  inter.bed <- left_join( inter.bed, RM.disease[, -c(1:8, 22:23, 26, 30)], by = 'identifier' )
  
  GSE.bed <- DEM_summary.bed[DEM_summary.bed$Dataset == GSE, ]
  GSE.bed$Peak_id <- gsub(' ', '-', GSE.bed$Peak_id)
  print(length(GSE.bed$Peak_id))
  print(dim(GSE.bed))
  condition = str_split(rownames(GSE.bed), '_vs', simplify = T)[,1]
  index = condition %in% i
  table(index)
  GSE.bed <- GSE.bed[index,]
  inter.bed <- inter.bed[inter.bed$Gene %in% GSE.bed$GeneSymbol,]
  
  bed <- left_join( inter.bed, GSE.bed[, 4:9], by = 'Peak_id' )
  colnames(bed)
  print(dim(bed))
  # print(unique(bed$SNPRNAName))
  bed$SNPTarSite <- paste0(bed$SNPChr,":", bed$SNPStart, "-", bed$SNPEnd, ":", bed$SNPStrand)
  bed$RNAmodSite <- paste0(bed$Chr,":", bed$Start, "-", bed$End)
  bed <- as.data.frame(bed)
  if(nrow(bed) !=0){
    write.table(bed[, c('SNPTarSite',colnames(inter.bed)[9:38], "GeneSymbol","Modification_type","Dataset","Enzyme","RNAmodSite")], bed.path, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
  }
}
##########################################
# changes <- unique(str_split(rownames(DEG_summary), '_vs', simplify = T)[,1])
# for (i in changes) {
#   GSE <- str_split(i, '_', simplify = T)[,1]
#   GSM <- remove_first_element(i)
#   bed.path <- file.path(plot_dir, GSE, GSM,  paste0( 'show_top_SNV.txt'))
#   DEM.path <-  file.path(plot_dir, GSE, GSM,  paste0( 'show_top_DEM.txt'))
#   GSE.bed <- read.table(DEM.path, sep = '\t', header = T)
#   print(dim(GSE.bed))
#   RM <- tolower(unique(GSE.bed$Modification_type))
#   inter.bed <- RM.disease[grepl(RM, RM.disease$type), ]
#   inter.bed <- inter.bed[inter.bed$Gene %in% GSE.bed$GeneSymbol,]
#   colnames(inter.bed)[30] <- 'GeneSymbol'
#   inter.bed$SNPTarSite <- paste0(inter.bed$seqnames,":", inter.bed$SNP_ChromStart, "-", inter.bed$SNP_ChromEnd, ":", inter.bed$SNP_Strand)
#   bed <- left_join( inter.bed, GSE.bed, by = 'GeneSymbol' )
#   # colnames(bed)
#   print(dim(bed))
#   # print(unique(bed$SNVRNAName))
#   bed$RNAmodSite <- paste0(bed$Chr,":", bed$Start, "-", bed$End)
#   # bed <- left_join( bed, SNV_merge.bed, by = 'SNVTarRand' )
#   if(nrow(bed) !=0){
#     write.table(bed[, c('SNPTarSite',colnames(inter.bed)[5:33], "Modification_type","Dataset","Enzyme","RNAmodSite")], bed.path, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
#   }
# }
# library(GenomicRanges)
# # Convert SNP region (columns 1-3) to a GRanges object
# snp_granges <- GRanges(
#   seqnames = bed$seqnames,
#   ranges = IRanges(start = bed$SNP_ChromStart, end = bed$SNP_ChromEnd)
# )
# # Convert RNA modification site (columns 34-36, i.e., Chr, Start, End) to a GRanges object
# rna_granges <- GRanges(
#   seqnames = bed$Chr,
#   ranges = IRanges(start = bed$Start, end = bed$End)
# )
# # Find overlaps between SNP site and RNA modification site
# overlaps <- findOverlaps(snp_granges, rna_granges)
# # Get the rows with overlaps
# overlap_rows <- unique(subjectHits(overlaps))
# # Filter the original data to retain only the rows where overlaps exist
# filtered_bed <- bed[overlap_rows, ]
# # # Merge the two data frames based on overlapping regions
# library(dplyr)
# # Define a function to check overlap
# is_overlap <- function(start1, end1, start2, end2) {
#   return(start1 <= end2 & end1 >= start2)
# }
# # Perform the join
# merged_df <- GSE.bed %>%
#   rowwise() %>%
#   mutate(
#     overlaps = list(inter.bed %>%
#                       filter(
#                         Chr == Chr,
#                         is_overlap(Start, End, Start, End)
#                       ) %>%
#                       pull(SNPTarSite))
#   ) %>%
#   unnest(overlaps)
# # View the result
# print(merged_df)
##########################################
# 9.4 #   RNA subcellular localization   #
##########################################
RNAloc <- read.table('/data/rluo4/EpiTrans/DataCollection/All RNA subcellular localization information.txt', fill = T, header = T, sep = '\t')
table(RNAloc$Species)
RNAloc <- RNAloc[RNAloc$Species %in% 'Homo sapiens', ]
multiomics_summary <- read.table('/data/rluo4/RPMfunc/Output/summary/multiomics_summary.txt', fill = T)

target_summary <- bulk_result[, c(1,4,5,6)]
target_summary <- target_summary[!duplicated(paste(target_summary$GeneSymbol, target_summary$Condition, target_summary$Dataset)),]
target_summary$Folder <- gsub('_vs_WT', '',target_summary$Condition)
index = paste( target_summary$Dataset, target_summary$Condition); length(unique(index))
multiomics_index = paste( multiomics_summary$project_id, multiomics_summary$patient_id )
length(unique(multiomics_index))
setdiff(unique(index), unique(multiomics_index))
table(index %in% multiomics_index)
target_summary$Tissue <- multiomics_summary$tissue[match(index, multiomics_index)]
id_mapped.hg19 <- read.table("/data/rluo4/hg38/gencode.v44lift37.annotation.gene.probeMap",header = T)
id_mapped.hg38 <- read.table("/data/rluo4/hg38/gencode.v36.annotation.gene.probeMap",header = T)
table(unique(target_summary$GeneSymbol) %in% unique(gene_info$geneSymbol))
table(unique(target_summary$GeneSymbol) %in% unique(id_mapped.hg38$gene))
table(unique(target_summary$GeneSymbol) %in% unique(id_mapped.hg19$gene))
table(unique(target_summary$GeneSymbol) %in% union(unique(id_mapped.hg38$gene), unique(gene_info$geneSymbol)))
all_version_genes <- union(union(unique(id_mapped.hg38$gene), unique(gene_info$geneSymbol)),
                           unique(id_mapped.hg19$gene))
table(unique(target_summary$GeneSymbol) %in% all_version_genes)
target_summary <- target_summary[target_summary$GeneSymbol %in% all_version_genes, ]
target_summary$uniprot <- gene_info$UniProtAcc[match(target_summary$GeneSymbol,gene_info$geneSymbol)]

uniprot <- read.csv('/data/rluo4/All/uniprot-hs.tsv',fill = T,header = T,sep = '\t')
library(org.Hs.eg.db)
table(target_summary$GeneSymbol %in% uniprot$Gene.Names)
table(target_summary$uniprot %in% uniprot$Entry)
uniprot$GeneSymbol <- str_split(uniprot$Gene.Names, ' ', simplify = T)[,1]
# table(target_summary$GeneSymbol %in% uniprot$GeneSymbol)
# target_loc <- left_join(target_summary, uniprot, by = 'GeneSymbol')
target_loc <- uniprot
target_loc$uniprotLoc <- target_loc$Subcellular.location..CC. #uniprot$Subcellular.location..CC.[match(gene_loc$geneSymbol,uniprot$Gene.Names)]
table(target_loc$uniprotLoc=='')#133382
target_loc$uniprotLoc <- gsub('SUBCELLULAR LOCATION: ','',target_loc$uniprotLoc)
target_loc$Nucleus <- ifelse(grepl('ucleus',target_loc$uniprotLoc),1,0)
target_loc$Secreted <- ifelse(grepl('ecreted',target_loc$uniprotLoc),1,0)
target_loc$Cytoplasm <- ifelse(grepl('ytoplasm',target_loc$uniprotLoc),1,0)
target_loc$Surface <- ifelse(grepl('ell membrane',target_loc$uniprotLoc),1,0)
multi.loc <- apply(target_loc[, c('Nucleus','Secreted','Cytoplasm','Surface')],1,function(target_loc){
  s <- sum(target_loc)
  return(s)
})
dim(target_loc);length(multi.loc)
target_loc$multi.loc <- multi.loc
table(target_loc$multi.loc)
table(target_loc$multi.loc %in% c(0,1))
# target_loc <- target_loc[target_loc$multi.loc %in% c(0,1),]
target_loc$uniprotLoc[is.na(target_loc$uniprotLoc)] = ''
table(target_loc$uniprotLoc=='')#3638
length(strsplit(target_loc$uniprotLoc[100], "[;]")[[1]])
target_loc$Nucleus <- gsub(1, 'Nucleus', target_loc$Nucleus)
target_loc$Secreted <- gsub(1, 'Secreted', target_loc$Secreted)
target_loc$Cytoplasm <- gsub(1, 'Cytoplasm', target_loc$Cytoplasm)
target_loc$Surface <- gsub(1, 'Surface', target_loc$Surface)
# if(target_loc$multi.loc>=2){
target_loc$proteinLoc <- paste(target_loc$Nucleus, target_loc$Secreted, target_loc$Cytoplasm, target_loc$Surface, sep = ';')
target_loc$proteinLoc <- gsub('0', '', gsub(';0', '', gsub('0;', '', target_loc$proteinLoc)))
# }
table(target_loc$proteinLoc)
# table(target_loc$proteinLoc[target_loc$proteinLoc %ni% c('Nucleus','Secreted','Cytoplasm','Surface')])
# target_loc$proteinLoc[target_loc$proteinLoc ==''] <- 'Others'
# table(target_loc$proteinLoc)
# target_summary <- target_loc[, -c(2,8:12)]
# colnames(target_summary)
# colnames(RNAloc) <- gsub('RNA_Symbol', 'GeneSymbol', colnames(RNAloc))
# target_summary <- left_join(target_summary, RNAloc[,-(1:2)], by = 'GeneSymbol')
write.table(target_summary, "/data/rluo4/RPMfunc/Output/summary/target_summary.txt", row.names = FALSE, sep = '\t', quote = F)
write.table(target_loc, "/data/rluo4/EpiTrans/DataCollection/uniprotLoc.txt", row.names = FALSE, sep = '\t', quote = F)

# colnames(RNAloc) <- gsub('GeneSymbol', 'geneSymbol', colnames(RNAloc))
# gene_summary <- left_join(gene_loc[, -c(1, 13:17)], RNAloc[,-(1:2)], by = 'geneSymbol')
write.table(gene_summary,'/data/rluo4/RPMfunc/Output/summary/gene_summary.txt',sep='\t',quote = F,row.names = F)

####################################################################################################
# 10) A tab separated table of chemicals and drugs interact with RMP-regulated genes from DGIdb 4.0 #
####################################################################################################
drug_interaction <- fread('/data/rluo4/EpiTrans/DataCollection/DGIdb_drug.tsv', fill = T)
dim(drug_interaction)
drug_info <-  fread('/data/rluo4/EpiTrans/DataCollection/interactions.tsv', fill = T)
dim(drug_info)
length(unique(drug_interaction$gene_claim_name))
table(drug_interaction$gene_name %in% drug_info$gene_name)
colnames(drug_interaction)
a = as.data.frame(drug_interaction[,-c(1,2,5,8,9)])
b = as.data.frame(drug_info[, c('gene_name', 'PMIDs')])
drug_interaction$PMID <-b$PMIDs[match(drug_interaction$gene_name, b$gene_name)] #left_join(a, b, by = 'gene_name')
# write.table(drug_interaction[,-c(1,2,5,8,9)], "/data/rluo4/RPMfunc/Output/summary/drug.DGIdb_summary.txt", row.names = FALSE, sep = '\t', quote = F)
table(unique(DEG_summary$GeneSymbol) %in% unique(drug_interaction$gene_name))
table(unique(DEM_summary$GeneSymbol) %in% unique(drug_interaction$gene_name))
table(inter_genes %in% unique(drug_interaction$gene_name)) #4302
# Patt <- read.table('/data/rluo4/All/Output/organ13-deg.txt',sep = '\t',fill = TRUE,header = TRUE) #PCTanno
# table(unique(Patt$Symbol) %in% unique(drug_interaction$gene_name))
# table(inter_genes %in% unique(Patt$Symbol)) #4302
# drug_info <- fread('/data/rluo4/EpiTrans/DataCollection/drugs.tsv')
# gene_drug <- fread('/data/rluo4/EpiTrans/DataCollection/genes.tsv')
# catogory <- fread('/data/rluo4/EpiTrans/DataCollection/categories.tsv')
# drug_summary <- fread('/data/rluo4/summary/drug_summary.txt', fill = T)
