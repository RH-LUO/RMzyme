#!/usr/local/bin/Rscript
# title: "Single-cell Landscape of Malignant Transition: Unraveling Cancer Cell-of-Origin and Heterogeneous Tissue Microenvironment"
# title:
# author: "Ruihan Luo"
# date: "May 29th,2024"
# rm(list=ls())
library(maftools)
library(data.table)
library(stringr)
library(dplyr)
options(bitmapType = 'cairo')
################################################################################
# 1) load pdata of PCTanno database
################################################################################
# load('~/lrh/All/Output/organ13_Epi.RData')# from database-sum.R in ts860
# data_path = '/home/lorihan/lrh/All/Output'
data_path = '/data/rluo4/All/Output'
setwd(data_path)
cell_summary <- read.csv('cell_summary.txt',sep = "\t",header = T)
tissue_summary <- read.table('/data/rluo4/RPMfunc/PCTanno_pdata/tissue_summary.txt')
extract_last_element <- function(x) {
  split_string <- strsplit(x, "_")
  last_element <- sapply(split_string, function(y) tail(y, n = 1))
  return(last_element)
}
################################################################################
# 2) arrange the SNV pdata
################################################################################
# Load Diseased Samples
Diseased_Samples <- read.table('/data/rluo4/RPMfunc/PCTanno_pdata/Diseased_Samples.txt')
# Load Normal Samples
normal_samples <- unique(cell_summary$orig.ident[cell_summary$Tissue %in% c('Healthy','ADJ')])
remove_samples <- setdiff(Diseased_Samples$orig.ident, tissue_summary$Replicates)# 4
normal_samples <- c(normal_samples, remove_samples) #417
# Tissue <- cell_summary[cell_summary$Organ==tis,]
# unique(Tissue$sample_name)
################################################################################
# index = ! cell_summary$Tissue %in% c('Healthy','ADJ')
# cell_summary <- cell_summary[index, ]
################################################################################
# # Method 1 (test on 57)
################################################################################
################################################################################
# # Method 2 (test on 36)
################################################################################
# # save(cell_summary, file = file.path(data_path, 'cell_summary_diseased.rds'))
################################################################################
# 3) arrange the SNV data
################################################################################
# CHROM  Start   End     REF     ALT     FILTER  Cell_types      Up_context      Down_context    N_ALT   
# Dp      Nc   -Bc      Cc      VAF     CCF     BCp     CCp     Cell_types_min_BC       Cell_types_min_CC       Rest_BC Rest_CC Fisher_p      Cell_type_Filter        INFO  
load(file.path(data_path, 'cell_summary_diseased.rds'))
# Helper function to load and filter files
load_and_filter_files <- function(dir, normal_samples, pattern) {
  solution <- list.files(dir, pattern = pattern)
  nor <- paste0(normal_samples, '.', pattern)  # Adjusted the separator
  solution <- solution[!solution %in% nor]
  file_paths <- file.path(dir, solution)
  # Use fread correctly within lapply
  solutions <- lapply(file_paths, function(fp) {
    fread(fp, sep = '\t', stringsAsFactors = FALSE)
  })
  names(solutions) <- sapply(solution, function(x) strsplit(x, '\\.variants.avinput')[[1]][1])
  return(solutions)
}
# Helper function to load and filter Annovar files
load_and_filter_annovarToMaf <- function(dir, normal_samples, pattern) {
  solution <- list.files(dir, pattern = pattern)
  nor <- paste(normal_samples, pattern, sep = '.')
  solution <- solution[!solution %in% nor]
  file_paths <- file.path(dir, solution)
  # solutions <- lapply(file_paths, annovarToMaf)
  solutions <- lapply(file_paths, function(fp) {
    annovarToMaf(fp, refBuild = "hg38")
  })
  names(solutions) <- sapply(solution, function(x) strsplit(x, '\\.hg38_multianno.txt')[[1]][1])
  return(solutions)
}
# Function to add VAF and CCF
add_vaf_ccf <- function(maf, solution, nameindex) {
  y <- as.data.frame(maf[[nameindex]])
  z <- gsub('.hg38_multianno.txt', '.variants.avinput', nameindex)
  w <- as.data.frame(solution[[z]]) 
  # Ensure w$V6 is split correctly and the lengths match
  split_v6 <- str_split(w$V6, '[-]', simplify = TRUE)  
  if (nrow(split_v6) == nrow(y)) {
    VAF <- split_v6[, 15]
    CCF <- split_v6[, 16]
    
    y$VAF <- VAF
    y$CCF <- CCF
    y$Tumor_Sample_Barcode <- nameindex
    return(y)
  } 
  # else {
  #   stop("Length of split V6 does not match the number of rows in y")
  # }
}

# Function to apply filters
apply_filters <- function(y) {
  filter_maf <- unique(c(
    which(y$X1000g2015aug_all > 0.05),
    which(y$ExAC_ALL > 0.05),
    which(y$gnomAD_genome_ALL > 0.05)
  ))
  y <- y[-filter_maf,]
  y <- y[y$VAF != '.' & y$VAF > 0.05, ]
  filter_repeat <- which(!is.na(y$rmsk) & !is.na(y$genomicSuperDups))
  y <- y[-filter_repeat, ]
  y <- y[!is.na(y$FATHMM_score) & y$FATHMM_score > 0.7, ]
  return(y)
}

# Process scMapping data
process_scMapping <- function(dir, annovar_mafs, cohort) {
  scMapping <- list.files(dir)
  scMapping <- scMapping[!scMapping %in% normal_samples]
  anno <- vector("list", length(scMapping))
  names(anno) <- scMapping
  
  for (i in scMapping) {
    path <- file.path(dir, i)
    solution <- list.files(path, pattern = 'single_cell_genotype.tsv')
    file_paths <- file.path(path, solution)
    solutions <- lapply(file_paths, fread, stringsAsFactors = FALSE)
    names(solutions) <- sapply(solution, function(x) strsplit(x, '\\.single_cell_genotype.tsv')[[1]][1])
    
    solutions <- lapply(names(solutions), function(nameindex) {
      b <- as.data.frame(annovar_mafs[[i]])
      c <- as.data.frame(solutions[[nameindex]])
      colnames(c)[1:2] <- colnames(b)[2:3]
      c$Start_Position <- as.character(c$Start_Position)
      y <- left_join(b, c[,-1], by = 'Start_Position')
      y <- y[!is.na(y$Cell_type_observed),]
      y <- y[y$Tumor_Sample_Barcode %in% Diseased_Samples$Replicates,]
      index <- match(y$Tumor_Sample_Barcode, Diseased_Samples[Diseased_Samples$Cohort==cohort,]$Replicates)
      y$Tumor_Sample_Barcode <- Diseased_Samples[Diseased_Samples$Cohort==cohort,]$orig.ident[index]
      
      Epi <- cell_summary[cell_summary$sample_name == cohort,]
      # Epi$CB <- extract_clean_barcode(Epi$barcode)
      Epi$CB <- paste(Epi$CB, Epi$orig.ident, sep = '--')
      y$CB <- paste(y$CB, y$Tumor_Sample_Barcode, sep = '--')
      Epi$CB <- paste(Epi$CB, Epi$Cell_Type, sep = '--')
      y$CB <- paste(y$CB, y$Cell_type_observed, sep = '--')
      
      y <- y[y$CB %in% Epi$CB, ]
      y <- apply_filters(y)
      return(y)
    })
    
    anno[[i]] <- data.table::rbindlist(solutions, fill = TRUE)
  }
  
  scRNA_annovar <- NULL
  scRNA_annovar <- do.call(rbind, anno)
  
  return(scRNA_annovar)
}

# Process scMapping data for Cervix
datasets <- unique(Diseased_Samples$Cohort)
# datasets = 'Losic'
for(cohort in datasets){
  # Define paths
  base_dir <- '/data/rluo4/RPMfunc/SCOMATIC'
  start_time <- Sys.time()
  tis <- Diseased_Samples$Tissue[match(cohort, Diseased_Samples$Cohort)]
  print(paste("start parse:", tis, '--', cohort, ' at', start_time))
  
  annovar_dir <- file.path(base_dir, tis, cohort, "anno.var")
  scMapping_dir <- file.path(base_dir, tis, cohort, "scMapping")
  out_dir <-  '/data/rluo4/RPMfunc/Output/scRNA/PCTanno_SNV'
  gse_acc <- Diseased_Samples$GSE[match(cohort, Diseased_Samples$Cohort)]
  SNV_data <- file.path(out_dir, paste0(gse_acc, '_SCOMATIC_annovar.RData'))
  # if(! tis %in% c('Breast', 'Cervix', 'CRC', 'Endometrium', 'Esophagus', 'GC')){
  #   print(paste0(tis, " run on ts860 !"))
  #   next;
  # }
  
  if( !file.exists(SNV_data)){
    # Load and filter Guo files
    test_solutions <- load_and_filter_files(annovar_dir, normal_samples, 'variants.avinput')
    test_annovar_mafs <- load_and_filter_annovarToMaf(annovar_dir, normal_samples, 'hg38_multianno.txt')
    test_annovar_mafs <- lapply(names(test_annovar_mafs), function(nameindex) {
      add_vaf_ccf(test_annovar_mafs, test_solutions, nameindex)
    })
    names(test_annovar_mafs) <- names(test_solutions)
    
    test_scRNA_annovar <- process_scMapping(scMapping_dir, test_annovar_mafs, cohort)
    # Combine Guo and Hua1 data
    # Save final data
    save(test_scRNA_annovar, file = SNV_data)
    end_time <- Sys.time()
    print(end_time)
    time_taken <- end_time - start_time
    print(paste("Time taken to process", tis, '--', cohort, ':', time_taken, 'mins' ))
    
  }
}