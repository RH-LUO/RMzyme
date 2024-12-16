#!/usr/local/bin/Rscript
# title: "RMzyme: Regulations of RNA Modifying Enzymes in Human"
# title:
# author: "Ruihan Luo"
# date: "July 15th,2024"
# rm(list=ls())
library(maftools)
library(data.table)
library(stringr)
library(dplyr)
options(bitmapType = 'cairo')
# Required libraries
library(future.apply)
options(future.globals.maxSize = 2 * 1024^3)  # Adjust if absolutely necessary
################################################################################
# 1) load pdata of PCTanno & disco database
################################################################################
# load('~/lrh/All/Output/organ13_Epi.RData')# from database-sum.R in ts860
data_path = '/data2/rluo4/All/Output'
setwd(data_path)
RMP_update <- readxl::read_excel('/data2/rluo4/EpiTrans/DataCollection/RMP_update.xlsx')
load(file = '/data2/rluo4/RPMfunc/Output/scRNA/sc_allsample_1122.RData', verbose = T)
# sc_allmeta <- readRDS(file = '/data2/rluo4/RPMfunc/Output/scRNA/sc_allmeta_1122.rds')
# RMzyme_alltissues <- sort(unique(sc_allmeta$tissue));
# RMzyme_alldatasets <- sort(unique(sc_allmeta$project_id));
# extract_last_element <- function(x) {
#   split_string <- strsplit(x, "_")
#   last_element <- sapply(split_string, function(y) tail(y, n = 1))
#   return(last_element)
# }
# # cell_summary <- read.csv('cell_summary.txt',sep = "\t",header = T)
# # sc_majortype <- readRDS(file = '/data2/rluo4/RPMfunc/Output/scRNA/sc_majortype_1122.rds')
# # tissue_summary <- read.table('/data2/rluo4/RPMfunc/PCTanno_pdata/tissue_summary.txt')
# ################################################################################
# # 2) arrange the SNV pdata
# ################################################################################
################################################################################
#           summary of sc_majortype for RMzyme database                        #
################################################################################
length(unique(sc_summary$project_id)) #325 --> 328 --> 332
sample.index <- c(grep("ormal", sc_summary$sample_type),grep("ealthy", sc_summary$sample_type),
                  grep("ontrol", sc_summary$sample_type))
sc_Diseased_Samples <- sc_summary[ -sample.index, c(2,3,5:8)] # 219 --> 222 --> 224 disease datasets
table(sc_Diseased_Samples$disease)
sc_Diseased_Samples <- sc_Diseased_Samples[! grepl('ealthy', sc_Diseased_Samples$disease), ]
# # unique(sc_majortype[! grepl('adjacent|ealthy', sc_majortype$disease),]$project_id)
# RMzyme_pheno <- as.data.frame(table(sc_summary$sample_type, sc_summary$disease))
# RMzyme_pheno <- RMzyme_pheno[RMzyme_pheno$Freq!=0,]
# ################################################################################
# #      organize the diseased data of disco named as sc_cell_summary            #
# ################################################################################
sc_cell_summary <- readRDS("/data2/rluo4/RPMfunc/Output/scRNA/sc_cell_summary.rds") # from SComatic_disco_GSE185381.R
################################################################################
# 3) arrange the SNV data
################################################################################
# # CHROM  Start   End     REF     ALT     FILTER  Cell_types      Up_context      Down_context    N_ALT   
# # Dp      Nc   -Bc      Cc      VAF     CCF     BCp     CCp     Cell_types_min_BC       Cell_types_min_CC       Rest_BC Rest_CC Fisher_p      Cell_type_Filter        INFO  
# Process scMapping data for disco
# datasets <-  readRDS('/data2/rluo4/RPMfunc/disco_pdata/sc_SNV_datasets.rds') #anno.var #unique(Diseased_Samples$Cohort)
# GSE192935 = disco_majortype[disco_majortype$project_id=='GSE192935', ]
# GSE222078 = disco_majortype[disco_majortype$project_id=='GSE222078', ]
# GSE163974 = disco_majortype[disco_majortype$project_id=='GSE163974', ]

datasets <- list.files('/data2/rluo4/EpiTrans/Disco/anno.var')
table(datasets %in% unique(sc_cell_summary$project_id))
setdiff(datasets, unique(sc_cell_summary$project_id))#"GSE121636" "GSE140042" "GSE156405"
# [1] "GSE158037" "GSE162498" "GSE166676" "GSE173468" "GSE184291"
# [6] "GSE201091" "GSE210066" "GSE211783" "GSE212038" "GSE212461"
# [11] "GSE213047" "GSE214207" "GSE217792" "GSE224273" "GSE227691"

# (base) rluo4@sbcsmadlp005:/data2/rluo4/EpiTrans/Disco/anno.var/GSE121636$ cd ..
# (base) rluo4@sbcsmadlp005:/data2/rluo4/EpiTrans/Disco/anno.var$ mv GSE121636 GSE121638
# (base) rluo4@sbcsmadlp005:/data2/rluo4/EpiTrans/Disco/anno.var/GSE158037$ cd ..
# (base) rluo4@sbcsmadlp005:/data2/rluo4/EpiTrans/Disco/anno.var$ mv GSE158037 GSE158038
# (base) rluo4@sbcsmadlp005:/data2/rluo4/EpiTrans/Disco/anno.var$ mv GSE162498 GSE162500
# datasets <- datasets[datasets %in% unique(sc_cell_summary$project_id)] #32 --> 33 --> 31
table(datasets %in% sc_cell_summary$ddt_id) #66
studies_running <- c('GSE185344', 'GSE162025', 'GSE221156', 'GSE211033', 'GSE207493', 
                     "GSE196638", "GSE203191", "GSE212447", "GSE185965")# , 'GSE185965'(2 samples undone!); 'GSE210543'(done in 11/30/2024)
asp_sc_file_done <- read.table(file="/data2/rluo4/EpiTrans/DataCollection/asp_sc_file_done.txt",sep = "\t")
table(datasets %in% unique(asp_sc_file_done$V1))
asp_sc_file_done_renew <- read.table(file="/data2/rluo4/EpiTrans/DataCollection/asp_sc_file_done_renew1.txt",sep = "\t")
setdiff(asp_sc_file_done_renew$V1, asp_sc_file_done$V1)
setdiff(asp_sc_file_done_renew$V2, asp_sc_file_done$V2)

table(datasets %in% unique(asp_sc_file_done_renew$V1))#58
table(datasets %in% studies_running)
datasets <- datasets[datasets %in% unique(asp_sc_file_done_renew$V1)] #58
table(datasets %in% studies_running)

library(data.table)
library(stringr)
library(dplyr)
library(future.apply)
library(parallel)

# Optimize future processing
n_cores <- min(detectCores() - 2, 4)
plan(multisession, workers = n_cores)
options(future.globals.maxSize = 8 * 1024^3)  # Increase memory limit

################################################################################
# Helper functions
################################################################################
# Helper function to load and filter files
load_and_filter_files <- function(dir,  pattern) {
  solution <- list.files(dir, pattern = pattern)
  file_paths <- file.path(dir, solution)
  # Use fread correctly within lapply
  solutions <- lapply(file_paths, function(fp) {
    fread(fp, sep = '\t', stringsAsFactors = FALSE)
  })
  names(solutions) <- sapply(solution, function(x) strsplit(x, '\\.variants.avinput')[[1]][1])
  return(solutions)
}
# Helper function to load and filter Annovar files
load_and_filter_annovarToMaf <- function(dir, pattern) {
  solution <- list.files(dir, pattern = pattern)
  file_paths <- file.path(dir, solution)
  # solutions <- lapply(file_paths, annovarToMaf)
  solutions <- lapply(file_paths, function(fp) {
    annovarToMaf(fp, refBuild = "hg38")
  })
  names(solutions) <- sapply(solution, function(x) strsplit(x, '\\.hg38_multianno.txt')[[1]][1])
  return(solutions)
}
# Add VAF and CCF columns to MAF data
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
# Apply variant filters
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
  # y <- y[!is.na(y$FATHMM_score) & y$FATHMM_score > 0.7, ]
  return(y)
}

################################################################################
# Process single-cell mapping data
################################################################################
# Helper function to handle the merging and filtering
process_cell_solution <- function(i, solutions, annovar_mafs, cohort) {
  b <- as.data.frame(annovar_mafs[[i]])
  
  # For each cell solution (cell type)
  results <- lapply(names(solutions), function(nameindex) {
    c <- as.data.frame(solutions[[nameindex]])
    colnames(c)[1:2] <- colnames(b)[2:3]
    c$Start_Position <- as.character(c$Start_Position)
    
    # Perform the left join and apply filtering
    y <- left_join(b, c[, -1], by = 'Start_Position')
    y <- y[!is.na(y$Cell_type_observed),]
    y <- y[y$Tumor_Sample_Barcode %in% sc_Diseased_Samples$patient_id,]
    
    index <- match(y$Tumor_Sample_Barcode, sc_Diseased_Samples[sc_Diseased_Samples$project_id == cohort,]$patient_id)
    y$Tumor_Sample_Barcode <- sc_Diseased_Samples[sc_Diseased_Samples$project_id == cohort,]$patient_id[index]
    
    # Combine barcode and cell type information
    Epi <- sc_cell_summary[sc_cell_summary$project_id == cohort, ]
    Epi$CB <- paste(Epi$CB, Epi$patient_id, sep = '--')
    y$CB <- paste(y$CB, y$Tumor_Sample_Barcode, sep = '--')
    Epi$CB <- paste(Epi$CB, Epi$Cell_Type, sep = '--')
    y$CB <- paste(y$CB, y$Cell_type_observed, sep = '--')
    
    # Filter based on barcode match
    y <- y[y$CB %in% Epi$CB, ]
    
    # Apply filters
    y <- apply_filters(y)
    
    return(y)
  })
  
  # Return the processed solutions for all cell types
  return(results)
}

# Main function to process scMapping data
process_scMapping <- function(dir, annovar_mafs, cohort) {
  scMapping_dirs <- list.files(dir, full.names = TRUE)
  scMapping_samples <- str_split(scMapping_dirs, '/', simplify = T)[,8]
  scMapping_dirs <- scMapping_dirs[scMapping_samples %in% names(annovar_mafs) ]
  anno <- vector("list", length(scMapping_dirs))
  names(anno) <- scMapping_dirs
  
  for (path in scMapping_dirs) {
    # Read single-cell genotype data
    solution_files <- list.files(path, pattern = 'single_cell_genotype.tsv', full.names = TRUE)
    solutions <- lapply(solution_files, fread, stringsAsFactors = FALSE)
    
    # Extract cell type names from solution file names
    cell_types <- sapply(solution_files, function(x) strsplit(basename(x), '\\.single_cell_genotype.tsv')[[1]][1])
    names(solutions) <- cell_types
    
    # Process each cell solution (cell type) for the current directory
    results <- process_cell_solution(basename(path), solutions, annovar_mafs, cohort)
    
    # Filter out empty or invalid results
    results <- results[!sapply(results, function(x) is.null(x) || nrow(x) == 0)]
    
    # Check column count, and filter cell types with fewer than 140 columns
    results <- lapply(results, function(x) {
      if (ncol(x) < 140) {
        print(paste0("Skipping cell type with fewer than 140 columns"))
        return(NULL)
      }
      return(x)
    })
    
    # Remove NULL elements from the list
    results <- results[!sapply(results, is.null)]
    
    # Combine all valid solutions for this directory into the `anno` list
    anno[[basename(path)]] <- data.table::rbindlist(results, fill = TRUE)
    rm(results)  # Free memory
    gc()
  }
  
  # Combine all processed data into a single data frame
  scRNA_annovar <- do.call(rbind, anno)
  rm(anno)
  gc()
  # Return the final combined data
  return(scRNA_annovar)
}


################################################################################
# Main Pipeline Execution
################################################################################
for(cohort in datasets) {
  # cohort = 'GSE185381'
  base_dir <- '/data2/rluo4/EpiTrans/Disco'
  start_time <- Sys.time()
  tis <- sc_Diseased_Samples$tissue[match(cohort, sc_Diseased_Samples$project_id)]
  print(paste("Starting parse:", tis, '--', cohort, 'at', start_time))
  
  annovar_dir <- file.path(base_dir, "anno.var", cohort)
  scMapping_dir <- file.path(base_dir, "scMapping", cohort)
  out_dir <- '/data2/rluo4/RPMfunc/Output/scRNA/disco_SNV'
  SNV_data <- file.path(out_dir, paste0(cohort, '_SCOMATIC_annovar.RData'))
  
  if (!file.exists(SNV_data)) {
    test_solutions <- load_and_filter_files(annovar_dir, 'variants.avinput')
    test_annovar_mafs <- load_and_filter_annovarToMaf(annovar_dir, 'hg38_multianno.txt')
    test_annovar_mafs <- future_lapply(names(test_annovar_mafs), function(nameindex) {
      add_vaf_ccf(test_annovar_mafs, test_solutions, nameindex)
    })
    names(test_annovar_mafs) <- names(test_solutions)
    test_annovar_mafs <- test_annovar_mafs[!sapply(test_annovar_mafs, is.null)]
    # GSE185381_allmeta <-  readRDS('/data2/rluo4/RPMfunc/Output/scRNA/GSE185381_allmeta.rds')#subset(sc_allmeta, ddt_id == 'GSE185381')
    # table(GSE185381_allmeta$split_ids)
    # GSE185381_allmeta <- subset(GSE185381_allmeta, split_ids == 'sc_adult_AML')
    # sc_adult_AML <- unique(GSE185381_allmeta$patient_id) #38
    # sc_adult_AML <- intersect(sc_adult_AML, sc_Diseased_Samples$patient_id[sc_Diseased_Samples$sample_type=='adult_AML'])
    
    # test_annovar_mafs <- readRDS('/data2/rluo4/RPMfunc/Output/scRNA/test_annovar_mafs.rds')
    # View(sc_Diseased_Samples[sc_Diseased_Samples$patient_id %in% sc_adult_AML,])
    # test_annovar_mafs <- test_annovar_mafs[names(test_annovar_mafs) %in% sc_adult_AML]
    # if(cohort %in% split_ddt_ids){
    #  split_ids <- unique(sc_allmeta$split_ids)
    # }
    test_scRNA_annovar <- process_scMapping(scMapping_dir, test_annovar_mafs, cohort)
    
    save(test_scRNA_annovar, file = SNV_data)
    end_time <- Sys.time()
    print(paste("Time taken to process", tis, '--', cohort, ':', end_time - start_time, 'mins'))
  }
}
