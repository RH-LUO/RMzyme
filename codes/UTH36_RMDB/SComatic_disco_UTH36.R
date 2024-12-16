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
# data_path = '/home/lorihan/lrh/All/Output'
data_path = '/data/rluo4/All/Output'
setwd(data_path)
RMP_update <- readxl::read_excel('/data/rluo4/EpiTrans/DataCollection/RMP_update.xlsx')
load(file = '/data/rluo4/RPMfunc/Output/scRNA/sc_allsample_1122.RData', verbose = T)
sc_allmeta <- readRDS(file = '/data/rluo4/RPMfunc/Output/scRNA/sc_allmeta_1122.rds')
RMzyme_alltissues <- sort(unique(sc_allmeta$tissue));
RMzyme_alldatasets <- sort(unique(sc_allmeta$project_id));
extract_last_element <- function(x) {
  split_string <- strsplit(x, "_")
  last_element <- sapply(split_string, function(y) tail(y, n = 1))
  return(last_element)
}
# cell_summary <- read.csv('cell_summary.txt',sep = "\t",header = T)
# sc_majortype <- readRDS(file = '/data/rluo4/RPMfunc/Output/scRNA/sc_majortype_1122.rds')
# tissue_summary <- read.table('/data/rluo4/RPMfunc/PCTanno_pdata/tissue_summary.txt')
################################################################################
# 2) arrange the SNV pdata
################################################################################
# Load Diseased Samples
# Diseased_Samples <- read.table('/data/rluo4/RPMfunc/PCTanno_pdata/Diseased_Samples.txt')
# # Load Normal Samples
# normal_samples <- unique(cell_summary$orig.ident[cell_summary$Tissue %in% c('Healthy','ADJ')])
# remove_samples <- setdiff(Diseased_Samples$orig.ident, tissue_summary$Replicates)# 4
# normal_samples <- c(normal_samples, remove_samples) #417

# sc_summary <- sc_majortype[! duplicated(paste(sc_majortype$patient_id, sc_majortype$tissue, 
#                                               sc_majortype$project_id)),]
# length(unique(sc_summary$project_id)) #396 --> 340
# sc_summary$omics <- 'scRNA'
# sc_Diseased_Samples <- sc_summary[ ! grepl('ealthy',sc_summary$disease), c(2,3,5,6)]
# # View(sc_Diseased_Samples[sc_Diseased_Samples$project_id %in% anno.var,])
# sc_cell_summary <- sc_majortype[! grepl('ealthy', sc_majortype$disease),]
# sc_cell_summary <- sc_cell_summary[ ! sc_cell_summary$project_id %in% tissue_summary$Dataset, ]
# sc_cell_summary$CB <- str_split(sc_cell_summary$barcode, '-', simplify = T)[,1]
# all_characters <- unique(unlist(strsplit(sc_cell_summary$ct, "")))
# potential_separators_and_whitespace <- all_characters[grep("[[:punct:][:space:]]", all_characters)]
# # Define the replacement character
# replacement <- "_"
# # Replace each potential separator with the replacement character
# modified_vector <- sc_cell_summary$ct
# for (sep in potential_separators_and_whitespace) {
#   modified_vector <- gsub(sep, replacement, modified_vector, fixed = TRUE)
# }
# # Replace consecutive occurrences of the replacement character with a single underscore
# sc_cell_summary$Cell_Type <- gsub("_+", replacement, modified_vector)

################################################################################
#           summary of sc_majortype for RMzyme database                        #
################################################################################
length(unique(sc_summary$project_id)) #325 --> 328 --> 332
sample.index <- c(grep("ormal", sc_summary$sample_type),grep("ealthy", sc_summary$sample_type),
                  grep("ontrol", sc_summary$sample_type))
sc_Diseased_Samples <- sc_summary[ -sample.index, c(2,3,5:8)] # 219 --> 222 --> 224 disease datasets
table(sc_Diseased_Samples$disease)
sc_Diseased_Samples <- sc_Diseased_Samples[! grepl('ealthy', sc_Diseased_Samples$disease), ]
# unique(sc_majortype[! grepl('adjacent|ealthy', sc_majortype$disease),]$project_id)
RMzyme_pheno <- as.data.frame(table(sc_summary$sample_type, sc_summary$disease))
RMzyme_pheno <- RMzyme_pheno[RMzyme_pheno$Freq!=0,]
################################################################################
#      organize the diseased data of disco named as sc_cell_summary            #
################################################################################
sc_meta_d <- sc_allmeta[! grepl('adjacent|ealthy', sc_allmeta$disease),] # 222 disease datasets
table(sc_meta_d$ddt_id == sc_meta_d$project_id)
index <- sc_meta_d$ddt_id == sc_meta_d$project_id
# View(sc_meta_d[index,])
# sc_cell_d <- (sc_meta_d[!index, ])
# unique(sc_cell_d$project_id) # GSE202813, GSE185991
# unique(sc_cell_d$ddt_id)

# index_pairs <- paste(sc_summary$project_id, sc_summary$patient_id)
# index_pair <- paste(sc_meta_d$project_id, sc_meta_d$patient_id)[! index]
# View(sc_summary[index_pairs %in% unique(index_pair),])
# View(sc_summary[sc_summary$project_id %in% unique(sc_meta_d$project_id[! index]),])

# table(sc_meta_d$rownames )
# index_GSM <- grepl('--', sc_majortype$barcode)
# View(sc_majortype[! index_GSM,])
# View(sc_summary[duplicated(sc_summary$patient_id),])

setdiff(tissue_summary$Dataset, sc_allsample$project_id)
sc_meta_d <- sc_meta_d%>% dplyr::select(ddt_tis, ddt_id, project_id, patient_id, split_ids,  tissue, disease,
                                        barcode,   Major_type, sub_type, ct, NMF_clusters = seurat_clusters)
sc_meta_d$RMP <- str_split(sc_meta_d$NMF_clusters, "[+]", simplify = T)[,1]
# sc_meta_d$disease <- gsub('adjacent tumor/disease section', 'healthy', sc_meta_d$disease)
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
sc_meta_d <- data.table(sc_meta_d)
sc_meta_d[, clusters := remove_first_element_fast(NMF_clusters)]
table(is.na(sc_meta_d$clusters))
sc_meta_d <- as.data.frame(sc_meta_d)
# View(sc_meta_d[sc_meta_d$sub_type=='Bcells' & sc_meta_d$clusters=='Tcells' ,])
sc_meta_d$cluster_type <- ifelse(sc_meta_d$clusters %in% sc_meta_d$sub_type, 'Major', 'Minor')
# sc_meta_d$disease_type <- sc_meta_d$disease
# sc_meta_d$disease_type[grepl('cancer|carcinoma|leukemia|lymphoma|blastoma|melanoma|myeloma', sc_meta_d$disease_type)] ='cancer'
# sc_meta_d$disease_type[! grepl('cancer|carcinoma|leukemia|lymphoma|blastoma|melanoma|myeloma|healthy', sc_meta_d$disease_type)] ='non-cancer'
# NMF_pheno <- as.data.frame(table(sc_meta_d$disease_type, sc_meta_d$disease))
# NMF_pheno <- NMF_pheno[NMF_pheno$Freq!=0,]
table(unique(sc_meta_d$RMP) %in% RMP_update$RMP)
index <- match(sc_meta_d$RMP, RMP_update$RMP)
sc_meta_d$RNA_modifcations <- RMP_update$RNA_modifications[index]
sc_meta_d$RMP_type <- RMP_update$RMP_type[index]
colnames(sc_meta_d)
################################################################################
extract_last_element <- function(x) {
  split_string <- strsplit(x, "_")
  last_element <- sapply(split_string, function(y) tail(y, n = 1))
  return(last_element)
}
extract_first_element <- function(x) {
  split_string <- strsplit(x, "[--]")
  first_element <- sapply(split_string, function(y) head(y, n = 1))
  return(first_element)
}
################################################################################
sc_cell_summary <- sc_meta_d[ ! sc_meta_d$project_id %in% tissue_summary$Dataset, ] # 169 --> 172 disease datasets
first_elements <- extract_first_element(sc_cell_summary$barcode)
CB = str_split(first_elements,'[-]',simplify = T)[,1]
sc_cell_summary$CB <- CB #str_split(sc_cell_summary$barcode, '-', simplify = T)[,1]
all_characters <- unique(unlist(strsplit(sc_cell_summary$ct, "")))
potential_separators_and_whitespace <- all_characters[grep("[[:punct:][:space:]]", all_characters)]
# Define the replacement character
replacement <- "_"
# Replace each potential separator with the replacement character
modified_vector <- sc_cell_summary$ct
for (sep in potential_separators_and_whitespace) {
  modified_vector <- gsub(sep, replacement, modified_vector, fixed = TRUE)
}
# Replace consecutive occurrences of the replacement character with a single underscore
sc_cell_summary$Cell_Type <- gsub("_+", replacement, modified_vector)
saveRDS(sc_cell_summary, file = '/data/rluo4/RPMfunc/Output/scRNA/sc_cell_summary.rds')
################################################################################
# 3) arrange the SNV data
################################################################################
# CHROM  Start   End     REF     ALT     FILTER  Cell_types      Up_context      Down_context    N_ALT   
# Dp      Nc   -Bc      Cc      VAF     CCF     BCp     CCp     Cell_types_min_BC       Cell_types_min_CC       Rest_BC Rest_CC Fisher_p      Cell_type_Filter        INFO  
# Process scMapping data for disco
# datasets <-  readRDS('/data/rluo4/RPMfunc/disco_pdata/sc_SNV_datasets.rds') #anno.var #unique(Diseased_Samples$Cohort)
# GSE192935 = disco_majortype[disco_majortype$project_id=='GSE192935', ]
# GSE222078 = disco_majortype[disco_majortype$project_id=='GSE222078', ]
# GSE163974 = disco_majortype[disco_majortype$project_id=='GSE163974', ]
# datasets <- list.files('/data/rluo4/EpiTrans/Disco/anno.var')
# table(datasets %in% unique(sc_cell_summary$project_id))
# setdiff(datasets, unique(sc_cell_summary$project_id))#"GSE121636" "GSE140042" "GSE156405"
# # (base) rluo4@sbcsmadlp005:/data/rluo4/EpiTrans/Disco/anno.var/GSE121636$ cd ..
# # (base) rluo4@sbcsmadlp005:/data/rluo4/EpiTrans/Disco/anno.var$ mv GSE121636 GSE121638
# # (base) rluo4@sbcsmadlp005:/data2/rluo4/EpiTrans/Disco/anno.var/GSE158037$ cd ..
# # (base) rluo4@sbcsmadlp005:/data2/rluo4/EpiTrans/Disco/anno.var$ mv GSE158037 GSE158038
# # (base) rluo4@sbcsmadlp005:/data2/rluo4/EpiTrans/Disco/anno.var$ mv GSE162498 GSE162500
# datasets <- datasets[datasets %in% unique(sc_cell_summary$project_id)] #32 --> 33 --> 31

sc_cell_summary <- readRDS("/data/rluo4/RPMfunc/Output/scRNA/sc_cell_summary.rds")
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
  scMapping_dirs <- scMapping_dirs[scMapping_samples %in% sc_adult_AML ]
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
# for(cohort in datasets) {
  cohort = 'GSE185381'
  base_dir <- '/data/rluo4/EpiTrans/Disco'
  start_time <- Sys.time()
  tis <- sc_Diseased_Samples$tissue[match(cohort, sc_Diseased_Samples$project_id)]
  print(paste("Starting parse:", tis, '--', cohort, 'at', start_time))
  
  annovar_dir <- file.path(base_dir, "anno.var", cohort)
  scMapping_dir <- file.path(base_dir, "scMapping", cohort)
  out_dir <- '/data/rluo4/RPMfunc/Output/scRNA/disco_SNV'
  SNV_data <- file.path(out_dir, paste0(cohort, '_SCOMATIC_annovar.RData'))
  
  if (!file.exists(SNV_data)) {
    # test_solutions <- load_and_filter_files(annovar_dir, 'variants.avinput')
    # test_annovar_mafs <- load_and_filter_annovarToMaf(annovar_dir, 'hg38_multianno.txt')
    # test_annovar_mafs <- future_lapply(names(test_annovar_mafs), function(nameindex) {
    #   add_vaf_ccf(test_annovar_mafs, test_solutions, nameindex)
    # })
    # names(test_annovar_mafs) <- names(test_solutions)
    GSE185381_allmeta <- subset(sc_allmeta, ddt_id == 'GSE185381')
    table(GSE185381_allmeta$split_ids)
    GSE185381_allmeta <- subset(GSE185381_allmeta, split_ids == 'sc_adult_AML')
    sc_adult_AML <- unique(GSE185381_allmeta$patient_id) #38
    sc_adult_AML <- intersect(sc_adult_AML, sc_Diseased_Samples$patient_id[sc_Diseased_Samples$sample_type=='adult_AML'])
    test_annovar_mafs <- readRDS('/data/rluo4/RPMfunc/Output/scRNA/test_annovar_mafs.rds')
    # View(sc_Diseased_Samples[sc_Diseased_Samples$patient_id %in% sc_adult_AML,])
    test_annovar_mafs <- test_annovar_mafs[names(test_annovar_mafs) %in% sc_adult_AML]
    test_scRNA_annovar <- process_scMapping(scMapping_dir, test_annovar_mafs, cohort)
    
    save(test_scRNA_annovar, file = SNV_data)
    end_time <- Sys.time()
    print(paste("Time taken to process", tis, '--', cohort, ':', end_time - start_time, 'mins'))
  }
# }

# # load(file.path(data_path, 'cell_summary_diseased.rds'))
# library(data.table)
# library(stringr)
# library(dplyr)
# library(future.apply)
# library(parallel)
# # n_cores <- min(detectCores() - 1, 20) # Set a lower limit, like 20 cores
# plan(multisession, workers = 120)#n_cores)  # Enable parallel processing
# 
# ################################################################################
# # Helper functions
# ################################################################################
# 
# # Load and filter variants files
# load_and_filter_files <- function(dir, pattern) {
#   solution_files <- list.files(dir, pattern = pattern, full.names = TRUE)
#   solutions <- future_lapply(solution_files, fread, sep = '\t', stringsAsFactors = FALSE)
#   names(solutions) <- sapply(solution_files, function(x) strsplit(basename(x), '\\.variants.avinput')[[1]][1])
#   return(solutions)
# }
# 
# # Load and convert Annovar MAF files
# load_and_filter_annovarToMaf <- function(dir, pattern) {
#   maf_files <- list.files(dir, pattern = pattern, full.names = TRUE)
#   solutions <- future_lapply(maf_files, function(fp) annovarToMaf(fp, refBuild = "hg38"))
#   names(solutions) <- sapply(maf_files, function(x) strsplit(basename(x), '\\.hg38_multianno.txt')[[1]][1])
#   return(solutions)
# }
# 
# # Add VAF and CCF columns to MAF data
# add_vaf_ccf <- function(maf, solution, nameindex) {
#   y <- as.data.frame(maf[[nameindex]])
#   z <- gsub('.hg38_multianno.txt', '.variants.avinput', nameindex)
#   w <- as.data.frame(solution[[z]]) 
#   split_v6 <- str_split(w$V6, '[-]', simplify = TRUE)
#   if (nrow(split_v6) == nrow(y)) {
#     y$VAF <- split_v6[, 15]
#     y$CCF <- split_v6[, 16]
#     y$Tumor_Sample_Barcode <- nameindex
#     return(y)
#   } 
# }
# 
# # Apply variant filters
# apply_filters <- function(y) {
#   filter_maf <- unique(c(
#     which(y$X1000g2015aug_all > 0.05),
#     which(y$ExAC_ALL > 0.05),
#     which(y$gnomAD_genome_ALL > 0.05)
#   ))
#   y <- y[-filter_maf,]
#   y <- y[y$VAF != '.' & y$VAF > 0.05, ]
#   filter_repeat <- which(!is.na(y$rmsk) & !is.na(y$genomicSuperDups))
#   y <- y[-filter_repeat, ]
#   return(y)
# }
# 
# ################################################################################
# # Process single-cell mapping data
# ################################################################################
# process_scMapping <- function(dir, annovar_mafs, cohort) {
#   scMapping_dirs <- list.dirs(dir, full.names = TRUE, recursive = FALSE)
#   anno <- vector("list", length(scMapping_dirs))
#   names(anno) <- basename(scMapping_dirs)
#   
#   future_lapply(scMapping_dirs, function(cell_dir) {
#     solution_files <- list.files(cell_dir, pattern = 'single_cell_genotype.tsv', full.names = TRUE)
#     solutions <- lapply(solution_files, fread, stringsAsFactors = FALSE)
#     cell_type_names <- sapply(solution_files, function(x) strsplit(basename(x), '\\.single_cell_genotype.tsv')[[1]][1])
#     names(solutions) <- cell_type_names
#     
#     cell_solutions <- lapply(names(solutions), function(nameindex) {
#       sample_id <- basename(cell_dir) #tail(strsplit(cell_dir, "/")[[1]], 1) 
#       maf_data <- as.data.frame(annovar_mafs[[sample_id]])
#       solution_data <- as.data.frame(solutions[[nameindex]])
#       if (nrow(solution_data) > 0) {
#         colnames(solution_data)[1:2] <- colnames(maf_data)[2:3]
#         solution_data$Start_Position <- as.character(solution_data$Start_Position)
#         
#         combined_data <- left_join(maf_data, solution_data[, -1], by = 'Start_Position')
#         combined_data <- combined_data[!is.na(combined_data$Cell_type_observed), ]
#         # Apply filters
#         combined_data <- apply_filters(combined_data)
#         return(combined_data)  # Return the processed data
#         } else {
#         return(data.frame(matrix(ncol = ncol(maf_data), nrow = 0, dimnames = list(NULL, colnames(maf_data)))))
#       }
#     })
#     
#     names(cell_solutions) <- cell_type_names
#     anno[[basename(cell_dir)]] <- data.table::rbindlist(cell_solutions, fill = TRUE)
#     rm(cell_solutions)  # Free memory
#     gc()
#   })
#   
#   scRNA_annovar <- do.call(rbind, anno)
#   rm(anno)
#   gc()
#   scRNA_annovar
# }

# ################################################################################
# # Main Pipeline Execution
# ################################################################################
# for(cohort in datasets) {
#   base_dir <- '/data/rluo4/EpiTrans/Disco'
#   start_time <- Sys.time()
#   tis <- sc_Diseased_Samples$tissue[match(cohort, sc_Diseased_Samples$project_id)]
#   print(paste("Starting parse:", tis, '--', cohort, 'at', start_time))
#   
#   annovar_dir <- file.path(base_dir, "anno.var", cohort)
#   scMapping_dir <- file.path(base_dir, "scMapping", cohort)
#   out_dir <- '/data/rluo4/RPMfunc/Output/scRNA/disco_SNV'
#   SNV_data <- file.path(out_dir, paste0(cohort, '_SCOMATIC_annovar.RData'))
#   
#   if (!file.exists(SNV_data)) {
#     # Load and process files
#     test_solutions <- load_and_filter_files(annovar_dir, 'variants.avinput')
#     test_annovar_mafs <- load_and_filter_annovarToMaf(annovar_dir, 'hg38_multianno.txt')
#     test_annovar_mafs <- future_lapply(names(test_annovar_mafs), function(nameindex) {
#       add_vaf_ccf(test_annovar_mafs, test_solutions, nameindex)
#     })
#     names(test_annovar_mafs) <- names(test_solutions)
#     
#     # Process single-cell mapping
#     test_scRNA_annovar <- process_scMapping(scMapping_dir, test_annovar_mafs, cohort)
#     
#     # Save final data
#     save(test_scRNA_annovar, file = SNV_data)
#     end_time <- Sys.time()
#     time_taken <- end_time - start_time
#     print(paste("Time taken to process", tis, '--', cohort, ':', time_taken, 'mins'))
#   }
# }
# # Required libraries
# library(future.apply)
# 
# # Set options for future package to handle larger data sizes if needed
# options(future.globals.maxSize = 2 * 1024^3)  # Increase if necessary
# 
# # Define main processing function for each cohort
# process_cohort <- function(cohort, sc_Diseased_Samples, datasets, base_dir, out_dir) {
#   start_time <- Sys.time()
#   tis <- sc_Diseased_Samples$tissue[match(cohort, sc_Diseased_Samples$project_id)]
#   print(paste("Starting parse:", tis, '--', cohort, 'at', start_time))
#   
#   # Define directories
#   annovar_dir <- file.path(base_dir, "anno.var", cohort)
#   scMapping_dir <- file.path(base_dir, "scMapping", cohort)
#   SNV_data <- file.path(out_dir, paste0(cohort, '_SCOMATIC_annovar.RData'))
#   
#   # Skip if data already exists
#   if (!file.exists(SNV_data)) {
#     
#     # Load initial data (filtered as needed)
#     test_solutions <- load_and_filter_files(annovar_dir, 'variants.avinput')
#     test_annovar_mafs <- load_and_filter_annovarToMaf(annovar_dir, 'hg38_multianno.txt')
#     
#     # Process in batches if the data is very large
#     batch_size <- 5  # Adjust based on available memory
#     name_batches <- split(names(test_annovar_mafs), ceiling(seq_along(names(test_annovar_mafs)) / batch_size))
#     
#     # Initialize empty list to store results
#     processed_annovar_mafs <- list()
#     
#     # Process each batch separately
#     for (batch in name_batches) {
#       # For each batch, apply function to each nameindex within the batch
#       batch_results <- future_lapply(batch, function(nameindex) {
#         solution_subset <- test_solutions[[nameindex]]
#         maf_subset <- test_annovar_mafs[[nameindex]]
#         add_vaf_ccf(maf_subset, solution_subset, nameindex)
#       })
#       processed_annovar_mafs <- c(processed_annovar_mafs, batch_results)
#     }
#     
#     # Set names for processed_annovar_mafs as per test_solutions
#     names(processed_annovar_mafs) <- names(test_solutions)
#     
#     # Process single-cell mapping for each cohort with processed results
#     test_scRNA_annovar <- process_scMapping(scMapping_dir, processed_annovar_mafs, cohort)
#     
#     # Save the final processed data to avoid re-processing
#     save(test_scRNA_annovar, file = SNV_data)
#     
#     # Calculate and print processing time
#     end_time <- Sys.time()
#     time_taken <- end_time - start_time
#     print(paste("Time taken to process", tis, '--', cohort, ':', time_taken, 'mins'))
#   }
# }
# 
# # Main pipeline execution
# base_dir <- '/data/rluo4/EpiTrans/Disco'
# out_dir <- '/data/rluo4/RPMfunc/Output/scRNA/disco_SNV'
# 
# # Loop over each cohort in datasets
# for (cohort in datasets) {
#   process_cohort(cohort, sc_Diseased_Samples, datasets, base_dir, out_dir)
# }