#!/usr/local/bin/Rscript
# title: "RMzyme: Regulations of RNA Modifying Enzymes in Human"
# author: "Ruihan Luo"
# date: "Dec 5th,2024"
# rm(list=ls())
options(bitmapType = 'cairo')

library(GenomicRanges)
library(data.table)
library(stringr)
library(dplyr)
library(future.apply)
library(parallel)

# Optimize future processing
n_cores <- min(detectCores() - 2, 4)
plan(multisession, workers = n_cores)
options(future.globals.maxSize = 8 * 1024^3)  # Increase memory limit

# Parse in required files
library(Seurat)
celltype_summary <- read.csv('/data/rluo4/summary/Annotation/celltype_summary.txt',sep = "\t",header = T, fill = T)
MisMatched_Samples_SNV <- read.table('/data/rluo4/RPMfunc/PCTanno_pdata/MisMatched_Samples_SNV.txt', header = T)
RMP_update <- readxl::read_excel('/data/rluo4/EpiTrans/DataCollection/RMP_update.xlsx')
cell_summary <- readRDS("/data/rluo4/RPMfunc/Output/scRNA/cell_summary.rds")
# sc_allmeta <- readRDS(file = '/data/rluo4/RPMfunc/Output/scRNA/sc_allmeta_1122.rds')#sc_allmeta_0925.rds') 
load(file = '/data/rluo4/RPMfunc/Output/scRNA/sc_allsample_1122.RData', verbose = T)
length(unique(sc_summary$project_id)) #325 --> 328 --> 332
sample.index <- c(grep("ormal", sc_summary$sample_type),grep("ealthy", sc_summary$sample_type),
                  grep("ontrol", sc_summary$sample_type))
sc_Diseased_Samples <- sc_summary[ -sample.index, c(2,3,5:8)] # 219 --> 222 --> 224 disease datasets
table(sc_Diseased_Samples$disease)
sc_Diseased_Samples <- sc_Diseased_Samples[! grepl('ealthy', sc_Diseased_Samples$disease), ]
###################################
#   Merge all scRNA SNV Results   #
###################################
# Get unique cohorts
base_dir <- '/public/lorihan/HPC'
SNV_dir <- '/data/rluo4/RPMfunc/Output/scRNA/PCTanno_SNV' # 89 datasets
combined_cohort <- gsub('_SCOMATIC_annovar.RData','',list.files(SNV_dir))
print(combined_cohort) # 32
# Subdirectory processing function
remove_last_element <- function(char_vector) {
  char_vector_parts <- strsplit(char_vector, "--")
  char_vector_parts <- lapply(char_vector_parts, function(parts) {
    if (length(parts) > 1) {
      return(paste(parts[-length(parts)], collapse = "--"))
    } else {
      return(parts[1]) # If there's only one part, return it as is.
    }
  })
  return(unlist(char_vector_parts))
}
extract_last_element <- function(x) {
  split_string <- strsplit(x, "--")
  last_element <- sapply(split_string, function(y) tail(y, n = 1))
  return(last_element)
}
# Initialize an empty data frame to store merged data (Done in terminal, skip)
maf_results <-  vector("list", length(combined_cohort))
names(maf_results) <- combined_cohort
maf <- data.frame()
# Loop through each cohort to load and combine data
for (cohort in combined_cohort) {
  SNV_data <- file.path(SNV_dir, paste0(cohort, '_SCOMATIC_annovar.RData'))
  
  # Load reference data for column validation
  load(file.path(SNV_dir, 'GSE181919_SCOMATIC_annovar.RData'))
  MAF_cols <- colnames(test_scRNA_annovar) # Columns of reference data (e.g., 140 columns)
  
  # Helper: Load and filter SNV data
  load_and_filter_snv <- function() {
    load(SNV_data)
    test_scRNA_annovar <- as.data.frame(test_scRNA_annovar)
    
    # Retain only valid columns
    valid_columns <- intersect(colnames(test_scRNA_annovar), MAF_cols)
    scRNA_annovar <- test_scRNA_annovar[, valid_columns]
    rm(test_scRNA_annovar); gc()
    
    # Apply FATHMM score filter
    high_FATHMM <- !is.na(scRNA_annovar$FATHMM_score) & scRNA_annovar$FATHMM_score > 0.7
    scRNA_annovar <- scRNA_annovar[high_FATHMM, ]
    print(paste0(cohort, ":", sum(high_FATHMM), " variants with FATHMM_score > 0.7"))
    
    if (nrow(scRNA_annovar) == 0) {
      print(paste0(cohort, ": no mutation data for Mut_celltypes!"))
      return(NULL)
    }
    short_ct <- scRNA_annovar$Cell_type_observed
    test_ct <- extract_last_element(scRNA_annovar$CB)
    print(table(short_ct==test_ct))
    full_ct <- celltype_summary$FullName[match(test_ct, celltype_summary$Minor.Celltype)]
    all_characters <- unique(unlist(strsplit(full_ct, "")))
    potential_separators_and_whitespace <- all_characters[grep("[[:punct:][:space:]]", all_characters)]
    # Define the replacement character
    replacement <- "_"
    # Replace each potential separator with the replacement character
    modified_vector <- full_ct
    for (sep in potential_separators_and_whitespace) {
      modified_vector <- gsub(sep, replacement, modified_vector, fixed = TRUE)
    }
    # Replace consecutive occurrences of the replacement character with a single underscore
    full_ct <- gsub("_+", replacement, modified_vector)
    scRNA_annovar$Cell_Type <- full_ct
    scRNA_annovar$CB <- remove_last_element(scRNA_annovar$CB)
    # Add cohort information to the CB column
    if(cohort %in% unique(MisMatched_Samples_SNV$project_id)){
      scRNA_annovar$CB <- remove_last_element(scRNA_annovar$CB)
      idx_match <- match(scRNA_annovar$Tumor_Sample_Barcode, MisMatched_Samples_SNV$patient_id)
      scRNA_annovar$Tumor_Sample_Barcode <- MisMatched_Samples_SNV$sample_id[idx_match]
      scRNA_annovar$CB <- paste(scRNA_annovar$CB, scRNA_annovar$Tumor_Sample_Barcode, full_ct, cohort, sep = '--')
    } else{
      scRNA_annovar$CB <- paste(scRNA_annovar$CB,full_ct, cohort, sep = '--')
    }
    return(scRNA_annovar)
  }
  
  # Default processing mode
    tryCatch({
      # Skip cohorts with no high FATHMM variants
      scRNA_annovar <- load_and_filter_snv()
    }
    )
  # Append filtered data to the master data frame
  maf <- rbind(maf, scRNA_annovar)
  maf_results[[cohort]] <- scRNA_annovar
  rm(scRNA_annovar); gc()
}
saveRDS(maf, file = file.path(SNV_dir, '../PCTanno_maf.rds'))
save(maf_results, file = file.path(SNV_dir, '../PCTanno_maf_results.RData'))

stats <- maf %>%
  filter(ExonicFunc.refGene %in% c('nonsynonymous SNV')) %>%
  group_by(CB) %>%
  summarise(AAChange.refGene = list(unique(AAChange.refGene))) %>%
  dplyr::mutate(Mut_count = lengths(AAChange.refGene),
                Mut_gene = sapply(AAChange.refGene, function(x) paste(unique(str_split(x, ':', simplify = TRUE)[, 1]), collapse = ', ')),
                AAChange = sapply(AAChange.refGene, function(x) paste(x, collapse = ', ')))
saveRDS(stats, file = file.path(SNV_dir, '../Variants_PCTanno.rds')) # this version skip all bone marrow datasets

