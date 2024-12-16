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
# RMBase V3.0
SNV <- fread('/data/rluo4/EpiTrans/DataCollection/hg38.modSNV.tar.gz', sep = '\t')
SNV$V1[1] <- 'SNV_site_1'
bed_df <- readRDS( '/data/rluo4/EpiTrans/DataCollection/RMBase_SNV.rds')

# Get unique cohorts
base_dir <- '/data/rluo4/EpiTrans/Disco/scMapping'
SNV_dir <- '/data/rluo4/RPMfunc/Output/scRNA/PCTanno_SNV' # 89 datasets
combined_cohort <- gsub('_SCOMATIC_annovar.RData','',list.files(SNV_dir))
print(combined_cohort) # 24
# Paths
indir <- '/data/rluo4/RPMfunc/Output/scRNA'
out_dir <- '/data/rluo4/RPMfunc/Output/scRNA/SComatic_PCTanno'
# Load data
large_split <- fread('/data/rluo4/EpiTrans/DataCollection/large_split.txt')
# Convert data.table columns to a vector of dataset IDs
split_ddt_ids <- unlist(large_split, use.names = FALSE)
# Remove empty strings
split_ddt_ids <- split_ddt_ids[split_ddt_ids != ""]


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
process_subdir <- function(cohort, subdir, allmeta, scRNA_annovar, use_split = TRUE) {
  # Subset metadata based on subdir only if use_split is TRUE
  if (use_split) {
    sub_allmeta <- subset(allmeta, split_ids == subdir)
  } else {
    sub_allmeta <- allmeta  # Use all rows for cohorts not in split_ddt_ids
  }
  
  if (nrow(sub_allmeta) == 0) {
    print(paste0(cohort, " (", subdir, "): No matching metadata!"))
    return(NULL)
  }
  
  samples <- unique(sub_allmeta$patient_id)
  samples <- intersect(samples, cell_summary$patient_id)#sc_Diseased_Samples$patient_id)
  
  print(paste0('Diseased Samples: ', paste(samples, collapse = ", ")))
  
  test <- scRNA_annovar[scRNA_annovar$Tumor_Sample_Barcode %in% samples, c(1:18, 130:140)]

  sc_meta_SNV <- cell_summary %>%
    filter(ddt_id == cohort) %>%
    filter(patient_id %in% samples) %>%
    {
      if (use_split) {
        filter(., split_ids == subdir)
      } else {
        .
      }
    } %>%
    mutate(
      CB = paste(CB, patient_id, sep = '--'),
      CB = paste(CB, Cell_Type, sep = '--'),
      CB = paste(CB, project_id, sep = '--')
    )
  
  sc_meta_SNV <- sc_meta_SNV[sc_meta_SNV$CB %in% test$CB, ]
  # Join sc_meta_SNV and stats data
  sc_meta_SNV <- left_join(sc_meta_SNV, test, by = 'CB')
  print(dim(sc_meta_SNV))
  
  # Ensure all chromosome names start with 'chr'
  if (!all(grepl('^chr', sc_meta_SNV$Chromosome))) {
    # Add 'chr' prefix to chromosomes that don't already have it
    sc_meta_SNV$Chromosome <- ifelse(grepl('^chr', sc_meta_SNV$Chromosome),
                                     sc_meta_SNV$Chromosome,
                                     paste0('chr', sc_meta_SNV$Chromosome))
  }
  # Load the GenomicRanges library
  library(GenomicRanges)
  # Print unique chromosomes for confirmation
  print(unique(sc_meta_SNV$Chromosome))
  
  # Create GRanges objects for SNVs and RNA modification sites
  snv_gr <- GRanges(
    seqnames = sc_meta_SNV$Chromosome,
    ranges = IRanges(start = as.numeric(sc_meta_SNV$Start_Position), end = as.numeric(sc_meta_SNV$End_Position)),
    id = sc_meta_SNV$CB
  ) # SComatic callings from scRNA-seq data 
  rna_mod_gr <- GRanges(
    seqnames = bed_df$chrom,
    ranges = IRanges(start = as.numeric(bed_df$chromStart), end = as.numeric(bed_df$chromEnd)),
    id = bed_df$name
  ) # RNA modification sites related to SNVs from RMBase v3.0
  # Find nearest RNA modification sites
  nearest_indices <- nearest(snv_gr, rna_mod_gr)
  # Extract corresponding RNA modification positions
  nearest_mod_sites <- rna_mod_gr[nearest_indices]
  # Calculate shiftPos
  shift_pos <- start(snv_gr) - start(nearest_mod_sites)
  # Combine results into a data frame
  SComatic_res <- data.frame(
    CB = mcols(snv_gr)$id,
    mod_id = mcols(nearest_mod_sites)$id,
    chrom = seqnames(snv_gr),
    snv_pos = start(snv_gr),
    mod_pos = start(nearest_mod_sites),
    shiftPos = shift_pos
  )
  select_cols <- c('ddt_id', 'patient_id', 'split_ids', 'tissue', 'disease',
                   'NMF_clusters', 'RMP', 'clusters', 'RNA_modifcations', 'RMP_type', #'CB', 
                   'Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele', 'Tumor_Seq_Allele2', 'Hugo_Symbol', 'Variant_Classification', 'tx', 'exon', 'txChange', 'aaChange', 'Variant_Type', 'Func.refGene', 'Gene.refGene', 'VAF', 'CCF')
  # join with nmf cluster information
  test_SNV <- sc_meta_SNV[, select_cols]
  SComatic_res <- cbind(SComatic_res, test_SNV)
  print(dim(SComatic_res))
  SComatic_res <- SComatic_res[abs(SComatic_res$shiftPos) <10, ]
  print(dim(SComatic_res))
  # join with RMBase SNV information
  test_RM <- SNV[SNV$V14 %in% SComatic_res$mod_id, c(1,8,9,10,14, 16,19)]
  colnames(test_RM) <- c('snv_id', 'RMBase_geneName', 'RMBase_disease', 'PMID', 'mod_id', 'mod_type', 'RMBase_geneSymbol')
  test_RM$mod_type <- gsub('A-I', 'A-to-I', gsub('Y', 'Psi(Pseudouridine)', test_RM$mod_type))
  
  print(unique(test_RM$mod_type[test_RM$mod_type %in% RMP_update$RNA_modifications]))
  print(setdiff(test_RM$mod_type, test_SNV$RNA_modifcations))
  test_RM$mod_type[ ! test_RM$mod_type %in% test_SNV$RNA_modifcations] <- 'Other'
  
  test_RM$RMBase_geneSymbol <- str_split(test_RM$RMBase_geneSymbol, ',', simplify = T)[,1]
  SComatic_res <- left_join(SComatic_res, test_RM, by = 'mod_id')
  print(dim(SComatic_res))
  
  outputFile <- file.path(out_dir, paste0(cohort, '_', subdir, '_SComatic_maf.rds'))
  saveRDS(SComatic_res, file = outputFile)
  return(SComatic_res)
}
process_cohort <- function(cohort) {
  # Load cohort-specific metadata
  allmeta <- subset(cell_summary, ddt_id == cohort) # sc_allmeta
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
  if (!(cohort %in% split_ddt_ids)) {
    outputFile <- file.path(out_dir, paste0(cohort, '_default_SComatic_maf.rds'))
    if (file.exists(outputFile)) {
      message(paste("Skipping cohort", cohort, "- output already exists."))
      return(NULL)
    }
    
    tryCatch({
      # Skip cohorts with no high FATHMM variants
      scRNA_annovar <- load_and_filter_snv()
      if (is.null(scRNA_annovar)) return(NULL)
      result <- process_subdir(cohort, "default", allmeta, scRNA_annovar, use_split = FALSE)
      if (!is.null(result)) {
        saveRDS(result, file = outputFile)
      }
      return(result)
    }, error = function(e) {
      message(paste("Error processing cohort", cohort, ":", e$message))
      return(NULL)
    })
  } else {
    # Split processing mode
    split_ids <- unique(allmeta$split_ids)
    results <- lapply(split_ids, function(subdir) {
      outputFile <- file.path(out_dir, paste0(cohort, '_', subdir, '_SComatic_maf.rds'))
      if (file.exists(outputFile)) {
        message(paste("Skipping subdir", subdir, "for cohort", cohort, "- output already exists."))
        return(NULL)
      }
      
      tryCatch({
        # Skip cohorts with no high FATHMM variants
        scRNA_annovar <- load_and_filter_snv()
        if (is.null(scRNA_annovar)) return(NULL)
        
        subdir_result <- process_subdir(cohort, subdir, allmeta, scRNA_annovar, use_split = TRUE)
        if (!is.null(subdir_result)) {
          saveRDS(subdir_result, file = outputFile)
        }
        return(subdir_result)
      }, error = function(e) {
        message(paste("Error processing subdir", subdir, "for cohort", cohort, ":", e$message))
        return(NULL)
      })
    })
    
    # Combine results if needed
    combined_result <- do.call(rbind, results)
    return(combined_result)
  }
}

# # test the function
# cohort = 'GSE211630'
# out_dir <- '/data/rluo4/RPMfunc/Output/scRNA/SComatic_PCTanno'
# res <- process_cohort(cohort)
# Run the process for all cohorts
results <- lapply(combined_cohort, process_cohort)
final_results <- do.call(rbind, results)
saveRDS(final_results, file = file.path(indir, 'SComatic_maf_PCTanno.rds'))

# # Split combined cohort into smaller batches
# batch_size <- 10
# cohort_batches <- split(combined_cohort, ceiling(seq_along(combined_cohort) / batch_size))
# 
# # Process each batch
# for (batch_idx in seq_along(cohort_batches)) {
#   message(paste("Processing batch", batch_idx, "of", length(cohort_batches)))
#   batch <- cohort_batches[[batch_idx]]
#   results <- lapply(batch, process_cohort)
#   
#   # Combine and save intermediate results
#   batch_results <- do.call(rbind, results)
#   saveRDS(batch_results, file = file.path(indir, paste0("SComatic_maf_batch_", batch_idx, ".rds")))
# }
# 
# # Combine all batch results
# all_batch_files <- list.files(indir, pattern = "SComatic_maf_batch_.*\\.rds", full.names = TRUE)
# final_results <- do.call(rbind, lapply(all_batch_files, readRDS))
# saveRDS(final_results, file = file.path(indir, "SComatic_maf.rds"))


