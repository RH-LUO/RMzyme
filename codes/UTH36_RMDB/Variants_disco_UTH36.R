#!/usr/local/bin/Rscript
# title: "RMzyme: Regulations of RNA Modifying Enzymes in Human"
# author: "Ruihan Luo"
# date: "Dec 5th,2024"
# rm(list=ls())
options(bitmapType = 'cairo')

library(GenomicRanges)
library(maftools)
library(data.table)
library(stringr)
library(dplyr)
library(future.apply)
library(parallel)

# Optimize future processing
n_cores <- min(detectCores() - 2, 4)
plan(multisession, workers = 24) # n_cores
options(future.globals.maxSize = 24 * 1024^3)  # Increase memory limit

################################################################################
# 1) load pdata of PCTanno & disco database
################################################################################
# load('~/lrh/All/Output/organ13_Epi.RData')# from database-sum.R in ts860
# data_path = '/home/lorihan/lrh/All/Output'
data_path = '/data/rluo4/All/Output'
setwd(data_path)
RMP_update <- readxl::read_excel('/data/rluo4/EpiTrans/DataCollection/RMP_update.xlsx')
# cell_summary <- read.csv('cell_summary.txt',sep = "\t",header = T)
# sc_majortype <- readRDS(file = '/data/rluo4/RPMfunc/Output/scRNA/sc_majortype_1122.rds')
# tissue_summary <- read.table('/data/rluo4/RPMfunc/PCTanno_pdata/tissue_summary.txt')
# sc_allmeta <- readRDS(file = '/data/rluo4/RPMfunc/Output/scRNA/sc_allmeta_1122.rds')
# RMzyme_alltissues <- sort(unique(sc_allmeta$tissue));
# RMzyme_alldatasets <- sort(unique(sc_allmeta$project_id));
load(file = '/data/rluo4/RPMfunc/Output/scRNA/sc_allsample_1122.RData', verbose = T)

################################################################################
# 2) arrange the SNV pdata
################################################################################
################################################################################
#           summary of sc_majortype for RMzyme database                        #
################################################################################
# Parse in required files
library(Seurat)
RMP_update <- readxl::read_excel('/data/rluo4/EpiTrans/DataCollection/RMP_update.xlsx')
sc_cell_summary <- readRDS("/data/rluo4/RPMfunc/Output/scRNA/sc_cell_summary.rds")
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
base_dir <- '/data/rluo4/EpiTrans/Disco/scMapping'
SNV_dir <- '/data/rluo4/RPMfunc/Output/scRNA/disco_SNV' # 89 datasets
combined_cohort <- gsub('_SCOMATIC_annovar.RData','',list.files(SNV_dir))
print(combined_cohort) # 89 -1
# (base) [rluo4@sbimrdplp001 Disco_SNV]$ ls GSE231535_SCOMATIC_annovar.RData -htrl
# -rw-rw-r-- 1 rluo4 rluo4 457K Jul 17 09:23 GSE231535_SCOMATIC_annovar.RData
# (base) [rluo4@sbimrdplp001 Disco_SNV]$ rm GSE231535_SCOMATIC_annovar.RData
# out_dir <- '/data/rluo4/RPMfunc/Output/scRNA/disco_SNV' # 33/55 datasets
# Aggregate mutation data
# SNV_cell_summary <-  sc_cell_summary %>%
#   filter(project_id %in% combined_cohort) 
# table(SNV_cell_summary$sub_type)
# # Epi_subtypes <- SNV_cell_summary$Cell_Type[SNV_cell_summary$sub_type %in% c('Tissue_Specific_Cells','Stem_Cells')]
# # print(unique(Epi_subtypes))
# # Imm_subtypes <- c(#SNV_cell_summary$Cell_Type[SNV_cell_summary$sub_type %in% c('Macrophages','Myeloids','Neutrophils')],
# #   SNV_cell_summary$Cell_Type[grepl('onocyte|Basophil',SNV_cell_summary$Cell_Type)])
# # print(unique(Imm_subtypes))
# # Str_subtypes <- SNV_cell_summary$Cell_Type[grepl('Adipocyte',SNV_cell_summary$Cell_Type)]
# # print(unique(Str_subtypes))
# # Other_subtypes <- SNV_cell_summary$Cell_Type[SNV_cell_summary$sub_type %in% c('Neuro_Cells')]
# # print(unique(Other_subtypes))
# Mut_celltypes <- SNV_cell_summary$Cell_Type[ ! SNV_cell_summary$sub_type %in% c('Bcells','Tcells')]#c(Epi_subtypes, Imm_subtypes, Str_subtypes, Other_subtypes)
# Initialize an empty data frame to store merged data (Done in terminal, skip)
load('/data/rluo4/RPMfunc/Output/scRNA/disco_SNV/GSE157376_SCOMATIC_annovar.RData')
MAF_cols <- colnames(test_scRNA_annovar) # UTH36_MAF_version: 140 columns

maf_results <-  vector("list", length(combined_cohort))
names(maf_results) <- combined_cohort
maf <- data.frame()
# Loop through each cohort to load and combine data
for (cohort in combined_cohort) {
  SNV_data <- file.path(SNV_dir, paste0(cohort, '_SCOMATIC_annovar.RData'))
  # Check if file exists before loading
  if (!file.exists(SNV_data)) next
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
    
    # Add cohort information to the CB column
    scRNA_annovar$CB <- paste(scRNA_annovar$CB, cohort, sep = '--')
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
saveRDS(maf, file = file.path(SNV_dir, '../disco_maf.rds'))
save(maf_results, file = file.path(SNV_dir, '../disco_maf_results.RData'))

# maf <- readRDS(file.path(out_dir, '../disco_maf.rds'))
# load(file = file.path(out_dir, '../disco_maf_results.RData'))
# 
# cohort = 'GSE185381'
# SNV_data <- file.path(out_dir, paste0(cohort, '_SCOMATIC_annovar.RData'))
# 
# # Check if file exists before loading
# # if (!file.exists(SNV_data)) next
# 
# # Load the data
# load(SNV_data)
# test_scRNA_annovar <- as.data.frame(test_scRNA_annovar)
# 
# # Filter by relevant columns and mutation cell types
# valid_columns <- intersect(colnames(test_scRNA_annovar), MAF_cols)
# scRNA_annovar <- test_scRNA_annovar[, valid_columns]
# rm(test_scRNA_annovar); gc()
# # scRNA_annovar <- scRNA_annovar[scRNA_annovar$Cell_type_observed %in% Mut_celltypes, ]
# 
# # Filter by FATHMM score
# high_FATHMM <- !is.na(scRNA_annovar$FATHMM_score) & scRNA_annovar$FATHMM_score > 0.7
# scRNA_annovar <- scRNA_annovar[high_FATHMM, ]
# 
# # Log filtered data information
# print(paste0(cohort, ":", sum(high_FATHMM), " variants with FATHMM_score >0.7"))
# if (nrow(scRNA_annovar) == 0) {
#   print(paste0(cohort, ": no mutation data for Mut_celltypes!"))
#   next
# }
# scRNA_annovar$CB <- paste(scRNA_annovar$CB,  cohort, sep = '--')
# # Append filtered data to the master data frame
# maf <- rbind(maf, scRNA_annovar)
# maf_results[[cohort]] <- scRNA_annovar
# rm(scRNA_annovar); gc()
# 
# saveRDS(maf, file = file.path(out_dir, '../disco_maf_32ddts.rds'))
# save(maf_results, file = file.path(out_dir, '../disco_maf_results_32ddts.RData'))

stats <- maf %>%
  filter(ExonicFunc.refGene %in% c('nonsynonymous SNV')) %>%
  group_by(CB) %>%
  summarise(AAChange.refGene = list(unique(AAChange.refGene))) %>%
  dplyr::mutate(Mut_count = lengths(AAChange.refGene),
                Mut_gene = sapply(AAChange.refGene, function(x) paste(unique(str_split(x, ':', simplify = TRUE)[, 1]), collapse = ', ')),
                AAChange = sapply(AAChange.refGene, function(x) paste(x, collapse = ', ')))
saveRDS(stats, file = file.path(SNV_dir, '../Variants_disco.rds')) # this version skip all bone marrow datasets

