#!/usr/local/bin/Rscript
# title: "RMzyme: Regulations of RNA Modifying Enzymes in Human"
# author: "Ruihan Luo"
# date: "December 9th,2024"
# rm(list=ls())
options(bitmapType = 'cairo')
################################################################################
# 0) visualize the SNV data of 13 organs in PCTanno
################################################################################
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
# 3) organize the SNV data of RMzyme
################################################################################
# On UTH145: nohup Rscript ../Variants_disco_UTH145.R >./Variants_disco_UTH145.R.o & 
# nohup Rscript ../Variants_RMzyme.R >./Variants_RMzyme.R.o &
# nohup Rscript ../Variants_RMzyme_sort.R >./Variants_RMzyme_sort.R.o  &

# On UTH36: nohup Rscript ../Variants_PCTanno_UTH36.R > ./Variants_PCTanno_UTH36.R.o &
# nohup Rscript ../Variants_PCTanno.R >> ./Variants_PCTanno.R.o &
################################################################################
# 4) Summarize the SNV data of 120 SNV datasets 
################################################################################
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(ggpubr)
library(ggsignif)
library(ggstatsplot)
library(ggplot2)
library(ggalluvial)
library(RColorBrewer)
library(ggthemes)
library(maftools)
###################################
#   Merge all scRNA SNV Results   #
##########################################
# 4.1 Base directory and output directory
##########################################
base_dir <- '/data/rluo4/RPMfunc/SCOMATIC'
out_dir <- '/data/rluo4/RPMfunc/Output/scRNA/disco_SNV' # 88 --> 87 datasets(exclude the AML--GSE231535)
# Get unique cohorts
combined_cohort <- gsub('_SCOMATIC_annovar.RData','',list.files(out_dir))
print(combined_cohort)
options(bitmapType = 'cairo')

RMP_update <- readxl::read_excel('/data/rluo4/EpiTrans/DataCollection/RMP_update.xlsx')
sc_cell_summary <- readRDS("/data/rluo4/RPMfunc/Output/scRNA/sc_cell_summary.rds")
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
load('/data/rluo4/EpiTrans/RMDatasets/EpiTrans.RData', verbose = T)
NMF_clusters <- readRDS('/data/rluo4/RPMfunc/Output/scRNA/NMF_clusters.rds')

##########################################
# 4.2 Filter and prepare SNV_cell_summary:
##########################################
# SNV_cell_summary <-  sc_cell_summary %>%
#   filter(project_id %in% combined_cohort) #%>% filter(!patient_id %in% unique(AML$patient_id))
# table(SNV_cell_summary$sub_type) # 16 sub_types (except for Germ_Cells: sc_majortype[sc_majortype$sub_type=='Germ_Cells',] all healthy)
# unique(SNV_cell_summary$project_id) #31
# unique(SNV_cell_summary$patient_id) #156
# 
# UTH36_scMapping <- list.files('/data/rluo4/EpiTrans/Disco/scMapping/')
# GSM_scMapping <- NULL
# for (i in UTH36_scMapping) {
#   j <- list.files(file.path('/data/rluo4/EpiTrans/Disco/scMapping/', i))
#   GSM_scMapping <<- c(GSM_scMapping, j)
# }
# table(GSM_scMapping %in% SNV_cell_summary$patient_id)
# AML <- readRDS(file=paste0('/data/rluo4/RPMfunc/disco_pdata/AML_meta.rds'))
# table(unique(SNV_cell_summary$patient_id) %in% AML$patient_id)

# Epi_subtypes <- SNV_cell_summary$Cell_Type[SNV_cell_summary$sub_type %in% c('Tissue_Specific_Cells','Stem_Cells')]
# print(unique(Epi_subtypes))
# Imm_subtypes <- c(#SNV_cell_summary$Cell_Type[SNV_cell_summary$sub_type %in% c('Macrophages','Myeloids','Neutrophils')],
#   SNV_cell_summary$Cell_Type[grepl('onocyte|Basophil',SNV_cell_summary$Cell_Type)])
# print(unique(Imm_subtypes))
# Str_subtypes <- SNV_cell_summary$Cell_Type[grepl('Adipocyte',SNV_cell_summary$Cell_Type)]
# print(unique(Str_subtypes))
# Other_subtypes <- SNV_cell_summary$Cell_Type[SNV_cell_summary$sub_type %in% c('Neuro_Cells')]
# print(unique(Other_subtypes))
# Mut_celltypes <- SNV_cell_summary$Cell_Type[ ! SNV_cell_summary$sub_type %in% c('Bcells','Tcells')]#c(Epi_subtypes, Imm_subtypes, Str_subtypes, Other_subtypes)

# Aggregate mutation data
## ( On Server UTH145: nohup Rscript ../Summary_disco_UTH36.R >./Summary_disco_UTH36.R.o & since Summary_disco.R is not working neither in 36 or 57)
# On Server UTH36: nohup Rscript ../Variants_disco_UTH36.R >./Variants_disco_UTH36.R.o &
# nohup Rscript ../Variants_PCTanno_UTH36.R >./Variants_PCTanno_UTH36.R.o &
# ...
# # Loop through each cohort to load and combine data
# ...
# saveRDS(maf, file = file.path(out_dir, '../disco_maf.rds'))
# save(maf_results, file = file.path(out_dir, '../disco_maf_results.RData'))

#########################################################################
# 4.3 summary of SNV samples for RMzyme database                        #
#########################################################################
# unique(sc_majortype[! grepl('adjacent|ealthy', sc_majortype$disease),]$project_id)
RMzyme_pheno <- as.data.frame(table(sc_summary$sample_type, sc_summary$disease))
RMzyme_pheno <- RMzyme_pheno[RMzyme_pheno$Freq!=0,]
Diseased_Samples_UTH36 <- read.table('/data/rluo4/RPMfunc/disco_pdata/Diseased_Samples_UTH36.txt') #32 = 31 + 1(GSE185381) cohorts completed in UTH36
Diseased_Samples <- read.table('/data/rluo4/RPMfunc/disco_pdata/Diseased_Samples_UTH57_renew2.txt')
table(Diseased_Samples$patient_id %in% Diseased_Samples_UTH36$patient_id )
Diseased_Samples_disco <- rbind(Diseased_Samples_UTH36, Diseased_Samples)
table( unique(Diseased_Samples_disco$project_id) %in% combined_cohort)
Diseased_Samples_disco <- Diseased_Samples_disco[Diseased_Samples_disco$project_id %in% combined_cohort, ]

setdiff(Diseased_Samples_disco$project_id, sc_cell_summary$project_id)
table(Diseased_Samples_disco$patient_id %in% sc_cell_summary$patient_id) #703 --> 766 samples from SNV calling --> 586
table(Diseased_Samples_disco$patient_id %in% cell_summary$patient_id) #
unique(Diseased_Samples_disco$project_id)#88
unique(Diseased_Samples_disco$patient_id)#610

Diseased_Samples_PCTanno <- read.csv('/data/rluo4/RPMfunc/PCTanno_pdata/Diseased_Samples_PCTanno.txt',sep = '\t', fill = T)
out_dir <- '/data/rluo4/RPMfunc/Output/scRNA/PCTanno_SNV/' # 55 datasets
# Get unique cohorts
combined_cohort_PCTanno <- gsub('_SCOMATIC_annovar.RData','',list.files(out_dir, pattern = '.RData'))
print(combined_cohort_PCTanno) # 32
table(unique(Diseased_Samples_PCTanno$project_id) %in% combined_cohort_PCTanno)
same_cols <- intersect(colnames(Diseased_Samples_PCTanno),colnames(Diseased_Samples_disco))

Diseased_Samples_SNV <- rbind(Diseased_Samples_PCTanno[,same_cols], 
                              Diseased_Samples_disco[,same_cols])
unique(Diseased_Samples_SNV$project_id)#120
unique( paste(Diseased_Samples_SNV$project_id, Diseased_Samples_SNV$patient_id))#1376
table(Diseased_Samples_SNV$patient_id %in% sc_allsample$patient_id)
setdiff(Diseased_Samples_SNV$patient_id, sc_allsample$patient_id)
# View(Diseased_Samples_SNV[ ! Diseased_Samples_SNV$patient_id %in% sc_allsample$patient_id, ])
MisMatched_Samples_SNV <- Diseased_Samples_SNV[ ! Diseased_Samples_SNV$patient_id %in% sc_allsample$patient_id, ]
MisMatched_Samples_SNV <- MisMatched_Samples_SNV[MisMatched_Samples_SNV$project_id %in% tissue_summary$Dataset,]
# MisMatched_Samples_SNV <- MisMatched_Samples_SNV[grepl('home', MisMatched_Samples_SNV$SNV_path),]
MisMatched_Samples_SNV$sample_id <- c('GSM5258385', 'GSM5258386', 'GSM5258387', 'GSM5258388', 'GSM5258389', 'GSM5258390',
                                      'P1', 'P2', 'P3', 'P4', 'P5', 'P6',
                                      'P1', 'P3', 'P5', 'P6', 'P7', 'P8', 'P9')
# write.table(MisMatched_Samples_SNV, file = '/data/rluo4/RPMfunc/PCTanno_pdata/MisMatched_Samples_SNV.txt',row.names = F, quote = F)

################################################################################
# 4.4 Arrange the SNV results                                                  #
################################################################################
# Paths
indir <- '/data/rluo4/RPMfunc/Output/scRNA'
out_dir <- '/data/rluo4/RPMfunc/Output/scRNA/SComatic_res'
# Load data
large_split <- fread('/data/rluo4/EpiTrans/DataCollection/large_split.txt')
# Convert data.table columns to a vector of dataset IDs
split_ddt_ids <- unlist(large_split, use.names = FALSE)
# Remove empty strings
split_ddt_ids <- split_ddt_ids[split_ddt_ids != ""]

# combine the maf results for RMzyme
maf_PCTanno <- readRDS(file.path(out_dir, '../PCTanno_maf.rds'))
unique(maf_PCTanno$Tumor_Sample_Barcode)# 646
maf_PCTanno$AAChange.refGene[is.na(maf_PCTanno$AAChange.refGene)] <- maf_PCTanno$Hugo_Symbol[is.na(maf_PCTanno$AAChange.refGene)]
table(maf_PCTanno$ExonicFunc.refGene)#外显子区的SNV or InDel变异类型
maf_disco <- readRDS(file.path(out_dir, '../disco_maf.rds'))
unique(maf_disco$Tumor_Sample_Barcode)# 572
# Handle NA values in AAChange.refGene
maf_disco$AAChange.refGene[is.na(maf_disco$AAChange.refGene)] <- maf_disco$Hugo_Symbol[is.na(maf_disco$AAChange.refGene)]
table(maf_disco$ExonicFunc.refGene)#外显子区的SNV or InDel变异类型
# nonsynonymous SNV          stopgain          stoploss
# 907474               585                25
# synonymous SNV           unknown
# 567              3744
maf <- rbind(maf_PCTanno[, -ncol(maf_PCTanno)], maf_disco)
maf$batch <- c(rep('PCTanno', nrow(maf_PCTanno)), rep('disco', nrow(maf_disco)))
extract_last_element <- function(x) {
  split_string <- strsplit(x, "--")
  last_element <- sapply(split_string, function(y) tail(y, n = 1))
  return(last_element)
}
saveRDS(maf, file = file.path(indir, "RMzyme_maf.rds"))                     
cohorts <- extract_last_element(maf$CB)
length(unique(cohorts)) # 120
identifier <- paste(cohorts, maf$Tumor_Sample_Barcode)
length(unique(identifier)) # 1227
SNV_samples <- unique( paste(Diseased_Samples_SNV$project_id, Diseased_Samples_SNV$patient_id))#1376
setdiff(SNV_samples, unique(identifier))

#################################################
# 4.5 RNA Modifications related to SComatic Var #
#################################################
# Combine all batch results for RM-related SNV
# Get unique PCTanno cohorts
SNV_dir <- '/data/rluo4/RPMfunc/Output/scRNA/PCTanno_SNV' # 89 datasets
combined_cohort <- gsub('_SCOMATIC_annovar.RData','',list.files(SNV_dir, pattern = 'RData'))
print(combined_cohort) # 32
out_dir <- '/data/rluo4/RPMfunc/Output/scRNA/SComatic_PCTanno'
all_batch_files <- list.files(out_dir, pattern = "*\\_SComatic_maf.rds", full.names = TRUE)
all_batch <- gsub('/data/rluo4/RPMfunc/Output/scRNA/SComatic_PCTanno/', '', all_batch_files)
all_batch <- str_split(all_batch, '_', simplify = T)[,1]
setdiff(combined_cohort,all_batch)
# [1] "GSE144236" "GSE163203" "GSE167297" "GSE172577" "GSE181254" 
# [6] "HTA11"  
# 1) no Epithelial cells in sc_allmeta for GSE172577 and GSE181254
# 2) no SNV for GSE167297 cells
PCTanno_RM_SNV <- do.call(rbind, lapply(all_batch_files, readRDS))
PCTanno_RM_SNV$batch <- 'PCTanno'
saveRDS(PCTanno_RM_SNV, file = file.path(indir, "SComatic_maf_PCTanno.rds"))

# Get unique disco cohorts
SNV_dir <- '/data/rluo4/RPMfunc/Output/scRNA/disco_SNV' # 89 datasets
combined_cohort <- gsub('_SCOMATIC_annovar.RData','',list.files(SNV_dir))
print(combined_cohort) # 88 = 89 -1
# -rw-rw-r-- 1 rluo4 rluo4  262 Dec  5 13:10 GSE181279_SCOMATIC_annovar.RData
# rm GSE181279_SCOMATIC_annovar.RData
combined_cohort <- sort(combined_cohort, decreasing = T)
print(combined_cohort) # 88
out_dir <- '/data/rluo4/RPMfunc/Output/scRNA/SComatic_res'
all_batch_files <- list.files(out_dir, pattern = "*\\_SComatic_maf.rds", full.names = TRUE)
all_batch <- gsub('/data/rluo4/RPMfunc/Output/scRNA/SComatic_res/', '', all_batch_files)
all_batch <- str_split(all_batch, '_', simplify = T)[,1]
all_batch <- unique(all_batch)
setdiff(combined_cohort,all_batch)
# character(0)
disco_RM_SNV <- do.call(rbind, lapply(all_batch_files, readRDS))
disco_RM_SNV$batch <- 'disco'
saveRDS(disco_RM_SNV, file = file.path(indir, "SComatic_maf_disco.rds"))
final_RM_SNV <- rbind(PCTanno_RM_SNV, disco_RM_SNV)
saveRDS(final_RM_SNV, file = file.path(indir, "final_RM_SNV.rds"))

################################################################################
# Parse in result files
################################################################################

#######################################
# 4.6 Summarize Coverage for samples  #
#######################################
# Coverage data as provided
# coverage_data <- read_delim('/data/rluo4/RPMfunc/PCTanno_pdata/HSIL_HPV_1.depth.chr.stat.gz')
# coverage_data <- coverage_data[!is.na(coverage_data$MeanDepth), ]
# coverage_data[, -1] <- apply(coverage_data[, -1], 2, function(x){
#   y <- as.numeric(x)
#   return(y)
# })
# indir='/data/rluo4/RPMfunc/disco_pdata/disco_SNV_Coverage/'
# tissues <- list.files(indir)
# load_files <- function(dir, normal_samples, pattern) {
#   solution <- list.files(dir, pattern = pattern)
#   # nor <- paste0(normal_samples, '.', pattern)  # Adjusted the separator
#   # solution <- solution[!solution %in% nor] #normal_samples is not fit here since it stems from cell_summry['orig.ident']
#   file_paths <- file.path(dir, solution)
#   # Use fread correctly within lapply
#   solutions <- lapply(file_paths, function(fp) {
#     fread(fp, sep = '\t', stringsAsFactors = FALSE)
#   })
#   names(solutions) <- sapply(solution, function(x) strsplit(x, '\\.depth.chr.stat.gz')[[1]][1])
#   return(solutions)
# }
# coverage_data <- list()
# for (tis in tissues) {
#   cov_dir <- file.path(indir, tis)
#   coverage_data[[tis]] <- load_files(cov_dir, normal_samples, 'depth.chr.stat.gz')
# }
# total_covered_bp <- data.frame()
# for (tis in tissues) {
#   cov_data <- coverage_data[[tis]]
#   df_cov <- lapply(cov_data, function(x){
#     return( sum(x$CoveredSite) )
#   })
#   df_cov <- data.frame(total_covered_mb = unlist(df_cov))
#   df_cov$Sample <- rownames(df_cov)
#   if(nrow(df_cov) == 0){
#     print(paste0('no coverage data for ', tis))
#     next;
#   }
#   # index  = tis=='Prostate' & rownames(df_cov) %in% c('P1','P2', 'P3', 'P4', 'P5', 'P6')
#   # rownames(df_cov)[index] <- paste0('Dong_',  rownames(df_cov)[index])
#   df_cov$Tissue <- tis
#   total_covered_bp <- rbind(total_covered_bp, df_cov)
# }
# # [1] "no coverage data for GSE121638"
# # [1] "no coverage data for GSE162500"
# table(total_covered_bp$Sample %in% Diseased_Samples_disco$patient_id)
# setdiff(total_covered_bp$Sample, Diseased_Samples_disco$patient_id)
# setdiff( Diseased_Samples_disco$patient_id, total_covered_bp$Sample)
# Diseased_Samples_disco$Total_Coverage <- total_covered_bp$total_covered_mb[match(Diseased_Samples_disco$patient_id, total_covered_bp$Sample)]
# Diseased_Samples_disco$Total_covered_mb <- Diseased_Samples_disco$Total_Coverage / 1e6
# # write.table(Diseased_Samples_disco, sep = "\t", file = '/data/rluo4/RPMfunc/disco_pdata/Diseased_Samples_disco.txt', quote = F)

outdir = '/data/rluo4/RPMfunc/Output/scRNA'
setwd(outdir)
# TMB per cell #
Diseased_Samples_PCTanno <- read.csv('/data/rluo4/RPMfunc/PCTanno_pdata/Diseased_Samples_PCTanno.txt',sep = '\t', fill = T)
Diseased_Samples_disco <- read.table('/data/rluo4/RPMfunc/disco_pdata/Diseased_Samples_disco.txt') # UTH36 + UTH57 from UTH57's DiscoDatasets.R
same_cols <- intersect(colnames(Diseased_Samples_PCTanno),colnames(Diseased_Samples_disco))

Diseased_Samples_SNV <- rbind(Diseased_Samples_PCTanno[,same_cols], 
                              Diseased_Samples_disco[,same_cols])
table(unique(maf$Tumor_Sample_Barcode) %in% Diseased_Samples_SNV$patient_id)
setdiff(unique(maf$Tumor_Sample_Barcode), Diseased_Samples_SNV$patient_id)
idx <- Diseased_Samples_SNV$project_id %in% unique(MisMatched_Samples_SNV$project_id)
idx_match <- match(Diseased_Samples_SNV$patient_id[idx], MisMatched_Samples_SNV$patient_id)
Diseased_Samples_SNV$patient_id[idx] <- MisMatched_Samples_SNV$sample_id[idx_match]
setdiff(unique(maf$Tumor_Sample_Barcode), Diseased_Samples_SNV$patient_id)
# write.table(Diseased_Samples_SNV, sep = "\t", file = '/data/rluo4/RPMfunc/Output/summary/Diseased_Samples_SNV.txt', quote = F)

# Extract the part between the first and second occurrence of "--"
# stats_disco <- readRDS('/data/rluo4/RPMfunc/Output/scRNA/Variants_disco.rds') # generated by Variant_disco_UTH36.R
# stats_PCTanno <- readRDS('/data/rluo4/RPMfunc/Output/scRNA/Variants_PCTanno.rds') # generated by Variant_disco_UTH36.R
# stats <- rbind(stats_PCTanno, stats_disco)
# stats$batch <- c(rep('PCTanno', nrow(stats_PCTanno)), rep('disco', nrow(stats_disco)))
# stats$Tumor_Sample_Barcode <- str_extract(stats$CB, "(?<=--)[^--]*(?=--)")
stats <- maf %>%
  # filter(ExonicFunc.refGene %in% c('nonsynonymous SNV')) %>%
  group_by(CB) %>%
  summarise(AAChange.refGene = list(unique(AAChange.refGene))) %>%
  dplyr::mutate(Mut_count = lengths(AAChange.refGene),
                Mut_gene = sapply(AAChange.refGene, function(x) paste(unique(str_split(x, ':', simplify = TRUE)[, 1]), collapse = ', ')),
                AAChange = sapply(AAChange.refGene, function(x) paste(x, collapse = ', ')))
saveRDS(stats, file = file.path(SNV_dir, '../Variants_stats.rds')) # this version skip all bone marrow datasets

library(stringr)
split_CB <- str_split(stats$CB, "--", simplify = TRUE)[, 2]
index = grepl("cd45---*?", stats$CB) #145 CB
split_CB[index]
# Further split the second part on the last occurrence of '_'
split_CB[index] <- paste0(split_CB[index], '-')#str_split_fixed(stats$CB[index], '---', 3)[, 2] #str_extract(stats$CB[index], "cd45-.*?")
# Assign the extracted part to Tumor_Sample_Barcode
stats$Tumor_Sample_Barcode <- split_CB
setdiff(unique(maf$Tumor_Sample_Barcode), unique(split_CB)) # "GSM3516662" "GSM3972020" "GSM4339781"
stats$project_id <- str_split(stats$CB, "--", simplify = TRUE)[, 4]
dim(stats)#1119322       8
# # Apply the function to the CB column
# stats$Tumor_Sample_Barcode <- sapply(stats$CB, extract_sample_barcode)
sc_meta_SNV <- rbind(sc_cell_summary, cell_summary) %>%
  mutate(CB = paste(CB, patient_id, sep = '--'),
         CB = paste(CB, Cell_Type, sep = '--'),
         CB = paste(CB, ddt_id, sep = '--'))

# Table of sc_meta_SNV CB matches
table(maf$CB %in% sc_meta_SNV$CB)
table(stats$CB %in% sc_meta_SNV$CB)
unmatched_stats <- stats[ ! stats$CB %in% sc_meta_SNV$CB, ]
dim(unmatched_stats) # [1] 253346      8
unique(unmatched_stats$project_id) # 32 PCTanno datasets
table(stats$batch)

table(stats$Tumor_Sample_Barcode %in% sc_meta_SNV$patient_id)
setdiff(stats$Tumor_Sample_Barcode, sc_meta_SNV$patient_id) # character(0) #"cirrhotic2_cd45" "cirrhotic3_cd45"
table(sc_meta_SNV$CB %in% stats$CB)
sc_meta_SNV <- sc_meta_SNV[sc_meta_SNV$CB %in% stats$CB,]
# Join sc_meta_SNV and stats data
sc_meta_SNV <- left_join(sc_meta_SNV[, -3], stats, by = 'CB')
setdiff(cohorts, unique(sc_meta_SNV$ddt_id))
# [1] "GSE172577" "GSE181254" from PCTanno which had no Epi cell in sc_allmeta

length(unique(maf$CB))
# [1] 1282466
length(unique(stats$CB))
# [1] 1119322
length(unique(sc_meta_SNV$CB))
# [1] 865976

# sc_sample_match <- sc_allsample %>% mutate(project_id = paste(project_id, patient_id))
# table(duplicated(sc_sample_match$project_id)) #35
# # View(sc_sample_match[duplicated(sc_sample_match$project_id),])
# table(stats$Tumor_Sample_Barcode %in% sc_sample_match$patient_id)
# stats <- stats %>% mutate(project_id = paste(project_id, Tumor_Sample_Barcode))
# stats <- left_join(stats, sc_sample_match, by = 'project_id')
# dim(stats) #  1135487      11
# length(unique(stats$project_id))
# mut <- stats %>% #scomatic_all %>%
#   dplyr::select(tissue, disease, project_id, Tumor_Sample_Barcode, barcode, ct, sub_type, CB, Mut_count, Mut_gene, AAChange) %>%
#   dplyr::select(tissue, disease, project_id, CB, Mut_count, Mut_gene, AAChange) %>%
# mutate(Sample = str_split(project_id, " ", simplify = T)[,2],
#        project_id = str_split(project_id, " ", simplify = T)[,1],
#        Mut_count = replace_na(Mut_count, 0))
sc_meta_SNV <- sc_meta_SNV %>%
  dplyr::select(ddt_tis, disease, ddt_id, split_ids, Tumor_Sample_Barcode, barcode, ct, sub_type, 
                CB, NMF_clusters, RMP, clusters, cluster_type, RNA_modifcations, RMP_type, Mut_count, Mut_gene, AAChange) %>%
  mutate(#Sample = Tumor_Sample_Barcode,
    Mut_count = replace_na(Mut_count, 0))
colnames(sc_meta_SNV) <- c('tissue', 'disease', 'project_id', 'subfolder', 'Sample', 'barcode', 'cell_type', 'sub_type',
                           'CB', 'NMF_clusters', 'RMP', 'clusters', 'cluster_type', 'RNA_modifcations', 'RMP_type', 'Mut_count', 'Mut_gene', 'AAChange')
table(sc_meta_SNV$Sample %in% Diseased_Samples_SNV$patient_id)
setdiff(unique(sc_meta_SNV$Sample), unique(Diseased_Samples_SNV$patient_id))
setdiff(unique(sc_meta_SNV$project_id), unique(Diseased_Samples_SNV$project_id))
sc_meta_SNV$cover_mb <- Diseased_Samples_SNV$Total_covered_mb[match(sc_meta_SNV$Sample, Diseased_Samples_SNV$patient_id)]
sc_meta_SNV$Mut_burden <- sc_meta_SNV$Mut_count/sc_meta_SNV$cover_mb
# sc_meta_SNV$Cell_Type <- str_split(sc_meta_SNV$CB, "--", simplify = TRUE)[, 3]
table(sc_meta_SNV$CB %in% maf$CB)
saveRDS(sc_meta_SNV, file = file.path('/data/rluo4/RPMfunc/Output/scRNA/sc_meta_SNV.rds')) # this version skip all bone marrow datasets

unique(maf_disco$Tumor_Sample_Barcode);unique(maf_PCTanno$Tumor_Sample_Barcode)
# 123 samples filtered
table(maf$Tumor_Sample_Barcode %in% sc_meta_SNV$Tumor_Sample_Barcode)
setdiff(unique(maf$Tumor_Sample_Barcode), sc_meta_SNV$Sample)
table(maf$CB %in% sc_meta_SNV$CB)
saveRDS(maf[maf$CB %in% sc_meta_SNV$CB, ], file = file.path('/data/rluo4/RPMfunc/Output/scRNA/sc_maf_SNV.rds')) # this version skip all bone marrow datasets
