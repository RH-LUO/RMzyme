#!/usr/local/bin/Rscript
# title: "RMzyme: Regulations of RNA Modifying Enzymes in Human"
# author: "Ruihan Luo"
# date: "July 19th,2024"
# rm(list=ls())
options(bitmapType = 'cairo')
################################################################################
# 0) visualize the SNV data of 13 organs in PCTanno
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
out_dir <-  '/data2/rluo4/RPMfunc/Output/scRNA/PCTanno_SNV'
data_path = '/data2/rluo4/All/Output'
# load('/data2/rluo4/All/Output/pheno_all.rds')#Cervix_Epi
# cell_summary <- read.csv('/data2/rluo4/All/Output/cell_summary.txt',sep = "\t",header = T)
# tissue_summary <- read.table('/data2/rluo4/RPMfunc/PCTanno_pdata/tissue_summary.txt')
# # Helper function to load and filter files
# # Load Diseased Samples
# Diseased_Samples <- read.table('/data2/rluo4/RPMfunc/PCTanno_pdata/Diseased_Samples.txt')
# # Load Normal Samples
# normal_samples <- unique(cell_summary$orig.ident[cell_summary$Tissue %in% c('Healthy','ADJ')])
# remove_samples <- setdiff(Diseased_Samples$orig.ident, tissue_summary$Replicates)# 4
# normal_samples <- c(normal_samples, remove_samples) #417
# Diseased_Samples <- Diseased_Samples[Diseased_Samples$orig.ident %in% tissue_summary$Replicates,]
# table(Diseased_Samples$orig.ident %in% cell_summary$orig.ident) #703 --> 766 samples from SNV calling
# unique(Diseased_Samples$Cohort)#32
# load(file.path(data_path, 'cell_summary_diseased.rds'))
# Diseased_Samples$Disease.Stage <- tissue_summary$Disease.Stage[
#   match(Diseased_Samples$orig.ident, tissue_summary$Replicates)]
# table(Diseased_Samples$Disease.Stage)
# ###################################
# # Summarize Coverage for samples  #
# ###################################
# indir='/data2/rluo4/RPMfunc/PCTanno_pdata/PCTanno_SNV_Coverage'
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
#   index  = tis=='Prostate' & rownames(df_cov) %in% c('P1','P2', 'P3', 'P4', 'P5', 'P6')
#   rownames(df_cov)[index] <- paste0('Dong_',  rownames(df_cov)[index])
#   df_cov$Tissue <- tis
#   total_covered_bp <- rbind(total_covered_bp, df_cov)
# }
# table(total_covered_bp$Sample %in% Diseased_Samples$Replicates)
# setdiff(total_covered_bp$Sample, Diseased_Samples$Replicates)
# setdiff( Diseased_Samples$Replicates, total_covered_bp$Sample)
# Diseased_Samples$Total_Coverage <- total_covered_bp$total_covered_mb[match(Diseased_Samples$Replicates, total_covered_bp$Sample)]
# Diseased_Samples$Total_covered_mb <- Diseased_Samples$Total_Coverage / 1e6
# colnames(Diseased_Samples)[c(5:6)] <- c('patient_id', 'project_id')
# write.table(Diseased_Samples, sep = "\t", file = '/data2/rluo4/RPMfunc/PCTanno_pdata/Diseased_Samples_PCTanno.txt', quote = F)
# Diseased_Samples <- Diseased_Samples_PCTanno <- read.csv('/data2/rluo4/RPMfunc/PCTanno_pdata/Diseased_Samples_PCTanno.txt',sep = '\t', fill = T)#Diseased_Samples
# ##################################
#   Merge all scRNA SNV Results   #
# ##################################
# Base directory and output directory
# out_dir <- '/data2/rluo4/RPMfunc/Output/scRNA/PCTanno_SNV'
# # Get unique cohorts
# combined_cohort <- unique(Diseased_Samples$Cohort)#[Diseased_Samples$Tissue == 'Prostate'])
# # Initialize an empty data frame to store merged data (Done in terminal, skip) # from Summary_PCTanno_UTH36.R
# # saveRDS(maf, file = file.path(out_dir, '../organ13_maf.rds'))
# # save(maf_results, file = file.path(out_dir, '../organ13_maf_results.RData')) # from Summary_PCTanno_UTH36.R
# # Display unique sample barcodes
# maf <- readRDS(file.path(out_dir, '../organ13_maf.rds'))
# unique(maf$Tumor_Sample_Barcode)# 624 --> 655
# # Handle NA values in AAChange.refGene
# maf$AAChange.refGene[is.na(maf$AAChange.refGene)] <- maf$Hugo_Symbol[is.na(maf$AAChange.refGene)]
# table(maf$ExonicFunc.refGene)#外显子区的SNV or InDel变异类型
# nonsynonymous SNV          stopgain          stoploss
# 4102404               162               665
# synonymous SNV           unknown
# 4810             25792
# Aggregate mutation data
# maf_PCTanno <- maf
# ###################################
# #    Summarize Mutation Burden    #
# ###################################
# stats <- maf %>%
#   filter(ExonicFunc.refGene %in% c('nonsynonymous SNV')) %>%
#   group_by(CB) %>%
#   summarise(AAChange.refGene = list(unique(AAChange.refGene))) %>%
#   mutate(Mut_count = lengths(AAChange.refGene),
#          Mut_gene = sapply(AAChange.refGene, function(x) paste(unique(str_split(x, ':', simplify = TRUE)[, 1]), collapse = ', ')),
#          AAChange = sapply(AAChange.refGene, function(x) paste(x, collapse = ', ')))
# # Extract the part between the first and second occurrence of "--"
# # stats$Tumor_Sample_Barcode <- str_extract(stats$CB, "(?<=--)[^--]*(?=--)")
# library(stringr)
# split_CB <- str_split(stats$CB, "--", simplify = TRUE)[, 2]
# index = grepl("cd45---.*?", stats$CB) #145 CB
# split_CB[index]
# # Further split the second part on the last occurrence of '_'
# split_CB[index] <- paste0(split_CB[index], '-')#str_split_fixed(stats$CB[index], '---', 3)[, 2] #str_extract(stats$CB[index], "cd45-.*?")
# # Assign the extracted part to Tumor_Sample_Barcode
# stats$Tumor_Sample_Barcode <- split_CB
# # # Apply the function to the CB column
# # stats$Tumor_Sample_Barcode <- sapply(stats$CB, extract_sample_barcode)
# # Filter and prepare Epi data
# Epi <- cell_summary %>%
#   filter(sample_name %in% combined_cohort) %>%
#   filter(!orig.ident %in% normal_samples) %>%
#   mutate(CB = paste(CB, orig.ident, sep = '--'),
#          CB = paste(CB, Cell_Type, sep = '--'))
# # Table of Epi CB matches
# table(maf$CB %in% Epi$CB)
# table(stats$CB %in% Epi$CB)
# table(stats$Tumor_Sample_Barcode %in% Epi$orig.ident)
# setdiff(stats$Tumor_Sample_Barcode, Epi$orig.ident)
# setdiff(unique(Epi$orig.ident), unique(stats$Tumor_Sample_Barcode))
# Epi <- Epi[Epi$CB %in% stats$CB,]
# # Join Epi and stats data
# Epi <- left_join(Epi, stats, by = 'CB')
# # Epi$GSE <- tissue_summary$Dataset[match(Epi$sample_name, tissue_summary$Cohort)]
# # saveRDS(Epi, file = file.path('/data2/rluo4/RPMfunc/Output/scRNA/SComatic_PCTanno.rds'))
# # system("mv SComatic_PCTanno.rds SComatic_pct.rds")
# ################################################################################
# # 1) load pdata of disco database
# # ################################################################################
# # # load('~/lrh/All/Output/organ13_Epi.RData')# from database-sum.R in ts860
# # # data_path = '/home/lorihan/lrh/All/Output'
# data_path = '/data2/rluo4/All/Output'
# setwd(data_path)
# RMP_update <- readxl::read_excel('/data2/rluo4/EpiTrans/DataCollection/RMP_update.xlsx')
# # cell_summary <- read.csv('cell_summary.txt',sep = "\t",header = T)
# # sc_majortype <- readRDS(file = '/data2/rluo4/RPMfunc/Output/scRNA/sc_majortype_1022.rds')
# # tissue_summary <- read.table('/data2/rluo4/RPMfunc/PCTanno_pdata/tissue_summary.txt')
# sc_allmeta <- readRDS(file = '/data2/rluo4/RPMfunc/Output/scRNA/sc_allmeta_1022.rds')
# RMzyme_alltissues <- sort(unique(sc_allmeta$tissue));
# RMzyme_alldatasets <- sort(unique(sc_allmeta$project_id));
# load(file = '/data2/rluo4/RPMfunc/Output/scRNA/sc_allsample_1022.RData', verbose = T)
# extract_last_element <- function(x) {
#   split_string <- strsplit(x, "_")
#   last_element <- sapply(split_string, function(y) tail(y, n = 1))
#   return(last_element)
# }
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
# Below done in UTH36 !
# ################################################################################
# # Renew on 12/07/2024:
# ################################################################################
# library(maftools)
# library(data.table)
# library(stringr)
# library(dplyr)
# options(bitmapType = 'cairo')
# sc_allmeta <- readRDS(file = '/data2/rluo4/RPMfunc/Output/scRNA/sc_allmeta_1122.rds')
# load(file = '/data2/rluo4/RPMfunc/Output/scRNA/sc_allsample_1122.RData', verbose = T)
# ##################################################
# # change 4                                       #
# ##################################################
# patient_id <- sc_allmeta$patient_id
# barcode_idx <- sc_allmeta$barcode
# idx9 <- grepl('GSE181254|PRJCA001063|GSE148673|GSE163203', sc_allmeta$project_id) 
# sc_allmeta$barcode[idx9][1:5]
# barcode_idx[idx9] <- extract_last_element(sc_allmeta$rownames[idx9])
# # View(sc_allmeta[grepl('[...]', sc_allmeta$rownames)  & idx9,])
# # table(sc_allmeta[grepl('[...]', sc_allmeta$rownames) & idx9, 'project_id'])
# barcode_ct <- str_split(barcode_idx[idx9], '[...]', simplify = T)[, 1]
# table(barcode_ct == barcode_idx[idx9])
# table(sc_allmeta$project_id[idx9])
# barcode_idx[idx9] <- barcode_ct
# barcode_idx[idx9] <- paste0(barcode_idx[idx9], '--', patient_id[idx9])
# 
# idx_dt <-  idx9
# dt <- unique(sc_allmeta[idx_dt, 'project_id'])
# print(dt)
# sc_allmeta$barcode[idx_dt] <- barcode_idx[idx_dt]
# 
# idx10 <- grepl('GSE144236', sc_allmeta$project_id) & grepl('P', sc_allmeta$patient_id)
# sc_allmeta$barcode[idx10][1:5]
# barcode_idx[idx10] <- extract_last_element(sc_allmeta$rownames[idx10])
# barcode_idx[idx10] <- paste0(barcode_idx[idx10], '--', patient_id[idx10])
# sc_allmeta$barcode[idx10] <- barcode_idx[idx10]
# sc_allmeta$disease[idx10] <- 'cutaneous squamous cell carcinoma'
# saveRDS(sc_allmeta, file = '/data2/rluo4/RPMfunc/Output/scRNA/sc_allmeta_1122.rds') # combined with diff meta from merge_nmf_meta and allmeta
# sc_allsample <- sc_allmeta %>% dplyr::select(patient_id, project_id, tissue, disease) %>%  distinct()
# save(disco, tissue_summary, sc_summary, sc_allsample, file = '/data2/rluo4/RPMfunc/Output/scRNA/sc_allsample_1122.RData')
# ################################################################################
# # Renew on 12/12/2024:
# ################################################################################
# idx11 <- grepl('GSE137829', sc_allmeta$project_id) 
# sc_allmeta$barcode[idx11][1:5]
# sc_allmeta$disease[idx11] <- 'neuroendocrine prostate cancer' # NEPC
# saveRDS(sc_allmeta, file = '/data2/rluo4/RPMfunc/Output/scRNA/sc_allmeta_1122.rds') # combined with diff meta from merge_nmf_meta and allmeta
# sc_allsample <- sc_allmeta %>% dplyr::select(patient_id, project_id, tissue, disease) %>%  distinct()
# save(disco, tissue_summary, sc_summary, sc_allsample, file = '/data2/rluo4/RPMfunc/Output/scRNA/sc_allsample_1122.RData')

##################################################
# 2) #        Statistics and Overview            #
##################################################
######################################################################
# 2.1. Summary of RNA modifications across single-cell transcriptome #
######################################################################
library(Seurat)
library(dplyr)
library(stringr)
RMP_update <- readxl::read_excel('/data2/rluo4/EpiTrans/DataCollection/RMP_update.xlsx')
sc_allmeta <- readRDS(file = '/data2/rluo4/RPMfunc/Output/scRNA/sc_allmeta_1122.rds')#sc_allmeta_0925.rds') 
load(file = '/data2/rluo4/RPMfunc/Output/scRNA/sc_allsample_1122.RData', verbose = T)
RMzyme_alltissues <- sort(unique(sc_allmeta$tissue));
RMzyme_alldatasets <- sort(unique(sc_allmeta$project_id));
table(sc_allsample$patient_id %in% sc_summary$patient_id)
diff_pt <- setdiff(sc_allsample$patient_id, sc_summary$patient_id)
# View(sc_allsample[sc_allsample$patient_id %in% diff_pt,])
sc_allmeta$omics <- 'scRNA'
sc_allsample$omics <- 'scRNA'
RMzyme_pheno <- as.data.frame(table(sc_allmeta$tissue, sc_allmeta$disease))
RMzyme_pheno <- RMzyme_pheno[RMzyme_pheno$Freq!=0,]
RMzyme_ct <- as.data.frame(table(sc_allmeta$sub_type, sc_allmeta$ct))
RMzyme_ct <- RMzyme_ct[RMzyme_ct$Freq!=0,]
######################################################################
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
write.table(tableS8, '/data2/rluo4/RPMfunc/Output/summary/NMF_clusters.txt', sep = '\t', row.names = F, quote = F)
saveRDS(NMF_clusters, file = '/data2/rluo4/RPMfunc/Output/scRNA/NMF_clusters.rds')

sc_meta_d <- sc_allmeta[! grepl('adjacent|ealthy', sc_allmeta$disease),] # 222 disease datasets
table(sc_meta_d$ddt_id == sc_meta_d$project_id)
index <- sc_meta_d$ddt_id == sc_meta_d$project_id

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
saveRDS(sc_cell_summary, file = '/data2/rluo4/RPMfunc/Output/scRNA/sc_cell_summary.rds')
# Disco cell summary done in UTH36 Rstudio Console

cell_summary <- sc_meta_d[ sc_meta_d$project_id %in% tissue_summary$Dataset, ] # 169 --> 172 disease datasets
first_elements <- extract_first_element(cell_summary$barcode)
CB = str_split(first_elements,'[-]',simplify = T)[,1]
cell_summary$CB <- CB #str_split(cell_summary$barcode, '-', simplify = T)[,1]
all_characters <- unique(unlist(strsplit(cell_summary$ct, "")))
potential_separators_and_whitespace <- all_characters[grep("[[:punct:][:space:]]", all_characters)]
# Define the replacement character
replacement <- "_"
# Replace each potential separator with the replacement character
modified_vector <- cell_summary$ct
for (sep in potential_separators_and_whitespace) {
  modified_vector <- gsub(sep, replacement, modified_vector, fixed = TRUE)
}
# Replace consecutive occurrences of the replacement character with a single underscore
cell_summary$Cell_Type <- gsub("_+", replacement, modified_vector)
saveRDS(cell_summary, file = '/data2/rluo4/RPMfunc/Output/scRNA/cell_summary.rds')
# PCTanno cell summary done in UTH36 Rstudio Console

################################################################################
# 3) organize the SNV data of RMzyme
################################################################################
# On UTH145: nohup Rscript ../Variants_disco_UTH145.R >./Variants_disco_UTH145.R.o & 
# nohup Rscript ../Variants_RMzyme.R >./Variants_RMzyme.R.o &
# nohup Rscript ../Variants_RMzyme_sort.R >./Variants_RMzyme_sort.R.o  &

# On UTH36: nohup Rscript ../Variants_PCTanno_UTH36.R > ./Variants_PCTanno_UTH36.R.o &
# nohup Rscript ../Variants_PCTanno.R >> ./Variants_PCTanno.R.o &
################################################################################
# 4) visualize the SNV data of 88 datasets 
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
base_dir <- '/data2/rluo4/RPMfunc/SCOMATIC'
out_dir <- '/data2/rluo4/RPMfunc/Output/scRNA/disco_SNV' # 88 --> 87 datasets(exclude the AML--GSE231535)
# Get unique cohorts
combined_cohort <- gsub('_SCOMATIC_annovar.RData','',list.files(out_dir))
print(combined_cohort)
options(bitmapType = 'cairo')

RMP_update <- readxl::read_excel('/data2/rluo4/EpiTrans/DataCollection/RMP_update.xlsx')
sc_cell_summary <- readRDS("/data2/rluo4/RPMfunc/Output/scRNA/sc_cell_summary.rds")
cell_summary <- readRDS("/data2/rluo4/RPMfunc/Output/scRNA/cell_summary.rds")
# sc_allmeta <- readRDS(file = '/data2/rluo4/RPMfunc/Output/scRNA/sc_allmeta_1122.rds')#sc_allmeta_0925.rds') 
load(file = '/data2/rluo4/RPMfunc/Output/scRNA/sc_allsample_1122.RData', verbose = T)
length(unique(sc_summary$project_id)) #325 --> 328 --> 332
sample.index <- c(grep("ormal", sc_summary$sample_type),grep("ealthy", sc_summary$sample_type),
                  grep("ontrol", sc_summary$sample_type))
sc_Diseased_Samples <- sc_summary[ -sample.index, c(2,3,5:8)] # 219 --> 222 --> 224 disease datasets
table(sc_Diseased_Samples$disease)
sc_Diseased_Samples <- sc_Diseased_Samples[! grepl('ealthy', sc_Diseased_Samples$disease), ]
# RMBase V3.0
SNV <- fread('/data2/rluo4/EpiTrans/DataCollection/hg38.modSNV.tar.gz', sep = '\t')
SNV$V1[1] <- 'SNV_site_1'
bed_df <- readRDS( '/data2/rluo4/EpiTrans/DataCollection/RMBase_SNV.rds')
load('/data2/rluo4/EpiTrans/RMDatasets/EpiTrans.RData', verbose = T)
NMF_clusters <- readRDS('/data2/rluo4/RPMfunc/Output/scRNA/NMF_clusters.rds')

##########################################
# 4.2 Filter and prepare SNV_cell_summary:
##########################################
# SNV_cell_summary <-  sc_cell_summary %>%
#   filter(project_id %in% combined_cohort) #%>% filter(!patient_id %in% unique(AML$patient_id))
# table(SNV_cell_summary$sub_type) # 16 sub_types (except for Germ_Cells: sc_majortype[sc_majortype$sub_type=='Germ_Cells',] all healthy)
# unique(SNV_cell_summary$project_id) #31
# unique(SNV_cell_summary$patient_id) #156
# 
# UTH36_scMapping <- list.files('/data2/rluo4/EpiTrans/Disco/scMapping/')
# GSM_scMapping <- NULL
# for (i in UTH36_scMapping) {
#   j <- list.files(file.path('/data2/rluo4/EpiTrans/Disco/scMapping/', i))
#   GSM_scMapping <<- c(GSM_scMapping, j)
# }
# table(GSM_scMapping %in% SNV_cell_summary$patient_id)
# AML <- readRDS(file=paste0('/data2/rluo4/RPMfunc/disco_pdata/AML_meta.rds'))
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
# write_xlsx(x = sc_allct, path = paste0('/data2/rluo4/RPMfunc/Output/summary/sc_allct_summary.xlsx'))
# Table S10
# sc_allct <- read_excel(paste0('/data2/rluo4/RPMfunc/Output/summary/sc_allct_summary.xlsx'))
# index <- sc_allct$Major_celltype %in% c('Bcells','Tcells')
# Mut_celltypes <- sc_allct$Minor_celltype[ ! index]#c(Epi_subtypes, Imm_subtypes, Str_subtypes, Other_subtypes)
# Mut_celltypes <- c(unique(Mut_celltypes), unique(sc_allct$Major_celltype[ ! index]))
index <- NMF_clusters$sub_type %in% c('Bcells','Tcells')
Mut_celltypes <- unique(NMF_clusters$clusters[! index])
# Aggregate mutation data
## ( On Server UTH145: nohup Rscript ../Summary_disco_UTH36.R >./Summary_disco_UTH36.R.o & since Summary_disco.R is not working neither in 36 or 57)
# On Server UTH36: nohup Rscript ../Variants_disco_UTH36.R >./Variants_disco_UTH36.R.o &
# nohup Rscript ../Variants_PCTanno_UTH36.R >./Variants_PCTanno_UTH36.R.o &
# ...
# # Loop through each cohort to load and combine data
# ...
# saveRDS(maf, file = file.path(out_dir, '../disco_maf.rds'))
# save(maf_results, file = file.path(out_dir, '../disco_maf_results.RData'))

#################################################
# # Extract the part between the first and second occurrence of "--"
# stats <- readRDS('/data2/rluo4/RPMfunc/Output/scRNA/Variants_disco.rds')
# # stats$Tumor_Sample_Barcode <- str_extract(stats$CB, "(?<=--)[^--]*(?=--)")
# library(stringr)
# split_CB <- str_split(stats$CB, "--", simplify = TRUE)[, 2]
# # index = grepl("cd45---.*?", stats$CB) #145 CB
# # split_CB[index]
# # Further split the second part on the last occurrence of '_'
# # split_CB[index] <- paste0(split_CB[index], '-')#str_split_fixed(stats$CB[index], '---', 3)[, 2] #str_extract(stats$CB[index], "cd45-.*?")
# # Assign the extracted part to Tumor_Sample_Barcode
# stats$Tumor_Sample_Barcode <- split_CB
# setdiff(unique(maf$Tumor_Sample_Barcode), unique(split_CB)) # "GSM3516662" "GSM3972020" "GSM4339781"
# stats$Cohort <- str_split(stats$CB, "--", simplify = TRUE)[, 4]
# # # Apply the function to the CB column
# for (ddt in combined_cohort) {
#   # stats$Tumor_Sample_Barcode <- sapply(stats$CB, extract_sample_barcode)
#   sc_meta_SNV <- sc_cell_summary %>%
#     # filter(project_id %in% 'GSE185381') %>%
#     filter(project_id %in% ddt) %>%
#     # filter(!orig.ident %in% normal_samples) %>%
#     mutate(CB = paste(CB, patient_id, sep = '--'),
#            CB = paste(CB, Cell_Type, sep = '--'),
#            CB = paste(CB, project_id, sep = '--'))
#   # Table of sc_meta_SNV CB matches
#   table(maf$CB %in% sc_meta_SNV$CB)
#   table(stats$CB %in% sc_meta_SNV$CB)
#   table(stats$Tumor_Sample_Barcode %in% sc_meta_SNV$patient_id)
#   setdiff(stats$Tumor_Sample_Barcode, sc_meta_SNV$orig.ident)
#   table(sc_meta_SNV$CB %in% stats$CB)
#   sc_meta_SNV <- sc_meta_SNV[sc_meta_SNV$CB %in% stats$CB,]
#   # Join sc_meta_SNV and stats data
#   stats_sub <- stats %>% filter(Cohort == ddt)
#   
#   sc_meta_SNV <- left_join(sc_meta_SNV, stats_sub, by = 'CB')
#   out_dir <- '/data2/rluo4/RPMfunc/Output/scRNA/SComatic_res'
#   saveRDS(sc_meta_SNV, file = file.path(out_dir, paste0(ddt,'_SComatic.rds'))) # this version skip all bone marrow datasets
# }

# A tab separated table of various of RNA modification sites related to SNP and SNV of Homo sapiens
# rm(list=ls())
library(maftools)
library(data.table)
library(stringr)
library(dplyr)
options(bitmapType = 'cairo')
# Required libraries
# library(future.apply)
# options(future.globals.maxSize = 2 * 1024^3)  # Adjust if absolutely necessary
library(future.apply)
library(parallel)
# Optimize future processing
n_cores <- min(detectCores() - 2, 4)
plan(multisession, workers = n_cores)
options(future.globals.maxSize = 8 * 1024^3)  # Increase memory limit

# Download the RM related SNV data from RMBase V3.0
SNP <- fread('/data2/rluo4/EpiTrans/DataCollection/hg38.modSNP.tar.gz', sep = '\t')
SNP$V1[1] <- 'SNP_site_1'
# table(SNP$V14)
SNV <- fread('/data2/rluo4/EpiTrans/DataCollection/hg38.modSNV.tar.gz', sep = '\t')
SNV$V1[1] <- 'SNV_site_1'
table(SNV$V16)
unique(SNV$V16[SNV$V16 %in% RMP_update$RNA_modifications])

# Input data
RMBase_SNV <- apply(SNV, 1, function(x) {
  info = paste0(x[c(14:16, 1)], collapse  = ' ')
  return(info)}
)
# Split the data
bed_data <- do.call(rbind, lapply(RMBase_SNV, function(x) {
  parts <- unlist(strsplit(x, " "))  # Split by space
  name <- parts[1]
  coords <- unlist(strsplit(parts[2], "[:-]"))  # Split chromosome and positions
  chrom <- coords[1]
  start <- as.numeric(coords[2]) - 1  # Convert to 0-based
  end <- as.numeric(coords[3])       # End position (exclusive)
  strand <- substr(coords[4], 1, 1)  # Extract strand
  annotation <- parts[3]
  id <- parts[4]
  c(chrom, start, end, name, 0, strand, annotation, id)  # Add score as 0
}))

# Convert to data frame
bed_df <- as.data.frame(bed_data, stringsAsFactors = FALSE)
colnames(bed_df) <- c("chrom", "chromStart", "chromEnd", "name", "score", "strand", "annotation", "modSNVid")
saveRDS(bed_df, file = '/data2/rluo4/EpiTrans/DataCollection/RMBase_SNV.rds')
# Save to a BED file (excluding the annotation column)
write.table(bed_df[, -(7:8)], "/data2/rluo4/EpiTrans/DataCollection/RMBase_SNV.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
table(bed_df$name %in% SNV$V14)

################################################################################
# Aggregate mutation data
################################################################################
# On Server UTH36:
################################################################################
# cohort_directory <- '/data2/rluo4/EpiTrans/Disco/scMapping'
# setwd(cohort_directory)
# scMapping_dir <- list.dirs(cohort_directory); print(scMapping_dir)
# Diseased_Samples <- scMapping_dir[grepl('scMapping', scMapping_dir)]
# Diseased_Samples <- Diseased_Samples[!grepl("/scMapping$", Diseased_Samples)]
# Diseased_Samples <- data.frame(SNV_path = c(Diseased_Samples) )
# extract_last_element <- function(x) {
#   split_string <- strsplit(x, "/")
#   last_element <- sapply(split_string, function(y) tail(y, n = 1))
#   return(last_element)
# }
# extract_tissue_element <- function(x) {
#   split_string <- strsplit(x, "/")[[1]]
#   return(split_string[6])
# }
# extract_cohort_element <- function(x) {
#   split_string <- strsplit(x, "/")[[1]]
#   return(split_string[7])
# }
# # Diseased_Samples$Tissue <- sapply(Diseased_Samples$SNV_path, extract_tissue_element)
# Diseased_Samples$project_id <- sapply(Diseased_Samples$SNV_path, extract_cohort_element)
# Diseased_Samples$patient_id <- extract_last_element(Diseased_Samples$SNV_path)
# table(GSM_scMapping %in% Diseased_Samples$patient_id)
# Diseased_Samples <- Diseased_Samples[ Diseased_Samples$patient_id %in% GSM_scMapping,]
# View(Diseased_Samples[ ! Diseased_Samples$patient_id %in% SNV_cell_summary$patient_id,])
# # Diseased_Samples <- Diseased_Samples[Diseased_Samples$patient_id %in% SNV_cell_summary$patient_id,]
# unique(Diseased_Samples$patient_id) # 144 --> 135
# # setdiff(Diseased_Samples$patient_id, tissue_summary$Replicates)
# # intersect(Diseased_Samples$Replicates, Chen$Samples)
# Diseased_Samples$SNV_path <- gsub('scMapping','Raw/cellranger_count', Diseased_Samples$SNV_path)
# Diseased_Samples$SNV_path <- gsub('GSE121638', 'GSE121636', Diseased_Samples$SNV_path)
# write.table(Diseased_Samples, sep = "\t", file = '/data2/rluo4/RPMfunc/disco_pdata/Diseased_Samples_UTH36.txt', quote = F)
################################################################################
# On Server UTH57:
################################################################################
library(data.table)
library(dplyr)
library(GenomicRanges)
# Paths
indir <- '/data2/rluo4/RPMfunc/Output/scRNA/disco_SNV'
out_dir <- '/data2/rluo4/RPMfunc/Output/scRNA/SComatic_res'
# Load data
large_split <- fread('/data2/rluo4/EpiTrans/DataCollection/large_split.txt')
# Convert data.table columns to a vector of dataset IDs
split_ddt_ids <- unlist(large_split, use.names = FALSE)
# Remove empty strings
split_ddt_ids <- split_ddt_ids[split_ddt_ids != ""]
# Parse in required files
SNV <- fread('/data2/rluo4/EpiTrans/DataCollection/hg38.modSNV.tar.gz', sep = '\t')
SNV$V1[1] <- 'SNV_site_1'
bed_df <- readRDS( '/data2/rluo4/EpiTrans/DataCollection/RMBase_SNV.rds')

################################################################################
# 4.3 Arrange the SNV pdata                                                    #
################################################################################
# cohort_directory <- '/data2/rluo4/EpiTrans/Disco/scMapping'
# setwd(cohort_directory)
# scMapping_dir <- list.dirs(cohort_directory); print(scMapping_dir)
# Diseased_Samples <- scMapping_dir[grepl('scMapping', scMapping_dir)]
# Diseased_Samples <- Diseased_Samples[!grepl("/scMapping$", Diseased_Samples)]
# Diseased_Samples <- data.frame(SNV_path = c(Diseased_Samples) )
# extract_last_element <- function(x) {
#   split_string <- strsplit(x, "/")
#   last_element <- sapply(split_string, function(y) tail(y, n = 1))
#   return(last_element)
# }
# extract_tissue_element <- function(x) {
#   split_string <- strsplit(x, "/")[[1]]
#   return(split_string[6])
# }
# extract_cohort_element <- function(x) {
#   split_string <- strsplit(x, "/")[[1]]
#   return(split_string[7])
# }
# # Diseased_Samples$Tissue <- sapply(Diseased_Samples$SNV_path, extract_tissue_element)
# Diseased_Samples$project_id <- sapply(Diseased_Samples$SNV_path, extract_cohort_element)
# Diseased_Samples$patient_id <- extract_last_element(Diseased_Samples$SNV_path)
# Diseased_Samples <- Diseased_Samples[ ! grepl('temp|GSE', Diseased_Samples$patient_id), ]
# # Diseased_Samples <- Diseased_Samples[ Diseased_Samples$patient_id %in% GSM_scMapping,]
# View(Diseased_Samples[ ! Diseased_Samples$patient_id %in% SNV_cell_summary$patient_id,])
# # Diseased_Samples <- Diseased_Samples[Diseased_Samples$patient_id %in% SNV_cell_summary$patient_id,]
# unique(Diseased_Samples$patient_id) # 493 -->542 -->553
# # setdiff(Diseased_Samples$patient_id, tissue_summary$Replicates)
# # intersect(Diseased_Samples$Replicates, Chen$Samples)
# Diseased_Samples$SNV_path <- gsub('scMapping','Raw/cellranger_count', Diseased_Samples$SNV_path)
# Diseased_Samples$SNV_path <- gsub('GSE162500', 'GSE162498', Diseased_Samples$SNV_path)
# table(unique(Diseased_Samples$project_id) %in% combined_cohort)
# setdiff(combined_cohort, unique(Diseased_Samples$project_id))
# UTH57_samples <- setdiff(combined_cohort,  unique(Diseased_Samples_UTH36$project_id))
# table(UTH57_samples %in% unique(Diseased_Samples$project_id))
# setdiff(UTH57_samples, unique(Diseased_Samples$project_id))
# # write.table(Diseased_Samples, sep = "\t", file = '/data2/rluo4/RPMfunc/disco_pdata/Diseased_Samples_UTH57.txt', quote = F)
# # write.table(Diseased_Samples, sep = "\t", file = '/data2/rluo4/RPMfunc/disco_pdata/Diseased_Samples_UTH57_renew1.txt', quote = F)
# write.table(Diseased_Samples, sep = "\t", file = '/data2/rluo4/RPMfunc/disco_pdata/Diseased_Samples_UTH57_renew2.txt', quote = F)
# # Diseased_Samples <- Diseased_Samples[Diseased_Samples$orig.ident %in% tissue_summary$Replicates,]
Diseased_Samples_UTH36 <- read.table('/data2/rluo4/RPMfunc/disco_pdata/Diseased_Samples_UTH36.txt') #32 = 31 + 1(GSE185381) cohorts completed in UTH36
Diseased_Samples <- read.table('/data2/rluo4/RPMfunc/disco_pdata/Diseased_Samples_UTH57_renew2.txt')
table(Diseased_Samples$patient_id %in% Diseased_Samples_UTH36$patient_id )
Diseased_Samples_disco <- rbind(Diseased_Samples_UTH36, Diseased_Samples)
table( unique(Diseased_Samples_disco$project_id) %in% combined_cohort)
Diseased_Samples_disco <- Diseased_Samples_disco[Diseased_Samples_disco$project_id %in% combined_cohort, ]

setdiff(Diseased_Samples_disco$project_id, sc_cell_summary$project_id)
table(Diseased_Samples_disco$patient_id %in% sc_cell_summary$patient_id) #703 --> 766 samples from SNV calling --> 586
table(Diseased_Samples_disco$patient_id %in% cell_summary$patient_id) #
unique(Diseased_Samples_disco$project_id)#88
unique(Diseased_Samples_disco$patient_id)#610

Diseased_Samples_PCTanno <- read.csv('/data2/rluo4/RPMfunc/PCTanno_pdata/Diseased_Samples_PCTanno.txt',sep = '\t', fill = T)
out_dir <- '/data2/rluo4/RPMfunc/Output/scRNA/PCTanno_SNV/' # 55 datasets
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
View(Diseased_Samples_SNV[ ! Diseased_Samples_SNV$patient_id %in% sc_allsample$patient_id, ])
MisMatched_Samples_SNV <- Diseased_Samples_SNV[ ! Diseased_Samples_SNV$patient_id %in% sc_allsample$patient_id, ]
MisMatched_Samples_SNV <- MisMatched_Samples_SNV[MisMatched_Samples_SNV$project_id %in% tissue_summary$Dataset,]
# MisMatched_Samples_SNV <- MisMatched_Samples_SNV[grepl('home', MisMatched_Samples_SNV$SNV_path),]
MisMatched_Samples_SNV$sample_id <- c('GSM5258385', 'GSM5258386', 'GSM5258387', 'GSM5258388', 'GSM5258389', 'GSM5258390',
                                      'P1', 'P2', 'P3', 'P4', 'P5', 'P6',
                                      'P1', 'P3', 'P5', 'P6', 'P7', 'P8', 'P9')
# write.table(MisMatched_Samples_SNV, file = '/data2/rluo4/RPMfunc/PCTanno_pdata/MisMatched_Samples_SNV.txt',row.names = F, quote = F)

################################################################################
# 4.4 Arrange the SNV results                                                  #
################################################################################
# Paths
indir <- '/data2/rluo4/RPMfunc/Output/scRNA'
out_dir <- '/data2/rluo4/RPMfunc/Output/scRNA/SComatic_res'
# Load data
large_split <- fread('/data2/rluo4/EpiTrans/DataCollection/large_split.txt')
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
SNV_dir <- '/data2/rluo4/RPMfunc/Output/scRNA/PCTanno_SNV' # 89 datasets
combined_cohort <- gsub('_SCOMATIC_annovar.RData','',list.files(SNV_dir, pattern = 'RData'))
print(combined_cohort) # 32
out_dir <- '/data2/rluo4/RPMfunc/Output/scRNA/SComatic_PCTanno'
all_batch_files <- list.files(out_dir, pattern = "*\\_SComatic_maf.rds", full.names = TRUE)
all_batch <- gsub('/data2/rluo4/RPMfunc/Output/scRNA/SComatic_PCTanno/', '', all_batch_files)
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
SNV_dir <- '/data2/rluo4/RPMfunc/Output/scRNA/disco_SNV' # 89 datasets
combined_cohort <- gsub('_SCOMATIC_annovar.RData','',list.files(SNV_dir))
print(combined_cohort) # 88 = 89 -1
# -rw-rw-r-- 1 rluo4 rluo4  262 Dec  5 13:10 GSE181279_SCOMATIC_annovar.RData
# rm GSE181279_SCOMATIC_annovar.RData
combined_cohort <- sort(combined_cohort, decreasing = T)
print(combined_cohort) # 88
out_dir <- '/data2/rluo4/RPMfunc/Output/scRNA/SComatic_res'
all_batch_files <- list.files(out_dir, pattern = "*\\_SComatic_maf.rds", full.names = TRUE)
all_batch <- gsub('/data2/rluo4/RPMfunc/Output/scRNA/SComatic_res/', '', all_batch_files)
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
# final_RM_SNV <- readRDS( file.path(indir, "final_RM_SNV.rds"))
load('/data2/rluo4/EpiTrans/RMDatasets/EpiTrans.RData', verbose = T)
cohort = 'GSE185381'
out_dir <- '/data/rluo4/RPMfunc/Output/scRNA/SComatic_res'
subdir = 'sc_adult_AML'
# outputFile <- file.path(out_dir, paste0(cohort, '_', subdir, '_SComatic_maf.rds'))
# SComatic_res <- readRDS(outputFile)
SComatic_res <- final_RM_SNV %>% filter(ddt_id == cohort) %>% filter(split_ids == subdir)
table(SComatic_res$mod_type)
table(SComatic_res$RNA_modifcations); table(is.na(SComatic_res$RNA_modifcations))

SComatic_res <- SComatic_res[!is.na(SComatic_res$RNA_modifcations),]
table(SComatic_res$mod_type %in% SComatic_res$RNA_modifcations)
table(SComatic_res$mod_type == SComatic_res$RNA_modifcations)
SComatic_res <- SComatic_res[SComatic_res$mod_type == SComatic_res$RNA_modifcations,]
table(SComatic_res$RNA_modifcations)
sort(table(SComatic_res$NMF_clusters))
table(SComatic_res$RMP)
table(SComatic_res$RMBase_geneSymbol)
summary_data <- as.data.frame(table(SComatic_res$RMP))
summary_data$RNAmodType <- SComatic_res$RNA_modifcations[match(summary_data$Var1, SComatic_res$RMP)]
colnames(summary_data)[1:2] <- c('RMP', 'Cell_count')
write_xlsx(x = summary_data, path = '/data2/rluo4/RPMfunc/Output/summary/SComatic_GSE185381_RMP_summary.xlsx')
# Table 2
summary_data <- as.data.frame(table(SComatic_res$NMF_clusters))
summary_data$RM <- SComatic_res$RNA_modifcations[match(summary_data$Var1, SComatic_res$NMF_clusters)]
summary_data$Percentage <- summary_data$Freq/sum(summary_data$Freq)*100
summary_data <- summary_data[order(summary_data$Percentage, decreasing = T),]
summary_data <- summary_data[1:20, ]
library(gcookbook)　# 为了使用数据
library(plyr)
ggplot(summary_data, aes(y= Var1, x = Percentage, fill = RM)) +
  geom_bar(stat="identity", colour="black") +
  guides(fill=guide_legend(reverse=TRUE)) +
  scale_fill_brewer(palette="Pastel1") +
  theme_minimal() +
  labs(
    title = "", #"Top 20 RNA Modification NMF clusters related to SNV",
    x = "Percentage (%)",
    y = "",
    fill = "RNA Modification"
  ) +  theme_minimal() +themes
# theme(axis.text.x = element_text(angle = 45, hjust = 1))
f = paste0('/data2/rluo4/RPMfunc/Output/scRNA/figures/SNV_chart/NMF_SNV_Percentage.png')
ggsave(height=6,width=5, filename=f, dpi = 500, device = "png")
# Fig. 5e
ggplot(summary_data, aes(y= Var1, x = Percentage, fill = RM)) +
  geom_bar(stat="identity", colour="black") +
  guides(fill=guide_legend(reverse=TRUE)) +
  scale_fill_brewer(palette="Pastel1") +
  theme_minimal() +
  labs(
    title = "", #"Top 20 RNA Modification NMF clusters related to SNV",
    x = "Percentage (%)",
    y = "",
    fill = "RNA Modification"
  ) 
colnames(SComatic_res)
paste(colnames(SComatic_res), collapse = ", ")

SComatic_maf <- final_RM_SNV %>% filter(!is.na(RNA_modifcations)) %>% filter(mod_type == RNA_modifcations) %>% 
  mutate(scSNV_id = paste0(Chromosome, ": ", Start_Position, "-", End_Position)) %>% 
  dplyr::select(scSNV_id, Hugo_Symbol, Variant_Classification, tx, exon, txChange, aaChange, Variant_Type, Func.refGene, Gene.refGene, VAF, CCF,
                tissue, disease, project_id = ddt_id , subset_id = split_ids, barcode = CB,  NMF_cluster =  NMF_clusters, RNA_modifcations, 
                RMBase_mod_id = mod_id, RMBase_mod_type = mod_type, RMBase_snv_id = snv_id, RMBase_disease, PMID, RMBase_geneSymbol)
table(SComatic_maf$RMBase_mod_id %in% SNV$V14)
SComatic_maf$RMBase_pos <- SNV$V15[match(SComatic_maf$RMBase_mod_id, SNV$V14)]
write.xlsx(SComatic_maf, file = "/data2/rluo4/RPMfunc/Output/summary/SComatic_maf_summary.xlsx", rowNames = FALSE)
# Table S15
SComatic_maf <- final_RM_SNV %>% filter(!is.na(RNA_modifcations)) %>% filter(mod_type == RNA_modifcations) %>% 
  mutate(scSNV_id = paste0(Chromosome, ": ", Start_Position, "-", End_Position)) %>% 
  dplyr::select(scSNV_id, Hugo_Symbol, Variant_Classification, tx, exon, txChange, aaChange, Variant_Type, Func.refGene, Gene.refGene, VAF, CCF,
                tissue, disease, project_id = ddt_id , subset_id = split_ids, barcode = CB,  NMF_cluster =  NMF_clusters,RMP, subtype = clusters, RNA_modifcations, 
                RMBase_mod_id = mod_id, RMBase_mod_type = mod_type, RMBase_snv_id = snv_id, RMBase_disease, PMID, RMBase_geneSymbol)
table(SComatic_maf$RMBase_mod_id %in% SNV$V14)
SComatic_maf$RMBase_pos <- SNV$V15[match(SComatic_maf$RMBase_mod_id, SNV$V14)]
saveRDS(SComatic_maf, file = "/data2/rluo4/RPMfunc/Output/summary/SComatic_maf.rds")
unique(SComatic_maf$tissue) # 41 
length(unique(SComatic_maf$scSNV_id)) # 9650
length(unique(SComatic_maf$barcode)) # 92246

# Group by scSNV_id, RMBase_mod_id, and tissue
library(dplyr)
# Count unique tissues and RMBase_mod_id for each scSNV_id
scSNV_summary <- SComatic_maf %>%
  dplyr::group_by(scSNV_id, RMBase_mod_id, RMBase_pos) %>%
  dplyr::summarize(
    unique_tissues = n_distinct(tissue),
    unique_RMBase_mod_id = n_distinct(RMBase_mod_id),
    .groups = "drop"
  )
library(dplyr)
library(tidyr)
# Step 1: Prepare the data
# Ensure data is clean and contains relevant columns
SComatic_maf_clean <- SComatic_maf %>%
  dplyr::filter(!is.na(RMBase_mod_id) & !is.na(tissue)) %>%
  dplyr::select(scSNV_id, RMBase_mod_id, RMBase_pos, tissue)
# Step 2: Tissue-specific scSNV_ids/sites
tissue_specific <- SComatic_maf_clean %>%
  dplyr::group_by(scSNV_id) %>%
  dplyr::summarize(tissue_count = n_distinct(tissue)) %>%
  dplyr::filter(tissue_count == 1) %>%
  left_join(SComatic_maf_clean, by = "scSNV_id") %>%
  dplyr::group_by(tissue) %>%
  dplyr::summarize(tissue_specific_count = n())
# Step 3: Common scSNV_ids/sites
common_scsnv <- SComatic_maf_clean %>%
  dplyr:: group_by(scSNV_id) %>%
  dplyr::summarize(tissue_count = n_distinct(tissue)) %>%
  dplyr::filter(tissue_count > 1) %>%
  left_join(SComatic_maf_clean, by = "scSNV_id") %>%
  group_by(tissue) %>%
  dplyr::summarize(common_scsnv_count = n())
# Step 4: Common RMBase_mod_id
common_rmbase <- SComatic_maf_clean %>%
  dplyr::group_by(RMBase_pos) %>%
  dplyr::summarize(tissue_count = n_distinct(tissue)) %>%
  dplyr::filter(tissue_count > 1) %>%
  left_join(SComatic_maf_clean, by = "RMBase_pos") %>%
  dplyr:: group_by(tissue) %>%
  dplyr::summarize(common_rmbase_count = n())
# Step 5: Combine and calculate percentages
total_count <- SComatic_maf_clean %>%
  dplyr::group_by(tissue) %>%
  dplyr::summarize(total_snv_count = n())

summary_data <- tissue_specific %>%
  left_join(common_scsnv, by = "tissue") %>%
  left_join(common_rmbase, by = "tissue") %>%
  left_join(total_count, by = "tissue") %>%
  dplyr:: mutate(
    Tissue_Specific = 100 * tissue_specific_count / total_snv_count,
    Common_scSNV_Sites = 100 * common_scsnv_count / total_snv_count,
    Common_RMBase_Mod_ID = 100 * common_rmbase_count / total_snv_count
  ) %>%
  dplyr::select(tissue, Tissue_Specific, Common_scSNV_Sites, Common_RMBase_Mod_ID)

# Rename for visualization
colnames(summary_data) <- c("Tissue", "Tissue_Specific", "Common_scSNV_Sites", "Common_RMBase_Mod_ID")
# Print summary data
print(summary_data)
# Required libraries
library(ggplot2)
# # Example data (replace this with your actual results)
# data <- data.frame(
#   Tissue = c("Liver", "Large Intestine", "Kidney"),
#   Tissue_Specific = c(50, 60, 70), # Percentages for tissue-specific scSNV_ids/sites
#   Common_scSNV_Sites = c(30, 20, 15), # Percentages for common scSNV_ids/sites
#   Common_RMBase_Mod_ID = c(20, 20, 15) # Percentages for common RMBase_mod_id
# )
# Reshape data for ggplot
data_long <- reshape2::melt(summary_data, id.vars = "Tissue", 
                            variable.name = "Category", 
                            value.name = "Percentage")

# Create grouped bar plot
ggplot(data_long, aes(x = Tissue, y = Percentage, fill = Category)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  scale_fill_brewer(palette = "Set2") +
  theme_minimal() +
  labs(
    title = "Tissue Specificity of SNVs and RNA Modification Sites",
    x = "Tissue",
    y = "Percentage",
    fill = "Category"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

library(gcookbook)　# 为了使用数据
library(plyr)
ggplot(data_long, aes(x = Tissue, y = Percentage, fill = Category)) +
  geom_bar(stat="identity", colour="black") +
  guides(fill=guide_legend(reverse=TRUE)) +
  scale_fill_brewer(palette="Pastel1") + 
  theme_minimal() +
  labs(
    title = "Tissue Specificity of SNVs and RNA Modification Sites",
    x = "Tissue",
    y = "Percentage",
    fill = "Category"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
################################################################################
library(dplyr)
library(tidyr)
# Step 1: Subset RM-related scSNV sites
rm_related_snvs <- SComatic_maf %>%
  dplyr::filter(!is.na(RMBase_mod_id)) %>%
  dplyr::select(scSNV_id, tissue) %>%
  distinct()

# Step 2: Create a list of scSNV_ids for each tissue
tissue_snvs <- rm_related_snvs %>%
  dplyr::group_by(tissue) %>%
  dplyr::summarize(scSNV_ids = list(unique(scSNV_id))) %>%
  deframe()

# Step 3: Compare all tissue pairs
# Create a matrix to store shared percentages
tissues <- names(tissue_snvs)
shared_matrix <- matrix(0, nrow = length(tissues), ncol = length(tissues), dimnames = list(tissues, tissues))

for (i in seq_along(tissues)) {
  for (j in seq_along(tissues)) {
    if (i != j) {
      # Get the intersection of scSNV_ids
      shared_snvs <- length(intersect(tissue_snvs[[tissues[i]]], tissue_snvs[[tissues[j]]]))
      # Calculate the percentage relative to the union of scSNV_ids
      total_snvs <- length(union(tissue_snvs[[tissues[i]]], tissue_snvs[[tissues[j]]]))
      shared_matrix[i, j] <- (shared_snvs / total_snvs) * 100
    }
  }
}

# Step 4: Convert the matrix to a data frame for visualization
shared_df <- as.data.frame(as.table(shared_matrix))
colnames(shared_df) <- c("Tissue1", "Tissue2", "Shared_Percentage")
# View the result
print(shared_df)
library(ggplot2)
themes <- theme(
  text = element_text(size = 12),
  panel.grid.major = element_line(colour = "grey90", size=0.2),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  # plot.title    = element_text(color = 'black', size   = 20, hjust = 0.5),
  plot.subtitle = element_text(color = 'black', size   = 20,hjust = 0.5),
  plot.caption  = element_text(color = 'black', size   = 20,face = 'italic', hjust = 1),
  axis.line = element_line(colour = "black"),
  axis.text.x= element_text(angle = 60, vjust = 1, hjust = 1, color = 'black', size = 12),
  # axis.text.y=element_text(angle = 0, vjust = 1, hjust = 1),
  axis.text.y   = element_text(color = 'black', size = 12, angle = 0),
  axis.title.x  = element_text(color = 'black', size = 14, angle = 0),
  axis.title.y  = element_text(color = 'black', size = 14, angle = 90),
  legend.position="none",
  legend.title  = element_text(color = 'black', size  = 20),
  legend.text   = element_text(color = 'black', size   = 20),
  axis.line.y = element_line(color = 'black', linetype = 'solid'), # y轴线特征
  axis.line.x = element_line (color = 'black',linetype = 'solid'), # x轴线特征
  # panel.border = element_rect(linetype = 'solid', size = 1.2,fill = NA) # 图四周框起来
)
rna_mod_colors <- c("#f89e81","#99a9cc","#dd9bc5","#84c7b3","#acd485","#1B9E77","#f6d573","skyblue",
                    '#A6CEE3','#1F78B4','#B2DF8A','#33A02C','#FB9A99','#E31A1C','#FDBF6F','#FF7F00','#CAB2D6','#6A3D9A','#FFFF99','#B15928' #RColorBrewer::brewer.pal(n = 12, name = 'Paired')
)
ggplot(shared_df, aes(x = Tissue1, y = Tissue2, fill = Shared_Percentage)) +
  geom_tile() + 
  scale_fill_gradient(low = "white", high = 'blue') +
  labs(
    title = "", #"Shared Percentage of RM-related SNV Sites Between Tissues",
    x = '', #"Tissue 1",
    y = '',#"Tissue 2",
    fill = "Shared (%)"
  ) +
  theme_minimal() +themes
  # theme(axis.text.x = element_text(angle = 45, hjust = 1))
f = paste0('/data2/rluo4/RPMfunc/Output/scRNA/figures/SNV_chart/SharedPercentage.png')
ggsave(height=8.5,width=8.5, filename=f, dpi = 500, device = "png")
ggsave(height=8.5,width=8.5, filename='/data2/rluo4/RPMfunc/Output/scRNA/figures/SNV_chart/SharedPercentage_legend.png', dpi = 500, device = "png")
# Fig. 5e
colnames(SComatic_maf)
table(DEM_summary$Enzyme)
table(SComatic_res$RMP)
table(DEM_summary$Enzyme %in% SComatic_res$RMP)
Enzymes <- str_split(DEM_summary$Enzyme, "_1|_2", simplify = T)[,1]
Enzymes <- gsub('FTO.2', 'FTO', Enzymes)
Enzymes <- gsub('DNMT2', 'TRDMT1', gsub('KIAA1429', 'VIRMA', 
                                        gsub('HAKAI','CBLL1', Enzymes)))
table(Enzymes %in% SComatic_res$RMP)
setdiff(unique(Enzymes), SComatic_res$RMP)
setdiff(unique(Enzymes), sc_cell_summary$RMP)
table(DEM_summary$Tissue_disease_type)

#######################################
# 4.6 Summarize Coverage for samples  #
#######################################
# # Coverage data as provided
# # coverage_data <- read_delim('/data2/rluo4/RPMfunc/PCTanno_pdata/HSIL_HPV_1.depth.chr.stat.gz')
# # coverage_data <- coverage_data[!is.na(coverage_data$MeanDepth), ]
# # coverage_data[, -1] <- apply(coverage_data[, -1], 2, function(x){
# #   y <- as.numeric(x)
# #   return(y)
# # })
# indir='/data2/rluo4/RPMfunc/disco_pdata/disco_SNV_Coverage/'
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
# write.table(Diseased_Samples_disco, sep = "\t", file = '/data2/rluo4/RPMfunc/disco_pdata/Diseased_Samples_disco.txt', quote = F)
# done in UTH36 !
outdir = '/data2/rluo4/RPMfunc/Output/scRNA'
setwd(outdir)
# TMB per cell #
Diseased_Samples_PCTanno <- read.csv('/data2/rluo4/RPMfunc/PCTanno_pdata/Diseased_Samples_PCTanno.txt',sep = '\t', fill = T)
Diseased_Samples_disco <- read.table('/data2/rluo4/RPMfunc/disco_pdata/Diseased_Samples_disco.txt') # UTH36 + UTH57 from UTH57's DiscoDatasets.R
same_cols <- intersect(colnames(Diseased_Samples_PCTanno),colnames(Diseased_Samples_disco))

Diseased_Samples_SNV <- rbind(Diseased_Samples_PCTanno[,same_cols], 
                              Diseased_Samples_disco[,same_cols])
table(unique(maf$Tumor_Sample_Barcode) %in% Diseased_Samples_SNV$patient_id)
setdiff(unique(maf$Tumor_Sample_Barcode), Diseased_Samples_SNV$patient_id)
idx <- Diseased_Samples_SNV$project_id %in% unique(MisMatched_Samples_SNV$project_id)
idx_match <- match(Diseased_Samples_SNV$patient_id[idx], MisMatched_Samples_SNV$patient_id)
Diseased_Samples_SNV$patient_id[idx] <- MisMatched_Samples_SNV$sample_id[idx_match]
setdiff(unique(maf$Tumor_Sample_Barcode), Diseased_Samples_SNV$patient_id)
# write.table(Diseased_Samples_SNV, sep = "\t", file = '/data2/rluo4/RPMfunc/Output/summary/Diseased_Samples_SNV.txt', quote = F)

#################################################
# Below was done in UTH36 by SNV_RMzyme_UTH36.R #
#################################################
# # Extract the part between the first and second occurrence of "--"
# # stats_disco <- readRDS('/data2/rluo4/RPMfunc/Output/scRNA/Variants_disco.rds') # generated by Variant_disco_UTH36.R
# # stats_PCTanno <- readRDS('/data2/rluo4/RPMfunc/Output/scRNA/Variants_PCTanno.rds') # generated by Variant_disco_UTH36.R
# # stats <- rbind(stats_PCTanno, stats_disco)
# # stats$batch <- c(rep('PCTanno', nrow(stats_PCTanno)), rep('disco', nrow(stats_disco)))
# # stats$Tumor_Sample_Barcode <- str_extract(stats$CB, "(?<=--)[^--]*(?=--)")
# stats <- maf %>%
#   # filter(ExonicFunc.refGene %in% c('nonsynonymous SNV')) %>%
#   group_by(CB) %>%
#   summarise(AAChange.refGene = list(unique(AAChange.refGene))) %>%
#   dplyr::mutate(Mut_count = lengths(AAChange.refGene),
#                 Mut_gene = sapply(AAChange.refGene, function(x) paste(unique(str_split(x, ':', simplify = TRUE)[, 1]), collapse = ', ')),
#                 AAChange = sapply(AAChange.refGene, function(x) paste(x, collapse = ', ')))
# saveRDS(stats, file = file.path(SNV_dir, '../Variants_stats.rds')) # this version skip all bone marrow datasets
# 
# library(stringr)
# split_CB <- str_split(stats$CB, "--", simplify = TRUE)[, 2]
# index = grepl("cd45---*?", stats$CB) #145 CB
# split_CB[index]
# # Further split the second part on the last occurrence of '_'
# split_CB[index] <- paste0(split_CB[index], '-')#str_split_fixed(stats$CB[index], '---', 3)[, 2] #str_extract(stats$CB[index], "cd45-.*?")
# # Assign the extracted part to Tumor_Sample_Barcode
# stats$Tumor_Sample_Barcode <- split_CB
# setdiff(unique(maf$Tumor_Sample_Barcode), unique(split_CB)) # "GSM3516662" "GSM3972020" "GSM4339781"
# stats$project_id <- str_split(stats$CB, "--", simplify = TRUE)[, 4]
# dim(stats)#1119322       8
# # # Apply the function to the CB column
# # stats$Tumor_Sample_Barcode <- sapply(stats$CB, extract_sample_barcode)
# sc_meta_SNV <- rbind(sc_cell_summary, cell_summary) %>%
#   mutate(CB = paste(CB, patient_id, sep = '--'),
#          CB = paste(CB, Cell_Type, sep = '--'),
#          CB = paste(CB, ddt_id, sep = '--'))
# 
# # Table of sc_meta_SNV CB matches
# table(maf$CB %in% sc_meta_SNV$CB)
# table(stats$CB %in% sc_meta_SNV$CB)
# unmatched_stats <- stats[ ! stats$CB %in% sc_meta_SNV$CB, ]
# dim(unmatched_stats) # [1] 253346      8
# unique(unmatched_stats$project_id) # 32 PCTanno datasets
# table(stats$batch)
# 
# table(stats$Tumor_Sample_Barcode %in% sc_meta_SNV$patient_id)
# setdiff(stats$Tumor_Sample_Barcode, sc_meta_SNV$patient_id) # character(0) #"cirrhotic2_cd45" "cirrhotic3_cd45"
# table(sc_meta_SNV$CB %in% stats$CB)
# sc_meta_SNV <- sc_meta_SNV[sc_meta_SNV$CB %in% stats$CB,]
# # Join sc_meta_SNV and stats data
# sc_meta_SNV <- left_join(sc_meta_SNV[, -3], stats, by = 'CB')
# setdiff(cohorts, unique(sc_meta_SNV$ddt_id))
# # [1] "GSE172577" "GSE181254" from PCTanno which had no Epi cell in sc_allmeta
# 
# length(unique(maf$CB))
# # [1] 1282466
# length(unique(stats$CB))
# # [1] 1119322
# length(unique(sc_meta_SNV$CB))
# # [1] 865976
# 
# # sc_sample_match <- sc_allsample %>% mutate(project_id = paste(project_id, patient_id))
# # table(duplicated(sc_sample_match$project_id)) #35
# # # View(sc_sample_match[duplicated(sc_sample_match$project_id),])
# # table(stats$Tumor_Sample_Barcode %in% sc_sample_match$patient_id)
# # stats <- stats %>% mutate(project_id = paste(project_id, Tumor_Sample_Barcode))
# # stats <- left_join(stats, sc_sample_match, by = 'project_id')
# # dim(stats) #  1135487      11
# # length(unique(stats$project_id))
# # mut <- stats %>% #scomatic_all %>%
# #   dplyr::select(tissue, disease, project_id, Tumor_Sample_Barcode, barcode, ct, sub_type, CB, Mut_count, Mut_gene, AAChange) %>%
# #   dplyr::select(tissue, disease, project_id, CB, Mut_count, Mut_gene, AAChange) %>%
# # mutate(Sample = str_split(project_id, " ", simplify = T)[,2],
# #        project_id = str_split(project_id, " ", simplify = T)[,1],
# #        Mut_count = replace_na(Mut_count, 0))
# sc_meta_SNV <- sc_meta_SNV %>% 
#   dplyr::select(ddt_tis, disease, ddt_id, split_ids, Tumor_Sample_Barcode, barcode, ct, sub_type, CB, NMF_clusters, RMP, RNA_modifcations, RMP_type, Mut_count, Mut_gene, AAChange) %>%
#   mutate(#Sample = Tumor_Sample_Barcode,
#     Mut_count = replace_na(Mut_count, 0))
# colnames(sc_meta_SNV) <- c('tissue', 'disease', 'project_id', 'subfolder', 'Sample', 'barcode', 'cell_type', 'sub_type',
#                            'CB', 'NMF_clusters', 'RMP', 'RNA_modifcations', 'RMP_type', 'Mut_count', 'Mut_gene', 'AAChange')
# table(sc_meta_SNV$Sample %in% Diseased_Samples_SNV$patient_id)
# setdiff(unique(sc_meta_SNV$Sample), unique(Diseased_Samples_SNV$patient_id))
# setdiff(unique(sc_meta_SNV$project_id), unique(Diseased_Samples_SNV$project_id))
# sc_meta_SNV$cover_mb <- Diseased_Samples_SNV$Total_covered_mb[match(sc_meta_SNV$Sample, Diseased_Samples_SNV$patient_id)]
# sc_meta_SNV$Mut_burden <- sc_meta_SNV$Mut_count/sc_meta_SNV$cover_mb
# # sc_meta_SNV$Cell_Type <- str_split(sc_meta_SNV$CB, "--", simplify = TRUE)[, 3]
# table(sc_meta_SNV$CB %in% maf$CB)
# # saveRDS(sc_meta_SNV, file = file.path('/data2/rluo4/RPMfunc/Output/scRNA/sc_meta_SNV.rds')) # this version skip all bone marrow datasets
# unique(maf_disco$Tumor_Sample_Barcode);unique(maf_PCTanno$Tumor_Sample_Barcode)
# # 123 samples filtered
# table(maf$Tumor_Sample_Barcode %in% sc_meta_SNV$Tumor_Sample_Barcode)
# setdiff(unique(maf$Tumor_Sample_Barcode), sc_meta_SNV$Sample)
# table(maf$CB %in% sc_meta_SNV$CB)
# # saveRDS(maf[maf$CB %in% sc_meta_SNV$CB, ], file = file.path('/data2/rluo4/RPMfunc/Output/scRNA/sc_maf_SNV.rds')) 

#################################################
sc_meta_SNV <- readRDS('/data2/rluo4/RPMfunc/Output/scRNA/sc_meta_SNV.rds')
# maf <- readRDS('/data2/rluo4/RPMfunc/Output/scRNA/RMzyme_maf.rds')
# sc_maf_SNV <- maf[maf$CB %in% sc_meta_SNV$CB, ]
sc_maf_SNV <- readRDS('/data2/rluo4/RPMfunc/Output/scRNA/sc_maf_SNV.rds')
table(sc_maf_SNV$ExonicFunc.refGene)#外显子区的SNV or InDel变异类型
# Data: Variants and their counts
variant_counts <- table(sc_maf_SNV$ExonicFunc.refGene)
names(variant_counts)[! grepl('synonymous',names(variant_counts))] <- 'others'
# table(sc_maf_SNV$Variant_Classification, sc_maf_SNV$ExonicFunc.refGene)
library(ggplot2)
variant_counts <- table(sc_maf_SNV$Func.refGene)#Variant_Classification)
# Convert data to a data frame
variant_df <- data.frame(
  Variant_class = names(variant_counts),
  Count = as.numeric(variant_counts)
)
# Create the pie chart
ggplot(variant_df, aes(x = "", y = Count, fill = Variant_class)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar(theta = "y") +
  labs(
    # title = "Functional Variants",
    x = NULL,
    y = NULL
  ) +
  theme_void() + # Remove axis elements
  scale_fill_brewer(palette = "Set3") # Add color palette
f = paste0('/data2/rluo4/RPMfunc/Output/scRNA/figures/SNV_chart/Func.refGene.png')
ggsave(height=4,width=4, filename=f, dpi = 500, device = "png")
# Fig. 5a
################################################################################
# Parse in test files
################################################################################
cohort = 'GSE185381'
out_dir <- '/data2/rluo4/RPMfunc/Output/scRNA/'
subdir = 'sc_adult_AML'
extract_last_element <- function(x) {
  split_string <- strsplit(x, "--")
  last_element <- sapply(split_string, function(y) tail(y, n = 1))
  return(last_element)
}
# cohorts <- extract_last_element(maf$CB)
test <-  sc_meta_SNV %>% 
  filter(project_id==cohort) %>% filter(subfolder==subdir) 
table(test$sub_type)
test_maf <-  sc_maf_SNV %>% 
  filter(CB %in% test$CB) %>%
  mutate(project_id = extract_last_element(CB)) 
unique(test_maf$Hugo_Symbol) # 490 known genes + unknown
# write.xlsx(test_maf[, -((ncol(test_maf)-1):ncol(test_maf))], file = paste0('/data2/rluo4/RPMfunc/Output/summary/test_maf_summary.xlsx'), rowNames = FALSE)
# Table S14
Mut_genes <- NULL
Mut_all <- NULL
# for (organ in organ_all) {
#   # organ = 'Cervix'
#   mut = Mut_scRNA[[organ]]
#   library(ggpubr)
#   table(mut$Tissue)
#   Mut_G <- mut[,c('cell_id','Cell_Type','Mut_count','Mut_gene','AAChange','Gene_count','Tissue')]
#   Mut_G$Organ <- organ
#   Mut_G <- Mut_G[! is.na(Mut_G$Mut_gene),]
#   Mut_G <- Mut_G[Mut_G$Mut_gene!='',]
#   Mut_all <<- rbind(Mut_all, Mut_G)
#   # Mut_G <- Mut_G[Mut_G$Cell_Type=='STM', ]
#   # Mut_G <- Mut_G[Mut_G$Tissue %in% precancer, ]
#   Mut_genes <- c(Mut_genes, paste0(Mut_G$Mut_gene))
# }
library(ggpubr)
table(test$tissue)
test$Gene_count <- apply(test,1,function(x){
  return(length(unique(strsplit(x["Mut_gene"],', ')[[1]])))
})
Mut_G <- test[,c('CB','cell_type', 'sub_type','Mut_count','Mut_gene','AAChange','Gene_count','tissue', 'disease')]
# Mut_G$cell_type <- test$sub_type
Mut_G <- Mut_G[! is.na(Mut_G$Mut_gene),]
Mut_G <- Mut_G[Mut_G$Mut_gene!='',]
Mut_all <<- rbind(Mut_all, Mut_G)
Mut_all$Index = 1:nrow(Mut_all)
Mut_all <- Mut_all[, c('Index', colnames(Mut_all)[1:(ncol(Mut_all)-1)])]
Mut_genes <- c(Mut_genes, paste0(Mut_G$Mut_gene))
n <- Mut_all$Mut_gene
m <- strsplit( Mut_all$Mut_gene,', ')
g <- lapply(m, function(x){
  index = unique(x %in% c('NPM1')) #
  # print(x)
  return(index)
})
# # View(Mut_all[g!='FALSE',])
table(Mut_all[g!='FALSE', 'tissue'])
table(Mut_all[g!='FALSE', 'cell_type'])
NPM1_mut <- test_maf[grepl('NPM1', test_maf$AAChange),]
NPM1 <- Mut_all[g!='FALSE',]
sort(table(NPM1$cell_type))
summary(as.numeric(NPM1_mut$FATHMM_score))
summary(as.numeric(NPM1_mut$DANN_score))
summary(as.numeric(test_maf$FATHMM_score))
# View(final_RM_SNV[final_RM_SNV$CB %in% NPM1$CB,])

ggplotdata <- as.data.frame(sort(table(NPM1$cell_type), decreasing = T))
colnames(ggplotdata)[1] <- 'cell_type'
ggplotdata <- ggplotdata[1:10, ]
bp <- ggbarplot(
  ggplotdata, x = "cell_type", y = 'Freq', fill = 'cell_type',  palette = "jco",#"palettetown", #View(paletteer::palettes_d_names)
  #palette = "quilava",
  # add = c('mean_se'), #facet.by = "TissueType"
) + xlab("") + ylab("") +
  # stat_pvalue_manual(stat.test,  label = "p.adj.signif", tip.length = 0.01, hide.ns = TRUE) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) + theme(
    text = element_text(size = 12),
    panel.grid.major = element_line(colour = "grey90", size=0.2),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    # axis.text.x = element_blank(),
    axis.text.x=element_text(angle = 60, vjust = 1, hjust = 1),
    axis.text.y=element_text(angle = 0, vjust = 1, hjust = 1),
    legend.position="none",
    legend.title=element_text(size=12))
bp
f = paste0('/data2/rluo4/RPMfunc/Output/scRNA/figures/SNV_chart/NPM1mut_Top10_celltypes.png')
ggsave(plot = bp, file = f, height =  3, width = 4, dpi = 500)

Mut_G <- paste0(unlist(strsplit(Mut_genes,', ')))
unique(Mut_G)#490 + UNKNOWN + Unknown = 492
getwd()
setwd('/data2/rluo4/RPMfunc/Output/scRNA/figures/SNV_chart')
Mut_CNV <- read.csv('/data2/rluo4/All/Mut_CNV_Tirosh.txt',sep = '\t')
# Mut_CNV <- Mut_CNV[!duplicated(Mut_CNV$Event),]
Mutants_TCGA <- unique(Mut_CNV$Event)
Mutants_TCGA <- unique(str_split(Mutants_TCGA, ' ',simplify = T)[,1])
Mutants_TCGA <- Mutants_TCGA[!grepl("^[0-9]+", Mutants_TCGA)]
Mutants_TCGA <- Mutants_TCGA[Mutants_TCGA %ni% c('Xp','Xq')]
table(Mutants_TCGA %in% unique(Mut_G))
Mutants_TCGA[Mutants_TCGA %in%  unique(Mut_G)]#28 for STM; 30 for all celltypes
# [1] "CASP8" "NPM1" 
# grep('IDH1',Mut_G)
# table(Mut_G=='TP53')#
# table(Mut_G=='RB1')#
# table(Mut_G=='TP53I3')#
# table(Mut_G=='CTNNB1')##SMAD4
# table(Mut_G=='RUNX1')##SMAD4
# 
# table(Mut_G=='DNMT3A')#IDH1
# table(Mut_G=='ASXL1')
# table(Mut_G=='KRAS')
# table(Mut_G=='NPM1')# high freq
# table(Mut_G=='CTCF')
Mut_df <- as.data.frame(table(Mut_G))
# View(Mut_df[Mut_df$Mut_G %in% Mutants_TCGA,])
# NCOR1 5066
# NPM1  4100
# CDKN2A 2705
high_mut <- Mut_df[Mut_df$Mut_G %ni% c('Unknown', 'UNKNOWN'),]#Mut_df[Mut_df$Mut_G %in% Mutants_TCGA,]
unique(high_mut$Mut_G)
high_mut <- high_mut$Mut_G[order(high_mut$Freq,decreasing = T)]
high_mut <- head(high_mut, 20)
g <- lapply(m, function(x){
  index = unique(x %in% high_mut)
  # print(x)
  return(index)
})
# View(Mut_all[g!='FALSE',])
table(Mut_all[g!='FALSE', 'tissue'])
table(Mut_all[g!='FALSE', 'cell_type'])
print(high_mut)
# high_mut_df <- as.data.frame(table(Mut_all[g!='FALSE', 'Cell_Type'], Mut_all[g!='FALSE', 'Tissue']))
top100 <- table(Mut_all[g!='FALSE', 'cell_type'])
top100 <- names(top100[top100>=100])
high_mut_df <- as.data.frame(table(Mut_all[g!='FALSE', 'cell_type'], Mut_all[g!='FALSE', 'tissue']))
colnames(high_mut_df)[1:2] <- c('Cell_Type', 'TissueType')
high_mut_df$Cell_Type <- as.character(high_mut_df$Cell_Type)
high_mut_df <- high_mut_df[high_mut_df$Cell_Type %in% top100,]
high_mut_df <- high_mut_df[high_mut_df$Freq>=1,]
high_mut_df <- high_mut_df[order(high_mut_df$Freq,decreasing = F),]
str(high_mut_df)
table(high_mut_df$TissueType)
palette=c("#f89e81","#99a9cc","#dd9bc5","#84c7b3","#acd485","#1B9E77","#f6d573","#66A61E","skyblue")
# color_cluster=c(RColorBrewer::brewer.pal(n = 9, name = 'Paired'))#,color_cluster)
library(ggpubr)
ggdotchart(high_mut_df, x = 'Cell_Type', y = 'Freq', color = 'TissueType', palette = c("#f89e81","#99a9cc"),
           sorting = 'asc', sort.by.groups = TRUE, add = 'segments', add.params = list(color = 'lightgray', size = 2),
           group = 'TissueType', dit.size = 6, ggtheme = theme_pubclean()) + font('x.text', size = 8, vjust = 0.5)
ggdotchart(high_mut_df, x = 'Cell_Type', y = 'Freq',
           color = "TissueType",                                # Color by groups
           palette = palette,#c("#00AFBB", "#E7B800", "#FC4E07"), # Custom color palette
           # sorting = "descending",                       # Sort value in descending order
           add = "segments",                             # Add segments from y = 0 to dots
           add.params = list(color = "lightgray", size = 3), # Change segment color and size
           group = "TissueType",                                # Order by groups
           dot.size = 6,                                 # Large dot size
           label = round(high_mut_df$Freq,1),                        # Add mpg values as dot labels
           font.label = list(color = "brown", size = 10,
                             vjust = 0.5),               # Adjust label parameters
           ggtheme = theme_pubr()                        # ggplot2 theme
) + geom_hline(yintercept = 0, linetype = 2, color = "lightgray") +  xlab('') + ylab("") + #ylab('Frequency of cells with top20 mutants') + 
  theme(
    plot.title    = element_text(color = 'black', size   = 16, hjust = 0.5),
    plot.subtitle = element_text(color = 'black', size   = 16,hjust = 0.5),
    plot.caption  = element_text(color = 'black', size   = 16,face = 'italic', hjust = 1),
    axis.title.x  = element_text(color = 'black', size = 16, angle = 0),
    axis.text.y   = element_text(color = 'black', size = 11, angle = 0),
    axis.title.y  = element_text(color = 'black', size = 12, angle = 90),
    # legend.title  = element_text(color = 'black', size  = 20),
    legend.text   = element_text(color = 'black', size   = 14),
    legend.title = element_blank()
    # axis.line.y = element_line(color = 'black', linetype = 'solid'), # y轴线特征
    # axis.line.x = element_line (color = 'black',linetype = 'solid'), # x轴线特征
    # panel.border = element_rect(linetype = 'solid', size = 1.2,fill = NA) # 图四周框起来
  ) 
# install.packages(c("Cairo", "ragg", "svglite"), lib = "/usr/local/lib/R/site-library/", dependencies = TRUE, repos = "http://cran.rstudio.com/")
f = paste0('/data2/rluo4/RPMfunc/Output/scRNA/figures/SNV_chart/Top20_mutated_celltype.png')
ggsave(height=4.5,width=6.5, filename=f, dpi = 500, device = "png")

# SComatic_maf <- SComatic_res[!is.na(SComatic_res$RNA_modifcations),]
# table(SComatic_maf$mod_type %in% SComatic_maf$RNA_modifcations)
# table(SComatic_maf$mod_type == SComatic_maf$RNA_modifcations)
# SComatic_maf <- SComatic_maf[SComatic_maf$mod_type == SComatic_maf$RNA_modifcations,]
# table(SComatic_maf$RNA_modifcations)
# # Apply the function to the CB column
# for (ddt in combined_cohort) {
#   # stats$Tumor_Sample_Barcode <- sapply(stats$CB, extract_sample_barcode)
# 
#   stats_sub <- sc_meta_SNV %>% filter(Cohort == ddt)
#   
#   out_dir <- '/data2/rluo4/RPMfunc/Output/scRNA/SComatic_RMzyme'
#   saveRDS(stats_sub, file = file.path(out_dir, paste0(ddt,'_SComatic.rds'))) # this version skip all bone marrow datasets
# }
# Table of mutations by tissue and cell type
# table(mut$tissue, mut$ct)
# table(maf$Tumor_Sample_Barcode %in% scomatic_disco$patient_id)
# # Method 1 #
# TMB_stats <- mut %>%
#   group_by(Sample) %>%
#   summarise(Mut_count = sum(Mut_count))
# TMB_stats$cover_mb <- Diseased_Samples_SNV$Total_covered_mb[match(TMB_stats$Sample, Diseased_Samples_SNV$patient_id)]
# TMB_stats$Mut_burden <- TMB_stats$Mut_count/TMB_stats$cover_mb
# # print(TMB_stats)
# print(unique(mut$Tumor_Sample_Barcode)) #819
# print(unique(mut$project_id)) #87
# library(stringr)
# split_CB <- str_split(test.tmb$Tumor_Sample_Barcode, "--", simplify = TRUE)[, 2]
# index = grepl("cd45---.*?", test.tmb$Tumor_Sample_Barcode) #145 CB
# split_CB[index] <- paste0(split_CB[index], '-')#
# # Assign the extracted part to Tumor_Sample_Barcode
# test.tmb$Sample <- split_CB

# Method 2 #
test.laml <- test_maf
test.laml$Sample <- test.laml$Tumor_Sample_Barcode
test.laml$Tumor_Sample_Barcode <- test.laml$CB # set identifiers as cell barcodes
table(test.laml$ExonicFunc.refGene)
# TMB per sample #
test.laml <- read.maf(maf = test.laml)
test.laml@data$Chromosome <- gsub('chr', '', test.laml@data$Chromosome)
# saveRDS(test.laml, file = '/data2/rluo4/RPMfunc/Output/scRNA/RMzyme_SNV_maf.rds')
# test.laml <- readRDS('/data2/rluo4/RPMfunc/Output/scRNA/RMzyme_SNV_maf.rds')
test.tmb <- tmb(test.laml) # the same with: TMB=as.data.frame(table(test.laml@data$Tumor_Sample_Barcode))
dim(test.laml@variants.per.sample)
table(test.tmb$total==0)
TMB=as.data.frame(table(test_maf$CB))
table(TMB$Freq==0)
colnames(TMB)[2] <- 'total'
test.tmb <- TMB
library(stringr)
# split_CB <- str_split(test.tmb$Tumor_Sample_Barcode, "--", simplify = TRUE)[, 2]
# index = grepl("cd45---.*?", test.tmb$Tumor_Sample_Barcode) #145 CB
# split_CB[index] <- paste0(split_CB[index], '-')#
# test.tmb$Sample <- split_CB
colnames(test.tmb)[1] <- 'CB'
test.tmb <- left_join(test, as.data.frame(test.tmb), by = 'CB')
# test.tmb <- test.tmb[!test.tmb$Dataset %in% exclusive_mastcell_basophil_projects$project_id,]
identifier = test.tmb$clusters
test.tmb_stats <- test.tmb %>% mutate(identifier = clusters) %>%
  dplyr::group_by(identifier) %>% dplyr::summarise(Mut_count = sum(total)) 
index = match(test.tmb_stats$identifier,identifier)
# Diseased_Samples$Total_Coverage <- total_covered_bp$total_covered_mb[match(Diseased_Samples$Replicates, total_covered_bp$Sample)]
# Diseased_Samples$Total_covered_mb <- Diseased_Samples$Total_Coverage / 1e6
test.tmb_stats$cover_mb <- test.tmb$cover_mb[index]
test.tmb_stats$Mut_burden <- test.tmb_stats$Mut_count/test.tmb_stats$cover_mb
test.tmb_stats <- test.tmb_stats[order(test.tmb_stats$Mut_burden, decreasing = F),]
# Load ggplot2 package
library(ggplot2)
# Convert Var1 to a factor to maintain the order
ylab<- "No. of mutations (/Mb)"#expression(paste("FCGR3A  ", "log"["2"], "(SCNA)"))
title = paste0(tis, '-related diseases')
# Create a bar plot using ggplot2
if(length(test.tmb_stats$Sample)>=3){
  axis.text.samples = element_blank()
} else{
  axis.text.samples = element_text(angle = 60, vjust = 1, hjust = 1, color = 'black', size = 12)
}
themes <- theme(
  text = element_text(size = 12),
  panel.grid.major = element_line(colour = "grey90", size=0.2),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  # plot.title    = element_text(color = 'black', size   = 20, hjust = 0.5),
  plot.subtitle = element_text(color = 'black', size   = 20,hjust = 0.5),
  plot.caption  = element_text(color = 'black', size   = 20,face = 'italic', hjust = 1),
  axis.line = element_line(colour = "black"),
  axis.text.x= axis.text.samples,
  # axis.text.y=element_text(angle = 0, vjust = 1, hjust = 1),
  axis.text.y   = element_text(color = 'black', size = 12, angle = 0),
  axis.title.x  = element_text(color = 'black', size = 14, angle = 0),
  axis.title.y  = element_text(color = 'black', size = 14, angle = 90),
  legend.position="none",
  legend.title  = element_text(color = 'black', size  = 20),
  legend.text   = element_text(color = 'black', size   = 20),
  axis.line.y = element_line(color = 'black', linetype = 'solid'), # y轴线特征
  axis.line.x = element_line (color = 'black',linetype = 'solid'), # x轴线特征
  panel.border = element_rect(linetype = 'solid', size = 1.2,fill = NA) # 图四周框起来
)
# # Ensure 'Sample' is a factor with levels ordered by 'Tissue'
# test.tmb_stats$Sample <- factor(test.tmb_stats$Sample, levels = test.tmb_stats$Sample[order(test.tmb_stats$Tissue)])
# # Sort the dataframe by 'Tissue' to maintain order in the plot
# test.tmb_stats <- test.tmb_stats[order(test.tmb_stats$Tissue, decreasing = F),]

ggplot(test.tmb_stats, aes(x = identifier, y = Mut_burden)) +
  geom_bar(stat = "identity", fill = "steelblue") +    # Create bars
  labs(title = "Mutational Burden Across Samples",     # Add titles and labels
       x = "Sample",
       y = ylab) + themes
# Create the faceted bar plot
p = ggplot(test.tmb_stats, aes(x = identifier, y = Mut_burden)) +
  geom_bar(stat = "identity", aes(fill = Tissue)) +  # Use 'fill' to color bars by 'Tissue'
  # geom_text(aes(label = round(Mut_burden, 2)),       # Add text labels on bars
  #           vjust = -0.3, size = 3.5) +
  labs(# title = title,#"Mutational Load Across Samples",   # Add titles and labels
    x = ifelse(length(test.tmb_stats$identifier)>=3, "Sample", ''),
    y = ylab) +
  facet_wrap(~ Tissue, scales = "free_x") + themes          # Facet by 'Tissue' with independent x-axis
f = paste0(file.path('/data/rluo4/RPMfunc/Output/scRNA/figures/SNV_chart/'), 'RMzyme_TMB.png')
sample_tis <- unique(test.tmb_stats$Tissue)
wid = ifelse(length(sample_tis)>=2, length(sample_tis)*1.5, length(sample_tis)*2)
if(length(sample_tis)==1){wid = length(sample_tis)*2 }
ggsave(plot=p,height=wid*0.75,width=wid*0.9, filename=f, dpi = 500, device = "png")
####################################
# 8.4 Clonality patterns analysis  #
####################################
# Load necessary libraries
library(dplyr)
library(ggplot2)
library(stringr)
library(ggpubr)
# Data Preparation
clonality <- test_maf %>% #test.laml@data %>%
  # mutate(identifier = Cell_type_observed) %>%
  dplyr::select(Tumor_Sample_Barcode, Chromosome, Start_Position, End_Position, Hugo_Symbol,
         aaChange, VAF, CCF, CB, Cell_type_observed, FATHMM_score)
clonality <- left_join(clonality, test, by = 'CB')
# clonality <- clonality[!clonality$Dataset %in% exclusive_mastcell_basophil_projects$project_id,]
# Clean and transform VAF and CCF columns
clonality <- clonality %>%
  mutate(VAF = as.numeric(str_split_fixed(VAF, '[,|]', 3)[, 1]),
         CCF = as.numeric(str_split_fixed(CCF, '[,|]', 3)[, 1]))
colnames(clonality)
# # Match Tissue information
# Handle missing values by removing rows with NA in VAF or CCF # Classify mutations as clonal or subclonal # Check the distribution of cell types and clonality
# Filter and prepare the data
clonality <- clonality %>%
  mutate(identifier = clusters) %>%
  filter( !is.na(identifier) & !is.na(VAF) & !is.na(CCF)) %>%
  mutate(clonality = ifelse(CCF >= 0.9, 'clonal', 'subclonal'))
colnames(clonality)
table(clonality$clusters == clonality$cell_type)
# Method1:
clonality_summary <- clonality %>%
  dplyr::group_by(identifier) %>%
  # dplyr::group_by(Tumor_Sample_Barcode, identifier) %>%
  dplyr::summarise(
    total_mutations = n(),
    clonal_mutations = sum(CCF >= 0.9, na.rm = TRUE),
    subclonal_mutations = total_mutations - clonal_mutations,
    clonal_proportion = clonal_mutations / total_mutations,
    subclonal_proportion = subclonal_mutations / total_mutations
  )
# # Merge the coverage data into the clonality data frame
clonality_summary <- test.tmb_stats %>%
  left_join(clonality_summary, by = c("identifier"))
clonality_summary <- clonality_summary %>%
  mutate(
    clonal_density = clonal_mutations / cover_mb,   # Clonal mutations per base pair
    subclonal_density = subclonal_mutations / cover_mb # Subclonal mutations per base pair
  )
clonality_summary <- clonality_summary[order(clonality_summary$clonal_proportion, decreasing = T),]
# clonality_summary <- clonality_summary[order(clonality_summary$clonal_density, decreasing = T),]
# # # Ensure 'Sample' is a factor with levels ordered by 'Tissue'
# clonality_summary$Sample <- factor(clonality_summary$Sample, levels = clonality_summary$Sample[order(clonality_summary$Tissue)])
# # Sort the dataframe by 'Tissue' to maintain order in the plot
# clonality_summary <- clonality_summary[order(clonality_summary$Tissue, decreasing = F),]
# print(table(clonality_summary$Tissue))
# print(summary(clonality_summary$clonal_density))
# print(summary(clonality_summary$clonal_proportion[clonality_summary$Stage == "Cancer"]))
# print(summary(clonality_summary$clonal_proportion[clonality_summary$Stage == "Precancer"]))
# clonality$Tissue <- clonality_summary$Tissue[match(clonality$Sample, clonality_summary$Sample)]
# clonality$Stage <- clonality_summary$Stage[match(clonality$Sample, clonality_summary$Sample)]

# clonality_summary_all <<- rbind(clonality_summary_all, clonality_summary)
library(reshape2)
library(ggplot2)
library(ggridges)
library(ArchR)
themes <- theme(
  plot.title    = element_text(color = 'black', size   = 10, hjust = 0.5),
  plot.subtitle = element_text(color = 'black', size   = 10,hjust = 0.5),
  plot.caption  = element_text(color = 'black', size   = 10,face = 'italic', hjust = 1),
  axis.text.x   = element_text(color = 'black', size = 12, angle = 0),
  axis.text.y   = element_text(color = 'black', size = 12, angle = 0),
  axis.title.x  = element_text(color = 'black', size = 16, angle = 0),
  axis.title.y  = element_text(color = 'black', size = 16, angle = 90),
  legend.title  = element_text(color = 'black', size  = 16),
  legend.text   = element_text(color = 'black', size   = 16),
  legend.position = 'none',
  axis.line.y = element_line(color = 'black', linetype = 'solid'), # y轴线特征
  axis.line.x = element_line (color = 'black',linetype = 'solid'), # x轴线特征
  # panel.border = element_rect(linetype = 'solid', size = 1.2,fill = NA) # 图四周框起来
)
set.seed(123)
unique(clonality_summary$identifier)
ggbarplotdata <- clonality_summary[clonality_summary$identifier %in% Mut_celltypes, ]
ggbarplotdata <- ggbarplotdata[ ! ggbarplotdata$identifier %in% 'Blood_Cells', ]
str(ggbarplotdata)
p1 <- ggbarplot(ggbarplotdata, y = 'clonal_proportion', x = 'identifier', fill = 'identifier',  palette = "palettetown", #View(paletteer::palettes_d_names)
  #palette = "quilava",
  # add = c('mean_se'), #facet.by = "TissueType"
) + xlab("") + ylab("") +
  # stat_pvalue_manual(stat.test,  label = "p.adj.signif", tip.length = 0.01, hide.ns = TRUE) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  theme(
    text = element_text(size = 12),
    panel.grid.major = element_line(colour = "grey90", size=0.2),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    # axis.text.x = element_blank(),
    axis.text.x=element_text(angle = 60, vjust = 1, hjust = 1),
    axis.text.y=element_text(angle = 0, vjust = 1, hjust = 1),
    legend.position="none",
    legend.title=element_text(size=12))
p1
f = paste0('/data2/rluo4/RPMfunc/Output/scRNA/figures/SNV_chart/ClonalProportion.png')
ggsave(plot=p1, height=4.5, width=5, filename=f, dpi = 500, device = "png")
p2 <- ggbarplot(ggbarplotdata, y = 'subclonal_proportion', x = 'identifier', fill = 'identifier',  palette = "bpalette", #View(paletteer::palettes_d_names)
                #palette = "quilava",
                # add = c('mean_se'), 
) + xlab("") + ylab("") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  theme(
    text = element_text(size = 12),
    panel.grid.major = element_line(colour = "grey90", size=0.2),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    # axis.text.x = element_blank(),
    axis.text.x=element_text(angle = 60, vjust = 1, hjust = 1),
    axis.text.y=element_text(angle = 0, vjust = 1, hjust = 1),
    legend.position="none",
    legend.title=element_text(size=12))
p2
f = paste0('/data2/rluo4/RPMfunc/Output/scRNA/figures/SNV_chart/SubclonalProportion.png')
ggsave(plot=p2, height=4.5, width=5, filename=f, dpi = 500, device = "png")

