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
library(tidyverse) 
library(GEOquery)
library(readxl)
options(bitmapType = 'cairo')
Sys.setenv("VROOM_CONNECTION_SIZE" = 999999999)
getOption('timeout')
options(timeout=100000000)
# path <- c('/home/rluo4/R/x86_64-pc-linux-gnu-library/4.3', '/opt/R/4.3.1/lib64/R/library', '/home/rluo4/R/x86_64-conda-linux-gnu-library/4.3', '/data2/rluo4/bin/miniconda3/lib/R/library')
# .libPaths(path)
library(future.apply)
options(future.globals.maxSize = 2 * 1024^3)  # Adjust if absolutely necessary

##################################################
# 7.2. Organization of merge nmf metadata
##################################################
large_split <- fread('/data/rluo4/EpiTrans/DataCollection/large_split.txt')
# Convert data.table columns to a vector of dataset IDs
split_ddt_ids <- unlist(large_split, use.names = FALSE)
# Remove empty strings
split_ddt_ids <- split_ddt_ids[split_ddt_ids != ""]

################################################################################
# 1) load pdata of PCTanno database
################################################################################
# load('~/lrh/All/Output/organ13_Epi.RData')# from database-sum.R in ts860
# data_path = '/home/lorihan/lrh/All/Output'
data_path = '/data/rluo4/All/Output'
setwd(data_path)
RMP_update <- readxl::read_excel('/data/rluo4/EpiTrans/DataCollection/RMP_update.xlsx')
sc_majortype <- readRDS(file = '/data/rluo4/RPMfunc/Output/scRNA/sc_majortype_1122.rds')
load(file = '/data/rluo4/RPMfunc/Output/scRNA/sc_allsample_1122.RData', verbose = T)
RMzyme_alltissues <- sort(unique(sc_summary$tissue));
RMzyme_alldatasets <- sort(unique(sc_summary$project_id));
extract_last_element <- function(x) {
  split_string <- strsplit(x, "_")
  last_element <- sapply(split_string, function(y) tail(y, n = 1))
  return(last_element)
}

library(Seurat)
library(dplyr)
library(stringr)

##################################################
# (2) merging allmeta by myself:
##################################################
allmeta_path <- "/data/rluo4/RPMfunc/Output/scRNA/allmeta"
allmeta_files <- list.files(path = allmeta_path, pattern = "allmeta.rds$", 
                            full.names = TRUE, recursive = T)
print(allmeta_files)
# View(RMzyme_all_meta[RMzyme_all_meta$project_id %in% na_barcode_projectid,])
disease_tissue_ids <- str_split(allmeta_files,'/', simplify = T)[,8]
table(unique(disease_tissue_ids) %in% sc_summary$tissue)
setdiff(RMzyme_alltissues, unique(disease_tissue_ids))

disease_project_ids <- str_split(allmeta_files,'/', simplify = T)[,9] # (255 + 1 GSE110686_allmeta.rds)
# missed_bc_files <- allmeta_files[disease_project_ids %in% paste0(na_barcode_projectid,'_allmeta.rds')]
# allmzeta_files <- allmeta_files[grepl('GSE148673|GSE163203|GSE181254|PRJCA001063',allmeta_files)]
# Function to extract patient_id
extract_patient_id <- function(barcodes) {
  patient_ids <- sapply(barcodes, function(bc) {
    if (str_detect(bc, "--")) {
      # Case 1: Extract after '--' and before any underscores
      str_extract(bc, "(?<=--)[^_]+")
    } else if (str_detect(bc, "^[^_]+_[^_]+_[^_]+_")) {
      # Case 2: If three underscores exist, capture the first three segments
      str_extract(bc, "^[^_]+_[^_]+_[^_]+")
    } else if (str_detect(bc, "^[^_]+_[^_]+_")) {
      # Case 3: If two underscores exist, capture the first two segments
      str_extract(bc, "^[^_]+_[^_]+")
    } else if (str_detect(bc, "_")) {
      # Case 4: Extract before the first underscore
      str_extract(bc, "^[^_]+")
    } else {
      # Case 5: Fallback to simple barcode
      bc
    }
  })
  return(patient_ids)
}
sc_allmeta<- NULL 
for (i in 1:length(allmeta_files)){
  # i=2
  print(allmeta_files[i])
  
  dataset <- gsub('_allmeta.rds', '', str_split(allmeta_files[i],'/', simplify = T)[,9])
  print(paste('disease dataset is: ', dataset))
  idx_dt <- sc_summary$project_id==dataset
  
  tissue <- gsub('_allmeta.rds', '', str_split(allmeta_files[i],'/', simplify = T)[,8])
  print(paste('disease tissue is: ', tissue))
  
  tissue_dt <- unique(sc_summary[idx_dt, 'tissue'])
  print(paste('analyzed tissues are: ', tissue_dt))
  idx_tis <- sc_summary$tissue %in% tissue_dt
  
  meta <- readRDS(allmeta_files[i])
  dim(meta)
  print('orig.idents are below: ')
  print(table(meta$orig.ident))
  meta$rownames <- rownames(meta)#str_split(rownames(meta), '[...]', simplify = T)[, 1]
  meta$tissue <- tissue
  
  table(meta$patient_id == meta$orig.ident)
  table(is.na(unique(meta$patient_id)))
  # View(meta[meta$patient_id != meta$orig.ident, ])
  # View(meta[is.na(meta$patient_id),])
  table(is.na(meta$tissue))
  print(table(meta$tissue))
  print(unique(meta$patient_id))
  
  if( NA %in% unique(meta$patient_id) | (! 'patient_id' %in% colnames(meta)) ){
    idx_patient1 = ( ! 'patient_id' %in% colnames(meta) )
    if( ! idx_patient1 ){
      meta$patient_id[is.na(meta$patient_id)] <- meta$orig.ident[is.na(meta$patient_id)]
    } else{
      meta$patient_id <- meta$orig.ident
    }
    idx_patient1 = idx_patient1 | ( ! 'project_id' %in% colnames(meta) ) | ( ! dataset %in% unique(meta$project_id) )
    idx_patient2 = length(unique(grep('Seurat', unique(meta$patient_id)))) != 0 #== 'SeuratObject|SeuratProject'
    idx_patient3 = NA %in% unique(meta$patient_id) 
    
    if( idx_patient1 | idx_patient2 | idx_patient3 ){
      # Regular expressions for different formats
      barcodes <- meta$rownames
      # Updated regex approach for all cases
      patient_ids <- extract_patient_id(barcodes)
      print(unique(patient_ids))
      idx_patient <- meta$rownames!=patient_ids
      meta$patient_id[idx_patient] <- str_split(patient_ids[idx_patient], '[...]', simplify = T)[, 1]
      
    }
  }
  table(unique(meta$patient_id) %in% sc_summary$patient_id[idx_dt])
  unique(meta$tissue[ ! meta$patient_id %in% sc_summary$patient_id[idx_dt] ])
  setdiff(unique(meta$patient_id), sc_summary$patient_id[idx_dt])
  setdiff(unique(meta$orig.ident), sc_summary$patient_id[idx_dt])
  
  if( (! 'project_id' %in% colnames(meta)) ){
    meta$project_id <- NA
  }
  if( NA %in% unique(meta$project_id) ){
    na_proj_idx <- is.na(meta$project_id)
    # View(meta[na_proj_idx,])
    print(table(na_proj_idx))
    
    index_h <- sc_summary$disease %in% c('adjacent tumor/disease section','healthy')
    meta_h <- meta[na_proj_idx, ]
    meta_h <- meta_h[meta_h$patient_id %in% sc_summary$patient_id[idx_tis & index_h ], ]
    idx_dt_match <- match(meta_h$patient_id, sc_summary$patient_id[idx_tis])
    table(is.na(idx_dt_match))
    meta_h$project_id <- sc_summary$project_id[idx_tis][idx_dt_match]
    
    meta_d <- meta[na_proj_idx, ]
    meta_d <- meta_d[ ! meta_d$patient_id %in% sc_summary$patient_id[idx_tis & index_h ], ]
    if(nrow(meta_d) !=0){
      meta_d$project_id <- dataset
      meta <- rbind(meta[! na_proj_idx, ], meta_h, meta_d)
    }else{
      meta <- rbind(meta[! na_proj_idx, ], meta_h)
    }
    
  }
  print(table(is.na(meta$project_id)))
  # if( NA %in% unique(meta$tissue) | (! 'tissue' %in% colnames(meta)) ){
  # meta$tissue <-  sc_summary$tissue[idx_tis][idx_dt_match]
  # }
  if( (! 'disease' %in% colnames(meta)) ){
    meta$disease <- NA
  }
  if( NA %in% unique(meta$disease) ){
    na_disea_idx <- is.na(meta$disease)
    print(table(na_disea_idx))
    
    index_h <- sc_summary$disease %in% c('adjacent tumor/disease section','healthy')
    meta_d <- meta[na_disea_idx, ]
    meta_d <- meta_d[ ! meta_d$patient_id %in% sc_summary$patient_id[idx_tis & index_h ], ]
    idx_dt_match <- match(meta_d$patient_id, sc_summary$patient_id[idx_dt])
    table(is.na(idx_dt_match))
    meta_d$disease <- sc_summary$disease[idx_dt][idx_dt_match]
    meta_h <- meta[na_disea_idx, ]
    meta_h <- meta_h[meta_h$patient_id %in% sc_summary$patient_id[idx_tis & index_h ], ]
    if(nrow(meta_h) !=0){
      meta_h$disease <- 'healthy'
      meta <- rbind(meta[!na_disea_idx, ], meta_h, meta_d)
    }else{
      meta <- rbind(meta[!na_disea_idx, ], meta_d)
    }
  }
  
  if( NA %in% unique(meta$project_id) | NA %in% unique(meta$disease) ){
    na_proj_idx <- is.na(meta$project_id) | is.na(meta$disease)
    # View(meta[na_proj_idx,])
    # meta$patient_id[na_proj_idx] <- meta$orig.ident[na_proj_idx]
    # Regular expressions for different formats
    barcodes <- meta$rownames[na_proj_idx]
    # Updated regex approach for all cases
    # Function to extract patient_id
    patient_ids <- extract_patient_id(barcodes)
    meta$patient_id[na_proj_idx] <- patient_ids
    # na_match_idx <- match(patient_ids, sc_summary$patient_id )
    # meta$project_id[na_proj_idx] <- sc_summary$project_id[na_match_idx]
    # meta$tissue[na_proj_idx] <- sc_summary$tissue[na_match_idx]
    # meta$disease[na_proj_idx] <- sc_summary$disease[na_match_idx]
    meta_h <- meta[na_proj_idx, ]
    index_h <- sc_summary$disease %in% c('adjacent tumor/disease section','healthy')
    meta_h <- meta_h[meta_h$patient_id %in% sc_summary$patient_id[idx_tis & index_h ], ]
    idx_dt_match <- match(meta_h$patient_id, sc_summary$patient_id[idx_tis])
    table(is.na(idx_dt_match))
    meta_h$project_id <- sc_summary$project_id[idx_tis][idx_dt_match]
    meta_h$tissue <- sc_summary$tissue[idx_tis][idx_dt_match]
    if(nrow(meta_h) !=0){
      meta_h$disease <- 'healthy'
    }
    meta_d <- meta[na_proj_idx, ]
    meta_d <- meta_d[! meta_d$patient_id %in% sc_summary$patient_id[idx_tis & index_h ], ]#[meta$orig.ident %in% sc_summary$patient_id[idx_dt], ]
    if(nrow(meta_d) !=0){
      meta_d$project_id <- dataset
    }
    idx_dt_match <- match(meta_d$patient_id, sc_summary$patient_id[idx_dt])
    table(is.na(idx_dt_match))
    meta_d$disease <- sc_summary$disease[idx_dt][idx_dt_match]
    meta_d$tissue <- sc_summary$tissue[idx_dt][idx_dt_match]
    
    meta <- rbind(meta[!na_proj_idx,],  meta_h, meta_d)
    
  }
  if(dataset %in% split_ddt_ids){
    meta$split_ids <- meta$subfolder
  }else{
    meta$split_ids <- NA
  }
  
  sc_meta <- meta %>%
    dplyr::select(barcode, rownames, patient_id, project_id, tissue, disease, sub_type, ct, seurat_clusters, split_ids) #sample_type, platform, Major_type, 
  rownames(sc_meta) <- 1:nrow(sc_meta)
  sc_meta$ddt_id <- dataset
  sc_meta$ddt_tis <- tissue
  print(table(sc_meta$project_id, sc_meta$ddt_id))
  print(table(sc_meta$project_id, sc_meta$patient_id))
  
  if(any(is.na(meta$project_id))){
    message(paste0("文件 ", allmeta_files[i], " 中 meta$project_id有NA"))
    break
  }
  sc_allmeta <<- rbind(sc_allmeta, sc_meta)
}
sc_allmeta$project_id <- gsub('https://ngdc.cncb.ac.cn/bioproject/browse/PRJCA001063', 'PRJCA001063', sc_allmeta$project_id)
# check the project_id
unique(sc_allmeta$project_id) # 269-->270
unique(sc_allmeta$tissue) # 63
setdiff(RMzyme_alltissues, unique(sc_allmeta$tissue))
# [1] "fallopian tube"  "rectum"          "skeletal muscle"
# [4] "thymus"
diff_proj <- setdiff(RMzyme_alldatasets, unique(sc_allmeta$project_id))
# View(sc_summary[sc_summary$project_id %in% diff_proj, ])

# check the patient_id
table(is.na(sc_allmeta$patient_id))
index <- sc_allmeta$patient_id==''
table(index)
unique(sc_allmeta$patient_id[sc_allmeta$project_id=='GSE160269'])
sc_allmeta$patient_id[sc_allmeta$project_id=='GSE160269'] <- gsub('-I', '-E', sc_allmeta$patient_id[sc_allmeta$project_id=='GSE160269'])
index = sc_allmeta$project_id=='GSE160269'
na_patientid <-  sc_allmeta$patient_id[index]
index.match <- match(na_patientid, sc_summary$patient_id)
sc_allmeta$tissue[index] <- sc_summary$tissue[index.match]
sc_allmeta$disease[index] <- sc_summary$disease[index.match]

# index = sc_allmeta$project_id=='GSE185991'
# index_scmajor = sc_majortype$project_id=='GSE185991'
# dts_barcode <- str_split(sc_majortype$barcode[index_scmajor], '[-1]', simplify = T)[, 1]
# dts_barcode <- paste0(dts_barcode, '--', sc_majortype$patient_id[index_scmajor])
# 
# index.match <- match(sc_allmeta$barcode[index], dts_barcode)
# sc_allmeta$tissue[index] <- sc_majortype$tissue[index_scmajor][index.match]
# sc_allmeta$disease[index] <- sc_majortype$disease[index_scmajor][index.match]

# na_disease_idx <- is.na(sc_allmeta$disease)
# table(unique(sc_allmeta$patient_id[na_disease_idx]) %in% sc_summary[duplicated(sc_summary$patient_id), 'patient_id'])
# View(sc_summary[duplicated(sc_summary$patient_id),])

##################################################
# change 1:
##################################################
diff_DT <- setdiff(unique(sc_allmeta$project_id), unique(sc_majortype$project_id)) 
# [1] "E-MTAB-8495" "GSE126030"  
# [3] "GSE150430"  -->:
# [1] "GSE159929" "GSE121080"
table(unique(sc_allmeta$patient_id) %in% unique(sc_majortype$patient_id))
setdiff(unique(sc_allmeta$patient_id), unique(sc_majortype$patient_id))
diff_pt <- setdiff(unique(sc_allmeta$patient_id), unique(sc_majortype$patient_id))#47

index = sc_allmeta$project_id %in% diff_DT
unique(sc_allmeta$project_id[index])
barcodes <- sc_allmeta$rownames[index]
# # Updated regex approach for all cases
# patient_ids <- sapply(barcodes, function(bc) {
#   if (str_detect(bc, "--")) {
#     # Case 1: Extract after '--'
#     str_extract(bc, "(?<=--)[^_]+")
#   } else if (str_detect(bc, "_")) {
#     # Case 2: Extract before the first underscore
#     str_extract(bc, "^[^_]+")
#   } else {
#     # Case 3: Fallback to simple barcode
#     bc
#   }
# })
patient_ids <- extract_patient_id(barcodes)
sc_allmeta$patient_id[index] <- patient_ids
patient.idx <- paste(sc_allmeta$patient_id, sc_allmeta$project_id)
match.idx <- match(patient.idx, paste(sc_summary$patient_id, sc_summary$project_id))
sc_allmeta$Major_type <- NA #sc_summary$Major_type[match.idx]
sc_allmeta$sample_type <- sc_summary$sample_type[match.idx]
sc_allmeta$platform <- sc_summary$platform[match.idx]
diff.new_majortype <- sc_allmeta[index,] %>% dplyr::select(colnames(sc_majortype))
unique(diff.new_majortype$tissue); unique(diff.new_majortype$disease)
diff.new_majortype$disease <- tolower(diff.new_majortype$disease)
# final version of sc_majortype
sc_majortype <- rbind(sc_majortype, diff.new_majortype)
sc_summary <- sc_majortype[! duplicated(paste(sc_majortype$patient_id, sc_majortype$tissue, 
                                              sc_majortype$project_id, sc_majortype$disease)),]
# saveRDS(sc_majortype, file = '/data/rluo4/RPMfunc/Output/scRNA/sc_majortype_1022.rds') # combined with diff meta from merge_nmf_meta and allmeta
sort(unique(sc_summary$disease)) #171
length(unique(sc_summary$project_id)) #396 --> 340 --> 327--> 332

# backup:
# sc_allmeta_copy <- sc_allmeta
sc_allsample <- sc_allmeta %>% dplyr::select(patient_id, project_id, tissue, disease) %>%  distinct()
table(sc_allsample$patient_id %in% sc_summary$patient_id)

##################################################
# change 2:
##################################################
bc_na_index <- is.na(sc_allmeta$barcode) | sc_allmeta$barcode=='NA'
table(is.na(sc_allmeta$barcode))
sc_allmeta$barcode[bc_na_index] <- sc_allmeta$rownames[bc_na_index]

table(sc_allmeta$disease[sc_allmeta$patient_id %in% diff_pt])  
# Healthy nasopharyngeal carcinoma              oral cancer 
# 143636                     3549                     8686 
# Acute myeloid leukemia          aplastic anemia                  healthy 
#                  41842                    11426                    48369 
#                Healthy nasopharyngeal carcinoma              oral cancer 
#                 144661                     3549                     8686 
index <- sc_allmeta$disease %in% c('Acute myeloid leukemia','aplastic anemia','nasopharyngeal carcinoma', 'oral cancer')
# View(sc_allmeta[index, ])
unique(sc_allmeta$project_id[index])
# [1] "GSE162025"  "GSE172577"  "CNP0001341" "GSE150825" 
# [1] "GSE162025"  "GSE181989"  "GSE185991"  "GSE211036" 
# [5] "GSE172577"  "CNP0001341" "GSE150825" 
table(sc_allmeta$barcode[index] == sc_allmeta$rownames[index])
table(grepl('--', sc_allmeta$rownames[index]))
# View(sc_allmeta[index,][sc_allmeta$barcode[index] != sc_allmeta$rownames[index], ])
sc_allmeta$barcode[index] <- sc_allmeta$rownames[index]
# View(sc_allmeta[grepl('[...]', sc_allmeta$rownames),])

##################################################
# change 3:
##################################################
diff_pt <- setdiff(unique(sc_allmeta$patient_id), unique(sc_majortype$patient_id))# 47 --> 106
index <- sc_allmeta$patient_id %in% diff_pt
# View(sc_allmeta[index, ])
unique(sc_allmeta[index, 'project_id'])
table(sc_allmeta$barcode[index] == sc_allmeta$rownames[index])
patient_id <- str_split(sc_allmeta$barcode[index], '--', simplify = T)[,2]
table(patient_id == '') # FALSE 156896
# View(sc_allmeta[index,][patient_id =='',])
unique(sc_allmeta[index,][patient_id =='', 'project_id'])
patient_ids <- sc_allmeta[index,][patient_id =='', 'rownames']
patient_id[patient_id==''] <- str_split(patient_ids, '_', simplify = T)[,1]
table(patient_id != sc_allmeta$patient_id[index])
# View(sc_allmeta[index,][patient_id != sc_allmeta$patient_id[index],])
table(patient_id %in% sc_summary$patient_id)
unique(sc_allmeta[index,'project_id'][ ! patient_id %in% sc_summary$patient_id])
# View(sc_allmeta[index,][! patient_id %in% sc_summary$patient_id,])
index1 <- patient_id != sc_allmeta$patient_id[index]
sc_allmeta[index,'patient_id'][index1] <- patient_id[index1]
table(sc_allmeta$project_id[index][index1])
# GSE156625 GSE161529 GSE162025 
# 1278     13571      3549 
# GSE156625 GSE161529 GSE162025 GSE181989 GSE185991 
# 1278     13571      3549     16545     80444 
# GSE188222 GSE211036 
# 2001      2647 
table(sc_allmeta$patient_id[index][index1])
table(patient_id)

diff_dt <- unique(sc_allmeta$project_id[index][index1])
table(sc_majortype$project_id %in% diff_dt)
unique(sc_majortype$patient_id[sc_majortype$project_id %in% diff_dt])
unique(sc_allmeta$patient_id[sc_allmeta$project_id %in% diff_dt])
# [1] "GSM5258386" "GSM5258387" "GSM5258390" "GSM5258388"
# [5] "GSM5258389" "GSM5258385"
# View(sc_allmeta[index,][index1,])
table(sc_allmeta$patient_id[index][index1])
table(sc_allmeta$orig.ident[index][index1])
# View(sc_allmeta[index,])
# unique(sc_allmeta[index, 'patient_id'])
diff_pt <- setdiff(unique(sc_allmeta$patient_id), unique(sc_majortype$patient_id))#37
table(sc_majortype$project_id %in% diff_dt)
diff_dt <- unique(sc_allmeta$project_id[sc_allmeta$patient_id %in% diff_pt])#10
# [1] "GSE164898" "GSE136103" "GSE156625" "GSE172577"
# [5] "GSE130973" 
diff_DT
setdiff(diff_dt, diff_DT)
# [1] "GSE164898" "GSE136103" "GSE156625" "GSE172577"
# [5] "GSE130973"
table(sc_allmeta$barcode %in% sc_majortype$barcode)
# View(sc_samples[sc_samples$project_id %in% diff_dt, ])
table(sc_allmeta$disease[sc_allmeta$project_id %in% diff_dt]) 
# Cirrhotic                  healthy                  Healthy 
# 12891                      308                   141880 
# hepatocellular carcinoma              oral cancer 
# 47534                     8686
table(is.na(sc_allmeta$barcode))
table(grepl('--', sc_allmeta$barcode))
index <- ! grepl('--', sc_allmeta$barcode)
index1 <- ! grepl('--', sc_allmeta$rownames)
sc_allmeta$barcode[index & (! index1)] <- sc_allmeta$rownames[index & (! index1)]

##################################################
# change 4:
##################################################
patient_id <- sc_allmeta$patient_id
barcode_idx <- sc_allmeta$barcode
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
# # barcode_idx[grepl('PRJCA001063', sc_allmeta$project_id)]
idx4 <- grepl('GSE164898|GSE172357', sc_allmeta$project_id) &  index
sc_allmeta$barcode[idx4][1:5]
barcode_idx[idx4] <- str_split(sc_allmeta$rownames[idx4], '[.]', simplify = T)[,2]
# barcode_idx[grepl('GSE158291', sc_allmeta$project_id)]
idx5 <- grepl('GSE158291', sc_allmeta$project_id)
sc_allmeta$barcode[idx5][1:5]
barcode_idx[idx5] <- remove_first_element(sc_allmeta$rownames[idx5])
# barcode_idx[grepl('GSE146115', sc_allmeta$project_id)]
idx6 <- grepl('GSE146115', sc_allmeta$project_id)
sc_allmeta$barcode[idx6][1:5]
barcode_idx[idx6] <- remove_first_element(sc_allmeta$rownames[idx6])
idx7 <- grepl('GSE203612', sc_allmeta$project_id) & ( ! grepl("_Vis", sc_allmeta$patient_id) ) # GSE149512
sc_allmeta$barcode[idx7][1:5]
barcode_idx[idx7] <- str_split(sc_allmeta$rownames[idx7], '.gz_', simplify = T)[,2]
idx8 <- grepl('HTA11', sc_allmeta$project_id)
sc_allmeta$barcode[idx8]
table(sc_allmeta$barcode[idx8] == sc_allmeta$rownames[idx8])
barcode_idx[idx8] <- str_split(sc_allmeta$rownames[idx8], '[-]', simplify = T)[,1]
# renew in 12/7/2024
idx9 <- grepl('GSE181254|PRJCA001063|GSE148673|GSE163203', sc_allmeta$project_id) 
sc_allmeta$barcode[idx9][1:5]
barcode_idx[idx9] <- extract_last_element(sc_allmeta$rownames[idx9])

idx_dt <-  idx4 | idx5 | idx6 | idx7 | idx8 | idx9
barcode_idx[idx_dt] <- paste0(barcode_idx[idx_dt], '--', patient_id[idx_dt])
# View(sc_allmeta[idx,])
dt <- unique(sc_allmeta[idx_dt, 'project_id'])
# setdiff( dt, unique(sc_allmeta$project_id[idx_diff]))
sc_allmeta$barcode[idx_dt] <- barcode_idx[idx_dt]

sc_allsample <- sc_allmeta[! duplicated(paste(sc_allmeta$patient_id, sc_allmeta$tissue, 
                                              sc_allmeta$project_id)),]
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
new_vector <- sc_allsample$barcode
# Check if any pattern matches for each element in the new_vector
logical_vector <- sapply(new_vector, function(x) {
  any(sapply(patterns, function(pattern) grepl(pattern, x)))
})
# View(sc_allsample[!logical_vector,])
unique(sc_allsample[!logical_vector, 'project_id'])
diff_bc <- paste(sc_allsample[!logical_vector, 'patient_id'], sc_allsample[!logical_vector, 'project_id'])
index_diff <- paste(sc_allmeta$patient_id, sc_allmeta$project_id) %in% diff_bc#sc_allmeta$patient_id %in% unique(diff_bc)
unique(sc_allmeta$project_id[index_diff])
# View(sc_allmeta[index_diff, ])
extract_clean_barcode <- function(barcodes) {
  library(stringr)
  # Define patterns to capture the barcode (ACGT sequence) more robustly
  patterns <- c(
    ".*_([ACGT]+)-\\d.*",                # General pattern ending in digits, common case
    ".*_([ACGT]+)\\.\\d.*",              # Pattern for "GSM5320155_scRNA_ctrl4_COH_Patients_CA1140_Ctrl_ACCAGTATCCTCTAGC.1"
    "^([ACGT]+)-\\d.*",                  # Pattern for "AAATCAGAAGTGATGC-0"
    ".*_([ACGT]+)_\\d.*",                # General pattern with underscores and digits
    ".*_([ACGT]+)-[A-Za-z0-9]+$",        # Pattern ending with alphanumeric characters
    ".*([ACGT]+)-[A-Za-z0-9]+$",         # Pattern ending with alphanumeric, from any location
    ".*-([ACGT]+)$",                     # Pattern for "P1T-E-AAACCTGCACTCTGTC", barcode near the end
    "([ACGT]+).*",                       # New pattern for cases like "AAAACCTCCAATGACCC-6034-YX-3-Cancer-2021"
    # "^.*_([ACGT]+)-.*$",               # Pattern for mixed content ending with barcode
    "^([ACGT]+)-[A-Za-z0-9]+$",           # Pattern for "TTTGTGTCTTATCTGT-MSS_A_2-2020"
    "^[A-Za-z0-9]*_([ACGT]+)-.*$",       # Pattern for mixed content ending with barcode
    # ".*\\.([ACGT]+)-\\d+$",             # Pattern for "GSM5252126_BPH283PrGF_Via.AAACCTGGTCCGAATT-1"
    ".*_([ACGT]+)-\\d.*",                # Pattern for barcodes ending with a hyphen and digits
    # Improved pattern to capture barcodes after a dot and before a hyphen (like "AAACCTGCAAGCCGTC-1")
    ".*\\.([ACGT]+)-\\d+$",              # General GSM-like pattern: barcode after dot and hyphen at the end
    "^[A-Za-z0-9]*_([ACGT]+)-.*$",       # Pattern for mixed content ending with barcode
    ".*_([ACGT]+)-\\d_\\d$",             # Specific pattern for barcodes with multiple digits
    ".*_([ACGT]+)-\\d+$",                # Case for "P4_S9_NORMAL_ATTCAGGTCGCTACGG-1"
    ".*_([ACGT]+)$"                      # Barcode at the very end
  )
  # Helper function to apply the patterns
  extract_barcode <- function(barcode) {
    for (pattern in patterns) {
      match <- str_match(barcode, pattern)
      if (!is.na(match[2])) {
        return(match[2])
      }
    }
    return(NA)  # Return NA if no pattern matches
  }
  # Apply the extraction function to the list of barcodes
  clean_barcodes <- sapply(barcodes, extract_barcode)
  return(clean_barcodes)
}
barcode.index <- extract_clean_barcode(sc_allmeta$barcode[index_diff])
barcode.index[grepl('GSM5022600_D2', sc_allmeta$rownames[index_diff])]
barcode.index[grepl('GSM5252460_', sc_allmeta$rownames[index_diff])]
# Check how many NAs were returned
table(is.na(barcode.index))
# See if any particular project_id has more missing barcodes
table(sc_allmeta$project_id[index_diff][is.na(barcode.index)]) #  barcode.index <- extract_clean_barcode(sc_allmeta$barcode[index])
# GSE149512 
# 20810
na_index <- is.na(barcode.index)
# View(sc_allmeta[index_diff,][na_index,])
barcode.index[na_index] <- sc_allmeta$barcode[index_diff][na_index]
sc_allmeta$barcode[index_diff] <-  paste0(barcode.index, '--', sc_allmeta$patient_id[index_diff])

##################################################
# change 5:
##################################################
index <- ! grepl('--', sc_allmeta$barcode)
index1 <- ! grepl('--', sc_allmeta$rownames)
# View( sc_allmeta[index, ] )
table(index); table(index & index1)
# index
# FALSE    TRUE 
# 7711961  216352 
# index
# FALSE    TRUE 
# 7628885  298161 
# 
# FALSE    TRUE 
# 7628885  298161 
index_dd <- sc_allmeta$barcode[index] != sc_allmeta$rownames[index]
table(index_dd)
# View(sc_allmeta[index,][index_dd, ])
patient_id <- sc_allmeta$patient_id[index&index1]
index2 <-  grepl('-1...', sc_allmeta$barcode[index & index1])
# View(sc_allmeta[index & index1, ][index2, ])
barcode.index <- extract_clean_barcode(sc_allmeta$barcode[index & index1])
# Check how many NAs were returned
table(is.na(barcode.index))
sc_allmeta$barcode[index & index1] <-  paste0(barcode.index, '--', sc_allmeta$patient_id[index & index1])
# idx_dot <- ! grepl('-1...', barcode.index)
# View(sc_allmeta[index & index1,][idx_dot, ])
table(is.na(sc_allmeta$barcode))
table(! grepl('--', sc_allmeta$barcode))
idx_dot <- grepl('[.]', sc_allmeta$barcode)
table(idx_dot); unique(sc_allmeta$patient_id[idx_dot])
# idx_dot
# FALSE    TRUE 
# 7915321   11725 
# [1] "Pt13.a" "Pt13.b" "Pt14.a" "Pt14.b" "Pt14.c"
# [6] "Pt14.d" "Pt13.c"
# View(sc_allmeta[grepl('[.]', sc_allmeta$barcode),])
table(grepl('[...]', sc_allmeta$rownames))
# FALSE    TRUE 
# 7062895  864151 
# View(sc_allmeta[grepl('[...]', sc_allmeta$rownames),])

##################################################
# change 6:
##################################################
# changes on ct of sc_allmeta
scMajortype_ct <- as.data.frame(table(sc_majortype$sub_type, sc_majortype$ct))
scMajortype_ct <- scMajortype_ct[scMajortype_ct$Freq!=0,]
index <- is.na(sc_allmeta$ct)
# View(sc_allmeta[index, ])
dts <- unique(sc_allmeta$project_id[index]) 
# [1] "GSE185381"    "GSE185991"    "E-MTAB-11948"
# [4] "GSE208653"    "E-MTAB-12305" "GSE168652"   
# [7] "SCP1950"-->
# [1] "GSE185381"    "GSE185991"    "GSE211036"   
# [4] "GSE208653"    "E-MTAB-11948" "E-MTAB-12305"
# [7] "GSE168652"    "S-BSST1035"   "SCP1950"  
table(unique(sc_majortype$project_id) %in% dts)

ct_dts <- subset(sc_majortype, project_id %in% dts)
table(sc_allmeta$barcode[index] %in% ct_dts$barcode)
table(sc_allmeta$patient_id[index] %in% ct_dts$patient_id)

index1 <- ! sc_allmeta$barcode[index] %in% ct_dts$barcode
dts <- unique(sc_allmeta[index,'project_id'][index1]) # 7 datasets
unique(sc_allmeta[index,'tissue'][index1])
unique(sc_allmeta[index,'patient_id'][index1])
# # View(sc_summary[sc_summary$project_id %in% unique(sc_allmeta[index,'project_id'][index1]),])
# ct_dts <- subset(sc_majortype, project_id %in% dts)
# unique(ct_dts$tissue) #[1] "bone marrow"
# # View(ct_dts[ct_dts$project_id %in% unique(sc_allmeta[index,'project_id'][index1]),])
# # View( sc_allmeta[index,][index1,] )
patient_id <- sc_allmeta$patient_id[index][index1]
barcode_ct <- str_split(sc_allmeta$barcode[index][index1], '--', simplify = T)[,1]
barcode_ct <- str_split(barcode_ct, '[...]', simplify = T)[, 1]
barcode_ct <- str_split(barcode_ct, '[-1]', simplify = T)[, 1]
barcode_ct <- paste0(barcode_ct, '--', patient_id)
table(barcode_ct %in% ct_dts$barcode)

barcode_ct <- c(barcode_ct, sc_allmeta$barcode[index][!index1])
table( sc_allmeta$barcode[index][index1] %in% ct_dts$barcode)

ct_dts_barcode <- ct_dts$barcode
ct_dts_barcode <- str_split(ct_dts_barcode, '-1_', simplify = T)[, 1]
index2 <-  grepl('--', ct_dts$barcode)
ct_dts_barcode[index2] <- str_split(ct_dts_barcode[index2] , '--', simplify = T)[,1]
index3 <-  grepl('_', ct_dts_barcode)
extract_last_element <- function(x) {
  split_string <- strsplit(x, "_")
  last_element <- sapply(split_string, function(y) tail(y, n = 1))
  return(last_element)
}
ct_dts_barcode[index3] <- extract_last_element(ct_dts_barcode[index3])
ct_dts_barcode <- str_split(ct_dts_barcode, '[-1]', simplify = T)[, 1]
ct_dts$barcode <- paste0(ct_dts_barcode, '--', ct_dts$patient_id)

# idx_match <- match(barcode_ct, ct_dts$barcode)
idx_match <- match( paste(sc_allmeta$project_id[index], barcode_ct), paste(ct_dts$project_id, ct_dts$barcode))
sc_allmeta$ct[index] <- ct_dts$ct[idx_match]
table(is.na(sc_allmeta$ct))
sc_allmeta$Major_type[index] <- ct_dts$Major_type[idx_match]
sc_allmeta$sub_type[index] <- ct_dts$sub_type[idx_match]
RMzyme_ct <- as.data.frame(table(sc_allmeta$sub_type, sc_allmeta$ct))
RMzyme_ct <- RMzyme_ct[RMzyme_ct$Freq!=0,] # ?Oligodendrocyte progenitor cell-Neuro, Mesothelial cell-epi/embryonic/mesenchymal

##################################################
# change 7:
##################################################
table(is.na(sc_allmeta$sub_type))
table(is.na(sc_allmeta$Major_type))
table(sc_allmeta$sub_type)
index <- sc_allmeta$sub_type=='MastCells' | sc_allmeta$sub_type=='Mast cell'
sc_allmeta$sub_type[index] <- 'MastCell'
index <- grepl('Stem-like cells', sc_allmeta$ct) & grepl('Tissue_Specific_Cells', sc_allmeta$sub_type)
sc_allmeta$sub_type[sc_allmeta$ct == 'Oligodendrocyte progenitor cell'] <- 'Neuro_Cells'
table(sc_allmeta$sub_type[index])
View(sc_allmeta[index,])
sc_allmeta$sub_type[index] <- 'Stem_Cells'

##################################################
index1 <- grepl('Bcells|Macrophages|Neutrophils|NKcells|OtherImmunecells|Tcells|MastCell|Myeloids',sc_allmeta$sub_type)
sc_allmeta$Major_type[index1] <- 'Immune_Cells'
index2 <- grepl('OtherMesenchymals|Fibroblasts',sc_allmeta$sub_type)
sc_allmeta$Major_type[index2] <- 'Mesenchymal_Cells'
index3 <- grepl('Endothelial_Cells',sc_allmeta$sub_type)
sc_allmeta$Major_type[index3] <- 'Endothelial_Cells'
index4 <- grepl('Tissue_Specific_Cells',sc_allmeta$sub_type)
sc_allmeta$Major_type[index4] <- 'Tissue_Specific_Cells'
index5 <- grepl('Stem_Cells',sc_allmeta$sub_type)
sc_allmeta$Major_type[index5] <- 'Stem_Cells'
index6 <- grepl('Blood_Cells',sc_allmeta$sub_type)
sc_allmeta$Major_type[index6] <- 'Blood_Cells'
table(sc_allmeta$Major_type)
table(is.na(sc_allmeta$Major_type))
# sc_allmeta$Major_type[sc_allmeta$Major_type=='Epithelial'] <- 'Tissue_Specific_Cells'
# sc_allmeta$Major_type[sc_allmeta$Major_type=='Stromal'] <- 'Mesenchymal_Cells'
# sc_allmeta$Major_type[sc_allmeta$Major_type=='Immune'] <- 'Immune_Cells'
index7 <- sc_allmeta$sub_type=='Stem_Cells'
table(sc_allmeta$Major_type[index7])
table(sc_allmeta$ct[index7 & sc_allmeta$Major_type=='Tissue_Specific_Cells'])
# # View(sc_allmeta[sc_allmeta$ct=='Metaplasia stem cells', ])
# sc_allmeta$Major_type[index7] <- gsub('Tissue_Specific_Cells','Stem_Cells', sc_allmeta$Major_type[index7])
index_na <- is.na(sc_allmeta$Major_type)
table(sc_allmeta$sub_type[index_na])
sc_allmeta$Major_type[index_na] <- sc_allmeta$sub_type[index_na] 
# sc_allmeta$sample_type[index_na] <- 'healthy'
# sc_allmeta$platform[index_na] <- '10x'

##################################################
# change 8:
##################################################
setdiff(unique(sc_allmeta$disease), unique(sc_majortype$disease))
setdiff(unique(sc_allmeta$tissue), unique(sc_majortype$tissue))
sort(unique(sc_allmeta$disease))
sort(unique(sc_majortype$disease))
setdiff(unique(sc_majortype$disease), unique(sc_allmeta$disease))
setdiff(unique(sc_majortype$tissue), unique(sc_allmeta$tissue))
setdiff(unique(sc_summary$disease), unique(sc_allmeta$disease))
setdiff(unique(sc_summary$tissue), unique(sc_allmeta$tissue))
setdiff(unique(sc_summary$patient_id), unique(sc_allmeta$patient_id))#
sc_allmeta$patient_id <- str_split(sc_allmeta$patient_id, '[...]', simplify = T)[, 1]
merge_patient <- unique(paste(sc_allmeta$project_id, sc_allmeta$patient_id)) #4463 (sc_majortype) --> 4048
table(merge_patient %in% sc_patient_id)
sc_patient_id <- paste(sc_summary$project_id, sc_summary$patient_id)
diff_pt <- setdiff(merge_patient, sc_patient_id)#setdiff(unique(sc_allmeta$patient_id),unique(sc_summary$patient_id))
index <- paste(sc_allmeta$project_id, sc_allmeta$patient_id) %in% diff_pt
diff_dt <- sc_allmeta[index,]
unique(diff_dt$project_id)
# [1] "GSE164898" "GSE171555" "GSE159929" "GSE134520"
# [5] "GSE136103" "GSE156625" "GSE112271" "GSE146409"
# [9] "GSE126030" "GSE121080" "GSE172577" "GSE137829"
# [13] "GSE130973" "GSE144236" "GSE167297"
# View(sc_summary[sc_summary$project_id %in% unique(diff_dt$project_id),])
diff_dt$patient_id <- str_split(diff_dt$barcode, '--', simplify = T)[,2]
table(unique(diff_dt$project_id) %in% tissue_summary$Dataset)
unique(diff_dt$disease)
sc_allmeta[index,'patient_id'] <- diff_dt$patient_id
merge_patient <- paste(sc_allmeta$project_id, sc_allmeta$patient_id)
sc_patient_id <- paste(sc_summary$project_id, sc_summary$patient_id)
index1 <- merge_patient %in% sc_patient_id
index2 <- match(merge_patient[index1], sc_patient_id)
sc_allmeta$disease[index1] <- sc_summary$disease[index2]
sc_allmeta$sample_type[index1] <- sc_summary$sample_type[index2]
sc_allmeta$tissue[index1] <- sc_summary$tissue[index2]
# View(sc_allmeta[ ! index1, ])
setdiff(unique(sc_summary$disease), unique(sc_allmeta$disease)) 
setdiff(unique(sc_summary$sample_type), unique(sc_allmeta$sample_type)) 
# test <- subset(sc_allmeta, project_id == 'GSE136103')
# unique(test$patient_id)
# test_major <- subset(sc_majortype, project_id == 'GSE136103')
# unique(test_major$patient_id)
# setdiff(unique(test$patient_id), unique(test_major$patient_id))
setdiff( unique(sc_allmeta$disease), unique(sc_summary$disease))
setdiff( unique(sc_allmeta$tissue), unique(sc_summary$tissue))
sc_allmeta$disease <- gsub('HCC', 'hepatocellular carcinoma', 
                           gsub('Healthy', 'healthy',
                                gsub('oral cancer', 'oral squamous cell carcinoma',
                                     gsub('PRAD', 'pancreatic ductal adenocarcinoma',
                                          gsub('cSCC', 'oral squamous cell carcinoma',
                                               gsub('AML', 'acute myeloid leukemia',
                                                    sc_allmeta$disease))))))
# unique(sc_allmeta$patient_id)

# sc_allmeta <- sc_allmeta %>% dplyr::select( c(colnames(sc_majortype), 'seurat_clusters') )
unique(sc_allmeta$project_id) #269
unique(sc_allmeta$patient_id) #3596 --> 3645
unique(sc_allmeta$platform)
index1 <- is.na(sc_allmeta$platform)
# table(sc_allmeta$patient_id %in% c( disco$patient_id, tissue_summary$Replicates) )
index2 <- match(sc_allmeta$patient_id[index1], disco$patient_id)
table(is.na(index2))
sc_allmeta$platform[index1] <- disco$platform[index2]
index1 <- is.na(sc_allmeta$platform)
unique(sc_allmeta$project_id[index1]) #[1] "CRA002497"   "GSE155673"   "covid_atlas" "GSE155224"
sc_allmeta$platform[sc_allmeta$project_id=='CRA002497'] = "10x5'v2"
sc_allmeta$platform[sc_allmeta$project_id=='GSE155673'] = "10x3'v3"
sc_allmeta$platform[sc_allmeta$project_id=='covid_atlas'] = "10x5'v2"
sc_allmeta$platform[sc_allmeta$project_id=='GSE155224'] = "10x5'"#"10x3'v3" #only pbmc? where are the blood tissues?
sc_allmeta$platform[sc_allmeta$platform=='10x3'] = "10x3'"
sc_allmeta$platform[sc_allmeta$project_id=='GSE161529'] = "10x3'"
sc_allmeta$platform[sc_allmeta$project_id=='GSE159929'] = "10x3'v3" 
sc_allmeta$platform[sc_allmeta$project_id=='GSE121080'] = "10x3'v2" 
sc_allmeta$platform[sc_allmeta$project_id=='GSE172577'] = "10x3'v3" 
unique(sc_allmeta$platform)
# [1] "10x3'v2"       "10x3'v3"       "10x3'"        
# [4] "10x5'"         "10x5'v2"       "10X Genomics" 
# [7] "10x"           "10x3'v1"       "inDrop V2"    
# [10] "BD Rhapsody"   "Microwell-seq" "Seq-well"     
# [13] "MARS-seq"
unique(sc_allmeta$project_id[sc_allmeta$platform=='10x'] )
sc_allmeta$platform[sc_allmeta$project_id=='GSE149689'] = "10x3'v3"
sc_allmeta$platform[sc_allmeta$project_id=='GSE168453'] = "10x5'"
sc_allmeta$platform[sc_allmeta$project_id=='GSE166992'] = "10x5'"
sc_allmeta$platform[sc_allmeta$project_id=='GSE171555'] = "10x5'v1"
sc_allmeta$platform[sc_allmeta$project_id=='GSE192391'] = "10x5'"
sc_allmeta$platform[sc_allmeta$project_id=='GSE154567'] = "10x3'v3"
sc_allmeta$platform[sc_allmeta$project_id=='GSE208337'] = "10x3'v3"
sc_allmeta$platform[sc_allmeta$project_id=='GSE150861'] = "10x3'v2"
unique(sc_allmeta$project_id[sc_allmeta$platform=='10X Genomics'] )
sc_allmeta$platform[sc_allmeta$platform=='10X Genomics'] = "10x"
NMF_clusters <- unique(sc_allmeta$seurat_clusters)
NMF_RMPs <- str_split(NMF_clusters, "[+]", simplify = T)[,1]
sort(unique(NMF_RMPs))
# index <- grepl("NA+", sc_allmeta$seurat_clusters)
# View(sc_allmeta[index,])
# NMF_clusters <- gsub('^NA\\+$', 'NoneRMP+', NMF_clusters)
sc_allmeta$seurat_clusters <- gsub('NA+', 'NoneRMP+', sc_allmeta$seurat_clusters, fixed = TRUE)
sc_allmeta$seurat_clusters <- gsub('MastCell', 'MastCells',  sc_allmeta$seurat_clusters, fixed = TRUE)
sc_allmeta$sub_type <- gsub('MastCell', 'MastCells',  sc_allmeta$sub_type, fixed = TRUE)
table(is.na(sc_allmeta$seurat_clusters))
saveRDS(sc_allmeta, file = '/data/rluo4/RPMfunc/Output/scRNA/sc_allmeta.rds')
# remove(sc_allmeta_copy);
# sc_samples <- sc_allmeta[! duplicated(paste(sc_allmeta$patient_id, sc_allmeta$tissue, 
#                                              sc_allmeta$project_id)),]
################################################################################
#      organize the diseased data of disco named as sc_cell_summary            #
################################################################################
sc_cell_summary <- sc_allmeta[! grepl('adjacent|ealthy', sc_allmeta$disease),] # 219 disease datasets
table(sc_cell_summary$ddt_id == sc_cell_summary$project_id)
index <- sc_cell_summary$ddt_id == sc_cell_summary$project_id
sc_cell_h <- (sc_cell_summary[!index, ])
sc_cell_h$rownames <- str_split(sc_cell_h$rownames, '[...]', simplify = T)[, 1]
table(sc_cell_h$rownames %in% sc_majortype$barcode)
index_pairs <- paste(sc_majortype$project_id, sc_majortype$patient_id)
index_pair <- paste(sc_cell_summary$project_id, sc_cell_summary$patient_id)[! index]
# sc_test <- sc_majortype[index_pairs %in% unique(index_pair),]
# sc_test$barcode <- str_split(sc_cell_h$rownames, '-1', simplify = T)[, 1]
idx_h_match <- match(sc_cell_h$rownames, sc_majortype$barcode)
print(table( is.na(idx_h_match)))
sc_cell_h$disease_sc <- sc_majortype$disease[idx_h_match]
idx_h_na <- is.na(sc_cell_h$disease_sc)
print(table(idx_h_na))
unique(sc_cell_h$project_id[ idx_h_na ])
idx_h_na_match <- match(sc_cell_h$barcode[idx_h_na], sc_majortype$barcode)
sc_cell_h$disease_sc[idx_h_na] <- sc_majortype$disease[idx_h_na_match]
idx_h_na <- is.na(sc_cell_h$disease_sc)
print(table(idx_h_na))
unique(sc_cell_h$project_id[ idx_h_na ])
# View( sc_cell_h[idx_h_na, ] )
idx_h_na_match <- match(sc_cell_h$barcode[idx_h_na], gsub("-1--", "--",sc_majortype$barcode) )
sc_cell_h$disease_sc[idx_h_na] <- sc_majortype$disease[idx_h_na_match]
idx_h_na <- is.na(sc_cell_h$disease_sc)
print(table(idx_h_na))
unique(sc_cell_h$project_id[ idx_h_na ])
# View( sc_cell_h[idx_h_na, ] )
# idx_h_na_match <- match( str_split(sc_cell_h$rownames[idx_h_na], '-1', simplify = T)[, 1],
#                          str_split(sc_majortype$barcode, '-1', simplify = T)[, 1],  )
# sc_cell_h$disease_sc[idx_h_na] <- sc_majortype$disease[idx_h_na_match]
# table(sc_cell_h$disease_sc)
# table(sc_cell_h$disease == sc_cell_h$disease_sc)
# index_equal <- sc_cell_h$disease == sc_cell_h$disease_sc
# # View(sc_cell_h[index_equal, ])
# table(sc_cell_h[index_equal, 'disease_sc'])
# table(sc_cell_h[!index_equal, 'disease_sc'])

idx_h_renew <- ! grepl('adjacent|ealthy', sc_allmeta$disease) & sc_allmeta$ddt_id != sc_allmeta$project_id
sc_allmeta[idx_h_renew, 'disease'] <- sc_cell_h$disease_sc

sc_allmeta <- readRDS(file = '/data/rluo4/RPMfunc/Output/scRNA/sc_allmeta_1022.rds')
sc_allmeta$omics <- 'scRNA'
sc_meta <- readRDS('/data/rluo4/RPMfunc/Output/scRNA/GSE185381_allmeta.rds')
sc_allmeta <- sc_allmeta[!sc_allmeta$ddt_id %in% 'GSE185381', ]
sc_allmeta <- rbind(sc_allmeta, sc_meta)

sc_allsample <- sc_allmeta %>% dplyr::select(patient_id, project_id, tissue, disease) %>%  distinct()
table(sc_allsample$patient_id %in% sc_summary$patient_id)
sc_summary <- sc_majortype[! duplicated(paste(sc_majortype$patient_id, sc_majortype$tissue, 
                                              sc_majortype$project_id, sc_majortype$disease)),]
saveRDS(sc_majortype, file = '/data/rluo4/RPMfunc/Output/scRNA/sc_majortype_1122.rds') # combined with diff meta from merge_nmf_meta and allmeta
saveRDS(sc_allmeta, file = '/data/rluo4/RPMfunc/Output/scRNA/sc_allmeta_1122.rds')
save(disco, tissue_summary, sc_summary, sc_allsample, file = '/data/rluo4/RPMfunc/Output/scRNA/sc_allsample_1122.RData')

