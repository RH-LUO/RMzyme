#!/usr/local/bin/Rscript
# title: "DISCO: a database of Deeply Integrated human Single-Cell Omics data"
# title:
# author: "Ruihan Luo"
# date: "March 28th,2024"
# rm(list=ls())
Sys.setenv("VROOM_CONNECTION_SIZE" = 999999999)
getOption('timeout')
options(timeout=100000000)
path <- c("/home/rluo4/R/x86_64-pc-linux-gnu-library/4.1" , '/home/rluo4/R/x86_64-pc-linux-gnu-library/4.3', '/opt/R/4.3.1/lib64/R/library', '/home/rluo4/R/x86_64-conda-linux-gnu-library/4.3', '/data2/rluo4/bin/miniconda3/lib/R/library')
path <- c('/home/rluo4/R/x86_64-pc-linux-gnu-library/4.3', '/opt/R/4.3.1/lib64/R/library', '/home/rluo4/R/x86_64-conda-linux-gnu-library/4.3', '/data2/rluo4/bin/miniconda3/lib/R/library')
.libPaths(path)
################################################################################
# 1. Adjust the sample_type & disease data
################################################################################
data_dir <- '/data2/rluo4/EpiTrans/DataCollection/disco_tsv'
setwd(data_dir)
disco_original <- readRDS('/data2/rluo4/RPMfunc/xuhaixia/disco_all_meta.rds')
disco <- disco_original
disco$patient_id[is.na(disco$patient_id)] <- str_split(rownames(disco)[is.na(disco$patient_id)], 
                                                       '--', simplify = T)[,2]
unique(disco$project_id)
unique(disco$sample_type)
unique(disco$disease)
table(disco$tissue)
table(is.na(disco$tissue))
table(disco$project_id[disco$tissue==''])
# E-MTAB-11536  E-MTAB-9389    GSE185381 
# 7407       105766        46013 
# 1) Check 20+ datasets with NA in tissue information
################################################################################
tissue_na.index = disco$tissue==''
# tissue_na <- disco[disco$tissue=='',]
#  Cross-tissue immune cell analysis reveals tissue-specific features in humans:
tissue_na1 <- read.table('/data2/rluo4/RPMfunc/disco_pdata/E-MTAB-11536.sdrf.txt', sep = '\t', header = T)
unique(tissue_na1$Comment.ENA_EXPERIMENT.)
table(tissue_na1$Characteristics.organism.part.)
table(unique(disco$patient_id) %in% tissue_na1$Comment.ENA_EXPERIMENT.)
tissue_na <- disco[tissue_na.index & disco$project_id=='E-MTAB-11536',]
unique(tissue_na$patient_id)
tissue_na$tissue <- tissue_na1$Characteristics.organism.part.[
  match(tissue_na$patient_id,tissue_na1$Comment.ENA_EXPERIMENT.)]
table(tissue_na$tissue)
disco[tissue_na.index & disco$project_id=='E-MTAB-11536','tissue'] <- tissue_na$tissue
table(disco[tissue_na.index & disco$project_id=='E-MTAB-11536','disease'])
disco[tissue_na.index & disco$project_id=='E-MTAB-11536','disease'] <- 'Healthy control'
disco[tissue_na.index & disco$project_id=='E-MTAB-11536','sample_type'] <- 'Healthy'
# Single cell analysis of emergent haematopoiesis in the human fetal bone marrow (10X)
tissue_na2 <- read.table('/data2/rluo4/RPMfunc/disco_pdata/E-MTAB-9389.sdrf.txt', sep = '\t', header = T)
unique(tissue_na2$Comment.ENA_EXPERIMENT.)
table(tissue_na2$Characteristics.organism.part.)
table(unique(disco$patient_id) %in% tissue_na2$Comment.ENA_EXPERIMENT.)
tissue_na <- disco[tissue_na.index & disco$project_id=='E-MTAB-9389',]
unique(tissue_na$patient_id)
tissue_na$tissue <- tissue_na2$Characteristics.organism.part.[
  match(tissue_na$patient_id,tissue_na2$Comment.ENA_EXPERIMENT.)]
unique(tissue_na$tissue)
tissue_na$disease <- tissue_na2$Characteristics.disease.[
  match(tissue_na$patient_id,tissue_na2$Comment.ENA_EXPERIMENT.)]
disco[tissue_na.index & disco$project_id=='E-MTAB-9389','tissue'] <- tissue_na$tissue
disco[tissue_na.index & disco$project_id=='E-MTAB-9389','disease'] <- tissue_na$disease
disco[tissue_na.index & disco$project_id=='E-MTAB-9389','sample_type'] <- 'normal'

tissue_na3 <- read.csv('/data2/rluo4/RPMfunc/disco_pdata/GSE185381_SRA_Table.txt',  header = T)
unique(tissue_na3$disease_state)
tissue_na <- disco[tissue_na.index & disco$project_id=='GSE185381',]
unique(tissue_na$patient_id)
tissue_na$tissue <- 'bone marrow'#tissue_na3$source_name
unique(tissue_na$tissue)
disco[tissue_na.index & disco$project_id=='GSE185381','tissue'] <- 'bone marrow'
disco[tissue_na.index & disco$project_id=='GSE185381','disease'] <- 'Healthy Control'
disco[tissue_na.index & disco$project_id=='GSE185381','sample_type'] <- 'Healthy'

tissue_na <- disco[disco$tissue == 'unknown',]
# tissue_na$patient_id <- str_split(rownames(tissue_na), '--', simplify = T)[,2]
disco$tissue[disco$tissue == 'unknown'] = 'head and neck'
disco$patient_id[disco$tissue == 'unknown'] = tissue_na$patient_id

table(disco$tissue)
table(is.na(disco$tissue))
# View(disco[is.na(disco$Major_type),])
# table(disco[is.na(disco$Major_type),]$ct)
table(disco$disease)
table(disco$sample_type)
table(is.na(disco$sample_type))
# disco$Major_type[is.na(disco$Major_type)] <- 'NonEpi'
# table(disco$Major_type)
# disco_epi <- disco[disco$Major_type == 'Epithelial',]
# str(disco_epi)
# table(is.na(disco_epi$sample_type))
# table(disco_epi$project_id)
# na.index <- which(is.na(disco_epi$sample_type))
# View(disco_epi[na.index,])
# table(disco_epi$project_id[na.index])
# disco_diseased_epi <-  disco_epi[-grep("ormal", disco_epi$sample_type),]
sample.index <- c(grep("ormal", disco$sample_type),grep("ealthy", disco$sample_type))
disco[sample.index, 'disease'] <- 'Healthy'
table(is.na(disco$disease))
disease_na.index <- is.na(disco$disease)
disease_na <- disco[disease_na.index,]
table(disease_na$project_id)
# E-MTAB-8884                      E-MTAB-9139 
# 1601                            67679 
# GSE120221                        GSE128639 
# 81752                             8501 
# GSE130430                        GSE133181 
# 5225                            26548 
# GSE135194                        GSE139369 
# 14610                            13254 
# GSE159624                        GSE159929 
# 21718                             2752 
# GSE163668                        GSE165645 
# 11930                             5454 
# GSE166895                        GSE169426 
# 22309                             3072 
# GSE175604                        GSE179346 
# 22389                             1845 
# GSE181989                        GSE184198 
# 5440                              798 
# GSE188222                        GSE193138 
# 1152                            10708 
# GSE202735 HCA_HematopoieticImmuneCellAtlas 
# 93                           198555 
# 2) Check 20+ datasets with NA in disease information
disease_na1 <- read.table('/data2/rluo4/RPMfunc/disco_pdata/E-MTAB-8884.sdrf.txt', sep = '\t', header = T)
unique(disease_na1$Comment.ENA_EXPERIMENT.)
table(disease_na1$Characteristics.disease.)
table(unique(disco$patient_id) %in% disease_na1$Comment.ENA_EXPERIMENT.)
disease_na <- disco[disease_na.index & disco$project_id=='E-MTAB-8884',]
unique(disease_na$patient_id)
unique(disease_na$tissue)
unique(disease_na$sample_type)
disease_na$disease <- disease_na1$Characteristics.disease.[
  match(disease_na$patient_id,disease_na1$Comment.ENA_EXPERIMENT.)]
unique(disease_na$disease)
disco[disease_na.index & disco$project_id=='E-MTAB-8884','disease'] <- disease_na$disease
disco[disease_na.index & disco$project_id=='E-MTAB-8884','sample_type'] <- 'normal'

disease_na2 <- read.table('/data2/rluo4/RPMfunc/disco_pdata/E-MTAB-9139.sdrf.txt', sep = '\t', header = T)
unique(disease_na2$Comment.ENA_EXPERIMENT.)
table(disease_na2$Characteristics.disease.)
table(unique(disco$patient_id) %in% disease_na2$Comment.ENA_EXPERIMENT.)
disease_na <- disco[disease_na.index & disco$project_id=='E-MTAB-9139',]
unique(disease_na$patient_id)
unique(disease_na$tissue)
unique(disease_na$sample_type)
disease_na$disease <- disease_na2$Characteristics.disease.[
  match(disease_na$patient_id,disease_na2$Comment.ENA_EXPERIMENT.)]
unique(disease_na$disease)
disco[disease_na.index & disco$project_id=='E-MTAB-9139','disease'] <- disease_na$disease
disco[disease_na.index & disco$project_id=='E-MTAB-9139','sample_type'] <- 'normal'

setwd('/data2/rluo4/RPMfunc/disco_pdata')
gse <- getGEO('GSE120221',destdir=".", AnnotGPL=F, getGPL=F)
disease_na3 <- pData(gse[[1]])
disco[disease_na.index & disco$project_id=='GSE120221','disease'] <- 'normal'
disco[disease_na.index & disco$project_id=='GSE120221','sample_type'] <- 'normal'

gse <- getGEO('GSE128639',destdir=".", AnnotGPL=F, getGPL=F)# scRNA + CITE-seq
disease_na3 <- pData(gse[[1]])
disco[disease_na.index & disco$project_id=='GSE128639','disease'] <- 'normal'
disco[disease_na.index & disco$project_id=='GSE128639','sample_type'] <- 'normal'

gse <- getGEO('GSE130430',destdir=".", AnnotGPL=F, getGPL=F)#
disease_na3 <- pData(gse[[1]])
disco[disease_na.index & disco$project_id=='GSE130430','disease'] <- 'normal'
disco[disease_na.index & disco$project_id=='GSE130430','sample_type'] <- 'normal'

gse <- getGEO('GSE133181',destdir=".", AnnotGPL=F, getGPL=F)#
disease_na3 <- pData(gse[[1]])
disease_na <- disco[disease_na.index & disco$project_id=='GSE133181',]
# View(disco[disco$project_id=='GSE133181',])
unique(disco$patient_id[disco$project_id=='GSE133181'])
unique(disease_na$patient_id)
unique(disco$disease[disco$project_id=='GSE133181'])
unique(disco$patient_id[disco$project_id=='GSE133181'])

unique(disease_na$tissue)
unique(disease_na$sample_type)
disease_na$disease <- disease_na3$`individual:ch1`[
  match(disease_na$patient_id,disease_na3$geo_accession)]
unique(disease_na$disease)
disease_na$disease <- gsub(' donor','', disease_na$disease)
disco[disease_na.index & disco$project_id=='GSE133181','disease'] <- disease_na$disease
disco[disease_na.index & disco$project_id=='GSE133181','sample_type'] <- disease_na$disease


gse <- getGEO('GSE135194',destdir=".", AnnotGPL=F, getGPL=F)#
disease_na3 <- pData(gse[[1]])
disease_na <- disco[disease_na.index & disco$project_id=='GSE135194',]
unique(disease_na$patient_id)
unique(disease_na$tissue)
unique(disease_na$sample_type)
unique(disco$patient_id[disco$project_id=='GSE135194'])
disease_na$disease <- disease_na3$title[
  match(disease_na$patient_id,disease_na3$geo_accession)]
unique(disease_na$disease)
disease_na$disease <- str_split(disease_na$disease, ' donor',simplify = T)[,1]
disco[disease_na.index & disco$project_id=='GSE135194','disease'] <- disease_na$disease
disco[disease_na.index & disco$project_id=='GSE135194','sample_type'] <- disease_na$disease

gse <- getGEO('GSE139369',destdir=".", AnnotGPL=F, getGPL=F)#
disease_na3 <- rbind(pData(gse[[1]]), pData(gse[[2]]) )
disease_na <- disco[disease_na.index & disco$project_id=='GSE139369',]
unique(disease_na$patient_id)
disease_na3$geo_accession
unique(disease_na$tissue)
unique(disease_na$sample_type)
unique(disco$patient_id[disco$project_id=='GSE139369'])
disease_na$disease <- disease_na3$characteristics_ch1[
  match(disease_na$patient_id,disease_na3$geo_accession)]
unique(disease_na$disease)
disease_na$disease <- str_split(disease_na$disease, 'disease state: ',simplify = T)[,2]
disco[disease_na.index & disco$project_id=='GSE139369','disease'] <- disease_na$disease
disco[disease_na.index & disco$project_id=='GSE139369','sample_type'] <- disease_na$disease

gse <- getGEO('GSE159624',destdir=".", AnnotGPL=F, getGPL=F)#
disease_na3 <- pData(gse[[1]])
disease_na <- disco[disease_na.index & disco$project_id=='GSE159624',]
unique(disease_na$patient_id)
disease_na3$geo_accession
table(disease_na$patient_id %in% disease_na3$geo_accession)
unique(disease_na$tissue)
unique(disease_na$sample_type)
unique(disco$patient_id[disco$project_id=='GSE159624'])
disease_na$disease <- disease_na3$`patient diagnosis:ch1`[
  match(disease_na$patient_id,disease_na3$geo_accession)]
unique(disease_na$disease)
# disease_na$disease <- str_split(disease_na$disease, 'disease state: ',simplify = T)[,2]
disco[disease_na.index & disco$project_id=='GSE159624','disease'] <- disease_na$disease
disco[disease_na.index & disco$project_id=='GSE159624','sample_type'] <- disease_na$disease


gse <- getGEO('GSE159929',destdir=".", AnnotGPL=F, getGPL=F)#
disease_na3 <- pData(gse[[1]])
disease_na <- disco[disease_na.index & disco$project_id=='GSE159929',]
unique(disease_na$patient_id)
disease_na3$geo_accession
table(disease_na$patient_id %in% disease_na3$geo_accession)
unique(disease_na$tissue)
unique(disease_na$sample_type)
unique(disco$patient_id[disco$project_id=='GSE159929'])
disco[disease_na.index & disco$project_id=='GSE159929','disease'] <- 'Healthy'
disco[disease_na.index & disco$project_id=='GSE159929','sample_type'] <- 'Healthy'

gse <- getGEO('GSE163668',destdir=".", AnnotGPL=F, getGPL=F)#
disease_na3 <- pData(gse[[1]])
disease_na <- disco[disease_na.index & disco$project_id=='GSE163668',]
unique(disease_na$patient_id)
disease_na3$geo_accession
setdiff(unique(disease_na$patient_id),disease_na3$geo_accession)
disease_na$patient_id <- substr(disease_na$patient_id, 1, 10)
table(disease_na$patient_id %in% disease_na3$geo_accession)
unique(disease_na$tissue)
unique(disease_na$sample_type)
unique(disease_na$disease)
unique(disco$patient_id[disco$project_id=='GSE163668'])
disease_na$sample_type <- disease_na3$`covid_status:ch1`[
  match(disease_na$patient_id,disease_na3$geo_accession)]
unique(disease_na$disease)
# disease_na$disease <- str_split(disease_na$disease, ' donor',simplify = T)[,1]
disco[disease_na.index & disco$project_id=='GSE163668','disease'] <- 'COVID-19'
disco[disease_na.index & disco$project_id=='GSE163668','sample_type'] <- disease_na$sample_type

gse <- getGEO('GSE165645',destdir=".", AnnotGPL=F, getGPL=F)#
disease_na3 <- pData(gse[[1]])
disease_na <- disco[disease_na.index & disco$project_id=='GSE165645',]
unique(disease_na$patient_id)
disease_na3$geo_accession
table(disease_na$patient_id %in% disease_na3$geo_accession)
unique(disease_na$tissue)
unique(disease_na$sample_type)
unique(disco$patient_id[disco$project_id=='GSE165645'])
disco[disease_na.index & disco$project_id=='GSE165645','disease'] <- 'Healthy'
disco[disease_na.index & disco$project_id=='GSE165645','sample_type'] <- 'Healthy'

gse <- getGEO('GSE166895',destdir=".", AnnotGPL=F, getGPL=F)#
disease_na3 <- pData(gse[[1]])
disease_na <- disco[disease_na.index & disco$project_id=='GSE166895',]
unique(disease_na$patient_id)
disease_na3$geo_accession
table(disease_na$patient_id %in% disease_na3$geo_accession)
unique(disease_na$tissue)
unique(disease_na$sample_type)
unique(disco$patient_id[disco$project_id=='GSE166895'])
unique(disco$disease[disco$project_id=='GSE166895'])
disease_na$sample_type <- disease_na3$`tissue:ch1`[
  match(disease_na$patient_id,disease_na3$geo_accession)]
unique(disease_na$sample_type)
disco[disease_na.index & disco$project_id=='GSE166895','disease'] <- 'Healthy'
disco[disease_na.index & disco$project_id=='GSE166895','sample_type'] <- 'Healthy'

gse <- getGEO('GSE169426',destdir=".", AnnotGPL=F, getGPL=F)#
disease_na3 <- pData(gse[[1]])
disease_na <- disco[disease_na.index & disco$project_id=='GSE169426',]
unique(disease_na$patient_id)
disease_na3$geo_accession
table(disease_na$patient_id %in% disease_na3$geo_accession)
unique(disease_na$tissue)
unique(disease_na$sample_type)
unique(disco$patient_id[disco$project_id=='GSE169426'])
unique(disco$disease[disco$project_id=='GSE169426'])
disease_na$disease <- disease_na3$title[
  match(disease_na$patient_id,disease_na3$geo_accession)]
unique(disease_na$disease)
disco[disease_na.index & disco$project_id=='GSE169426','disease'] <- 'Healthy'
disco[disease_na.index & disco$project_id=='GSE169426','sample_type'] <- 'Healthy'

gse <- getGEO('GSE175604',destdir=".", AnnotGPL=F, getGPL=F)#
disease_na3 <- pData(gse[[1]])
disease_na <- disco[disease_na.index & disco$project_id=='GSE175604',]
unique(disease_na$patient_id)
disease_na3$geo_accession
table(disease_na$patient_id %in% disease_na3$geo_accession)
unique(disease_na$tissue)
unique(disease_na$sample_type)
unique(disco$patient_id[disco$project_id=='GSE175604'])
unique(disco$disease[disco$project_id=='GSE175604'])
disease_na$disease <- disease_na3$title[
  match(disease_na$patient_id,disease_na3$geo_accession)]
unique(disease_na$disease)
disco[disease_na.index & disco$project_id=='GSE175604','disease'] <- 'Healthy'
disco[disease_na.index & disco$project_id=='GSE175604','sample_type'] <- 'Healthy'

gse <- getGEO('GSE179346',destdir=".", AnnotGPL=F, getGPL=F)#
disease_na3 <- pData(gse[[1]])
disease_na <- disco[disease_na.index & disco$project_id=='GSE179346',]
unique(disease_na$patient_id)
disease_na3$geo_accession
table(disease_na$patient_id %in% disease_na3$geo_accession)
unique(disease_na$tissue)
unique(disease_na$sample_type)
unique(disco$patient_id[disco$project_id=='GSE179346'])
unique(disco$disease[disco$project_id=='GSE179346'])
disease_na$disease <- disease_na3$`disease state:ch1`[
  match(disease_na$patient_id,disease_na3$geo_accession)]
unique(disease_na$disease)
disco[disease_na.index & disco$project_id=='GSE179346','disease'] <- disease_na$disease
disco[disease_na.index & disco$project_id=='GSE179346','sample_type'] <- disease_na$disease

gse <- getGEO('GSE181989',destdir=".", AnnotGPL=F, getGPL=F)#
disease_na3 <- pData(gse[[1]])
disease_na <- disco[disease_na.index & disco$project_id=='GSE181989',]
unique(disease_na$patient_id)
disease_na3$geo_accession
table(disease_na$patient_id %in% disease_na3$geo_accession)
unique(disease_na$tissue)
unique(disease_na$sample_type)
unique(disco$patient_id[disco$project_id=='GSE181989'])
unique(disco$disease[disco$project_id=='GSE181989'])
disease_na$disease <- disease_na3$`disease state:ch1`[
  match(disease_na$patient_id,disease_na3$geo_accession)]
unique(disease_na$disease)
disco[disease_na.index & disco$project_id=='GSE181989','disease'] <- 'Healthy'
disco[disease_na.index & disco$project_id=='GSE181989','sample_type'] <- 'Healthy'

gse <- getGEO('GSE184198',destdir=".", AnnotGPL=F, getGPL=F)#
disease_na3 <- pData(gse[[1]])
disease_na <- disco[disease_na.index & disco$project_id=='GSE184198',]
unique(disease_na$patient_id)
disease_na3$geo_accession
table(disease_na$patient_id %in% disease_na3$geo_accession)
unique(disease_na$tissue)
unique(disease_na$sample_type)
unique(disco$patient_id[disco$project_id=='GSE184198'])
unique(disco$disease[disco$project_id=='GSE184198'])
disease_na$disease <- disease_na3$`tissue:ch1`[
  match(disease_na$patient_id,disease_na3$geo_accession)]
unique(disease_na$disease)
disco[disease_na.index & disco$project_id=='GSE184198','disease'] <- 'gastric cancer'

gse <- getGEO('GSE188222',destdir=".", AnnotGPL=F, getGPL=F)#
disease_na3 <- pData(gse[[1]])
disease_na <- disco[disease_na.index & disco$project_id=='GSE188222',]
unique(disease_na$patient_id)
disease_na3$geo_accession
table(disease_na$patient_id %in% disease_na3$geo_accession)
unique(disease_na$tissue)
unique(disease_na$sample_type)
unique(disco$patient_id[disco$project_id=='GSE188222'])
unique(disco$disease[disco$project_id=='GSE188222'])
disease_na$disease <- disease_na3$`disease:ch1`[
  match(disease_na$patient_id,disease_na3$geo_accession)]
unique(disease_na$disease)
disco[disease_na.index & disco$project_id=='GSE188222','disease'] <- 'Healthy'
disco[disease_na.index & disco$project_id=='GSE188222','sample_type'] <- 'Healthy'

gse <- getGEO('GSE193138',destdir=".", AnnotGPL=F, getGPL=F)#
disease_na3 <- pData(gse[[1]])
disease_na <- disco[disease_na.index & disco$project_id=='GSE193138',]
unique(disease_na$patient_id)
disease_na3$geo_accession
table(disease_na$patient_id %in% disease_na3$geo_accession)
unique(disease_na$tissue)
unique(disease_na$sample_type)
unique(disco$patient_id[disco$project_id=='GSE193138'])
unique(disco$disease[disco$project_id=='GSE193138'])
disease_na$disease <- disease_na3$`disease state:ch1`[
  match(disease_na$patient_id,disease_na3$geo_accession)]
unique(disease_na$disease)
disco[disease_na.index & disco$project_id=='GSE193138','disease'] <- 'Healthy'
disco[disease_na.index & disco$project_id=='GSE193138','sample_type'] <- 'Healthy'

gse <- getGEO('GSE193138',destdir=".", AnnotGPL=F, getGPL=F)#
disease_na3 <- pData(gse[[1]])
disease_na <- disco[disease_na.index & disco$project_id=='GSE193138',]
unique(disease_na$patient_id)
disease_na3$geo_accession
table(disease_na$patient_id %in% disease_na3$geo_accession)
unique(disease_na$tissue)
unique(disease_na$sample_type)
unique(disco$patient_id[disco$project_id=='GSE193138'])
unique(disco$disease[disco$project_id=='GSE193138'])
disco[disease_na.index & disco$project_id=='GSE193138','disease'] <- 'Healthy'
disco[disease_na.index & disco$project_id=='GSE193138','sample_type'] <- 'Healthy'

gse <- getGEO('GSE202735',destdir=".", AnnotGPL=F, getGPL=F)#
disease_na3 <- pData(gse[[1]])
disease_na <- disco[disease_na.index & disco$project_id=='GSE202735',]
unique(disease_na$patient_id)
disease_na3$geo_accession
table(disease_na$patient_id %in% disease_na3$geo_accession)
unique(disease_na$tissue)
unique(disease_na$sample_type)
unique(disco$patient_id[disco$project_id=='GSE202735'])
unique(disco$disease[disco$project_id=='GSE202735'])
disease_na$disease <- disease_na3$title[
  match(disease_na$patient_id,disease_na3$geo_accession)]
unique(disease_na$disease)
disco[disease_na.index & disco$project_id=='GSE202735','disease'] <- 'mitochondrial encephalopathy, lactic acidosis, and stroke-like episodes'
disco[disease_na.index & disco$project_id=='GSE202735','sample_type'] <- 'MELAS'


disease_na <- disco[disease_na.index & disco$project_id=='HCA_HematopoieticImmuneCellAtlas',]
unique(disease_na$patient_id)
disco[disease_na.index & disco$project_id=='HCA_HematopoieticImmuneCellAtlas','disease'] <- 'Healthy'
disco[disease_na.index & disco$project_id=='HCA_HematopoieticImmuneCellAtlas','sample_type'] <- 'Healthy'
table(is.na(disco$disease))

na.index <- which(is.na(disco$sample_type))
table(disco$tissue[na.index])
table(disco$project_id[na.index])
table(disco$disease[na.index])
table(disco$disease[na.index], disco$ct[na.index])
table(disco$disease[na.index], disco$project_id[na.index])
disco$sample_type[na.index & disco$disease=='COVID-19'] <- 'Infection'
unique(disco$sample_type)
table(is.na(disco$sample_type))

# 3) check the ascp links for diseased samples
sample.index <- c(grep("ormal", disco$sample_type),grep("ealthy", disco$sample_type))
disco_diseased <-  disco[-sample.index,]
table(disco_diseased$sample_type)
disco_meta <- data.frame(table(disco_diseased$sample_type, disco_diseased$disease))
table(is.na(disco_diseased$patient_id))
length(unique(disco_diseased$patient_id)) #1332 diseased samples
disco_sample_summary <- disco_diseased[! duplicated(disco_diseased$patient_id),]
disco_sample_summary$disease <- gsub('inflamed','ulcerative colitis', disco_sample_summary$disease) #
PCTanno_sample_summary <- tissue_summary[! tissue_summary$Disease.Stage %in% c('ADJ','Healthy'),]
length(unique(PCTanno_sample_summary$Patients)) #912
disco_sample_summary$tissue <- gsub('colon', 'colorectum', disco_sample_summary$tissue)
tissue_all <- unique(c(tolower(PCTanno_sample_summary$Tissue), disco_sample_summary$tissue))
length(unique(disco_diseased$project_id))
c(unique(disco_diseased$disease), unique(PCTanno_sample_summary$Disease.Stage))
## first 120 datasets
setwd(data_dir)
disco_GEO <- list.files('./', pattern='.txt')
disco_GEO.acc <- str_split(gsub('.txt','',disco_GEO), '_', simplify = T)[,1]
unique(disco_GEO.acc)
table(unique(disco_GEO.acc) %in% unique(disco_diseased$project_id))
setdiff(unique(disco_GEO.acc), unique(disco_diseased$project_id))

tissue_summary <- read.csv('/data2/rluo4/All/Output/tissue_summary.txt',sep = "\t",header = T)
PCTanno_GEO <- unique(tissue_summary$Dataset)
table(disco_GEO.acc %in% PCTanno_GEO)
PCTanno_disco <- intersect(disco_GEO.acc, PCTanno_GEO) # 9 --> 17 intersected datasets
PCTanno_disco
# View(tissue_summary[tissue_summary$Dataset %in% PCTanno_disco,])
disco_GEO <- disco_GEO[!disco_GEO.acc %in% PCTanno_GEO]#111
disco_GEO.acc <- str_split(gsub('.txt','',disco_GEO), '_', simplify = T)[,1]
table(disco_GEO.acc %in% unique(disco_diseased$project_id))
disco_GEO <- disco_GEO[disco_GEO.acc %in% unique(disco_diseased$project_id)]
# cell_summary <- read.csv('cell_summary.txt',sep = "\t",header = T)

filePath <- sapply(disco_GEO, function(x){ 
  paste(data_dir,x,sep='/')})  
disco_GEO <- lapply(filePath, function(x){
  fread(x)})
disco_allDat <- NULL#data.frame()
unique(disco_diseased$patient_id)
disco_sum <- lapply(disco_GEO,function(y){
  disco_allDat <<- rbind(disco_allDat, y)
  return(disco_allDat)  
})

unique(disco_allDat$sample_alias)
table(disco_allDat$sample_alias %in% disco_diseased$patient_id)
disco_allDat <- disco_allDat[disco_allDat$sample_alias %in% disco_diseased$patient_id,]
unique(disco_allDat$submitted_aspera)
unique(disco_allDat$sample_alias)# 132 --> 413 samples for raw data 
colnames(disco_allDat)
disco_allDat <- disco_allDat[, c('study_alias','sample_alias','fastq_aspera','submitted_aspera','sra_aspera','run_alias')]
disco_allDat$submitted_aspera[is.na(disco_allDat$submitted_aspera)] <- disco_allDat$fastq_aspera[is.na(disco_allDat$submitted_aspera)]
disco_allDat$submitted_aspera
disco_allDat$submitted_aspera[disco_allDat$submitted_aspera==''] <- disco_allDat$fastq_aspera[disco_allDat$submitted_aspera=='']
disco_allDat$submitted_aspera
disco_allDat$datatype[grepl('.bam', disco_allDat$submitted_aspera)] <- 'Bam'
disco_allDat$datatype[grepl('_1.fastq', disco_allDat$submitted_aspera)] <- 'FastQ'
disco_allDat$datatype[is.na(disco_allDat$datatype)] <- 'Sra'
index = disco_allDat$datatype=='Sra'
disco_allDat$submitted_aspera[index] <- disco_allDat$sra_aspera[index]
disco_allDat_GSE <- paste(disco_allDat$study_alias, disco_allDat$study_alias, sep = ';')
disco_allDat_GSE[index] <- disco_allDat$study_alias[index]
disco_allDat_GSE
disco_allDat_ascp <- unlist(str_split(disco_allDat$submitted_aspera,";"))#submitted_aspera#fastq_ftp#
disco_allDat_GSM <- paste(disco_allDat$sample_alias, disco_allDat$sample_alias, sep = ';')
disco_allDat_GSM[index] <- disco_allDat$sample_alias[index]
disco_allDat_GSM
disco_allDat_dt <- paste(disco_allDat$datatype, disco_allDat$datatype, sep = ';')
disco_allDat_dt[index] <- disco_allDat$datatype[index]
disco_allDat_dt
asp_sc <- data.frame(study_alias = unlist(str_split(disco_allDat_GSE,";")), 
                     sample_alias = unlist(str_split(disco_allDat_GSM,";")), 
                     sra_aspera = disco_allDat_ascp, datatype = unlist(str_split(disco_allDat_dt,";"))
)
# write.table(asp_sc,file="/data2/rluo4/EpiTrans/DataCollection/asp_disco_scRNA.link.txt",sep = "\t",quote = F, col.names = F, row.names = F)
# # saveRDS(disco,  file=paste0('/data2/rluo4/RPMfunc/disco_pdata/disco_all_meta.rds'))
################################################################################
# 4) Check the lineage for celltypes 
################################################################################
disco <- readRDS(paste0('/data2/rluo4/RPMfunc/disco_pdata/disco_all_meta.rds'))
length(unique(disco$ct))
Lineage <- data.frame(Major_type = disco$ct, Minor_type = disco$ct, Cell_type = disco$ct)
Lineage <- Lineage[!duplicated(Lineage$Cell_type),]
Lineage <- Lineage[order(Lineage$Cell_type), ]

table_ct <- read_excel('/data2/rluo4/RPMfunc/disco_pdata/table_all_ct.xlsx', sheet = 1)
x <- unlist(table_ct)
str(x)
x <- na.omit(x)
x <- as.character(x)
x = setdiff(Lineage$Cell_type, x)
as.factor(x)
Lineage$Lineage[Lineage$Cell_type %in% table_ct$Embryonic_Cells] <- 'Embryonic_Cells'
Lineage$Lineage[Lineage$Cell_type %in% table_ct$Stem_Cells] <- 'Stem_Cells'
Lineage$Lineage[Lineage$Cell_type %in% table_ct$Germ_Cells] <- 'Germ_Cells'
Lineage$Lineage[Lineage$Cell_type %in% table_ct$Immune_Cells] <- 'Immune_Cells'
Lineage$Lineage[Lineage$Cell_type %in% table_ct$Mesenchymal_Cells] <- 'Mesenchymal_Cells'
Lineage$Lineage[Lineage$Cell_type %in% table_ct$Endothelial_Cells] <- 'Endothelial_Cells'
Lineage$Lineage[Lineage$Cell_type %in% table_ct$Blood_Cells] <- 'Blood_Cells'
Lineage$Lineage[Lineage$Cell_type %in% table_ct$Neuro_Cells] <- 'Neuro_Cells'
Lineage$Lineage[Lineage$Cell_type %in% table_ct$Tissue_Specific_Cells] <- 'Tissue_Specific_Cells'
patterns <- c("lymphoid","gdT","ILC"," B", "B cell","lasma", "T ", " T", "CD4", "CD8", "reg", "LymP", # LymP, lymphoid progenitor
              "yeloid", "onocyte", "acrophage", "eutrophi", "ranulocyte", "Tfh", "NK", "DC", "endritic", "mast","MAIT", "Mast")
# "LZ B cells" refers to B cells that reside in the light zone (LZ) of the germinal centers (GCs) within lymphoid tissues such as lymph nodes and the spleen. In the germinal center, B cells undergo somatic hypermutation and affinity maturation, processes crucial for the generation of high-affinity antibodies during an immune response. The light zone (LZ) is where B cells interact with follicular dendritic cells (FDCs) and T follicular helper (Tfh) cells, facilitating these processes.
# Create regular expression pattern to match any of the patterns
pattern <- paste(patterns, collapse = "|")
# Extract strings matching the pattern
matches <- grep(pattern, Lineage$Cell_type, value = TRUE)
Lineage$Lineage[Lineage$Cell_type %in% matches] <-  'Immune_Cells'

patterns <- c(" EC","ndothelia")
# Create regular expression pattern to match any of the patterns
pattern <- paste(patterns, collapse = "|")
# Extract strings matching the pattern
matches <- grep(pattern, Lineage$Cell_type, value = TRUE)
Lineage$Lineage[Lineage$Cell_type %in% matches] <-  'Endothelial_Cells'

patterns <- c("ibro","muscl", 'peri', 'vascula', 'mura') #, 'HSC' #hepatocellular stella cells
# Create regular expression pattern to match any of the patterns
pattern <- paste(patterns, collapse = "|")
# Extract strings matching the pattern
matches <- grep(pattern, Lineage$Cell_type, value = TRUE)
Lineage$Lineage[Lineage$Cell_type %in% matches] <-  'Mesenchymal_Cells'

patterns <- c("astro","glia", "Glia",'euro')
# Create regular expression pattern to match any of the patterns
pattern <- paste(patterns, collapse = "|")
# Extract strings matching the pattern
matches <- grep(pattern, Lineage$Cell_type, value = TRUE)
Lineage$Lineage[Lineage$Cell_type %in% matches] <-  'Neuro_Cells'

patterns <- c("oocyte",'germ','spermatid','spermatocyte')
# Create regular expression pattern to match any of the patterns
pattern <- paste(patterns, collapse = "|")
# Extract strings matching the pattern
matches <- grep(pattern, Lineage$Cell_type, value = TRUE)
Lineage$Lineage[Lineage$Cell_type %in% matches] <-  'Germ_Cells'

patterns <- c("rythroblast",'Megakaryocyte','erythroid','erythrocyte', 'MEP')#MEP, megakaryocyte progenitor;
# Create regular expression pattern to match any of the patterns
pattern <- paste(patterns, collapse = "|")
# Extract strings matching the pattern
matches <- grep(pattern, Lineage$Cell_type, value = TRUE)
Lineage$Lineage[Lineage$Cell_type %in% matches] <-  'Blood_Cells'

patterns <- c("pithelia","asal","uminal", "duct", "keratino","ndocrine", "ilia",'olonocyte','Leydig','trophoblast','actocyte','ranulosa') 
#cytotrophoblast:  cytotrophoblasts and syncytiotrophoblasts, which are both types of trophoblasts, can be considered epithelial cells. They form the outer layer of the blastocyst during early embryonic development and later contribute to the formation of the placenta.
# The mammalian Leydig cell is a polyhedral epithelioid cell with a single eccentrically located ovoid nucleus. The nucleus contains one to three prominent nucleoli and large amounts of dark-staining peripheral heterochromatin.# Create regular expression pattern to match any of the patterns
# Based on the naming convention "Fibroblast/Granulosa doublet like cell," it seems that these cells are a hybrid or mixture of fibroblasts and granulosa cells. Granulosa cells are a type of somatic cell found in the ovarian follicle, which supports oocyte development.# Given this hybrid nature, it's challenging to definitively categorize them as strictly mesenchymal or tissue-specific. They may exhibit characteristics of both cell types, depending on their specific properties and functions. It's possible that they share features with mesenchymal cells due to the presence of fibroblasts, but they may also have tissue-specific functions related to granulosa cells. Further context or information about their origin and function would be helpful for a more precise classification
pattern <- paste(patterns, collapse = "|")
# Extract strings matching the pattern
matches <- grep(pattern, Lineage$Cell_type, value = TRUE)
Lineage$Lineage[Lineage$Cell_type %in% matches] <-  'Tissue_Specific_Cells'

Lineage$Cell_type[is.na(Lineage$Lineage)]
Lineage$Lineage[is.na(Lineage$Lineage)] <- 'Tissue_Specific_Cells'

table(disco$tissue[disco$ct == 'HSC'] )
Lineage$Lineage[Lineage$Cell_type == 'HSC'] <- 'Stem_Cells'
table(is.na(match(disco$ct, Lineage$Cell_type)))
disco$Major_type <- Lineage$Lineage[match(disco$ct, Lineage$Cell_type)]
saveRDS(disco,  file=paste0('/data2/rluo4/RPMfunc/disco_pdata/disco_all_meta.rds'))

################################################################################
# 2. Arrange the clinical pdata
################################################################################
# sc_majortype <- readRDS(file = '/data2/rluo4/RPMfunc/Output/scRNA/sc_majortype_0725.rds')
# pancreas <- sc_majortype[sc_majortype$tissue == 'pancreas',]
# table( pancreas$sub_type, pancreas$project_id)
# library(DropletUtils)
# pdata.Lin <- read.table('/data2/rluo4/EpiTrans/DataCollection/Lin_pdata.txt', fill = T, sep = '\t', header = T)
# GSE160269 <- sc_majortype[sc_majortype$project_id == 'GSE160269',]
# GSE160269  <- GSE160269 [GSE160269 $disease!='Adjacent tumor/disease section',]
# index <- match(GSE160269$patient_id, paste0(pdata.Lin$Sample.ID, 'T-E'))
# pdata.Lin$patient_id <- paste0(pdata.Lin$Sample.ID, 'T-E')
# GSE160269 <- left_join(GSE160269,pdata.Lin, by = 'patient_id' )
# # GSE160269$stage <- pdata.Lin$Pathologic.stage[index]
# saveRDS(GSE160269, file = '/data2/rluo4/EpiTrans/DataCollection/GSE160269_Lin_meta.rds')
# library(GEOquery)
# metadata <- read.csv('/data2/rluo4/EpiTrans/DataCollection/manifest_20230407_193024.csv')
# HTA12.fq <- metadata[grep('HTA12_',metadata$sample_id),]
# # raw_files <- read.table('/public/home/lorihan/lrh/database/Pancreas/Daniel/raw_files')
# # table(raw_files$V1 %in% HTA12.fq$name)
# # setdiff(raw_files$V1, HTA12.fq$name)
# rawid_Daniel <- HTA12.fq[HTA12.fq$Library.source %in% c('Single Cell','Single Nucleus'),]
# colnames(rawid_Daniel[,c(2,10,18,37)])
# pdata.Daniel <- read.table("/data2/rluo4/EpiTrans/DataCollection/Daniel_pdata.txt",header = F,sep = ';')
# head(pdata.Daniel)
# colnames(pdata.Daniel) <- c(colnames(rawid_Daniel[,c(2,10,18,37)]),'files','sample_name')
# table(pdata.Daniel$Library.source)
# table(pdata.Daniel$files %in% rawid_Daniel$name)#Pancreatobiliary-type carcinoma
# colnames(rawid_Daniel)[2] <- 'files'
# dim(rawid_Daniel)
# pdata.Daniel <- left_join(pdata.Daniel[,-4], rawid_Daniel[, c(2,20:65)], by = 'files')
# HTA12 <- sc_majortype[sc_majortype$project_id == 'HTA12',]
# table(unique(HTA12$patient_id) %in% pdata.Daniel$sample_id)
# HTA12$grade <- pdata.Daniel$Tumor.grade[match(HTA12$patient_id, pdata.Daniel$sample_id)]
# HTA12$grade[HTA12$grade=='Not Applicable'] <- 'G1'
# table(HTA12$grade, HTA12$ct)
# table(pdata.Daniel$Tumor.tissue.type)
# pdata.Daniel$Disease.State <- gsub('Not reported','Unknown tumor status',pdata.Daniel$Last.known.disease.status)
# pdata.Daniel$Disease.State[grepl('progression', pdata.Daniel$Disease.State )] <- 'Recurrence'
# pdata.Daniel$Site[!grepl('pancreas', pdata.Daniel$Site )] <- 'Unknown'
# table(pdata.Daniel$Disease.State)
# table(pdata.Daniel$experimental_strategy)
# # pdata.Daniel$Tumor.tissue.type[pdata.Daniel$Tumor.tissue.type != 'Primary'] <- 'NonPrimary'
# table(pdata.Daniel$Tumor.tissue.type, pdata.Daniel$Library.source)
# HTA12$Tumor.tissue.type <- pdata.Daniel$Disease.State[match(HTA12$patient_id, pdata.Daniel$sample_id)]
# HTA12$Library.source <- pdata.Daniel$Library.source[match(HTA12$patient_id, pdata.Daniel$sample_id)]
# HTA12$Site <- pdata.Daniel$Site[match(HTA12$patient_id, pdata.Daniel$sample_id)]
# table(HTA12$Site)
# table(HTA12$Library.source)
# HTA12$class <- paste(HTA12$Library.source, HTA12$grade, HTA12$Tumor.tissue.type)
# HTA12$class <- paste(HTA12$Library.source, HTA12$grade)#, HTA12$Site)
# table(HTA12$class, HTA12$ct)
# table(HTA12$Tumor.tissue.type, HTA12$ct)
# saveRDS(HTA12, file = '/data2/rluo4/EpiTrans/DataCollection/HTA12_Daniel_meta.rds')
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
# 3. Organize the cell annotation data
################################################################################
# # Load Diseased Samples
# sc_majortype <- readRDS(file = '/data2/rluo4/RPMfunc/Output/scRNA/sc_majortype_0725.rds') # 340 --> 312 datasets
# sc_summary <- sc_majortype[! duplicated(paste(sc_majortype$patient_id, sc_majortype$tissue, 
#                                               sc_majortype$project_id)),]
# length(unique(sc_summary$project_id)) #396 --> 340
# sc_summary$omics <- 'scRNA'
# sample.index <- c(grep("ormal", sc_summary$sample_type),grep("ealthy", sc_summary$sample_type),
#                   grep("control", sc_summary$sample_type))
# sc_Diseased_Samples <- sc_summary[ -sample.index, c(2,3,5,6)]
# sample.index <- c(grep("ormal", sc_majortype$sample_type),grep("ealthy", sc_majortype$sample_type),
#                   grep("Norm", sc_majortype$sample_type), grep("control", sc_majortype$sample_type))
# sc_cell_summary <- sc_majortype[ -sample.index,]
# sc_cell_summary <- sc_cell_summary[ ! sc_cell_summary$project_id %in% tissue_summary$Dataset, ]
# first_elements <- extract_first_element(sc_cell_summary$barcode)
# CB = str_split(first_elements,'[-]',simplify = T)[,1]
# sc_cell_summary$CB <- CB#str_split(sc_cell_summary$barcode, '-', simplify = T)[,1]
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
# disco_pheno <- as.data.frame(table(sc_cell_summary$sample_type, sc_cell_summary$disease))
# disco_pheno <- disco_pheno[disco_pheno$Freq!=0,]
# # sample_type inconsistent with disease: GSE168453 (Critical/Moderate -- Healthy)
# sample.index <- c(grep("ormal", disco_pheno$Var1),grep("ealthy", disco_pheno$Var1),
#                   grep("Norm", disco_pheno$Var1), grep("control", disco_pheno$Var1))
# # View(disco_pheno[sample.index,])
################################################################################
# 4. Find the missed BM&MAST datasets
################################################################################
ct_all <- unique(sc_cell_summary$ct)
# Filter for project_id values with ct as 'Mast cell' or 'Basophil'
mastcell_basophil_projects <- sc_cell_summary %>%
  filter(ct %in% c("Mast cell", "Mast cells", "Basophil")) %>%
  distinct(project_id)
unique(mastcell_basophil_projects$project_id)
non_mastcell_basophil_projects <- sc_cell_summary %>%
  # filter(! ct %in% c("Mast cell", "Mast cells", "Basophil")) 
  filter(ct %in% setdiff(ct_all,  c("Mast cell", "Mast cells", "Basophil")) ) %>%
  distinct(project_id)
exclusive_mastcell_basophil_projects <- setdiff(unique(mastcell_basophil_projects$project_id), 
                                                unique(non_mastcell_basophil_projects$project_id))
sort(exclusive_mastcell_basophil_projects)
exclusive_mastcell_basophil_projects <- sc_cell_summary %>%
  filter(project_id %in% exclusive_mastcell_basophil_projects) %>%
  distinct(project_id)
# exclusive_mastcell_basophil_projects <- sc_cell_summary %>%
#   filter( ! ct %in% setdiff(ct_all,  c("Mast cell", "Mast cells", "Basophil")) ) %>%
#   distinct(project_id)
print(exclusive_mastcell_basophil_projects)
index = match(exclusive_mastcell_basophil_projects$project_id, sc_cell_summary$project_id)
exclusive_mastcell_basophil_projects$tissue <- sc_cell_summary$tissue[index]
exclusive_mastcell_basophil_projects$disease <- sc_cell_summary$disease[index]

disco <- readRDS(paste0('/data2/rluo4/RPMfunc/disco_pdata/disco_all_meta.rds')) # made by Ruihan + Haixia
BM <- disco[disco$tissue == 'bone marrow',]
pancreas <- disco[disco$tissue == 'pancreas',]
leukemia <- disco[grepl('eukemia', disco$disease),]
unique(leukemia$project_id)
#[1] "GSE211036" "GSE185991" "GSE211033"
#[4] "GSE185381"

# 4.1 disco changes on AML
GSE185991 <- readRDS('/data2/rluo4/EpiTrans/DataCollection/leukemia/all_meta_GSE185991.rds')#sub_NA_GSE185991_metadata.rds')
# GSE185991 <- readRDS('/data2/rluo4/EpiTrans/DataCollection/leukemia/add_GSE185991_meta.rds')
# GSE185991$BC <- paste(GSE185991$orig.ident, rownames(GSE185991), sep = '_')
GSE185991$sub_type <- as.character(GSE185991$ct)
GSE185991$sub_type <- gsub('Adipocytes','OtherMesenchymals',GSE185991$sub_type)
GSE185991$sub_type <- gsub('B-cells','Bcells',GSE185991$sub_type)
GSE185991$sub_type[grepl('CD4|CD8',GSE185991$sub_type)] <- 'Tcells'
GSE185991$sub_type <- gsub('NK cells','NKcells',GSE185991$sub_type)
GSE185991$sub_type <- gsub('HSC','Stem_Cells',GSE185991$sub_type)
table(GSE185991$sub_type)
# GSE185991 <- rbind(GSE185991, GSE185991[,colnames(GSE185991)])
# HSC – hematopoietic stem cells, MPP – multipotent progenitors, GMP – granulocyte-monocyte progenitors, MEP – megakaryocyte progenitors, LymP – lymphoid progenitors, DC – dendritic cells, Ery – erythrocytes
GSE185381 <- readRDS('/data2/rluo4/EpiTrans/DataCollection/leukemia/sub_NA_GSE185381_metadata.rds')
GSE185381 <- readRDS('/data2/rluo4/EpiTrans/DataCollection/GSE2185381_merged_data.rds')
GSE185381$ct <- GSE185381$Identity
GSE185381$ct[is.na(GSE185381$ct)] <- GSE185381$Broad_cell_identity[is.na(GSE185381$ct)]
table(GSE185381$ct %in% sc_majortype$ct)
GSE185381$sub_type <- sc_majortype$sub_type[match(GSE185381$ct, sc_majortype$ct)]
table(GSE185381$ct[is.na(GSE185381$sub_type)])
GSE185381$sub_type[is.na(GSE185381$sub_type)] <- 'Tcells'
table(GSE185381$sub_type)
library(Matrix, lib.loc = "/data2/rluo4/bin/miniconda3/lib/R/library")
library(Seurat)
library(readr)
library(hdf5r)
GSE211036 <- Read10X_h5("/data2/rluo4/EpiTrans/DataCollection/leukemia/GSM6447735_ID216_filtered_feature_bc_matrix.h5")
dim(GSE211036)
cell_type_df <- read_delim('/data2/rluo4/EpiTrans/DataCollection/leukemia/GSM6447735_ID216_scRNAseq_Barcode_cluster_assignments.tsv.gz')
rownames(cell_type_df) <- cell_type_df$Cell_Barcode
colnames(cell_type_df) <- c('barcode', 'ClusterID', 'Cell_Type')
GSE211036 <- CreateSeuratObject(GSE211036, project = "GSM6447735")#, meta.data = cell_type_df) #后面就可以单细胞处理的标准流程啦
GSE211036$barcode <- colnames(GSE211036)
library(dplyr)
GSE211036 <- GSE211036[, GSE211036$barcode %in% cell_type_df$barcode]
dim(GSE211036)
GSE211036@meta.data <- left_join(GSE211036@meta.data, cell_type_df, by = 'barcode')
# GSE211036 <- AddMetaData(GSE211036, cell_type_df$CellType, col.name = "Cell_Type")
table(GSE211036$Cell_Type)
GSE211036$sub_type <- gsub('LSC', 'Stem_Cells', gsub('Blast', 'Blood_Cells', GSE211036$Cell_Type))
GSE211036$ct <- gsub('LSC', 'Leukemic stem cells', gsub('Blast', 'Leukemic blast cells', GSE211036$Cell_Type))

library(stringr)
GSE185381$sample_id <- str_split(rownames(GSE185381), '_', simplify = T)[,1]
index1 <- grepl('CD11c+|monocyte|DC|LymP|GMP',GSE185381$ct)
GSE185381$sub_type[index1] <- 'OtherImmunecells'
index2 <- grepl('egakaryocyte|Ery|MEP',GSE185381$ct)
GSE185381$sub_type[index2] <- 'Blood_Cells'
CT <- as.data.frame(table(GSE185381$ct, GSE185381$sub_type))
CT <- CT[CT$Freq!=0,]

GSE185991$sample_id <- str_split(GSE185991$orig.ident, '_', simplify = T)[,1]#str_split(rownames(GSE185991), '_', simplify = T)[,1]
GSE185991$barcode <- rownames(GSE185991)#str_split(rownames(GSE185991), 'rds_', simplify = T)[,2]
index1 <- grepl('DC|Monocytes|Eosinophils|Neutrophils',GSE185991$sub_type)
GSE185991$sub_type[index1] <- 'OtherImmunecells'
index2 <- grepl('Erythrocytes',GSE185991$ct)
GSE185991$sub_type[index2] <- 'Blood_Cells'
table(GSE185991$sub_type)
CT <- as.data.frame(table(GSE185991$ct, GSE185991$sub_type))
CT <- CT[CT$Freq!=0,]
AML1 <- data.frame(project_id = 'GSE185381', patient_id = GSE185381$patient_id, ct = GSE185381$ct,
                  tissue = 'bone marrow', platform = '10x3', sample_type = GSE185381$ap_aml_age,
                  disease = GSE185381$ap_aml_age, sub_type = GSE185381$sub_type)
AML1$disease <- ifelse(grepl('AML', AML1$disease), 'AML', 'healthy')
rownames(AML1) <- rownames(GSE185381) #paste(GSE185381$barcode, GSE185381$sample_id, GSE185381$ap_aml_age, sep = '--') 
saveRDS(AML1,  file=paste0('/data2/rluo4/RPMfunc/disco_pdata/GSE185381_meta.rds'))

AML2 <- data.frame(project_id = 'GSE185991', patient_id = GSE185991$sample_id, ct = GSE185991$ct,
                  tissue = 'bone marrow', platform = '10x3', sample_type = GSE185991$Timepoint,
                  disease = 'AML', sub_type = GSE185991$sub_type)
rownames(AML2) <- paste(GSE185991$barcode, GSE185991$sample_id, sep = '--')

AML3 <- data.frame(project_id = 'GSE211036', patient_id = GSE211036$orig.ident, ct = GSE211036$ct,
                  tissue = 'bone marrow', platform = '10x3', sample_type = 'AML',
                  disease = 'AML', sub_type = GSE211036$sub_type)
rownames(AML3) <- paste(GSE211036$barcode, GSE211036$orig.ident, sep = '--')

AML <- rbind(AML1, AML2, AML3)
# saveRDS(AML,  file=paste0('/data2/rluo4/RPMfunc/disco_pdata/AML_meta.rds'))

disco <- disco[!disco$project_id %in% BM$project_id,] #
table(AML$sub_type)
colnames(AML)[8] <- 'Major_type'
index1 <- grepl('Bcells|Macrophages|Neutrophils|NKcells|OtherImmunecells|Tcells',AML$Major_type)
AML$Major_type[index1] <- 'Immune_Cells'
index2 <- grepl('OtherMesenchymals',AML$Major_type)
AML$Major_type[index2] <- 'Mesenchymal_Cells'
disco <- rbind(disco, AML)

# 4.2 disco changes on pancreas
GSE205049 <- readRDS('/data2/rluo4/EpiTrans/DataCollection/PDAC/GSE205049_meta.rds')
GSE205049$sub_type <- as.character(GSE205049$type)
GSE205049$sub_type <- gsub('B cells','Bcells',GSE205049$sub_type)
GSE205049$sub_type[grepl('CD4|CD8| T',GSE205049$sub_type)] <- 'Tcells'
GSE205049$sub_type[grepl('M2',GSE205049$sub_type)] <- 'Macrophages'
GSE205049$sub_type[grepl('NK',GSE205049$sub_type)] <- 'NKcells'
GSE205049$sub_type <- gsub('Plasma cells','Bcells',GSE205049$sub_type)
GSE205049$sub_type <- gsub('Mast cells','MastCell',GSE205049$sub_type)
GSE205049$sub_type[grepl('onocyte|DC',GSE205049$sub_type)] <- 'OtherImmunecells'
table(GSE205049$sub_type)
# GSE205049$sample_type <-  gsub('PDAC','Primary Tumor',
#                                gsub('Adjacent normal','AdjNorm',GSE205049$DiseaseState) )
GSE205049$sample_type <-  gsub('PDAC','Primary Tumor',
                               gsub('AdjNorm','Adjacent normal',GSE205049$DiseaseState) )
PDAC1 <- data.frame(project_id = 'GSE205049', patient_id = GSE205049$orig.ident, ct = GSE205049$type,
                   tissue = 'pancreas', platform = '10x3', sample_type = GSE205049$sample_type, #'Primary Tumor'
                   disease = GSE205049$DiseaseState, sub_type = GSE205049$sub_type)
rownames(PDAC1) <- paste(GSE205049$barcode, GSE205049$gsm_name, sep = '--')
# saveRDS(PDAC1,  file=paste0('/data2/rluo4/RPMfunc/disco_pdata/PDAC1_meta.rds'))

disco <- disco[!disco$project_id %in% c('GSE205049','GSE217837','GSE231535'),]#pancreas$project_id,]
colnames(PDAC1)[8] <- 'Major_type'
index1 <- grepl('Bcells|Macrophages|Neutrophils|NKcells|OtherImmunecells|Tcells|MastCell',PDAC1$Major_type)
PDAC1$Major_type[index1] <- 'Immune_Cells'
disco <- rbind(disco, PDAC1)
# saveRDS(disco,  file=paste0('/data2/rluo4/RPMfunc/disco_pdata/disco_adjustmeta_0725.rds')) #AdjNorm for PDAC1
# saveRDS(GSE211036,  file=paste0('/data2/rluo4/RPMfunc/disco_pdata/GSE211036.rds'))
disco <- disco[! duplicated(paste(disco$project_id, disco$patient_id)),]
saveRDS(disco,  file=paste0('/data2/rluo4/RPMfunc/disco_pdata/disco_adjustmeta_1122.rds')) #AdjNorm for PDAC1

# 4.3 disco changes on bone marrow/mast
# adjust the disco meta info for SCOMATIC input:
disco <- readRDS(paste0('/data2/rluo4/RPMfunc/disco_pdata/disco_adjustmeta_0725.rds'))
sc_majortype <- readRDS(file = '/data2/rluo4/RPMfunc/Output/scRNA/sc_majortype_0725.rds') # 340 --> 312 datasets made from Summary_RMzyme.R
error_datasets <- read.csv('/data2/rluo4/EpiTrans/DataCollection/error_ct.csv')
# BM_MAST39 <- readRDS('/data2/rluo4/EpiTrans/DataCollection/all_table_of_disco_39.rds')
# BM_MAST39 <- BM_MAST39[!BM_MAST39$Project.ID %in% c('dgcMatrix(.rds)', '10Xh5(.h5)', 'Cell type(.txt)','Project ID', ''),]
# colnames(BM_MAST39) <- c('PatientID', colnames(BM_MAST39)[1:(ncol(BM_MAST39)-1)])
# setdiff(unique(BM_MAST39$Project.ID), unique(BM_MAST_meta_p1$project_id)) #11
# BM_MAST39 <- data.frame(project_id = BM_MAST39$Project.ID, patient_id = BM_MAST39$PatientID, ct = 'UnknownYet',
#                       tissue = BM_MAST39$Tissue, platform = BM_MAST39$Platform, sample_type = BM_MAST39$Sample.type, 
#                       disease = BM_MAST39$Disease, sub_type = 'UnknownYet')
# Below is the newly added BM & MAST datasets:
# all:
disco_BM_MAST <- readRDS('/data2/rluo4/EpiTrans/DataCollection/all_table_of_disco.rds') # 72 new datasets
disco_BM_MAST <- disco_BM_MAST[!disco_BM_MAST$Project.ID %in% c('dgcMatrix(.rds)', '10Xh5(.h5)', 'Cell type(.txt)','Project ID', ''),]
colnames(disco_BM_MAST) <- c('PatientID', colnames(disco_BM_MAST)[1:(ncol(disco_BM_MAST)-1)])
# # rownames(disco_BM_MAST) <- paste(disco_BM_MAST$barcode, disco_BM_MAST$gsm_name, sep = '--')
disco_BM_MAST <- data.frame(project_id = disco_BM_MAST$Project.ID, patient_id = disco_BM_MAST$PatientID, ct = 'UnknownYet',
                            tissue = disco_BM_MAST$Tissue, platform = disco_BM_MAST$Platform, sample_type = disco_BM_MAST$Sample.type, 
                            disease = disco_BM_MAST$Disease, Major_type = 'UnknownYet')
unique(disco_BM_MAST$project_id)
table(duplicated(disco_BM_MAST$patient_id))
table(is.na(disco_BM_MAST$tissue))
table(is.na(disco_BM_MAST$disease))
table(disco_BM_MAST$disease=='')
disco_BM_MAST[(disco_BM_MAST$disease==''),'disease'] <- 'Healthy'
# View(disco_BM_MAST[duplicated(disco_BM_MAST$patient_id),])
# write.table(disco_BM_MAST, file = '/data2/rluo4/EpiTrans/DataCollection/disco_BM_MAST.txt', quote = F)
# part1:
BM_MAST_meta_p1 <- readRDS('/data2/rluo4/EpiTrans/DataCollection/add_disco_part1_40datasets.rds')#add_disco_part1.rds')#BM_merge_meta.rds') # GSE214693 CITE-seq (2T + 2N)
table(is.na(BM_MAST_meta_p1$ct))
BM_MAST_meta_p1 <- BM_MAST_meta_p1[!is.na(BM_MAST_meta_p1$ct),]
table(is.na(BM_MAST_meta_p1$tissue))
index = is.na(BM_MAST_meta_p1$tissue)
table(BM_MAST_meta_p1$project_id[index])
GSE149512_meta <- disco_BM_MAST[disco_BM_MAST$project_id=='GSE149512',]
GSE149512_meta <- GSE149512_meta[, ! colnames(GSE149512_meta) %in% c('ct','Major_type')]
GSE149512 <- BM_MAST_meta_p1[index,]
GSE149512 <- GSE149512[, ! colnames(GSE149512) %in% colnames(GSE149512_meta)]
table(GSE149512$sample_id %in% GSE149512_meta$patient_id)
colnames(GSE149512_meta) <- gsub('patient_id', 'sample_id', colnames(GSE149512_meta))
GSE149512 <- left_join(GSE149512, GSE149512_meta, by = 'sample_id')
rownames(GSE149512) <- GSE149512$barcode
BM_MAST_meta_p1 <- BM_MAST_meta_p1[! index, ]
index = colnames(GSE149512) %in% colnames(BM_MAST_meta_p1)
BM_MAST_meta_p1 <- rbind(BM_MAST_meta_p1, GSE149512[, index ])
table(BM_MAST_meta_p1$sample_type, BM_MAST_meta_p1$disease)
sample.index <- grepl("ormal|ontrol",BM_MAST_meta_p1$sample_type) 
table(BM_MAST_meta_p1$disease[sample.index])
BM_MAST_meta_p1$disease[sample.index] <- 'AdjNorm'
CT <- as.data.frame(table(BM_MAST_meta_p1$ct, BM_MAST_meta_p1$sub_type))
CT <- CT[CT$Freq!=0,] #As resident macrophages of the central nervous system (CNS), microglia are associated with diverse functions essential to the developing and adult brain during homeostasis and disease. They are aided in their tasks by intricate bidirectional communication with other brain cells under steady-state conditions as well as with infiltrating peripheral immune cells during perturbations. Harmonious cell-cell communication involving microglia are considered crucial to maintain the healthy state of the tissue environment and to overcome pathology such as neuroinflammation. Analyses of such intercellular pathways have contributed to our understanding of the heterogeneous but context-associated microglial responses to environmental cues across neuropathology, including inflammatory conditions such as infections and autoimmunity, as well as immunosuppressive states as seen in brain tumors. Here, we summarize the latest evidence demonstrating how these interactions drive microglia immune and non-immune functions, which coordinate the transition from homeostatic to disease-related cellular states.
table(unique(BM_MAST_meta_p1$project_id) %in% c('GSE211033','GSE196052','GSE196676'))

# part2:
BM_MAST_meta_p2 <- readRDS('/data2/rluo4/EpiTrans/DataCollection/add_disco_part2_30datasets.rds')
table(is.na(BM_MAST_meta_p2$ct))
# BM_MAST_meta_p2 <- BM_MAST_meta_p2[!is.na(BM_MAST_meta_p2$ct),]
table(is.na(BM_MAST_meta_p2$tissue))
index = is.na(BM_MAST_meta_p2$tissue)
table(BM_MAST_meta_p2$project_id[index])
GSE222078_meta <- disco_BM_MAST[disco_BM_MAST$project_id=='GSE222078',]
GSE222078_meta <- GSE222078_meta[, ! colnames(GSE222078_meta) %in% c('ct','Major_type')]
GSE222078 <- BM_MAST_meta_p2[index,]
GSE222078 <- GSE222078[, ! colnames(GSE222078) %in% colnames(GSE222078_meta)]
table(GSE222078$sample_id %in% GSE222078_meta$patient_id)
colnames(GSE222078_meta) <- gsub('patient_id', 'sample_id', colnames(GSE222078_meta))
GSE222078 <- left_join(GSE222078, GSE222078_meta, by = 'sample_id')
BM_MAST_meta_p2 <- BM_MAST_meta_p2[! index, ]
index = colnames(GSE222078) %in% colnames(BM_MAST_meta_p2)
BM_MAST_meta_p2 <- rbind(BM_MAST_meta_p2, GSE222078[, index ])
table(BM_MAST_meta_p2$sample_type, BM_MAST_meta_p2$disease)
table(BM_MAST_meta_p2$disease=='')
index = BM_MAST_meta_p2$disease==''
table(BM_MAST_meta_p2$project_id[index]) #E-MTAB-8142
table(BM_MAST_meta_p2$sample_type[index] ) #control
BM_MAST_meta_p2[index,'disease'] <- 'Healthy'
sample.index <- grepl("ormal|ontrol",BM_MAST_meta_p2$sample_type) 
table(BM_MAST_meta_p2$disease[sample.index])
BM_MAST_meta_p2$disease[sample.index] <- 'AdjNorm'

# part2 -- missed 1 + 2:
# first:
# setwd('/data2/rluo4/RPMfunc/disco_pdata')
# gse <- getGEO('GSE181294',destdir=".", AnnotGPL=F, getGPL=F)
# GSE181294_pdata <- pData(gse[[1]])
# saveRDS(GSE181294_pdata, file = '/data2/rluo4/EpiTrans/DataCollection/GSE181294_pdata.rds')#,sep = '\t', quote = F)
GSE181294_pdata <- readRDS('/data2/rluo4/EpiTrans/DataCollection/GSE181294_pdata.rds')#read.table('/data2/rluo4/EpiTrans/DataCollection/GSE181294_pdata.txt',sep = '\t',header = T, fill = T)
colnames(GSE181294_pdata)
GSE181294_pdata <- data.frame(project_id = "GSE181294", patient_id = GSE181294_pdata$title, 
                              ct = "UnknownYet", tissue = 'prostate', platform = "10x3'", 
                              sample_type = GSE181294_pdata[,'grade:ch1'], 
                              disease = GSE181294_pdata[,'condition:ch1'], Major_type = "UnknownYet")
GSE181294_pdata$sample_type[GSE181294_pdata$sample_type=='NA'] <- 'normal'
GSE181294_pdata$disease <- gsub('Tumor', 'Prostate cancer', GSE181294_pdata$disease)
GSE181294_pdata$sample_type[GSE181294_pdata$disease=='Normal'] <- 'normal'
GSE181294_meta <- disco_BM_MAST[disco_BM_MAST$project_id == 'GSE181294',]
# table(unique(GSE181294_pdata$patient_id) %in% unique(GSE181294_meta$patient_id))
# View(GSE181294_pdata[ ! GSE181294_pdata$sample_id %in% unique(GSE181294_meta$patient_id),])
# disco_BM_MAST <- disco_BM_MAST[ ! disco_BM_MAST$project_id %in% c('PRJNA434002'),] #PRJNA434002 is Single-nucleus RNA sequencing of post-mortem brain tissue from Autism Spectrum Disorder patients
table(unique(disco_BM_MAST$project_id) %in% c('GSE211033','GSE196052','GSE196676'))

# disco_BM_MAST$project_id <- gsub('PRJNA434002', 'GSE185965', disco_BM_MAST$project_id)
missed_BM_MAST <- readRDS('/data2/rluo4/EpiTrans/DataCollection/GSE965GSE677_table_of_disco.rds')
missed_BM_MAST <- missed_BM_MAST[, 2:7] # should be: GSE185965, GSE197677
missed_BM_MAST_1 <- BM_MAST_meta_p1 %>% filter(project_id=='GSE185965') %>% 
  select(project_id, patient_id = sample_id, ct, tissue, platform, 
         sample_type, disease) #,Major_type
missed_BM_MAST_1 <- missed_BM_MAST_1[!duplicated(missed_BM_MAST_1$patient_id),]
missed_BM_MAST_2 <- BM_MAST_meta_p2 %>% filter(project_id=='GSE197677') %>% 
  select(project_id, patient_id = sample_id, ct, tissue, platform, 
         sample_type, disease)
missed_BM_MAST_2 <- missed_BM_MAST_2[!duplicated(missed_BM_MAST_2$patient_id),]

# dt 2:
missed_data2 <- readRDS('/data2/rluo4/EpiTrans/DataCollection/add_70_to_72_meta.rds')
missed_BM_MAST_3 <-  missed_data2  %>% 
  select(project_id, patient_id = sample_id, ct, tissue, platform, 
         sample_type, disease) #,Major_type
missed_BM_MAST_3 <- missed_BM_MAST_3[!duplicated(missed_BM_MAST_3$patient_id),]
# dt 1:
missed_data1 <- readRDS('/data2/rluo4/EpiTrans/DataCollection/bile_duct_healthy_meta.rds')
missed_BM_MAST_4 <-  missed_data1  %>% 
  select(project_id, patient_id = sample_id, ct, tissue, platform, 
         sample_type, disease) #,Major_type
missed_BM_MAST_4 <- missed_BM_MAST_4[!duplicated(missed_BM_MAST_4$patient_id),]
missed_BM_MAST <- rbind(missed_BM_MAST_1, missed_BM_MAST_2, missed_BM_MAST_3, missed_BM_MAST_4)
missed_BM_MAST$Major_type <- 'UnknownYet'

table(disco_BM_MAST$project_id %in% missed_BM_MAST$project_id)
disco_BM_MAST <- rbind(disco_BM_MAST, missed_BM_MAST)
unique(disco_BM_MAST$project_id) # 74 --> 76 -->77
disco_BM_MAST <- disco_BM_MAST[disco_BM_MAST$project_id != 'GSE181294',]
disco_BM_MAST <- rbind(disco_BM_MAST, GSE181294_pdata)
unique(disco_BM_MAST$project_id) # 74 --> 76 -->77
# saveRDS(disco_BM_MAST,  file=paste0('/data2/rluo4/RPMfunc/disco_pdata/BM_MAST_meta.rds')) # 76-->77 new datasets
table(unique(disco_BM_MAST$project_id) %in% c('GSE211033','GSE196052','GSE196676'))

# wrong 3:
wrong_3 <- readRDS('/data2/rluo4/EpiTrans/DataCollection/add_disco_part1.rds')#BM_merge_meta.rds') # GSE214693 CITE-seq (2T + 2N)
table(unique(wrong_3$project_id) %in% c('GSE211033','GSE196052','GSE196676'))
wrong_3 <- subset(wrong_3, project_id %in% c('GSE211033','GSE196052','GSE196676'))

#
GSE181294_meta <- disco_BM_MAST[disco_BM_MAST$project_id=='GSE181294',]
GSE181294_meta <- GSE181294_meta[, ! colnames(GSE181294_meta) %in% c('ct','Major_type')]
colnames(GSE181294_meta) <- gsub('patient_id', 'sample_id', colnames(GSE181294_meta))
GSE181294 <- read.csv('/data2/rluo4/EpiTrans/DataCollection/GSE181294_scRNAseq.ano.csv')
colnames(GSE181294) <- c('barcode','ct','sample_id')
GSE181294 <- left_join(GSE181294, GSE181294_meta, by = 'sample_id')
index <-  ! colnames(BM_MAST_meta_p2) %in% colnames(GSE181294)
GSE181294_fake <- BM_MAST_meta_p2[1:nrow(GSE181294), index]
rownames(GSE181294_fake) <- paste(str_split(GSE181294$barcode,'_',simplify = T)[,2],
                                  GSE181294$sample_id, sep = '--')
# GSE181294_fake$sample_id <- GSE181294$sample_id
GSE181294_fake$orig.ident <- GSE181294$sample_id
GSE181294_fake <- cbind(GSE181294_fake, GSE181294)#left_join(GSE181294_fake, GSE181294, by = 'sample_id')
BM_MAST_meta <- unique(c(BM_MAST_meta_p1$project_id, BM_MAST_meta_p2$project_id)) #69
setdiff(unique(disco_BM_MAST$project_id), BM_MAST_meta) # "GSE181294" 
# [1] "GSE211033"   "GSE196676"   "GSE196052"   "PRJNA434002" "GSE181294"  
index = match(colnames(BM_MAST_meta_p2), colnames(GSE181294_fake))
BM_MAST_meta_p2 <- rbind(BM_MAST_meta_p2, GSE181294_fake[, index ])
# setdiff(unique(disco_BM_MAST$project_id), unique(BM_MAST_meta_p1$project_id)) 
# setdiff(unique(disco_BM_MAST$project_id), unique(BM_MAST_meta_p2$project_id)) 
inter_sets <- sort(intersect(disco_BM_MAST$project_id, sc_majortype$project_id))#17 --> 55-->57
inter_sets
# [1] "GSE130560" "GSE149512" "GSE150825" "GSE151192" "GSE152042"
# [6] "GSE155468" "GSE156285" "GSE156625" "GSE157703" "GSE159677"
# [11] "GSE166352" "GSE173706" "GSE173896" "GSE179159" "GSE181294"
# [16] "GSE181688" "GSE184198"
inter_majortype <- sc_majortype[sc_majortype$project_id %in% disco_BM_MAST$project_id,]
table(inter_majortype$disease)
unique(inter_majortype$ct[! inter_majortype$disease %in% 'healthy']) # [1] "Mast cell" "Basophil" 
# inter_summary <- sc_cell_summary[sc_cell_summary$project_id %in% disco_BM_MAST$project_id,]
table(disco$project_id %in% disco_BM_MAST$project_id)
################################################################################
disco <- disco[!disco$project_id %in% disco_BM_MAST$project_id,] #disco_BM_MAST72
disco <- rbind(disco, disco_BM_MAST)
disco_idents <- paste(disco$project_id, disco$patient_id)
table(duplicated(disco_idents))
unique(disco$project_id[!duplicated(disco_idents)])
disco <- disco[!duplicated(disco_idents),]
# saveRDS(disco,  file=paste0('/data2/rluo4/RPMfunc/disco_pdata/disco_adjustmeta_0827.rds'))#0821.rds')) # added 76-->77 renew disco_BM_MAST datasets
# saveRDS(disco,  file=paste0('/data2/rluo4/RPMfunc/disco_pdata/disco_adjustmeta_0809.rds')) # added 74 new disco_BM_MAST datasets
# 0804.rds'))
unique(disco$project_id) #340
# setdiff(unique(BM_MAST39$project_id), unique(BM_MAST_meta_p1$project_id)) 
# intersect(unique(BM_MAST39$project_id), unique(BM_MAST_meta_p1$project_id)) 
################################################################################
# 5. Adjust disco's sc_majortype:
################################################################################
disco <- readRDS(paste0('/data2/rluo4/RPMfunc/disco_pdata/disco_adjustmeta_0827.rds')) # 340 datasets
disco_majortype <- readRDS('/data2/rluo4/RPMfunc/Output/scRNA/disco_all_adjustmeta_0515.rds') # made by Haixia
################################################################################
sc_majortype <- readRDS(file = '/data2/rluo4/RPMfunc/Output/scRNA/sc_majortype_0828.rds')#08.rds')
sc_majortype$Major_type[sc_majortype$Major_type=='Epithelial'] <- 'Tissue_Specific_Cells'
sc_majortype$Major_type[sc_majortype$Major_type=='Stromal'] <- 'Mesenchymal_Cells'
sc_majortype$Major_type[sc_majortype$Major_type=='Immune'] <- 'Immune_Cells'
sc_majortype[sc_majortype$ct=='HSC','Major_type'] <- 'Stem_Cells'

index_ct_subtype  <- sc_majortype %>% select(ct, sub_type, Major_type, tissue, disease)
index_ct_subtype <- index_ct_subtype[!duplicated(index_ct_subtype$ct),]
table(is.na(match( BM_MAST_meta_p2$ct, index_ct_subtype$ct )))

index <- match( BM_MAST_meta_p2$ct, index_ct_subtype$ct )
BM_MAST_meta_p2$sub_type <- index_ct_subtype$sub_type[index]
index <- is.na(index)
table(BM_MAST_meta_p2$ct[index])
subtype  <- BM_MAST_meta_p2[index, c('ct','sub_type')]
patterns <- c(" EC","ndothelia")
pattern <- paste(patterns, collapse = "|")
subtype[grepl(pattern,subtype$ct), 'sub_type'] <- 'Endothelial_Cells'
patterns <- c("muscl", 'peri', 'vascula', 'mura', 'stromal')
pattern <- paste(patterns, collapse = "|")
table(subtype[grepl(pattern,subtype$ct),'ct'])
subtype[grepl(pattern,subtype$ct), 'sub_type'] <-  'OtherMesenchymals'
subtype[grepl('ibro',subtype$ct),'sub_type'] <- 'Fibroblasts'
subtype[grepl('B cell',subtype$ct),'sub_type'] <- 'Bcells'
subtype[grepl('T/NK cell',subtype$ct),'sub_type'] <- 'NKcells'
subtype[grepl('CD8',subtype$ct),'sub_type'] <- 'Tcells'
subtype[grepl('acrophage',subtype$ct),'sub_type'] <- 'Macrophages'
# subtype[grepl('onocyte|DC',subtype)] <- 'OtherImmunecells'
table(subtype$ct[is.na(subtype$sub_type)])
patterns <- c("pithelia","asal","uminal", "duct", "keratino","ndocrine", "ilia",'olonocyte','Leydig','trophoblast','actocyte','ranulosa') 
pattern <- paste(patterns, collapse = "|")
subtype[grepl(pattern,subtype$ct), 'sub_type'] <-  'Tissue_Specific_Cells'
table(subtype$ct, subtype$sub_type)
BM_MAST_meta_p2[index, c('ct','sub_type')] <- subtype
table(BM_MAST_meta_p2$sub_type)
CT <- as.data.frame(table(BM_MAST_meta_p2$ct, BM_MAST_meta_p2$sub_type))
CT <- CT[CT$Freq!=0,] 
intersect(unique(BM_MAST_meta_p1$project_id), unique(BM_MAST_meta_p2$project_id))
# "GSE185965"
BM_MAST_meta_p2 <- BM_MAST_meta_p2[BM_MAST_meta_p2$project_id != "GSE185965",]
CT <- as.data.frame(table(wrong_3$ct, wrong_3$sub_type))
CT <- CT[CT$Freq!=0,] 
# index1 <- grepl('eutrophil',wrong_3$sub_type)
# wrong_3$sub_type[index1] <- 'Neutrophils'
unique(disco_majortype$project_id) # 301 --> 275
BM <- disco_majortype[disco_majortype$tissue == 'bone marrow',] # 29 datasets
disco_BM_MAST <- readRDS(file=paste0('/data2/rluo4/RPMfunc/disco_pdata/BM_MAST_meta.rds'))
# leukemia <- disco[grepl('eukemia', disco$disease),]
table(unique(BM$project_id) %in% unique(disco_majortype$project_id))#29
table(unique(BM$project_id) %in% unique(BM_MAST_meta_p1$project_id))#8-->10
table(unique(BM$project_id) %in% unique(disco$project_id))#16-->19
table(unique(BM$project_id) %in% unique(disco_BM_MAST$project_id))#14-->16
#
AML <- readRDS(file=paste0('/data2/rluo4/RPMfunc/disco_pdata/AML_meta.rds')) # made through /data2/rluo4/Rcode/RMDB/DiscoDatasets.R
#
PDAC1 <- readRDS(file=paste0('/data2/rluo4/RPMfunc/disco_pdata/PDAC1_meta.rds'))
# delete the wrong datasets (mast/basophils) from disco
pancreas_pid <- c('GSE205049','GSE217837','GSE231535')
excluded_data <- c(unique(AML$project_id), pancreas_pid, unique(BM$project_id), 
                   unique(disco_BM_MAST$project_id)) #107 excluded datasets -->115
sort(excluded_data)
excluded_data <- data.frame(excluded_datasets = excluded_data, source = c(rep('AML', length(unique(AML$project_id))),
                                                                          rep('pancreas', length(unique(pancreas_pid))), rep('BM', length(unique(BM$project_id))), 
                                                                          rep('BM_MAST', length(unique(disco_BM_MAST$project_id)))) # exclude 3 BM_MAST datasets c('GSE211033','GSE196052','GSE196676')
) 
disco_majortype <- disco_majortype[!disco_majortype$project_id %in% excluded_data$excluded_datasets,]
table(index_ct_subtype$Major_type)
table(AML$sub_type)
AML$Major_type <- AML$sub_type
index1 <- grepl('Bcells|Macrophages|Neutrophils|NKcells|OtherImmunecells|Tcells',AML$Major_type)
AML$Major_type[index1] <- 'Immune_Cells'
index2 <- grepl('OtherMesenchymals',AML$Major_type)
AML$Major_type[index2] <- 'Mesenchymal_Cells'
# AML$Major_type <- index_ct_subtype$Major_type[match( AML$sub_type, index_ct_subtype$sub_type)]
table(AML$sub_type, AML$Major_type)
head(AML)
AML$barcode <- rownames(AML)
disco_majortype <- rbind(disco_majortype, AML[, colnames(disco_majortype)])

PDAC1$Major_type <- PDAC1$sub_type
index1 <- grepl('Bcells|Macrophages|Neutrophils|NKcells|OtherImmunecells|Tcells|MastCell',PDAC1$Major_type)
PDAC1$Major_type[index1] <- 'Immune_Cells'
PDAC1$barcode <- rownames(PDAC1)
disco_majortype <- rbind(disco_majortype, PDAC1[, colnames(disco_majortype)])

BM_MAST_meta_p1$Major_type <- BM_MAST_meta_p1$sub_type
index1 <- grepl('Bcells|Macrophages|Neutrophils|NKcells|MastCell|Myeloids|OtherImmunecells|Tcells',BM_MAST_meta_p1$Major_type)
BM_MAST_meta_p1$Major_type[index1] <- 'Immune_Cells'
index2 <- grepl('OtherMesenchymals|Fibroblasts',BM_MAST_meta_p1$Major_type)
BM_MAST_meta_p1$Major_type[index2] <- 'Mesenchymal_Cells'
table(BM_MAST_meta_p1$Major_type)

BM_MAST_meta_p2$Major_type <- BM_MAST_meta_p2$sub_type
index1 <- grepl('Bcells|Macrophages|Neutrophils|NKcells|MastCell|Myeloids|OtherImmunecells|Tcells',BM_MAST_meta_p2$Major_type)
BM_MAST_meta_p2$Major_type[index1] <- 'Immune_Cells'
index2 <- grepl('OtherMesenchymals|Fibroblasts',BM_MAST_meta_p2$Major_type)
BM_MAST_meta_p2$Major_type[index2] <- 'Mesenchymal_Cells'
table(BM_MAST_meta_p2$Major_type)

missed_data2$Major_type <- missed_data2$sub_type
index1 <- grepl('Bcells|Macrophages|Neutrophils|NKcells|MastCell|Myeloids|OtherImmunecells|Tcells',missed_data2$Major_type)
missed_data2$Major_type[index1] <- 'Immune_Cells'
index2 <- grepl('OtherMesenchymals|Fibroblasts',missed_data2$Major_type)
missed_data2$Major_type[index2] <- 'Mesenchymal_Cells'
table(missed_data2$Major_type)

missed_data1$Major_type <- missed_data1$sub_type
index1 <- grepl('Bcells|Macrophages|Neutrophils|NKcells|MastCell|Myeloids|OtherImmunecells|Tcells',missed_data1$Major_type)
missed_data1$Major_type[index1] <- 'Immune_Cells'
index2 <- grepl('OtherMesenchymals|Fibroblasts',missed_data1$Major_type)
missed_data1$Major_type[index2] <- 'Mesenchymal_Cells'
table(missed_data1$Major_type)

wrong_3$Major_type <- wrong_3$sub_type
index1 <- grepl('Bcells|Macrophages|Neutrophils|NKcells|MastCell|Myeloids|OtherImmunecells|Tcells',wrong_3$Major_type)
wrong_3$Major_type[index1] <- 'Immune_Cells'
index2 <- grepl('OtherMesenchymals|Fibroblasts',wrong_3$Major_type)
wrong_3$Major_type[index2] <- 'Mesenchymal_Cells'
table(wrong_3$Major_type)

disco_majortype_add_p1 <- BM_MAST_meta_p1 %>% 
  select(project_id, patient_id = sample_id, ct, tissue, platform, 
         sample_type, disease, Major_type, sub_type)
disco_majortype_add_p2 <- BM_MAST_meta_p2 %>% 
  select(project_id, patient_id = sample_id, ct, tissue, platform, 
         sample_type, disease, Major_type, sub_type)
disco_majortype_add_p3 <- missed_data2 %>%
  select(project_id, patient_id = sample_id, ct, tissue, platform,
         sample_type, disease, Major_type, sub_type)
disco_majortype_add_p4 <- wrong_3 %>%
  select(project_id, patient_id = sample_id, ct, tissue, platform,
         sample_type, disease, Major_type, sub_type)
disco_majortype_add_p5 <- missed_data1 %>%
  select(project_id, patient_id = sample_id, ct, tissue, platform,
         sample_type, disease, Major_type, sub_type)

table(disco_majortype_add_p1$sample_type)
GSE201425 <- subset(disco_majortype_add_p1, project_id=='GSE201425') 
unique(GSE201425$disease) # cholangiocarcinoma, but include other tissue from tumor patient
table(GSE201425$tissue) # blood include the PBMC and peripheral blood
# bile duct      blood      liver lymph node peritoneum 
# 25321      28696       7894      11960       4443 
disco_pheno <- as.data.frame(table(disco_majortype_add_p1$sample_type, disco_majortype_add_p1$disease))
disco_pheno <- disco_pheno[disco_pheno$Freq!=0,]
sample.index <- c(grep("ormal", disco_pheno$Var1),grep("ealthy", disco_pheno$Var1),
                  grep("control", disco_pheno$Var1))
disco_pheno[sample.index,]
table(unique(disco_majortype_add_p1$project_id) %in% unique(sc_majortype$project_id))#30-->40
table(unique(disco_majortype_add_p1$project_id) %in% unique(disco_majortype$project_id))#40
# sort(intersect(disco_BM_MAST$project_id, sc_cell_summary$project_id))

disco_majortype_add_p1$barcode <- rownames(disco_majortype_add_p1)
disco_majortype <- rbind(disco_majortype, disco_majortype_add_p1[, colnames(disco_majortype)])
disco_majortype_add_p2$barcode <- rownames(disco_majortype_add_p2)
disco_majortype <- rbind(disco_majortype, disco_majortype_add_p2[, colnames(disco_majortype)])
table(unique(disco_majortype$project_id) %in% c('GSE211033','GSE196052','GSE196676'))
table(unique(disco_majortype$project_id) %in% BM_MAST_newdatasets)
disco_majortype_add_p3$barcode <- rownames(disco_majortype_add_p3)
disco_majortype <- rbind(disco_majortype, disco_majortype_add_p3[, colnames(disco_majortype)])
disco_majortype_add_p4$barcode <- rownames(disco_majortype_add_p4)
disco_majortype <- rbind(disco_majortype, disco_majortype_add_p4[, colnames(disco_majortype)])
# saveRDS(disco_majortype,  file=paste0('/data2/rluo4/RPMfunc/Output/scRNA/disco_all_adjustmeta_0827.rds'))
# All Disco project id #
all_tissue_project <- read_excel('/data2/rluo4/EpiTrans/DataCollection/all_tissue_project.xlsx', sheet = 1)
all_tissue_project <- apply(all_tissue_project, 2, function(x){
  y <- str_split(x, '[ （]', simplify = T)[,1]
  return(y)
})
# Convert the matrix to a data frame
all_tissue_project_df <- as.data.frame(all_tissue_project)
# Gather the data into a long format with "tissues" and "project_id"
library(tidyr)
long_df <- gather(all_tissue_project_df, key = "tissues", value = "project_id")
# Remove rows where "project_id" is NA
long_df <- long_df[!is.na(long_df$project_id), ]
# View the new data frame
print(long_df)
table(unique(sc_majortype$project_id) %in% long_df$project_id)
setdiff(unique(sc_majortype$project_id),long_df$project_id)
table(long_df$project_id %in% unique(sc_majortype$project_id))
setdiff(long_df$project_id, unique(sc_majortype$project_id))
# 325 = 129 healthy + 196
################################################################################
# run on UTH36
################################################################################
# disco_supp <- list.files('/data2/rluo4/EpiTrans/DataCollection/disco_supplement_meta')
# datlist <- lapply(disco_supp, function(x){
#   readRDS(file = file.path('/data2/rluo4/EpiTrans/DataCollection/disco_supplement_meta/',x)) 
# })
# total_supp <- data.frame()
# for (GSE in 1:length(disco_supp)) {
#   names(datlist)[[GSE]] <- strsplit(disco_supp[GSE],'_')[[1]][1]
#   dataL <- datlist[[GSE]]
#   dataL$project_id <- names(datlist)[[GSE]]
#   total_supp <<- rbind(total_supp, dataL)
# }
# table(is.na(total_supp$ct))
# total_supp <- total_supp[!is.na(total_supp$ct),]
# total_supp$sub_type <- as.character(total_supp$ct)
# table(total_supp$sub_type)
# total_supp$sub_type[grepl('acrophage',total_supp$sub_type)] <- 'Macrophages'
# total_supp$sub_type[grepl('ibroblast',total_supp$sub_type)] <- 'Fibroblasts'
# total_supp$sub_type[grepl('B cell|B cells',total_supp$sub_type)] <- 'Bcells'
# total_supp$sub_type[grepl('CD4|CD8| T|MAIT',total_supp$sub_type)] <- 'Tcells'
# total_supp$sub_type[grepl('M2',total_supp$sub_type)] <- 'Macrophages'
# total_supp$sub_type[grepl('NK',total_supp$sub_type)] <- 'NKcells'
# total_supp$sub_type <- gsub('Plasma cell','Bcells',total_supp$sub_type)
# total_supp$sub_type <- gsub('Hematopoietic stem cell','Stem_Cells',total_supp$sub_type)
# total_supp$sub_type <- gsub('Mast cells','MastCell',total_supp$sub_type)
# total_supp$sub_type[grepl('onocyte|DC|eutrophil|myelo|lymphoid|ILC',total_supp$sub_type)] <- 'OtherImmunecells'
# index1 <- grepl('CD11c+|monocyte|endritic|DC|LymP|GMP',total_supp$ct)
# total_supp$sub_type[index1] <- 'OtherImmunecells'
# index2 <- grepl('Red blood|erythrocyte|erythroblast|egakaryocyte|Ery|MEP',total_supp$ct)
# total_supp$sub_type[index2] <- 'Blood_Cells'
# index3 <- grepl('EC',total_supp$ct)
# total_supp$sub_type[index3] <- 'Endothelial_Cells'
# CT <- as.data.frame(table(total_supp$ct, total_supp$sub_type))
# table(total_supp$sub_type)
# unique(total_supp$sub_type)
# unique(total_supp$project_id)
# saveRDS(total_supp, file = '/data2/rluo4/EpiTrans/DataCollection/GSE_5datasets.rds')
################################################################################
# 6. Check the last ascp links
################################################################################
## last 100+ datasets on UTH57
UTH57_Bam <- list.files('/data2/rluo4/EpiTrans/Disco/Raw/Bam/') #70
UTH57_FastQ <- list.files('/data2/rluo4/EpiTrans/Disco/Raw/FastQ/') #70
UTH57_Sra <- list.files('/data2/rluo4/EpiTrans/Disco/Raw/Sra/') #70
UTH57_Raw <- unique(c(UTH57_Bam, UTH57_FastQ, UTH57_Sra)) #26
# Raw_done <- intersect(UTH57_Raw, BM_MAST39$project_id)
Raw_done <- intersect(UTH57_Raw, BM_MAST_meta_p1$project_id)
# for (i in UTH57_Raw ) {
#   GSM <-  list.files(file.path("/data2/rluo4/EpiTrans/Disco/scMapping/",i ))
#   print(GSM)
#   for (j in GSM) {
#     print(paste('to be removed dataset: ', i))
#     print(paste('to be removed GSM: ', j))
#     scMapping_files <-  list.files(file.path("/data2/rluo4/EpiTrans/Disco/scMapping/",i,j ))
#     print(scMapping_files)
#   }
# }
################################################################################
# 5.1 Renew the last ascp links (last means the samples downloaded in the UTH57 not UTH36)
################################################################################
disco <- readRDS(paste0('/data2/rluo4/RPMfunc/disco_pdata/disco_adjustmeta_0827.rds'))#0809.rds'))#0725.rds'))
table(disco$sample_type)
sample.index <- c(grep("ormal", disco$sample_type),grep("ealthy", disco$sample_type),
                  grep("control", disco$sample_type))#sample.index <- c(grep("ormal", disco$sample_type),grep("ealthy", disco$sample_type))
disco_diseased <-  disco[-sample.index,]
table(disco$sample_type)
sample.index <- c(grep("ormal", disco$sample_type),grep("ealthy", disco$sample_type),
                  grep("control", disco$sample_type))#sample.index <- c(grep("ormal", disco$sample_type),grep("ealthy", disco$sample_type))
disco_diseased <-  disco[-sample.index,]
table(disco_diseased$sample_type)
disco_meta <- data.frame(table(disco_diseased$sample_type, disco_diseased$disease))
table(is.na(disco_diseased$patient_id))
length(unique(disco_diseased$project_id)) #204 diseased datasets -->197 -->196-->214-->213(1 covid dataset:GSE168453)
length(unique(disco_diseased$patient_id)) #1332 diseased samples -->1359 -->1360-->1525-->1838-->1822-->1876
disco_sample_summary <- disco_diseased[! duplicated(disco_diseased$patient_id),]

data_dir <- '/data2/rluo4/EpiTrans/DataCollection/disco_tsv'
setwd(data_dir)
library(stringr)
disco_GEO <- list.files('./', pattern='.txt') #/data2/rluo4/EpiTrans/DataCollection/disco_tsv
disco_GEO.acc <- str_split(gsub('.txt','',disco_GEO), '_', simplify = T)[,1]
unique(disco_GEO.acc)#117 --> 294
table(unique(disco_GEO.acc) %in% unique(disco_diseased$project_id))
setdiff(unique(disco_GEO.acc), unique(disco_diseased$project_id))

tissue_summary <- read.table('/data2/rluo4/RPMfunc/PCTanno_pdata/tissue_summary.txt')#read.csv('/data2/rluo4/All/Output/tissue_summary.txt',sep = "\t",header = T)
PCTanno_GEO <- unique(tissue_summary$Dataset)
table(disco_GEO.acc %in% PCTanno_GEO)
PCTanno_disco <- intersect(disco_GEO.acc, PCTanno_GEO) # 9 datasets --> 17 intersection
PCTanno_disco
# View(tissue_summary[tissue_summary$Dataset %in% PCTanno_disco,])
disco_GEO <- disco_GEO[!disco_GEO.acc %in% PCTanno_GEO]
disco_GEO.acc <- str_split(gsub('.txt','',disco_GEO), '_', simplify = T)[,1]#281 datasets with raw data
table(disco_GEO.acc %in% unique(disco_diseased$project_id))
disco_GEO <- disco_GEO[disco_GEO.acc %in% unique(disco_diseased$project_id)]#165-->170-->168
# cell_summary <- read.csv('cell_summary.txt',sep = "\t",header = T)
library(data.table)
filePath <- sapply(disco_GEO, function(x){ 
  paste(data_dir,x,sep='/')})  
disco_GEO <- lapply(filePath, function(x){
  fread(x)})
disco_allDat <- NULL#data.frame()
unique(disco_diseased$patient_id)
disco_sum <- lapply(disco_GEO,function(y){
  disco_allDat <<- rbind(disco_allDat, y)
  return(disco_allDat)  
})
# GSE162498 --> GSE162500 (UTH57), GSE121636 --> GSE121638 (UTH36)
table(disco_allDat$sample_alias %in% disco_diseased$patient_id)
table(disco_diseased$project_id %in% c('GSE162498','GSE121636'))
table(disco_diseased$project_id %in% c('GSE162500','GSE121638'))
disco_diseased$project_id <- gsub('GSE121638', 'GSE121636',
                                  gsub('GSE162500', 'GSE162498', disco_diseased$project_id))
table(disco_allDat$sample_alias %in% disco_diseased$patient_id)
disco_allDat$sample_alias[disco_allDat$study_alias=='GSE181294'] <- disco_allDat$sample_title[disco_allDat$study_alias=='GSE181294']
disco_allDat <- disco_allDat[disco_allDat$sample_alias %in% disco_diseased$patient_id,]
unique(disco_allDat$submitted_aspera)
unique(disco_allDat$sample_alias)# 413 sample for raw data --> 466 --> 475 --> 466 -->826-->810 -->882
colnames(disco_allDat)
disco_allDat <- disco_allDat[, c('study_alias','sample_alias','fastq_aspera','submitted_aspera','sra_aspera','run_alias')]
disco_allDat$submitted_aspera[is.na(disco_allDat$submitted_aspera)] <- disco_allDat$fastq_aspera[is.na(disco_allDat$submitted_aspera)]
disco_allDat$submitted_aspera
disco_allDat$submitted_aspera[disco_allDat$submitted_aspera==''] <- disco_allDat$fastq_aspera[disco_allDat$submitted_aspera=='']
disco_allDat$submitted_aspera
disco_allDat$datatype[grepl('.bam', disco_allDat$submitted_aspera)] <- 'Bam'
disco_allDat$datatype[grepl('_1.fastq', disco_allDat$submitted_aspera)] <- 'FastQ'
disco_allDat$datatype[is.na(disco_allDat$datatype)] <- 'Sra'
index = disco_allDat$datatype=='Sra'
disco_allDat$submitted_aspera[index] <- disco_allDat$sra_aspera[index]
disco_allDat_GSE <- paste(disco_allDat$study_alias, disco_allDat$study_alias, sep = ';')
disco_allDat_GSE[index] <- disco_allDat$study_alias[index]
disco_allDat_GSE
disco_allDat_ascp <- unlist(str_split(disco_allDat$submitted_aspera,";"))#submitted_aspera#fastq_ftp#
disco_allDat_GSM <- paste(disco_allDat$sample_alias, disco_allDat$sample_alias, sep = ';')
disco_allDat_GSM[index] <- disco_allDat$sample_alias[index]
disco_allDat_GSM
disco_allDat_dt <- paste(disco_allDat$datatype, disco_allDat$datatype, sep = ';')
disco_allDat_dt[index] <- disco_allDat$datatype[index]
disco_allDat_dt
asp_sc <- data.frame(study_alias = unlist(str_split(disco_allDat_GSE,";")), 
                     sample_alias = unlist(str_split(disco_allDat_GSM,";")), 
                     sra_aspera = disco_allDat_ascp, datatype = unlist(str_split(disco_allDat_dt,";"))
)
asp_sc_last <- read.table(file="/data2/rluo4/EpiTrans/DataCollection/asp_disco_scRNA_last_renew76.link.txt",sep = "\t")

asp_sc_first <- read.table(file="/data2/rluo4/EpiTrans/DataCollection/asp_disco_scRNA.link.txt",sep = "\t")
table(asp_sc_first$V3 %in% asp_sc$sra_aspera) 
table(asp_sc$sra_aspera %in% asp_sc_first$V3) 
asp_sc_last <- asp_sc[! asp_sc$sra_aspera %in% asp_sc_first$V3,]
intersect(asp_sc_first$V1, asp_sc_last$study_alias)
# [1] "GSE130560" "GSE150825" "GSE151192" "GSE155468"
# split 4 renewed datasets into 2 parts in 2 servers:
# UTH57: "GSE130560" "GSE150825"
# UTH36: "GSE151192" "GSE155468"

# write.table(asp_sc_last,file="/data2/rluo4/EpiTrans/DataCollection/asp_disco_scRNA_last.link.txt",sep = "\t",quote = F, col.names = F, row.names = F)
# write.table(asp_sc_last,file="/data2/rluo4/EpiTrans/DataCollection/asp_disco_scRNA_AML.link.txt",sep = "\t",quote = F, col.names = F, row.names = F)
# paste0(unique(AML$patient_id[AML$disease=='Control']), collapse = "|")
# cat $infolder/asp_disco_scRNA_AML.link.txt|grep -E 'GSE185381|GSE185991|GSE211036'|grep -Ev 'GSM6412447|GSM6412451|GSM6412453|GSM6412454|GSM6447735'|grep -E 'GSM5613747|GSM5613748|GSM5613749|GSM5613750|GSM5613751|GSM5613752|GSM5613756|GSM5613757|GSM5613758|GSM5613759|GSM5613769|GSM5613770|GSM5613771|GSM5613774|GSM5613775|GSM5613776|GSM5613781|GSM5613787'
asp_sc_last_AML <- read.table(file="/data2/rluo4/EpiTrans/DataCollection/asp_disco_scRNA_AML.link.txt",sep = "\t") # made through /data2/rluo4/Rcode/RMDB/DiscoDatasets.R
# f5 <- c('GSM6412447', 'GSM6412451', 'GSM6412453', 'GSM6412454', 'GSM6447735')
# asp_sc_last <- asp_sc_last[ (asp_sc_last$V1 %in% unique(AML$project_id)) & (!asp_sc_last$V2 %in% f5),] # 170
# which(asp_sc_last$V2 %in% unique(AML$patient_id[AML$disease=='Control']))
# write.table(asp_sc_last,file="/data2/rluo4/EpiTrans/DataCollection/asp_disco_scRNA_disease.link.txt",sep = "\t",quote = F, col.names = F, row.names = F)
# write.table(asp_sc_last,file="/data2/rluo4/EpiTrans/DataCollection/asp_disco_scRNA_last_renew39.link.txt",sep = "\t",quote = F, col.names = F, row.names = F)
asp_sc_last_disease <- read.table(file="/data2/rluo4/EpiTrans/DataCollection/asp_disco_scRNA_disease.link.txt",sep = "\t") # made through /data2/rluo4/Rcode/RMDB/DiscoDatasets.R
asp_sc_last_renew39 <- read.table(file="/data2/rluo4/EpiTrans/DataCollection/asp_disco_scRNA_last_renew39.link.txt",sep = "\t") # made through /data2/rluo4/Rcode/RMDB/DiscoDatasets.R
# write.table(asp_sc_last,file="/data2/rluo4/EpiTrans/DataCollection/asp_disco_scRNA_last_renew72.link.txt",sep = "\t",quote = F, col.names = F, row.names = F)
asp_sc_last_renew72 <- read.table(file="/data2/rluo4/EpiTrans/DataCollection/asp_disco_scRNA_last_renew72.link.txt",sep = "\t") # made through /data2/rluo4/Rcode/RMDB/DiscoDatasets.R
# write.table(asp_sc_last,file="/data2/rluo4/EpiTrans/DataCollection/asp_disco_scRNA_last_renew74.link.txt",sep = "\t",quote = F, col.names = F, row.names = F)
asp_sc_last_renew74 <- read.table(file="/data2/rluo4/EpiTrans/DataCollection/asp_disco_scRNA_last_renew74.link.txt",sep = "\t") # made through /data2/rluo4/Rcode/RMDB/DiscoDatasets.R
# write.table(asp_sc_last,file="/data2/rluo4/EpiTrans/DataCollection/asp_disco_scRNA_last_renew76.link.txt",sep = "\t",quote = F, col.names = F, row.names = F)

################################################################################
# 5.2 Renew the Disco sample metadata for SComatic analysis
################################################################################
# Renew on 11/22/2024:
AML1 <- readRDS(file=paste0('/data2/rluo4/RPMfunc/disco_pdata/GSE185381_meta.rds'))
disco <- readRDS(paste0('/data2/rluo4/RPMfunc/disco_pdata/disco_adjustmeta_1122.rds'))
RMP_update <- readxl::read_excel('/data2/rluo4/EpiTrans/DataCollection/RMP_update.xlsx')
# cell_summary <- read.csv('cell_summary.txt',sep = "\t",header = T)
sc_majortype <- readRDS(file = '/data2/rluo4/RPMfunc/Output/scRNA/sc_majortype_1022.rds')
# tissue_summary <- read.table('/data2/rluo4/RPMfunc/PCTanno_pdata/tissue_summary.txt')
# sc_allmeta <- readRDS(file = '/data2/rluo4/RPMfunc/Output/scRNA/sc_allmeta_1022.rds')
RMzyme_alltissues <- sort(unique(sc_allmeta$tissue));
RMzyme_alldatasets <- sort(unique(sc_allmeta$project_id));
load(file = '/data2/rluo4/RPMfunc/Output/scRNA/sc_allsample_1022.RData', verbose = T)
extract_last_element <- function(x) {
  split_string <- strsplit(x, "_")
  last_element <- sapply(split_string, function(y) tail(y, n = 1))
  return(last_element)
}
GSE185381_majortype <- AML1
GSE185381_majortype$barcode <- rownames(GSE185381_majortype)
GSE185381_majortype$Major_type <- NA#sc_majortype$Major_type[match(AML1$sub_type, sc_majortype$sub_type)]
GSE185381_majortype <- GSE185381_majortype %>% dplyr::select(colnames(sc_majortype))
sc_majortype <- sc_majortype[!sc_majortype$project_id %in% 'GSE185381', ]
sc_majortype <- rbind(sc_majortype, GSE185381_majortype)
unique(sc_majortype$project_id[grepl('[...]', sc_majortype$patient_id)])
# [1] "GSE112271" "GSE171555" "GSE121080"
sc_majortype$patient_id[sc_majortype$project_id=='GSE121080'] <- str_split(sc_majortype$patient_id[sc_majortype$project_id=='GSE121080'],
                                                                           '[...]', simplify = T)[, 1]
rownames(sc_majortype) <- 1:nrow(sc_majortype)
sc_summary <- sc_majortype[! duplicated(paste(sc_majortype$patient_id, sc_majortype$tissue, 
                                              sc_majortype$project_id, sc_majortype$disease)),]
# saveRDS(sc_majortype, file = '/data2/rluo4/RPMfunc/Output/scRNA/sc_majortype_1122.rds') # combined with diff meta from merge_nmf_meta and allmeta

# Renew on 12/02/2024:
# sc_meta generated by testing on allmeta_merge_renew.R using GSE185381 as an example
library(maftools)
library(data.table)
library(stringr)
library(dplyr)
options(bitmapType = 'cairo')
load(file = '/data2/rluo4/RPMfunc/Output/scRNA/sc_allsample_1122.RData', verbose = T)
##################################################
#  Re-organization of merge nmf metadata
##################################################
large_split <- fread('/data2/rluo4/EpiTrans/DataCollection/large_split.txt')
# Convert data.table columns to a vector of dataset IDs
split_ddt_ids <- unlist(large_split, use.names = FALSE)
# Remove empty strings
split_ddt_ids <- split_ddt_ids[split_ddt_ids != ""]

allmeta_path <- "/data2/rluo4/RPMfunc/Output/scRNA/allmeta"
allmeta_files <- list.files(path = allmeta_path, pattern = "allmeta.rds$", 
                            full.names = TRUE, recursive = T)
print(allmeta_files)
# View(RMzyme_all_meta[RMzyme_all_meta$project_id %in% na_barcode_projectid,])
disease_tissue_ids <- str_split(allmeta_files,'/', simplify = T)[,8]
table(unique(disease_tissue_ids) %in% sc_summary$tissue)
setdiff(RMzyme_alltissues, unique(disease_tissue_ids))

disease_project_ids <- str_split(allmeta_files,'/', simplify = T)[,9] # (255 + 1 GSE110686_allmeta.rds)
# missed_bc_files <- allmeta_files[disease_project_ids %in% paste0(na_barcode_projectid,'_allmeta.rds')]
sc_meta<- NULL 
# for (i in 1:length(allmeta_files)){
i=44
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
  # Regular expressions for different formats
  barcodes <- meta$rownames[na_proj_idx]
  # Updated regex approach for all cases
  # Function to extract patient_id
  patient_ids <- extract_patient_id(barcodes)
  meta$patient_id[na_proj_idx] <- patient_ids
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
# check the patient_id
patient.idx <- paste(sc_meta$patient_id, sc_meta$project_id)
match.idx <- match(patient.idx, paste(sc_summary$patient_id, sc_summary$project_id))
sc_meta$Major_type <- NA #sc_summary$Major_type[match.idx]
sc_meta$sample_type <- sc_summary$sample_type[match.idx]
sc_meta$platform <- sc_summary$platform[match.idx]
##################################################
index1 <- grepl('Bcells|Macrophages|Neutrophils|NKcells|OtherImmunecells|Tcells|MastCell|Myeloids',sc_meta$sub_type)
sc_meta$Major_type[index1] <- 'Immune_Cells'
index2 <- grepl('OtherMesenchymals|Fibroblasts',sc_meta$sub_type)
sc_meta$Major_type[index2] <- 'Mesenchymal_Cells'
index3 <- grepl('Endothelial_Cells',sc_meta$sub_type)
sc_meta$Major_type[index3] <- 'Endothelial_Cells'
index4 <- grepl('Tissue_Specific_Cells',sc_meta$sub_type)
sc_meta$Major_type[index4] <- 'Tissue_Specific_Cells'
index5 <- grepl('Stem_Cells',sc_meta$sub_type)
sc_meta$Major_type[index5] <- 'Stem_Cells'
index6 <- grepl('Blood_Cells',sc_meta$sub_type)
sc_meta$Major_type[index6] <- 'Blood_Cells'
table(sc_meta$Major_type)
table(is.na(sc_meta$Major_type))
sc_meta$patient_id <- str_split(sc_meta$patient_id, '[...]', simplify = T)[, 1]

setdiff( unique(sc_meta$disease), unique(sc_summary$disease))
setdiff( unique(sc_meta$tissue), unique(sc_summary$tissue))
sc_meta$disease <- gsub('HCC', 'hepatocellular carcinoma', 
                        gsub('Healthy', 'healthy',
                             gsub('oral cancer', 'oral squamous cell carcinoma',
                                  gsub('PRAD', 'pancreatic ductal adenocarcinoma',
                                       gsub('cSCC', 'oral squamous cell carcinoma',
                                            gsub('AML', 'acute myeloid leukemia',
                                                 sc_meta$disease))))))
# sc_meta generated by testing on allmeta_merge_renew.R using GSE185381 as an example
ct_dts_barcode <- str_split(sc_meta$barcode, '[-1]', simplify = T)[, 1]
sc_meta$barcode <- paste0(ct_dts_barcode, '--', sc_meta$patient_id)
saveRDS(sc_meta, file = '/data2/rluo4/RPMfunc/disco_pdata/GSE185381_allmeta.rds')
# sc_meta <- readRDS('/data2/rluo4/RPMfunc/disco_pdata/GSE185381_allmeta.rds')
sc_allmeta <- sc_allmeta[!sc_allmeta$ddt_id %in% 'GSE185381', ]
sc_allmeta <- rbind(sc_allmeta, sc_meta)
# saveRDS(sc_allmeta, file = '/data2/rluo4/RPMfunc/Output/scRNA/sc_allmeta_1122.rds') # combined with diff meta from merge_nmf_meta and allmeta
# sc_allsample <- sc_allmeta %>% dplyr::select(patient_id, project_id, tissue, disease) %>%  distinct()
# save(disco, tissue_summary, sc_summary, sc_allsample, file = '/data2/rluo4/RPMfunc/Output/scRNA/sc_allsample_1122.RData')
################################################################################
# Renew on 12/07/2024:
################################################################################
library(maftools)
library(data.table)
library(stringr)
library(dplyr)
options(bitmapType = 'cairo')
sc_allmeta <- readRDS(file = '/data/rluo4/RPMfunc/Output/scRNA/sc_allmeta_1122.rds')
load(file = '/data/rluo4/RPMfunc/Output/scRNA/sc_allsample_1122.RData', verbose = T)
##################################################
# change 4                                       #
##################################################
patient_id <- sc_allmeta$patient_id
barcode_idx <- sc_allmeta$barcode
idx9 <- grepl('GSE181254|PRJCA001063|GSE148673|GSE163203', sc_allmeta$project_id) 
sc_allmeta$barcode[idx9][1:5]
barcode_idx[idx9] <- extract_last_element(sc_allmeta$rownames[idx9])
# View(sc_allmeta[grepl('[...]', sc_allmeta$rownames)  & idx9,])
# table(sc_allmeta[grepl('[...]', sc_allmeta$rownames) & idx9, 'project_id'])
barcode_ct <- str_split(barcode_idx[idx9], '[...]', simplify = T)[, 1]
table(barcode_ct == barcode_idx[idx9])
table(sc_allmeta$project_id[idx9])
barcode_idx[idx9] <- barcode_ct
barcode_idx[idx9] <- paste0(barcode_idx[idx9], '--', patient_id[idx9])

idx_dt <-  idx9
dt <- unique(sc_allmeta[idx_dt, 'project_id'])
print(dt)
sc_allmeta$barcode[idx_dt] <- barcode_idx[idx_dt]
saveRDS(sc_allmeta, file = '/data/rluo4/RPMfunc/Output/scRNA/sc_allmeta_1122.rds') # combined with diff meta from merge_nmf_meta and allmeta
sc_allsample <- sc_allmeta %>% dplyr::select(patient_id, project_id, tissue, disease) %>%  distinct()
save(disco, tissue_summary, sc_summary, sc_allsample, file = '/data/rluo4/RPMfunc/Output/scRNA/sc_allsample_1122.RData')

################################################################################
sc_allmeta$omics <- 'scRNA'
sample.index <- c(grep("ormal", sc_summary$sample_type),grep("ealthy", sc_summary$sample_type),
                  grep("ontrol", sc_summary$sample_type))
sc_Diseased_Samples <- sc_summary[ -sample.index, c(2,3,5:8)] # 219 --> 222 --> 224 disease datasets
disco_diseased <- sc_Diseased_Samples # checked GSM5613751 that have both AML and control cells
table(disco_diseased$sample_type)
disco_meta <- data.frame(table(disco_diseased$sample_type, disco_diseased$disease))
table(is.na(disco_diseased$patient_id))
length(unique(disco_diseased$project_id)) #204 diseased datasets -->197 -->196-->214-->213(1 covid dataset:GSE168453)
length(unique(disco_diseased$patient_id)) #1332 diseased samples -->1359 -->1360-->1525-->1838-->1822-->1876
disco_sample_summary <- disco_diseased[! duplicated(disco_diseased$patient_id),]

data_dir <- '/data2/rluo4/EpiTrans/DataCollection/disco_tsv'
setwd(data_dir)
library(stringr)
disco_GEO <- list.files('./', pattern='.txt') #/data2/rluo4/EpiTrans/DataCollection/disco_tsv
disco_GEO.acc <- str_split(gsub('.txt','',disco_GEO), '_', simplify = T)[,1]
unique(disco_GEO.acc)#117 --> 294
table(unique(disco_GEO.acc) %in% unique(disco_diseased$project_id))
setdiff(unique(disco_GEO.acc), unique(disco_diseased$project_id))

tissue_summary <- read.table('/data2/rluo4/RPMfunc/PCTanno_pdata/tissue_summary.txt')#read.csv('/data2/rluo4/All/Output/tissue_summary.txt',sep = "\t",header = T)
PCTanno_GEO <- unique(tissue_summary$Dataset)
table(disco_GEO.acc %in% PCTanno_GEO)
PCTanno_disco <- intersect(disco_GEO.acc, PCTanno_GEO) # 9 datasets --> 17 intersection
PCTanno_disco
# View(tissue_summary[tissue_summary$Dataset %in% PCTanno_disco,])
disco_GEO <- disco_GEO[!disco_GEO.acc %in% PCTanno_GEO]
disco_GEO.acc <- str_split(gsub('.txt','',disco_GEO), '_', simplify = T)[,1]#281 datasets with raw data
table(disco_GEO.acc %in% unique(disco_diseased$project_id))
disco_GEO <- disco_GEO[disco_GEO.acc %in% unique(disco_diseased$project_id)]#165-->170-->168
# cell_summary <- read.csv('cell_summary.txt',sep = "\t",header = T)
library(data.table)
filePath <- sapply(disco_GEO, function(x){ 
  paste(data_dir,x,sep='/')})  
disco_GEO <- lapply(filePath, function(x){
  fread(x)})
disco_allDat <- NULL#data.frame()
unique(disco_diseased$patient_id)
disco_sum <- lapply(disco_GEO,function(y){
  disco_allDat <<- rbind(disco_allDat, y)
  return(disco_allDat)  
})
# GSE162498 --> GSE162500 (UTH57), GSE121636 --> GSE121638 (UTH36)
table(disco_allDat$sample_alias %in% disco_diseased$patient_id)
table(disco_diseased$project_id %in% c('GSE162498','GSE121636'))
table(disco_diseased$project_id %in% c('GSE162500','GSE121638'))
disco_diseased$project_id <- gsub('GSE121638', 'GSE121636',
                                  gsub('GSE162500', 'GSE162498', disco_diseased$project_id))
table(disco_allDat$sample_alias %in% disco_diseased$patient_id)
disco_allDat$sample_alias[disco_allDat$study_alias=='GSE181294'] <- disco_allDat$sample_title[disco_allDat$study_alias=='GSE181294']
# disco_allDat <- disco_allDat[disco_allDat$sample_alias %in% disco_diseased$patient_id,]
# disco_allDat <- disco_allDat[disco_allDat$study_alias=='GSE185381',]
updated_disease_samples <- c(disco_diseased$patient_id, AML1$patient_id[AML1$disease=='AML'])#GSE185381$patient_id)
disco_allDat <- disco_allDat[disco_allDat$sample_alias %in% unique(updated_disease_samples),]

unique(disco_allDat$submitted_aspera)
unique(disco_allDat$sample_alias)# 413 sample for raw data --> 466 --> 475 --> 466 -->826-->810 -->882
colnames(disco_allDat)
disco_allDat <- disco_allDat[, c('study_alias','sample_alias','fastq_aspera','submitted_aspera','sra_aspera','run_alias')]
disco_allDat$submitted_aspera[is.na(disco_allDat$submitted_aspera)] <- disco_allDat$fastq_aspera[is.na(disco_allDat$submitted_aspera)]
disco_allDat$submitted_aspera
disco_allDat$submitted_aspera[disco_allDat$submitted_aspera==''] <- disco_allDat$fastq_aspera[disco_allDat$submitted_aspera=='']
disco_allDat$submitted_aspera
disco_allDat$datatype[grepl('.bam', disco_allDat$submitted_aspera)] <- 'Bam'
disco_allDat$datatype[grepl('_1.fastq', disco_allDat$submitted_aspera)] <- 'FastQ'
disco_allDat$datatype[is.na(disco_allDat$datatype)] <- 'Sra'
index = disco_allDat$datatype=='Sra'
disco_allDat$submitted_aspera[index] <- disco_allDat$sra_aspera[index]
disco_allDat_GSE <- paste(disco_allDat$study_alias, disco_allDat$study_alias, sep = ';')
disco_allDat_GSE[index] <- disco_allDat$study_alias[index]
disco_allDat_GSE
disco_allDat_ascp <- unlist(str_split(disco_allDat$submitted_aspera,";"))#submitted_aspera#fastq_ftp#
disco_allDat_GSM <- paste(disco_allDat$sample_alias, disco_allDat$sample_alias, sep = ';')
disco_allDat_GSM[index] <- disco_allDat$sample_alias[index]
disco_allDat_GSM
disco_allDat_dt <- paste(disco_allDat$datatype, disco_allDat$datatype, sep = ';')
disco_allDat_dt[index] <- disco_allDat$datatype[index]
disco_allDat_dt
asp_sc <- data.frame(study_alias = unlist(str_split(disco_allDat_GSE,";")), 
                     sample_alias = unlist(str_split(disco_allDat_GSM,";")), 
                     sra_aspera = disco_allDat_ascp, datatype = unlist(str_split(disco_allDat_dt,";"))
)
index <- asp_sc$sample_alias %in% AML1$patient_id[AML1$disease=='AML']
table(index)
table(asp_sc$study_alias=='GSE185381')
asp_sc_GSE185381 <- asp_sc[index,]
# write.table(asp_sc_GSE185381,file="/data2/rluo4/EpiTrans/DataCollection/asp_disco_scRNA_GSE185381.link.txt",sep = "\t",quote = F, col.names = F, row.names = F)

# f5 <- c('GSM6412447', 'GSM6412451', 'GSM6412453', 'GSM6412454', 'GSM6447735')
# asp_sc_last <- asp_sc_last[ (asp_sc_last$V1 %in% unique(AML$project_id)) & (!asp_sc_last$V2 %in% f5),] # 170
BM_MAST_newdatasets <- setdiff(unique(asp_sc_last$study_alias), unique(asp_sc_last_disease$V1))
BM_MAST_newdatasets <- setdiff(unique(asp_sc_last$study_alias), unique(asp_sc_last_renew39$V1))
BM_MAST_newdatasets <- setdiff(unique(asp_sc_last$study_alias), unique(asp_sc_last_renew74$V1))
# [1] "GSE139324" "GSE162025"
raw_AML1 <- list.files('/data2/rluo4/EpiTrans/Disco/Raw/cellranger_count/GSE185381/')
raw_AML2 <- list.files('/data2/rluo4/EpiTrans/Disco/Raw/cellranger_count/GSE185991/')
raw_AML3 <- list.files('/data2/rluo4/EpiTrans/Disco/Raw/cellranger_count/GSE211036/')
raw_done_AML <- c(raw_AML1, raw_AML2, raw_AML3)
paste0(raw_done_AML,collapse = '|')
uncommon <- readxl::read_excel('/data2/rluo4/RPMfunc/xuhaixia/disco_tsv/uncommon.xlsx')
snRNAGSE <- readxl::read_excel('/data2/rluo4/RPMfunc/xuhaixia/disco_tsv/snRNAGSE.xlsx',sheet = 2)
table(snRNAGSE$GSE %in% asp_sc$study_alias)
table(snRNAGSE$GSE %in% disco_GEO.acc)
table(snRNAGSE$GSE %in% disco_diseased$project_id)
setdiff(snRNAGSE$GSE, asp_sc$study_alias)
intersect(snRNAGSE$GSE, asp_sc$study_alias)#"GSE141552" "GSE144136" "GSE157783" "GSE174367"
intersect(snRNAGSE$GSE, asp_sc_last$study_alias)
unique(asp_sc$study_alias)#107-->103-->124
unique(disco_allDat_GSE)#112-->108-->129
# https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR22957234&display=data-access
# datasets with _1,_2,_3,_4.fastq.gz-->I1,I2,R1,R2_001.fastq.gz
# GSE134809: GSM3972009;
# GSE222078: GSM6914629;GSM6914627;GSM6914631;
path <- "/data2/rluo4/EpiTrans/Disco/Raw/"
setwd('/data2/rluo4/EpiTrans/Disco/Raw/Bam/')
bam_files <- list.dirs('/data2/rluo4/EpiTrans/Disco/Raw/Bam', full.names = TRUE)
bam_files <- bam_files[grep('GSE', bam_files)]
types <- c('Bam','FastQ','Sra')
file_info_allDat <- NULL
for (type in types) {
  raw_files <- list.dirs(paste0('/data2/rluo4/EpiTrans/Disco/Raw/',type), full.names = TRUE)
  raw_files <- raw_files[grep('GSE', raw_files)]
  
  filePath <- sapply(raw_files, function(x){
    y <- list.files(x)#, pattern = '.bam')
    return(paste(x,y, sep='/'))
  }
  )
  file_info <- lapply(filePath, function(x){
    y <- file.info(x)
    return(y)  
  })
  disco_sum <- lapply(file_info, function(y){
    # print(y)
    file_info_allDat <<- rbind(file_info_allDat, y)
    return(file_info_allDat)
  })
}
extract_last_element <- function(x) {
  split_string <- strsplit(x, "/")
  last_element <- sapply(split_string, function(y) tail(y, n = 1))
  return(last_element)
}
file_info_allDat$filename = extract_last_element(rownames(file_info_allDat))
table(file_info_allDat$filename %in% extract_last_element(asp_sc$sra_aspera))
asp_sc$filename <- extract_last_element(asp_sc$sra_aspera)
# asp_sc_file <- asp_sc_GSE185381#[-grep('_1.fastq.gz', asp_sc$filename), ]
# asp_sc_file <- asp_sc_file[ ! grepl('bai', asp_sc_file$sra_aspera),]
# index <- asp_sc_file$datatype=='Sra'
# asp_sc_file$filename[index] <- paste0(asp_sc_file$filename[index], '_2.fastq.gz')
# asp_sc_file$filesize <- file_info_allDat$size[match(asp_sc_file$filename, file_info_allDat$filename)]
# max(asp_sc_file$filesize[!is.na(asp_sc_file$filesize)]) #25301137945
# unique(asp_sc_file$study_alias)
# write.table(asp_sc_file, file="/data2/rluo4/EpiTrans/DataCollection/asp_sc_file_GSE185381.txt",sep = "\t",quote = F, col.names = F, row.names = F)

asp_sc_file <- asp_sc[-grep('_1.fastq.gz', asp_sc$filename), ]
# asp_sc_last$filename <- extract_last_element(asp_sc_last$sra_aspera)
# asp_sc_file <- asp_sc_last[-grep('_1.fastq.gz', asp_sc_last$filename), ]
asp_sc_file <- asp_sc_file[ ! grepl('bai', asp_sc_file$sra_aspera),]
# asp_sc_file$filename <- gsub('.bai','',asp_sc_file$filename)
index <- asp_sc_file$datatype=='Sra'
asp_sc_file$filename[index] <- paste0(asp_sc_file$filename[index], '_2.fastq.gz')
asp_sc_file$filesize <- file_info_allDat$size[match(asp_sc_file$filename, file_info_allDat$filename)]
max(asp_sc_file$filesize[!is.na(asp_sc_file$filesize)]) #25301137945
# View(asp_sc_file[is.na(asp_sc_file$filesize),])
na.asp_sc_file <- asp_sc_file[is.na(asp_sc_file$filesize),]
unique(asp_sc_file$study_alias)
# write.table(asp_sc_file, file="/data2/rluo4/EpiTrans/DataCollection/asp_sc_file_all.txt",sep = "\t",quote = F, col.names = F, row.names = F)
# saveRDS(asp_sc_file, file="/data2/rluo4/EpiTrans/DataCollection/asp_sc_file.rds")

intersect(asp_sc_file$study_alias, sc_allmeta$project_id) #77-->108
intersect(asp_sc_file$sample_alias, sc_allmeta$patient_id) #646-->823
asp_sc_file_GSE185381 <- asp_sc_file[asp_sc_file$study_alias=='GSE185381',]
# asp_sc_file_GSE185381 <- asp_sc_file[asp_sc_file$sample_alias %in% AML1$patient_id[AML1$disease=='AML'],]
# write.table(asp_sc_file_GSE185381, file="/data2/rluo4/EpiTrans/DataCollection/asp_sc_file_GSE185381.txt",sep = "\t",quote = F, col.names = F, row.names = F)
updated_GSE185381_samples <- setdiff(asp_sc_file_GSE185381$sample_alias, sc_summary$patient_id[sc_summary$project_id=='GSE185381'])
paste(updated_GSE185381_samples, collapse = '|')
# save(disco, tissue_summary, sc_summary, sc_allsample, file = '/data2/rluo4/RPMfunc/Output/scRNA/sc_allsample_1122.RData')


################################################################################
# 5.3 Organize the files for disco's SComatic: filtering the files that were done
################################################################################
GSE.first <- list.files('/data2/rluo4/EpiTrans/Disco/Raw/cellranger_count/')
Bam.first <- NULL
for (gse in GSE.first) {
  setwd('/data2/rluo4/EpiTrans/Disco/Raw/cellranger_count/')
  f <- list.files(gse)
  bam <- list.files(paste0('/data2/rluo4/EpiTrans/Disco/Raw/cellranger_count/',gse, '/',f, '/outs/'))
  bam <- bam[grepl('possorted_genome_bam.bam', bam) & !grepl('.bai', bam)]
  bam <- paste0(gse, '_', f, '_',bam)
  print(f)
  Bam.first <<- c(Bam.first, bam)
}
Bam.first <- as.data.frame(Bam.first)
Bam.first$GSE <- str_split(Bam.first$Bam.first, '_', simplify = T)[,1]
Bam.first$GSM <- str_split(Bam.first$Bam.first, '_', simplify = T)[,2]
raw_done <- unique(Bam.first$GSM)
table(unique(asp_sc_last$sample_alias) %in% raw_done)
# table(unique(BM_MAST$patient_id) %in% raw_done)
diff <- setdiff(unique(asp_sc$study_alias), unique(Bam.first$GSE))
raw_undone <- setdiff(diff, asp_sc_first$V1)
# [1] "GSE135194" "GSE139369" "GSE159624" "GSE159677" "GSE162025"
# [6] "GSE165645" "GSE166352" "GSE169426" "GSE181294" "GSE181989"
# [11] "GSE196756" "GSE201425" "GSE203115" "GSE219210" "GSE220116"
# [16] "GSE220243"
View(sc_Diseased_Samples[sc_Diseased_Samples$project_id %in% raw_undone, ])
asp_sc_last <- asp_sc_last[asp_sc_last$study_alias %in% disco_BM_MAST$project_id, ]
unique(asp_sc_last$sample_alias) #217-->440-->469--541
asp_sc_last <- asp_sc_last[ ! asp_sc_last$sample_alias %in% raw_done, ]
unique(asp_sc_last$sample_alias) #180-->309-->334
table(unique(asp_sc_last$sample_alias) %in% raw_done_AML)
table(unique(asp_sc_last$sample_alias) %in% unique(asp_sc_last_disease$V2))
# table(unique(asp_sc_last$sample_alias) %in% unique(sc_cell_summary$patient_id))
table(unique(asp_sc_last$sample_alias) %in% unique(disco$patient_id))
unique(asp_sc_last$study_alias) #22-->42-->43
setdiff(BM_MAST_newdatasets, unique(asp_sc_last$study_alias))
# [1] "GSE231535" --> character(0) --> "GSE185381"
setdiff(unique(asp_sc_last$study_alias), BM_MAST_newdatasets)
# [1] "GSE173706" "GSE181688" -->35 -- >43
# [1] "GSE181294" "GSE185965"
# write.table(asp_sc_last,file="/data2/rluo4/EpiTrans/DataCollection/asp_disco_scRNA_BM_MAST.link.txt",sep = "\t",quote = F, col.names = F, row.names = F)
# summary(1:689)
asp_sc_last_renew39 <- read.table("/data2/rluo4/EpiTrans/DataCollection/asp_disco_scRNA_BM_MAST.link.txt") # part1: new samples from 39 datasets
table(asp_sc_last$sample_alias %in% asp_sc_last_renew39$V2)
asp_sc_last <- asp_sc_last[ ! asp_sc_last$sample_alias %in% asp_sc_last_renew39$V2, ] # part1: new samples from 72-39=33 datasets
# write.table(asp_sc_last,file="/data2/rluo4/EpiTrans/DataCollection/asp_disco_scRNA_BM_MAST_part2.link.txt",sep = "\t",quote = F, col.names = F, row.names = F)
asp_sc_last_renew72 <- read.table("/data2/rluo4/EpiTrans/DataCollection/asp_disco_scRNA_BM_MAST_part2.link.txt") # part2: new samples from 72 datasets
table(asp_sc_last$sample_alias %in% asp_sc_last_renew72$V2)
asp_sc_last <- asp_sc_last[ ! asp_sc_last$sample_alias %in% asp_sc_last_renew72$V2, ] # part1: new samples from 72-39=33 datasets
# write.table(asp_sc_last,file="/data2/rluo4/EpiTrans/DataCollection/asp_disco_scRNA_BM_MAST_part3.link.txt",sep = "\t",quote = F, col.names = F, row.names = F)
asp_sc_last_renew74 <- read.table("/data2/rluo4/EpiTrans/DataCollection/asp_disco_scRNA_BM_MAST_part3.link.txt") # part2: new samples from 72 datasets
table(asp_sc_last$sample_alias %in% asp_sc_last_renew74$V2)
asp_sc_last <- asp_sc_last[ ! asp_sc_last$sample_alias %in% asp_sc_last_renew74$V2, ] # 
unique(asp_sc_last$study_alias)
# [1] "GSE139324" "GSE162025"
# write.table(asp_sc_last,file="/data2/rluo4/EpiTrans/DataCollection/asp_disco_scRNA_BM_MAST_part4.link.txt",sep = "\t",quote = F, col.names = F, row.names = F)
unique(asp_sc_last$study_alias)
unique(asp_sc_last$sample_alias)
summary(1:nrow(asp_sc_last))
# table(unique(asp_sc_last$study_alias) %in% BM_MAST$project_id)
table(unique(asp_sc_last$study_alias) %in% asp_sc_last_renew39$V1)
table(unique(asp_sc_last$study_alias) %in% sc_majortype$project_id)
table(unique(asp_sc_last_renew39$V1) %in% sc_majortype$project_id)


################################################################################
# 7. Write down cell_barcode_annotations.tsv for each disco dataset
################################################################################
################################################################################
# 7.1. Organization and adjustment of scRNA datasets: disco + PCTanno
# Skip below since it's time-consuming:
################################################################################
# (1) Adjust sub_type in sc_majortype of diso & PCTanno
################################################################################
disco <- readRDS(paste0('/data2/rluo4/RPMfunc/disco_pdata/disco_adjustmeta_0827.rds'))#0809.rds'))
disco_majortype <- readRDS('/data2/rluo4/RPMfunc/Output/scRNA/disco_all_adjustmeta_0827.rds')#0809.rds') # made by DiscoDatasets.R:
disco_pheno <- as.data.frame(table(disco_majortype$sample_type, disco_majortype$disease))
disco_pheno <- disco_pheno[disco_pheno$Freq!=0,]
# sample_type inconsistent with disease: GSE168453 (Critical/Moderate -- Healthy)
# sample.index <- c(grep("ormal", disco_pheno$Var1),grep("ealthy", disco_pheno$Var1),
#                   grep("control", disco_pheno$Var1))
# View(disco_pheno[sample.index,])
# disco changes of disease column
sample.index <- grepl("ormal|ealthy|control", disco_majortype$sample_type)
table(disco_majortype$disease[sample.index])
disco_majortype$disease[sample.index] <- 'Healthy'
covid.index <- disco_majortype$project_id == 'GSE168453'
table(disco_majortype$disease[covid.index])
table(disco_majortype$disease[covid.index], disco_majortype$sample_type[covid.index])
disco_majortype$disease[covid.index] <- disco_majortype$sample_type[covid.index]# critical, moderate and severe are COVID negative with NIH clinical symptomes
GSE168453_meta <- read.csv('/data2/rluo4/EpiTrans/DataCollection/GSE168453_covars.csv')
# GSE168453$NIH_clinical <- GSE168453_meta$NIH_clinical[match(GSE168453$patient_id, GSE168453_meta$donor)]
# GSE168453$COVID_status <- GSE168453_meta$COVID_status[match(GSE168453$patient_id, GSE168453_meta$donor)]
# GSE168453$pulmonary_infection <- GSE168453_meta$pulmonary_infection[match(GSE168453$patient_id, GSE168453_meta$donor)]
# GSE168453_healty <- read.table('/data2/rluo4/EpiTrans/DataCollection/COVID_healthy.txt')
# GSE168453_healty$COVID_status <- GSE168453_meta$COVID_status[match(GSE168453_healty$V1, GSE168453_meta$donor)]
# GSE168453_healty$NIH_clinical <- GSE168453_meta$NIH_clinical[match(GSE168453_healty$V1, GSE168453_meta$donor)]
# GSE168453_healty$sample_type <- GSE168453$sample_type[match(GSE168453_healty$V1, GSE168453$patient_id)]
# # write.table(GSE168453_disease$COVID_disease,file = '/data2/rluo4/EpiTrans/DataCollection/COVID_disease.txt',
# #             col.names = F, row.names = F, quote = F)
# GSE168453_disease <- read.table('/data2/rluo4/EpiTrans/DataCollection/COVID_disease.txt')
# GSE168453_disease$COVID_status <- GSE168453_meta$COVID_status[match(GSE168453_disease$V1, GSE168453_meta$donor)]
# GSE168453_disease$NIH_clinical <- GSE168453_meta$NIH_clinical[match(GSE168453_disease$V1, GSE168453_meta$donor)]
# GSE168453_disease$sample_type <- GSE168453$sample_type[match(GSE168453_disease$V1, GSE168453$patient_id)]
CT <- as.data.frame(table(disco_majortype$ct, disco_majortype$sub_type))
CT <- CT[CT$Freq!=0,]
# View(CT[duplicated(CT$Var1),])
index1 <- grepl('eutrophil',disco_majortype$ct)
disco_majortype$sub_type[index1] <- 'Neutrophils'
index2 <- grepl('EMT',disco_majortype$ct)
disco_majortype$sub_type[index2] <- 'Tissue_Specific_Cells'
index3 <- grepl('NKT',disco_majortype$ct)
disco_majortype$sub_type[index3] <- 'Tcells'
index4 <- grepl(' EC',disco_majortype$ct) # from BM_MAST_meta_p1
disco_majortype$sub_type[index4] <- 'Endothelial_Cells'
index5 <- grepl('Granulocyte|Myelocyte',disco_majortype$ct)
disco_majortype$sub_type[index5] <- 'OtherImmunecells'
# Granulocyte <- disco_majortype[disco_majortype$ct=='Granulocyte',]
# table(Granulocyte$project_id, Granulocyte$sub_type)#GSE185381--BM
# table(Granulocyte$project_id, Granulocyte$tissue)#GSE185381--BM
# HSC <- disco_majortype[disco_majortype$ct=='HSC',]
# table(HSC$project_id, HSC$sub_type)#GSE181989--BM
# table(HSC$project_id, HSC$tissue)#GSE181989--BM
# Placental mesenchymal stem/stromal cells (PMSCs) are fibroblast-like, self-renewing, and multipotent cells that reside in the chorionic plate, umbilical cord, and amnion of the placenta
index6 <- grepl('HSC',disco_majortype$ct) & disco_majortype$tissue %in% c('blood','bone marrow', 'Buffy coat','PBMC','thymus')
disco_majortype$sub_type[index6] <- 'Stem_Cells'
index7 <- grepl('myeloid',disco_majortype$ct)
disco_majortype$sub_type[index7] <- 'Myeloids'
index8 <- grepl(' T/NK',disco_majortype$ct)
disco_majortype$sub_type[index8] <- 'NKcells'
index9 <- grepl('Mast',disco_majortype$ct)
disco_majortype$sub_type[index9] <- 'MastCell'
index10 <- grepl('Oligodendrocyte',disco_majortype$ct)
disco_majortype$sub_type[index10] <- 'Neuro_Cells'
index11 <- grepl('epithelia',disco_majortype$ct)
disco_majortype$sub_type[index11] <- 'Tissue_Specific_Cells'
# saveRDS(disco_majortype,  file=paste0('/data2/rluo4/RPMfunc/Output/scRNA/disco_majortype_0828.rds'))

PCT_majortype <- readRDS('/data2/rluo4/RPMfunc/Output/scRNA/PCT_all_meta.rds') # made by Haixia
table(unique(PCT_majortype$sample_name) %in% tissue_summary$Cohort)
PCT_majortype$project_id <- tissue_summary$Dataset[match(PCT_majortype$sample_name, tissue_summary$Cohort)]
celltype_summary <- read_excel('/data2/rluo4/All/Output/table_xlsx/TableS1.xlsx')
table(PCT_majortype$Cell_Type %in% celltype_summary$Minor.Celltype)
setdiff(PCT_majortype$Cell_Type, celltype_summary$Minor.Celltype)#  "MSC.ADIPO" "MSC.MVA"   "MSC.SEC"   "MSC.PVA"
index <- match(PCT_majortype$Cell_Type, celltype_summary$Minor.Celltype)
PCT_majortype$FullName <- celltype_summary$FullName[index]
table(is.na(PCT_majortype$FullName))
unique(PCT_majortype$Cell_Type[is.na(PCT_majortype$FullName)])
PCT_majortype$ct[ is.na(PCT_majortype$FullName) ] <- 'Mesenchymal stem cells'
PCT_majortype$ct[! is.na(PCT_majortype$FullName) ] <- PCT_majortype$FullName[! is.na(PCT_majortype$FullName) ]
unique(PCT_majortype$ct)
unique(PCT_majortype$sub_type)
unique(disco_majortype$sub_type)
index1 <- PCT_majortype$Major_type=='Epithelial'
PCT_majortype$sub_type[index1] <- 'Tissue_Specific_Cells'
index2 <- PCT_majortype$Cell_Type=='MDSC'
PCT_majortype$sub_type[index2] <- 'Myeloids'
index3 <- PCT_majortype$Cell_Type=='STM'
PCT_majortype$sub_type[index3] <- 'Stem_Cells'
index4 <- grepl('Mast',PCT_majortype$ct)
PCT_majortype$sub_type[index4] <- 'MastCell'
index5 <- grepl('endritic',PCT_majortype$ct)
PCT_majortype$sub_type[index5] <- 'OtherImmunecells'
index6 <- grepl('uppfer',PCT_majortype$ct)
PCT_majortype$sub_type[index6] <- 'Macrophages'
index7 <- grepl('Megakaryocytes|Erythrocytes',PCT_majortype$ct)
PCT_majortype$sub_type[index7] <- 'Blood_Cells'
CT <- as.data.frame(table(PCT_majortype$ct, PCT_majortype$sub_type))
CT <- CT[CT$Freq!=0,]
# View(CT[duplicated(CT$Var1),])
PCT_majortype$Major_type <- gsub('Epithelial','Tissue_Specific_Cells',
                                 gsub('Stromal','Mesenchymal_Cells',
                                      gsub('Immune','Immune_Cells',PCT_majortype$Major_type)))
table(PCT_majortype$project_id %in% tissue_summary$Dataset)
PCT_majortype$platform <- tissue_summary$Platform[match(PCT_majortype$project_id, tissue_summary$Dataset)]
# saveRDS(PCT_majortype[, -10], file = '/data2/rluo4/RPMfunc/Output/scRNA/PCT_all_adjustmeta_0715.rds')
# write_rds() will save a big file.
#Langerhans cell histiocytosis (LCH) cells are a type of dendritic cell that normally helps the body fight infection. Sometimes mutations (changes) develop in genes that control how dendritic cells function. These include mutations of the BRAF, MAP2K1, RAS, and ARAF genes. These mutations may cause too many LCH cells to grow and build up in certain parts of the body, where they can damage tissue or form lesions.
colnames(disco_majortype)
# disco_majortype <- disco_majortype[, c(1:4,5,7:10)]
sc_majortype <- PCT_majortype %>% dplyr::select(barcode, project_id, patient_id = orig.ident, ct ,
                                                tissue = Organ, platform, sample_type = Tissue, disease = Tissue, Major_type, sub_type)
# sc_majortype$project_id <- PCTanno$project_id[match(sc_majortype$patient_id, PCTanno$patient_id)]
sc_majortype <- rbind(sc_majortype, disco_majortype[! disco_majortype$project_id %in% sc_majortype$project_id,])
unique(sc_majortype$sub_type) # 17
unique(sc_majortype$project_id) # 340 --> 312 --> 296
# sc_majortype$sub_type[sc_majortype$sub_type=='Fibroblast'] = 'Fibroblasts'
# sc_majortype$sub_type[sc_majortype$sub_type=='Macrophage'] = 'Macrophages'
table(sc_majortype$sub_type)
table(sc_majortype$tissue)
unique(sc_majortype$tissue) #79
CT <- as.data.frame(table(sc_majortype$ct, sc_majortype$sub_type))
CT <- CT[CT$Freq!=0,]
# View(CT[duplicated(CT$Var1),])
################################################################################
# (2) Adjust tissue and disease in sc_majortype of diso & PCTanno
################################################################################
sc_majortype$tissue <-  gsub('Acute meyloid leukemia|Acute leukemia|Acute monocytic leukemia|Acute myeloid leukemia|Chronic lymphocytic leukemia|Chronic myeloid leukemia|Leukemia|Myelogenous leukemia', 'bone marrow',
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
                                                                                                                                                     sc_majortype$tissue,
                                                                                                                                                )))))))))))))))))))))))))
sc_majortype$tissue <- tolower(sc_majortype$tissue)
table(sc_majortype$tissue) # 72 kinds of tissues in scRNA datasets
sc_majortype$tissue <- gsub('head and neck', "mouth",sc_majortype$tissue)# all oral cavity --> head and neck --> mouth, sc_majortype$tissue <- gsub("oral cavity", 'mouth',sc_majortype$tissue)
sc_majortype$tissue <- gsub("head", 'head and neck',sc_majortype$tissue) # sc_majortype$tissue <- gsub("teeth", 'gingiva',sc_majortype$tissue)

library(dplyr)
# Create a named vector with abbreviations and their corresponding full names (PCTanno)
abbreviations <- c(
  ADJ = "Adjacent tumor/disease section",
  IDC = "Invasive ductal carcinoma",
  DCIS = "Ductal carcinoma in situ",
  `BRCA1-mut` = "Precancerous lesion from BRCA1 mutation carriers",
  CC = "Cervix cancer",
  N_HPV = "HPV-infected normal cervix",
  HSIL_HPV = "HPV-infected high-grade squamous intraepithelial lesions",
  FAP = "Familial adenomatous polyposis",
  CRC = "Colorectal cancer",
  AD = "Adenomas",
  SER = "Sessile serrated lesions",
  `MSI-H` = "Microsatellite-high colorectal cancer",
  MSS = "Microsatellite stable colorectal cancer",
  AEH = "Atypical endometrial hyperplasia",
  EEC = "Endometrioid cancer",
  HGIN = "High-grade intraepithelial neoplasias",
  LGIN = "High-grade intraepithelial neoplasias",
  ESCC = "Esophageal squamous cell carcinoma",
  HCC = "Hepatocellular carcinoma",
  NAFLD = "Non-alcoholic fatty liver disease",
  AAH = "Atypical adenomatous hyperplasia",
  AIS = "Adenocarcinoma in situ",
  MIAC = "Minimally invasive adenocarcinoma",
  IAC = "Invasive lung adenocarcinoma",
  LP = "Leukoplakia",
  NEOLP = "Non-erosive oral lichen planus",
  EOLP = "Erosive oral lichen planus",
  OSCC = "Oral squamous cell carcinoma",
  PDAC = 'Pancreatic ductal adenocarcinoma',
  PanIN = 'Pancreatic Intraepithelial Neoplasia',
  BPH = "Benign prostatic hyperplasia",
  PRAD = "Prostate adenocarcinoma",
  AK = "Actinic keratosis",
  SCCIS = "Squamous cell carcinoma in situ",
  cSCC = "Cutaneous squamous cell carcinoma",
  CAG = "Chronic atrophic gastritis",
  `CAG with IM` = "Chronic atrophic gastritis with intestinal metaplasia",
  CSG = "Chronic superficial gastritis",
  GC = "Gastric cancer",
  SIM = "Severe intestinal metaplasia",
  WIM = "Wild intestinal metaplasia",
  PTC = "Papillary thyroid cancer",
  ATC = "Anaplastic thyroid cancer",
  HT = "Hashimoto's thyroiditis"
)
# Replace the abbreviations with full names
sc_majortype$disease <- recode(sc_majortype$disease, !!!abbreviations)
sc_majortype$disease <- gsub("AdjNorm", 'Adjacent tumor/disease section',
                             gsub('AML', 'Acute myeloid leukemia', gsub("’s", "'s", sc_majortype$disease)) )
sc_majortype$disease[sc_majortype$disease == "Infection"] <- 'COVID-19'
sc_majortype$disease[sc_majortype$project_id=='GSE145927'] <- 'acute antibody-mediated rejection' # kidney transplant
sc_majortype$disease <- tolower(sc_majortype$disease)
# Check the result
head(sc_majortype$disease)
sort(unique(sc_majortype$disease)) # ~ 154 kinds --> 142 -->173
# saveRDS(sc_majortype, file = '/data2/rluo4/RPMfunc/Output/scRNA/sc_majortype_0725.rds')#sc_majortype.rds')
# saveRDS(sc_majortype, file = '/data2/rluo4/RPMfunc/Output/scRNA/sc_majortype_0808.rds')# BM_MAST_p1
# saveRDS(sc_majortype, file = '/data2/rluo4/RPMfunc/Output/scRNA/sc_majortype_0809.rds')# BM_MAST_p1+p2+p3(GSE181294)
table(sc_majortype$project_id[is.na(sc_majortype$patient_id)])
# GSE181989 <- subset(sc_majortype, project_id == 'GSE181989')
# GSM <- str_split(rownames(GSE181989), '--', simplify = T)[,2]
index <- is.na(sc_majortype$patient_id)
sc_majortype$patient_id[index] <- str_split(rownames(sc_majortype[index, ]), '--', simplify = T)[,2]
table(sc_majortype$project_id[is.na(sc_majortype$patient_id)])
saveRDS(sc_majortype, file = '/data2/rluo4/RPMfunc/Output/scRNA/sc_majortype_0828.rds')# BM_MAST_p1+p2+p3(GSE181294) + last2 + wrong3 + last1

################################################################################
# sc_majortype <- readRDS(file = '/data2/rluo4/RPMfunc/Output/scRNA/sc_majortype_0828.rds')
table(unique(disco$project_id) %in% c('GSE211033','GSE196052','GSE196676')) # 3 yes
table(unique(sc_majortype$project_id) %in% c('GSE211033','GSE196052','GSE196676')) # all no
table(unique(BM_MAST_meta_p1$project_id) %in% c('GSE211033','GSE196052','GSE196676')) # all no
table(unique(disco_BM_MAST$project_id) %in% c('GSE211033','GSE196052','GSE196676')) # 3 yes
table(unique(BM_MAST_meta_p1$project_id) %in% unique(sc_majortype$project_id)) # all yes
table(unique(disco_majortype_add_p1$project_id) %in% unique(sc_majortype$project_id)) # all yes

table(na.asp_sc_file$sample_alias %in% Bam.first$GSM)
setdiff(na.asp_sc_file$sample_alias, Bam.first$GSM) #9 single-end sequencing samples?
################################################################################
# [1] "GSM4952959" "GSM4952960" "GSM6290522"
# [4] "GSM6422822" "GSM6422825" "GSM6422824"
# [7] "GSM6422823" "GSM6840122" "GSM6840117"
# (base) [rluo4@sbimrdplp001 GSE162498]$ ls -hltr
# total 38G
# -rw-r--r-- 1 rluo4 rluo4  11G Mar 30 06:29 SRR13181489
# -rw-r--r-- 1 rluo4 rluo4 9.3G Mar 30 09:49 SRR13181490
# -rw-rw-r-- 1 rluo4 rluo4 7.5G Apr  6 23:17 SRR13181490_1.fastq.gz
# -rw-rw-r-- 1 rluo4 rluo4  11G Apr  7 01:10 SRR13181489_1.fastq.gz
# (base) [rluo4@sbimrdplp001 GSE162498]$ cd ..
# (base) [rluo4@sbimrdplp001 Sra]$ cd GSE207493
# (base) [rluo4@sbimrdplp001 GSE207493]$ ls -hltr
# total 8.3G
# -rw-r--r-- 1 rluo4 rluo4 4.7G Apr  4 11:24 SRR19987215
# -rw-rw-r-- 1 rluo4 rluo4 3.7G Apr  7 08:30 SRR19987215_1.fastq.gz
# (base) [rluo4@sbimrdplp001 GSE207493]$ cd ..
# (base) [rluo4@sbimrdplp001 Sra]$ cd GSE210152
# (base) [rluo4@sbimrdplp001 GSE210152]$ ls -hltr
# total 44G
# -rw-r--r-- 1 rluo4 rluo4 4.5G Apr  3 19:17 SRR20723945
# -rw-r--r-- 1 rluo4 rluo4 2.1G Apr  3 20:19 SRR20723946
# -rw-r--r-- 1 rluo4 rluo4 1.8G Apr  3 20:49 SRR20723947
# -rw-rw-r-- 1 rluo4 rluo4 1.7G Apr  7 08:57 SRR20723947_1.fastq.gz
# -rw-rw-r-- 1 rluo4 rluo4 1.7G Apr  7 09:21 SRR20723944_1.fastq.gz
# -rw-rw-r-- 1 rluo4 rluo4 4.0G Apr  7 10:24 SRR20723945_1.fastq.gz
# -rw-rw-r-- 1 rluo4 rluo4 1.9G Apr  7 10:53 SRR20723946_1.fastq.gz
# -rw-rw-r-- 1 rluo4 rluo4  26G Apr  8 17:47 SRR20723944
# (base) [rluo4@sbimrdplp001 GSE210152]$ cd ..
# (base) [rluo4@sbimrdplp001 Sra]$ cd GSE220116
# (base) [rluo4@sbimrdplp001 GSE220116]$ ls -hltr
# total 20G
# -rw-r--r-- 1 rluo4 rluo4 6.9G Apr  2 12:41 SRR22746967
# -rw-r--r-- 1 rluo4 rluo4 3.4G Apr  2 14:14 SRR22746972
# -rw-rw-r-- 1 rluo4 rluo4 6.3G Apr  7 13:43 SRR22746967_1.fastq.gz
# -rw-rw-r-- 1 rluo4 rluo4 2.9G Apr  7 14:29 SRR22746972_1.fastq.gz
# (base) [rluo4@sbimrdplp001 GSE220116]$
# load('DiscoMut.R')
################################################################################
index = (asp_sc_file$sample_alias %in% Bam.first$GSM) & is.na(asp_sc_file$filesize)
# View(asp_sc_file[index,])
transformed_files <- asp_sc_file[index, ]
asp_sc_file$filesize[index] <- 25301137945*2 
index = index & (asp_sc_file$datatype!='Bam')
asp_sc_file$filename[index] <- paste0(asp_sc_file$sample_alias[index],'_S1_L001_R2_001.fastq.gz')
  
table(duplicated(asp_sc_file$sample_alias))
rownames(asp_sc_file) = 1:nrow(asp_sc_file)
tmp=by(asp_sc_file,asp_sc_file$sample_alias, function(x) rownames(x)[which.max(x$filesize)])
probes = as.character(tmp)
table(is.na(probes))
table(rownames(asp_sc_file) %in% probes)
asp_sc_file=asp_sc_file[rownames(asp_sc_file) %in% probes ,]
asp_sc_file$path = rownames(file_info_allDat)[match(asp_sc_file$filename, file_info_allDat$filename)]
asp_sc_file$path[is.na(asp_sc_file$path)] = asp_sc_file$filename[is.na(asp_sc_file$path)]

table(asp_sc_file$sample_alias %in% disco_allDat$sample_alias)
table(asp_sc_file$sample_alias %in% unique(disco_majortype_add_p1$patient_id))
table(unique(asp_sc_file$study_alias) %in% unique(disco_majortype_add_p1$project_id))
table(asp_sc_file$sample_alias %in% unique(BM_MAST_meta_p1$sample_id))
table(unique(asp_sc_file$study_alias) %in% unique(BM_MAST_meta_p1$project_id))
# FALSE  TRUE 
# 60    29 
# FALSE  TRUE 
# 64    30
# write.table(asp_sc_file[asp_sc_file$study_alias %in% unique(BM_MAST_meta_p1$project_id),], 
#             file="/data2/rluo4/EpiTrans/DataCollection/asp_sc_file_transformed_BM_MAST_p1.txt",
#             sep = "\t", quote = F, col.names = F, row.names = F)
table(unique(asp_sc_file$study_alias) %in% unique(BM_MAST_meta_p2$project_id))
# FALSE  TRUE 
# 74    18 
asp_sc_file_BM_MAST_p2 <- asp_sc_file[asp_sc_file$study_alias %in% unique(BM_MAST_meta_p2$project_id),]
# write.table(asp_sc_file_BM_MAST_p2, 
#             file="/data2/rluo4/EpiTrans/DataCollection/asp_sc_file_transformed_BM_MAST_p2.txt",
#             sep = "\t", quote = F, col.names = F, row.names = F)
asp_sc_file_BM_MAST_p4 <- asp_sc_file[asp_sc_file$study_alias %in% unique(wrong_3$project_id),]
# write.table(asp_sc_file_BM_MAST_p4, 
#             file="/data2/rluo4/EpiTrans/DataCollection/asp_sc_file_transformed_BM_MAST_p4.txt",
#             sep = "\t", quote = F, col.names = F, row.names = F)

asp_sc_file_BM_MAST_p3 <- asp_sc_file[asp_sc_file$study_alias %in% unique(missed_data2$project_id),]
# write.table(asp_sc_file_BM_MAST_p3,
#             file="/data2/rluo4/EpiTrans/DataCollection/asp_sc_file_transformed_BM_MAST_p3.txt",
#             sep = "\t", quote = F, col.names = F, row.names = F)
# # write.table(asp_sc_file, file="/data2/rluo4/EpiTrans/DataCollection/asp_sc_file_transformed_disease_rest.txt",sep = "\t",quote = F, col.names = F, row.names = F)
# # write.table(asp_sc_file, file="/data2/rluo4/EpiTrans/DataCollection/asp_sc_file_transformed_disease.txt",sep = "\t",quote = F, col.names = F, row.names = F)
# # 68 datasets -- 297 samples (-rw-r--r-- 1 rluo4 rluo4  64K Jul 25 16:03 asp_sc_file_transformed_disease.txt)
# # asp_sc_file <- read.table(file="/data2/rluo4/EpiTrans/DataCollection/asp_sc_file_transformed_disease.txt", header = F)
# # Disco_SComatic_disease <- asp_sc_file[asp_sc_file[,1] %in% AML$project_id,]
# # paste0(Disco_SComatic_disease[,2], collapse = '|')
asp_sc_file_BM_MAST_p1 <- read.table(file="/data2/rluo4/EpiTrans/DataCollection/asp_sc_file_transformed_BM_MAST_p1.txt", header = F)
unique(asp_sc_file_BM_MAST_p1$V1)
asp_sc_file_BM_MAST_p1_missed <- asp_sc_file[asp_sc_file$study_alias %in% unique(BM_MAST_meta_p1$project_id),]
unique(asp_sc_file_BM_MAST_p1_missed$study_alias)
setdiff(unique(asp_sc_file_BM_MAST_p1_missed$study_alias), unique(asp_sc_file_BM_MAST_p1$V1))
#  "GSE135194" "GSE139369" "GSE203115"
index <- ! asp_sc_file_BM_MAST_p1_missed$study_alias %in% unique(asp_sc_file_BM_MAST_p1$V1)
write.table(asp_sc_file_BM_MAST_p1_missed[index, ],
            file="/data2/rluo4/EpiTrans/DataCollection/asp_sc_file_transformed_BM_MAST_p1_missed3.txt",
            sep = "\t", quote = F, col.names = F, row.names = F)

table(asp_sc_file_BM_MAST_p2$study_alias %in% asp_sc_file_BM_MAST_p1$V1)
asp_sc_file_BM_MAST_GSE130560 <- asp_sc_file[asp_sc_file$study_alias %in% 'GSE130560',]
# write.table(asp_sc_file_BM_MAST_GSE130560,
#             file="/data2/rluo4/EpiTrans/DataCollection/asp_sc_file_transformed_BM_MAST_GSE130560.txt",
#             sep = "\t", quote = F, col.names = F, row.names = F)
asp_sc_file_BM_MAST_p2 <- read.table(file="/data2/rluo4/EpiTrans/DataCollection/asp_sc_file_transformed_BM_MAST_p2.txt", header = F)
asp_sc_file_BM_MAST_p3 <- read.table(file="/data2/rluo4/EpiTrans/DataCollection/asp_sc_file_transformed_BM_MAST_p3.txt", header = F)
asp_sc_file_BM_MAST_p4 <- read.table(file="/data2/rluo4/EpiTrans/DataCollection/asp_sc_file_transformed_BM_MAST_p4.txt", header = F)
table(raw_undone %in% asp_sc_file_BM_MAST_p1$V1)
table(raw_undone %in% asp_sc_file_BM_MAST_p1_missed$study_alias)
table(raw_undone %in% asp_sc_file_BM_MAST_p2$V1)
table(raw_undone %in% asp_sc_file_BM_MAST_p3$V1)
table(raw_undone %in% asp_sc_file_BM_MAST_p4$V1)

# delete bone marrow samples except for those from 3 AML datasets
# write.table(asp_sc_file, file="/data2/rluo4/EpiTrans/DataCollection/asp_sc_file_transformed_AML.txt",sep = "\t",quote = F, col.names = F, row.names = F)
# 69 datasets -- 254 samples (-rw-r--r-- 1 rluo4 rluo4 54K Jul 21 16:06 asp_sc_file_transformed_AML.txt)
# asp_sc_file <- read.table(file="/data2/rluo4/EpiTrans/DataCollection/asp_sc_file_transformed_AML.txt", header = F)
# Disco_SComatic_AML <- asp_sc_file[asp_sc_file$V1 %in% AML$project_id,]
# paste0(Disco_SComatic_AML$V2, collapse = '|')
# table(asp_sc_file$V2 %in% disco_allDat$sample_alias)
# write.table(asp_sc_file, file="/data2/rluo4/EpiTrans/DataCollection/asp_sc_file_transformed.txt",sep = "\t",quote = F, col.names = F, row.names = F)
# 72 datasets -- 265 samples
# write.table(asp_sc_file, file="/data2/rluo4/EpiTrans/DataCollection/asp_sc_file.txt",sep = "\t",quote = F, col.names = F, row.names = F)
extract_first_element <- function(x) {
  split_string <- strsplit(x, "[--]")
  first_element <- sapply(split_string, function(y) head(y, n = 1))
  return(first_element)
}
datasets <- unique(disco$project_id)
datasets_AML <- unique(disco$project_id[disco$tissue=='bone marrow'])
datasets_PDAC <- unique(PDAC1$project_id)
datasets_BM_MAST <- unique(disco_majortype_add_p1$project_id)
table(datasets_BM_MAST %in% UTH57_wrong_mapping)
datasets_BM_MAST <- unique(disco_majortype_add_p2$project_id)
datasets_BM_MAST <- unique(disco_majortype_add_p3$project_id)
# datasets_BM_MAST <- unique(disco_majortype_add_p4$project_id)
for(i in datasets_BM_MAST){#datasets_PDAC){#datasets_AML # datasets){
  # disco_fisrt <- disco[disco$patient_id %in% asp_sc_file$sample_alias,]
  disco_meta <- disco[disco$project_id %in% i,] # for datasets, datasets_AML, datasets_PDAC
  disco_meta <- disco_majortype_add_p1[disco_majortype_add_p1$project_id %in% i,]

  disco_meta <- disco_majortype_add_p2[disco_majortype_add_p2$project_id %in% i,]
  disco_meta <- disco_majortype_add_p3[disco_majortype_add_p3$project_id %in% i,]
  # disco_meta <- disco_majortype_add_p4[disco_majortype_add_p4$project_id %in% i,]
  disco_meta <- AML1
  all_characters <- unique(unlist(strsplit(disco_meta$ct, "")))
  potential_separators_and_whitespace <- all_characters[grep("[[:punct:][:space:]]", all_characters)]
  # Define the replacement character
  replacement <- "_"
  # Replace each potential separator with the replacement character
  modified_vector <- disco_meta$ct
  for (sep in potential_separators_and_whitespace) {
    modified_vector <- gsub(sep, replacement, modified_vector, fixed = TRUE)
  }
  # Replace consecutive occurrences of the replacement character with a single underscore
  disco_meta$celltype <- gsub("_+", replacement, modified_vector)

  # Call the function to extract the last element from each string
  first_elements <- extract_first_element(rownames(disco_meta))
  i = 'GSE185381'
  first_elements <- extract_first_element(GSE185381$barcode)
  
  Index = first_elements
  Index = str_split(Index,'[-]',simplify = T)[,1]
  meta <- data.frame('Index' = Index,
                     'Cell_type'= disco_meta[,c('celltype')])
  print(table(meta$Index==''))
  print(table(is.na(meta$Index)))
  write.table(meta, paste0('/data2/rluo4/EpiTrans/Disco/MetaData/',i,'.cell_barcode_annotations.tsv'),sep = '\t',
              row.names = F,quote = F)
}
# for disco_majortype_add_p2
# disco_majortype <- readRDS(file=paste0('/data2/rluo4/RPMfunc/Output/scRNA/disco_all_adjustmeta_0809.rds'))
UTH57_scMapping <- list.files('/data2/rluo4/EpiTrans/Disco/scMapping/') #45
# UTH57_scMapping <- list.files('/data2/rluo4/EpiTrans/Disco/anno.var/') #57
UTH57_scMapping <- list.files('/data2/rluo4/EpiTrans/Disco/SCOMATIC/') #45
# save(UTH57_scMapping,  file=paste0('/data2/rluo4/RPMfunc/disco_pdata/UTH57_scMapping_p2.rds'))
UTH57_wrong_mapping <- intersect(UTH57_scMapping , inter_sets) #2 [1] "GSE130560" "GSE185965"
# rm /data2/rluo4/EpiTrans/Disco/scMapping/GSE185965 only
# for disco_majortype_add_p1
# disco_majortype <- readRDS(file=paste0('/data2/rluo4/RPMfunc/Output/scRNA/disco_all_adjustmeta_0808.rds'))
# the wrong metadata files was removed, skip !
# UTH57_scMapping <- list.files('/data2/rluo4/EpiTrans/Disco/scMapping/') #70
# UTH57_scMapping <- list.files('/data2/rluo4/EpiTrans/Disco/anno.var/') #70
# save(UTH57_scMapping,  file=paste0('/data2/rluo4/RPMfunc/disco_pdata/UTH57_scMapping.rds'))
# UTH57_wrong_mapping <- intersect(UTH57_scMapping , inter_sets) #26
UTH57_wrong_mapping <- unique(wrong_3$project_id)
# for (i in UTH57_wrong_mapping ) {
#   GSM <-  list.files(file.path("/data2/rluo4/EpiTrans/Disco/scMapping/",i ))
#   print(GSM)
#   for (j in GSM) {
#     print(paste('to be removed dataset: ', i))
#     print(paste('to be removed GSM: ', j))
#     scMapping_files <-  list.files(file.path("/data2/rluo4/EpiTrans/Disco/scMapping/",i,j ))
#     print(scMapping_files)
#   }
#   system(paste0('rm -r ', '/data2/rluo4/EpiTrans/Disco/scMapping/', i))
#   system(paste0('rm ', '/data2/rluo4/EpiTrans/Disco/MetaData/', i, '.cell_barcode_annotations.tsv'))
#   system(paste0('rm -r ', '/data2/rluo4/EpiTrans/Disco/SCOMATIC/', i))
#   system(paste0('rm -r ', '/data2/rluo4/EpiTrans/Disco/anno.var/', i))
# }
# sc_summary <- sc_majortype[! duplicated(paste(sc_majortype$patient_id, sc_majortype$tissue, 
#                                               sc_majortype$project_id)),]
# length(unique(sc_summary$project_id)) #325 --> 328
# sample.index <- c(grep("ormal", sc_summary$sample_type),grep("ealthy", sc_summary$sample_type),
#                   grep("ontrol", sc_summary$sample_type))
# sc_Diseased_Samples <- sc_summary[ -sample.index, c(2,3,5:8)] # 219 --> 222 disease datasets
# # unique(sc_majortype[! grepl('adjacent|ealthy', sc_majortype$disease),]$project_id)
# RMzyme_pheno <- as.data.frame(table(sc_majortype$sample_type, sc_majortype$disease))
# RMzyme_pheno <- RMzyme_pheno[RMzyme_pheno$Freq!=0,]
# # organize the diseased data of disco named as sc_cell_summary 
# sc_cell_summary <- sc_majortype[! grepl('adjacent|ealthy', sc_majortype$disease),] # 219 disease datasets
# sc_cell_summary <- sc_cell_summary[ ! sc_cell_summary$project_id %in% tissue_summary$Dataset, ] # 169 --> 172 disease datasets
# first_elements <- extract_first_element(sc_cell_summary$barcode)
# CB = str_split(first_elements,'[-]',simplify = T)[,1]
# sc_cell_summary$CB <- CB #str_split(sc_cell_summary$barcode, '-', simplify = T)[,1]
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
disco_BM_MAST <- readRDS(file=paste0('/data2/rluo4/RPMfunc/disco_pdata/BM_MAST_meta.rds'))
#
AML <- readRDS(file=paste0('/data2/rluo4/RPMfunc/disco_pdata/AML_meta.rds')) # made through /data2/rluo4/Rcode/RMDB/DiscoDatasets.R
#
PDAC1 <- readRDS(file=paste0('/data2/rluo4/RPMfunc/disco_pdata/PDAC1_meta.rds'))
# delete the wrong datasets (mast/basophils) from disco
pancreas_pid <- c('GSE205049','GSE217837','GSE231535')
excluded_data <- c(unique(AML$project_id), pancreas_pid, #unique(BM$project_id), 
                   unique(disco_BM_MAST$project_id)) #107 excluded datasets -->115
sort(excluded_data)
excluded_data <- data.frame(excluded_datasets = excluded_data, source = c(rep('AML', length(unique(AML$project_id))),
                                                                          rep('pancreas', length(unique(pancreas_pid))), #rep('BM', length(unique(BM$project_id))), 
                                                                          rep('BM_MAST', length(unique(disco_BM_MAST$project_id)))) # exclude 3 BM_MAST datasets c('GSE211033','GSE196052','GSE196676')
) 
renew_datasets <- intersect(excluded_data$excluded_datasets, asp_sc_last$study_alias)# all FastQ
intersect(asp_sc_last$study_alias, asp_sc_first$V1)
# [1] "GSE155468" "GSE152042" "GSE151192" "GSE150825" "GSE130560"
table(renew_datasets %in% disco_BM_MAST$project_id)
setdiff(renew_datasets, disco_BM_MAST$project_id)
# [1] "GSE185381" "GSE185991" "GSE211036"
asp_sc_renew <- asp_sc[ asp_sc$study_alias %in% renew_datasets,]
################################################################################
# 5.3 Organize the files for disco's SComatic: filtering the files that were done
################################################################################
asp_sc_file <- readRDS('/data2/rluo4/EpiTrans/DataCollection/asp_sc_file.rds')
length(unique(asp_sc_file$study_alias)) #127
# asp_sc_file <- read.table("/data2/rluo4/EpiTrans/DataCollection/asp_sc_file_all.txt", sep = '\t', header = F)
GSE.first <- list.files('/data2/rluo4/EpiTrans/Disco/Raw/cellranger_count/')
Bam.first <- NULL
for (gse in GSE.first) {
  setwd('/data2/rluo4/EpiTrans/Disco/Raw/cellranger_count/')
  f <- list.files(gse)
  bam <- list.files(paste0('/data2/rluo4/EpiTrans/Disco/Raw/cellranger_count/',gse, '/',f, '/outs/'))
  bam <- bam[grepl('possorted_genome_bam.bam', bam) & !grepl('.bai', bam)]
  bam <- paste0(gse, '_', f, '_',bam)
  print(f)
  Bam.first <<- c(Bam.first, bam)
}
Bam.first <- as.data.frame(Bam.first)
Bam.first$GSE <- str_split(Bam.first$Bam.first, '_', simplify = T)[,1]
Bam.first$GSM <- str_split(Bam.first$Bam.first, '_', simplify = T)[,2]
raw_done <- unique(Bam.first$GSM)
diff <- setdiff(unique(asp_sc$study_alias), unique(Bam.first$GSE))
raw_undone <- setdiff(diff, asp_sc_last$study_alias)
table(raw_undone %in% asp_sc_last$study_alias)

table(unique(asp_sc_renew$sample_alias) %in% raw_done)
# table(unique(BM_MAST$patient_id) %in% raw_done)
table(asp_sc_renew$study_alias %in% disco_BM_MAST$project_id)
# asp_sc_renew <- asp_sc_renew[asp_sc_renew$study_alias %in% disco_BM_MAST$project_id, ]
unique(asp_sc_renew$sample_alias) #
table( asp_sc_renew$sample_alias %in% raw_done)
# asp_sc_renew <- asp_sc_renew[ ! asp_sc_renew$sample_alias %in% raw_done, ]
unique(asp_sc_renew$sample_alias) #823
table(unique(asp_sc_renew$sample_alias) %in% unique(disco$patient_id))
unique(asp_sc_renew$study_alias) #59

index = (asp_sc_file$sample_alias %in% Bam.first$GSM) & is.na(asp_sc_file$filesize)
# View(asp_sc_file[index,])
transformed_files <- asp_sc_file[index, ]
length(unique(transformed_files$study_alias)) #89
asp_sc_file$filesize[index] <- 25301137945*2 
index = index & (asp_sc_file$datatype!='Bam')
asp_sc_file$filename[index] <- paste0(asp_sc_file$sample_alias[index],'_S1_L001_R2_001.fastq.gz')

table(duplicated(asp_sc_file$sample_alias))
rownames(asp_sc_file) = 1:nrow(asp_sc_file)
tmp=by(asp_sc_file,asp_sc_file$sample_alias, function(x) rownames(x)[which.max(x$filesize)])
probes = as.character(tmp)
table(is.na(probes))
table(rownames(asp_sc_file) %in% probes)
asp_sc_file=asp_sc_file[rownames(asp_sc_file) %in% probes ,]
asp_sc_file$path = rownames(file_info_allDat)[match(asp_sc_file$filename, file_info_allDat$filename)]
asp_sc_file$path[is.na(asp_sc_file$path)] = asp_sc_file$filename[is.na(asp_sc_file$path)]
unique(asp_sc_file$study_alias) # 94 --> 79
# write.table(asp_sc_file, file="/data2/rluo4/EpiTrans/DataCollection/asp_sc_file_renew.txt",sep = "\t",quote = F, col.names = F, row.names = F)

sc_cell_summary <- readRDS("/data2/rluo4/RPMfunc/Output/scRNA/sc_cell_summary.rds")
Diseased_Samples_UTH36 <- read.table('/data2/rluo4/RPMfunc/disco_pdata/Diseased_Samples_UTH36.txt')
# Base directory and output directory
base_dir <- '/data2/rluo4/RPMfunc/SCOMATIC'
out_dir <- '/data2/rluo4/RPMfunc/Output/scRNA/disco_SNV' # 88-1 datasets 
# Get unique cohorts
combined_cohort <- gsub('_SCOMATIC_annovar.RData','',list.files(out_dir))
setdiff( unique(Diseased_Samples_UTH36$project_id), combined_cohort)
# [1] "GSE140042" "GSE156405" were excluded in SComatic analysis

UTH57_scMapping <- list.files('/data2/rluo4/EpiTrans/Disco/scMapping/')
GSM_scMapping <- NULL
for (i in UTH57_scMapping) {
  j <- list.files(file.path('/data2/rluo4/EpiTrans/Disco/scMapping/', i))
  GSM_scMapping <<- c(GSM_scMapping, j)
}
GSM_scMapping <- c(GSM_scMapping, Diseased_Samples_UTH36$patient_id)
length(unique(GSM_scMapping))
raw_undone <- setdiff(asp_sc_file$sample_alias, GSM_scMapping)
# View(asp_sc_file[asp_sc_file$sample_alias %in% raw_undone,])
table(raw_undone %in% unique(sc_cell_summary$patient_id))

scMapping_undone <- sc_Diseased_Samples[sc_Diseased_Samples$patient_id %in% raw_undone,]
table(scMapping_undone$patient_id %in% sc_cell_summary$patient_id)

index <- asp_sc_file_BM_MAST_p1$V2 %in% raw_undone
asp_sc_file_BM_MAST_p1_undone <- asp_sc_file_BM_MAST_p1[index, ]
# write.table(asp_sc_file_BM_MAST_p1_undone,
#             file="/data2/rluo4/EpiTrans/DataCollection/asp_sc_file_transformed_BM_MAST_p1_undone.txt",
#             sep = "\t", quote = F, col.names = F, row.names = F)
index <- asp_sc_file_BM_MAST_p2$V2 %in% raw_undone
asp_sc_file_BM_MAST_p2_undone <- asp_sc_file_BM_MAST_p2[index, ]
# write.table(asp_sc_file_BM_MAST_p2_undone,
#             file="/data2/rluo4/EpiTrans/DataCollection/asp_sc_file_transformed_BM_MAST_p2_undone.txt",
#             sep = "\t", quote = F, col.names = F, row.names = F)
index <- asp_sc_file_BM_MAST_p3$V2 %in% raw_undone
asp_sc_file_BM_MAST_p3_undone <- asp_sc_file_BM_MAST_p3[index, ]
# write.table(asp_sc_file_BM_MAST_p3_undone,
#             file="/data2/rluo4/EpiTrans/DataCollection/asp_sc_file_transformed_BM_MAST_p3_undone.txt",
#             sep = "\t", quote = F, col.names = F, row.names = F)
index <- asp_sc_file_BM_MAST_p4$V2 %in% raw_undone
asp_sc_file_BM_MAST_p4_undone <- asp_sc_file_BM_MAST_p4[index, ]
# write.table(asp_sc_file_BM_MAST_p4_undone,
#             file="/data2/rluo4/EpiTrans/DataCollection/asp_sc_file_transformed_BM_MAST_p4_undone.txt",
#             sep = "\t", quote = F, col.names = F, row.names = F)

asp_sc_file_BM_MAST <- rbind(asp_sc_file_BM_MAST_p1, asp_sc_file_BM_MAST_p2, asp_sc_file_BM_MAST_p3, asp_sc_file_BM_MAST_p4)
table(duplicated(asp_sc_file_BM_MAST$V1))
# write.table(asp_sc_file_BM_MAST,
#             file="/data2/rluo4/EpiTrans/DataCollection/asp_sc_file_transformed_BM_MAST_undone.txt",
#             sep = "\t", quote = F, col.names = F, row.names = F)
asp_sc_file_undone_v1 <- asp_sc_file[asp_sc_file$sample_alias %in% scMapping_undone$patient_id,]
table(is.na(asp_sc_file_undone_v1$filesize))
table(grepl('/data2/', asp_sc_file_undone_v1$path))
asp_sc_file_undone_v1 <- asp_sc_file_undone_v1[grepl('/data2/', asp_sc_file_undone_v1$path),]
write.table(asp_sc_file_undone_v1,
            file="/data2/rluo4/EpiTrans/DataCollection/asp_sc_file_undone.txt",
            sep = "\t", quote = F, col.names = F, row.names = F)
# asp_sc_file <- fread("/data2/rluo4/EpiTrans/DataCollection/asp_sc_file_transformed.txt")
table(asp_sc_file_BM_MAST$V2 %in% asp_sc_file$sample_alias)
setdiff(asp_sc_file_BM_MAST$V2, asp_sc_file$sample_alias)

asp_sc_file <- readRDS('/data2/rluo4/EpiTrans/DataCollection/asp_sc_file.rds')
# asp_sc_file_undone <- asp_sc_file[asp_sc_file$sample_alias %in% scMapping_undone$patient_id,]
# asp_sc_file_undone$path = rownames(file_info_allDat)[match(asp_sc_file_undone$filename, file_info_allDat$filename)]
# asp_sc_file_undone$path[is.na(asp_sc_file_undone$path)] = asp_sc_file_undone$filename[is.na(asp_sc_file_undone$path)]
# table(is.na(asp_sc_file_undone$filesize))
# table(grepl('/data2/', asp_sc_file_undone$path))
# asp_sc_file_undone <- asp_sc_file_undone[grepl('/data2/', asp_sc_file_undone$path),]
# tmp=by(asp_sc_file_undone,asp_sc_file_undone$sample_alias, function(x) rownames(x)[which.max(x$filesize)])
# probes = as.character(tmp)
# table(is.na(probes))
# table(rownames(asp_sc_file_undone) %in% probes)
# asp_sc_file_undone=asp_sc_file_undone[rownames(asp_sc_file_undone) %in% probes ,]
# unique(asp_sc_file$study_alias)

# update in 11/27/2024
studies_done <- c(unique(Diseased_Samples_UTH36$project_id), UTH57_scMapping) # UTH57_scMapping: 84 cohorts
setdiff(studies_done, sc_Diseased_Samples$project_id)
studies_done <- studies_done[studies_done %in% sc_Diseased_Samples$project_id]
studies_running <- c('GSE185344', 'GSE162025', 'GSE221156', 'GSE211033', 'GSE207493', 
                     "GSE196638", "GSE203191", "GSE212447", "GSE185965")#'GSE162025' (7 sample undone!) , 'GSE185965'(2 samples undone!); 'GSE210543'(done in 11/30/2024)
studies_done <- studies_done[ ! studies_done %in% studies_running]
setdiff(studies_done, asp_sc_file$study_alias) # [1] "GSE121638" "GSE131882" "GSE162500"
# asp_sc_file_done <- asp_sc_file[asp_sc_file$study_alias %in% studies_done,]
asp_sc_file_done_renew <- asp_sc_file[asp_sc_file$study_alias %in% studies_done,]
unique(asp_sc_file_done_renew$study_alias) # 87 --> 89
# asp_sc_file_done <- asp_sc_file_done[!asp_sc_file_done$study_alias %in% unique(Diseased_Samples_UTH36$project_id),]
write.table(asp_sc_file_done_renew, file="/data2/rluo4/EpiTrans/DataCollection/asp_sc_file_done_renew1.txt",sep = "\t",quote = F, col.names = F, row.names = F)


