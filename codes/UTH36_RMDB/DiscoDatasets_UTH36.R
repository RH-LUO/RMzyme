#!/usr/local/bin/Rscript
# title: "DISCO: a database of Deeply Integrated human Single-Cell Omics data"
# title:
# author: "Ruihan Luo"
# date: "March 25th,2024"
# rm(list=ls())
data_dir <- '/data/rluo4/EpiTrans/DataCollection/disco_tsv'
setwd(data_dir)
disco_original <- readRDS('/data/rluo4/RPMfunc/xuhaixia/disco_all_meta.rds')
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
tissue_na.index = disco$tissue==''
# tissue_na <- disco[disco$tissue=='',]
#  Cross-tissue immune cell analysis reveals tissue-specific features in humans:
tissue_na1 <- read.table('/data/rluo4/RPMfunc/disco_pdata/E-MTAB-11536.sdrf.txt', sep = '\t', header = T)
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
tissue_na2 <- read.table('/data/rluo4/RPMfunc/disco_pdata/E-MTAB-9389.sdrf.txt', sep = '\t', header = T)
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

tissue_na3 <- read.csv('/data/rluo4/RPMfunc/disco_pdata/GSE185381_SRA_Table.txt',  header = T)
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
disease_na1 <- read.table('/data/rluo4/RPMfunc/disco_pdata/E-MTAB-8884.sdrf.txt', sep = '\t', header = T)
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

disease_na2 <- read.table('/data/rluo4/RPMfunc/disco_pdata/E-MTAB-9139.sdrf.txt', sep = '\t', header = T)
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

setwd('/data/rluo4/RPMfunc/disco_pdata')
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

tissue_summary <- read.csv('/data/rluo4/All/Output/tissue_summary.txt',sep = "\t",header = T)
PCTanno_GEO <- unique(tissue_summary$Dataset)
table(disco_GEO.acc %in% PCTanno_GEO)
PCTanno_disco <- intersect(disco_GEO.acc, PCTanno_GEO) # 9 datasets
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
unique(disco_allDat$sample_alias)# 132 sample for raw data 
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
# write.table(asp_sc,file="/data/rluo4/EpiTrans/DataCollection/asp_disco_scRNA.link.txt",sep = "\t",quote = F, col.names = F, row.names = F)
# # saveRDS(disco,  file=paste0('/data/rluo4/RPMfunc/disco_pdata/disco_all_meta.rds'))
# 4) check the lineage for celltypes 
disco <- readRDS(paste0('/data/rluo4/RPMfunc/disco_pdata/disco_all_meta.rds'))
length(unique(disco$ct))
Lineage <- data.frame(Major_type = disco$ct, Minor_type = disco$ct, Cell_type = disco$ct)
Lineage <- Lineage[!duplicated(Lineage$Cell_type),]
Lineage <- Lineage[order(Lineage$Cell_type), ]

table_ct <- read_excel('/data/rluo4/RPMfunc/disco_pdata/table_all_ct.xlsx', sheet = 1)
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
patterns <- c("lymphoid","gdT","ILC"," B", "B cell","lasma", "T ", " T", "CD4", "CD8", "reg", 
              "yeloid", "onocyte", "acrophage", "eutrophi", "ranulocyte", "Tfh", "NK", "DC", "endritic", "mast","MAIT")
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

patterns <- c("ibro","muscl", 'peri', 'vascula', 'mura', 'HSC')
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

patterns <- c("rythroblast",'Megakaryocyte','erythroid','erythrocyte')
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

table(is.na(match(disco$ct, Lineage$Cell_type)))
disco$Major_type <- Lineage$Lineage[match(disco$ct, Lineage$Cell_type)]
saveRDS(disco,  file=paste0('/data/rluo4/RPMfunc/disco_pdata/disco_all_meta.rds'))

# 5) Check the last ascp links 
## last 100+ datasets
disco <- readRDS(paste0('/data/rluo4/RPMfunc/disco_pdata/disco_all_meta.rds'))
sample.index <- c(grep("ormal", disco$sample_type),grep("ealthy", disco$sample_type))
disco_diseased <-  disco[-sample.index,]
table(disco_diseased$sample_type)
disco_meta <- data.frame(table(disco_diseased$sample_type, disco_diseased$disease))
table(is.na(disco_diseased$patient_id))
length(unique(disco_diseased$patient_id)) #1332 diseased samples
disco_sample_summary <- disco_diseased[! duplicated(disco_diseased$patient_id),]

setwd(data_dir)
library(stringr)
disco_GEO <- list.files('./', pattern='.txt')
disco_GEO.acc <- str_split(gsub('.txt','',disco_GEO), '_', simplify = T)[,1]
unique(disco_GEO.acc)#117 --> 294
table(unique(disco_GEO.acc) %in% unique(disco_diseased$project_id))
setdiff(unique(disco_GEO.acc), unique(disco_diseased$project_id))

tissue_summary <- read.csv('/data/rluo4/All/Output/tissue_summary.txt',sep = "\t",header = T)
PCTanno_GEO <- unique(tissue_summary$Dataset)
table(disco_GEO.acc %in% PCTanno_GEO)
PCTanno_disco <- intersect(disco_GEO.acc, PCTanno_GEO) # 9 datasets --> 17 intersection
PCTanno_disco
# View(tissue_summary[tissue_summary$Dataset %in% PCTanno_disco,])
disco_GEO <- disco_GEO[!disco_GEO.acc %in% PCTanno_GEO]#281
disco_GEO.acc <- str_split(gsub('.txt','',disco_GEO), '_', simplify = T)[,1]
table(disco_GEO.acc %in% unique(disco_diseased$project_id))
disco_GEO <- disco_GEO[disco_GEO.acc %in% unique(disco_diseased$project_id)]#165
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

table(disco_allDat$sample_alias %in% disco_diseased$patient_id)
disco_allDat <- disco_allDat[disco_allDat$sample_alias %in% disco_diseased$patient_id,]
unique(disco_allDat$submitted_aspera)
unique(disco_allDat$sample_alias)# 413 sample for raw data 
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
asp_sc_first <- read.table(file="/data/rluo4/EpiTrans/DataCollection/asp_disco_scRNA.link.txt",sep = "\t")
table(asp_sc_first$V3 %in% asp_sc$sra_aspera) 
table(asp_sc$sra_aspera %in% asp_sc_first$V3) 
asp_sc_last <- asp_sc[! asp_sc$sra_aspera %in% asp_sc_first$V3,]
# write.table(asp_sc_last,file="/data2/rluo4/EpiTrans/DataCollection/asp_disco_scRNA_last.link.txt",sep = "\t",quote = F, col.names = F, row.names = F)
uncommon <- read_excel('/data/rluo4/RPMfunc/xuhaixia/disco_tsv/uncommon.xlsx')
snGSE <- read_excel('/data/rluo4/RPMfunc/xuhaixia/disco_tsv/snRNAGSE.xlsx',sheet = 2)
table(snRNAGSE$GSE %in% asp_sc$study_alias)
table(snRNAGSE$GSE %in% disco_GEO.acc)
table(snRNAGSE$GSE %in% disco_diseased$project_id)
setdiff(snRNAGSE$GSE, asp_sc$study_alias)
intersect(snRNAGSE$GSE, asp_sc$study_alias)#"GSE141552" "GSE144136" "GSE157783" "GSE174367"
intersect(snRNAGSE$GSE, asp_sc_last$study_alias)
unique(asp_sc$study_alias)#107
unique(disco_allDat_GSE)#112
# 6)  Write down cell_barcode_annotations.tsv for each dataset
path <- "/data/rluo4/EpiTrans/Disco/Raw/"
setwd('/data/rluo4/EpiTrans/Disco/Raw/Bam/')
bam_files <- list.dirs('/data/rluo4/EpiTrans/Disco/Raw/Bam', full.names = TRUE)
bam_files <- bam_files[grep('GSE', bam_files)]
types <- c('Bam','FastQ','Sra')
file_info_allDat <- NULL
for (type in types) {
  raw_files <- list.dirs(paste0('/data/rluo4/EpiTrans/Disco/Raw/',type), full.names = TRUE)
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
asp_sc_file <- asp_sc[-grep('_1.fastq.gz', asp_sc$filename), ]
asp_sc_file$filename <- gsub('.bai','',asp_sc_file$filename)
index <- asp_sc_file$datatype=='Sra'
asp_sc_file$filename[index] <- paste0(asp_sc_file$filename[index], '_2.fastq.gz')
asp_sc_file$filesize <- file_info_allDat$size[match(asp_sc_file$filename, file_info_allDat$filename)]

table(duplicated(asp_sc_file$sample_alias))
rownames(asp_sc_file) = 1:nrow(asp_sc_file)
tmp=by(asp_sc_file,asp_sc_file$sample_alias,function(x) rownames(x)[which.max(x$filesize)])
probes = as.character(tmp)
table(is.na(probes))
table(rownames(asp_sc_file) %in% probes)
asp_sc_file=asp_sc_file[rownames(asp_sc_file) %in% probes ,]
asp_sc_file$path = rownames(file_info_allDat)[match(asp_sc_file$filename, file_info_allDat$filename)]
table(asp_sc_file$sample_alias %in% disco_allDat$sample_alias)
# write.table(asp_sc_file, file="/data/rluo4/EpiTrans/DataCollection/asp_sc_file.txt",sep = "\t",quote = F, col.names = F, row.names = F)
asp_sc_file <- read.table("/data/rluo4/EpiTrans/DataCollection/asp_sc_file.txt")
extract_first_element <- function(x) {
  split_string <- strsplit(x, "[--]")
  first_element <- sapply(split_string, function(y) head(y, n = 1))
  return(first_element)
}
datasets <- unique(disco$project_id)
for(i in datasets){
  # disco_fisrt <- disco[disco$patient_id %in% asp_sc_file$sample_alias,]
  disco_meta <- disco[disco$project_id %in% i,]
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
  Index = first_elements
  Index = str_split(Index,'[-]',simplify = T)[,1]
  meta <- data.frame('Index' = Index,
                     'Cell_type'= disco_meta[,c('celltype')])
  print(table(meta$Index==''))
  print(table(is.na(meta$Index)))
  write.table(meta, paste0('/data/rluo4/EpiTrans/Disco/MetaData/',i,'.cell_barcode_annotations.tsv'),sep = '\t',
              row.names = F,quote = F)
}

# all_separators <- unique(unlist(strsplit(disco_meta$ct, "[^;|]+")))[-1]
# Filter for characters that could be potential separators
# all_characters <- unique(unlist(strsplit(disco_meta$ct, "")))
# potential_separators <- all_characters[grep("[[:punct:]]", all_characters)]
# potential_separators_and_whitespace <- all_characters[grep("[[:punct:][:space:]]", all_characters)]
# # Print the potential separators
# print(potential_separators)
# # Define the replacement character
# replacement <- "_"
# # Replace potential separators with the replacement character
# modified_vector <- gsub(paste0("[", paste(potential_separators, collapse = ""), "]"), replacement, disco_meta$ct)
# # Replace potential separators and whitespace with the replacement character
# modified_vector <- gsub(paste0("[", paste(potential_separators_and_whitespace, collapse = ""), "]"), replacement, disco_meta$ct)
# # Replace potential separators and whitespace with the replacement character
# modified_vector <- gsub(paste0("[", paste(potential_separators_and_whitespace, collapse = ""), "]+"), replacement, disco_meta$ct)
# # Print the modified vector
# print(modified_vector)








#!/usr/local/bin/Rscript
# title: "RBM33 is a unique m(6)A RNA-binding protein that regulates ALKBH5 demethylase activity and substrate selectivity"
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
out_dir <-  '/data/rluo4/RPMfunc/Output/scRNA/PCTanno_SNV'
data_path = '/data/rluo4/All/Output'
# load('/data/rluo4/All/Output/pheno_all.rds')#Cervix_Epi
# cell_summary <- read.csv('/data/rluo4/All/Output/cell_summary.txt',sep = "\t",header = T)
# tissue_summary <- read.table('/data/rluo4/RPMfunc/PCTanno_pdata/tissue_summary.txt')
# # Helper function to load and filter files
# # Load Diseased Samples
# Diseased_Samples <- read.table('/data/rluo4/RPMfunc/PCTanno_pdata/Diseased_Samples.txt')
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
# indir='/data/rluo4/RPMfunc/PCTanno_pdata/PCTanno_SNV_Coverage'
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
# write.table(Diseased_Samples, sep = "\t", file = '/data/rluo4/RPMfunc/PCTanno_pdata/Diseased_Samples_PCTanno.txt', quote = F)
Diseased_Samples <- Diseased_Samples_PCTanno <- read.csv('/data/rluo4/RPMfunc/PCTanno_pdata/Diseased_Samples_PCTanno.txt',sep = '\t', fill = T)#Diseased_Samples
# ##################################
#   Merge all scRNA SNV Results   #
# ##################################
# Base directory and output directory
base_dir <- '/data/rluo4/RPMfunc/SCOMATIC'
out_dir <- '/data/rluo4/RPMfunc/Output/scRNA/PCTanno_SNV'
# Get unique cohorts
combined_cohort <- unique(Diseased_Samples$Cohort)#[Diseased_Samples$Tissue == 'Prostate'])
# Initialize an empty data frame to store merged data (Done in terminal, skip) # from Summary_PCTanno_UTH36.R
# saveRDS(maf, file = file.path(out_dir, '../organ13_maf.rds'))
# save(maf_results, file = file.path(out_dir, '../organ13_maf_results.RData')) # from Summary_PCTanno_UTH36.R
# Display unique sample barcodes
maf <- readRDS(file.path(out_dir, '../organ13_maf.rds'))
unique(maf$Tumor_Sample_Barcode)# 624 --> 655
# Handle NA values in AAChange.refGene
maf$AAChange.refGene[is.na(maf$AAChange.refGene)] <- maf$Hugo_Symbol[is.na(maf$AAChange.refGene)]
table(maf$ExonicFunc.refGene)#外显子区的SNV or InDel变异类型
# nonsynonymous SNV          stopgain          stoploss
# 4102404               162               665
# synonymous SNV           unknown
# 4810             25792
# Aggregate mutation data
maf_PCTanno <- maf
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
# # saveRDS(Epi, file = file.path('/data/rluo4/RPMfunc/Output/scRNA/SComatic_PCTanno.rds'))
# # system("mv SComatic_PCTanno.rds SComatic_pct.rds")
# ################################################################################
# # 1) load pdata of disco database
# # ################################################################################
# # # load('~/lrh/All/Output/organ13_Epi.RData')# from database-sum.R in ts860
# # # data_path = '/home/lorihan/lrh/All/Output'
# data_path = '/data/rluo4/All/Output'
# setwd(data_path)
# RMP_update <- readxl::read_excel('/data/rluo4/EpiTrans/DataCollection/RMP_update.xlsx')
# # cell_summary <- read.csv('cell_summary.txt',sep = "\t",header = T)
# # sc_majortype <- readRDS(file = '/data/rluo4/RPMfunc/Output/scRNA/sc_majortype_1022.rds')
# # tissue_summary <- read.table('/data/rluo4/RPMfunc/PCTanno_pdata/tissue_summary.txt')
# sc_allmeta <- readRDS(file = '/data/rluo4/RPMfunc/Output/scRNA/sc_allmeta_1022.rds')
# RMzyme_alltissues <- sort(unique(sc_allmeta$tissue));
# RMzyme_alldatasets <- sort(unique(sc_allmeta$project_id));
# load(file = '/data/rluo4/RPMfunc/Output/scRNA/sc_allsample_1022.RData', verbose = T)
# extract_last_element <- function(x) {
#   split_string <- strsplit(x, "_")
#   last_element <- sapply(split_string, function(y) tail(y, n = 1))
#   return(last_element)
# }
################################################################################
# 1.1. Check the first ascp links
################################################################################
## last 100+ datasets on UTH57
UTH57_Bam <- list.files('/data2/rluo4/EpiTrans/Disco/Raw/Bam/') #70
UTH57_FastQ <- list.files('/data2/rluo4/EpiTrans/Disco/Raw/FastQ/') #70
UTH57_Sra <- list.files('/data2/rluo4/EpiTrans/Disco/Raw/Sra/') #70
UTH57_Raw <- unique(c(UTH57_Bam, UTH57_FastQ, UTH57_Sra)) #26
# Raw_done <- intersect(UTH57_Raw, BM_MAST39$project_id)
# Raw_done <- intersect(UTH57_Raw, BM_MAST_meta_p1$project_id)
# for (i in UTH57_Raw ) {
#   GSM <-  list.files(file.path("/data2/rluo4/EpiTrans/Disco/scMapping/",i ))
#   print(GSM)
#   for (j in GSM) {
#     print(paste('to be removed dataset: ', i))
#     print(paste('to be removed GSM: ', j))
#     scMapping_files <-  list.files(file.path("/data/rluo4/EpiTrans/Disco/scMapping/",i,j ))
#     print(scMapping_files)
#   }
# }
################################################################################
# 1.2 Renew the last ascp links (last means the samples downloaded in the UTH57 not UTH36)
################################################################################
# disco <- readRDS(paste0('/data/rluo4/RPMfunc/disco_pdata/disco_adjustmeta_0809.rds'))#0725.rds'))
table(disco$sample_type)
sample.index <- c(grep("ormal", disco$sample_type),grep("ealthy", disco$sample_type),
                  grep("control", disco$sample_type))#sample.index <- c(grep("ormal", disco$sample_type),grep("ealthy", disco$sample_type))
disco_diseased <-  disco[-sample.index,]
table(disco_diseased$sample_type)
disco_meta <- data.frame(table(disco_diseased$sample_type, disco_diseased$disease))
table(is.na(disco_diseased$patient_id))
length(unique(disco_diseased$project_id)) #204 diseased datasets -->197 -->196-->214-->213(1 covid dataset:GSE168453)
length(unique(disco_diseased$patient_id)) #1332 diseased samples -->1359 -->1360-->1525-->1838-->1822
disco_sample_summary <- disco_diseased[! duplicated(disco_diseased$patient_id),]

data_dir <- '/data/rluo4/EpiTrans/DataCollection/disco_tsv'
setwd(data_dir)
library(stringr)
disco_GEO <- list.files('./', pattern='.txt') #/data/rluo4/EpiTrans/DataCollection/disco_tsv
disco_GEO.acc <- str_split(gsub('.txt','',disco_GEO), '_', simplify = T)[,1]
unique(disco_GEO.acc)#117 --> 294 --> 277
table(unique(disco_GEO.acc) %in% unique(disco_diseased$project_id))
setdiff(unique(disco_GEO.acc), unique(disco_diseased$project_id))
tissue_summary <- read.table('/data/rluo4/RPMfunc/PCTanno_pdata/tissue_summary.txt')#read.csv('/data/rluo4/All/Output/tissue_summary.txt',sep = "\t",header = T)
PCTanno_GEO <- unique(tissue_summary$Dataset)
table(disco_GEO.acc %in% PCTanno_GEO)
PCTanno_disco <- intersect(disco_GEO.acc, PCTanno_GEO) # 9 datasets --> 17 intersection
PCTanno_disco
# View(tissue_summary[tissue_summary$Dataset %in% PCTanno_disco,])
disco_GEO <- disco_GEO[!disco_GEO.acc %in% PCTanno_GEO]
disco_GEO.acc <- str_split(gsub('.txt','',disco_GEO), '_', simplify = T)[,1]#281 datasets with raw data
table(disco_GEO.acc %in% unique(disco_diseased$project_id))
disco_GEO <- disco_GEO[disco_GEO.acc %in% unique(disco_diseased$project_id)]#60
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
unique(disco_allDat$sample_alias)# 255 sample for raw data 
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
)# 42 studies
# asp_sc_last <- asp_sc[! asp_sc$sra_aspera %in% asp_sc_first$V3,]
asp_sc_first <- read.table(file="/data/rluo4/EpiTrans/DataCollection/asp_disco_scRNA.link.txt",sep = "\t")
table(asp_sc_first$V3 %in% asp_sc$sra_aspera) 
table(asp_sc$sra_aspera %in% asp_sc_first$V3) 
asp_sc_last <- read.table(file="/data/rluo4/EpiTrans/DataCollection/asp_disco_scRNA_last_renew76.link.txt",sep = "\t")
table(asp_sc$study_alias %in% asp_sc_last$V1)
table(asp_sc$study_alias %in% asp_sc_first$V1)
# intersect(asp_sc$study_alias, asp_sc_last$V1)
# intersect(asp_sc$study_alias, asp_sc_first$V1)
#
disco_BM_MAST <- readRDS(file=paste0('/data/rluo4/RPMfunc/disco_pdata/BM_MAST_meta.rds'))
#
AML <- readRDS(file=paste0('/data/rluo4/RPMfunc/disco_pdata/AML_meta.rds')) # made through /data2/rluo4/Rcode/RMDB/DiscoDatasets.R
#
PDAC1 <- readRDS(file=paste0('/data/rluo4/RPMfunc/disco_pdata/PDAC1_meta.rds'))
# delete the wrong datasets (mast/basophils) from disco
pancreas_pid <- c('GSE205049','GSE217837','GSE231535')
excluded_data <- c(unique(AML$project_id), pancreas_pid, #unique(BM$project_id), 
                   unique(disco_BM_MAST$project_id)) #107 excluded datasets -->115
sort(excluded_data)
excluded_data <- data.frame(excluded_datasets = excluded_data, source = c(rep('AML', length(unique(AML$project_id))),
                                                                          rep('pancreas', length(unique(pancreas_pid))), #rep('BM', length(unique(BM$project_id))), 
                                                                          rep('BM_MAST', length(unique(disco_BM_MAST$project_id)))) # exclude 3 BM_MAST datasets c('GSE211033','GSE196052','GSE196676')
) 
renew_datasets <- intersect(excluded_data$excluded_datasets, asp_sc_first$V1)# all FastQ
# [1] "GSE155468" "GSE152042" "GSE151192" "GSE150825" "GSE130560"
table(renew_datasets %in% disco_BM_MAST$project_id)
# View(asp_sc_first[asp_sc_first$V1 %in% renew_datasets,])
# View(disco[disco$project_id %in% renew_datasets,])
intersect(renew_datasets, asp_sc_last$V1)
# [1] "GSE155468" "GSE151192" "GSE150825" "GSE130560"
setdiff(renew_datasets, asp_sc_last$V1)
# [1] "GSE152042"
# UTH36_wrong_mapping <- renew_datasets
# for (i in UTH36_wrong_mapping ) {
#   GSM <-  list.files(file.path("/data/rluo4/EpiTrans/Disco/scMapping/",i ))
#   print(GSM)
#   for (j in GSM) {
#     print(paste('to be removed dataset: ', i))
#     print(paste('to be removed GSM: ', j))
#     scMapping_files <-  list.files(file.path("/data/rluo4/EpiTrans/Disco/scMapping/",i,j ))
#     print(scMapping_files)
#   }
#   system(paste0('rm -r ', '/data/rluo4/EpiTrans/Disco/scMapping/', i))
#   system(paste0('rm ', '/data/rluo4/EpiTrans/Disco/MetaData/', i, '.cell_barcode_annotations.tsv'))
#   system(paste0('rm -r ', '/data/rluo4/EpiTrans/Disco/SCOMATIC/', i))
#   system(paste0('rm -r ', '/data/rluo4/EpiTrans/Disco/anno.var/', i))
# }
intersect(asp_sc_first$V1, asp_sc_last$V1)
# [1] "GSE130560" "GSE150825" "GSE151192" "GSE155468"
# split 4 renewed datasets into 2 parts in 2 servers:
# 1)  cd /data2/rluo4/EpiTrans/Disco/Raw/FastQ
# rsync -avz -P -e "ssh -o 'ProxyJump rluo4@129.106.31.39'" GSE155468/ rluo4@129.106.31.36:/data/rluo4/EpiTrans/Disco/Raw/FastQ/GSE150825/
# rsync -avz -P -e "ssh -o 'ProxyJump rluo4@129.106.31.39'" GSE151192/ rluo4@129.106.31.36:/data/rluo4/EpiTrans/Disco/Raw/FastQ/GSE151192/
# 2) cd /data2/rluo4/EpiTrans/Disco/Raw/cellranger_count
# rsync UTH36:GSE150825 --> UTH57:GSE150825
# rsync -avz -P -e "ssh -o 'ProxyJump rluo4@129.106.31.39'" rluo4@129.106.31.36:/data/rluo4/EpiTrans/Disco/Raw/SCOMATIC/GSE150825/ /data2/rluo4/EpiTrans/Disco/Raw/SCOMATIC/GSE150825/
# rsync -avz -P -e "ssh -o 'ProxyJump rluo4@129.106.31.39'" rluo4@129.106.31.36:/data/rluo4/EpiTrans/Disco/Raw/scMapping/GSE150825/ /data2/rluo4/EpiTrans/Disco/Raw/scMapping/GSE150825/
# drwxrwxr-x 2 rluo4 rluo4 4.0K Aug 12 15:53 GSM4559281
# drwxrwxr-x 2 rluo4 rluo4 4.0K Aug 13 04:59 GSM4586332
# drwxrwxr-x 2 rluo4 rluo4 4.0K Aug 13 09:51 GSM4586335
# drwxrwxr-x 2 rluo4 rluo4 4.0K Aug 13 19:39 GSM4586339
# rsync UTH36:GSE130560 --> UTH57:GSE130560
# rsync -avz -P -e "ssh -o 'ProxyJump rluo4@129.106.31.39'" rluo4@129.106.31.36:/data/rluo4/EpiTrans/Disco/Raw/cellranger_count/GSE130560/ ./GSE130560/

asp_sc_renew <- asp_sc[ asp_sc$study_alias %in% renew_datasets,]
extract_first_element <- function(x) {
  split_string <- strsplit(x, "[--]")
  first_element <- sapply(split_string, function(y) head(y, n = 1))
  return(first_element)
}
################################################################################
# for(i in renew_datasets){#datasets_PDAC){#datasets_AML # datasets){
#   # disco_fisrt <- disco[disco$patient_id %in% asp_sc_file$sample_alias,]
#   disco_meta <- sc_majortype[sc_majortype$project_id %in% i,]
#   all_characters <- unique(unlist(strsplit(disco_meta$ct, "")))
#   potential_separators_and_whitespace <- all_characters[grep("[[:punct:][:space:]]", all_characters)]
#   # Define the replacement character
#   replacement <- "_"
#   # Replace each potential separator with the replacement character
#   modified_vector <- disco_meta$ct
#   for (sep in potential_separators_and_whitespace) {
#     modified_vector <- gsub(sep, replacement, modified_vector, fixed = TRUE)
#   }
#   # Replace consecutive occurrences of the replacement character with a single underscore
#   disco_meta$celltype <- gsub("_+", replacement, modified_vector)
# 
#   # Call the function to extract the last element from each string
#   first_elements <- extract_first_element(rownames(disco_meta))
#   Index = first_elements
#   Index = str_split(Index,'[-]',simplify = T)[,1]
#   meta <- data.frame('Index' = Index,
#                      'Cell_type'= disco_meta[,c('celltype')])
#   print(table(meta$Index==''))
#   print(table(is.na(meta$Index)))
#   write.table(meta, paste0('/data/rluo4/EpiTrans/Disco/MetaData/',i,'.cell_barcode_annotations.tsv'),sep = '\t',
#               row.names = F,quote = F)
# }
################################################################################
# 1.3 Renew the Disco sample metadata for SComatic analysis
################################################################################
path <- "/data/rluo4/EpiTrans/Disco/Raw/"
setwd('/data/rluo4/EpiTrans/Disco/Raw/Bam/')
bam_files <- list.dirs('/data/rluo4/EpiTrans/Disco/Raw/Bam', full.names = TRUE)
bam_files <- bam_files[grep('GSE', bam_files)]
types <- c('Bam','FastQ','Sra')
file_info_allDat <- NULL
for (type in types) {
  raw_files <- list.dirs(paste0('/data/rluo4/EpiTrans/Disco/Raw/',type), full.names = TRUE)
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
asp_sc_file <- asp_sc[-grep('_1.fastq.gz', asp_sc$filename), ]

asp_sc_file <- asp_sc_file[ ! grepl('bai', asp_sc_file$sra_aspera),]
# asp_sc_file$filename <- gsub('.bai','',asp_sc_file$filename)
index <- asp_sc_file$datatype=='Sra'
asp_sc_file$filename[index] <- paste0(asp_sc_file$filename[index], '_2.fastq.gz')
asp_sc_file$filesize <- file_info_allDat$size[match(asp_sc_file$filename, file_info_allDat$filename)]
max(asp_sc_file$filesize[!is.na(asp_sc_file$filesize)]) #25301137945
# View(asp_sc_file[is.na(asp_sc_file$filesize),])
na.asp_sc_file <- asp_sc_file[is.na(asp_sc_file$filesize),]
unique(asp_sc_file$study_alias)
################################################################################
# 1.4 Organize the files for disco's SComatic: filtering the files that were done
################################################################################
GSE.first <- list.files('/data/rluo4/EpiTrans/Disco/Raw/cellranger_count/')
Bam.first <- NULL
for (gse in GSE.first) {
  setwd('/data/rluo4/EpiTrans/Disco/Raw/cellranger_count/')
  f <- list.files(gse)
  bam <- list.files(paste0('/data/rluo4/EpiTrans/Disco/Raw/cellranger_count/',gse, '/',f, '/outs/'))
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
raw_undone <- setdiff(diff, asp_sc_first$V1)
table(raw_undone %in% asp_sc_last$V1)

table(unique(asp_sc_renew$sample_alias) %in% raw_done)
# table(unique(BM_MAST$patient_id) %in% raw_done)
table(asp_sc_renew$study_alias %in% disco_BM_MAST$project_id)
# asp_sc_renew <- asp_sc_renew[asp_sc_renew$study_alias %in% disco_BM_MAST$project_id, ]
unique(asp_sc_renew$sample_alias) #
table( asp_sc_renew$sample_alias %in% raw_done)
# asp_sc_renew <- asp_sc_renew[ ! asp_sc_renew$sample_alias %in% raw_done, ]
unique(asp_sc_renew$sample_alias) #180-->309-->334
table(unique(asp_sc_renew$sample_alias) %in% unique(disco$patient_id))
unique(asp_sc_renew$study_alias) #22-->42-->43

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

UTH36_scMapping <- list.files('/data/rluo4/EpiTrans/Disco/scMapping/')
GSM_scMapping <- NULL
for (i in UTH36_scMapping) {
  j <- list.files(file.path('/data/rluo4/EpiTrans/Disco/scMapping/', i))
  GSM_scMapping <<- c(GSM_scMapping, j)
}
raw_undone <- setdiff(asp_sc_file$sample_alias, GSM_scMapping)
View(asp_sc_file[asp_sc_file$sample_alias %in% raw_undone,])
unique(asp_sc_file[asp_sc_file$sample_alias %in% raw_undone,]$study_alias)
# [1] "GSE130560" "GSE151192" "GSE152042"
# [4] "GSE155468"

table(asp_sc_file$sample_alias %in% disco_allDat$sample_alias)
unique(asp_sc_file$study_alias)
unique(asp_sc_file$sample_alias)

asp_sc_file_BM_MAST <- asp_sc_file[asp_sc_file$study_alias %in% renew_datasets[-5],] # "GSE130560" done on UTH57
asp_sc_file_BM_MAST <- asp_sc_file[asp_sc_file$study_alias %in% renew_datasets[-(4:5)],] # mv GSE150825 and GSE130560 to UTH57
write.table(asp_sc_file_BM_MAST, 
            file="/data/rluo4/EpiTrans/DataCollection/asp_sc_file_transformed_BM_MAST_new.txt", #asp_sc_file_transformed_BM_MAST.txt",
            sep = "\t", quote = F, col.names = F, row.names = F)
unique(asp_sc_file_BM_MAST$sample_alias)

table(asp_sc_file_BM_MAST$sample_alias %in% asp_sc_file$sample_alias)

################################################################################
# Epi$project_id <- tissue_summary$Dataset[match(Epi$sample_name, tissue_summary$Cohort)]
# # write_rds(Epi, file.path('/data/rluo4/RPMfunc/Output/scRNA/SComatic_PCTanno.rds'))
# table(Epi$barcode %in% sc_majortype$barcode)
# Epi$ct <- sc_majortype$ct[match(Epi$barcode, sc_majortype$barcode)]
# Epi$cell_type <- str_split(Epi$CB, '--', simplify = T)[,3]#Epi$Cell_Type
# index = grepl("cd45---.*?", Epi$CB) #145 CB
# Epi$CB[index];Epi$cell_type[index]
# # Further split the second part on the last occurrence of '_'
# Epi$cell_type[index] <- gsub('-', '', Epi$cell_type[index])
# table(Epi$cell_type==Epi$Cell_Type)
# all_characters <- unique(unlist(strsplit(Epi$ct, "")))
# potential_separators_and_whitespace <- all_characters[grep("[[:punct:][:space:]]", all_characters)]
# # Define the replacement character
# replacement <- "_"
# # Replace each potential separator with the replacement character
# modified_vector <- Epi$ct
# for (sep in potential_separators_and_whitespace) {
#   modified_vector <- gsub(sep, replacement, modified_vector, fixed = TRUE)
# }
# # Replace consecutive occurrences of the replacement character with a single underscore
# Epi$Cell_Type <- gsub("_+", replacement, modified_vector)
# Epi$disease <- sc_majortype$disease[match(Epi$barcode, sc_majortype$barcode)]
# Epi$sub_type <- sc_majortype$sub_type[match(Epi$barcode, sc_majortype$barcode)]
# Epi$Organ <- sc_majortype$tissue[match(Epi$barcode, sc_majortype$barcode)]
# saveRDS(Epi, file = file.path('/data/rluo4/RPMfunc/Output/scRNA/SComatic_PCTanno.rds'))
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

# Renew on 12/07/2024:
library(maftools)
library(data.table)
library(stringr)
library(dplyr)
options(bitmapType = 'cairo')
load(file = '/data/rluo4/RPMfunc/Output/scRNA/sc_allsample_1122.RData', verbose = T)
##################################################
# change 1:
##################################################
patient_id <- sc_allmeta$patient_id
barcode_idx <- sc_allmeta$barcode
idx9 <- grepl('GSE181254|PRJCA001063|GSE148673|GSE163203', sc_allmeta$project_id) 
sc_allmeta$barcode[idx9][1:5]
barcode_idx[idx9] <- extract_last_element(sc_allmeta$rownames[idx9])

idx_dt <-  idx9
barcode_idx[idx_dt] <- paste0(barcode_idx[idx_dt], '--', patient_id[idx_dt])
dt <- unique(sc_allmeta[idx_dt, 'project_id'])
# setdiff( dt, unique(sc_allmeta$project_id[idx_diff]))
sc_allmeta$barcode[idx_dt] <- barcode_idx[idx_dt]
saveRDS(sc_allmeta, file = '/data/rluo4/RPMfunc/Output/scRNA/sc_allmeta_1122.rds') # combined with diff meta from merge_nmf_meta and allmeta
sc_allsample <- sc_allmeta %>% dplyr::select(patient_id, project_id, tissue, disease) %>%  distinct()
save(disco, tissue_summary, sc_summary, sc_allsample, file = '/data/rluo4/RPMfunc/Output/scRNA/sc_allsample_1122.RData')

################################################################################
# Renew on 12/12/2024:
################################################################################
idx11 <- grepl('GSE137829', sc_allmeta$project_id) 
sc_allmeta$barcode[idx11][1:5]
sc_allmeta$disease[idx11] <- 'neuroendocrine prostate cancer' # NEPC
saveRDS(sc_allmeta, file = '/data/rluo4/RPMfunc/Output/scRNA/sc_allmeta_1122.rds') # combined with diff meta from merge_nmf_meta and allmeta
sc_allsample <- sc_allmeta %>% dplyr::select(patient_id, project_id, tissue, disease) %>%  distinct()
save(disco, tissue_summary, sc_summary, sc_allsample, file = '/data/rluo4/RPMfunc/Output/scRNA/sc_allsample_1122.RData')
