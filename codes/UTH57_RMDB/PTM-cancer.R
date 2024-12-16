#!/usr/local/bin/Rscript
# title: "RMzyme: Regulations of RNA Modifying Enzymes in Human"
# author: "Ruihan Luo"
# date: "April 19th,2024"
# rm(list=ls())
library(data.table)
library(stringr)
library(tidyverse) 
library(GEOquery)
library(readxl)
Sys.setenv("VROOM_CONNECTION_SIZE" = 999999999)
getOption('timeout')
options(timeout=100000000)
path <- c('/home/rluo4/R/x86_64-pc-linux-gnu-library/4.3', '/opt/R/4.3.1/lib64/R/library', '/home/rluo4/R/x86_64-conda-linux-gnu-library/4.3', '/data2/rluo4/bin/miniconda3/lib/R/library')
.libPaths(path)
raw_directory <- '/data2/rluo4/EpiTrans/RMDatasets/Raw'
data_directory <- '/data2/rluo4/EpiTrans/RMDatasets/GEO'
setwd(data_directory)
plot_dir <- '/data2/rluo4/RPMfunc/Output/GEO'


##################################
# 1) load the Rdata from RMdeg.R #
##################################
load(paste0('/data2/rluo4/EpiTrans/RMDatasets/allDat.RData')) # pdata from RMDatasets.R
load(paste0('/data2/rluo4/EpiTrans/RMDatasets/RMP_allDat.RData')) # asp metadata from RMDatasets_UTH36.R
load(paste0('/data2/rluo4/EpiTrans/RMDatasets/RMdeg_allDat.RData')) # from RMdeg.R
RMP_update <- read_excel('/data2/rluo4/EpiTrans/DataCollection/RMP_update.xlsx') # from RM.R
asp_sc <- fread(file="/data2/rluo4/EpiTrans/DataCollection/asp_RMP_epitrans.link.txt", header = F)

##################################################
# 1) PTM of PDC data at https://pdc.cancer.gov/  #
##################################################
library(org.Hs.eg.db)
library(clusterProfiler)
# clinical data Cancer Cell(Article) ----------------------------------------------------
cancercellcli <- read_xlsx("/data2/rluo4/EpiTrans/PTM/CPTAC/inputdata/cancer_cell_clinical.xlsx",sheet = 2)
PTM_var <- fread('/data2/rluo4/EpiTrans/PTM/CELL/ref/var_map_full_v4.tsv')
table(PTM_var$feature)
PTM_var <- PTM_var[! PTM_var$feature %in% c('proteome', 'transcriptome'),]
uniprot <- read.csv('/data2/rluo4/All/uniprot-hs.tsv',fill = T,header = T,sep = '\t')
uniprot_loc <- read.table('/data2/rluo4/RPMfunc/Output/summary/uniprotLoc.txt',fill = T,header = T,sep = '\t')
# 1.1.batch acety data processing -------------------------------------------------
acetfil <- list.files("/data2/rluo4/EpiTrans/PTM/CPTAC/inputdata/Acetylome_UMich_GENCODE34_v1/Acetylome_UMich_GENCODE34_v1/",pattern = "_single-site_MD",recursive = T)
datlist <- lapply(acetfil, function(x){
  read.csv(file = file.path('/data2/rluo4/EpiTrans/PTM/CPTAC/inputdata/Acetylome_UMich_GENCODE34_v1/Acetylome_UMich_GENCODE34_v1/',x),
           sep = "\t",check.names = F) 
})
for (tumor in 1:length(acetfil)) {
  names(datlist)[[tumor]] <- strsplit(acetfil[tumor],'/')[[1]][1]
}

# acetylist <- list()
acetylome <- data.frame()
for (i in names(datlist)) {
  # i = names(datlist)[2]
  message(paste0("Running ",i,"..."))
  acety <- datlist[[i]][,1:6] %>% 
    mutate(ENSEMBL = str_sub(Gene,1,15)) %>% 
    mutate(location = map_chr(Index, ~ strsplit(.x, "_")[[1]][2]))
  set.seed(123)
  acety_imid <- bitr(str_sub(acety$Gene,1,15),fromType = "ENSEMBL",
                     toType = c("SYMBOL","UNIPROT"),OrgDb ="org.Hs.eg.db") 
  acety_im_symbol <- merge(acety_imid,acety,by.x = "ENSEMBL",by.y = "ENSEMBL",all = F) %>% 
    mutate(GENE_locat = paste(SYMBOL,location)) %>% 
    mutate(Disease = i) %>%
    select(ENSEMBL,SYMBOL, UNIPROT,location,GENE_locat,Disease,ProteinID,Peptide) #%>% t() 
  
  acetylome <<- rbind(acetylome, acety_im_symbol)
}
acetylome$feature <- 'acetylome' # including tumor + normal
dim(acetylome)
rm_duplicates <- paste(acetylome$ENSEMBL, acetylome$UNIPROT, acetylome$GENE_locat, acetylome$Disease, acetylome$Peptide)
table(duplicated(rm_duplicates))
# acetylome <- acetylome[! duplicated(rm_duplicates), ]
# save(acetylome,file = "/data2/rluo4/EpiTrans/PTM/CELL/acetylome.rdata")
load("/data2/rluo4/EpiTrans/PTM/CPTAC/inputdata/acetylist.rdata")
otherCancers <- list.files("/data2/rluo4/EpiTrans/PTM/CELL/PDCdata/",pattern = "Acetylome.acetylsite")
# otherCancers <- otherCancers[-grep('Breast', otherCancers)]
otherT <- lapply(otherCancers, function(x) {
  file_path <- file.path('/data2/rluo4/EpiTrans/PTM/CELL/PDCdata/', x)
  data <- read.table(file = file_path, header = TRUE, sep = "\t", check.names = FALSE)
  return(data)
})
names(otherT) <- c('BT')
BrainCancer <- str_split(colnames(otherT$BT)[grepl('CPT', colnames(otherT$BT))], ' ', simplify = T)[,1]
table(BrainCancer %in% rownames(acetylist$GBM))
otheracetylome <- data.frame()
for (i in names(otherT)) {
  message(paste0("Running ",i,"..."))
  col.index <- c(1, (ncol(otherT[[i]])-2):ncol(otherT[[i]]))
  # print(colnames(otherT[[i]])[col.index])
  print(head(otherT[[i]][1:5, col.index]))
  otherphos <- otherT[[i]] %>%  as.data.frame() %>% 
    # round(.,1) %>% 
    mutate(SYMBOL = Gene )%>% #
    mutate(ProteinID = map_chr(Acetylsite, ~ strsplit(.x, ":")[[1]][1])) %>% 
    mutate(location = map_chr(Acetylsite, ~ strsplit(.x, ":")[[1]][2])) #map_chr(Index, ~ strsplit(.x, "_")[[1]][2]))
  set.seed(123)
  # otherphos_imid <- bitr(otherphos$Gene,fromType = "SYMBOL",
  #                       toType = c("ENSEMBL","UNIPROT","ENSEMBLPROT"),OrgDb ="org.Hs.eg.db") 
  # colnames(otherphos_imid) <- gsub('ENSEMBLPROT','ProteinID',colnames(otherphos_imid))
  otherphos_imid <- bitr(otherphos$Gene,fromType = "SYMBOL",
                         toType = c("ENSEMBL","UNIPROT"),OrgDb ="org.Hs.eg.db") 
  otherphos_im_symbol <- merge(otherphos_imid,otherphos,by.x = "SYMBOL",by.y = "SYMBOL",all = F) %>% 
    mutate(GENE_locat = paste(SYMBOL,location)) %>% 
    mutate(Disease = i) %>%
    select(ENSEMBL,SYMBOL, UNIPROT,location,GENE_locat,Disease,ProteinID,Peptide) #%>% t() 
  otheracetylome <<- rbind(otheracetylome, otherphos_im_symbol)
}
otheracetylome$feature <- 'acetylome' # including tumor + normal
dim(otherphos)
rm_duplicates <- paste(otheracetylome$ENSEMBL, otheracetylome$UNIPROT, otheracetylome$GENE_locat, otheracetylome$Disease, otheracetylome$Peptide)
table(duplicated(rm_duplicates))
# View(otheracetylome[duplicated(rm_duplicates),])
otheracetylome <- otheracetylome[! duplicated(rm_duplicates), ]

# 1.2.batch phos data processing -------------------------------------------------
# phos Tumor
profilT <- list.files("/data2/rluo4/EpiTrans/PTM/CPTAC/inputdata/Phosphoproteome_BCM_GENCODE_v34_harmonized_v1/Phosphoproteome_BCM_GENCODE_v34_harmonized_v1/",pattern = "_Tumor")
datlistT <- lapply(profilT, function(x){
  read.table(file = file.path('/data2/rluo4/EpiTrans/PTM/CPTAC/inputdata/Phosphoproteome_BCM_GENCODE_v34_harmonized_v1/Phosphoproteome_BCM_GENCODE_v34_harmonized_v1/',x),
             header = T,sep = "\t",check.names = F) %>% 
    remove_rownames() %>% 
    column_to_rownames("idx") %>% as.matrix()
})

for (tumor in 1:length(profilT)) {
  names(datlistT)[[tumor]] <- strsplit(profilT[tumor],'_')[[1]][1]
}
# phossymbolistT <- list()
phosphoproteome <- data.frame()
for (i in names(datlistT)) {
  message(paste0("Running ",i,"..."))
  phos <- datlistT[[i]]
  # phos <- phos[rowSums(is.nan(phos) | is.infinite(phos)) == 0,]#error外接函数调用时不能有NA/NaN/Inf(arg1)
  # phos_im <- impute.knn(phos[rowSums(is.na(phos))/ncol(phos) <= 0.5,])
  phos <- phos %>% as.data.frame() %>% 
    # round(.,1) %>% 
    mutate(Gene = map_chr(rownames(.), ~ strsplit(.x, "\\|")[[1]][1])) %>% 
    mutate(ProteinID = map_chr(rownames(.), ~ strsplit(.x, "\\|")[[1]][2])) %>%
    mutate(Peptide = map_chr(rownames(.), ~ strsplit(.x, "\\|")[[1]][4]))
  phos_im <- phos %>% 
    mutate(ENSEMBL = str_sub(rownames(.),1,15)) %>% 
    mutate(location = map_chr(rownames(.), ~ strsplit(.x, "\\|")[[1]][3]))
  set.seed(123)
  phos_imid <- bitr(str_sub(rownames(phos_im),1,15),fromType = "ENSEMBL",
                    toType = c("SYMBOL","UNIPROT"),OrgDb ="org.Hs.eg.db") 
  phos_im_symbol <- merge(phos_im,phos_imid,by.x = "ENSEMBL",by.y = "ENSEMBL",all = F) %>% 
    mutate(GENE_locat = paste(SYMBOL,location)) %>% 
    mutate(Disease = i) %>%
    select(ENSEMBL,SYMBOL, UNIPROT,location,GENE_locat,Disease,ProteinID,Peptide) %>% as.data.frame() 
  phosphoproteome <<- rbind(phosphoproteome, phos_im_symbol)
}
# phos Normal
profilN <- list.files("/data2/rluo4/EpiTrans/PTM/CPTAC/inputdata/Phosphoproteome_BCM_GENCODE_v34_harmonized_v1/Phosphoproteome_BCM_GENCODE_v34_harmonized_v1/",pattern = "_Normal")
datlistN <- lapply(profilN, function(x){
  read.table(file = file.path('/data2/rluo4/EpiTrans/PTM/CPTAC/inputdata/Phosphoproteome_BCM_GENCODE_v34_harmonized_v1/Phosphoproteome_BCM_GENCODE_v34_harmonized_v1/',x),
             header = T,sep = "\t",check.names = F) %>% 
    remove_rownames() %>% 
    column_to_rownames("idx") %>% as.matrix()
})

for (tumor in 1:length(profilN)) {
  names(datlistN)[[tumor]] <- strsplit(profilN[tumor],'_')[[1]][1]
}
# phossymbolistN <- list()
for (i in names(datlistN)) {
  message(paste0("Running ",i,"..."))
  phos <- datlistN[[i]]
  phos <- phos %>% as.data.frame() %>% 
    mutate(Gene = map_chr(rownames(.), ~ strsplit(.x, "\\|")[[1]][1])) %>% 
    mutate(ProteinID = map_chr(rownames(.), ~ strsplit(.x, "\\|")[[1]][2])) %>%
    mutate(Peptide = map_chr(rownames(.), ~ strsplit(.x, "\\|")[[1]][4]))
  phos_im <- phos %>% 
    mutate(ENSEMBL = str_sub(rownames(.),1,15)) %>% 
    mutate(location = map_chr(rownames(.), ~ strsplit(.x, "\\|")[[1]][3]))
  set.seed(123)
  phos_imid <- bitr(str_sub(rownames(phos_im),1,15),fromType = "ENSEMBL",
                    toType = c("SYMBOL","UNIPROT"),OrgDb ="org.Hs.eg.db") 
  phos_im_symbol <- merge(phos_im,phos_imid,by.x = "ENSEMBL",by.y = "ENSEMBL",all = F) %>% 
    mutate(GENE_locat = paste(SYMBOL,location)) %>% 
    mutate(Disease = i) %>%
    select(ENSEMBL,SYMBOL, UNIPROT,location,GENE_locat,Disease,ProteinID,Peptide) %>% as.data.frame() 
  
  phosphoproteome <<- rbind(phosphoproteome, phos_im_symbol)
}
dim(phosphoproteome)
phosphoproteome$feature <- 'phosphoproteome' # including tumor + normal
rm_duplicates <- paste(phosphoproteome$ENSEMBL, phosphoproteome$UNIPROT, phosphoproteome$GENE_locat, phosphoproteome$Disease, phosphoproteome$Peptide)
table(duplicated(rm_duplicates))
# View(phosphoproteome[duplicated(rm_duplicates),])
phosphoproteome <- phosphoproteome[! duplicated(rm_duplicates), ]
# save(phosphoproteome,file = "/data2/rluo4/EpiTrans/PTM/CELL/phosphoproteome.rdata")
# load("/data2/rluo4/EpiTrans/PTM/CPTAC/inputdata/phossymbolistT.rdata")
# load("/data2/rluo4/EpiTrans/PTM/CPTAC/inputdata/phossymbolistN.rdata")
# View(phosphoproteome)
otherCancers <- list.files("/data2/rluo4/EpiTrans/PTM/CELL/PDCdata/",pattern = "Phosphoproteome.phosphosite")
otherCancers <- otherCancers[-grep('Breast', otherCancers)]
otherT <- lapply(otherCancers, function(x) {
  file_path <- file.path('/data2/rluo4/EpiTrans/PTM/CELL/PDCdata/', x)
  data <- read.table(file = file_path, header = TRUE, sep = "\t", check.names = FALSE)
  return(data)
})
# names(otherT) <- c('Glioblastoma;Gliomas;Meningiomas','AML','GC', 'Pediatric/AYA Brain Tumors','ICC','HCC') , of 218 tumors across 7 histological types of childhood brain cancer: low-grade glioma (n = 93), ependymoma (32), high-grade glioma (25), medulloblastoma (22), ganglioglioma (18), craniopharyngioma (16), and atypical teratoid rhabdoid tumor (12). 
names(otherT) <- c('BT','AML','GC', 'PBT','ICC','HCC') #MB stands for Medulloblastoma,
# BrainCancer <- str_split(colnames(otherT[[4]])[grepl('CPT', colnames(otherT[[4]]))], ' ', simplify = T)[,1]
# table(BrainCancer %in% rownames(acetylist$GBM))
# table(BrainCancer %in% rownames(phossymbolistT$GBM))
otherphosphoproteome <- data.frame()
for (i in names(otherT)) {
  message(paste0("Running ",i,"..."))
  col.index <- c(1, (ncol(otherT[[i]])-2):ncol(otherT[[i]]))
  # print(colnames(otherT[[i]])[col.index])
  print(head(otherT[[i]][1:5, col.index]))
  otherphos <- otherT[[i]] %>%  as.data.frame() %>% 
    # round(.,1) %>% 
    mutate(SYMBOL = Gene )%>% #
    mutate(ProteinID = map_chr(Phosphosite, ~ strsplit(.x, ":")[[1]][1])) %>% 
    mutate(location = map_chr(Phosphosite, ~ strsplit(.x, ":")[[1]][2])) #map_chr(Index, ~ strsplit(.x, "_")[[1]][2]))
  set.seed(123)
  # otherphos_imid <- bitr(otherphos$Gene,fromType = "SYMBOL",
  #                       toType = c("ENSEMBL","UNIPROT","ENSEMBLPROT"),OrgDb ="org.Hs.eg.db") 
  # colnames(otherphos_imid) <- gsub('ENSEMBLPROT','ProteinID',colnames(otherphos_imid))
  otherphos_imid <- bitr(otherphos$Gene,fromType = "SYMBOL",
                         toType = c("ENSEMBL","UNIPROT"),OrgDb ="org.Hs.eg.db") 
  otherphos_im_symbol <- merge(otherphos_imid,otherphos,by.x = "SYMBOL",by.y = "SYMBOL",all = F) %>% 
    mutate(GENE_locat = paste(SYMBOL,location)) %>% 
    mutate(Disease = i) %>%
    select(ENSEMBL,SYMBOL, UNIPROT,location,GENE_locat,Disease,ProteinID,Peptide) #%>% t() 
  otherphosphoproteome <<- rbind(otherphosphoproteome, otherphos_im_symbol)
}
otherphosphoproteome$feature <- 'phosphoproteome' # including tumor + normal
dim(otherphos)
rm_duplicates <- paste(otherphosphoproteome$ENSEMBL, otherphosphoproteome$UNIPROT, otherphosphoproteome$GENE_locat, otherphosphoproteome$Disease, otherphosphoproteome$Peptide)
table(duplicated(rm_duplicates))
# View(otherphosphoproteome[duplicated(rm_duplicates),])
otherphosphoproteome <- otherphosphoproteome[! duplicated(rm_duplicates), ]
# save(otherphosphoproteome,file = "/data2/rluo4/EpiTrans/PTM/CELL/otherphosphoproteome.rdata")

# 1.3.batch glyco data processing -------------------------------------------------
# glyco Tumor + Normal
glyfilT <- list.files("/data2/rluo4/EpiTrans/PTM/CPTAC/inputdata/Glycoproteome_JHU_v1/N-linked_Glycoproteome-GENCODE_v34/",pattern = "abundances-MD_norm")
glyfilT <- glyfilT[-6]
datlistT <- lapply(glyfilT, function(x) {
  file_path <- file.path('/data2/rluo4/EpiTrans/PTM/CPTAC/inputdata/Glycoproteome_JHU_v1/N-linked_Glycoproteome-GENCODE_v34/', x)
  data <- read.table(file = file_path, header = TRUE, sep = "\t", check.names = FALSE)
  data <- data %>%   # Replace 'None' with NA 
    mutate(across(everything(), ~ gsub('None', NA, .))) %>% mutate(across(-(1:3), as.numeric)) %>%  # Convert columns excluding the first three to numeric
    remove_rownames() %>%  as.matrix()
  return(data)
})
for (tumor in 1:length(glyfilT)) {
  names(datlistT)[[tumor]] <- strsplit(glyfilT[tumor],'_')[[1]][1]
}

glycoproteome <- data.frame()
for (i in names(datlistT)) {
  message(paste0("Running ",i,"..."))
  glyco <- datlistT[[i]] %>%  as.data.frame() %>% 
    mutate(SYMBOL = map_chr(Gene, ~ strsplit(.x, "\\|")[[1]][1])) %>%
    mutate(Gene = map_chr(Gene, ~ strsplit(.x, "\\|")[[1]][2])) %>% 
    mutate(ENSEMBL = substr(Gene, 1, 15)) %>%
    mutate(Peptide = map_chr(Sequence, ~ strsplit(.x, "-")[[1]][1])) %>%
    mutate(location = map_chr(Sequence, ~ strsplit(.x, "-")[[1]][2])) #map_chr(Index, ~ strsplit(.x, "_")[[1]][2]))
  set.seed(123)
  glyco_imid <- bitr(glyco$ENSEMBL,fromType = "ENSEMBL",
                     toType = c("UNIPROT","ENSEMBLPROT"),OrgDb ="org.Hs.eg.db") 
  colnames(glyco_imid) <- gsub('ENSEMBLPROT','ProteinID',colnames(glyco_imid))
  glyco_im_symbol <- merge(glyco_imid,glyco,by.x = "ENSEMBL",by.y = "ENSEMBL",all = F) %>% 
    mutate(GENE_locat = paste(SYMBOL,location)) %>% 
    mutate(Disease = i) %>%
    select(ENSEMBL,SYMBOL, UNIPROT,location,GENE_locat,Disease,ProteinID,Peptide) #%>% t() 
  glycoproteome <<- rbind(glycoproteome, glyco_im_symbol)
}
glycoproteome$feature <- 'glycoproteome' # including tumor + normal
dim(glycoproteome)
rm_duplicates <- paste(glycoproteome$ENSEMBL, glycoproteome$UNIPROT, glycoproteome$GENE_locat, glycoproteome$Disease, glycoproteome$Peptide)
table(duplicated(rm_duplicates))
glycoproteome <- glycoproteome[! duplicated(rm_duplicates), ]
# save(glycoproteome,file = "/data2/rluo4/EpiTrans/PTM/CELL/glycoproteome.rdata")

# 1.4.batch ubiquity data processing -------------------------------------------------
# ubiquity Tumor + Normal
ubifilT <- list.files("/data2/rluo4/EpiTrans/PTM/CPTAC/inputdata/Ubiquitylome_CDAP_v1/",pattern = "CPTAC3_Lung_Squamous_Cell_Carcinoma_Ubiquitylome.ubiquitylsite.tmt11.tsv")
# ubifilT <- ubifilT[-6]
datlistT <- lapply(ubifilT, function(x) {
  file_path <- file.path('/data2/rluo4/EpiTrans/PTM/CPTAC/inputdata/Ubiquitylome_CDAP_v1/', x)
  read.table(file = file_path, header = TRUE, sep = "\t", check.names = FALSE)
})
names(datlistT) = 'LSCC'#[[tumor]] <- strsplit(ubifilT[tumor],'_')[[1]][1]

ubiquitylome <- data.frame()
for (i in names(datlistT)) {
  message(paste0("Running ",i,"..."))
  ubiquity <- datlistT[[i]] %>%  as.data.frame() %>% 
    # round(.,1) %>% 
    mutate(SYMBOL = Gene )%>% #str_sub(Gene,1,15)) %>% 
    mutate(location = map_chr(Ubiquitylsite, ~ strsplit(.x, ":")[[1]][2])) #map_chr(Index, ~ strsplit(.x, "_")[[1]][2]))
  set.seed(123)
  ubiquity_imid <- bitr(ubiquity$Gene,fromType = "SYMBOL",
                        toType = c("ENSEMBL","UNIPROT","ENSEMBLPROT"),OrgDb ="org.Hs.eg.db") 
  colnames(ubiquity_imid) <- gsub('ENSEMBLPROT','ProteinID',colnames(ubiquity_imid))
  ubiquity_im_symbol <- merge(ubiquity_imid,ubiquity,by.x = "SYMBOL",by.y = "SYMBOL",all = F) %>% 
    mutate(GENE_locat = paste(SYMBOL,location)) %>% 
    mutate(Disease = i) %>%
    select(ENSEMBL,SYMBOL, UNIPROT,location,GENE_locat,Disease,ProteinID,Peptide) #%>% t() 
  
  ubiquitylome <<- rbind(ubiquitylome, ubiquity_im_symbol)
}
ubiquitylome$feature <- 'ubiquitylome' # including tumor + normal
dim(ubiquitylome)
rm_duplicates <- paste(ubiquitylome$ENSEMBL, ubiquitylome$UNIPROT, ubiquitylome$GENE_locat, ubiquitylome$Disease, ubiquitylome$Peptide)
table(duplicated(rm_duplicates))
# save(ubiquitylome,file = "/data2/rluo4/EpiTrans/PTM/CELL/ubiquitylome.rdata")
# load("/data2/rluo4/EpiTrans/PTM/CPTAC/inputdata/ubiquitylist.rdata")

# 1.5.Lacty data processing -------------------------------------------------
# https://www.biosino.org/node/browse?keyword=OEZ00008344
# datlistGC <- read_xlsx('/data2/rluo4/EpiTrans/PTM/CELL/PDCdata/PTM-Kla-iScience-GC.xlsx')
datlistHCC <- read_xlsx('/data2/rluo4/EpiTrans/PTM/CELL/PDCdata/PTM-Kla-NM-HCC.xlsx')
# Surgically resected samples were obtained from 52 patients with HBV-related HCC from Zhongshan Hospital, Fudan University. Totally, 58 tumor and 52 adjacent liver tissue samples were collected. All sample identities were renamed with codes (such as 78T, 78P and so on) instead of the patient’s private information. Written informed consent was provided by all patients. All patients had an HBV but no HCV infection background and had not undergone prior chemotherapy or radiotherapy. 
# datlistmHCC <- read_xlsx('/data2/rluo4/EpiTrans/PTM/CELL/PDCdata/PTM-Kla-Proteomics-HCC.xlsx')
datlistLung <- read_xlsx('/data2/rluo4/EpiTrans/PTM/CELL/PDCdata/PTM-Kla-Proteomics-NormalLung.xlsx')
datlistLungReport <- read_xlsx('/data2/rluo4/EpiTrans/PTM/CELL/PDCdata/PTM-lung-Kla-literature.xlsx')
table(datlistLung$`Protein accession` %in% datlistLungReport$Access)
table(datlistLung$`Protein accession` %in% uniprot_loc$Entry)
table(datlistLung$`Protein accession` %in% gene_info$UniProtAcc)
setdiff(datlistLung$`Protein accession`, gene_info$UniProtAcc)
table(datlistLung$`Protein accession` %in% datlistHCC$`Protein accession`)

lacfillT <- list.files("/data2/rluo4/EpiTrans/PTM/CELL/PDCdata/",pattern = "PTM-Kla-")
datlistT <- lapply(lacfillT, function(x) {
  # Read data using fread for faster reading
  data <- read_xlsx(file.path('/data2/rluo4/EpiTrans/PTM/CELL/PDCdata/', x))
  # Extract required columns
  df_new <- data[, grepl('Protein accession|Position|KlaSite|Gene name|Modified sequence|Unique peptide', colnames(data) )]
  # Set column names same as file names
  colnames(df_new) <- c('Protein.accession', 'Position', 'Gene.name', 'Modified.sequence')
  print(table(is.na(df_new$Gene.name)))
  print(table(df_new$Gene.name=='--'))
  df_new$Gene.name[df_new$Gene.name=='--'] <- NA
  return(df_new)
})
str(datlistT)
names(datlistT) <- c('GC', 'HCC', 'mHCC', 'LUNG')
set.seed(123)
lacty_imid <- bitr(datlistT$LUNG$Protein.accession, fromType = "UNIPROT",#"SYMBOL",
                   toType = c("SYMBOL"),OrgDb ="org.Hs.eg.db") 
datlistT$LUNG$Gene.name <- datlistT$HCC$Gene.name[match(datlistT$LUNG$Protein.accession, lacty_imid$UNIPROT)]
table(is.na(datlistT$LUNG$Gene.name))
index <- is.na(datlistT$LUNG$Gene.name)
unique(datlistT$LUNG$Protein.accession[index])
table(datlistT$LUNG$Protein.accession[index] %in% datlistT$HCC$Protein.accession)
table(datlistT$LUNG$Protein.accession[index] %in% datlistT$mHCC$Protein.accession)
table(datlistT$LUNG$Protein.accession[index] %in% gene_info$UniProtAcc)
datlistT$LUNG$Gene.name[index] <- datlistT$HCC$Gene.name[match(datlistT$LUNG$Protein.accession[index], datlistT$HCC$Protein.accession)]
index <- is.na(datlistT$LUNG$Gene.name)
datlistT$LUNG$Gene.name[index] <- datlistT$mHCC$Gene.name[match(datlistT$LUNG$Protein.accession[index], datlistT$mHCC$Protein.accession)]
index <- is.na(datlistT$LUNG$Gene.name)
datlistT$LUNG$Gene.name[index] <- gene_info$geneSymbol[match(datlistT$LUNG$Protein.accession[index], gene_info$UniProtAcc)]
table(is.na(datlistT$LUNG$Gene.name))

lactylome <- data.frame()
for (i in names(datlistT)) {
  message(paste0("Running ",i,"..."))
  lacty <- datlistT[[i]] %>%  as.data.frame() %>% 
    mutate(UNIPROT = Protein.accession )%>% 
    mutate(Peptide = Modified.sequence )%>% 
    # mutate(SYMBOL = Gene.name )%>% 
    mutate(location = paste0("K", Position)) 
  set.seed(123)
  lacty_imid <- bitr(lacty$UNIPROT, fromType = "UNIPROT",#"SYMBOL",
                     toType = c("ENSEMBL", 'SYMBOL',"ENSEMBLPROT"),OrgDb ="org.Hs.eg.db") 
  colnames(lacty_imid) <- gsub('ENSEMBLPROT','ProteinID',colnames(lacty_imid))
  lacty_im_symbol <- merge(lacty_imid,lacty,by.x = "UNIPROT",by.y = "UNIPROT",all = F) %>% 
    mutate(GENE_locat = paste(SYMBOL,location)) %>% 
    mutate(Disease = i) %>%
    dplyr::select(ENSEMBL,SYMBOL, UNIPROT,location,GENE_locat,Disease,ProteinID,Peptide) #%>% t()
  lactylome <<- rbind(lactylome, lacty_im_symbol)
}
lactylome$feature <- 'lactylome' # including tumor + normal
dim(lactylome)
table(is.na(lactylome$SYMBOL))
rm_duplicates <- paste(lactylome$ENSEMBL, lactylome$UNIPROT, lactylome$GENE_locat, lactylome$Disease, lactylome$Peptide)
table(duplicated(rm_duplicates))
lactylome <- lactylome[! duplicated(rm_duplicates), ]
# save(lactylome,file = "/data2/rluo4/EpiTrans/PTM/CELL/lactylome.rdata")

# 1.6.PTM data summarizing -------------------------------------------------
save(acetylome, otheracetylome, phosphoproteome, otherphosphoproteome, 
     glycoproteome, lactylome, ubiquitylome, file = '/data2/rluo4/EpiTrans/PTM/CELL/PDCdata/PTM_PDCdata.rdata')



