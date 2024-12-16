#!/usr/local/bin/Rscript
# title: "RMBase v3.0: decode the landscape, mechanisms and functions of RNA modifications"
# title:
# author: "Ruihan Luo"
# date: "March 15th,2024"
# rm(list=ls())
datasetlists <- read.table('/data2/rluo4/EpiTrans/DataCollection/HsRMdatasets.txt', fill = T, header = T, sep = '\t')
RMDatasets <- read_excel('/data2/rluo4/EpiTrans/DataCollection/Regulatory-enzymes-RNA-modifications.xlsx', sheet = 1)
unique(RMDatasets$GSE)
library(data.table)
library(stringr)
library(tidyverse) 
library(GEOquery)
Sys.setenv("VROOM_CONNECTION_SIZE" = 999999999)
getOption('timeout')
options(timeout=100000000)
data_directory <- '/data2/rluo4/EpiTrans/RMDatasets/GEO/'
setwd(data_directory)
GSE.datasets <- grep("GSE", datasetlists$GSE, value = TRUE)

################################################################################
# 1) below is done in R of terminal
################################################################################
# allDat <- lapply(GSE.datasets, function(x){
#   # pdata <- NULL
#   # gse <- getGEO("GSE70486", GSEMatrix = TRUE)  # show(x)
#   gse <- getGEO(x,destdir=".", AnnotGPL=F, getGPL=F)
#   print(paste0('getGEO  for ', x,  " !"))
#   pdata <- pData(gse[[1]])
#   outfile <- paste0(data_directory, x)
#   if ( dir.exists(outfile)) {
#     print(paste0('getGEOSuppFiles  for ', x,  " is done, so skip !"))
#     # next;
#   } else{
#     # gSupp <- getGEOSuppFiles(x)
#   }
#   return(pdata)
#   })
# names(allDat) <- GSE.datasets
# save(allDat, file = paste0('/data2/rluo4/EpiTrans/RMDatasets/RM_pdata_GEO.RData'))
################################################################################
# 2) below steps of getGEOSuppFiles() is done in rstudio one-by-one 
################################################################################
# 3) continue to download the main RMDatasets in GEO #
################################################################################
load(paste0('/data2/rluo4/EpiTrans/RMDatasets/RM_pdata_GEO.RData'))
gse <- getGEO('GSE90963',destdir=".", AnnotGPL=F, getGPL=F)
pdata <- pData(gse[[2]])
View(pData(gse[[2]]))
allDat$GSE90963 <- pdata

gse <- getGEO('GSE97909',destdir=".", AnnotGPL=F, getGPL=F)
pdata <- pData(gse[[1]])
# gSupp <- getGEOSuppFiles('GSE97909')
allDat$GSE97908 <- pdata
names(allDat)[names(allDat)=='GSE97908'] <- 'GSE97909'

setwd(data_directory)
gse <- getGEO('GSE124509',destdir=".", AnnotGPL=F, getGPL=F)
pdata <- pData(gse[[2]])
allDat$GSE124509 <- pdata

setwd(data_directory)
gse <- getGEO('GSE38957',destdir=".", AnnotGPL=F, getGPL=F)
pdata <- pData(gse[[1]])
# gSupp <- getGEOSuppFiles('GSE38957')
allDat$GSE38958 <- pdata
names(allDat)[names(allDat)=='GSE38958'] <- 'GSE38957'

setwd(data_directory)
gse <- getGEO('GSE155413',destdir=".", AnnotGPL=F, getGPL=F)
pdata <- pData(gse[[1]])
View(pData(gse[[2]]))
# gSupp <- getGEOSuppFiles('GSE155413')
setwd('GSE155413')
# gSupp <- getGEOSuppFiles('GSE155447')
allDat$GSE155448 <- pdata
names(allDat)[names(allDat)=='GSE155448'] <- 'GSE155413'

setwd(data_directory)
gse <- getGEO('GSE93750',destdir=".", AnnotGPL=F, getGPL=F)
pdata <- pData(gse[[1]])
# gSupp <- getGEOSuppFiles('GSE93750')
allDat$GSE93751 <- pdata
names(allDat)[names(allDat)=='GSE93751'] <- 'GSE93750'

setwd(data_directory)
gse <- getGEO('GSE149989',destdir=".", AnnotGPL=F, getGPL=F)
pdata <- pData(gse[[1]])
# gSupp <- getGEOSuppFiles('GSE149989')
allDat$GSE150076 <- pdata
names(allDat)[names(allDat)=='GSE150076'] <- 'GSE149989'

setwd(data_directory)
gse <- getGEO('GSE90684',destdir=".", AnnotGPL=F, getGPL=F)
pdata <- pData(gse[[1]])
setwd('GSE90642')
# gSupp <- getGEOSuppFiles('GSE90684')
# gSupp <- getGEOSuppFiles('GSE90639')
allDat$GSE90642 <- pdata
names(allDat)[names(allDat)=='GSE90642'] <- 'GSE90684'

setwd(data_directory)
gse <- getGEO('GSE240674',destdir=".", AnnotGPL=F, getGPL=F)
pdata <- pData(gse[[1]])
# gSupp <- getGEOSuppFiles('GSE240674')
allDat$GSE240674 <- pdata

setwd(data_directory)
gse <- getGEO('GSE171497',destdir=".", AnnotGPL=F, getGPL=F)
pdata <- pData(gse[[1]])
# gSupp <- getGEOSuppFiles('GSE171497')
allDat$GSE171497 <- pdata
setwd(data_directory)
gse <- getGEO('GSE171227',destdir=".", AnnotGPL=F, getGPL=F)
pdata <- pData(gse[[1]])
setwd('GSE171497')
# gSupp <- getGEOSuppFiles('GSE171227')
View(pData(gse[[2]]))

setwd(data_directory)
gse <- getGEO('GSE254232',destdir=".", AnnotGPL=F, getGPL=F)
pdata <- pData(gse[[1]])
# gSupp <- getGEOSuppFiles('GSE254232')
allDat$GSE254232 <- pdata

setwd(data_directory)
gse <- getGEO('GSE252752',destdir=".", AnnotGPL=F, getGPL=F)
pdata <- pData(gse[[1]])
# gSupp <- getGEOSuppFiles('GSE252752')
allDat$GSE252752 <- pdata

setwd(data_directory)
gse <- getGEO('GSE242276',destdir=".", AnnotGPL=F, getGPL=F)
pdata <- pData(gse[[1]])
# gSupp <- getGEOSuppFiles('GSE242276')
allDat$GSE242276 <- pdata

setwd(data_directory)
gse <- getGEO('GSE223731',destdir=".", AnnotGPL=F, getGPL=F)
pdata <- pData(gse[[2]])
View(pData(gse[[1]]))
# gSupp <- getGEOSuppFiles('GSE223731')
allDat$GSE223731 <- pdata
setwd(data_directory)
gse <- getGEO('GSE223728',destdir=".", AnnotGPL=F, getGPL=F)
pdata <- pData(gse[[1]])
setwd('GSE223731')
# gSupp <- getGEOSuppFiles('GSE223728')
setwd(data_directory)
gse <- getGEO('GSE223730',destdir=".", AnnotGPL=F, getGPL=F)
pdata <- pData(gse[[1]])
setwd('GSE223731')
# gSupp <- getGEOSuppFiles('GSE223730')

setwd(data_directory)
gse <- getGEO('GSE226129',destdir=".", AnnotGPL=F, getGPL=F)
pdata <- pData(gse[[1]])
# gSupp <- getGEOSuppFiles('GSE226129')
allDat$GSE226129 <- pdata

setwd(data_directory)
gse <- getGEO('GSE163310',destdir=".", AnnotGPL=F, getGPL=F)
pdata <- pData(gse[[1]])
# gSupp <- getGEOSuppFiles('GSE163310')
allDat$GSE163310 <- pdata

setwd(data_directory)
gse <- getGEO('GSE195637',destdir=".", AnnotGPL=F, getGPL=F)
pdata <- pData(gse[[1]])
# gSupp <- getGEOSuppFiles('GSE195637')
allDat$GSE195637 <- pdata
setwd(data_directory)
gse <- getGEO('GSE195703',destdir=".", AnnotGPL=F, getGPL=F)
pdata <- pData(gse[[1]])
setwd('GSE195637')
# gSupp <- getGEOSuppFiles('GSE195703')

setwd(data_directory)
gse <- getGEO('GSE144620',destdir=".", AnnotGPL=F, getGPL=F)
pdata <- pData(gse[[1]])
# gSupp <- getGEOSuppFiles('GSE144620')
allDat$GSE144620 <- pdata

setwd(data_directory)
gse <- getGEO('GSE141994',destdir=".", AnnotGPL=F, getGPL=F)
pdata <- pData(gse[[1]])
# gSupp <- getGEOSuppFiles('GSE141994')
allDat$GSE141994 <- pdata
setwd(data_directory)
gse <- getGEO('GSE159551',destdir=".", AnnotGPL=F, getGPL=F)
pdata <- pData(gse[[1]])
setwd('GSE141994')
# gSupp <- getGEOSuppFiles('GSE159551')
allDat <- allDat[-which(names(allDat) == "GSE141993")]

setwd(data_directory)
gse <- getGEO('GSE133517',destdir=".", AnnotGPL=F, getGPL=F)
pdata <- pData(gse[[1]])
# gSupp <- getGEOSuppFiles('GSE133517')
allDat$GSE133517 <- pdata

setwd(data_directory)
gse <- getGEO('GSE142386',destdir=".", AnnotGPL=F, getGPL=F)
pdata <- pData(gse[[1]])
# gSupp <- getGEOSuppFiles('GSE142386')
allDat$GSE142386 <- pdata

setwd(data_directory)
gse <- getGEO('GSE79577',destdir=".", AnnotGPL=F, getGPL=F)
pdata <- pData(gse[[1]])
# gSupp <- getGEOSuppFiles('GSE79577')
allDat$GSE79577 <- pdata

setwd(data_directory)
gse <- getGEO('GSE86214',destdir=".", AnnotGPL=F, getGPL=F)
pdata <- pData(gse[[1]])
# gSupp <- getGEOSuppFiles('GSE86214')
allDat$GSE86214 <- pdata
View(pData(gse[[2]]))

setwd(data_directory)
gse <- getGEO('GSE63591',destdir=".", AnnotGPL=F, getGPL=F)
pdata <- pData(gse[[1]])
# gSupp <- getGEOSuppFiles('GSE63591')
allDat$GSE63591 <- pdata
View(pData(gse[[2]]))

setwd(data_directory)
gse <- getGEO('GSE49339',destdir=".", AnnotGPL=F, getGPL=F)
pdata <- pData(gse[[1]])
# gSupp <- getGEOSuppFiles('GSE49339')
allDat$GSE49339 <- pdata

setwd(data_directory)
gse <- getGEO('GSE40132',destdir=".", AnnotGPL=F, getGPL=F)
pdata <- pData(gse[[1]])
# gSupp <- getGEOSuppFiles('GSE40132')
allDat$GSE40132 <- pdata

# setwd(data_directory)
# gse <- getGEO('GSE180400',destdir=".", AnnotGPL=F, getGPL=F)
# pdata <- pData(gse[[1]])
# # gSupp <- getGEOSuppFiles('GSE180400')
# allDat$GSE180400 <- pdata
# allDat <- allDat[-which(names(allDat) == "GSE180400")]

setwd(data_directory)
gse <- getGEO('GSE174374',destdir=".", AnnotGPL=F, getGPL=F)
pdata <- pData(gse[[1]])
# gSupp <- getGEOSuppFiles('GSE174374')
allDat$GSE174374 <- pdata

setwd(data_directory)
gse <- getGEO('GSE174492',destdir=".", AnnotGPL=F, getGPL=F)
pdata <- pData(gse[[1]])
# gSupp <- getGEOSuppFiles('GSE174492')
allDat$GSE174492 <- pdata
View(pData(gse[[3]]))

setwd(data_directory)
gse <- getGEO('GSE71096',destdir=".", AnnotGPL=F, getGPL=F)
pdata <- pData(gse[[1]])
# gSupp <- getGEOSuppFiles('GSE71096')
allDat$GSE71096 <- pdata
setwd(data_directory)
gse <- getGEO('GSE71095',destdir=".", AnnotGPL=F, getGPL=F)
pdata <- pData(gse[[1]])
setwd('GSE71096')
# gSupp <- getGEOSuppFiles('GSE71095')

setwd(data_directory)
gse <- getGEO('GSE207643',destdir=".", AnnotGPL=F, getGPL=F)
pdata <- pData(gse[[1]])
# gSupp <- getGEOSuppFiles('GSE207643')
allDat$GSE207643 <- pdata

setwd(data_directory)
gse <- getGEO('GSE202815',destdir=".", AnnotGPL=F, getGPL=F)
pdata <- pData(gse[[1]])
setwd('GSE202815')
# gSupp <- getGEOSuppFiles('GSE202815')
allDat$GSE202815 <- pdata

getwd()
gse <- getGEO('GSE110323',destdir=".", AnnotGPL=F, getGPL=F)
pdata <- pData(gse[[1]])
setwd('GSE121942')
# gSupp <- getGEOSuppFiles('GSE121949')
# gSupp <- getGEOSuppFiles('GSE121952')
setwd(data_directory)
gse <- getGEO('GSE110320',destdir=".", AnnotGPL=F, getGPL=F)
pdata <- pData(gse[[1]])
setwd('GSE121942')
# gSupp <- getGEOSuppFiles('GSE110320')
# gSupp <- getGEOSuppFiles('GSE121948')
allDat$GSE121942 <- pdata

################################################################################
GSE.download <- list.dirs('/data2/rluo4/EpiTrans/RMDatasets/GEO')
GSE.download <- gsub("/data2/rluo4/EpiTrans/RMDatasets/GEO/", "", GSE.download)
missed.download <- setdiff(GSE.datasets, GSE.download)
# print(missed.download)
# missedDat <- lapply(missed.download, function(x){
#   # pdata <- NULL
#   # gse <- getGEO("GSE70486", GSEMatrix = TRUE)  # show(x)
#   gse <- getGEO(x,destdir=".", AnnotGPL=F, getGPL=F)
#   print(paste0('getGEO  for ', x,  " !"))
#   pdata <- pData(gse[[1]])
#   outfile <- paste0(data_directory, x)
#   if ( dir.exists(outfile)) {
#     print(paste0('getGEOSuppFiles  for ', x,  " is done, so skip !"))
#     # next;
#   } else{
#     # gSupp <- getGEOSuppFiles(x)
#   }
#   return(pdata)
# })
# missed.download <- setdiff(GSE.datasets, names(allDat))
# missedDat <- lapply(missed.download, function(x){
#   # pdata <- NULL
#   # gse <- getGEO("GSE70486", GSEMatrix = TRUE)  # show(x)
#   gse <- getGEO(x,destdir=".", AnnotGPL=F, getGPL=F)
#   print(paste0('getGEO  for ', x,  " !"))
#   pdata <- pData(gse[[1]])
#   outfile <- paste0(data_directory, x)
#   if ( dir.exists(outfile)) {
#     print(paste0('getGEOSuppFiles  for ', x,  " is done, so skip !"))
#     # next;
#   } else{
#     # gSupp <- getGEOSuppFiles(x)
#   }
#   return(pdata)
# })
# names(missedDat) <- missed.download
# allDat <- c(allDat, missedDat)
# allDat <- allDat[na.omit(match(datasetlists$GSE, names(allDat)))]
missed.download <- setdiff(RMDatasets$GSE, GSE.download)
missed.download #PRJNA498900
names(allDat)
################################################################################
# remove other species: c('GSE145686', 'GSE133138', 'GSE66090')[3,4,5]
allDat <- allDat[ ! names(allDat) %in% c('GSE145686', 'GSE133138', 'GSE66090')]
save(allDat, file = paste0('/data2/rluo4/EpiTrans/RMDatasets/allDat.RData'))
# gse <- getGEO('GSE130011',destdir=".", AnnotGPL=F, getGPL=F)
# pdata <- pData(gse[[1]])
# # gSupp <- getGEOSuppFiles('GSE130011')
# allDat$GSE145686 <- pdata

# cd /data2/rluo4/EpiTrans/RMDatasets/GEO/GSE148764
# gse <- getGEO('GSE148763',destdir=".", AnnotGPL=F, getGPL=F)
# pdata <- pData(gse[[1]])
# setwd('GSE148764')
# # gSupp <- getGEOSuppFiles('GSE148763')

# getwd()
# gse <- getGEO('GSE211076',destdir=".", AnnotGPL=F, getGPL=F)
# pdata <- pData(gse[[1]])
# setwd('GSE211442')
# # gSupp <- getGEOSuppFiles('GSE211076')
load(paste0('/data2/rluo4/EpiTrans/RMDatasets/allDat.RData'))
################################################################################
# 4) download the supplementary RMDatasets in GEO #
################################################################################
getwd()
gse <- getGEO('GSE37001',destdir=".", AnnotGPL=F, getGPL=F)
pdata <- pData(gse[[1]])
setwd('GSE37002')
# gSupp <- getGEOSuppFiles('GSE37001')

getwd()
setwd(data_directory)
gse <- getGEO('GSE87187',destdir=".", AnnotGPL=F, getGPL=F)
pdata <- pData(gse[[1]])
gse <- getGEO('GSE87189',destdir=".", AnnotGPL=F, getGPL=F)
pdata <- pData(gse[[1]])
setwd('GSE87190')
# gSupp <- getGEOSuppFiles('GSE87187')
# gSupp <- getGEOSuppFiles('GSE87189')

getwd()
setwd(data_directory)
gse <- getGEO('GSE93054',destdir=".", AnnotGPL=F, getGPL=F)
pdata <- pData(gse[[1]])
setwd('GSE87515')
# gSupp <- getGEOSuppFiles('GSE93054')

getwd()
setwd(data_directory)
gse <- getGEO('GSE84944',destdir=".", AnnotGPL=F, getGPL=F)
pdata <- pData(gse[[1]])
setwd('GSE76414')
# gSupp <- getGEOSuppFiles('GSE84944')
setwd(data_directory)
gse <- getGEO('GSE85008',destdir=".", AnnotGPL=F, getGPL=F)
pdata <- pData(gse[[1]])
setwd('GSE76414')
# gSupp <- getGEOSuppFiles('GSE85008')

getwd()
setwd(data_directory)
gse <- getGEO('GSE107956',destdir=".", AnnotGPL=F, getGPL=F)
pdata <- pData(gse[[1]])
setwd('GSE107954')
# gSupp <- getGEOSuppFiles('GSE107956')

setwd(data_directory)
gse <- getGEO('GSE112182',destdir=".", AnnotGPL=F, getGPL=F)
pdata <- pData(gse[[1]])
setwd('GSE112181')
# gSupp <- getGEOSuppFiles('GSE112182')
setwd(data_directory)
gse <- getGEO('GSE112180',destdir=".", AnnotGPL=F, getGPL=F)
pdata <- pData(gse[[1]])
setwd('GSE112181')
# gSupp <- getGEOSuppFiles('GSE112180')

setwd(data_directory)
gse <- getGEO('GSE122800',destdir=".", AnnotGPL=F, getGPL=F)
pdata <- pData(gse[[1]])
setwd('GSE122803')
# gSupp <- getGEOSuppFiles('GSE122800')
setwd(data_directory)
gse <- getGEO('GSE122802',destdir=".", AnnotGPL=F, getGPL=F)
pdata <- pData(gse[[1]])
setwd('GSE122803')
# gSupp <- getGEOSuppFiles('GSE122802')

setwd(data_directory)
gse <- getGEO('GSE132306',destdir=".", AnnotGPL=F, getGPL=F)
pdata <- pData(gse[[1]])
setwd('GSE120024')
# gSupp <- getGEOSuppFiles('GSE132306')

setwd(data_directory)
gse <- getGEO('GSE95372',destdir=".", AnnotGPL=F, getGPL=F)
pdata <- pData(gse[[1]])
setwd('GSE106122')
# gSupp <- getGEOSuppFiles('GSE95372')
setwd(data_directory)
gse <- getGEO('GSE105782',destdir=".", AnnotGPL=F, getGPL=F)
pdata <- pData(gse[[1]])
setwd('GSE106122')
# gSupp <- getGEOSuppFiles('GSE105782')

setwd(data_directory)
gse <- getGEO('GSE120659',destdir=".", AnnotGPL=F, getGPL=F)
pdata <- pData(gse[[1]])
setwd('GSE120860')
# gSupp <- getGEOSuppFiles('GSE120659')

setwd(data_directory)
gse <- getGEO('GSE129945',destdir=".", AnnotGPL=F, getGPL=F)
pdata <- pData(gse[[1]])
setwd('GSE130172')
# gSupp <- getGEOSuppFiles('GSE129945')

setwd(data_directory)
gse <- getGEO('GSE141991',destdir=".", AnnotGPL=F, getGPL=F)
pdata <- pData(gse[[1]])
setwd('GSE141993')
# gSupp <- getGEOSuppFiles('GSE141991')

setwd(data_directory)
gse <- getGEO('GSE128581',destdir=".", AnnotGPL=F, getGPL=F)
pdata <- pData(gse[[1]])
setwd('GSE128580')
# gSupp <- getGEOSuppFiles('GSE128581')

setwd(data_directory)
gse <- getGEO('GSE144959',destdir=".", AnnotGPL=F, getGPL=F)
pdata <- pData(gse[[1]])
setwd('GSE144984')
# gSupp <- getGEOSuppFiles('GSE144959')
setwd(data_directory)
gse <- getGEO('GSE144968',destdir=".", AnnotGPL=F, getGPL=F)
pdata <- pData(gse[[1]])
setwd('GSE144984')
# gSupp <- getGEOSuppFiles('GSE144968')

setwd(data_directory)
gse <- getGEO('GSE158740',destdir=".", AnnotGPL=F, getGPL=F)
pdata <- pData(gse[[1]])
setwd('GSE158742')
# gSupp <- getGEOSuppFiles('GSE158740')


setwd(data_directory)
gse <- getGEO('GSE161789',destdir=".", AnnotGPL=F, getGPL=F)
pdata <- pData(gse[[1]])
setwd('GSE161792')
# gSupp <- getGEOSuppFiles('GSE161789')

setwd(data_directory)
gse <- getGEO('GSE147884',destdir=".", AnnotGPL=F, getGPL=F)
pdata <- pData(gse[[1]])
setwd('GSE147891')
# gSupp <- getGEOSuppFiles('GSE147884')

setwd(data_directory)
gse <- getGEO('GSE161301',destdir=".", AnnotGPL=F, getGPL=F)
pdata <- pData(gse[[1]])
setwd('GSE161304')
# gSupp <- getGEOSuppFiles('GSE161301')
setwd(data_directory)
gse <- getGEO('GSE161302',destdir=".", AnnotGPL=F, getGPL=F)
pdata <- pData(gse[[1]])
setwd('GSE161304')
# gSupp <- getGEOSuppFiles('GSE161302')
setwd(data_directory)
gse <- getGEO('GSE161303',destdir=".", AnnotGPL=F, getGPL=F)
pdata <- pData(gse[[1]])
setwd('GSE161304')
# gSupp <- getGEOSuppFiles('GSE161303')

setwd(data_directory)
gse <- getGEO('GSE134101',destdir=".", AnnotGPL=F, getGPL=F)
pdata <- pData(gse[[1]])
setwd('GSE134103')
# gSupp <- getGEOSuppFiles('GSE134101')
setwd(data_directory)
gse <- getGEO('GSE134102',destdir=".", AnnotGPL=F, getGPL=F)
pdata <- pData(gse[[1]])
setwd('GSE134103')
# gSupp <- getGEOSuppFiles('GSE134102')

setwd(data_directory)
gse <- getGEO('GSE163500',destdir=".", AnnotGPL=F, getGPL=F)
pdata <- pData(gse[[1]])
setwd('GSE163491')
# gSupp <- getGEOSuppFiles('GSE163500')

setwd(data_directory)
gse <- getGEO('GSE66011',destdir=".", AnnotGPL=F, getGPL=F)
pdata <- pData(gse[[1]])
setwd('GSE66012')
# gSupp <- getGEOSuppFiles('GSE66011')

setwd(data_directory)
gse <- getGEO('GSE93749',destdir=".", AnnotGPL=F, getGPL=F)
pdata <- pData(gse[[1]])
setwd('GSE93750')
# gSupp <- getGEOSuppFiles('GSE93749')

setwd(data_directory)
gse <- getGEO('GSE210865',destdir=".", AnnotGPL=F, getGPL=F)
pdata <- pData(gse[[1]])
setwd('GSE210867')
# gSupp <- getGEOSuppFiles('GSE210865')
setwd(data_directory)
gse <- getGEO('GSE224671',destdir=".", AnnotGPL=F, getGPL=F)
pdata <- pData(gse[[1]])
setwd('GSE210867')
# gSupp <- getGEOSuppFiles('GSE224671')

setwd(data_directory)
gse <- getGEO('GSE217254',destdir=".", AnnotGPL=F, getGPL=F)
pdata <- pData(gse[[1]])
setwd('GSE217256')
# gSupp <- getGEOSuppFiles('GSE217254')
setwd(data_directory)
gse <- getGEO('GSE217255',destdir=".", AnnotGPL=F, getGPL=F)
pdata <- pData(gse[[1]])
setwd('GSE217256')
# gSupp <- getGEOSuppFiles('GSE217255')

setwd(data_directory)
gse <- getGEO('GSE130011',destdir=".", AnnotGPL=F, getGPL=F)
pdata <- pData(gse[[1]])
setwd('GSE213207')
# gSupp <- getGEOSuppFiles('GSE130011')

setwd(data_directory)
gse <- getGEO('GSE202848',destdir=".", AnnotGPL=F, getGPL=F)
pdata <- pData(gse[[1]])
View(pData(gse[[2]]))
View(pData(gse[[3]]))
setwd('GSE202848')
# gSupp <- getGEOSuppFiles('GSE202844')
# gSupp <- getGEOSuppFiles('GSE202846')
# gSupp <- getGEOSuppFiles('GSE228733')
# Read in the matrix file #

setwd(data_directory)
gse <- getGEO('GSE122254',destdir=".", AnnotGPL=F, getGPL=F)
pdata <- pData(gse[[1]])
setwd('GSE122254')
gSupp <- getGEOSuppFiles('GSE122260')
setwd(data_directory)
gse <- getGEO('GSE122259',destdir=".", AnnotGPL=F, getGPL=F)
pdata <- pData(gse[[2]])

setwd(data_directory)
gse <- getGEO('GSE125046',destdir=".", AnnotGPL=F, getGPL=F)
pdata <- pData(gse[[2]]) # batch 2,3,4 - 	Illumina HiSeq 4000
allDat[[i]] <- rbind(allDat[[i]], pdata)

setwd(data_directory)
gse <- getGEO('GSE83438',destdir=".", AnnotGPL=F, getGPL=F)
pdata <- pData(gse[[2]])
allDat[[i]] <- rbind(allDat[[i]], pdata)

setwd(data_directory)
gse <- getGEO('GSE55572',destdir=".", AnnotGPL=F, getGPL=F)
pdata <- pData(gse[[2]])

setwd(data_directory)
gse <- getGEO('GSE112181',destdir=".", AnnotGPL=F, getGPL=F)
pdata <- pData(gse[[2]])
gse <- getGEO('GSE120455',destdir=".", AnnotGPL=F, getGPL=F)
pdata <- pData(gse[[1]])

setwd(data_directory)
gse <- getGEO('GSE90684',destdir=".", AnnotGPL=F, getGPL=F)
pdata <- pData(gse[[2]])
gse <- getGEO('GSE90642',destdir=".", AnnotGPL=F, getGPL=F)
pdata <- pData(gse[[1]])

setwd(data_directory)
# gse <- getGEO('GSE106122',destdir=".", AnnotGPL=F, getGPL=F)
# pdata <- pData(gse[[2]])
gse <- getGEO('GSE117314',destdir=".", AnnotGPL=F, getGPL=F)
pdata <- pData(gse[[1]])
setwd('GSE106122')
gSupp <- getGEOSuppFiles('GSE117314')

setwd(data_directory)
# gse <- getGEO('GSE106122',destdir=".", AnnotGPL=F, getGPL=F)
# pdata <- pData(gse[[2]])
gse <- getGEO('GSE120659',destdir=".", AnnotGPL=F, getGPL=F)
setwd('GSE120860')
pdata <- pData(gse[[1]])
setwd(data_directory)
gSupp <- getGEOSuppFiles('GSE120659')
setwd(data_directory)
gse <- getGEO('GSE120611',destdir=".", AnnotGPL=F, getGPL=F)
setwd('GSE120860')
pdata <- pData(gse[[1]])
gSupp <- getGEOSuppFiles('GSE120611')


setwd(data_directory)
gse <- getGEO('GSE158020',destdir=".", AnnotGPL=F, getGPL=F)
setwd('GSE158742')
pdata <- pData(gse[[1]])
gSupp <- getGEOSuppFiles('GSE158020')

setwd(data_directory)
gse <- getGEO('GSE169589',destdir=".", AnnotGPL=F, getGPL=F)
setwd('GSE169589')
pdata <- pData(gse[[3]])

setwd(data_directory)
gse <- getGEO('GSE86214',destdir=".", AnnotGPL=F, getGPL=F)
setwd('GSE86214')
pdata <- pData(gse[[2]])

#n = 113
setwd(data_directory)
gse <- getGEO('GSE253795',destdir=".", AnnotGPL=F, getGPL=F)
pdata <- pData(gse[[1]])
allDat$GSE253795 <- pdata

setwd(data_directory)
gse <- getGEO('GSE211442',destdir=".", AnnotGPL=F, getGPL=F)
pdata <- pData(gse[[1]])
allDat$GSE211442 <- pdata

setwd(data_directory)
gse <- getGEO('GSE44386',destdir=".", AnnotGPL=F, getGPL=F)
pdata <- pData(gse[[1]])
allDat$GSE44386 <- pdata


setwd(data_directory)
gse <- getGEO('GSE73941',destdir=".", AnnotGPL=F, getGPL=F)
pdata <- pData(gse[[1]])
allDat$GSE73941 <- pdata

setwd(data_directory)
gse <- getGEO('GSE198643',destdir=".", AnnotGPL=F, getGPL=F)
pdata <- pData(gse[[1]])
allDat$GSE198643 <- pdata

setwd(data_directory)
gse <- getGEO('GSE78509',destdir=".", AnnotGPL=F, getGPL=F)
pdata <- pData(gse[[1]])
allDat$GSE78509 <- pdata
pdata <- pData(gse[[2]])

setwd(data_directory)
gse <- getGEO('GSE149989',destdir=".", AnnotGPL=F, getGPL=F)
pdata <- pData(gse[[1]])
allDat$GSE149989 <- pdata

setwd(data_directory)
gse <- getGEO('GSE210563',destdir=".", AnnotGPL=F, getGPL=F)
pdata <- pData(gse[[1]])
allDat$GSE210563 <- pdata
pdata <- pData(gse[[3]])

################################################################################
# 5) parse in supplementary RMDatasets in GEO #
################################################################################
x <- list.files(paste0(data_directory,'GSE106122/'))
x <- str_split(x[grep('cpm', x)], '[.]', simplify = T)[,1]
x
x <- paste(str_split(x, 'cpm_', simplify = T)[,2], collapse = ';')
x <- gsub('Sample_', '', gsub('data_', '', x))
x
extract_last_element <- function(x) {
  split_string <- strsplit(x, "_")
  last_element <- sapply(split_string, function(y) tail(y, n = 1))
  return(last_element)
}

x <- list.files(paste0(data_directory,'GSE106122/'))
x <- str_split(x[grep('HTSeqCts', x)], '[.]', simplify = T)[,1]
x

x <- read.table(paste0(data_directory,'GSE253795/GSE253795_A375_shALKBH5_inputVSA375_shCtrl_input_Gene_differential_expression.txt.gz'), fill =T, sep = '\t')
library(readxl)
x <-  read_excel(paste0(data_directory,'GSE93749/GSE93750_all-samples-reads-count.xls'),sheet = 1)
y <-  read_excel(paste0(data_directory,'GSE93749/GSE93750_all-samples-reads-count.xls'),sheet = 2)
x <-  read_excel(paste0(data_directory,'GSE98623/GSE98623_rnaseq.counts.xls'),sheet = 1)
x <-  read.table(paste0(data_directory,'GSE90914/GSE90914_input_RNAseq_fpkm.xls'), header = T)
colnames(x)
x <- read.xlsx(paste0(data_directory,'GSE87516/GSE87516_ZIKV-processed.xlsx'),sheet = 1)

x <- read.table(paste0(data_directory,'GSE87190/GSE87189_Drug1_vs_PBS1.all.gene.result.xls'), fill = T, sep = '\t', header = T)
x <- read.table(paste0(data_directory,'GSE87190/GSE87189_Drug2_vs_PBS2.all.gene.result.xls'), fill = T, sep = '\t', header = T)
x <- read.table(paste0(data_directory,'GSE87190/GSE87187_Group1_vs_Group2.all.gene.result.txt.gz'), sep = '\t', header = T)
x <- read.table(paste0(data_directory,'GSE87190/GSE87187_Group2_vs_Group3.all.gene.result.txt.gz'), fill = T, sep = '\t')
x[1,]
x <- read_excel(paste0(data_directory,'GSE38957/GSE38957_04-NSUN2_Aza-IP_VarScan_signature_analysis.xls'),sheet = 1)
x <- read_xlsx(paste0(data_directory,'GSE104947/GSE104947_Krogh_RMS_data_snRNA.xlsx'),sheet = 1)
x <- read_xlsx(paste0(data_directory,'GSE113798/pheno_data.xlsx'),sheet = 1)
x <- read_xlsx(paste0(data_directory,'GSE102336/GSE102336_ZCCHC4_polysome_profiling.xlsx'),sheet = 1)
x <- import("/data2/rluo4/EpiTrans/RMDatasets/GEO/GSE122803/GSM3485692_hg38-MEL624-Input-KO-RNA.bw")
x <- read_xlsx(paste0(data_directory,'GSE117299/GSE117299_RNA-seq_Polysome_profiling_20160428.xlsx'),sheet = 1)
x <- read_xlsx(paste0(data_directory,'GSE129716/GSE129716_Treat_inputVSNC_input_Transcript_differential_expression-PROCESSING_DATA.xlsx'),sheet = 1)
x <- read_xlsx(paste0(data_directory,'GSE137675/GSE137675_Reads_count_featureCounts.xlsx'),sheet = 1)
x <- read_xlsx(paste0(data_directory,'GSE137675/GSE137675_OCM1.OE.HINT2.and.OCM1.OE.NC.RNA-seq.processed.data.xlsx'),sheet = 1)
x <- read_xlsx(paste0(data_directory,'GSE137675/GSE137675_OCM1.sgALKBH5.and.OCM1.sgNC.RNA-seq.processed.data.xlsx'),sheet = 1)
x <- read_excel(paste0(data_directory,'GSE137675/GSE137675_PIG1.shMETTL3.and.PIG1.shNC.RNA-seq.processed.data.xls'),sheet = 1)
x <- read_xlsx(paste0(data_directory,'GSE171373/GSE171373_Methylated_RNA_sites.circRNA.xlsx'),sheet = 1)
x <- read_xlsx(paste0(data_directory,'GSE171472/GSE171472_All_Comparison_mRNA_.xlsx'), sheet = 1)
x <- read_xlsx(paste0(data_directory,'GSE154555/GSE154555_Methylated_RNA_sites.mRNA.xlsx'), sheet = 1)
x <- read_xlsx(paste0(data_directory,'GSE161792/GSE161789_raw_counts_and_FPKM_for_RNA-seq.xlsx'), sheet = 1)
x <- read_xlsx(paste0(data_directory,'GSE179042/GSE179042_4_genes_fpkm_expression.xlsx'), sheet = 1)
# x <- read_xlsx(paste0(data_directory,'GSE181540/GSE181540_mRNA_Expression_Profiling.xlsx')
x <- read_xlsx(paste0(data_directory,'GSE144032/GSE144032_m6A_young_called_peaks.xlsx'), sheet = 1)
x <- read_xlsx(paste0(data_directory,'GSE167075/GSE167075_shM3_vs_NTC_human_.xlsx'), sheet = 1)
x <- read_excel(paste0(data_directory,'GSE131316/GSE131316_HeLa.DEseq2.results.xlsx'), sheet = 1)
x <- read_excel(paste0(data_directory,'GSE38957/GSE38957_02-DNMT2_Aza-IP_All_Genes.xls'), sheet = 1)
x <- read_excel(paste0(data_directory,'GSE38957/GSE38957_03-NSUN2_Aza-IP_All_Genes.xls'), sheet = 2)
x <- read_excel(paste0(data_directory,'GSE38957/GSE38957_04-NSUN2_Aza-IP_VarScan_signature_analysis.xls'), sheet = 2)
x <- read_excel(paste0(data_directory,'GSE122413/GSE122413_Methylated_RNA_sites.xlsx'), sheet = 3)
x <- read.table('GSE169589_Polyribosome_mRNA_expression.txt.gz')
x <- read_excel(paste0(data_directory,'GSE240879/GSE240879_Gene_Expression_Profiling.xlsx'), sheet = 1)
x <- read_excel(paste0(data_directory,'GSE149989/GSE149989_LNZ308_METTL1_KD_Trac-seq.xlsx'), sheet = 1)
# x <- read.table(paste0(data_directory,'GSE240674/GSE240674_RNAseq_all_TPM.txt.gz')
x <- read_excel(paste0(data_directory,'GSE171497/GSE171227_ALKBH5-RNA-seq-FPKM.xls'), sheet = 1)
y <- read_excel(paste0(data_directory,'GSE171497/GSE171227_ALKBH5-RNA-seq-FPKM.xls'), sheet = 2)
x <- read.table(paste0(data_directory,'GSE254232/GSE254232_gene_count.xls'), fill = T, header = T, sep = '\t')
x <- read_xlsx(paste0(data_directory,'GSE252752/GSE252752_processed_data_file.xlsx'), sheet = 3)
paste0(colnames(x), collapse = ';')
x <- read.table(paste0(data_directory,'GSE242276/GSE242276_gene_input_all_cpm.txt.gz'), fill = T, header = T, sep = '\t')
x <- read.table(paste0(data_directory,'GSE223731/GSE223728_Hg3_METTL3KO_diff.txt.gz'), fill = T, header = T, sep = '\t')
x <- read_xlsx(paste0(data_directory,'GSE226129/GSE226129_Methylated_RNA_sites.mRNA.xlsx'),sheet = 3)
x <- read.table(paste0(data_directory,'GSE195637/GSE195703_mRNA-seq_All.fpkm.anno.xls'), fill = T, header = T, sep = '\t')
x <- read_xlsx(paste0(data_directory,'GSE142386/GSE142386_mRNA_Expression_Profiling.xlsx'), sheet = 1)
x <- read_xlsx(paste0(data_directory,'GSE86214/GSE86214_siYTHDF3_RiboProfile.xlsx'), sheet = 1)
y <- read_xlsx(paste0(data_directory,'GSE86214/GSE86214_siYTHDF3_YTHDF12RIP.xlsx'), sheet = 1)
x <- read_xlsx(paste0(data_directory,'GSE174374/GSE174374_dsRNA-seq_WT_vs_NSUN2-KO.xlsx'), sheet = 1)
y <- read_xlsx(paste0(data_directory,'GSE174374/GSE174374_RNA-seq_A549_WT_vs_NSUN2-KO.xlsx'), sheet = 1)
x <- read_xlsx(paste0(data_directory,'GSE213207/GSE130011_M2_count_RNAseq.xlsx'), sheet = 1)

# # To access the phenotypic information about the samples, the best way is to use getGEO() function to obtain the GSE object and then extract the phenoData object from that. Unfortunately this means downloadint the entire GSE Matrix file.
# dim(pData(gse[[1]]))
# head(pData(gse[[1]])[, 1:3])
# # Sometimes GSEs are include separate data tables with the sample information. If these exist, you can uuse the getGSEDataTables() function. For example here is the phenoData object from a different GSE accession GSE3494 with a Data Table.
# df1 <- getGSEDataTables("GSE102113")
# lapply(df1, head)
# # save(gset,pdata,file = "Quan_anno.RData")
# # load('Quan_anno.RData')
# # pdata.Quan <- pData(gset[[1]])
# # 
# # pdata.Guo$title <- gsub(', replicate ','_',gsub(',scRNAseq','',pdata.Guo$title))
# # pdata.Guo$SRA <- ascp_ena_link$run_accession[match(pdata.Guo$geo_accession,
# #                                                    ascp_ena_link$sample_alias) ]
# # write.table(pdata.Guo[,c(1,2, 44, 47)],"Guo_pdata.txt",
# #             quote = F,sep = ";",row.names = F,col.names = F)
# path = '/public/home/lorihan/lrh/database/Cervix/Guo/GSE208653'
RMDatasets <- read_excel('/data2/rluo4/EpiTrans/DataCollection/Regulatory-enzymes-RNA-modifications.xlsx', sheet = 1)
unique(RMDatasets$Modifcation)
RMDatasets$GSE[duplicated(RMDatasets$GSE)]
table(RMDatasets$GSE %in% names(allDat))
setdiff(RMDatasets$GSE, names(allDat))
table(RMDatasets$Layout)
table(RMDatasets$Modifcation)
table(RMDatasets$Treatment.Condition)
length(unique(RMDatasets$GSE))
# RMDatasets$GSE[duplicated(RMDatasets$GSE)]
RMDatasets$SRP[duplicated(RMDatasets$SRP)]
unique(RMDatasets$PMID)#168
unique(RMDatasets$GSE)#177
RMP <- read.csv(file = '/data2/rluo4/EpiTrans/DataCollection/RMP.csv')
# Identify the rows with specific characters in the "Gene" column
matching_rows <- grepl(paste(RMP$RMP_update, collapse = "|"), RMDatasets$Treatment.Condition)
# Print the dataframe with only the matching rows
# View(RMDatasets[matching_rows, ])
# View(RMDatasets[! matching_rows, ])
matching_rows1 <- grepl(paste(RMP$RMP_update, collapse = "|"), RMDatasets$Sample)
# Print the dataframe with only the matching rows
# View(RMDatasets[! (matching_rows|matching_rows1), ])
RMP_datasets <- RMDatasets[ matching_rows|matching_rows1, ]
unique(RMP_datasets$GSE)
table(RMP_datasets$Layout)
Layout <- gsub('count|FPKM|RPKM|result','matrix',RMP_datasets$Layout)
table(Layout)
grep("GSE|and|PAIRED matrix", Layout, value = TRUE)
index = grep("GSE|and|PAIRED matrix", Layout)#, value = TRUE)
length(index) #77
# View(RMP_datasets[-index,])

setwd(data_directory)
GEO_dir <- unique(RMP_datasets$GSE)
# GEO_dir <- gsub('./', '',list.dirs('.'))
# GEO_dir <- RMP_datasets$GSE[RMP_datasets$GSE %in% GEO_dir]

# RMP_datasets are GEO datasets with peak callings
RMP_GEO <- vector("list", length(GEO_dir))
names(RMP_GEO) <- GEO_dir
for (i in GEO_dir) {
  RMP_GEO[[i]] <- list.files(i)
}
RMP_datasets$peak <- apply(RMP_datasets, 1, function(x){
  i <- x['GSE']
  level3 <- grep("eak|bed|bw|wig", RMP_GEO[[i]], value = TRUE)
  level3 <- ifelse(length(level3)>=1, 'peakcalling','nopeak')
  return(level3)
})
table(RMP_datasets$peak)

################################################################################
# 6) ascp RMDatasets raw data in ENA #
################################################################################
setwd(data_directory)
GEO.dt <- list.dirs(data_directory)
GEO.dt <- gsub('/data2/rluo4/EpiTrans/RMDatasets/GEO//', '', GEO.dt)
GEO.dt <- GEO.dt[grep('GSE', GEO.dt)]
table(RMDatasets$GSE %in% GEO.dt)

setwd('/data2/rluo4/RPMfunc/xuhaixia/mismatch_RRM/')
a <- list.files('./')
a <- a[! a %in% 'RMP_nomatched.xlsx']
b <- NULL
for(i in a){
  setwd('/data2/rluo4/RPMfunc/xuhaixia/mismatch_RRM/')
  setwd(i)
  c <- list.files('./', pattern = '.txt')
  b <- c(b, c)
}
b[duplicated(b)]
setwd('/data2/rluo4/EpiTrans/DataCollection/RPMfunc_tsv')
d <- list.files('./', pattern = '.txt')
setdiff(b, d)
table(RMDatasets$Layout)
Layout <- gsub('count|FPKM|RPKM|result','matrix',RMDatasets$Layout)
table(Layout)
grep("GSE|and|PAIRED matrix", Layout, value = TRUE)
index = grep("GSE|and|PAIRED matrix", Layout)#, value = TRUE)
length(index) #99
RMDatasets$SRP[-index]
setdiff(RMDatasets$SRP[-index], gsub('.txt','',b))
setdiff(RMDatasets$SRP, gsub('.txt','',b)) # "SRP313533": 2024 PMID:?
# GSE171497	NA	NA	NA	U87;Glioblastoma	m6A	shALKBH5;shScr;Normoxia;Hypoxia	MeRIP-seq	shALKBH5-H;shScr-H;U87 WT-N;U87 WT-H	SRP313533
data_dir = '/data2/rluo4/EpiTrans/DataCollection/RPMfunc_tsv'
setwd(data_dir)
library(stringr)
RMP_GEO <- list.files('./', pattern='.txt')
RMP_GEO.acc <- str_split(gsub('.txt','',RMP_GEO), '_', simplify = T)[,1]
unique(RMP_GEO.acc)#79 --> 178
library(data.table)
filePath <- sapply(RMP_GEO, function(x){ 
  paste(data_dir,x,sep='/')})  
RMP_GEO <- lapply(filePath, function(x){
  fread(x)})
RMP_GEO$SRP064176.txt$fastq_aspera[RMP_GEO$SRP064176.txt$sample_alias=='GSM2011452'] <- "fasp.sra.ebi.ac.uk:/vol1/fastq/SRR306/000/SRR3066590/SRR3066590.fastq.gz"
RMP_GEO$SRP121512.txt$fastq_aspera[RMP_GEO$SRP121512.txt$sample_alias=='GSM2830063'] <- "fasp.sra.ebi.ac.uk:/vol1/fastq/SRR621/004/SRR6211494/SRR6211494_1.fastq.gz;fasp.sra.ebi.ac.uk:/vol1/fastq/SRR621/004/SRR6211494/SRR6211494_2.fastq.gz"

RMP_allDat <- NULL#data.frame()
RMP_sum <- lapply(RMP_GEO,function(y){
  RMP_allDat <<- rbind(RMP_allDat, y)
  return(RMP_allDat)  
})
unique(RMP_allDat$submitted_aspera)
unique(RMP_allDat$sample_alias)# ? sample for raw data 
colnames(RMP_allDat)
RMP_ENA <- RMP_allDat

table(RMP_ENA$scientific_name)
RMP_allDat <- RMP_ENA[RMP_ENA$scientific_name=='Homo sapiens',]
table(unique(RMP_allDat$study_alias) %in% RMDatasets$GSE)
table(gsub('.txt','',names(RMP_GEO)) %in% RMDatasets$SRP)
setdiff(unique(RMP_allDat$study_alias), RMDatasets$GSE)
index = duplicated(RMP_allDat$run_accession)
# View(RMP_allDat[index,])
RMP_allDat <- RMP_allDat[! index, ]

table(is.na(RMP_allDat$sample_alias))
table(is.na(RMP_allDat$experiment_alias))
index = is.na(RMP_allDat$sample_alias)
unique(RMP_allDat$study_alias[index])
RMP_allDat$sample_alias[index] <- str_split(RMP_allDat$experiment_alias[index], '_', simplify = T)[,1] 
table(is.na(RMP_allDat$sample_alias))
length(unique(RMP_allDat$sample_alias)) # 5351
table(duplicated(RMP_allDat$run_accession))
index = duplicated(RMP_allDat$sample_alias)
# View(RMP_allDat[index,])
RMP_allDat <- RMP_allDat[! index, ]
table(is.na(RMP_allDat$sample_title))
index = is.na(RMP_allDat$sample_title)
unique(RMP_allDat$study_alias[index])
# View(RMP_allDat[index,])
RMP_allDat$sample_title[index] <- str_split(RMP_allDat$experiment_alias[index], '_', simplify = T)[,1] 
table(duplicated(RMP_allDat$sample_title))
# View(RMP_allDat[duplicated(RMP_allDat$sample_title),])

RMP_allDat <-RMP_allDat[, c('study_alias','sample_alias','fastq_aspera','submitted_aspera','sra_aspera','sample_title')]
RMP_allDat$submitted_aspera[is.na(RMP_allDat$submitted_aspera)] <- RMP_allDat$fastq_aspera[is.na(RMP_allDat$submitted_aspera)]
RMP_allDat$submitted_aspera
RMP_allDat$submitted_aspera[RMP_allDat$submitted_aspera==''] <- RMP_allDat$fastq_aspera[RMP_allDat$submitted_aspera=='']
RMP_allDat$submitted_aspera
RMP_allDat$datatype[grepl('.bam', RMP_allDat$submitted_aspera)] <- 'Bam'
RMP_allDat$datatype[grepl('_1.fastq', RMP_allDat$submitted_aspera)] <- 'FastQ'
RMP_allDat$datatype[is.na(RMP_allDat$datatype)] <- 'Sra'
index = RMP_allDat$datatype=='Sra'
RMP_allDat$submitted_aspera[index] <- RMP_allDat$sra_aspera[index]
RMP_allDat_GSE <- paste(RMP_allDat$study_alias, RMP_allDat$study_alias, sep = ';')
RMP_allDat_GSE[index] <- RMP_allDat$study_alias[index]
RMP_allDat_GSE 
unique(RMP_allDat$study_alias) #80 -->198 -->200-->206
table(unique(RMP_allDat$study_alias) %in% RMDatasets$GSE)
setdiff(unique(RMP_allDat$study_alias), RMDatasets$GSE) # 52-->54-->57
# PRJNA865373.txt: "GSE66011"    "GSE66010"    "GSE129842"   "GSE189261"   "GSE202814"   "PRJNA865373"
# SRP093982.txt: "PRJNA355164" 
# SRP168947.txt: "GSE122600"
RMP_allDat_ascp <- unlist(str_split(RMP_allDat$submitted_aspera,";"))#submitted_aspera#fastq_ftp#
table(is.na(RMP_allDat_ascp))
RMP_allDat_GSM <- paste(RMP_allDat$sample_alias, RMP_allDat$sample_alias, sep = ';')
RMP_allDat_GSM[index] <- RMP_allDat$sample_alias[index]
RMP_allDat_GSM #2875 -->5455
RMP_allDat_dt <- paste(RMP_allDat$datatype, RMP_allDat$datatype, sep = ';')
RMP_allDat_dt[index] <- RMP_allDat$datatype[index]
RMP_allDat_dt
asp_sc <- data.frame(study_alias = RMP_allDat_GSE,
                     sample_alias = RMP_allDat_GSM, 
                     sra_aspera = RMP_allDat$submitted_aspera, datatype = RMP_allDat_dt
)
asp_sc$n1 <-  apply(asp_sc, 1, function(x){
  # n1 = sum(grepl(";", x['study_alias']))
  n1 <- str_count(x['study_alias'], ";")
  return(n1)
})
asp_sc$n2 <-  apply(asp_sc, 1, function(x){
  # n1 = sum(grepl(";", x['sample_alias']))
  n1 <- str_count(x['sample_alias'], ";")
  return(n1)
})
asp_sc$n3 <-  apply(asp_sc, 1, function(x){
  # n1 = sum(grepl(";", x['sra_aspera']))
  n1 <- str_count(x['sra_aspera'], ";")
  return(n1)
}) # I found the row 'GSM2011452' column 'fastq_aspera' in RMP_GEO["SRP064176.txt"] is: 
# fasp.sra.ebi.ac.uk:/vol1/fastq/SRR306/000/SRR3066590/SRR3066590.fastq.gz;fasp.sra.ebi.ac.uk:/vol1/fastq/SRR306/000/SRR3066590/SRR3066590_1.fastq.gz;fasp.sra.ebi.ac.uk:/vol1/fastq/SRR306/000/SRR3066590/SRR3066590_2.fastq.gz
# and found the row 'GSM2830063' column 'fastq_aspera' in RMP_GEO["SRP121512.txt"] is: 
# fasp.sra.ebi.ac.uk:/vol1/fastq/SRR621/004/SRR6211494/SRR6211494.fastq.gz;fasp.sra.ebi.ac.uk:/vol1/fastq/SRR621/004/SRR6211494/SRR6211494_1.fastq.gz;fasp.sra.ebi.ac.uk:/vol1/fastq/SRR621/004/SRR6211494/SRR6211494_2.fastq.gz
max(asp_sc$n3)
asp_sc <- data.frame(study_alias = unlist(str_split(RMP_allDat_GSE,";")), 
                     sample_alias = unlist(str_split(RMP_allDat_GSM,";")), 
                     sra_aspera = RMP_allDat_ascp, datatype = unlist(str_split(RMP_allDat_dt,";"))
)
write.table(asp_sc,file="/data2/rluo4/EpiTrans/DataCollection/asp_RMP_epitrans.link.txt",sep = "\t",quote = F, col.names = F, row.names = F)

