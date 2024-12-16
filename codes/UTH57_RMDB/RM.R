#!/usr/local/bin/Rscript
# title: "RMDisease V2.0: an updated database of genetic variants that affect RNA modifications with disease and trait implication"
# title:
# author: "Ruihan Luo"
# date: "March 15th,2024"
# rm(list=ls())
library(openxlsx)
# 1) RMDisease database
ac4c.disease <- read.csv('/data2/rluo4/EpiTrans/DataCollection/RMDisease/ac4c_human_associatedSNPs.csv')
am.disease <- read.csv('/data2/rluo4/EpiTrans/DataCollection/RMDisease/am_human_associatedSNPs.csv')
atoi.disease <- read.csv('/data2/rluo4/EpiTrans/DataCollection/RMDisease/a-to-i_human_associatedSNPs.csv')
cm.disease <- read.csv('/data2/rluo4/EpiTrans/DataCollection/RMDisease/cm_human_associatedSNPs.csv')
d.disease <- read.csv('/data2/rluo4/EpiTrans/DataCollection/RMDisease/d_pombe_associatedSNPs.csv')
f5c.disease <- read.csv('/data2/rluo4/EpiTrans/DataCollection/RMDisease/f5c_yeast_associatedSNPs.csv')
gm.disease <- read.csv('/data2/rluo4/EpiTrans/DataCollection/RMDisease/gm_human_associatedSNPs.csv')
hm5c.disease <- read.csv('/data2/rluo4/EpiTrans/DataCollection/RMDisease/hm5c_fly_associatedSNPs.csv')
m1a.disease <- read.csv('/data2/rluo4/EpiTrans/DataCollection/RMDisease/m1a_human_associatedSNPs.csv')
m5c.disease <- read.csv('/data2/rluo4/EpiTrans/DataCollection/RMDisease/m5c_human_associatedSNPs.csv')
m5u.disease <- read.csv('/data2/rluo4/EpiTrans/DataCollection/RMDisease/m5u_human_associatedSNPs.csv')
m6a.disease <- read.csv('/data2/rluo4/EpiTrans/DataCollection/RMDisease/m6a_human_associatedSNPs.csv')
m6am.disease <- read.csv('/data2/rluo4/EpiTrans/DataCollection/RMDisease/m6am_human_associatedSNPs.csv')
m7g.disease <- read.csv('/data2/rluo4/EpiTrans/DataCollection/RMDisease/m7g_human_associatedSNPs.csv')
psi.disease <- read.csv('/data2/rluo4/EpiTrans/DataCollection/RMDisease/psi_human_associatedSNPs.csv')
um.disease <- read.csv('/data2/rluo4/EpiTrans/DataCollection/RMDisease/um_human_associatedSNPs.csv')

ls()
d <- ls()[grep('.disease',ls())]
d
RM.disease <- bind_rows(ac4c.disease, am.disease, atoi.disease, cm.disease, d.disease, f5c.disease,
                                 gm.disease, hm5c.disease, m1a.disease, m5c.disease, m5u.disease, m6a.disease,
                                 m6am.disease, m7g.disease, psi.disease, um.disease)
extract_second_char <- function(x) {
  sapply(strsplit(x, ";"), function(y) substr(y[2], 1, 1))
}
extract_second_word <- function(x) {
  sapply(strsplit(x, ";"), function(y) if(length(y) >= 2) y[2] else NA)
}
# vector <- c("GSE1234;CRA5432", "CRA9876_XYZ;ABC_GSE5678")
# Function to extract words starting with "GSE" or "CRA" from each string
extract_GSE_CRA <- function(x) {
  words <- unlist(strsplit(x, "[;_]+"))  # Split by semicolons, underscores, or other delimiters
  matches <- grep("^GSE|^CRA", words, value = TRUE)
  if (length(matches) == 0) {
    return(NA)
  } else {
    return(matches)
  }
}
extract_GSE_CRA <- function(x) {
  words <- unlist(strsplit(x, "[;_]+"))  # Split by semicolons, underscores, or other delimiters
  matches <- grep("GSE|CRA", words, value = TRUE)
  if (length(matches) == 0) {
    return(NA)
  } else {
    return(matches)
  }
}
# Apply the function to each element of the vector
RMD.datasets <- sapply(RM.disease$MD_Source, extract_GSE_CRA)
# RMD.datasets <- extract_second_word(RM.disease$MD_Source)
# unique(RMD.datasets)
dd <- unique(unlist(RMD.datasets))
dd <- unique(unlist(strsplit(dd, "[|]")))
dd
dd <- grep("GSE|CRA", dd, value = TRUE)
RMD.datasets <- gsub('GSA:','', dd)
# d <- str_split(RM.disease$MD_Source, "[;_]+", simplify = T)
# dd <- unlist(strsplit(RM.disease$MD_Source, "[;_]+"))
# dd <- unlist(strsplit(dd, "|"))
# unique(dd)
# dim(d)
# d0 <- unique(d[,1]);d1 <- unique(d[,2]);d2 <- unique(d[,3]);d3 <- unique(d[,4]);d4 <- unique(d[,5]);d5 <- unique(d[,6]);d6 <- unique(d[,7]);d7 <- unique(d[,8]);d8 <- unique(d[,9])
test <- RM.disease[grep("GSE90164", RM.disease$MD_Source),]
length(unique(test$Gene))
test <- RMBase[grep("GSE90164", RMBase$supportList),]
length(unique(test$geneName))
table(test$cellList)

RMD.HighRes <- read.xlsx('/data2/rluo4/EpiTrans/DataCollection/RMDisease/Experimental-Validated.xlsx',sheet = 1)
RMD.LowRes <- read.xlsx('/data2/rluo4/EpiTrans/DataCollection/RMDisease/Experimental-Validated.xlsx',sheet = 2)
RMD.COVID <- read.xlsx('/data2/rluo4/EpiTrans/DataCollection/RMDisease/Experimental-Validated.xlsx',sheet = 3)
RMD.validation <- unique(c(RMD.HighRes$GSE[ grep("human", RMD.HighRes$Experiment.ID) ],
                           RMD.LowRes$GSE[ grep("human", RMD.LowRes$Experiment.ID) ], RMD.COVID$GSE))
RMD.validation <- na.omit(RMD.validation)

table(RMD.validation %in% RMD.datasets)
setdiff(RMD.validation, RMD.datasets)
setdiff(RMD.datasets, RMD.validation)

# 2) RMBase database
Header <- read.table('/data2/rluo4/EpiTrans/DataCollection/RMBase/RMBaseSNVHeader.txt', fill = T, sep = ':')
Header <- str_split(Header$V2, " -> ", simplify = T)[,1]
Header <- gsub(' ', '', Header)

# Ψ#(Pseudouridine) 
# RMBaseV3 <- read.table('/data2/rluo4/EpiTrans/DataCollection/RMBase/trans.hg38.tar.gz', sep = '\t', fill = T)
m1A <- read.table('/data2/rluo4/EpiTrans/DataCollection/RMBase/hg38.m1A.tar.gz', sep = '\t', fill = T)
m5C <- read.table('/data2/rluo4/EpiTrans/DataCollection/RMBase/hg38.m5C.tar.gz', sep = '\t', fill = T)
Pseudo <- read.table('/data2/rluo4/EpiTrans/DataCollection/RMBase/hg38.Pseudo.tar.gz', sep = '\t', fill = T)
m7G <- read.table('/data2/rluo4/EpiTrans/DataCollection/RMBase/hg38.m7G.tar.gz', sep = '\t', fill = T)
Nm <- read.table('/data2/rluo4/EpiTrans/DataCollection/RMBase/hg38.Nm.tar.gz', sep = '\t', fill = T)
otherMod <- read.table('/data2/rluo4/EpiTrans/DataCollection/RMBase/hg38.otherMod.tar.gz', sep = '\t', fill = T)
m6A <- read.table('/data2/rluo4/EpiTrans/DataCollection/RMBase/hg38.m6A.tar.gz', sep = '\t', fill = T)
RNAediting <- read.table('/data2/rluo4/EpiTrans/DataCollection/RMBase/hg38.RNA-editing.tar.gz', sep = '\t', fill = T)

Header <- read.table('/data2/rluo4/EpiTrans/DataCollection/RMBase/RMBaseV3Header.txt', fill = T, sep = ':')
Header <- str_split(Header$V2, " -> ", simplify = T)[,1]
Header <- gsub(' ', '', Header)
m6A <- m6A[-1,]
colnames(m6A) <- Header
m6A$motifScore <- as.character(m6A$motifScore)

RNAediting <- RNAediting[-1,]
colnames(RNAediting) <- Header
m1A <- m1A[-1,]
colnames(m1A) <- Header
m5C <- m5C[-1,]
colnames(m5C) <- Header
Pseudo <- Pseudo[-1,]
colnames(Pseudo) <- Header
m7G <- m7G[-1,]
colnames(m7G) <- Header
Nm <- Nm[-1,]
colnames(Nm) <- Header
otherMod <- otherMod[-1,]
colnames(otherMod) <- Header

RMBase <- bind_rows(m6A, m1A, m5C, Pseudo, m7G, Nm, otherMod)
RMBase.datasets <- unique(unlist(strsplit(RMBase$supportList, "[,_]+")))
RMBase.datasets <- grep("GSE|CRA", RMBase.datasets, value = TRUE)
intersect(RMBase.datasets, RMD.datasets)
# 3) m6AAtlas database
# m6A.datasets <- read.table('/data2/rluo4/EpiTrans/DataCollection/m6AAtlas/m6A-seq-GEO.txt', fill = T, header = T)
m6A.datasets1 <- read.xlsx('/data2/rluo4/EpiTrans/DataCollection/m6AAtlas/m6AAtlas.xlsx',sheet = 1)
m6A.datasets2 <- read.xlsx('/data2/rluo4/EpiTrans/DataCollection/m6AAtlas/m6AAtlas.xlsx',sheet = 2)
m6A.datasets2$GSE <- 'CRA001315'
m6A.datasets3 <- read.xlsx('/data2/rluo4/EpiTrans/DataCollection/m6AAtlas/m6AAtlas.xlsx',sheet = 3)
m6A.datasets <- unique(c(m6A.datasets1$GSE[grep('Homo', m6A.datasets1$Species)]
                         ,  m6A.datasets2$GSE[grep('Homo', m6A.datasets2$Species)]
                         , m6A.datasets3$GSE[grep('Homo', m6A.datasets3$Species)]
                         ))
intersect(RMBase.datasets, m6A.datasets)
intersect(RMD.datasets, m6A.datasets)
m5C.Atlas <- read.csv('/data2/rluo4/EpiTrans/DataCollection/m6AAtlas/m5CAtlas.csv', fill = T, header = T)
m5C.datasets <- unique(m5C.Atlas$GSE)

RM.datasets <- union(union(union(RMD.datasets, RMBase.datasets), m6A.datasets), m5C.datasets)                        
table(m5C.datasets %in% RM.datasets)
RM.datasets <- RM.datasets[RM.datasets !='NO'] # 161 datasets of human 
table(RMD.validation %in% RM.datasets)
setdiff(RMD.validation, RM.datasets)
RM.datasets <- union(RM.datasets,RMD.validation ) # 163 datasets of human 
# write.table(RM.datasets, file = '/data2/rluo4/EpiTrans/DataCollection/HsRMdatasets.txt', quote = F, row.names = F, col.names = F)
# + GSE169589;GSE215095;GSE210867;GSE228267;GSE240879
datasetlists <- read.table('/data2/rluo4/EpiTrans/DataCollection/HsRMdatasets.txt')
colnames(datasetlists) <- 'GSE'
index <- datasetlists$GSE %in% RMD.HighRes$GSE
datasetlists$PMID[index] <- RMD.HighRes$Pubmed.ID[match(datasetlists$GSE[index], RMD.HighRes$GSE)]
datasetlists$PaperName[index] <- RMD.HighRes$Paper.Name[match(datasetlists$GSE[index], RMD.HighRes$GSE)]
datasetlists$Journal[index] <- RMD.HighRes$Journal[match(datasetlists$GSE[index], RMD.HighRes$GSE)]
datasetlists$Cell.Line.Tissue[index] <- RMD.HighRes$Cell.Line[match(datasetlists$GSE[index], RMD.HighRes$GSE)]
datasetlists$Modification[index] <- RMD.HighRes$Modification[match(datasetlists$GSE[index], RMD.HighRes$GSE)]
datasetlists$Treatment.Condition[index] <- RMD.HighRes$Treatment[match(datasetlists$GSE[index], RMD.HighRes$GSE)]
datasetlists$Technique[index] <- RMD.HighRes$Technique[match(datasetlists$GSE[index], RMD.HighRes$GSE)]

index <- datasetlists$GSE %in% RMD.LowRes$GSE
datasetlists$PMID[index] <- RMD.LowRes$Pubmed.ID[match(datasetlists$GSE[index], RMD.LowRes$GSE)]
datasetlists$PaperName[index] <- RMD.LowRes$Paper.Name[match(datasetlists$GSE[index], RMD.LowRes$GSE)]
datasetlists$Journal[index] <- RMD.LowRes$Journal[match(datasetlists$GSE[index], RMD.LowRes$GSE)]
datasetlists$Cell.Line.Tissue[index] <- RMD.LowRes$`Cell_line/tissues`[match(datasetlists$GSE[index], RMD.LowRes$GSE)]
datasetlists$Modification[index] <- RMD.LowRes$modeifcatioin[match(datasetlists$GSE[index], RMD.LowRes$GSE)]
datasetlists$Treatment.Condition[index] <- RMD.LowRes$`Condition/Treatments`[match(datasetlists$GSE[index], RMD.LowRes$GSE)]
datasetlists$Technique[index] <- 'MeRIP-seq'

m6A.datasets <- m6A.datasets1[grep('Homo', m6A.datasets1$Species), ]
index <- datasetlists$GSE %in% m6A.datasets$GSE
# View(datasetlists[index,] )
datasetlists$PMID[index] <- m6A.datasets$Pubmed_ID[match(datasetlists$GSE[index], m6A.datasets$GSE)]
datasetlists$PaperName[index] <- m6A.datasets$Paper_Name[match(datasetlists$GSE[index], m6A.datasets$GSE)]
datasetlists$Journal[index] <- m6A.datasets$Journal[match(datasetlists$GSE[index], m6A.datasets$GSE)]
datasetlists$Cell.Line.Tissue[index] <- m6A.datasets$Cell_line[match(datasetlists$GSE[index], m6A.datasets$GSE)]
datasetlists$Modification[index] <- m6A.datasets$modeifcatioin[match(datasetlists$GSE[index], m6A.datasets$GSE)]
datasetlists$Treatment.Condition[index] <- m6A.datasets$Treatment[match(datasetlists$GSE[index], m6A.datasets$GSE)]
# datasetlists$Technique[index] <- m6A.datasets$Technique[match(datasetlists$GSE[index], m6A.datasets$GSE)]
# datasetlists$Sample[index] <- m6A.datasets$Sample[match(datasetlists$GSE[index], m6A.datasets$GSE)]
# datasetlists$SRP[index] <- m6A.datasets$Treatment[match(datasetlists$GSE[index], m6A.datasets$GSE)]

m6A.datasets <- m6A.datasets3[grep('Homo', m6A.datasets3$Species), ]
index <- datasetlists$GSE %in% m6A.datasets$GSE
# View(datasetlists[index,] )
datasetlists$PMID[index] <- m6A.datasets$Pubmed_ID[match(datasetlists$GSE[index], m6A.datasets$GSE)]
datasetlists$PaperName[index] <- m6A.datasets$Paper_Name[match(datasetlists$GSE[index], m6A.datasets$GSE)]
datasetlists$Journal[index] <- m6A.datasets$Journal[match(datasetlists$GSE[index], m6A.datasets$GSE)]
datasetlists$Cell.Line.Tissue[index] <- m6A.datasets$Cell_line[match(datasetlists$GSE[index], m6A.datasets$GSE)]
datasetlists$Modification[index] <- m6A.datasets$modeifcatioin[match(datasetlists$GSE[index], m6A.datasets$GSE)]
datasetlists$Treatment.Condition[index] <- m6A.datasets$Treatment[match(datasetlists$GSE[index], m6A.datasets$GSE)]
# datasetlists$Technique[index] <- m6A.datasets$Technique[match(datasetlists$GSE[index], m6A.datasets$GSE)]
# datasetlists$Sample[index] <- m6A.datasets$Sample[match(datasetlists$GSE[index], m6A.datasets$GSE)]
# datasetlists$SRP[index] <- m6A.datasets$Treatment[match(datasetlists$GSE[index], m6A.datasets$GSE)]
write.csv(datasetlists, file = '/data2/rluo4/EpiTrans/DataCollection/HsRMdatasets.csv', quote = F, row.names = F, col.names = F)

# 4) RMPs
modSNV <- read.table('/data2/rluo4/EpiTrans/DataCollection/RMBase/hg38.modSNV.tar.gz', sep = '\t', fill = T)
modSNV <- modSNV[-1,]
Header <- read.table('/data2/rluo4/EpiTrans/DataCollection/RMBase/RMBaseSNVHeader.txt', fill = T, sep = ':')
Header <- str_split(Header$V2, " -> ", simplify = T)[,1]
Header <- gsub(' ', '', Header)

colnames(modSNV) <- Header
RMP <-  str_split(modSNV$geneName, ";", simplify = T)[,1] 
RMP <-  str_split(RMP, ",", simplify = T)[,1] 
unique(RMP)

#######intergrate the RMPs#####
RMP_regulator <- read.csv('/data2/rluo4/EpiTrans/DataCollection/RMP_regulator.csv') # from RMBaseV3.0 http://bioinformaticsscience.cn/rmbase/rbpCancer.php
RM.reg <- read.xlsx('/data2/rluo4/EpiTrans/DataCollection/RM-regulators.xlsx')#,1)
RM.reg$RNA.modifications <- zoo::na.locf(RM.reg$RNA.modifications)#, fromLast = TRUE)
colnames(RM.reg)
RMP_vector <- tidyr::gather(RM.reg, key = "Key", value = "Value", 3:5)
RMP_vector <- unlist(strsplit(na.omit(RMP_vector$Value), ', '))
RMP_vector <- unique(RMP_vector)

# m6A <- RM.reg[RM.reg$RNA.modifications=='m6A', 3:5]
# Convert the dataframe to a long vector using gather
# m6A_vector <- tidyr::gather(RM.reg[RM.reg$RNA.modifications=='m6A', ], key = "Key", value = "Value", -RNA.modifications)
m6A_vector <- tidyr::gather(RM.reg[RM.reg$RNA.modifications=='m6A', ], key = "Key", value = "Value", 3:5)
m6A_vector <- unlist(strsplit(na.omit(m6A_vector$Value), ', '))
m6A_vector <- unique(m6A_vector)
# Print the long vector
# print("\nLong Vector:")
print(m6A_vector) # 32 unique m6A enzymes + 2

colnames(m6A)
m6A.gene <- unique(str_split(m6A$geneName, ',', simplify = T)[,1])
m6A.RMP <- unique(str_split(m6A$writerID, ':', simplify = T)[,1])
# DEG <- unique(Patt$Symbol[Patt$Organ=='Liver'])
# table(DEG %in% m6A.gene)
m6A.cancer <- read.table('/data2/rluo4/EpiTrans/DataCollection/m6A-regulator-cancer.txt')
m6A.23 <- read.table('/data2/rluo4/EpiTrans/DataCollection/m6A-23.txt')
m6A.cancer <- unique(m6A.cancer$V1)
m6A.23 <- unique(m6A.23$V1) # cancer-related: 24 unique m6A enzymes
setdiff(m6A.cancer, m6A.23)
setdiff(m6A.23, m6A.cancer)

# table(rownames(Liver_Epi) %in% m6A.23)
# table(DEG %in% m6A.23)
# setdiff(m6A.23, unique(Patt$Symbol))
setdiff(m6A_vector,m6A.23)
#  [1] "VIRMA"  "HAKAI"  "ZCCHC4" "METTL5" "HNRNPG" "EIF3A"  "EIF3B" 
# [8] "EIF3C"  "EIF3H"  "PRRC2A" "SND1"   "SAFB"   "SAFB2" 
setdiff(m6A.23, m6A_vector)
# [1] "CBLL1"    "ELAVL1"   "KIAA1429"
# alias:  HAKAI=CBLL1;HUR=ELAVL1;VIRMA=KIAA1429;
m6A.reg <- toupper(unique(c(m6A.23[! m6A.23 %in% c('CBLL1','KIAA1429')], # 2 duplicates: HAKAI=CBLL1 and VIRMA=KIAA1429
                            m6A_vector))) # 35 m6A regulatory enzymes in the literatures.

RMP <- toupper(unique(c(m6A.23[! m6A.23 %in% c('CBLL1','KIAA1429')], 
                        RMP_vector)))
# X29 and its human homologue H29K (Nudt16) are nuclear nucleoside diphosphatase proteins localized within foci in the nucleolus and nucleoplasm. These proteins can remove m(7)G and m(227)G caps from RNAs, rendering them substrates for 5'-3' exonucleases for degradation in vivo.
RMP <- gsub('H29K','NUDT16',RMP)
# Component of the cap-binding complex (CBC), which binds co-transcriptionally to the 5' cap of pre-mRNAs and is involved in various processes such as pre-mRNA splicing, translation regulation, nonsense-mediated mRNA decay, RNA-mediated gene silencing (RNAi) by microRNAs (miRNAs) and mRNA export. The CBC complex is involved in mRNA export from the nucleus via its interaction with ALYREF/THOC4/ALY, leading to the recruitment of the mRNA export machinery to the 5' end of mRNA and to mRNA export in a 5' to 3' direction through the nuclear pore. The CBC complex is also involved in mediating U snRNA and intronless mRNAs export from the nucleus. The CBC complex is essential for a pioneer round of mRNA translation, before steady state translation when the CBC complex is replaced by cytoplasmic cap-binding protein eIF4E. The pioneer round of mRNA translation mediated by the CBC complex plays a central role in nonsense-mediated mRNA decay (NMD), NMD only taking place in mRNAs bound to the CBC complex, but not on eIF4E-bound mRNAs. The CBC complex enhances NMD in mRNAs containing at least one exon-junction complex (EJC) via its interaction with UPF1, promoting the interaction between UPF1 and UPF2. The CBC complex is also involved in 'failsafe' NMD, which is independent of the EJC complex, while it does not participate in Staufen-mediated mRNA decay (SMD). During cell proliferation, the CBC complex is also involved in microRNAs (miRNAs) biogenesis via its interaction with SRRT/ARS2, thereby being required for miRNA-mediated RNA interference. The CBC complex also acts as a negative regulator of PARN, thereby acting as an inhibitor of mRNA deadenylation. In the CBC complex, NCBP2/CBP20 recognizes and binds capped RNAs (m7GpppG-capped RNA) but requires NCBP1/CBP80 to stabilize the movement of its N-terminal loop and lock the CBC into a high affinity cap-binding state with the cap structure. The conventional cap-binding complex with NCBP2 binds both small nuclear RNA (snRNA) and messenger (mRNA) and is involved in their export from the nucleus 
RMP <- gsub('CBC','NCBP2',RMP)
RMP <- unique(c(RMP,"NCBP1")) 
# H/ACA snoRNPs are the complexes of box H/ACA snoRNA and four proteins, dyskerin (NAP57/DKC1), NOP10, NHP2 and GAR1 to target rRNA and snRNA for pseudouridylation. Only DKC1 is catalytically active and catalyses U-to-Ψ isomerization, and belongs to the TruB family. NOP10 and GAR1 enhance the catalytic activity of DKC1 [176]. In addition, GAR1 contributes to substrate recruitment and product release [177]. NHP2 and Box H/ACA help target RNAs to enter the catalytic centre of H/ACA RNPs [178]. scaRNA is small Cajal body-specific RNA, a specialised snoRNA in Cajal bodies (CBs), which can also guide Ψ modification in conjunction with the above four proteins [179].
# 2′-O-Me writers include two types of enzymes: independent proteases and the C/D-box snoRNPs complex (SNORNP). The independent proteases in humans are CCDC76 (hTrmt13), FTSJ1, TARBP1, MRM1, MRM2 (FTSJ2), MRM3 (RNMTL1), FTSJ3, HENMT1 and CMTR1/2. But, no reader or eraser for 2′-O-Me has yet been reported.
RMP <- unique(c(RMP, 'NAP57', 'DKC1', 'NOP10', 'NHP2', 'GAR1', #H/ACA snoRNPs
         'CCDC76', 'FTSJ1', 'TARBP1', 'MRM1', 'MRM2', 'MRM3', 'FTSJ3', 'HENMT1', 'CMTR1', 'CMTR2'# independent proteases
         ))
RMP <- RMP[! RMP %in% c('H/ACA SNORNPS', 'C/D-BOX SNORNPS')]
RMP <- data.frame(RMP = RMP)
unique(RMP$RMP)# 121 RMPs in the literatures that I reviewed

###########change the gene symbols of RMP##########
# setdiff(m6A.reg, unique(Patt$Symbol))
# "RBM15"  "HAKAI" "ZCCHC4" "HNRNPG"
setdiff(m6A.reg, unique(gene_info$symbols))
# [1] "HAKAI"    "HNRNPG"
intersect(RMP$RMP, RMP_regulator$Rmbase_regulator) #73
setdiff(RMP$RMP, RMP_regulator$Rmbase_regulator) #48
setdiff(RMP_regulator$Rmbase_regulator, RMP$RMP) #72
# [1] "ADAR"     "ADARB1"   "ADAT1"    "ADAT2"   
# [5] "ADAT3"    "ALKBH8"   "BCDIN3D"  "BUD23"   
# [9] "CBLL1"    "CDK5RAP1" "CDKAL1"   "CTU1"...
setdiff(RMP$RMP, unique(gene_info$symbols))
alias <- sort(setdiff(RMP$RMP, unique(gene_info$geneSymbol))) #18
# View(gene_info[grepl(paste0(strings, collapse = '|'), gene_info$geneAlias),  c('symbols','geneAlias')])
alias_df <- data.frame(symbols = alias, geneSymbol = c('ADAR','ADARB1','ADARB2','TRMT13','TRDMT1','FMR1','CBLL1','RBMX','TRMT13',
                                                    'MARS1','DKC1','RRP8','NOP2','DDX46','RET1', 'RET2', 'HSD17B10','BUD23'))# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6232591/
setdiff(RMP_regulator$Rmbase_regulator, unique(gene_info$symbols))
# [1] "HuR"
RMP_regulator$Rmbase_regulator <- gsub('HuR','ELAVL1',RMP_regulator$Rmbase_regulator)
# alias: ADAR=ADAR1; ADARB1=AVAR2; CBLL1=HAKAI; 
# alias: KIAA1429=VIRMA;HAKAI=CBLL1;HNRNPG=RBMX;
RMP$RMP[match(alias_df$symbols,RMP$RMP)]<- alias_df$geneSymbol[alias_df$symbols %in% RMP$RMP]
RMP <- union(RMP_regulator$Rmbase_regulator, RMP$RMP)
RMP <- data.frame(RMP_update = RMP)
setdiff(RMP$RMP_update, unique(gene_info$geneSymbol))
write.csv(RMP, file = '/data2/rluo4/EpiTrans/DataCollection/RMP.csv', quote = F, row.names = F)

cell_RM <- read.table('/data2/rluo4/EpiTrans/DataCollection/cell-RM.txt', sep = '\t', header = T)
disease_RM <- read.table('/data2/rluo4/EpiTrans/DataCollection/disease-RM.txt',sep = '\t', header = T)

RMP_classified <- read_excel('/data2/rluo4/EpiTrans/DataCollection/RMP.xlsx', sheet = 2)
# colnames(RMP_classified) <- c('unclassified', 'modifications', 'type', '4','5','6')
table(RMP_classified$RMP %in% RMP$RMP_update)
intersect(RMP_classified$RMP, alias_df$symbols)
setdiff(RMP_classified$RMP,RMP$RMP_update)
setdiff(RMP$RMP_update, RMP_classified$RMP)
# [1] "ADAT1"    "ADAT2"    "ADAT3"    "ALKBH8"  
# [5] "BCDIN3D"  "CDK5RAP1" "CDKAL1"   "CTU1"    
# [9] "CTU2"     "DCPS"     "DGCR8"    "DIMT1"   
# [13] "DUS1L"    "DUS2"     "DUS3L"    "DUS4L"   
# [17] "ELP1"     "ELP2"     "ELP3"     "ELP4"    
# [21] "EMG1"     "FBL"      "FTSJ3"    "FXR1"    
# [25] "GAR1"     "GTPBP3"   "HENMT1"   "ELAVL1"  
# [29] "LCMT2"    "MEPCE"    "METTL2A"  "METTL2B" 
# [33] "METTL6"   "METTL8"   "MOCS2"    "MOCS3"   
# [37] "MTERF4"   "MTO1"     "NHP2"     "NOP10"   
# [41] "NOP56"    "NOP58"    "NXF1"     "PUS7L"   
# [45] "PUSL1"    "QTRT1"    "QTRT2"    "RNGTT"   
# [49] "RPUSD1"   "SNU13"    "TET3"     "TFB1M"   
# [53] "THG1L"    "TRIT1"    "TRMO"     "TRMT1"   
# [57] "TRMT10A"  "TRMT12"   "TRMT1L"   "TRMT2B"  
# [61] "TRMT44"   "TRMT5"    "TRMT9B"   "TRMU"    
# [65] "TYW1"     "TYW3"     "TYW5"     "URM1"    
# [69] "DCP2"     "NCBP2"    "MARS1"   
modifications = c(rep('A-to-I', 3), 'cm5U', 'Capping', 'ms2i6A', 'ms2t6A', rep('mcm5U',2), 
                  'm7GMP','m6A','Other', 'Other', 'Other', 'Other', 'Other', rep('mcm5U',4),
                  'Psi(Pseudouridine)', rep("2'-O-Me",2), 'm6A', 'Psi(Pseudouridine)', 'Other', "2'-O-Me", 'm6A/m5C',
                  'Other', 'Capping', rep('Other', 4), rep('Other', 5), 'Psi(Pseudouridine)',
                  rep('Other', 2), 'm6A', rep('Psi(Pseudouridine)',2), rep('Other', 3),
                  'Psi(Pseudouridine)', rep('Other', 20),"2'-O-Me", 'Other'
                  )
types = c(rep('Writer', 3), 'Eraser', 'Other', 'Writer', 'Writer',rep('Reader', 2),
          'Other', 'Writer', 'Other', 'Other', 'Other', 'Other', 'Other', rep('Writer', 4),
          'Writer', rep('Writer', 2), 'Reader', 'Writer', 'Other', 'Writer', 'Reader',
          'Other', 'Other', rep('Writer', 4), rep('Other', 5), 'Writer', 
          rep('Other', 2), 'Other', rep('Writer', 2), rep('Other', 3),
          'Writer', rep('Other', 20), 'Writer', 'Other'
          )
RMP_unclassified <- data.frame( "RNA_modifications" = modifications, 
 "RMP" = setdiff(RMP$RMP_update, RMP_classified$RMP), "RMP_type" = types)
colnames(RMP_classified) <- colnames(RMP_unclassified)
RMP_update <- rbind(RMP_classified,RMP_unclassified)
table(RMP_update$RMP %in% RMP$RMP_update)
RMP_update<- RMP_update[order(RMP_update$RNA_modifications),]
length(unique(RMP_update$RMP)) #179
write.xlsx(RMP_update, '/data2/rluo4/EpiTrans/DataCollection/RMP_update.xlsx')
write.csv(RMP_update, '/data2/rluo4/EpiTrans/DataCollection/RMP_update.csv')
# write.xlsx(x = PID[[i]]@result[PID[[i]]@result$pvalue<cutoff,-(1:2)], 
#           file = paste0(organ,'_MP_all',".xlsx"),#paste0(organ,'_MP_',i,".xlsx"),
#           sheetName = sheet_name[[i]],
#           append = TRUE)

############nohup above##################

###############################################
# 4.3 Guitar visualization for each condition #
###############################################
# Load GTF
# gtf_v36 <- rtracklayer::import.gff(con = '/data2/rluo4/hg38/gencode.v36.annotation.gtf.gz', format = "gtf")
# gtf_v19 <- rtracklayer::import.gff(con = '/data2/rluo4/hg38/gencode.v19.annotation.gtf.gz', format = "gtf")
load( paste0('/data2/rluo4/EpiTrans/RMDatasets/txdb_Guitar.RData') )
txdb_hg19 <- makeTxDbFromGFF(file = '/data2/rluo4/hg38/gencode.v19.annotation.gtf.gz',
                             format="gtf",  dataSource="Ensembl", organism="Homo sapiens")
txdb_hg38 <- makeTxDbFromGFF(file = '/data2/rluo4/hg38/gencode.v36.annotation.gtf.gz',
                             format="gtf",  dataSource="Ensembl", organism="Homo sapiens")
# save(gtf_v19, gtf_v36, file = paste0('/data2/rluo4/EpiTrans/RMDatasets/txdb_Guitar.RData') )
options(stringsAsFactors = F)
# load package
library(Guitar)
package.version("Guitar")
# [1] "2.10.0"
getwd()
setwd(data_directory)
# i='GSE122948'
# setwd(i)
summit_bed_files <- data.frame()
for (n in 1:nrow(RMPdeg_datasets)) {
  # for (n in bw.rank) {
  print(n)
  print(RMPdeg_datasets[n, 'GSE'])
  i <- RMPdeg_datasets[n, 'GSE']
  files <- list.files(file.path(data_directory, i, 'metaPlotR'))
  # files <- files[grep('summits.bed', files)]
  files <- files[grep('summits.bed', files)]
  if(length(files)==0){
    print('no summit.bed !')
    next;
  }
  summit_bed <- data.frame(rank = n, GSE = i, summit.bed = files)
  summit_bed_files <<- rbind(summit_bed_files, summit_bed)
}
RMP_pattern <- paste(c(RMP_update$RMP,'WT'), collapse = "|")
# Extract the matching words
matches <- str_extract(summit_bed_files$summit.bed, RMP_pattern)
# Add the matches as a new column in the dataframe
summit_bed_files$RMP <- matches
summit_bed_files$RMP[grepl('TGIRT',summit_bed_files$summit.bed)] <- 'TRMT61A'
summit_bed_files$RMP[grepl('KIAA1429',summit_bed_files$summit.bed)] <- 'KIAA1429' #VIRMA=KIAA1429;
summit_bed_files$RMP[grepl('HAKAI',summit_bed_files$summit.bed)] <- 'HAKAI' #HAKAI=CBLL1
summit_bed_files$RMP[grepl('SETD2',summit_bed_files$summit.bed)] <- 'SETD2' # not RMP
summit_bed_files$RMP[grepl('DNMT2',summit_bed_files$summit.bed)] <- 'DNMT2' # not RMP
summit_bed_files$RMP[grepl('ontrol|GFP|m6A_rep|METTL14WT',summit_bed_files$summit.bed)] <- 'WT'
summit_bed_files[is.na(summit_bed_files$RMP),]
# summit_bed_files <- summit_bed_files[!grepl('GSM2324292|GSM2324294|GSM2324302|GSM2324298|GSM2324310|GSM2324306',summit_bed_files$summit.bed),]
summit_bed_files <- summit_bed_files[!grepl('GSM2324292|GSM2324294',summit_bed_files$summit.bed),] # remove the files of NOMO-1 cell line
# GSM2324291	NOMO-1 PBS input
# GSM2324292	NOMO-1 PBS m6A-seq
# GSM2324293	NOMO-1 R-2HG input
# GSM2324294	NOMO-1 R-2HG m6A-seq
# summit_bed_files$RMP[is.na(summit_bed_files$RMP)] <- c('FTO','WT','FTO','WT')
summit_bed_files$RMP[summit_bed_files$rank==33] <- c('WT','WT','FTO','FTO', 'WT','WT','FTO','FTO')
flatten_condition_IP_list <- function(condition_list) {
  unlist(unique(unlist(condition_list)))
}
conditions <- flatten_condition_IP_list(condition_IP_list)
conditions <- conditions[conditions!='WT']
conditions <- c(conditions, 'WT')
conditions_pattern <- paste(conditions, collapse = "|")
matches <- str_extract(summit_bed_files$summit.bed, conditions_pattern)
# Add the matches as a new column in the dataframe
summit_bed_files$condition <- matches
summit_bed_files$condition[grepl('TGIRT',summit_bed_files$summit.bed)] <- 'TRMT61AOE'
summit_bed_files$condition[grepl('siNSUN2',summit_bed_files$summit.bed)] <- 'NSUN2KD'
summit_bed_files$condition[grepl('siMETTL3|shMETTL3',summit_bed_files$summit.bed)] <- 'METTL3KD'
summit_bed_files$condition[grepl('siMETTL14|shMETTL14',summit_bed_files$summit.bed)] <- 'METTL14KD'
summit_bed_files$condition[grepl('siMETTL3andsiMETTL14',summit_bed_files$summit.bed)] <- 'METTL3_14KD'

summit_bed_files$condition[grepl('siWTAP|shWTAP',summit_bed_files$summit.bed)] <- 'WTAPKD'
summit_bed_files$condition[grepl('FTOIB',summit_bed_files$summit.bed)] <- 'FTO.2IB' 
summit_bed_files$condition[grepl('siKIAA1429',summit_bed_files$summit.bed)] <- 'KIAA1429KD'#VIRMA=KIAA1429;
summit_bed_files$condition[grepl('SETD2_2KD',summit_bed_files$summit.bed)] <- 'SETD2_2KD' # not RMP condition
summit_bed_files$condition[grepl('DNMT2',summit_bed_files$summit.bed)] <- 'DNMT2PD' # not RMP condition
summit_bed_files$condition[grepl('ontrol|GFP|m6A_rep|METTL14WT',summit_bed_files$summit.bed)] <- 'WT'
summit_bed_files$condition[grepl('YTHDF3KD_YTHDF',summit_bed_files$summit.bed)] <- c('YTHDF3KD_YTHDF1PD', 'YTHDF3KD_YTHDF1PD', 'YTHDF3KD_YTHDF2PD', 'YTHDF3KD_YTHDF2PD')

summit_bed_files[is.na(summit_bed_files$condition),]
summit_bed_files$condition[summit_bed_files$rank==33] <- c('WT','WT','FTOKD','FTOKD','WT','WT','FTOOE','FTOOE')
# summit_bed_files$condition[grepl('METTL3KO_N_',summit_bed_files$summit.bed)] <- 'METTL3KO_N' 
# summit_bed_files$condition[grepl('METTL3KO_NL_',summit_bed_files$summit.bed)] <- 'METTL3KO_NL' 
# summit_bed_files$condition[grepl('METTL3KO_TL_',summit_bed_files$summit.bed)] <- 'METTL3KO_TL' 
summit_bed_files <- summit_bed_files[!grepl('GSE132306',summit_bed_files$GSE),] # rm GSE132306(n = 48) whose result is the same with GSE120024(n = 47)
# summit_bed_files <- summit_bed_files[!grepl('GSE130172',summit_bed_files$GSE),]
summit_bed_files <- summit_bed_files[!grepl('GSE171497',summit_bed_files$GSE),] # rm GSE171497(n = 107) whose result is not good for now
setdiff(names(target_results), unique(summit_bed_files$GSE))#[1] "GSE125046" "GSE109183" "GSE122961" "GSE122948" "GSE128699"
# with RIP/miCLIP without Input, so no summit.bed
summit_bed_files$path = file.path(data_directory, summit_bed_files$GSE, 'metaPlotR', summit_bed_files$summit.bed)
write.table(summit_bed_files, quote = F, row.names = F, sep = '\t', col.names = F,
            file = paste0('/data2/rluo4/EpiTrans/RMDatasets/summit_bed_files.txt'))

# http://localhost:8888/edit/Shell/UTH_RPMfunc/RPMdb-bw.sh
# original <- c('GSE97419', 'GSE55572', 'GSE87190', 'GSE103497', 'GSE141994', 'GSE145924')
# 'GSE97419', 'GSE55572', 'GSE87190', 'GSE76414',
# 'GSE103497', 'GSE141994', 'GSE145924', 'GSE124509', 
# 'GSE210867',  'GSE207643', 'GSE171497', 'GSE198643'
paste0(original, collapse =  "' | '")
# 'GSE97419' | 'GSE55572' | 'GSE87190' | 'GSE103497' | 'GSE141994'| 'GSE145924' | 'GSE207643' | 'GSE198643'
summit_bed_files$newpath <- file.path('/data2/rluo4/RPMfunc/GEO_metaPlotR',summit_bed_files$GSE,
                                      summit_bed_files$summit.bed)
# summit_bed_guitar <- vector("list", length(target_results))
# names(summit_bed_guitar) <- names(target_results)
summit_bed_guitar = target_results
for (n in 1:nrow(RMPdeg_datasets)) {# 1, 7, 13, 19, 25, 28, 43, 45, 48, c(95, 108)){#
  # for (n in bw.rank) {112?
  print(n)
  print(RMPdeg_datasets[n, 'GSE'])
  i <- RMPdeg_datasets[n, 'GSE']
  # if( genome_version$skip[genome_version$GSE==i] =='NoData' ){
  #   print(paste0(i, ": no data finally, skip  !"))
  #   cmd = paste0('rm -r ',  file.path(plot_dir, i))
  #   system(cmd)
  #   next;
  # }
  target_res <- target_results[[i]]
  if( n %in% c(1,7, 13, 19, 25, 28,  43, 45,  58, 65, 68, 90, 99, 100, 102, 108, 116)){
    print(paste0(i, ": error in summit.bed !"))
    next;
  }
  # if( n %in% c(1,7, 13, 19, 25, 28, 36, 43, 45, 48, 58, 65, 68, 90, 99, 100, 102, 107, 108, 116)){
  # Guitar for below still empty
  # [1] "start figure plotting for mrna ..."
  # Error in density.default(siteID, adjust = adjust, from = 0, to = 1, n = 256,  :
  #                            'x' contains missing values
  #                          In addition: There were 13 warnings (use warnings() to see them)
  # n = 1, GSE102113/Hela_NAT10KO
  # n = 7, GSE93749/Hela_ALYREFPD 
  # n = 13, GSE121942/HepG2.2_SETD2KD + HepG2_SETD2KD 
  # n = 19, GSE102493/Hela_METTL14MU + HeLa_ZC3H13KD + HeLa_VIRMAKD
  # n = 25, GSE120229/HEK293T_PCIF1KD
  # n = 28, GSE55572/HEK293_WTAPKD + HEK293_METTL3KD + A549_METTL14KD 
  # n = 43, GSE90684/HepG2_METTL14KD
  # n = 45, GSE103497/NB4_FTO.2IB
  # n = 58, GSE128575/THP1_ALKBH5KD
  # n = 65, GSE144984/NOMO1_ALKBH5KD
  # n = 68, GSE129469/SU-DHL-8_WTAPKD
  # n = 90, GSE202815/HEK293_ALKBH1KO
  # n = 99, GSE144620/Hela_ALKBH5MU
  # n = 100, GSE195637/Hela.m1A_ALKBH3KO + Hela.m6A_ALKBH3KO
  # n = 102, GSE226129/HCT15_NSUN2KD
  # n = 108, GSE240674/THP1.N_METTL3KO + THP1.NL_METTL3KO + THP1.TL_METTL3KO
  if( ! i %in% summit_bed_files$GSE ){
    print(paste0(i, ": no summit.bed in /data2/rluo4/EpiTrans/RMDatasets !"))
    next;
  }
  outdir <- file.path(plot_dir, i)
  if(! dir.exists( outdir )) {
    dir.create(outdir)
  }
  resG <- RMdeg_list[[i]]
  resM <- RMdem_list[[i]]
  resP <- RMdep_list[[i]]
  # table(target_res)
  cl <- unique(str_split( names(target_res), '_', simplify = TRUE)[, 1])
  print(cl)
  print(names(target_res))
  # for ( contrast in names(target_res)[-(1:4)] ) {
  for ( contrast in names(target_res) ) {
    # contrast <- names(target_res)[1]
    # contrast <- names(target_res)[2]
    outdir <- file.path(plot_dir, i, contrast)
    if(! dir.exists( outdir )) {
      dir.create(outdir)
    }
    outdir <- file.path(plot_dir, i, contrast, "Guitar/")
    if( file.exists( file.path(outdir, 'guitarPlot_mRNA.png') )){
      print(paste0(i, ": GuitarPlot done !"))
      next;
    }
    if( file.exists( file.path(outdir, 'guitarPlot_ncRNA.png') )){
      print(paste0(i, ": GuitarPlot done !"))
      next;
    }
    if(! dir.exists( outdir )) {dir.create(outdir)}
    contrast_IP_list[[i]][[contrast]]
    condition_IP_list[[i]][[contrast]]
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
    RMP =  remove_first_element(contrast)#unique(str_split(contrast, '_', simplify = TRUE)[, 2])
    RMP <- gsub('VIRMA_1MU', 'VIRMAMU', RMP)
    RMP <- gsub('VIRMA_2MU', 'VIRMAMU', RMP)
    if(contrast %in% c('Hela.YTHDF3KD_YTHDF1PD','Hela.YTHDF3KD_YTHDF2PD')){
      RMP <- gsub('Hela.', '', contrast)
    }
    c <- unique(str_split(contrast, '_', simplify = TRUE)[, 2])
    if( length(cl) >1 ){
      c <- contrast
    }
    fullname <- paste0(contrast, '_vs_WT')
    # if( i %in% c('GSE97419', 'GSE55572', 'GSE87190', 'GSE76414', 'GSE103497', 'GSE141994')){
    format_gse <- 'bw'
    track_folder <- file.path(data_directory, i, "metaPlotR")
    index1 = summit_bed_files$GSE==i
    bedfiles <- summit_bed_files$summit.bed[index1]
    bedfiles
    condition <- RMP#paste0('WT|', RMP)
    index2 = summit_bed_files$condition[index1]==RMP#grepl(RMP, summit_bed_files$condition[index1])
    index3 = summit_bed_files$condition[index1]=='WT'
    bedfiles <- c(bedfiles[index2], bedfiles[index3])
    # bedfiles <- bedfiles[! grepl('annot|MEF', bedfiles)] #mouse embryonic fibroblasts (MEFs)
    summit_bed <- rbind(summit_bed_files[index1,][index2,], summit_bed_files[index1,][index3,])
    summit_bed_guitar[[i]][[contrast]] = summit_bed
    #   }
    # }
    stBedFiles <- list()
    # Loop through each file path and add it to the list
    # for (file in  file.path(track_folder, summit_bed$summit.bed)) {
    for (file in  summit_bed$newpath) {
      stBedFiles <- c(stBedFiles, list(file))
    }
    stBedFiles
    
    if(n == 13 & contrast %in% c('HepG2_METTL3KD', 'HepG2_METTL14KD', 'HepG2_SETD2KD','HepG2_WTAPKD')){
      stBedFiles <- stBedFiles[ ! grepl('WT_2', stBedFiles)]
    }
    if(n == 13 & contrast %in% c('HepG2.2_SETD2KD')){
      stBedFiles <- stBedFiles[ ! grepl('WT_rep', stBedFiles)]
      stBedFiles <- gsub('SETD2', 'SETD2_2', stBedFiles)
    }
    
    if(n == 17 & contrast %in% c('Hela_METTL1KD')){
      stBedFiles <- stBedFiles[ grepl('Hela_', stBedFiles)]
    }
    if(n == 17 & contrast %in% c('HepG2_METTL1KD')){
      stBedFiles <- stBedFiles[ grepl('HepG2_', stBedFiles)]
    }
    
    if(n == 19 & contrast %in% c('HeLa_VIRMAPD', 'HeLa_METTL3MU','HeLa_METTL14MU','HeLa_VIRMA_2MU')){
      stBedFiles <- stBedFiles[ ! grepl('WT_rep3|WT_rep4|WT_rep5|WT_rep6|WT_rep7|WT_rep8|WT_rep11|WT_rep12', stBedFiles) ]
    }
    if(n == 19 & contrast %in% c('HeLa_ZC3H13KD', 'HeLa_HAKAIKD', 'HeLa_VIRMAKD')){
      stBedFiles <- stBedFiles[ ! grepl('WT_rep1_|WT_rep2_|WT_rep7|WT_rep8', stBedFiles) ]
    }
    if(n == 19 & contrast %in% c('HeLa_VIRMA_1MU')){
      stBedFiles <- stBedFiles[ grepl('VIRMAMU_rep1|VIRMAMU_rep2|WT_rep7|WT_rep8', stBedFiles) ]
    }
    
    if(n == 21 & contrast %in% c('HEC-1-A_METTL14MU')){
      stBedFiles <- stBedFiles[ ! grepl('/WT_summits.bed', stBedFiles)]
    }
    if(n == 21 & contrast %in% c('HEC-1-A_METTL3KD')){
      stBedFiles <- stBedFiles[ ! grepl('METTL14WT', stBedFiles)]
    }
    
    if(n == 24 & contrast %in% c('NB4_METTL14KD')){
      stBedFiles <- stBedFiles[ ! grepl('MM6_', stBedFiles)]
    }
    if(n == 24 & contrast %in% c('MM6_METTL14KD')){
      stBedFiles <- stBedFiles[ ! grepl('NB4_', stBedFiles)]
    }
    
    if(n == 28 & contrast %in% c('A549_METTL14KD', 'A549_METTL3KD', 'A549_KIAA1429KD', 'A549_METTL3_14KD')){
      stBedFiles <- stBedFiles[ ! grepl('Exp1si', stBedFiles) ]
    }
    if(n == 28 & contrast %in% c('HEK293_WTAPKD', 'HEK293_METTL3KD')){
      stBedFiles <- stBedFiles[  grepl('Exp1si', stBedFiles) ]
    }
    
    if(n == 33 & contrast %in% c('MA9.3ITD_FTOKD')){
      stBedFiles <- stBedFiles[ ! grepl('GSM2324304|GSM2324306', stBedFiles)]
    }
    if(n == 33 & contrast %in% c('MA9.3RAS_FTOOE')){
      stBedFiles <- stBedFiles[ ! grepl('GSM2324296|GSM2324298', stBedFiles)]
    }
    
    if(n == 40 & contrast %in% c('Hela_NSUN2PD')){
      print(paste0(i, ": no summit.bed for this contrast -- Hela_NSUN2PD vs WT !"))
      next;
    }
    
    if(n == 46 & contrast %in% c('MEL624.1_PCIF1KO')){
      stBedFiles <- stBedFiles[ ! grepl('PCIF1KO.2_', stBedFiles)]
    }
    if(n == 46 & contrast %in% c('MEL624.2_PCIF1KO')){
      stBedFiles <- stBedFiles[ ! grepl('PCIF1KO_', stBedFiles)]
    }
    
    if(n == 63 & contrast %in% c('PC-3_FTOKD')){
      stBedFiles <- stBedFiles[ ! grepl('_rep', stBedFiles)]
    }
    if(n == 63 & contrast %in% c('SK-BR-3_FTOKD')){
      stBedFiles <- stBedFiles[  grepl('_rep', stBedFiles)]
    }
    
    if(n == 84 & contrast %in% c('Hela_YTHDF3PD')){
      stBedFiles <- stBedFiles[ ! grepl('WT_YTHDF', stBedFiles) ]
      file <- c(list("/data2/rluo4/RPMfunc/GEO_metaPlotR/GSE46705/WT_rep1_summits.bed"), 
                list("/data2/rluo4/RPMfunc/GEO_metaPlotR/GSE46705/WT_rep2_summits.bed"),
                list("/data2/rluo4/RPMfunc/GEO_metaPlotR/GSE46705/WT_rep3_summits.bed"),
                list("/data2/rluo4/RPMfunc/GEO_metaPlotR/GSE46705/WT_rep4_summits.bed")
      )
      stBedFiles <<- c(stBedFiles, file)
    }
    
    
    if(n == 95 & contrast %in% c('Hela_YTHDF3PD')){
      stBedFiles <- stBedFiles[ ! grepl('WT_YTHDF', stBedFiles) ]
      file <- c(list("/data2/rluo4/RPMfunc/GEO_metaPlotR/GSE46705/WT_rep1_summits.bed"), 
                list("/data2/rluo4/RPMfunc/GEO_metaPlotR/GSE46705/WT_rep2_summits.bed"),
                list("/data2/rluo4/RPMfunc/GEO_metaPlotR/GSE46705/WT_rep3_summits.bed"),
                list("/data2/rluo4/RPMfunc/GEO_metaPlotR/GSE46705/WT_rep4_summits.bed")
      )
      stBedFiles <<- c(stBedFiles, file)
    }
    if(n == 95 & contrast %in% c('Hela.YTHDF3KD_YTHDF1PD')){
      stBedFiles <- stBedFiles[ ! grepl('WT_YTHDF2', stBedFiles) ]
    }
    if(n == 95 & contrast %in% c('Hela.YTHDF3KD_YTHDF2PD')){
      stBedFiles <- stBedFiles[ ! grepl('WT_YTHDF1', stBedFiles) ]
    }
    
    if(n == 100 & contrast %in% c('Hela.m1A_ALKBH3KO')){
      stBedFiles <- stBedFiles[  grepl('m1A', stBedFiles) ]
    }
    if(n == 100 & contrast %in% c('Hela.m6A_ALKBH3KO')){
      stBedFiles <- stBedFiles[ !  grepl('m1A', stBedFiles) ]
    }
    
    if(n == 112 & contrast %in% c("A549.RSV_NSUN2KO")){
      stBedFiles <- stBedFiles[ ! grepl('VSV', stBedFiles)]
    }
    if(n == 112 & contrast %in% c("A549.VSV_NSUN2KO")){
      stBedFiles <- stBedFiles[ ! grepl('RSV', stBedFiles)]
    }
    
    
    # if(n == 108 & contrast %in% c('THP1.N_METTL3KO')){
    #   stBedFiles <- stBedFiles[  grepl('_N_', stBedFiles) ]
    # }
    # if(n == 108 & contrast %in% c('THP1.NL_METTL3KO')){
    #   stBedFiles <- stBedFiles[  grepl('_NL_', stBedFiles) ]
    # }
    # if(n == 108 & contrast %in% c('THP1.TL_METTL3KO')){
    #   stBedFiles <- stBedFiles[  grepl('_TL_', stBedFiles) ]
    # }
    
    # condition <- condition[1:length(stBedFiles)]
    # conditions_pattern <- paste(unique(summit_bed_files$condition), collapse = "|")
    # condition <- str_extract(stBedFiles, conditions_pattern)
    # if( ! n %in% c(13, 19)){
    condition_index <- match(unlist(stBedFiles), summit_bed_files$newpath)
    condition = summit_bed_files$condition[condition_index] #condition = summit_bed$condition#
    # }
    print(condition)
    print(unlist(stBedFiles))
    ############
    # gtf_file <- '/data2/rluo4/lorihan/hg38/gencode.v36.annotation.gtf.gz'
    txdb <- if (genome_version$HG[n] == 'hg19') txdb_hg19 else txdb_hg38
    # makeTxDbFromGFF(file = "/data2/rluo4/lorihan/hg38/Homo_sapiens.GRCh38.102.gtf.gz",
    #                       format="gtf",  dataSource="Ensembl", organism="Homo sapiens") #loadTaxonomyDb()
    start_time <- Sys.time()
    print(paste("start visualizing GuitarPlot:", contrast, 'from', i, 'at', start_time))
    pM <- GuitarPlot(txTxdb = txdb,
                     stBedFiles = stBedFiles,
                     headOrtail = TRUE,#FALSE,
                     enableCI = FALSE,
                     mapFilterTranscript = TRUE,
                     pltTxType = c("mrna"),
                     stGroupName = condition)
    end_time <- Sys.time()
    print(end_time)
    time_taken <- end_time - start_time
    print(paste("Time taken of GuitarPlot:", i , '--', contrast,  'as below: '))
    print(time_taken)
    pM <- pM + ggtitle(label = "")  +  xlab("Position") + ylab("Peak density") +
      theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_blank(),
            axis.title.x = element_text(color = 'black', size = 14),
            axis.title.y = element_text(color = 'black', size = 14, angle = 90),
            legend.position = 'top',
            legend.text = element_text(color = 'black', size = 14),
            legend.title = element_blank()) #+ theme_bw()
    ggsave(plot = pM, file = file.path(outdir, 'guitarPlot_mRNA.png'), height =  4, width = 5.6, dpi = 500)
    
    pN <- GuitarPlot(txTxdb = txdb,
                     stBedFiles = stBedFiles,
                     headOrtail = TRUE,#FALSE,
                     enableCI = FALSE,
                     mapFilterTranscript = TRUE,
                     pltTxType = c("ncrna"),
                     stGroupName = condition)
    pN <- pN  +  xlab("Position") + ylab("Peak density")  + ggtitle(label = "") +
      theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_blank(),
            axis.title.x = element_text(color = 'black', size = 14),
            axis.title.y = element_text(color = 'black', size = 14, angle = 90),
            legend.position = 'top',
            legend.text = element_text(color = 'black', size = 14),
            legend.title = element_blank()) #+ theme_bw()
    ggsave(plot = pN, file = file.path(outdir, 'guitarPlot_ncRNA.png'), height =  4, width = 5.6, dpi = 500)
    
  }
}


# solutions <- lapply(summit_bed$newpath, function(x){
#   fread(x) # read.table(paste0(data_directory,i,x),header = T)
# })
# 
# # Merge dataframes in countGTF list by 'gene_id'
# df_in <- Reduce(function(x, y) rbind(x, y), solutions)
# 
# p <- GuitarPlot(txTxdb = txdb,
#                 stBedFiles = stBedFiles,
#                 headOrtail = TRUE,#FALSE,
#                 enableCI = FALSE,
#                 mapFilterTranscript = TRUE,
#                 pltTxType = c("mrna"),
#                 stGroupName = c("Pooled_peaks"))
# # setwd('metaPlotR/')
# 


# genomic features imported into named list
# stBedFiles <- list(system.file("extdata", "m6A_mm10_exomePeak_1000peaks_bed12.bed", package="Guitar"))
# With the following script, we may generate the transcriptomic distribution of genomic features to be tested, and the result will be automatically saved into a PDF file under the working directory with prefix “example”. With the GuitarPlot function, the gene annotation can be downloaded from internet automatically with a genome assembly number provided; however, this feature requires working internet and might take a longer time. The toy Guitar coordinates generated internally should never be re-used in other real data analysis.
BiocManager::install('TxDb.Hsapiens.UCSC.hg38.knownGene')
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
# txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
count <- GuitarPlot(txGenomeVer = "hg38",
                    stBedFiles = stBedFiles,
                    miscOutFilePrefix = NA)
# In a more efficent protocol, in order to re-use the gene annotation and Guitar coordinates, you will have to build Guitar Coordiantes from a txdb object in a separate step. The transcriptDb contains the gene annotation information and can be obtained in a number of ways, .e.g, download the complete gene annotation of species from UCSC automatically, which might takes a few minutes. In the following analysis, we load the Txdb object from a toy dataset provided with the Guitar package. Please note that this is only a very small part of the complete hg19 transcriptome, and the Txdb object provided with Guitar package should not be used in real data analysis. With a TxDb object that contains gene annotation information, we in the next build Guitar coordiantes, which is essentially a bridge connects the transcriptomic landmarks and genomic coordinates.
txdb_file <- system.file("extdata", "mm10_toy.sqlite", package="Guitar")
txdb <- loadDb(txdb_file)
guitarTxdb <- makeGuitarTxdb(txdb = txdb, txPrimaryOnly = FALSE)
############
# You may now generate the Guitar plot from the named list of genome-based features.
GuitarPlot(txTxdb = txdb, stBedFiles = stBedFiles, miscOutFilePrefix = i)#"example")
# Alternatively, you may also optionally include the promoter DNA region and tail DNA region on the 5’ and 3’ side of a transcript in the plot with parameter headOrtail =TRUE.
GuitarPlot(txTxdb = txdb, stBedFiles = stBedFiles, headOrtail = TRUE)
GuitarPlot(txTxdb = txdb, stBedFiles = stBedFiles, headOrtail = TRUE, enableCI = FALSE)
# import different data formats into a named list object.
# These genomic features are using hg38 genome assembly
# Build Guitar Coordinates
GuitarPlot(txTxdb = txdb,  stBedFiles = stBedFiles,   headOrtail = TRUE,    enableCI = FALSE,
           mapFilterTranscript = TRUE, pltTxType = c("mrna"),   stGroupName = c("WT","PCIF1KO1","PCIF1KO2"))
############
# 3 Processing of sampling sites information
# We can select parameters for site sampling.
stGRangeLists = vector("list", length(stBedFiles))
sitesPoints <- list()
for (i in seq_len(length(stBedFiles))) {
  stGRangeLists[[i]] <- blocks(import(stBedFiles[[i]]))#, format = 'bed'
}

for (i in seq_len(length(stGRangeLists))) {
  sitesPoints[[i]] <- samplePoints(stGRangeLists[i], stSampleNum = 10,stAmblguity = 5,
                                   pltTxType = c("mrna"),stSampleModle = "Equidistance",mapFilterTranscript = FALSE,guitarTxdb = guitarTxdb)
}

############
# 4 Guitar Coordinates – Transcriptomic Landmarks Projected on Genome
# The guitarTxdb object contains the genome-projected transcriptome coordinates, which can be valuable for evaluating transcriptomic information related applications, such as checking the quality of MeRIP-Seq data. The Guitar coordinates are essentially the genomic projection of standardized transcript-based coordiantes, making a viable bridge beween the landmarks on transcript and genome-based coordinates. It is based on the txdb object input, extracts the transcript information in txdb, selects the transcripts that match the parameters according to the component parameters set by the user, and saves according to the transcript type (tx, mrna, ncrna).
guitarTxdb <- makeGuitarTxdb(txdb = txdb,   txAmblguity = 5,
                             txMrnaComponentProp = c(0.1,0.15,0.6,0.05,0.1),
                             txLncrnaComponentProp = c(0.2,0.6,0.2),
                             pltTxType = c("tx","mrna","ncrna"),
                             txPrimaryOnly = FALSE)
# 5 Check the Overlapping between Different Components
# We can also check the distribution of the Guitar coordinates built.
gcl <- list(guitarTxdb$tx$tx)
p = GuitarPlot(txTxdb = txdb,  stGRangeLists = gcl, #stGroupName = condition,
           stSampleNum = 200,  enableCI = TRUE, pltTxType = c("tx"),  txPrimaryOnly = FALSE
)
GuitarPlot(txTxdb = txdb,
           stGRangeLists = gcl,
           stSampleNum = 200,
           enableCI = TRUE,
           pltTxType = c("tx"),
           txPrimaryOnly = FALSE
)
# Alternatively, we can extract the RNA components, check the distribution of tx components in the transcriptome.
GuitarCoords <- guitarTxdb$tx$txComponentGRange
type <- paste(mcols(GuitarCoords)$componentType,mcols(GuitarCoords)$txType)
key <- unique(type)
landmark <- list(1,2,3,4,5,6,7,8,9,10,11)
names(landmark) <- key
for (i in 1:length(key)) {
  landmark[[i]] <- GuitarCoords[type==key[i]]
}
GuitarPlot(txTxdb = txdb, stGRangeLists = landmark[1:3], pltTxType = c("tx"),
           enableCI = FALSE)
# Check the distribution of mRNA components in the transcriptome
GuitarPlot(txTxdb = txdb,   stGRangeLists = landmark[4:8],
           pltTxType = c("mrna"),
           enableCI = FALSE)
# Check the distribution of lncRNA components in the transcriptome
GuitarPlot(txTxdb = txdb ,
           stGRangeLists = landmark[9:11],
           pltTxType = c("ncrna"),
           enableCI = FALSE)


############
# gtf_file <- '/data2/rluo4/lorihan/hg38/gencode.v36.annotation.gtf.gz'
txdb <- if (genome_version$HG[n] == 'hg19') txdb_hg19 else txdb_hg38
# makeTxDbFromGFF(file = "/data2/rluo4/lorihan/hg38/Homo_sapiens.GRCh38.102.gtf.gz",
#                       format="gtf",  dataSource="Ensembl", organism="Homo sapiens") #loadTaxonomyDb()
# condition_index <- match(bedfiles, summit_bed_files$summit.bed)
condition = summit_bed$condition#summit_bed_files$condition[condition_index]
start_time <- Sys.time()
print(paste("start visualizing GuitarPlot:", contrast, 'from', i, 'at', start_time))
pM <- GuitarPlot(txTxdb = txdb,
                 stBedFiles = stBedFiles,
                 headOrtail = TRUE,#FALSE,
                 enableCI = FALSE,
                 mapFilterTranscript = TRUE,
                 pltTxType = c("mrna"),
                 stGroupName = condition)
end_time <- Sys.time()
print(end_time)
time_taken <- end_time - start_time
print(paste("Time taken of GuitarPlot:", i , '--', contrast,  'as below: '))
print(time_taken)
pM <- pM + ggtitle(label = "")  +  xlab("Position") + ylab("Peak density") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_blank(),
        axis.title.x = element_text(color = 'black', size = 14),
        axis.title.y = element_text(color = 'black', size = 14, angle = 90),
        legend.position = 'top',
        legend.text = element_text(color = 'black', size = 14),
        legend.title = element_blank()) #+ theme_bw()
ggsave(plot = pM, file = file.path(outdir, 'guitarPlot_mRNA.png'), height =  4, width = 5.6, dpi = 500)

start_time <- Sys.time()
print(paste("start visualizing GuitarPlot:", contrast, 'from', i, 'at', start_time))
# p = GuitarPlot(txTxdb = txdb,  stBedFiles = stBedFiles, stGRangeLists = gcl, stGroupName = condition,
#                stSampleNum = 200,  enableCI = TRUE, pltTxType =  c("mrna"),  txPrimaryOnly = FALSE
# )
pN <- GuitarPlot(txTxdb = txdb,
                 stBedFiles = stBedFiles,
                 headOrtail = TRUE,#FALSE,
                 enableCI = FALSE,
                 mapFilterTranscript = TRUE,
                 pltTxType = c("ncrna"),
                 stGroupName = condition)
end_time <- Sys.time()
print(end_time)
time_taken <- end_time - start_time
print(paste("Time taken of GuitarPlot:", i , '--', contrast,  'as below: '))
print(time_taken)
pN <- pN + ggtitle(label = "Distribution on ncRNA")  +  xlab("Position") + ylab("Peak density") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_blank(),
        axis.title.x = element_text(color = 'black', size = 14),
        axis.title.y = element_text(color = 'black', size = 14, angle = 90),
        legend.position = 'top',
        legend.text = element_text(color = 'black', size = 14),
        legend.title = element_blank()) #+ theme_bw()
ggsave(plot = pN, file = file.path(outdir, 'guitarPlot_ncRNA.png'), height =  3.6, width = 7.2, dpi = 500)



# Load local cytoband data
load_cytoband_data <- function(file_path) {
  tryCatch({
    cyto_data <- read.table(file_path, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    colnames(cyto_data) <- c("chrom", "chromStart", "chromEnd", "name", "gieStain")
    return(cyto_data)
  }, error = function(e) {
    message("Error loading cytoband data: ", e$message)
    return(NULL)
  })
}
# Example of loading local cytoband data for hg19
cytoband_hg19 <- load_cytoband_data("/data2/rluo4/hg38/hg19_cytoBand.txt.gz")

# Function to plot ideogram with error handling and using local cytoband data
plot_ideogram <- function(gene, rank, cytoband_data = NULL) {
  tryCatch({
    # Check if the rank is within the valid range of genome_version$HG
    if(rank > length(genome_version$HG) || rank < 1) {
      stop("Invalid rank: out of range.")
    }
    # Get the genome version based on rank
    genome_ver <- genome_version$HG[rank]
    # Ensure genome_ver is not NULL or empty
    if(is.null(genome_ver) || genome_ver == "") {
      stop("Invalid genome version.")
    }
    # Check if cytoband data is provided
    if(is.null(cytoband_data)) {
      # Attempt to get cytoband data from the function (might hit rate limits)
      ideogram_layer <- geom_ideogram(genome = genome_ver, plot.space = 0, highlight.centromere = TRUE)
    } else {
      # Use local cytoband data
      ideogram_layer <- geom_ideogram(genome = genome_ver, plot.space = 0, highlight.centromere = TRUE, cytoband = cytoband_data)
    }
    
    # Print the ideogram data for debugging
    print(ideogram_layer)
    
    # Combine the basic coverage plot with the ideogram
    plot <- basic_coverage +
      geom_gc(bs.fa.seq = BS) +
      geom_gene(gtf.gr = target_annotation) +
      ideogram_layer # Added ideogram layer
    
    # Print the combined plot
    print(plot)
    
    # Return the plot object
    return(plot)
  }, error = function(e) {
    message("Error occurred: ", e$message)
    return(NULL)  # Return NULL in case of error
  })
}

# Example call for plotting
if (grepl('MT-', g)) {
  plotx <- basic_coverage +
    geom_gc(bs.fa.seq = BS) +
    geom_gene(gtf.gr = target_annotation)
} else {
  plotx <- plot_ideogram(gene = g, rank = n, cytoband_data = cytoband_hg19)
}

# Check if plotx is not NULL before attempting to save
if (!is.null(plotx)) {
  ggsave(plot = plotx, file = paste0(outdir, g, '_peak.png'), height = 3.6, width = 5.2, dpi = 500)
} else {
  message("Plot could not be generated for gene: ", g)
}
##############deeptools###############
# BiocManager::install('CAGEfightR')
# BiocManager::install("Gviz", force = TRUE)
BiocManager::install("BiocStyle")#  on base env of server 36
remotes::install_github("ivanek/Gviz") # on base env by:.libPaths(c("/data2/rluo4/bin/miniconda3/lib/R/library", .libPaths()[-6]))
# and conda install Pandoc, then on R terminal of base env:
devtools::install_github("MalteThodberg/CAGEfightR", build_vignettes = TRUE)
# on Rstudio:
devtools::install_git("https://github.com/koustav-pal/HiCBricks")
remotes::install_github("showteeth/ggcoverage")

library(Gviz)
library(CAGEfightR)
library(rtracklayer)
library("ggcoverage")
library("ggpattern")
# 1) Load GTF
# To add gene annotation, the gtf file should contain gene_type and gene_name attributes in column 9; to add transcript annotation, the gtf file should contain transcript_name attribute in column 9.
# gtf_file <- system.file("extdata", "used_hg19.gtf", package = "ggcoverage")
gtf_v36 <- rtracklayer::import.gff(con = '/data2/rluo4/lorihan/hg38/Homo_sapiens.GRCh38.102.gtf.gz', format = "gtf")
gtf_v36 <- rtracklayer::import.gff(con = '/data2/rluo4/lorihan/hg38/hg38.refGene.gtf.gz', format = "gtf")
gtf_v19 <- rtracklayer::import.gff(con = '/data2/rluo4/lorihan/hg38/hg19.refGene.gtf.gz', format = "gtf")
# gtf_gr <- rtracklayer::import.gff(con = gtf_file, format = "gtf")

test.bed <- as.data.frame(gtf_v36[match(target_res$geneSymbol[1],gtf_v36$gene_name),c(1:6)])
hg38.bed <- read.delim("/data2/rluo4/lorihan/hg38/hg38.txt",sep ="",header = FALSE)
rownames(hg38.bed)<-hg38.bed$V4
# .libPaths(c(.libPaths(), '/data2/rluo4/bin/miniconda3/lib/R/library'))
library(GenomicFeatures)
library(RMariaDB)
# BiocManager::install('RMariaDB')
# hg38.refseq.db <- makeTxDbFromUCSC(genome="hg38", table="refGene")
ensembl_ids <- c("ENSG00000000003", "ENSG00000000005", "ENSG00000000419",
                 "ENSG00000000457", "ENSG00000000460", "ENSG00000000938",
                 "ENSG00000000971", "ENSG00000001036", "ENSG00000001084",
                 "ENSG00000001167")
ensembl_ids <- c('ENSG00000168209')
test.bed <- hg38.bed[match(ensembl_ids,rownames(hg38.bed)),c(1:6)]
test.bed$V4 <- rownames(test.bed)
# write.table(test.bed,file ="test.bed",sep ="\t", quote = FALSE, col.names = FALSE,row.names=FALSE)
# Write the Ensembl IDs to a file named '10.txt'
# writeLines(ensembl_ids, "10.txt")
# Subset the GTF file for the FBF1 gene
target_annotation <- subset(gtf_v36, gene_name == "DDIT4")#subset(gtf_v19, gene_id == "FBF1")

# 2) Load track files and metadata

# pdf(file=paste0(i, "_guitarPlot_mRNA.pdf"), width=8,height=6)
# print(p)
# dev.off()
# p <- GuitarPlot(txTxdb = txdb,
#                 stBedFiles = stBedFiles,
#                 headOrtail = TRUE,#FALSE,
#                 enableCI = FALSE,
#                 mapFilterTranscript = TRUE,
#                 pltTxType = c("ncrna"),
#                 stGroupName = c('WT', 'WT',"TRMT61AOE", "TRMT61AOE" ))
# p <- p + ggtitle(label = "Distribution on ncRNA") +
#   theme(plot.title = element_text(hjust = 0.5)) + theme_bw()
# png(filename=paste0(i, "_guitarPlot_ncRNA.png"),width = 1000,height = 800, res = 180)
# print(p)
# dev.off()

# Function to filter data frames in the list
# filter_significant <- function(df) {
#   # df[df$label == "sig", ]
#   df[ ! df$logFC %in% c(-Inf, Inf), ]
# }
# # Apply the function to each data frame in the list while maintaining the structure
# new_list <- lapply(RMdeg_list, function(sublist) {
#   lapply(sublist, filter_significant)
# })
# remove_chr_prefix <- function(BS) {
#   # Get the current sequence names
#   seq_names <- seqnames(BS)
#   new_seq_names <- sub("^chr", "", seq_names)
#   seqnames(BS) <- new_seq_names
#
#   seq_names <- names(BS@user_seqnames)
#   new_seq_names <- sub("^chr", "", seq_names)
#   names(BS@user_seqnames) <- new_seq_names
#
#   seq_names <- BS@seqinfo@seqnames
#   new_seq_names <- sub("^chr", "", seq_names)
#   BS@seqinfo@seqnames <- new_seq_names
#   return(BS)
# }
# library(GenomicRanges)
# library(diffloop)
# regA <- GRanges(c('1'),IRanges(c(36200000),c(36300000)))
# regA <- addchr(regA)
# regA
# regA <- rmchr(regA)
# regA
# if( ! i %in% c('GSE97419', 'GSE125046')){
#   gtf <- rmchr(gtf)
#   # gtf@seqnames <- gsub('chr', '', gtf@seqnames)
#   # gtf$remap_original_location <- gsub('chr', '', gtf$remap_original_location)
#   BS <- remove_chr_prefix(BS)
# }
target_tss <- start(target_annotation[target_annotation$type == "start_codon"])
# Define the promoter region based on the TSS
promoter_start <- target_tss - 1000  # 1000 bp upstream of the TSS
promoter_end <- target_tss - 1       # 1 bp before the TSS
# Extract the sequence of the promoter region from the genome # Assuming you have the genome sequence loaded into the variable 'genome_sequence'
# target_promoter_sequence <- genome_sequence[Promoter(start = promoter_start, end = promoter_end)]
promoter_df <- data.frame(seqnames = unique(seqnames(target_annotation)),
                          start = unique(promoter_start), end = unique(promoter_end), strand = unique(strand(target_annotation)),
                          gene_name = i, label = "Promoter")
print(promoter_df)
mark_region <- promoter_df[, c('start', 'end', 'label')]
ggcoverage(
  data = track_df,
  plot.type = "joint",
  facet.key = "Group",
  group.key = "Type",
  # mark.region = mark_region,
  range.position = "out"
)
ggcoverage(
  data = track_df,
  plot.type = "facet",
  facet.key = "Type",
  group.key = "Group",
  # mark.region = mark_region,
  range.position = "out"
)
basic_coverage <- ggcoverage(data = track_df,
                             range.position = "in",
                             group.key = "Group",
                             facet.key = "Type",
                             plot.type = "joint",
                             facet.y.scale = "fixed",
                             #mark.region = mark_region,
                             show.mark.label = TRUE)
plotx <-  basic_coverage +
  geom_gc(bs.fa.seq = BS) +
  geom_gene(gtf.gr = target_annotation) +
  # geom_peak(bed.file = peak_file) +
  geom_ideogram(genome = genome_version$HG[n], plot.space = 0, highlight.centromere = TRUE
  )

# get consensus peak file
peak_file <- 'metaPlotR/newPeakFile.txt'
basic_coverage +
  geom_gene(gtf.gr = target_annotation) +
  # geom_peak(bed.file = peak_file) +
  geom_ideogram(genome = "hg19", plot.space = 0)

head(track_df)
ggcoverage(
  data = track_df,
  plot.type = "joint",
  facet.key = "Type",
  group.key = "Type",
  # mark.region = mark_region,
  range.position = "out"
)
ggcoverage(
  data = track_df, color = "auto", 
  facet.y.scale = "fixed",
  #mark.region = mark_region,
  show.mark.label = TRUE,
  
  plot.type = "joint",
  facet.key = "Group",
  group.key = "Type",
  joint.avg = TRUE,
  range.position = "in"
)
############
# metaPlotR #
############
library(getopt)
spec <- matrix(c(
  'help',  'h',  0,  "logical",
  'dist',  'd',  1,  "character",
  'name',  'n',  1,  "character",
  'od',    'o',  1,  "character"),
  byrow = TRUE, ncol = 4)
opt <- getopt(spec)

## usages
print_usage <- function(spec=NULL){
  cat("
  Usage example:
  Rscript metaPlot.R  --dist m6a.dist.measures.txt --name m6a --od  ./
  Options:
    --help   -h      NULL get this help
    --dist   character  m6a.dist.measures.txt [forced]
    --name   character  sample name
    --od     character  outdir [forced]
      \n")
  q(status=1)
}

## default paramter setting
if (!is.null(opt$help) |is.null(opt$dist) |is.null(opt$od)) { print_usage(spec) }
if (!file.exists(opt$od))  dir.create(opt$od)

# 加载包
library(ggplot2)
library(scales)

# 读取数据
m6a.dist <- read.delim(opt$dist, header = T)
# m6a.dist <- read.delim('metaPlotR/GSM3489019_PCIF1KO1.m6A.DRACH.FDR01.bed.gz.m6a.dist.measures.txt', header = T)
dim(m6a.dist)
head(m6a.dist)

# Determine longest length transcript for each gene
trx_len <- m6a.dist$utr5_size + m6a.dist$cds_size + m6a.dist$utr3_size
temp <- data.frame(m6a.dist$gene_name, m6a.dist$refseqID, trx_len)
colnames(temp) <- c("gene_name", "gid", "trx_len")
temp.df <- temp[order(temp$gene_name, temp$gid, -temp$trx_len),]
temp.df <- temp[!duplicated(temp$gene_name),]

# limit m6a data to one transcript per gene (longest)
m6a.dist <- m6a.dist[m6a.dist$refseqID %in% temp.df$gid,]

# View size of our dataset (rows, columns)
dim(m6a.dist)

# re-scale the widths of the 5'UTR and 3'UTR relative to the CDS
utr5.SF <- median(m6a.dist$utr5_size, na.rm = T)/median(m6a.dist$cds_size, na.rm = T)
utr3.SF <- median(m6a.dist$utr3_size, na.rm = T)/median(m6a.dist$cds_size, na.rm = T)

# assign the regions to new dataframes
utr5.m6a.dist <- m6a.dist[m6a.dist$rel_location < 1, ]
cds.m6a.dist <- m6a.dist [m6a.dist$rel_location < 2 & m6a.dist$rel_location >= 1, ]
utr3.m6a.dist <- m6a.dist[m6a.dist$rel_location >= 2, ]

# rescale 5'UTR and 3'UTR
utr5.m6a.dist$rel_location <- rescale(utr5.m6a.dist$rel_location, to = c(1-utr5.SF, 1), from = c(0,1))
utr3.m6a.dist$rel_location <- rescale(utr3.m6a.dist$rel_location, to = c(2, 2+utr3.SF), from = c(2,3))

m6a.metagene.coord <- c(utr5.m6a.dist$rel_location, cds.m6a.dist$rel_location, utr3.m6a.dist$rel_location)
df <- data.frame(m6a.metagene.coord = m6a.metagene.coord)

##---------------- Visualizing the metagene
# A smooth density plot
setwd(opt$od)
p <- ggplot(df,aes(x =m6a.metagene.coord)) + geom_density(aes(colour="red")) + geom_vline(xintercept = 1:2, col = "grey") +
  theme_bw() + xlim(0, 3) + theme(legend.position="none")
p

png(filename = paste0(opt$name,".png"),width = 800,height = 600,res=200)
print(p)
dev.off()

write.table(m6a.metagene.coord, file = paste0(opt$name,"_m6a.metagene.coord.txt"),row.names = F,sep = "\t",quote = F)



library(devtools)
# install_github("scottzijiezhang/MeRIPtools")
library("MeRIPtools")

### THIS PORTION OF THIS SCRIPT IS RUN IN R USING A .R SCRIPT ###
#############################################################################
## LIBRARIES TO USE
library(ggfortify)
library(cluster)
library(DESeq2)
library(apeglm)
library(limma)
library(edgeR)

# rm(list=ls())
options(stringsAsFactors = F)
#BiocManager::install("ChIPseeker",ask = F)
library(ChIPseeker)
library(org.Hs.eg.db)
library(clusterProfiler)
library(GenomicFeatures)
# annotation file gtf
txdb <- makeTxDbFromGFF(file = "data/KO1/Mus_musculus.GRCm38.101.gtf.gz", format="gtf",
                        dataSource="Ensembl", organism="Mus musculus")
# overlap
peak_file <- readPeakFile("data/KO1/Mod.bed")
## peak annotation
#region select："Promoter", "5UTR", "3UTR", "Exon", "Intron", "Downstream", "Intergenic"
peak_anno <- annotatePeak(peak_file,
                          tssRegion = c(-3000, 3000),
                          TxDb = txdb,
                          assignGenomicAnnotation = TRUE,
                          genomicAnnotationPriority = c("5UTR", "3UTR", "Exon","Intron","Intergenic"),
                          addFlankGeneInfo = TRUE,
                          flankDistance = 5000)

pdf(file = "KO1.Anno.Pie.pdf",width = 6,height = 5)
plotAnnoPie(peak_anno)
dev.off()
png(filename="KO1.Anno.Pie.png" ,width=1000, height=850, res=200)
plotAnnoPie(peak_anno)
dev.off()
# # genomic features imported into named list
# stBedFiles <- list(system.file("extdata", "m6A_mm10_exomePeak_1000peaks_bed12.bed", package="Guitar"))
# # With the following script, we may generate the transcriptomic distribution of genomic features to be tested, and the result will be automatically saved into a PDF file under the working directory with prefix “example”. With the GuitarPlot function, the gene annotation can be downloaded from internet automatically with a genome assembly number provided; however, this feature requires working internet and might take a longer time. The toy Guitar coordinates generated internally should never be re-used in other real data analysis.
# count <- GuitarPlot(txGenomeVer = "mm10",
#                     stBedFiles = stBedFiles,
#                     miscOutFilePrefix = NA)
# # In a more efficent protocol, in order to re-use the gene annotation and Guitar coordinates, you will have to build Guitar Coordiantes from a txdb object in a separate step. The transcriptDb contains the gene annotation information and can be obtained in a number of ways, .e.g, download the complete gene annotation of species from UCSC automatically, which might takes a few minutes. In the following analysis, we load the Txdb object from a toy dataset provided with the Guitar package. Please note that this is only a very small part of the complete hg19 transcriptome, and the Txdb object provided with Guitar package should not be used in real data analysis. With a TxDb object that contains gene annotation information, we in the next build Guitar coordiantes, which is essentially a bridge connects the transcriptomic landmarks and genomic coordinates.
# txdb_file <- system.file("extdata", "mm10_toy.sqlite", package="Guitar")
# txdb <- loadDb(txdb_file)
# guitarTxdb <- makeGuitarTxdb(txdb = txdb, txPrimaryOnly = FALSE)
# # You may now generate the Guitar plot from the named list of genome-based features.
# GuitarPlot(txTxdb = txdb,stBedFiles = stBedFiles,miscOutFilePrefix = "example")
# # Alternatively, you may also optionally include the promoter DNA region and tail DNA region on the 5’ and 3’ side of a transcript in the plot with parameter headOrtail =TRUE.
# GuitarPlot(txTxdb = txdb, stBedFiles = stBedFiles, headOrtail = TRUE)
# GuitarPlot(txTxdb = txdb, stBedFiles = stBedFiles, headOrtail = TRUE, enableCI = FALSE)
# # import different data formats into a named list object.
# # These genomic features are using mm10 genome assembly
# stBedFiles <- list(system.file("extdata", "m6A_mm10_exomePeak_1000peaks_bed12.bed",
#                                package="Guitar"),
#                    system.file("extdata", "m6A_mm10_exomePeak_1000peaks_bed6.bed",
#                                package="Guitar"))
# # Build Guitar Coordinates
# txdb_file <- system.file("extdata", "mm10_toy.sqlite", package="Guitar")
# txdb <- loadDb(txdb_file)
# # Guitar Plot
# GuitarPlot(txTxdb = txdb,  stBedFiles = stBedFiles,   headOrtail = TRUE,    enableCI = FALSE,
#            mapFilterTranscript = TRUE, pltTxType = c("mrna"),   stGroupName = c("BED12","BED6"))
# # 3 Processing of sampling sites information
# # We can select parameters for site sampling.
# stGRangeLists = vector("list", length(stBedFiles))
# sitesPoints <- list()
# for (i in seq_len(length(stBedFiles))) {
#   stGRangeLists[[i]] <- blocks(import(stBedFiles[[i]]))#, format = 'bed'
# }
#
# for (i in seq_len(length(stGRangeLists))) {
#   sitesPoints[[i]] <- samplePoints(stGRangeLists[i], stSampleNum = 10,stAmblguity = 5,
#                                    pltTxType = c("mrna"),stSampleModle = "Equidistance",mapFilterTranscript = FALSE,guitarTxdb = guitarTxdb)
# }
# # 4 Guitar Coordinates – Transcriptomic Landmarks Projected on Genome
# # The guitarTxdb object contains the genome-projected transcriptome coordinates, which can be valuable for evaluating transcriptomic information related applications, such as checking the quality of MeRIP-Seq data. The Guitar coordinates are essentially the genomic projection of standardized transcript-based coordiantes, making a viable bridge beween the landmarks on transcript and genome-based coordinates. It is based on the txdb object input, extracts the transcript information in txdb, selects the transcripts that match the parameters according to the component parameters set by the user, and saves according to the transcript type (tx, mrna, ncrna).
# guitarTxdb <- makeGuitarTxdb(txdb = txdb,   txAmblguity = 5,
#                              txMrnaComponentProp = c(0.1,0.15,0.6,0.05,0.1),
#                              txLncrnaComponentProp = c(0.2,0.6,0.2),
#                              pltTxType = c("tx","mrna","ncrna"),
#                              txPrimaryOnly = FALSE)
# # 5 Check the Overlapping between Different Components
# # We can also check the distribution of the Guitar coordinates built.
# gcl <- list(guitarTxdb$tx$tx)
# GuitarPlot(txTxdb = txdb,  stGRangeLists = gcl,
#            stSampleNum = 200,  enableCI = TRUE, pltTxType = c("tx"),  txPrimaryOnly = FALSE
# )
# # Alternatively, we can extract the RNA components, check the distribution of tx components in the transcriptome.
# GuitarCoords <- guitarTxdb$tx$txComponentGRange
# type <- paste(mcols(GuitarCoords)$componentType,mcols(GuitarCoords)$txType)
# key <- unique(type)
# landmark <- list(1,2,3,4,5,6,7,8,9,10,11)
# names(landmark) <- key
# for (i in 1:length(key)) {
#   landmark[[i]] <- GuitarCoords[type==key[i]]
# }
# GuitarPlot(txTxdb = txdb,   stGRangeLists = landmark[1:3], pltTxType = c("tx"),
#            enableCI = FALSE)
# # Check the distribution of mRNA components in the transcriptome
# GuitarPlot(txTxdb = txdb,   stGRangeLists = landmark[4:8],
#            pltTxType = c("mrna"),
#            enableCI = FALSE)
# # Check the distribution of lncRNA components in the transcriptome
# GuitarPlot(txTxdb = txdb ,
#            stGRangeLists = landmark[9:11],
#            pltTxType = c("ncrna"),
#            enableCI = FALSE)



# #################### Replicability of processed data ########
# ## Define a function to draw a scatter plot for a pair of variables (samples) with density colors
# plotFun <- function(x,y){
#   dns <- densCols(x,y);
#   points(x,y, col="black", pch=".",cex=4,panel.first=grid());
#   abline(a=0, b=1, col="red")
# }
# 
# vsd_count=as.data.frame(assay(vst(dds,blind=F)))
# count.table2=as.matrix(subset(vsd_count,rownames(vsd_count)%in%rownames(res1))[1:3])
# 
# colnames(count.table2)=factor(c("mo replicate1","mo replicate2","mo replicate3"),levels =c("mo replicate3","mo replicate2","mo replicate1") )
# ## Plot the scatter plot for a few pairs of variables selected at random
# set.seed(123) # forces the random number generator to produce fixed results. Should generally not be used, except for the sake of demonstration with a particular selection.
# pairs(log2(count.table2[,sample(ncol(count.table2), 3)] + 1),
#       panel=plotFun)
# 
# ######### Gene Ontology enrichment of top 500 highly expressed genes #############
# library(clusterProfiler)
# library(org.Hs.eg.db)
# ####### M0 vs mo ########
# 
# ## read the differential genes between mo and M0 ##
# DEG_M0vsmo=read.csv("DEG_M0vsmo.csv")
# rownames(DEG_M0vsmo)=DEG_M0vsmo[,1]
# mo_M0_GO<- enrichGO(gene       = rownames(DEG_M0vsmo),
#                     OrgDb         = org.Hs.eg.db,
#                     universe      = rownames(vsd_count), 
#                     keyType       = 'SYMBOL',
#                     ont           = "BP",
#                     pAdjustMethod = "BH",
#                     pvalueCutoff  = 0.05,
#                     qvalueCutoff  = 0.05)
# ## read the differential genes between M0 and M1 ##
# DEG_M1vsM0=read.csv("DEG_M1vsM0.csv")
# rownames(DEG_M1vsM0)=DEG_M1vsM0[,1]
# M0_M1_GO<- enrichGO(gene       = rownames(DEG_M1vsM0),
#                     OrgDb         = org.Hs.eg.db,
#                     universe      = rownames(vsd_count), 
#                     keyType       = 'SYMBOL',
#                     ont           = "BP",
#                     pAdjustMethod = "BH",
#                     pvalueCutoff  = 0.05,
#                     qvalueCutoff  = 0.05)
# 
# ## read the differential genes between M0 and M2 ##
# DEG_M2vsM0=read.csv("DEG_M2vsM0.csv")
# rownames(DEG_M2vsM0)=DEG_M2vsM0[,1]
# M0_M2_GO<- enrichGO(gene       = rownames(DEG_M2vsM0),
#                     OrgDb         = org.Hs.eg.db,
#                     universe      = rownames(vsd_count), 
#                     keyType       = 'SYMBOL',
#                     ont           = "BP",
#                     pAdjustMethod = "BH",
#                     pvalueCutoff  = 0.05,
#                     qvalueCutoff  = 0.05)
# 
# ### Convert to DataFrame #####
# GO_mo_M0_RNAseq=as.data.frame(mo_M0_GO)
# GO_M0_M1_RNAseq=as.data.frame(M0_M1_GO)
# GO_M0_M2_RNAseq=as.data.frame(M0_M2_GO)
plot_ideogram <- function(gene, rank) {
  tryCatch({
    # Check if the rank is within the valid range of genome_version$HG
    if(rank > length(genome_version$HG) || rank < 1) {
      stop("Invalid rank: out of range.")
    }
    # Get the genome version based on rank
    genome_ver <- genome_version$HG[rank]
    # Ensure genome_ver is not NULL or empty
    if(is.null(genome_ver) || genome_ver == "") {
      stop("Invalid genome version.")
    }
    # Create the ideogram layer
    ideogram_layer <- geom_ideogram(genome = genome_ver, plot.space = 0, highlight.centromere = TRUE)
    # Print the ideogram data for debugging
    print(ideogram_layer)
    # Combine the basic coverage plot with the ideogram
    plot <- basic_coverage +
      geom_gc(bs.fa.seq = BS) +
      geom_gene(gtf.gr = target_annotation) +
      ideogram_layer # Added ideogram layer
    # Print the combined plot
    # print(plot)
    # Return the plot object
    return(plot)
  }, error = function(e) {
    message("Error occurred: ", e$message)
    return(NULL)  # Return NULL in case of error
  })
}

# Load necessary libraries
library(dplyr)
library(readr)
setwd('/data2/rluo4/EpiTrans/DataCollection')
# Specify the file path
file_path <- "/data2/rluo4/EpiTrans/DataCollection/hg38.m5C_reader.modSite.tar.gz"
# Read the compressed file and skip the first line
con <- gzfile(file_path, "rt")
lines <- readLines(con)
close(con)
# # Remove the unwanted string from each line/ Skip the first line which contains the unwanted string
data_lines <-  cleaned_lines <- gsub("human.hg38.modrbp.m5C.reader.bed\\^@\\^@\\^@\\^@\\^@\\^@\\^@\\^@\\^@\\^@\\^@\\^@\\^@\\^@\\^@\\^@\\^@\\^@\\^@\\^@\\^@\\^@\\^@\\^@\\^@\\^@\\^@\\^@\\^@\\^@\\^@\\^@\\^@\\^@\\^@\\^@\\^@\\^@\\^@\\^@\\^@\\^@\\^@\\^@\\^@\\^@\\^@\\^@\\^@\\^@\\^@", "", lines)
#lines[-1]
# Create a temporary file to store the cleaned data (without the first line)
temp_file <- tempfile()
# Write the cleaned data to the temporary file
writeLines(data_lines, temp_file)
# Read the cleaned data from the temporary file
data_lines <- fread(temp_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
# Clean up the temporary file
unlink(temp_file)
# Assign column names
colnames(data_lines) <- c("rbpSiteID", "rbpName", "narrowPeak", "broadPeak", "rbpRand", 
                          "clipExpNum", "dataSetIdList", "RNAmodID", "RNAmodLoc", 
                          "RNAmodType", "geneID", "transcriptID", "geneName", 
                          "geneType", "Region", "Conservation", "rbpType")# Display the first few rows with column names
head(data_lines)
