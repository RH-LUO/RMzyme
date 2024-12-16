#!/usr/local/bin/Rscript
# title: "RMzyme, a comprehensive platform of regulations of RNA modifying enzymes in human based on multiomics data collections"
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

library(org.Hs.eg.db)
library(clusterProfiler)
# PTM data from databases -------------------------------------------------
PTM_info <- read_excel(paste0('/data2/rluo4/RPMfunc/Output/summary/ALL_RMP_PTM1015.xlsx'))
# clinical data Cancer Cell(Article) ----------------------------------------------------
cancercellcli <- read_xlsx("/data2/rluo4/EpiTrans/PTM/CPTAC/inputdata/cancer_cell_clinical.xlsx",sheet = 2)
uniprot <- read.csv('/data2/rluo4/All/uniprot-hs.tsv',fill = T,header = T,sep = '\t')
uniprot$geneSymbol <- str_split(uniprot$Gene.Names, ' ', simplify = T)[,1]
PTM_mutation <- read_excel('/data2/rluo4/EpiTrans/DataCollection/Table_S3_Mutations_Near_PTMs.tsv.xlsx',sheet = 1)
PTM_mutation$geneSymbol <- uniprot$Gene.Names[match(PTM_mutation$Uniprot_ID, uniprot$Entry)]
PTM_mutation$geneSymbol <- str_split(PTM_mutation$geneSymbol, " ", simplify = T)[,1]
length(unique(PTM_mutation$geneSymbol))#[1] 2284
table(unique(PTM_mutation$geneSymbol) %in% RMP_update$RMP)
# View(PTM_mutation[PTM_mutation$Mutation_AA_Pos=='700' & PTM_mutation$geneSymbol=='SF3B1',])
# PTM_full <- fread('/data2/rluo4/EpiTrans/PTM/CELL/ref/PTM_full.tsv')
# package -----------------------------------------------------------------
library(dplyr)
library(tidyverse)
library(stringr)
library(maftools)
load("/data2/rluo4/EpiTrans/PTM/CPTAC/inputdata/RBP.rdata")
RMP_update <- read_excel('/data2/rluo4/EpiTrans/DataCollection/RMP_update.xlsx') # from RM.R
# GSVA contrustion by R 4.3.3 -------------------------------------------
RBPgene <- dplyr::pull(RBP,1)
mufl <-  list.files("/data2/rluo4/EpiTrans/PTM/CPTAC/inputdata/Mutation_BCM_v1/Mutation_BCM_v1/",pattern = ".maf")
mulistT <- lapply(paste0("/data2/rluo4/EpiTrans/PTM/CPTAC/inputdata/Mutation_BCM_v1/Mutation_BCM_v1/",mufl), function(x){
  print(x)
  read.maf(x) 
})
for (tumor in 1:length(mufl)) {
  names(mulistT)[[tumor]] <- strsplit(mufl[tumor],'_')[[1]][1]
}
names(mulistT)[[2]] <- "CCRCC"
PTM_PDC_full <- read.table(file = '/data2/rluo4/EpiTrans/PTM/CELL/ref/PTM_PDC_full.tsv',sep = '\t', header = T)
table(names(mulistT) %in% unique(PTM_PDC_full$Disease))
# Extract the amino acid position from HGVSp_Short (e.g., p.V600E -> 600)
# maf_data <- mulistT$BRCA@data %>%
#   mutate(aa_position = gsub("p\\.[A-Za-z]+(\\d+)[A-Za-z]+", "\\1", aaChange))
PDC_PTM_mutation <- NULL
library(purrr)
extract_positions <- function(hgvsp) {
  # Frame Shift Mutations (e.g., p.A379Gfs*18)
  if (grepl("fs", hgvsp)) {
    match <- regmatches(hgvsp, regexec("p\\.[A-Za-z](\\d+)[A-Za-z]*fs\\*\\d+", hgvsp))
    if (length(match[[1]]) > 1) {
      return(c(as.numeric(match[[1]][2]), as.numeric(match[[1]][2])))  # Single position for frameshift
    }
  }
  # Deletion Mutations with Range (e.g., p.L503_D507del)
  if (grepl("del", hgvsp)) {
    match <- regmatches(hgvsp, regexec("p\\.[A-Za-z](\\d+)_?[A-Za-z]*(\\d+)?del", hgvsp))
    if (length(match[[1]]) > 1) {
      pos_start <- as.numeric(match[[1]][2])
      pos_end <- ifelse(is.na(match[[1]][3]), pos_start, as.numeric(match[[1]][3]))  # Handle single/multiple sites
      return(c(pos_start, pos_end))
    }
  }
  # Insertion and Delins Mutations (e.g., p.S57_A58insPSS, p.Q219_Y221delinsH)
  if (grepl("ins|delins", hgvsp)) {
    match <- regmatches(hgvsp, regexec("p\\.[A-Za-z](\\d+)_?[A-Za-z]*(\\d+)?[A-Za-z]*ins[A-Za-z]*", hgvsp))
    if (length(match[[1]]) > 1) {
      pos_start <- as.numeric(match[[1]][2])
      pos_end <- ifelse(is.na(match[[1]][3]), pos_start, as.numeric(match[[1]][3]))  # Handle single/multiple sites
      return(c(pos_start, pos_end))
    }
  }
  # Stop Codon Mutations (e.g., p.C492*, p.M1*)
  if (grepl("\\*", hgvsp)) {
    match <- regmatches(hgvsp, regexec("p\\.[A-Za-z](\\d+)\\*", hgvsp))  # Handle nonsense mutations
    if (length(match[[1]]) > 1) {
      return(c(as.numeric(match[[1]][2]), as.numeric(match[[1]][2])))  # Single position for nonsense mutation
    }
  }
  # Special case for Translation Start Site
  if (grepl("^p\\.M1\\*?$", hgvsp)) {
    return(c(1, 1))  # Explicitly return position 1 for translation start site
  }
  # Missense Mutations (e.g., p.M102T)
  match <- regmatches(hgvsp, regexec("p\\.[A-Za-z](\\d+)[A-Za-z]", hgvsp))  # Handle missense mutations
  if (length(match[[1]]) > 1) {
    return(c(as.numeric(match[[1]][2]), as.numeric(match[[1]][2])))  # Single position for missense mutation
  }
  # Handle Stop Codon Delins Mutations (e.g., p.*1910delinsNRDSLEC*)
  if (grepl("delins", hgvsp)) {
    match <- regmatches(hgvsp, regexec("p\\.\\*(\\d+)delins[A-Za-z]+\\*", hgvsp))
    if (length(match[[1]]) > 1) {
      return(c(as.numeric(match[[1]][2]), as.numeric(match[[1]][2])))  # Position for stop codon with delins
    }
  }
  # Translation Start Site Mutations (e.g., p.M1? or p.M1*)
  if (grepl("p\\.M1", hgvsp)) {
    return(c(1, 1))  # Translation start site at position 1
  }
  # If no match found, return NA
  return(c(NA, NA))
}
for (tumor in names(mulistT)) {
  # Function to extract positions based on different mutation types
  # Updated mutation classification and position extraction
  print(tumor)
  print(table(mulistT[[tumor]]@data$Variant_Classification))
  maf_data <- mulistT[[tumor]]@data %>%
    dplyr::select(gene_name, Chromosome, Start_Position, End_Position, Strand,
                  Tumor_Sample_Barcode, Hugo_Symbol, n_ref_count, n_alt_count,
                  t_ref_count, t_alt_count, TumorVAF, AAChange.refGene, txChange, aaChange, Variant_Classification) %>%
    mutate(
      # Classify mutation types
      mutation_type = case_when(
        grepl("fs", aaChange) ~ "frameshift",
        grepl("ins", aaChange) ~ "insertion",
        grepl("delins", aaChange) ~ "delins",
        grepl("del", aaChange) ~ "deletion",
        grepl("\\*", aaChange) & !grepl("fs|delins", aaChange) ~ "nonsense",  # Nonsense mutation if it ends with "*"
        grepl("delins", aaChange) & grepl("\\*", aaChange) ~ "nonsense",  # Handle stop codon delins as nonsense
        grepl("^p\\.M1\\*?$", aaChange) ~ "translation_start_site",  # Only p.M1* or p.M1 for start site
        TRUE ~ "missense"  # All other changes with amino acid substitutions
      ),
      # Extract start and end positions using the improved extract_positions function
      position_info = purrr::map(aaChange, extract_positions),
      position_start = purrr::map_dbl(position_info, 1),
      position_end = purrr::map_dbl(position_info, 2)
    ) %>%
    # Ensure that missing positions are filled with aa_position where applicable
    mutate(
      aa_position = gsub("p\\.[A-Za-z]+(\\d+)[A-Za-z]+", "\\1", aaChange),
      position_start = coalesce(as.numeric(position_start), as.numeric(aa_position)),
      position_end = coalesce(as.numeric(position_end), as.numeric(aa_position))
    )
  # Output the result
  # head(maf_data)
  index <- is.na(maf_data$position_start)
  table(index)
  maf_data <- maf_data[!index,] # delete the rows with Splice_Site
  
  maf_data$Disease <- tumor
  maf_data$Uniprot_ID <- uniprot$Entry[match(maf_data$Hugo_Symbol, uniprot$geneSymbol)]
  library(dplyr)
  # Step 1: Prepare maf_data and PTM_data
  # Extract relevant columns from maf_data
  maf_filtered <- as.data.frame(maf_data) %>%
    dplyr::select(CPTAC_Mutation_ID = Tumor_Sample_Barcode, Disease,
                  Uniprot_ID, Mutation_AA_Pos = position_start,  # Use extracted start position for mutation
                  # PTM_AA_Pos = position_end,          # Assume position_end matches PTM position
                  geneSymbol = Hugo_Symbol)
  # Prepare PTM_data
  PTM_filtered <- PTM_PDC_full %>%
    filter(Disease == tumor) %>%
    dplyr::select(Uniprot_ID = UNIPROT, PTM_AA_Pos = location, PTM_Type = feature, geneSymbol = SYMBOL) #as.numeric(gsub("K", "", PTM_AA_Pos)))  # Convert location to numeric if needed
  # Step 2: Ensure Mutation_AA_Pos in maf_filtered is also numeric (if needed)
  maf_filtered <- maf_filtered %>%
    mutate(Mutation_AA_Pos = as.numeric(Mutation_AA_Pos))
  PTM_filtered <- PTM_filtered %>%
    mutate(PTM_AA_Pos = as.numeric(substr(PTM_AA_Pos, 2, nchar(PTM_AA_Pos))))
  # coinciding_mutations <- maf_filtered %>%
  #   inner_join(PTM_filtered, by = c("Uniprot_ID", "geneSymbol", "PTM_AA_Pos")) %>%
  #   mutate(Mutation_Distance_from_PTM = Mutation_AA_Pos - PTM_AA_Pos)
  # Step 3: Merge DataFrames to identify coinciding mutations and surrounding mutations
  coinciding_mutations <- maf_filtered %>%
    inner_join(PTM_filtered, by = c("Uniprot_ID", "geneSymbol")) %>%
    # Calculate the distance between mutation and PTM sites
    mutate(Mutation_Distance_from_PTM = (Mutation_AA_Pos - PTM_AA_Pos)) %>%#abs(Mutation_AA_Pos - PTM_AA_Pos)) %>%
    # Filter for mutations that are either coinciding or within 10 residues of PTM
    filter(abs(Mutation_Distance_from_PTM) <= 10)
  # Step 4: Add a new column indicating the nature of the mutation relative to the PTM site
  coinciding_mutations <- coinciding_mutations %>%
    mutate(Mutation_Type_Relative_to_PTM = case_when(
      Mutation_Distance_from_PTM == 0 ~ "Coinciding",
      Mutation_Distance_from_PTM <= 10 ~ "Surrounding",
      TRUE ~ "Not near PTM"
    ))
  PDC_PTM_mutation <<- rbind(PDC_PTM_mutation, coinciding_mutations)
}
PDC_PTM_mutation$Tissue <-  gsub('Acute meyloid leukemia|Acute meyloid leukemia|Acute monocytic leukemia|Acute myeloid leukemia|Chronic lymphocytic leukemia|Chronic myeloid leukemia|Leukemia|Myelogenous leukemia', 'bone marrow',
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
                                                                                                                                               PDC_PTM_mutation$Disease,         
                                                                                                                                          )))))))))))))))))))))))))

write.table(PDC_PTM_mutation, file = '/data2/rluo4/RPMfunc/Output/summary/PTM_all_mutation.txt', row.names = F, quote = F, sep = '\t')
PDC_PTM_mutation <- read.table('/data2/rluo4/RPMfunc/Output/summary/PTM_all_mutation.txt',header = T, sep = '\t')
write_xlsx(x = PDC_PTM_mutation, path = paste0('/data2/rluo4/RPMfunc/Output/summary/PDC_PTM_mutation_summary.xlsx'))
# Table S18
summary_data <- PDC_PTM_mutation %>%
  group_by(Disease, PTM_Type) %>%
  summarise(Mutation_Count = n()) %>%
  ungroup()
summary_data <- summary_data[order(summary_data$Mutation_Count, decreasing = F),]
# Reorder PTMs factor based on Percent_multiple_sites
summary_data$Disease <- factor(summary_data$Disease, levels = unique(summary_data$Disease))
library(ggprism)
ggplot(summary_data, aes(x = Mutation_Count, y = Disease)) +
  geom_col(aes(fill = PTM_Type), width = 0.7) +  # Bar plot with PTM types as fill color
  labs(x = "PTM Mutation Counts",
       y = "", # title = "Percentage of Regulators with Multiple PTM Sites by PTM Type"
  ) + scale_fill_brewer(palette = "Set3") +  # Use a color palette for PTM types
  # scale_fill_manual(values = color_cluster)+ #scale_fill_brewer(palette = "Blues")+
  theme(axis.text=element_text(face = "bold", color="black")) + #labs(title = "Predicted neopeptides in TCGA-SKCM")+ 
  # coord_flip() +
  theme_prism(#palette = "winter_bright",
    base_fontface = "plain", # ×ÖÌåÑùÊ½£¬¿ÉÑ¡ bold, plain, italic
    base_family = "serif", # ×ÖÌå¸ñÊ½£¬¿ÉÑ¡ serif, sans, mono, ArialµÈ
    base_size = 16,  # Í¼ÐÎµÄ×ÖÌå´óÐ¡
    base_line_size = 0.8, # ×ø±êÖáµÄ´ÖÏ¸
    axis_text_angle = 45) #+ theme(legend.position='none')
ggsave('/data2/rluo4/RPMfunc/WebSite/PTM_Mutation.png', height = 5.5, width = 7, dpi = 500)
# Fig. 6d
RMP_PTM_mutation <- subset(PDC_PTM_mutation, geneSymbol %in% RMP_update$RMP)
table(RMP_PTM_mutation$PTM_Type)
sort(table(RMP_PTM_mutation$geneSymbol))
unique(RMP_PTM_mutation$geneSymbol)
# write.table(RMP_PTM_mutation, file = '')
# Display the result
# print(RMP_PTM_mutation)
# package -----------------------------------------------------------------
# BiocManager::install("GSVA")
library(GSVA)
library(GSEABase)
library(dplyr)
library(tidyverse)
library(stringr)
library(maftools)
plot_cols = c( '#96dce0','#FB9A99',"#66A61E",'#CAB2D6', '#6A3D9A','#FDBF6F','#FF7F00',"#85cded","#FAEFD1")
vc_cols = c(
  'Frame_Shift_Del' = '#e4e4e4',   'Missense_Mutation' = "#c2a5cf",   'Nonsense_Mutation' = '#FFB6C1',
  'Multi_Hit' = "#9fd7c6", 'Frame_Shift_Ins' =  "#85cded" ,  'In_Frame_Ins' = '#E6E6FA',      
  'Splice_Site' = '#96dce0', 'In_Frame_Del' = '#f2938c',  'Translation_Start_Site' = '#E31A1C'
)
names(plot_cols) <- names(vc_cols)
plot_cols
metabolic_pathways_to_plot <- str_split("HALLMARK_XENOBIOTIC_METABOLISM
HALLMARK_FATTY_ACID_METABOLISM
HALLMARK_OXIDATIVE_PHOSPHORYLATION
KEGG_GLYCOLYSIS_GLUCONEOGENESIS
KEGG_PYRUVATE_METABOLISM
KEGG_PROPANOATE_METABOLISM
KEGG_BUTANOATE_METABOLISM
KEGG_BIOSYNTHESIS_OF_UNSATURATED_FATTY_ACIDS","\n")[[1]]

immune_pathways_to_plot <- strsplit("HALLMARK_IL6_JAK_STAT3_SIGNALING
HALLMARK_INTERFERON_ALPHA_RESPONSE
HALLMARK_INTERFERON_GAMMA_RESPONSE
HALLMARK_COMPLEMENT
HALLMARK_INFLAMMATORY_RESPONSE
HALLMARK_IL2_STAT5_SIGNALING","\n")[[1]]
#   dir.create(file.path("/data2/rluo4/EpiTrans/PTM/CPTAC/fig/mutation/", fold))
# }
kcdf="Poisson"
require(clusterProfiler)
library(limma)
options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 120000 * 1024^2)
library(future)
my_genesets <- read.csv('/data2/rluo4/All/my_genesets.csv')
my_genesets <- my_genesets[, -1]#[-grep("KEGG", my_genesets$gs_name),-1]
# my_genesets <- rbind(Senescense, Metaplasia, H_genesets, cellcycle.gmt, reactome.gmt)
colnames(my_genesets) <- c("term","gene")
my_genesets$term <- as.factor(my_genesets$term)
unique(my_genesets$term)
table(unique(my_genesets$term) %in% metabolic_pathways_to_plot)
# my.gmt <- split(pathways_to_plot$gene, pathways_to_plot$term)
# pathways_to_plot <- pathways_to_plot[!is.na(pathways_to_plot$phos),]
# GSVA + Mutation: 179 RBP group -----------------------------------------------------
load("/data2/rluo4/EpiTrans/PTM/CPTAC/inputdata/phossymbolistT.rdata")
load("/data2/rluo4/EpiTrans/PTM/CPTAC/inputdata/acetylist.rdata")
load("/data2/rluo4/EpiTrans/PTM/CPTAC/inputdata/glycolist.rdata")
load("/data2/rluo4/EpiTrans/PTM/CPTAC/inputdata/ubiquitylist.rdata")
# Helper function to prepare the GeneSetCollection for metabolic and immune pathways with unique set names
create_gene_sets <- function(ptm_type_data, gene_column, metabolic_pathways_to_plot, immune_pathways_to_plot, my_genesets, ptm_type) {
  gene_df <- data.frame(ptm_type = rownames(ptm_type_data), gene = str_split(rownames(ptm_type_data), ' ', simplify = TRUE)[, 1])
  # Metabolic Pathways
  metabolic_genes <- left_join(my_genesets[my_genesets$term %in% metabolic_pathways_to_plot, ], gene_df, by = 'gene') %>%
    na.omit() %>%
    pull(ptm_type) %>%
    unique()
  # Immune Pathways
  immune_genes <- left_join(my_genesets[my_genesets$term %in% immune_pathways_to_plot, ], gene_df, by = 'gene') %>%
    na.omit() %>%
    pull(ptm_type) %>%
    unique()
  
  # # Create GeneSetCollections with unique set names
  # metabolic_geneSet <- GeneSetCollection(lapply(seq_along(metabolic_genes), function(i) {
  #   gene <- metabolic_genes[i]
  #   setName <- paste0(gene_column, "_", ptm_type, "_metabolic_", i)  # Unique set name
  #   GeneSet(gene, geneIdType = EntrezIdentifier(), setName = setName)
  # }))
  # 
  # immune_geneSet <- GeneSetCollection(lapply(seq_along(immune_genes), function(i) {
  #   gene <- immune_genes[i]
  #   setName <- paste0(gene_column, "_", ptm_type, "_immune_", i)  # Unique set name
  #   GeneSet(gene, geneIdType = EntrezIdentifier(), setName = setName)
  # }))
  pathway_genes <- list( metabolic_genes); names(pathway_genes) <- paste0('Metabolic_',ptm_type)
  metabolic_geneSet <- GeneSetCollection(mapply(function(geneIds, keggId){
    GeneSet(geneIds, geneIdType=EntrezIdentifier(),
            collectionType=KEGGCollection(keggId),
            setName=keggId)
  }, pathway_genes, names(pathway_genes)))
  
  pathway_genes <- list( immune_genes); names(pathway_genes) <- paste0('Immune_',ptm_type)
  immune_geneSet <- GeneSetCollection(mapply(function(geneIds, keggId){
    GeneSet(geneIds, geneIdType=EntrezIdentifier(),
            collectionType=KEGGCollection(keggId),
            setName=keggId)
  }, pathway_genes, names(pathway_genes)))
  
  return(list(metabolic_geneSet = metabolic_geneSet, immune_geneSet = immune_geneSet))
}
# Updated GSVA processing function
process_PTM_data_and_GSVA <- function(tumor, phossymbolistT, acetylist, glycolist, ubiquitylist, my_genesets, metabolic_pathways_to_plot, immune_pathways_to_plot) {
  
  # Initialize Data Container
  data_PTM <- vector("list", 4)
  names(data_PTM) <- c('phos', 'acety', 'glyco', 'ubiquity')
  
  # Containers to hold the GeneSetCollection for each PTM type
  metabolic_geneSets <- list()
  immune_geneSets <- list()
  
  # Helper function to ensure numeric matrices for GSVA input
  ensure_numeric <- function(mat) {
    mat <- as.matrix(mat) # Ensure it's a matrix
    mode(mat) <- "numeric" # Convert to numeric mode
    return(mat)
  }
  
  # Process Phosphorylation Data
  if (tumor %in% names(phossymbolistT)) {
    data_PTM[['phos']] <- t(phossymbolistT[[tumor]]) %>% as.matrix()
    data_PTM[['phos']] <- ensure_numeric(data_PTM[['phos']])
    gs <- create_gene_sets(data_PTM[['phos']], "Metabolic_phos", metabolic_pathways_to_plot, immune_pathways_to_plot, my_genesets, "phos")
    metabolic_geneSets[['phos']] <- gs$metabolic_geneSet
    immune_geneSets[['phos']] <- gs$immune_geneSet
  }
  
  # Process Acetylation Data
  if (tumor %in% names(acetylist)) {
    data_PTM[['acety']] <- t(acetylist[[tumor]]) %>% as.matrix()
    data_PTM[['acety']] <- ensure_numeric(data_PTM[['acety']])
    gs <- create_gene_sets(data_PTM[['acety']], "Metabolic_acety", metabolic_pathways_to_plot, immune_pathways_to_plot, my_genesets, "acety")
    metabolic_geneSets[['acety']] <- gs$metabolic_geneSet
    immune_geneSets[['acety']] <- gs$immune_geneSet
  }
  
  # Process Glycosylation Data
  if (tumor %in% names(glycolist)) {
    data_PTM[['glyco']] <- t(glycolist[[tumor]]) %>% as.matrix()
    data_PTM[['glyco']] <- ensure_numeric(data_PTM[['glyco']])
    gs <- create_gene_sets(data_PTM[['glyco']], "Metabolic_glyco", metabolic_pathways_to_plot, immune_pathways_to_plot, my_genesets, "glyco")
    metabolic_geneSets[['glyco']] <- gs$metabolic_geneSet
    immune_geneSets[['glyco']] <- gs$immune_geneSet
  }
  
  # Process Ubiquitylation Data
  if (tumor %in% names(ubiquitylist)) {
    data_PTM[['ubiquity']] <- t(ubiquitylist[[tumor]]) %>% as.matrix()
    data_PTM[['ubiquity']] <- ensure_numeric(data_PTM[['ubiquity']])
    gs <- create_gene_sets(data_PTM[['ubiquity']], "Metabolic_ubiquity", metabolic_pathways_to_plot, immune_pathways_to_plot, my_genesets, "ubiquity")
    metabolic_geneSets[['ubiquity']] <- gs$metabolic_geneSet
    immune_geneSets[['ubiquity']] <- gs$immune_geneSet
  }
  
  # # Run GSVA for all PTM data types
  # gsva_results <- lapply(names(data_PTM), function(ptm_type) {
  #   if (!is.null(data_PTM[[ptm_type]])) {
  #     list(
  #       metabolic = gsva(data_PTM[[ptm_type]], metabolic_geneSets[[ptm_type]], method = 'ssgsea', kcdf = 'Poisson', parallel.sz = 32),
  #       immune = gsva(data_PTM[[ptm_type]], immune_geneSets[[ptm_type]], method = 'ssgsea', kcdf = 'Poisson', parallel.sz = 32)
  #     )
  #   }
  # })
  # Run GSVA for all PTM data types and organize results
  gsva_results <- list()  # Initialize an empty list to store GSVA results
  for (ptm_type in names(data_PTM)) {
    if (!is.null(data_PTM[[ptm_type]])) {
      
      # Run GSVA for metabolic pathways
      gsva_metabolic <- gsva(
        data_PTM[[ptm_type]],
        metabolic_geneSets[[ptm_type]],
        method = 'ssgsea',
        kcdf = 'Poisson',
        parallel.sz = 32
      )
      # Run GSVA for immune pathways
      gsva_immune <- gsva(
        data_PTM[[ptm_type]],
        immune_geneSets[[ptm_type]],
        method = 'ssgsea',
        kcdf = 'Poisson',
        parallel.sz = 32
      )
      # Store the results in the gsva_results list with clear naming
      gsva_results[[ptm_type]] <- list(
        metabolic = gsva_metabolic,
        immune = gsva_immune
      )
    }
  }
  # Return the final organized results
  return(gsva_results)
}
col1 = c('#FDBF6F', "#9fd7c6");names(col1) = c('High','Low')
col2 = c( 'brown', 'beige' );names(col2) = c('High','Low')

for (tumor in names(mulistT)) {
  print(tumor)
  # Run the GSVA analysis
  gsva_results <- process_PTM_data_and_GSVA(
    tumor, 
    phossymbolistT, 
    acetylist, 
    glycolist, 
    ubiquitylist, 
    my_genesets, 
    metabolic_pathways_to_plot, 
    immune_pathways_to_plot
  )
  # Access GSVA results for each PTM type with error handling
  Metabolic_phos_gsva <- if ("phos" %in% names(gsva_results)) gsva_results$phos$metabolic %>% t() else NULL
  Metabolic_acety_gsva <- if ("acety" %in% names(gsva_results)) gsva_results$acety$metabolic %>% t() else NULL
  Metabolic_glyco_gsva <- if ("glyco" %in% names(gsva_results)) gsva_results$glyco$metabolic %>% t() else NULL
  Metabolic_ubiquity_gsva <- if ("ubiquity" %in% names(gsva_results)) gsva_results$ubiquity$metabolic %>% t() else NULL
  
  immune_phos_gsva <- if ("phos" %in% names(gsva_results)) gsva_results$phos$immune %>% t() else NULL
  immune_acety_gsva <- if ("acety" %in% names(gsva_results)) gsva_results$acety$immune %>% t() else NULL
  immune_glyco_gsva <- if ("glyco" %in% names(gsva_results)) gsva_results$glyco$immune %>% t() else NULL
  immune_ubiquity_gsva <- if ("ubiquity" %in% names(gsva_results)) gsva_results$ubiquity$immune %>% t() else NULL
  # Mutation Data and Clinical Info
  mudata <- mulistT[[tumor]]
  # cli <- mudata@clinical.data
  cli <- cancercellcli[cancercellcli$tumor_code %in% tumor & is.na(cancercellcli$reason_for_exclusion),]
  # Merge clinical data with GSVA results
  cligsva <- mudata@clinical.data
  # Add Metabolic and Immune scores if they are available
  if (!is.null(Metabolic_phos_gsva) && !is.null(immune_phos_gsva)) {
    table(rownames(Metabolic_phos_gsva) %in% cligsva$Tumor_Sample_Barcode)
    cligsva <- merge(cligsva, data.frame(Tumor_Sample_Barcode = rownames(Metabolic_phos_gsva), 
                                         Metabolic_phos = Metabolic_phos_gsva,
                                         Immune_phos = immune_phos_gsva), 
                     by = "Tumor_Sample_Barcode", all.x = TRUE)
  }
  if (!is.null(Metabolic_acety_gsva) && !is.null(immune_acety_gsva)) {
    index <- match(rownames(Metabolic_acety_gsva), cli$`specimen/aliquout_id_protein_tumor`) %>% na.omit()
    cligsva <- merge(cligsva, data.frame(Tumor_Sample_Barcode = cli$Case_ID[index], 
                                         Metabolic_acety = Metabolic_acety_gsva[index],
                                         Immune_acety = immune_acety_gsva[index]), 
                     by = "Tumor_Sample_Barcode", all.x = TRUE)
  }
  if (!is.null(Metabolic_glyco_gsva) && !is.null(immune_glyco_gsva)) {
    t_strings <- grepl("\\.T", rownames(Metabolic_glyco_gsva))
    index <- match(str_split(rownames(Metabolic_glyco_gsva)[t_strings], '[.]', simplify = T)[,2], 
                   cli$Case_ID) %>% na.omit() %>% unique()
    # index <- str_split(rownames(Metabolic_glyco_gsva)[t_strings], '[.]', simplify = T)[,2]
    cligsva <- merge(cligsva, data.frame(Tumor_Sample_Barcode = cli$Case_ID[index], 
                                         Metabolic_glyco = Metabolic_glyco_gsva[t_strings][index],
                                         Immune_glyco = immune_glyco_gsva[t_strings][index]), 
                     by = "Tumor_Sample_Barcode", all.x = TRUE)
  }
  
  if (!is.null(Metabolic_ubiquity_gsva) && !is.null(immune_ubiquity_gsva)) {
    index <- match(rownames(Metabolic_acety_gsva), cli$`specimen/aliquout_id_protein_tumor`) %>% na.omit()
    cligsva <- merge(cligsva, data.frame(Tumor_Sample_Barcode = cli$Case_ID[index], 
                                         Metabolic_ubiquity = Metabolic_ubiquity_gsva[index],
                                         Immune_ubiquity = immune_ubiquity_gsva[index]), 
                     by = "Tumor_Sample_Barcode", all.x = TRUE)
  }
  # Sort by both Metabolic and Immune scores and create new group columns
  mudata@clinical.data <- cligsva #%>% arrange(desc(Metabolic_phos), desc(Immune_phos)) 
  # Filter PTM features that actually exist in the clinical data
  PTM_features <-  c('Metabolic_phos', 'Metabolic_acety', 'Metabolic_glyco', 'Metabolic_ubiquity',
                     'Immune_phos', 'Immune_acety', 'Immune_glyco', 'Immune_ubiquity')
  valid_PTM_features <- PTM_features[PTM_features %in% colnames(mudata@clinical.data)]
  # Define the list of possible PTM features and their respective grouping rules
  PTM_features_groups <- list(
    Metabolic_phos = "Metabolic_phosphoproteome",
    Metabolic_acety = "Metabolic_acetylome",
    Metabolic_glyco = "Metabolic_glycoproteome",
    Metabolic_ubiquity = "Metabolic_ubiquitylome",
    Immune_phos = "Immune_phosphoproteome",
    Immune_acety = "Immune_acetylome",
    Immune_glyco = "Immune_glycoproteome",
    Immune_ubiquity = "Immune_ubiquitylome"
  )
  # Check which features in PTM_features_groups exist in the clinical data
  valid_PTM_features <- names(PTM_features_groups)[names(PTM_features_groups) %in% colnames(mudata@clinical.data)]
  # Apply mutate using across only on the valid PTM features
  mudata@clinical.data <- mudata@clinical.data %>%
    mutate(across(all_of(valid_PTM_features), ~ ifelse(. > median(., na.rm = TRUE), "High", "Low"), 
                  .names = "{PTM_features_groups[.col]}"))  #%>% arrange(desc(Metabolic_acety))
  # Get valid PTM features based on those that exist in clinical data
  valid_PTM_features_group <- PTM_features_groups[names(PTM_features_groups) %in% colnames(mudata@clinical.data)]
  # Define color mappings for each valid feature group
  color_map <- list(
    Metabolic_phosphoproteome = col1,#c('High' = 'darkred', 'Low' = 'lightpink'),
    Metabolic_acetylome = col1,#c('High' = 'darkgreen', 'Low' = 'lightgreen'),
    Metabolic_glycoproteome = col1,#c('High' = 'purple', 'Low' = 'lavender'),
    Metabolic_ubiquitylome = col1,#c('High' = 'navy', 'Low' = 'skyblue'),
    Immune_phosphoproteome = col2,#c('High' = 'darkblue', 'Low' = 'lightblue'),
    Immune_acetylome = col2,#c('High' = 'darkorange', 'Low' = 'lightyellow'),
    Immune_glycoproteome = col2,#c('High' = 'brown', 'Low' = 'beige'),
    Immune_ubiquitylome = col2#c('High' = 'darkgoldenrod', 'Low' = 'wheat')
  )
  # Filter the color map for only valid PTM feature groups
  annotationColor <- list()
  for (feature in valid_PTM_features_group) {
    if (feature %in% names(color_map)) {
      annotationColor[[feature]] <- color_map[[feature]]
    }
  }
  # Generate Oncoplots with Expanded Clinical Features
  genesummary <- getGeneSummary(mudata)
  index = match(RBPgene, genesummary$Hugo_Symbol)
  mutagene <- genesummary[na.omit(index), 'Hugo_Symbol']
  num_plots <- ceiling(length(mutagene$Hugo_Symbol) / 30)
  print(paste('figure numbers: ', num_plots))
  for (j in 1:num_plots) {
    start_index <- (j - 1) * 25 + 1
    end_index <- min(j * 25, length(mutagene$Hugo_Symbol))
    genes_subset <- mutagene[start_index:end_index]
     # Output Oncoplot as PNG
    png_filename <- paste0("/data2/rluo4/EpiTrans/PTM/CPTAC/fig/mutation/", tumor, "/RBP", start_index, "_to_", end_index, ".png")
    # Clinical Features now include both metabolic and immune pathway scores
    # oncoplot(maf = mudata, genes = genes_subset$Hugo_Symbol,
    #          clinicalFeatures =  names(annotationColor),
    #          sortByAnnotation = TRUE, showTumorSampleBarcodes = TRUE,
    #          draw_titv = TRUE, top = 30, fontSize = 0.75, SampleNamefontSize = 0.7, legend_height = 5, gene_mar = 7)
    # mudata@clinical.data <- cligsva %>% arrange(desc(Metabolic_phos), desc(Immune_phos)) %>%  # Sort by both Metabolic and Immune scores
    # mutate(Metabolic_phos_group = ifelse(Metabolic_phos > median(Metabolic_phos), 'High', 'Low'),
    #        Immune_phos_group = ifelse(Immune_phos > median(Immune_phos), 'High', 'Low')
    # )
    # Generate the oncoplot with the updated group features and colors
    if(tumor == 'LSCC'){
      png(png_filename, width = 13, height = 10, res = 500, units = 'in')
      oncoplot(maf = mudata,
               genes = genes_subset$Hugo_Symbol,
               clinicalFeatures = names(annotationColor),  # Only valid PTM group features will be used
               colors = plot_cols, draw_titv = TRUE,
               anno_height = 2,  # Annotation bar height
               bgCol = "#dbdbdb",  # Background color
               annoBorderCol = "white",  # Border color for annotations
               annotationColor = annotationColor,  # Dynamic annotation colors
               legendFontSize = 1,
               annotationFontSize = 1,
               legend_height = 2,
               sortByAnnotation = TRUE)
    } else{
      png(png_filename, width = 12, height = 10, res = 500, units = 'in')
      oncoplot(maf = mudata,
               genes = genes_subset$Hugo_Symbol,
               clinicalFeatures = names(annotationColor),  # Only valid PTM group features will be used
               colors = plot_cols, draw_titv = TRUE,
               anno_height = 2,  # Annotation bar height
               bgCol = "#dbdbdb",  # Background color
               annoBorderCol = "white",  # Border color for annotations
               annotationColor = annotationColor,  # Dynamic annotation colors
               legendFontSize = 1.2, 
               annotationFontSize = 1.1,
               legend_height = 2,
               sortByAnnotation = TRUE)
    }


    
    dev.off()  # Close the graphics device
  }
}




# # for (i in grep("OV", names(phossymbolistT), invert = TRUE, value = TRUE)) {
# for (tumor in names(mulistT)) {
#   # 
#   pathways_to_plot <- my_genesets[my_genesets$term %in% metabolic_pathways_to_plot,]
#   data_PTM <- vector("list", 4)
#   names(data_PTM) <- c('phos', 'acety', 'glyco', 'ubiquity')
#   if( tumor %in% names(phossymbolistT) ){
#     data_exp <- phossymbolistT[[tumor]] %>% t() %>% as.matrix(); class(data_exp) <- "numeric"
#     phos_protein <- data.frame(phos = rownames(data_exp), gene = str_split(rownames(data_exp), ' ', simplify = T)[,1])
#     pathways_to_plot <- left_join(pathways_to_plot, phos_protein, by = 'gene')
#     data_PTM[['phos']] <- data_exp
#   }
#   if( tumor %in% names(acetylist) ){
#     data_exp <- acetylist[[tumor]] %>% t() %>% as.matrix();class(data_exp) <- "numeric"
#     acety_protein <- data.frame(acety = rownames(data_exp), gene = str_split(rownames(data_exp), ' ', simplify = T)[,1])
#     pathways_to_plot <- left_join(pathways_to_plot, acety_protein, by = 'gene')
#     data_PTM[['acety']] <- data_exp
#   }
#   if( tumor %in% names(glycolist) ){
#     data_exp <- glycolist[[tumor]] %>% t() %>% as.matrix();class(data_exp) <- "numeric"
#     glyco_protein <- data.frame(glyco = rownames(data_exp), gene = str_split(rownames(data_exp), ' ', simplify = T)[,1])
#     pathways_to_plot <- left_join(pathways_to_plot, glyco_protein, by = 'gene')
#     data_PTM[['glyco']] <- data_exp
#   }
#   if( tumor %in% names(ubiquitylist) ){
#     data_exp <- ubiquitylist[[tumor]] %>% t() %>% as.matrix();class(data_exp) <- "numeric"
#     ubiquity_protein <- data.frame(ubiquity = rownames(data_exp), gene = str_split(rownames(data_exp), ' ', simplify = T)[,1])
#     pathways_to_plot <- left_join(pathways_to_plot, ubiquity_protein, by = 'gene')
#     data_PTM[['ubiquity']] <- data_exp
#   }
#    
#   # RBPphoset <- rownames(data_exp)[grep(paste(RBPgene, collapse="|"), rownames(data_exp))]
#   # RBPphoset <- list(RBPphoset = RBPphoset)
#   Metabolic.gene <- na.omit(unique(pathways_to_plot$phos))#unique(dplyr::pull(pathways_to_plot,3))
#   Metabolic_phos <- list(Metabolic_phos = Metabolic.gene)
#   geneSet <- GeneSetCollection(mapply(function(geneIds, keggId){
#     GeneSet(geneIds, geneIdType=EntrezIdentifier(),
#             collectionType=KEGGCollection(keggId),
#             setName=keggId)
#   }, Metabolic_phos, names(Metabolic_phos)))
#   
#   set.seed(123)
#   gsva_res <- gsva(data_PTM[['phos']], geneSet, kcdf=kcdf, method='ssgsea',parallel.sz=32)
#   gsva_res <- gsva_res %>% t() %>% as.data.frame()
#   gsva_res <- gsva_res %>% 
#     mutate("Metabolic_phos" = ifelse(gsva_res$Metabolic_phos > median(gsva_res$Metabolic_phos),'High', 'Low')) %>% 
#     as.data.frame() %>% 
#     rownames_to_column("Tumor_Sample_Barcode") %>% arrange(desc(Metabolic_phos))
#   mudata <- mulistT[[i]]
#   cli <- mudata@clinical.data
#   cligsva <- merge(cli,gsva_res, by = "Tumor_Sample_Barcode") %>% arrange(desc(Metabolic_phos))
#   mudata@clinical.data <- cligsva
#   mutagene <- getGeneSummary(mudata)
#   # mutagene <- genesummary[genesummary$total > 0,"Hugo_Symbol"]
#   # mutagene <- genesummary[genesummary$Hugo_Symbol %in% RBPgene, "Hugo_Symbol"] %>% dplyr::pull(.,1)
#   index = match(RBPgene, genesummary$Hugo_Symbol)
#   mutagene <- genesummary[ na.omit(index), 'Hugo_Symbol']
#   num_plots <- ceiling(length(mutagene) / 30)
#   for (j in 1:num_plots) {
#     start_index <- (j - 1) * 25 + 1
#     end_index <- min(j * 25, length(mutagene))
#     genes_subset <- mutagene[start_index:end_index]
#     # genes_subset <- RMP_update[match(genes_subset, RMP_update$RMP),]
#     png_filename <- paste("/data2/rluo4/EpiTrans/PTM/CPTAC/fig/mutation/",i,"/","RBP",start_index, "_to_", end_index,".png",sep = "")
#     png(file = png_filename, width = 800, height = 700)  
#     oncoplot(maf = mudata,genes = genes_subset$RMP,
#              clinicalFeatures = c('Metabolic_phos', 'Metabolic_acety', 'Metabolic_glyco', 'Metabolic_ubiquity'), 
#              colors = plot_cols,anno_height = 1,# 设置临床信息注释条的高度
#              bgCol = "#dbdbdb",annoBorderCol = "white",
#              annotationColor = list(group = col1),
#              sortByAnnotation = TRUE)
#     dev.off()
#   }
# }
