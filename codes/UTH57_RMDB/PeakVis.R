#!/usr/local/bin/Rscript
# title: "RBM33 is a unique m(6)A RNA-binding protein that regulates ALKBH5 demethylase activity and substrate selectivity"
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
# 1) load the rdata from RMdeg.R
load(paste0('/data2/rluo4/EpiTrans/RMDatasets/allDat.RData'))
load(paste0('/data2/rluo4/EpiTrans/RMDatasets/RMP_allDat.RData'))
load(paste0('/data2/rluo4/EpiTrans/RMDatasets/RMdeg_allDat.RData'))
RMP_update <- read_excel('/data2/rluo4/EpiTrans/DataCollection/RMP_update.xlsx')
asp_sc <- fread(file="/data2/rluo4/EpiTrans/DataCollection/asp_RMP_epitrans.link.txt", header = F)
# 2) load the functions
extract_names <- function(df) {
  colnames(df)
}
# Function to transform the strings
transform_path <- function(x) {
  x <- gsub("^X\\.", "/", x)
  x <- gsub("\\.(?=.*_)", "/", x, perl = TRUE)
  return(x)
}
# # Replace 'X.' with '/'
# modified_string <- gsub("^X\\.", "/", original_string)
# # Replace the last '_' with '.'
# modified_string <- gsub("_(?=.*\\.bam)", ".", modified_string, perl = TRUE)
# Define a function to extract gene_name from the attribute string
extract_gene_name <- function(attribute_string) {
  # Use a regular expression to match the gene_name
  gene_name <- str_extract(attribute_string, "gene_name [^;]+")
  # Remove the "gene_name " prefix to get the actual gene name
  gene_name <- str_replace(gene_name, "gene_name ", "")
  gene_name <- str_replace_all(gene_name, '"', '')
  return(gene_name)
}
# 3) the details about RMPdeg_datasets
genome_version <- data.frame(GSE = RMPdeg_datasets$GSE, HG = 'hg19', skip = '')
hg38_list <- c(9,12,13,15,20,21,22,24,30,34,35,36,41,44,45,46,47,53,58,59,65,84,90,95,103,104,108,112,114)
skip_list <- c(23, 26, 27, 31, 32, 52, 54, 55, 60, 61, 70, 72, 83, 91 )
genome_version$HG[hg38_list] <- 'hg38'
genome_version$skip[skip_list] <- 'NoData'
genome_version$macs2_dir <- genome_version$GSE
genome_version$macs2_dir[genome_version$GSE=='GSE90963'] = 'PRJNA355164'
genome_version$macs2_dir[genome_version$GSE=='GSE121942'] = 'GSE110320'
genome_version$macs2_dir[genome_version$GSE=='GSE103497'] = 'GSE103495'
genome_version$macs2_dir[genome_version$GSE=='GSE122803'] = 'GSE122801'
genome_version$macs2_dir[genome_version$GSE=='GSE120024'] = 'GSE132306'
genome_version$macs2_dir[genome_version$GSE=='GSE128575'] = 'GSE128520'
genome_version$macs2_dir[genome_version$GSE=='GSE202815'] = 'GSE202814'
genome_version$macs2_dir[genome_version$GSE=='GSE223731'] = 'GSE223729'
genome_version$macs2_dir[genome_version$GSE=='GSE90684'] = 'GSE90642'

library(stringr)
library(dplyr)
###################################
#  Cell number-normalized (CNN)   #
#     MeRIP RIP visualization     #
###################################
library(GenomicFeatures)
library(Guitar)
library(ggcoverage)
# # 1) Load GTF
# gtf_v36 <- rtracklayer::import.gff(con = '/data2/rluo4/hg38/gencode.v36.annotation.gtf.gz', format = "gtf")
# gtf_v19 <- rtracklayer::import.gff(con = '/data2/rluo4/hg38/gencode.v19.annotation.gtf.gz', format = "gtf")
load( paste0('/data2/rluo4/EpiTrans/RMDatasets/txdb_Guitar.RData') )
# To add gene annotation, the gtf file should contain gene_type and gene_name attributes in column 9; to add transcript annotation, the gtf file should contain transcript_name attribute in column 9.
library(Gviz)
library(CAGEfightR)
library(rtracklayer)
library("ggcoverage")
library("ggpattern")
library("BSgenome.Hsapiens.UCSC.hg19")
library("BSgenome.Hsapiens.UCSC.hg38")
BSv19 = BSgenome.Hsapiens.UCSC.hg19
BSv38 = BSgenome.Hsapiens.UCSC.hg38
# load the target annotation results #
load(file.path(data_directory,'../target_results.RData'))
original <- c('GSE97419', 'GSE55572', 'GSE87190', 'GSE76414',
              'GSE103497', 'GSE141994', 'GSE145924', 'GSE124509',
              'GSE210867',  'GSE207643', 'GSE171497', 'GSE198643')
bw.rank <- match(original, RMPdeg_datasets$GSE)
final_datasets <- read_excel('/data2/rluo4/RPMfunc/Output/summary/dataset_summary.xlsx')#,sep='\t',quote = F,row.names = F)
# RMPdeg_datasets <- final_datasets #106 datasets
# for (n in 1:30) { #/data2/rluo4/Rcode/RMDB/Rs_outs/PeakVis-Q1st.R.o
# for (n in 31:60) { #/data2/rluo4/Rcode/RMDB/Rs_outs/PeakVis-median.R.o
# for (n in 61:90) { #/data2/rluo4/Rcode/RMDB/Rs_outs/PeakVis-Q3rd.R.o
# for (n in 91:120) { #/data2/rluo4/Rcode/RMDB/Rs_outs/PeakVis-max.R.o
# for (n in 13:nrow(RMPdeg_datasets)) { #/data2/rluo4/Rcode/RMDB/Rs_outs/PeakVis-new.R
for (n in 21:nrow(RMPdeg_datasets)) { #/data2/rluo4/Rcode/RMDB/Rs_outs/PeakVis.R
# for (n in rev(1:nrow(RMPdeg_datasets))) { #/data2/rluo4/Rcode/RMDB/Rs_outs/PeakVis-rev.R.o
### for (n in 44:71) { #/data2/rluo4/Rcode/RMDB/Rs_outs/PeakVis-mid.R.o: ATP9B of n=42 PCDH1 of n=43
  if( n %in% bw.rank ){
    print(paste0(i, ": run on UTH36 !"))
    next;
  }
  print(n)
  print(RMPdeg_datasets[n, 'GSE'])
  i <- RMPdeg_datasets[n, 'GSE']
  # result <- process_samples(n, RMPdeg_datasets, data_directory, RMdem_list, RMdep_list, genome_version)
  # target_results[[i]] <- result
  ###################################
  #  Cell number-normalized (CNN)   #
  #     MeRIP RIP Visualization     #
  ###################################
  target_res <- target_results[[i]]
  length(target_res)
  if( length(target_res) ==0 ){
    print(paste0(i, ": no data in target_results !"))
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
  
  for (contrast in names(target_res)) {
    # contrast <- names(target_res)[1]
    start_time <- Sys.time()
    print(paste("start visualizing:", contrast, 'from', i, 'at', start_time))
    
    outdir <- file.path(plot_dir, i, contrast)
    if(! dir.exists( outdir )) {
      dir.create(outdir)
    }
    c <- unique(str_split(contrast, '_', simplify = TRUE)[, 2])
    if( length(cl) >1 | n == 111){
      c <- contrast
    }
    fullname <- paste0(contrast, '_vs_WT')
    # bwfiles <- list.files('./metaPlotR/')
    # bwfiles <- bwfiles[grep('.bw', bwfiles)]
    # bwfiles
    target <- target_res[[contrast]]
    target_genes = unique(target$geneSymbol)
    print(paste0('no. of target peaks: ',length(target_genes)))
    DEG <- resG[[fullname]]
    DEM <- resM[[paste0(fullname, '_DESeq2')]]
    index = grepl('DESeq2', names(resM))
    if (length(unique(index)) > 1 & !is.null(DEM) ) {
      DEM <- resM[[paste0(fullname, '_DESeq2')]]
    } else{
      DEM <- resM[[paste0(fullname, '_limma')]]
    }
    if (length(unique(index)) > 1 & !is.null(DEM) & nrow(DEM) == 0) {
      DEM <- resM[[paste0(fullname, '_limma')]]
    }
    DEM$peak_id <- rownames(DEM)
    DEM <- left_join(target, DEM, by = 'peak_id')
    DEM <- DEM[!is.na(DEM$change),]
    intersect(DEG$geneSymbol, target$geneSymbol)
    print(length(intersect(DEG$geneSymbol, target$geneSymbol)))
    write.table(DEG, file = paste0(outdir, '/',i, '_DEG.txt'), sep = '\t', quote = F)
    write.table(DEM, file = paste0(outdir, '/',i, '_DEM.txt'), sep = '\t', quote = F)
    gtf <- if (genome_version$HG[n] == 'hg19') gtf_v19 else gtf_v36
    # 1) Load track files and metadata
    # Load track files and metadata
    sample_meta <-  sample_meta_InP <-  sample_meta_IP <- NULL
    if (!is.null( contrast_input_list[[i]] )) {
      # RNA-seq:
      bedfiles <- contrast_input_list[[i]][[c]]
      condition <- condition_input_list[[i]][[c]]
      # table(gsub('.bam','.bw',bedfiles) %in% bwfiles)
      
      if (grepl('X.data2', bedfiles[1])) {
        bedfiles <- as.character(sapply(bedfiles, transform_path))
        bedfiles <- str_split(bedfiles, '/', simplify = TRUE)[, 9]
      }
      sample_meta_InP <- data.frame(SampleName = c( gsub('.bam','',bedfiles)),
                                    Type = condition, 
                                    Group = 'Input')
    }
    # RIP-seq:
    if (!is.null(contrast_IP_list[[i]])) {
      bedfiles <- contrast_IP_list[[i]][[c]]
      condition <- condition_IP_list[[i]][[c]]
      if (grepl('X.data2', bedfiles[1])) {
        bedfiles <- as.character(sapply(bedfiles, transform_path))
        bedfiles <- str_split(bedfiles, '/', simplify = TRUE)[, 9]
      }
      sample_meta_IP <- data.frame(SampleName = c( gsub('.bam','',bedfiles)),
                                   Type = condition, 
                                   Group = 'IP')
    }
    # Combine the metadata for RNA-seq and RIP-seq if both are available
    if (exists("sample_meta_InP") & exists("sample_meta_IP")) {
      sample_meta <- rbind(sample_meta_IP, sample_meta_InP)
    } else if (exists("sample_meta_InP")) {
      sample_meta <- sample_meta_InP
    } else if (exists("sample_meta_IP")) {
      sample_meta <- sample_meta_IP
    } else {
      sample_meta <- NULL
    }
    
    # Check and print the combined sample metadata
    # if (!is.null(sample_meta) & g==target_genes[1]){
    if (!is.null(sample_meta)){
      print(sample_meta)
    }
    # format_gse <- 'bw'
    if( i %in% original){
      format_gse <- 'bw'
      track_folder <- file.path(data_directory, i, "metaPlotR")
    } else{
      format_gse <- 'bedgraph'
      i_macs2 <- genome_version$macs2_dir[match(i, genome_version$GSE)]
      track_folder <- file.path(raw_directory, "Out/macs2", i_macs2)
    }
    # 2) Visualize and annotate omics coverage with ggplot2
    target_genes = unique(DEM$geneSymbol)
    print(paste0('no. of DEM peaks: ',length(target_genes)))
    for (g in target_genes) {
      # if(i == 'GSE198643' & g == 'RP1-63G5.8'){
      #   print(paste0(i, ': ', g, " error in LoadTrackFile: TrackFile, value = GSM5954204_M-KD2_Input.bw !"))
      #   next;
      # }
      outdir <- file.path(plot_dir, i, contrast, "PeakVis/")
      if(! dir.exists( outdir )) {
        dir.create(outdir)
      }
      if( file.exists(paste0(outdir, g, '_peak.png')) ){
        print(paste0(i, ': ', g, "_peak.png visualization done !"))
        next;
      }
      # g = intersect(DEG$geneSymbol, target$geneSymbol)[1]
      print(g)
      print(genome_version$HG[n])
      BS <- if (genome_version$HG[n] == 'hg19') BSv19 else BSv38
      
      target.bed <- as.data.frame(gtf[match(g, gtf$gene_name),])
      # region <- paste0(test.bed$V1,":", test.bed$V5)
      # region <- paste0('chr',test.bed$seqnames,":", test.bed$start,"-", test.bed$end)
      # region <- paste0(target_annotation$seqnames, ":", target_annotation$start, "-", target_annotation$end)
      region <- paste0(target.bed$seqnames,":", target.bed$start,"-", target.bed$end)
      # Subset the GTF file for the FBF1 gene
      target_annotation <- subset(gtf, gene_name == g)#subset(gtf, gene_id == "FBF1")
      # # View(as.data.frame(target_annotation))
      # # Identify the transcription start site (TSS) of FBF1 # Assuming FBF1 is a protein-coding gene, you can find the TSS by looking for the "start_codon" annotation
      # target_tss <- start(target_annotation[target_annotation$type == "start_codon"])
      # # Define the promoter region based on the TSS
      # promoter_start <- target_tss - 1000  # 1000 bp upstream of the TSS
      # promoter_end <- target_tss - 1       # 1 bp before the TSS
      # # Extract the sequence of the promoter region from the genome # Assuming you have the genome sequence loaded into the variable 'genome_sequence'
      # # target_promoter_sequence <- genome_sequence[Promoter(start = promoter_start, end = promoter_end)]
      # promoter_df <- data.frame(seqnames = unique(seqnames(target_annotation)),
      #                           start = unique(promoter_start), end = unique(promoter_end), strand = unique(strand(target_annotation)),
      #                           gene_name = i, label = "Promoter")
      # print(promoter_df)
      # mark_region <- promoter_df[, c('start', 'end', 'label')]
      track_df <- LoadTrackFile(
        track.folder = track_folder,
        format = format_gse,
        region = region,#"chr18:76822285-76900000",
        meta.info = sample_meta
      )
      # # create plot
      basic_coverage <- ggcoverage(data = track_df,
                                   range.position = "in",
                                   # plot.type = "joint",
                                   facet.y.scale = "fixed",
                                   #mark.region = mark_region,
                                   show.mark.label = TRUE)
      # Example call for plotting
      if (n == 42 & g=='ATP9B' | n ==43 & g == 'PCDH1') {
        plotx <- basic_coverage +
          # geom_gc(bs.fa.seq = BS) +
          geom_gene(gtf.gr = target_annotation)
      } else {
        plotx <-  basic_coverage +
          geom_gc(bs.fa.seq = BS) +
          geom_gene(gtf.gr = target_annotation)
        # plotx <- plot_ideogram(gene = g, rank = n)
      }
      
      # Check if plotx is not NULL before attempting to save
      if (!is.null(plotx)) {
        ggsave(plot = plotx, file = paste0(outdir, g, '_peak.png'), height = 3.6, width = 5.2, dpi = 500)
      } else {
        message("Plot could not be generated for gene: ", g)
      }
      
      # if(grepl('MT-', g)){
      #   plotx <-  basic_coverage +
      #     geom_gc(bs.fa.seq = BS) +
      #     geom_gene(gtf.gr = target_annotation)
      # } else{
      #   plotx <-  basic_coverage +
      #     geom_gc(bs.fa.seq = BS) +
      #     geom_gene(gtf.gr = target_annotation) +
      #     # geom_peak(bed.file = peak_file) +
      #     geom_ideogram(genome = genome_version$HG[n], plot.space = 0, highlight.centromere = TRUE
      #     ) 
      # }
      # ggsave(plot = plotx, file = paste0(outdir, g, '_peak.png'), height =  3.6, width = 5.2, dpi = 500)
      # Function to plot ideogram with detailed error handling and skipping warnings
      # plot_ideogram <- function(gene, rank) {
      #   tryCatch({
      #     # Check if the rank is within the valid range of genome_version$HG
      #     if(rank > length(genome_version$HG) || rank < 1) {
      #       stop("Invalid rank: out of range.")
      #     }
      #     # Get the genome version based on rank
      #     genome_ver <- genome_version$HG[rank]
      #     # Ensure genome_ver is not NULL or empty
      #     if(is.null(genome_ver) || genome_ver == "") {
      #       stop("Invalid genome version.")
      #     }
      #     # Create the ideogram layer
      #     ideogram_layer <- geom_ideogram(genome = genome_ver, plot.space = 0, highlight.centromere = TRUE)
      #     # Print the ideogram data for debugging
      #     print(ideogram_layer)
      #     # Combine the basic coverage plot with the ideogram
      #     plot <- basic_coverage +
      #       geom_gc(bs.fa.seq = BS) +
      #       geom_gene(gtf.gr = target_annotation) +
      #       ideogram_layer # Added ideogram layer
      #     # Print the combined plot
      #     # print(plot)
      #     # Return the plot object
      #     return(plot)
      #   }, error = function(e) {
      #     message("Error occurred: ", e$message)
      #     return(NULL)  # Return NULL in case of error
      #   })
      # }

      }
  }
}