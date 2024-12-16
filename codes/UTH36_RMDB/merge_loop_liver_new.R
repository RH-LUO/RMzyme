#!/usr/local/bin/Rscript
# title: "RMzyme: Regulations of RNA Modifying Enzymes in Human"
# title:
# author: "Ruihan Luo"
# date: "July 15th,2024"
# rm(list=ls())
path = .libPaths(); .libPaths(path[c(6,2:5)])# (env_nf) rluo4@sbcsmadlp005:/data/rluo4/bin/miniconda3/envs/env_nf/lib/R/bin$ pwd
# path <- c('/home/rluo4/R/x86_64-pc-linux-gnu-library/4.3', '/opt/R/4.3.1/lib64/R/library', '/home/rluo4/R/x86_64-conda-linux-gnu-library/4.3', '/data/rluo4/bin/miniconda3/lib/R/library')
# path = c( "/home/rluo4/R/x86_64-pc-linux-gnu-library/4.1" ,  "/usr/local/lib/R/site-library"  ,   "/usr/local/lib/R/library" ,    "/data/rluo4/lorihan/R/site-library" , "/data/rluo4/lorihan/R/x86_64-pc-linux-gnu-library/4.1")
# .libPaths(path)
library(ggplot2)
library(Seurat)
library(tidyverse)
library(NMF)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)
options(stringsAsFactors = F)
# rm(list = ls())
################### 需要修改绝对路径为当前工作路径
# 1. 读入当前所有数据集编号
loop_path <- "/data/rluo4/EpiTrans/DataCollection/liver_NEW"#loop_path <- "/home/jjguo/xuhaixia/RPMfunc/disco/adjust_disco_split_tissue/skin"
subfolders <- list.dirs(path = loop_path, full.names = TRUE, recursive = FALSE)
subfolder_names <- basename(subfolders); subfolder_names <- subfolder_names[!subfolder_names %in% c('healthy_split', '.ipynb_checkpoints')]
print(subfolder_names)
# 2. 循环数据集
for (j in 1:length(subfolder_names) ){#1:length(subfolder_names)) {
  print(subfolder_names[j])
  # obj/NK
  analysis_path <- file.path(loop_path,subfolder_names[j])
  healthy_path <- paste0(loop_path, "/healthy_split"); setwd(analysis_path)
  # obj <- readRDS(paste0(subfolder_names[j],"_SCT.rds"))#readRDS(paste0(analysis_path,".rds"))
  dataset_tissue <- "liver_NEW"
  dataset_name <- subfolder_names[j]
  # 检查有没有细胞数小于30超过15000的
  # low_count_types <- names(table(obj$sub_type)[table(obj$sub_type) < 30])
  # high_count_types <- names(table(obj$sub_type)[table(obj$sub_type) > 15000])
  # # 从原始数据中删除这些细胞类型
  # obj <- obj[, !(obj$sub_type %in% low_count_types)]
  # obj <- obj[, !(obj$sub_type %in% high_count_types)]
  # # 打印过滤后的细胞类型表格
  # print(table(obj$sub_type))
  # 
  # disease_list <- SplitObject(object = obj, split.by = "sub_type")
  # print(names(disease_list))
  disease_list <- list.files(file.path(analysis_path, 'disease_split'))
  names(disease_list) <- gsub('.rds','',disease_list) 
  for (i in 1:length(names(disease_list))) {
    analysis_celltype <- names(disease_list)[[i]]
    save_path=paste0(file.path(analysis_path,analysis_celltype))
    if( dir.exists(save_path) ){
      print(paste0(analysis_celltype, " already done !"))
      next;
    }else{
      dir.create(save_path)
    }
    setwd(save_path)
    print(save_path)
    print("————————————————————————————start————————————————————————————")
    print(analysis_celltype)
    
    analysis_celltype_gse <- readRDS(paste0(file.path(analysis_path,"disease_split",analysis_celltype),".rds"))
    healthy_analysis_celltype_gse <- readRDS(paste0(file.path(healthy_path,analysis_celltype),".rds"))
    
    # 合并disease的Tcells和healthy的Tcells
    merge_seurat <- merge(analysis_celltype_gse,healthy_analysis_celltype_gse)
    analysis_celltype_gse <- merge_seurat
    
    print(paste0("————————————————", dataset_name, analysis_celltype, "——————————————————"))
    print(analysis_celltype_gse)
    analysis_celltype_gse <- subset(analysis_celltype_gse, subset = nFeature_RNA > 200 & nFeature_RNA < 7000)
    # 过滤线粒体
    analysis_celltype_gse[["percent.mt"]] <- PercentageFeatureSet(analysis_celltype_gse, pattern = "^MT-")
    analysis_celltype_gse <- subset(analysis_celltype_gse, subset = percent.mt < 10)
    print(analysis_celltype_gse)
    if( dim(analysis_celltype_gse)[2] > 16000 ){
      print(paste0(analysis_celltype, " cannot be coped with in UTH36, skip !"))
      next;
    }
    # SCTransform
    analysis_celltype_gse <- SCTransform(analysis_celltype_gse, return.only.var.genes = TRUE,variable.features.n = 2000)
    # 保存SCT后的Tcells
    saveRDS(analysis_celltype_gse,paste0(file.path(save_path,analysis_celltype),"_gse.rds"))
    
    start_time <- Sys.time()
    print(paste("start NMF:", subfolder_names[j], '--', analysis_celltype, ' at', start_time)) #print(paste0("————————————————START", analysis_celltype, ":NMF ——————————————————"))
    # NMF
    range = 4:12 # 迭代次数
    gmin = 5 #最低关联数量
    ncores = 50 #运行核数量
    
    # data <- GetAssayData(analysis_celltype_gse, layer = 'scale.data') %>% as.matrix()
    highly_variable_genes <- VariableFeatures(analysis_celltype_gse)
    highly_variable_genes_matrix <- GetAssayData(analysis_celltype_gse, layer = 'scale.data') %>% as.matrix()
    dim(highly_variable_genes_matrix)
    class(highly_variable_genes_matrix)
    
    valid_genes <- highly_variable_genes[highly_variable_genes %in% rownames(highly_variable_genes_matrix)]
    data <- highly_variable_genes_matrix[valid_genes, ]
    print(dim(data))
    
    colnames(data) <- colnames(analysis_celltype_gse)
    # rownames(data) <- VariableFeatures(analysis_celltype_gse)
    rownames(data) <- valid_genes
    print(dim(data))
    data[data < 0] <- 0
    data <- data[apply(data, 1, var) > 0, ]
    res.list <- parallel::mclapply(range, function(r){
      result <- nmf(data, rank = r, nrun = 1, seed = 'ica', method = 'nsNMF')
      return(result)
    }, mc.cores = ncores)
    names(res.list) = range
    
    saveRDS(res.list,paste0(file.path(save_path,analysis_celltype),"_gse_res.list.rds"))
    
    # NMFToModules function
    NMFToModules = function(
      res,
      gmin = 5
    ){
      scores = basis(res)
      coefs = coefficients(res)
      
      # Remove if fewer than gmin genes
      ranks_x = t(apply(-t(t(scores) / apply(scores, 2, mean)), 1, rank))
      ranks_y = apply(-t(t(scores) / apply(scores, 2, mean)), 2, rank)
      for (i in 1:ncol(scores)){
        ranks_y[ranks_x[,i] > 1,i] = Inf
      }
      modules = apply(ranks_y, 2, function(m){
        a = sort(m[is.finite(m)])
        a = a[a == 1:length(a)]
        names(a)
      })
      l = sapply(modules, length)
      keep = (l >= gmin)
      scores = scores[, keep]
      coefs = coefs[keep, ]
      
      # Find modules
      ranks_x = t(apply(-t(t(scores) / apply(scores, 2, mean)), 1, rank))
      ranks_y = apply(-t(t(scores) / apply(scores, 2, mean)), 2, rank)
      for (i in 1:ncol(scores)){
        ranks_y[ranks_x[,i] > 1,i] = Inf
      }
      modules = apply(ranks_y, 2, function(m){
        a = sort(m[is.finite(m)])
        a = a[a == 1:length(a)]
        names(a)
      })
      
      names(modules) = sapply(modules, '[', 1)
      # names(modules) = seq_along(modules)
      names(modules) = paste('m', names(modules), sep = '_')
      names(modules) = gsub('-','_',names(modules))
      
      return(modules)
    }
    
    # GSEXXXXXX.nmf.rds
    modules.list = lapply(res.list, NMFToModules, gmin = gmin)
    print(sapply(modules.list,length))
    comp = as.numeric(names(modules.list)) - sapply(modules.list, length)
    mi = min(comp)
    r = names(which(comp == mi))
    r = r[length(r)]
    res = res.list[[r]]
    modules = NMFToModules(res, gmin = gmin)
    print(names(modules))
    
    analysis_celltype_gse <- RunPCA(analysis_celltype_gse)
    analysis_celltype_gse@reductions$nmf <- analysis_celltype_gse@reductions$pca
    analysis_celltype_gse@reductions$nmf@cell.embeddings <- t(coef(res))    
    analysis_celltype_gse@reductions$nmf@feature.loadings <- basis(res) 
    
    analysis_celltype_gse.nmf <- RunUMAP(analysis_celltype_gse, reduction = 'nmf', dims = 1:r) %>%
      FindNeighbors(reduction = 'nmf', dims = 1:r) %>% FindClusters()
    
    end_time <- Sys.time()
    print(end_time)
    time_taken <- end_time - start_time
    print(paste("Time taken to process NMF on", subfolder_names[j], '--', dim(analysis_celltype_gse.nmf)[2], analysis_celltype,  'as below:'))
    print(time_taken)
    
    
    print(paste0("————————————————START", analysis_celltype, ":RMP+sub_type ——————————————————"))
    RMP <- readxl::read_excel('/data/rluo4/EpiTrans/DataCollection/RMP_update.xlsx')#read.csv("/home/jjguo/xuhaixia/RPMfunc/RMP_12types.csv",header = T)
    # add RMP+celltype
    # 先FindAllMarkers计算每个seurat_clusters的marker
    Idents(analysis_celltype_gse.nmf) <- analysis_celltype_gse.nmf$seurat_clusters
    # 不同簇之间文库大小不一致
    analysis_celltype_gse.nmf <- PrepSCTFindMarkers(analysis_celltype_gse.nmf)
    markers <- FindAllMarkers(analysis_celltype_gse.nmf,
                              only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.25, 
                              max.cells.per.ident = 250)
    
    matching_markers <- markers[markers$gene %in% RMP$RMP, ]
    grouped_markers <- split(matching_markers, matching_markers$cluster)
    # 定义一个函数,找到logfc最大的RMP
    find_max_abs_logFC <- function(df) {
      max_idx <- which.max(abs(df$avg_log2FC))
      return(df[max_idx, ])
    }
    max_logFC_genes <- do.call(rbind, lapply(grouped_markers, find_max_abs_logFC))
    print(max_logFC_genes)
    # 将 cluster 列设置为行名
    rownames(max_logFC_genes) <- max_logFC_genes$cluster
    name_to_RMP <- data.frame(cluster_name = character(),
                              RMP = character(),
                              stringsAsFactors = FALSE)
    # 定义一个函数，通过 cluster_name 获取对应的 gene
    get_gene_by_cluster <- function(cluster_name) {
      if (cluster_name %in% rownames(max_logFC_genes)) {
        return(max_logFC_genes[as.character(cluster_name), "gene"])
      } else {
        return(NA)
      }
    }
    
    # 获取每个cluster中的RMP
    test_list <- list()
    for (cluster_id in unique(markers$cluster)) {
      nmf_names <- subset(markers, cluster == cluster_id)
      nmf_names <- rownames(nmf_names)
      test_list[[as.character(cluster_id)]] <- nmf_names[nmf_names %in% RMP$RMP]
    }
    # 获取每个cluster中的logfc绝对值最高的RMP
    for (cluster_name in names(test_list)) {
      first_gene <- test_list[[cluster_name]][1]
      if (length(first_gene) > 0) {
        new_name <- get_gene_by_cluster(cluster_name)
      } else {
        new_name <- "NoneRMP"
      }
      name_to_RMP <- rbind(name_to_RMP, data.frame(cluster_name = cluster_name, RMP = new_name))
    }
    name_to_RMP$RMP[is.na(name_to_RMP$RMP)] <- "NoneRMP"
    print(name_to_RMP)
    new_cluster_names <- name_to_RMP$RMP[match(as.character(analysis_celltype_gse.nmf$seurat_clusters), 
                                               as.character(name_to_RMP$cluster_name))]
    analysis_celltype_gse.nmf$seurat_clusters <- paste(new_cluster_names, "+", analysis_celltype, sep = "")
    # 保存analysis_celltype_gse.nmf
    saveRDS(analysis_celltype_gse.nmf,paste0(file.path(save_path,analysis_celltype),"_gse.nmf.rds"))
    
    print(paste0("————————————————START", analysis_celltype, ":RMP_boxplot stop ——————————————————"))
    # ################ 收集特有RMP,也就是用来画box的RMP,也就是这个GSE用来搜索的RMP
    # RMP_plot <- max_logFC_genes$gene
    # # 写入GSEXXXX对应的RMP
    # max_logFC_genes$tissue <- dataset_tissue # 添加tissue
    # max_logFC_genes$dataset <- dataset_name
    # last_three_columns <- max_logFC_genes[, c("dataset", "gene","tissue")]
    # print(head(last_three_columns))
    # RMP_file_path <- file.path(save_path, analysis_celltype, "boxplot_RMP.txt")
    # # 检查路径
    # dir_path <- dirname(RMP_file_path)
    # if (!dir.exists(dir_path)) {
    #   dir.create(dir_path, recursive = TRUE)
    # }
    # write.table(last_three_columns, file = RMP_file_path, sep = "\t",
    #             quote = FALSE, row.names = FALSE)
    # print(head(last_three_columns))
    # rawCounts <- GetAssayData(object = analysis_celltype_gse, layer = "counts")
    # rmp_counts <- as(log2(t(t(rawCounts)/Matrix::colSums(rawCounts)) * 10000 + 1), "matrix")
    # print(dim(rmp_counts))
    # 
    # ##### 如果只有一个RMP+sub_type
    # if(class(rmp_counts[RMP_plot, ])[1] != "numeric") {
    #   sub_rmp_counts <- as.data.frame(t(rmp_counts[RMP_plot, ]))
    # } else {    
    #   sub_rmp_counts <- as.data.frame(rmp_counts[RMP_plot, ])
    #   print("只有一个RMP+sub_type")
    #   
    # }
    # 
    # colnames(sub_rmp_counts) <- RMP_plot
    # print(dim(sub_rmp_counts))
    # print(head(sub_rmp_counts))
    # for (p in colnames(sub_rmp_counts)) {
    #   RMP_data <- as.data.frame(sub_rmp_counts[, p])
    #   rownames(RMP_data) <- rownames(sub_rmp_counts)
    #   colnames(RMP_data) <- p
    #   
    #   RMP_data$disease <- analysis_celltype_gse$disease
    #   
    #   RMP_boxplot <-
    #     ggplot(RMP_data, aes(x = disease, y = !!sym(colnames(RMP_data)[1]))) +
    #     stat_boxplot(geom="errorbar")+
    #     geom_boxplot(aes(fill=disease),lwd = 0.1) +
    #     labs(x = "Disease", y = paste(colnames(RMP_data)[1], "Expression")) +
    #     theme_bw()
    #   ggsave(paste0(file.path(save_path,colnames(RMP_data)[1]),"_RMP_boxplot.png"),
    #          RMP_boxplot,
    #          width = 5, height = 5)
    # }
    # 
    # # 计算nono_celltype 与 RMP_celltype的DEGs
    # Idents(analysis_celltype_gse.nmf) <- analysis_celltype_gse.nmf$seurat_clusters
    # RMP_DEGs <- FindAllMarkers(analysis_celltype_gse.nmf, 
    #                            only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.25, 
    #                            max.cells.per.ident = 250)
    # write.csv(RMP_DEGs,paste0(file.path(save_path,analysis_celltype),"_gse.nmf_RMP_celltype_DEGs.csv")
    #           ,row.names = F)
    # saveRDS(RMP_DEGs,paste0(file.path(save_path,analysis_celltype),"_gse.nmf_RMP_celltype_DEGs.rds"))
    # 
    # # 画analysis_celltype_gse.nmf umap
    # celltype_gse.nmf_umap <- DimPlot(analysis_celltype_gse.nmf, reduction = 'umap',pt.size = 0.3,
    #                                  group.by = 'seurat_clusters') + 
    #   ggtitle("NMF cluster")
    # ggsave(paste0(file.path(save_path,analysis_celltype),"_gse.nmf_UMAP.png"), 
    #        celltype_gse.nmf_umap, width = 6, height = 6)
    # 
    # # 计算healthy_disease_DEGs
    # Idents(analysis_celltype_gse.nmf) <- analysis_celltype_gse.nmf$disease
    # DEGs <- FindMarkers(analysis_celltype_gse.nmf, 
    #                     ident.1 = names(table(Idents(analysis_celltype_gse.nmf)))[1], 
    #                     ident.2 = names(table(Idents(analysis_celltype_gse.nmf)))[2],
    #                     only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.25
    # )
    # DEGs$gene <- rownames(DEGs)
    # write.csv(DEGs,paste0(file.path(save_path,analysis_celltype),"_gse.nmf_healthy_disease_DEGs.csv")
    #           ,row.names = F)
    # saveRDS(DEGs,paste0(file.path(save_path,analysis_celltype),"_gse.nmf_healthy_disease_DEGs.rds"))
    # 
    # ########## healthy vs disease nmf proportion
    # library(dittoSeq)
    # head(analysis_celltype_gse.nmf@meta.data)
    # ditto_plot <- dittoBarPlot(
    #   object = analysis_celltype_gse.nmf,
    #   var = "seurat_clusters",
    #   group.by = "disease")
    # ggsave(paste0(file.path(save_path,analysis_celltype),"_healthy_vs_disease_proportion.png"), 
    #        ditto_plot, width = 6, height = 6)
  }
  
}
