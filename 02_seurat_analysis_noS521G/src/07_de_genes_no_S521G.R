#!/usr/bin/env Rscript

#load required packages
library(Seurat)
library(dplyr)
library(ggplot2)
library(gridExtra)

#set all directory paths
basedir <- paste0("/gladstone/bioinformatics/projects/",
                  "mn-1318-maxine-nelson-yadong-huang-snrnaseq-mm10-aug-2022/",
                  "results/data")
indir <- paste0(basedir,
                "/02_clustering_no_S521G/")
outdir <- paste0(basedir,
                 "/07_de_genes_no_S521G/")
setwd(outdir)

#load the Seurat object
pca_dim <- 15
cluster_res <- 0.7
dat <- readRDS(paste0(indir,
                      "sct_data_no_S521G_post_cluster_pcadims_",
                      pca_dim,
                      "_res_",
                      cluster_res,
                      ".rds"))

#1. Clusters 3 is likely degenerating oligodendrocytes, which are reduced in...
# PS19-E4-R/S and PS19-E4-S/S mice. Perform DE genes analysis of cluster 3 vs 2.
de_cluster_3_vs_2 <- FindMarkers(object = dat,
                                 ident.1 = 3,
                                 ident.2 = 2,
                                 assay = "SCT",
                                 slot = "data",
                                 test.use = "wilcox",
                                 logfc.threshold = 0.1)
de_cluster_3_vs_2 <- cbind(gene=rownames(de_cluster_3_vs_2),
                           de_cluster_3_vs_2)
write.csv(de_cluster_3_vs_2,
          file = paste0("de_genes_cluster_3_vs_2_",
                        "post_clustering_sct_no_S521G_pcadims_",
                        pca_dim,
                        "_res_",
                        cluster_res,
                        ".csv"),
          row.names = FALSE)

#2. Clusters 9 is likely degenerating oligodendrocytes, which are reduced in...
# PS19-E4-R/S and PS19-E4-S/S mice. Perform DE genes analysis of cluster 9 vs 2.
de_cluster_9_vs_2 <- FindMarkers(object = dat,
                                 ident.1 = 9,
                                 ident.2 = 2,
                                 assay = "SCT",
                                 slot = "data",
                                 test.use = "wilcox",
                                 logfc.threshold = 0.1)
de_cluster_9_vs_2 <- cbind(gene=rownames(de_cluster_9_vs_2),
                           de_cluster_9_vs_2)
write.csv(de_cluster_9_vs_2,
          file = paste0("de_genes_cluster_9_vs_2_",
                        "post_clustering_sct_no_S521G_pcadims_",
                        pca_dim,
                        "_res_",
                        cluster_res,
                        ".csv"),
          row.names = FALSE)

#3. Clusters 1, 4, 6, 7, 8, 22, 26, 28 are likely normal excitatory neurons,...
# which are increased in PS19-E4-R/S and PS19-E4-S/S mice. Therefore, we want...
# to know the DE genes for each of the clusters, comparing...
# PS19-E4-R/S vs PS19-E4 and PS19-E4-S/S vs PS19-E4
for(cluster_of_interest in c(1, 4, 6, 7, 8, 22, 26, 28)){
  de_cluster_rs_vs_e4 <- FindMarkers(object = dat,
                                     ident.1 = "fE4R_S",
                                     ident.2 = "fE4",
                                     verbose = TRUE,
                                     group.by="genotype",
                                     subset.ident = cluster_of_interest,
                                     assay = "SCT",
                                     slot = "data",
                                     test.use = "wilcox",
                                     logfc.threshold = 0.1)
  de_cluster_ss_vs_e4 <- FindMarkers(object = dat,
                                     ident.1 = "fE4S_S",
                                     ident.2 = "fE4",
                                     verbose = TRUE,
                                     group.by="genotype",
                                     subset.ident = cluster_of_interest,
                                     assay = "SCT",
                                     slot = "data",
                                     test.use = "wilcox",
                                     logfc.threshold = 0.1)
  #add gene symbols as a column
  de_cluster_rs_vs_e4 <- cbind(gene=rownames(de_cluster_rs_vs_e4),
                               de_cluster_rs_vs_e4)
  de_cluster_ss_vs_e4 <- cbind(gene=rownames(de_cluster_ss_vs_e4),
                               de_cluster_ss_vs_e4)

  #write out the results to a csv file
  write.csv(de_cluster_rs_vs_e4,
            file = paste0("de_genes_cluster_",
                          cluster_of_interest,
                          "_PS19-E4-R_S_vs_PS19-E4_",
                          "post_clustering_sct_no_S521G_pcadims_",
                          pca_dim,
                          "_res_",
                          cluster_res,
                          ".csv"),
            row.names = FALSE)
  write.csv(de_cluster_ss_vs_e4,
            file = paste0("de_genes_cluster_",
                          cluster_of_interest,
                          "_PS19-E4-S_S_vs_PS19-E4_",
                          "post_clustering_sct_no_S521G_pcadims_",
                          pca_dim,
                          "_res_",
                          cluster_res,
                          ".csv"),
            row.names = FALSE)
}

#4. DE genes of PS19-E4-S/S versus PS19-E4 for cluster 9
de_cluster_9_ss_vs_e4 <- FindMarkers(object = dat,
                                     ident.1 = "fE4S_S", 
                                     ident.2 = "fE4",
                                     verbose = TRUE, 
                                     group.by="genotype",
                                     subset.ident = 9,
                                     assay = "SCT", 
                                     slot = "data", 
                                     test.use = "wilcox",
                                     logfc.threshold = 0.1)
de_cluster_9_ss_vs_e4 <- cbind(gene=rownames(de_cluster_9_ss_vs_e4),
                               de_cluster_9_ss_vs_e4)
write.csv(de_cluster_9_ss_vs_e4,
          file = paste0("de_genes_cluster_9",
                        "_PS19-E4-S_S_vs_PS19-E4_",
                        "post_clustering_sct_no_S521G_pcadims_",
                        pca_dim,
                        "_res_",
                        cluster_res,
                        ".csv"),
          row.names = FALSE)


#5. find the list of background genes
all_data <- GetAssayData(dat, assay = "SCT")
background_genes <- rownames(all_data[rowSums(all_data[])>0,])
write.csv(as.data.frame(background_genes),
          file = paste0("nonzero_bakcground_gene_list_",
                        "post_clustering_sct_no_S521G_pcadims_",
                        pca_dim,
                        "_res_",
                        cluster_res,
                        ".csv"),
          row.names = FALSE)


print("********** Script completed! **********")


################## END ################## 
