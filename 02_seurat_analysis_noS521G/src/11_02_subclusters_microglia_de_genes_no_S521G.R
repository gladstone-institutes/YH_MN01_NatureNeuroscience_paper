#!/usr/bin/env Rscript

#load required packages
library(Seurat)
library(dplyr)

#set all directory paths
basedir <- paste0("/gladstone/bioinformatics/projects/",
                  "mn-1318-maxine-nelson-yadong-huang-snrnaseq-mm10-aug-2022/",
                  "results/data")
indir <- paste0(basedir,
                "/06_subcluster_microglia_clusters_17_19/")
outdir <- paste0(basedir,
                 "/11_subclusters_de_genes_no_S521G/microglia/")
setwd(outdir)

#load the Seurat object
pca_dim_micro <- 15
cluster_res_micro <- 0.9
dat_micro <- readRDS(paste0(indir,
                            "microglia_data_post_subcluster_no_S521G_pcadims_",
                            pca_dim_micro,
                            "_res_",
                            cluster_res_micro,
                            ".rds"))

# 1. DE gene for:
## a. subcluster 2 vs other microglia
## b. subcluster 8 vs other microglia
## c. subcluster 11 vs other microglia
## d. For subcluster 2, PS19-E4-R/S vs PS19-E4 and PS19-E4-S/S vs PS19-E4.
## e. For subcluster 8, PS19-E4-R/S vs PS19-E4 and PS19-E4-S/S vs PS19-E4.
## f. For subcluster 11, PS19-E4-R/S vs PS19-E4 and PS19-E4-S/S vs PS19-E4.
clusters_for_de <- c(2,8,11)
other_clusters <- seq(1, max(as.numeric(dat_micro$seurat_clusters)))
other_clusters <- other_clusters[!(other_clusters %in% clusters_for_de)]
print(other_clusters)

for(de_clus in clusters_for_de){
  de_list <- FindMarkers(object = dat_micro, 
                         ident.1 = de_clus,
                         ident.2 = other_clusters,
                         assay = "SCT", 
                         slot = "data", 
                         test.use = "wilcox",
                         logfc.threshold = 0.1)
  
  de_cluster_rs_vs_e4 <- FindMarkers(object = dat_micro,
                                     ident.1 = "fE4R_S", 
                                     ident.2 = "fE4",
                                     verbose = TRUE, 
                                     group.by="genotype",
                                     subset.ident = de_clus,
                                     assay = "SCT", 
                                     slot = "data", 
                                     test.use = "wilcox",
                                     logfc.threshold = 0.1)
  
  de_cluster_ss_vs_e4 <- FindMarkers(object = dat_micro,
                                     ident.1 = "fE4S_S", 
                                     ident.2 = "fE4",
                                     verbose = TRUE, 
                                     group.by="genotype",
                                     subset.ident = de_clus,
                                     assay = "SCT", 
                                     slot = "data", 
                                     test.use = "wilcox",
                                     logfc.threshold = 0.1)
  
  #add gene symbols as a column
  de_list <- cbind(gene=rownames(de_list), 
                   de_list)
  de_cluster_rs_vs_e4 <- cbind(gene=rownames(de_cluster_rs_vs_e4), 
                               de_cluster_rs_vs_e4)
  de_cluster_ss_vs_e4 <- cbind(gene=rownames(de_cluster_ss_vs_e4), 
                               de_cluster_ss_vs_e4)
  
  #write out the results to a csv file
  write.csv(de_list,
            file = paste0("de_genes_microglia_subcluster_",
                          de_clus,
                          "_vs_other_microglia_",
                          "no_S521G_pcadims_",
                          pca_dim_micro,
                          "_res_",
                          cluster_res_micro,
                          ".csv"),
            row.names = FALSE)
  write.csv(de_cluster_rs_vs_e4,
            file = paste0("de_genes_microglia_subcluster_",
                          de_clus,
                          "_PS19-E4-R_S_vs_PS19-E4_",
                          "no_S521G_pcadims_",
                          pca_dim_micro,
                          "_res_",
                          cluster_res_micro,
                          ".csv"),
            row.names = FALSE)
  write.csv(de_cluster_ss_vs_e4,
            file = paste0("de_genes_microglia_subcluster_",
                          de_clus,
                          "_PS19-E4-S_S_vs_PS19-E4_",
                          "no_S521G_pcadims_",
                          pca_dim_micro,
                          "_res_",
                          cluster_res_micro,
                          ".csv"),
            row.names = FALSE)
}



#2. find the list of background genes
all_data <- GetAssayData(dat_micro, assay = "SCT")
background_genes <- rownames(all_data[rowSums(all_data[])>0,])
write.csv(as.data.frame(background_genes),
          file = paste0("nonzero_bakcground_gene_list_",
                        "microglia_data_no_S521G_pcadims_",
                        pca_dim_micro,
                        "_res_",
                        cluster_res_micro,
                        ".csv"),
          row.names = FALSE)

print("********** Script completed! **********")

################## END ################## 

