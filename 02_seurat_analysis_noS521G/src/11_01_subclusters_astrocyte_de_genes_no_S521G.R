#!/usr/bin/env Rscript

#load required packages
library(Seurat)
library(dplyr)

#set all directory paths
basedir <- paste0("/gladstone/bioinformatics/projects/",
                  "mn-1318-maxine-nelson-yadong-huang-snrnaseq-mm10-aug-2022/",
                  "results/data")
indir <- paste0(basedir,
                "/05_subcluster_astrocyte_clusters_13_36/")
outdir <- paste0(basedir,
                 "/11_subclusters_de_genes_no_S521G/astrocyte/")
setwd(outdir)

#load the Seurat object
pca_dim_astro <- 15
cluster_res_astro <- 0.9
dat_astro <- readRDS(paste0(indir,
                            "astrocyte_data_post_subcluster_no_S521G_pcadims_",
                            pca_dim_astro,
                            "_res_",
                            cluster_res_astro,
                            ".rds"))

# 1. DE gene for:
## a. subcluster 3 vs other astrocyte
## b. subcluster 5 vs other astrocyte
## c. subcluster 7 vs other astrocyte
## d. For subcluster 3, PS19-E4-R/S vs PS19-E4 and PS19-E4-S/S vs PS19-E4.
## e. For subcluster 5, PS19-E4-R/S vs PS19-E4 and PS19-E4-S/S vs PS19-E4.
## f. For subcluster 7, PS19-E4-R/S vs PS19-E4 and PS19-E4-S/S vs PS19-E4.
clusters_for_de <- c(3,5,7)
other_clusters <- seq(1, max(as.numeric(dat_astro$seurat_clusters)))
other_clusters <- other_clusters[!(other_clusters %in% clusters_for_de)]
print(other_clusters)

for(de_clus in clusters_for_de){
  de_list <- FindMarkers(object = dat_astro, 
                         ident.1 = de_clus,
                         ident.2 = other_clusters,
                         assay = "SCT", 
                         slot = "data", 
                         test.use = "wilcox",
                         logfc.threshold = 0.1)
  
  de_cluster_rs_vs_e4 <- FindMarkers(object = dat_astro,
                                     ident.1 = "fE4R_S", 
                                     ident.2 = "fE4",
                                     verbose = TRUE, 
                                     group.by="genotype",
                                     subset.ident = de_clus,
                                     assay = "SCT", 
                                     slot = "data", 
                                     test.use = "wilcox",
                                     logfc.threshold = 0.1)
  
  de_cluster_ss_vs_e4 <- FindMarkers(object = dat_astro,
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
            file = paste0("de_genes_astrocyte_subcluster_",
                          de_clus,
                          "_vs_other_astrocyte_",
                          "no_S521G_pcadims_",
                          pca_dim_astro,
                          "_res_",
                          cluster_res_astro,
                          ".csv"),
            row.names = FALSE)
  write.csv(de_cluster_rs_vs_e4,
            file = paste0("de_genes_astrocyte_subcluster_",
                          de_clus,
                          "_PS19-E4-R_S_vs_PS19-E4_",
                          "no_S521G_pcadims_",
                          pca_dim_astro,
                          "_res_",
                          cluster_res_astro,
                          ".csv"),
            row.names = FALSE)
  write.csv(de_cluster_ss_vs_e4,
            file = paste0("de_genes_astrocyte_subcluster_",
                          de_clus,
                          "_PS19-E4-S_S_vs_PS19-E4_",
                          "no_S521G_pcadims_",
                          pca_dim_astro,
                          "_res_",
                          cluster_res_astro,
                          ".csv"),
            row.names = FALSE)
}



#2. find the list of background genes
all_data <- GetAssayData(dat_astro, assay = "SCT")
background_genes <- rownames(all_data[rowSums(all_data[])>0,])
write.csv(as.data.frame(background_genes),
          file = paste0("nonzero_bakcground_gene_list_",
                        "astrocyte_data_no_S521G_pcadims_",
                        pca_dim_astro,
                        "_res_",
                        cluster_res_astro,
                        ".csv"),
          row.names = FALSE)

print("********** Script completed! **********")

################## END ################## 

