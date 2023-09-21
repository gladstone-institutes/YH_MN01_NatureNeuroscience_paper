#!/usr/bin/env Rscript

#load required packages
library(Seurat)
library(dplyr)
library(ggplot2)
library(gridExtra)

#set all directory paths
outdir <- paste0("/gladstone/bioinformatics/projects/",
                 "mn-1318-maxine-nelson-yadong-huang-snrnaseq-mm10-aug-2022/",
                 "results")
indir <- paste0(outdir,
                "/data/02_clustering_no_S521G/")
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

#sub cluster the astrocytes clusters 13 and 36
dat_astrocytes <- subset(x = dat, idents = c(13,36))

#normalization
dat_astrocytes <- SCTransform(dat_astrocytes, 
                              method="glmGamPoi")
#PCA
dat_astrocytes <- RunPCA(dat_astrocytes, 
                         assay = "SCT", 
                         verbose = FALSE, 
                         approx=FALSE)
pdf(paste0(outdir,
           "/plot/05_subcluster_astrocyte_clusters_13_36/",
           "astrocyte_data_no_S521G_pcaplot.pdf"))
print(ElbowPlot(dat_astrocytes))
print(ElbowPlot(dat_astrocytes, ndims = 50))
dev.off()

#identify the list of marker genes for astrocytes 
marker_genes_astrocytes <- read.csv(paste0("/gladstone/bioinformatics/projects/",
                                           "mn-1318-maxine-nelson-yadong-huang-",
                                           "snrnaseq-mm10-aug-2022/data/",
                                           "Hippocampus_cellcluster_",
                                           "MarkerGenes.csv"
)
)
marker_genes_astrocytes <- 
  marker_genes_astrocytes$marker_gene_symbol[
    marker_genes_astrocytes$hippocampal_region == "Astrocytes"]
marker_genes_astrocytes <- c(marker_genes_astrocytes,
                             "Apoe","hapoE-transgene",
                             "Mapt","Human-MAPT")

#set the number of PCs and resolution to be used for clustering as recommended by Yadong
astro_pca_dim <- 15
astro_cluster_res <- 0.9

#perform clustering and UMAP 
dat_astrocytes_processed <- FindNeighbors(dat_astrocytes, 
                                          dims = 1:astro_pca_dim)
dat_astrocytes_processed <- RunUMAP(dat_astrocytes_processed, 
                                    dims = 1:astro_pca_dim)
dat_astrocytes_processed <- FindClusters(dat_astrocytes_processed, 
                                         resolution =  astro_cluster_res)

#rename the cluster ids to start from 1 instead of 0
dat_astrocytes_processed$orig_seurat_clusters <- dat_astrocytes_processed$seurat_clusters
dat_astrocytes_processed$seurat_clusters <- factor(as.numeric(
  as.character(dat_astrocytes_processed$orig_seurat_clusters))+1)
Idents(dat_astrocytes_processed) <- dat_astrocytes_processed$seurat_clusters

#visualizations
#i. #Generate UMAP without cluster labels
DimPlot(dat_astrocytes_processed, 
        pt.size = 1,
        raster = FALSE, 
        order = TRUE, 
        label = FALSE, 
        reduction = "umap") %>%
  ggsave(file = file.path(outdir, 
                          paste0("plot/05_subcluster_astrocyte_clusters_13_36/",
                                 "astrocyte_data_no_S521G", 
                                 "_pcadims_",
                                 astro_pca_dim,
                                 "_res_",
                                 astro_cluster_res,
                                 "_umap.pdf")),
         plot = .,
         width = 12,
         height = 7)

#ii. Generate UMAP with cluster labels
DimPlot(dat_astrocytes_processed, 
        pt.size = 1,
        raster = FALSE, 
        order = TRUE, 
        label = TRUE, 
        reduction = "umap") %>%
  ggsave(file = file.path(outdir, 
                          paste0("plot/05_subcluster_astrocyte_clusters_13_36/",
                                 "astrocyte_data_no_S521G", 
                                 "_pcadims_",
                                 astro_pca_dim,
                                 "_res_",
                                 astro_cluster_res,
                                 "_labeled_umap.pdf")),
         plot = .,
         width = 12,
         height = 7)

#iii. Generate UMAP split by genotypes
pdf(file = paste0(outdir,
                  "/plot/05_subcluster_astrocyte_clusters_13_36/",
                  "astrocyte_data_no_S521G", 
                  "_pcadims_",
                  astro_pca_dim,
                  "_res_",
                  astro_cluster_res,
                  "_split_genotype_umap.pdf"),
    width = 15,
    height = 7)
print(DimPlot(dat_astrocytes_processed, 
              group.by = "genotype",
              pt.size = 1,
              raster=FALSE,
              order = TRUE))
print(DimPlot(dat_astrocytes_processed, 
              split.by = "genotype", 
              group.by = "genotype",
              pt.size = 1,
              raster=FALSE,
              order = TRUE))
print(DimPlot(dat_astrocytes_processed, 
              split.by = "genotype",
              pt.size = 1,
              raster=FALSE,
              order = TRUE))
print(DimPlot(dat_astrocytes_processed, 
              split.by = "genotype",
              pt.size = 1,
              raster=FALSE,
              order = TRUE,
              label = TRUE))
dev.off()

#iv. Generate featureplot of all astrocyte marker genes
pdf(file = paste0(outdir,
                  "/plot/05_subcluster_astrocyte_clusters_13_36/",
                  "astrocyte_data_no_S521G", 
                  "_pcadims_",
                  astro_pca_dim,
                  "_res_",
                  astro_cluster_res,
                  "_featureplot.pdf"),
    width = 14,
    height = 7)
for(f in 1:length(marker_genes_astrocytes)){
  print(FeaturePlot(dat_astrocytes_processed, 
                    features = marker_genes_astrocytes[f],
                    pt.size = 1.5,
                    raster=FALSE,
                    order = TRUE))
}
dev.off()
pdf(file = paste0(outdir,
                  "/plot/05_subcluster_astrocyte_clusters_13_36/",
                  "astrocyte_data_no_S521G", 
                  "_pcadims_",
                  astro_pca_dim,
                  "_res_",
                  astro_cluster_res,
                  "_labeled_featureplot.pdf"),
    width = 14,
    height = 7)
for(f in 1:length(marker_genes_astrocytes)){
  print(FeaturePlot(dat_astrocytes_processed, 
                    features = marker_genes_astrocytes[f],
                    pt.size = 1.5,
                    raster=FALSE,
                    order = TRUE,
                    label = TRUE))
}
dev.off()
pdf(file = paste0(outdir,
                  "/plot/05_subcluster_astrocyte_clusters_13_36/",
                  "astrocyte_data_no_S521G", 
                  "_pcadims_",
                  astro_pca_dim,
                  "_res_",
                  astro_cluster_res,
                  "_split_genotype_featureplot.pdf"),
    width = 25,
    height = 7)
for(f in 1:length(marker_genes_astrocytes)){
  print(FeaturePlot(dat_astrocytes_processed, 
                    features= marker_genes_astrocytes[f], 
                    pt.size = 1.5,
                    split.by = "genotype",
                    raster = FALSE, 
                    order = TRUE, 
                    label = FALSE))
  print(FeaturePlot(dat_astrocytes_processed, 
                    features= marker_genes_astrocytes[f],
                    pt.size = 1.5,
                    split.by = "genotype",
                    raster = FALSE, 
                    order = TRUE, 
                    label = TRUE))
}
dev.off()

#v. Generate dotplot of all astrocyte marker genes
pdf(file = paste0(outdir,
                  "/plot/05_subcluster_astrocyte_clusters_13_36/",
                  "astrocyte_data_no_S521G", 
                  "_pcadims_",
                  astro_pca_dim,
                  "_res_",
                  astro_cluster_res,
                  "_dotplot.pdf"),
    width = 14,
    height = 7)
print(DotPlot(dat_astrocytes_processed, 
              features = marker_genes_astrocytes
))
print(DotPlot(dat_astrocytes_processed, 
              features = marker_genes_astrocytes, 
              group.by = "genotype"
))
dev.off()

#vi. #Generate visualizations using the various metadata
meta_features <- c("new_samplenames", "genotype", "date_of_birth",
                   "date_perfused", "age_at_perfusion", "sex", 
                   "date_of_nuclear_isolation")
for(mf in meta_features){
  DimPlot(dat_astrocytes_processed, 
          pt.size = 1,
          raster = FALSE, 
          order = TRUE, 
          label = FALSE,
          group.by = mf, 
          reduction = "umap") %>%
    ggsave(file = file.path(outdir, 
                            paste0("plot/05_subcluster_astrocyte_clusters_13_36/",
                                   "astrocyte_data_no_S521G", 
                                   "_pcadims_",
                                   astro_pca_dim,
                                   "_res_",
                                   astro_cluster_res, "_",
                                   mf, "_umap.pdf")),
           plot = .,
           width = 12,
           height = 7)
  DimPlot(dat_astrocytes_processed, 
          pt.size = 1,
          raster = FALSE, 
          order = TRUE, 
          label = FALSE, 
          group.by = mf,
          split.by = mf, 
          reduction = "umap") %>%
    ggsave(file = file.path(outdir, 
                            paste0("plot/05_subcluster_astrocyte_clusters_13_36/",
                                   "astrocyte_data_no_S521G", 
                                   "_pcadims_",
                                   astro_pca_dim,
                                   "_res_",
                                   astro_cluster_res, "_",
                                   mf, "_split_umap.pdf")),
           plot = .,
           width = 30,
           height = 10)
}


#save the astrocyte processed sub-clustered data
saveRDS(dat_astrocytes_processed, 
        file = paste0(outdir,
                      "/data/05_subcluster_astrocyte_clusters_13_36/",
                      "astrocyte_data_post_subcluster_no_S521G_pcadims_",
                      astro_pca_dim,
                      "_res_",
                      astro_cluster_res,
                      ".rds"))

#Find differentially expressed genes for ID'ing the clusters
MarkersRes <- FindAllMarkers(dat_astrocytes_processed, 
                             assay = "SCT", 
                             slot = "data", 
                             test.use = "wilcox",
                             only.pos = TRUE,
                             min.pct = 0.25,
                             logfc.threshold = 0.5)
write.csv(MarkersRes,
          file = file.path(outdir, 
                           paste0("/data/05_subcluster_astrocyte_clusters_13_36/",
                                  "marker_genes_astrocyte_data_post_subcluster",
                                  "_no_S521G_pcadims_",
                                  astro_pca_dim,
                                  "_res_",
                                  astro_cluster_res,
                                  ".csv")),
          row.names = FALSE)


print("*********** Script completed! ***********")

############### END ###############

