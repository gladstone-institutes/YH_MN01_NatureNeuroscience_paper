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

#sub cluster the microglia clusters 17 and 19
dat_microglia <- subset(x = dat, idents = c(17,19))

#normalization
dat_microglia <- SCTransform(dat_microglia, 
                              method="glmGamPoi")
#PCA
dat_microglia <- RunPCA(dat_microglia, 
                         assay = "SCT", 
                         verbose = FALSE, 
                         approx=FALSE)
pdf(paste0(outdir,
           "/plot/06_subcluster_microglia_clusters_17_19/",
           "microglia_data_no_S521G_pcaplot.pdf"))
print(ElbowPlot(dat_microglia))
print(ElbowPlot(dat_microglia, ndims = 50))
dev.off()

#identify the list of marker genes for microglia 
marker_genes_microglia <- read.csv(paste0("/gladstone/bioinformatics/projects/",
                                           "mn-1318-maxine-nelson-yadong-huang-",
                                           "snrnaseq-mm10-aug-2022/data/",
                                           "Hippocampus_cellcluster_",
                                           "MarkerGenes.csv"))
marker_genes_microglia <- 
  marker_genes_microglia$marker_gene_symbol[
    marker_genes_microglia$hippocampal_region == "Microglia"]
marker_genes_microglia <- c(marker_genes_microglia,
                             "Apoe","hapoE-transgene",
                             "Mapt","Human-MAPT")
marker_genes_homeostatic_microglia <- c("P2ry12","Csf1r","Hexb","Cst3",
                                        "Cx3cr1","Siglech","Tgfbr1","Selplg",
                                        "Mef2a","Serinc3")
marker_genes_dam <- c("Cd9","Fth1","Plp1")
marker_genes_erbb_pathway <- c("Nrg1","Nrg2","Nrg3","Camk2a",
                               "Akt3","Ptk2","Mapk8","Erbb4")
marker_genes_alzheimers_disease <- c("Grin1","Grin2a","Grin2b","Itpr1",
                                     "Itpr2","Ryr2","Ryr3","Plcb1",
                                     "Plcb4","Cacna1c","Ppp3ca")
marker_genes_microglia <- unique(c(marker_genes_microglia,
                                   marker_genes_homeostatic_microglia,
                                   marker_genes_dam,
                                   marker_genes_erbb_pathway,
                                   marker_genes_alzheimers_disease))

#set the number of PCs and resolution to be used for clustering as recommended by Yadong
micro_pca_dim <- 15
micro_cluster_res <- 0.9

#perform clustering and UMAP 
dat_microglia_processed <- FindNeighbors(dat_microglia, 
                                          dims = 1:micro_pca_dim)
dat_microglia_processed <- RunUMAP(dat_microglia_processed, 
                                    dims = 1:micro_pca_dim)
dat_microglia_processed <- FindClusters(dat_microglia_processed, 
                                         resolution =  micro_cluster_res)

#rename the cluster ids to start from 1 instead of 0
dat_microglia_processed$orig_seurat_clusters <- dat_microglia_processed$seurat_clusters
dat_microglia_processed$seurat_clusters <- factor(as.numeric(
  as.character(dat_microglia_processed$orig_seurat_clusters))+1)
Idents(dat_microglia_processed) <- dat_microglia_processed$seurat_clusters

#visualizations
#i. #Generate UMAP without cluster labels
DimPlot(dat_microglia_processed, 
        pt.size = 1,
        raster = FALSE, 
        order = TRUE, 
        label = FALSE, 
        reduction = "umap") %>%
  ggsave(file = file.path(outdir, 
                          paste0("plot/06_subcluster_microglia_clusters_17_19/",
                                 "microglia_data_no_S521G", 
                                 "_pcadims_",
                                 micro_pca_dim,
                                 "_res_",
                                 micro_cluster_res,
                                 "_umap.pdf")),
         plot = .,
         width = 12,
         height = 7)

#ii. Generate UMAP with cluster labels
DimPlot(dat_microglia_processed, 
        pt.size = 1,
        raster = FALSE, 
        order = TRUE, 
        label = TRUE, 
        reduction = "umap") %>%
  ggsave(file = file.path(outdir, 
                          paste0("plot/06_subcluster_microglia_clusters_17_19/",
                                 "microglia_data_no_S521G", 
                                 "_pcadims_",
                                 micro_pca_dim,
                                 "_res_",
                                 micro_cluster_res,
                                 "_labeled_umap.pdf")),
         plot = .,
         width = 12,
         height = 7)

#iii. Generate UMAP split by genotypes
pdf(file = paste0(outdir,
                  "/plot/06_subcluster_microglia_clusters_17_19/",
                  "microglia_data_no_S521G", 
                  "_pcadims_",
                  micro_pca_dim,
                  "_res_",
                  micro_cluster_res,
                  "_split_genotype_umap.pdf"),
    width = 15,
    height = 7)
print(DimPlot(dat_microglia_processed, 
              group.by = "genotype",
              pt.size = 1,
              raster=FALSE,
              order = TRUE))
print(DimPlot(dat_microglia_processed, 
              split.by = "genotype", 
              group.by = "genotype",
              pt.size = 1,
              raster=FALSE,
              order = TRUE))
print(DimPlot(dat_microglia_processed, 
              split.by = "genotype",
              pt.size = 1,
              raster=FALSE,
              order = TRUE))
print(DimPlot(dat_microglia_processed, 
              split.by = "genotype",
              pt.size = 1,
              raster=FALSE,
              order = TRUE,
              label = TRUE))
dev.off()

#iv. Generate featureplot of all microglia marker genes
pdf(file = paste0(outdir,
                  "/plot/06_subcluster_microglia_clusters_17_19/",
                  "microglia_data_no_S521G", 
                  "_pcadims_",
                  micro_pca_dim,
                  "_res_",
                  micro_cluster_res,
                  "_featureplot.pdf"),
    width = 14,
    height = 7)
for(f in 1:length(marker_genes_microglia)){
  print(FeaturePlot(dat_microglia_processed, 
                    features = marker_genes_microglia[f],
                    pt.size = 1.5,
                    raster=FALSE,
                    order = TRUE))
}
dev.off()
pdf(file = paste0(outdir,
                  "/plot/06_subcluster_microglia_clusters_17_19/",
                  "microglia_data_no_S521G", 
                  "_pcadims_",
                  micro_pca_dim,
                  "_res_",
                  micro_cluster_res,
                  "_labeled_featureplot.pdf"),
    width = 14,
    height = 7)
for(f in 1:length(marker_genes_microglia)){
  print(FeaturePlot(dat_microglia_processed, 
                    features = marker_genes_microglia[f],
                    pt.size = 1.5,
                    raster=FALSE,
                    order = TRUE,
                    label = TRUE))
}
dev.off()
pdf(file = paste0(outdir,
                  "/plot/06_subcluster_microglia_clusters_17_19/",
                  "microglia_data_no_S521G", 
                  "_pcadims_",
                  micro_pca_dim,
                  "_res_",
                  micro_cluster_res,
                  "_split_genotype_featureplot.pdf"),
    width = 25,
    height = 7)
for(f in 1:length(marker_genes_microglia)){
  print(FeaturePlot(dat_microglia_processed, 
                    features= marker_genes_microglia[f], 
                    pt.size = 1.5,
                    split.by = "genotype",
                    raster = FALSE, 
                    order = TRUE, 
                    label = FALSE))
  print(FeaturePlot(dat_microglia_processed, 
                    features= marker_genes_microglia[f],
                    pt.size = 1.5,
                    split.by = "genotype",
                    raster = FALSE, 
                    order = TRUE, 
                    label = TRUE))
}
dev.off()

#v. Generate dotplot of all microglia marker genes
pdf(file = paste0(outdir,
                  "/plot/06_subcluster_microglia_clusters_17_19/",
                  "microglia_data_no_S521G", 
                  "_pcadims_",
                  micro_pca_dim,
                  "_res_",
                  micro_cluster_res,
                  "_dotplot.pdf"),
    width = 24,
    height = 7)
print(DotPlot(dat_microglia_processed, 
              features = marker_genes_microglia
))
print(DotPlot(dat_microglia_processed, 
              features = marker_genes_microglia, 
              group.by = "genotype"
))
dev.off()

#vi. #Generate visualizations using the various metadata
meta_features <- c("new_samplenames", "genotype", "date_of_birth",
                   "date_perfused", "age_at_perfusion", "sex", 
                   "date_of_nuclear_isolation")
for(mf in meta_features){
  DimPlot(dat_microglia_processed, 
          pt.size = 1,
          raster = FALSE, 
          order = TRUE, 
          label = FALSE,
          group.by = mf, 
          reduction = "umap") %>%
    ggsave(file = file.path(outdir, 
                            paste0("plot/06_subcluster_microglia_clusters_17_19/",
                                   "microglia_data_no_S521G", 
                                   "_pcadims_",
                                   micro_pca_dim,
                                   "_res_",
                                   micro_cluster_res, "_",
                                   mf, "_umap.pdf")),
           plot = .,
           width = 12,
           height = 7)
  DimPlot(dat_microglia_processed, 
          pt.size = 1,
          raster = FALSE, 
          order = TRUE, 
          label = FALSE, 
          group.by = mf,
          split.by = mf, 
          reduction = "umap") %>%
    ggsave(file = file.path(outdir, 
                            paste0("plot/06_subcluster_microglia_clusters_17_19/",
                                   "microglia_data_no_S521G", 
                                   "_pcadims_",
                                   micro_pca_dim,
                                   "_res_",
                                   micro_cluster_res, "_",
                                   mf, "_split_umap.pdf")),
           plot = .,
           width = 30,
           height = 10)
}

#vii. #generate different UMAPs taht are split by sex
DimPlot(dat_microglia_processed, 
        pt.size = 1,
        raster = FALSE, 
        order = TRUE, 
        label = TRUE, 
        split.by = "sex",
        reduction = "umap") %>%
  ggsave(file = file.path(outdir, 
                          paste0("plot/06_subcluster_microglia_clusters_17_19/",
                                 "microglia_data_no_S521G", 
                                 "_pcadims_",
                                 micro_pca_dim,
                                 "_res_",
                                 micro_cluster_res,
                                 "_sex_split_labeled_umap.pdf")),
         plot = .,
         width = 30,
         height = 10)
DimPlot(dat_microglia_processed, 
        pt.size = 1,
        raster = FALSE, 
        order = TRUE, 
        label = FALSE,
        group.by = "new_samplenames",
        split.by = "sex",
        reduction = "umap") %>%
  ggsave(file = file.path(outdir, 
                          paste0("plot/06_subcluster_microglia_clusters_17_19/",
                                 "microglia_data_no_S521G", 
                                 "_pcadims_",
                                 micro_pca_dim,
                                 "_res_",
                                 micro_cluster_res,
                                 "_sex_split_samplename_labeled_umap.pdf")),
         plot = .,
         width = 30,
         height = 10)

#save the microglia processed sub-clustered data
saveRDS(dat_microglia_processed, 
        file = paste0(outdir,
                      "/data/06_subcluster_microglia_clusters_17_19/",
                      "microglia_data_post_subcluster_no_S521G_pcadims_",
                      micro_pca_dim,
                      "_res_",
                      micro_cluster_res,
                      ".rds"))

#Find differentially expressed genes for ID'ing the clusters
MarkersRes <- FindAllMarkers(dat_microglia_processed, 
                             assay = "SCT", 
                             slot = "data", 
                             test.use = "wilcox",
                             only.pos = TRUE,
                             min.pct = 0.25,
                             logfc.threshold = 0.5)
write.csv(MarkersRes,
          file = file.path(outdir, 
                           paste0("/data/06_subcluster_microglia_clusters_17_19/",
                                  "marker_genes_microglia_data_post_subcluster",
                                  "_no_S521G_pcadims_",
                                  micro_pca_dim,
                                  "_res_",
                                  micro_cluster_res,
                                  ".csv")),
          row.names = FALSE)


#Find differentially expressed genes for ID'ing the clusters
#no marker genes detected for subcluster 4
#re-run FindAllMarkers with default only.pos=FALSE
MarkersRes <- FindAllMarkers(dat_microglia_processed,
                             assay = "SCT",
                             slot = "data",
                             test.use = "wilcox",
                             #only.pos = TRUE,
                             min.pct = 0.25,
                             logfc.threshold = 0.25)
write.csv(MarkersRes,
          file = file.path(outdir,
                           paste0("/data/06_subcluster_microglia_clusters_17_19/",
                                  "marker_genes_only_pos_false_microglia_data_post_subcluster",
                                  "_no_S521G_pcadims_",
                                  micro_pca_dim,
                                  "_res_",
                                  micro_cluster_res,
                                  ".csv")),
          row.names = FALSE)


print("*********** Script completed! ***********")

############### END ###############

