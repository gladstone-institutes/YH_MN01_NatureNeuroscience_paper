#run locally on Ayushi's laptop

#load required packages
library(Seurat)
library(dplyr)
library(ggplot2)
library(scales)

#set all directory paths
basedir <- paste0("~/Dropbox (Gladstone)/YH_MN02/",
                  "seurat_analysis_no_S521G/")
outdir <- paste0("~/Dropbox (Gladstone)/YH_MN02/",
                 "seurat_analysis_no_S521G/14_paper_figures/")
if(!dir.exists(outdir)){
  dir.create(outdir)
}
setwd(outdir)


#load the Seurat object
pca_dim <- 15
cluster_res <- 0.7
dat <- readRDS(paste0(basedir,
                      "02_clustering_no_S521G/",
                      "sct_data_no_S521G_post_cluster_pcadims_",
                      pca_dim,
                      "_res_",
                      cluster_res,
                      ".rds"))

#load the astrocytes Seurat object
pca_dim_astro <- 15
cluster_res_astro <- 0.9
dat_astrocytes <- readRDS(paste0(basedir,
                                 "05_subcluster_astrocyte_clusters_13_36/",
                                 "astrocyte_data_post_subcluster_no_S521G_pcadims_",
                                 pca_dim_astro,
                                 "_res_",
                                 cluster_res_astro,
                                 ".rds"))

#load the microglia Seurat object
pca_dim_micro <- 15
cluster_res_micro <- 0.9
dat_microglia <- readRDS(paste0(basedir,
                                "06_subcluster_microglia_clusters_17_19/",
                                "microglia_data_post_subcluster_no_S521G_pcadims_",
                                pca_dim_micro,
                                "_res_",
                                cluster_res_micro,
                                ".rds"))

##############################################################
##1.A dot plot of all 38 clusters for the following genes:
#   Syn1, Dgkh, Trhde, Slc17a7, Man1a, Galntl6, Meis2, Slit2, 
#   Mbp, Vcan, Phkg1, Aldh1l1, Gfap, Csf1r, Trem2, Cx3cr1, Folr1.
##############################################################
genes_to_plot <- c("Syn1", "Dgkh", "Trhde", "Slc17a7", "Man1a", "Galntl6", 
                   "Meis2", "Slit2", "Mbp", "Vcan", "Phkg1", "Aldh1l1", "Gfap", 
                   "Csf1r", "Trem2", "Cx3cr1", "Folr1")
(DotPlot(dat, 
         features = genes_to_plot, 
         cols = c("blue", "red")) + 
    RotatedAxis())%>%
  ggsave(file = "dotplot_with_selected_genes.pdf",
         plot = .,
         width = 12,
         height = 9)


##############################################################
##2. Make a genotype-split human-APOE feature plot (with 
#   labels) for all 38 clusters. Please make 3 versions with 
#   different APOE expression scales (1.5, 2.0, 2.5?)
##############################################################
#option #1
(FeaturePlot(dat, 
             features = "hapoE-transgene",
             split.by = "genotype",
             raster = FALSE, 
             order = TRUE,
             label = TRUE,
             max.cutoff = 2.5) + 
   theme(legend.position = "right")) %>%
  ggsave(file = "featureplot_genotype_split_hAPOE_maxcutoff_2.5.pdf",
         plot = .,
         width = 20,
         height = 7)
#option #2
(FeaturePlot(dat, 
             features = "hapoE-transgene",
             split.by = "genotype",
             raster = FALSE, 
             order = TRUE,
             label = TRUE,
             max.cutoff = 2) + 
    theme(legend.position = "right")) %>%
  ggsave(file = "featureplot_genotype_split_hAPOE_maxcutoff_2.0.pdf",
         plot = .,
         width = 20,
         height = 7)
#option #3
(FeaturePlot(dat, 
             features = "hapoE-transgene",
             split.by = "genotype",
             raster = FALSE, 
             order = TRUE,
             label = TRUE,
             max.cutoff = 1.5) +  
   theme(legend.position = "right")) %>%
  ggsave(file = "featureplot_genotype_split_hAPOE_maxcutoff_1.5.pdf",
         plot = .,
         width = 20,
         height = 7)
#option #4
(FeaturePlot(dat, 
             features = "hapoE-transgene",
             split.by = "genotype",
             raster = FALSE, 
             order = TRUE,
             label = TRUE,
             max.cutoff = 1) + 
    theme(legend.position = "right")) %>%
  ggsave(file = "featureplot_genotype_split_hAPOE_maxcutoff_1.0.pdf",
         plot = .,
         width = 20,
         height = 7)


##############################################################
##3. Make a genotype-split UMAP for all 38 clusters, with the 
#    following clusters highlighted in different colors: 
#    1, 6, 7, 9, 28. All other clusters in gray color.
##############################################################
dat$clusters_of_interest <- ifelse(dat$seurat_clusters %in% c(1, 6, 7, 9, 28),
                                   dat$seurat_clusters,
                                   0)
#option 1
(DimPlot(dat, 
         split.by = "genotype",
         group.by = "clusters_of_interest",
         raster = FALSE, 
         order = TRUE, 
         cols = c("grey",hue_pal()(5)),
         reduction = "umap") + 
    scale_color_manual(labels = c("Other", "Cluster1", 
                                  "Cluster6", "Cluster7",
                                  "Cluster9","Cluster28"),
                       values = c("grey",hue_pal()(5)))) %>%
  ggsave(file = "umap_genotype_split_clusters_of_interest.pdf",
         plot = .,
         width = 20,
         height = 7)
#option 2
(DimPlot(dat, 
         split.by = "genotype",
         group.by = "clusters_of_interest",
         raster = TRUE, 
         order = TRUE, 
         cols = c("grey",hue_pal()(5)),
         reduction = "umap") + 
    scale_color_manual(labels = c("Other", "Cluster1", 
                                  "Cluster6", "Cluster7",
                                  "Cluster9","Cluster28"),
                       values = c("grey",hue_pal()(5)))) %>%
  ggsave(file = "umap_genotype_split_clusters_of_interest_rasterized.pdf",
         plot = .,
         width = 20,
         height = 7)
#option 3
(DimPlot(dat, 
         split.by = "genotype",
         group.by = "clusters_of_interest",
         raster = FALSE, 
         order = TRUE, 
         cols = c("grey",hue_pal()(5)),
         reduction = "umap") + 
    scale_color_manual(labels = c("Other", "Cluster1", 
                                  "Cluster6", "Cluster7",
                                  "Cluster9","Cluster28"),
                       values = c("grey",hue_pal()(5)))) %>%
  ggsave(file = "umap_genotype_split_clusters_of_interest_smaller_width.pdf",
         plot = .,
         width = 12,
         height = 7)
#option 4
(DimPlot(dat, 
         split.by = "genotype",
         group.by = "clusters_of_interest",
         raster = TRUE, 
         order = TRUE, 
         cols = c("grey",hue_pal()(5)),
         reduction = "umap") + 
    scale_color_manual(labels = c("Other", "Cluster1", 
                                  "Cluster6", "Cluster7",
                                  "Cluster9","Cluster28"),
                       values = c("grey",hue_pal()(5)))) %>%
  ggsave(file = "umap_genotype_split_clusters_of_interest_rasterized_and_smaller_width.pdf",
         plot = .,
         width = 12,
         height = 7)
#option 5
(DimPlot(dat, 
         split.by = "genotype",
         group.by = "clusters_of_interest",
         raster = FALSE, 
         order = TRUE, 
         cols = c("grey",hue_pal()(5)),
         reduction = "umap") + 
    scale_color_manual(labels = c("Other", "Cluster1", 
                                  "Cluster6", "Cluster7",
                                  "Cluster9","Cluster28"),
                       values = c("grey",hue_pal()(5)))+
    theme(text = element_text(size = 12))) %>%
  ggsave(file = "umap_genotype_split_clusters_of_interest_smaller_fontsize.pdf",
         plot = .,
         width = 10,
         height = 7)


##############################################################
##4. Make a genotype-split human APOE feature plot (with labels)
#    for all 12 astrocyte subclusters. Please make 3 versions 
#    with different APOE expression scales (1.5, 2.0, 2.5?).
##############################################################
#option #1
(FeaturePlot(dat_astrocytes, 
             features = "hapoE-transgene",
             split.by = "genotype",
             raster = FALSE, 
             order = TRUE,
             label = TRUE,
             max.cutoff = 2.5) + 
    theme(legend.position = "right")) %>%
  ggsave(file = "astrocyte_featureplot_genotype_split_hAPOE_maxcutoff_2.5.pdf",
         plot = .,
         width = 15,
         height = 7)
#option #2
(FeaturePlot(dat_astrocytes, 
             features = "hapoE-transgene",
             split.by = "genotype",
             raster = FALSE, 
             order = TRUE,
             label = TRUE,
             max.cutoff = 2) + 
    theme(legend.position = "right")) %>%
  ggsave(file = "astrocyte_featureplot_genotype_split_hAPOE_maxcutoff_2.0.pdf",
         plot = .,
         width = 15,
         height = 7)
#option #3
(FeaturePlot(dat_astrocytes, 
             features = "hapoE-transgene",
             split.by = "genotype",
             raster = FALSE, 
             order = TRUE,
             label = TRUE,
             max.cutoff = 1.5) +  
    theme(legend.position = "right")) %>%
  ggsave(file = "astrocyte_featureplot_genotype_split_hAPOE_maxcutoff_1.5.pdf",
         plot = .,
         width = 15,
         height = 7)
#option #4
(FeaturePlot(dat_astrocytes, 
             features = "hapoE-transgene",
             split.by = "genotype",
             raster = FALSE, 
             order = TRUE,
             label = TRUE,
             max.cutoff = 1) + 
    theme(legend.position = "right")) %>%
  ggsave(file = "astrocyte_featureplot_genotype_split_hAPOE_maxcutoff_1.0.pdf",
         plot = .,
         width = 15,
         height = 7)


##############################################################
##5. Make a genotype-split UMAP for all 12 astrocyte subclusters, 
#    with the following clusters highlighted in different 
#    colors: 3, 5, 7. All other clusters in gray color.
##############################################################
dat_astrocytes$clusters_of_interest <- 
  ifelse(dat_astrocytes$seurat_clusters %in% c(3, 5, 7),
         dat_astrocytes$seurat_clusters,
         0)
(DimPlot(dat_astrocytes, 
         split.by = "genotype",
         group.by = "clusters_of_interest",
         raster = FALSE, 
         order = TRUE, 
         cols = c("grey",hue_pal()(3)),
         reduction = "umap") + 
    scale_color_manual(labels = c("Other", "Cluster3", "Cluster5", "Cluster7"),
                       values = c("grey",hue_pal()(3)))) %>%
  ggsave(file = "astrocyte_umap_genotype_split_clusters_of_interest.pdf",
         plot = .,
         width = 15,
         height = 7)


##############################################################
##6. Make a genotype-split human APOE feature plot (with labels) 
#   for all 15 microglia subclusters. Please make 3 versions 
#   with different APOE expression scales (1.0, 1.5, 2.0?).
##############################################################
#range of human-APOE is 0-1.386294 for microglia
#so, generate featureplot for max.cutoff=1
(FeaturePlot(dat_microglia, 
             features = "hapoE-transgene",
             split.by = "genotype",
             raster = FALSE, 
             order = TRUE,
             label = TRUE,
             max.cutoff = 1) + 
    theme(legend.position = "right")) %>%
  ggsave(file = "microglia_featureplot_genotype_split_hAPOE_maxcutoff_1.0.pdf",
         plot = .,
         width = 15,
         height = 7)


##############################################################
##7. Make a genotype-split UMAP for all 15 microglia subclusters, 
#   with the following clusters highlighted in different 
#   colors: 2, 8, 11. All other clusters in gray color.
##############################################################
dat_microglia$clusters_of_interest <- 
  ifelse(dat_microglia$seurat_clusters %in% c(2, 8, 11),
         dat_microglia$seurat_clusters,
         0)
(DimPlot(dat_microglia, 
         split.by = "genotype",
         group.by = "clusters_of_interest",
         raster = FALSE, 
         order = TRUE, 
         cols = c("grey",hue_pal()(3)),
         reduction = "umap") + 
    scale_color_manual(labels = c("Other", "Cluster2", "Cluster8", "Cluster11"),
                       values = c("grey",hue_pal()(3)))) %>%
  ggsave(file = "microglia_umap_genotype_split_clusters_of_interest.pdf",
         plot = .,
         width = 15,
         height = 7)


##############################################################
##8. UMAP for all Cell Clusters with labels (square dimension if possible)
##############################################################
#option 1
(DimPlot(dat, 
         raster = FALSE, 
         order = TRUE, 
         label = TRUE, 
         reduction = "umap")  +  
   theme(legend.key.size = unit(0.36, "cm"),
         legend.text = element_text(size = 8)) +
   guides(color = guide_legend(override.aes = list(size=1),
                               ncol=1) )) %>%
  ggsave(file = paste0(outdir, 
                       "all_cell_clusters_labeled_umap.pdf"),
         plot = .,
         width = 7,
         height = 7)
#option 2
(DimPlot(dat, 
         raster = TRUE, 
         order = TRUE, 
         label = TRUE, 
         reduction = "umap")  +  
    theme(legend.key.size = unit(0.36, "cm"),
          legend.text = element_text(size = 8)) +
    guides(color = guide_legend(override.aes = list(size=1),
                                ncol=1) )) %>%
  ggsave(file = paste0(outdir, 
                       "all_cell_clusters_labeled_umap_rasterized.pdf"),
         plot = .,
         width = 7,
         height = 7)


##############################################################
##9. UMAP of Astrocyte Subclusters with labels (square dimension if possible)
##############################################################
(DimPlot(dat_astrocytes, 
         raster = FALSE, 
         order = TRUE, 
         label = TRUE, 
         reduction = "umap")  +  
   theme(legend.key.size = unit(0.36, "cm"),
         legend.text = element_text(size = 8)) +
   guides(color = guide_legend(override.aes = list(size=1),
                               ncol=1) )) %>%
  ggsave(file = paste0(outdir, 
                       "astrocyte_labeled_umap.pdf"),
         plot = .,
         width = 7,
         height = 7)


##############################################################
##10. UMAP of Microglia Subclusters with labels (square dimension if possible)
##############################################################
(DimPlot(dat_microglia, 
         raster = FALSE, 
         order = TRUE, 
         label = TRUE, 
         reduction = "umap")  +  
   theme(legend.key.size = unit(0.36, "cm"),
         legend.text = element_text(size = 8)) +
   guides(color = guide_legend(override.aes = list(size=1),
                               ncol=1) )) %>%
  ggsave(file = paste0(outdir, 
                       "microglia_labeled_umap.pdf"),
         plot = .,
         width = 7,
         height = 7)


##############################################################
##11.Dotplot of selected genes for all 12 astrocyte subclusters.
##############################################################
genes_to_plot_astro <- c("Luzp2","Slc7a10","Mfge8","Gfap","Id3","Aqp4","Myoc","Id1",
                   "Fabp7","Ctsb","Vim","Osmr","Serpina3n","Gsn","Ggta1","Trpm3",
                   "Csmd1","C4b","Cd9","Sparcl1","Plce1","Sgcd","Fos","S100a6",
                   "hapoE-transgene","Clu")
(DotPlot(dat_astrocytes, 
         features = genes_to_plot_astro, 
         cols = c("blue", "red")) + 
    RotatedAxis())%>%
  ggsave(file = "astrocyte_dotplot_with_selected_genes.pdf",
         plot = .,
         width = 12,
         height = 9)


##############################################################
##12.Dotplot of selected genes for all 12 astrocyte subclusters.
##############################################################
genes_to_plot_micro <- c("Hexb","Cst3","Cx3cr1","Ctsd","Csf1r","Ctss","Sparc",
                         "Tmsb4x","P2ry12","C1qa","C1qb","Tmem119","Tyrobp","Ctsb",
                         "hapoE-transgene","B2m","Fth1","Lyz2","Trem2","Axl",
                         "Cst7","Ctsl","Lpl","Cd9","Csf1","Ccl6","Itgax","Timp2")
(DotPlot(dat_microglia, 
         features = genes_to_plot_micro, 
         cols = c("blue", "red")) + 
    RotatedAxis())%>%
  ggsave(file = "microglia_dotplot_with_selected_genes.pdf",
         plot = .,
         width = 12,
         height = 9)


################### END ###################
