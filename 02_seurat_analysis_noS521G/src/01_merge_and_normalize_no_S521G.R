#!/usr/bin/env Rscript

#load required packages
library(Seurat)
library(dplyr)
library(ggplot2)


#set all directory paths
indir <- paste0("/gladstone/bioinformatics/projects/",
                "mn-1318-maxine-nelson-yadong-huang-snrnaseq-mm10-aug-2022/",
                "data")
outdir <- paste0("/gladstone/bioinformatics/projects/",
                 "mn-1318-maxine-nelson-yadong-huang-snrnaseq-mm10-aug-2022/",
                 "results")
setwd(indir)

#load the seurat processed object
load("YH_MN02_Seurat.rdata")

#read in the sample metadata
mouse.data <- read.csv("DataonMice_snRNAseq_MetaData_02_AA.csv")
#modify the mouse number to match the new_samplenames in Seurat metadata
mouse.data$Mouse.. <- gsub("-", "_", mouse.data$Mouse..)
mouse.data$Mouse.. <- gsub("/", "_", mouse.data$Mouse..)
mouse.data <- mouse.data[ , !(names(mouse.data) == "Genotype")]

#add the metadata to the Seurat object
meta.data <- merged[[]] 
meta.data$date_of_birth <- mouse.data$DOB[match(meta.data$new_samplenames, 
                                                mouse.data$Mouse..)]
meta.data$date_perfused <- mouse.data$Date.Perfused[match(meta.data$new_samplenames, 
                                                          mouse.data$Mouse..)]
meta.data$age_at_perfusion <- mouse.data$Age.at.Perfusion[match(meta.data$new_samplenames, 
                                                                mouse.data$Mouse..)]
meta.data$sex <- mouse.data$Sex[match(meta.data$new_samplenames, 
                                      mouse.data$Mouse..)]
meta.data$date_of_nuclear_isolation <- mouse.data$Date.of.Nuc.Isolation[match(meta.data$new_samplenames, 
                                                                              mouse.data$Mouse..)]
merged <- AddMetaData(merged, meta.data)

#based on the QC results, remove sample "PS19_fE3_521G" from the merged data
DefaultAssay(merged) <- "RNA"
merged <- subset(merged, subset = new_samplenames != 'PS19_fE3_521G')
print(merged)

saveRDS(merged, 
        file = paste0(outdir,
                      "/data/01_merge_and_normalize_no_S521G/",
                      "merged_data_no_S521G.rds"))
rm(merged)

print("***** Dataset merge completed! *****")


#read in the merged dataset
data <- readRDS(paste0(outdir,
                       "/data/01_merge_and_normalize_no_S521G/",
                       "merged_data_no_S521G.rds"))

#normalization using SCTransform
data_sct <- SCTransform(data, 
                        method="glmGamPoi")
saveRDS(data_sct, 
        file = paste0(outdir,
                      "/data/01_merge_and_normalize_no_S521G/",
                      "merged_data_no_S521G_post_sct.rds"))

#Perform and PCA and UMAP on the SCTransformed data
data_sct_processed <- RunPCA(object = data_sct, 
                             assay = "SCT")
pdf(file.path(outdir, 
              paste0("/plot/01_merge_and_normalize_no_S521G/",
                     "sctransform_merged_data_no_S521G", 
                     "_pca_plot.pdf")),
    width = 12,
    height = 7)
ElbowPlot(data_sct_processed, ndims = 50, reduction = "pca")
ElbowPlot(data_sct_processed, ndims = 20, reduction = "pca")
ElbowPlot(data_sct_processed, ndims = 16, reduction = "pca")
dev.off()

data_sct_processed <- data_sct_processed %>% 
  RunUMAP(., dims = 1:15, verbose = TRUE)

write.csv(data_sct_processed[[]],
          file = paste0(outdir,
                        "/data/01_merge_and_normalize_no_S521G/",
                        "merged_data_no_S521G_",
                        "post_sct_processed_metadata.csv"))
saveRDS(data_sct_processed, 
        file = paste0(outdir,
                      "/data/01_merge_and_normalize_no_S521G/",
                      "merged_data_no_S521G_post_sct_processed.rds"))

#Generate visualizations using the various metadata
DimPlot(data_sct_processed, 
        raster = FALSE, 
        order = TRUE, 
        label = FALSE, 
        group.by="new_samplenames", 
        reduction = "umap") %>%
  ggsave(file = file.path(outdir, 
                          paste0("plot/01_merge_and_normalize_no_S521G/",
                                 "sctransform_merged_data_no_S521G", 
                                 "_sample_name_umap.pdf")),
         plot = .,
         width = 12,
         height = 7)
DimPlot(data_sct_processed, 
        raster = FALSE, 
        order = TRUE, 
        label = FALSE, 
        group.by="new_samplenames",
        split.by="new_samplenames", 
        reduction = "umap",
        ncol = 2) %>%
  ggsave(file = file.path(outdir, 
                          paste0("plot/01_merge_and_normalize_no_S521G/",
                                 "sctransform_merged_data_no_S521G", 
                                 "_sample_name_split_umap.pdf")),
         plot = .,
         width = 12,
         height = 48)

DimPlot(data_sct_processed, 
        raster = FALSE,
        order = TRUE, 
        label = FALSE, 
        group.by="genotype", 
        reduction = "umap") %>%
  ggsave(file = file.path(outdir, 
                          paste0("plot/01_merge_and_normalize_no_S521G/",
                                 "sctransform_merged_data_no_S521G", 
                                 "_genotype_umap.pdf")),
         plot = .,
         width = 12,
         height = 7)
DimPlot(data_sct_processed, 
        raster = FALSE, 
        order = TRUE, 
        label = FALSE, 
        group.by="genotype",
        split.by="genotype", 
        reduction = "umap") %>%
  ggsave(file = file.path(outdir, 
                          paste0("plot/01_merge_and_normalize_no_S521G/",
                                 "sctransform_merged_data_no_S521G", 
                                 "_genotype_split_umap.pdf")),
         plot = .,
         width = 12,
         height = 7)

DimPlot(data_sct_processed, 
        raster = FALSE, 
        order = TRUE, 
        label = FALSE, 
        group.by="sex", 
        reduction = "umap") %>%
  ggsave(file = file.path(outdir, 
                          paste0("plot/01_merge_and_normalize_no_S521G/",
                                 "sctransform_merged_data_no_S521G", 
                                 "_sex_umap.pdf")),
         plot = .,
         width = 12,
         height = 7)
DimPlot(data_sct_processed, 
        raster = FALSE, 
        order = TRUE, 
        label = FALSE, 
        group.by="sex",
        split.by="sex", 
        reduction = "umap",
        ncol = 2) %>%
  ggsave(file = file.path(outdir, 
                          paste0("plot/01_merge_and_normalize_no_S521G/",
                                 "sctransform_merged_data_no_S521G", 
                                 "_sex_split_umap.pdf")),
         plot = .,
         width = 12,
         height = 7)

DimPlot(data_sct_processed, 
        raster = FALSE, 
        order = TRUE, 
        label = FALSE, 
        group.by="date_of_birth", 
        reduction = "umap") %>%
  ggsave(file = file.path(outdir, 
                          paste0("plot/01_merge_and_normalize_no_S521G/",
                                 "sctransform_merged_data_no_S521G", 
                                 "_date_of_birth_umap.pdf")),
         plot = .,
         width = 12,
         height = 7)
DimPlot(data_sct_processed, 
        raster = FALSE, 
        order = TRUE, 
        label = FALSE, 
        group.by="date_of_birth",
        split.by="date_of_birth", 
        reduction = "umap",
        ncol = 2) %>%
  ggsave(file = file.path(outdir, 
                          paste0("plot/01_merge_and_normalize_no_S521G/",
                                 "sctransform_merged_data_no_S521G", 
                                 "_date_of_birth_split_umap.pdf")),
         plot = .,
         width = 12,
         height = 48)

DimPlot(data_sct_processed, 
        raster = FALSE, 
        order = TRUE, 
        label = FALSE, 
        group.by="age_at_perfusion", 
        reduction = "umap") %>%
  ggsave(file = file.path(outdir, 
                          paste0("plot/01_merge_and_normalize_no_S521G/",
                                 "sctransform_merged_data_no_S521G", 
                                 "_age_at_perfusion_umap.pdf")),
         plot = .,
         width = 12,
         height = 7)
DimPlot(data_sct_processed, 
        raster = FALSE, 
        order = TRUE, 
        label = FALSE, 
        group.by="age_at_perfusion",
        split.by="age_at_perfusion", 
        reduction = "umap",
        ncol = 2) %>%
  ggsave(file = file.path(outdir, 
                          paste0("plot/01_merge_and_normalize_no_S521G/",
                                 "sctransform_merged_data_no_S521G", 
                                 "_age_at_perfusion_split_umap.pdf")),
         plot = .,
         width = 12,
         height = 29)

DimPlot(data_sct_processed, 
        raster = FALSE, 
        order = TRUE, 
        label = FALSE, 
        group.by="date_perfused",
        reduction = "umap") %>%
  ggsave(file = file.path(outdir, 
                          paste0("plot/01_merge_and_normalize_no_S521G/",
                                 "sctransform_merged_data_no_S521G", 
                                 "_date_perfused_umap.pdf")),
         plot = .,
         width = 12,
         height = 7)
DimPlot(data_sct_processed, 
        raster = FALSE, 
        order = TRUE, 
        label = FALSE, 
        group.by="date_perfused",
        split.by="date_perfused", 
        reduction = "umap",
        ncol = 2) %>%
  ggsave(file = file.path(outdir, 
                          paste0("plot/01_merge_and_normalize_no_S521G/",
                                 "sctransform_merged_data_no_S521G", 
                                 "_date_perfused_split_umap.pdf")),
         plot = .,
         width = 12,
         height = 36)

DimPlot(data_sct_processed, 
        raster = FALSE, 
        order = TRUE, 
        label = FALSE, 
        group.by="date_of_nuclear_isolation", 
        reduction = "umap") %>%
  ggsave(file = file.path(outdir, 
                          paste0("plot/01_merge_and_normalize_no_S521G/",
                                 "sctransform_merged_data_no_S521G", 
                                 "_date_of_nuclear_isolation_umap.pdf")),
         plot = .,
         width = 12,
         height = 7)
DimPlot(data_sct_processed, 
        raster = FALSE, 
        order = TRUE, 
        label = FALSE, 
        group.by="date_of_nuclear_isolation",
        split.by="date_of_nuclear_isolation", 
        reduction = "umap",
        ncol = 2) %>%
  ggsave(file = file.path(outdir, 
                          paste0("plot/01_merge_and_normalize_no_S521G/",
                                 "sctransform_merged_data_no_S521G", 
                                 "_date_of_nuclear_isolation_split_umap.pdf")),
         plot = .,
         width = 12,
         height = 15)

FeaturePlot(data_sct_processed, 
            raster = FALSE, 
            order = TRUE, 
            label = FALSE, 
            features="nCount_RNA", 
            reduction = "umap") %>%
  ggsave(file = file.path(outdir, 
                          paste0("plot/01_merge_and_normalize_no_S521G/",
                                 "sctransform_merged_data_no_S521G", 
                                 "_nCount_RNA_umap.pdf")),
         plot = .,
         width = 12,
         height = 7)

FeaturePlot(data_sct_processed, 
            raster = FALSE, 
            order = TRUE, 
            label = FALSE, 
            features="nFeature_RNA", 
            reduction = "umap") %>%
  ggsave(file = file.path(outdir, 
                          paste0("plot/01_merge_and_normalize_no_S521G/",
                                 "sctransform_merged_data_no_S521G", 
                                 "_nFeature_RNA_umap.pdf")),
         plot = .,
         width = 12,
         height = 7)

FeaturePlot(data_sct_processed, 
            raster = FALSE, 
            order = TRUE, 
            label = FALSE, 
            features="percent.mt", 
            reduction = "umap") %>%
  ggsave(file = file.path(outdir, 
                          paste0("plot/01_merge_and_normalize_no_S521G/",
                                 "sctransform_merged_data_no_S521G", 
                                 "_percent.mt_umap.pdf")),
         plot = .,
         width = 12,
         height = 7)

FeaturePlot(data_sct_processed, 
            raster = FALSE, 
            order = TRUE, 
            label = FALSE, 
            features="nCount_SCT", 
            reduction = "umap") %>%
  ggsave(file = file.path(outdir, 
                          paste0("plot/01_merge_and_normalize_no_S521G/",
                                 "sctransform_merged_data_no_S521G", 
                                 "_nCount_SCT_umap.pdf")),
         plot = .,
         width = 12,
         height = 7)

FeaturePlot(data_sct_processed, 
            raster = FALSE, 
            order = TRUE, 
            label = FALSE, 
            features="nFeature_SCT", 
            reduction = "umap") %>%
  ggsave(file = file.path(outdir, paste0("plot/01_merge_and_normalize_no_S521G/",
                                         "sctransform_merged_data_",
                                         "no_S521G", 
                                         "_nFeature_SCT_umap.pdf")),
         plot = .,
         width = 12,
         height = 7)

print("***** Normalization completed! *****")


########################## END ########################## 

