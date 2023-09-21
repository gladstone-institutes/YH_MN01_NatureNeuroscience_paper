library(Seurat)
library(ggplot2)
library(dplyr)
require(emmeans)

sc_dat <- readRDS(paste0("~/Dropbox (Gladstone)/YH_MN02/seurat_analysis_no_S521G/",
                         "02_clustering_no_S521G/sct_data_no_S521G_post_cluster_pcadims_15_res_0.7.rds"))

outdir <- paste0("~/Dropbox (Gladstone)/YH_MN02/seurat_analysis_no_S521G/",
                 "revision_plots")
dir.create(outdir)
setwd(outdir)

pdf("violin_plot_apoe_SCTassay_v1.pdf", width = 20)
print(VlnPlot(sc_dat, features = "hapoE-transgene",
              group.by = "seurat_clusters",
              split.by = "genotype", 
              cols = c("yellow","pink","green","red"),
              raster = F))
dev.off()

pdf("violin_plot_apoe_RNAassay_v2.pdf", width = 20)
print(VlnPlot(sc_dat, features = "hapoE-transgene", assay = "RNA",
              group.by = "seurat_clusters",
              split.by = "genotype", 
              cols = c("yellow","pink","green","red"),
              raster = F))
dev.off()

pdf("violin_plot_apoe_SCTassay_neurons_only.pdf", width = 20)
print(VlnPlot(sc_dat, features = "hapoE-transgene", assay = "SCT",
              idents = c(1, 4, 5, 6, 7, 8, 10, 11, 12, 15, 18, 20, 22, 23, 24, 26, 28, 30, 31, 33, 34),
              group.by = "seurat_clusters",
              split.by = "genotype", 
              cols = c("yellow","pink","green","red"),
              raster = F))
dev.off()

pdf("violin_plot_apoe_SCTassay_glia_cell_clusters.pdf", width = 10)
print(VlnPlot(sc_dat, features = "hapoE-transgene", assay = "SCT",
              idents = c(13, 17, 19, 36, 38),
              group.by = "seurat_clusters",
              split.by = "genotype",
              cols = c("yellow","pink","green","red"),
              raster = F))
dev.off()

#ridge plot
sc_dat$genotype_label <- ifelse(sc_dat$genotype == "fE3", "PS19-E3",
                                ifelse(sc_dat$genotype == "fE4", "PS19-E4",
                                       ifelse(sc_dat$genotype == "fE4R_S", "PS19-E4-R/S", "PS19-E4-S/S")))
cluster_genotype_levels <- expand.grid(c("PS19-E4", "PS19-E3", "PS19-E4-R/S", "PS19-E4-S/S"), levels(sc_dat$seurat_clusters)) %>% 
  mutate(reorder_levels = paste0(Var2, "_", Var1))
sc_dat$cluster_genotype <- paste0(sc_dat$seurat_clusters, "_", sc_dat$genotype_label)
sc_dat$cluster_genotype <- factor(sc_dat$cluster_genotype,
                                  levels = cluster_genotype_levels$reorder_levels)
pdf("ridge_plot_apoe_SCTassay_v1.pdf", width = 10, height = 25)
print(RidgePlot(sc_dat, features = "hapoE-transgene", 
                group.by = "cluster_genotype") + theme(legend.position = 'none'))
dev.off()

pdf("ridge_plot_apoe_SCTassay_neurons_only.pdf", width = 10, height = 15)
print(RidgePlot(sc_dat, features = "hapoE-transgene", 
                idents = c(1, 4, 5, 6, 7, 8, 10, 11, 12, 15, 18, 20, 22, 23, 24, 26, 28, 30, 31, 33, 34),
                group.by = "cluster_genotype") + theme(legend.position = 'none'))
dev.off()

pdf("ridge_plot_apoe_SCTassay_glia_cell_clusters.pdf", width = 10)
print(RidgePlot(sc_dat, features = "hapoE-transgene", idents = c(13, 17, 19, 36, 38), 
                group.by = "cluster_genotype") + theme(legend.position = 'none'))
dev.off()

#For the clusters 13, 16, 26, 36, 38, is it possible to compare the APOE expression levels 
#among the four genotypes to see whether the level in the S/S group is significantly different from other groups?
x <- data.frame(apoe_exp = FetchData(sc_dat, "hapoE-transgene", assay = "SCT", slot = "data")$`hapoE-transgene`,
                sample_names = sc_dat$new_samplenames,
                cluster = sc_dat$seurat_clusters,
                genotype = sc_dat$genotype)
x <- x[x$cluster %in% c(13, 16, 26, 36, 38),]
x <- x %>% group_by(cluster, sample_names, genotype) %>% summarise(mean_apoe_exp = mean(apoe_exp)) %>% as.data.frame()
x$cluster <- factor(x$cluster)
x$genotype <- factor(x$genotype)
x$genotype <- relevel(x$genotype, ref="fE4")

#fit a linear model
lmfit <- lm(mean_apoe_exp ~ cluster * genotype, data=x)
print(summary(lmfit))
contrast_means <- emmeans(lmfit, pairwise ~ genotype|cluster, adjust="none")
results <- as.data.frame(contrast_means$contrasts) %>% 
  mutate(adjusted.p.value=p.adjust(p.value, method = "BH"))
write.csv(results, "human_apoe_genotype_association_analysis_clusters_13_16_26_36_38.csv", row.names = F)

pdf("human_apoe_expression_across_genotypes_clusters_13_16_26_36_38_boxplot.pdf", width = 10)
print(ggplot(x, aes(x=cluster, y=mean_apoe_exp, fill=genotype)) + geom_boxplot(outlier.shape = NA) +
        geom_point(position=position_jitterdodge(jitter.width = 0.2,dodge.width = 0.75)) +
        theme_classic() + theme(text = element_text(size = 16)))
dev.off()


########################## END ########################## 
