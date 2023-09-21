#Quality control metrics per cluster
#supplementary figure

#run locally on Ayushi's laptop

#load required packages
library(Seurat)
library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(cowplot)

#set all directory paths
basedir <- "~/Dropbox (Gladstone)/YH_MN02/"
indir <- paste0(basedir, 
                "seurat_analysis_no_S521G/02_clustering_no_S521G/")
outdir <- paste0(basedir,
                 "extended_data_plots/")

#load the Seurat object
dat <- readRDS(paste0(indir,
                      "sct_data_no_S521G_post_cluster_pcadims_15_res_0.7.rds"))

##the number of cells in each mouse in each cluster
all_metadata <- dat[[]]

#make a table of the required data
avg_qc_features_per_cluster <- all_metadata %>% group_by(seurat_clusters) %>% 
  summarise(total_cells = n(),
            avg_ngene = mean(nFeature_RNA),
            sd_ngene = sd(nFeature_RNA),
            se_ngene = sd(nFeature_RNA) / sqrt(length(nFeature_RNA)),
            avg_nUMI = mean(nCount_RNA),
            sd_nUMI = sd(nCount_RNA),
            se_nUMI = sd(nCount_RNA) / sqrt(length(nCount_RNA)),
            avg_mt =  mean(percent.mt),
            sd_mt =  sd(percent.mt),
            se_mt = sd(percent.mt) / sqrt(length(percent.mt))) %>% 
  as.data.frame()

##########
#panel b
##########
# Create plot with legend
ggp1_legend <- ggplot(avg_qc_features_per_cluster, 
                      aes(x=seurat_clusters, y=total_cells, 
                          fill=seurat_clusters)) +
  geom_bar(position=position_dodge(), stat="identity") + 
  guides(fill=guide_legend(ncol=4)) + 
  scale_fill_discrete(name = "Cluster Identity", 
                      labels = c("1 - Ex Neuron (Granule Cell)",
                                 "2 - Oligodendrocyte",
                                 "3 - Oligodendrocyte",
                                 "4 - Ex Neuron (CA1 Neuron)",
                                 "5 - In Neuron",
                                 "6 - Ex Neuron (CA2/CA3)",
                                 "7 - Ex Neuron (CA2/CA3)",
                                 "8 - Ex Neuron",
                                 "9 - Oligodendrocyte",
                                 "10 - In Neuron",
                                 "11 - Ex Neuron",
                                 "12 - In Neuron",
                                 "13 - Astrocyte",
                                 "14 - OPC",
                                 "15 - Ex Neuron",
                                 "16 - Choroid Plexus",
                                 "17 - Microglia",
                                 "18 - Ex Neuron",
                                 "19 - Microglia",
                                 "20 - Ex Neuron",
                                 "21 - Unknown",
                                 "22 - Ex Neuron",
                                 "23 - Ex Neuron",
                                 "24 - Ex Neuron",
                                 "25 - Unknown",
                                 "26 - Ex Neuron",
                                 "27 - Unknown",
                                 "28 - Ex Neuron",
                                 "29 - Unknown",
                                 "30 - Ex Neuron",
                                 "31 - In Neuron",
                                 "32 - Unknown",
                                 "33 - Ex Neuron",
                                 "34 - In Neuron",
                                 "35 - Unknown",
                                 "36 - Astrocyte",
                                 "37 - Oligodendrocyte",
                                 "38 - Unknownn")) +
  theme(text = element_text(size = 26))
# Create user-defined function, which extracts legends from ggplots
extract_legend <- function(my_ggp) {
  step1 <- ggplot_gtable(ggplot_build(my_ggp))
  step2 <- which(sapply(step1$grobs, function(x) x$name) == "guide-box")
  step3 <- step1$grobs[[step2]]
  return(step3)
}

# Apply user-defined function to extract legend
shared_legend <- extract_legend(ggp1_legend)

##########
#panel c
##########
p1 <- ggplot(avg_qc_features_per_cluster, 
             aes(x=seurat_clusters, y=total_cells, fill=seurat_clusters)) +
  geom_bar(position=position_dodge(), stat="identity") + 
  scale_x_discrete(expand = c(0.04,0.04))+
  theme_classic() + theme(
    panel.border = element_rect(colour = "black", fill=NA, size=2),
    text = element_text(size = 20),
    legend.position = "none") + 
  ggtitle("Number of Cells per Cluster") +
  xlab("Cell Cluster") + 
  ylab("Number of Cells")

##########
#fig c
##########
# p2 <- ggplot(avg_qc_features_per_cluster, 
#              aes(x=seurat_clusters, y=avg_ngene, fill=seurat_clusters)) + 
#   geom_bar(stat="identity", position=position_dodge()) +
#   geom_errorbar(aes(ymin=avg_ngene-se_ngene, ymax=avg_ngene+se_ngene), 
#                 width=.2,
#                 position=position_dodge(.9)) + 
#   scale_x_discrete(expand = c(0.04,0.04))+
#   theme_classic() + theme(
#     panel.border = element_rect(colour = "black", fill=NA, size=2),
#     text = element_text(size = 20),
#     legend.position = "none") + 
#   ggtitle("Average Genes per Cells by Cluster") +
#   xlab("Cell Cluster") + 
#   ylab("Average Number of Genes") 
p2 <- ggplot(all_metadata, 
             aes(x=seurat_clusters, y=nFeature_RNA, fill=seurat_clusters)) + 
  geom_boxplot(outlier.color = NA) + stat_boxplot(geom='errorbar') +
  geom_point(alpha =0.2, size = 0.3) +
  scale_x_discrete(expand = c(0.04,0.04))+
  theme_classic() + theme(
    panel.border = element_rect(colour = "black", fill=NA, size=2),
    text = element_text(size = 20),
    legend.position = "none") + 
  ggtitle("Number of Genes per Cell by Cluster") +
  xlab("Cell Cluster") + 
  ylab("Number of Genes per Cell") 

##########
#panel d
##########
# p3 <- ggplot(avg_qc_features_per_cluster, 
#              aes(x=seurat_clusters, y=avg_nUMI, fill=seurat_clusters)) + 
#   geom_bar(stat="identity", position=position_dodge()) +
#   geom_errorbar(aes(ymin=avg_nUMI-se_nUMI, ymax=avg_nUMI+se_nUMI), 
#                 width=.2,
#                 position=position_dodge(.9)) + 
#   scale_x_discrete(expand = c(0.04,0.04))+
#   theme_classic() + theme(
#     panel.border = element_rect(colour = "black", fill=NA, size=2),
#     text = element_text(size = 20),
#     legend.position = "none") + 
#   ggtitle("Average nUMI per Cells by Cluster") +
#   xlab("Cell Cluster") + 
#   ylab("Average nUMI")
p3 <- ggplot(all_metadata, 
             aes(x=seurat_clusters, y=nCount_RNA, fill=seurat_clusters)) + 
  geom_boxplot(outlier.color = NA) + stat_boxplot(geom='errorbar') +
  geom_point(alpha =0.2, size = 0.3) +
  scale_x_discrete(expand = c(0.04,0.04))+
  theme_classic() + theme(
    panel.border = element_rect(colour = "black", fill=NA, size=2),
    text = element_text(size = 20),
    legend.position = "none") + 
  ggtitle("nUMI per Cell by Cluster") +
  xlab("Cell Cluster") + 
  ylab("nUMI per Cell") 


##########
#panel e
##########
# p4 <- ggplot(avg_qc_features_per_cluster, 
#              aes(x=seurat_clusters, y=avg_mt, fill=seurat_clusters)) + 
#   geom_bar(stat="identity", position=position_dodge()) +
#   geom_errorbar(aes(ymin=avg_mt-se_mt, ymax=avg_mt+se_mt), 
#                 width=.2,
#                 position=position_dodge(.9)) + 
#   scale_x_discrete(expand = c(0.04,0.04))+
#   theme_classic() + theme(
#     panel.border = element_rect(colour = "black", fill=NA, size=2),
#     text = element_text(size = 20),
#     legend.position = "none") + 
#   ggtitle("Average % Mitochondrial Genes per Cells by Cluster") +
#   xlab("Cell Cluster") + 
#   ylab("Average % Mito Genes")
p4 <- ggplot(all_metadata, 
             aes(x=seurat_clusters, y=percent.mt, fill=seurat_clusters)) + 
  geom_boxplot(outlier.color = NA) + stat_boxplot(geom='errorbar') +
  geom_point(alpha =0.2, size = 0.3) +
  scale_x_discrete(expand = c(0.04,0.04))+
  theme_classic() + theme(
    panel.border = element_rect(colour = "black", fill=NA, size=2),
    text = element_text(size = 20),
    legend.position = "none") + 
  ggtitle("% Mitochondrial Genes per Cells by Cluster") +
  xlab("Cell Cluster") + 
  ylab("% Mito Genes per Cell")


# Draw plots with shared legend
#with standard error
pdf(paste0(outdir, "quality_control_supplemental_figure_boxplots.pdf"),
    width = 25, height = 26)
grid.arrange(shared_legend,
             arrangeGrob(p1, p2, p3,p4, ncol = 2),
             heights=c(2, 10)) 
dev.off()

print("********** Script completed! **********")

################## END ################## 
