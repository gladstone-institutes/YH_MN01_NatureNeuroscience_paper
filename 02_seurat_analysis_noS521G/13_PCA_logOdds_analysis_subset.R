rm(list=ls())
logOddsdir <- "~/Dropbox (Gladstone)/YH_MN02/seurat_analysis_no_S521G/04_log_odds_histopathology_calculation/"
logOdds_astro_dir <- "~/Dropbox (Gladstone)/YH_MN02/seurat_analysis_no_S521G/10_subclusters_log_odds_histopathology_calculation/"
logOdds_micro_dir <- "~/Dropbox (Gladstone)/YH_MN02/seurat_analysis_no_S521G/10_subclusters_log_odds_histopathology_calculation/"

histodir <- "~/Dropbox (Gladstone)/MT/Maxine/input/"

outdir <- "~/Dropbox (Gladstone)/YH_MN02/seurat_analysis_no_S521G/13_logOdds_PCA/"
setwd(outdir)
library(factoextra)
library(FactoMineR)
library(dplyr)
library(magrittr)
library(RColorBrewer)




##logOdds histopathology all clusters
logOdds_histo=read.csv(paste0(logOddsdir,"Maxine_log_odds_ratio_per_unit_per_histopathology_per_cluster_by_sample.csv"))
logOdds_histo <- logOdds_histo[,c("X", "Cluster1", "Cluster6","Cluster7","Cluster9","Cluster28")]
##logOdds histopathology all astrocyte subclusters
logOdds_histo_astro=read.csv(paste0(logOdds_astro_dir,"astrocyte_log_odds_ratio_per_unit_per_histopathology_per_cluster_by_sample.csv"))
colnames(logOdds_histo_astro) <- paste("AS", colnames(logOdds_histo_astro), sep = "_")
logOdds_histo_astro <- logOdds_histo_astro[,c("AS_X", "AS_Cluster3", "AS_Cluster5","AS_Cluster7")]
##logOdds histopathology all microglia subclusters
logOdds_histo_micro=read.csv(paste0(logOdds_micro_dir,"microglia_log_odds_ratio_per_unit_per_histopathology_per_cluster_by_sample.csv"))
colnames(logOdds_histo_micro) <- paste("MG", colnames(logOdds_histo_micro), sep = "_")
logOdds_histo_micro <- logOdds_histo_micro[,c("MG_X", "MG_Cluster2", "MG_Cluster8","MG_Cluster11")]


#merge and extract all subclusters clusters - logOdds only
logOdds_histo_clusters_of_interest_all_sub <- merge(merge(logOdds_histo, logOdds_histo_astro, by.x="X", by.y="AS_X") , logOdds_histo_micro,  by.x="X",by.y="MG_X") 
#all pheno
logOdds_histo_clusters_of_interest_all_sub %<>% slice_head(n = 6)
colnames(logOdds_histo_clusters_of_interest_all_sub)[1] <- "logOddsHistopathology"






#all subclusters

t_histo<-t(logOdds_histo_clusters_of_interest_all_sub)

histo_log<-data.frame(cluster=colnames(logOdds_histo_clusters_of_interest_all_sub)[-1], as.numeric(t_histo[-1,1]),as.numeric(t_histo[-1,2]), as.numeric(t_histo[-1,3]), as.numeric(t_histo[-1,4]), as.numeric(t_histo[-1,5]), as.numeric(t_histo[-1,6]))
histo_log$cluster<-factor(histo_log$cluster, levels = histo_log$cluster)
colnames(histo_log)<-c("cluster",logOdds_histo_clusters_of_interest_all_sub$logOddsHistopathology)
histo_log_PCA<-histo_log[,-1]
histo_log_PCA<-data.frame(histo_log_PCA)
colnames(histo_log_PCA) <- c("logOddsRatio_for_unit_Hippocampus_Vol_mm3" ,     "logOddsRatio_for_unit_prop_AT8_Coverage_Area" , "logOddsRatio_for_unit_prop_CD68_Coverage_Area" ,  "logOddsRatio_for_unit_prop_GFAP_Quant_Coverage_Area", "logOddsRatio_for_unit_prop_IBA1_Coverage_Area","logOddsRatio_for_unit_prop_S100beta_Coverage_Area" )

res.pca.norm.histo<- prcomp(histo_log_PCA,  scale=T, center = T, rank.=2)

pca_histo_log<-data.frame(cluster=histo_log$cluster,res.pca.norm.histo$x)

#
#1) For astrocytes without gender effect, heatmap of subclusters 3, 5, and 7 with all pathologies.
#2) For microglia without gender effect, heatmap of subclusters 2, 8, 11 with all pathologies.
#3) For microglia with gender effect, heatmap of subclusters 2, 8, 11 with all pathologies.


library(RColorBrewer)
#histopatohlogy
getPalette = colorRampPalette(brewer.pal(8, "Set3"))
nb.cols <-6


#rownames(res.pca.norm.histo) <- pca_histo_log$cluster
#pdf(paste0(outdir,"PCA_logodds_histopathology_cl17_25_all_subclusters_refined.pdf"))


p <- fviz_pca_ind(res.pca.norm.histo,
             geom.ind = c("point"),
                pointshape = 21,repel = T,
                pointsize = 5.5, title = "", fill.ind = pca_histo_log$cluster, mean.point = FALSE, addEllipses =F, legend.title = "") + 
                guides(x.sec = "axis", y.sec = "axis")+
                scale_fill_manual(values = c("gold", "pink2","forestgreen","limegreen","navy",getPalette(nb.cols)))+
                #geom_text(aes(label=ifelse(data.frame(res.pca.norm.histo$x)$PC1>2 & data.frame(res.pca.norm.histo$x)$PC2>0,as.character(pca_histo_log$cluster),'')),position = position_dodge(0.9))
                annotate("text", label = "MG cluster 11", x = -2.8 ,y = -1.15, size = 4) +
                annotate("text", label = "MG cluster 8", x = 2.55 ,y = -0.3, size = 4) +
                annotate("text", label = "MG cluster 2", x = -2.9 ,y = 0.36, size = 4) +
                annotate("text", label = "AS cluster 3",x = 0.7 ,y = 0.24, size = 4) +
                annotate("text", label = "AS cluster 5", x = 2.2 ,y = 0.14, size = 4) +
                annotate("text", label = "AS cluster 7", x = 1.6 ,y = -0.38, size = 4) +
                 annotate("text", label = "Cluster 1", x = 0.3 ,y = 0.37, size = 4) +
                 annotate("text", label = "Cluster 6", x = -1.8 ,y = 0.23, size = 4) +
                 annotate("text", label = "Cluster 7", x = -0.5 ,y = 0.21, size = 4)+
                 annotate("text", label = "Cluster 9", x = 2.5 ,y = -0.082, size = 4)+
                 annotate("text", label = "Cluster 28", x = -1.4 ,y = 0.33, size = 4)
  ggsave(file = file.path(outdir, paste0("PCA_logodds_histopathology_cl1_6_7_9_28_astro_3_5_7_micro_2_8_11.pdf")),
         plot = p,
         width = 10,
         height = 10)








