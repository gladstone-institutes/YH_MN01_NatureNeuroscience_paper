 #to identify clusters/cell-types whose membership changes with histopathological parameters
rm(list=ls())
library(lme4)
library(dplyr)
require(magrittr)
library(Seurat)
library(ggplot2)

#set all directory paths
sourcedir<-"~/Dropbox (Gladstone)/YH_MN02/seurat_analysis_no_S521G/03_log_odds_calculation/"
histodir <- "~/Dropbox (Gladstone)/MT/Maxine/input/"
outdir <- "~/Dropbox (Gladstone)/YH_MN02/seurat_analysis_no_S521G/04_log_odds_histopathology_calculation/"

setwd("~/Dropbox (Gladstone)/YH_MN02/seurat_analysis_no_S521G/04_log_odds_histopathology_calculation/")


histop <- read.csv(paste0(histodir,"Maxine_Pathological_Data_of_Mice_snRNAseq_07.14.22.csv"), header = T, check.names=F)
counts <- read.csv(paste0(sourcedir,"counts_per_sample_per_cluster.csv"), header = T)
pheno <- counts %>%
  select(sample_id, animal_model, total_numbers_of_cells_per_sample) %>%
  unique()
Clusters <- unique(counts$cluster_id)

histopheno <- merge(histop, pheno, by="sample_id") 
head(histopheno)
dim(histopheno)

##function to estimate the change in the odds of cluster membership from the E4 to the other genotypes
estimateCellStateChange <- function(k, counts, histopheno, optimizer) {
  require(lme4)
  require(gdata)
  print(paste("Cluster", k))
  cluster_counts <- counts %>% 
    filter(cluster_id == k)
  cluster_counts %<>% merge(., histopheno, all.y=TRUE)
  cluster_counts$number_of_cells_per_sample_in_cluster[is.na(cluster_counts$number_of_cells_per_sample_in_cluster)] <- 0
  
  cluster_counts %<>% arrange(animal_model) %>% mutate(proportion=number_of_cells_per_sample_in_cluster/total_numbers_of_cells_per_sample)



  colnames(cluster_counts)[8:13] <- c("Hippocampus_Vol_mm3",  "prop_IBA1_Coverage_Area", "prop_CD68_Coverage_Area", "prop_GFAP_Quant_Coverage_Area", "prop_S100beta_Coverage_Area", "prop_AT8_Coverage_Area")
 cluster_counts %<>% mutate(animal_model = as.factor(animal_model))
  cluster_counts$animal_model <- relevel(cluster_counts$animal_model, ref="fE4")
 
TempRes<-NULL

for (ph in 8:(ncol(cluster_counts)-2)){
pdf(paste0(outdir,"proportion_of_cells_per_sample_cluster_",k,"_per_pheno_",colnames(cluster_counts[ph]),".pdf"))

p <- ggplot(cluster_counts, aes(((number_of_cells_per_sample_in_cluster+0.01)/total_numbers_of_cells_per_sample), cluster_counts[,ph]))
print(p + geom_point(aes(colour = factor(animal_model)) ) + scale_x_log10() + geom_smooth(aes(color = animal_model), method = "lm", se = FALSE, size=1) + xlab(paste0("proportion of cells per sample in cluster ", k)) + ylab(paste0(colnames(cluster_counts[ph]))))
dev.off()

 
    formula1=as.formula(paste0("cbind(number_of_cells_per_sample_in_cluster, ",
                             "total_numbers_of_cells_per_sample - ",
                             "number_of_cells_per_sample_in_cluster) ~ ",
                             "(1 | animal_model/sample_id) + cluster_counts[,ph] "))
  glmerFit <- glmer(formula1 ,data = cluster_counts, family = binomial, nAGQ=1)
  
  sglmerFit1 <- summary(glmerFit)
  TempRes1 <- (sglmerFit1$coefficients[-1,])
  
  #print(TempRes1)
      TempRes <- rbind(TempRes, TempRes1)

             
}  

  #rownames(TempRes) <- c(`Hippocampus_Vol_mm3`, `prop_IBA1_Coverage_Area`, `prop_CD68_Coverage_Area`,`prop_GFAP_Quant_Coverage_Area`,`prop_S100beta_Coverage_Area`,`prop_AT8_Coverage_Area`)

print(TempRes)

}


#run the log odds function for all clusters for each histopathological variable
ClusterRes <- sapply(Clusters, estimateCellStateChange, counts, histopheno)
#reformat the results table
ClusterRes %<>% 
  as.data.frame() %>% 
  t() 
row.names(ClusterRes) <-  paste0("Cluster", Clusters)
ClusterRes <- data.frame(ClusterRes)
colnames(ClusterRes)[c(1:12, 19:24)] <- c("logOddsRatio_for_unit_Hippocampus_Vol_mm3",
										  "logOddsRatio_for_unit_prop_IBA1_Coverage_Area",
										  "logOddsRatio_for_unit_prop_CD68_Coverage_Area",
										  "logOddsRatio_for_unit_prop_GFAP_Quant_Coverage_Area",
                                          "logOddsRatio_for_unit_prop_S100beta_Coverage_Area",
                                          "logOddsRatio_for_unit_prop_AT8_Coverage_Area",
                                          "StdErr_for_unit_Hippocampus_Vol_mm3",
                                          "StdErr_for_unit_prop_IBA1_Coverage_Area",
										  "StdErr_for_unit_prop_CD68_Coverage_Area",
										  "StdErr_for_unit_prop_GFAP_Quant_Coverage_Area",
                                          "StdErr_for_unit_prop_S100beta_Coverage_Area",
                                          "StdErr_for_unit_prop_AT8_Coverage_Area",
                                          "pvalue_for_unit_Hippocampus_Vol_mm3",
										  "pvalue_for_unit_prop_IBA1_Coverage_Area",
										  "pvalue_for_unit_prop_CD68_Coverage_Area",
										  "pvalue_for_unit_prop_GFAP_Quant_Coverage_Area",
                                          "pvalue_for_unit_prop_S100beta_Coverage_Area",
                                          "pvalue_for_unit_prop_AT8_Coverage_Area")
                                         
                                         
# make a vector of all p-values and p.adjust for all p-values together
p.adjust_all <- p.adjust(c(ClusterRes$`pvalue_for_unit_Hippocampus_Vol_mm3`, 
                           ClusterRes$`pvalue_for_unit_prop_IBA1_Coverage_Area`,
                           ClusterRes$`pvalue_for_unit_prop_CD68_Coverage_Area`,
                           ClusterRes$`pvalue_for_unit_prop_GFAP_Quant_Coverage_Area`,
                           ClusterRes$`pvalue_for_unit_prop_S100beta_Coverage_Area`,
                           ClusterRes$`pvalue_for_unit_prop_AT8_Coverage_Area`), method = "BH")
##perform multiple-testing correction
ClusterRes[,"p.adjust_for_unit_Hippocampus_Vol_mm3"] = p.adjust_all[1:nrow(ClusterRes)]
ClusterRes[,"p.adjust_for_unit_prop_IBA1_Coverage_Area"] = p.adjust_all[(nrow(ClusterRes) + 1):(nrow(ClusterRes)*2)]
ClusterRes[,"p.adjust_for_unit_prop_CD68_Coverage_Area"] = p.adjust_all[(nrow(ClusterRes)*2 + 1):(nrow(ClusterRes)*3)]
ClusterRes[,"p.adjust_for_unit_prop_GFAP_Quant_Coverage_Area"] = p.adjust_all[(nrow(ClusterRes)*3 + 1):(nrow(ClusterRes)*4)]
ClusterRes[,"p.adjust_for_unit_prop_S100beta_Coverage_Area"] = p.adjust_all[(nrow(ClusterRes)*4 + 1):(nrow(ClusterRes)*5)]
ClusterRes[,"p.adjust_for_unit_prop_AT8_Coverage_Area"] = p.adjust_all[(nrow(ClusterRes)*5 + 1):(length(p.adjust_all))]



##output the results
ClusterRes <- ClusterRes[,!(colnames(ClusterRes) %in% c("X13","X14","X15","X16","X17","X18"))]
print(ClusterRes)
write.csv(t(ClusterRes), file = paste0(outdir,"Maxine_log_odds_ratio_per_unit_per_histopathology_per_cluster_by_sample.csv"))





sourcedir <- "~/Dropbox (Gladstone)/YH_MN02/seurat_analysis_no_S521G/04_log_odds_histopathology_calculation/"
outdir <-  "~/Dropbox (Gladstone)/YH_MN02/seurat_analysis_no_S521G/04_log_odds_histopathology_calculation/"
setwd(sourcedir)
library(dplyr)
library(ggplot2)
library(dendsort)
library(pheatmap)
library(RColorBrewer)

histo_all<-read.csv(paste0(sourcedir,"Maxine_log_odds_ratio_per_unit_per_histopathology_per_cluster_by_sample.csv"))

histo<- histo_all  %>% slice_head(n = 6)
#hippocampal volume, AT8 area, Iba1 area, CD68 area, GFAP area, and S100b area
histo_ordered<- rbind(histo[c(1,6,2,3,4,5),])

#without rowmeans across clusters normalization
histo_no_norm <- data.frame(histo_parameters=histo_ordered[1], histo_ordered[2:ncol(histo_ordered)])
colnames(histo_no_norm)[1]<-"histo_parameters"


#with rowmeans across clusters normalization
histo_rowMeans_norm <- data.frame(histo_parameters=histo_ordered[1], histo_ordered[2:ncol(histo_ordered)] - rowMeans(histo_ordered[2:ncol(histo_ordered)]))
colnames(histo_rowMeans_norm)[1]<-"histo_parameters"

#subset with cluster of interest 1, 3, 4, 6, 7, 8, 9, 22, 26, 28.
histo_no_norm_cluster_of_interest <-  histo_no_norm %>% select(Cluster1,Cluster3,Cluster4,Cluster6,Cluster7,Cluster8, Cluster9, Cluster22,Cluster26,Cluster28)
histo_rowMeans_norm_cluster_of_interest <-  histo_rowMeans_norm %>% select(Cluster1,Cluster3,Cluster4,Cluster6,Cluster7,Cluster8, Cluster9, Cluster22,Cluster26,Cluster28)

#clustering and heatmap
#Without rowmeans normalization

#Elbow Method for finding the optimal number of clusters
set.seed(123)
# Compute and plot wss for k = 2 to k = 15.
k.max <- 5
data <- histo_no_norm_cluster_of_interest
set.seed(123)
wss <- sapply(1:k.max, 
              function(k){kmeans(data, k, nstart=50,iter.max = 15 )$tot.withinss})
wss

#plot(1:k.max, wss,
#     type="b", pch = 19, frame = FALSE, 
#      xlab="Number of clusters K",
#      ylab="Total within-clusters sum of squares")





##perform k-means clustering with k=3
kNClust <- 3
kmeanClust <- kmeans(histo_no_norm_cluster_of_interest,kNClust)

row.names(histo_no_norm_cluster_of_interest) <- histo_no_norm$histo_parameters
no_norm_histo_kmeans <- histo_no_norm_cluster_of_interest #[order(kmeanClust$cluster),]


kGaps <- vector(mode = "numeric")
kGaps[1] <- kmeanClust$size[1] + 1
for(i in 2:(kNClust-1)) {
  kGaps[i] <- kGaps[i-1] + kmeanClust$size[i]
}

simpleredbluecols = colorRampPalette(c("blue","white","red"))(400)

df <- data.frame(histo_parameters=histo_no_norm$histo_parameters)
row.names(df) <- histo_no_norm$histo_parameters





paletteLength <- 400
##center the color scale so that while represents 0 in the row-centered matrix
myBreaks <- c(seq(min(histo_no_norm_cluster_of_interest, na.rm = TRUE), 0, length.out=ceiling(paletteLength/2) + 1),
              seq(max(histo_no_norm_cluster_of_interest, na.rm = TRUE)/paletteLength, max(histo_no_norm_cluster_of_interest, na.rm=TRUE), length.out=floor(paletteLength/2)))


#for dendrogram
#mat_cluster_cols <- hclust(dist(t(no_norm_histo_kmeans)))
#sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))
#cmat_cluster_cols <- sort_hclust(mat_cluster_cols)
#mat_cluster_rows <- sort_hclust(hclust(dist(no_norm_histo_kmeans)))

cl<- data.frame(colnames(no_norm_histo_kmeans))
row.names(cl) <- colnames(no_norm_histo_kmeans)
no_norm_histo_kmeans



pdf(paste0(outdir,"Fig.5h_pathology_cl1_3_4_6_7_8_9_22_26_28.pdf"), height=6, width=12)
print(pheatmap(no_norm_histo_kmeans, cluster_rows=F, show_rownames=T, show_colnames=T, cluster_cols=F, treeheight_row=T,  gaps_row = F, annot_cols=df,annot_row=cl, color = simpleredbluecols, breaks = myBreaks, display_numbers = TRUE, number_format="%.3f",
         number_color = "black", 
         fontsize_number = 14))
dev.off()




