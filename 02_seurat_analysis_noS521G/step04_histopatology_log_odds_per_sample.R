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





