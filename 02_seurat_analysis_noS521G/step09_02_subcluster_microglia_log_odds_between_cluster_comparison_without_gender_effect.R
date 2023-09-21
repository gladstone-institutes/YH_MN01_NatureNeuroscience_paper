### Between -cluster comparison
#For this analyses, we will use the **lme4** package in R to fit generalized 
#linear mixed effects models. We are going the model the change in the chance 
#(or more formally the odds) of cells from a given mouse belonging to a given 
#cluster from the 3 genotypes to the fE4 genotype. The random effects part
#of these models captures the inherent correlation between the cells coming 
#from the same mouse.

#run locally on Ayushi's laptop

library(lme4)
library(dplyr)
require(magrittr)
library(Seurat)
library(ggplot2)

#set all directory paths
basedir <- paste0("~/Dropbox (Gladstone)/YH_MN02/",
                  "seurat_analysis_no_S521G/")
indir <- paste0(basedir, 
                "06_subcluster_microglia_clusters_17_19/")
outdir <- paste0(basedir, 
                 "09_subclusters_log_odds_calculation/microglia_without_gender_effect/")
if(!(dir.exists(outdir))){
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE) 
}
setwd(outdir)

#load the Seurat object
pca_dim <- 15
cluster_res <- 0.9
dat <- readRDS(paste0(indir,
                      "microglia_data_post_subcluster_no_S521G_pcadims_",
                      pca_dim,
                      "_res_",
                      cluster_res,
                      ".rds"))

##the number of cells in each mouse in each cluster
all_metadata <- dat[[]]

#make a table of the required data
smp_cluster_counts <- unique(all_metadata %>%
                               group_by(new_samplenames) %>%
                               mutate(total_numbers_of_cells_per_sample = 
                                        n()) %>%
                               group_by(seurat_clusters, .add=TRUE) %>%
                               mutate(number_of_cells_per_sample_in_cluster = 
                                        n()) %>%
                               select(new_samplenames,
                                      genotype,
                                      seurat_clusters,
                                      total_numbers_of_cells_per_sample,
                                      number_of_cells_per_sample_in_cluster,
                                      sex))
colnames(smp_cluster_counts)[1:3] <- c("sample_id","animal_model","cluster_id")

smp_cluster_counts <- smp_cluster_counts[order(smp_cluster_counts$cluster_id),]
write.csv(smp_cluster_counts, 
          file = "counts_per_sample_per_microglia_subcluster.csv",
          row.names = FALSE)

counts <- read.csv("counts_per_sample_per_microglia_subcluster.csv", header = T)

pheno <- counts %>%
  select(sample_id, animal_model, total_numbers_of_cells_per_sample, sex) %>%
  unique()
rownames(pheno) <- seq(1,nrow(pheno))
Clusters <- unique(counts$cluster_id)

##function to estimate the change in the odds of cluster membership from the E4 to the other genotypes
estimateCellStateChange <- function(k, counts, pheno) {
  require(lme4)
  require(gdata)
  print(paste("Cluster", k))
  cluster_counts <- counts %>% 
    filter(cluster_id == k)
  cluster_counts %<>% merge(., pheno, all.y=TRUE)
  cluster_counts[is.na(cluster_counts$number_of_cells_per_sample_in_cluster),
                 "number_of_cells_per_sample_in_cluster"] <- 0
  
  cluster_counts %<>% 
    arrange(animal_model) %>% 
    mutate(proportion=
             number_of_cells_per_sample_in_cluster/total_numbers_of_cells_per_sample)
  ##plot proportion of cells per genotype
  pdf(paste0("proportion_of_cells_per_genotype_microglia_subcluster_",
             k,
             ".pdf"))
  print(ggplot(cluster_counts, 
               aes(x=animal_model, 
                   y=((number_of_cells_per_sample_in_cluster+0.01)/total_numbers_of_cells_per_sample))) +
          geom_boxplot(outlier.color = NA) +
          geom_jitter(position=position_jitter(0.2)) +
          scale_y_log10() +
          ylab("proportion of cells per genotype")+
          ggtitle(paste0("Microglia sub cluster ",k)))
  dev.off()
  
  ##plot proportion of cells per genotype and grouped by gender
  pdf(paste0("proportion_of_cells_per_genotype_groupby_sex_microglia_subcluster_",
             k,
             ".pdf"))
  print(ggplot(cluster_counts, 
               aes(x=animal_model, 
                   y=((number_of_cells_per_sample_in_cluster+0.01)/total_numbers_of_cells_per_sample),
                   color=sex)) +
          geom_boxplot(outlier.shape=NA) +
          geom_point(position=position_jitterdodge()) +
          scale_y_log10() +
          ylab("proportion of cells per genotype")+
          ggtitle(paste0("Microglia sub cluster ",k)))
  dev.off()
  
  cluster_counts %<>% mutate(animal_model = as.factor(animal_model))
  cluster_counts$animal_model <- relevel(cluster_counts$animal_model, ref="fE4")
  
  formula1=as.formula(paste0("cbind(number_of_cells_per_sample_in_cluster, ",
                             "total_numbers_of_cells_per_sample - ",
                             "number_of_cells_per_sample_in_cluster) ~ ",
                             "(1 | sample_id) + animal_model"))
  glmerFit <- glmer(formula1 ,data = cluster_counts, family = binomial, nAGQ=10)
  
  sglmerFit1 <- summary(glmerFit)
  TempRes1 <- (sglmerFit1$coefficients[-1,])
  
  return(TempRes1)
}

#run the log odds function for all clusters
ClusterRes <- sapply(Clusters, estimateCellStateChange, counts, pheno)
#reformat the results table
ClusterRes %<>% 
  as.data.frame() %>% 
  t() 
row.names(ClusterRes) <-  paste0("Cluster", Clusters)
ClusterRes <- data.frame(ClusterRes)
colnames(ClusterRes)[c(1:6, 10:12)] <- c("logOddsRatio_fE3_vs_fE4",
                                         "logOddsRatio_fE4R_S_vs_fE4",
                                         "logOddsRatio_fE4S_S_vs_fE4",
                                         "standardError_fE3_vs_fE4",
                                         "standardError_fE4R_S_vs_fE4",
                                         "standardError_fE4S_S_vs_fE4",
                                         "pvalue-fE3",
                                         "pvalue-fE4R_S",
                                         "pvalue-fE4S_S")
# make a vector of all p-values and p.adjust for all p-values together
p.adjust_all <- p.adjust(c(ClusterRes$`pvalue-fE3`, 
                           ClusterRes$`pvalue-fE4R_S`,
                           ClusterRes$`pvalue-fE4S_S`), method = "BH")
##perform multiple-testing correction
ClusterRes[,"p.adjust-fE3"] = p.adjust_all[1:nrow(ClusterRes)]
ClusterRes[,"p.adjust-fE4R_S"] = p.adjust_all[(nrow(ClusterRes) + 1):(nrow(ClusterRes)*2)]
ClusterRes[,"p.adjust-fE4S_S"] = p.adjust_all[(nrow(ClusterRes)*2 + 1):length(p.adjust_all)]

##output the results
ClusterRes <- ClusterRes[,!(colnames(ClusterRes) %in% c("X7","X8","X9"))]
print(ClusterRes)
write.csv(ClusterRes, 
          file = "log_odds_ratio_without_gender_effect_per_genotype_per_microglia_subcluster.csv")

#make boxplots for each cluster
#get data in long form
library(ggplot2)
ClusterRes$cluster_id <- rownames(ClusterRes)
x1 <- ClusterRes[,c("cluster_id", 
                    "logOddsRatio_fE3_vs_fE4", 
                    "standardError_fE3_vs_fE4",
                    "pvalue-fE3",
                    "p.adjust-fE3")]
x1$genotype <- "fE3"
colnames(x1) <- c("cluster_id","logOddsRatio","standardError","pvalue","p.adjust","genotype")
x2 <- ClusterRes[,c("cluster_id", 
                    "logOddsRatio_fE4R_S_vs_fE4", 
                    "standardError_fE4R_S_vs_fE4",
                    "pvalue-fE4R_S",
                    "p.adjust-fE4R_S")]
x2$genotype <- "fE4R_S"
colnames(x2) <- c("cluster_id","logOddsRatio","standardError","pvalue","p.adjust","genotype")
x3 <- ClusterRes[,c("cluster_id", 
                    "logOddsRatio_fE4S_S_vs_fE4", 
                    "standardError_fE4S_S_vs_fE4",
                    "pvalue-fE4S_S",
                    "p.adjust-fE4S_S")]
x3$genotype <- "fE4S_S"
colnames(x3) <- c("cluster_id","logOddsRatio","standardError","pvalue","p.adjust","genotype")

ClusterRes_plot <- rbind(x1,x2,x3)
rownames(ClusterRes_plot) <- seq(1,nrow(ClusterRes_plot))
ClusterRes_plot$cluster_id <- factor(ClusterRes_plot$cluster_id,
                                     levels = paste0("Cluster",seq(1,max(Clusters))))

#make the box plot for each cluster
#add label for p.adjusted < 0.05
label.df <- ClusterRes_plot[ClusterRes_plot$p.adjust < 0.05,]
label.df$logOddsRatio <- 4.3
label.df2 <- ClusterRes_plot[ClusterRes_plot$pvalue < 0.05,]
label.df2$logOddsRatio <- 3.9
pdf("log_odds_boxplot_without_gender_effect_microglia_subcluster_with_padj.pdf",
    height = 20,
    width = 35
)
ggplot(ClusterRes_plot, 
       aes(x=genotype, y=logOddsRatio, fill=genotype)) +
  geom_hline(yintercept=0) +
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=logOddsRatio-standardError, ymax=logOddsRatio+standardError), width=.2,
                position=position_dodge(.9)) +
  #scale_x_discrete(labels= c("fE3","fE4R_S","fE4S_S")) +
  facet_wrap(~ cluster_id) +
  theme(legend.position = "none",
        text = element_text(size = 28)
  ) +
  geom_text(data = label.df, label = "*",
            size = 15) +
  xlab("") +
  ylab("log odds ratio: animal model vs fE4")
dev.off()

pdf("log_odds_boxplot_without_gender_effect_microglia_subcluster_with_pvalue_and_padj.pdf",
    height = 20,
    width = 35
)
ggplot(ClusterRes_plot, 
       aes(x=genotype, y=logOddsRatio, fill=genotype)) +
  geom_hline(yintercept=0) +
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=logOddsRatio-standardError, ymax=logOddsRatio+standardError), width=.2,
                position=position_dodge(.9)) +
  #scale_x_discrete(labels= c("fE3","fE4R_S","fE4S_S")) +
  facet_wrap(~ cluster_id) +
  theme(legend.position = "none",
        text = element_text(size = 28)
  ) +
  geom_text(data = label.df, label = "*",
            size = 15) +
  geom_text(data = label.df2[label.df2$pvalue < 0.001,], label = "***",
            col = "gray", #hjust = -0.5,
            size = 15) +
  geom_text(data = label.df2[label.df2$pvalue > 0.001 & label.df2$pvalue < 0.01,], label = "**",
            col = "gray", #hjust = -0.5,
            size = 15) +
  geom_text(data = label.df2[label.df2$pvalue > 0.01 & label.df2$pvalue < 0.05,], label = "*",
            col = "gray", #hjust = -0.5,
            size = 15) +
  xlab("") +
  ylab("log odds ratio: animal model vs fE4")
dev.off()

#add session info
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")


#################### END ####################
