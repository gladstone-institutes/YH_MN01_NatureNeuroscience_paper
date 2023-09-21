#run locally on Ayushi's laptop
#Goal: volcano plots of the enriched genes using EnhancedVolcano package

#load required packages
library(EnhancedVolcano)
library(tools)
library(dplyr)

#set all directory paths
indir <- paste0("~/Dropbox (Gladstone)/YH_MN02/seurat_analysis_no_S521G")
outdir <- file.path(indir,"15_volcano_plots/")

if(!dir.exists(outdir)){
  dir.create(outdir)
}
setwd(indir)

############
#load the data
############
#get the list of all the de gene files
de_files <- c("07_de_genes_no_S521G/de_genes_cluster_9_vs_2_post_clustering_sct_no_S521G_pcadims_15_res_0.7.csv",
              "07_de_genes_no_S521G/de_genes_cluster_9_PS19-E4-S_S_vs_PS19-E4_post_clustering_sct_no_S521G_pcadims_15_res_0.7.csv",
              "11_subclusters_de_genes_no_S521G/astrocyte/de_genes_astrocyte_subcluster_5_vs_other_astrocyte_no_S521G_pcadims_15_res_0.9.csv",
              "11_subclusters_de_genes_no_S521G/astrocyte/de_genes_astrocyte_subcluster_5_PS19-E4-S_S_vs_PS19-E4_no_S521G_pcadims_15_res_0.9.csv",
              "11_subclusters_de_genes_no_S521G/astrocyte/de_genes_astrocyte_subcluster_7_vs_other_astrocyte_no_S521G_pcadims_15_res_0.9.csv",
              "11_subclusters_de_genes_no_S521G/astrocyte/de_genes_astrocyte_subcluster_7_PS19-E4-S_S_vs_PS19-E4_no_S521G_pcadims_15_res_0.9.csv",
              "11_subclusters_de_genes_no_S521G/microglia/de_genes_microglia_subcluster_8_vs_other_microglia_no_S521G_pcadims_15_res_0.9.csv",
              "11_subclusters_de_genes_no_S521G/microglia/de_genes_microglia_subcluster_8_PS19-E4-S_S_vs_PS19-E4_no_S521G_pcadims_15_res_0.9.csv")


############
#volcano plots
############
for(i in 1:length(de_files)){
  #these files have DE genes that are significant
  signif_res <- read.table(de_files[i], 
                           sep = ",",header = TRUE)

  #label the top 30 DE genes with the highest and lowest Log2FC and p<0.05 (not adj-p).
  signif_res_to_label <- signif_res %>% filter(p_val < 0.05) %>% arrange(desc(abs(avg_log2FC))) %>% head(30)
  write.csv(signif_res_to_label, 
            file = file.path(outdir, paste0("top30genes_",basename(de_files[i]))))
  genes_to_label <- signif_res_to_label$gene
  
  #set the height based on the input
  height_selected <- 12
  
  #set the limits for x-axis
  xlim_selected <- c(-3,3)
  if(max(abs(signif_res$avg_log2FC)) < 2){
    xlim_selected <- c(-2,2)
  }
  
  pdf(paste0(outdir,
             "volcano_plot_",
             tools::file_path_sans_ext(basename(de_files[i])),
             ".pdf"),
      height = height_selected,
      width = 20)
  print(EnhancedVolcano(signif_res,
                        lab = signif_res$gene,
                        title = tools::file_path_sans_ext(basename(de_files[i])),
                        selectLab = genes_to_label,
                        x = 'avg_log2FC',
                        y = 'p_val',
                        FCcutoff = 0.4,
                        pCutoff = 0.05,
                        max.overlaps = Inf,
                        xlim = xlim_selected,
                        drawConnectors = TRUE,
                        legendLabSize = 10,
                        ylab = bquote(~-Log[10] ~ italic(p-value)),
                        legendLabels = c("NS", expression(Log[2] ~ FC), 
                                         "p-value", 
                                         expression(p-value ~ and
                                                    ~ log[2] ~ FC))
  ))
  dev.off()
}

#save the session info
writeLines(capture.output(sessionInfo()), 
           file.path(outdir,"sessionInfo.txt"))

print("********** Script completed! **********")

################## END ################## 
