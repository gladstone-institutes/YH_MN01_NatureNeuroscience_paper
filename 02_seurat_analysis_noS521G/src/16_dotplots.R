#!/usr/bin/env Rscript

#load required packages
library(Seurat)
library(dplyr)
library(ggplot2)

#set all directory paths
basedir <- paste0("/gladstone/bioinformatics/projects/",
                  "mn-1318-maxine-nelson-yadong-huang-snrnaseq-mm10-aug-2022/",
                  "results/data/")
outdir <- paste0("/gladstone/bioinformatics/projects/",
                 "mn-1318-maxine-nelson-yadong-huang-snrnaseq-mm10-aug-2022/",
                 "results/plot/16_dotplots/")
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
##1.A dot plot of the following genes for the four 
#   oligodendrocyte clusters 2, 3, 9, 37
##############################################################
genes_to_plot <- "Mbp, Plp1, Mag, St18, Pde4b, Rnf220, Prr5l, Neat1, Plcl1, Lsamp, Snhg11, Kalrn, Lrrc7, Meg3, Gria1, Grin2a, Grin2b, Epha6, Frmpd4, Cntnap2, Dab1, Celf2, Lrrtm4, Gria2, Serpina3n, H2-D1, H2-K1, B2M, C4b, Klk6, Cd63, Il33, Sgk1, Cd9, Fxyd1, Plvap, Fxyd7, Rnase4, hapoE-transgene, Gstm1, Snca, Calm1, Calm2, Ubb, Hsp90aa1, Hsp90ab1"
genes_to_plot <- unique(as.character(unlist(strsplit(as.character(genes_to_plot), ", "))))
(DotPlot(dat, 
         features = genes_to_plot,
         cols = c("blue", "red")) + 
    RotatedAxis())%>%
  ggsave(file = "dotplot_with_selected_genes_oligo_clusters_2_3_9_37.pdf",
         plot = .,
         width = 20,
         height = 9)


##############################################################
##2.A dot plot of the select genes for all 12 astrocyte subclusters
##############################################################
genes_to_plot <- "Luzp2, Slc7a10, Mfge8, Gfap, Id3, Aqp4, Myoc, Id1, Fabp7, Ctsb, Vim, Osmr, Serpina3n, Gsn, Ggta1, Trpm3, Csmd1, C4b, Cd9, Sparcl1, Plce1, Sgcd, Fos, S100a6, hapoE-transgene, Clu, Nkain2, Pde4b, Slc24a2, St18, Tmeff2, Frmd5, Edil3, Neat1, Plp1, Pex5l, Ptprd, Plcl1, Magi2, Zfp536, Rnf220, Pcdh9, Ano4, Mast4, St6galnac3, Aspa, Dgki, Ank3, Csmd3, Lrp1b, Calm1, Calm2, Ubb, Hsp90aa1, Hsp90ab1"
genes_to_plot <- unique(as.character(unlist(strsplit(as.character(genes_to_plot), ", "))))
(DotPlot(dat_astrocytes, 
         features = genes_to_plot, 
         cols = c("blue", "red")) + 
    RotatedAxis())%>%
  ggsave(file = "dotplot_with_selected_genes_all_astrocyte_subclusters.pdf",
         plot = .,
         width = 20,
         height = 9)


##############################################################
##3.A dot plot of the select genes for all 15 microglia subclusters
##############################################################
genes_to_plot <- "Hexb, Cst3, Cx3cr1, Ctsd, Csf1r, Ctss, Sparc, Tmsb4x, P2ry12, C1qa, C1qb, Tmem119, Tyrobp, Ctsb, hapoE-transgene, B2m, Fth1, Lyz2, Trem2, Axl, Cst7, Ctsl, Lpl, Cd9, Csf1, Ccl6, Itgax, Clec7a, Lilrb4a, Timp2, Gpnmb, Atp6v0d2, Apobec1, Igf1, Lyst, Arhgap24, Mamdc2, Fnip2, Colec12, Apbb2, Myo1e, Pid1, Myo5a, Spp1, Ddhd1, Osbpl8, Fmn1, Neat1, Prkg1, Rftn1, Fgf13, Fat3, Sash1, Ctsb, Tanc2, Nav3, Csmd3, Calm1, Calm2, Ubb, Hsp90aa1, Hsp90ab1"
genes_to_plot <- unique(as.character(unlist(strsplit(as.character(genes_to_plot), ", "))))
(DotPlot(dat_microglia, 
         features = genes_to_plot, 
         cols = c("blue", "red")) + 
    RotatedAxis())%>%
  ggsave(file = "dotplot_with_selected_genes_all_microglia_subclusters.pdf",
         plot = .,
         width = 20,
         height = 9)


#save the session info
writeLines(capture.output(sessionInfo()), 
           file.path(outdir,"sessionInfo.txt"))

print("********** Script completed! **********")

################## END ################## 
