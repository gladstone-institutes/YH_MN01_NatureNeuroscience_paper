#KEGG Enrichment Analysis of the enriched genes
#this script uses clusterProfiler::enrichKEGG() function

#run locally on Ayushi's laptop

#load the required packages
library(org.Mm.eg.db)
library(clusterProfiler)
library(tidyr)

#set the working directories
indir <- paste0("~/Dropbox (Gladstone)/YH_MN02/",
                "seurat_analysis_no_S521G/11_subclusters_de_genes_no_S521G/")
setwd(indir)
outdir_astro <- "../12_02_subclusters_KEGG_enriched_pathways_no_S521G_significant_DE_genes_only/astrocyte/"
outdir_micro <- "../12_02_subclusters_KEGG_enriched_pathways_no_S521G_significant_DE_genes_only/microglia/"
dir.create(outdir_astro, showWarnings = FALSE, recursive = TRUE)
dir.create(outdir_micro, showWarnings = FALSE, recursive = TRUE)


de_files <- c(list.files(path = "./astrocyte", pattern = "de_genes"),
              list.files(path = "./microglia", pattern = "de_genes"))
background_genes_astrocytes <- read.table(
  "astrocyte/nonzero_bakcground_gene_list_astrocyte_data_no_S521G_pcadims_15_res_0.9.csv",
  sep = ",",
  header = TRUE)
background_genes_microglia <- read.table(
  "microglia/nonzero_bakcground_gene_list_microglia_data_no_S521G_pcadims_15_res_0.9.csv",
  sep = ",",
  header = TRUE)

for(i in 1:length(de_files)){
  
  if(grepl("astrocyte", de_files[i], fixed = TRUE)){
    #read in the DE gene list
    #filter the data for DE genes that are significant
    signif_res <- read.table(paste0("./astrocyte/", de_files[i]), 
                             sep = ",",
                             header = TRUE) %>% 
      dplyr::filter(p_val_adj <0.05) %>%
      as.data.frame()
    if(nrow(signif_res) < 100){
      print(paste0(nrow(signif_res) ,
                   " significant DE genes for: ", de_files[i]))
      next
    }
    #select background set of genes
    background_genes <- background_genes_astrocytes
    #set the output directory
    outdir <- outdir_astro
  } else if(grepl("microglia", de_files[i], fixed = TRUE)){
    #read in the DE gene list
    #filter the data for DE genes that are significant
    signif_res <- read.table(paste0("./microglia/", de_files[i]), 
                             sep = ",",
                             header = TRUE)%>% 
      dplyr::filter(p_val_adj <0.05) %>%
      as.data.frame()
    if(nrow(signif_res) < 100){
      print(paste0(nrow(signif_res) ,
                   " significant DE genes for: ", de_files[i]))
      next
    }
    #select background set of genes
    background_genes <- background_genes_microglia
    #set the output directory
    outdir <- outdir_micro
  }else{
    print("ERROR: Input file is NOT an astrocyte or microglia DE list!")
  }
  
  #rename the first column of the DE gene list
  colnames(signif_res)[1] <- "Geneid"
  
  #get the ENTREZ IDs
  entrezIDs <- AnnotationDbi::select(org.Mm.eg.db, 
                                     keys = background_genes$background_genes%>% 
                                       as.character(),
                                     columns = c("ENTREZID", "SYMBOL"),
                                     keytype = "SYMBOL")
  entrezIDs_signif <- AnnotationDbi::select(org.Mm.eg.db, 
                                            keys = signif_res$Geneid %>% 
                                              as.character(),
                                            columns = c("ENTREZID", "SYMBOL"),
                                            keytype = "SYMBOL")
  
  #-----------------------------------
  #KEGG analysis
  #------------------------------
  ekeggbp <- enrichKEGG(
    gene     = entrezIDs_signif$ENTREZID %>% subset(., !is.na(.)),
    universe = entrezIDs$ENTREZID %>% subset(., !is.na(.)),
    organism    = "mmu",
    minGSSize = 10,
    pvalueCutoff = 0.8,
    keyType = "ncbi-geneid"
  )
  
  #translating gene IDs to human readable symbols
  ekeggbp <- setReadable(ekeggbp, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
  
  #Visualize
  ## save images
  pdf(paste0(outdir,
             tools::file_path_sans_ext(de_files[i]),
             "_ekegg-dot.pdf"),
      height = 8)
  print(dotplot(ekeggbp, showCategory = 20, orderBy="GeneRatio") )
  dev.off()
  
  x2 <- enrichplot::pairwise_termsim(ekeggbp)
  pdf(paste0(outdir,
             tools::file_path_sans_ext(de_files[i]),
             "_ekegg-emap.pdf"),
      height = 8)
  print(emapplot(x2, showCategory = 40) )
  dev.off()
  
  #save the list of enriched pathways
  write.csv(ekeggbp,file = paste0(outdir,
                                  tools::file_path_sans_ext(de_files[i]),
                                  "_ekegg.csv"))
  
}

#save the session info
writeLines(capture.output(sessionInfo()), 
           paste0("../12_02_subclusters_KEGG_enriched_pathways_no_S521G_significant_DE_genes_only/",
                  "sessionInfo.txt"))

print("********** Script completed! **********")

################## END ################## 
