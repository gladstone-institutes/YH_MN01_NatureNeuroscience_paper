# Dockerfile for the Seurat 4.1.1
FROM satijalab/seurat:4.1.0

#install specific versions of R packages
RUN R -e 'install.packages("remotes")'
RUN R -e "library('remotes');install_version('ggplot2', '3.3.6')"
RUN R -e "library('remotes');install_version('dplyr', '1.0.9')"
RUN R -e "library('remotes');install_version('gdata', '2.18.0.1')"
RUN R -e 'library("remotes");install_version("harmony", "0.1.0")'
RUN R -e "library('remotes');install_version('gplots', '3.1.3')"
RUN R -e "library('remotes');install_version('RColorBrewer', '1.1-3')"
RUN R -e 'library("remotes");install_version("pheatmap", "1.0.12")'
RUN R -e 'BiocManager::install(version = "3.14",ask = FALSE)'
RUN R -e 'BiocManager::install("glmGamPoi", version = "3.14")'
