# Install the R packages used in this repository (R 4.3.3, Bioconductor 3.18)

cran <- c("vegan", "fossil", "ggplot2", "ggpubr", "rstatix", "dplyr", "tidyr",
          "openxlsx", "zCompositions", "RVAideMemoire", "glmnet", "pROC", "plotROC",
          "verification", "VIM", "Hmisc", "pheatmap", "ggvenn", "gtsummary",
          "MicrobiomeStat", "BiocManager")
bioc <- c("phyloseq", "ANCOMBC")

install.packages(setdiff(cran, rownames(installed.packages())))
BiocManager::install(setdiff(bioc, rownames(installed.packages())), update = FALSE)
