###
### Expression Analysis, Configuration Commands
###


# load libraries


## core utilities

# install.packages("BiocManager")
# install.packages("devtools")
# install.packages("https://cran.r-project.org/src/contrib/Archive/curl/curl_4.0.tar.gz",repo=NULL,type="source")
BiocManager::install()

## plotting

# install.packages("ggplot2")
# install.packages("ggfortify")
# install.packages("pheatmap")
# install.packages("RColorBrewer")
# install.packages("Rtsne")
# install.packages("tsne")
# install.packages("pROC")
# install.packages("VennDiagram")
# install.packages("cowplot")
# install.packages("gridGraphics")
library("ggplot2")
# library("ggfortify")
library("pheatmap")
library("RColorBrewer")
library("Rtsne")
library("tsne")
# library("pROC")
library("VennDiagram")
library("cowplot")
library("gridGraphics")


## data parsing

# BiocManager::install("statmod")
# BiocManager::install("apeglm")
# BiocManager::install("BiocParallel")
# install.packages("tidyverse")
# install.packages("dplyr")
# install.packages("readr")
library("statmod")
library("apeglm")
library("BiocParallel")
# library("tidyverse")
library("readr")
library("dplyr")

## annotation

# BiocManager::install("AnnotationDbi")
# BiocManager::install("org.Hs.eg.db")
# BiocManager::install("org.Mm.eg.db")
# BiocManager::install("ensembldb")
# BiocManager::install("EnsDb.Hsapiens.v86")
# BiocManager::install("EnsDb.Mmusculus.v79")
library("AnnotationDbi")
library("org.Hs.eg.db")
library("ensembldb")
library("EnsDb.Hsapiens.v86")

## gene/txp expression

# BiocManager::install("EDASeq")
# BiocManager::install("limma")
# BiocManager::install("edgeR")
# BiocManager::install("DESeq2")
# BiocManager::install("sva")
# BiocManager::install("RUVSeq")
# BiocManager::install("WGCNA")
# BiocManager::install("tximport")
library("EDASeq")
library("limma")
library("edgeR")
library("DESeq2")
library("sva")
library("RUVSeq")
# library("WGCNA"); allowWGCNAThreads()
library("tximport")

## gene set enrichment

# BiocManager::install("goseq")
# BiocManager::install("GO.db")
# BiocManager::install("fgsea")
# BiocManager::install("reactome.db")
# BiocManager::install("GSEABase")
# BiocManager::install("GSVA")
# devtools::install_github('dviraran/xCell')
library("goseq")
library("GO.db")
library("fgsea")
library("reactome.db")
library("GSEABase")            # must be loaded for xCell
library("GSVA")                # must be loaded for xCell
# library("xCell")

## survival

# install.packages("survival")
# library("survival")

## data

# BiocManager::install("TCGAbiolinks") #for 3.5.1 required archived version of cmprsk `install_version("cmprsk", version="2.2-7", repos="http://cran.us.r-project.org")`
# library('TCGAbiolinks')
# BiocManager::install("curatedCRCData")
# library("curatedCRCData")


# define data objects

set_group_colors <- function(pheno,goi) {
  ## function to define group colors
  colors <- list()
  for (i in 1:length(goi)) {
    if (class(unlist(pheno[, goi[i]])) == "character") {
      gn <- length(levels(as.factor(unlist(pheno[, goi[i]]))))
      ## get hcl color wheel coordinates for various cohorts (ggplot starts at 15 and increases to 360 by 360/gn)
      inc <- 360/gn
      cords <- vector()
      for (j in seq(15, 360, inc)) {
        if (j==15) {
          cords <- j
        } else {
          cords <- append(cords, j)
        }
      }
      ## get hcl hexidecimal color codes
      hex <- hcl(h=cords, c=100, l=65)
      names(hex) <- levels(as.factor(unlist(pheno[,goi[i]])))
      colors[[i]] <- hex
    } else {
      colors[[i]] <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
    }
  }
  names(colors) <- goi
  return(colors)
}


## datasets
quant_tools <- c("sal","rsem")
#filtered_sets <- c("UqSamp")
#analysis_groups <- c("NatCrcPar", "IndSmpAll", "IndSmpTrn", "IndSmpTst")
mol_types <- c("ge","tx")

## directories
work_dir <- "/scratch/chd5n/glioma/"
ann_dir <- paste0(work_dir,"annotations/")
sal_qnt_dir <- paste0(work_dir,"quant/")
rsem_qnt_dir <- paste0(work_dir,"rsem/")
r_dir <- paste0(work_dir,"rdata/")
res_dir <- paste0(work_dir,"results/")
plot_dir <- paste0(res_dir,"plots/")
table_dir <- paste0(res_dir,"tables/")
scripts_dir <- "/home/chd5n/projects/glioma/scripts/"

## colors
setwd(ann_dir)
filename <- "ann.tsv"
pheno <- read_tsv(filename)
goi <- c("cell_line","treatment_group","rin")  ## groups of interest
colors <- set_group_colors(pheno,goi)  ## get colors for each group

# conclusion
cat("glioma R configuration loaded\nREADY TO GO!\n\n")
