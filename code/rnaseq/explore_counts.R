###
### Expression Analysis, Preliminary Exploration
###

  ## load pheno data and txi counts and explore global variation

  ## provides exploratory gene/txp expression analysis
  ## relies on: DESeq2
  ## output: PCA, heatmap

## load config
args <- commandArgs(trailingOnly=TRUE)
batch <- args[2]
cell_line <- args[3]
drop_prnt <- args[4]
quant_tool <- args[5]
mol_type <- args[6]
run_type <- args[7]
cat("batch",batch,"\n");cat("cell_line",cell_line,"\n");cat("drop_prnt",drop_prnt,"\n")
cat("quant_tool",quant_tool,"\n");cat("mol_type",mol_type,"\n");cat("run_type",run_type,"\n\n")
source("~/projects/glioma/code/rnaseq/config.R")

## set functions
make_null_dds <- function(txi,pheno) {
  ## set up dds data object with null model
  dds <- DESeqDataSetFromTximport(txi, pheno, ~1)
  ## pre-filter lowly expressed genes to speed up downstream analysis
      ## want 5-10 counts in a library to prove gene is expressed in that library
      ## want minimum expression in at least the number of libraries in the smallest group of interest or 1/3 of samples
      ## minimum count threshold for a gene across all samples should be 5-10*N
      ## strict filtering performed internally later; this is simply for size reduction and speed
  count_lim <- 5
  sample_lim <- (1/3) * ncol(dds)
  keep <- rowSums( counts(dds) > count_lim ) >= sample_lim    ## sum logical (T/F) output for each row (i.e. gene) to get total number samples with that many counts for that gene
  dds <- dds[keep,]
  ## calculate normalization factors for tx length, library size
  dds <- estimateSizeFactors(dds)
  return(dds)
}

get_naive_vsd <- function(dds) {
  ## preliminary transformations for visualization (e.g. heatmaps, pca plots)
      ## blind=TRUE if exploring data for quality assurance
      ## blind=FALSE when analyzing design
  setwd(r_dir)
  if (cell_line!="none") {
    filename <- paste(batch,cell_line,drop_prnt,quant_tool,mol_type,"vsd.Rdata",sep="_")
  } else {
    filename <- paste(batch,drop_prnt,quant_tool,mol_type,"vsd.Rdata",sep="_")
  }
  if (file.exists(filename)) {
    load_name <- load(filename)
  } else {
    vsd <- vst(dds, blind=TRUE)  #*********************rate limiting
    save(vsd, file=filename)
  }
  return(vsd)
}

calc_sample_dist <- function(vsd) {
  ## calculate sample distances
  setwd(r_dir)
  if (cell_line!="none") {
    filename <- paste(batch,cell_line,drop_prnt,quant_tool,mol_type,"sdm.Rdata",sep="_")
  } else {
    filename <- paste(batch,drop_prnt,quant_tool,mol_type,"sdm.Rdata",sep="_")
  }
  if (file.exists(filename)) {
    load_name <- load(filename)
  } else {
    sdm <- list()
    sdm[["s_dist"]] <- dist(t(assay(vsd)))
    sdm[["sdm"]] <- as.matrix(sdm[["s_dist"]])
    rownames(sdm[["sdm"]]) <- colData(vsd)[,"sample_id"]
    colnames(sdm[["sdm"]]) <- NULL
    save(sdm, file=filename)
  }
  return(sdm)
}

calc_pcs <- function(vsd) {
  ## pca plot calculations
  pcs <- list()
  single_factors <- goi  ## define groups of interest
  pcs[["intgroup"]] <- list()
  for ( i in 1:length(single_factors) ) {
    pcs[["intgroup"]][[i]] <- single_factors[i]
  }
  names(pcs[["intgroup"]]) <- c("Cell Line","Group","RIN")
  pcs[["ggp_title_names"]] <- c("all samples")  ## set dataset names for title
  ## calculate pcs
  pcs[["pca_data"]] <- list()
  for ( i in 1:length(pcs[["intgroup"]]) ) {
    pcs[["pca_data"]][[i]] <- DESeq2::plotPCA(vsd, intgroup=pcs[["intgroup"]][[i]], returnData=TRUE)
    pcs[["pca_data"]][[i]][,"rin"] <- as.factor(colData(vsd)[,"rin"])
    # pcs[["pca_data"]][[i]][,"num_map"] <- as.factor(round(colData(vsd)[,"num_map"]/10^6,1))
    pcs[["pca_data"]][[i]][,"sample_id"] <- as.factor(colData(vsd)[,"sample_id"])
  }
  pcs[["perc_var"]] <- round( 100 * attr(pcs[["pca_data"]][[1]], "percentVar") )
  return(pcs)
}

pca_single <- function(pca_data, perc_var, base_name, aes_name, ig_name) {
  color <- sym(aes_name)
  sz=1.5
  ggp_title <- paste("PCA plot for",base_name,"by",ig_name)
  color_values <- colors[[aes_name]]
  ggp <- ggplot(pca_data, aes(PC1, PC2, color=!!color, label=name)) +
    geom_point(size=sz) +
    xlab(paste0("PC1: ",perc_var[1],"% variance")) +
    ylab(paste0("PC2: ",perc_var[2],"% variance")) +
    coord_fixed() + theme_classic() +
    theme(plot.title=element_text(hjust=0.5)) +
    ggtitle(ggp_title) +
    scale_color_manual(values=color_values)
  return(ggp)
}

pca_single_label <- function(pca_data, perc_var, base_name, aes_name, ig_name, plab) {
  ## single factor, with label
  color <- sym(aes_name)
  sz=1
  ggp_title <- paste("PCA plot for",base_name,"by",ig_name,"with",plab)
  color_values <- colors[[aes_name]]
  ggp <- ggplot(pca_data, aes(PC1, PC2, color=!!color)) +
    geom_point(size=sz) +
    xlab(paste0("PC1: ",perc_var[1],"% variance")) +
    ylab(paste0("PC2: ",perc_var[2],"% variance")) +
    coord_fixed() + theme_classic() +
    theme(plot.title=element_text(hjust=0.5)) +
    ggtitle(ggp_title) +
    scale_color_manual(values=color_values) +
    geom_text(aes(label=!!sym(plab)), hjust=0, vjust=0, size=3)
  return(ggp)
}


## load counts
setwd(r_dir)
if (cell_line!="none") {
  filename <- paste(batch,cell_line,drop_prnt,quant_tool,mol_type,"txi.Rdata",sep="_")
} else {
  filename <- paste(batch,drop_prnt,quant_tool,mol_type,"txi.Rdata",sep="_")
}
load_name <- load(filename)

## load annotations
setwd(ann_dir)
filename <- "ann.tsv"
pheno <- read_tsv(filename)
pheno <- pheno[pheno$chrom==batch,]
if (cell_line!="none") {
  pheno <- pheno[pheno$cell_line==cell_line,]
}
if (drop_prnt=="prdp") {
  pheno <- pheno[pheno$treatment_group!="parental",]
}


## re-order levels for select categories
# cat <- c("location","uva_loc")
# custom_levels <- c("right","transverse","left")
# for (each in cat) {
#   pheno[,each] <- factor(as.factor(unlist(pheno[,each])), levels=custom_levels)
# }

## calculations
cat("::: starting DESeq2 :::\n")
dds <- make_null_dds(txi,pheno)
vsd <- get_naive_vsd(dds)
sdm <- calc_sample_dist(vsd)
pcs <- calc_pcs(vsd)

## main routine
if (run_type == "live") {
  ## figures
  setwd(plot_dir)
  if (cell_line!="none") {
    filename <- paste(batch,cell_line,drop_prnt,quant_tool,mol_type,"des","prm.pdf",sep="_")
  } else {
    filename <- paste(batch,drop_prnt,quant_tool,mol_type,"des","prm.pdf",sep="_")
  }
  cat("plotting to file",paste(getwd(),filename,sep="/"),"\n")
  ## open device
  pdf(file=filename)
  colrs <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  ## map between sample distances, cluster by distances
  pheatmap(sdm[["sdm"]], clustering_distance_rows=sdm[["s_dist"]],
    clustering_distance_cols=sdm[["s_dist"]], col=colrs,
    main=c("Heatmap of sample to sample distances, cluster by dist")
  )
  ## single frame pcas
  for ( i in 1:length(pcs[["intgroup"]]) ) {
    ggp <- list()
    ggp[[1]] <- pca_single(pcs[["pca_data"]][[i]], pcs[["perc_var"]],
                  pcs[["ggp_title_names"]], pcs[["intgroup"]][[i]],
                  names(pcs[["intgroup"]])[i])
    ggp[[2]] <- pca_single_label(pcs[["pca_data"]][[i]], pcs[["perc_var"]],
                  pcs[["ggp_title_names"]], pcs[["intgroup"]][[i]],
                  names(pcs[["intgroup"]])[i], "rin")
    ggp[[3]] <- pca_single_label(pcs[["pca_data"]][[i]], pcs[["perc_var"]],
                  pcs[["ggp_title_names"]], pcs[["intgroup"]][[i]],
                  names(pcs[["intgroup"]])[i], "rin")
    ggp[[4]] <- pca_single_label(pcs[["pca_data"]][[i]], pcs[["perc_var"]],
                  pcs[["ggp_title_names"]], pcs[["intgroup"]][[i]],
                  names(pcs[["intgroup"]])[i], "sample_id")
    for (i in 1:length(ggp)) {
      print(ggp[[i]])
    }
  }
  ## close device to save plots
  dev.off()
} else {
  cat("test complete\n")
}
cat("DESeq2 end\n\n")
