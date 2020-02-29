###
### Full R-based analysis for CRISPR deletion
###

  ## input dds for chrom20, use deseq2 to test for DGE, fgsea to perform GSEA


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
source("~/projects/glioma/scripts/config.R")


## set functions
load_counts <- function() {
  ## load counts
  setwd(r_dir)
  if (cell_line!="none") {
    filename <- paste(batch,cell_line,drop_prnt,quant_tool,mol_type,"txi.Rdata",sep="_")
  } else {
    filename <- paste(batch,drop_prnt,quant_tool,mol_type,"txi.Rdata",sep="_")
  }
  load_name <- load(filename)
  return(txi)
}


load_annotations <- function() {
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
  pheno$pair_id <- factor(pheno$pair_id)
  pheno$treatment_group <- factor(pheno$treatment_group)
  levels(pheno$treatment_group)[levels(pheno$treatment_group)=="cas"] <- "DEL"
  levels(pheno$treatment_group)[levels(pheno$treatment_group)=="NC"] <- "CTL"
  pheno$treatment_group <- factor(pheno$treatment_group, levels=c("CTL","DEL"))
  return(pheno)
}


make_paired_dds <- function(txi,pheno) {
  ## set up dds data object with model to test effect of treatment within pairs
  cat(paste("making paired dds for",batch,cell_line,drop_prnt,quant_tool,mol_type),"...\n")
  dds <- DESeqDataSetFromTximport(txi, pheno, ~pair_id+treatment_group)
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
  cat(paste("paired dds made for",batch,cell_line,drop_prnt,quant_tool,mol_type),"\n\n")
  return(dds)
}


make_pooled_dds <- function(txi,pheno) {
  ## set up dds data object with model to test effect of treatment group across all days
  cat(paste("making pooled dds for",batch,cell_line,drop_prnt,quant_tool,mol_type),"...\n")
  dds <- DESeqDataSetFromTximport(txi, pheno, ~treatment_group)
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
  cat(paste("pooled dds made for",batch,cell_line,drop_prnt,quant_tool,mol_type),"\n\n")
  return(dds)
}


run_deseq <- function(dds,design) {
  ## run differential expression tests
  setwd(r_dir)
  if (cell_line!="none") {
    filename <- paste(batch,cell_line,drop_prnt,quant_tool,mol_type,design,"dds.Rdata",sep="_")
  } else {
    filename <- paste(batch,drop_prnt,quant_tool,mol_type,design,"dds.Rdata",sep="_")
  }
  if (file.exists(filename)) {
    cat(filename,"already exists; loading prior...\n\n")
    lname <- load(filename)
  } else {
    cat(filename,"not found; running DESeq() now...\n\n")
    dds <- DESeq(dds)  ## estimate size factors and dispersions, fit negative binomial model, test
    save(dds, file=filename)
  }
  return(dds)
}


get_deseq_res <- function(dds,design,fc=1,alpha=0.05) {
  ## test for DE vs fold change threshold, FDR 5%
  if (!"results" %in% mcols(mcols(dds))$type) {  ## test for presence of results, taken from Mike Love's code at https://github.com/mikelove/DESeq2/blob/106c19064bec437091b50496a827fd0b1823375e/R/results.R
    cat("Total DEG for",paste(batch,cell_line,drop_prnt,quant_tool,mol_type,paste0(design,":")),"none\n")
  } else {
    setwd(r_dir)
    if (cell_line!="none") {
      filename <- paste(batch,cell_line,drop_prnt,quant_tool,mol_type,design,"wald_res.Rdata",sep="_")
    } else {
      filename <- paste(batch,drop_prnt,quant_tool,mol_type,design,"wald_res.Rdata",sep="_")
    }
    if (file.exists(filename)) {
      cat(filename,"already exists; loading prior...\n")
      lname <- load(filename)  ## load wald_res
      wald_res_rank <- wald_res[order(wald_res$padj), ]
      wald_res_rank$rank <- seq(1,nrow(wald_res_rank))
      return(list(wald_res,wald_res_rank))
    } else {
      wald_res <- results(dds, lfcThreshold=log2(fc), altHypothesis="greaterAbs", alpha=alpha)
      wald_res$symbol <- mapIds(org.Hs.eg.db, keys=base::substr(row.names(wald_res),start=1,stop=15), column=c("SYMBOL"), keytype=c("ENSEMBL"), multiVals=c("first"))
      wald_res$entrez <- mapIds(org.Hs.eg.db, keys=base::substr(row.names(wald_res),start=1,stop=15), column=c("ENTREZID"), keytype=c("ENSEMBL"), multiVals=c("first"))
      wald_res <- wald_res[ , c(7,8,1:6)]
      wald_res_rank <- wald_res[order(wald_res$padj), ]
      wald_res_rank$rank <- seq(1,nrow(wald_res_rank))
      summary(wald_res)
      cat(paste("Total DEG for",paste(batch,cell_line,drop_prnt,quant_tool,mol_type,paste0(design,":")),sum(wald_res$padj < alpha, na.rm=TRUE),"\n\n"))
      print(wald_res_rank[1:5,c("symbol","padj")])
      save(wald_res, file=filename)
      setwd(table_dir)
      filename <- sub("wald_res.Rdata","wald_res_rank.tsv",filename)
      write.table(wald_res_rank, file=filename, quote=FALSE, sep="\t", row.names=TRUE)
      return(list(wald_res,wald_res_rank))
    }
  }
}


expr_from_dds_norm <- function(ensembl_ids,dds) {
  ## this function generates a data matrix of gene counts for selected genes
  ## based on ensembl ids in dds rownames
  count_data <- data.frame()
  for (id_x in ensembl_ids) {
    if (length(grep(id_x,rownames(counts(dds))))>0) {
      target_row <- grep(id_x,rownames(counts(dds)))
      target_counts <- log2(counts(dds, normalized=TRUE)[target_row,] + 1)
      target_df <- data.frame(
                      gene_id=rep(id_x,ncol(dds)),
                      count=target_counts,
                      tx_group=colData(dds)$treatment_group,
                      gene_name=rep(names(ensembl_ids[ensembl_ids==id_x]),ncol(dds))
                    )
      count_data <- rbind.data.frame(count_data, target_df)
    }
  }
  return(count_data)
}


plot_deseq_res <- function(dds,wald_res) {
  if (cell_line!="none") {
    filename <- paste(batch,cell_line,drop_prnt,quant_tool,mol_type,design,"des_res.pdf",sep="_")
  } else {
    filename <- paste(batch,drop_prnt,quant_tool,mol_type,design,"des_res.pdf",sep="_")
  }
  setwd(plot_dir)
  pdf(file=filename, height=7, width=7)
  DESeq2::plotMA(wald_res, ylim=c(-5.5,5.5), main="TRT vs CTL", alpha=0.05)
  text(x=150000, y=5.5, labels=paste(sum(!is.na(wald_res$padj) & wald_res$padj<0.05 & wald_res$log2FoldChange>0), "up"), pos=2, cex=0.95)
  text(x=150000, y=(-5.5), labels=paste(sum(!is.na(wald_res$padj) & wald_res$padj<0.05 & wald_res$log2FoldChange<0), "down"), pos=2, cex=0.95)
  top_id <- rownames(wald_res[which.min(wald_res$padj),])
  names(top_id) <- wald_res[which.min(wald_res$padj),"symbol"]
  count_data <- expr_from_dds_norm(top_id,dds)
  ggp_ylab <- "Log2(counts)"
  ggp <- ggplot(count_data, aes(x=gene_name, y=count, colour=tx_group)) +
            geom_point(size=5,position=position_jitterdodge()) +
            # geom_boxplot(outlier.size=-1,width=0.5) +
            # geom_dotplot(binaxis="y", stackdir="center", position="dodge") +
            # geom_jitter(width=0.1, size=1) +
            labs(
              title=paste("Expression of",names(top_id)),
              x=element_blank(),
              y=ggp_ylab
            ) +
            theme(
              panel.background=element_rect(fill="white", colour="black", size=0.5, linetype="solid"),
              panel.grid.major=element_line(color="white"),
              panel.grid.minor=element_line(color="white"),
              plot.title=element_text(size=14,face="bold",hjust=0.5),
              axis.title.x=element_text(size=12, face="plain"),
              axis.title.y=element_text(size=12, face="plain"),
              axis.text.x=element_text(size=12, face="plain", angle=45, hjust=1),
              axis.text.y=element_text(size=10, face="plain"),
              legend.key=element_rect(fill="white")
            )
  print(ggp)
  dev.off()
  return(ggp)
}


run_fgsea <- function(wald_res,design) {
  ## run subramanian-style GSEA with hallmark pathways
  if (cell_line!="none") {
    filename <- paste(batch,cell_line,drop_prnt,quant_tool,mol_type,design,"fgsea.Rdata",sep="_")
  } else {
    filename <- paste(batch,drop_prnt,quant_tool,mol_type,design,"fgsea.Rdata",sep="_")
  }
  ## load msigdb gsea pathways
  pathways <- gmtPathways("https://data.broadinstitute.org/gsea-msigdb/msigdb/release/6.2/h.all.v6.2.entrez.gmt")
  ## clean dge data
  res <- wald_res[(!is.na(wald_res[,"entrez"]) & !duplicated(wald_res[,"entrez"])),c("entrez","stat")]
  ## run fgsea
  set.seed(1)
  stats <- sort(setNames(as.numeric(res[,"stat"]), res[,"entrez"]), decreasing=TRUE)
  fgseaRes <- fgsea(pathways=pathways, stats=stats, nperm=10000)
  ## format names
  temp <- fgseaRes  ## this is hack to get around unpredictable list behavior
  temp$pathName <- gsub("_", " ", substr(as.character(temp$pathway), (regexpr("_", as.character(temp$pathway))+1), nchar(temp$pathway)))
  fgseaRes[,"pathName"] <- temp$pathName
  fgseaRes <- fgseaRes[order(fgseaRes$padj,decreasing=FALSE),]
  setwd(r_dir)
  save(fgseaRes, file=filename)
  return(fgseaRes)
}


plot_fgsea <- function(fgseaRes,design) {
  if (cell_line!="none") {
    filename <- paste(batch,cell_line,drop_prnt,quant_tool,mol_type,design,"fgsea.pdf",sep="_")
  } else {
    filename <- paste(batch,drop_prnt,quant_tool,mol_type,design,"fgsea.pdf",sep="_")
  }
  x_lab <- "Gene set"
  y_lab <- "Normalized Enrichment Score"
  title_lab <- "Hallmark sets"
  subtitle_lab <- paste(cell_line,"TRT vs CTL")
  enr_colors <- c("grey","red")
  ggp <- ggplot(fgseaRes[abs(NES)>1.2], aes(reorder(pathName, NES), NES)) +
    geom_col(aes(fill=padj<0.05)) + coord_flip() +
    labs(x=x_lab, y=y_lab, title=title_lab, subtitle=subtitle_lab) + theme_minimal() +
    scale_fill_manual(values=enr_colors) + theme(legend.position="none")
  setwd(plot_dir)
  ggsave(filename=filename, plot=ggp, device="pdf",
    width=7, height=7, units=c("in"))  ## to plot without age_strat, 7x7 is good
  return(ggp)
}


plot_select <- function(dds,wald_res_rank,selection) {
  if (selection=="top") {
    select_ids <- head(rownames(wald_res_rank),6)
    names(select_ids) <- head(wald_res_rank$symbol,6)
    ggp_title <- "Expression Levels of Top Ranked Genes by p-value"
  } else if (selection=="qpcr") {
    select_ids <- c("ENSG00000258366", "ENSG00000026036", "ENSG00000125508", "ENSG00000101216", "ENSG00000197457", "ENSG00000101213", "ENSG00000101246", "ENSG00000125520", "ENSG00000203896")
    names(select_ids) <- c("RTEL1", "RTEL1-TNFRSF6B", "SRMS", "GMEB2", "STMN3", "PTK6", "ARFRP1", "SLC2A4RG", "LIME1")
    ggp_title <- "Expression Levels of qPCR Target Genes"
  } else if (selection=="bio") {
    stop("selection type",selection,"not yet supported")
    ggp_title <- "Expression Levels of Biologically Relevant Genes"
  }
  if (cell_line!="none") {
    filename <- paste(batch,cell_line,drop_prnt,quant_tool,mol_type,design,selection,"slx.pdf",sep="_")
  } else {
    filename <- paste(batch,drop_prnt,quant_tool,mol_type,design,selection,"slx.pdf",sep="_")
  }
  count_data <- expr_from_dds_norm(select_ids,dds)
  dat <- dplyr::group_by(count_data, gene_name) %>% summarise(mean=mean(count, na.rm=TRUE)) %>% arrange(desc(mean))
  ord <- as.character(dat$gene_name) ## put them in descending order by mean
  count_data$gene_name <- factor(count_data$gene_name, levels=ord)
  ggp_ylab <- "Log2 Expression Level"
  ggp_legend <- "Treatment"
  ggp <- ggplot(count_data, aes(x=gene_name, y=count, colour=tx_group)) +
            geom_point(size=2,position=position_jitterdodge()) +
            # geom_boxplot(outlier.size=-1,width=0.75) +
            # geom_boxplot(outlier.size=-1,position=position_dodge2()) +
            # geom_dotplot(binaxis="y", stackdir="center", position="dodge") +
            # geom_jitter(width=0.1, size=1.5) +
            stat_summary(fun.y=median,fun.ymin=median,fun.ymax=median,
              geom="crossbar",size=0.5,position=position_dodge2()) +
            labs(
              title=ggp_title,
              x=element_blank(),
              y=ggp_ylab,
              colour=ggp_legend
            ) +
            # scale_colour_discrete(labels=c("DEL","CTL")) +
            scale_colour_manual(values=c("black","grey")) +   ## ,labels=c("DEL","CTL")
            theme(
              panel.background=element_rect(fill="white", colour="black", size=0.5, linetype="solid"),
              panel.grid.major=element_line(color="white",linetype="solid"),
              panel.grid.minor.y=element_line(color="white",linetype="dotted"),
              panel.grid.minor.x=element_line(color="black",linetype="dotted"),
              plot.title=element_text(size=14,face="bold",hjust=0.5),
              axis.title.x=element_text(size=12, face="plain"),
              axis.title.y=element_text(size=12, face="plain"),
              axis.text.x=element_text(size=12, face="plain", angle=45, hjust=1),
              axis.text.y=element_text(size=10, face="plain"),
              legend.key=element_rect(fill="white")
            )
  ggsave(filename=filename,plot=ggp,device="pdf",width=10,height=6,units=c("in"))
}


examine_select <- function(wald_res_rank,selection) {
  if (selection=="top") {
    select_ids <- head(rownames(wald_res_rank),6)
  } else if (selection=="qpcr") {
    select_ids <- c("ENSG00000258366", "ENSG00000026036", "ENSG00000125508", "ENSG00000101216", "ENSG00000197457", "ENSG00000101213", "ENSG00000101246", "ENSG00000125520", "ENSG00000203896")
  } else if (selection=="bio") {
    stop("selection type",selection,"not yet supported")
  }
  if (cell_line!="none") {
    filename <- paste(batch,cell_line,drop_prnt,quant_tool,mol_type,design,selection,"res.tsv",sep="_")
  } else {
    filename <- paste(batch,drop_prnt,quant_tool,mol_type,design,selection,"res.tsv",sep="_")
  }
  focus <- rownames(wald_res_rank) %in% select_ids
  cols <- c("symbol","stat","pvalue","padj")
  print(wald_res_rank[focus,cols])
  select_res <- wald_res_rank[focus,]
  setwd(table_dir)
  write.table(select_res, file=filename, quote=FALSE, sep="\t", row.names=TRUE)
  return(select_res)
}


## calculations
cat("::: starting DESeq2 :::\n")

txi <- load_counts()
pheno <- load_annotations()


design <- "pair"
dds <- make_paired_dds(txi,pheno)
dds <- run_deseq(dds,design)
payload <- get_deseq_res(dds,design,fc=1,alpha=0.05)
wald_res <- payload[[1]]
wald_res_rank <- payload[[2]]
ggp <- plot_deseq_res(dds,wald_res)
fgseaRes <- run_fgsea(wald_res,design)
ggp <- plot_fgsea(fgseaRes,design)
ggp <- plot_select(dds,wald_res_rank,"top")
ggp <- plot_select(dds,wald_res_rank,"qpcr")
# ggp <- plot_select(dds,wald_res_rank,"bio")
select_res <- examine_select(wald_res_rank,"qpcr")


design <- "pool"
dds <- make_pooled_dds(txi,pheno)
dds <- run_deseq(dds,design)
payload <- get_deseq_res(dds,design,fc=1,alpha=0.05)
wald_res <- payload[[1]]
wald_res_rank <- payload[[2]]
ggp <- plot_deseq_res(dds,wald_res)
fgseaRes <- run_fgsea(wald_res,design)
ggp <- plot_fgsea(fgseaRes,design)
ggp <- plot_select(dds,wald_res_rank,"top")
ggp <- plot_select(dds,wald_res_rank,"qpcr")
# ggp <- plot_select(dds,wald_res_rank,"bio")
select_res <- examine_select(wald_res_rank,"qpcr")



################################################################################
################################################################################
################################################################################



if(FALSE) {

ebl_ids <- get_topgenes(wald_res,fgseaRes)
count_data <- expr_from_dds_norm(dds,ebl_ids=ebl_ids)
plot_topgenes(dds,wald_res,count_data,genes_per_fig=4)


plot_topgenes <- function(dds,wald_res,count_data,genes_per_fig=4) {
  gpf <- genes_per_fig
  color_values <- c(brewer.pal(9,"Greys")[9],brewer.pal(9,"Greys")[5])
  y_lab <- "Log transformed counts"
  setwd("/scratch/chd5n/organoids/results/plots/calcium/")
  filename <- sub("dds.Rdata","des_topgenes.pdf",set_dds_filename(bso))
  pdf(file=filename, height=7, width=7)
  r1 <- 1
  for (r2 in seq(from=ncol(dds)*gpf,to=nrow(count_data),by=ncol(dds)*gpf)) {
    count_subset <- count_data[r1:r2,]
    r1 <- r2 + 1
    ggp <- ggplot(count_subset, aes(x=gene_name, y=count, colour=tx_group)) +
              # geom_boxplot(outlier.size=-1,width=0.5) +
              # geom_jitter(width=0.1, size=1) +
              geom_boxplot(outlier.size=-1,position=position_dodge2()) +
              # geom_jitter(size=1,position=position_dodge2(0.8,preserve="single")) +
              geom_point(size=1,position=position_jitterdodge()) +
              labs(title="Differential gene expression", x=element_blank(), y=y_lab) +
              scale_colour_manual(values=color_values) +
              theme(
                panel.background=element_rect(fill="white", colour="black", size=0.5, linetype="solid"),
                panel.grid.major=element_line(color="white"),
                panel.grid.minor=element_line(color="white"),
                plot.title=element_text(size=14,face="bold",hjust=0.5),
                axis.title.x=element_text(size=12, face="plain"),
                axis.title.y=element_text(size=12, face="plain"),
                axis.text.x=element_text(size=12, face="plain", angle=45, hjust=1),
                axis.text.y=element_text(size=10, face="plain"),
                legend.key=element_rect(fill="white")
              )
      print(ggp)
  }
  dev.off()
}


filename <- sub("dds","fgsea",set_dds_filename(bso))
fgseaRes <- run_fgsea(wald_res,filename)
filename <- sub("dds.Rdata","fgsea.pdf",set_dds_filename(bso))
ggp <- plot_fgsea(fgseaRes,filename)
ebl_ids <- get_topgenes(wald_res,fgseaRes)
count_data <- expr_from_dds_norm(dds,ebl_ids=ebl_ids)
plot_topgenes(dds,wald_res,count_data,genes_per_fig=4)


make_boxplot <- function(count_data,data_type,bso,expr_type,y_lab="Transformed counts") {
  ## for all organoids
  target_col <- colData(dds)$tx_group == "CTL"
  dds <- dds[,target_col]
  ## for barcelona
  target_col <- rep(TRUE,ncol(dds))

  y_lab <- "Log transformed counts"
  setwd("/scratch/chd5n/organoids/results/plots/calcium/")
  ggp <- ggplot(count_data, aes(x=tx_group, y=count)) +

  ggp <- ggplot(count_data, aes(x=gene_name, y=count)) +
            geom_boxplot(outlier.size=-1,width=0.5) +
            geom_jitter(width=0.1, size=1) +
            labs(title="Expression of PRKCA", x=element_blank(), y=y_lab) +
            theme(
              panel.background=element_rect(fill="white", colour="black", size=0.5, linetype="solid"),
              panel.grid.major=element_line(color="white"),
              panel.grid.minor=element_line(color="white"),
              plot.title=element_text(size=14,face="bold",hjust=0.5),
              axis.title.x=element_text(size=10, face="plain"),
              axis.title.y=element_text(size=10, face="plain"),
              axis.text.x=element_text(size=10, face="plain", angle=45, hjust=1),
              axis.text.y=element_text(size=8, face="plain"),
              legend.key=element_rect(fill="white")
            )

  ## barcelona
  filename <- paste("bc","ge","des","norm","sx.pdf",sep="_")
  ## organoids
  filename <- paste("all_ctl",tolower(expr_type),"des","norm","sx.pdf",sep="_")
  ## calcium results

  ggsave(filename=filename, plot=ggp, device="pdf", height=6, width=5, units="in")

}


get_topgenes <- function(wald_res,fgseaRes) {
  ## set up gene-pathway relationships for results
  library(tidyverse)
  pathways <- gmtPathways("https://data.broadinstitute.org/gsea-msigdb/msigdb/release/6.2/h.all.v6.2.entrez.gmt")
  wald_df <- base::as.data.frame(wald_res)
  wald_df$ensembl <- rownames(base::as.data.frame(wald_res))
  pway_nested <- tibble::enframe(pathways,name="pathway",value="entrez")
  pway_unnest <- tidyr::unnest(pway_nested, cols="entrez")
  res_enrich <- dplyr::inner_join(pway_unnest, wald_df, by="entrez")
  ## set up biologically meaningful genes
  keep <- abs(fgseaRes$NES)>2.5 & fgseaRes$padj<0.05  ## pathways with highest enrichment
  key_pways <- data.frame(fgseaRes[,c("pathway","ES","NES","padj","size")])[keep,]
  key_pways <- key_pways[order(key_pways[,"NES"], decreasing=TRUE),"pathway"]
  bio_genes <- vector()
  for (pway in key_pways) {
    bio_genes <- append(bio_genes, res_enrich[res_enrich$pathway==pway,"entrez"])
  }
  bio_genes <- unlist(bio_genes, use.names=FALSE)
  bio_genes <- unique(bio_genes)
  bio_dge_ebl <- vector()
  for (enz_x in bio_genes) {
    if (all(!is.na(res_enrich[res_enrich$entrez==enz_x,"padj"]) & res_enrich[res_enrich$entrez==enz_x,"padj"]<0.05)) {
      bio_dge_ebl <- append(bio_dge_ebl, unique(res_enrich[res_enrich$entrez==enz_x,"ensembl"]))
    }
  }
  bio_dge_ebl <- unlist(bio_dge_ebl, use.names=FALSE)
  ## identify topgenes
  keep <- res_enrich$ensembl %in% bio_dge_ebl
  topgene_df <- res_enrich[keep,]
  topgene_df <- topgene_df[order(topgene_df$stat,decreasing=TRUE),]
  topgene_df <- topgene_df %>% nest(pathways=c(pathway))
  ebl_ids <- setNames(topgene_df$ensembl, topgene_df$symbol)
  return(ebl_ids)
}
}
