###
### expression analysis count prep
###


## README


  ## samples section reads pheno file from tsv into R
  ## in theory should rely on output from metaFilters.R

  ## files section sets paths to Salmon/RSEM quant files; paths generated
  ## from information in pheno list and hard-coded text

  ## txi section imports quant files (TPM estimates) into lists of count
  ## (gene and transcript) matrices, saves txi data objects as Rdata objects for
  ## downstream GE/TE analyses

  ## takes 4m, 4GB


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


## define functions

set_file_names <- function(batch,cell_line,drop_prnt,quant_tool) {
  ## load pheno
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
  ## make filesets storing path to quant file in concat column
  if (quant_tool == "sal") {
    fileset <- data.frame(
      in_path_top=rep(sal_qnt_dir,nrow(pheno)),
      in_path_mid=paste(batch,pheno$fq_id,sep="/"),
      filename=rep("quant.sf",nrow(pheno))
    )
    fileset <- within(fileset, concat <- paste(paste0(in_path_top,in_path_mid),filename,sep="/"))
  } else if (quant_tool == "rsem") {
    fileset <- data.frame(
      in_path_top=rep(rsem_qnt_dir,nrow(pheno)),
      in_path_mid=paste(pheno$chrom,pheno$fq_id,sep="/"),
      id=pheno$fq_id,
      filename=rep(".rsem.genes.results",nrow(pheno))  ## change later in script if want tx
    )
    fileset <- within(fileset, concat <- paste(paste0(in_path_top,in_path_mid),paste0(id,filename),sep="/"))
  } else {
    stop("Please indicate quantification tool. Supported codes include `sal` and `rsem`.")
  }
  if (all(file.exists(fileset$concat))) {
    cat("all files exist\n")
  } else {
    stop("not all files exist\n")
  }
  return(fileset)
}


get_tx2gene <- function() {
  ## create tx2gene (2 columns: txID:geneID)
  txdb <- transcripts(EnsDb.Hsapiens.v86, columns=c("tx_id", "tx_name", "gene_id"))
  txdf <- as.data.frame(txdb)
  tx2gene <- as.data.frame( cbind(txdf$tx_id, txdf$gene_id) )
  colnames(tx2gene) <- c("txName", "geneID")
  return(tx2gene)
}


import_data <- function(file_paths,quant_tool,mol_type) {
  ## import the transcript abundances as original count estimates with offset
  ## read length default is 75 bp, although unclear how it affects output
  ## basic import with length matrix as offset for downstream analysis
  ## alternatives include library scaledTPM, length and library lengthScaledTPM
  ## uses transcript counts and gene-level summarized counts
  if (quant_tool == "sal") {
    if (mol_type == "ge") {
      tx2gene <- get_tx2gene()
      txi <- tximport(file_paths, type=quant_tool, tx2gene=tx2gene, ignoreTxVersion=TRUE)
    } else if (mol_type == "tx") {
      txi <- tximport(file_paths, type=quant_tool, txOut=TRUE, ignoreTxVersion=FALSE)
    }
  } else if (quant_tool == "rsem") {
    if (mol_type == "ge") {
      txi <- tximport(file_paths, type=quant_tool, txIn=FALSE, txOut=FALSE, ignoreTxVersion=TRUE,
        geneIdCol="gene_id", txIdCol="transcript_id(s)", abundanceCol="TPM", countsCol="expected_count",
        lengthCol="effective_length")
      txi$length[txi$length < 1] <- 1
    } else if (mol_type == "tx") {
      file_paths <- gsub(".genes.", ".isoforms.", file_paths)
      txi <- tximport(file_paths, type=quant_tool, txIn=TRUE, txOut=TRUE, ignoreTxVersion=FALSE,
        txIdCol="transcript_id", abundanceCol="TPM", countsCol="expected_count",
        lengthCol="effective_length")
    }
  }
  return(txi)
}


write_txi <- function(txi, batch, cell_line, drop_prnt, quant_tool, mol_type) {
  ## write txi to file for downstream analyses
  setwd(r_dir)
  if (cell_line!="none") {
    filename <- paste(batch,cell_line,drop_prnt,quant_tool,mol_type,"txi.Rdata",sep="_")
  } else {
    filename <- paste(batch,drop_prnt,quant_tool,mol_type,"txi.Rdata",sep="_")
  }
  save(txi, file=filename)
}


## main script
fileset <- set_file_names(batch,cell_line,drop_prnt,quant_tool)
if (run_type == "live") {
  txi <- import_data(fileset$concat,quant_tool,mol_type)
  write_txi(txi,batch,cell_line,drop_prnt,quant_tool,mol_type)
} else {
  cat("test complete\n")
}
