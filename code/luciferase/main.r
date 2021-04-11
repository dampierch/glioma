## usage: Rscript main.r --args <data_file> <out_prefix> <trans>
##        Rscript main.r --args positives.tsv pos log2
##        Rscript main.r

library(readr)
library(dplyr)
library(geepack)
library(ggplot2)
library(cowplot)
library(scales)

source("input.r")
source("gee.r")
source("ggp_themes.r")
source("plot.r")

main <- function(args) {

    ## set project variables
    pvars <- list()
    pvars$proj_dir <- "../../"
    pvars$data_dir <- paste0(pvars$proj_dir, "data/")
    pvars$data_file <- ifelse(is.na(args[2]), "negatives.tsv", args[2])
    pvars$out_prefix <- ifelse(is.na(args[3]), "neg", args[3])
    pvars$data_name <- paste0(pvars$data_dir, pvars$data_file)
    pvars$trans <- ifelse(is.na(args[4]), NA, args[4])

    ## parse data
    dfs <- get_data(pvars$data_name)
    fits <- lapply(dfs, fit_gee)
    sums <- lapply(fits, parse_gee)
    dfs2 <- lapply(dfs, summarise_data)

    ## make plots
    if (pvars$out_prefix == "neg") {
        ggp <- lapply(dfs, ggp_box, 0.1)
        ggp2 <- lapply(dfs2, ggp_box)
    } else {
        ggp <- lapply(dfs, ggp_box, 0.1, 0, pvars$trans)
        ggp2 <- lapply(dfs2, ggp_box, 0.1, 0, pvars$trans)
    }
    target <- paste0(pvars$proj_dir, pvars$out_prefix, "_wells.pdf")
    if (pvars$out_prefix == "neg") {
        ggp_write(ggp[5:12], target, nrow=2, ncol=4, ht=6, wd=12)
        ggp_write(
            ggp[c(1:4, 13:22)], sub(".pdf", "2.pdf", target),
            nrow=4, ncol=4, ht=12, wd=12
        )
    } else {
        ggp_write(ggp, target, nrow=1, ncol=2, ht=3, wd=6)
    }
    target <- paste0(pvars$proj_dir, pvars$out_prefix, "_means.pdf")
    if (pvars$out_prefix == "neg") {
        ggp_write(ggp2[5:12], target, nrow=2, ncol=4, ht=6, wd=12)
        ggp_write(
            ggp2[c(1:4, 13:22)], sub(".pdf", "2.pdf", target),
            nrow=4, ncol=4, ht=12, wd=12
        )
    } else {
        ggp_write(ggp2, target, nrow=1, ncol=2, ht=3, wd=6)
    }
    invisible(
        lapply(seq.int(from=1, to=length(ggp), by=2), FUN <- function(i) {
            title <- sub(" ", "_", tolower(ggp[[i]]$labels$title))
            target <- paste0(
                pvars$proj_dir, pvars$out_prefix, "_wells_", title, ".pdf"
            )
            if (pvars$out_prefix == "neg") {
                ggp_write(ggp[i:(i + 1)], target)
            } else {
                ggp_write(ggp[i:(i + 1)], target, nrow=1, ncol=2, ht=3.5, wd=7)
            }
        })
    )
    invisible(
        lapply(seq.int(from=1, to=length(ggp2), by=2), FUN <- function(i) {
            title <- sub(" ", "_", tolower(ggp2[[i]]$labels$title))
            target <- paste0(
                pvars$proj_dir, pvars$out_prefix, "_means_", title, ".pdf"
            )
            if (pvars$out_prefix == "neg") {
                ggp_write(ggp2[i:(i + 1)], target)
            } else {
                ggp_write(ggp2[i:(i + 1)], target, nrow=1, ncol=2, ht=3.5, wd=7)
            }
        })
    )

    ## make table
    target <- paste0(pvars$proj_dir, pvars$out_prefix, "_pvals.tsv")
    write_gee(sums, target)
}

args <- commandArgs(trailingOnly=TRUE)
main(args)
