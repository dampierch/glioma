## to be sourced from main analysis script

ggp_box <- function(data, jwd=0, jht=0, trans=NA) {
    ## relative luminescence is lum of putative element divided by lum of neg
    ## y-axis should include interval from 0 to 1 to show inactive
    ## y-axis is easier to interpret if not log transformed
    ggp_title <- paste(data$locus[1], data$orientation[1])
    ggp_subtitle <- data$cell_line[1]
    ggp_xlab <- gsub("_", ", ", data$rsid[1])
    ggp_ylab <- "Relative Luminescence"
    ggp_legend <- "Experiment"
    ggp_shape <- c(16, 15, 17)
    ggp <- ggplot(data, aes(x=allele, y=rel_lum)) +
        geom_boxplot(outlier.size=-1, width=0.3) +
        stat_summary(
            fun=mean, geom="point", shape=23, size=5, color="grey40",
            stroke=1.5
        ) +
        geom_point(
            aes(shape=factor(plate_id)), size=2, alpha=0.5, stroke=1,
            position=position_jitter(width=jwd, height=jht)
        ) +
        labs(title=ggp_title, subtitle=ggp_subtitle, x=ggp_xlab, y=ggp_ylab) +
        scale_shape_manual(values=ggp_shape) +
        ggp_theme_default
    if (is.na(trans)) {
        ggp <- ggp +
            scale_y_continuous(limits=c(0, max(data$rel_lum) * 1.1))
    } else {
        ggp_ylab <- expression("Relative Luminescence "*"(Log"["2"]*" Scaled)")
        ggp <- ggp +
            labs(
                title=ggp_title, subtitle=ggp_subtitle, x=ggp_xlab, y=ggp_ylab
            ) +
            scale_y_continuous(
                trans=trans, limits=c(0.9, max(data$rel_lum) * 1.1)
            )
    }
    if (length(unlist(strsplit(data$rsid[1], "_"))) > 3) {
        ## when rsIDs do not fit in x label space
        ggp <- ggp +
            theme(axis.title.x=element_text(size=8, face="plain"))
    }
    return(ggp)
}

ggp_annotate <- function(ggp, pval) {
    ## not utilized in current workflow
    ptxt <- format(round(pval[1, 1], 2), digits=2)
    ggp_ann <- bquote(italic("p") == .(ptxt))
    if (nrow(pval) == 2) {
        ptxt2 <- format(round(pval[2, 1], 2), digits=2)
        ggp_ann2 <- bquote(italic("p") == .(ptxt2))
    }
    if (nrow(pval) == 1) {
        ggp <- ggp +
            annotate(
                "text", x=1.5, y=1.1 * max(data$rel_lum), label=ggp_ann, size=3
            )
    } else if (nrow(pval) == 2) {
        ggp <- ggp +
            annotate(
                "text", x=1.5, y=1.1 * max(data$rel_lum), label=ggp_ann, size=3
            ) +
            annotate(
                "text", x=2.5, y=1.1 * max(data$rel_lum), label=ggp_ann2, size=3
            )
    }
}

ggp_write <- function(ggp, target, nrow=1, ncol=2, ht=3, wd=6) {
    ## default is for a pair of 3x3in plots
    cwp <- plot_grid(
        plotlist=ggp, nrow=nrow, ncol=ncol, labels=LETTERS[seq_len(length(ggp))]
    )
    pdf(target, height=ht, width=wd)
    print(cwp)
    dev.off()
    cat("plot written to", target, "\n")
}
