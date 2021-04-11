## to be sourced from main analysis script

gee_multi_plate <- function(data, levels) {
    fit <- geeglm(
        rel_lum ~ factor(allele, levels=levels) + factor(plate_id),
        data=data,
        id=factor(plate_id):factor(clone_id),
        family=gaussian(link="log"),
        corstr="exchangeable",
        std.err="san.se",
        zcor=NULL
    )
    return(fit)
}

gee_single_plate <- function(data, levels) {
    fit <- geeglm(
        rel_lum ~ factor(allele, levels=levels),
        data=data,
        id=factor(clone_id),
        family=gaussian(link="log"),
        corstr="exchangeable",
        std.err="san.se",
        zcor=NULL
    )
    return(fit)
}

fit_gee <- function(data) {
    if (length(levels(factor(data$plate_id))) > 1) {
        cat("testing multi plate experiment\n")
        if (length(levels(factor(data$allele))) == 3) {
            cat("testing tri alleleic site\n")
            levs <- levels(factor(data$allele))
            fit1 <- gee_multi_plate(data, levs)
            idx <- data$allele %in% levs[2:3]
            fit2 <- gee_multi_plate(data[idx, ], levs[2:3])
            fit <- setNames(list(fit1, fit2), c("default", "extra"))
        } else {
            cat("testing bi alleleic site\n")
            levs <- levels(factor(data$allele))
            fit <- gee_multi_plate(data, levs)
        }
    } else {
        cat("testing single plate experiment\n")
        if (length(levels(factor(data$allele))) == 3) {
            cat("testing tri alleleic site\n")
            levs <- levels(factor(data$allele))
            fit1 <- gee_single_plate(data, levs)
            idx <- data$allele %in% levs[2:3]
            fit2 <- gee_single_plate(data[idx, ], levs[2:3])
            fit <- setNames(list(fit1, fit2), c("default", "extra"))
        } else {
            cat("testing bi alleleic site\n")
            levs <- levels(factor(data$allele))
            fit <- gee_single_plate(data, levs)
        }
    }
    return(fit)
}

fit_parser <- function(fit) {
    df <- summary(fit)$coefficients
    df <- df[2:nrow(df), ]
    colnames(df) <- c("Effect", "SE", "Wald_Stat", "pval")
    df$Locus_Cell <- paste(unname(unlist(
        fit$data[1, c("locus", "orientation", "cell_line")]
    )), collapse="_")
    df$A1 <- levels(fit$model[ , "factor(allele, levels = levels)"])[1]
    df$A2 <- unlist(lapply(strsplit(rownames(df), ")"), "[[", 2))
    df$Label <- paste(df$Locus_Cell, paste(df$A1, df$A2, sep="_"), sep="::")
    col_order <- c(
        "Label", "Locus_Cell", "A1", "A2", "Effect", "SE", "Wald_Stat", "pval"
    )
    df <- df[col_order]
    i <- suppressWarnings(which(is.na(as.numeric(df$A2))))
    df <- df[i, ]
    return(df)
}

parse_gee <- function(fit) {
    if (class(fit)[1] == "list") {
        cat("extracting coefficient summary for tri alleleic site\n")
        df <- lapply(fit, fit_parser)
    } else {
        cat("extracting coefficient summary for bi alleleic site\n")
        df <- fit_parser(fit)
    }
    return(df)
}

write_gee <- function(sums, target) {
    obj <- lapply(sums, FUN <- function(sum) {
        if (class(sum) == "list") {
            df <- do.call(rbind, sum)
        } else {
            df <- sum
        }
        return(df)
    })
    df <- do.call(rbind, obj)
    write.table(
        df, file=target, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE
    )
    cat("table written to", target, "\n")
}
