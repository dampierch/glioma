## to be sourced from main analysis script

get_data <- function(target) {
    ## -- reads luciferase luminescence data from tab-separated values text
    ## -- normalizes luminescence with negative control
    ## -- returns list of tbls, one for each experiment/dataset

    ## task 1: read data
    if (file.exists(target)) {
        tbl <- read_tsv(target)
    } else {
        stop(target, "does not exist")
    }

    ## task 2: normalize intensity to negative control
    dfs <- lapply(unique(tbl$data_batch), FUN <- function(batch) {
        df <- tbl[tbl$data_batch == batch, ]
        dfs <- lapply(unique(df$plate_id), FUN2 <- function(plate) {
            df2 <- df[df$plate_id == plate, ]
            df2$rel_lum <- df2$intensity / df2$plate_neg_cntrl_mean
            return(df2)
        })
        df <- do.call(rbind, dfs)
        return(df)
    })
    return(dfs)
}

summarise_data <- function(df) {
    ## use dplyr to summarize data by clone_id
    df <- df %>% group_by(plate_id, clone_id) %>% summarise(
        locus=first(locus),
        rsid=first(rsid),
        orientation=first(orientation),
        gene=first(gene),
        cell_line=first(cell_line),
        allele=first(allele),
        rel_lum=mean(rel_lum),
        .groups="keep"
    )
    return(df)
}
