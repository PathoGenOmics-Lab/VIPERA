#!/usr/bin/env Rscript

# Write stdout and stderr to log file
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, type = "output")

library(tidyverse)

read_tsv(snakemake@input$tsv) %>%
    separate_rows(contains("[*]"), sep = snakemake@params$sep) %>%
    rename(all_of(unlist(snakemake@params$colnames_mapping))) %>%
    filter(
        !!!map2(
            names(snakemake@params$filter_include),
            snakemake@params$filter_include,
            ~ expr(.data[[!!.x]] %in% !!.y)
        ),
        !!!map2(
            names(snakemake@params$filter_exclude),
            snakemake@params$filter_exclude,
            ~ expr(!(.data[[!!.x]] %in% !!.y))
        )
    ) %>%
    distinct() %>%
    mutate(VARIANT_NAME = str_glue(snakemake@params$variant_name_pattern)) %>%
    write_tsv(snakemake@output$tsv)
