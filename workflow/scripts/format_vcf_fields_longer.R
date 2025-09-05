#!/usr/bin/env Rscript

# Write stdout and stderr to log file
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, type = "output")

library(tidyverse)

read_tsv(snakemake@input$tsv) %>%
    separate_rows(contains("[*]"), sep = snakemake@params$sep) %>%
    write_tsv(snakemake@output$tsv)
