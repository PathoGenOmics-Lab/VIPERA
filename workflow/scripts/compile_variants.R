#!/usr/bin/env Rscript

# Write stdout and stderr to log file
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, type = "output")

library(readr)
library(dplyr)
library(logger)

log_threshold(INFO)

log_info("Reading input partial variant tables and adding sample column")
tables <- lapply(
  snakemake@input,
  function(path) {
    read_tsv(path) %>%
      mutate(SAMPLE = sub("\\.variants.tsv$", "", basename(path)))
  }
)

log_info("Binding and writing table")
bind_rows(tables) %>%
  write_tsv(snakemake@output$tsv)
