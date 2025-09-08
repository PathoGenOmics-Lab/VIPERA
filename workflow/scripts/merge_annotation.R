#!/usr/bin/env Rscript

# Write stdout and stderr to log file
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, type = "output")

library(tidyverse)
library(logger)
log_threshold(INFO)

log_info("Reading tables")
variants <- read_tsv(snakemake@input$tsv)
annotation <- read_tsv(
  snakemake@input$annot,
  col_select = c("CHROM", "POS", "REF", "ALT", "VARIANT_NAME")
)

log_info("Merging tables")
merged <- left_join(
  variants,
  annotation,
  by = c("CHROM", "POS", "REF", "ALT")
)

log_info("Saving results")
write_tsv(
  tsv,
  snakemake@output$tsv
)
