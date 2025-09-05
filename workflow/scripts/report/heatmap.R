#!/usr/bin/env Rscript

# Write stdout and stderr to log file
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, type = "output")

library(tidyverse)
library(logger)
log_threshold(INFO)

vcf <- read_tsv(snakemake@input[["vcf"]])

# Obtain sample names ordered by CollectionDate
date_order <- read_csv(snakemake@input[["metadata"]]) %>%
  arrange(CollectionDate) %>%
  pull(ID) %>%
  unique()

# Create SNP variable and select useful variables from vcf
vcf <- vcf %>%
  dplyr::select(variant, SAMPLE, ALT_FREQ)

vcf <- vcf %>%
  pivot_wider(
    names_from = variant,
    values_from = ALT_FREQ,
    values_fill = 0,
    values_fn = sum
  ) %>%
  arrange(factor(SAMPLE, levels = date_order)) %>%
  column_to_rownames(var = "SAMPLE")

log_info("Saving table to create heatmap")
write.csv(vcf, snakemake@output[["table"]])
