#!/usr/bin/env Rscript

library(tidyverse)
library(logger)
log_threshold(INFO)


vcf <- read_tsv(snakemake@input[["vcf"]])
metadata <- read.csv(snakemake@params[["metadata"]])

# Obtain sample names ordered by CollectionDate
date_order <- read_csv(snakemake@params[["metadata"]]) %>%
  arrange(CollectionDate) %>%
  pull(ID) %>%
  unique()

# Create SNP variable and select useful variables from vcf
vcf <- vcf %>%
  dplyr::select(variant, REGION, ALT_FREQ)

vcf <- vcf %>%
  pivot_wider(
    names_from = variant,
    values_from = ALT_FREQ,
    values_fill = 0,
    values_fn = sum
  ) %>%
  arrange(factor(REGION, levels = date_order)) %>%
  column_to_rownames(var = "REGION")

log_info("Saving table to create heatmap")
write.csv(vcf, snakemake@output[["table"]])
