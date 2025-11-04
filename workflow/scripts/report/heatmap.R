#!/usr/bin/env Rscript

# Write stdout and stderr to log file
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, type = "output")

library(dplyr)
library(readr)
library(tidyr)
library(tibble)
library(stringr)
library(logger)

log_threshold(INFO)

variants <- read_tsv(snakemake@input[["variants"]])

# Obtain sample names ordered by CollectionDate
date_order <- read_csv(snakemake@input[["metadata"]]) %>%
  arrange(CollectionDate) %>%
  pull(ID) %>%
  unique()

# Create SNP variable and select useful variables
variants <- variants %>% select(VARIANT_NAME, SAMPLE, ALT_FREQ)

variants <- variants %>%
  pivot_wider(
    names_from = VARIANT_NAME,
    values_from = ALT_FREQ,
    values_fill = 0,
    values_fn = sum
  ) %>%
  arrange(factor(SAMPLE, levels = date_order)) %>%
  # Removes "|"-separated annotations, keeping the first one + ellipsis (clarifies heatmap)
  rename_with(~ str_replace(., "^([^|]+)\\|.*$", "\\1(...)"), -SAMPLE) %>%
  column_to_rownames(var = "SAMPLE")

log_info("Saving table to create heatmap")
write.csv(variants, snakemake@output[["table"]])
