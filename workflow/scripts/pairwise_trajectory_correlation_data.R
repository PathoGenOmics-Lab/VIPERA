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

log_info("Reading variants")
variants <- read_tsv(snakemake@input[["variants"]])

# Obtain sample names ordered by CollectionDate
date_order <- read_csv(snakemake@input[["metadata"]]) %>%
  arrange(CollectionDate) %>%
  pull(ID) %>%
  unique()

all_variants_wider <- variants %>%
  select(SAMPLE, VARIANT_NAME, ALT_FREQ) %>%
  pivot_wider(
    names_from = VARIANT_NAME,
    values_from = ALT_FREQ
  ) %>%
  # Apply chronological ordering
  arrange(factor(SAMPLE, levels = date_order)) %>%
  # Removes "|"-separated annotations, keeping the first one + ellipsis
  rename_with(~ str_replace(., "^([^|]+)\\|.*$", "\\1(...)"), -SAMPLE) %>%
  column_to_rownames(var = "SAMPLE")

log_info("Saving table")
write.csv(all_variants_wider, snakemake@output[["table"]])
