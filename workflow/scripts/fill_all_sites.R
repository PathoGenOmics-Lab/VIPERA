#!/usr/bin/env Rscript

# Write stdout and stderr to log file
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, type = "output")

library(dplyr)
library(readr)
library(tidyr)
library(logger)

log_threshold(INFO)

log_info("Reading variants")
variants <- read_tsv(snakemake@input[["variants"]])

# Create a mapping of variant names to their genomic position
variant_coords <- variants %>%
  select(VARIANT_NAME, REGION, POS) %>%
  distinct()

log_info("Reading filtered sites")
sites <- read_tsv(snakemake@input[["sites"]]) %>%
  select(SAMPLE, POS) %>%  # TODO: consider region/chrom
  distinct() %>%
  mutate(FILTER_PASS = TRUE)

log_info("Processing variants")
all_variants <- variants %>%
  # Select minimal columns
  select(VARIANT_NAME, REGION, SAMPLE, ALT_FREQ) %>%
  # Handle duplicates
  group_by(SAMPLE, VARIANT_NAME, REGION) %>%
  summarise(ALT_FREQ = sum(ALT_FREQ, na.rm = TRUE), .groups = "drop") %>%
  # Complete with NA
  complete(SAMPLE, VARIANT_NAME, REGION) %>%
  # Assign genomic positions for all combinations
  left_join(variant_coords, by = c("REGION", "VARIANT_NAME")) %>%
  # Merge filtered sites
  # TODO: consider region/chrom
  left_join(sites, by = c("SAMPLE", "POS")) %>%
  replace_na(list(FILTER_PASS = FALSE)) %>%
  # Fill missing frequencies conditionally
  mutate(
    ALT_FREQ = case_when(
      !is.na(ALT_FREQ) ~ ALT_FREQ,  # variant was called, keep the frequency
      FILTER_PASS ~ 0,              # missing but site is covered: 0 (reference)
      TRUE ~ NA                     # missing and not covered: NA (unknown)
    )
  )

log_info("Saving table")
write_tsv(all_variants, snakemake@output[["variants"]])
