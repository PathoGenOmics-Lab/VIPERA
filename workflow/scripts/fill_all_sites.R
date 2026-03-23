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

# Select variables that should be constant across samples
variant_meta <- variants %>%
  select(VARIANT_NAME, CHROM, POS, REF, ALT,
    EFFECT, IMPACT, BIOTYPE, GENE, GENEID,
    FEATURE, FEATUREID, HGVS_P, HGVS_C, ERRORS) %>%
  distinct()

log_info("Reading filtered sites")
sites <- read_tsv(snakemake@input[["sites"]]) %>%
  distinct(SAMPLE, POS) %>%  # TODO: consider region/chrom
  mutate(FILTER_PASS = TRUE)

log_info("Processing variants")
all_variants <- variants %>%
  # Keep sample-level measurement columns alongside keys
  distinct(VARIANT_NAME, CHROM, SAMPLE, ALT_FREQ, RAW_DEPTH, STRAND_BIAS) %>%
  # Handle duplicates
  group_by(SAMPLE, VARIANT_NAME, CHROM) %>%
  summarise(
    ALT_FREQ = sum(ALT_FREQ, na.rm = TRUE),
    RAW_DEPTH = sum(RAW_DEPTH, na.rm = TRUE),
    STRAND_BIAS = first(STRAND_BIAS),
    .groups = "drop"
  ) %>%
  # Complete with NA
  complete(SAMPLE, nesting(VARIANT_NAME, CHROM)) %>%
  # Restore variant-level annotations (also brings POS back in)
  left_join(variant_meta, by = c("VARIANT_NAME", "CHROM")) %>%
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
