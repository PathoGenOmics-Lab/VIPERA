#!/usr/bin/env Rscript
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

# Check if each sample, variant and position have 1 frequency
variants %>%
  distinct(VARIANT_NAME, CHROM, POS, SAMPLE, ALT_FREQ) %>%
  group_by(VARIANT_NAME, CHROM, POS, SAMPLE) %>%
  filter(n() > 1) %>%
  {
    if (nrow(.) > 0) {
      log_warn(
        "Found {nrow(.)} ambiguous (SAMPLE, VARIANT_NAME, POS) combinations"
      )
    }
  }

log_info("Reading filtered sites")
sites <- read_tsv(snakemake@input[["sites"]]) %>%
  distinct(SAMPLE, POS) %>%  # TODO: consider region/chrom
  mutate(FILTER_PASS = TRUE)

log_info("Processing variants")
all_variants <- variants %>%
  # Select minimal columns
  distinct(VARIANT_NAME, CHROM, POS, SAMPLE, ALT_FREQ) %>%
  # Complete with NA
  complete(SAMPLE, nesting(VARIANT_NAME, CHROM, POS)) %>%
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
