#!/usr/bin/env Rscript
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

log_info("Sorting dates")
date_order <- read_csv(snakemake@input[["metadata"]]) %>%
  arrange(CollectionDate) %>%
  pull(ID) %>%
  unique()

log_info("Formatting variants")
all_variants_wider <- variants %>%
  # Collapse positions (treat as uncertain if filter is inconsistent)
  group_by(SAMPLE, VARIANT_NAME) %>%
  summarise(
    ALT_FREQ = ifelse(n_distinct(ALT_FREQ, na.rm = TRUE) > 1, NA, first(ALT_FREQ)),
    FILTER_PASS = ifelse(n_distinct(FILTER_PASS) > 1, NA, first(FILTER_PASS)),
    .groups = "drop"
  ) %>%
  mutate(ALT_FREQ = if_else(is.na(FILTER_PASS), NA, ALT_FREQ)) %>%
  distinct(SAMPLE, VARIANT_NAME, ALT_FREQ) %>%
  pivot_wider(
    names_from = VARIANT_NAME,
    values_from = ALT_FREQ
  ) %>%
  # Apply chronological ordering
  arrange(factor(SAMPLE, levels = date_order)) %>%
  # Removes "|"-separated annotations, keeping the first one + ellipsis
  rename_with(~ str_replace(., "^([^|]+)\\|.*$", "\\1(...)"), -SAMPLE) %>%
  column_to_rownames(var = "SAMPLE")

log_info("Saving table of frequencies")
write.csv(all_variants_wider, snakemake@output[["table"]])

log_info("Filtering variant columns by unique values and amplitude thresholds")
filtered <- all_variants_wider %>%
  select(
    where(
      function(col) {
        clean_col <- na.omit(col)
        # Drop empty
        if (length(clean_col) == 0) return(FALSE)
        # Unique values check
        pass_unique <- length(unique(clean_col)) >= snakemake@params$min_unique_n
        # Amplitude check
        amplitude <- max(clean_col) - min(clean_col)
        pass_amp <- amplitude >= snakemake@params$min_freq_amplitude
        return(pass_unique && pass_amp)
      }
    )
  )

log_info("Calculating correlations")
cor.mat <- cor(
  filtered,
  method = snakemake@params$cor_method,
  use = snakemake@params$cor_use
)

log_info("Saving correlation matrix")
write.csv(cor.mat, snakemake@output[["matrix"]])
