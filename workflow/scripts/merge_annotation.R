#!/usr/bin/env Rscript

# Write stdout and stderr to log file
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, type = "output")

library(tidyverse)
library(logger)
log_threshold(INFO)

log_info("Reading variants table, replacing REGION with the reference name")
variants <- read_tsv(
  snakemake@input$tsv,
  col_types = list(
    CHROM = col_character(),
    POS = col_integer(),
    REF = col_character(),
    ALT = col_character()
  )
) %>%
  mutate(REGION = snakemake@params$ref_name)

log_info("Reading annotation table")
annotation <- read_tsv(
  snakemake@input$annot,
  col_select = c("CHROM", "POS", "REF", "ALT", "VARIANT_NAME"),
  col_types = list(
    CHROM = col_character(),
    POS = col_integer(),
    REF = col_character(),
    ALT = col_character(),
    VARIANT_NAME = col_character()
  )
) %>%
  distinct() %>%
  group_by(CHROM, POS, REF, ALT) %>%
  mutate(VARIANT_NAME = paste(VARIANT_NAME, collapse = "|")) %>%
  ungroup()

log_info("Merging tables")
merged <- left_join(
  variants,
  annotation,
  by = c("REGION" = "CHROM", "POS", "REF", "ALT")
)

log_info("Saving results")
write_tsv(
  merged,
  snakemake@output$tsv
)
