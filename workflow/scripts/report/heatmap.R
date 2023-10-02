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
  mutate(
    SNP = case_when(
      !is.na(REF_AA) ~ paste(
        GFF_FEATURE,
        ":",
        REF_AA,
        POS_AA,
        ALT_AA,
        sep = ""
      ),
      TRUE ~ paste(REF, POS, ALT, sep = "")
    )
  ) %>%
  unique() %>%
  dplyr::select(SNP, REGION, ALT_FREQ)

vcf <- vcf %>%
  pivot_wider(
    names_from = SNP,
    values_from = ALT_FREQ,
    values_fill = 0
  ) %>%
  arrange(factor(REGION, levels = date_order)) %>%
  column_to_rownames(var = "REGION")

log_info("Saving table to create heatmap")
write.csv(vcf, snakemake@output[["table"]])
