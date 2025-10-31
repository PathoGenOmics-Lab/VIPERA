#!/usr/bin/env Rscript

# Write stdout and stderr to log file
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, type = "output")

library(tidyverse)
library(jsonlite)
library(logger)
log_threshold(INFO)

log_info("Reading variants")
variants <- read_delim(
  snakemake@input$variants,
  col_select = c(
    "SAMPLE",
    "REGION",
    "POS",
    "ALT",
    "VARIANT_NAME",
    "ALT_FREQ",
    "EFFECT",
    "SYNONYMOUS"
  )
)
log_debug("Read {nrow(variants)} rows")

log_info("Reading dates from metadata")
dates <- read_delim(
  snakemake@input$metadata,
  col_select = c("ID", "CollectionDate")
)

log_info("Dating variants")
variants <- left_join(
  variants, dates,
  by = c("SAMPLE" = "ID")
)

log_info("Classifying variants")
variants <- variants %>%
  mutate(
    NV_class = ifelse(
      str_detect(ALT, fixed("-")) | str_detect(ALT, fixed("+")),
      "INDEL",
      "SNP"
    ),
    Class = ifelse(EFFECT == "intergenic_region", "Intergenic", SYNONYMOUS),
    POS = as.numeric(POS)
  ) %>%
  # rowwise() %>%
  mutate(
    indel_len = ifelse(NV_class == "INDEL", str_length(ALT) - 1, NA),
    indel_class = case_when(
      EFFECT == "intergenic_region" ~ "Intergenic",
      (NV_class == "INDEL") & (indel_len %% 3 == 0) ~ "In frame",
      (NV_class == "INDEL") & (indel_len %% 3 != 0) ~ "Frameshift"
    )
  ) %>%
  ungroup() %>%
  mutate(
    group = case_when(
      EFFECT == "intergenic_region" ~ "Intergenic",
      NV_class == "SNP" ~ Class,
      NV_class == "INDEL" ~ indel_class
    )
  ) %>%
  filter(ALT_FREQ > 0)
log_debug("Processed {nrow(variants)} rows")

if (nrow(variants) == 0) {
  log_warning("No variants found, using an empty table")
  variants <- tibble(
    SAMPLE = date_order,
    REGION = as.character(NA),
    VARIANT_NAME = as.character(NA),
    ALT_FREQ = as.numeric(NA),
    EFFECT = as.character(NA),
    SYNONYMOUS = as.character(NA),
    POS = as.numeric(NA),
    ALT = as.character(NA),
    NV_class = as.character(NA),
    group = as.character(NA)
  )
}

log_info("Writing processed table")
write_csv(variants, snakemake@output$table)

log_info("Writing NV summary")
nv_counts <- variants %>%
  distinct(VARIANT_NAME, NV_class) %>%
  count(NV_class) %>%
  deframe()
list(
  "INDELS" = nv_counts["INDEL"],
  "SNV" = nv_counts["SNP"]
) %>%
  write_json(
    snakemake@output$json,
    auto_unbox = TRUE,
    digits = NA
  )
