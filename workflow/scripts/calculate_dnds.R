#!/usr/bin/env Rscript

# Write stdout and stderr to log file
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, type = "output")

library(dplyr)
library(readr)
library(tidyr)
library(purrr)
library(logger)

log_threshold(INFO)

# Read inputs
log_info("Reading variants table")
variants <- read_delim(
  snakemake@input[["variants"]],
  col_select = c(
    "SAMPLE",
    "VARIANT_NAME",
    "ALT_FREQ",
    "SYNONYMOUS",
    "POS"
  )
)

log_info("Reading metadata table")
metadata <- read_delim(snakemake@input[["metadata"]]) %>%
  mutate(
    interval = difftime(
      as.Date(CollectionDate),
      min(as.Date(CollectionDate)),
      units = "days"
    ) |> as.numeric()
  ) %>%
  select(ID, interval) %>%
  rename(SAMPLE = ID)

log_info("Reading masked sites BED")
masked <- read_tsv(
  snakemake@input$masked,
  col_names = c("chrom", "start", "end"),
  col_types = "cii",
  comment = "#"
) %>%
  mutate(POS = map2(start + 1, end, seq.int)) %>%
  unnest(POS) %>%
  pull(POS)

log_info("Filtering variants on masked sites")
variants <- variants %>%
  filter(!POS %in% masked, !is.na(SYNONYMOUS))

log_debug("Adding metadata to variants table")
variants <- left_join(variants, metadata)

log_info("Reading N/S sites")
n_s_position <- read_delim(snakemake@input[["n_s_sites"]])

log_info("Computing dN/dS over time (NG86)")
dn.ds <- variants %>%
  group_by(SAMPLE, SYNONYMOUS) %>%
  summarise(
    Freq = sum(ALT_FREQ, na.rm = TRUE)
  ) %>%
  pivot_wider(
    names_from = SYNONYMOUS,
    values_from = Freq,
    values_fill = 0
  ) %>%
  transmute(
    dn = No / sum(n_s_position$N),
    ds = Yes / sum(n_s_position$S)
  ) %>%
  ungroup() %>%
  left_join(unique(select(variants, SAMPLE, interval))) %>%
  transmute(
    sample = SAMPLE,
    day = interval,
    dN = dn,
    dS = ds,
    w = dn / ds
  )

log_info("Writing results")
write_csv(dn.ds, snakemake@output[["table"]])
