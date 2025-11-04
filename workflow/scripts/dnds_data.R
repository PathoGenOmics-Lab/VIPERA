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

# Read inputs
log_info("Reading variants table")
variants <- read_delim(
  snakemake@input[["variants"]],
  col_select = c(
    "SAMPLE",
    "VARIANT_NAME",
    "REGION",
    "ALT_FREQ",
    "SYNONYMOUS",
    "POS"
  )
)

log_info("Reading metadata table")
metadata <- read_delim(snakemake@input[["metadata"]]) %>%
  mutate(
    interval = as.numeric(
      as.Date(CollectionDate) - min(as.Date(CollectionDate))
    )
  ) %>%
  select(ID, interval) %>%
  rename(SAMPLE = ID)

log_debug("Adding metadata to variants table")
variants <- left_join(variants, metadata)

log_info("Reading N/S sites")
N_S_position <- read_delim(snakemake@input[["n_s_sites"]])

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
    dn = No / sum(N_S_position$N),
    ds = Yes / sum(N_S_position$S)
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
