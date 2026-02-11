#!/usr/bin/env Rscript

# Write stdout and stderr to log file
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, type = "output")

library(readr)
library(dplyr)
library(jsonlite)
library(logger)

log_threshold(INFO)

log_info("Reading variants")
variants <- read_delim(snakemake@input$variants)

log_info("Reading metadata")
metadata <- read_delim(snakemake@input$metadata)

log_info("Calculating heterozygous sites")
sites <- variants %>%
  filter(ALT_FREQ <= snakemake@params$max_alt_freq) %>%
  left_join(
    metadata,
    by = c("SAMPLE" = "ID")
  ) %>%
  group_by(SAMPLE) %>%
  summarise(
    CollectionDate = min(as.Date(CollectionDate)),
    n = n_distinct(POS)
  ) %>%
  ungroup() %>%
  arrange(CollectionDate) %>%
  mutate(
    Day = as.numeric(
      difftime(CollectionDate, min(CollectionDate), units = "days")
    )
  )

if (nrow(sites) == 0) {
  log_warn("There are none, using an empty table and no linear regression")
  sites <- tibble(
    SAMPLE = date_order,
    CHROM = as.character(NA),
    VARIANT_NAME = as.character(NA),
    ALT_FREQ = as.numeric(NA),
    EFFECT = as.character(NA),
    SYNONYMOUS = as.character(NA),
    POS = as.numeric(NA),
    ALT = as.character(NA),
    NV_class = as.character(NA),
    group = as.character(NA)
  )
  r_squared <- "none"
  p_value_string <- "none"
} else if (nrow(sites) > 2) {
  log_info("Calculating linear regression")
  model <- lm(n ~ CollectionDate, data = sites)
  r_squared <- summary(model)$r.squared[[1]]
  p_value <- summary(model)$coefficients[2, 4]
  p_value_string <- ifelse(p_value < 0.001, "< 0.001", p_value)
} else {
  log_warn("Not enough data points for a linear regression")
  r_squared <- "none"
  p_value_string <- "none"
}

log_info("Writing JSON summary")
list(
  "r2" = r_squared,
  "value" = p_value_string
) %>%
  write_json(
    snakemake@output$json,
    auto_unbox = TRUE,
    digits = NA
  )

log_info("Writing processed table")
write_csv(sites, snakemake@output$table)
