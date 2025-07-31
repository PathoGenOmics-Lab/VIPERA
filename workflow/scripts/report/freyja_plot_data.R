#!/usr/bin/env Rscript

# Write stdout and stderr to log file
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, type = "output")

library(tidyverse)
library(logger)
log_threshold(INFO)

# Read inputs
demix <- read_csv(snakemake@input[["summary_demixing"]])

# Data processing
log_info("Obtaining main lineages")
main_lineages <- demix %>%
  group_by(sample) %>%
  top_n(1, abundance) %>%
  ungroup() %>%
  pull(lineage) %>%
  unique()

# Obtain sample names ordered by CollectionDate
log_info("Sorting dates")
metadata <- read_csv(snakemake@input[["metadata"]])
date_order <- metadata %>%
  arrange(CollectionDate) %>%
  filter(ID %in% demix$sample) %>%
  pull(ID) %>%
  unique()

# Build plot data
log_info("Building plot data")
plot.data <- demix %>%
  group_by(lineage, sample) %>%
  summarize(abundance = sum(abundance)) %>%
  ungroup() %>%
  left_join(
    select(
      metadata,
      ID,
      CollectionDate
    ),
    by = c("sample" = "ID")
  ) %>%
  mutate(
    is_main = lineage %in% main_lineages
  )

# Save plot data
log_info("Saving plot data")
write_csv(
  plot.data,
  snakemake@output[["data"]]
)
