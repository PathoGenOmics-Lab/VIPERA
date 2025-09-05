#!/usr/bin/env Rscript

# Write stdout and stderr to log file
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, type = "output")

library(tidyverse)
library(logger)
log_threshold(INFO)

# Read inputs
data <- read_tsv(snakemake@input[["tsv"]])

# Filtering
is.deletion <- str_detect(
  data$ALT,
  "^[A-Z]",
  negate = TRUE
)
strand.mask <- data$ALT_RV > snakemake@params$min_alt_rv &
  data$ALT_DP > snakemake@params$min_alt_dp
depth.mask <- (data$ALT_RV + data$ALT_DP) >= snakemake@params$min_depth

log_info("Filtering variants")
data <- filter(
  data,
  as.logical(PASS),
  depth.mask,
  strand.mask | is.deletion
)

log_info("Finding synonymous and non synonymous variants")
# Adding synonymous variable
data <- mutate(
  data,
  synonimous = case_when(
    REF_AA == ALT_AA ~ "Yes",
    TRUE ~ "No"
  )
)

# Remove duplicated features
data <- distinct(data, pick(!GFF_FEATURE), .keep_all = TRUE)

# Change annotation to gb2seq annotation
features <- read_csv(snakemake@input[["annotation"]])

data <- data %>%
  select(!GFF_FEATURE) %>%
  left_join(features) %>%
  rename(GFF_FEATURE = GEN)

log_info("Saving results")
write_tsv(
  data,
  snakemake@output[["filtered_tsv"]]
)
