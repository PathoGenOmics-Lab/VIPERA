#!/usr/bin/env Rscript

# Write stdout and stderr to log file
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, type = "output")

library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(rlang)
library(glue)
library(logger)

empty.to.na <- function(x) {
  x[x == ""] <- NA
  x
}

# Replace "" values with NA in R filter list
# Snakemake passes filters like: list(ERRORS = c(""))
log_info("Preprocessing data")
filter.include <- lapply(snakemake@params$filter_include, empty.to.na)
filter.exclude <- lapply(snakemake@params$filter_exclude, empty.to.na)

# Process input table
log_info("Applying filters and writing results")
read_tsv(snakemake@input$tsv) %>%

  # Separate <sep>-delimited "...[*]..." columns (e.g. ANN[*].EFFECT)
  separate_rows(
    contains("[*]"),
    sep = snakemake@params$sep,
    convert = TRUE
  ) %>%

  # Rename "...[*]..." columns using the provided lookup via Snakemake config
  rename(all_of(unlist(snakemake@params$colnames_mapping))) %>%

  # Separate &-delimited error column (more than one error/warning/info message per row is possible)
  mutate(split_errors = strsplit(ERRORS, "&")) %>%
  # Keep rows with none of the excluded ERRORS terms, if any
  filter(map_lgl(split_errors, ~ !any(. %in% filter.exclude[["ERRORS"]]))) %>%
  select(-split_errors) %>%

  # Apply filters
  filter(
    # Keep variants that include required values in each field
    !!!map2(
      names(filter.include),
      filter.include,
      ~ expr(.data[[!!.x]] %in% !!.y)
    ),
    # Keep variants that exclude required values in each field
    !!!map2(
      names(filter.exclude),
      filter.exclude,
      ~ expr(!(.data[[!!.x]] %in% !!.y))
    )
  ) %>%

  # Keep unique rows
  distinct() %>%

  # Assign variant name using the pattern defined via Snakemake config
  mutate(VARIANT_NAME = str_glue(snakemake@params$variant_name_pattern)) %>%

  # Write output file
  write_tsv(snakemake@output$tsv)
