#!/usr/bin/env Rscript

# Write stdout and stderr to log file
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, type = "output")

library(tidyverse)
library(logger)
log_threshold(INFO)

# Empty dataframe to be filled with data
demix <- data.frame(
  "lineage" = NA,
  "abundance" = NA,
  "sample" = NA
) %>%
  filter(!is.na(sample))

log_info("Summarizing demixing tables")
lapply(
  snakemake@input[["tables"]],
  function(tsv_file) {
    read_tsv(
      tsv_file,
      col_names = c("variable", "value"),
      show_col_types = FALSE
    ) %>%
      filter(
        variable %in% c("lineages", "abundances")
      ) %>%
      pivot_wider(
        names_from = variable,
        values_from = value
      ) %>%
      separate_rows(
        lineages,
        abundances,
        sep = " "
      ) %>%
      rename(
        lineage = lineages,
        abundance = abundances
      ) %>%
      mutate(
        sample = str_extract(
          basename(tsv_file),
          "(.+)_demixed.tsv",
          group = 1
        )
      )
  }
) %>%
  bind_rows %>%
  write_csv(snakemake@output[["summary_df"]])
