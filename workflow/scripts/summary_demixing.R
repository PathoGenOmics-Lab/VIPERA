#!/usr/bin/env Rscript

# Write stdout and stderr to log file
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, type = "output")

library(tidyverse)
library(logger)
log_threshold(INFO)

demix <- data.frame( # Empty dataframe to be fille with data
  "lineages" = NA,
   "abundances" = NA,
   "sample" = NA
   ) %>%
  filter(!is.na(sample))

log_info("Summarizing demixing tables")
lapply(
  snakemake@input[["tables"]],
  function(tsv_file) {
    read_tsv(
        tsv_file,
        col_names = c("variable", "valor"),
        show_col_types = FALSE
      ) %>%
      filter(
        row_number() %in% c(3, 4)
      ) %>%
      pivot_wider(
        names_from = variable,
        values_from = valor
      ) %>%
      separate_rows(
        lineages,
        abundances,
        sep = " "
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
