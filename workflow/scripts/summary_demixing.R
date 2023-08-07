#!/usr/bin/env Rscript
# Jordi Sevilla

library(tidyverse)

# Write stdout and stderr to log file
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, type = "output")

# Crear dataframe vacío donde contener los datos
demix <- data.frame("lineages" = NA, "abundances" = NA, "sample" = NA) %>%
  filter(!is.na(sample))

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
