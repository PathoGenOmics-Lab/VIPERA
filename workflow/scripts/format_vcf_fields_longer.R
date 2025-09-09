#!/usr/bin/env Rscript

# Write stdout and stderr to log file
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, type = "output")

library(tidyverse)

empty.to.na <- function(x) {
    x[x == ""] <- NA
    x
}

# Replace "" values with NA in R filter list
# Snakemake passes filters like: list(ERRORS = c(""))
filter.include <- lapply(snakemake@params$filter_include, empty.to.na)
filter.exclude <- lapply(snakemake@params$filter_exclude, empty.to.na)

# Process input table
read_tsv(snakemake@input$tsv) %>%
    # Separate <sep>-delimited "...[*]..." columns (e.g. ANN[*].EFFECT)
    separate_rows(
        contains("[*]"),
        sep = snakemake@params$sep,
        convert = TRUE
    ) %>%
    # Separate &-delimited error column (more than one error/warning/info message per row is possible)
    separate_rows("ANN[*].ERRORS", sep = "&") %>%
    # Rename "...[*]..." columns using the provided lookup via Snakemake config
    rename(all_of(unlist(snakemake@params$colnames_mapping))) %>%
    # Apply dynamic filters from the Snakemake config:
    # map2 pairs column names (.x) with value vectors (.y) and builds boolean expressions.
    # Inside the expr call, !! injects a single value into each expression.
    # The resulting list of expressions is spliced with !!! so each becomes its
    # own condition as if written directly inside the filter call.
    filter(
        # Keep variants that include the required values in each defined field (e.g. empty ERRORS)
        !!!map2(
            names(filter.include),
            filter.include,
            ~ expr(.data[[!!.x]] %in% !!.y)
        ),
        # Keep variants that exclude the required values in each defined field  (e.g. EFFECT != "upstream_gene_variant")
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
