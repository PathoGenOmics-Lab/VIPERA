#!/usr/bin/env Rscript

library(tidyverse)
library(logger)
log_threshold(INFO)


# Write stdout and stderr to log file
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, type = "output")


metadata <- read.csv(snakemake@input[["metadata"]])
pango_report <- read.csv(snakemake@input[["report"]])


# Obtain sample names ordered by CollectionDate
date_order <- metadata %>%
    arrange(CollectionDate) %>%
    filter(
        ID %in% pango_report$taxon
    ) %>%
pull(ID) %>%
unique()

# Create a temporal index for samples
index  <- data.frame(
    Sample = date_order,
    Index = seq(1, length(date_order), 1)
    )

metadata <- select(
        metadata,
        ID,
        CollectionDate
    ) %>%
    filter(
        ID %in% pango_report$taxon
    ) %>%
    left_join(
        select(pango_report, taxon, lineage),
        by = c("ID" = "taxon")
    ) %>%
    rename(
        Sample = ID,
        Collection_Date = CollectionDate,
        Lineage = lineage
    )

log_info("Creating summary table")
table <- left_join(index, metadata)

log_info("Saving table")
write.csv(table, snakemake@output[["table"]], row.names = FALSE)
