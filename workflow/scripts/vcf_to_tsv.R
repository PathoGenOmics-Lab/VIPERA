#!/usr/bin/env Rscript

library(tidyverse)
library(logger)
log_threshold(INFO)

# Write stdout and stderr to log file
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, type = "output")

# read data
log_info("Reading data")
vcf <- read_tsv(snakemake@input[["ann_vcf"]], comment = "##")
tsv <- read_tsv(snakemake@input[["pre_tsv"]])

tsv["variant"] <- str_extract(vcf$INFO, "p\\.[^|]*")

tsv <- tsv %>%
    mutate(
        variant = case_when(
            is.na(variant) ~ paste(POS, REF, ">", ALT, sep = ""),
            TRUE ~ paste(GFF_FEATURE, ":", variant, sep = "")
        )
    )

log_info("Saving results")
write_tsv(
    tsv,
    snakemake@output[["tsv"]]
)
