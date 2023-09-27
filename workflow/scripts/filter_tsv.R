#!/usr/bin/env Rscript

library(tidyverse)
library(logger)
log_threshold(INFO)

# Write stdout and stderr to log file
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, type = "output")


data <- read_tsv(snakemake@input[["tsv"]])


# Filtering criteria:
# - P-value < 0.05
# - Depth >= 20
# - For SNP more than 2 reads in each strand
is.deletion <- str_detect(
                        data$ALT,
                        "^[A-Z]",
                        negate = TRUE
                        )
inBothStrands <- data$ALT_RV > 2 & data$ALT_DP > 2
Depth <- (data$ALT_RV + data$ALT_DP) >= 20


log_info("Filtering variants")
data <- filter(
    data,
    as.logical(PASS),
    Depth,
    inBothStrands | is.deletion
    )

log_info("Finding synonymous and non synonymous variants")
# Adding synonymous variable
data <- mutate(
    data,
    synonimous = case_when(
        REF_AA == ALT_AA ~ "yes",
        TRUE             ~  "No"
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