#!/usr/bin/env Rscript

# Write stdout and stderr to log file
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, type = "output")

library(tidyverse)
library(logger)
log_threshold(INFO)

# Read data
log_info("Reading variants table")
variants <- read_delim(
  snakemake@input[["variants"]],
  col_select = c(
    "VARIANT_NAME",
    "SAMPLE",
    "REGION",
    "ALT_FREQ",
    "POS"
  )
) %>%
  # Fill positions without alt frequency with 0
  complete(
    nesting(REGION, VARIANT_NAME, POS),
    SAMPLE,
    fill = list(ALT_FREQ = 0)
  )

log_info("Reading metadata")
metadata <- read_csv(snakemake@input[["metadata"]]) %>%
  filter(
    ID %in% variants$SAMPLE
  ) %>%
  select(
    ID,
    CollectionDate
  )

log_info("Dating variants")
variants <- left_join(variants, metadata, by = c("SAMPLE" = "ID")) %>%
  # Calculate days since first sample
  arrange(
    CollectionDate
  ) %>%
  mutate(
    interval = as.numeric(CollectionDate - min(CollectionDate))
  )

# Save processed input
log_info("Writing dated and frequency-filled variants")
write_csv(variants, snakemake@output$fmt_variants)

log_info("Calculating correlations")
log_debug("Calculating unique SNPs")
# Get list with all different polymorphisms
unique.snps <- pull(
  variants,
  VARIANT_NAME
) %>%
  unique()

# Create an empty dataframe to be filled
cor.df <- data.frame(
  variant = "",
  coefficient = 0,
  p.value = 0,
  p.value.adj = 0
) %>%
  filter(p.value != 0)

log_debug("Calculating correlation using method = {snakemake@params$cor_method} and exact p-value = {snakemake@params$cor_exact}")
correlations <- lapply(
  unique.snps,
  function(snp) {
    # Select SNP
    df <- filter(
      variants,
      VARIANT_NAME == snp
    )
    # Perform calculation
    test <- cor.test(
      df$ALT_FREQ,
      df$interval,
      method = snakemake@params$cor_method,
      exact = snakemake@params$cor_exact
    )
    # Adjust p-value
    p.value.adj <- p.adjust(
      test$p.value,
      method = "BH",
      n = length(unique.snps)
    )
    # Add row to dataframe
    add_row(
      cor.df,
      variant = snp,
      coefficient = test$estimate,
      p.value = test$p.value,
      p.value.adj = p.value.adj
    )
  }
) %>%
  bind_rows()

log_info("Writing correlations table")
write_csv(
  correlations,
  snakemake@output[["correlations"]]
)
