#!/usr/bin/env Rscript

# Write stdout and stderr to log file
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, type = "output")

library(dplyr)
library(readr)
library(tidyr)
library(tibble)
library(logger)
log_threshold(INFO)

# Read data
log_info("Reading variants table")
variants <- read_delim(
  snakemake@input[["variants"]],
  col_select = c(
    "VARIANT_NAME",
    "SAMPLE",
    "CHROM",
    "ALT_FREQ",
    "POS"
  )
) %>%
  # Fill positions without alt frequency with NA
  complete(
    nesting(CHROM, VARIANT_NAME, POS),
    SAMPLE,
    fill = list(ALT_FREQ = NA)
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
    interval = as.numeric(difftime(CollectionDate, min(CollectionDate), units = "days"))
  )

# Save processed input
log_info("Writing dated and frequency-filled variants")
write_csv(variants, snakemake@output$fmt_variants)

log_info("Calculating correlations")
log_debug("Calculating unique SNPs")
# Get list with all different polymorphisms
unique.snps <- unique(variants$VARIANT_NAME)

# Create an empty dataframe to be filled
cor.df <- data.frame(
  variant = "",
  min_af = 0,
  max_af = 0,
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
      min_af = min(df$ALT_FREQ, na.rm = TRUE),
      max_af = max(df$ALT_FREQ, na.rm = TRUE),
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

log_info("Selecting variants whose allele frequency is significantly correlated with time")
significant.variants <- correlations %>%
  filter(
    p.value.adj <= snakemake@params$max_p_adj_threshold,
    abs(coefficient) >= snakemake@params$min_abs_cor_threshold,
    (max_af - min_af) >= snakemake@params$min_diff_af_threshold
  ) %>%
  pull(variant) %>%
  unique()

log_info("Significant: {significant.variants}")

log_info("Selecting variants in positions with more than one alternative allele")
mult.alt.variants <- variants %>%
  select(
    VARIANT_NAME,
    POS
  ) %>%
  distinct() %>%
  group_by(POS) %>%
  filter(n() > 1) %>%
  ungroup() %>%
  pull(VARIANT_NAME) %>%
  unique()

log_info("Mult all: {mult.alt.variants}")

# Build selected subset to represent
variant.selection <- unique(c(significant.variants, mult.alt.variants))

log_info("Selection: {variant.selection}")

log_info("Writing selected variants subset")
write_lines(variant.selection, snakemake@output$subset)
