#!/usr/bin/env Rscript

# Write stdout and stderr to log file
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, type = "output")

library(tidyverse)
library(logger)
log_threshold(INFO)

# Import file with plots style
source(snakemake@params[["design"]])

# Read inputs
vcf <- read_delim(snakemake@input[["vcf"]])
metadata <- read_delim(snakemake@input[["metadata"]])
N_S_position <- read_delim(snakemake@input[["N_S"]])

# DATA PROCESSING
# Create SNP variable and select useful variables
vcf <- vcf %>%
  dplyr::select(
    SAMPLE,
    VARIANT_NAME,
    REGION,
    ALT_FREQ,
    SYNONYMOUS,
    POS
  )

# Create variable for days sins first sample in metadata
metadata <- metadata %>%
  mutate(
    interval = as.numeric(
      as.Date(CollectionDate) - min(as.Date(CollectionDate))
    )
  ) %>%
  select(ID, interval) %>%
  rename(SAMPLE = ID)

vcf <- left_join(vcf, metadata)

# PLOT
log_info("Ploting dN and dS over time")
plot_df <- vcf %>%
  group_by(SAMPLE, SYNONYMOUS) %>%
  summarise(
    Freq = sum(ALT_FREQ, na.rm = TRUE)
  ) %>%
  pivot_wider(
    names_from = SYNONYMOUS,
    values_from = Freq,
    values_fill = 0
  ) %>%
  transmute(
    dn = No / sum(N_S_position$N),
    ds = Yes / sum(N_S_position$S)
  ) %>%
  ungroup() %>%
  mutate(
    w = dn / ds,
  ) %>%
  filter(w != Inf) %>%
  pivot_longer(
    c("dn", "ds", "w"),
    values_to = "value",
    names_to = "d"
  ) %>%
  left_join(unique(select(vcf, SAMPLE, interval)))

plot <- plot_df %>%
  filter(d != "w") %>%
  ggplot() +
  aes(
    x = interval,
    y = value,
    color = d,
    shape = d
  ) +
  geom_point(size = 3) +
  geom_line() +
  scale_color_manual(
    name = "Parameter",
    labels = DNDS_LABELS,
    values = DNDS_COLORS
  ) +
  scale_shape_manual(
    name = "Parameter",
    values = DNDS_SHAPES,
    labels = DNDS_LABELS
  ) +
  labs(
    y = "Substitution rate",
    x = "Days since the initial sampling",
    color = "Parameter"
  )

ggsave(
  filename = snakemake@output[["plot"]],
  plot = plot,
  width = 159.2,
  height = 119.4,
  units = "mm",
  dpi = 250
)

# Plot for omega

plot_omega <- plot_df %>%
  filter(d == "w") %>%
  ggplot() +
  aes(
    x = interval,
    y = value,
  ) +
  geom_point(color = "black") +
  geom_line(color = "black") +
  labs(
    y = "w (dN/dS)",
    x = "Days since the initial sampling",
    color = "Parameter"
  )

ggsave(
  filename = snakemake@output[["plot_omega"]],
  plot = plot_omega,
  width = 159.2,
  height = 119.4,
  units = "mm",
  dpi = 250
)

# PLOT TABLES
log_info("Saving plot table")
vcf %>%
  group_by(SAMPLE, SYNONYMOUS) %>%
  summarise(
    Freq = sum(ALT_FREQ, na.rm = TRUE)
  ) %>%
  pivot_wider(
    names_from = SYNONYMOUS,
    values_from = Freq,
    values_fill = 0
  ) %>%
  transmute(
    dn = No / sum(N_S_position$N),
    ds = Yes / sum(N_S_position$S)
  ) %>%
  ungroup() %>%
  left_join(unique(select(vcf, SAMPLE, interval))) %>%
  transmute(
    sample = SAMPLE,
    DaysSinceFirst = interval,
    dN = dn,
    dS = ds,
    w = dn / ds
  ) %>%
  write.csv(snakemake@output[["table"]], row.names = FALSE)
