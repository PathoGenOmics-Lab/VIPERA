#!/usr/bin/env Rscript

library(tidyverse)
library(logger)
log_threshold(INFO)

# Write stdout and stderr to log file
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, type = "output")


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
    variant,
    REGION,
    ALT_FREQ,
    GFF_FEATURE,
    synonimous,
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
  rename(REGION = ID)

vcf <- left_join(vcf, metadata)

# PLOT
log_info("Ploting dN and dS over time")
plot_df <- vcf %>%
  group_by(REGION, synonimous) %>%
  summarise(
    Freq = sum(ALT_FREQ, na.rm = TRUE)
  ) %>%
  pivot_wider(
    names_from = synonimous,
    values_from = Freq,
    values_fill = 0
  )  %>%
  transmute(
    dn = No / sum(N_S_position$N),
    ds = Yes / sum(N_S_position$S)
  ) %>%
  ungroup() %>%
  mutate(
    w_raw = dn / ds,
    w = w_raw * (mean(c(dn, ds)) / mean(w_raw[w_raw != Inf])),
    mean_value = mean(c(dn, ds))
    ) %>%
  filter(w_raw != Inf) %>%
  pivot_longer(
    c("dn", "ds", "w"),
    values_to = "value",
    names_to = "d"
  ) %>%
  left_join(unique(select(vcf, REGION, interval)))

  plot <- plot_df %>%
  ggplot() +
  aes(
    x = interval,
    y = value,
    color = d,
    linetype = d
  ) +
  geom_point() +
  geom_line() +
  scale_color_manual(
    labels = dnds.labels,
    values = dnds.colors
    ) +
  scale_linetype_manual(
    values = c(dn = 2, ds = 2, w = 1),
    guide = "none"
    ) +
  scale_y_continuous(
    sec.axis = sec_axis(
      ~ . * (mean(plot_df$w_raw)) / as.numeric(plot_df$mean_value[1]),
      name = "w (dN/dS)")) +
  labs(
    y = "",
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


# PLOT TABLES
log_info("Saving plot table")
vcf %>%
  group_by(REGION, synonimous) %>%
  summarise(
    Freq = sum(ALT_FREQ, na.rm = TRUE)
  ) %>%
  pivot_wider(
    names_from = synonimous,
    values_from = Freq,
    values_fill = 0
  ) %>%
  transmute(
    dn = No / sum(N_S_position$N),
    ds = Yes / sum(N_S_position$S)
  ) %>%
  ungroup() %>%
  left_join(unique(select(vcf, REGION, interval))) %>%
  transmute(
    sample = REGION,
    DaysSinceFirst = interval,
    dN = dn,
    dS = ds,
    w = dn / ds
  ) %>%
  write.csv(snakemake@output[["table"]], row.names = FALSE)
