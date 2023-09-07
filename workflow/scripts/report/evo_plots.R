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


vcf <- read_delim(snakemake@input[["vcf"]])
metadata <- read_delim(snakemake@params[["metadata"]])
N_S_position <- read_delim(snakemake@input[["N_S"]])

# DATA PROCESSING

# Create SNP variable and select useful variables
vcf <- vcf %>%
  mutate(
    SNP = paste(REF, POS, ALT, sep = "-")) %>%
  dplyr::select(
    SNP,
    REGION,
    ALT_FREQ,
    GFF_FEATURE,
    synonimous
  ) %>%
  rowwise() %>%
  mutate(POS = strsplit(SNP, "-")[[1]][2]) %>%
  ungroup()

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
plot <- vcf %>%
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
    dn = No / sum(N_S_position$S),
    ds = yes / sum(N_S_position$S)
  ) %>%
  ungroup() %>%
  pivot_longer(
    c("dn", "ds"),
    values_to = "value",
    names_to = "d"
  ) %>%
  left_join(unique(select(vcf, REGION, interval))) %>%
  ggplot() +
  aes(
    x = interval,
    y = value,
    color = d
  ) +
  geom_point() +
  geom_line() +
  scale_color_hue(labels = c("dN", "dS")) +
  labs(
    y = "",
    x = "Time since first sample",
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
    dn = No / sum(N_S_position$S),
    ds = yes / sum(N_S_position$S)
  ) %>%
  ungroup() %>%
  left_join(unique(select(vcf, REGION, interval))) %>%
  transmute(
    sample = REGION,
    DaysSinceFirst = interval,
    dN = dn,
    dS = ds
  ) %>%
  write.csv(snakemake@output[["table"]], row.names = FALSE)
