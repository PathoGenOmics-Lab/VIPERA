#!/usr/bin/env Rscript

# Write stdout and stderr to log file
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, type = "output")

library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(logger)

log_threshold(INFO)

# Import file with plots style
source(snakemake@params[["design"]])

# Read inputs
log_info("Reading plot data")
plot.df <- read_delim(snakemake@input$table) %>%
  filter(w != Inf) %>%
  pivot_longer(
    c("dN", "dS", "w"),
    values_to = "value",
    names_to = "ratio"
  )

log_info("Plotting dN and dS")
p.dn.ds <- plot.df %>%
  filter(ratio != "w") %>%
  ggplot() +
  aes(
    x = day,
    y = value,
    color = ratio,
    shape = ratio
  ) +
  geom_point(size = 3) +
  geom_line() +
  scale_color_manual(
    name = "",
    labels = DNDS_LABELS,
    values = DNDS_COLORS
  ) +
  scale_shape_manual(
    name = "",
    values = DNDS_SHAPES,
    labels = DNDS_LABELS
  ) +
  labs(
    y = "Substitution rate",
    x = "Days since the initial sampling",
    color = ""
  )

ggsave(
  filename = snakemake@output$plot_dn_ds,
  plot = p.dn.ds,
  width = snakemake@params$plot_width_mm,
  height = snakemake@params$plot_height_mm,
  units = "mm",
  dpi = 250
)

log_info("Plotting w (dN/dS)")
p.omega <- plot.df %>%
  filter(ratio == "w") %>%
  ggplot() +
  aes(
    x = day,
    y = value,
  ) +
  geom_point(color = "black") +
  geom_line(color = "black") +
  labs(
    y = "w (dN/dS)",
    x = "Days since the initial sampling",
    color = ""
  )

ggsave(
  filename = snakemake@output$plot_omega,
  plot = p.omega,
  width = snakemake@params$plot_width_mm,
  height = snakemake@params$plot_height_mm,
  units = "mm",
  dpi = 250
)
