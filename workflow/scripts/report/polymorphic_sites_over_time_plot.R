#!/usr/bin/env Rscript

# Write stdout and stderr to log file
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, type = "output")

library(readr)
library(dplyr)
library(ggplot2)
library(logger)

log_threshold(INFO)

# Import file with plots style
source(snakemake@params[["design"]])

log_info("Reading plot data")
df <- read_csv(snakemake@input$table)

log_info("Plotting")
p <- df %>%
  ggplot() +
  aes(x = Day, y = n) +
  geom_smooth(
    method = "lm",
    fill = "gray95",
    alpha = 0.6,
    colour = "orange"
  ) +
  geom_point(size = 2, shape = 1) +
  labs(
    x = "Days since the initial sampling",
    y = "No. of polimorphic sites"
  )

log_info("Saving plot")
ggsave(
  filename = snakemake@output$plot,
  plot = p,
  width = snakemake@params$plot_width_mm,
  height = snakemake@params$plot_height_mm,
  units = "mm",
  dpi = 250
)
