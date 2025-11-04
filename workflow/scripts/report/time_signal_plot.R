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

log_info("Reading time signal data")
time.signal <- read_csv(snakemake@input$table)

# PLOTS ####
# TempEst
log_info("Plotting time signal data")
p <- time.signal %>%
  ggplot() +
  aes(
    x = date_interval,
    y = distance
  ) +
  geom_smooth(
    method = "lm",
    fill = "gray95",
    alpha = 0.6,
    color = "orange"
  ) +
  geom_point(size = 2, shape = 1) +
  labs(
    y = "Root-to-tip distance",
    x = "Days since the initial sampling"
  )

ggsave(
  filename = snakemake@output$plot,
  plot = p,
  width = snakemake@params[["plot_width_mm"]],
  height = snakemake@params[["plot_height_mm"]],
  units = "mm",
  dpi = 250
)
