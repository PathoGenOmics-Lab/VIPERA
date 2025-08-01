#!/usr/bin/env Rscript

# Write stdout and stderr to log file
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, type = "output")

library(tidyverse)
library(jsonlite)
library(logger)
log_threshold(INFO)

# Import file with plots style
source(snakemake@params[["design"]])

# Read data
log_info("Reading diversity results")
divs <- read_lines(snakemake@input$divs) %>%
  as.numeric()
json <- read_json(snakemake@input$json)

# Plot and save
log_info("Plotting")
p <- data.frame(pi = divs) %>%
  ggplot() +
  geom_density(
    aes(x = pi),
    fill = DIVERSITY_PALETTE["density_fill"],
    alpha = 0.7,
    bw = 0.000001,
    color = DIVERSITY_PALETTE["density_color"]
  ) +
  geom_vline(
    xintercept = json$diversity,
    color = DIVERSITY_PALETTE["value_color"]
  ) +
  stat_function(
    fun = dnorm,
    args = list(mean = mean(divs), sd = sd(divs)),
    color = DIVERSITY_PALETTE["dnorm_color"]
  ) +
  labs(
    x = "π",
    y = "Density"
  )

ggsave(
  filename = snakemake@output[["plot"]],
  plot = p,
  width = snakemake@params[["plot_width_mm"]],
  height = snakemake@params[["plot_height_mm"]],
  units = "mm",
  dpi = 250
)
