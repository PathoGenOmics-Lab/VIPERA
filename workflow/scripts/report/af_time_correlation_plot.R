#!/usr/bin/env Rscript

# Write stdout and stderr to log file
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, type = "output")

library(dplyr)
library(readr)
library(ggplot2)
library(ggrepel)
library(logger)
log_threshold(INFO)

# Import file with plots style
source(snakemake@params[["design"]])

log_info("Reading correlation data")
correlations <- read_csv(snakemake@input$correlations)

log_info("Plotting coorrelation coefficients and p-values of each variants")
p <- correlations %>%
  mutate(
    trans.p = -log10(p.value.adj),
    label = ifelse(p.value.adj < 0.05, variant, NA)
  ) %>%
  ggplot() +
  aes(
    x = coefficient,
    y = trans.p
  ) +
  geom_text_repel(aes(label = label), max.overlaps = 10000, direction = "x") +
  geom_point(
    data = function(x) subset(x, !is.na(label)),
    color = "orange",
    size = 2
  ) +
  geom_point(size = 2, shape = 1) +
  xlim(c(-1, 1)) +
  geom_hline(
    aes(yintercept = -log10(0.05)),
    linetype = 2,
    color = "orange"
  ) +
  labs(
    x = "Correlation coefficient",
    y = "-log10(p-value)"
  )

ggsave(
  filename = snakemake@output$plot,
  plot = p,
  width = snakemake@params$plot_width_mm,
  height = snakemake@params$plot_height_mm,
  units = "mm",
  dpi = 250
)
