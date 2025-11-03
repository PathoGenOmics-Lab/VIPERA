#!/usr/bin/env Rscript

# Write stdout and stderr to log file
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, type = "output")

library(tidyverse)
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
  geom_point() +
  geom_text_repel(aes(label = label), max.overlaps = 1000, direction = "x") +
  xlim(c(-1, 1)) +
  geom_hline(
    aes(yintercept = -log10(0.05)),
    linetype = 2,
    color = "red"
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
