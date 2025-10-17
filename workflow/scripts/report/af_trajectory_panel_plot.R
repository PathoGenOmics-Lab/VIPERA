#!/usr/bin/env Rscript

# Write stdout and stderr to log file
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, type = "output")

set.seed(snakemake@params$random_color_seed)  # seed for sampling colors

library(tidyverse)
library(logger)
log_threshold(INFO)

# Import file with plots style
source(snakemake@params[["design"]])

log_info("Reading formatted variants table")
variants <- read_csv(snakemake@input$fmt_variants)

log_info("Reading subset of variants to represent")
selected.variants <- read_lines(snakemake@input$subset)

# Set plot height depending on the number of SNPs assuming 4 columns in the plot
log_debug("Calculating plot height")
plot.height <- ceiling(length(selected.variants) / 4) * 42

log_info("Plotting {length(selected.variants)} SNPs allele frequency trajectories in time")
selected.colors <- sample(TRAJECTORY.PANEL.COLORS, length(selected.variants))
log_debug("Selected color: {selected.colors}")
p <- variants %>%
  filter(VARIANT_NAME %in% selected.variants) %>%
  ggplot() +
  aes(
    x = interval,
    y = ALT_FREQ,
    color = VARIANT_NAME
  ) +
  scale_color_manual(values = selected.colors) +
  geom_point() +
  geom_line() +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 9)
  ) +
  labs(
    x = "Days since first sample",
    y = "Frequency",
    color = "NV"
  ) +
  guides(color = guide_legend(ncol = 3))

if (length(selected.variants) > 1) {
  p <- p +
    facet_wrap(
      vars(POS),
      nrow = ceiling(length(selected.variants) / 4),
      ncol = 4
    )
}

ggsave(
  filename = snakemake@output[["plot"]],
  plot = p,
  width = 159.2,
  height = max(100, plot.height),
  units = "mm",
  dpi = 250
)
