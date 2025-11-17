#!/usr/bin/env Rscript

# Write stdout and stderr to log file
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, type = "output")

set.seed(snakemake@params$random_color_seed)  # seed for sampling colors

library(dplyr)
library(readr)
library(ggplot2)
library(glue)
library(logger)
library(logger)
log_threshold(INFO)

# Import file with plots style
source(snakemake@params[["design"]])

log_info("Reading formatted variants table")
variants <- read_csv(snakemake@input$fmt_variants)

log_info("Reading subset of variants to represent")
selected.variants <- read_lines(snakemake@input$subset)

# Set plot height depending on the number of SNPs assuming 4 columns in the plot
plot.rows <- ceiling(
  length(selected.variants) / snakemake@params$n_plot_columns
)
plot.height <- max(100, plot.rows * snakemake@params$plot_row_height_mm)
log_debug("Setting total plot height to {plot.height} mm with {plot.rows} rows")

log_info("Plotting {length(selected.variants)} SNPs allele frequency trajectories in time")
selected.colors <- sample(TRAJECTORY.PANEL.COLORS, length(selected.variants))
log_debug("Selected color: {selected.colors}")
p <- variants %>%
  filter(VARIANT_NAME %in% selected.variants) %>%
  mutate(
    gPOS = paste0("g.", POS),
    gPOS = reorder(gPOS, POS)
  ) %>%
  ggplot() +
  aes(
    x = interval,
    y = ALT_FREQ,
    color = reorder(glue("{gPOS}\n{VARIANT_NAME}"), POS)
  ) +
  scale_color_manual(values = selected.colors) +
  geom_point() +
  geom_line() +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 9),
    legend.title = element_blank(),
    legend.spacing.y = unit(3, "mm")
  ) +
  labs(
    x = "Days since first sample",
    y = "Frequency"
  ) +
  guides(color = guide_legend(ncol = 3))

if (length(selected.variants) > 1) {
  p <- p +
    facet_wrap(
      vars(gPOS),
      ncol = snakemake@params$n_plot_columns
    )
}

ggsave(
  filename = snakemake@output[["plot"]],
  plot = p,
  width = snakemake@params$plot_width_mm,
  height = plot.height,
  units = "mm",
  dpi = 250
)
