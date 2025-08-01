#!/usr/bin/env Rscript

# Write stdout and stderr to log file
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, type = "output")

library(tidyverse)
library(logger)
log_threshold(INFO)

# Import file with plots style
source(snakemake@params[["design"]])

# Read inputs
log_info("Plotting")
p <- read_csv(snakemake@input$data) %>%
  mutate(
    sample = reorder(sample, CollectionDate),
    lineage = ifelse(is_main, lineage, NA)
  ) %>%
  ggplot() +
  aes(
    x = sample,
    y = abundance,
    fill = lineage
  ) +
  scale_fill_viridis_d(
    na.value = "gray50",
    labels = function(x) {
      ifelse(is.na(x), "Other", x)
    }
  ) +
  geom_col() +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1),
    legend.position = "bottom"
  ) +
  labs(
    x = "Sample",
    y = "Relative abundance",
    fill = "Lineage"
  )

log_info("Saving plot")
ggsave(
  filename = snakemake@output[["plot"]],
  plot = p,
  width = snakemake@params$plot_width_mm,
  height = snakemake@params$plot_height_mm,
  units = "mm",
  dpi = 250
)
