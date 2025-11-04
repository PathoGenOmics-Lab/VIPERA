#!/usr/bin/env Rscript

# Write stdout and stderr to log file
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, type = "output")

library(dplyr)
library(readr)
library(ggplot2)
library(ggpubr)
library(jsonlite)
library(tibble)
library(purrr)
library(logger)

log_threshold(INFO)

# Import file with plots style
source(snakemake@params[["design"]])

# Subplot A: windows ============================
log_info("Reading variant windows")
window <- read_csv(snakemake@input$window)

log_debug("Calculating X axis limits")
xlim_values <- c(
  max(0, min(window$position) - snakemake@params$window_step),
  max(window$position) + snakemake@params$window_step
)

log_info("Plotting variant windows")
p1 <- window %>%
  ggplot() +
  aes(
    x = position,
    y = fraction,
    color = feature
  ) +
  geom_point() +
  geom_line(
    aes(group = 1),
    colour = "black",
    alpha = 0.3
  ) +
  scale_y_continuous(
    label = scales::percent,
    limits = c(0, max(window$fraction) + mean(window$fraction))
  ) +
  xlim(xlim_values) +
  # TODO: change GENE_PALETTE to selection of TRAJECTORY.PANEL.COLORS ?
  scale_color_manual(values = GENE_PALETTE) +
  labs(
    y = "Proportion of\nsites with NV",
    x = "",
    color = "Gen"
  )

log_info("Reading regions to highlight")
highlight <- read_csv(snakemake@input$highlight_window_regions)

if (nrow(highlight) != 0) {
  log_info("Highlighting {nrow(highlight)} regions")
  highlighted_sites <- unique(c(highlight$start, highlight$end))
  p1 <- p1 +
    geom_vline(xintercept = highlighted_sites, color = "orange") +
    geom_text(
      data = highlight,
      aes(
        x = (start + end) / 2,
        y = max(window$fraction + mean(window$fraction) / 2),
        label = region
      ),
      inherit.aes = FALSE,
      size = 3,
      angle = 45
    )
} else {
  log_info("No highlighted regions")
}

# Subplot B: windows ============================
log_info("Reading panel data")
panel <- read_csv(snakemake@input$panel)

log_info("Calculating sample chronological order")
date_order <- panel %>%
  arrange(CollectionDate) %>%
  pull(SAMPLE) %>%
  unique()

log_info("Reading regions")
regions <- read_json(snakemake@input$regions, simplifyVector = TRUE) %>%
  enframe(name = "region", value = "coords") %>%
  mutate(
    start = map_dbl(coords, 1),
    end = map_dbl(coords, 2),
    length = end - start + 1
  ) %>%
  select(-coords)

log_info("Plotting variants")
p2 <- panel %>%
  ggplot() +
  aes(
    x = POS,
    y = factor(SAMPLE, date_order),
    shape = factor(NV_class, c("SNP", "INDEL")),
    color = group,
    alpha = ALT_FREQ
  ) +
  geom_point(size = 3) +
  geom_col(
    data = regions,
    aes(
      x = length,
      y = 0.3,
      # TODO: legend shows as NA those without match in `regions`
      fill = factor(region, rev(regions$region))
    ),
    inherit.aes = FALSE,
    width = 0.3
  ) +
  # TODO: change GENE_PALETTE to selection of TRAJECTORY.PANEL.COLORS ?
  scale_fill_manual(values = GENE_PALETTE) +
  xlim(xlim_values) +
  scale_color_manual(
    labels = NV_TYPE_NAMES,
    values = NV_TYPE_PALETTE
  ) +
  scale_y_discrete(drop = FALSE) +
  labs(
    x = "Genome position",
    y = "Sample",
    shape = "Variant class",
    color = "Classification",
    alpha = "Frequency",
    fill = "Region"
  ) +
  guides(
    fill = guide_legend(reverse = TRUE)
  )

# Plot arrangement ============================
log_info("Arranging plots")
p <- ggarrange(
  p1,
  p2,
  nrow = 2,
  align = "v",
  legend.grob = get_legend(p2),
  heights = c(2, 6),
  legend = "right",
  labels = c("A", "B")
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
