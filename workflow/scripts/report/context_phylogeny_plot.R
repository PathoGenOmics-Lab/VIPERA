#!/usr/bin/env Rscript

# Write stdout and stderr to log file
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, type = "output")

library(dplyr)
library(readr)
library(jsonlite)
library(ape)
library(ggplot2)
library(ggtree)
library(logger)

log_threshold(INFO)

# Import file with plots style
source(snakemake@params[["design"]])

# Format tree label for well supported nodes
TREE_LEGEND_NAMES["boot_alrt_pass"] <- sprintf(
  TREE_LEGEND_NAMES["boot_alrt_pass"],
  100 * snakemake@params[["boot_th"]], "%",
  100 * snakemake@params[["alrt_th"]], "%"
)

log_info("Reading JSON data")
json <- read_json(snakemake@input$json)

log_info("Reading tree")
tree_ml <- read.tree(snakemake@input[["tree"]]) %>%
  root(
    json$root,
    resolve.root = TRUE
  )

log_info("Reading tree annotation")
annotation <- read_csv(snakemake@input$annotation)

log_info("Plotting M-L tree with context samples")
p <- ggtree(tree_ml, layout = "circular") %<+%
  annotation +
  geom_highlight(node = json$target_mrca, colour = "orange", alpha = 0) +
  geom_point(aes(color = Class, size = Class)) +
  geom_treescale(1.05 * json$max_tip_length) +
  geom_rootedge(0.05 * json$max_tip_length) +
  xlim(-json$max_tip_length / 3, NA) +
  scale_color_manual(
    values = setNames(
      TREE_PALETTE[names(TREE_LEGEND_NAMES)],
      TREE_LEGEND_NAMES
    ),
    na.translate = FALSE
  ) +
  scale_size_manual(
    values = setNames(
      TREE_NODE_SIZE[names(TREE_LEGEND_NAMES)],
      TREE_LEGEND_NAMES
    ),
    na.translate = FALSE
  )

ggsave(
  filename = snakemake@output[["plot"]],
  plot = p,
  width = snakemake@params[["plot_width_mm"]],
  height = snakemake@params[["plot_height_mm"]],
  units = "mm",
  dpi = 250
)
