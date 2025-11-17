#!/usr/bin/env Rscript

# Write stdout and stderr to log file
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, type = "output")

library(dplyr)
library(readr)
library(ape)
library(ggplot2)
library(ggtree)
library(logger)

log_threshold(INFO)

# Import file with plots style
source(snakemake@params[["design"]])

log_info("Building tree annotation")

# Read metadata
log_debug("Reading metadata")
metadata <- read_csv(snakemake@input$metadata)

# Extract names from records
log_debug("Reading FASTA")
study_records <- read.dna(
  snakemake@input$study_fasta,
  format = "fasta",
  as.matrix = FALSE,
)
log_debug("Extracting names")
study_records <- study_records[
  !startsWith(names(study_records), snakemake@params$outgroup_id)
]
study_names <- names(study_records)

# Obtain sample names ordered by CollectionDate
log_debug("Sorting names by collection date")
date_order <- metadata %>%
  arrange(CollectionDate) %>%
  filter(ID %in% study_names) %>%
  pull(ID) %>%
  unique()

log_debug("Building annotation dataframe")
tree_tiplab <- data.frame(
  ID = date_order,
  order = seq(1, length(date_order), 1)
) %>%
  rowwise() %>%
  mutate(
    tip_label = sprintf(
      "(%s)-%s",
      order,
      ID
    )
  ) %>%
  ungroup() %>%
  add_row(
    ID = snakemake@params$outgroup_id,
    order = 0,
    tip_label = snakemake@params$outgroup_id
  )

# Read distance tree and root
log_info("Reading tree")
tree <- read.tree(snakemake@input$tree) %>%
  root(snakemake@params$outgroup_id, resolve.root = TRUE)

# Plot
log_info("Plotting distance tree")
max.tip.length <- max(
  node.depth.edgelength(tree)[seq_along(length(tree$tip.label))]
)
p <- ggtree(tree) %<+% tree_tiplab +
  geom_tiplab(aes(label = tip_label)) +
  geom_treescale(1.1 * max.tip.length) +
  geom_rootedge(0.05 * max.tip.length)

ggsave(
  filename = snakemake@output$plot,
  plot = p,
  width = snakemake@params$plot_width_mm,
  height = snakemake@params$plot_height_mm,
  units = "mm",
  dpi = 250
)
