#!/usr/bin/env Rscript

# Write stdout and stderr to log file
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, type = "output")

library(tidyverse)
library(stringi)
library(ape)
library(ggtree)
library(data.table)
library(ggpubr)
library(pegas)
library(jsonlite)
library(logger)

log_threshold(INFO)

# Import file with plots style
source(snakemake@params[["design"]])

matrix <- read_csv(snakemake@input[["dist"]])
metadata <- read_csv(snakemake@input[["metadata"]])

study_names <- read.dna(
  snakemake@input[["study_fasta"]],
  format = "fasta",
  as.matrix = FALSE,
)

study_names <- study_names[
  !startsWith(names(study_names), snakemake@config[["ALIGNMENT_REFERENCE"]])
]
study_names <- names(study_names)

# Obtain sample names ordered by CollectionDate
date_order <- metadata %>%
  arrange(CollectionDate) %>%
  filter(ID %in% study_names) %>%
  pull(ID) %>%
  unique()

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
    ID = snakemake@params[["ref_name"]],
    order = 0,
    tip_label = snakemake@params[["ref_name"]]
  )


# DATA PROCESSING ####
# Resolve possible negative edge lengths in n-j tree
fix_negative_edge_length <- function(nj.tree) {
  edge_infos <- cbind(
    nj.tree$edge,
    nj.tree$edge.length
  ) %>%
    as.data.table
  colnames(edge_infos) <- c("from", "to", "length")

  nega_froms <- edge_infos[length < 0, sort(unique(from))]

  for (nega_from in nega_froms) {
    minus_length <- edge_infos[from == nega_from, ][order(length)][1, length]
    edge_infos[from == nega_from, length := length - minus_length]
    edge_infos[to == nega_from, length := length + minus_length]
  }
  nj.tree$edge.length <- edge_infos$length
  nj.tree
}

# Tree construction
if (snakemake@params$use_bionj) {
  log_info("Constructing BIONJ tree based on distances")
  tree_method = bionj
} else {
  log_info("Constructing NJ tree based on distances")
  tree_method = nj
}
tree <- matrix %>%
  column_to_rownames(var = "...1") %>%
  as.matrix() %>%
  as.dist() %>%
  tree_method() %>%
  root(snakemake@params[["ref_name"]], resolve.root = TRUE)

tree <- fix_negative_edge_length(tree)

# Get patristic distances to ancestor from n-j tree
log_info("Getting patristic distances to ancestor from n-j tree")
tempest <- adephylo::distRoot(
  tree,
  "all",
  method = "patristic"
) %>%
  as.data.frame() %>%
  rownames_to_column(var = "ID") %>%
  filter(ID != snakemake@params[["ref_name"]]) %>%
  rename(distance = ".") %>%
  left_join(
    select(
      metadata,
      ID,
      CollectionDate
    )
  ) %>%
  mutate(
    date_interval = as.numeric(
      as.Date(CollectionDate) - min(as.Date(CollectionDate))
    )
  )


# PLOTS ####
# (BIO)NJ tree
log_info("Plotting distance tree")
tree_plot <- ggtree(tree) %<+%
  tree_tiplab +
  geom_tiplab(aes(label = tip_label)) +
  geom_treescale() +
  geom_rootedge(0.5) +
  xlim(0, max(tempest$distance) + 3.5)

ggsave(
  filename = snakemake@output[["tree"]],
  plot = tree_plot,
  width = snakemake@params[["plot_width_mm"]],
  height = snakemake@params[["plot_height_mm"]],
  units = "mm",
  dpi = 250
)

# TempEst
log_info("Plotting temporal signal analysis")
tempest_fig <- tempest %>%
  ggplot() +
  aes(
    x = date_interval,
    y = distance
  ) +
  geom_smooth(
    method = "lm",
    fill = "gray95",
    alpha = 0.6,
    color = "red"
  ) +
  geom_point() +
  labs(
    y = "Root-to-tip distance",
    x = "Days since the initial sampling"
  )

ggsave(
  filename = snakemake@output[["temest"]],
  plot = tempest_fig,
  width = snakemake@params[["plot_width_mm"]],
  height = snakemake@params[["plot_height_mm"]],
  units = "mm",
  dpi = 250
)

# PLOT TABLES
log_info("Saving temporal signal table")
tempest %>%
  transmute(
    sample = ID,
    Distance = distance,
    CollectionDate = CollectionDate,
    DaysSinceFirst = date_interval
  ) %>%
  write.csv(snakemake@output[["table"]], row.names = FALSE)

# TEMPEST STATS
model <- lm(distance ~ date_interval, data = tempest)
p.value <- summary(model)$coefficients[2, 4]

# TREE STATS
list(
  "sub_rate" = model$coefficients[[2]] * 365,
  "r2" = summary(model)$r.squared[[1]],
  "pvalue" = ifelse(p.value < 0.001, "< 0.001", p.value)
) %>%
  toJSON() %>%
  write(snakemake@output[["json"]])
