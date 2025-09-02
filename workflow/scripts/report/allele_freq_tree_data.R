#!/usr/bin/env Rscript

# Write stdout and stderr to log file
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, type = "output")

library(tidyverse)
library(data.table)
library(ape)
library(logger)

fix_negative_edge_length <- function(nj.tree) {
  edge_infos <- cbind(
    nj.tree$edge,
    nj.tree$edge.length
  ) %>%
    as.data.table()
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

log_threshold(INFO)

matrix <- read_csv(snakemake@input$dist)

# Tree construction
if (snakemake@params$use_bionj) {
  log_info("Constructing BIONJ tree based on distances")
  tree_method <- bionj
} else {
  log_info("Constructing NJ tree based on distances")
  tree_method <- nj
}
tree <- matrix %>%
  column_to_rownames(var = "...1") %>%
  as.matrix() %>%
  as.dist() %>%
  tree_method() %>%
  root(snakemake@params$ref_name, resolve.root = TRUE)

# Resolve possible negative edge lengths in tree
log_info("Resolving negative edge lengths")
tree <- fix_negative_edge_length(tree)

log_info("Saving tree")
write.tree(tree, snakemake@output$tree)
