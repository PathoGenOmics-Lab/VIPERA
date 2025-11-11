#!/usr/bin/env Rscript

# Write stdout and stderr to log file
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, type = "output")

library(dplyr)
library(readr)
library(tibble)
library(jsonlite)
library(ape)
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

log_info("Reading tree")
tree_ml <- read.tree(snakemake@input[["tree"]]) %>%
  root(
    snakemake@params[["ref_name"]],
    resolve.root = TRUE
  )
log_debug("Read tree with {length(tree_ml$node.label)} labels")

log_info("Reading target names from FASTA")
target_names <- read.dna(
  snakemake@input[["target_fasta"]],
  format = "fasta",
  as.matrix = FALSE,
)
log_debug("Read {length(target_names)} records")

log_info("Processing target names")
target_names <- target_names[
  !startsWith(names(target_names), snakemake@config[["ALIGNMENT_REFERENCE"]])
]
target_names <- names(target_names)
log_debug("{length(target_names)} records remaining after processing")

# ML tree with context data
# Internal nodes color
# Node labels contain SH-aLRT/UFboot values
log_info("Reading support values from labels")
labels <- strsplit(tree_ml$node.label, "/")
aLRT.values <- sapply(
  labels,
  function(x) as.numeric(x[1])
)
bootstrap.values <- sapply(
  labels,
  function(x) as.numeric(x[2])
)

log_info("Calculating support mask for the given thresholds")
aLRT.mask <- aLRT.values >= snakemake@params[["alrt_th"]]
boot.mask <- bootstrap.values >= snakemake@params[["boot_th"]]

# MRCA for target samples
log_info("Calculating MRCA of target samples")
target.mrca <- getMRCA(tree_ml, target_names)
target.mrca.clade <- extract.clade(tree_ml, target.mrca)
target.mrca.clade.ntips <- Ntip(target.mrca.clade)

log_info("Building M-L tree annotation")
tree_ml.ntips <- length(tree_ml$tip.label)
tree_ml.nodes <- tree_ml$Nnode + tree_ml.ntips
tree_ml.labels <- tree_ml$tip.label[1:tree_ml.nodes]
tree_ml.node.pass <- c(rep(FALSE, tree_ml.ntips), aLRT.mask & boot.mask)

ml.tree.annot <- tibble(
  node = 1:tree_ml.nodes,
) %>%
  mutate(
    Class = case_when(
      tree_ml.labels %in% target_names ~ TREE_LEGEND_NAMES["tip_label"],
      tree_ml.node.pass ~ TREE_LEGEND_NAMES["boot_alrt_pass"],
      TRUE ~ NA
    )
  )

# Write output files
log_info("Writing tree annotation")
write_csv(ml.tree.annot, snakemake@output$annotation)

log_info("Writing JSON data")
target.node <- tree_ml$node.label[target.mrca - length(tree_ml$tip.label)]
list(
  "boot" = strsplit(target.node, "/")[[1]][2] %>% as.numeric(),
  "alrt" = strsplit(target.node, "/")[[1]][1] %>% as.numeric(),
  "monophyly" = ifelse(
    is.monophyletic(tree_ml, target_names),
    "are",
    "are not"
  ),
  "target_mrca" = target.mrca,
  "clade_tips" = target.mrca.clade.ntips,
  "max_tip_length" = max(node.depth.edgelength(tree_ml)[
    1:length(tree_ml$tip.label)
  ]),
  "root" = snakemake@params[["ref_name"]]
) %>%
  write_json(
    snakemake@output$json,
    auto_unbox = TRUE,
    digits = NA
  )
