#!/usr/bin/env Rscript

# Write stdout and stderr to log file
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, type = "output")

library(readr)
library(dplyr)
library(tibble)
library(jsonlite)
library(ape)
library(adephylo)
library(logger)

log_threshold(INFO)

# Read distance tree and root
log_info("Reading tree")
tree <- read.tree(snakemake@input$tree) %>%
  root(snakemake@params$outgroup_id, resolve.root = TRUE)

# Read metadata
log_info("Reading metadata")
metadata <- read_csv(snakemake@input$metadata)

# Get patristic distances to ancestor in distance tree
log_info("Calculating patristic distances to ancestor in distance tree")
time.signal <- distRoot(
  tree,
  "all",
  method = "patristic"
) %>%
  as.data.frame() %>%
  rownames_to_column(var = "ID") %>%
  filter(ID != snakemake@params$outgroup_id) %>%
  rename(distance = ".") %>%
  left_join(
    select(
      metadata,
      ID,
      CollectionDate
    ),
    by = "ID"
  ) %>%
  mutate(
    date_interval = as.numeric(
      as.Date(CollectionDate) - min(as.Date(CollectionDate))
    )
  )

# Save table
log_info("Saving time signal data")
log_debug("Saving table")
write_csv(time.signal, snakemake@output$table)

# Time signal stats
log_debug("Building linear model")
model <- lm(distance ~ date_interval, data = time.signal)
p.value <- summary(model)$coefficients[2, 4]

# TREE STATS
log_debug("Saving linear model")
list(
  "sub_rate" = model$coefficients[[2]] * 365,
  "r2" = summary(model)$r.squared[[1]],
  "pvalue" = ifelse(p.value < 0.001, "< 0.001", p.value)
) %>%
  toJSON() %>%
  write(snakemake@output$json)
