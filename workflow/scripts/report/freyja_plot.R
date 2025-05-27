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
demix <- read_csv(snakemake@input[["summary_demixing"]])

# DATA PROCESSING
log_info("Obtaining main lineages")
main_lineages <- demix %>%
  group_by(sample) %>%
  top_n(1, abundances) %>%
  ungroup() %>%
  pull(lineages) %>%
  unique()

# Obtain sample names ordered by CollectionDate
metadata <- read_csv(snakemake@input[["metadata"]])
date_order <- metadata %>%
  arrange(CollectionDate) %>%
  filter(ID %in% demix$sample) %>%
  pull(ID) %>%
  unique()

# PLOT
log_info("Plotting summary demixing")
demix_plot <- demix %>%
  mutate(
    lineages = case_when(
      lineages %in% main_lineages ~ lineages
    )
  ) %>%
  group_by(lineages, sample) %>%
  mutate(
    abundances = sum(abundances)
  ) %>%
  unique() %>%
  ggplot() +
  aes(
    x = factor(sample, date_order),
    y = as.numeric(abundances),
    fill = lineages
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

ggsave(
  filename = snakemake@output[["fig"]],
  plot = demix_plot,
  width = 159.2,
  height = 119.4,
  units = "mm",
  dpi = 250
)


# PLOT TABLES
log_info("Saving plot table")
demix %>%
  mutate(
    lineages = case_when(
      lineages %in% main_lineages ~ lineages,
      TRUE ~ "Other"
    )
  ) %>%
  group_by(sample, lineages) %>%
  summarise(abundances = sum(abundances)) %>%
  ungroup() %>%
  left_join(
    select(
      metadata,
      ID,
      CollectionDate
    ),
    by = c("sample" = "ID")
  ) %>%
  write.csv(snakemake@output[["table"]], row.names = FALSE)
