#!/usr/bin/env Rscript

library(tidyverse)
library(logger)
log_threshold(INFO)

# Write stdout and stderr to log file
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, type = "output")

# Import file with plots style
source(snakemake@params[["design"]])


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
date_order <- read_csv(snakemake@params[["metadata"]]) %>%
  arrange(CollectionDate) %>%
  filter(ID %in% demix$sample) %>%
  pull(ID) %>%
  unique()

# PLOT
log_info("Plotting summary demixing")
demix_plot <- demix %>%
  mutate(
    lin_2 = case_when(
      lineages %in% main_lineages ~ lineages,
      TRUE ~ "Other"
    )
  ) %>%
  ggplot() +
  aes(
    x = factor(sample, date_order),
    y = as.numeric(abundances),
    fill = lin_2
  ) +
  scale_fill_viridis_d(option = "Mako") +
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
      read_csv(snakemake@params[["metadata"]]),
      ID,
      CollectionDate
    ),
    by = c("sample" = "ID")
  ) %>%
  write.csv(snakemake@output[["table"]], row.names = FALSE)
