#!/usr/bin/env Rscript

# Write stdout and stderr to log file
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, type = "output")

library(tidyverse)
library(stringi)
library(ggpubr)
library(jsonlite)
library(logger)
log_threshold(INFO)


# Import file with plots style
source(snakemake@params[["design"]])

log_info("Reading inputs")

# Anotation for genome scheme
coordinates <- read_json(snakemake@input$coordinates)

vcf <- read_delim(snakemake@input[["vcf"]])
vcf_snp <- vcf
window <- read_csv(snakemake@input[["window"]])
metadata <- read_csv(snakemake@input[["metadata"]])

# DATA PROCESSING
# Obtain sample names ordered by CollectionDate
date_order <- metadata %>%
  filter(ID %in% snakemake@params[["samples"]]) %>%
  arrange(CollectionDate) %>%
  pull(ID) %>%
  unique()

empty_vcf <- tibble(
  REGION = date_order,
  variant = as.character(NA),
  ALT_FREQ = as.numeric(NA),
  GFF_FEATURE = as.character(NA),
  synonimous = as.character(NA),
  POS = as.numeric(NA),
  ALT = as.character(NA),
  NV_class = as.character(NA),
  group = as.character(NA)
)

# Create SNP variable and select useful variables
vcf <- vcf %>%
  dplyr::select(
    variant,
    REGION,
    ALT_FREQ,
    GFF_FEATURE,
    synonimous,
    POS,
    ALT
  )

# Create dataframe with gene length for scheme
notation_empty <- data.frame(
  gene = "",
  len = 0
) %>%
  filter(len != 0)

notation <- lapply(
  names(coordinates),
  function(x) {
    add_row(
      notation_empty,
      gene = x,
      len = length(coordinates[[x]])
    )
  }
) %>%
  bind_rows()

# Classifying variants
log_info("Classifying variants")
vcf <- vcf %>%
  mutate(
    NV_class = case_when(
      str_detect(ALT, fixed("-")) |
        str_detect(ALT, fixed("+")) ~
        "INDEL",
      TRUE ~ "SNP"
    ),
    Class = case_when(
      GFF_FEATURE == "Intergenic" ~ "Intergenic",
      TRUE ~ synonimous
    ),
    POS = as.numeric(POS)
  ) %>%
  rowwise() %>%
  mutate(
    indel_len = case_when(
      NV_class == "INDEL" ~ str_length(ALT) - 1
    ),
    indel_class = case_when(
      GFF_FEATURE == "Intergenic" ~ "Intergenic",
      NV_class == "INDEL" &
        indel_len %% 3 == 0 ~
        "In frame",
      NV_class == "INDEL" &
        indel_len %% 3 > 0 ~
        "Frameshift"
    )
  ) %>%
  ungroup() %>%
  mutate(
    group = case_when(
      GFF_FEATURE == "Intergenic" ~ "Intergenic",
      NV_class == "SNP" ~ Class,
      NV_class == "INDEL" ~ indel_class
    )
  )

# NSP data
npc <- read_csv(snakemake@input[["regions"]]) %>%
  mutate(
    summary_nsp = case_when(
      region %in% paste("nsp", seq(4, 12, 1), sep = "") ~ "nsp4-12",
      region %in% paste("nsp", seq(14, 16, 1), sep = "") ~ "nsp14-16",
      TRUE ~ region
    ),
    summary_start = case_when(
      region %in% paste("nsp", seq(4, 12, 1), sep = "") ~ 8555,
      region %in% paste("nsp", seq(14, 16, 1), sep = "") ~ 18040,
      TRUE ~ start
    ),
    summary_end = case_when(
      region %in% paste("nsp", seq(4, 12, 1), sep = "") ~ 16236,
      region %in% paste("nsp", seq(14, 16, 1), sep = "") ~ 21552,
      TRUE ~ end
    )
  ) %>%
  filter(region != "nsp1")


# PLOTS
## SUMMARY FIGURE FOR WHOLE GENOME
log_info("Plotting summary figure for whole genome")

if (nrow(vcf) == 0) {
  log_warn("Whole-genome VCF has no rows")
  vcf <- empty_vcf
}

variants <- vcf %>%
  filter(ALT_FREQ > 0) %>%
  ggplot() +
  aes(
    x = POS,
    y = factor(REGION, date_order),
    shape = factor(NV_class, c("SNP", "INDEL")),
    color = group,
    alpha = ALT_FREQ
  ) +
  geom_point(size = 3) +
  geom_col(
    data = notation,
    aes(
      x = len,
      y = 0.3,
      fill = factor(gene, rev(names(coordinates)))
    ),
    inherit.aes = FALSE,
    width = 0.3
  ) +
  scale_fill_manual(values = GENE_PALETTE) +
  xlim(c(0, 29903)) +
  scale_color_manual(
    labels = NV_TYPE_NAMES,
    values = NV_TYPE_PALETTE
  ) +
  scale_y_discrete(drop = FALSE) +
  labs(
    x = "SARS-CoV-2 genome position",
    y = "Sample",
    shape = "Variant class",
    color = "Classification",
    alpha = "Frequency",
    fill = "Region"
  ) +
  guides(
    fill = guide_legend(reverse = TRUE)
  )

# Window plot
window_plot <- window %>%
  ggplot() +
  aes(
    x = position,
    y = fractions,
    color = gen
  ) +
  geom_point() +
  geom_line(
    aes(group = 1),
    colour = "black",
    alpha = 0.3
  ) +
  scale_y_continuous(
    label = scales::percent,
    limits = c(0, max(window$fractions) + 0.005)
  ) +
  xlim(c(0, 29903)) +
  scale_color_manual(values = GENE_PALETTE) +
  labs(
    y = "Proportion of \n sites with NV",
    x = "",
    color = "Gen"
  )

# Add nsp info
window_plot_nsp <- window_plot +
  geom_vline(
    data = npc,
    aes(xintercept = summary_start),
    color = "red"
  ) +
  geom_vline(
    data = npc,
    aes(xintercept = summary_end),
    color = "red"
  ) +
  geom_text(
    data = npc,
    aes(
      x = (summary_start + summary_end) / 2,
      y = max(window$fractions) + 0.002,
      label = summary_nsp
    ),
    inherit.aes = FALSE,
    size = 5,
    angle = 60
  )

figura <- ggarrange(
  window_plot_nsp,
  variants,
  nrow = 2,
  align = "v",
  legend.grob = get_legend(variants),
  heights = c(2, 6),
  legend = "right",
  labels = c("A", "B")
)

ggsave(
  filename = snakemake@output[["fig"]],
  plot = figura,
  width = 250,
  height = 240,
  units = "mm",
  dpi = 250
)


# Zoom in in spike
log_info("Plotting summary for variants in the spike")

spike_pos <- window %>%
  filter(gen == "S") %>%
  pull(position)

vcf_spike <- vcf %>%
  filter(
    ALT_FREQ > 0,
    POS %in% c(min(spike_pos):max(spike_pos))
  )

window_plot_spike <- window %>%
  filter(gen == "S") %>%
  ggplot() +
  aes(
    x = position,
    y = fractions,
    color = gen
  ) +
  geom_point() +
  geom_line(
    aes(group = 1),
    colour = "black",
    alpha = 0.3
  ) +
  scale_y_continuous(
    label = scales::percent,
    limits = c(0, max(window$fractions) + 0.005)
  ) +
  xlim(c(min(spike_pos), max(spike_pos))) +
  scale_color_manual(
    values = GENE_PALETTE,
    guide = "none"
  ) +
  labs(
    y = "Proportion of \n sites with NV",
    x = ""
  )

if (nrow(vcf_spike) == 0) {
  log_warn("Spike VCF has no rows")
  vcf_spike <- empty_vcf
}

variants_spike <- vcf_spike %>%
  ggplot() +
  aes(
    x = POS,
    y = factor(REGION, date_order),
    shape = factor(NV_class, c("SNP", "INDEL")),
    color = group,
    alpha = ALT_FREQ
  ) +
  geom_point(size = 3) +
  xlim(c(min(spike_pos), max(spike_pos))) +
  scale_color_manual(
    labels = NV_TYPE_NAMES,
    values = NV_TYPE_PALETTE
  ) +
  scale_y_discrete(drop = FALSE) +
  labs(
    x = "SARS-CoV-2 genome position",
    y = "Sample",
    shape = "Variant class",
    color = "Classification",
    alpha = "Frequency",
    fill = "Region"
  ) +
  guides(
    fill = guide_legend(reverse = TRUE)
  )

figura_spike <- ggarrange(
  window_plot_spike,
  variants_spike,
  nrow = 2,
  align = "v",
  legend.grob = get_legend(variants),
  heights = c(2, 6),
  legend = "right",
  labels = c("A", "B")
)

ggsave(
  filename = snakemake@output[["fig_s"]],
  plot = figura_spike,
  width = 250,
  height = 240,
  units = "mm",
  dpi = 250
)


# Figure for no of heterozygous sites for each sample
log_info("Plotting no. of heterozygous sites for each sample")
figur_SNP_table <- vcf_snp %>%
  filter(ALT_FREQ <= snakemake@params$max_alt_freq) %>%
  left_join(
    metadata,
    by = c("REGION" = "ID")
  ) %>%
  group_by(REGION) %>%
  summarise(
    CollectionDate = min(as.Date(CollectionDate)),
    n = n_distinct(POS)
  ) %>%
  ungroup()

if (nrow(figur_SNP_table) == 0) {
  log_warn("Filtered SNP table has no rows")
  figur_SNP_table <- empty_vcf
}

figur_SNP_time <- figur_SNP_table %>%
  ggplot() +
  aes(
    x = CollectionDate,
    y = n
  ) +
  geom_smooth(
    method = "lm",
    fill = "gray95",
    alpha = 0.6
  ) +
  geom_point() +
  labs(
    x = "Collection date",
    y = "No. of polimorphic sites"
  )

ggsave(
  filename = snakemake@output[["fig_cor"]],
  plot = figur_SNP_time,
  width = 250,
  height = 119.4,
  units = "mm",
  dpi = 250
)


# PLOT TABLES
log_info("Saving plot tables")

# Variants summary table
vcf %>%
  select(
    REGION,
    POS,
    variant,
    ALT_FREQ,
    NV_class,
    group
  ) %>%
  rename(
    sample = REGION,
    Variant = variant,
    Class = group
  ) %>%
  filter(ALT_FREQ > 0) %>%
  mutate(
    Class = case_when(
      Class == "Yes" ~ "synonymous",
      Class == "No" ~ "non_synonymous",
      TRUE ~ Class
    )
  ) %>%
  write.csv(snakemake@output[["table_2"]], row.names = FALSE)

# Window plot table
window %>%
  transmute(
    POS = position,
    feature = gen,
    prop_PolymorphicSites = fractions
  ) %>%
  write.csv(snakemake@output[["table_1"]], row.names = FALSE)

# Heterozygous sites per sample table
vcf_snp %>%
  filter(ALT_FREQ <= snakemake@params$max_alt_freq) %>%
  select(!GFF_FEATURE) %>%
  left_join(
    metadata,
    by = c("REGION" = "ID")
  ) %>%
  group_by(REGION) %>%
  summarise(
    CollectionDate = min(as.Date(CollectionDate)),
    n_PolymorphicSites = n_distinct(POS)
  ) %>%
  ungroup() %>%
  rename(sample = REGION) %>%
  unique() %>%
  write.csv(snakemake@output[["table_3"]], row.names = FALSE)

# Stats for reporting
n_indels <- vcf %>%
  filter(NV_class == "INDEL") %>%
  pull(variant) %>%
  unique() %>%
  length()

n_snv <- length(unique(vcf$variant)) - n_indels
model <- lm(n ~ CollectionDate, data = figur_SNP_table)

# Calculate correlation, if possible
if (nrow(figur_SNP_table) > 2) {
  p.value <- summary(model)$coefficients[2,4]
} else {
  p.value <- NA
}

list(
  "INDELS" = n_indels,
  "SNV" = n_snv,
  "window" = snakemake@params[["window"]],
  "step" = snakemake@params[["step"]],
  "r2" = summary(model)$r.squared[[1]],
  "value" = ifelse(p.value < 0.001, "< 0.001", p.value)
) %>%
  toJSON() %>%
  write(snakemake@output[["json"]])
