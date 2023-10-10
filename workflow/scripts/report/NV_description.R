#!/usr/bin/env Rscript

library(tidyverse)
library(stringi)
library(ggpubr)
library(jsonlite)
library(logger)
log_threshold(INFO)


# Import file with plots style
source(snakemake@params[["design"]])

# Write stdout and stderr to log file
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, type = "output")


# SARS-CoV-2 anotation for genome scheme
SCov2_annotation <- list(
  "five_prime_UTR"  = c(1:265),
  "orf1ab"          = c(266:21555),
  "Intergenic_1"    = c(21556:21562),
  "S"               = c(21563:25384),
  "Intergenic_2"    = c(25385:25392),
  "ORF3a"           = c(25393:26220),
  "Intergenic_3"    = c(26221:26244),
  "E"               = c(26245:26472),
  "Intergenic_4"    = c(26473:26522),
  "M"               = c(26523:27191),
  "Intergenic_5"    = c(27192:27201),
  "ORF6"            = c(27202:27387),
  "Intergenic_6"    = c(27388:27393),
  "ORF7"            = c(27394:27759),
  "Intergenic_7"    = c(27760:27893),
  "ORF8"            = c(27894:28259),
  "Intergenic_8"    = c(28260:28273),
  "N"               = c(28274:29533),
  "Intergenic_9"    = c(29534:29557),
  "ORF10"           = c(29558:29674),
  "three_prime_UTR" = c(29675:29903))



vcf <- read_delim(snakemake@input[["vcf"]])
vcf_snp <- vcf
window <- read_csv(snakemake@input[["window"]])

# DATA PROCESSING

# Obtain sample names ordered by CollectionDate
date_order <- read_csv(snakemake@params[["metadata"]]) %>%
  filter(ID %in% snakemake@params[["samples"]]) %>%
  arrange(CollectionDate) %>%
  pull(ID) %>%
  unique()

# Create SNP variable and select useful variables
vcf <- vcf %>%
  mutate(
    SNP = paste(REF, POS, ALT, sep = "-")
    ) %>%
  dplyr::select(
    SNP,
    REGION,
    ALT_FREQ,
    GFF_FEATURE,
    synonimous
    ) %>%
  rowwise() %>%
  mutate(POS = strsplit(SNP, "-")[[1]][2]) %>%
  ungroup()


# Df with gene length for scheme
notation_empty <- data.frame(
  gene = "",
  len = 0
  ) %>%
  filter(len != 0)

notation <- lapply(
  names(SCov2_annotation),
  function(x) {
   add_row(
    notation_empty,
    gene = x,
    len = length(SCov2_annotation[[x]])
   )
  }
) %>%
bind_rows()


# Classifying variants
log_info("Classifying variants")
vcf <- vcf %>%
  mutate(
    NV_class = case_when(
      str_detect(SNP, fixed("--")) |
      str_detect(SNP, fixed("+")) ~ "INDEL",
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
      NV_class == "INDEL" &
        str_detect(SNP, fixed("--")) ~
        str_length(strsplit(SNP, "--")[[1]][2]),
      NV_class == "INDEL" &
        str_detect(SNP, fixed("-+")) ~
        str_length(strsplit(SNP, "-+")[[1]][2])
      ),
    indel_class = case_when(
      GFF_FEATURE == "Intergenic" ~ "Intergenic",
      NV_class == "INDEL" &
        indel_len %% 3 == 0 ~ "In frame",
      NV_class == "INDEL" &
        indel_len %% 3 > 0 ~ "Frameshift"
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
npc <- read_csv(snakemake@params[["nsp"]]) %>%
  mutate(
    summary_nsp = case_when(
      NSP %in% paste("nsp", seq(4, 12, 1), sep = "") ~ "nsp4-12",
      NSP %in% paste("nsp", seq(14, 16, 1), sep = "") ~ "nsp14-16",
      TRUE ~ NSP
      ),
    summaary_start = case_when(
      NSP %in% paste("nsp", seq(4, 12, 1), sep = "") ~ 8555,
      NSP %in% paste("nsp", seq(14, 16, 1), sep = "") ~ 18040,
      TRUE ~ POS_i
      ),
    summaary_end = case_when(
      NSP %in% paste("nsp", seq(4, 12, 1), sep = "") ~ 16236,
      NSP %in% paste("nsp", seq(14, 16, 1), sep = "") ~ 21552,
      TRUE ~ POS_f
      )
    ) %>%
filter(NSP != "nsp1")


# PLOTS

## SUMMARY FIGURE FOR WHOLE GENOME

log_info("Plotting summary figure for whole genome")
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
      fill = factor(gene, rev(names(SCov2_annotation)))
      ),
      inherit.aes = FALSE,
      width = 0.3
      ) +
  scale_fill_manual(values = gene_colors) +
  xlim(c(0, 29903)) +
  scale_color_manual(
    labels = NV_names,
    values = NV_colors
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
    limits = c(0, max(window$fractions) + 0.005)) +
  xlim(c(0, 29903)) +
  scale_color_manual(values = gene_colors) +
  labs(
    y = "Proportion of \n sites with SNV",
    x = "",
    color = "Gen")

# Add nsp info
window_plot_nsp <- window_plot +
  geom_vline(
    data = npc,
    aes(xintercept = summaary_start),
    color = "red"
    ) +
  geom_vline(
    data = npc,
    aes(xintercept = summaary_end),
    color = "red"
    ) +
  geom_text(
    data = npc,
    aes(
      x = (summaary_start + summaary_end) / 2,
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
  dpi = 250)


# Zoom in in spike
log_info("Plotting summary for variants in the spike")
spike.pos <- window %>%
  filter(gen == "S") %>%
  pull(position)


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
  xlim(c(min(spike.pos), max(spike.pos))) +
  scale_color_manual(
    values = gene_colors,
    guide = "none"
    ) +
  labs(
    y = "Proportion of \n sites with SNV",
    x = ""
    )

variants_spike <- vcf %>%
  filter(
    ALT_FREQ > 0,
      POS %in% c(min(spike.pos):max(spike.pos))
      ) %>%
  ggplot() +
  aes(
    x = POS,
    y = factor(REGION, date_order),
    shape = factor(NV_class, c("SNP", "INDEL")),
    color = group,
    alpha = ALT_FREQ
    ) +
  geom_point(size = 3) +
  xlim(c(min(spike.pos), max(spike.pos))) +
  scale_color_manual(
    labels = NV_names,
    values  = NV_colors
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


# Figure for nº of heterozygus sites for each sample
log_info("Plotting nº of heterozygus sites for each sample")
figur_SNP_table <- vcf_snp %>%
  filter(ALT_FREQ <= 0.95) %>%
  select(!GFF_FEATURE) %>%
  unique() %>%
  left_join(
    read_csv(snakemake@params[["metadata"]]),
    by = c("REGION" = "ID")
    ) %>%
  group_by(REGION) %>%
  summarise(
    CollectionDate = min(as.Date(CollectionDate)),
    n = n_distinct(POS)
    ) %>%
  ungroup()

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
    x = "Date",
    y = "Nº of polimorphic sites"
    )

ggsave(
  filename = snakemake@output[["fig_cor"]],
  plot = figur_SNP_time,
  width = 250,
  height = 119.4,
  units = "mm",
  dpi = 250)


# PLOT TABLES
log_info("Saving plot tables")


# Variants summary table
vcf %>%
  select(
    REGION,
    POS,
    SNP,
    ALT_FREQ,
    NV_class,
    group
  ) %>%
  rename(
    sample = REGION,
    NV = SNP,
    Class = group
  ) %>%
  filter(ALT_FREQ > 0) %>%
  mutate(
    Class = case_when(
      Class == "yes" ~ "synonymous",
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

# Heterzygus sites per sample table
vcf_snp %>%
  filter(ALT_FREQ <= 0.95) %>%
  select(!GFF_FEATURE) %>%
  left_join(
    read_csv(snakemake@params[["metadata"]]),
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


# STATS FOR REPORTING

n_indels <- vcf %>%
  filter(NV_class == "INDEL") %>%
  pull(SNP) %>%
  unique() %>%
  length()

n_snv <- length(unique(vcf$SNP)) - n_indels
model <- lm(n ~ CollectionDate, data = figur_SNP_table)
cortest <- cor.test(figur_SNP_table$n, as.numeric(figur_SNP_table$CollectionDate))
list(
  "INDELS" = n_indels,
  "SNV" = n_snv,
  "window" = snakemake@params[["window"]],
  "step"   = snakemake@params[["step"]],
  "r2"     = summary(model)$r.squared[[1]],
  "value"  = ifelse(cortest$p.value < 0.001, "< 0.001", cortest$p.value)
) %>%
toJSON() %>%
write(snakemake@output[["json"]])
