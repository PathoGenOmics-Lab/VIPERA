#!/usr/bin/env Rscript

library(tidyverse)
library(stringi)
library(ggrepel)
library(logger)
log_threshold(INFO)

# Get colors
color <- grDevices::colors()[grep("gr(a|e)y", grDevices::colors(), invert = TRUE)]
color <- color[grep("white", color, invert = TRUE)]
# Import file with plots style
source(snakemake@params[["design"]])

# Write stdout and stderr to log file
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, type = "output")

# DATA PREPROCESSING #####
metadata <- read_csv(snakemake@input[["metadata"]])

vcf <- read_delim(snakemake@input[["vcf"]])
data <- metadata %>%
  filter(
    ID %in% vcf$REGION
  ) %>%
  dplyr::select(
    ID,
    CollectionDate
  )

# Obtain sample names ordered by CollectionDate
date_order <- metadata %>%
  arrange(CollectionDate) %>%
  pull(ID) %>%
  unique()

# Simplify features names and create SNP variable
vcf <- vcf %>%
  dplyr::select(
    variant,
    REGION,
    ALT_FREQ,
    POS
    )

# Get a list with studied samples ID
IDs <- pull(
    vcf,
    REGION
    ) %>%
  unique()

vcf <- vcf %>%
  pivot_wider( # Obtain 0 in positions without freq
    names_from = REGION,
    values_from = ALT_FREQ,
    values_fill = 0
    ) %>%
  pivot_longer(
    IDs,
    names_to = "REGION",
    values_to = "ALT_FREQ"
    ) %>%
  ungroup()

# Join variants file and metadata file
vcf <- left_join(vcf, data, by = c("REGION" = "ID"))

# Calculate days since first sample
vcf <- arrange(
    vcf,
    CollectionDate
    ) %>%
  mutate(
    interval = as.numeric(CollectionDate - min(CollectionDate))
  )

# PLOTS ####
## VOLCANO PLOT ####

# Get list with all different polymorphisms
SNPs <- pull(
    vcf,
    variant
    ) %>%
  unique()

# Create an empty dataframe to be filled
cor.df <- data.frame(
    snp = "",
    cor = 0,
    p.value = 0
    ) %>%
  filter(p.value != 0)

cor.df.fill <- lapply(
    SNPs,
    function(snp) {
      df <- filter(
        vcf,
        variant == snp
      )

      test <- cor.test(
        df$ALT_FREQ,
        df$interval
      )

      pvalue <- p.adjust(
        test$p.value,
        method = "BH",
        n = length(SNPs)
      )
      add_row(
        cor.df,
        snp = snp,
        cor = test$estimate,
        p.value = pvalue
      )
    }
  ) %>%
  bind_rows()

# Plot a pseudo volcano plot with coorrelation index and p-value
log_info("Plotting pseudovolcano figure")
volcano <- cor.df.fill %>%
  mutate(
    trans.p = -log10(p.value),
    label = case_when(p.value < 0.05 ~ snp)
  ) %>%
  ggplot() +
  aes(
    x = cor,
    y = trans.p
  ) +
  geom_point() +
  geom_label_repel(aes(label = label)) +
  xlim(c(-1, 1)) +
  geom_hline(
    aes(yintercept = -log10(0.05)),
    linetype = 2,
    color = "red"
  ) +
  labs(
    x = "Correlation",
    y = "-log10(p-value)"
  )

ggsave(
  filename = snakemake@output[["pseudovolcano"]],
  plot = volcano,
  width = 159.2,
  height = 119.4,
  units = "mm",
  dpi = 250
)


## SNP PANEL ####
# SNPs significantly correlated with time
sign <- filter(
            cor.df.fill,
            p.value < 0.05
            ) %>%
  pull(snp) %>%
  unique()

# SNPs which are in positions with more than one alternative allele
dup <- vcf %>%
  select(
    variant,
    POS
  ) %>%
  unique() %>%
  group_by(POS) %>%
  filter(n() > 1) %>%
  ungroup() %>%
  pull(variant) %>%
  unique()

subset <- c(sign, dup) %>%
  unique()

# Plot height depending on the number of SNPs assuming 4 columns in the plot
plot.height <- ceiling(length(subset) / 4) * 42

log_info("PLotting SNPs trends in time")
panel <- vcf %>%
  filter(variant %in% subset) %>%
  ggplot() +
  aes(
    x = interval,
    y = ALT_FREQ,
    color = variant
  ) +
  scale_color_manual(values = sample(color, length(subset))) +
  geom_point() +
  geom_line() +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 9)
  ) +
  labs(
    x = "Days since first sample",
    y = "Frequency",
    color = "NV") +
  guides(color = guide_legend(ncol = 3))

if (length(subset) > 1) {
  panel <- panel +
    facet_wrap(
      vars(POS),
      nrow = ceiling(length(subset) / 4),
      ncol = 4
    )
}

ggsave(
  filename = snakemake@output[["snp_panel"]],
  plot = panel, width = 159.2,
  height = max(100, plot.height),
  units = "mm",
  dpi = 250
)


# PLOT TABLES ####
log_info("Saving coorelation table")
cor.df.fill %>%
  transmute(
    NV = snp,
    PearsonCor = cor,
    adj_pvalue = p.value
  ) %>%
  write.csv(
    snakemake@output[["table_1"]],
    row.names = FALSE
    )

log_info("Saving SNPs trends table")
vcf %>%
  filter(variant %in% subset) %>%
  transmute(
    sample = REGION,
    POS = POS,
    NV = variant,
    ALT_FREQ = ALT_FREQ,
    DaysSinceFirst = interval
    ) %>%
    write.csv(snakemake@output[["table_2"]], row.names = FALSE)
