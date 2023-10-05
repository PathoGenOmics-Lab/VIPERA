#!/usr/bin/env Rscript

library(ape)
library(pegas)
library(future.apply)
library(tidyverse)
library(jsonlite)
library(logger)
log_threshold(INFO)


# Write stdout and stderr to log file
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, type = "output")

# Pi calculation
nucleotide.diversity <- function(dna_object, record.names, sample.size) {

  sample <- sample(record.names, sample.size, replace = FALSE)
  dna_subset <- dna_object[record.names %in% sample]
  nuc.div(dna_subset)
}

# Parallel bootstrapping
boot.nd.parallel <- function(aln, sample.size = 12, reps = 100) {

  record.names <- names(aln)
  future_sapply(
    1:reps,
    function(x) {
      nucleotide.diversity(aln, record.names, sample.size)
    },
    future.seed = TRUE
  )
}

# Import file with plots style
source(snakemake@params[["design"]])

# Outgroup/context alignment
gene_ex <- read.dna(
  snakemake@input[["context_fasta"]],
  format = "fasta",
  as.matrix = FALSE
)
gene_ex <- gene_ex[!startsWith(names(gene_ex), snakemake@config[["ALIGNMENT_REFERENCE"]])]

# Study alignment
study_aln <- read.dna(
  snakemake@input[["study_fasta"]],
  format = "fasta",
  as.matrix = FALSE
)
study_aln <- study_aln[!startsWith(names(study_aln), snakemake@config[["ALIGNMENT_REFERENCE"]])]

# Diversity value for our samples
log_info("Calculating diversity value for studied samples")
diversity <- nuc.div(study_aln)


# Perform bootstrap
log_info("Performing bootstraped calculation for nucleotide diversity in oontext samples")
plan(multisession, workers = snakemake@threads)
divs <- boot.nd.parallel(gene_ex, length(study_aln), snakemake@params[["bootstrap_reps"]])
plan(sequential)

# Test normality
log_info("Normality test for nucleotide diversity values")
st <- shapiro.test(divs)

# Calculate p-value (assuming normal distribution)

log_info("Calculating p-value (assuming normal distribution)")

test <- t.test(
  divs,
  alternative = "greater",
  mu = diversity,
  conf.level = 0.95
)

pvalue.norm <- test$p.value


# Estimate p-value empirically
log_info("Estimating p-value empirically")
empirical.probs <- ecdf(divs)
pvalue.emp <- empirical.probs(diversity)

# Plot and save
log_info("Plotting diversity plot")
p <- data.frame(pi = divs) %>%
  ggplot() +
  geom_density(
    aes(x = pi),
    fill = "#fcbf49",
    alpha = 0.7,
    bw = 0.000001,
    color = "#eae2b7"
  ) +
  geom_vline(
    aes(xintercept = diversity),
    color = "#d62828"
  ) +
  stat_function(
    fun = dnorm,
    args = list(mean = mean(divs), sd = sd(divs)),
    color = "#f77f00") +
  labs(
    x = "π",
    y = "Density")


ggsave(
  filename = snakemake@output[["fig"]],
  plot = p,
  width = snakemake@params[["plot_width"]],
  height = snakemake@params[["plot_height"]],
  units = "mm",
  dpi = 250
)

# DATA JSON #####
p.value <- ifelse(st$p.value >= 0.05, pvalue.norm,pvalue.emp)

list.div <- list(
                "diversity" = diversity,
                "p.value" = ifelse(p.value >= 0.001, p.value, "< 0.001"),
                "normal.pvalue" = ifelse(st$p.value >= 0.001, p.value, "< 0.001"),
                "norm.text" = ifelse(st$p.value >= 0.05, "", "not"),
                "type.test" = ifelse(st$p.value >= 0.05, "", "empirical"),
                "boot.reps" = snakemake@params[["bootstrap_reps"]],
                "sample.size" = length(study_aln)
)

json <- toJSON(list.div)

write(json, snakemake@output[["json"]])


# PLOT TABLES

data.frame(
  pi = divs,
  prop.value = diversity
  ) %>%
  write.csv(snakemake@output[["table"]], row.names = FALSE)
