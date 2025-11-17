#!/usr/bin/env Rscript

# Write stdout and stderr to log file
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, type = "output")

set.seed(snakemake@params$seed)

library(ape)
library(pegas)
library(future)
library(future.apply)
library(readr)
library(jsonlite)
library(logger)

log_threshold(INFO)

# Pi calculation
nucleotide.diversity <- function(dna_object, record.names, sample.size) {
  sample <- sample(record.names, sample.size, replace = FALSE)
  dna_subset <- dna_object[record.names %in% sample]
  nuc.div(dna_subset)
}

# Parallel bootstrapping
boot.nd.parallel <- function(aln, sample.size, reps) {
  record.names <- names(aln)
  future_sapply(
    1:reps,
    function(x) {
      nucleotide.diversity(aln, record.names, sample.size)
    },
    future.seed = TRUE
  )
}

# Read outgroup/context alignment
log_info("Reading context")
gene_ex <- read.dna(
  snakemake@input[["context_fasta"]],
  format = "fasta",
  as.matrix = FALSE
)
gene_ex <- gene_ex[
  !startsWith(names(gene_ex), snakemake@params$aln_reference)
]

# Read target (study) alignment
log_info("Reading target alignment")
study_aln <- read.dna(
  snakemake@input[["study_fasta"]],
  format = "fasta",
  as.matrix = FALSE
)
study_aln <- study_aln[
  !startsWith(names(study_aln), snakemake@params$aln_reference)
]

# Diversity value for our samples
log_info("Calculating diversity value for studied samples")
diversity <- nuc.div(study_aln)

# Perform bootstrap
log_info("Performing calculation for nucleotide diversity in context samples")
plan(multisession, workers = snakemake@threads)
divs <- boot.nd.parallel(
  gene_ex,
  length(study_aln),
  snakemake@params[["bootstrap_reps"]]
)
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

# Data for JSON file
log_info("Building JSON data")
p.value <- ifelse(st$p.value >= 0.05, pvalue.norm, pvalue.emp)
list.div <- list(
  "diversity" = diversity,
  "p.value" = ifelse(p.value >= 0.001, p.value, "< 0.001"),
  "normal.pvalue" = ifelse(st$p.value >= 0.001, p.value, "< 0.001"),
  "norm.text" = ifelse(st$p.value >= 0.05, "", "not"),
  "type.test" = ifelse(st$p.value >= 0.05, "", "empirical"),
  "boot.reps" = snakemake@params[["bootstrap_reps"]],
  "sample.size" = length(study_aln)
)

log_info("Writing diversity distribution")
write_lines(divs, snakemake@output[["divs"]])

log_info("Writing results JSON")
write_json(
  list.div,
  snakemake@output[["json"]],
  auto_unbox = TRUE,
  digits = NA  # maximum precision
)
