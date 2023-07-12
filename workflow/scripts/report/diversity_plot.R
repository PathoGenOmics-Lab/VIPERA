# LIBRERIAS #######

library(ape)
library(pegas)
library(future.apply)
library(tidyverse)


# Pi calculation
nucleotide.diversity <- function(dna_object, record.names, sample.size){
  sample <- sample(record.names, sample.size, replace = F)
  dna_subset <- dna_object[record.names %in% sample]
  nuc.div(dna_subset)
}

# Parallel bootstrapping
boot.nd.parallel <- function(aln, sample.size = 12, reps = 100) {
  record.names <- names(aln)
  future_sapply(
    1:reps,
    function(x) nucleotide.diversity(aln, record.names, sample.size),
    future.seed = TRUE
  )
}

# Plot design
source(snakemake@params[["design"]])

# Outgroup/context alignment
gene_ex <- read.dna(snakemake@input[["context_fasta"]], format = "fasta", as.matrix = F)
gene_ex <- gene_ex[!startsWith(names(gene_ex), snakemake@config[["ALIGNMENT_REFERENCE"]])]

# Study alignment
study_aln <- read.dna(snakemake@input[["study_fasta"]],format = "fasta", as.matrix = F)
study_aln <- study_aln[!startsWith(names(study_aln), snakemake@config[["ALIGNMENT_REFERENCE"]])]

# Diversity value for our samples
diversity <- nuc.div(study_aln)
write.table(data.frame(div = diversity), snakemake@output[["value"]], row.names = F)

# Perform bootstrap
plan(multisession, workers = snakemake@threads)
divs <- boot.nd.parallel(gene_ex, length(study_aln), snakemake@params[["bootstrap_reps"]])
plan(sequential)

# Test normality
st <- shapiro.test(divs)

# Calculate p-value (assuming normal distribution)
standardized.value <- (diversity - mean(divs)) / sd(divs)
pvalue.norm <- pnorm(standardized.value)

# Estimate p-value empirically
empirical.probs <- ecdf(divs)
pvalue.emp <- empirical.probs(diversity)

# Plot and save
p <- data.frame(pi = divs) %>%
  ggplot() +
  geom_density(aes(x = pi), fill = "#fcbf49", alpha = 0.7, bw = 0.000001, color = "#eae2b7") +
  geom_vline(aes(xintercept = diversity), color = "#d62828") +
  stat_function(fun = dnorm, args = list(mean = mean(divs), sd = sd(divs)), color = "#f77f00") +
  labs(x = "π", y = "Density") +
  ggtitle(
    "Diversity distribution",
    paste0(
      "Study π: ", prettyNum(diversity), "\n",
      "eCDF-estimated p =", prettyNum(pvalue.emp), " (left-tailed)", "\n",
      "Assumed normal p = ", prettyNum(pvalue.norm), " (left-tailed)", "\n",
      "Normal distribution: ", st$p.value >= 0.05, " (Shapiro-Wilk test, p = ", prettyNum(st$p.value), ")", "\n",
      snakemake@params[["bootstrap_reps"]], " reps with size=", length(study_aln)
    )
  ) +
  theme(
    plot.title = element_text(size = 10),
    plot.subtitle = element_text(size = 8)
  )

ggsave(
  filename = snakemake@output[["fig"]],
  plot = p,
  width = snakemake@params[["plot_width"]],
  height = snakemake@params[["plot_height"]],
  units = "mm",
  dpi = 250
)
