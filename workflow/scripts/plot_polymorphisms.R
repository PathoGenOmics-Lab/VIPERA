library(tidyverse)
library(ggrepel)
library(logger)

# Write stdout and stderr to log file
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, type = "output")

log_threshold(INFO)
NUCLEOTIDES <- c("A", "C", "G", "T")

log_info("Reading tables")
polymorphisms <- read_delim(snakemake@input[["table"]], show_col_types = FALSE)
metadata <- read_delim(
    snakemake@input[["metadata"]],
    show_col_types = FALSE,
    col_select = c("Virus name", "CollectionDate")
)

log_info("Joining tables")
sites.with.metadata <- polymorphisms %>%
    left_join(metadata, by = c("SampleID" = "Virus name"))

log_info("Plotting")
p <- sites.with.metadata %>%
    pivot_longer(
        -c(SampleID, CollectionDate),
        names_to = "position",
        values_to = "allele") %>%
    mutate(
        snv = paste0(position, allele),
        position = as.numeric(position)) %>%  # TODO: use Wuhan reference to build it like A123T
    ggplot(aes(x = CollectionDate, y = position, color = SampleID)) +
        geom_point(size = 1) +
        geom_label_repel(aes(label = snv)) +
        xlab("Collection date") +
        ylab("Position in genome") +
        theme_bw()

log_info("Saving plot")
ggsave(
    snakemake@output[["plot"]],
    plot = p,
    width = snakemake@params[["width"]],
    height = snakemake@params[["height"]],
    units = snakemake@params[["units"]]
)

log_info("Done")
