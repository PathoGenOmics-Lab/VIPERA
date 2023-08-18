#!/usr/bin/env Rscript


library(GISAIDR)
library(tidyverse)
library(glue)
library(logger)
log_threshold(INFO)

CHKPT.ERROR.MSG <- paste(
    "Please provide a sequence dataset through the CONTEXT_FASTA parameter",
    "by editing the config/targets.yaml or via the command line (with --config)"
)

chunk <- function(x, chunk_length = snakemake@params[["chunk_length"]]) {
    split(x, ceiling(seq_along(x) / chunk_length))
}


# ==== #
# MAIN #
# ==== #

# Write stdout and stderr to log file
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, type = "output")

# Read original sample metadata
sample.metadata <- read_delim(snakemake@input[["metadata"]], show_col_types = FALSE)

# Checkpoint: needed columns exist
needed.columns <- c(
    snakemake@params[["date_column"]],
    snakemake@params[["location_column"]],
    snakemake@params[["samples_gisaid_accession_column"]]
)
needed.columns.mask <- needed.columns %in% colnames(sample.metadata)
if (!all(needed.columns.mask)) {
    print(glue("Please ensure column '{needed.columns[!needed.columns.mask]}' is present"))
    stop(glue("Missing columns in '{snakemake@input[['metadata']]}'. Alternatively:\n{CHKPT.ERROR.MSG}"))
}

log_info("Getting time windows for context samles")
# Get time windows
dates <- sample.metadata %>%
    pull(snakemake@params[["date_column"]]) %>%
    as.numeric
window.quantile.offset <- (1 - snakemake@params[["date_window_span"]]) / 2
min.date <- as_date(quantile(dates, window.quantile.offset))
max.date <- as_date(quantile(dates, 1 - window.quantile.offset))
padded.min.date <- min.date - snakemake@params[["date_window_paddding_days"]]
padded.max.date <- max.date + snakemake@params[["date_window_paddding_days"]]
print(glue("Time window (span={snakemake@params[['date_window_span']]}): {round(interval(min.date, max.date) / days(1))} days (from {min.date} to {max.date})"))
print(glue("Padded time window (padding={snakemake@params[['date_window_paddding_days']]} days): {round(interval(padded.min.date, padded.max.date) / days(1))} days (from {padded.min.date} to {padded.max.date})"))

log_info("Getting locations for context samples")
# Get locations (if there are multiple, sample from all of them)
locations <- sample.metadata %>%
    pull(snakemake@params[["location_column"]]) %>%
    unique

# GISAID login
creds.list <- yaml::read_yaml(snakemake@params[["gisaid_creds"]])
credentials <- login(
    username = creds.list[["USERNAME"]],
    password = creds.list[["PASSWORD"]]
)

# Get accession IDs (EPI_ISL codes)
dataframes <- lapply(
    locations,
    function(location) {
        query(
            credentials = credentials,
            from = as.character(strftime(padded.min.date, format = "%Y-%m-%d")),
            to = as.character(strftime(padded.max.date, format = "%Y-%m-%d")),
            location = location,
            fast = TRUE,
            low_coverage_excl = snakemake@params[["exclude_low_coverage"]],
            complete = snakemake@params[["complete"]],
            collection_date_complete = snakemake@params[["collection_date_complete"]],
            high_coverage = snakemake@params[["high_coverage"]]
        )
    }
)

# Join results
metadata <- bind_rows(dataframes)

log_info("Removeing overlapping sequences")
# Checkpoint: remove samples that overlap with target samples according to GISAID ID
samples.accids <- sample.metadata %>%
    pull(snakemake@params[["samples_gisaid_accession_column"]])
filtered.accids <- metadata %>% filter(accession_id %in% samples.accids)
metadata <- metadata %>% filter(!accession_id %in% samples.accids)
print(glue("{nrow(metadata)} accession_ids remaining after GISAID ID filter"))

# Checkpoint: enforce a minimum number of samples to have at least
# as many possible combinations as bootstrap replicates.
# This is done by calculating the root of a function based on the
# formula for calculating combinations with replacement
# for n ≥ r ≥ 0: combinations with replacement = n! / (r! (n-r)!)
r <- nrow(sample.metadata)
min.comb <- snakemake@params[["min_theoretical_combinations"]]
solution <- uniroot(
    function(n) {
        factorial(n) / (factorial(r) * factorial(n - r)) - min.comb
    },
    lower = r,   # determined by sample size (n ≥ r)
    upper = 170  # determined by default number precision
)
if (nrow(metadata) < floor(solution$root)) {
    stop(glue("Too few available samples (n={nrow(metadata)}).\n{CHKPT.ERROR.MSG}"))
}

# Split sequences in chunks of sequences
chunks <- chunk(metadata$accession_id)

# Download chunks from GISAID and get metadata
downloads <- do.call(
    rbind,
    lapply(
        seq_along(chunks),
        function(x) {
            res <- download(
                credentials = login(
                    username = creds.list[["USERNAME"]],
                    password = creds.list[["PASSWORD"]]
                ),
                list_of_accession_ids = chunks[[x]]
            )
            time.sleep <- runif(
                1,
                min = snakemake@params[["min_sleep"]],
                max = snakemake@params[["max_sleep"]]
            )
            Sys.sleep(time.sleep)
            res
        }
    )
)

# Filter host
downloads <- filter(downloads, host == snakemake@params[["host"]])

log_info("Saving fasta file")
# Export sequences in FASTA format
cat("", file = snakemake@output[["fasta"]])
for (i in seq_len(nrow(downloads))) {
    # Write header
    cat(
        ">", downloads$accession_id[i], "\n",
        sep = "",
        file = snakemake@output[["fasta"]],
        append = TRUE
    )
    # Write sequence
    cat(
        downloads$sequence[i], "\n",
        sep = "",
        file = snakemake@output[["fasta"]],
        append = TRUE
    )
}

# Remove sequence data from metadata dataframe
downloads$sequence <- NULL

# Write metadata
# TODO: select columns
downloads %>%
    write_csv(
        snakemake@output[["metadata"]],
        na = "",
        col_names = TRUE,
        progress = FALSE
    )

# Write filtered IDs
write_lines(filtered.accids, snakemake@output[["duplicate_accids"]])
