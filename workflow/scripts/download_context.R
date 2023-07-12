#!/usr/bin/env Rscript


library(GISAIDR)
library(tidyverse)
library(glue)

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

# Read original sample metadata
sample.metadata <- read_delim(snakemake@input[["metadata"]], show_col_types = FALSE)

# Get time windows
dates <- sample.metadata %>%
    pull(snakemake@params[["date_column"]])
min.date <- min(dates, na.rm = TRUE)
max.date <- max(dates, na.rm = TRUE)
print(glue("Time window: {difftime(max.date, min.date, units = 'days')} days (from {min.date} to {max.date})"))

# Checkpoint: time window cannot be zero
if (min.date == max.date) {
    stop(glue("Time window is too short.\n{CHKPT.ERROR.MSG}"))
}

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
            from = as.character(min.date),
            to = as.character(max.date),
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

# begin test
index <- sample(seq_len(nrow(metadata)), nrow(sample.metadata))
metadata <- metadata %>% slice(index)
# end test

# Checkpoint: at least as many context samples as our dataset
if (nrow(metadata) < nrow(sample.metadata)) {
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
