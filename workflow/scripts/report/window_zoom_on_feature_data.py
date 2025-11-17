#!/usr/bin/env python3

import sys
import logging
import json

import pandas as pd


if __name__ == "__main__":

    logging.basicConfig(filename=snakemake.log[0], format=snakemake.config["LOG_PY_FMT"], level=logging.INFO)

    logging.info("Reading regions")
    with open(snakemake.input.regions) as f:
        regions = json.load(f)

    if snakemake.wildcards.region_name not in regions:
        logging.error(f"Region {snakemake.wildcards.region_name} is absent in {snakemake.input.regions}")
        sys.exit(1)
    
    start, end = regions[snakemake.wildcards.region_name]

    logging.info("Reading input table")
    df = pd.read_csv(snakemake.input.table)

    logging.info("Filtering sites and writing results")
    df[(df.position >= start) & (df.position <= end)].to_csv(snakemake.output.table, index=False)
