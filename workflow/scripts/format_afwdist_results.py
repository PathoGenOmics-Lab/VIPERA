#!/usr/bin/env python3

import logging
import pandas as pd


if __name__ == "__main__":

    logging.basicConfig(filename=snakemake.log[0], format=snakemake.config["LOG_PY_FMT"], level=logging.INFO)

    logging.info("Read pairwise distances")
    df = pd.read_csv(snakemake.input.distances)

    logging.info("Initializing formatted output")
    output = pd.DataFrame(
        data=0.0,
        columns=snakemake.params.samples,
        index=snakemake.params.samples,
        dtype="float64"
    )

    logging.info("Filling table")
    for i, row in df.iterrows():
        output.loc[row.sample_m, row.sample_n] = row.distance
        output.loc[row.sample_n, row.sample_m] = row.distance
    
    logging.info("Writing formatted results")
    output.to_csv(snakemake.output.distances)
