#!/usr/bin/env python3

import logging
import json
import pandas as pd
from gb2seq.features import Features


def get_polimorphic_sites(df: pd.DataFrame) -> set:
    return set(df.POS)


def window_calculation(sites: set, window: int, step: int, genome_size: int, coord: str) -> pd.DataFrame:
    ft = Features(coord) # Create Features object to obtain annotations
    positions = []
    pct = []
    genes = []
    lim_sup = genome_size + 1
    for position in range(1, lim_sup, step):
        if len(ft.getFeatureNames(position)) == 0:
            genes.append("Intergenic")
        else:
            genes.append(list(ft.getFeatureNames(position))[0])
        # Add percent (excluding initial and final positions)
        if position - window not in range(1, lim_sup):
            pct.append(0)
        else:
            # Calculate nº of polimorphisms in the window
            num_snp = len([x for x in sites if x in range(position - window, position + 1)])
            pct.append(num_snp / window)
        positions.append(position)
    return pd.DataFrame({"position": positions, "fractions": pct, "gen": genes})


def main():

    logging.basicConfig(filename=snakemake.log[0], format=snakemake.config["LOG_PY_FMT"], level=logging.INFO)

    logging.info("Reading input VCF")
    df = pd.read_table(snakemake.input.vcf)

    logging.info("Reading features")
    with open(snakemake.input.features) as f:
        features_key = json.load(f)

    logging.info("Getting polimorphic sites")
    sites = get_polimorphic_sites(df)

    logging.info("Sliding window calculation of proportion of polimorphic sites")

    frame = window_calculation(sites, snakemake.params.window, snakemake.params.step, 29903, snakemake.input.gb)
    
    logging.info("Saving results")
    frame.replace(features_key, inplace = True)
    frame.to_csv(snakemake.output.window_df, index= False)


if __name__ == "__main__":
    main()
