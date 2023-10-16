#!/usr/bin/env python3

import logging
import json
import pandas as pd
from gb2seq.features import Features


def window_calculation(sites: set, window: int, step: int, coord: str) -> pd.DataFrame:
    ft = Features(coord) # Create Features object to obtain annotations
    positions = []
    pct = []
    genes = []
    lim_sup = len(ft.reference.sequence) + 1
    for position in range(1, lim_sup, step):
        if len(ft.getFeatureNames(position)) == 0:
            genes.append("Intergenic")
        else:
            genes.append(list(ft.getFeatureNames(position))[0])
        # Add percent (excluding initial and final positions)
        if position - window not in range(1, lim_sup):
            pct.append(0)
        else:
            # Calculate no. of polimorphisms in the window
            num_snp = len([x for x in sites if x in range(position - window, position + 1)])
            pct.append(num_snp / window)
        positions.append(position)
    return pd.DataFrame({"position": positions, "fractions": pct, "gen": genes})


def main():

    logging.basicConfig(filename=snakemake.log[0], format=snakemake.config["LOG_PY_FMT"], level=logging.INFO)

    logging.info("Reading input VCF")
    df = pd.read_table(snakemake.input.vcf)
    sites = set(df.POS)

    logging.info("Reading features")
    with open(snakemake.input.features) as f:
        features_key = json.load(f)

    logging.info("Sliding window calculation of proportion of polimorphic sites")
    frame = window_calculation(sites, snakemake.params.window, snakemake.params.step, snakemake.input.gb)
    
    logging.info("Saving results")
    frame.replace(features_key, inplace = True)
    frame.to_csv(snakemake.output.window_df, index= False)


if __name__ == "__main__":
    main()
