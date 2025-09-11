#!/usr/bin/env python3

import logging
import json
import pandas as pd
from gb2seq.features import Features


def window_calculation(sites: set, window: int, step: int, gb_features: Features) -> pd.DataFrame:
    positions, fractions, features = [], [], []
    lim_sup = len(gb_features.reference.sequence) + 1
    for position in range(1, lim_sup, step):
        if len(gb_features.getFeatureNames(position)) == 0:
            features.append("Intergenic")
        else:
            # Include all features on site
            features.append(gb_features.getFeatureNames(position))
        # Add percent (excluding initial and final positions)
        if position - window not in range(1, lim_sup):
            fractions.append(0.0)
        else:
            # Calculate no. of polimorphisms in the window
            num_snp = sum(
                1 for x in sites if x in range(position - window, position + 1)
            )
            fractions.append(num_snp / window)
        positions.append(position)
    return pd.DataFrame({"position": positions, "fraction": fractions, "feature": features})


def select_and_format_features(window_features: set) -> str | None:
    selected_features = sorted(
        snakemake.params.select_gb_features[feat]
        for feat in window_features
        if feat in snakemake.params.select_gb_features
    )
    if len(selected_features) != 0:
        return "|".join(selected_features)
    else:
        return None


def format_features(window_features: set) -> str:
    if len(window_features) != 0:
        return "|".join(sorted(window_features))
    else:
        return None


def main():

    logging.basicConfig(
        filename=snakemake.log[0], format=snakemake.config["LOG_PY_FMT"],
        level=logging.INFO
    )

    logging.info("Reading input VCF")
    df = pd.read_table(snakemake.input.vcf)
    sites = set(df.POS)

    logging.info("Reading genbank features")
    features = Features(snakemake.input.gb)

    logging.info("Calculating polimorphic sites sliding window")
    windows = window_calculation(
        sites,
        snakemake.params.window,
        snakemake.params.step,
        features
    )

    if len(snakemake.params.select_gb_features) != 0:
        logging.info("Selecting and formatting genbank features")
        windows.feature = windows.feature.map(select_and_format_features)
    else:
        logging.info("Formatting genbank features")
        windows.feature = windows.feature.map(format_features)

    logging.info("Saving results")
    windows.to_csv(snakemake.output.window_df, index=False)


if __name__ == "__main__":
    main()
