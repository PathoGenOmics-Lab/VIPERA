#!/usr/bin/env python3

import logging
import pandas as pd
from typing import NewType, Dict, Iterable, List

from Bio import SeqIO
from Bio.SeqFeature import SeqFeature
from Bio.SeqFeature import CompoundLocation, FeatureLocation


FeatureIndex = NewType("FeatureIndex", Dict[str, FeatureLocation | CompoundLocation])


def index_features(features: Iterable[SeqFeature]) -> FeatureIndex:
    index = FeatureIndex({})
    for feature in features:
        identifier = "|".join(feature.qualifiers.get(snakemake.params.gb_qualifier_display, []))
        if identifier == "":
            logging.error(f"Feature at {feature.location} has no qualifier '{snakemake.params.gb_qualifier_display}' and will be skipped")
        elif feature.location is None:
            logging.warning(f"Feature '{identifier}' has no location and will be skipped")
        elif identifier in index:
            logging.warning(f"Identifier '{identifier}' is already in feature index and will not be replaced by feature at {feature.location}")
        else:
            index[identifier] = feature.location
    return index


def feature_names_at(position: int, feature_index: FeatureIndex) -> List[str]:
    return [name for name, location in feature_index.items() if position in location]


def window_calculation(sites: set, window: int, step: int, size: int, feature_index: FeatureIndex) -> pd.DataFrame:
    positions, fractions, features = [], [], []
    lim_sup = size + 1
    for position in range(1, lim_sup, step):
        feature_names = ";".join(feature_names_at(position, feature_index))
        if len(feature_names) == 0:
            features.append("Intergenic")
        else:
            # Include all features on site
            features.append(feature_names)
        # Add percent (excluding initial and final positions)
        if position - window not in range(1, lim_sup):
            fractions.append(0.0)
        else:
            # Calculate number of polymorphisms in the window
            num_snp = sum(
                1 for x in sites if x in range(position - window, position + 1)
            )
            fractions.append(num_snp / window)
        positions.append(position)
    return pd.DataFrame({"position": positions, "fraction": fractions, "feature": features})


def main():

    logging.basicConfig(
        filename=snakemake.log[0], format=snakemake.config["LOG_PY_FMT"],
        level=logging.INFO
    )

    # Read input files
    logging.info("Reading input variants file")
    df = pd.read_table(snakemake.input.variants)
    sites = set(df.POS)

    logging.info("Reading GenBank file")
    gb = SeqIO.read(snakemake.input.gb, format="gb")

    # Build feature iterator given the provided filters
    if len(snakemake.params.features) == 0:
        logging.debug("Selecting all features")
        feature_iterator = iter(gb.features)
    else:
        included = snakemake.params.features.get("include", {})
        excluded = snakemake.params.features.get("exclude", {})
        logging.debug(f"Selecting features including any of {included} and excluding all of {excluded}")
        feature_iterator = (
            feature for feature in gb.features
            if any(
                (qualifier_value in included.get(qualifier_key, []))
                for qualifier_key in included.keys()
                for qualifier_value in feature.qualifiers.get(qualifier_key, [])
            ) and all(
                (qualifier_value not in excluded.get(qualifier_key, []))
                for qualifier_key in excluded.keys()
                for qualifier_value in feature.qualifiers.get(qualifier_key, [])
            )
        )

    logging.info("Indexing feature locations")
    feature_index = index_features(feature_iterator)

    logging.info("Calculating polimorphic sites sliding window")
    windows = window_calculation(
        sites,
        snakemake.params.window,
        snakemake.params.step,
        len(gb),
        feature_index
    )

    logging.info("Saving results")
    windows.to_csv(snakemake.output.window_df, index=False)


if __name__ == "__main__":
    main()
