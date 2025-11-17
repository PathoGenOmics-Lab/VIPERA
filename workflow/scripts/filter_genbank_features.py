#!/usr/bin/env python3

import logging
from typing import Dict, Iterable

from Bio import SeqIO
from Bio.SeqFeature import SeqFeature


def iter_features_filtering(features: Iterable[SeqFeature], included: Dict[str, str], excluded: Dict[str, str]) -> Iterable[SeqFeature]:
    # No filters
    if len(included) == 0 and len(excluded) == 0:
        logging.debug("Selecting all features")
        return iter(features)
    # Only inclusion filter
    elif len(included) == 0 and len(excluded) != 0:
        logging.debug(f"Selecting features excluding all of {excluded}")
        return (
            feature for feature in features
            if all(
                (qualifier_value not in excluded.get(qualifier_key, []))
                for qualifier_key in excluded.keys()
                for qualifier_value in feature.qualifiers.get(qualifier_key, [])
            )
        )
    # Only exclusion filter
    elif len(included) != 0 and len(excluded) == 0:
        logging.debug(f"Selecting features including any of {included}")
        return (
            feature for feature in features
            if any(
                (qualifier_value in included.get(qualifier_key, []))
                for qualifier_key in included.keys()
                for qualifier_value in feature.qualifiers.get(qualifier_key, [])
            )
        )
    # Inclusion then exclusion filter
    else:
        logging.debug(f"Selecting features including any of {included} and then excluding all of {excluded}")
        included_features = (
            feature for feature in features
            if any(
                (qualifier_value in included.get(qualifier_key, []))
                for qualifier_key in included.keys()
                for qualifier_value in feature.qualifiers.get(qualifier_key, [])
            )
        )
        return (
            feature for feature in included_features
            if all(
                (qualifier_value not in excluded.get(qualifier_key, []))
                for qualifier_key in excluded.keys()
                for qualifier_value in feature.qualifiers.get(qualifier_key, [])
            )
        )


if __name__ == "__main__":
    logging.basicConfig(
        filename=snakemake.log[0], format=snakemake.config["LOG_PY_FMT"],
        level=logging.INFO
    )

    logging.info("Reading GenBank file")
    gb = SeqIO.read(snakemake.input.gb, format="gb")

    logging.info("Extracting CDS")
    features = []
    for feature in iter_features_filtering(gb.features, snakemake.params.included, snakemake.params.excluded):
        logging.debug(
            "Processing SeqFeature: "
            f"ID={feature.id} type={feature.type} location={feature.location} "
            f"gene={feature.qualifiers.get('gene', [])} "
            f"locus_tag={feature.qualifiers.get('locus_tag', [])} "
            f"product={feature.qualifiers.get('product', [])}"
        )
        features.append(feature)
    
    logging.info("Replacing features")
    gb.features = features

    logging.info("Writing filtered GenBank files")
    SeqIO.write(gb, snakemake.output.gb, "gb")
