#!/usr/bin/env python3

import logging
import json
import itertools as it
from typing import Dict

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


NTS = ("A", "C", "G", "T")


def read_monofasta(path: str) -> SeqRecord:
    return SeqIO.read(path, format="fasta")


def split_into_codons(seq: Seq) -> list:
    """Split the complete CDS feature in to a list of codons"""
    return [
        seq[i:i + 3] for i in range(0, len(seq)-2, 3) if all(char in NTS for char in seq[i:i + 3])
    ]


def calculate_potential_changes(genetic_code: dict) -> dict:
    """Generate a dictionary with S and N pre-calculated for all possible codons"""
    # Initialize structure
    potential_changes = {"S": {}, "N": {}}
    for codon in it.product(NTS, repeat=3):
        potential_changes["S"]["".join(codon)] = 0.0
        potential_changes["N"]["".join(codon)] = 0.0
    # Mutate (substitutions) all possible codons in the given genetic code
    # and count proportions of mutations that are synonymous and non-synonmyous
    for codon in genetic_code.keys():
        for codon_p in range(0, 3):
            nts = list(NTS)
            # Do not consider self substitutions, e.g. A->A
            nts.remove(codon[codon_p])
            for nt in nts:
                codon_mutated = list(codon)
                # Mutate the basepair
                codon_mutated[codon_p] = nt
                # Count how many of them are synonymous
                if genetic_code[codon] == genetic_code["".join(codon_mutated)]:
                    potential_changes["S"][codon] += 1 / 3.0
                else:
                    potential_changes["N"][codon] += 1 / 3.0
    return potential_changes


def get_feature_codons(coding_records: Dict[str, SeqRecord]) -> dict:
    return {name: split_into_codons(record.seq) for name, record in coding_records.items()}


def calculate_ns_sites(codons: dict, genetic_code: dict) -> pd.DataFrame:
    features = []
    N_sites, S_sites = [], []
    values = calculate_potential_changes(genetic_code)
    for key, item in codons.items():
        features.append(key)
        N = sum([values["N"][x] for x in item if x in values["N"].keys()])
        S = sum([values["S"][x] for x in item if x in values["S"].keys()])
        N_sites.append(N)
        S_sites.append(S)
    return pd.DataFrame({"feature": features, "N": N_sites, "S": S_sites})


def main():

    logging.basicConfig(
        filename=snakemake.log[0], format=snakemake.config["LOG_PY_FMT"],
        level=logging.INFO
    )

    logging.info("Reading genetic code")
    with open(snakemake.input.genetic_code) as f:
        genetic_code = json.load(f)

    logging.info("Reading GenBank file")
    gb = SeqIO.read(snakemake.input.gb, format="gb")

    logging.info("Reading input FASTA")
    record = read_monofasta(snakemake.input.fasta)

    if len(snakemake.params.features) == 0:
        logging.debug("Selecting all features")
        feature_iterator = iter(gb.features)
    else:
        included = snakemake.params.features.get("INCLUDE", {})
        excluded = snakemake.params.features.get("EXCLUDE", {})
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

    logging.info("Extracting CDS")
    coding_records = {}
    for feature in feature_iterator:
        logging.debug(
            "Processing SeqFeature: "
            f"ID={feature.id} type={feature.type} location={feature.location} "
            f"gene={feature.qualifiers.get('gene', [])} "
            f"locus_tag={feature.qualifiers.get('locus_tag', [])} "
            f"product={feature.qualifiers.get('product', [])}"
        )
        identifier = "|".join(feature.qualifiers.get(snakemake.params.gb_qualifier_display, []))
        if identifier == "":
            logging.error(f"Feature at {feature.location} has no qualifier '{snakemake.params.gb_qualifier_display}'")
        elif identifier in coding_records:
            logging.warning(f"Identifier '{identifier}' is already among coding records and will not be replaced by feature at {feature.location}")
        else:
            coding_records[identifier] = feature.extract(record)

    logging.info(f"Splitting {len(coding_records)} records into codons")
    codons = get_feature_codons(coding_records)

    logging.info("Calculating synonymous and non synonymous sites")
    df = calculate_ns_sites(codons, genetic_code)

    logging.info("Saving results")
    df.to_csv(snakemake.output.csv, index=False)


if __name__ == "__main__":
    main()
