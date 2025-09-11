#!/usr/bin/env python3

import logging
import json
import itertools as it
import pandas as pd
from gb2seq.alignment import Gb2Alignment
from gb2seq.features import Features
from dark.fasta import FastaReads


def split_into_codons(seq: str) -> list:
    """Split the complete CDS feature in to a list of codons"""
    codons = [
        seq[i:i + 3] for i in range(0, len(seq), 3) if "N" not in seq[i:i + 3]
    ]
    return codons


def calculate_potential_changes(genetic_code: dict) -> dict:
    """Generate a dictionary with S and N pre-calculated for all possible codons"""
    # Initialize structure
    nts = set(["A", "G", "T", "C"])
    potential_changes = {"S": {}, "N": {}}
    for codon in it.product(nts, repeat=3):
        potential_changes["S"]["".join(codon)] = 0.
        potential_changes["N"]["".join(codon)] = 0.
    # Mutate (substitutions) all possible codons in the given genetic code
    # and count proportions of mutations that are synonymous and non-synonmyous
    for codon in genetic_code.keys():
        for codon_p in range(0, 3):
            nts = ["A", "G", "T", "C"]
            # Do not consider self substitutions, e.g. A->A
            nts.remove(codon[codon_p])
            for nt in nts:
                codon_mutated = list(codon)
                # Mutate the basepair
                codon_mutated[codon_p] = nt
                # Count how many of them are synonymous
                if genetic_code[codon] == genetic_code["".join(codon_mutated)]:
                    potential_changes["S"][codon] += 1/3.
                else:
                    potential_changes["N"][codon] += 1/3.
    return potential_changes


def get_feature_codons(alignment: Gb2Alignment, annotation: list) -> dict:
    dct = {key: alignment.ntSequences(key)[1].sequence for key in annotation}
    return {key: split_into_codons(item) for key, item in dct.items()}


def calculate_ns_sites(codons: dict, genetic_code: dict) -> pd.DataFrame:
    features = []
    N_sites = []
    S_sites = []
    values = calculate_potential_changes(genetic_code)
    for key, item in codons.items():
        features.append(key)
        N = sum([values["N"][x] for x in item if x in values["N"].keys()])
        S = sum([values["S"][x] for x in item if x in values["S"].keys()])
        N_sites.append(N)
        S_sites.append(S)
    return pd.DataFrame({"gene": features, "N": N_sites, "S": S_sites})


def main():

    logging.basicConfig(
        filename=snakemake.log[0], format=snakemake.config["LOG_PY_FMT"],
        level=logging.INFO
    )

    logging.info("Reading features")
    feature_list = list(snakemake.params.gb_features.keys())

    logging.info("Reading genetic code")
    with open(snakemake.input.genetic_code) as f:
        genetic_code = json.load(f)

    logging.info("Create alignment object")
    features = Features(snakemake.input.gb)
    fasta_reads = list(FastaReads(snakemake.input.fasta))
    if len(fasta_reads) > 1:
        logging.warning(
            f"More than one record found in {snakemake.input.fasta}, selecting the first one")
    seq = fasta_reads[0]
    aln = Gb2Alignment(seq, features)

    logging.info("Splitting input sequence into codons")
    codons = get_feature_codons(aln, feature_list)

    logging.info("Calculating synonymous and non synonymous sites")
    df = calculate_ns_sites(codons, genetic_code)

    logging.info("Saving results")
    df.to_csv(snakemake.output.csv, index=False)


if __name__ == "__main__":
    main()
