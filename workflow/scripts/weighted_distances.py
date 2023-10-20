#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Adapted script from https://github.com/PathoGenOmics-Lab/genetic-distances

import logging
from typing import List, Tuple

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


def read_monofasta(path: str) -> SeqRecord:
    fasta = SeqIO.parse(path, "fasta")
    record = next(fasta)
    if next(fasta, None) is not None:
        logging.warning(f"There are unread records left in '{path}'")
    return record


def read_masked_sites(vcf_path: str, mask_classes: List[str]) -> List[int]:
    """
    Parse a VCF containing positions for masking. Assumes the VCF file is
    formatted as in:
    github.com/W-L/ProblematicSites_SARS-CoV2/blob/master/problematic_sites_sarsCov2.vcf
    with a "mask" or "caution" recommendation in column 7.
    Masked sites are specified with params.
    """
    vcf = pd.read_csv(
        vcf_path,
        delim_whitespace=True,
        comment="#",
        names=("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")
    )
    return vcf.loc[vcf.FILTER.isin(mask_classes), "POS"].tolist()


def build_ancestor_variant_table(ancestor: Seq, reference: Seq, reference_name: str, masked_positions: List[int]) -> pd.DataFrame:
    pos = []
    alt = []
    for i in range(1, len(ancestor) + 1):
        if i not in masked_positions and ancestor[i-1] != reference[i-1]:
            pos.append(i)
            alt.append(reference[i-1])
    df = pd.DataFrame({"POS": pos, "ALT": alt})
    df["ALT_FREQ"] = 1  # As a reference genome, we assume all positions are monomorphic
    df["REGION"] = reference_name
    return df


def get_frequencies_in_position(variant_table: pd.DataFrame, sample_name: str, position: int, reference: Seq) -> Tuple[float]:
    frequencies = {alt: 0. for alt in variant_table["ALT"].unique()}
    alt_keys = sorted(frequencies.keys())
    sample_variant_table = variant_table[variant_table["REGION"] == sample_name]
    # If the position has polimorphisims, allele frequencies are captured
    if position in sample_variant_table["POS"].values:
        alleles = sample_variant_table[sample_variant_table["POS"] == position]
        for base in alt_keys:
            if base in alleles["ALT"].values:
                frequencies[base] = float(alleles["ALT_FREQ"][alleles["ALT"] == base].iloc[0])
    # Obtain frequency for reference allele
    reference_allele = reference[position-1]
    reference_freq = 1 - sum(frequencies.values())
    if reference_allele in frequencies:
        frequencies[reference[position-1]] += reference_freq
    else:
        frequencies[reference[position-1]] = reference_freq
    return tuple(frequencies[base] for base in alt_keys)


def heterozygosity(freqs: Tuple[float]) -> float:
    return 1 - sum([f ** 2 for f in freqs])


def calc_fst_weir_cockerham(hs: float, ht: float) -> float:
    return (ht - hs) / ht if ht != 0 else 0


def build_cache(variant_table: pd.DataFrame, reference: Seq):
    cache = {"freq": {}, "hz": {}}
    for sample_name in variant_table["REGION"].unique():
        for position in variant_table["POS"].unique().astype("Int64"):
            if sample_name not in cache["freq"]:
                cache["freq"][sample_name] = {}
            if sample_name not in cache["hz"]:
                cache["hz"][sample_name] = {}
            cache["freq"][sample_name][position] = get_frequencies_in_position(variant_table, sample_name, position, reference)
            cache["hz"][sample_name][position] = heterozygosity(cache["freq"][sample_name][position])
    return cache


def calc_heterozygosities(sample1_name: str, sample2_name: str, pos: int, cache: dict):
    # Retrieve pre-computed values
    freqs1 = cache["freq"][sample1_name][pos]
    freqs2 = cache["freq"][sample2_name][pos]
    hs1 = cache["hz"][sample1_name][pos]
    hs2 = cache["hz"][sample2_name][pos]
    # Calculate pairwise values
    total_freqs = [(f1 + f2) / 2 for f1, f2 in zip(freqs1, freqs2)]
    ht = heterozygosity(total_freqs)
    hs = (hs1 + hs2) / 2
    return hs, ht


def calculate_pairwise_distance(positions: List[int], sample1_name: str, sample2_name: str, cache: dict) -> float:
    if len(positions) == 0:
        return 0
    else:
        return sum([calc_fst_weir_cockerham(*calc_heterozygosities(sample1_name, sample2_name, pos, cache)) for pos in positions])


def calculate_sample_distances(positions: List[int], sample_name: str, samples: List[str], cache: dict) -> List[float]:
    return [calculate_pairwise_distance(positions, sample_name, other_sample_name, cache) for other_sample_name in samples]


def calculate_distance_matrix(variant_table: pd.DataFrame, samples: List[str], reference: Seq) -> pd.DataFrame:
    positions = variant_table["POS"].astype("Int64").unique().tolist()
    cache = build_cache(variant_table, reference)
    distance_matrix = {}
    for sample_name in samples:
        distance_matrix[sample_name] = calculate_sample_distances(positions, sample_name, samples, cache)
    for i in range(len(samples)):
        for j in range(i+1, len(samples)):
            distance_matrix[samples[j]][i] = distance_matrix[samples[i]][j]
    return pd.DataFrame(distance_matrix, index=samples)


def main():

    logging.basicConfig(filename=snakemake.log[0], format=snakemake.config["LOG_PY_FMT"], level=logging.INFO)

    logging.info("Reading input FASTA files")
    ancestor = read_monofasta(snakemake.input.ancestor)
    reference = read_monofasta(snakemake.input.reference)

    logging.info("Reading input tables")
    masked_positions = read_masked_sites(snakemake.input.vcf, snakemake.params.mask_class)
    input_table = pd.read_table(snakemake.input.tsv, sep="\t")
    ancestor_table = build_ancestor_variant_table(ancestor.seq, reference.seq, reference.id, masked_positions)
    variant_table = pd.concat([input_table, ancestor_table], ignore_index=True)

    logging.info(f"Calculating distance matrix")
    sample_names = snakemake.params.samples + [reference.id]
    distances = calculate_distance_matrix(variant_table, sample_names, ancestor.seq)

    logging.info("Writing results")
    distances.to_csv(snakemake.output.distances)


if __name__ == "__main__":
    main()
