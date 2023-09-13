#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Adapted script from https://github.com/PathoGenOmics-Lab/genetic-distances

import logging
import multiprocessing as mp

import pandas as pd
from Bio import SeqIO


def parse_vcf() -> tuple:
    """
    Parse a VCF containing positions for masking. Assumes the VCF file is
    formatted as:
    github.com/W-L/ProblematicSites_SARS-CoV2/blob/master/problematic_sites_sarsCov2.vcf
    with a "mask" or "caution" recommendation in column 7.
    Masked sites are specified with params.
    """
    vcf = pd.read_csv(
        snakemake.input["vcf"],
        delim_whitespace=True,
        comment="#",
        names=("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")
    )
    positions = tuple(vcf.loc[ vcf.FILTER.isin(snakemake.params.mask_class), "POS" ])

    return positions


def create_freq_dict(input_file:str) -> dict:

    alt_values = set()  # Set to avoid duplicates

    df = pd.read_table(input_file, sep = "\t") # Get all possible alternative alleles
    alt_values.update(df['ALT'].unique())
    
    freq = { alt_value: 0 for alt_value in alt_values } # Initialize a dict for allele frequencies

    return freq


def tsv_from_seq(tsv_reference:str ,reference:str , reference_name:str) -> pd.DataFrame:

    mask = parse_vcf()
    
    pos = []
    alt = []
    for i in range(1,len(tsv_reference) + 1):
        if i not in mask and tsv_reference[i -1] != reference[i-1]:
            pos.append(i)
            alt.append(reference[i -1])

    df = pd.DataFrame({"POS":pos,"ALT": alt})

    df["ALT_FREQ"] = 1  # As is a reference genome, it is assumed for all positions to be monomorphic
    df["REGION"] = reference_name

    return df


def get_pos_tup(df:pd.DataFrame, pos:int, reference:str, freq:dict) -> tuple:

    freq = freq.copy()  
    alt_keys = sorted(list(freq.keys()))  

    # If studied position has polimorphisims, allele frequencies are captured
    if pos + 1 in df["POS"].values:

        df_ = df[df["POS"] == pos+1]
        for base in alt_keys:
            if base in df_["ALT"].values:
                freq[base] = float(df_["ALT_FREQ"][df_["ALT"] == base].iloc[0])

    
    ref = 1 - sum(freq.values()) # Obtain frequency for reference allele
    freq[reference[pos]] += ref
    
    return tuple(freq[base] for base in alt_keys)


def calc_heterozygosities(df1:pd.DataFrame, df2:pd.DataFrame, pos:int, reference:str, freq:dict):

    freqs1 = get_pos_tup(df1, pos, reference, freq)
    freqs2 = get_pos_tup(df2, pos, reference, freq)

    hs1 = heterozygosity(freqs1)
    hs2 = heterozygosity(freqs2)
    hs = (hs1 + hs2) / 2

    total_freqs = [ (f1 + f2) / 2 for f1, f2 in zip(freqs1, freqs2) ]
    ht = heterozygosity(total_freqs)

    return hs, ht


def heterozygosity(freqs:tuple) -> float:
    return 1 - sum([ f ** 2 for f in freqs ])


def calc_fst_weir_cockerham(hs:float, ht:float) -> float:
    return (ht - hs) / ht if ht != 0 else 0


def get_dif_n(df:pd.DataFrame, COV1:str, COV2:str, reference:str, freq:dict) -> float:

    positions = df["POS"].unique().tolist()
    if len(positions) == 0:
        return 0
    
    df1 = df[df["REGION"] == COV1]
    df2 = df[df["REGION"] == COV2]

    return sum([calc_fst_weir_cockerham(*calc_heterozygosities(df1, df2, i-1, reference, freq)) 
                for i in positions])


def _calculate_distance(df:pd.DataFrame, sample:str,reference:str, freq:dict, cov_list:list) -> list:
    return [get_dif_n(df, sample, cov, reference, freq) for cov in cov_list]


def get_matrix(df:pd.DataFrame, cov_list:list, reference:str, freq:dict, num_jobs:int) -> pd.DataFrame:

    distance_matrix = {}

    with mp.Pool(num_jobs) as pool:

        results = pool.starmap(
            _calculate_distance,
            [ (df, sample, reference, freq, cov_list) for sample in cov_list ]
        )

    for i, sample in enumerate(cov_list):
        distance_matrix[sample] = results[i]
    
    for i in range(len(cov_list)):
        for j in range(i+1, len(cov_list)):
            distance_matrix[cov_list[j]][i] = distance_matrix[cov_list[i]][j]
    
    return pd.DataFrame(distance_matrix, index=cov_list)


def read_and_concatenate_tsvs(input:str, tsv_reference:str, reference:str, reference_name:str) -> pd.DataFrame:

    df_1 = pd.read_table(input, sep = "\t")

    return pd.concat([ df_1, tsv_from_seq(tsv_reference,reference,reference_name) ], ignore_index=True)


def main():

    logging.basicConfig(filename=snakemake.log[0], format=snakemake.config["LOG_PY_FMT"], level=logging.INFO)

    logging.info("Reading input FASTA files")
    reference = str(next(SeqIO.parse(snakemake.params.tsv_reference, "fasta")).seq)
    outgroup = str(next(SeqIO.parse(snakemake.params.reference, "fasta")).seq)
    outgroup_name = str(next(SeqIO.parse(snakemake.params.reference, "fasta")).id)

    logging.info("Reading input tables")
    df = read_and_concatenate_tsvs(snakemake.input.tsv, reference, outgroup, outgroup_name)
    cov_list = df["REGION"].unique().tolist()

    logging.info(f"Parallelizing the calculation with {snakemake.threads} jobs")
    freq = create_freq_dict(snakemake.input.tsv)
    df = get_matrix(df, cov_list, reference, freq, snakemake.threads)

    logging.info("Writing results")
    df.to_csv(snakemake.output.distances)


if __name__ == "__main__":
    main()
