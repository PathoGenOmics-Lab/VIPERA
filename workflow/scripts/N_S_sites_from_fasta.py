#!/usr/bin/env python3

import logging
import json
import copy
import pandas as pd
from gb2seq.alignment import Gb2Alignment
from gb2seq.features import Features
from dark.fasta import FastaReads


def split_into_codons(seq: str) -> list:
    """Split the complete CDS feature in to a list of codons"""
    codons = [seq[i:i + 3] for i in range(0, len(seq), 3) if "N" not in seq[i:i + 3]]
    return codons


def potential_changes_dict(genetic_code: dict) -> dict:
    """Generate a dictionary with S and N pre-calculated for all possible codons (key:codon, value: (S,N)"""
    potential_changes = {   "S": {  "AAA":0.0,"AAC":0.0,"AAG":0.0,"AAT":0.0,"ACA":0.0,"ACC":0.0,"ACG":0.0,"ACT":0.0,"AGA":0.0,"AGC":0.0, \
                                    "AGG":0.0,"AGT":0.0,"ATA":0.0,"ATC":0.0,"ATG":0.0,"ATT":0.0,"CAA":0.0,"CAC":0.0,"CAG":0.0,"CAT":0.0, \
                                    "CCA":0.0,"CCC":0.0,"CCG":0.0,"CCT":0.0,"CGA":0.0,"CGC":0.0,"CGG":0.0,"CGT":0.0,"CTA":0.0,"CTC":0.0,"CTG":0.0, \
                                    "CTT":0.0,"GAA":0.0,"GAC":0.0,"GAG":0.0,"GAT":0.0,"GCA":0.0,"GCC":0.0,"GCG":0.0,"GCT":0.0,"GGA":0.0,"GGC":0.0, \
                                    "GGG":0.0,"GGT":0.0,"GTA":0.0,"GTC":0.0,"GTG":0.0,"GTT":0.0,"TAA":0.0,"TAC":0.0,"TAG":0.0,"TAT":0.0,"TCA":0.0, \
                                    "TCC":0.0,"TCG":0.0,"TCT":0.0,"TGA":0.0,"TGC":0.0,"TGG":0.0,"TGT":0.0,"TTA":0.0,"TTC":0.0,"TTG":0.0,"TTT":0.0},
                            "N": {  "AAA":0.0,"AAC":0.0,"AAG":0.0,"AAT":0.0,"ACA":0.0,"ACC":0.0,"ACG":0.0,"ACT":0.0,"AGA":0.0,"AGC":0.0,"AGG":0.0, \
                                    "AGT":0.0,"ATA":0.0,"ATC":0.0,"ATG":0.0,"ATT":0.0,"CAA":0.0,"CAC":0.0,"CAG":0.0,"CAT":0.0,"CCA":0.0,"CCC":0.0,"CCG":0.0, \
                                    "CCT":0.0,"CGA":0.0,"CGC":0.0,"CGG":0.0,"CGT":0.0,"CTA":0.0,"CTC":0.0,"CTG":0.0,"CTT":0.0,"GAA":0.0,"GAC":0.0,"GAG":0.0, \
                                    "GAT":0.0,"GCA":0.0,"GCC":0.0,"GCG":0.0,"GCT":0.0,"GGA":0.0,"GGC":0.0,"GGG":0.0,"GGT":0.0,"GTA":0.0,"GTC":0.0,"GTG":0.0, \
                                    "GTT":0.0,"TAA":0.0,"TAC":0.0,"TAG":0.0,"TAT":0.0,"TCA":0.0,"TCC":0.0,"TCG":0.0,"TCT":0.0,"TGA":0.0,"TGC":0.0,"TGG":0.0, \
                                    "TGT":0.0,"TTA":0.0,"TTC":0.0,"TTG":0.0,"TTT":0.0}}
    # Mutate (substitutions) all possible codons in the given genetic code, and count proportions of mutations that are synonymous and non-synonmyous
    for codon in genetic_code.keys():
        for codon_p in range(0, 3):
            nts = ["A", "G", "T", "C"]
            # Do not consider self substitutions, e.g. A->A
            nts.remove(codon[codon_p]) 
            # ...and for each nucleotide that the bp can change 
            for nt in nts:
                codon_mutated = list(copy.deepcopy(codon))
                codon_mutated[codon_p] = nt  # mutate the basepair
                codon_mutated = "".join(codon_mutated)
                # ...count how many of them are synonymous.
                if genetic_code[codon] == genetic_code[codon_mutated]:
                    potential_changes["S"][codon]+=1.0/3.0 
                else:
                    potential_changes["N"][codon]+=1.0/3.0
    return potential_changes


def get_feature_codons(alignment: Gb2Alignment, annotation: list) -> dict:
    dct = {key:alignment.ntSequences(key)[1].sequence for key in annotation}
    return {key:split_into_codons(item) for key,item in dct.items()}


def get_df(codons: dict, genetic_code: dict) -> pd.DataFrame:
    keys = []
    N_sites = []
    S_sites = []
    values = potential_changes_dict(genetic_code)
    for key,item in codons.items():
        keys.append(key)
        N = sum([values["N"][x] for x in item if x in values["N"].keys()])
        S = sum([values["S"][x] for x in item if x in values["S"].keys()])
        N_sites.append(N)
        S_sites.append(S) 
    
    return pd.DataFrame({ "gene":keys, "N":N_sites, "S":S_sites })


def main():

    logging.basicConfig(filename=snakemake.log[0], format=snakemake.config["LOG_PY_FMT"], level=logging.INFO)
    
    logging.info("Reading features")
    with open(snakemake.input.features) as f:
        feature_list = list(json.load(f).keys())

    logging.info("Reading genetic code")
    with open(snakemake.input.genetic_code) as f:
        genetic_code = json.load(f)
    
    logging.info("Create alignment object")
    features = Features(snakemake.input.gb)
    seq = list(FastaReads(snakemake.input.fasta))[0]
    aln = Gb2Alignment(seq, features)

    logging.info("Splitting ancestral sequence into codons")
    codons_dict = get_feature_codons(aln, feature_list)

    logging.info("Calculating synonymous and non synonymous sites")
    df = get_df(codons_dict, genetic_code)

    logging.info("Saving results")
    df.to_csv(snakemake.output.csv,index= False)


if __name__ == "__main__":
    main()
