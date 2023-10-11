#!/usr/bin/env python3

import json
import logging
import pandas as pd
from gb2seq.features import Features
from Bio import SeqIO


def main():
    logging.basicConfig(filename=snakemake.log[0], format=snakemake.config["LOG_PY_FMT"], level=logging.INFO)
    
    logging.info("Reading features")
    ft = Features(snakemake.input.gb)
    with open(snakemake.input.features) as f:
        feature_key = list(json.load(f))

    logging.info("Reading reference")
    reference = str(next(SeqIO.parse(snakemake.input.ref, "fasta")).seq)
    positions = [x for x in range(len(reference))]

    logging.info("Building annotation")
    genes = []
    for pos in positions:
        if len(ft.getFeatureNames(pos)) == 0:
            genes.append("Intergenic")
        else:
            genes.append(feature_key[list(ft.getFeatureNames(pos))[0]])

    logging.info("Writing table")
    df = pd.DataFrame({"POS": positions, "GEN": genes})
    df.to_csv(snakemake.output.df, index= False)


if __name__ == "__main__":
    main()
