#!/usr/bin/env python3

import logging
import pandas as pd
from gb2seq.features import Features
from Bio import SeqIO

replace_terry = {
    "ORF1ab polyprotein": "orf1ab",
    "ORF1a polyprotein": "orf1ab",
    "surface glycoprotein": "S",
    "ORF3a protein": "ORF3a",
    "envelope protein": "E",
    "membrane glycoprotein": "M",
    "ORF6 protein": "ORF6",
    "ORF7a protein": "ORF7",
    "ORF7b": "ORF7",
    "ORF8 protein": "ORF8",
    "nucleocapsid phosphoprotein": "N",
    "ORF10 protein": "ORF10"
}

def main():
    logging.basicConfig(filename=snakemake.log[0], format=snakemake.config["LOG_PY_FMT"], level=logging.INFO)
    ft = Features(snakemake.input.gb)

    reference = str(next(SeqIO.parse(snakemake.input.ref, "fasta")).seq)

    positions = [x for x in range(len(reference))]

    genes = []
    for pos in positions:

        if len(ft.getFeatureNames(pos)) == 0:
            genes.append("Intergenic")
        else:
            genes.append(replace_terry[list(ft.getFeatureNames(pos))[0]])

    df = pd.DataFrame({"POS":positions, "GEN": genes})

    df.to_csv(snakemake.output.df,index= False)


if __name__ == "__main__":
    main()

    