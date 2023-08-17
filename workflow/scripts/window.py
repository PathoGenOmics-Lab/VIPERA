#!/usr/bin/env python3



import logging
import pandas as pd
from gb2seq.features import Features


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


def get_polimorphic_sites(df):
    return set(df.POS)


def window_calculation(sites, step, genome_size, coord):

    ft = Features(coord)
    
    positions = []
    pct = []
    genes = []
    lim_sup = genome_size + 1
    for position in range(1, lim_sup):
        # Add gene
        if len(ft.getFeatureNames(position)) == 0:
            genes.append("Intergenic")
        else:
            genes.append(list(ft.getFeatureNames(position))[0])
    
        # Add percent (excluding initial and final positions)
        if position - step not in range(1, lim_sup):
            pct.append(0)
        else:
            # Calculate SNPs
            num_snp = 0
            for x in sites:
                if x in range(position - step, position + 1):
                    num_snp += 1
            pct.append(num_snp / step)
        # Add positions
        positions.append(position)
    return pd.DataFrame({"position": positions, "fractions": pct, "gen": genes})


def main():
    logging.basicConfig(filename=snakemake.log[0], format=snakemake.config["LOG_PY_FMT"], level=logging.INFO)

    # Process
    df = pd.read_table(snakemake.input.vcf)    
    sites = get_polimorphic_sites(df)
    frame = window_calculation(sites, snakemake.params.step, 29903, snakemake.input.gb)
    frame.replace(replace_terry, inplace = True)
    frame.to_csv(snakemake.output.window_df,index= False)


if __name__ == "__main__":
    main()
