#!/usr/bin/env python3

import logging
import pandas as pd
 

def read_problematic_positions(path: str, mask_classes: list) -> tuple:
    """
    Parse a VCF containing positions for masking. Assumes the VCF file is
    formatted as:
    github.com/W-L/ProblematicSites_SARS-CoV2/blob/master/problematic_sites_sarsCov2.vcf
    with a "mask" or "caution" recommendation in column 7.
    Masked sites are specified with params.
    """
    vcf = pd.read_csv(
        path,
        sep="\s+",
        comment="#",
        names=("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")
    )
    return tuple(vcf.loc[vcf.FILTER.isin(mask_classes), "POS"])


def main():
    logging.basicConfig(filename=snakemake.log[0], format=snakemake.config["LOG_PY_FMT"], level=logging.INFO)
    
    logging.info("Getting sites to mask")
    positions = read_problematic_positions(snakemake.input.vcf, snakemake.params.mask_class)

    # Open unmasked VCF and masked VCF
    f = open(snakemake.input.tsv,"r") 
    mv = open(snakemake.output.masked_tsv,"w")
    logging.info("Masking VCF")
    for line in f:
        if line[0] == "R":
            mv.write(line)
            continue
        elif int(line.split("\t")[1]) not in positions:
            mv.write(line)
    mv.close()
    f.close()

    
if __name__ == "__main__":
    main()
