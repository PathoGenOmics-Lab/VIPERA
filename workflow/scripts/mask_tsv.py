#!/usr/bin/env python3

import logging
import pandas as pd
 
def parse_vcf():
    """
    Parse a VCF containing positions for masking. Assumes the VCF file is
    formatted as:
    github.com/W-L/ProblematicSites_SARS-CoV2/blob/master/problematic_sites_sarsCov2.vcf
    with a "mask" or "caution" recommendation in column 7.
    Masked sites are specified with params.
    """
    vcf = pd.read_csv(
        snakemake.config["PROBLEMATIC_VCF"],
        delim_whitespace=True,
        comment="#",
        names=("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")
    )
    positions = tuple(vcf.loc[vcf.FILTER.isin(snakemake.params.mask_class), "POS"])
    return positions

def main():
    logging.basicConfig(filename=snakemake.log[0], format=snakemake.config["LOG_PY_FMT"], level=logging.INFO)
    
    # parse VCF to get list of sites to mask
    iffy_sites = parse_vcf()

    # open unmasked VCF 
    f = open(snakemake.input.tsv,"r")

    # open masked vcf
    mv = open(snakemake.output.masked_tsv,"w")

    for line in f:
        
        if line[0] == "R":
            mv.write(line)
            continue
        elif int(line.split("\t")[1]) not in iffy_sites:

            mv.write(line)

    mv.close()
    f.close()




if __name__ == "__main__":
    main()
