#!/usr/bin/env python3

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
    # parse VCF to get list of sites to mask
    iffy_sites = parse_vcf()

    # open unmasked VCF 
    f = open(snakemake.input.vcf,"r")

    # open masked vcf
    mv = open(snakemake.output.masked_vcf,"w")

    for line in f:
        
        if line[0] == "#":
            mv.write(line)
            continue
        elif line.split("\t")[1] not in iffy_sites:
            mv.write(line)

    mv.close()
    f.close()




if __name__ == "__main__":
    main()
