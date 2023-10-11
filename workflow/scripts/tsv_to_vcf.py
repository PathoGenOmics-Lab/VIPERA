#!/usr/bin/env python3

import pandas as pd
import logging

def tsv_to_vcf(tsv_file, vcf_file):
    # Read the TSV file
    tsv_df = pd.read_csv(tsv_file, sep='\t')

    # Open a new VCF file for writing
    with open(vcf_file, 'w') as vcf:
        # Write the VCF header
        vcf.write('##fileformat=VCFv4.2\n')
        vcf.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')

        # Process each row in the TSV file
        for index, row in tsv_df.iterrows():
            # Extract fields from the TSV row
            chrom = 'NC_045512.2'
            pos = row['POS']
            ref = row['REF']
            alt = row['ALT']
            # Add other necessary fields for the VCF file
            
            # Handle INDELs
            if alt[0] == "+":
                alt = ref + alt[1:]
            elif alt[0] == "-":
                alt2 = ref
                ref += alt[1:]
                alt = alt2  # Adjust this line to use the correct value for ALT after a deletion

            # Write the VCF row
            vcf.write(f'{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\n')

def main():
    logging.basicConfig(filename=snakemake.log[0], format=snakemake.config["LOG_PY_FMT"], level=logging.INFO)

    logging.info("Converting from TSV to VCF format")
    input_tsv_file = snakemake.input.tsv
    output_vcf_file = snakemake.output.vcf
    tsv_to_vcf(input_tsv_file, output_vcf_file)


if __name__ == '__main__':
    main()