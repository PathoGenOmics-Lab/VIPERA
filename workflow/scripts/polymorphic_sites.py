#!/usr/bin/env python3

import pandas as pd
from Bio import SeqIO


def encode_nucleotide(nt):
    if nt == snakemake.params.indeterminate_char:
        return None
    else:
        return nt


def is_different(char1, char2):
    return char1 != char2


def is_different_indet(char1, char2):
    return (char1 != char2) and \
        (snakemake.params.indeterminate_char != char1) and \
        (snakemake.params.indeterminate_char != char2)


if __name__ == "__main__":

    # Set dif function
    if snakemake.params.include_indeterminations:
        dif = is_different
    else:
        dif = is_different_indet

    # Read reference and init table
    reference_record = SeqIO.read(snakemake.input.reference_fasta, "fasta")
    
    # Get polymorphic sites
    sites_0based = set()
    for record in SeqIO.parse(snakemake.input.alignment_fasta, "fasta"):
        # Check sequence length
        if len(record.seq) != len(reference_record.seq):
            raise ValueError("Length of record '{record.id}' is different from reference '{reference_record.id}'")
        # Fill sites
        if record.id not in snakemake.params.omit_record_IDs:
            for i in range(len(reference_record.seq)):
                if dif(record.seq[i], reference_record.seq[i]):
                    sites_0based.add(i)
    sites_0based = sorted(sites_0based)

    # Init dataframe
    df = pd.DataFrame(columns=[i+1 for i in sites_0based])

    # Fill reference row
    df.loc[reference_record.description] = [encode_nucleotide(reference_record.seq[i]) for i in sites_0based]
    
    # Fill alignment rows
    for record in SeqIO.parse(snakemake.input.alignment_fasta, "fasta"):
        df.loc[record.description] = [encode_nucleotide(record.seq[i]) for i in sites_0based]

    # Write table
    df.index.name = "SampleID"
    df.to_csv(snakemake.output.table, sep=",")
