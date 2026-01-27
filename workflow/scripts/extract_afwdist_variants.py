#!/usr/bin/env python3

import logging
from typing import List, Set

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


def read_monofasta(path: str) -> SeqRecord:
    fasta = SeqIO.parse(path, "fasta")
    record = next(fasta)
    if next(fasta, None) is not None:
        logging.warning(f"There are unread records left in '{path}'")
    return record


def read_masked_sites(vcf_path: str, mask_classes: List[str]) -> Set[int]:
    """
    Parse a VCF containing positions for masking. Assumes the VCF file is
    formatted as in:
    github.com/W-L/ProblematicSites_SARS-CoV2/blob/master/problematic_sites_sarsCov2.vcf
    with a "mask" or "caution" recommendation in column 7.
    Masked sites are specified with params.
    """
    vcf = pd.read_csv(
        vcf_path,
        sep="\s+",
        comment="#",
        names=("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")
    )
    return set(vcf.loc[vcf.FILTER.isin(mask_classes), "POS"].unique())


def build_ancestor_variant_table(ancestor: Seq, reference: Seq, reference_name: str, masked_positions: Set[int]) -> pd.DataFrame:
    pos = []
    alt = []
    for i in range(1, len(ancestor) + 1):
        if i not in masked_positions and ancestor[i-1] != reference[i-1]:
            pos.append(i)
            alt.append(reference[i-1])
    df = pd.DataFrame({snakemake.params.position_col: pos, snakemake.params.sequence_col: alt})
    df[snakemake.params.frequency_col] = 1  # As a reference genome, we assume all positions have fixed alleles
    df[snakemake.params.sample_col] = reference_name
    return df


DTYPES = {
    "sample": "object",
    "position": "int64",
    "sequence": "object",
    "frequency": "float64"
}


if __name__ == "__main__":

    logging.basicConfig(filename=snakemake.log[0], format=snakemake.config["LOG_PY_FMT"], level=logging.INFO)

    colnames = {
        snakemake.params.sample_col: "sample",
        snakemake.params.position_col: "position",
        snakemake.params.sequence_col: "sequence",
        snakemake.params.frequency_col: "frequency"
    }

    logging.info("Reading input tables")
    # Variants
    variants = pd.read_table(snakemake.input.variants, sep="\t")
    logging.info(f"Read {len(variants)} variant records")
    # VCF with sites to mask
    masked_sites = read_masked_sites(snakemake.input.mask_vcf, snakemake.params.mask_class)
    logging.info(f"Read {len(masked_sites)} masked positions")

    logging.info("Reading input FASTA files")
    # Case ancestor
    ancestor = read_monofasta(snakemake.input.ancestor)
    logging.info(f"Ancestor: '{ancestor.description}', length={len(ancestor)}")
    # Alignment reference
    reference = read_monofasta(snakemake.input.reference)
    logging.info(f"Reference: '{reference.description}', length={len(reference)}")
    assert len(ancestor) == len(reference)

    logging.info("Processing ancestor variants")
    ancestor_table = build_ancestor_variant_table(ancestor.seq, reference.seq, reference.id, masked_sites)
    logging.info(f"Ancestor has {len(ancestor_table)} variants")
    all_variants = pd.concat([variants, ancestor_table], ignore_index=True)
    logging.info(f"Combined table has {len(all_variants)} variants")

    logging.info("Renaming and selecting columns")
    output = (
        all_variants
        .rename(columns=colnames)
        [list(colnames.values())]
        .astype(DTYPES)
    )
    logging.info("Filtering sites")
    output = output[~output.position.isin(masked_sites)]
    logging.info(f"There are {len(output)} rows left")
    
    logging.info("Writing results")
    output.to_csv(snakemake.output.variants, index=False)
