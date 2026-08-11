#!/usr/bin/env python3

import logging

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


def read_masked_sites(vcf_path: str, mask_classes: list[str]) -> set[int]:
    """
    Parse a VCF containing positions for masking. Assumes the VCF file is
    formatted as in:
    github.com/W-L/ProblematicSites_SARS-CoV2/blob/master/problematic_sites_sarsCov2.vcf
    with a "mask" or "caution" recommendation in column 7.
    Masked sites are specified with params.
    """
    vcf = pd.read_csv(
        vcf_path,
        sep="\\s+",
        comment="#",
        names=("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")
    )
    return set(vcf.loc[vcf.FILTER.isin(mask_classes), "POS"].unique())


def remove_ref_matching_variants(variants: pd.DataFrame, reference: Seq, position_col: str, sequence_col: str) -> pd.DataFrame:
    """Filters out variants whose ALT matches the base at the reference sequence (i.e. REF[1] == ALT)."""
    reference_bases = variants[position_col].map(lambda p: reference[p - 1]).astype(str)
    mask = reference_bases.str.upper() == variants[sequence_col].str.upper()
    if (n_dropped := mask.sum()) > 0:
        logging.warning(f"Dropping {n_dropped} rows where REF[1] equals ALT")
    return variants[~mask]


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
    variants = pd.read_table(
        snakemake.input.variants,
        sep="\t",
        usecols=list(colnames.keys())
    )

    logging.info(f"Read {len(variants)} variant records")
    # VCF with sites to mask
    masked_sites = read_masked_sites(snakemake.input.mask_vcf, snakemake.params.mask_class)
    logging.info(f"Read {len(masked_sites)} masked positions")

    logging.info("Reading input FASTA files")
    # Case ancestor (variant calling reference)
    ancestor = read_monofasta(snakemake.input.ancestor)
    logging.info(f"Ancestor: '{ancestor.description}', length={len(ancestor)}")
    
    # Alignment reference
    # TODO: clean unused variables
    reference = read_monofasta(snakemake.input.reference)
    logging.info(f"Reference: '{reference.description}', length={len(reference)}")
    assert len(ancestor) == len(reference)

    logging.info("Removing ancestor-matching variants")
    variants = remove_ref_matching_variants(
        variants,
        ancestor.seq,
        snakemake.params.position_col,
        snakemake.params.sequence_col
    ).drop_duplicates(keep="first")

    logging.info("Renaming and selecting columns")
    output = (
        variants
        .rename(columns=colnames)
        [list(colnames.values())]
        .astype(DTYPES)
    )
    logging.info("Filtering sites")
    output = output[~output.position.isin(masked_sites)]
    logging.info(f"There are {len(output)} rows left")
    
    logging.info("Writing results")
    output.to_csv(snakemake.output.variants, index=False)
