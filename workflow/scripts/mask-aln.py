#!/usr/bin/env python3

# ---------------------------------------------------------------------------------------
# Copyright (C) 2020 EMBL - European Bioinformatics Institute
# Contact: goldman@ebi.ac.uk, cwalker@ebi.ac.uk
#
# This program is free software: you can redistribute it and/or modify it under the terms
# of the GNU General Public License as published by the Free Software Foundation,
# either version 3 of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
# more details.
#
# You should have received a copy of the GNU General Public License along with this
# program. If not, see <http://www.gnu.org/licenses/>.
# ---------------------------------------------------------------------------------------

"""
Source: https://github.com/W-L/ProblematicSites_SARS-CoV2
Modified by Miguel Alvarez <m.alvarez.herrera@csic.es>
"""

import logging
from Bio import SeqIO
import pandas as pd


def read_fasta_keep_name(file:str):
    sample_headers = []
    sample_sequences = []
    for record in SeqIO.parse(file, "fasta"):
        sample_headers.append(str(record.description))
        sample_sequences.append(str(record.seq))
    return sample_headers, sample_sequences


def ref_coords_to_align_coords(ref_align_seq) -> dict:
    """
    Generate a dictionary of reference sequence coordinates mapped to MSA
    coordinates. Used to update the VCF positions in case the reference
    sequence contains gaps.
    """
    ref_coord_dic = {}
    ref_align_seq = "".join(ref_align_seq).strip()
    for i in range(len(ref_align_seq.replace("-",""))):
        ref_coord_dic[i] = 0
    seq_count = 0
    align_count = 0
    for c in ref_align_seq:
        if c != "-":
            ref_coord_dic[seq_count] = align_count
            seq_count += 1
            align_count += 1
        else:
            align_count += 1
    return ref_coord_dic


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
        delim_whitespace=True,
        comment="#",
        names=("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")
    )
    return tuple(vcf.loc[vcf.FILTER.isin(mask_classes), "POS"])


def main():
    logging.basicConfig(filename=snakemake.log[0], format=snakemake.config["LOG_PY_FMT"], level=logging.INFO)
    
    logging.info("getting sites to mask")
    # parse VCF to get list of sites to mask
    positions = read_problematic_positions(snakemake.input.vcf, snakemake.params.mask_class)

    logging.info("Parsing Fasta file")
    # parse existing MSA FASTA file
    with open(snakemake.input.fasta, "r") as fasta_fi:
        headers, sequences = read_fasta_keep_name(fasta_fi)

    # get reference sequence
    with open(snakemake.input.ref_fasta, "r") as fasta_fi:
        _, ref_sequence = read_fasta_keep_name(fasta_fi)
        ref_sequence = ref_sequence[0]  # only 1 record -hopefully

    # convert reference sequence coords to corresponding alignment coords
    ref_coord_dic = ref_coords_to_align_coords(ref_sequence)

    # convert list of iffy sites to corresponding alignment positions
    # using zero-based indexing
    # Mask problematic sites from VCF and leading and trailing positions
    zb_sites = tuple(ref_coord_dic[int(i)-1] for i in positions)

    logging.info("Writing masked alignment")
    # using list of iffy alignment sites, mask corresponding positions
    # in the input MSA
    # Also: only mask sites with letters - no gaps!
    if snakemake.params.remove_sites:
        mask_sub_char = ""
    else:
        mask_sub_char = snakemake.params.mask_character

    f = open(snakemake.output.fasta, "w")
    for i in range(len(headers)):
        seq = list(sequences[i])
        for p in zb_sites:
            if seq[p] != "-":
                seq[p] = mask_sub_char
        seq = "".join(seq)
        f.write(">" + headers[i] + "\n" + seq + "\n")
    f.close()


if __name__ == "__main__":
    main()
