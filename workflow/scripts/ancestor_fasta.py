#!/usr/bin/env python3


import logging
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


def filter_node_lines(file:str, node:int):
    with open(file) as f:
        for line in f:
            if not line.startswith("#"):
                fields = line.strip().split("\t")
                if fields[0] == node:
                    yield fields


if __name__ == "__main__":

    logging.basicConfig(filename=snakemake.log[0], format=snakemake.config["LOG_PY_FMT"], level=logging.INFO)

    logging.info("Reconstricting ancestral sequence for studied samples")
    sequence = []
    for fields in filter_node_lines(snakemake.input.state_file, snakemake.params.node_id):
        char = fields[2]
        # Replace uncertainties encoded by "-"
        if char == "-":
            sequence.append(snakemake.params.indeterminate_char)
        else:
            sequence.append(char)
    
    logging.info("Saving ancestral sequence")
    SeqIO.write(
        SeqRecord(
            Seq("".join(sequence)),
            snakemake.params.name,
            description=f"[Ancestral sequence reconstruction of node '{snakemake.params.node_id}']"
        ),
        snakemake.output.fasta,
        "fasta"
    )
