#!/usr/bin/env bash

# BAM files
while read -r path; do
    scp garnatxa:${path} data/bam
done < data/bam_paths.txt

# FASTA files
while read -r path; do
    scp garnatxa:${path} data/fasta
done < data/fasta_paths.txt
