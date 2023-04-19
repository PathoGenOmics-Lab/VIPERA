#!/usr/bin/env bash

# Crear carpetas
mkdir -p data/bam
mkdir -p data/fasta

# BAM files
while read -r path; do
    scp garnatxa:${path} data/bam
done < data/bam_paths.txt

# FASTA files
while read -r path; do
    scp garnatxa:${path} data/fasta
done < data/fasta_paths.txt
