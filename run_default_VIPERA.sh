#!/usr/bin/env bash

set -e

logthis () {
    echo $(date) "|" $@
}

NCPU=1
ZENODO_URL=""

tmpdir=$(mktemp -d)

logthis "Downloading compressed data"
curl -O ${tmpdir}/zenodo.zip ${ZENODO_URL}

logthis "Creating directories"
mkdir -p data/bam data/fasta

logthis "Decompressing"
unzip -d ${tmpdir} ${tmpdir}/zenodo.zip
rm ${tmpdir}/zenodo.zip

logthis "Organizing files"
mv ${tmpdir}/*.bam data/bam
mv ${tmpdir}/*.fa data/fasta

logthis "Running VIPERA"
snakemake --use-conda -c ${NCPU}

logthis "Cleaning up"
rmdir ${tmpdir}

logthis "Done!"
