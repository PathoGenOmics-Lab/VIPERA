#!/usr/bin/env bash

# run_default_VIPERA.sh
# 
# This script is used to download and analyze intra-patient SARS-CoV-2 data
# using VIPERA, a tool for the automated detection of serially sampled infections
# and the identification of evolutionary patterns within the same viral infection.
# 
# Access the full data record via DOI: 10.20350/digitalCSIC/15648

set -e

logthis () {
    echo $(date) "|" $@
}

NCPU=1
DATA_URL="https://digital.csic.es/bitstream/10261/337461/1/data.zip"
MD5_SUM="07447bdae794a6c82adbf79423b79c30"

tmpdir=$(mktemp -d)

logthis "Downloading compressed data from '$DATA_URL'"
curl -o ${tmpdir}/data.zip ${DATA_URL}

logthis "Validating file"
md5_sum_dwld="$(md5sum ${tmpdir}/data.zip | cut -d' ' -f1)"
if [ "$md5_sum_dwld" != "$MD5_SUM" ]; then
    logthis "ERROR: MD5 checksum does not match"
    exit 1
fi

logthis "Creating data directories"
mkdir -p data/bam data/fasta

logthis "Decompressing"
unzip -d ${tmpdir} ${tmpdir}/data.zip
rm ${tmpdir}/data.zip

logthis "Organizing files"
mv ${tmpdir}/**/*.bam data/bam
mv ${tmpdir}/**/*.fa data/fasta
mv ${tmpdir}/**/*.csv data

logthis "Running VIPERA"
snakemake --use-conda -c ${NCPU}

logthis "Cleaning up"
rm -r ${tmpdir}

logthis "Done!"
