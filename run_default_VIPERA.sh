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
MD5_SUM="ff59d513309a3af47f2c9248f7f3518d"

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
sm_version=$(snakemake --version | grep -oE "^[0-9]+")
if [ $sm_version -ge 8 ]; then
    logthis "Launching with snakemake >= 8"
    snakemake --software-deployment-method conda -c ${NCPU}
else
    logthis "Launching with snakemake < 8"
    snakemake --use-conda -c ${NCPU}
fi

logthis "Cleaning up"
rm -r ${tmpdir}

logthis "Done!"
