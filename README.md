# VIPERA

<p align="center">
  <img src="logo.jpg" title="VIPERA logo">
</p>

[![PGO badge](https://img.shields.io/badge/PathoGenOmics-Lab-yellow.svg)](https://pathogenomics.github.io/)
[![Release](https://img.shields.io/github/v/release/PathoGenOmics-Lab/VIPERA)](https://github.com/PathoGenOmics-Lab/VIPERA/releases)
[![Snakemake](https://img.shields.io/badge/snakemake-≥7.19-brightgreen.svg?style=flat)](https://snakemake.readthedocs.io)
![Install workflow](https://github.com/PathoGenOmics-Lab/VIPERA/actions/workflows/install.yml/badge.svg)
![Test workflow](https://github.com/PathoGenOmics-Lab/VIPERA/actions/workflows/test.yml/badge.svg)

A pipeline for SARS-CoV-2 Virus Intra-Patient Evolution Reporting and Analysis.

## Instructions

To run VIPERA locally with the default configuration, you only need one line of code after
[installing Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)
and [configuring the inputs and outputs](config/README.md#inputs-and-outputs):

```shell
snakemake --use-conda -c4  # runs VIPERA on 4 cores (reccomended)
```

Alternatively, you can use a simple script ([`run_default_VIPERA.sh`](run_default_VIPERA.sh))
that downloads the data from our study and performs the analysis in one step.

This Snakemake workflow can be also [executed in an HPC environment](config/README.md#run-modes).

Please refer to [config/README.md](config/README.md) for detailed setup instructions.

## Contributors

[![Contributors figure](https://contrib.rocks/image?repo=PathoGenOmics-Lab/VIPERA)](https://github.com/PathoGenOmics-Lab/VIPERA/graphs/contributors)

## Citation

Manuscript in progress.
