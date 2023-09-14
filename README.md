# VIPERA

<p align="center">
  <img src="logo.jpg" title="VIPERA logo">
</p>

[![PGO badge](https://img.shields.io/badge/PathoGenOmics-Lab-yellow.svg)](https://pathogenomics.github.io/)
[![Release](https://img.shields.io/github/v/release/PathoGenOmics-Lab/Case-study-SARS-CoV-2)](https://github.com/PathoGenOmics-Lab/Case-study-SARS-CoV-2/releases)
[![Snakemake](https://img.shields.io/badge/snakemake-≥7.19-brightgreen.svg?style=flat)](https://snakemake.readthedocs.io)
![Install workflow](https://github.com/PathoGenOmics-Lab/Case-study-SARS-CoV-2/actions/workflows/install.yml/badge.svg)
![Test workflow](https://github.com/PathoGenOmics-Lab/Case-study-SARS-CoV-2/actions/workflows/test.yml/badge.svg)

A pipeline for SARS-CoV-2 Virus Intra-Patient Evolution Reporting and Analysis.

## Instructions

You can run VIPERA locally with the default configuration with one line of code after
[installing Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)
and [configuring inputs and outputs](config/README.md#inputs-and-outputs):

```shell
snakemake --use-conda -c4
```

The pipeline also allows for its [execution in an HPC environment](config/README.md#run-modes).
Please check [config/README.md](config/README.md) for in-detail setup instructions.

## Methodology

### Pairwise distances between samples

In order to describe in a better way the relationship between samples, distances beween them are calculated tacking into account allele frequencies in sequencing data. Our aproach to compute the distance between sets of allele frequencies is based on FST formula and define the distance between two samples as:

```math
d(M,N)=\sum\limits_{i=1}^I \frac {\sum\limits_{j=1}^J (M_{ij} -N_{ij})^2 } {4 - \sum\limits_{j=1}^J (M_{ij} +N_{ij})^2 }
```

where:

$M$ y $N$: Two sequences.

$i = 1... I :$ Index over polymorphic sites.

$j = 1... J :$ Index over alleles.

$M_{ij}$ : Frequency of allel $j$ in position $i$ for sequence $M$.

## Contributors

[![Contributors figure](https://contrib.rocks/image?repo=PathoGenOmics-Lab/Case-study-SARS-CoV-2)](https://github.com/PathoGenOmics-Lab/Case-study-SARS-CoV-2/graphs/contributors)

## Citation

Manuscript in progress.
