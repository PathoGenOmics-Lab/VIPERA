# VIPERA

<p align="center">
  <img src="logo.jpg" title="VIPERA logo">
</p>

[![PGO badge](https://img.shields.io/badge/PathoGenOmics-Lab-yellow.svg)](https://pathogenomics.github.io/)
[![DOI](https://img.shields.io/badge/Virus_Evolution-10.1093/ve/veae018-387088.svg)](https://doi.org/10.1093/ve/veae018)
[![Release](https://img.shields.io/github/v/release/PathoGenOmics-Lab/VIPERA)](https://github.com/PathoGenOmics-Lab/VIPERA/releases)
[![Snakemake](https://img.shields.io/badge/Snakemake-≥7.19-brightgreen.svg?style=flat)](https://snakemake.readthedocs.io)
![Install workflow](https://github.com/PathoGenOmics-Lab/VIPERA/actions/workflows/install.yml/badge.svg)
![Test workflow](https://github.com/PathoGenOmics-Lab/VIPERA/actions/workflows/test.yml/badge.svg)

A pipeline for SARS-CoV-2 Viral Intra-Patient Evolution Reporting and Analysis.

## First steps

To run VIPERA locally with the default configuration, you only need one line of code after
[installing Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html),
configuring [the inputs and outputs](config/README.md#inputs-and-outputs) and
[the context dataset](config/README.md#automated-construction-of-a-context-dataset):

```shell
snakemake --use-conda -c4  # runs VIPERA on 4 cores
```

We provide a simple script that downloads the [data](https://doi.org/10.20350/digitalCSIC/15648) from [our study](https://doi.org/10.1093/ve/veae018)
and performs the analysis in a single step:

```shell
./run_default_VIPERA.sh
```

This Snakemake workflow is compatible with both local execution and HPC environments utilizing SLURM. It supports dependency management through either conda or Singularity, as detailed in the [run modes documentation](config/README.md#run-modes).

Please refer to the [full workflow documentation](config/README.md) for detailed setup instructions.

> The documentation in this repository provides instructions for running VIPERA
> with Snakemake <8. We recommend using Snakemake 7.32.
> However, using Snakemake 8 is possible with minimal modifications (see the
> [migration documentation](https://snakemake.readthedocs.io/en/stable/getting_started/migration.html)).
> For example, `--use-conda` and `--use-singularity` are deprecated in Snakemake 8,
> and `--software-deployment-method conda apptainer` is the preferred way to provide the options. Additionally, SLURM support is only available after installing an
> [executor plugin for slurm](https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/slurm.html).

## Contributors

[![Contributors figure](https://contrib.rocks/image?repo=PathoGenOmics-Lab/VIPERA)](https://github.com/PathoGenOmics-Lab/VIPERA/graphs/contributors)

## Citation

> Álvarez-Herrera, M. & Sevilla, J., Ruiz-Rodriguez, P., Vergara, A., Vila, J., Cano-Jiménez, P., González-Candelas, F., Comas, I., & Coscollá, M. (2024). VIPERA: Viral Intra-Patient Evolution Reporting and Analysis. Virus Evolution, 10(1), veae018. https://doi.org/10.1093/ve/veae018

```bibtex
@article{ahs2024,
  title = {{VIPERA}: {Viral} {Intra}-{Patient} {Evolution} {Reporting} and {Analysis}},
  volume = {10},
  issn = {2057-1577},
  shorttitle = {{VIPERA}},
  url = {https://doi.org/10.1093/ve/veae018},
  doi = {10.1093/ve/veae018},
  abstract = {Viral mutations within patients nurture the adaptive potential of severe acute respiratory syndrome coronavirus 2 (SARS-CoV-2) during chronic infections, which are a potential source of variants of concern. However, there is no integrated framework for the evolutionary analysis of intra-patient SARS-CoV-2 serial samples. Herein, we describe Viral Intra-Patient Evolution Reporting and Analysis (VIPERA), a new software that integrates the evaluation of the intra-patient ancestry of SARS-CoV-2 sequences with the analysis of evolutionary trajectories of serial sequences from the same viral infection. We have validated it using positive and negative control datasets and have successfully applied it to a new case, which revealed population dynamics and evidence of adaptive evolution. VIPERA is available under a free software license at https://github.com/PathoGenOmics-Lab/VIPERA.},
  number = {1},
  journal = {Virus Evolution},
  author = {Álvarez-Herrera$^*$, Miguel and Sevilla$^*$, Jordi and Ruiz-Rodriguez, Paula and Vergara, Andrea and Vila, Jordi and Cano-Jiménez, Pablo and González-Candelas, Fernando and Comas, Iñaki and Coscollá, Mireia},
  month = jan,
  year = {2024},
  pages = {veae018},
  note = {$^*$ indicates equal contribution}
}
```
