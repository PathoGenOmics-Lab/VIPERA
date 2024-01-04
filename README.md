# VIPERA

<p align="center">
  <img src="logo.jpg" title="VIPERA logo">
</p>

[![PGO badge](https://img.shields.io/badge/PathoGenOmics-Lab-yellow.svg)](https://pathogenomics.github.io/)
[![DOI:10.1101/2023.10.24.561010](https://img.shields.io/badge/DOI-10.1101/2023.10.24.561010-blue.svg)](https://doi.org/10.1101/2023.10.24.561010)
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

We provide a simple script that downloads the [data](https://doi.org/10.20350/digitalCSIC/15648) from [our study](https://doi.org/10.1101/2023.10.24.561010)
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

Álvarez-Herrera M. & Sevilla, J., Ruiz-Rodriguez, P., Vergara, A., Vila, J., Cano-Jiménez, P., González-Candelas, F., Comas, I., & Coscolla, M. (2023). VIPERA: Viral Intra-Patient Evolution Reporting and Analysis. bioRxiv. https://doi.org/10.1101/2023.10.24.561010

```bibtex
@misc{AHS_VIPERA_2023,
  title = {{VIPERA}: {Viral} {Intra}-{Patient} {Evolution} {Reporting} and {Analysis}},
  shorttitle = {{VIPERA}},
  author = {Álvarez-Herrera$^*$, Miguel and Sevilla$^*$, Jordi and Ruiz-Rodriguez, Paula and Vergara, Andrea and Vila, Jordi and Cano-Jiménez, Pablo and González-Candelas, Fernando and Comas, Iñaki and Coscolla, Mireia},
  url = {https://www.biorxiv.org/content/10.1101/2023.10.24.561010},
  doi = {10.1101/2023.10.24.561010},
  language = {en},
  urldate = {2023-10-25},
  publisher = {bioRxiv},
  note = {$^*$ indicates equal contribution}
}
```
