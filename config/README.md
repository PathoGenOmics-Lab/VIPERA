# Instructions

To run VIPERA, an environment with Snakemake version 7.19 or later is needed
(see [the Snakemake docs](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)
for setup instructions).

> This guide provides command-line instructions for running VIPERA with Snakemake versions prior to 8.
All configuration parameters are fully cross-compatible. The original publication used Snakemake 7.32,
but newer versions can also be used with only minor changes. For details, see the Snakemake
[migration guide](https://snakemake.readthedocs.io/en/stable/getting_started/migration.html).
For example, existing [profiles](https://snakemake.readthedocs.io/en/stable/getting_started/migration.html#profiles)
are cross-compatible as well, but note that the `--use-conda` flag is deprecated starting with Snakemake 8.
Instead, use `--software-deployment-method conda`.

## Inputs and outputs

The workflow requires a set of FASTA files (one per target sample), a corresponding set of
BAM files (also one per target sample), and a metadata table in CSV format with one row per
sample. The metadata must include the following columns: unique sample identifier (default column `ID`,
used to match sequencing files with metadata), the date the sample was collected (default `CollectionDate`),
the location where the sample was collected (default `ResidenceCity`), and GISAID accession (default `GISAIDEPI`).
The default column names but can be customized if needed via the workflow parameters.

These parameters are set in two configuration files in YAML format:
[config.yaml](/config/config.yaml) (for general workflow settings) and
[targets.yaml](/config/targets.yaml) (for specific dataset-related settings).
The latter must be modified by the user to point the `SAMPLES` and `METADATA`
parameters to your data. The `OUTPUT_DIRECTORY` parameter should point to your
desired results directory.

The script [`build_targets.py`](/build_targets.py) simplifies the process of creating
the targets configuration file. To run this script, you need to have PyYAML installed. It
takes a list of sample names, a directory with BAM and FASTA files, the path to
the metadata table and the name of your dataset as required inputs. Then, it searches the
directory for files that have the appropriate extensions and sample names and adds them
to the configuration file.

An example file could look like this:

```yaml
OUTPUT_NAME:
  "your-dataset-name"
SAMPLES:
  sample1:
    bam: "path/to/sorted/bam1.bam"
    fasta: "path/to/sequence1.fasta"
  sample2:
    bam: "path/to/sorted/bam2.bam"
    fasta: "path/to/sequence2.fasta"
  ...
METADATA:
  "path/to/metadata.csv"
OUTPUT_DIRECTORY:
  "output"
CONTEXT_FASTA:
  null
MAPPING_REFERENCES_FASTA:
  null
```

This information may also be provided through the `--config` parameter.

## Automated construction of a context dataset

Setting the `CONTEXT_FASTA` parameter to `null` (default) will enable
the automatic download of sequences from the GISAID EpiCoV SARS-CoV-2 database.
An unset parameter has the same effect.
To enable this, you must also [sign up to the GISAID platform](https://gisaid.org/register/)
and provide your credentials by creating and filling an additional configuration
file (default: `config/gisaid.yaml`) as follows:

```yaml
USERNAME: "your-username"
PASSWORD: "your-password"
```

A set of samples that meet the spatial, temporal and phylogenetic criteria
set through the [`download_context` rule](/workflow/rules/context.smk#L1)
will be retrieved automatically from GISAID. These criteria are:

- Location matching the place(s) of sampling of the target samples
- Collection date within the time window that includes 95% of the date distribution of the
target samples (2.5% is trimmed at each end to account for extreme values) ± 2 weeks
- Pango lineage matching that of the target samples

Then, a series of checkpoint steps are executed for quality assurance:

- Remove context samples whose GISAID ID match any of the target samples
- Enforce a minimum number of samples to have at least as many possible
  combinations as random subsample replicates for the diversity assessment
  (set in [config.yaml](/config/config.yaml))

The workflow will continue its execution until completion if the obtained
context dataset passes these checkpoints. Otherwise, the execution will be
terminated and, to continue the analysis, an external context dataset must
be provided through the `CONTEXT_FASTA` parameter. This can be done
by editing [targets.yaml](/config/targets.yaml) or via the command line:

```shell
snakemake --config CONTEXT_FASTA="path/to/fasta"
```

> [!IMPORTANT]
> The GISAID EpiCoV database is proprietary and not openly accessible. For details, refer
> to the [GISAID Terms of Use](https://gisaid.org/terms-of-use/).
> VIPERA uses [GISAIDR](https://github.com/Wytamma/GISAIDR) to automate access to GISAID data.
> However, this access can be unstable or occasionally fail due to changes in the platform.
> Possible workarounds are documented in the GISAIDR repository (e.g. issues
> [#55](https://github.com/Wytamma/GISAIDR/issues/55) and [#58](https://github.com/Wytamma/GISAIDR/issues/58)).
> If programmatic access fails, a suitable context dataset must be manually provided by setting the
> `CONTEXT_FASTA` parameter to the path of a FASTA file.
> As a last resort, some of the analyses can be allowed to run even if context-dependent rules
> fail by passing the `--keep-going` flag to Snakemake.
> To replicate [our work](https://doi.org/10.1093/ve/veae018), the automatic context dataset is available via
> DOI: [10.55876/gis8.250718er](https://doi.org/10.55876/gis8.250718er) (EPI_SET_250718er).
> Read more about EPI_SETs [here](https://gisaid.org/episet/).

## Mapping reference sequence

Setting `MAPPING_REFERENCES_FASTA` to `null` (default) will enable the automatic download of the
reference sequence(s) that were used to map the reads and generate the BAM files.
An unset parameter has the same effect.
If the required sequence is not available publically or the user already has it
at your disposal, it can be provided manually by setting the parameter to the
path of the reference FASTA file.

## Workflow configuration variables

All of the following variables are pre-defined in [config.yaml](/config/config.yaml):

- `ALIGNMENT_REFERENCE`: NCBI accession number of the reference record for sequence alignment.
- `PROBLEMATIC_VCF`: URL or path of a VCF file containing problematic genome positions for masking.
- `GENETIC_CODE_JSON`: path of a JSON file containing a genetic code for gene translation.
- `TREE_MODEL`: substitution model used by IQTREE (see [docs](http://www.iqtree.org/doc/Substitution-Models)).
- `UFBOOT` & `SHALRT`: settings for ultrafast bootstrap (see [UFBoot](https://doi.org/10.1093/molbev/msx281))
  and Shimodaira–Hasegawa approximate likelihood ratio test bootstrap
  (see [SH-aLRT](https://doi.org/10.1093/sysbio/syq010)) in IQTREE runs:
  - `REPS`: number of replicates.
  - `THRESHOLD`: value cutoff for visualization.
- `VC`: variant calling configuration:
  - `MAX_DEPTH`: maximum depth at a position for `samtools mpileup` (option `-d`).
  - `MIN_QUALITY`: minimum base quality for `samtools mpileup` (option `-Q`).
  - `IVAR_QUALITY`: minimum base quality for `ivar variants` (option `-q`).
  - `IVAR_FREQ`: minimum frequency threshold for `ivar variants` (option `-t`).
  - `IVAR_DEPTH`: minimum read depth for `ivar variants` (option `-m`).
- `DEMIX`: demixing configuration (uses [Freyja](https://github.com/andersen-lab/Freyja), see also [its docs](https://andersen-lab.github.io/Freyja/index.html)):
  - `PATHOGEN`: pathogen of interest for `freyja update` (option `--pathogen`); must be 'SARS-CoV-2'.
  - `MIN_QUALITY`: minimum quality for `freyja variants` (option `--minq`).
  - `COV_CUTOFF`: minimum depth to calculate the reported "coverage" (percent of sites with that depth) for `freyja demix` (option `--covcut`).
  - `MIN_ABUNDANCE`: minimum lineage estimated abundance for `freyja demix` (option `--eps`).
  - `CONFIRMED_ONLY`: exclude unconfirmed lineages in `freyja demix` (option `--confirmedonly`).
  - `DEPTH_CUTOFF`: minimum depth on each site for `freyja demix` (option `--depthcutoff`).
  - `RELAXED_MRCA`: assign clusters using relaxed (as opposed to strict) MRCA, used with `DEPTH_CUTOFF`, for `freyja demix` (option `--relaxedmrca`).
  - `RELAXED_MRCA_THRESH`: `RELAXED_MRCA` threshold for `freyja demix` (option `--relaxedthresh`).
  - `AUTO_ADAPT`: use error profile to set adaptive lasso penalty parameter for `freyja demix` (option `--autoadapt`).
- `WINDOW`: sliding window of nucleotide variants per site configuration:
  - `WIDTH`: number of sites within windows.
  - `STEP`: number of sites between windows.
- `GB_FEATURES`: optional mapping to filter which features from the GenBank file are used
  by some analyses (e.g. rules `window` and `N_S_sites`).
  If `GB_FEATURES` is empty or unset, all features are used.
  Filtering is applied in order:
  - `INCLUDE`: mapping of qualifier names to sequences of values.
    If present, only features that match at least one key/value pair in `INCLUDE`
    are included in the analyses. For example, having `INCLUDE: {gene: [S, N]}` keeps features whose `gene`
    qualifier equals `S` or `N`.
  - `EXCLUDE`: mapping of qualifier names to sequences of values.
    After `INCLUDE` is applied, any feature that matches any key/value pair in `EXCLUDE`
    is omitted. For example, having `EXCLUDE: {gene: [S, N]}` removes features whose
    `gene` qualifier equals `S` or `N`.
- `GISAID`: automatic context download configuration.
  - `CREDENTIALS`: path of the GISAID credentials in YAML format.
  - `DATE_COLUMN`: name of the column that contains sampling dates (YYYY-MM-DD) in the input target metadata.
  - `LOCATION_COLUMN`: name of the column that contains sampling locations (e.g. city names) in the input target metadata.
  - `ACCESSION_COLUMN`: name of the column that contains GISAID EPI identifiers in the input target metadata.
- `DIVERSITY_REPS`: number of random sample subsets of the context dataset for the nucleotide diversity comparison.
- `USE_BIONJ`: use the BIONJ algorithm ([Gascuel, 1997](https://doi.org/10.1093/oxfordjournals.molbev.a025808)) instead of NJ (neighbor-joining; [Saitou & Nei, 1987](https://doi.org/10.1093/oxfordjournals.molbev.a040454)) to reconstruct phylogenetic trees from pairwise distances.
- `COR`: configuration for correlation analyses of allele frequency data over time and between variants. This parameter controls how correlation tests are performed using R's `cor.test` and `cor` functions (see [R documentation](https://search.r-project.org/CRAN/refmans/correlation/html/cor_test.html)).
  - `METHOD`: correlation method to use. Valid options are "pearson" (default), "kendall", or "spearman".
  - `EXACT`: boolean flag indicating whether to compute an exact p-value when possible. This option applies only to certain methods and may be set to `null` (default) to let R decide automatically.
- `LOG_PY_FMT`: logging format string for Python scripts.
- `PLOTS`: path of the R script that sets the design and style of data visualizations.
- `PLOT_GENOME_REGIONS`: path of a CSV file containing genome regions, e.g. SARS-CoV-2 non-structural protein (NSP) coordinates, for data visualization.
- `REPORT_QMD`: path of the report template in Quarto markdown (QMD) format.

TODO: describe new entries in config.yaml
TODO: describe issues with many cores (rules freyja and snps_to_ancestor)

## Workflow graphs

To generate a simplified rule graph, run:

```shell
snakemake --rulegraph | dot -Tpng > .rulegraph.png
```

![Snakemake rule graph](/.rulegraph.png)

To generate the directed acyclic graph (DAG) of all rules
to be executed, run:

```shell
snakemake --forceall --dag | dot -Tpng > .dag.png
```

![Snakemake rule graph](/.dag.png)

## Run modes

To run the analysis with the default configuration, run the following command
(change the `-c/--cores` argument to use a different number of CPUs):

```shell
snakemake --use-conda -c4
```

To run the analysis in an HPC environment using SLURM, we provide a
default [SLURM profile configuration](/profile/slurm) as an example that
should be modified to fit your needs. Read more about Snakemake profiles
[here](https://snakemake.readthedocs.io/en/stable/executing/cli.html#executing-profiles).
To use the profile, install the
[Snakemake executor plugin for SLURM](https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/slurm.html)
and run one of the following commands:

```shell
snakemake --slurm --profile profile/slurm  # Snakemake v7
snakemake --profile profile/slurm          # Snakemake v8+
```

Additionally, we offer the option of running the workflow within a containerized
environment using a [pre-built Docker image](https://hub.docker.com/r/ahmig/vipera),
provided that [Apptainer/Singularity](https://en.wikipedia.org/wiki/Singularity_(software))
is available on the system. This eliminates the need for further conda package
downloads and environment configuration.
To do that, simply add the option `--use-apptainer` to any of the previous commands.

Using Apptainer for running VIPERA in the Windows Subsystem for Linux (WSL)
may encounter errors due to the default file permissions configuration, which
conflicts with Snakemake's containerized conda environment activation mechanism.
Thus, running the containerized VIPERA workflow on the WSL is not advised.
Additionally, certain known issues arise when utilizing non-default temporary
directories and Snakemake shadow directories. To address this issue, use the
default temporary directory (e.g. `export TMPDIR=/tmp` in Linux machines) and
specify the shadow prefix (`--shadow-prefix /tmp`) before executing the containerized workflow.
