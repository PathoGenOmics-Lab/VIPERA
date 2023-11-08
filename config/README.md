# Instructions

To run the pipeline, you will need an environment with `snakemake`
(check [the Snakemake docs](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)).

## Inputs and outputs

Modify the [target configuration file](config/targets.yaml)
to point the `SAMPLES` and `METADATA` parameters to your data. The `OUTPUT_DIRECTORY`
parameter should point to your desired results directory.

The script [`build_targets.py`](build_targets.py) simplifies the process of creating
the configuration file. To run this script, you need to have PyYAML installed. It
takes a list of sample names, a directory with BAM and FASTA files, the path to
the metadata table and the name of your dataset as required inputs. Then, it searches the
directory for files that have the appropriate extensions and sample names and adds them
to the configuration file.

The file should look like this:

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
```

You may also provide this information through the `--config` parameter.

## Automated construction of a context dataset

Setting the `CONTEXT_FASTA` parameter to `null` (default) will enable
the automatic download of sequences from the GISAID SARS-CoV-2 database
(see [the following section](README.md#context-checkpoints) for further details).
To enable this, you must also [sign up to the GISAID platform](https://gisaid.org/register/)
and provide your credentials by creating and filling an additional configuration
file `config/gisaid.yaml` as follows:

```yaml
USERNAME: "your-username"
PASSWORD: "your-password"
```

## Mapping reference sequence

Setting `MAPPING_REFERENCES_FASTA` to `null` (default) will enable the automatic download of the
reference sequence(s) that were used to map the reads and generate the BAM files.
If the required sequence is not available publically or you already have it
at your disposal, it may be provided manually by setting the parameter to the
path of the reference FASTA file.

## Context checkpoints

By default, a dataset of samples that meet the spatial
and temporal criteria set through the [`download_context` rule](workflow/rules/context.smk):

- Location matching the place(s) of sampling of the target samples
- Collection date within the time window that includes 95% of the date distribution of the
target samples (2.5% is trimmed at each end to account for extreme values) ± 2 weeks
- Pango lineage matching that of the target samples

Then, a series of checkpoints are enforced:

- Remove context samples whose GISAID ID match any of the target samples
- Enforce a minimum number of samples to have at least as many possible combinations as bootstrap replicates for the diversity assessment (set in [the configuration file](config/config.yaml))

If these requirements are not met, a custom sequence dataset must be
provided through the `CONTEXT_FASTA` parameter by editing [the target configuration file](config/targets.yaml)
or via the command line:

```shell
snakemake --config CONTEXT_FASTA="path/to/fasta"
```

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

To run the analysis with the default configuration, just run the following command
(change the `-c/--cores` argument to use a different number of CPUs):

```shell
snakemake --use-conda -c8
```

To run the analysis in an HPC environment using SLURM, we provide a
[default profile configuration](profile/default/config.yaml) as an example that
should be modified to fit your needs. To use it, run the following command:

```shell
snakemake --use-conda --slurm --profile profile/default
```
