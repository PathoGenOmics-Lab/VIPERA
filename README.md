# Case-study-SARS-CoV-2

## Instructions
To run the pipeline, first fetch the data (you may need to modify the script to include your credentials):

```bash
./fetch_data.sh
```

Then, within an environment with `snakemake>6.0` (see [the Snakemake docs](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)), run the following:

```bash
snakemake --use-conda -c8
```

You may change the `-c` argument to use a different number of CPUs.

## TO DO
- [ ] get reference with `samtools view -H some.bam` (BAM files may not be aligned to same reference)
- [ ] calculate BAM stats with samtools and summarize
- [ ] make a config table with infection groups to do ASR separately
- [ ] do Freyja demixing and summarize results
