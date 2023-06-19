# Case-study-SARS-CoV-2

![](.rulegraph.png)


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


## Generating a rule graph

```bash
snakemake --rulegraph | dot -Tpng > .rulegraph.png
```
