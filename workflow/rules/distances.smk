rule weighted_distances:
    threads: 4
    conda: "../envs/biopython.yaml"
    params:
        mask_class = ["mask"],
        tsv_reference = OUTDIR/f"{OUTPUT_NAME}.ancestor.fasta",
        reference = OUTDIR/"reference.fasta"
    input:
        tsv = OUTDIR/f"{OUTPUT_NAME}.masked.filtered.tsv",
        vcf = OUTDIR/"problematic_sites.vcf"
    output:
        distances = OUTDIR/f"{OUTPUT_NAME}.weighted_distances.csv"
    log:
        LOGDIR / "weighted_distances" / "log.txt"
    script:
        "../scripts/weighted_distances.py"
