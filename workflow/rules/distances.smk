rule weighted_distances:
    threads: 1
    conda: "../envs/biopython.yaml"
    params:
        samples = expand("{sample}", sample = iter_samples()),
        mask_class = ["mask"]
    input:
        tsv = OUTDIR/f"{OUTPUT_NAME}.masked.filtered.tsv",
        vcf = lambda wildcards: select_problematic_vcf(),
        ancestor = OUTDIR/f"{OUTPUT_NAME}.ancestor.fasta",
        reference = OUTDIR/"reference.fasta"
    output:
        distances = REPORT_DIR_TABLES/f"distances.csv"
    log:
        LOGDIR / "weighted_distances" / "log.txt"
    script:
        "../scripts/weighted_distances.py"
