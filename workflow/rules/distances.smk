rule weighted_distances:
    threads: 4
    conda: "../envs/biopython.yaml"
    params:
        samples = expand("{sample}", sample = iter_samples()),
        mask_class = ["mask"]
    input:
        tsv = OUTDIR/f"{OUTPUT_NAME}.masked.filtered.tsv",
        vcf = OUTDIR/"problematic_sites.vcf",
        ancestor = OUTDIR/f"{OUTPUT_NAME}.ancestor.fasta",
        reference = OUTDIR/"reference.fasta"
    output:
        distances = REPORT_DIR_TABLES/f"figure_4.csv"
    log:
        LOGDIR / "weighted_distances" / "log.txt"
    script:
        "../scripts/weighted_distances.py"
