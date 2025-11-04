rule extract_afwdist_variants:
    conda: "../envs/biopython.yaml"
    params:
        sample_col = "SAMPLE",
        position_col = "POS",
        sequence_col = "ALT",
        frequency_col = "ALT_FREQ",
        mask_class = ["mask"],
    input:
        variants = OUTDIR/f"{OUTPUT_NAME}.variants.tsv",
        mask_vcf = lambda wildcards: select_problematic_vcf(),
        ancestor = OUTDIR/f"{OUTPUT_NAME}.ancestor.fasta",
        reference = OUTDIR/"reference.fasta",
    output:
        variants = temp(OUTDIR/f"{OUTPUT_NAME}.variants.afwdist.csv"),
    log:
        LOGDIR/"extract_afwdist_variants"/"log.txt"
    script:
        "../scripts/extract_afwdist_variants.py"


rule afwdist_weighted_distances:
    conda: "../envs/afwdist.yaml"
    params:
        extra_args = ""
    input:
        variants = OUTDIR/f"{OUTPUT_NAME}.variants.afwdist.csv",
        reference = OUTDIR/f"{OUTPUT_NAME}.ancestor.fasta",
    output:
        distances = temp(OUTDIR/"distances.raw.csv"),
    log:
        LOGDIR/"afwdist_weighted_distances"/"log.txt"
    shell:
        "afwdist "
        "-i {input.variants:q} "
        "-r {input.reference:q} "
        "-o {output.distances:q} "
        "{params.extra_args} >{log:q} 2>&1"


rule format_afwdist_results:
    conda: "../envs/biopython.yaml"
    params:
        samples = sorted(iter_samples()),
    input:
        distances = OUTDIR/"distances.raw.csv",
    output:
        distances = OUTDIR/"distances.csv",
    log:
        LOGDIR/"format_afwdist_results"/"log.txt"
    script:
        "../scripts/format_afwdist_results.py"
