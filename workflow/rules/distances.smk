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
        mask_vcf = OUTDIR / "all_mask_sites.vcf",
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
        extra_args = "",
    input:
        variants = OUTDIR/f"{OUTPUT_NAME}.variants.afwdist.csv",
        reference = OUTDIR/f"{OUTPUT_NAME}.ancestor.fasta",
    output:
        distances = temp(OUTDIR/f"{OUTPUT_NAME}.distances.raw.csv"),
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
        samples = sorted(iter_samples()) + [config["ALIGNMENT_REFERENCE"]],
    input:
        distances = OUTDIR/f"{OUTPUT_NAME}.distances.raw.csv",
    output:
        distances = OUTDIR/f"{OUTPUT_NAME}.distances.csv",
    log:
        LOGDIR/"format_afwdist_results"/"log.txt"
    script:
        "../scripts/format_afwdist_results.py"


rule allele_freq_tree_data:
    conda: "../envs/renv.yaml"
    params:
        use_bionj = config["USE_BIONJ"],
        outgroup_id = config["ALIGNMENT_REFERENCE"],
    input:
        dist = OUTDIR/f"{OUTPUT_NAME}.distances.csv",
    output:
        tree = REPORT_DIR_TABLES/"allele_freq_tree.nwk",
    log:
        LOGDIR / "allele_freq_tree_data" / "log.txt"
    script:
        "../scripts/report/allele_freq_tree_data.R"


rule time_signal_data:
    conda: "../envs/renv.yaml"
    params:
        outgroup_id = config["ALIGNMENT_REFERENCE"],
    input:
        tree = report(REPORT_DIR_TABLES/"allele_freq_tree.nwk"),
        metadata = config["METADATA"],
    output:
        table = report(REPORT_DIR_TABLES/"time_signal.csv"),
        json = REPORT_DIR_TABLES/"time_signal.json",
    log:
        LOGDIR / "time_signal_data" / "log.txt"
    script:
        "../scripts/report/time_signal_data.R"
