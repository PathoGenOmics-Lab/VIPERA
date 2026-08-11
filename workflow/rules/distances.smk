rule extract_afwdist_variants:
    conda: "../envs/biopython.yaml"
    params:
        sample_col = "SAMPLE",
        position_col = "POS",
        sequence_col = "ALT",
        frequency_col = "ALT_FREQ",
        mask_class = ["mask"],
    input:
        variants = "<results>/<dataset>/variants.tsv",
        mask_vcf = "<results>/<dataset>/all_mask_sites.vcf",
        ancestor = "<results>/<dataset>/ancestor.fasta",
        reference = "<results>/<dataset>/reference.fasta",
    output:
        variants = temp("<results>/<dataset>/variants.afwdist.csv"),
    log:
        "<logs>/<dataset>/extract_afwdist_variants/log.txt"
    script:
        "../scripts/extract_afwdist_variants.py"


rule afwdist_weighted_distances:
    conda: "../envs/afwdist.yaml"
    params:
        extra_args = "-s",
    input:
        variants = "<results>/<dataset>/variants.afwdist.csv",
        reference = "<results>/<dataset>/ancestor.fasta",
    output:
        distances = temp("<results>/<dataset>/distances.raw.csv"),
    log:
        "<logs>/<dataset>/afwdist_weighted_distances/log.txt"
    shell:
        "afwdist "
        "-i {input.variants:q} "
        "-r {input.reference:q} "
        "-o {output.distances:q} "
        "{params.extra_args} >{log:q} 2>&1"


rule format_afwdist_results:
    conda: "../envs/biopython.yaml"
    params:
        samples = sorted(iter_samples()) + ["case_ancestor"],
    input:
        distances = "<results>/<dataset>/distances.raw.csv",
    output:
        distances = "<results>/<dataset>/distances.csv",
    log:
        "<logs>/<dataset>/format_afwdist_results/log.txt"
    script:
        "../scripts/format_afwdist_results.py"


rule allele_freq_tree_data:
    conda: "../envs/renv.yaml"
    params:
        use_bionj = config["USE_BIONJ"],
        outgroup_id = "case_ancestor",
    input:
        dist = "<results>/<dataset>/distances.csv",
    output:
        tree = "<results>/<dataset>/report/tables/allele_freq_tree.nwk",
    log:
        "<logs>/<dataset>/allele_freq_tree_data/log.txt"
    script:
        "../scripts/report/allele_freq_tree_data.R"


rule time_signal_data:
    conda: "../envs/renv.yaml"
    params:
        outgroup_id = "case_ancestor",
        confidence_interval = 0.95,
    input:
        tree = report("<results>/<dataset>/report/tables/allele_freq_tree.nwk"),
        metadata = config["METADATA"],
    output:
        table = report("<results>/<dataset>/report/tables/time_signal.csv"),
        json = "<results>/<dataset>/report/tables/time_signal.json",
    log:
        "<logs>/<dataset>/time_signal_data/log.txt"
    script:
        "../scripts/report/time_signal_data.R"
