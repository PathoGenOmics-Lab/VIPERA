rule filter_genbank_features:
    threads: 1
    conda: "../envs/biopython.yaml"
    params:
        included = config.get("GB_FEATURES", {}).get("INCLUDE", {}),
        excluded = config.get("GB_FEATURES", {}).get("EXCLUDE", {}),
    input:
        gb = OUTDIR/"reference.gb",
    output:
        gb = OUTDIR/"reference.cds.gb",
    log:
        LOGDIR / "filter_genbank_features" / "log.txt"
    script:
        "../scripts/filter_genbank_features.py"


rule n_s_sites:
    threads: 1
    conda: "../envs/biopython.yaml"
    params:
        gb_qualifier_display = "gene",
    input:
        fasta = OUTDIR/f"{OUTPUT_NAME}.ancestor.fasta",
        gb = OUTDIR/"reference.cds.gb",
        genetic_code = Path(config["GENETIC_CODE_JSON"]).resolve(),
    output:
        csv = temp(OUTDIR/f"{OUTPUT_NAME}.ancestor.n_s.sites.csv"),
    log:
        LOGDIR / "n_s_sites" / "log.txt"
    script:
        "../scripts/n_s_sites_from_fasta.py"


rule calculate_dnds:
    conda: "../envs/renv.yaml"
    input: 
        n_s_sites = OUTDIR/f"{OUTPUT_NAME}.ancestor.n_s.sites.csv",
        variants =  OUTDIR/f"{OUTPUT_NAME}.variants.tsv",
        metadata = config["METADATA"]
    output:
        table = OUTDIR/f"{OUTPUT_NAME}.dnds.csv",
    log:
        LOGDIR / "calculate_dnds" / "log.txt"
    script:
        "../scripts/calculate_dnds.R"
