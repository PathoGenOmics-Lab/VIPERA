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


rule N_S_sites:
    threads: 1
    conda: "../envs/biopython.yaml"
    params:
        gb_qualifier_display = "gene",
    input:
        fasta = OUTDIR/f"{OUTPUT_NAME}.ancestor.fasta",
        gb = OUTDIR/"reference.cds.gb",
        genetic_code = Path(config["GENETIC_CODE_JSON"]).resolve(),
    output:
        csv = temp(OUTDIR/f"{OUTPUT_NAME}.ancestor.N_S.sites.csv"),
    log:
        LOGDIR / "N_S_sites" / "log.txt"
    script:
        "../scripts/N_S_sites_from_fasta.py"


rule dnds_data:
    conda: "../envs/renv.yaml"
    input: 
        n_s_sites = OUTDIR/f"{OUTPUT_NAME}.ancestor.N_S.sites.csv",
        variants =  OUTDIR/f"{OUTPUT_NAME}.variants.tsv",
        metadata = config["METADATA"]
    output:
        table = report(REPORT_DIR_TABLES/"dnds.csv")
    log:
        LOGDIR / "dnds_data" / "log.txt"
    script:
        "../scripts/dnds_data.R"
