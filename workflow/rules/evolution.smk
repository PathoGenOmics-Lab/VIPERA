rule N_S_sites:
    threads: 1
    conda: "../envs/biopython.yaml"
    params:
        features = config.get("GB_FEATURES", {}),
        gb_qualifier_display = "gene",
    input:
        fasta = OUTDIR/f"{OUTPUT_NAME}.ancestor.fasta",
        gb = OUTDIR/"reference.gb",
        genetic_code = Path(config["GENETIC_CODE_JSON"]).resolve(),
    output:
        csv = temp(OUTDIR/f"{OUTPUT_NAME}.ancestor.N_S.sites.csv"),
    log:
        LOGDIR / "N_S_sites" / "log.txt"
    script:
        "../scripts/N_S_sites_from_fasta.py"
