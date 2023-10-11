rule N_S_sites:
    threads: 1
    conda: "../envs/biopython.yaml"
    input:
        fasta = OUTDIR/f"{OUTPUT_NAME}.ancestor.fasta",
        gb = OUTDIR/"reference.gb",
        features = config["FEATURES_JSON"],
        genetic_code = config["GENETIC_CODE_JSON"]
    output:
        csv = temp(OUTDIR/f"{OUTPUT_NAME}.ancestor.N_S.sites.csv")
    log:
        LOGDIR / "N_S_sites" / "log.txt"
    script:
        "../scripts/N_S_sites_from_fasta.py"
