rule N_S_sites:
    threads: 1
    conda: "../envs/biopython.yaml"
    params:
        genetic_code = "standard"
    input:
        fasta = OUTDIR/f"{OUTPUT_NAME}.ancestor.fasta"
    output:
        csv = OUTDIR/f"{OUTPUT_NAME}.ancestor.N_S.sites.csv"
    script:
        "../scripts/N_S_sites_from_fasta.py"
