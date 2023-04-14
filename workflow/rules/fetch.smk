rule fetch_reference:
    threads: 1
    conda: "../envs/fetch.yaml"
    output:
        fasta = OUTDIR/"reference.fasta"
    shell:
        "esearch -db nucleotide -query {config[REFSEQ_REFERENCE]} | efetch -format fasta > {output.fasta}"
