rule fetch_alignment_reference:
    threads: 1
    conda: "../envs/fetch.yaml"
    params:
        ref = config["ALIGNMENT_REFERENCE"]
    output:
        fasta = OUTDIR/"reference.fasta"
    log:
        LOGDIR / "fetch_alignment_reference" / "log.txt"
    shell:
        "esearch -db nucleotide -query {params.ref} | efetch -format fasta > {output.fasta} 2> {log}"


rule fetch_mapping_references:
    threads: 1
    conda: "../envs/fetch.yaml"
    params:
        refs = config["MAPPING_REFERENCES"]
    output:
        fasta = OUTDIR/"mapping_references.fasta"
    log:
        LOGDIR / "fetch_mapping_references" / "log.txt"
    shell:
        """
        # Create empty file
        echo -n > {output.fasta}
        # Download and concatenate each reference
        echo {params.refs} | while read -r ref_id; do
            esearch -db nucleotide -query "$ref_id" | efetch -format fasta >> {output.fasta} 2>> {log}
        done
        """


rule fetch_alignment_annotation:
    threads: 1
    params:
        ref = config["ALIGNMENT_REFERENCE"]
    output:
        OUTDIR/"reference.gff3"
    log:
        LOGDIR / "fetch_alignment_annotation" / "log.txt"
    shell:
        "curl 'https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=gff3&id={params.ref}' -o {output} -s 2>{log}"
