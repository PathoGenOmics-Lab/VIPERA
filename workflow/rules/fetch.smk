rule fetch_alignment_reference:
    threads: 1
    conda: "../envs/fetch.yaml"
    output:
        fasta = OUTDIR/"reference.fasta"
    log:
        LOGDIR / "fetch_alignment_reference" / "log.txt"
    shell:
        "esearch -db nucleotide -query {config[ALIGNMENT_REFERENCE]} | efetch -format fasta > {output.fasta} 2> {log}"


rule fetch_mapping_references:
    threads: 1
    conda: "../envs/fetch.yaml"
    output:
        fasta = OUTDIR/"mapping_references.fasta"
    log:
        LOGDIR / "fetch_mapping_references" / "log.txt"
    shell:
        """
        # Create empty file
        echo -n > {output.fasta}
        # Download and concatenate each reference
        echo {config[MAPPING_REFERENCES]} | while read -r ref_id; do
            esearch -db nucleotide -query "$ref_id" | efetch -format fasta >> {output.fasta} 2>> {log}
        done
        """


rule fetch_alignment_annotation:
    threads: 1
    output:
        OUTDIR/"reference.gff3"
    log:
        LOGDIR / "fetch_alignment_annotation" / "log.txt"
    shell:
        "curl 'https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=gff3&id={config[ALIGNMENT_REFERENCE]}' -o {output} -s 2>{log}"
