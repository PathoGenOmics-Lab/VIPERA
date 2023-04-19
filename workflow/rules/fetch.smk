rule fetch_alignment_reference:
    threads: 1
    conda: "../envs/fetch.yaml"
    output:
        fasta = OUTDIR/"reference.fasta"
    shell:
        "esearch -db nucleotide -query {config[ALIGNMENT_REFERENCE]} | efetch -format fasta > {output.fasta}"


rule fetch_mapping_references:
    threads: 1
    conda: "../envs/fetch.yaml"
    output:
        fasta = OUTDIR/"mapping_references.fasta"
    shell:
        """
        # Create empty file
        echo -n > {output.fasta}
        # Download and concatenate each reference
        echo {config[MAPPING_REFERENCES]} | while read -r ref_id; do
            esearch -db nucleotide -query "$ref_id" | efetch -format fasta >> {output.fasta}
        done
        """
