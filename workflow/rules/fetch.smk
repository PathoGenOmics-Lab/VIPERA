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


rule fetch_reference_gb:
    threads: 1
    conda: "../envs/fetch.yaml"
    params:
        ref = config["ALIGNMENT_REFERENCE"],
        database = "nucleotide",
        format = "gb"
    output:
        fasta = OUTDIR/"reference.gb"
    log:
        LOGDIR / "fetch_reference_gb" / "log.txt"
    shell:
        "esearch -db {params.database} -query {params.ref} | efetch -format {params.format} > {output.fasta} 2> {log}"



rule fetch_mapping_references:
    threads: 1
    conda: "../envs/fetch.yaml"
    input:
        OUTDIR / "bam_ids.txt"
    output:
        fasta = select_mapping_references_fasta()
    log:
        LOGDIR / "fetch_mapping_references" / "log.txt"
    shell:
        """
        cat {input} | while read ref_id || [[ -n $ref_id ]]; do
            esearch -db nucleotide -query "$ref_id" | efetch -format fasta >> {output.fasta} 2>> {log}
        done
        """


rule fetch_alignment_annotation:
    threads: 1
    params:
        ref = config["ALIGNMENT_REFERENCE"]
    output:
        temp(OUTDIR/"reference.gff3")
    log:
        LOGDIR / "fetch_alignment_annotation" / "log.txt"
    shell:
        "curl 'https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=gff3&id={params.ref}' -o {output} -s 2>{log}"


rule fetch_problematic_vcf:
    threads: 1
    params:
        url = config["PROBLEMATIC_VCF"]
    log:
        LOGDIR / "fetch_problematic_vcf" / "log.txt"
    output:
        select_problematic_vcf()
    shell:
        "curl {params.url} -o {output} -s 2> {log}"
