rule download_context:
    threads: 1
    shadow: "shallow"
    conda: "../envs/gisaidr.yaml"
    input:
        metadata = config["METADATA"]
    params:
        gisaid_creds = config["GISAID_YAML"],
        date_window_span = 0.95,
        date_window_paddding_days = 14,
        date_column = "CollectionDate",
        location_column = "ResidenceCity",
        samples_gisaid_accession_column = "GISAIDEPI",
        context_gisaid_accession_column = "accession_id",
        host = "Human",
        exclude_low_coverage = True,
        complete = True,
        collection_date_complete = True,
        high_coverage = False,
        min_sleep = 1,
        max_sleep = 3,
        chunk_length = 3000,
        min_theoretical_combinations = config["DIVERSITY_BOOTSTRAP_REPS"]
    output:
        fasta = OUTDIR/"context"/"sequences.fasta",
        metadata = OUTDIR/"context"/"metadata.csv"
    script:
        "../scripts/download_context.R"


rule align_context:
    threads: 32
    shadow: "shallow"
    conda: "../envs/nextalign.yaml"
    params:
        name = "context_sequences"
    input:
        ref_fasta = OUTDIR/"reference.fasta",
        fasta = CONTEXT_FASTA
    output:
        folder = directory(OUTDIR/"context"/"nextalign"),
        fasta = OUTDIR/"context"/"nextalign"/"context_sequences.aligned.fasta"
    shell:
        "nextalign run -j {threads} -O {output.folder} -o {output.fasta} -n {params.name} --include-reference -r {input.ref_fasta} {input.fasta}"


rule mask_context:
    threads: 1
    conda: "../envs/biopython.yaml"
    params:
        remove_sites = False,
        mask_character = "N",
        mask_class = ["mask"]
    input:
        fasta = OUTDIR/"context"/"nextalign"/"context_sequences.aligned.fasta",
        ref_fasta = OUTDIR/"reference.fasta"
    output:
        fasta = OUTDIR/"context"/"nextalign"/"context_sequences.aligned.masked.fasta"
    script:
        "../scripts/mask-aln.py"
