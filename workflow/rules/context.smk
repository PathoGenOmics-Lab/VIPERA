rule download_context:
    threads: 1
    shadow: "shallow"
    conda: "../envs/gisaidr.yaml"
    input:
        metadata = config["METADATA"],
        pango_report = OUTDIR/f"{OUTPUT_NAME}.lineage_report.csv"
    params:
        gisaid_creds = config["GISAID"]["CREDENTIALS"],
        date_window_span = 0.95,
        date_window_paddding_days = 14,
        date_column = config["GISAID"]["DATE_COLUMN"],
        location_column = config["GISAID"]["LOCATION_COLUMN"],
        samples_gisaid_accession_column = config["GISAID"]["ACCESSION_COLUMN"],
        context_gisaid_accession_column = "accession_id",
        host = "Human",
        exclude_low_coverage = True,
        complete = True,
        collection_date_complete = True,
        high_coverage = False,
        min_sleep = 1,
        max_sleep = 3,
        chunk_length = 3000,
        min_theoretical_combinations = config["DIVERSITY_REPS"]
    output:
        fasta = temp(OUTDIR/"context"/"sequences.fasta"),
        metadata = temp(OUTDIR/"context"/"metadata.csv"),
        duplicate_accids = OUTDIR/"context"/"duplicate_accession_ids.txt"
    log:
        LOGDIR / "download_context" / "log.txt"
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
    log:
        LOGDIR / "align_context" / "log.txt"
    shell:
        "nextalign run -j {threads} -O {output.folder} -o {output.fasta} -n {params.name} --include-reference -r {input.ref_fasta} {input.fasta} >{log} 2>&1"


rule mask_context:
    threads: 1
    conda: "../envs/biopython.yaml"
    params:
        remove_sites = False,
        mask_character = "N",
        mask_class = ["mask"]
    input:
        fasta = OUTDIR/"context"/"nextalign"/"context_sequences.aligned.fasta",
        ref_fasta = OUTDIR/"reference.fasta",
        vcf = OUTDIR/"problematic_sites.vcf"
    output:
        fasta = OUTDIR/"context"/"nextalign"/"context_sequences.aligned.masked.fasta"
    log:
        LOGDIR / "mask_context" / "log.txt"
    script:
        "../scripts/mask-aln.py"


rule ml_context_tree:
    threads: 4
    conda: "../envs/iqtree.yaml"
    shadow: "shallow"
    params:
        seqtype = "DNA",
        name = OUTPUT_NAME,
        ufboot = config["UFBOOT_REPS"],
        alrt = config["SHALRT_REPS"],
        outgroup = config["ALIGNMENT_REFERENCE"],
        model = config["TREE_MODEL"]
    input:
        fasta = OUTDIR/"nextalign"/f"{OUTPUT_NAME}.aligned.masked.fasta",
        outgroup_aln = OUTDIR/"context"/"nextalign"/"context_sequences.aligned.masked.fasta"
    output:
        folder = directory(OUTDIR/"tree_context"),
        ml = OUTDIR/f"tree_context/{OUTPUT_NAME}.treefile"
    log:
        LOGDIR / "ml_context_tree" / "log.txt"
    shell:
        """
        exec >{log}                                                                    
        exec 2>&1
        
        awk '/^>/{{p=seen[$0]++}}!p' {input.fasta} {input.outgroup_aln} > aln.fasta
        mkdir -p {output.folder}
        iqtree2 \
            -B {params.ufboot} -alrt {params.alrt} \
            -o {params.outgroup} -T AUTO --threads-max {threads} -s aln.fasta \
            --seqtype {params.seqtype} -m {params.model} --prefix {output.folder}/{params.name}
        """
