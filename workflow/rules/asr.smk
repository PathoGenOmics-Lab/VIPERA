rule reconstruct_ancestral_sequence:
    threads: 4
    conda: "../envs/iqtree.yaml"
    params:
        seed = 7291,
        seqtype = "DNA",
        name = OUTPUT_NAME,
        outgroup = config["ALIGNMENT_REFERENCE"],
        model = config["TREE_MODEL"]
    input:
        fasta = OUTDIR/"nextalign"/f"{OUTPUT_NAME}.aligned.masked.fasta"
    output:
        folder = directory(OUTDIR/"tree"),
        state_file = OUTDIR/"tree"/f"{OUTPUT_NAME}.state"
    log:
        LOGDIR / "reconstruct_ancestral_sequence" / "log.txt"
    shell:
        "mkdir -p {output.folder} && "
        "iqtree2 -seed {params.seed} "
            "-asr "
            "-o {params.outgroup} -T AUTO --threads-max {threads} -s {input.fasta} "
            "--seqtype {params.seqtype} -m {params.model} --prefix {output.folder}/{params.name} >{log} 2>&1"


rule ancestor_fasta:
    threads: 1
    conda: "../envs/biopython.yaml"
    params:
        node_id = "Node1",
        indeterminate_char = "N",
        name = "case_ancestor",
    input:
        state_file = OUTDIR/"tree"/f"{OUTPUT_NAME}.state"
    output:
        fasta = report(OUTDIR/f"{OUTPUT_NAME}.ancestor.fasta")
    log:
        LOGDIR / "ancestor_fasta" / "log.txt"
    script:
        "../scripts/ancestor_fasta.py"
