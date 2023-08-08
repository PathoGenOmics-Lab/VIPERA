rule reconstruct_ancestral_sequence:
    threads: 4
    conda: "../envs/iqtree.yaml"
    params:
        seqtype = "DNA",
        name = OUTPUT_NAME,
        etc = ETC_TREE_PARAMS
    input:
        fasta = OUTDIR/"nextalign"/f"{OUTPUT_NAME}.aligned.masked.fasta"
    output:
        folder = directory(OUTDIR/"tree"),
        state_file = OUTDIR/"tree"/f"{OUTPUT_NAME}.state"
    log:
        LOGDIR / "reconstruct_ancestral_sequence" / "log.txt"
    shell:
        """
        mkdir -p {output.folder}
        iqtree2 \
            {params.etc} -asr \
            -o {config[ALIGNMENT_REFERENCE]} -T AUTO --threads-max {threads} -s {input.fasta} \
            --seqtype {params.seqtype} -m {config[TREE_MODEL]} --prefix {output.folder}/{params.name} >{log} 2>&1
        """


rule ancestor_fasta:
    threads: 1
    conda: "../envs/biopython.yaml"
    params:
        node_id = "Node1",
        indeterminate_char = "N",
        name = OUTPUT_NAME
    input:
        state_file = OUTDIR/"tree"/f"{OUTPUT_NAME}.state"
    output:
        fasta = report(OUTDIR/f"{OUTPUT_NAME}.ancestor.fasta")
    log:
        LOGDIR / "ancestor_fasta" / "log.txt"
    script:
        "../scripts/ancestor_fasta.py"
