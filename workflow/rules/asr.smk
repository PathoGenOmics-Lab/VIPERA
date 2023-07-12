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
    shell:
        """
        mkdir -p {output.folder}
        iqtree2 \
            {params.etc} -asr \
            -o {config[ALIGNMENT_REFERENCE]} -T AUTO --threads-max {threads} -s {input.fasta} \
            --seqtype {params.seqtype} -m {config[TREE_MODEL]} --prefix {output.folder}/{params.name}
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
    script:
        "../scripts/ancestor_fasta.py"


rule ml_context_tree:
    threads: 4
    conda: "../envs/iqtree.yaml"
    shadow: "shallow"
    params:
        seqtype = "DNA",
        name = OUTPUT_NAME,
        etc = ETC_TREE_PARAMS
    input:
        fasta = OUTDIR/"nextalign"/f"{OUTPUT_NAME}.aligned.masked.fasta",
        outgroup_aln = CONTEXT_FASTA
    output:
        folder = directory(OUTDIR/"tree_context"),
        state_file = OUTDIR/"tree_context"/f"{OUTPUT_NAME}.state",
        treefile = OUTDIR/"tree_context"/f"{OUTPUT_NAME}.treefile"
    shell:
        """
        awk '/^>/{{p=seen[$0]++}}!p' {input.fasta} {input.outgroup_aln} > aln.fasta
        mkdir -p {output.folder}
        iqtree2 \
            {params.etc} -asr \
            -o {config[ALIGNMENT_REFERENCE]} -T AUTO --threads-max {threads} -s aln.fasta \
            --seqtype {params.seqtype} -m {config[TREE_MODEL]} --prefix {output.folder}/{params.name}
        """
