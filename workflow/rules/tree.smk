rule infer_tree:
    threads: 4
    conda: "../envs/iqtree.yaml"
    params:
        seqtype = "DNA",
        name = OUTPUT_NAME,
        etc = ETC_TREE_PARAMS
    input:
        fasta = OUTDIR/"nextalign"/f"{OUTPUT_NAME}.aligned.masked.fasta"
    output:
        touch(OUTDIR/"tree"/".iqtree_done"),
        folder = directory(OUTDIR/"tree")
    shell:
        """
        mkdir -p {output.folder}
        iqtree2 \
            {params.etc} \
            -o {config[REFSEQ_REFERENCE]} -T AUTO --threads-max {threads} -s {input.fasta} \
            --seqtype {params.seqtype} -m {config[TREE_MODEL]} --prefix {output.folder}/{params.name}
        """
