rule reconstruct_ancestral_sequence:
    threads: 4
    conda: "../envs/iqtree.yaml"
    params:
        seed = 7291,
        seqtype = "DNA",
        outgroup = config["ALIGNMENT_REFERENCE"],
        model = config["TREE_MODEL"]
    input:
        fasta = "<results>/<dataset>/aligned.masked.fasta"
    output:
        folder = directory("<results>/<dataset>/tree"),
        state_file = "<results>/<dataset>/tree/asr.state"
    log:
        "<logs>/<dataset>/reconstruct_ancestral_sequence/log.txt"
    shell:
        "mkdir -p {output.folder} && "
        "iqtree2 -seed {params.seed} "
        "-asr "
        "-o {params.outgroup} -T AUTO --threads-max {threads} -s {input.fasta} "
        "--seqtype {params.seqtype} -m {params.model} --prefix {output.folder}/asr >{log} 2>&1"


rule ancestor_fasta:
    threads: 1
    conda: "../envs/biopython.yaml"
    params:
        node_id = "Node1",
        indeterminate_char = "N",
        name = "case_ancestor",
    input:
        state_file = "<results>/<dataset>/tree/asr.state"
    output:
        fasta = report("<results>/<dataset>/ancestor.fasta")
    log:
        "<logs>/<dataset>/ancestor_fasta/log.txt"
    script:
        "../scripts/ancestor_fasta.py"
