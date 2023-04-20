rule concat_fasta:
    threads: 1
    shadow: "shallow"
    input:
        iter_files_in_path(FASTA_FOLDER)
    output:
        fasta = OUTDIR/f"{OUTPUT_NAME}.fasta"
    shell:
        "cat {input} > {output.fasta}"


rule align_fasta:
    threads: 32
    shadow: "shallow"
    conda: "../envs/nextalign.yaml"
    params:
        name = OUTPUT_NAME
    input:
        ref_fasta = OUTDIR/"reference.fasta",
        fasta = OUTDIR/f"{OUTPUT_NAME}.fasta"
    output:
        folder = directory(OUTDIR/"nextalign"),
        fasta = OUTDIR/"nextalign"/f"{OUTPUT_NAME}.aligned.fasta"
    shell:
        "nextalign run -j {threads} -O {output.folder} -o {output.fasta} -n {params.name} --include-reference -r {input.ref_fasta} {input.fasta}"


rule mask_alignment:
    threads: 1
    conda: "../envs/biopython.yaml"
    params:
        remove_sites = False,
        mask_character = "N",
        mask_class = ["mask"]
    input:
        fasta = OUTDIR/"nextalign"/f"{OUTPUT_NAME}.aligned.fasta",
        ref_fasta = OUTDIR/"reference.fasta"
    output:
        fasta = OUTDIR/"nextalign"/f"{OUTPUT_NAME}.aligned.masked.fasta"
    script:
        "../scripts/mask-aln.py"
