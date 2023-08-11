rule rename_fastas:
    input:
        fasta = get_input_fasta
    output:
        renamed = temp(OUTDIR/"renamed.{sample}.fasta")
    log:
        LOGDIR / "rename_fastas" / "{sample}.log.txt"
    shell:
        "sed 's/>.*/>'{wildcards.sample}'/g' {input.fasta} > {output.renamed} 2> {log}"


rule concat_fasta:
    threads: 1
    shadow: "shallow"
    input:
        expand(OUTDIR/"renamed.{sample}.fasta", sample = iter_samples())
    output:
        fasta = OUTDIR/f"{OUTPUT_NAME}.fasta"
    log:
        LOGDIR / "concat_fasta" / "log.txt"
    shell:
        "cat {input} > {output.fasta} 2> {log}"


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
    log:
        LOGDIR / "align_fasta" / "log.txt"
    shell:
        "nextalign run -j {threads} -O {output.folder} -o {output.fasta} -n {params.name} --include-reference -r {input.ref_fasta} {input.fasta} >{log} 2>&1"


rule mask_alignment:
    threads: 1
    conda: "../envs/biopython.yaml"
    params:
        remove_sites = False,
        mask_character = "N",
        mask_class = ["mask"]
    input:
        fasta = OUTDIR/"nextalign"/f"{OUTPUT_NAME}.aligned.fasta",
        ref_fasta = OUTDIR/"reference.fasta",
        vcf = OUTDIR/"problematic_sites.vcf"
    output:
        fasta = OUTDIR/"nextalign"/f"{OUTPUT_NAME}.aligned.masked.fasta"
    log:
        LOGDIR / "mask_alignment" / "log.txt"
    script:
        "../scripts/mask-aln.py"
