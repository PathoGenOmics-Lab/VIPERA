rule read_bam_refs:
    threads: 1
    shadow: "minimal"
    conda: "../envs/var_calling.yaml"
    input:
        iter_files("bam")
    output:
        temp("<results>/<dataset>/bam_ids.txt")
    log:
        LOGDIR / "read_bam_refs" / "log.txt"
    shell:
        """
        for bam_file in {input:q}; do
            samtools view -H "$bam_file" | grep ^@SQ | cut -d"\t" -f2 | sed 's/SN://g' >> ids.txt 2>> {log}
        done
        sort ids.txt | uniq > {output}
        """


rule rename_fastas:
    input:
        fasta = get_input_fasta
    output:
        renamed = temp("<results>/<dataset>/renamed/{sample}.fasta")
    log:
        LOGDIR / "rename_fastas" / "{sample}.log.txt"
    shell:
        "sed 's/>.*/>'{wildcards.sample}'/g' {input.fasta} > {output.renamed} 2> {log}"


rule concat_fasta:
    threads: 1
    shadow: "shallow"
    input:
        expand("<results>/<dataset>/renamed/{sample}.fasta", sample = iter_samples())
    output:
        fasta = "<results>/<dataset>/sequences.fasta"
    log:
        LOGDIR / "concat_fasta" / "log.txt"
    shell:
        "cat {input} > {output.fasta} 2> {log}"


rule align_fasta:
    threads: 32
    shadow: "shallow"
    conda: "../envs/nextalign.yaml"
    input:
        ref_fasta = "<results>/<dataset>/reference.fasta",
        fasta = "<results>/<dataset>/sequences.fasta"
    output:
        folder = directory("<results>/<dataset>/nextalign"),
        fasta = "<results>/<dataset>/nextalign/sequences.aligned.fasta"
    log:
        LOGDIR / "align_fasta" / "log.txt"
    shell:
        "nextalign run -j {threads} -O {output.folder} -o {output.fasta} -n sequences --include-reference -r {input.ref_fasta} {input.fasta} >{log} 2>&1"


rule mask_alignment:
    threads: 1
    conda: "../envs/biopython.yaml"
    params:
        remove_sites = False,
        mask_character = "N",
        mask_class = ["mask"]
    input:
        fasta = "<results>/<dataset>/nextalign/sequences.aligned.fasta",
        ref_fasta = "<results>/<dataset>/reference.fasta",
        vcf = lambda wildcards: select_problematic_vcf()
    output:
        fasta = "<results>/<dataset>/aligned.masked.fasta"
    log:
        LOGDIR / "mask_alignment" / "log.txt"
    script:
        "../scripts/mask_aln.py"
