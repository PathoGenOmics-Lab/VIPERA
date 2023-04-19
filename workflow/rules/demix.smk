rule calculate_mixtures:
    threads: 1
    conda: "../envs/freyja.yaml"
    shadow: "shallow"
    input:
        bam = BAM_FOLDER/"{sample}.trim.sort.bam",
        ref_fasta = OUTDIR/"mapping_references.fasta"
    params:
        coverage_cutoff = 30,
        minimum_abundance = 0.0001
    output:
        folder = directory(OUTDIR/"demixing"/"{sample}")  # TODO: add explicit tsv name
    shell:
        """
        freyja variants \
            "{input.bam}" \
            --variants "{output.folder}/{wildcards.sample}_variants" \
            --depths "{output.folder}/{wildcards.sample}_depth.txt" \
            --ref {input.ref_fasta}

        freyja demix \
            "{output.folder}/{wildcards.sample}_variants.tsv" \
            "{output.folder}/{wildcards.sample}_depth.txt" \
            --eps {params.minimum_abundance} \
            --covcut {params.coverage_cutoff} \
            --output "{output.folder}/{wildcards.sample}_demixed.tsv"
        """
