rule calculate_bam_variants:
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
        # Create output folder (why is this a necessary step?)
        mkdir -p {output.folder}

        echo Calculating variants of sample '{wildcards.sample}'
        freyja variants \
            "{input.bam}" \
            --variants "{output.folder}/{wildcards.sample}_variants.tsv" \
            --depths "{output.folder}/{wildcards.sample}_depth.txt" \
            --ref {input.ref_fasta}

        echo Demixing sample '{wildcards.sample}'
        freyja demix \
            "{output.folder}/{wildcards.sample}_variants.tsv" \
            "{output.folder}/{wildcards.sample}_depth.txt" \
            --eps {params.minimum_abundance} \
            --covcut {params.coverage_cutoff} \
            --output "{output.folder}/{wildcards.sample}_demixed.tsv"
        """
