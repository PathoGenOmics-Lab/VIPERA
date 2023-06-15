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
        folder = directory(OUTDIR/"demixing"/"{sample}/"),  # TODO: add explicit tsv name
        depth_file = OUTDIR/"demixing"/"{sample}/{sample}_depth.txt",
        variants_file = OUTDIR/"demixing"/"{sample}/{sample}_variants.tsv"
    shell:
        """
        # Create output folder (why is this a necessary step?)
        mkdir -p {output.folder}

        echo Calculating variants of sample '{wildcards.sample}'
        freyja variants \
            "{input.bam}" \
            --variants {output.variants_file} \
            --depths {output.depth_file} \
            --ref {input.ref_fasta}

        echo Demixing sample '{wildcards.sample}'
        freyja demix \
            {output.variants_file} \
            {output.depth_file} \
            --eps {params.minimum_abundance} \
            --covcut {params.coverage_cutoff} \
            --output "{output.folder}/{wildcards.sample}_demixed.tsv"
        """


rule summarise_demixing:
    threads: 1
    conda: "../envs/renv.yaml"
    shadow: "shallow"
    params:
        demix_folder = Path(OUTDIR/"demixing")
    input:
        directories = expand(OUTDIR/"demixing"/"{sample}", sample=iter_samples())
    output:
        summary_df = OUTDIR/"summary_freyja_demixing.csv"
    script: 
        "../scripts/summary_demixing.R"
