rule demix_preprocessing:
    threads: 1
    conda: "../envs/freyja.yaml"
    shadow: "minimal"
    input:
        bam = get_input_bam,
        ref_fasta = MAPPING_REFERENCES_FASTA
    params:
        minq = config["DEMIX"]["MIN_QUALITY"]
    output:
        depth_file = OUTDIR/"demixing"/"{sample}/{sample}_depth.txt",
        variants_file = OUTDIR/"demixing"/"{sample}/{sample}_variants.tsv"
    log:
        LOGDIR / "demix_preprocessing" / "{sample}.log.txt"
    shell:
        """
        freyja variants \
            "{input.bam}" \
            --variants {output.variants_file} \
            --depths {output.depth_file} \
            --minq {params.minq} \
            --ref {input.ref_fasta} >{log} 2>&1
        """


rule demix:
    threads: 1
    conda: "../envs/freyja.yaml"
    shadow: "minimal"
    input:
        depth_file = OUTDIR/"demixing"/"{sample}/{sample}_depth.txt",
        variants_file = OUTDIR/"demixing"/"{sample}/{sample}_variants.tsv"
    params:
        coverage_cutoff = config["DEMIX"]["COV_CUTOFF"],
        minimum_abundance = config["DEMIX"]["MIN_ABUNDANCE"]
    output:
        demix_file = OUTDIR/"demixing"/"{sample}/{sample}_demixed.tsv"
    log:
        LOGDIR / "demix" / "{sample}.log.txt"
    shell:
        """
        freyja demix \
            "{input.variants_file}" \
            "{input.depth_file}" \
            --eps {params.minimum_abundance} \
            --covcut {params.coverage_cutoff} \
            --output {output.demix_file} >{log} 2>&1
        """


rule summarise_demixing:
    threads: 1
    conda: "../envs/renv.yaml"
    shadow: "shallow"
    input:
        tables = expand(OUTDIR/"demixing"/"{sample}/{sample}_demixed.tsv", sample=iter_samples())
    output:
        summary_df = report(OUTDIR/"summary_freyja_demixing.csv")
    log:
        LOGDIR / "summarise_demixing" / "log.txt"
    script: 
        "../scripts/summary_demixing.R"
