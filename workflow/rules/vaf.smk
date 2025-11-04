rule snps_to_ancestor:
    threads: 1
    retries: 3
    shadow: "minimal"
    conda: "../envs/var_calling.yaml"
    params:
        mpileup_depth = config["VC"]["MAX_DEPTH"],
        mpileup_quality = 0,
        ivar_quality = config["VC"]["MIN_QUALITY"],
        ivar_freq = config["VC"]["MIN_FREQ"],
        ivar_depth = config["VC"]["MIN_DEPTH"],
    input:
        reference_fasta = OUTDIR/f"{OUTPUT_NAME}.ancestor.fasta",
        bam = get_input_bam,
        gff = OUTDIR/"reference.gff3"
    output:
        tsv = temp(OUTDIR/"vaf"/"{sample}.tsv")
    log:
        LOGDIR / "snps_to_ancestor" / "{sample}.log.txt"
    shell:
        """
        set -e
        exec >{log}
        exec 2>&1

        ref=`samtools view -H {input.bam} | grep ^@SQ | cut -d"\t" -f2 | sed 's/SN://g'`
        echo Reference: $ref
        echo FASTA before:
        grep ">" {input.reference_fasta}
        sed 's/>.*/>'$ref'/g' {input.reference_fasta} > renamed_reference.fasta
        echo FASTA after:
        grep ">" renamed_reference.fasta
        
        echo Starting VC
        samtools mpileup \
            -aa \
            --ignore-overlaps \
            -d {params.mpileup_depth} \
            --count-orphans \
            --no-BAQ \
            -Q {params.mpileup_quality} \
            -f renamed_reference.fasta \
            {input.bam} \
            | ivar variants \
                -p variants \
                -q {params.ivar_quality} \
                -t {params.ivar_freq} \
                -m {params.ivar_depth} \
                -g {input.gff} \
                -r renamed_reference.fasta
        mv variants.tsv {output.tsv:q}
        """


rule mask_tsv:
    threads: 1
    conda: "../envs/biopython.yaml"
    params:
         mask_class = ["mask"]
    input: 
        tsv = OUTDIR/"vaf"/"{sample}.tsv",
        vcf = lambda wildcards: select_problematic_vcf()
    output:
        masked_tsv = temp(OUTDIR/"vaf"/"{sample}.masked.tsv")
    log:
        LOGDIR / "mask_tsv" / "{sample}.log.txt"
    script:
        "../scripts/mask_tsv.py"


rule filter_tsv:
    threads: 1
    conda: "../envs/renv.yaml"
    params:
        min_depth = 20,
        min_alt_rv = 2,
        min_alt_dp = 2,
    input: 
        tsv = OUTDIR/"vaf"/"{sample}.masked.tsv"
    output:
        filtered_tsv = temp(OUTDIR/"vaf"/"{sample}.masked.prefiltered.tsv")
    log:
        LOGDIR / "filter_tsv" / "{sample}.log.txt"
    script:
        "../scripts/filter_tsv.R"


rule tsv_to_vcf:
    threads: 1
    conda: "../envs/biopython.yaml"
    params:
        ref_name = config["ALIGNMENT_REFERENCE"],
    input: 
        tsv = OUTDIR/"vaf"/"{sample}.masked.prefiltered.tsv",
    output:
        vcf = temp(OUTDIR/"vaf"/"{sample}.vcf")
    log:
        LOGDIR / "tsv_to_vcf" / "{sample}.log.txt"
    script:
        "../scripts/tsv_to_vcf.py"


rule variants_effect:
    threads: 1
    shadow: "minimal"
    conda: "../envs/snpeff.yaml"
    params:
        ref_name = config["ALIGNMENT_REFERENCE"],
        snpeff_data_dir = (BASE_PATH / "config" / "snpeff").resolve()
    input:
        vcf = OUTDIR/"vaf"/"{sample}.vcf"
    output:
        ann_vcf = OUTDIR/"vaf"/"{sample}.annotated.vcf"
    log:
        LOGDIR / "variants_effect" / "{sample}.log.txt"
    retries: 2
    shell:
        """
        exec >{log}
        exec 2>&1

        # Check if snpEff database is available
        if [ -d "{params.snpeff_data_dir}/{params.ref_name}" ]; then
            echo "Using local database at '{params.snpeff_data_dir}'"
        else
            echo "Local database not found at '{params.snpeff_data_dir}', downloading from repository"
        fi

        snpEff eff -dataDir {params.snpeff_data_dir} -noStats {params.ref_name} {input.vcf} >{output.ann_vcf}
        """


rule extract_vcf_fields:
    threads: 1
    conda: "../envs/snpeff.yaml"
    params:
        extract_columns = [f"'{col}'" for col in config["ANNOTATION"]["SNPEFF_COLS"].values()],
        sep = ",",
    input:
        vcf = OUTDIR/"vaf"/"{sample}.annotated.vcf"
    output:
        tsv = OUTDIR/"vaf"/"{sample}.vcf_fields.tsv"
    log:
        LOGDIR / "tsv_to_vcf" / "{sample}.log.txt"
    shell:
        "SnpSift extractFields -e 'NA' -s {params.sep:q} {input.vcf:q} {params.extract_columns} >{output.tsv:q} 2>{log:q}"


rule format_vcf_fields_longer:
    conda: "../envs/renv.yaml"
    params:
        sample = "{sample}",
        colnames_mapping = config["ANNOTATION"]["SNPEFF_COLS"],
        filter_include = config["ANNOTATION"]["FILTER_INCLUDE"],
        filter_exclude = config["ANNOTATION"]["FILTER_EXCLUDE"],
        variant_name_pattern = lambda wildcards: config["ANNOTATION"]["VARIANT_NAME_PATTERN"],  # lambda to deactivate automatic wildcard expansion in pattern
        sep = ",",
    input:
        tsv = OUTDIR/"vaf"/"{sample}.vcf_fields.tsv",
    output:
        tsv = OUTDIR/"vaf"/"{sample}.vcf_fields.longer.tsv",
    log:
        LOGDIR / "format_vcf_fields_longer" / "{sample}.log.txt"
    script:
        "../scripts/format_vcf_fields_longer.R"


rule concat_vcf_fields:
    params:
        sep = "\t",
    input:
        expand(OUTDIR/"vaf"/"{sample}.vcf_fields.longer.tsv", sample=iter_samples()),
    output:
        OUTDIR/f"{OUTPUT_NAME}.vcf_fields.longer.tsv",
    run:
        import pandas as pd
        from functools import reduce
        reduce(
            lambda a, b: pd.concat((a, b), axis="rows", ignore_index=True),
            (pd.read_csv(path, sep=params.sep) for path in input)
        ).to_csv(output[0], sep=params.sep, index=False)


rule merge_annotation:
    threads: 1
    conda: "../envs/renv.yaml"
    params:
        sample = "{sample}",
        ref_name = config["ALIGNMENT_REFERENCE"],
    input:
        tsv = OUTDIR/"vaf"/"{sample}.masked.prefiltered.tsv",
        annot = OUTDIR/"vaf"/"{sample}.vcf_fields.longer.tsv",
    output:
        tsv = OUTDIR/"vaf"/"{sample}.variants.tsv"
    log:
        LOGDIR / "merge_annotation" / "{sample}.log.txt"
    script:
        "../scripts/merge_annotation.R"


use rule concat_vcf_fields as concat_variants with:
    input:
        expand(OUTDIR/"vaf"/"{sample}.variants.tsv", sample=iter_samples()),
    output:
        OUTDIR/f"{OUTPUT_NAME}.variants.tsv",


rule pairwise_trajectory_correlation:
    conda: "../envs/renv.yaml"
    input:
        variants =  OUTDIR/f"{OUTPUT_NAME}.variants.tsv",
        metadata = config["METADATA"],
    output:
        table = report(OUTDIR/"vaf"/"pairwise_trajectory_correlation.csv"),
    log:
        LOGDIR / "pairwise_trajectory_correlation" / "log.txt"
    script:
        "../scripts/pairwise_trajectory_correlation.R"
