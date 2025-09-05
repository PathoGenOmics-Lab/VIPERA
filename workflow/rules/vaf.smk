rule snps_to_ancestor:
    threads: 1
    retries: 3
    shadow: "minimal"
    conda: "../envs/var_calling.yaml"
    params:
        max_depth = config["VC"]["MAX_DEPTH"],
        min_quality = config["VC"]["MIN_QUALITY"],
        ivar_quality = config["VC"]["IVAR_QUALITY"],
        ivar_freq = config["VC"]["IVAR_FREQ"],
        ivar_depth = config["VC"]["IVAR_DEPTH"]
    input:
        reference_fasta = OUTDIR/f"{OUTPUT_NAME}.ancestor.fasta",
        bam = get_input_bam,
        gff = OUTDIR/"reference.gff3"
    output:
        tsv = temp(OUTDIR/"{sample}.tsv")
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
            -d {params.max_depth} \
            --count-orphans \
            --no-BAQ \
            -Q {params.min_quality} \
            -f renamed_reference.fasta \
            {input.bam} \
            | ivar variants \
                -p {wildcards.sample} \
                -q {params.ivar_quality} \
                -t {params.ivar_freq} \
                -m {params.ivar_depth} \
                -g {input.gff} \
                -r renamed_reference.fasta

        sed 's/'$ref'/'{wildcards.sample}'/g' {wildcards.sample}.tsv | cat > {output.tsv}
        """


rule annotation:
    threads:1
    conda: "../envs/biopython.yaml"
    shadow: "shallow"
    input:
        gb = OUTDIR/"reference.gb",
        ref = OUTDIR/"reference.fasta",
        features = Path(config["FEATURES_JSON"]).resolve()
    output:
        df = temp(OUTDIR/"annotation.csv")
    log:
        LOGDIR / "annotation" / "log.txt"
    script:
        "../scripts/report/get_annotation.py"


rule format_tsv:
    threads:1
    shadow: "shallow"
    input:
        expand(OUTDIR/"{sample}.tsv", sample = iter_samples())
    output:
        tsv = OUTDIR/f"{OUTPUT_NAME}.tsv"
    log:
        LOGDIR / "format_tsv" / "log.txt"
    shell:
        """
        path=`echo {input} | awk '{{print $1}}'`
        grep "^REGION" "$path" > header
        for tsv in {input}; do
            tail -n +2 "$tsv"  >> body
        done
        cat header body > "{output.tsv}"
        rm header
        """


rule mask_tsv:
    threads: 1
    conda: "../envs/biopython.yaml"
    params:
         mask_class = ["mask"]
    input: 
        tsv = OUTDIR/f"{OUTPUT_NAME}.tsv",
        vcf = lambda wildcards: select_problematic_vcf()
    output:
        masked_tsv = temp(OUTDIR/f"{OUTPUT_NAME}.masked.tsv")
    log:
        LOGDIR / "mask_tsv" / "log.txt"
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
        tsv = OUTDIR/f"{OUTPUT_NAME}.masked.tsv",
        annotation = OUTDIR/"annotation.csv"
    output:
        filtered_tsv = temp(OUTDIR/f"{OUTPUT_NAME}.masked.prefiltered.tsv")
    log:
        LOGDIR / "filter_tsv" / "log.txt"
    script:
        "../scripts/filter_tsv.R"


rule tsv_to_vcf:
    threads: 1
    conda: "../envs/biopython.yaml"
    input: 
        tsv = OUTDIR/f"{OUTPUT_NAME}.masked.prefiltered.tsv",
    output:
        vcf = temp(OUTDIR/f"{OUTPUT_NAME}.vcf")
    log:
        LOGDIR / "tsv_to_vcf" / "log.txt"
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
        vcf = OUTDIR/f"{OUTPUT_NAME}.vcf"
    output:
        ann_vcf = temp(OUTDIR/f"{OUTPUT_NAME}.annotated.vcf")
    log:
        LOGDIR / "variants_effect" / "log.txt"
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

        snpEff eff -dataDir {params.snpeff_data_dir} -noStats {params.ref_name} {input.vcf} > {output.ann_vcf}
        """


rule extract_vcf_fields:
    threads: 1
    conda: "../envs/snpeff.yaml"
    params:
        extract_columns = [
            "CHROM", "POS", "REF", "ALT",
            '"ANN[*].IMPACT"', '"ANN[*].BIOTYPE"',
            '"ANN[*].GENE"', '"ANN[*].GENEID"', '"ANN[*].FEATURE"', '"ANN[*].HGVS_P"', '"ANN[*].HGVS_C"'
        ],
        sep = ","
    input:
        vcf = OUTDIR/f"{OUTPUT_NAME}.annotated.vcf"
    output:
        tsv = OUTDIR/f"{OUTPUT_NAME}.vcf_fields.tsv"
    log:
        LOGDIR / "tsv_to_vcf" / "log.txt"
    shell:
        'SnpSift extractFields -s {params.sep:q} {input.vcf:q} {params.extract_columns} >{output.tsv:q} 2>{log:q}'


rule vcf_to_tsv:
    threads: 1
    conda: "../envs/renv.yaml"
    input:
        ann_vcf = OUTDIR/f"{OUTPUT_NAME}.annotated.vcf",
        pre_tsv = OUTDIR/f"{OUTPUT_NAME}.masked.prefiltered.tsv"
    output:
        tsv = OUTDIR/f"{OUTPUT_NAME}.masked.filtered.tsv"
    log:
        LOGDIR / "vcf_to_tsv" / "log.txt"
    script:
        "../scripts/vcf_to_tsv.R"
