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
        features = config["FEATURES_JSON"]
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
        vcf = OUTDIR / "problematic_sites.vcf"
    output:
        masked_tsv = temp(OUTDIR/f"{OUTPUT_NAME}.masked.tsv")
    log:
        LOGDIR / "mask_tsv" / "log.txt"
    script:
        "../scripts/mask_tsv.py"


rule filter_tsv:
    threads: 1
    conda: "../envs/renv.yaml"
    input: 
        tsv = OUTDIR/f"{OUTPUT_NAME}.masked.tsv",
        annotation = OUTDIR/"annotation.csv"
    output:
        filtered_tsv = OUTDIR/f"{OUTPUT_NAME}.masked.prefiltered.tsv"
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
        vcf = OUTDIR/f"{OUTPUT_NAME}.vcf"
    log:
        LOGDIR / "tsv_to_vcf" / "log.txt"
    script:
        "../scripts/tsv_to_vcf.py"

rule variants_effect:
    threads: 1
    conda: "../envs/snpeff.yaml"
    params:
        ref_name = config["ALIGNMENT_REFERENCE"]
    input:
        vcf = OUTDIR/f"{OUTPUT_NAME}.vcf"
    output:
        ann_vcf = OUTDIR/f"{OUTPUT_NAME}.annotated.vcf"
    log:
        LOGDIR / "variants_effect" / "log.txt"
    shell:
        """
        exec >{log}                                                                    
        exec 2>&1
        
        snpEff eff {params.ref_name} {input.vcf} > {output.ann_vcf} || true
        rm snpEff_genes.txt snpEff_summary.html
        """

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


