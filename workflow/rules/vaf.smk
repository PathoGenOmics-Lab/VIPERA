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
        tsv = OUTDIR/f"{OUTPUT_NAME}.masked.tsv"
    output:
        filtered_tsv = OUTDIR/f"{OUTPUT_NAME}.masked.filtered.tsv"
    log:
        LOGDIR / "filter_tsv" / "log.txt"
    script:
        "../scripts/filter_tsv.R"
