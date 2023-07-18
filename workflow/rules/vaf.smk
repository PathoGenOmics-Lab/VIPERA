rule snps_to_ancestor:
    threads: 1
    shadow: "full"
    conda: "../envs/var_calling.yaml"
    params:
        max_depth = config["VC"]["MAX_DEPTH"],
        min_quality = config["VC"]["MIN_QUALITY"],
        ivar_quality = config["VC"]["IVAR_QUALITY"],
        ivar_freq = config["VC"]["IVAR_FREQ"],
        ivar_depth = config["VC"]["IVAR_DEPTH"],
        gff = config["ANNOTATION_GFF"],
        ref_id = config["ANNOTATION_GFF_SEQNAME"]
    input:
        reference_fasta = OUTDIR/f"{OUTPUT_NAME}.ancestor.fasta",
        bam = get_input_bam
    output:
        tsv = temp(OUTDIR/"{sample}.tsv")
    shell:
        """
        sed 's/^>.*/>{params.ref_id}/g' {input.reference_fasta} > renamed_reference.fasta

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
                -g {params.gff} \
                -r renamed_reference.fasta
        
        sed 's/{params.ref_id}/'{wildcards.sample}'/g' {wildcards.sample}.tsv | cat > {output.tsv}
        """


rule format_tsv:
    threads:1
    shadow: "shallow"
    input:
        expand(OUTDIR/"{sample}.tsv", sample = iter_samples())
    output:
        tsv = OUTDIR/f"{OUTPUT_NAME}.tsv"
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
        tsv = OUTDIR/f"{OUTPUT_NAME}.tsv"
    output:
        masked_tsv = OUTDIR/f"{OUTPUT_NAME}.masked.tsv"
    script:
        "../scripts/mask_tsv.py"


rule filter_tsv:
    threads: 1
    conda: "../envs/renv.yaml"
    params:
         
    input: 
        tsv = OUTDIR/f"{OUTPUT_NAME}.masked.tsv"
    output:
        filtered_tsv = OUTDIR/f"{OUTPUT_NAME}.masked.filtered.tsv"
    script:
        "../scripts/filter_tsv.R"
