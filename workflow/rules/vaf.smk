

rule snps_to_ancestor:
    threads: 1
    shadow: "shallow"
    conda: "../envs/var_calling.yaml"
    params:
        max_depth = 0,
        min_quality = 0,
        ivar_quality = 20,
        ivar_freq = 0.05,
        ivar_depth = 30,
        gff = config["ANNOTATION_GFF"]
    input:
        reference_fasta = OUTDIR/f"{OUTPUT_NAME}.ancestor.fasta",
        bam = BAM_FOLDER/"{sample}.trim.sort.bam"
    output:
        tsv = OUTDIR/"{sample}.tsv"
    shell:
        """
        set +o pipefail

        sed 's/>.*/>MN908947.3/g' {input.reference_fasta} > renamed_reference.fasta
        
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

        
         sed 's/MN908947.3/'{wildcards.sample}'/g' {wildcards.sample}.tsv | cat > {output.tsv}
        
        
        """
rule format_tsv:
    threads:1
    shadow: "shallow"
    input:
        expand(OUTDIR/"{sample}.tsv", sample = iter_samples_in_path(BAM_FOLDER))
    output:
        tsv = OUTDIR/f"{OUTPUT_NAME}.tsv"
    shell:
        """
    

        path=`echo {input} | awk '{{print $1}}'`
        grep REGION $path > header
    
        for tsv in {input}; do
            tail -n +2 $tsv  >> body
            rm $tsv
        done

        cat header body > {output.tsv}
        rm header
        exit 0
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

