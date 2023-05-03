# pipeline para variant calling de https://github.com/grubaughlab/2022_paper_chronic_infection/blob/main/variant_calls.sh 
# No se porque no es --ploidy 1 

rule snps_to_ancestor:
    threads: 1
    shadow: "shallow"
    conda: "../envs/var_calling.yaml"
    params:
        max_depth = 1000000,
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
        
        ID=`echo {input.bam} | grep -o -E COV[0-9]{{6}}`
        samtools mpileup \
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

        
         sed 's/MN908947.3/'$ID'/g' {wildcards.sample}.tsv > {output.tsv}
        
        
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

