# pipeline para variant calling de https://github.com/grubaughlab/2022_paper_chronic_infection/blob/main/variant_calls.sh 
# No se porque no es --ploidy 1 

rule snps_to_ancestor:
    threads: 1
    shadow: "shallow"
    conda: "../envs/bcftools.yaml"
    params:
        max_depth = 1000000,
        min_quality = 0
    input:
        reference_fasta = OUTDIR/f"{OUTPUT_NAME}.ancestor.fasta",
        bam = BAM_FOLDER/"{sample}.trim.sort.bam"
    output:
        vcf = OUTDIR/"{sample}.vcf"
    shell:
        """
        set +o pipefail

        sed 's/>.*/>MN908947.3/g' {input.reference_fasta} > renamed_reference.fasta
        
        ID=`echo {input.bam} | grep -o -E COV[0-9]{{6}}`
        bcftools mpileup \
         -Ou -d {params.max_depth} \
         --max-idepth {params.max_depth} \
         --per-sample-mF \
         --count-orphans \
          --no-BAQ \
          --annotate AD,ADF,ADR,DP,SP,AD,ADF,ADR \
          -Q {params.min_quality} \
          -f renamed_reference.fasta \
           {input.bam} \
           | bcftools +fill-tags - \
            | bcftools call \
            --prior-freqs AN,AC \
            --variants-only \
            --multiallelic-caller -Ov \
            | sed 's/$/\t'$ID'/g' > {output.vcf}
        
        
        """
rule format_vcf:
    threads:1
    shadow: "shallow"
    input:
        expand(OUTDIR/"{sample}.vcf", sample = iter_samples_in_path(BAM_FOLDER))
    output:
        vcf = OUTDIR/f"{OUTPUT_NAME}.vcf"
    shell:
        """
        set +o pipefail

        path=`echo {input} | awk '{{print $1}}'`
        grep "##" $path | sed 's/\tCOV......//g' > header_1
        grep "#" $path | grep -v "##" | sed 's/[^\t][^\t]*$/COV/g' > header_2
        cat header_1 header_2 > header
        for vcf in {input}; do
            grep -v "#" $vcf | cat >> body
            rm $vcf
        done

        cat header body > {output.vcf}
        rm header
        exit 0
        """

rule mask_vcf:
    threads: 1
    conda: "../envs/biopython.yaml"
    params:
         mask_class = ["mask"]
    input: 
        vcf = OUTDIR/f"{OUTPUT_NAME}.vcf"
    output:
        masked_vcf = OUTDIR/f"{OUTPUT_NAME}.masked.vcf"
    script:
        "../scripts/mask_vcf.py"

rule annotate_vcf:
    threads: 1
    conda: "../envs/vcf-annotator.yaml"
    params:
        gb = config["ANNOTATION_GB"]
    input:
        vcf =  OUTDIR/f"{OUTPUT_NAME}.masked.vcf"
    output:
        annotated = OUTDIR/f"{OUTPUT_NAME}.masked.annotated.vcf"
    shell:
        """
        vcf-annotator --output {output.annotated} {input.vcf} {params.gb}
        """
rule compute_vaf:
    threads: 1
    conda: "../envs/renv.yaml"
    input:
        vcf = OUTDIR/f"{OUTPUT_NAME}.masked.annotated.vcf"
    output:
        vaf_file = OUTDIR/f"{OUTPUT_NAME}.masked.annotated.vaf.tsv"
    script:
        "../scripts/vaf_from_vcf.R"
