rule snps_to_ancestor:
    threads: 1
    shadow: "shallow"
    conda: "../envs/bcftools.yaml"
    input:
        reference_fasta = OUTDIR/f"{OUTPUT_NAME}.ancestor.fasta",
        bams = expand("{bam}", bam = iter_files_in_path(BAM_FOLDER))
    output:
        vcf = OUTDIR/f"{OUTPUT_NAME}.vcf"
    shell:
        """
        set +o pipefail

        sed 's/>.*/>MN908947.3/g' {input.reference_fasta} > renamed_reference.fasta
        for path in {input.bams}; do 
            ID=`echo $path | grep -o COV[0-9]{6}`
            bcftools mpileup -Ou -d 2000 -f renamed_reference.fasta $path  | bcftools call --ploidy 1 -mv -Ov | grep -v "#" | sed 's/MN908947.3/'$ID'/g' > $ID.pre.vcf
        done
        cat *pre.vcf > {output.vcf}
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
