rule bcftools_mpileup_all_sites:
    threads: 1
    conda: "../envs/var_calling.yaml"
    params:
        min_mq = 0,
        min_bq = config["VC"]["MIN_QUALITY"],
        mpileup_extra = "--no-BAQ"
    input:
        bam = get_input_bam,
        reference = OUTDIR/"vaf"/"{sample}.reference.fasta",
    output:
        mpileup = temp(OUTDIR / "all_sites" / "{sample}.mpileup.vcf"),
        query = temp(OUTDIR / "all_sites" / "{sample}.query.tsv"),
    log:
        mpileup = LOGDIR / "bcftools_mpileup_all_sites" / "{sample}.mpileup.txt",
        query = LOGDIR / "bcftools_mpileup_all_sites" / "{sample}.query.txt",
    shell:
        "bcftools mpileup {params.mpileup_extra} -a AD,ADF,ADR --fasta-ref {input.reference:q} --threads {threads} -Q {params.min_bq} -q {params.min_mq} -Ov -o {output.mpileup:q} {input.bam:q} >{log.mpileup:q} 2>&1 && "
        "echo 'CHROM\tPOS\tREF\tALT\tDP\tAD\tADF\tADR' >{output.query:q} && "
        "bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%DP\t[ %AD]\t[ %ADF]\t[ %ADR]\n' {output.mpileup:q} >>{output.query:q} 2>{log.query:q}"


rule filter_mpileup_all_sites:
    threads: 1
    params:
        min_total_AD = config["VC"]["MIN_DEPTH"],
        min_total_ADF = 0,
        min_total_ADR = 0,
    input:
        OUTDIR / "all_sites" / "{sample}.query.tsv",
    output:
        sites_pass = temp(OUTDIR / "all_sites" / "{sample}.filtered_sites.tsv"),
        sites_fail = temp(OUTDIR / "all_sites" / "{sample}.fail_sites.tsv"),
    run:
        import pandas as pd
        df = pd.read_csv(input[0], sep="\t")
        df["SAMPLE"] = wildcards.sample
        df["REF_AD"] = df.AD.str.split(",").apply(lambda values: int(values[0]))
        df["TOTAL_AD"] = df.AD.str.split(",").apply(lambda values: sum(int(n) for n in values))
        df["TOTAL_ADF"] = df.ADF.str.split(",").apply(lambda values: sum(int(n) for n in values))
        df["TOTAL_ADR"] = df.ADR.str.split(",").apply(lambda values: sum(int(n) for n in values))
        mask = (
            (df.TOTAL_AD >= params.min_total_AD) &
            (df.TOTAL_ADF >= params.min_total_ADF) &
            (df.TOTAL_ADR >= params.min_total_ADR)
        )
        df[mask].to_csv(output.sites_pass, sep="\t", index=False)
        df[~mask].to_csv(output.sites_fail, sep="\t", index=False)


use rule concat_vcf_fields as merge_filtered_mpileup_all_sites with:
    input:
        expand(OUTDIR / "all_sites" / "{sample}.filtered_sites.tsv", sample=iter_samples()),
    output:
        OUTDIR / f"{OUTPUT_NAME}.filtered_sites.tsv",


use rule concat_vcf_fields as merge_fail_mpileup_all_sites with:
    input:
        expand(OUTDIR / "all_sites" / "{sample}.fail_sites.tsv", sample=iter_samples()),
    output:
        OUTDIR / f"{OUTPUT_NAME}.fail_sites.tsv",


rule fill_all_sites:
    conda: "../envs/renv.yaml"
    input:
        variants = OUTDIR/f"{OUTPUT_NAME}.variants.tsv",
        sites = OUTDIR / f"{OUTPUT_NAME}.filtered_sites.tsv",
    output:
        variants = OUTDIR/f"{OUTPUT_NAME}.variants.all_sites.tsv",
    log:
        LOGDIR / "fill_all_sites" / "log.txt"
    script:
        "../scripts/fill_all_sites.R"


rule compile_fail_sites_vcf:
    params:
        filter_text = "mask",
        sub_text = "NA",
        exc_text = "site_qual",
    input:
        sites = OUTDIR / f"{OUTPUT_NAME}.fail_sites.tsv",
    output:
        sites = temp(OUTDIR / f"{OUTPUT_NAME}.fail_sites.vcf"),
    run:
        import pandas as pd
        HEADER = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
        sites = (
            pd.read_table(input.sites, sep="\t")
                .drop_duplicates(subset=("CHROM", "POS", "REF"))
                .rename(columns={"CHROM": "#CHROM"})
        )
        sites["ID"] = "."
        sites["ALT"] = "."
        sites["QUAL"] = "."
        sites["FILTER"] = params.filter_text
        sites["INFO"] = f"SUB={params.sub_text};EXC={params.exc_text}"
        sites[HEADER].to_csv(output.sites, sep="\t", index=False)


rule merge_mask_sites_vcf:
    input:
        lambda wildcards: select_problematic_vcf(),
        OUTDIR / f"{OUTPUT_NAME}.fail_sites.vcf",
    output:
        sites = temp(OUTDIR / "all_mask_sites.vcf"),
    run:
        import pandas as pd
        HEADER = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
        (
            pd.concat(
                [pd.read_table(path, sep="\t", comment="#", names=HEADER, dtype={"POS": "int64"}) for path in input],
                axis="rows",
                ignore_index=True
            )
            .drop_duplicates(subset=("#CHROM", "POS", "FILTER"), keep="first")
            .sort_values(by=["#CHROM", "POS"])
            .to_csv(output.sites, sep="\t", index=False)
        )
