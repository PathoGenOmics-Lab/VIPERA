rule compile_fail_sites_vcf:
    params:
        header = ("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"),
        filter_text = "mask",
        sub_text = "NA",
        exc_text = "site_qual"
    input:
        sites = OUTDIR / f"{OUTPUT_NAME}.fail_sites.tsv",
    output:
        sites = temp(OUTDIR / f"{OUTPUT_NAME}.fail_sites.vcf"),
    run:
        import pandas as pd
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
        sites[list(params.header)].to_csv(output.sites, sep="\t", index=False)


rule merge_sites:
    params:
        header = ("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")
    input:
        lambda wildcards: select_problematic_vcf(),
        OUTDIR / f"{OUTPUT_NAME}.fail_sites.vcf",
    output:
        sites = temp(OUTDIR / "all_mask_sites.vcf"),
    run:
        import pandas as pd
        (
            pd.concat(
                [pd.read_table(path, sep="\t", comment="#", names=params.header) for path in input],
                axis="rows",
                ignore_index=True
            )
                .drop_duplicates(subset=("#CHROM", "POS", "FILTER"), keep="first")
                .sort_values(list(params.header))
                .to_csv(output.sites, sep="\t", index=False)
        )


rule extract_afwdist_variants:
    conda: "../envs/biopython.yaml"
    params:
        sample_col = "SAMPLE",
        position_col = "POS",
        sequence_col = "ALT",
        frequency_col = "ALT_FREQ",
        mask_class = ["mask"],
    input:
        variants = OUTDIR/f"{OUTPUT_NAME}.variants.tsv",
        mask_vcf = OUTDIR / "all_mask_sites.vcf",
        ancestor = OUTDIR/f"{OUTPUT_NAME}.ancestor.fasta",
        reference = OUTDIR/"reference.fasta",
    output:
        variants = temp(OUTDIR/f"{OUTPUT_NAME}.variants.afwdist.csv"),
    log:
        LOGDIR/"extract_afwdist_variants"/"log.txt"
    script:
        "../scripts/extract_afwdist_variants.py"


rule afwdist_weighted_distances:
    conda: "../envs/afwdist.yaml"
    params:
        extra_args = "",
    input:
        variants = OUTDIR/f"{OUTPUT_NAME}.variants.afwdist.csv",
        reference = OUTDIR/f"{OUTPUT_NAME}.ancestor.fasta",
    output:
        distances = temp(OUTDIR/f"{OUTPUT_NAME}.distances.raw.csv"),
    log:
        LOGDIR/"afwdist_weighted_distances"/"log.txt"
    shell:
        "afwdist "
        "-i {input.variants:q} "
        "-r {input.reference:q} "
        "-o {output.distances:q} "
        "{params.extra_args} >{log:q} 2>&1"


rule format_afwdist_results:
    conda: "../envs/biopython.yaml"
    params:
        samples = sorted(iter_samples()) + [config["ALIGNMENT_REFERENCE"]],
    input:
        distances = OUTDIR/f"{OUTPUT_NAME}.distances.raw.csv",
    output:
        distances = OUTDIR/f"{OUTPUT_NAME}.distances.csv",
    log:
        LOGDIR/"format_afwdist_results"/"log.txt"
    script:
        "../scripts/format_afwdist_results.py"
