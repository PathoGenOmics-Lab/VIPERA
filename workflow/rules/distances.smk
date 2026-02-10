rule compile_fail_sites_vcf:
    params:
        filter_text = "fail_site",
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


rule extract_afwdist_variants:
    conda: "../envs/biopython.yaml"
    params:
        sample_col = "SAMPLE",
        position_col = "POS",
        sequence_col = "ALT",
        frequency_col = "ALT_FREQ",
        mask_class = ["mask", "fail_site"],
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
