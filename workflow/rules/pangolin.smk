rule pangolin_report:
    threads: 1
    conda: "../envs/pangolin.yaml"
    input:
        fastas = OUTDIR/f"{OUTPUT_NAME}.fasta"
    output:
        report = OUTDIR/f"{OUTPUT_NAME}.lineage_report.csv"
    shell:
        """
        pangolin {input.fastas} --outfile {output.report}
        """
