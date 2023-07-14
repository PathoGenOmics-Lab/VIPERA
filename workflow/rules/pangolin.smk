rule pangolin_report:
    threads: 8
    conda: "../envs/pangolin.yaml"
    input:
        fastas = OUTDIR/f"{OUTPUT_NAME}.fasta"
    output:
        report = report(OUTDIR/f"{OUTPUT_NAME}.lineage_report.csv")
    shell:
        """
        pangolin {input.fastas} --outfile {output.report} -t {threads}
        """
