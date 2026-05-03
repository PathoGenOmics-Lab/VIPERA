rule pangolin_report:
    threads: 8
    conda: "../envs/pangolin.yaml"
    input:
        fastas = "<results>/<dataset>/sequences.fasta"
    output:
        report = report("<results>/<dataset>/lineage_report.csv")
    log:
        LOGDIR / "pangolin_report" / "log.txt"
    shell:
        """
        pangolin {input.fastas} --outfile {output.report} -t {threads} >{log} 2>&1
        """
