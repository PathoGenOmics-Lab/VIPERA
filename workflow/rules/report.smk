rule window:
    conda: "../envs/biopython.yaml"
    input:
        vcf = OUTDIR/f"{OUTPUT_NAME}.masked.filtered.tsv",
    params:
        step = 1000
    output:
        window_df = OUTDIR/f"{OUTPUT_NAME}.window.csv",
    script:
        "../scripts/window.py"


rule diversity:
    threads: 4
    conda: "../envs/renv.yaml"
    params:
        design = config["PLOTS"],
        bootstrap_reps = config["DIVERSITY_BOOTSTRAP_REPS"],
        plot_width = 159.2,
        plot_height = 119.4
    input:
        study_fasta = OUTDIR/"nextalign"/f"{OUTPUT_NAME}.aligned.masked.fasta",
        context_fasta = OUTDIR/"context"/"nextalign"/"context_sequences.aligned.masked.fasta"
    output:
        fig = report((REPORT_DIR/"div.plot.png").resolve()),
        value = temp((OUTDIR/"our_diversity.txt").resolve())
    script:
        "../scripts/report/diversity_plot.R"


rule freyja_plot:
    conda: "../envs/renv.yaml"
    params:
        design = config["PLOTS"],
        metadata = config["METADATA"]
    input:
        summary_demixing =  OUTDIR/"summary_freyja_demixing.csv"
    output:
        fig = report((REPORT_DIR/"freyja.plot.png").resolve())
    script:
        "../scripts/report/freyja_plot.R"


rule general_NV_description:
    conda: "../envs/renv.yaml"
    params:
        design = config["PLOTS"],
        nsp = config["NSP"],
        metadata = config["METADATA"]
    input:
        window = OUTDIR/f"{OUTPUT_NAME}.window.csv",
        vcf =  OUTDIR/f"{OUTPUT_NAME}.masked.filtered.tsv"
    output:
        fig = report((REPORT_DIR/"NV.description.png").resolve()),
        fig_cor = report((REPORT_DIR/"cor_snp_time.png").resolve()),
        summary_nv = temp((OUTDIR/"summary_nv.csv").resolve())
    script:
        "../scripts/report/NV_description.R"


rule phylo_plots:
    conda: "../envs/renv.yaml"
    params: 
        design = config["PLOTS"],
        metadata = config["METADATA"],
        ref_name = config["ALIGNMENT_REFERENCE"],
        boot_th = 85,
        boot_color = "green"
    input: 
        dist = OUTDIR/f"{OUTPUT_NAME}.weighted_distances.csv",
        study_fasta = OUTDIR/f"{OUTPUT_NAME}.fasta",
        ml = OUTDIR/f"tree_context/{OUTPUT_NAME}.treefile"
    output:
        temest = report((REPORT_DIR/"temp_est.png").resolve()),
        tree = report((REPORT_DIR/"tree.png").resolve()),
        tree_ml = report((REPORT_DIR/"tree_ml.png").resolve()),
        stats_lm = temp((OUTDIR/"stats.lm.csv").resolve())
    script:
        "../scripts/report/pylo.R"


rule evo_plots:
    conda: "../envs/renv.yaml"
    params: 
        design = config["PLOTS"],
        metadata = config["METADATA"]
    input: 
        N_S = OUTDIR/f"{OUTPUT_NAME}.ancestor.N_S.sites.csv",
        vcf =  OUTDIR/f"{OUTPUT_NAME}.masked.filtered.tsv"
    output:
        plot = report((REPORT_DIR/"dn_ds.png").resolve())
    script:
        "../scripts/report/evo_plots.R"


rule snp_plots:
    conda: "../envs/renv.yaml"
    params:
        design = config["PLOTS"],
        metadata = config["METADATA"],
    input:
         vcf =  OUTDIR/f"{OUTPUT_NAME}.masked.filtered.tsv"
    output:
        pseudovolcano = report((REPORT_DIR/"volcano.png").resolve()),
        snp_panel = report((REPORT_DIR/"panel.png").resolve())
    script:
        "../scripts/report/snp_plots.R"

rule summary_table:
    conda: "../envs/renv.yaml"
    params:
        metadata = config["METADATA"],
    input:
        report = report(OUTDIR/f"{OUTPUT_NAME}.lineage_report.csv")
    output:
        table = temp((OUTDIR/"summary_table.csv").resolve())
    script:
        "../scripts/report/summary_table.R"

rule report:
    conda: "../envs/renv.yaml"
    shadow: "shallow"
    input:
        qmd       = Path(config["REPORT_QMD"]).resolve(),
        diversity = report(rules.diversity.output.fig),
        freyja    = report(rules.freyja_plot.output.fig),
        tree      = report(rules.phylo_plots.output.tree),
        temest    = report(rules.phylo_plots.output.temest),
        SNV       = report(rules.general_NV_description.output.fig),
        evo       = report(rules.evo_plots.output.plot),
        value     = rules.diversity.output.value,
        panel     = report(rules.snp_plots.output.pseudovolcano),
        volcano   = report(rules.snp_plots.output.snp_panel),
        tree_ml   = report(rules.phylo_plots.output.tree_ml),
        fig_cor   = report(rules.general_NV_description.output.fig_cor),
        stats_lm  = rules.phylo_plots.output.stats_lm,
        table     = rules.summary_table.output.table,
        sum_nv    = rules.general_NV_description.output.summary_nv
    output:
        html = report(OUTDIR/f"{OUTPUT_NAME}.report.html")
    shell:
        """
        set +o pipefail
        Rscript -e 'library(quarto)' -e \"quarto_render(input = '{input.qmd}',\
                                           execute_params=list(div='{input.diversity}',\
                                                       freyja ='{input.freyja}',\
                                                       tree = '{input.tree}',\
                                                       tempest = '{input.temest}',\
                                                       SNV = '{input.SNV}',\
                                                       evo = '{input.evo}',
                                                       div_value = '{input.value}',
                                                       panel = '{input.panel}',
                                                       volcano = '{input.volcano}',
                                                       tree_ml = '{input.tree_ml}',
                                                       fig_cor_snp = '{input.fig_cor}',
                                                       stats_lm = '{input.stats_lm}',
                                                       table = '{input.table}',
                                                       sum_nv = '{input.sum_nv}'))\"    
        mv "$(dirname {input.qmd:q})/report.html" {output.html}
        """
