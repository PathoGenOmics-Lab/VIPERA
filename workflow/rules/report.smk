rule demix_plot_data:
    conda: "../envs/renv.yaml"
    input:
        summary_demixing = OUTDIR/"demixing"/"summary.csv",
        metadata = config["METADATA"]
    output:
        data = report(REPORT_DIR_TABLES/"demix.csv")
    log:
        LOGDIR / "demix_plot_data" / "log.txt"
    script:
        "../scripts/report/demix_plot_data.R"


rule demix_plot:
    conda: "../envs/renv.yaml"
    params:
        design = config["PLOTS"],
        plot_width_mm = 159.2,
        plot_height_mm = 119.4,
    input:
        data = REPORT_DIR_TABLES/"demix.csv"
    output:
        plot = report(REPORT_DIR_PLOTS/"demix.png")
    log:
        LOGDIR / "demix_plot" / "log.txt"
    script:
        "../scripts/report/demix_plot.R"


rule heatmap_plot_data:
    conda: "../envs/renv.yaml"
    input:
        vcf =  OUTDIR/f"{OUTPUT_NAME}.variants.tsv",
        metadata = config["METADATA"]
    output:
        table = report(REPORT_DIR_TABLES/"heatmap.csv")
    log:
        LOGDIR / "heatmap" / "log.txt"
    script:
        "../scripts/report/heatmap.R"


rule window:
    conda: "../envs/biopython.yaml"
    params:
        window = config["WINDOW"]["WIDTH"],
        step = config["WINDOW"]["STEP"],
        features = config.get("GB_FEATURES", {}),
        gb_qualifier_display = "gene"
    input:
        variants = OUTDIR/f"{OUTPUT_NAME}.variants.tsv",
        gb = OUTDIR/"reference.gb",
    output:
        window_df = temp(OUTDIR/f"{OUTPUT_NAME}.window.csv"),
    log:
        LOGDIR / "window" / "log.txt"
    script:
        "../scripts/window.py"


rule diversity_data:
    threads: 4
    conda: "../envs/renv.yaml"
    params:
        bootstrap_reps = config["DIVERSITY_REPS"],
        aln_reference = config["ALIGNMENT_REFERENCE"],
    input:
        study_fasta = OUTDIR/"nextalign"/f"{OUTPUT_NAME}.aligned.masked.fasta",
        context_fasta = OUTDIR/"context"/"nextalign"/"context_sequences.aligned.masked.fasta",
    output:
        divs = REPORT_DIR_TABLES/"diversity.txt",
        json = REPORT_DIR_TABLES/"diversity.json",
    log:
        LOGDIR / "diversity_data" / "log.txt"
    script:
        "../scripts/report/diversity_data.R"


rule diversity_plot:
    threads: 1
    conda: "../envs/renv.yaml"
    params:
        design = config["PLOTS"],
        plot_width_mm = 159.2,
        plot_height_mm = 119.4,
    input:
        divs = REPORT_DIR_TABLES/"diversity.txt",
        json = REPORT_DIR_TABLES/"diversity.json",
    output:
        plot = report(REPORT_DIR_PLOTS/"diversity.png"),
    log:
        LOGDIR / "diversity_plot" / "log.txt"
    script:
        "../scripts/report/diversity_plot.R"


rule general_NV_description:
    conda: "../envs/renv.yaml"
    params:
        samples = expand("{sample}", sample = iter_samples()),
        design = config["PLOTS"],
        window = config["WINDOW"]["WIDTH"],
        step = config["WINDOW"]["STEP"],
        max_alt_freq = 1.0 - config["VC"]["IVAR_FREQ"]
    input:
        coordinates = config["COORDINATES_JSON"],
        regions = config["PLOT_GENOME_REGIONS"],
        window = OUTDIR/f"{OUTPUT_NAME}.window.csv",
        vcf =  OUTDIR/f"{OUTPUT_NAME}.variants.tsv",
        metadata = config["METADATA"]
    output:
        fig = report(REPORT_DIR_PLOTS/"figure_5a.png"),
        fig_s = report(REPORT_DIR_PLOTS/"figure_5b.png"),
        fig_cor = report(REPORT_DIR_PLOTS/"figure_4.png"),
        json = temp(OUTDIR/"summary_nv.json"),
        table_1 = report(REPORT_DIR_TABLES/"figure_5a.csv"),
        table_2 = report(REPORT_DIR_TABLES/"figure_5b.csv"),
        table_3 = report(REPORT_DIR_TABLES/"figure_4.csv")
    log:
        LOGDIR / "general_NV_description" / "log.txt"
    script:
        "../scripts/report/NV_description.R"


rule context_phylogeny_data:
    conda: "../envs/renv.yaml"
    params:
        design = config["PLOTS"],
        ref_name = config["ALIGNMENT_REFERENCE"],
        boot_th = 95.0,
        alrt_th = 80.0,
    input:
        target_fasta = OUTDIR/f"{OUTPUT_NAME}.fasta",
        tree = OUTDIR/f"tree_context/{OUTPUT_NAME}.treefile",
    output:
        json = REPORT_DIR_TABLES/"context_phylogeny.json",
        annotation = REPORT_DIR_TABLES/"context_phylogeny.csv",
    log:
        LOGDIR / "context_phylogeny_data" / "log.txt"
    script:
        "../scripts/report/context_phylogeny_data.R"


rule context_phylogeny_plot:
    conda: "../envs/renv.yaml"
    params:
        design = config["PLOTS"],
        plot_height_mm = 119.4,
        plot_width_mm = 159.2,
    input:
        tree = OUTDIR/f"tree_context/{OUTPUT_NAME}.treefile",
        json = REPORT_DIR_TABLES/"context_phylogeny.json",
        annotation = REPORT_DIR_TABLES/"context_phylogeny.csv"
    output:
        plot = report(REPORT_DIR_PLOTS/"context_phylogeny.png"),
    log:
        LOGDIR / "context_phylogeny_plot" / "log.txt"
    script:
        "../scripts/report/context_phylogeny_plot.R"


rule allele_freq_tree_data:
    conda: "../envs/renv.yaml"
    params:
        use_bionj = config["USE_BIONJ"],
        ref_name = config["ALIGNMENT_REFERENCE"],
    input:
        dist = REPORT_DIR_TABLES/"distances.csv",
    output:
        tree = report(REPORT_DIR_TABLES/"allele_freq_tree.nwk"),
    log:
        LOGDIR / "allele_freq_tree_data" / "log.txt"
    script:
        "../scripts/report/allele_freq_tree_data.R"


rule allele_freq_tree_plot:
    conda: "../envs/renv.yaml"
    params:
        design = config["PLOTS"],
        ref_name = config["ALIGNMENT_REFERENCE"],
        plot_height_mm = 119.4,
        plot_width_mm = 159.2,
    input:
        tree = report(REPORT_DIR_TABLES/"allele_freq_tree.nwk"),
        study_fasta = OUTDIR/f"{OUTPUT_NAME}.fasta",
        metadata = config["METADATA"],
    output:
        plot = report(REPORT_DIR_PLOTS/"allele_freq_tree.png"),
    log:
        LOGDIR / "allele_freq_tree_plot" / "log.txt"
    script:
        "../scripts/report/allele_freq_tree_plot.R"


rule time_signal_data:
    conda: "../envs/renv.yaml"
    params:
        ref_name = config["ALIGNMENT_REFERENCE"],
    input:
        tree = report(REPORT_DIR_TABLES/"allele_freq_tree.nwk"),
        metadata = config["METADATA"],
    output:
        table = report(REPORT_DIR_TABLES/"time_signal.csv"),
        json = REPORT_DIR_TABLES/"time_signal.json"
    log:
        LOGDIR / "time_signal_data" / "log.txt"
    script:
        "../scripts/report/time_signal_data.R"


rule time_signal_plot:
    conda: "../envs/renv.yaml"
    params:
        design = config["PLOTS"],
        plot_height_mm = 119.4,
        plot_width_mm = 159.2,
    input:
        table = report(REPORT_DIR_TABLES/"time_signal.csv"),
    output:
        plot = report(REPORT_DIR_PLOTS/"time_signal.png"),
    log:
        LOGDIR / "time_signal_plot" / "log.txt"
    script:
        "../scripts/report/time_signal_plot.R"


rule dnds_plots:
    conda: "../envs/renv.yaml"
    params: 
        design = config["PLOTS"]
    input: 
        table = REPORT_DIR_TABLES/"dnds.csv",
    output:
        plot_dn_ds = report(REPORT_DIR_PLOTS/"dn_and_ds.png"),
        plot_omega = report(REPORT_DIR_PLOTS/"dnds.png"),
    log:
        LOGDIR / "evo_plots" / "log.txt"
    script:
        "../scripts/report/dnds_plots.R"


rule af_time_correlation_per_snp_data:
    conda: "../envs/renv.yaml"
    params:
        cor_method = config["COR"]["METHOD"],
        cor_exact = config["COR"]["EXACT"],
    input:
        variants = OUTDIR/f"{OUTPUT_NAME}.variants.tsv",
        metadata = config["METADATA"],
    output:
        fmt_variants = temp(OUTDIR/f"{OUTPUT_NAME}.variants.filled.dated.tsv"),
        correlations = report(REPORT_DIR_TABLES/"af_time_correlation_per_snp.csv"),  # fig 6
    log:
        LOGDIR / "af_time_correlation_per_snp_data" / "log.txt"
    script:
        "../scripts/report/af_time_correlation_per_snp_data.R"


rule af_time_correlation_per_snp_plot:
    conda: "../envs/renv.yaml"
    params:
        design = config["PLOTS"],
    input:
        correlations = REPORT_DIR_TABLES/"af_time_correlation_per_snp.csv",
    output:
        plot = report(REPORT_DIR_PLOTS/"af_time_correlation_per_snp.png"),  # fig 6
    log:
        LOGDIR / "af_time_correlation_per_snp_plot" / "log.txt"
    script:
        "../scripts/report/af_time_correlation_per_snp_plot.R"


rule af_trajectory_panel_data:  # panel with AF trajectories in time
    # TODO
    conda: "../envs/renv.yaml"
    input:
        fmt_variants = OUTDIR/f"{OUTPUT_NAME}.variants.filled.dated.tsv",
        correlations = REPORT_DIR_TABLES/"af_time_correlation_per_snp.csv",
    output:
        table = report(REPORT_DIR_TABLES/"af_trajectory_panel.csv"),  # fig 7
    log:
        LOGDIR / "af_trajectory_panel_data" / "log.txt"
    script:
        "../scripts/report/af_trajectory_panel_data.R"


rule af_trajectory_panel_plot:  # panel with AF trajectories in time
    # TODO
    conda: "../envs/renv.yaml"
    params:
        design = config["PLOTS"],
    input:
        table = REPORT_DIR_TABLES/"af_trajectory_panel.csv",
    output:
        plot = report(REPORT_DIR_PLOTS/"af_trajectory_panel.png"),  # fig 7
    log:
        LOGDIR / "af_trajectory_panel_plot" / "log.txt"
    script:
        "../scripts/report/af_trajectory_panel_plot.R"


# rule snp_plots:
#     conda: "../envs/renv.yaml"
#     params:
#         design = config["PLOTS"],
#         cor_method = config["COR"]["METHOD"],
#         cor_exact = config["COR"]["EXACT"]
#     input:
#         vcf =  OUTDIR/f"{OUTPUT_NAME}.variants.tsv",
#         metadata = config["METADATA"]
#     output:
#         pseudovolcano = report(REPORT_DIR_PLOTS/"figure_6.png"),
#         snp_panel = report(REPORT_DIR_PLOTS/"figure_7.png"),
#         table_1 = report(REPORT_DIR_TABLES/"figure_6.csv"),
#         table_2 = report(REPORT_DIR_TABLES/"figure_7.csv")
#     log:
#         LOGDIR / "snp_plots" / "log.txt"
#     script:
#         "../scripts/report/snp_plots.R"


rule summary_table:
    conda: "../envs/renv.yaml"
    input:
        report = report(OUTDIR/f"{OUTPUT_NAME}.lineage_report.csv"),
        metadata = config["METADATA"]
    output:
        table = temp(OUTDIR/"summary_table.csv")
    log:
        LOGDIR / "summary_table" / "log.txt"
    script:
        "../scripts/report/summary_table.R"


rule report:
    conda: "../envs/quarto_render.yaml"
    shadow: "shallow"
    input:
        qmd        = Path(config["REPORT_QMD"]).resolve(),
        demix     = report(REPORT_DIR_PLOTS/"demix.png"),
        tree_ml    = report(REPORT_DIR_PLOTS/"context_phylogeny.png"),
        diversity  = report(REPORT_DIR_PLOTS/"diversity.png"),
        fig_cor    = report(REPORT_DIR_PLOTS/"figure_4.png"),
        SNV        = report(REPORT_DIR_PLOTS/"figure_5a.png"),
        SNV_spike  = report(REPORT_DIR_PLOTS/"figure_5b.png"),
        volcano    = report(REPORT_DIR_TABLES/"af_time_correlation_per_snp.png"),
        panel      = report(REPORT_DIR_TABLES/"af_trajectory_panel.png"),
        tree       = report(REPORT_DIR_PLOTS/"allele_freq_tree.png"),
        temest     = report(REPORT_DIR_PLOTS/"time_signal.png"),
        heat_table = report(REPORT_DIR_TABLES/"heatmap.csv"),
        evo        = report(REPORT_DIR_PLOTS/"dn_and_ds.png"),
        omega_plot = report(REPORT_DIR_PLOTS/"dnds.png"),
        freyja_ts  = OUTDIR/"demixing"/"freyja_data"/"last_barcode_update.txt",
        value      = REPORT_DIR_TABLES/"diversity.json",
        stats_lm   = REPORT_DIR_TABLES/"time_signal.json",
        stats_ml   = REPORT_DIR_TABLES/"context_phylogeny.json",
        table      = OUTDIR/"summary_table.csv",
        sum_nv     = OUTDIR/"summary_nv.json",
    params:
        workflow_version = get_repo_version(BASE_PATH.as_posix(), __version__),
        min_ivar_freq = config["VC"]["IVAR_FREQ"],
        ufboot_reps = config["UFBOOT_REPS"],
        shalrt_reps = config["SHALRT_REPS"],
        name = config["OUTPUT_NAME"],
        use_bionj = config["USE_BIONJ"],
        cor_method = config["COR"]["METHOD"],
    output:
        html = report(OUTDIR/f"{OUTPUT_NAME}.report.html")
    log:
        LOGDIR / "report" / "log.txt"
    shell:
        "set +o pipefail; "
        "Rscript -e 'library(quarto)' "
        "-e \"quarto_render("
            "input = '{input.qmd}', "
            "execute_params=list("
                "ufboot_reps='{params.ufboot_reps}', "
                "shalrt_reps='{params.shalrt_reps}', "
                "min_ivar_freq='{params.min_ivar_freq}', "
                "workflow_version='{params.workflow_version}', "
                "use_bionj='{params.use_bionj}', "
                "cor_method='{params.cor_method}', "
                "div='{input.diversity}', "
                "demix ='{input.demix}', "
                "tree = '{input.tree}', "
                "tempest = '{input.temest}', "
                "SNV = '{input.SNV}', "
                "SNV_s = '{input.SNV_spike}', "
                "evo = '{input.evo}', "
                "div_value = '{input.value}', "
                "panel = '{input.panel}', "
                "volcano = '{input.volcano}', "
                "tree_ml = '{input.tree_ml}', "
                "fig_cor_snp = '{input.fig_cor}', "
                "stats_lm = '{input.stats_lm}', "
                "stats_ml = '{input.stats_ml}', "
                "table = '{input.table}', "
                "sum_nv = '{input.sum_nv}', "
                "heat_tab = '{input.heat_table}', "
                "omega_plot = '{input.omega_plot}', "
                "freyja_ts = '{input.freyja_ts}', "
                "name = '{params.name}'))\" "
        ">{log:q} 2>&1 && "
        'mv "$(dirname {input.qmd:q})/report.html" {output.html:q}'
