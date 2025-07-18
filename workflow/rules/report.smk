rule heatmap:
    conda: "../envs/renv.yaml"
    input:
        vcf =  OUTDIR/f"{OUTPUT_NAME}.masked.filtered.tsv",
        metadata = config["METADATA"]
    output:
        table = report(REPORT_DIR_TABLES/"figure_10.csv")
    log:
        LOGDIR / "heatmap" / "log.txt"
    script:
        "../scripts/report/heatmap.R"


rule window:
    conda: "../envs/biopython.yaml"
    params:
        window = config["WINDOW"]["WIDTH"],
        step = config["WINDOW"]["STEP"]
    input:
        vcf = OUTDIR/f"{OUTPUT_NAME}.masked.filtered.tsv",
        gb = OUTDIR/"reference.gb",
        features = config["FEATURES_JSON"]
    output:
        window_df = temp(OUTDIR/f"{OUTPUT_NAME}.window.csv"),
    log:
        LOGDIR / "window" / "log.txt"
    script:
        "../scripts/window.py"


rule diversity:
    threads: 4
    conda: "../envs/renv.yaml"
    params:
        design = config["PLOTS"],
        bootstrap_reps = config["DIVERSITY_REPS"],
        plot_width = 159.2,
        plot_height = 119.4
    input:
        study_fasta = OUTDIR/"nextalign"/f"{OUTPUT_NAME}.aligned.masked.fasta",
        context_fasta = OUTDIR/"context"/"nextalign"/"context_sequences.aligned.masked.fasta"
    output:
        fig = report(REPORT_DIR_PLOTS/"figure_3.png"),
        json = temp(OUTDIR/"diversity.json"),
        table = REPORT_DIR_TABLES/"figure_3.csv"
    log:
        LOGDIR / "diversity" / "log.txt"
    script:
        "../scripts/report/diversity_plot.R"


rule freyja_plot:
    conda: "../envs/renv.yaml"
    params:
        design = config["PLOTS"]
    input:
        summary_demixing =  OUTDIR/"summary_freyja_demixing.csv",
        metadata = config["METADATA"]
    output:
        fig = report(REPORT_DIR_PLOTS/"figure_1.png"),
        table = report(REPORT_DIR_TABLES/"figure_1.csv")
    log:
        LOGDIR / "freyja_plot" / "log.txt"
    script:
        "../scripts/report/freyja_plot.R"


rule general_NV_description:
    conda: "../envs/renv.yaml"
    params:
        samples = expand("{sample}", sample = iter_samples()),
        design = config["PLOTS"],
        regions = config["PLOT_GENOME_REGIONS"],
        window = config["WINDOW"]["WIDTH"],
        step = config["WINDOW"]["STEP"],
        max_alt_freq = 1.0 - config["VC"]["IVAR_FREQ"]
    input:
        window = OUTDIR/f"{OUTPUT_NAME}.window.csv",
        vcf =  OUTDIR/f"{OUTPUT_NAME}.masked.filtered.tsv",
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


rule phylo_plots:
    conda: "../envs/renv.yaml"
    params:
        design = config["PLOTS"],
        ref_name = config["ALIGNMENT_REFERENCE"],
        boot_th = 95,
        alrt_th = 80,
        plot_height_mm = 119.4,
        plot_width_mm = 159.2,
        use_bionj = config["USE_BIONJ"]
    input:
        dist = REPORT_DIR_TABLES/f"figure_8.csv",
        study_fasta = OUTDIR/f"{OUTPUT_NAME}.fasta",
        ml = OUTDIR/f"tree_context/{OUTPUT_NAME}.treefile",
        metadata = config["METADATA"]
    output:
        temest = report(REPORT_DIR_PLOTS/"figure_9.png"),
        tree = report(REPORT_DIR_PLOTS/"figure_8.png"),
        tree_ml = report(REPORT_DIR_PLOTS/"figure_2.png"),
        table = report(REPORT_DIR_TABLES/"figure_9.csv"),
        json = temp(OUTDIR/"stats.lm.json")
    log:
        LOGDIR / "phylo_plots" / "log.txt"
    script:
        "../scripts/report/phylo_plots.R"


rule evo_plots:
    conda: "../envs/renv.yaml"
    params: 
        design = config["PLOTS"]
    input: 
        N_S = OUTDIR/f"{OUTPUT_NAME}.ancestor.N_S.sites.csv",
        vcf =  OUTDIR/f"{OUTPUT_NAME}.masked.filtered.tsv",
        metadata = config["METADATA"]
    output:
        plot = report(REPORT_DIR_PLOTS/"figure_11.png"),
        plot_omega = report(REPORT_DIR_PLOTS/"figure_12.png"),
        table = report(REPORT_DIR_TABLES/"figure_11.csv")
    log:
        LOGDIR / "evo_plots" / "log.txt"
    script:
        "../scripts/report/evo_plots.R"


rule snp_plots:
    conda: "../envs/renv.yaml"
    params:
        design = config["PLOTS"],
        cor_method = config["COR"]["METHOD"],
        cor_exact = config["COR"]["EXACT"]
    input:
        vcf =  OUTDIR/f"{OUTPUT_NAME}.masked.filtered.tsv",
        metadata = config["METADATA"]
    output:
        pseudovolcano = report(REPORT_DIR_PLOTS/"figure_6.png"),
        snp_panel = report(REPORT_DIR_PLOTS/"figure_7.png"),
        table_1 = report(REPORT_DIR_TABLES/"figure_6.csv"),
        table_2 = report(REPORT_DIR_TABLES/"figure_7.csv")
    log:
        LOGDIR / "snp_plots" / "log.txt"
    script:
        "../scripts/report/snp_plots.R"


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
        freyja     = report(REPORT_DIR_PLOTS/"figure_1.png"),
        tree_ml    = report(REPORT_DIR_PLOTS/"figure_2.png"),
        diversity  = report(REPORT_DIR_PLOTS/"figure_3.png"),
        fig_cor    = report(REPORT_DIR_PLOTS/"figure_4.png"),
        SNV        = report(REPORT_DIR_PLOTS/"figure_5a.png"),
        SNV_spike  = report(REPORT_DIR_PLOTS/"figure_5b.png"),
        volcano    = report(REPORT_DIR_PLOTS/"figure_6.png"),
        panel      = report(REPORT_DIR_PLOTS/"figure_7.png"),
        tree       = report(REPORT_DIR_PLOTS/"figure_8.png"),
        temest     = report(REPORT_DIR_PLOTS/"figure_9.png"),
        heat_table = report(REPORT_DIR_TABLES/"figure_10.csv"),
        evo        = report(REPORT_DIR_PLOTS/"figure_11.png"),
        omega_plot = report(REPORT_DIR_PLOTS/"figure_12.png"),
        value      = OUTDIR/"diversity.json",
        stats_lm   = OUTDIR/"stats.lm.json",
        table      = OUTDIR/"summary_table.csv",
        sum_nv     = OUTDIR/"summary_nv.json"
    params:
        workflow_version = get_repo_version(BASE_PATH.as_posix(), __version__),
        min_ivar_freq = config["VC"]["IVAR_FREQ"],
        ufboot_reps = config["UFBOOT_REPS"],
        shalrt_reps = config["SHALRT_REPS"],
        name = config["OUTPUT_NAME"],
        use_bionj = config["USE_BIONJ"],
        cor_method = config["COR"]["METHOD"]
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
                "freyja ='{input.freyja}', "
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
                "table = '{input.table}', "
                "sum_nv = '{input.sum_nv}', "
                "heat_tab = '{input.heat_table}', "
                "omega_plot = '{input.omega_plot}', "
                "name = '{params.name}'))\" "
        ">{log} 2>&1 && "
        'mv "$(dirname {input.qmd:q})/report.html" {output.html:q}'
