rule demix_plot_data:
    conda: "../envs/renv.yaml"
    input:
        summary_demixing = OUTDIR/"demixing"/"summary.csv",
        metadata = config["METADATA"]
    output:
        data = REPORT_DIR_TABLES/"demix.csv"
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


rule diversity_data:
    threads: 4
    conda: "../envs/renv.yaml"
    params:
        bootstrap_reps = config["DIVERSITY_REPS"],
        aln_reference = config["ALIGNMENT_REFERENCE"],
        seed = 7291,
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


rule extract_genbank_regions:
    conda: "../envs/biopython.yaml"
    params:
        gb_qualifier = "gene",
    input:
        gb = OUTDIR/"reference.cds.gb",
    output:
        regions = temp(REPORT_DIR_TABLES/"genbank_regions.json"),
    log:
        LOGDIR / "extract_genbank_regions" / "log.txt"
    script:
        "../scripts/report/extract_genbank_regions.py"


rule polymorphic_sites_over_time_data:
    conda: "../envs/renv.yaml"
    params:
        max_alt_freq = 1.0 - config["VC"]["MIN_FREQ"],
    input:
        variants = OUTDIR/f"{OUTPUT_NAME}.variants.tsv",
        metadata = config["METADATA"],
    output:
        table = REPORT_DIR_PLOTS/"polymorphic_sites_over_time.csv",
        json = temp(REPORT_DIR_TABLES/"polymorphic_sites_over_time.json"),
    log:
        LOGDIR / "polymorphic_sites_over_time_data" / "log.txt"
    script:
        "../scripts/report/polymorphic_sites_over_time_data.R"


rule polymorphic_sites_over_time_plot:
    conda: "../envs/renv.yaml"
    params:
        design = config["PLOTS"],
        plot_width_mm = 159.2,
        plot_height_mm = 119.4,
    input:
        table = REPORT_DIR_PLOTS/"polymorphic_sites_over_time.csv",
    output:
        plot = report(REPORT_DIR_PLOTS/"polymorphic_sites_over_time.png"),
    log:
        LOGDIR / "polymorphic_sites_over_time_plot" / "log.txt"
    script:
        "../scripts/report/polymorphic_sites_over_time_plot.R"


rule window_data:
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
        window_df = REPORT_DIR_TABLES/"window.csv",
        json = temp(REPORT_DIR_TABLES/"window.json"),
    log:
        LOGDIR / "window_data" / "log.txt"
    script:
        "../scripts/report/window_data.py"


rule nv_panel_data:
    conda: "../envs/renv.yaml"
    input:
        variants = OUTDIR/f"{OUTPUT_NAME}.variants.tsv",
        metadata = config["METADATA"],
    output:
        table = REPORT_DIR_TABLES/"nv_panel.csv",
        json = temp(REPORT_DIR_TABLES/"nv_panel.json"),
    log:
        LOGDIR / "nv_panel_data" / "log.txt"
    script:
        "../scripts/report/nv_panel_data.R"


rule nv_panel_zoom_on_feature_data:
    input:
        table = REPORT_DIR_TABLES/"nv_panel.csv",
        regions = REPORT_DIR_TABLES/"genbank_regions.json",
    output:
        table = temp(REPORT_DIR_TABLES/"nv_panel.{region_name}.csv"),
    log:
        LOGDIR / "nv_panel_zoom_on_feature_data" / "{region_name}.log.txt"
    script:
        "../scripts/report/nv_panel_zoom_on_feature_data.py"


rule window_zoom_on_feature_data:
    input:
        table = REPORT_DIR_TABLES/"window.csv",
        regions = REPORT_DIR_TABLES/"genbank_regions.json",
    output:
        table = temp(REPORT_DIR_TABLES/"window.{region_name}.csv"),
    log:
        LOGDIR / "window_zoom_on_feature_data" / "{region_name}.log.txt"
    script:
        "../scripts/report/window_zoom_on_feature_data.py"


rule nv_panel_plot:
    conda: "../envs/renv.yaml"
    params:
        design = config["PLOTS"],
        window_step = config["WINDOW"]["STEP"],
        plot_height_mm = 250.0,
        plot_width_mm = 240.0,
    input:
        panel = REPORT_DIR_TABLES/"nv_panel.csv",
        window = REPORT_DIR_TABLES/"window.csv",
        regions = REPORT_DIR_TABLES/"genbank_regions.json",
        highlight_window_regions = config["PLOT_GENOME_REGIONS"],
    output:
        plot = report(REPORT_DIR_PLOTS/"nv_panel.png"),
    log:
        LOGDIR / "nv_panel_plot" / "log.txt"
    script:
        "../scripts/report/nv_panel_plot.R"


use rule nv_panel_plot as nv_panel_plot_S with:
    input:
        panel = REPORT_DIR_TABLES/"nv_panel.S.csv",
        window = REPORT_DIR_TABLES/"window.S.csv",
        regions = REPORT_DIR_TABLES/"genbank_regions.json",
        highlight_window_regions = OUTDIR/"empty.txt",
    output:
        plot = report(REPORT_DIR_PLOTS/"nv_panel.S.png"),
    log:
        LOGDIR / "nv_panel_plot_S" / "log.txt"


rule merge_json_files:
    input:
        REPORT_DIR_TABLES/"nv_panel.json",
        REPORT_DIR_TABLES/"polymorphic_sites_over_time.json",
        REPORT_DIR_TABLES/"window.json",
    output:
        json = REPORT_DIR_TABLES/"nv_panel_summary.json",
    run:
        import json
        result = {}
        for path in input:
            with open(path) as f:
                d = json.load(f)
            result |= d  # will replace existing keys
        with open(output.json, "w") as fw:
            json.dump(result, fw, indent=2)


rule context_phylogeny_data:
    conda: "../envs/renv.yaml"
    params:
        design = config["PLOTS"],
        ref_name = config["ALIGNMENT_REFERENCE"],
        boot_th = config["UFBOOT"]["THRESHOLD"],
        alrt_th = config["SHALRT"]["THRESHOLD"],
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
        boot_th = config["UFBOOT"]["THRESHOLD"],
        alrt_th = config["SHALRT"]["THRESHOLD"],
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
        outgroup_id = config["ALIGNMENT_REFERENCE"],
    input:
        dist = OUTDIR/f"{OUTPUT_NAME}.distances.csv",
    output:
        tree = REPORT_DIR_TABLES/"allele_freq_tree.nwk",
    log:
        LOGDIR / "allele_freq_tree_data" / "log.txt"
    script:
        "../scripts/report/allele_freq_tree_data.R"


rule allele_freq_tree_plot:
    conda: "../envs/renv.yaml"
    params:
        design = config["PLOTS"],
        outgroup_id = config["ALIGNMENT_REFERENCE"],
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
        outgroup_id = config["ALIGNMENT_REFERENCE"],
    input:
        tree = report(REPORT_DIR_TABLES/"allele_freq_tree.nwk"),
        metadata = config["METADATA"],
    output:
        table = report(REPORT_DIR_TABLES/"time_signal.csv"),
        json = REPORT_DIR_TABLES/"time_signal.json",
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
        design = config["PLOTS"],
        plot_height_mm = 119.4,
        plot_width_mm = 159.2,
    input: 
        table = OUTDIR/f"{OUTPUT_NAME}.dnds.csv",
    output:
        plot_dn_ds = report(REPORT_DIR_PLOTS/"dn_and_ds.png"),
        plot_omega = report(REPORT_DIR_PLOTS/"dnds.png"),
    log:
        LOGDIR / "evo_plots" / "log.txt"
    script:
        "../scripts/report/dnds_plots.R"


rule af_time_correlation_data:
    conda: "../envs/renv.yaml"
    params:
        cor_method = config["COR"]["METHOD"],
        cor_exact = config["COR"]["EXACT"],
    input:
        variants = OUTDIR/f"{OUTPUT_NAME}.variants.tsv",
        metadata = config["METADATA"],
    output:
        fmt_variants = temp(REPORT_DIR_TABLES/"variants.filled.dated.tsv"),
        correlations = report(REPORT_DIR_TABLES/"af_time_correlation.csv"),
        subset = REPORT_DIR_TABLES/"af_time_correlation.subset.txt",
    log:
        LOGDIR / "af_time_correlation_data" / "log.txt"
    script:
        "../scripts/report/af_time_correlation_data.R"


rule af_time_correlation_plot:
    conda: "../envs/renv.yaml"
    params:
        design = config["PLOTS"],
        plot_height_mm = 119.4,
        plot_width_mm = 159.2,
    input:
        correlations = REPORT_DIR_TABLES/"af_time_correlation.csv",
    output:
        plot = report(REPORT_DIR_PLOTS/"af_time_correlation.png"),
    log:
        LOGDIR / "af_time_correlation_plot" / "log.txt"
    script:
        "../scripts/report/af_time_correlation_plot.R"


rule af_trajectory_panel_plot:
    conda: "../envs/renv.yaml"
    params:
        design = config["PLOTS"],
        n_plot_columns = 3,
        plot_row_height_mm = 35,
        plot_width_mm = 159.2,
        random_color_seed = 7291,
    input:
        fmt_variants = REPORT_DIR_TABLES/"variants.filled.dated.tsv",
        subset = REPORT_DIR_TABLES/"af_time_correlation.subset.txt"
    output:
        plot = report(REPORT_DIR_PLOTS/"af_trajectory_panel.png"),
    log:
        LOGDIR / "af_trajectory_panel_plot" / "log.txt"
    script:
        "../scripts/report/af_trajectory_panel_plot.R"


rule summary_table:
    conda: "../envs/renv.yaml"
    input:
        report = report(OUTDIR/f"{OUTPUT_NAME}.lineage_report.csv"),
        metadata = config["METADATA"]
    output:
        table = REPORT_DIR_TABLES/"summary_table.csv"
    log:
        LOGDIR / "summary_table" / "log.txt"
    script:
        "../scripts/report/summary_table.R"


rule report:
    conda: "../envs/quarto_render.yaml"
    shadow: "shallow"
    input:
        qmd        = Path(config["REPORT_QMD"]).resolve(),
        css        = Path(config["REPORT_CSS"]).resolve(),
        demix      = report(REPORT_DIR_PLOTS/"demix.png"),
        tree_ml    = report(REPORT_DIR_PLOTS/"context_phylogeny.png"),
        diversity  = report(REPORT_DIR_PLOTS/"diversity.png"),
        fig_cor    = report(REPORT_DIR_PLOTS/"polymorphic_sites_over_time.png"),
        SNV        = report(REPORT_DIR_PLOTS/"nv_panel.png"),
        SNV_spike  = report(REPORT_DIR_PLOTS/"nv_panel.S.png"),
        volcano    = report(REPORT_DIR_PLOTS/"af_time_correlation.png"),
        panel      = report(REPORT_DIR_PLOTS/"af_trajectory_panel.png"),
        tree       = report(REPORT_DIR_PLOTS/"allele_freq_tree.png"),
        temest     = report(REPORT_DIR_PLOTS/"time_signal.png"),
        evo        = report(REPORT_DIR_PLOTS/"dn_and_ds.png"),
        omega_plot = report(REPORT_DIR_PLOTS/"dnds.png"),
        heat_table = report(OUTDIR/"vaf"/"pairwise_trajectory_correlation.csv"),
        freyja_ts  = OUTDIR/"demixing"/"freyja_data"/"last_barcode_update.txt",
        value      = REPORT_DIR_TABLES/"diversity.json",
        stats_lm   = REPORT_DIR_TABLES/"time_signal.json",
        stats_ml   = REPORT_DIR_TABLES/"context_phylogeny.json",
        table      = REPORT_DIR_TABLES/"summary_table.csv",
        sum_nv     = REPORT_DIR_TABLES/"nv_panel_summary.json",
    params:
        workflow_version = get_repo_version(BASE_PATH.as_posix(), __version__),
        min_ivar_freq = config["VC"]["MIN_FREQ"],
        ufboot_reps = config["UFBOOT"]["REPS"],
        shalrt_reps = config["SHALRT"]["REPS"],
        name = config["OUTPUT_NAME"],
        use_bionj = config["USE_BIONJ"],
        cor_method = config["COR"]["METHOD"],
    output:
        html = report(OUTDIR/f"{OUTPUT_NAME}.report.html")
    log:
        LOGDIR / "report" / "log.txt"
    shell:
        "set +o pipefail; "
        "Rscript -e \"quarto::quarto_render("
            "input = '{input.qmd:q}', "
            "execute_params=list("
                "css='{input.css:q}', "
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
