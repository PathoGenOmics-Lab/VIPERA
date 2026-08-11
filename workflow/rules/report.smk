rule demix_plot_data:
    conda: "../envs/renv.yaml"
    input:
        summary_demixing = "<results>/<dataset>/demixing/summary.csv",
        metadata = config["METADATA"]
    output:
        data = "<results>/<dataset>/report/tables/demix.csv"
    log:
        "<logs>/<dataset>/demix_plot_data/log.txt"
    script:
        "../scripts/report/demix_plot_data.R"


rule demix_plot:
    conda: "../envs/renv.yaml"
    params:
        design = config["PLOTS"],
        plot_width_mm = 159.2,
        plot_height_mm = 119.4,
    input:
        data = "<results>/<dataset>/report/tables/demix.csv"
    output:
        plot = report("<results>/<dataset>/report/plots/demix.png")
    log:
        "<logs>/<dataset>/demix_plot/log.txt"
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
        study_fasta = "<results>/<dataset>/aligned.masked.fasta",
        context_fasta = "<results>/<dataset>/context/nextalign/context_sequences.aligned.masked.fasta",
    output:
        divs = "<results>/<dataset>/report/tables/diversity.txt",
        json = "<results>/<dataset>/report/tables/diversity.json",
    log:
        "<logs>/<dataset>/diversity_data/log.txt"
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
        divs = "<results>/<dataset>/report/tables/diversity.txt",
        json = "<results>/<dataset>/report/tables/diversity.json",
    output:
        plot = report("<results>/<dataset>/report/plots/diversity.png"),
    log:
        "<logs>/<dataset>/diversity_plot/log.txt"
    script:
        "../scripts/report/diversity_plot.R"


rule extract_genbank_regions:
    conda: "../envs/biopython.yaml"
    params:
        gb_qualifier = "gene",
    input:
        gb = "<results>/<dataset>/reference.cds.gb",
    output:
        regions = temp("<results>/<dataset>/report/tables/genbank_regions.json"),
    log:
        "<logs>/<dataset>/extract_genbank_regions/log.txt"
    script:
        "../scripts/report/extract_genbank_regions.py"


rule polymorphic_sites_over_time_plot:
    conda: "../envs/renv.yaml"
    params:
        design = config["PLOTS"],
        plot_width_mm = 159.2,
        plot_height_mm = 119.4,
    input:
        table = "<results>/<dataset>/report/tables/polymorphic_sites_over_time.csv",
    output:
        plot = report("<results>/<dataset>/report/plots/polymorphic_sites_over_time.png"),
    log:
        "<logs>/<dataset>/polymorphic_sites_over_time_plot/log.txt"
    script:
        "../scripts/report/polymorphic_sites_over_time_plot.R"


rule nv_panel_plot:
    conda: "../envs/renv.yaml"
    params:
        design = config["PLOTS"],
        window_step = config["WINDOW"]["STEP"],
        plot_height_mm = 250.0,
        plot_width_mm = 240.0,
    input:
        panel = "<results>/<dataset>/report/tables/nv_panel.csv",
        window = "<results>/<dataset>/report/tables/window.csv",
        regions = "<results>/<dataset>/report/tables/genbank_regions.json",
        highlight_window_regions = config["PLOT_GENOME_REGIONS"],
    output:
        plot = report("<results>/<dataset>/report/plots/nv_panel.png"),
    log:
        "<logs>/<dataset>/nv_panel_plot/log.txt"
    script:
        "../scripts/report/nv_panel_plot.R"


use rule nv_panel_plot as nv_panel_plot_S with:
    input:
        panel = "<results>/<dataset>/report/tables/nv_panel.S.csv",
        window = "<results>/<dataset>/report/tables/window.S.csv",
        regions = "<results>/<dataset>/report/tables/genbank_regions.json",
        highlight_window_regions = "<results>/<dataset>/empty.txt",
    output:
        plot = report("<results>/<dataset>/report/plots/nv_panel.S.png"),
    log:
        "<logs>/<dataset>/nv_panel_plot_S/log.txt"


rule merge_json_files:
    input:
        "<results>/<dataset>/report/tables/nv_panel.json",
        "<results>/<dataset>/report/tables/polymorphic_sites_over_time.json",
        "<results>/<dataset>/report/tables/window.json",
    output:
        json = "<results>/<dataset>/report/tables/nv_panel_summary.json",
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
        target_fasta = "<results>/<dataset>/sequences.fasta",
        tree = "<results>/<dataset>/tree_context/context.treefile",
    output:
        json = "<results>/<dataset>/report/tables/context_phylogeny.json",
        annotation = "<results>/<dataset>/report/tables/context_phylogeny.csv",
    log:
        "<logs>/<dataset>/context_phylogeny_data/log.txt"
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
        tree = "<results>/<dataset>/tree_context/context.treefile",
        json = "<results>/<dataset>/report/tables/context_phylogeny.json",
        annotation = "<results>/<dataset>/report/tables/context_phylogeny.csv"
    output:
        plot = report("<results>/<dataset>/report/plots/context_phylogeny.png"),
    log:
        "<logs>/<dataset>/context_phylogeny_plot/log.txt"
    script:
        "../scripts/report/context_phylogeny_plot.R"


rule allele_freq_tree_plot:
    conda: "../envs/renv.yaml"
    params:
        design = config["PLOTS"],
        outgroup_id = "case_ancestor",
        plot_height_mm = 119.4,
        plot_width_mm = 159.2,
    input:
        tree = report("<results>/<dataset>/report/tables/allele_freq_tree.nwk"),
        study_fasta = "<results>/<dataset>/sequences.fasta",
        metadata = config["METADATA"],
    output:
        plot = report("<results>/<dataset>/report/plots/allele_freq_tree.png"),
    log:
        "<logs>/<dataset>/allele_freq_tree_plot/log.txt"
    script:
        "../scripts/report/allele_freq_tree_plot.R"


rule time_signal_plot:
    conda: "../envs/renv.yaml"
    params:
        design = config["PLOTS"],
        plot_height_mm = 119.4,
        plot_width_mm = 159.2,
    input:
        table = report("<results>/<dataset>/report/tables/time_signal.csv"),
    output:
        plot = report("<results>/<dataset>/report/plots/time_signal.png"),
    log:
        "<logs>/<dataset>/time_signal_plot/log.txt"
    script:
        "../scripts/report/time_signal_plot.R"


rule dnds_plots:
    conda: "../envs/renv.yaml"
    params: 
        design = config["PLOTS"],
        plot_height_mm = 119.4,
        plot_width_mm = 159.2,
    input: 
        table = "<results>/<dataset>/dnds.csv",
    output:
        plot_dn_ds = report("<results>/<dataset>/report/plots/dn_and_ds.png"),
        plot_omega = report("<results>/<dataset>/report/plots/dnds.png"),
    log:
        "<logs>/<dataset>/evo_plots/log.txt"
    script:
        "../scripts/report/dnds_plots.R"


rule af_time_correlation_plot:
    conda: "../envs/renv.yaml"
    params:
        design = config["PLOTS"],
        plot_height_mm = 119.4,
        plot_width_mm = 159.2,
    input:
        correlations = "<results>/<dataset>/report/tables/af_time_correlation.csv",
    output:
        plot = report("<results>/<dataset>/report/plots/af_time_correlation.png"),
    log:
        "<logs>/<dataset>/af_time_correlation_plot/log.txt"
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
        fmt_variants = "<results>/<dataset>/report/tables/variants.filled.dated.tsv",
        subset = "<results>/<dataset>/report/tables/af_time_correlation.subset.txt"
    output:
        plot = report("<results>/<dataset>/report/plots/af_trajectory_panel.png"),
    log:
        "<logs>/<dataset>/af_trajectory_panel_plot/log.txt"
    script:
        "../scripts/report/af_trajectory_panel_plot.R"


rule summary_table:
    conda: "../envs/renv.yaml"
    input:
        report = report("<results>/<dataset>/lineage_report.csv"),
        metadata = config["METADATA"]
    output:
        table = "<results>/<dataset>/report/tables/summary_table.csv"
    log:
        "<logs>/<dataset>/summary_table/log.txt"
    script:
        "../scripts/report/summary_table.R"


rule report:
    conda: "../envs/quarto_render.yaml"
    shadow: "minimal"
    input:
        qmd        = Path(config["REPORT_QMD"]).resolve(),
        css        = Path(config["REPORT_CSS"]).resolve(),
        demix      = report("<results>/<dataset>/report/plots/demix.png"),
        tree_ml    = report("<results>/<dataset>/report/plots/context_phylogeny.png"),
        diversity  = report("<results>/<dataset>/report/plots/diversity.png"),
        fig_cor    = report("<results>/<dataset>/report/plots/polymorphic_sites_over_time.png"),
        SNV        = report("<results>/<dataset>/report/plots/nv_panel.png"),
        SNV_spike  = report("<results>/<dataset>/report/plots/nv_panel.S.png"),
        volcano    = report("<results>/<dataset>/report/plots/af_time_correlation.png"),
        panel      = report("<results>/<dataset>/report/plots/af_trajectory_panel.png"),
        tree       = report("<results>/<dataset>/report/plots/allele_freq_tree.png"),
        temest     = report("<results>/<dataset>/report/plots/time_signal.png"),
        evo        = report("<results>/<dataset>/report/plots/dn_and_ds.png"),
        omega_plot = report("<results>/<dataset>/report/plots/dnds.png"),
        heat_table = report("<results>/<dataset>/report/tables/pairwise_trajectory_correlation_matrix.csv"),
        freyja_ts  = "<results>/<dataset>/demixing/freyja_data/last_barcode_update.txt",
        value      = "<results>/<dataset>/report/tables/diversity.json",
        stats_lm   = "<results>/<dataset>/report/tables/time_signal.json",
        stats_ml   = "<results>/<dataset>/report/tables/context_phylogeny.json",
        table      = "<results>/<dataset>/report/tables/summary_table.csv",
        sum_nv     = "<results>/<dataset>/report/tables/nv_panel_summary.json",
    params:
        workflow_version = get_repo_version(Path(workflow.basedir).parent, __version__),
        min_ivar_freq = config["VC"]["MIN_FREQ"],
        ufboot_reps = config["UFBOOT"]["REPS"],
        shalrt_reps = config["SHALRT"]["REPS"],
        name = config["OUTPUT_NAME"],
        use_bionj = config["USE_BIONJ"],
        cor_method = config["COR"]["METHOD"],
    output:
        html = report("<results>/<dataset>/report.html")
    log:
        "<logs>/<dataset>/report/log.txt"
    shell:
        """
        set +o pipefail
        exec >{log} && exec 2>&1

        printf "%s\n" \
            "ufboot_reps: '{params.ufboot_reps}'" \
            "shalrt_reps: '{params.shalrt_reps}'" \
            "min_ivar_freq: '{params.min_ivar_freq}'" \
            "workflow_version: '{params.workflow_version}'" \
            "use_bionj: '{params.use_bionj}'" \
            "cor_method: '{params.cor_method}'" \
            "div: '{input.diversity}'" \
            "demix: '{input.demix}'" \
            "tree: '{input.tree}'" \
            "tempest: '{input.temest}'" \
            "SNV: '{input.SNV}'" \
            "SNV_s: '{input.SNV_spike}'" \
            "evo: '{input.evo}'" \
            "div_value: '{input.value}'" \
            "panel: '{input.panel}'" \
            "volcano: '{input.volcano}'" \
            "tree_ml: '{input.tree_ml}'" \
            "fig_cor_snp: '{input.fig_cor}'" \
            "stats_lm: '{input.stats_lm}'" \
            "stats_ml: '{input.stats_ml}'" \
            "table: '{input.table}'" \
            "sum_nv: '{input.sum_nv}'" \
            "heat_tab: '{input.heat_table}'" \
            "omega_plot: '{input.omega_plot}'" \
            "freyja_ts: '{input.freyja_ts}'" \
            "name: '{params.name}'" \
            >params.yaml

        sed "s|__CSSPLACEHOLDER__|{input.css}|g" {input.qmd:q} >report.qmd

        quarto render report.qmd --execute-params params.yaml
        mv report.html {output.html:q}
        """
