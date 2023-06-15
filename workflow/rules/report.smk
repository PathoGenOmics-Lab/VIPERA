rule report:
    conda: "../envs/renv.yaml"
    params:
        id = expand("{sample2}",sample2 = iter_samples()),
        metadata = config["METADATA"],
        context_tree = config["ML_TREE"], 
        design = config["PLOTS"],
        diversity_alignment = config["DIVERSITY"] ,
        nsp = config["NSP"],
        qmd = "case_study.report.qmd"
    input:
        vcf = OUTDIR/f"{OUTPUT_NAME}.masked.filtered.tsv",
        depth_files = expand(OUTDIR/"demixing/{sample}/{sample}_depth.txt", sample = iter_samples()),
        summary_demixing =  OUTDIR/"summary_freyja_demixing.csv",
        window = OUTDIR/f"{OUTPUT_NAME}.window.csv", 
        distancias = OUTDIR/f"{OUTPUT_NAME}.weighted_distances.csv",
        N_S = OUTDIR/f"{OUTPUT_NAME}.ancestor.N_S.sites.csv"
    output:
        html = f"{OUTPUT_NAME}.report.html"
    shell:         
        """
         set +o pipefail
        Rscript -e 'library(quarto)' -e \"quarto_render(input = '{params.qmd}',\
                                           output_file = '{output.html}',\
                                           execute_params=list(id='{params.id}',\
                                                       metadata='{params.metadata}',\
                                                       context_tree = '{params.context_tree}',\
                                                       design = '{params.design}',\
                                                       diversity_alignment = '{params.diversity_alignment}',\
                                                       nsp = '{params.nsp}',
                                                       vcf = '{input.vcf}',\
                                                       depth_files = '{input.depth_files}',\
                                                       summary_demixing = '{input.summary_demixing}',\
                                                       window = '{input.window}',\
                                                       distancias = '{input.distancias}',\
                                                       N_S = '{input.N_S}'))\"    
        """


rule window:
    conda: "../envs/biopython.yaml"
    input:
        vcf = OUTDIR/f"{OUTPUT_NAME}.masked.filtered.tsv",
    output:
        window_df = OUTDIR/f"{OUTPUT_NAME}.window.csv",
    script:
        "../scripts/window.py"
