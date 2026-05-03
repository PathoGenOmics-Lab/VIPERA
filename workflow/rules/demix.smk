rule demix_barcode_update:
    threads: 1
    shadow: "shallow"
    conda:
        "../envs/freyja.yaml"
    params:
        pathogen = config["DEMIX"]["PATHOGEN"]
    output:
        folder = directory("<results>/<dataset>/demixing/freyja_data"),
        curated_lineages = "<results>/<dataset>/demixing/freyja_data/curated_lineages.json",
        last_barcode_update = "<results>/<dataset>/demixing/freyja_data/last_barcode_update.txt",
        lineage_mutations = "<results>/<dataset>/demixing/freyja_data/lineage_mutations.json",
        lineage_yml = "<results>/<dataset>/demixing/freyja_data/lineages.yml",
        pathogens = "<results>/<dataset>/demixing/freyja_data/pathogen_config.yml",
        usher_barcodes = "<results>/<dataset>/demixing/freyja_data/usher_barcodes.feather"
    log:
        LOGDIR / "demix_barcode_update" / "log.txt"
    shell:
        "mkdir -p {output.folder:q} && "
        "freyja update --outdir {output.folder:q} --pathogen {params.pathogen:q} >{log} 2>&1"


rule demix_preprocessing:
    threads: 1
    conda: "../envs/var_calling.yaml"
    shadow: "minimal"
    input:
        bam = get_input_bam,
        ref_fasta = lambda wildcards: select_mapping_references_fasta()
    params:
        max_depth = config["DEMIX"]["MAX_DEPTH"],
        minq = config["DEMIX"]["MIN_QUALITY"],
    output:
        depth_file = "<results>/<dataset>/demixing/{sample}/{sample}_depth.txt",
        variants_file = "<results>/<dataset>/demixing/{sample}/{sample}_variants.tsv",
    log:
        pileup = LOGDIR / "demix_preprocessing" / "{sample}_pileup.log.txt",
        ivar = LOGDIR / "demix_preprocessing" / "{sample}_ivar.log.txt",
    shell:
        "set -euo pipefail && "
        "samtools mpileup -aa -A -d {params.max_depth} -Q {params.minq} -q 0 -B -f {input.ref_fasta:q} {input.bam:q} >sample.pileup 2>{log.pileup:q} && "
        "ivar variants -p variants -q {params.minq} -r {input.ref_fasta:q} >{log.ivar:q} 2>&1 <sample.pileup && "
        "cut -f1-4 sample.pileup >{output.depth_file:q} && "
        "mv variants.tsv {output.variants_file:q}"


rule demix:
    threads: 1
    conda: "../envs/freyja.yaml"
    shadow: "minimal"
    input:
        depth_file = "<results>/<dataset>/demixing/{sample}/{sample}_depth.txt",
        variants_file = "<results>/<dataset>/demixing/{sample}/{sample}_variants.tsv",
        barcodes = "<results>/<dataset>/demixing/freyja_data/usher_barcodes.feather",
        curated_lineages = "<results>/<dataset>/demixing/freyja_data/curated_lineages.json",
        lineage_yml = "<results>/<dataset>/demixing/freyja_data/lineages.yml",
    params:
        coverage_cutoff = config["DEMIX"]["COV_CUTOFF"],
        minimum_abundance = config["DEMIX"]["MIN_ABUNDANCE"],
        confirmed_only = "--confirmedonly " if config["DEMIX"]["CONFIRMED_ONLY"] else "",
        depth_cutoff = config["DEMIX"]["DEPTH_CUTOFF"],
        auto_adapt = "--autoadapt " if config["DEMIX"]["AUTO_ADAPT"] else "",
        relaxed_mrca = "--relaxedmrca " if config["DEMIX"]["RELAXED_MRCA"] else "",
        relaxed_mrca_thresh = config["DEMIX"]["RELAXED_MRCA_THRESH"],
        solver = config["DEMIX"]["SOLVER"],
    output:
        demix_file = "<results>/<dataset>/demixing/samples/{sample}/{sample}_demixed.tsv"
    log:
        LOGDIR / "demix" / "{sample}.log.txt"
    shell:
        "freyja demix "
        "{input.variants_file:q} "
        "{input.depth_file:q} "
        "--barcodes {input.barcodes:q} "
        "--meta {input.curated_lineages:q} "
        "--lineageyml {input.lineage_yml:q} "
        "--eps {params.minimum_abundance} "
        "--covcut {params.coverage_cutoff} "
        "--depthcutoff {params.depth_cutoff} "
        "{params.confirmed_only}"
        "{params.auto_adapt}"
        "{params.relaxed_mrca}"
        "--relaxedthresh {params.relaxed_mrca_thresh} "
        "--solver {params.solver} "
        "--output {output.demix_file} "
        "--max-solver-threads {threads} "
        ">{log} 2>&1"


rule summarise_demix:
    threads: 1
    conda: "../envs/renv.yaml"
    shadow: "shallow"
    input:
        tables = expand("<results>/<dataset>/demixing/samples/{sample}/{sample}_demixed.tsv", sample=iter_samples())
    output:
        summary_df = report("<results>/<dataset>/demixing/summary.csv")
    log:
        LOGDIR / "summarise_demix" / "log.txt"
    script: 
        "../scripts/summarise_demix.R"
