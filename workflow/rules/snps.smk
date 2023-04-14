rule polymorphic_sites:
    threads: 1
    conda: "../envs/biopython.yaml"
    params:
        indeterminate_char = "N",
        include_indeterminations = False,
        omit_record_IDs = [config["REFSEQ_REFERENCE"]]  # Omit sites from these records
    input:
        reference_fasta = OUTDIR/f"{OUTPUT_NAME}.ancestor.fasta",
        alignment_fasta = OUTDIR/"nextalign"/f"{OUTPUT_NAME}.aligned.fasta"
    output:
        table = OUTDIR/f"{OUTPUT_NAME}.polymorphic_sites.csv"
    script:
        "scripts/polymorphic_sites.py"


rule plot_polymorphisms:
    threads: 1
    conda: "../envs/renv.yaml"
    params:
        width = 30,
        height = 15,
        units = "in"
    input:
        table = OUTDIR/f"{OUTPUT_NAME}.polymorphic_sites.csv",
        metadata = METADATA
    output:
        plot = OUTDIR/f"{OUTPUT_NAME}.polymorphic_sites.pdf"
    script:
        "scripts/plot_polymorphisms.R"
