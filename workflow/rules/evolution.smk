rule filter_genbank_features:
    threads: 1
    conda: "../envs/biopython.yaml"
    params:
        included = config.get("GB_FEATURES", {}).get("INCLUDE", {}),
        excluded = config.get("GB_FEATURES", {}).get("EXCLUDE", {}),
    input:
        gb = "<results>/<dataset>/reference.gb",
    output:
        gb = "<results>/<dataset>/reference.cds.gb",
    log:
        "<logs>/<dataset>/filter_genbank_features/log.txt"
    script:
        "../scripts/filter_genbank_features.py"


rule n_s_sites:
    threads: 1
    conda: "../envs/biopython.yaml"
    params:
        gb_qualifier_display = "gene",
    input:
        fasta = "<results>/<dataset>/ancestor.fasta",
        masked = "<results>/<dataset>/sites_masked.bed",
        gb = "<results>/<dataset>/reference.cds.gb",
        genetic_code = Path(config["GENETIC_CODE_JSON"]).resolve(),
    output:
        csv = temp("<results>/<dataset>/ancestor.n_s.sites.csv"),
    log:
        "<logs>/<dataset>/n_s_sites/log.txt"
    script:
        "../scripts/n_s_sites_from_fasta.py"


rule calculate_dnds:
    conda: "../envs/renv.yaml"
    input: 
        n_s_sites = "<results>/<dataset>/ancestor.n_s.sites.csv",
        masked = "<results>/<dataset>/sites_masked.bed",
        variants =  "<results>/<dataset>/variants.tsv",
        metadata = config["METADATA"]
    output:
        table = "<results>/<dataset>/dnds.csv",
    log:
        "<logs>/<dataset>/calculate_dnds/log.txt"
    script:
        "../scripts/calculate_dnds.R"
