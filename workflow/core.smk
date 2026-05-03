BASE_PATH = Path(workflow.basedir).parent.resolve()

include: "rules/common.smk"

# Logging
LOGDIR = Path("<results>/<logs>")

# Report
REPORT_DIR_PLOTS = Path("<results>/<dataset>/report/plots")
REPORT_DIR_TABLES = Path("<results>/<dataset>/report/tables")

include: "rules/fetch.smk"
include: "rules/fasta.smk"
include: "rules/asr.smk"
include: "rules/vaf.smk"
include: "rules/sites.smk"
include: "rules/distances.smk"
include: "rules/evolution.smk"
