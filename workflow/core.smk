BASE_PATH = Path(workflow.basedir).parent.resolve()

include: "rules/common.smk"

# Outputs
OUTPUT_NAME = config["OUTPUT_NAME"]
OUTDIR = Path(config["OUTPUT_DIRECTORY"])

# Logging
LOGDIR = OUTDIR / "logs"

# Report
REPORT_DIR_PLOTS = Path(OUTDIR/"report/plots")
REPORT_DIR_TABLES = Path(OUTDIR/"report/tables")

include: "rules/fetch.smk"
include: "rules/fasta.smk"
include: "rules/asr.smk"
include: "rules/demix.smk"
include: "rules/vaf.smk"
include: "rules/pangolin.smk"
include: "rules/distances.smk"
include: "rules/evolution.smk"
