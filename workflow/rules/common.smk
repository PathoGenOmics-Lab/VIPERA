def iter_files(file_type: str):
    for sample_name in config["SAMPLES"].keys():
        yield config["SAMPLES"][sample_name][file_type]


def iter_samples():
    for sample_name in config["SAMPLES"].keys():
        yield sample_name


def get_input_bam(wildcards):
    return config["SAMPLES"][wildcards.sample]["bam"]


def get_input_fasta(wildcards):
    return config["SAMPLES"][wildcards.sample]["fasta"]


def get_version_str(base_dir: str) -> str:
    if not Path(f"{base_dir}/.git").is_dir():
        print("Not a git repository!")
        return "no version available"
    last_tag_description = subprocess.check_output(
        ["git", f"--git-dir={base_dir}/.git", "describe", "--always"]
    ).strip().decode("utf-8")
    if last_tag_description.startswith("v"):
        return last_tag_description
    else:
        return f"commit {last_tag_description}"
