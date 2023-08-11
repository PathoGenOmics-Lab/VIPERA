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
    try:
        last_tag_description = subprocess.check_output(
            ["git", f"--git-dir={base_dir}/.git", "describe", "--always"],
            shell=True
        ).strip().decode("utf-8")
    except subprocess.CalledProcessError as e:
        print(f"Version not found. Error: '{e}'")
        return "N/A"
    if last_tag_description.startswith("v"):
        return last_tag_description
    else:
        return f"commit {last_tag_description}"
