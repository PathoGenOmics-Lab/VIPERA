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


def get_repo_version(base_dir: str, default: str, warn=False) -> str:
    """Determines the workflow version from the git repository status.
    It first checks if the workflow is a git repository by looking for
    a .git folder. If it is, it reads the latest tag from the .git folder
    and compares it with the default version provided as an argument.
    It returns one of the following strings:
    - 'v{version from git tag}', if the default version matches the latest git tag
    - 'v{default version} ({git description})', if the default version differs from the latest git tag
    - 'commit {commit hash}', if the repository has no tags
    - 'v{default version} (no git)', if there is no git repository
    """
    try:
        last_tag_description = subprocess.check_output(
            f"git --git-dir={base_dir}/.git describe --always",
            shell=True
        ).strip().decode("utf-8")
        if last_tag_description.startswith("v"):
            if last_tag_description[1:] == default:
                return last_tag_description
            else:
                return f"v{default} ({last_tag_description})"
        else:
            return f"commit {last_tag_description}"
    except subprocess.CalledProcessError as e:
        if warn:
            print(f"Repository tag not found: '{e}'")
        return f"v{default} (no git)"
