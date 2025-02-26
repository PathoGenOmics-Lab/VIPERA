def iter_files(file_type: str):
    for sample_name in config.get("SAMPLES", {}).keys():
        yield config.get("SAMPLES", {})[sample_name][file_type]


def iter_samples():
    for sample_name in config.get("SAMPLES", {}).keys():
        yield sample_name


def get_input_bam(wildcards):
    return config.get("SAMPLES", {})[wildcards.sample]["bam"]


def get_input_fasta(wildcards):
    return config.get("SAMPLES", {})[wildcards.sample]["fasta"]


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
            shell=True,
            stderr=subprocess.DEVNULL
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


def select_context_fasta():
    """Set context to be fetched automatically if CONTEXT_FASTA=null"""
    if "CONTEXT_FASTA" not in config.keys() or config["CONTEXT_FASTA"] is None:
        return OUTDIR/"context"/"sequences.fasta"
        if not Path(config["GISAID"]["CREDENTIALS"]).is_file():
            raise FileNotFoundError(f"Tried to download a context dataset, but no GISAID credentials were found at '{config['GISAID']['CREDENTIALS']}' (see README.md).")
    elif Path(config["CONTEXT_FASTA"]).is_file():
        return config["CONTEXT_FASTA"]
    else:
        raise FileNotFoundError(f"No context FASTA file was found at '{config['CONTEXT_FASTA']}' (see README.md).")


def select_mapping_references_fasta():
    """Set mapping references to be fetched automatically if MAPPING_REFERENCES_FASTA=null"""
    if "MAPPING_REFERENCES_FASTA" not in config.keys() or config["MAPPING_REFERENCES_FASTA"] is None:
        return OUTDIR/"mapping_references.fasta"
    elif Path(config["MAPPING_REFERENCES_FASTA"]).is_file():
        return config["MAPPING_REFERENCES_FASTA"]
    else:
        raise FileNotFoundError(f"No mapping references FASTA file was found at '{config['MAPPING_REFERENCES_FASTA']}' (see README.md).")
