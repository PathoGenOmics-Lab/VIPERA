#!/usr/bin/env python

"""build_targets.py

This script simplifies the process of creating the configuration file.
It takes a list of sample names, a directory with BAM and FASTA files, the path
to the metadata table and the name of your dataset as required inputs. Then,
it searches the directory for files that have the appropriate extensions
and sample names and adds them to the YAML configuration file.
"""

import sys
import argparse
from pathlib import Path
from typing import List

import yaml


def find_file_with_extension(directory: Path, prefix: str, extensions: List[str]) -> str:
    candidate_files = []
    for path in directory.rglob(f"*"):
        if any(path.name.endswith(ext) for ext in extensions) and path.name.startswith(prefix):
            candidate_files.append(path)
    if len(candidate_files) == 1:
        return candidate_files[0].as_posix()
    else:
        sys.exit(f"ERROR: {len(candidate_files)} candidates found in '{directory}' for prefix '{prefix}' with extensions {extensions}: {candidate_files}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "-i", "--sample-directory",
        help="Directory containing sample sequencing data",
        required=True,
        type=Path
    )
    parser.add_argument(
        "-s", "--sample-names",
        nargs="+",
        help="Sample names to look for in the sample directory",
        required=True
    )
    parser.add_argument(
        "-b", "--bam-extensions",
        nargs="+",
        help="File extensions for BAM files",
        required=False,
        default=[".trim.sort.bam"]
    )
    parser.add_argument(
        "-f", "--fasta-extensions",
        nargs="+",
        help="File extensions for FASTA files",
        required=False,
        default=[".fa", ".fasta"]
    )
    parser.add_argument(
        "-n", "--output-name",
        help="Dataset name",
        required=False,
        default="case_study"
    )
    parser.add_argument(
        "-m", "--metadata-csv",
        help="Dataset metadata CSV file",
        required=True,
        type=Path
    )
    parser.add_argument(
        "-o", "--output-yaml",
        help="Output YAML target configuration file for VIPERA",
        required=False,
        default="targets.yaml"
    )
    args = parser.parse_args()
    
    # Build targets
    targets = {"SAMPLES": {}}
    for sample_name in args.sample_names:
        targets["SAMPLES"][sample_name] = {}
        targets["SAMPLES"][sample_name]["bam"] = find_file_with_extension(args.sample_directory, sample_name, args.bam_extensions)
        targets["SAMPLES"][sample_name]["fasta"] = find_file_with_extension(args.sample_directory, sample_name, args.fasta_extensions)
    
    # Write empty fields
    if args.metadata_csv.is_file():
        targets["METADATA"] = args.metadata_csv.as_posix()
    else:
        sys.exit(f"ERROR: metadata file '{args.metadata_csv}' does not exist")
    targets["OUTPUT_NAME"] = args.output_name
    targets["OUTPUT_DIRECTORY"] = "output"
    targets["CONTEXT_FASTA"] = None
    targets["MAPPING_REFERENCES_FASTA"] = None

    # Write output
    with open(args.output_yaml, "w") as fw:
        yaml.dump(targets, fw, indent=2)
