#!/usr/bin/env python3
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from argparse import ArgumentParser
import shutil
import sys

@dataclass
class PrepareGreasyArrayJobArgs:
    basedir: str
    generate_config_cmd: str
    reference_genome: str
    krakendb: str
    reads_directory: str

    ...
    @classmethod
    def get_arguments(cls, args=sys.argv[1:]) -> PrepareGreasyArrayJobArgs:
        parser = ArgumentParser(
            description="Prepares the greasy array job to submit multiple meTAline jobs in the cluster."
        )
        parser.add_argument(
            "--basedir",
            type=str,
            required=True,
            help="Base directory in which the pipeline will create output"
        )
        parser.add_argument(
            "--generate_config_cmd",
            type=str,
            default="singularity run --cleanenv ./metaline.sif metaline-generate-config",
            help="The external command to generate the config files.",
        )
        parser.add_argument(
            "--reference_genome",
            type=str,
            required=True,
            default="singularity run --cleanenv ./metaline.sif metaline-generate-config",
            help="Path to the indexed reference genome.",
        )
        parser.add_argument(
            "--krakendb",
            type=str,
            required=True,
            default="singularity run --cleanenv ./metaline.sif metaline-generate-config",
            help="Path to the decompressed Kraken database.",
        )
        parser.add_argument(
            "--reads_directory",
            type=str,
            required=True,
            default="singularity run --cleanenv ./metaline.sif metaline-generate-config",
            help="Path to the directory containing the reads.",
        )
        return PrepareGreasyArrayJobArgs(**vars(parser.parse_args(args)))

def extract_reads_prefixes(reads_directory: Path | str) -> list[str]:
    reads_prefixes: list[str] = []
    # Define the directory containing the files
    reads_dir = Path(reads_directory)  # Update with your actual path
    # Find and sort matching files
    files = sorted(reads_dir.glob("*[!a-zA-Z0-9]1.f*q.gz"))
    # Process each file
    for file_path in files:
        file_name = file_path.name  # Extract filename
        file_prefix = file_name.rsplit(".", 2)[0].rsplit("1", 1)[0].rstrip("._-")  # Extract prefix
        reads_prefixes.append(file_prefix)
    return reads_prefixes


def main() -> int:
    args = PrepareGreasyArrayJobArgs.get_arguments()
    # 1. make config files (./_scripts/make_configs.sh)
    config_file_output_dir = Path(args.basedir) / "configs"
    shutil.rmtree(str(config_file_output_dir), ignore_errors=True)
    config_file_output_dir.mkdir(parents=True, exist_ok=True)
    reads_prefixes = extract_reads_prefixes(args.reads_directory)

    # 2. make job list and greasy list (./_scripts/make_joblist_and_greasy_list.sh)
    # 3. [OPTIONAL] submit directly the .job file
    return 0

if __name__ == "__main__":
    raise SystemExit(main())