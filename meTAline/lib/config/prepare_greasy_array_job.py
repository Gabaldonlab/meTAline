#!/usr/bin/env python3
from __future__ import annotations

import shutil
import sys
from argparse import ArgumentParser
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from subprocess import getstatusoutput
from typing import Iterable, TypedDict


@dataclass
class PrepareGreasyArrayJobArgs:
    metaline_rule: str
    prefix_fr: list[str]
    basedir: str
    generate_config_cmd: str
    metaline_cmd: str
    greasy_cmd: str
    reference_genome: str
    krakendb: str
    reads_directory: str
    fastq_extension: str
    metaphlan_db: str
    metaphlan_index: str
    n_db: str
    protein_db: str
    max_workers: int
    joblist_size: int
    cluster_node_cores: int
    trimmo_cores: int
    hisat2_cores: int
    kraken2_cores: int
    snakemake_cores: int

    @classmethod
    def get_arguments(cls, args=sys.argv[1:]) -> PrepareGreasyArrayJobArgs:
        parser = ArgumentParser(
            description="Prepares the greasy array job to submit multiple meTAline jobs in the cluster."
        )
        parser.add_argument(
            "--metaline-rule",
            type=str,
            choices=[
                "all",
                "trimming",
                "host_depletion",
                "kmer_taxonomy",
                "RAnalysis",
                "BioBakery",
            ],
            default="all",
            help="The meTAline rule to use for the subcommands of the greasy array.",
        )
        parser.add_argument(
            "--prefix-fr",
            nargs="+",
            default=("_",),
            help=(
                "List of prefixes of the forward-reverse annotations in the "
                'filenames. (E.G.: \'--prefix-fr ".R" "_"\' for '
                "'PONSJIB_191.H100_1' and 'PONSJIB_191.H200.R1')"
            ),
        )
        parser.add_argument(
            "--basedir",
            type=str,
            required=True,
            help="Base directory in which the pipeline will create output",
        )
        parser.add_argument(
            "--generate_config_cmd",
            type=str,
            default="metaline-generate-config",  # Simply this way, because it will try to get it from inside the Singularity image!
            help="The external command to generate the config files.",
        )
        parser.add_argument(
            "--metaline-cmd",
            type=str,
            default="singularity run --cleanenv ./metaline.sif metaline",
            help="The external command to run metaline with the config files.",
        )
        parser.add_argument(
            "--greasy-cmd",
            type=str,
            default="module load greasy && greasy/2.2.4.1",
            help="The external command to run greasy with the job files.",
        )
        parser.add_argument(
            "--reference-genome",
            type=str,
            required=True,
            help="Path to the indexed reference genome.",
        )
        parser.add_argument(
            "--krakendb",
            type=str,
            required=True,
            help="Path to the decompressed Kraken database.",
        )
        parser.add_argument(
            "--reads-directory",
            type=str,
            required=True,
            help="Path to the directory containing the reads.",
        )
        parser.add_argument(
            "--fastq-extension",
            type=str,
            default="fq.gz",
            help="Extension of the fastq files.",
        )
        parser.add_argument(
            "--metaphlan-db",
            type=str,
            required=True,
            help="Path to the Metaphlan database directory.",
        )
        parser.add_argument(
            "--metaphlan-index",
            type=str,
            required=True,
            help="Name of the index in the given Metaphlan database directory.",
        )
        parser.add_argument(
            "--n-db",
            type=str,
            required=True,
            help="Path to the database to do the nucleotide search. (*.tar.gz)",
        )
        parser.add_argument(
            "--protein-db",
            type=str,
            required=True,
            help="Path to the database to do the translation search.",
        )
        parser.add_argument(
            "--max_workers",
            type=int,
            default=4,
            help="Number of max. workers to use to generate the configuration files for meTAline.",
        )
        parser.add_argument(
            "--joblist-size",
            type=int,
            default=2,
            help=(
                "Max. number of commands to be contained in"
                " a single joblist file. (Greasy array will "
                "contain multiple of these files to be launched parallel.)"
            ),
        )
        parser.add_argument(
            "--cluster-node-cores",
            type=int,
            default=40,
            help=("Total cores to assign for the cluster node."),
        )
        parser.add_argument(
            "--trimmo-cores",
            type=int,
            default=10,
            help=("Total cores to assign for trimmomatic."),
        )
        parser.add_argument(
            "--hisat2-cores",
            type=int,
            default=10,
            help=("Total cores to assign for hisat2."),
        )
        parser.add_argument(
            "--kraken2-cores",
            type=int,
            default=10,
            help=("Total cores to assign for kraken2."),
        )
        parser.add_argument(
            "--snakemake-cores",
            type=int,
            default=5,
            help=("Total cores to assign for snakemake."),
        )
        return PrepareGreasyArrayJobArgs(**vars(parser.parse_args(args)))


class PreparedConfigGenCmd(TypedDict):
    config_output_file: Path
    config_gen_cmd: str


class WriteJoblistFileArgs(TypedDict):
    idx: int
    basedir: Path
    cmds_split: list[str]


def get_common_prefixes_for_fasta_pairs(
    paths: tuple[Path, ...], fr_annotation_prefixes: Iterable
) -> set[str]:
    longest_unique_prefixes: set[str] = set()
    path_names: list[str] = [path.name for path in paths]
    combined_path_names: str = "\t".join(path_names)
    for path_name in path_names:
        previous_substring: str = ""
        longest_substring: str = ""
        for char in path_name:
            previous_substring = longest_substring
            longest_substring += char
            subcount = combined_path_names.count(longest_substring)
            if (
                subcount < 2
                and previous_substring != ""
                and previous_substring.endswith(tuple(fr_annotation_prefixes))
            ):
                longest_unique_prefixes.add(
                    previous_substring.rstrip("".join(fr_annotation_prefixes))
                )
                break
    return longest_unique_prefixes


def extract_reads_prefixes(reads_directory: Path | str, fr_annotation_prefixes: Iterable) -> set[str]:
    reads_dir = Path(reads_directory)
    files = tuple(reads_dir.rglob("*.f*q.gz"))
    found_prefixes = get_common_prefixes_for_fasta_pairs(files, fr_annotation_prefixes)
    return found_prefixes


def prepare_config_generation_commands(
    args: PrepareGreasyArrayJobArgs,
    reads_prefixes: Iterable[str],
    config_file_output_dir: Path,
) -> list[PreparedConfigGenCmd]:
    config_gen_cmds: list[PreparedConfigGenCmd] = []
    for prefix in reads_prefixes:
        config_output_file = config_file_output_dir / f"config.{prefix}.json"
        intermediate_files_output_dir = Path(args.basedir) / prefix

        metaphlan_db_abs_path = Path(args.metaphlan_db).absolute()
        n_db_abs_path = Path(args.n_db).absolute()
        protein_db_abs_path = Path(args.protein_db).absolute()
        alignment_out_abs_path = (intermediate_files_output_dir / "BAM").absolute()
        trimmomatic_out_abs_path = (
            intermediate_files_output_dir / "TRIMMOMATIC"
        ).absolute()
        kraken_out_abs_path = (
            intermediate_files_output_dir / "KRAKEN_ASSIGN"
        ).absolute()
        krona_out_abs_path = (intermediate_files_output_dir / "KRONA_HTML").absolute()
        extracted_fa_out_abs_path = (
            intermediate_files_output_dir / "EXTRACTED_FASTA"
        ).absolute()
        ranalysis_out_abs_path = (
            intermediate_files_output_dir / "R_ANALYSIS"
        ).absolute()
        metaphlan4_out_abs_path = (
            intermediate_files_output_dir / "METAPHLAN4"
        ).absolute()
        fr_annotation_prefixes_formatted = " ".join([f'"{prefix}"' for prefix in args.prefix_fr])

        cmd: str = " ".join(
            (
                f"{args.generate_config_cmd}",
                f"--config-file {config_output_file}",
                f"--extension {args.fastq_extension}",
                f"--basedir {args.basedir}",
                f"--prefix-fr {fr_annotation_prefixes_formatted}",
                f"--reads-directory {args.reads_directory}",
                f"--reference-genome {args.reference_genome}",
                f"--krakendb {args.krakendb}",
                f"--sample-barcode {prefix}",
                f"--fastq-prefix {prefix}",
                f"--metaphlan-db {metaphlan_db_abs_path}",
                f"--metaphlan-index {args.metaphlan_index}",
                f"--n-db {n_db_abs_path}",
                f"--protein-db {protein_db_abs_path}",
                f"--alignment-out {alignment_out_abs_path}",
                f"--trimmomatic-out {trimmomatic_out_abs_path}",
                f"--kraken-out {kraken_out_abs_path}",
                f"--krona-out {krona_out_abs_path}",
                f"--extracted-fa-out {extracted_fa_out_abs_path}",
                f"--ranalysis-out {ranalysis_out_abs_path}",
                f"--metaphlan4-out {metaphlan4_out_abs_path}",
                f"--trimmo-cores {args.trimmo_cores}",
                f"--hisat2-cores {args.hisat2_cores}",
                f"--kraken2-cores {args.kraken2_cores}",
            )
        )
        config_gen_cmds.append(
            {
                "config_output_file": config_output_file,
                "config_gen_cmd": cmd,
            }
        )
    return config_gen_cmds


def _exec_shell_cmd(cmd: str) -> str:
    exit_code, output = getstatusoutput(cmd)
    if exit_code != 0:
        raise ChildProcessError(output)
    return output


def _run_shell_cmds_parallel(cmds: list[str], max_workers: int = 4) -> tuple[str, ...]:
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        results = tuple(executor.map(_exec_shell_cmd, cmds))
    return results


def _write_joblist_file(args: WriteJoblistFileArgs) -> Path:
    joblist_output_file = args["basedir"] / f"joblist{args['idx']}.txt"
    with open(joblist_output_file, "w", encoding="UTF-8") as file:
        file.writelines(tuple(map(lambda line: f"{line}\n", args["cmds_split"])))
    return joblist_output_file


def write_joblist_files(
    max_workers: int, metaline_cmds_splits: list[list[str]], basedir: Path
) -> tuple[Path, ...]:
    write_joblist_args: list[WriteJoblistFileArgs] = [
        {"idx": idx, "cmds_split": cmds_split, "basedir": basedir}
        for idx, cmds_split in enumerate(metaline_cmds_splits)
    ]
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        res = tuple(executor.map(_write_joblist_file, write_joblist_args))
        return res


def write_greasy_list(
    basedir: Path,
    joblist_output_files: tuple[Path, ...],
    greasy_cmd: str,
) -> Path:
    output_greasy_list_file = basedir / "list_greasy.txt"
    greasy_list_lines = [f"{greasy_cmd} {joblist}" for joblist in joblist_output_files]
    with open(output_greasy_list_file, "w", encoding="UTF-8") as file:
        file.writelines(tuple(map(lambda line: f"{line}\n", greasy_list_lines)))
    return output_greasy_list_file


def write_slurm_job_file_template(
    basedir: Path,
    greasy_list_path,
    joblist_num: int,
    jobslist_max_size: int,
    cpus_per_task: int,
) -> Path:
    output_job_file = basedir / "launch_metaline.job"
    log_output_directory = str(basedir).rstrip("/")
    array_range = f"1-{joblist_num}"
    if joblist_num == 1:
        array_range = "1"
    template = f"""\
#!/bin/bash

#SBATCH --job-name=metaline
#SBATCH --qos=gp_bscls
#SBATCH --account=bsc40
#SBATCH --output={log_output_directory}/jobnm_%j.out # Output file, where
#SBATCH --error={log_output_directory}/jobnm_%j.err # File where the error is written

#SBATCH --array={array_range} # The number of jobs in the array
#SBATCH --ntasks={joblist_num} # The number of parallel tasks
#SBATCH --tasks-per-node={jobslist_max_size}
#SBATCH --cpus-per-task={cpus_per_task} # Number of CPUs per run task
#SBATCH --constraint=highmem  # Use this only if your the database is too big and this constrainst is available in your cluster.
#SBATCH --time=24:00:00  # 1 day.

module load greasy
# Print the task id
$(sed -n "${{SLURM_ARRAY_TASK_ID}}p" {greasy_list_path})
"""
    output_job_file.write_text(template)
    return output_job_file


def main() -> int:
    args = PrepareGreasyArrayJobArgs.get_arguments()
    basedir = Path(args.basedir)
    config_file_output_dir = basedir / "configs"
    shutil.rmtree(str(config_file_output_dir), ignore_errors=True)
    config_file_output_dir.mkdir(parents=True, exist_ok=True)
    reads_prefixes = extract_reads_prefixes(args.reads_directory, args.prefix_fr)
    gen_config_cmds = prepare_config_generation_commands(
        args,
        reads_prefixes,
        config_file_output_dir,
    )
    cmd_entries = [entry["config_gen_cmd"] for entry in gen_config_cmds]
    _ = _run_shell_cmds_parallel(cmd_entries, args.max_workers)
    generated_config_files = [entry["config_output_file"] for entry in gen_config_cmds]
    metaline_cmds = [
        f"{args.metaline_cmd} {args.metaline_rule} -j {args.snakemake_cores} --configfile {file}"
        for file in generated_config_files
    ]
    metaline_cmds_splits = [
        metaline_cmds[i : i + args.joblist_size]
        for i in range(0, len(metaline_cmds), args.joblist_size)
    ]
    joblist_output_files = write_joblist_files(
        max_workers=args.max_workers,
        basedir=basedir,
        metaline_cmds_splits=metaline_cmds_splits,
    )
    greasy_list = write_greasy_list(
        basedir=basedir,
        greasy_cmd=args.greasy_cmd,
        joblist_output_files=joblist_output_files,
    )
    write_slurm_job_file_template(
        basedir=basedir,
        greasy_list_path=greasy_list,
        joblist_num=len(joblist_output_files),
        jobslist_max_size=args.joblist_size,
        cpus_per_task=args.cluster_node_cores,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
