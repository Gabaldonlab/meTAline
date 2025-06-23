#!/usr/bin/env python3

from __future__ import annotations

import json
import argparse
import sys
from dataclasses import dataclass
from io import StringIO
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from subprocess import getstatusoutput

@dataclass
class GenerateHtmlReport:
    output_html_path: Path
    config_file: Path

    @classmethod
    def get_arguments(cls, args=sys.argv[1:]) -> GenerateHtmlReport:
        parser = argparse.ArgumentParser(
            description="Generates benchmark report in html format based on the latest 'benchmark' outputs of an already finished job."
        )
        parser.add_argument(
            "--output-html-path",
            type=Path,
            default=Path("./report.html"),
            help="Output path of the html report file.",
        )
        parser.add_argument(
            "--config-file",
            type=Path,
            required=True,
            help="Path to the config file of the job to create the report",
        )
        return GenerateHtmlReport(**vars(parser.parse_args(args)))

def _exec_shell_cmd(cmd: str) -> str:
    exit_code, output = getstatusoutput(cmd)
    if exit_code != 0:
        raise ChildProcessError(f"->cmd: {cmd}\n->stderr: {output}")
    return output


def main(argv: list[str] = sys.argv[1:]) -> int:
    args = GenerateHtmlReport.get_arguments(argv)
    config_blob: dict[Any, Any] = json.loads(args.config_file.read_text())

    base_dir: Path = Path(config_blob["Parameters"]["basedir"])
    if not base_dir.exists():
        raise FileNotFoundError(f"{base_dir} does not exists!")

    benchmark_dir = base_dir / "Benchmark"
    if not benchmark_dir.exists():
        raise FileNotFoundError(f"{benchmark_dir} does not exists!")


    output_dir = benchmark_dir / "plots"
    output_dir.mkdir(parents=True, exist_ok=True)

    benchmark_files = tuple(benchmark_dir.glob("*.benchmark.txt"))
    benchmark_lines: list[str] = []
    for idx, file in enumerate(benchmark_files):
        file_content = file.read_text().strip().splitlines()
        rule_name = file.name.replace(".benchmark.txt", "").split(".")[-1]
        if idx == 0:
            benchmark_lines.append(f"rule_name\t{file_content[0].strip()}")
        benchmark_lines.append(f"{rule_name}\t{file_content[1].strip()}")
    benchmark_data = "\n".join(benchmark_lines)

    column_renames = {
        "s": "Runtime (s)",
        "h:m:s": "Runtime (h:m:s)",
        "max_rss": "Max Memory (RSS, MB)",
        "max_vms": "Max Virtual Memory (VMS, MB)",
        "max_uss": "Max Unique Memory (USS, MB)",
        "max_pss": "Max Proportional Memory (PSS, MB)",
        "io_in": "I/O Read (MB)",
        "io_out": "I/O Write (MB)",
        "mean_load": "Mean CPU Load",
        "cpu_time": "CPU Time (s)",
    }

    df = pd.read_csv(StringIO(benchmark_data), sep="\t")
    df.rename(columns=column_renames, inplace=True)
    df.set_index("rule_name", inplace=True)

    for column in df.columns:
        if column in ["rule_name", "Runtime (h:m:s)"]:
            continue

        plt.figure(figsize=(10, 6))
        ax = sns.barplot(
            data=df,
            x="rule_name",
            y=column,
            hue="rule_name",
            palette="viridis",
            legend=False,
        )

        # Add value labels
        for container in ax.containers:
            ax.bar_label(container, fmt="%.2f", label_type="edge", fontsize=9)  # type: ignore

        plt.title(column, fontsize=14, fontweight="bold")
        plt.xlabel("Rule Name", fontsize=12)
        plt.ylabel(column, fontsize=12)
        plt.xticks(rotation=45, ha="right")
        plt.grid(axis="y", linestyle="--", alpha=0.5)
        plt.tight_layout()
        output_plot_name = (
            column.lower()
            .replace(" ", "_")
            .replace("(", "")
            .replace(")", "")
            .replace("/", "")
            .replace(",", "")
            + ".svg"
        )

        plot_output_path: Path = output_dir / output_plot_name
        plt.savefig(plot_output_path, format="svg")
        plt.close()
        print(f"-> Created '{column}' plot on path {plot_output_path}")

    _exec_shell_cmd(f"metaline --configfile {args.config_file} --report {args.output_html_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
