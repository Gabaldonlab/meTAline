#!/usr/bin/env python3

from __future__ import annotations

import sys
import json
import argparse
from pathlib import Path
from io import StringIO
from dataclasses import dataclass
from typing import Any
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

@dataclass
class GenerateHtmlReport:
    job_config: Path

    @classmethod
    def get_arguments(cls, args=sys.argv[1:]) -> GenerateHtmlReport:
        parser = argparse.ArgumentParser(
            description="Generates benchmark report in html format based on the 'benchmark' outputs of an already finished job."
        )
        parser.add_argument(
            "--job-config",
            type=Path,
            required=True,
            help="Path to the configuration file of the job.",
        )
        return GenerateHtmlReport(**vars(parser.parse_args(args)))


def main(argv: list[str] = sys.argv[1:]) -> int:
    args = GenerateHtmlReport.get_arguments(argv)
    config_content = args.job_config.read_text()
    config_payload: dict[Any, Any] = json.loads(config_content)
    job_basedir = Path(config_payload["Parameters"]["basedir"])
    benchmark_dir = job_basedir / "Benchmark"
    benchmark_files = tuple(benchmark_dir.glob("*_test_run.*.benchmark.txt"))
    run_date_ids = [file.name.split("_", maxsplit=1)[0] for file in benchmark_files]
    latest_date_id = max(set(run_date_ids))
    latest_benchmark_files = [file for file in benchmark_files if file.name.startswith(latest_date_id)]
    benchmark_lines: list[str] = []
    for idx, file in enumerate(latest_benchmark_files):
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
        "cpu_time": "CPU Time (s)"
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
            legend=False
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
        output_plot_name = column.lower().replace(" ", "_").replace("(", "").replace(")", "").replace("/", "") + ".png"
        plot_output_path = benchmark_dir.absolute() / output_plot_name
        plt.savefig(plot_output_path)
        plt.close()
        print(f"-> Created '{column}' plot on path {plot_output_path}")
    return 0

if __name__ == "__main__":
    raise SystemExit(main(["--job-config", "./updated_test_output/test_run.json"]))
