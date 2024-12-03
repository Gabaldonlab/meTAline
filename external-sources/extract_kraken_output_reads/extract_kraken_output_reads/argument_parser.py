#!/usr/bin/env python3

from __future__ import annotations
import argparse
from dataclasses import dataclass
import sys


@dataclass
class ExtractKrakenOutputReadsArgs:
    kraken_file: str
    seq_file1: str
    seq_file2: str
    taxid: list[str]
    output_file: str
    output_file2: str
    append: bool
    max_reads: int
    report_file: str
    parents: bool
    children: bool
    exclude: bool
    fastq_out: bool

    @classmethod
    def get_arguments(
        cls, args: list[str] = sys.argv[1:]
    ) -> ExtractKrakenOutputReadsArgs:
        parser = argparse.ArgumentParser(
            description=(
                "Extracts reads classified by Kraken as a specified taxonomy ID."
                " Those reads are extracted into a new FASTA file."
            )
        )
        parser.add_argument(
            "-k",
            dest="kraken_file",
            required=True,
            help="Kraken output file to parse",
        )
        parser.add_argument(
            "-s",
            "-s1",
            "-1",
            "-U",
            dest="seq_file1",
            required=True,
            help="FASTA/FASTQ File containing the raw sequence letters.",
        )
        parser.add_argument(
            "-s2",
            "-2",
            dest="seq_file2",
            default="",
            help="2nd FASTA/FASTQ File containing the raw sequence letters (paired).",
        )
        parser.add_argument(
            "-t",
            "--taxid",
            dest="taxid",
            required=True,
            nargs="+",
            help="Taxonomy ID[s] of reads to extract (space-delimited)",
        )
        parser.add_argument(
            "-o",
            "--output",
            dest="output_file",
            required=True,
            help="Output FASTA/Q file containing the reads and sample IDs",
        )
        parser.add_argument(
            "-o2",
            "--output2",
            dest="output_file2",
            required=False,
            default="",
            help="Output FASTA/Q file containig the second pair of reads [required for paired input]",
        )
        parser.add_argument(
            "--append",
            dest="append",
            action="store_true",
            help="Append the sequences to the end of the output FASTA file specified.",
        )
        parser.add_argument(
            "--noappend",
            dest="append",
            action="store_false",
            help="Create a new FASTA file containing sample sequences and IDs \
                (rewrite if existing) [default].",
        )
        parser.add_argument(
            "--max",
            dest="max_reads",
            required=False,
            default=100000000,
            type=int,
            help="Maximum number of reads to save [default: 100,000,000]",
        )
        parser.add_argument(
            "-r",
            "--report",
            dest="report_file",
            required=False,
            default="",
            help="Kraken report file. [required only if --include-parents/children \
            is specified]",
        )
        parser.add_argument(
            "--include-parents",
            dest="parents",
            required=False,
            action="store_true",
            default=False,
            help="Include reads classified at parent levels of the specified taxids",
        )
        parser.add_argument(
            "--include-children",
            dest="children",
            required=False,
            action="store_true",
            default=False,
            help="Include reads classified more specifically than the specified taxids",
        )
        parser.add_argument(
            "--exclude",
            dest="exclude",
            required=False,
            action="store_true",
            default=False,
            help="Instead of finding reads matching specified taxids, finds all reads NOT matching specified taxids",
        )
        parser.add_argument(
            "--fastq-output",
            dest="fastq_out",
            required=False,
            action="store_true",
            default=False,
            help="Print output FASTQ reads [requires input FASTQ, default: output is FASTA]",
        )
        parser.set_defaults(append=False)

        return ExtractKrakenOutputReadsArgs(**vars(parser.parse_args(args)))
