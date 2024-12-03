#!/usr/bin/env python
import argparse
import gzip
from functools import wraps
from time import gmtime
from time import strftime
from typing import Any
from typing import IO
from typing import Literal
from typing import TypedDict

from Bio import SeqIO  # type: ignore
from extract_kraken_output_reads.argument_parser import (
    ExtractKrakenOutputReadsArgs,
)


class Tree(object):
    """
    usage: Tree node used in constructing taxonomy tree.
           Includes only taxonomy levels and genomes identified in the Kraken report.
    """

    def __init__(self, taxid, level_num, level_id, children=None, parent=None):
        self.taxid = taxid
        self.level_num = level_num
        self.level_id = level_id
        self.children = []
        self.parent = parent
        if children is not None:
            for child in children:
                self.add_child(child)

    def add_child(self, node):
        if not isinstance(node, Tree):
            raise ValueError(
                "Given node is not an instance of the Tree class."
            )
        self.children.append(node)


class EvaulatedKrakenSamplesResult(TypedDict):
    count_kraken: int
    saved_read_ids_1: set[str]
    saved_read_ids_2: set[str]


def timing(f):
    @wraps(f)
    def wrap(*args, **kw):
        # Start Program
        start_time = strftime("%m-%d-%Y %H:%M:%S", gmtime())
        print(f"PROGRAM START TIME: {start_time}")
        result = f(*args, **kw)
        # End of program
        end_time = strftime("%m-%d-%Y %H:%M:%S", gmtime())
        print(f"PROGRAM END TIME: {end_time}")
        return result

    return wrap


def evaluate_kraken_file_samples(
    args: ExtractKrakenOutputReadsArgs,
    save_taxids: dict[int, int],
    exclude_taxids: dict[int, int],
) -> EvaulatedKrakenSamplesResult:
    """
    Process Kraken file for classified read IDs.
    Evaluate each sample in the Kraken file.

    return:
        - count_kraken: int
        - saved_read_ids_1: set[str]
        - saved_read_ids_2: set[str]

    """
    opened_kraken_file: IO[Any] = open(args.kraken_file, "r", encoding="UTF-8")
    result: EvaulatedKrakenSamplesResult = {
        "count_kraken": 0,
        "saved_read_ids_1": set(),
        "saved_read_ids_2": set(),
    }
    for line in opened_kraken_file:
        result["count_kraken"] += 1
        if result["count_kraken"] % 10000 == 0:
            print(
                f"\t{result['count_kraken'] / 1000000.0:.2f} million reads processed"
            )

        # Parse line for results
        [tax_id, read_id] = process_kraken_output(line)
        if tax_id == -1:
            continue

        # Skip if reads are human/artificial/synthetic
        if (tax_id in save_taxids) and not args.exclude:
            save_taxids[tax_id] += 1
            result["saved_read_ids_1"].add(read_id)
            result["saved_read_ids_2"].add(read_id)
        elif (tax_id not in exclude_taxids) and args.exclude:
            save_taxids.setdefault(tax_id, 1)
            if tax_id in save_taxids:
                save_taxids[tax_id] += 1

            result["saved_read_ids_1"].add(read_id)
            result["saved_read_ids_2"].add(read_id)

        if len(result["saved_read_ids_1"]) >= args.max_reads:
            break

    opened_kraken_file.close()
    return result


def process_kraken_output(kraken_line: str) -> tuple[int, str]:
    """
    usage: Parses single line from kraken output and returns taxonomy ID and readID.
    input: Kraken output file with readid and taxid in the second and third tab-delimited columns.
    returns:
        - taxonomy ID
        - read ID
    """
    l_vals = kraken_line.split("\t")
    if len(l_vals) < 5:
        return (-1, "")
    if "taxid" in l_vals[2]:
        temp = l_vals[2].split("taxid ")[-1]
        tax_id = temp[:-1]
    else:
        tax_id = l_vals[2]

    read_id = l_vals[1]
    default_tax_id: int = 81077
    if tax_id != "A":
        default_tax_id = int(tax_id)
    return (default_tax_id, read_id)


def process_kraken_report(report_line: str) -> tuple[int, int, str]:
    """
    usage: Parses single line from report output and returns taxID, levelID.
    input: Kraken report file with the following tab delimited lines:
        - percent of total reads
        - number of reads (including at lower levels)
        - number of reads (only at this level)
        - taxonomy classification of level
            (U, - (root), - (cellular org), D, P, C, O, F, G, S)
        - taxonomy ID (0 = unclassified, 1 = root, 2 = Bacteria...etc)
        - spaces + name
    returns:
    - taxonomy ID
    - level number (number of spaces before name)
    - level_type (type of taxonomy level - U, R, D, P, C, O, F, G, S, etc)
    """
    l_vals: list[str] = report_line.strip().split("\t")
    if len(l_vals) < 5:
        raise ValueError(
            "Not enough fields in report line. Expected more than 5."
        )
    int(l_vals[1])

    # Extract relevant information
    try:
        taxid = int(l_vals[-3])
        level_type = l_vals[-2]
        map_kuniq = {
            "species": "S",
            "genus": "G",
            "family": "F",
            "order": "O",
            "class": "C",
            "phylum": "P",
            "superkingdom": "D",
            "kingdom": "K",
        }
        if level_type not in map_kuniq:
            level_type = "-"
        else:
            level_type = map_kuniq[level_type]
    except ValueError:
        taxid = int(l_vals[-2])
        level_type = l_vals[-3]
    # Get spaces to determine level num
    spaces = 0
    for char in l_vals[-1]:
        if char == " ":
            spaces += 1
        else:
            break
    level_num = int(spaces / 2)
    return (taxid, level_num, level_type)


def get_file_type(first_line: str) -> Literal["fasta", "fastq"]:
    """
    Determine the file type based on the first character of the first line.

    return: "fasta" or "fastq"
    """
    if first_line == ">":
        return "fasta"
    elif first_line == "@":
        return "fastq"
    raise ValueError("ERROR: sequence file must be FASTA or FASTQ.")


def get_output_file_io_wrappers(
    args: ExtractKrakenOutputReadsArgs,
) -> tuple[IO[Any], IO[Any] | None]:
    """
    Returns the file IO wrapper objects for the input files.

    return:
        (output_file1: TextIO, output_file2: TextIO)
    """
    output_files_io_mode: str = "w"
    if args.append:
        output_files_io_mode = "a"

    output_file1: IO[Any] = open(args.output_file, output_files_io_mode)
    output_file2: IO[Any] | None = None
    if args.output_file2:
        output_file2 = open(args.output_file2, output_files_io_mode)
    return output_file1, output_file2


def process_sequence_file(
    args: Any,
    sequence_file: IO[Any],
    save_readids: set[str],
    output_file: IO[Any],
    filetype: Literal["fastq", "fasta"],
) -> tuple[int, int]:
    """
    Parses the first input sequence file.
    return:
        (sequence_count: int, output_count: int)
    """
    # Process SEQUENCE 1 file
    seqs_count: int = 0
    output_count: int = 0
    for record in SeqIO.parse(sequence_file, filetype):
        seqs_count += 1
        # Print update
        if seqs_count % 1000 == 0:
            print(
                f"\t{output_count} read IDs found ({seqs_count / 1000000.0:.2f} mill reads processed)"
            )

        test_id, test_id2 = get_test_ids(record.id)
        if test_id in save_readids or test_id2 in save_readids:
            output_count += 1

            # Print update
            print(
                f"\t{output_count} read IDs found ({seqs_count / 1000000.0:.2f} mill reads processed)"
            )
            save_seqio_record_to_file(args, record, output_file)

        # If no more reads to find
        if len(save_readids) == output_count:
            break

    # Close files
    sequence_file.close()
    output_file.close()

    return seqs_count, output_count


def get_test_ids(record_id: str | int) -> tuple[str, str]:
    """
    Reformats the test ID from a Bio.SeqIO based record ID.
    """
    test_id: str = str(record_id)
    test_id2: str = test_id
    if "/1" in test_id or "/2" in test_id:
        test_id2 = test_id[:-2]
    return test_id, test_id2


def save_seqio_record_to_file(args, record: Any, output_file: IO[Any]) -> int:
    if args.fastq_out:
        return SeqIO.write(record, output_file, "fastq")

    return SeqIO.write(record, output_file, "fasta")


def get_input_sequence_file_io(input_sequence_file_path: str) -> IO[Any]:
    """
    Returns the IOWrapper object for the given file path based on,
    if it's gz compressed or regular text file.
    """
    if input_sequence_file_path.endswith(".gz"):
        return gzip.open(input_sequence_file_path, "rt")
    return open(input_sequence_file_path, "rt")


def check_input_file_type(
    args: ExtractKrakenOutputReadsArgs, input_sequence_file: str
) -> Literal["fasta", "fastq"]:
    """
    Checks and assures that the given input sequence file is fasta or fastq.
    return:
        - file type (fasta / fastq)
    """
    with get_input_sequence_file_io(input_sequence_file) as s_file1:
        first_line: str = s_file1.readline().strip()
        if not first_line:
            raise ValueError("ERROR: sequence file's first line is blank.")

        file_type: Literal["fasta", "fastq"] = get_file_type(first_line[0])
        if file_type != "fastq" and args.fastq_out:
            raise argparse.ArgumentTypeError(
                "ERROR: for FASTQ output, input file must be FASTQ"
            )
        return file_type


@timing
def main() -> int:
    args = ExtractKrakenOutputReadsArgs.get_arguments()

    # Check input
    if (len(args.output_file2) == 0) and (len(args.seq_file2) > 0):
        raise ValueError(
            "Must specify second output file -o2 for paired input"
        )

    # Initialize taxids
    save_taxids: dict[int, int] = {int(tid): 0 for tid in args.taxid}
    main_lvls: tuple[str, ...] = ("R", "K", "D", "P", "C", "O", "F", "G", "S")

    # STEP 0: READ IN REPORT FILE AND GET ALL TAXIDS
    if args.parents or args.children:
        # check that report file exists
        if not args.report_file.strip():
            raise argparse.ArgumentTypeError(
                ">> ERROR: --report not specified."
            )
        print(f">> STEP 0: PARSING REPORT FILE {args.report_file}")

        # create tree and save nodes with taxids in the list
        base_nodes: dict[int, Tree] = {}
        r_file = open(args.report_file, "r")
        prev_node: Tree | None = None
        for line in r_file:
            # extract values
            report_vals = process_kraken_report(line)
            if len(report_vals) == 0:
                continue
            (taxid, level_num, level_id) = report_vals
            if taxid == 0:
                continue
            # tree root
            if taxid == 1:
                level_id = "R"
                root_node = Tree(taxid, level_num, level_id)
                prev_node = root_node
                # save if needed
                if taxid in save_taxids:
                    base_nodes[taxid] = root_node
                continue

            if prev_node is None:
                raise ValueError("Could not assign the previous node in Tree.")

            # move to correct parent
            while prev_node is not None and level_num != (
                prev_node.level_num + 1
            ):
                prev_node = prev_node.parent

            if prev_node is None:
                raise ValueError("Could not assign the previous node in Tree.")

            # determine correct level ID
            if prev_node and (level_id == "-" or len(level_id) > 1):
                if prev_node.level_id in main_lvls:
                    level_id = prev_node.level_id + "1"
                else:
                    num = int(prev_node.level_id[-1]) + 1
                    level_id = prev_node.level_id[:-1] + str(num)

            # make node
            curr_node = Tree(taxid, level_num, level_id, None, prev_node)
            prev_node.add_child(curr_node)
            prev_node = curr_node

            # save if taxid matches
            if taxid in save_taxids:
                base_nodes[taxid] = curr_node

        r_file.close()

        # FOR SAVING PARENTS
        if args.parents:
            # For each node saved, traverse up the tree and save each taxid
            for tid in base_nodes:
                curr_node = base_nodes[tid]
                while curr_node.parent != None:
                    curr_node = curr_node.parent
                    save_taxids[curr_node.taxid] = 0

        # FOR SAVING CHILDREN
        if args.children:
            for tid in base_nodes:
                curr_nodes = base_nodes[tid].children
                while len(curr_nodes) > 0:
                    # For this node
                    curr_n = curr_nodes.pop()
                    if curr_n.taxid not in save_taxids:
                        save_taxids[curr_n.taxid] = 0
                    # Add all children
                    if curr_n.children != None:
                        for child in curr_n.children:
                            curr_nodes.append(child)

    ##############################################################################
    print(f"\t{len(save_taxids)} taxonomy IDs to parse")
    print(f">> STEP 1: PARSING KRAKEN FILE FOR READIDS {args.kraken_file}")

    exclude_taxids = {}
    if args.exclude:
        exclude_taxids = save_taxids
        save_taxids = {}

    kraken_samples_result: EvaulatedKrakenSamplesResult = (
        evaluate_kraken_file_samples(args, save_taxids, exclude_taxids)
    )

    print(
        f"\t{kraken_samples_result['count_kraken'] / 1000000.0:.2f} million reads processed\n"
    )
    print(f"\t{len(kraken_samples_result['saved_read_ids_1'])} read IDs saved")

    ####TEST IF INPUT IS FASTA OR FASTQ

    file_type: Literal["fasta", "fastq"] = check_input_file_type(
        args, args.seq_file1
    )
    s_file1: IO[Any] = get_input_sequence_file_io(args.seq_file1)
    s_file2: IO[Any] = get_input_sequence_file_io(args.seq_file1)

    # PROCESS INPUT FILE AND SAVE FASTA FILE
    print(">> STEP 2: READING SEQUENCE FILES AND WRITING READS")
    print("\t0 read IDs found (0 mill reads processed)")

    # Open output file
    output_file1, output_file2 = get_output_file_io_wrappers(args)
    seq_file1_seq_count, seq_file1_output_count = process_sequence_file(
        args,
        s_file1,
        kraken_samples_result["saved_read_ids_1"],
        output_file1,
        file_type,
    )

    print(
        f"\r\t{seq_file1_output_count} read IDs found"
        f" ({seq_file1_seq_count / 1000000.0:.2f} mill reads processed)\n"
    )

    if args.seq_file2 and output_file2 is not None:
        seq_file2_seq_count, seq_file2_output_count = process_sequence_file(
            args,
            s_file2,
            kraken_samples_result["saved_read_ids_2"],
            output_file2,
            file_type,
        )
        print(
            f"\r\t{seq_file2_output_count} read IDs found"
            f" ({seq_file2_seq_count / 1000000.0:.2f} mill reads processed)"
        )

    total_output_count: int = seq_file1_output_count + seq_file2_output_count
    print(f"\t{total_output_count} reads printed to file")
    print(f"\tGenerated file: {args.output_file}.")
    if args.output_file2:
        print("\tGenerated file: {args.output_file2}.")

    return 0
