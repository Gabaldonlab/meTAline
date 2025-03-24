#!/usr/bin/env python3

"""
Author: Diego Fuentes and Olfat Khannous
Contact email: olfat.khannous@bsc.es
Date:2024-07-23
"""

import os
import json
import argparse
import sys
import re


class CreateConfigurationFile:
    """Class which manages Configuration file Manager"""

    def __init__(self):
        """Class constructor of the pipeline"""

        # GENERAL PARAMETERS
        self.configFile = (
            "config.json"  # Name of the json configuration file to be created.
        )
        self.version = 1  # Pipeline version
        self.logs_dir = "WGS_logs"  # Directory to keep all the log files
        self.extension = "fastq.gz"  # Extension of the raw read files
        self.sample_barcode = None  # Sample barcode
        self.basedir = self.sample_barcode  # Base directory for the pipeline run

        # INPUT PARAMETERS
        self.reads_directory = None  # Directory where the illumina fastqs are stored
        self.reference_genome = None  # Indexed reference genome (Hisat2 index) path
        self.krakendb = None  # Kraken2 DB path that is going to be used
        self.taxid = 0  # Default taxid used to extract reads, based on the kraken report. By default 0 -> unclassified
        self.metaphlan_db = None  # Metaphlan4 database path
        self.metaphlan_Index = None  # Metaphlan4 index (last release of the database)
        self.n_db = None  # Humann database for nucleotide search
        self.protein_db = None  # Humann database for translation search

        # OUTPUT PARAMETERS
        self.alignment_out = "BAM"  # Out directory of the alignment step
        self.trimmomatic_out = "TRIMMOMATIC"  # Out directory of the sv calls
        self.kraken_out = "KRAKEN_ASSIGN"  # Out directory of the Kraken2 reports
        self.krona_out = "KRONA_HTML"  # Out directory of the Krona visualization tool
        self.extracted_fa_out = (
            "EXTRACTED_FASTA"  # Out directory of the extracted fasta
        )
        self.ranalysis_out = "R_ANALYSIS"  # Out directory of R analysis
        self.metaphlan4_out = "METAPHLAN4"  # Metaphlan4 and Humann directory

        # WILDCARD PARAMETER
        self.fastq_prefix = None  # List with basename of the fastq-prefix if any

        # TRIMMOMATIC PARAMETERS
        self.trimmo_cores = 4  # Number of threads to run the trimmomatic
        self.leading = 3  # leading parameter
        self.trailing = 3  # trailing parameter
        self.slidingwindow = "4:15"  # sliding window parameter (two values, ie: 4:15)
        self.minlen = 35  # min read lenght paramter
        self.illuminaclip = "Illumina_MGI_adapters.fa:2:30:10"  # Fasta file with the combination of both Illumina and MGI/BGI adapters

        # FASTQC parameters
        self.adapters = "adapter_list_new.txt"  # This is the extended list of sequences to be checked as adapters when making the quality assessment

        # HISAT PARAMETERS (set of default parameters for local-sensitive)
        self.hisat2_cores = 8  # Number of threads to run the ngmlr aligner

        # KRAKEN2 PARAMETERS
        self.kraken2_cores = 16

        # DICTIONARIES
        self.allParameters = {}
        self.generalParameters = {}
        self.inputParameters = {}
        self.outputParameters = {}
        self.wildcardParameters = {}
        self.TrimmomaticParameters = {}
        self.fastqcParameters = {}
        self.hisat2Parameters = {}
        self.Kraken2Parameters = {}

    def register_parameter(self, parser: argparse.ArgumentParser):
        """Register all parameters with the given
        argparse parser"""
        self.register_general(parser)
        self.register_input(parser)
        self.register_output(parser)
        self.register_wildcards(parser)
        self.register_trimmomatic(parser)
        self.register_fastqc(parser)
        self.register_hisat2(parser)
        self.register_kraken2(parser)

    def register_general(self, parser: argparse.ArgumentParser):
        """Register all general parameters with the given
        argparse parser
        parser -- the argparse parser
        """
        general_group = parser.add_argument_group("General Parameters")
        general_group.add_argument(
            "--configFile",
            dest="configFile",
            metavar="configFile",
            help="Configuration JSON to be generated.",
            default=self.configFile,
        )
        general_group.add_argument(
            "--version",
            type=int,
            dest="version",
            metavar="version",
            default=self.version,
            help="Pipeline run version.",
        )
        general_group.add_argument(
            "--logs-dir",
            dest="logs_dir",
            metavar="logs_dir",
            default=self.logs_dir,
            help="Directory to keep all the log files.",
        )
        general_group.add_argument(
            "--extension",
            dest="extension",
            metavar="extension",
            default=self.extension,
            help="Extension of the illumina raw read files.",
        )
        general_group.add_argument(
            "--sample-barcode",
            dest="sample_barcode",
            metavar="sample_barcode",
            default=self.sample_barcode,
            help="Sample barcode.",
        )
        general_group.add_argument(
            "--basedir",
            dest="basedir",
            metavar="basedir",
            default=self.basedir,
            help="Base directory for the pipeline run.",
        )

    def register_input(self, parser: argparse.ArgumentParser):
        """Register all input parameters with the given
        argparse parser
        parser -- the argparse parser
        """
        input_group = parser.add_argument_group("Inputs")
        input_group.add_argument(
            "--reads-directory",
            dest="reads_directory",
            metavar="reads_directory",
            default=self.reads_directory,
            help="Directory where the Illumina fastqs are stored.",
        )
        input_group.add_argument(
            "--reference-genome",
            dest="reference_genome",
            metavar="reference_genome",
            default=self.reference_genome,
            help="Indexed reference genome path (ie: reference/Human_index).",
        )
        input_group.add_argument(
            "--krakendb",
            dest="krakendb",
            required=True,
            metavar="krakendb",
            default=self.krakendb,
            help="Kraken2 database used for the taxonomic assignment.",
        )
        input_group.add_argument(
            "--taxid",
            dest="taxid",
            metavar="taxid",
            default=self.taxid,
            help="Selects the taxid used for exctracting fasta reads. You can add more than one using space as the separator.",
        )
        input_group.add_argument(
            "--metaphlan_db",
            dest="metaphlan_db",
            metavar="metaphlan_db",
            default=self.metaphlan_db,
            help="Path of the Metaphlan4 database.",
        )
        input_group.add_argument(
            "--metaphlan_Index",
            dest="metaphlan_Index",
            metavar="metaphlan_Index",
            default=self.metaphlan_Index,
            help="Index of the metaphlan4 database.",
        )
        input_group.add_argument(
            "--protein_db",
            dest="protein_db",
            metavar="protein_db",
            default=self.protein_db,
            help="Humann database to do the translation search (by default this is by-passed).",
        )
        input_group.add_argument(
            "--n_db",
            dest="n_db",
            metavar="n_db",
            default=self.n_db,
            help="Humann database to do the nucleotide search (based on already built annotations).",
        )

    def register_output(self, parser: argparse.ArgumentParser):
        """Register all output parameters with the given
        argparse parser
        parser -- the argparse parser
        """

        output_group = parser.add_argument_group("Outputs")
        output_group.add_argument(
            "--alignment-out",
            dest="alignment_out",
            default=self.alignment_out,
            help="Out directory of the alignment step.",
        )
        output_group.add_argument(
            "--trimmomatic-out",
            dest="trimmomatic_out",
            default=self.trimmomatic_out,
            help="Out directory of the trimmomatic output.",
        )
        output_group.add_argument(
            "--kraken-out",
            dest="kraken_out",
            default=self.kraken_out,
            help="Out directory of the Kraken2 taxonomic assignation step.",
        )
        output_group.add_argument(
            "--krona-out",
            dest="krona_out",
            default=self.krona_out,
            help="Out directory of the Krona visualization tool.",
        )
        output_group.add_argument(
            "--extracted-fa-out",
            dest="extracted_fa_out",
            default=self.extracted_fa_out,
            help="Out directory of the extracted fasta reads.",
        )
        output_group.add_argument(
            "--ranalysis-out",
            dest="ranalysis_out",
            default=self.ranalysis_out,
            help="Out directory of the R analysis.",
        )
        output_group.add_argument(
            "--metaphlan4-out",
            dest="metaphlan4_out",
            default=self.metaphlan4_out,
            help="Out directory of the Metaphlan4 and Humann analysis.",
        )

    def register_wildcards(self, parser: argparse.ArgumentParser):
        """Register all wildcards parameters with the given
        argparse parser
        parser -- the argparse parser
        """
        wildcards_group = parser.add_argument_group("Wildcards")
        wildcards_group.add_argument(
            "--fastq-prefix",
            dest="fastq_prefix",
            metavar="fastq_prefix",
            required=True,
            help="List with basename of the fastq-prefix.",
        )

    def register_trimmomatic(self, parser: argparse.ArgumentParser):
        """Register all trimmomatic parameters with the given
        argparse parser
        parser -- the argparse parser
        """
        trimmomatic_group = parser.add_argument_group("Trimmomatic parameters")
        trimmomatic_group.add_argument(
            "--trimmo-cores",
            type=int,
            dest="trimmo_cores",
            metavar="trimmo_cores",
            default=self.trimmo_cores,
            help="Number of threads to run trimmomatic.",
        )
        trimmomatic_group.add_argument(
            "--leading",
            type=int,
            dest="leading",
            metavar="leading",
            default=self.leading,
            help="Leading parameter.",
        )
        trimmomatic_group.add_argument(
            "--trailing",
            type=int,
            dest="trailing",
            metavar="trailing",
            default=self.trailing,
            help="Trailing parameter.",
        )
        trimmomatic_group.add_argument(
            "--slidingwindow",
            dest="slidingwindow",
            metavar="slidingwindow",
            default=self.slidingwindow,
            help="Sliding window parameter.",
        )
        trimmomatic_group.add_argument(
            "--minlen",
            type=int,
            dest="minlen",
            metavar="minlen",
            default=self.minlen,
            help="Minimum read length parameter.",
        )
        trimmomatic_group.add_argument(
            "--illuminaclip",
            dest="illuminaclip",
            metavar="illuminaclip",
            default=self.illuminaclip,
            help="Illumina clip information.",
        )

    def register_fastqc(self, parser: argparse.ArgumentParser):
        """Register all fastqc parameters with the given
        argparse parser
        parser -- the argparse parser
        """
        fastqc_group = parser.add_argument_group("fastqc parameters")
        fastqc_group.add_argument(
            "--adapters",
            dest="adapters",
            metavar="adapters",
            default=self.adapters,
            help="Adapter sequences to be checked in the quality assessment.",
        )

    def register_hisat2(self, parser: argparse.ArgumentParser):
        """Register all hisat2 aligner parameters with the given
        argparse parser
        parser -- the argparse parser
        """
        hisat2_group = parser.add_argument_group("hisat2 parameters")
        hisat2_group.add_argument(
            "--hisat2-cores",
            type=int,
            dest="hisat2_cores",
            metavar="hisat2_cores",
            default=self.hisat2_cores,
            help="Number of threads to run the hisat2 aligner.",
        )

    def register_kraken2(self, parser: argparse.ArgumentParser):
        """Register all kraken2 parameters with the given
        argparse parser
        parser -- the argparse parser
        """
        kraken2_group = parser.add_argument_group("Kraken2 parameters")
        kraken2_group.add_argument(
            "--kraken2-cores",
            type=int,
            dest="kraken2_cores",
            metavar="kraken2_cores",
            default=self.kraken2_cores,
            help="Number of threads to run the Kraken2 taxonomic assignment step.",
        )

    def check_parameters(self, args, parser: argparse.ArgumentParser):
        """Check parameters consistency

        args -- set of parsed arguments

        Note that by using os.path.exist it only checks if the file/directory exist but function may return False,
        if permission is not granted to execute os.stat() on the requested file, even if the path exists.
        """

        if args.configFile is None:
            parser.print_help()
            sys.exit(-1)

        working_dir = os.getcwd() + "/"

        if args.sample_barcode is None:
            print(
                "No sample_barcode specified. "
                "A barcode or identification is required",
                file=sys.stderr,
            )
            parser.print_help()
            sys.exit(1)

        if args.basedir:
            args.basedir = os.path.abspath(args.basedir) + "/"
        else:
            args.basedir = (
                f"{working_dir}v{args.version}/{args.sample_barcode}/"
            )

        if args.logs_dir:
            args.logs_dir = os.path.abspath(args.logs_dir) + "/"
        else:
            args.logs_dir = f"{args.basedir}{self.logs_dir}/"

        if args.reads_directory:
            args.reads_directory = os.path.abspath(args.reads_directory) + "/"
        else:
            args.reads_directory = (
                f"{working_dir}reads/Illumina/{args.sample_barcode}/"
            )
        if not os.path.exists(args.reads_directory):
            print(
                f"{args.reads_directory} not found. "
                "The directory where the reads are located is required. "
                "Exiting now.",
                file=sys.stderr,
            )
            parser.print_help()
            sys.exit(1)

        if args.reference_genome != None:
            args.reference_genome = os.path.abspath(args.reference_genome)

            if not os.path.exists(args.reference_genome):
                if not re.search("index", args.reference_genome):
                    print(
                        f"A reference genome index (including index as suffix, ie:'HUMAN_index') has been not provided or it has not been found in "
                        f"{args.reference_genome}. Enabling the taxonomic assignment of a environmental sample."
                    , file=sys.stderr)
                    args.reference_genome = None

        if args.krakendb:
            args.krakendb = os.path.abspath(args.krakendb)
        if not os.path.exists(args.krakendb):
            print(
                f"The Kraken2 DB has not been provided. Check the your provided path {args.krakendb}. Note that this step is mandatory"
            , file=sys.stderr)
            parser.print_help()
            sys.exit(1)

        if args.alignment_out:
            args.alignment_out = os.path.abspath(args.alignment_out) + "/"
        else:
            args.alignment_out = f"{args.basedir}{self.alignment_out}/"

        if args.trimmomatic_out:
            args.trimmomatic_out = os.path.abspath(args.trimmomatic_out) + "/"
        else:
            args.trimmomatic_out = f"{args.basedir}{self.trimmomatic_out}/"

        if args.kraken_out:
            args.kraken_out = os.path.abspath(args.kraken_out) + "/"
        else:
            args.kraken_out = f"{args.basedir}{self.kraken_out}/"

        if args.krona_out:
            args.krona_out = os.path.abspath(args.krona_out) + "/"
        else:
            args.krona_out = f"{args.basedir}{self.krona_out}/"

        if args.extracted_fa_out:
            args.extracted_fa_out = os.path.abspath(args.extracted_fa_out) + "/"
        else:
            args.extracted_fa_out = f"{args.basedir}{self.extracted_fa_out}/"

        if args.ranalysis_out:
            args.ranalysis_out = os.path.abspath(args.ranalysis_out) + "/"
        else:
            args.ranalysis_out = f"{args.basedir}{self.ranalysis_out}/"

        if args.metaphlan4_out:
            args.metaphlan4_out = os.path.abspath(args.metaphlan4_out) + "/"
        else:
            args.metaphlan4_out = f"{args.basedir}{self.metaphlan4_out}/"

    def store_general_parameters(self, args):
        """Updates general parameters to the map of parameters to be store in a JSON file
        args -- set of parsed arguments
        """
        self.generalParameters["configFile"] = args.configFile
        self.generalParameters["version"] = args.version
        self.generalParameters["basedir"] = args.basedir
        self.generalParameters["logs_dir"] = args.logs_dir
        self.generalParameters["extension"] = args.extension
        self.generalParameters["sample_barcode"] = args.sample_barcode
        self.allParameters["Parameters"] = self.generalParameters

    def store_input_parameters(self, args):
        """Updates input parameters to the map of parameters to be store in a JSON file
        args -- set of parsed arguments
        """

        self.inputParameters["reads_directory"] = args.reads_directory
        self.inputParameters["reference_genome"] = args.reference_genome
        self.inputParameters["krakendb"] = args.krakendb
        self.inputParameters["taxid"] = args.taxid
        self.inputParameters["metaphlan_db"] = args.metaphlan_db
        self.inputParameters["metaphlan_Index"] = args.metaphlan_Index
        self.inputParameters["n_db"] = args.n_db
        self.inputParameters["protein_db"] = args.protein_db
        self.allParameters["Inputs"] = self.inputParameters

    def store_output_parameters(self, args):
        """Updates output parameters to the map of parameters to be store in a JSON file
        args -- set of parsed arguments
        """
        self.outputParameters["alignment_out"] = args.alignment_out
        self.outputParameters["trimmomatic_out"] = args.trimmomatic_out
        self.outputParameters["kraken_out"] = args.kraken_out
        self.outputParameters["krona_out"] = args.krona_out
        self.outputParameters["extracted_fa_out"] = args.extracted_fa_out
        self.outputParameters["ranalysis_out"] = args.ranalysis_out
        self.outputParameters["metaphlan4_out"] = args.metaphlan4_out
        self.allParameters["Outputs"] = self.outputParameters

    def store_wildcard_parameters(self, args):
        """Updates wildcard parameters to the map of parameters to be store in a JSON file
        args -- set of parsed arguments
        """
        self.wildcardParameters["fastq_prefix"] = args.fastq_prefix
        self.allParameters["Wildcards"] = self.wildcardParameters

    def store_trimmomatic_parameters(self, args):
        """Updates the trimmomatic parameters to the map of parameters to be store in a JSON file
        args -- set of parsed arguments
        """
        self.TrimmomaticParameters["trimmo_cores"] = args.trimmo_cores
        self.TrimmomaticParameters["leading"] = args.leading
        self.TrimmomaticParameters["trailing"] = args.trailing
        self.TrimmomaticParameters["slidingwindow"] = args.slidingwindow
        self.TrimmomaticParameters["minlen"] = args.minlen
        self.TrimmomaticParameters["illuminaclip"] = args.illuminaclip
        self.allParameters["Trimmomatic"] = self.TrimmomaticParameters

    def store_fastqc_parameters(self, args):
        """Updates the Fastqc parameters to the map of parameters to be store in a JSON file
        args -- set of parsed arguments
        """
        self.fastqcParameters["adapters"] = args.adapters
        self.allParameters["fastqc"] = self.fastqcParameters

    def store_hisat2_parameters(self, args):
        """Updates hisat2 aligner parameters to the map of parameters to be store in a JSON file
        args -- set of parsed arguments
        """
        self.hisat2Parameters["hisat2_cores"] = args.hisat2_cores
        self.allParameters["hisat2"] = self.hisat2Parameters

    def store_kraken2_parameters(self, args):
        """Updates kraken2 to the map of parameters to be store in a JSON file
        args -- set of parsed arguments
        """
        self.Kraken2Parameters["kraken2_cores"] = args.kraken2_cores
        self.allParameters["Kraken2"] = self.Kraken2Parameters


def main() -> int:
    # 1.Create object class Configuration File
    configManager = CreateConfigurationFile()

    # 2.Create object for argument parsinng
    parser = argparse.ArgumentParser(
        prog="metaline-generate-config",
        description="Create a configuration json file for the meTAline pipeline.",
    )

    # 2.1 Updates arguments and parsing
    configManager.register_parameter(parser)

    args = parser.parse_args()

    # 2.2 Check Parameters
    configManager.check_parameters(args, parser)

    # 3. store arguments to super map structure
    configManager.store_general_parameters(args)
    configManager.store_input_parameters(args)
    configManager.store_output_parameters(args)
    configManager.store_wildcard_parameters(args)
    configManager.store_trimmomatic_parameters(args)
    configManager.store_fastqc_parameters(args)
    configManager.store_hisat2_parameters(args)
    configManager.store_kraken2_parameters(args)

    # 4. Store JSON fileMGI_adapters.fa
    with open(args.configFile, "w") as output_file:
        json.dump(configManager.allParameters, output_file, indent=2)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
