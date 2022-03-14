#!/usr/bin/env python3
import os
import json
import argparse
import sys
import re

#Author: Diego Fuentes
#Contact email: diegofupa@gmail.com
#Date:2021-11-01

barcodes = []    #List to store the barcodes

#######################
###CONFIG FILE CLASS###
#######################
class CreateConfigurationFile(object):
    """Class which manages Configuration file Manager"""
      
    def __init__(self):
        """Class constructor of the pipeline"""
        #GENERAL PARAMETERS

        self.configFile = "config.json"                      #Name of the json configuration file to be created.
        self.version = 1                                     #Pipeline version
        self.logs_dir = "WGS_logs"                           #Directory to keep all the log files
        self.extension = "fastq.gz"                          #Extension of the illumina raw read files
        self.sample_barcode = None                           #Sample barcode 
        self.basedir = self.sample_barcode                   #Base directory for the pipeline run

        #INPUT PARAMETERS

        self.reads_directory = None              #Directory where the illumina fastqs are stored
        self.reference_genome = None             #Indexed reference genome (bowtie2 index) path
        self.krakendb = None                     #Kraken2 DB path that is going to be used
        self.kmer_dist = None                    #Kraken2 kmer distribution path that is going to be used for Bracken
        self.taxid = 0                           #Default taxid used to extract reads, based on the kraken report. By default 0 -> unclassified

        #OUTPUT PARAMETERS

        self.alignment_out = "BAM"                     #Out directory of the alignment step
        self.trimmomatic_out =  "TRIMMOMATIC"          #Out directory of the sv calls
        self.kraken_out = "KRAKEN_ASSIGN"              #Out directory of the Kraken2 reports
        self.krona_out = "KRONA_HTML"                  #Out directory of the Krona visualization tool
        self.extracted_fa_out = "EXTRACTED_FASTA"      #Out directory of the extracted fasta
    
        #WILDCARD PARAMETER

        self.fastqs = None                       #List with basename of the fastqs if any

        #TRIMMOMATIC PARAMETERS

        self.trimmo_cores = 4                  #Number of threads to run the trimmomatic
        self.leading= 3                        #leading parameter
        self.trailing= 3                       #trailing parameter
        self.slidingwindow = "4:15"              #sliding window parameter (two values, ie: 4:15)
        self.minlen = 35                       #min read lenght paramter
        self.illuminaclip = "TruSeq3-PE-2.fa:2:30:10"   #Illuminaclip

        #BOWTIE PARAMETERS (set of default parameters for local-sensitive)

        self.bowtie2_cores = 8                   #Number of threads to run the ngmlr aligner
        self.bowtie2_D_param = 15                #D parameter: give up extending after <int> failed extends in a row 
        self.bowtie2_R_param = 2                 #R parameter: for reads w/ repetitive seeds, try <int> sets of seeds
        self.bowtie2_N_param = 0                 #N parameter: max # mismatches in seed alignment
        self.bowtie2_L_param = 20                #L parameter: length of seed substrings; must be >3, <32
        self.bowtie2_i_param = "S,1,0.75"        #i parameter: interval between seed substrings w/r/t read len (S,1,1.15)
        self.bowtie2_score_min = "G,20,8"        #score min: min acceptable alignment score w/r/t read length(G,20,8 for local, L,-0.6,-0.6 for end-to-end)

        #KRAKEN2 PARAMETERS

        self.kraken2_cores = 16

###
        #DICTIONARIES
        self.allParameters = {}
        self.generalParameters = {}
        self.inputParameters = {}
        self.outputParameters = {}
        self.wildcardParameters = {}
        self.TrimmomaticParameters = {}
        self.Bowtie2Parameters = {}
        self.Kraken2Parameters = {}
        
####

    def register_parameter(self, parser):
        """Register all parameters with the given
        argparse parser"""
        self.register_general(parser)
        self.register_input(parser)
        self.register_output(parser)
        self.register_wildcards(parser)
        self.register_trimmomatic(parser)
        self.register_bowtie2(parser)
        self.register_kraken2(parser)

    def register_general(self, parser):
        """Register all general parameters with the given
        argparse parser
        parser -- the argparse parser
        """
        general_group = parser.add_argument_group('General Parameters')
        general_group.add_argument('--configFile', dest="configFile", metavar="configFile", help='Configuration JSON to be generated. Default %s.' % self.configFile)
        general_group.add_argument('--version', type=int, dest="version", metavar="version", default=self.version, help='Pipeline run version. Default %s.' % self.version)
        general_group.add_argument('--logs-dir', dest="logs_dir", metavar="logs_dir", help='Directory to keep all the log files. Default WGS_logs.')
        general_group.add_argument('--extension', dest="extension", metavar="extension", default=self.extension, help='Extension of the illumina raw read files. Default %s.' % self.extension)
        general_group.add_argument('--sample-barcode', dest="sample_barcode", metavar="sample_barcode", help='Sample barcode. Default %s.' % self.sample_barcode)
        general_group.add_argument('--basedir', dest="basedir", metavar="basedir", help='Base directory for the pipeline run. Default %s.' % self.basedir)

    def register_input(self, parser):
        """Register all input parameters with the given
        argparse parser
        parser -- the argparse parser
        """
        input_group = parser.add_argument_group('Inputs')
        input_group.add_argument('--reads-directory', dest="reads_directory", metavar="reads_directory", help='Directory where the Illumina fastqs are stored. Default %s.' % self.reads_directory)
        input_group.add_argument('--reference-genome', dest="reference_genome", metavar="reference_genome", help='Indexed reference genome path (ie: reference/Human_index). Your path is  %s.' % self.reference_genome)
        input_group.add_argument('--krakendb', dest="krakendb", metavar="krakendb", help='Kraken2 database used for the taxonomic assignation. Your path is  %s.' % self.krakendb) 
        input_group.add_argument('--kmer_dist', dest="kmer_dist", metavar="kmer_dist", help='Kraken2 kmer distribution path used for Bracken. Your path is  %s.' % self.kmer_dist) 
        input_group.add_argument('--taxid', dest="taxid", metavar="taxid", default=self.taxid, help='Selects the taxid used for exctracting fasta reads. You can add more than one using space as the separator. Default \"%s\".' % self.taxid)

    def register_output(self, parser):
        """Register all output parameters with the given
        argparse parser
        parser -- the argparse parser
        """

        output_group = parser.add_argument_group('Outputs')
        output_group.add_argument('--alignment-out', dest="alignment_out", help='Out directory of the alignment step. Default "/%s"' % self.alignment_out)
        output_group.add_argument('--trimmomatic-out', dest="trimmomatic_out", help='Out directory of the trimmomatic output. Default "/%s"' % self.trimmomatic_out)
        output_group.add_argument('--kraken-out', dest="kraken_out", help='Out directory of the Kraken2 taxonomic assignation step. Default "/%s"' % self.kraken_out)
        output_group.add_argument('--krona-out', dest="krona_out", help='Out directory of the Krona visualization tool. Default "/%s"' % self.krona_out)
        output_group.add_argument('--extracted-fa-out', dest="extracted_fa_out", help='Out directory of the extracted fasta reads. Default "/%s"' % self.extracted_fa_out)


    def register_wildcards(self, parser):
        """Register all wildcards parameters with the given
        argparse parser
        parser -- the argparse parser
        """
        wildcards_group = parser.add_argument_group('Wildcards')
        wildcards_group.add_argument('--fastqs', dest="fastqs", metavar="fastqs", help='List with basename of the fastqs. Default %s' % self.fastqs)

    def register_trimmomatic(self, parser):
        """Register all trimmomatic parameters with the given
        argparse parser
        parser -- the argparse parser
        """
        trimmomatic_group = parser.add_argument_group('Trimmomatic parameters')
        trimmomatic_group.add_argument('--trimmo-cores', type=int, dest="trimmo_cores", metavar="trimmo_cores", default=self.trimmo_cores, help='Number of threads to run trimmomatic. Default %s.' % self.trimmo_cores)
        trimmomatic_group.add_argument('--leading', type=int, dest="leading", metavar="leading", default=self.leading, help='Leading parameter. Default %s.' % self.leading)
        trimmomatic_group.add_argument('--trailing', type=int, dest="trailing", metavar="trailing", default=self.trailing, help='Trailing parameter. Default %s.' % self.trailing)
        trimmomatic_group.add_argument('--slidingwindow', dest="slidingwindow", metavar="slidingwindow", default=self.slidingwindow, help='Sliding window parmeter. Default %s.' % self.slidingwindow)
        trimmomatic_group.add_argument('--minlen', type = int, dest="minlen", metavar="minlen", default=self.minlen, help='Minimum readl length parameter. Default %s.' % self.minlen)
        trimmomatic_group.add_argument('--illuminaclip', dest="illuminaclip", metavar="illuminaclip", default=self.illuminaclip, help='Illumina clip information. Default %s.' % self.illuminaclip)

    def register_bowtie2(self, parser):
        """Register all bowtie2 aligner parameters with the given
        argparse parser
        parser -- the argparse parser
        """
        bowtie2_group = parser.add_argument_group('Bowtie2 parameters')
        bowtie2_group.add_argument('--bowtie2-cores', type=int, dest="bowtie2_cores", metavar="bowtie2_cores", default=self.bowtie2_cores, help='Number of threads to run the bowtie2 aligner. Default %s.' % self.bowtie2_cores)
        bowtie2_group.add_argument('--bowtie2-D-param', type=float, dest="bowtie2_D_param", metavar="bowtie2_D_param", default=self.bowtie2_D_param, help='D parameter: give up extending after <int> failed extends in a row. Default %s.' % self.bowtie2_D_param)
        bowtie2_group.add_argument('--bowtie2-R-param', type=int, dest="bowtie2_R_param", metavar="bowtie2_R_param", default=self.bowtie2_R_param, help='R parameter: for reads w/ repetitive seeds, try <int> sets of seeds. Default %s.' % self.bowtie2_R_param)
        bowtie2_group.add_argument('--bowtie2-N-param', type=int, dest="bowtie2_N_param", metavar="bowtie2_N_param", default=self.bowtie2_N_param, help='N parameter: max # mismatches in seed alignment. Default %s.' % self.bowtie2_N_param)
        bowtie2_group.add_argument('--bowtie2-L-param', type = int, dest="bowtie2_L_param", metavar="bowtie2_L_param", default=self.bowtie2_L_param, help='L parameter: length of seed substrings; must be >3, <32. Default %s.' % self.bowtie2_L_param)
        bowtie2_group.add_argument('--bowtie2-i-param', dest="bowtie2_i_param", metavar="bowtie2_i_param", default=self.bowtie2_i_param, help='i parameter: interval between seed substrings w/r/t read len (S,1,1.15). Default %s.' % self.bowtie2_i_param)
        bowtie2_group.add_argument('--bowtie2-score-min', dest="bowtie2_score_min", metavar="bowtie2_score_min", default=self.bowtie2_score_min, help='score min: min acceptable alignment score w/r/t read length(G,20,8 for local, L,-0.6,-0.6 for end-to-end). Default %s.' % self.bowtie2_score_min)

    def register_kraken2(self, parser):
        """Register all sniffles sv caller parameters with the given
        argparse parser
        parser -- the argparse parser
        """
        kraken2_group = parser.add_argument_group('Kraken2 parameters')
        kraken2_group.add_argument('--kraken2-cores', type=int, dest="kraken2_cores", metavar="kraken2_cores", default=self.kraken2_cores, help='Number of threads to run the Kraken2 taxonomic assignation step. Default %s.' % self.kraken2_cores)
####

    def check_parameters(self,args):
        """Check parameters consistency
            
        args -- set of parsed arguments
        
        Note that by using os.path.exist it only checks if the file/directory exist but function may return False, 
        if permission is not granted to execute os.stat() on the requested file, even if the path exists. """

        if args.configFile==None:
            parser.print_help()
            sys.exit(-1)

        working_dir = os.getcwd() + "/"

        if args.sample_barcode == None:
            print("No sample_barcode specified. A barcode or identification is required")
            parser.print_help()
            sys.exit(-1)

        if args.basedir:
            args.basedir = os.path.abspath(args.basedir) + "/"
        else: 
            args.basedir = working_dir + "v" + str(args.version) + "/" + args.sample_barcode + "/"

        if args.logs_dir:
            args.logs_dir = os.path.abspath(args.logs_dir) + "/"
        else:
            args.logs_dir = args.basedir + self.logs_dir + "/"
            

        if args.reads_directory:
            args.reads_directory = os.path.abspath(args.reads_directory) + "/"
        else:
            args.reads_directory =  working_dir + "reads/Illumina/" + args.sample_barcode + "/"
        if not os.path.exists(args.reads_directory):
            print(args.reads_directory + " not found. The directory where the reads are located is required. Exiting now.")
            parser.print_help()
            sys.exit(-1)
            
        if args.reference_genome != None:
            args.reference_genome = os.path.abspath(args.reference_genome)
        #else:
        #    args.reference_genome = working_dir + "reference/genome.fa"

            if not os.path.exists(args.reference_genome):
                #print("The reference genome has been not provided or it has not been found in "+args.reference_genome+". Enabling the taxonomic assignation of a environmental sample")
                #parser.print_help()
                #sys.exit(-1)
                if not re.search('index', args.reference_genome):
                    print("A reference genome index (including index as suffix, ie:'HUMAN_index') has been not provided or it has not been found in "+args.reference_genome+". Enabling the taxonomic assignation of a environmental sample")
                    args.reference_genome = None

        if args.krakendb:
            args.krakendb = os.path.abspath(args.krakendb)
        if not os.path.exists(args.krakendb):
            print("The Kraken2 DB has not been provided. Check the your provided path "+args.krakendb+" . Note that this step is mandatory")
            parser.print_help()
            sys.exit(-1)

        if args.kmer_dist:
            args.kmer_dist = os.path.abspath(args.kmer_dist)
        if not os.path.exists(args.kmer_dist):
            print("The Kraken2 kmer distribution has not been provided. Check the your provided path "+args.kmer_dist+" . Note that this step is mandatory")
            parser.print_help()
            sys.exit(-1)

        if args.alignment_out:
            args.alignment_out = os.path.abspath(args.alignment_out) + "/"
        else:
            args.alignment_out = args.basedir + self.alignment_out + "/"

        if args.trimmomatic_out:
            args.trimmomatic_out = os.path.abspath(args.trimmomatic_out) + "/"
        else:
            args.trimmomatic_out = args.basedir  + self.trimmomatic_out + "/"

        if args.kraken_out:
            args.kraken_out = os.path.abspath(args.kraken_out) + "/"
        else:
            args.kraken_out = args.basedir  + self.kraken_out + "/"

        if args.krona_out:
            args.krona_out = os.path.abspath(args.krona_out) + "/"
        else:
            args.krona_out = args.basedir  + self.krona_out + "/"

        if args.extracted_fa_out:
            args.extracted_fa_out = os.path.abspath(args.extracted_fa_out) + "/"
        else:
            args.extracted_fa_out = args.basedir  + self.extracted_fa_out + "/"

        ##Assign wildcards

        if args.fastqs == None:
            for r, d, f in os.walk(args.reads_directory):
                for file in f:
                    if re.search('.fastq.gz', file) or re.search('.fq.gz', file):
                        if file.endswith('_1.fastq.gz') or file.endswith('_1.fq.gz'):
                            if file.endswith('_1.fastq.gz'):
                                a = file.replace('_1.fastq.gz','')
                            elif file.endswith('_1.fq.gz'):
                                a = file.replace('_1.fq.gz','')
                        elif file.endswith('_2.fastq.gz') or file.endswith('_2.fastq.gz'):
                            if file.endswith('_2.fastq.gz'):
                                a = file.replace('_2.fastq.gz','')
                            elif file.endswith('_2.fq.gz'):
                                a = file.replace('_2.fq.gz','')
                        barcodes.append(a)
                        if args.fastqs == None:
                            args.fastqs = a
                        else:
                            args.fastqs += "," + a
                    elif re.search('.fastq', file) or re.search('.fq', file):
                        if file.endswith('_1.fastq') or file.endswith('_1.fq'):
                            if file.endswith('_1.fastq'):
                                a = file.replace('_1.fastq','')
                            elif file.endswith('_1.fq'):
                                a = file.replace('_1.fq.','')
                        elif file.endswith('_2.fastq') or file.endswith('_2.fq'):
                            if file.endswith('_2.fastq'):
                                a = file.replace('_2.fastq','')
                            elif file.endswith('_2.fq'):
                                a = file.replace('_2.fq','')
                        barcodes.append(a)
                        if args.fastqs == None:
                            args.fastqs = a
                        else:
                            args.fastqs += "," + a

###

    def storeGeneralParameters(self,args):
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

    def storeInputParameters(self,args):
        """Updates input parameters to the map of parameters to be store in a JSON file
        args -- set of parsed arguments
        """

        self.inputParameters["reads_directory"] = args.reads_directory
        self.inputParameters["reference_genome"] = args.reference_genome
        self.inputParameters["krakendb"] = args.krakendb
        self.inputParameters["kmer_dist"] = args.kmer_dist
        self.inputParameters["taxid"] = args.taxid
        self.allParameters ["Inputs"] = self.inputParameters

    def storeOutputParameters(self,args):
        """Updates output parameters to the map of parameters to be store in a JSON file
        args -- set of parsed arguments
        """
        self.outputParameters["alignment_out"] = args.alignment_out
        self.outputParameters["trimmomatic_out"] = args.trimmomatic_out
        self.outputParameters["kraken_out"] = args.kraken_out
        self.outputParameters["krona_out"] = args.krona_out
        self.outputParameters["extracted_fa_out"] = args.extracted_fa_out
        self.allParameters ["Outputs"] = self.outputParameters

    def storeWildcardParameters(self,args):
        """Updates wildcard parameters to the map of parameters to be store in a JSON file
        args -- set of parsed arguments
        """
        self.wildcardParameters["fastqs"] = args.fastqs
        self.allParameters ["Wildcards"] = self.wildcardParameters

    def storeTrimmomaticParameters(self,args):
        """Updates minimap2 aligner parameters to the map of parameters to be store in a JSON file
        args -- set of parsed arguments
        """
        self.TrimmomaticParameters["trimmo_cores"] = args.trimmo_cores
        self.TrimmomaticParameters["leading"] = args.leading
        self.TrimmomaticParameters["trailing"] = args.trailing
        self.TrimmomaticParameters["slidingwindow"] = args.slidingwindow
        self.TrimmomaticParameters["minlen"] = args.minlen
        self.TrimmomaticParameters["illuminaclip"] = args.illuminaclip
        self.allParameters ["Trimmomatic"] = self.TrimmomaticParameters

    def storeBowtie2Parameters(self,args):
        """Updates Ngmlr aligner parameters to the map of parameters to be store in a JSON file
        args -- set of parsed arguments
        """
        self.Bowtie2Parameters["bowtie2_cores"] = args.bowtie2_cores
        self.Bowtie2Parameters["bowtie2_D_param"] = args.bowtie2_D_param
        self.Bowtie2Parameters["bowtie2_R_param"] = args.bowtie2_R_param
        self.Bowtie2Parameters["bowtie2_N_param"] = args.bowtie2_N_param
        self.Bowtie2Parameters["bowtie2_L_param"] = args.bowtie2_L_param
        self.Bowtie2Parameters["bowtie2_i_param"] = args.bowtie2_i_param
        self.Bowtie2Parameters["bowtie2_score_min"] = args.bowtie2_score_min
        self.allParameters ["Bowtie2"] = self.Bowtie2Parameters


    def storeKraken2Parameters(self,args):
        """Updates Sniffles SV caller parameters to the map of parameters to be store in a JSON file
        args -- set of parsed arguments
        """
        self.Kraken2Parameters["kraken2_cores"] = args.kraken2_cores
        self.allParameters ["Kraken2"] = self.Kraken2Parameters



#####

#1.Create object class Configuration File
configManager = CreateConfigurationFile()

#2.Create object for argument parsinng
parser = argparse.ArgumentParser(prog="create_configuration_file",
                description="Create a configuration json file for the MASV pipeline."
                )     

#2.1 Updates arguments and parsing
configManager.register_parameter(parser)

args = parser.parse_args()

#2.2 Check Parameters
configManager.check_parameters(args)

#3. store arguments to super map structure
configManager.storeGeneralParameters(args)
configManager.storeInputParameters(args)
configManager.storeOutputParameters(args)
configManager.storeWildcardParameters(args)
configManager.storeTrimmomaticParameters(args)
configManager.storeBowtie2Parameters(args)
configManager.storeKraken2Parameters(args)


        


###
#4. Store JSON file
with open(args.configFile, 'w') as of:
    json.dump(configManager.allParameters, of, indent=2)
