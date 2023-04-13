#If you want to set up your environment PATH in the bashrc you can use the line below to source it
#shell.prefix("source ~/.bashrc; ")


#Author: Diego Fuentes
#Contact email: diego.fuentes@bsc.es
#Barcelona
#Date:2022-10-27

import os
from datetime import datetime
import sys

date1 = str(datetime.now())
tmp = str.replace(date1," ",".") 
tmp2 = str.replace(tmp,":","")
date = str.replace(tmp2,"-","")

#First determine the config file, if changed by the user, need to either specified here or through snakemake commands (if using commands, comment the line below)
#configfile: os.path.join(workflow.basedir, "config.json")

##############
# PARAMETERS #
##############

#Global parameters

#Sample barcode or identifier data. Can be used as name of the project if not required
sample = config["Parameters"]["sample_barcode"]

#Usefull when updating versions of the pipeline
pipeline_version = "v" + str(config["Parameters"]["version"])

workingdir = config["Parameters"]["basedir"]
alignment_out = os.path.join(workingdir, config["Outputs"]["alignment_out"])
trimmomatic_out = os.path.join(workingdir, config["Outputs"]["trimmomatic_out"])
kraken_out = os.path.join(workingdir, config["Outputs"]["kraken_out"])
krona_out = os.path.join(workingdir, config["Outputs"]["krona_out"])
extracted_fa_out = os.path.join(workingdir, config["Outputs"]["extracted_fa_out"])
ranalysis_out = os.path.join(workingdir, config["Outputs"]["ranalysis_out"])

reference_genome = config["Inputs"]["reference_genome"]

#Benchmark directory
benchmark_dir = os.path.join(workingdir, "Benchmark")
if not os.path.exists(benchmark_dir):
    os.makedirs(benchmark_dir)

#logs directory
logs_dir = os.path.join(workingdir, config["Parameters"]["logs_dir"])
if not os.path.exists(logs_dir):
    os.makedirs(logs_dir)

#Create the necessary dir. Do not create alignment_out directory if there is no reference in the config
if reference_genome != None or reference_genome != "null":
    if not os.path.exists(alignment_out):
        os.makedirs(alignment_out)
if not os.path.exists(trimmomatic_out):
    os.makedirs(trimmomatic_out)
if not os.path.exists(kraken_out):
    os.makedirs(kraken_out)
if not os.path.exists(krona_out):
    os.makedirs(krona_out)
if not os.path.exists(extracted_fa_out):
    os.makedirs(extracted_fa_out)



# file wildcard
files = config["Wildcards"]["fastqs"] 


################
# RULE MODULES #
################

include: "lib/rules/trimming.smk"
include: "lib/rules/alignment.smk"
include: "lib/rules/taxonomy.smk"
include: "lib/rules/RAnalysis.smk"


##############
# MAIN RULES #
##############

read_processing = list()

if reference_genome == None or reference_genome == "null":
    read_processing.append(trimmomatic_out +sample+"_qc")
else:
    read_processing.append(trimmomatic_out +sample+"_qc")
    read_processing.append(alignment_out + sample +".hisat2.bam")
    read_processing.append(alignment_out + sample +".hisat2.bam.bai")
    read_processing.append(alignment_out + sample +".unmapped.fastq.gz")

#Rule by default, runs the complete pipeline
rule all:
    input:
        read_processing,
        kraken_out= kraken_out + sample + ".kraken2.txt",
        kraken_report= kraken_out + sample + ".kraken2.report",
        bracken_report= kraken_out + sample + ".bracken_abundance.txt",
        bracken_species_report = kraken_out + sample + ".kraken2_bracken_species.report",
        krona_file= krona_out + sample + ".krona",
        krona_html= krona_out + sample + ".html",
        ext1= extracted_fa_out + sample + ".1.fastq.gz",
        ext2= extracted_fa_out + sample + ".2.fastq.gz",
        biom_file= kraken_out + sample + ".Braken.biom",
        phyloseq_object =  ranalysis_out + sample + ".rds",
        rich_object =  ranalysis_out + sample + "_alpha_div.csv",
        out_plot =  ranalysis_out + sample + "_Barplot_phyla.jpeg"

    log:
        logs_dir + str(date) + "_" + sample +".all.rule.log"


if reference_genome != None or reference_genome != "null":
    rule alignment:
        input:
            read_processing
        log:
            logs_dir + str(date) + "_" + sample +".alignment.rule.log"


#This rule is used if you already have the concatenated/unmapped reads and you want to perform the assignation.
rule taxonomy_assignation:
    input:
        kraken_out= kraken_out + sample + ".kraken2.txt",
        kraken_report= kraken_out + sample + ".kraken2.report",
        bracken_report= kraken_out + sample + ".bracken_abundance.txt",
        bracken_species_report = kraken_out + sample + ".kraken2_bracken_species.report",
        krona_file= krona_out + sample + ".krona",
        krona_html= krona_out + sample + ".html",
        #ext1= extracted_fa_out + sample + ".1.fastq.gz",
        #ext2= extracted_fa_out + sample + ".2.fastq.gz",
        biom_file= kraken_out + sample + ".Braken.biom",
        phyloseq_object =  ranalysis_out + sample + ".rds",
        rich_object =  ranalysis_out + sample + "_alpha_div.csv",
        out_plot =  ranalysis_out + sample + "_Barplot_phyla.jpeg"
    log:
        logs_dir + str(date) + "_" + sample +".taxonomy_assignation.rule.log"

#This rule is used if you already have the kraken2 taxonomic assignation and you want to generate the krona images and extract reads
rule krona_and_reads:
    input:
        bracken_report= kraken_out + sample + ".bracken_abundance.txt",
        krona_file= krona_out + sample + ".krona",
        krona_html= krona_out + sample + ".html",
        ext1= extracted_fa_out + sample + ".1.fastq.gz",
        ext2= extracted_fa_out + sample + ".2.fastq.gz"
    log:
        logs_dir + str(date) + "_" + sample +".krona_and_reads.rule.log"

#This rule is used if you already have the Kraken2 and bracken assignment and you want to pass from their reports to biom format for further analysis
rule biom_format:
    input:
        bracken_species_report = kraken_out + sample + ".kraken2_bracken_species.report",
        biom_file= kraken_out + sample + ".Braken.biom"
    log:
        logs_dir + str(date) + "_" + sample +".biom_format.rule.log"
#This rule is used to make a basic R analysis
rule R_Analysis:
    input:
        biom_file= kraken_out + sample + ".Braken.biom",
        biom= rules.convert_biom.output.biom_file,
        phyloseq_object =  ranalysis_out + sample + ".rds",
        rich_object =  ranalysis_out + sample + "_alpha_div.csv",
        out_plot =  ranalysis_out + sample + "_Barplot_phyla.jpeg"
    log:
        logs_dir + str(date) + "_" + sample +".R_analysis.rule.log"
