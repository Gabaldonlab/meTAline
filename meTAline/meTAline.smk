#If you want to set up your environment PATH in the bashrc you can use the line below to source it
#shell.prefix("source ~/.bashrc; ")


#Author: Diego Fuentes, Dániel Majer and Olfat Khannous Lleiffe
#Contact email: olfat.khannous@bsc.es
#Barcelona
#Date:2025-03-27

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

#Useful when updating versions of the pipeline
pipeline_version = "v" + str(config["Parameters"]["version"])

workingdir = config["Parameters"]["basedir"]
alignment_out = os.path.join(workingdir, config["Outputs"]["alignment_out"])
trimmomatic_out = os.path.join(workingdir, config["Outputs"]["trimmomatic_out"])
kraken_out = os.path.join(workingdir, config["Outputs"]["kraken_out"])
krona_out = os.path.join(workingdir, config["Outputs"]["krona_out"])
extracted_fa_out = os.path.join(workingdir, config["Outputs"]["extracted_fa_out"])
ranalysis_out = os.path.join(workingdir, config["Outputs"]["ranalysis_out"])
metaphlan4_out = os.path.join(workingdir, config["Outputs"]["metaphlan4_out"])

#Three optional arguments, that we are assessing here:
reference_genome = config["Inputs"]["reference_genome"]


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
if not os.path.exists(ranalysis_out):
    os.makedirs(ranalysis_out)
if not os.path.exists(metaphlan4_out):
    os.makedirs(metaphlan4_out)


# file wildcard
files = config["Wildcards"]["fastq_prefix"]


################
# RULE MODULES #
################

include: "lib/rules/trimming.smk"
include: "lib/rules/host_depletion.smk"
include: "lib/rules/kmer_taxonomy.smk"
include: "lib/rules/RAnalysis.smk"
include: "lib/rules/Biobakery.smk"



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
        krona_file= krona_out + sample + ".krona",
        krona_html= krona_out + sample + ".html",
        ext1= extracted_fa_out + sample + ".1.fastq.gz",
        ext2= extracted_fa_out + sample + ".2.fastq.gz",
        phyloseq_object =  ranalysis_out + sample + ".rds",
        rich_object =  ranalysis_out + sample + "_alpha_div.csv",
        out_plot =  ranalysis_out + sample + "_Barplot_phyla.jpeg",
        out_profile = metaphlan4_out + sample + "_profiled.txt", 
        outdir_h = config["Outputs"]["metaphlan4_out"] + sample


#Rule in case you want to only perform trimming and qualiy assessment of the data
rule Trimming:
    input:
        read_processing
     

#Rule to deplete host reads, if a reference genome is provied
if reference_genome != None or reference_genome != "null":
    rule host_depletion:
        input:
            fastq = alignment_out + sample +".unmapped.fastq.gz"


#This rule is used to perform the taxonomy assignment based on the k-mer based approach
rule kmer_taxonomy:
    input:
        kraken_out= kraken_out + sample + ".kraken2.txt",
        kraken_report= kraken_out + sample + ".kraken2.report",
        krona_file= krona_out + sample + ".krona",
        krona_html= krona_out + sample + ".html"


#This rule is used to make a basic R analysis
rule R_Analysis:
    input:
        biom_file= kraken_out + sample + ".Kraken2.biom",
        biom= rules.convert_biom.output.biom_file,
        phyloseq_object =  ranalysis_out + sample + ".rds",
        rich_object =  ranalysis_out + sample + "_alpha_div.csv",
        out_plot =  ranalysis_out + sample + "_Barplot_phyla.jpeg"


#Rule to run the BioBakery rules: Taxonomy assignment by MetaPhlAn4 and functional profiling by HUMAnN
rule BioBakery:
    input: 
        out_profile = metaphlan4_out + sample + "_profiled.txt",
        out_vsc = metaphlan4_out + sample + ".vsc.txt",
        outdir_h = config["Outputs"]["metaphlan4_out"] + sample