#If you want to set up your environment PATH in the bashrc you can use the line below to source it
#shell.prefix("source ~/.bashrc; ")

#Author: Laura McCluskey, Diego Fuentes, Dániel Majer and Olfat Khannous Lleiffe
#Contact email: olfat.khannous@bsc.es
#Barcelona
#Date:2026-05-12

import os
import sys

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
ranalysis_out = os.path.join(workingdir, config["Outputs"]["ranalysis_out"])
metaphlan4_out = os.path.join(workingdir, config["Outputs"]["metaphlan4_out"])
sylph_out = os.path.join(workingdir, config["Outputs"]["sylph_out"])
filtered_euk = os.path.join(workingdir, config["Outputs"]["filtered_euk"])
filtered_bact = os.path.join(workingdir, config["Outputs"]["filtered_bact"])
sample_dir=os.path.join(workingdir, sample)
reference_phix = config["Inputs"]["reference_phix"]

#Benchmark directory
benchmark_dir = os.path.join(workingdir, "Benchmark")

os.makedirs(sample_dir, exist_ok=True)
os.makedirs(trimmomatic_out, exist_ok=True)
os.makedirs(kraken_out, exist_ok=True)
os.makedirs(krona_out, exist_ok=True)
os.makedirs(ranalysis_out, exist_ok=True)
os.makedirs(metaphlan4_out, exist_ok=True)
os.makedirs(sylph_out, exist_ok=True)
os.makedirs(filtered_euk, exist_ok=True)
os.makedirs(filtered_bact, exist_ok=True)


#Create the necessary dir. Do not create alignment_out directory if there is no reference in the config
reference_genome = config["Inputs"]["reference_genome"]
if reference_genome != None or reference_genome != "null":
    os.makedirs(alignment_out, exist_ok=True)


# file wildcard
files = config["Wildcards"]["fastq_prefix"]


################
# RULE MODULES #
################

include: "lib/rules/trimming.smk"
include: "lib/rules/host_depletion.smk"
include: "lib/rules/PhiX_removal.smk"
include: "lib/rules/kmer_taxonomy.smk"
include: "lib/rules/RAnalysis.smk"
include: "lib/rules/Biobakery.smk"
include: "lib/rules/sylph.smk"


##############
# MAIN RULES #
##############

read_processing = list()

if reference_genome == None or reference_genome == "null":
    read_processing.append(os.path.join(trimmomatic_out, f"{sample}_qc"))
else:
    read_processing.append(os.path.join(trimmomatic_out, f"{sample}_qc"))
    read_processing.append(os.path.join(alignment_out, f"{sample}.hisat2.bam"))
    read_processing.append(os.path.join(alignment_out, f"{sample}.hisat2.bam.bai"))
    read_processing.append(os.path.join(alignment_out, f"{sample}.unmapped.fastq.gz"))

if reference_genome not in [None, "null"] and reference_phix in [None, "null"]:
    sylph_taxonomy_suffix="unmapped.fastq.gz.sylphmpa"
if reference_genome not in [None, "null"] and reference_phix not in [None, "null"]:
    sylph_taxonomy_suffix="clean.fastq.gz.sylphmpa"
if reference_genome in [None, "null"] and reference_phix in [None, "null"]:
    sylph_taxonomy_suffix="1.fastq.gz.sylphmpa"

#Rule by default, runs the complete pipeline
rule all:
    input:
        read_processing,
        kraken_output = os.path.join(kraken_out, f"{sample}.kraken2.txt.gz"),
        kraken_report= os.path.join(kraken_out, f"{sample}.kraken2.report"),
        krona_file= os.path.join(krona_out, f"{sample}.krona"),
        krona_html= os.path.join(krona_out, f"{sample}.html"),
        phyloseq_object = os.path.join( ranalysis_out, f"{sample}.rds"),
        rich_object = os.path.join( ranalysis_out, f"{sample}_alpha_div.csv"),
        out_plot = os.path.join( ranalysis_out, f"{sample}_Barplot_phyla.jpeg"),
        out_profile = os.path.join(metaphlan4_out, f"{sample}_profiled.txt"),
        outdir_h = os.path.join(config["Outputs"]["metaphlan4_out"], sample),
        out_tsv = os.path.join(sylph_out, f"{sample}.tsv"),
        complete_taxonomy = os.path.join(config["Outputs"]["filtered_euk"], f"{sample}.{sylph_taxonomy_suffix}"),
        filtered_taxonomy_fungi = os.path.join(config["Outputs"]["filtered_euk"], f"{sample}.tsv"),
        filtered_taxonomy_bact = os.path.join(config["Outputs"]["filtered_bact"], f"{sample}.tsv"),

#Rule in case you want to only perform trimming and qualiy assessment of the data
rule trimming:
    input:
        read_processing


#Rule to deplete host reads, if a reference genome is provided
if reference_genome not in [None, "null"] and reference_phix in [None, "null"]:
    rule host_depletion:
        input:
            fastq = os.path.join(alignment_out, f"{sample}.unmapped.fastq.gz")

elif reference_genome not in [None, "null"] and reference_phix not in [None, "null"]:
    rule host_depletion:
        input:
            fastq_phix = os.path.join(alignment_out, f"{sample}.clean.fastq.gz")

elif reference_genome in [None, "null"] and reference_phix not in [None, "null"]:
    rule PhiX_removal:
        input:
            fastq_phix= os.path.join(alignment_out, f"{sample}.unmapped.fastq.gz")


#This rule is used to perform the taxonomy assignment based on the k-mer based approach
rule kmer_taxonomy:
    input:
        kraken_output = os.path.join(kraken_out, f"{sample}.kraken2.txt.gz"),
        kraken_report = os.path.join(kraken_out, f"{sample}.kraken2.report"),
        krona_file = os.path.join(krona_out, f"{sample}.krona"),
        krona_html = os.path.join(krona_out, f"{sample}.html")


#This rule is used to make a basic R analysis
rule R_Analysis:
    input:
        biom_file = os.path.join(kraken_out, f"{sample}.Kraken2.biom"),
        biom = rules.convert_biom.output.biom_file,
        phyloseq_object = os.path.join( ranalysis_out, f"{sample}.rds"),
        rich_object = os.path.join( ranalysis_out, f"{sample}_alpha_div.csv"),
        out_plot = os.path.join( ranalysis_out, f"{sample}_Barplot_phyla.jpeg")


#Rule to run the BioBakery rules: Taxonomy assignment by MetaPhlAn4 and functional profiling by HUMAnN
rule BioBakery:
    input:
        out_profile = os.path.join(metaphlan4_out, f"{sample}_profiled.txt"),
        out_vsc = os.path.join(metaphlan4_out, f"{sample}.vsc.txt"),
        outdir_h = os.path.join(config["Outputs"]["metaphlan4_out"], sample)


#Rule to run the sylph taxonomic assignment
rule sylph:
    input:
        out_tsv = os.path.join(sylph_out, f"{sample}.tsv"),
        complete_taxonomy = os.path.join(config["Outputs"]["filtered_euk"], f"{sample}.{sylph_taxonomy_suffix}"),
        filtered_taxonomy_fungi = os.path.join(config["Outputs"]["filtered_euk"], f"{sample}.tsv"),
        filtered_taxonomy_bact = os.path.join(config["Outputs"]["filtered_bact"], f"{sample}.tsv")

