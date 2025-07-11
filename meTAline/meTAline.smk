#If you want to set up your environment PATH in the bashrc you can use the line below to source it
#shell.prefix("source ~/.bashrc; ")

#Author: Diego Fuentes, Dániel Majer and Olfat Khannous Lleiffe
#Contact email: olfat.khannous@bsc.es
#Barcelona
#Date:2025-06-13

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
extracted_fa_out = os.path.join(workingdir, config["Outputs"]["extracted_fa_out"])
ranalysis_out = os.path.join(workingdir, config["Outputs"]["ranalysis_out"])
metaphlan4_out = os.path.join(workingdir, config["Outputs"]["metaphlan4_out"])

#Three optional arguments, that we are assessing here:
reference_genome = config["Inputs"]["reference_genome"]

#Benchmark directory
benchmark_dir = os.path.join(workingdir, "Benchmark")
if not os.path.exists(benchmark_dir) :
    os.makedirs(benchmark_dir)

benchmark_plots_dir = os.path.join(benchmark_dir, "plots")
if not os.path.exists(benchmark_plots_dir) :
    os.makedirs(benchmark_plots_dir)


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
    read_processing.append(os.path.join(trimmomatic_out, f"{sample}_qc"))
else:
    read_processing.append(os.path.join(trimmomatic_out, f"{sample}_qc"))
    read_processing.append(os.path.join(alignment_out, f"{sample}.hisat2.bam"))
    read_processing.append(os.path.join(alignment_out, f"{sample}.hisat2.bam.bai"))
    read_processing.append(os.path.join(alignment_out, f"{sample}.unmapped.fastq.gz"))

#Rule by default, runs the complete pipeline
rule all:
    input:
        read_processing,
        kraken_out = os.path.join(kraken_out, f"{sample}.kraken2.txt"),
        kraken_report= os.path.join(kraken_out, f"{sample}.kraken2.report"),
        krona_file= os.path.join(krona_out, f"{sample}.krona"),
        krona_html= os.path.join(krona_out, f"{sample}.html"),
        ext1= os.path.join(extracted_fa_out, f"{sample}.1.fastq.gz"),
        ext2= os.path.join(extracted_fa_out, f"{sample}.2.fastq.gz"),
        phyloseq_object = os.path.join( ranalysis_out, f"{sample}.rds"),
        rich_object = os.path.join( ranalysis_out, f"{sample}_alpha_div.csv"),
        out_plot = os.path.join( ranalysis_out, f"{sample}_Barplot_phyla.jpeg"),
        out_profile = os.path.join(metaphlan4_out, f"{sample}_profiled.txt"),
        outdir_h = os.path.join(config["Outputs"]["metaphlan4_out"], sample),
        benchmark_plots = benchmark_plots_dir



#Rule in case you want to only perform trimming and qualiy assessment of the data
rule Trimming:
    input:
        read_processing

#Rule to deplete host reads, if a reference genome is provied
if reference_genome != None or reference_genome != "null":
    rule host_depletion:
        input:
            fastq = os.path.join(alignment_out, f"{sample}.unmapped.fastq.gz")


#This rule is used to perform the taxonomy assignment based on the k-mer based approach
rule kmer_taxonomy:
    input:
        kraken_out = os.path.join(kraken_out, f"{sample}.kraken2.txt"),
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

# Private rule. It exists only to include the custom plots of the jobs in the html report.
rule __benchmark_report:
    input:
        benchmark_dir
    output:
       report(
            directory(benchmark_plots_dir),
            patterns=["{name}.svg"],
            category="Resource benchmark"
       )
    shell:
       "sleep 0.1"
