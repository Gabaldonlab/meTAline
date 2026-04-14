#If you want to set up your environment PATH in the bashrc you can use the line below to source it
#shell.prefix("source ~/.bashrc; ")

#Author: Diego Fuentes, Dániel Majer and Olfat Khannous Lleiffe
#Contact email: olfat.khannous@bsc.es
#Barcelona
#Date:2026-03-04

import os
import sys
import json
import glob

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
sylph_out = os.path.join(workingdir, config["Outputs"]["sylph_out"])
filtered_euk = os.path.join(workingdir, config["Outputs"]["filtered_euk"])
filtered_bact = os.path.join(workingdir, config["Outputs"]["filtered_bact"])
merged_out = os.path.join(workingdir, config["Outputs"]["merged_out"])
completed_samples_folder = os.path.join(merged_out, "completed_samples")



#Three optional arguments, that we are assessing here:
reference_genome = config["Inputs"]["reference_genome"]

#Benchmark directory
benchmark_dir = os.path.join(workingdir, "Benchmark")
os.makedirs(benchmark_dir, exist_ok=True)

benchmark_plots_dir = os.path.join(benchmark_dir, "plots")
os.makedirs(benchmark_plots_dir, exist_ok=True)


#Create the necessary dir. Do not create alignment_out directory if there is no reference in the config
if reference_genome != None or reference_genome != "null":
    os.makedirs(alignment_out, exist_ok=True)
os.makedirs(trimmomatic_out, exist_ok=True)
os.makedirs(kraken_out, exist_ok=True)
os.makedirs(krona_out, exist_ok=True)
os.makedirs(extracted_fa_out, exist_ok=True)
os.makedirs(ranalysis_out, exist_ok=True)
os.makedirs(metaphlan4_out, exist_ok=True)
os.makedirs(sylph_out, exist_ok=True)
os.makedirs(filtered_euk, exist_ok=True)
os.makedirs(filtered_bact, exist_ok=True)
os.makedirs(merged_out, exist_ok=True)
os.makedirs(completed_samples_folder, exist_ok=True)

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
include: "lib/rules/sylph.smk"
include: "lib/rules/merge_taxonomy_outputs.smk"


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

# need to include the outputs of all smaples in the rule all

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
        out_tsv = os.path.join(sylph_out, f"{sample}.tsv"),
        complete_taxonomy = os.path.join(config["Outputs"]["filtered_euk"], f"{sample}.unmapped.fastq.gz.sylphmpa"),
        filtered_taxonomy_fungi = os.path.join(config["Outputs"]["filtered_euk"], f"{sample}.tsv"),
        filtered_taxonomy_bact = os.path.join(config["Outputs"]["filtered_bact"], f"{sample}.tsv"),
        merged_kraken = os.path.join(completed_samples_folder, f"{sample}_kraken_merged_flag.txt"),
        merged_metaphlan = os.path.join(completed_samples_folder, f"{sample}_metaphlan4_merged_flag.txt"),
        merged_sylph_fungi = os.path.join(completed_samples_folder, f"{sample}_sylph_fungi_merged_flag.txt"),
        merged_sylph_bacteria = os.path.join(completed_samples_folder, f"{sample}_sylph_bacteria_merged_flag.txt"),
        benchmark_plots = benchmark_plots_dir



#Rule in case you want to only perform trimming and qualiy assessment of the data
rule Trimming:
    input:
        read_processing

#Rule to deplete host reads, if a reference genome is provided
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

#Rule to run the sylph taxonomic assignment
rule sylph:
    input:
        out_tsv = os.path.join(sylph_out, f"{sample}.tsv"),
        complete_taxonomy = os.path.join(config["Outputs"]["filtered_euk"], f"{sample}.unmapped.fastq.gz.sylphmpa"),
        filtered_taxonomy_fungi = os.path.join(config["Outputs"]["filtered_euk"], f"{sample}.tsv"),
        filtered_taxonomy_bact = os.path.join(config["Outputs"]["filtered_bact"], f"{sample}.tsv")


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
