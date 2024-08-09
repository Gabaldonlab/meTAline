#Author: Diego Fuentes and Olfat Khannous
#Contact email: olfat.khannous@bsc.es
#Barcelona
#Date:2024-07-23

###################
# BIOBAKERY_TOOLS (Metaphlan4 and Humann) #
###################

#This is a rule for both taxonomy and functional profiling, based on gene markers. It will only be run if there is previously a host removal (indicated by reference genome) and the user indicates a metaphlan database

rule metaphlan4:
    input:
        unmapped_in = rules.extracted_unmapped_reads.output.fastq
    output:
        out_bz2 = metaphlan4_out + sample +".bz2",
        out_sam = metaphlan4_out + sample + "sam.bz2",
        out_profile = metaphlan4_out + sample + "_profiled.txt",
        out_vsc = metaphlan4_out + sample + ".vsc.txt"
    params:
        outdir = config["Outputs"]["metaphlan4_out"],
        metaphlan_db = config["Inputs"]["metaphlan_db"], #/gpfs/projects/bsc40/current/okhannous/Metaphlan4/db
        metaphlan_index = config["Inputs"]["metaphlan_Index"] #mpa_vJun23_CHOCOPhlAnSGB_202307

    shell:
        """
        metaphlan --bowtie2db {params.metaphlan_db} --index {params.metaphlan_index} {input.unmapped_in} --input_type fastq --bowtie2out {output.out_bz2} -s {output.out_sam} --profile_vsc -o {output.out_profile} --nproc 4 --vsc_out {output.out_vsc};
        """

rule humann:
    input: 
        unmapped_f = rules.extracted_unmapped_reads.output.fastq,
        profiled_sample = rules.metaphlan4.output.out_profile
    output:
        outdir_h = directory(config["Outputs"]["metaphlan4_out"] + "{sample}")
    params:
        protein_db = config["Inputs"]["protein_db"],
        metaphlan_index = config["Inputs"]["metaphlan_Index"],
        metaphlan_db = config["Inputs"]["metaphlan_db"],
        n_db = config["Inputs"]["n_db"]

    shell:
        """
        humann --taxonomic-profile {input.profiled_sample} --input {input.unmapped_f} --output {output.outdir_h} --nucleotide-database {params.n_db} --protein-database {params.protein_db} --metaphlan-options "--bowtie2db {params.metaphlan_db}" --metaphlan-options "--index {params.metaphlan_index}" --bypass-translated-search;
        """