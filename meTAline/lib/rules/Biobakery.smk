#Author: Diego Fuentes, Dániel Majer and Olfat Khannous Lleiffe
#Contact email: olfat.khannous@bsc.es
#Barcelona
#Date:2025-03-27

###################
# BIOBAKERY_TOOLS (Metaphlan4 and Humann) #Taxonomy and functional profiling based on marker genes.
###################

#This is a rule for both taxonomy and functional profiling, based on gene markers.


#Check whether we use the filtered reads or the unmapped reads, which happens when the host genome is provided. In case for instance of an environmental sample the trimmed reads are the ones that will be used. s
if reference_genome == None or reference_genome == "null":
    read1_selected = rules.Concat_reads.output.concat1,
    read2_selected = rules.Concat_reads.output.concat2
else:
    read1_selected = rules.extracted_unmapped_reads.output.fastq,
    read2_selected = rules.extracted_unmapped_reads.output.fastq

rule metaphlan4:
    input:
        read1 = read1_selected,
        read2 = read2_selected
    output:
        out_bz2 = os.path.join(metaphlan4_out, f"{sample}.bz2"),
        out_sam = os.path.join(metaphlan4_out, f"{sample}sam.bz2"),
        out_profile = os.path.join(metaphlan4_out, f"{sample}_profiled.txt"),
        out_vsc = os.path.join(metaphlan4_out, f"{sample}.vsc.txt")
    params:
        outdir = config["Outputs"]["metaphlan4_out"],
        metaphlan_db = config["Inputs"]["metaphlan_db"],
        metaphlan_index = config["Inputs"]["metaphlan_Index"]

    #Run interpretes the following block as python code, keep python synthax
    run:
        if reference_genome == None or reference_genome == "null":
            #Running kraken2 for paired data
            shell("metaphlan --bowtie2db {params.metaphlan_db} --index {params.metaphlan_index} {input.read1},{input.read2} --input_type fastq --bowtie2out {output.out_bz2} -s {output.out_sam} --profile_vsc -o {output.out_profile} --nproc 4 --vsc_out {output.out_vsc}")

        else:
            #Running metaphlan for reads extracted from BAM file
            shell("metaphlan --bowtie2db {params.metaphlan_db} --index {params.metaphlan_index} {input.read1} --input_type fastq --bowtie2out {output.out_bz2} -s {output.out_sam} --profile_vsc -o {output.out_profile} --nproc 4 --vsc_out {output.out_vsc}")



# Check again if reference genome is provided
if reference_genome == None or reference_genome == "null":
    read1_selected = rules.Concat_reads.output.concat1,
    read2_selected = rules.Concat_reads.output.concat2
else:
    read1_selected = rules.extracted_unmapped_reads.output.fastq,
    read2_selected = rules.extracted_unmapped_reads.output.fastq

rule humann:
    input:
        read1 = read1_selected,
        read2 = read2_selected,
        profiled_sample = rules.metaphlan4.output.out_profile
    output:
        combined_reads = os.path.join(config["Outputs"]["metaphlan4_out"], f"{sample}.fastq"),
        outdir_h = directory(os.path.join(config["Outputs"]["metaphlan4_out"], f"{sample}"))
    params:
        protein_db = config["Inputs"]["protein_db"],
        metaphlan_index = config["Inputs"]["metaphlan_Index"],
        metaphlan_db = config["Inputs"]["metaphlan_db"],
        n_db = config["Inputs"]["n_db"]

    run:
        if reference_genome == None or reference_genome == "null":
            # Concatenate paired-end reads using zcat if gzipped
            shell("zcat {input.read1} {input.read2} > {output.combined_reads}")
            # Run humann on the concatenated reads
            shell("humann --taxonomic-profile {input.profiled_sample} --input {output.combined_reads} --output {output.outdir_h} --nucleotide-database {params.n_db} --protein-database {params.protein_db} --metaphlan-options '--bowtie2db {params.metaphlan_db} --index {params.metaphlan_index}' --bypass-translated-search")
        else:
            # Running humann for unmapped reads to the human genome
            #shell("humann --taxonomic-profile {input.profiled_sample} --input {input.read1} --output {output.outdir_h} --nucleotide-database {params.n_db} --protein-database {params.protein_db} --metaphlan-options '--bowtie2db {params.metaphlan_db} --index {params.metaphlan_index}' --bypass-translated-search")
            shell("zcat {input.read1} > {output.combined_reads}")
            shell("humann --taxonomic-profile {input.profiled_sample} --input {output.combined_reads} --output {output.outdir_h} --nucleotide-database {params.n_db} --protein-database {params.protein_db} --metaphlan-options '--bowtie2db {params.metaphlan_db}' --metaphlan-options '--index {params.metaphlan_index}' --bypass-translated-search;")
