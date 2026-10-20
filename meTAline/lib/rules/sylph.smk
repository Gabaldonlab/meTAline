#Author: Laura McCluskey, Diego Fuentes, Dániel Majer and Olfat Khannous Lleiffe
#Contact email: olfat.khannous@bsc.es
#Barcelona
#Date:2026-05-12

###################
# SYLPH  #Metagenomic profiling with ANIs
###################

#This is a rule for taxonomy assignment using sylph


#Check whether we use the filtered reads or the unmapped reads, which happens when the host genome is provided. In case for instance of an environmental sample the trimmed reads are the ones that will be used.
if reference_genome in [None, "null"] and reference_phix in [None, "null"]:
    read1_selected = rules.Concat_reads.output.concat1,
    read2_selected = rules.Concat_reads.output.concat2,
    sylph_taxonomy_suffix="1.fastq.gz.sylphmpa"

elif reference_genome not in [None, "null"] and reference_phix in [None, "null"]:
    read1_selected = rules.extracted_unmapped_reads.output.fastq,
    read2_selected = rules.extracted_unmapped_reads.output.fastq,
    sylph_taxonomy_suffix="unmapped.fastq.gz.sylphmpa"

elif reference_genome not in [None, "null"] and reference_phix not in [None, "null"]:
    read1_selected = rules.extract_phix_unmapped_reads.output.fastq_phix,
    read2_selected = rules.extract_phix_unmapped_reads.output.fastq_phix,
    sylph_taxonomy_suffix="clean.fastq.gz.sylphmpa"

elif reference_genome in [None, "null"] and reference_phix not in [None, "null"]:
    read1_selected = rules.remove_phix_reads.output.fastq_phix,
    read2_selected = rules.remove_phix_reads.output.fastq_phix,
    sylph_taxonomy_suffix="unmapped.fastq.gz.sylphmpa"

rule sylph_profiling:
    input:
        read1 = read1_selected,
        read2 = read2_selected
    output:
        out_tsv = os.path.join(config["Outputs"]["sylph_out"], f"{sample}.tsv")
    params:
        sylph_db = config["Inputs"]["sylph_db"],
        ani_score_sylph = config["Inputs"]["ani_score_sylph"]
    benchmark:
        os.path.join(benchmark_dir, f"{sample}.sylph.benchmark.txt")
    run:
        if reference_genome in [None, "null"] and reference_phix in [None, "null"]:
            shell("sylph profile {params.sylph_db} -1 {input.read1} -2 {input.read2} --estimate-read-counts -m {params.ani_score_sylph} -t 48 -o {output.out_tsv}")

        else:
            shell("sylph profile {params.sylph_db} {input.read1} --estimate-read-counts -m {params.ani_score_sylph} -t 48 -o {output.out_tsv}")


#This will take the sylph output and filter only Fungi that are greater than or equal to the ANI score given
rule add_taxonomy_to_sylph:
    input:
        profiled_sample = rules.sylph_profiling.output.out_tsv
    output:
        complete_taxonomy = os.path.join(config["Outputs"]["filtered_euk"], f"{sample}.{sylph_taxonomy_suffix}")
    params:
        sylph_taxonomy = config["Inputs"]["sylph_taxonomy"],
        out_prefix = os.path.join(config["Outputs"]["filtered_euk"], "")
    benchmark:
        os.path.join(benchmark_dir, f"{sample}.sylph_taxonomy.benchmark.txt")
    run:
       #Incorporation of taxonomy into sylph assignment
       shell("sylph-tax taxprof {input.profiled_sample} -t {params.sylph_taxonomy} -o {params.out_prefix}" )

#This will take the sylph result, in which we already added taxonomy, and filter only to have hits of Fungi at greater than or equal to the ani score given
rule filter_fungi:
    input:
        profile_taxo_sylph = rules.add_taxonomy_to_sylph.output.complete_taxonomy
    output:
        filtered_taxonomy_fungi = os.path.join(config["Outputs"]["filtered_euk"], f"{sample}.tsv")
    params:
        ani_score_fungi_filtering = config["Inputs"]["ani_score_fungi_filtering"]
    benchmark:
        os.path.join(benchmark_dir, f"{sample}.filter_fungi_sylph.benchmark.txt")
    run:
        #Filter the results, considering only Fungi hits
        shell("awk -F'[[:space:]]+' 'NR==2 || ($0 ~ /k__Fungi/ && $(NF-1)+0 >= {params.ani_score_fungi_filtering})' \
        \"{input.profile_taxo_sylph}\" > \"{output.filtered_taxonomy_fungi}\"")
        
rule filter_bacteria:
    input:
        profile_taxo_sylph = rules.add_taxonomy_to_sylph.output.complete_taxonomy
    output:
        filtered_taxonomy_bact = os.path.join(config["Outputs"]["filtered_bact"], f"{sample}.tsv")
    params:
        ani_score_bacteria_filtering = config["Inputs"]["ani_score_bacteria_filtering"]
    benchmark:
        os.path.join(benchmark_dir,  f"{sample}.filter_bact_sylph.benchmark.txt")
    run:
        #Filter the results, considering only bacteria hits
        shell("awk -F'[[:space:]]+' 'NR==2 || ($(NF-1)+0 >= {params.ani_score_bacteria_filtering})' \
        \"{input.profile_taxo_sylph}\" > \"{output.filtered_taxonomy_bact}\"")
