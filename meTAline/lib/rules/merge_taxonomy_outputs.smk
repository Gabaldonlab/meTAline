#Author: Diego Fuentes, Dániel Majer and Olfat Khannous Lleiffe
#Contact email: olfat.khannous@bsc.es
#Barcelona
#Date:2026-02-03

###################
# Merging taxonomy assignment of different samples
###################

#Kraken2 reports to biom format that will be exported to phyloseq
rule merge_kraken2: 
    input:
        report= rules.Kraken2.output.report ##Take all the kraken2 reports
    output:
        biom_file = protected(os.path.join(merged_out, ".Kraken2.biom"))
    benchmark:
        os.path.join(benchmark_dir, ("conversion_bio.benchmark.txt"))
    threads: config["Kraken2"]["kraken2_cores"]
    shell:
        "kraken-biom {input.report} -o {output.biom_file} --fmt json"

#Metaphlan4 reports to txt table
rule merge_metaphlan4:
    input:
        metaphlan_report = rules.metaphlan4.output.out_profile ##Take all the metaphlan4 reports
    output:
        metaphlan_merged = protected(os.path.join(merged_out, "_abundance_table.txt"))
    benchmark:
        os.path.join(benchmark_dir, ("merging_metaphlan4.benchmark.txt"))
    shell:
        "merge_metaphlan_tables.py {input.metaphlan_report} > {output.metaphlan_merged}"


#Sylph, merge the fungal identification by sylph to tsv table
rule sylph_fungi_merge:
    input:
        fungi_identified = #All the tsv that are create in filter_fungi --> for all the samples 
    output:
        out_fungi = protected(os.path.join(merged_out, ".tsv")) ##Take all the samples
    benchmark:
        os.path.join(benchmark_dir, ("merging_sylph_fungi.benchmark.txt"))
    shell:
        "python3 merge_sylph_taxprof.py {input.fungi_identified} --column sequence_abundance -o {output.out_fungi}"

#Sylph, merge the bacterial identification by sylph to tsv table
rule sylph_bacteria_merge:
    input:
        bacteria_identified = #All the tsv that are created in filter_bacteria 
    output:
        bacteria_identified = protected(os.path.join(merged_out, "bacteria.tsv")) ##Take all the samples
    benchmark:
        os.path.join(benchmark_dir, ("merging_sylph_bacteria.benchmark.txt"))