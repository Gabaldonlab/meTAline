#Author: Laura McCluskey, Diego Fuentes, Dániel Majer and Olfat Khannous Lleiffe
#Contact email: olfat.khannous@bsc.es
#Barcelona
#Date:2026-05-12

######################
# TAXONOMY ABUNDANCE #
######################

#Check whether we use the filtered reads or the unmapped reads when the host genome is provided

if reference_genome in [None, "null"] and reference_phix in [None, "null"]:
    read1_selected = rules.Concat_reads.output.concat1,
    read2_selected = rules.Concat_reads.output.concat2

elif reference_genome not in [None, "null"] and reference_phix in [None, "null"]:
    read1_selected = rules.extracted_unmapped_reads.output.fastq,
    read2_selected = rules.extracted_unmapped_reads.output.fastq

elif reference_genome not in [None, "null"] and reference_phix not in [None, "null"]:
    read1_selected = rules.extract_phix_unmapped_reads.output.fastq_phix,
    read2_selected = rules.extract_phix_unmapped_reads.output.fastq_phix

elif reference_genome in [None, "null"] and reference_phix not in [None, "null"]:
    read1_selected = rules.remove_phix_reads.output.fastq_phix,
    read2_selected = rules.remove_phix_reads.output.fastq_phix


#Taxonomy assignment with kraken2 (k-mer based approach)
rule Kraken2:
    input:
        read1 = read1_selected,
        read2 = read2_selected
    output:
        report = protected(os.path.join(kraken_out, f"{sample}.kraken2.report")),
        out  =  protected(os.path.join(kraken_out, f"{sample}.kraken2.txt.gz"))
    params:
        database = config["Inputs"]["krakendb"],
        uncompressed_file  = (os.path.join(kraken_out, f"{sample}.kraken2.txt"))
    benchmark:
        os.path.join(benchmark_dir, (sample +".kraken2.benchmark.txt"))
    threads:
        config["Kraken2"]["kraken2_cores"]
    #Run interpretes the following block as python code, keep python synthaxes
    run:
        if reference_genome in [None, "null"] and reference_phix in [None, "null"]:
            #Running kraken2 for paired data
            shell(
            """
            kraken2 --threads {threads} --db {params.database}  --paired {input.read1} {input.read2} --use-names --report {output.report} --output {params.uncompressed_file}
            gzip {params.uncompressed_file}
            """)
        else:
            #Running kraken2 for read extracted from BAM file
            shell(
            """
            kraken2 --threads {threads} --db {params.database} {input.read1} --use-names --report {output.report} --output {params.uncompressed_file}
            gzip {params.uncompressed_file}
            """)

#Krona representation of the kraken2 assignment
rule Krona:
    input:
        report = rules.Kraken2.output.report
    output:
        krona_file = protected(os.path.join(krona_out, f"{sample}.krona")),
        krona_html = protected(os.path.join(krona_out, f"{sample}.html"))
    threads: config["Kraken2"]["kraken2_cores"]
    benchmark:
        os.path.join(benchmark_dir, (sample +".krona.benchmark.txt"))
    shell:
        "kreport2krona.py -r  {input.report} -o {output.krona_file}; ktImportText {output.krona_file} -o {output.krona_html}"
