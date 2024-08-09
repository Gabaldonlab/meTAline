#Author: Diego Fuentes and Olfat Khannous
#Contact email: olfat.khannous@bsc.es
#Barcelona
#Date:2024-07-17

######################
# TAXONOMY ABUNDANCE #
######################
#Check whether we use the filtered reads or the unmapped reads when the host genome is provided
if reference_genome == None or reference_genome == "null":
    read1_selected = rules.Concat_reads.output.concat1,
    read2_selected = rules.Concat_reads.output.concat2
else:
    read1_selected = rules.extracted_unmapped_reads.output.fastq,
    read2_selected = rules.extracted_unmapped_reads.output.fastq

rule Kraken2:
    input:
        read1 = read1_selected,
        read2 = read2_selected
    output:
        report = protected(kraken_out + sample + ".kraken2.report"),
        out  =  protected(kraken_out + sample + ".kraken2.txt")
    params:
        database = config["Inputs"]["krakendb"],
        outdir = config["Outputs"]["kraken_out"]
    benchmark:
        benchmark_dir + "/" + sample +".kraken2.benchmark.txt"
    log:
        logs_dir + str(date) + ".illumina_kraken2.log"
    threads: config["Kraken2"]["kraken2_cores"]
    #Run interpretes the following block as python code, keep python synthax
    run:
        if reference_genome == None or reference_genome == "null":
            #Running kraken2 for paired data
            shell("kraken2 --threads {threads} --db {params.database}  --paired {input.read1} {input.read2} --use-names --report {output.report} --output {output.out}")

        else:
            #Running kraken2 for read extracted from BAM file
            shell("kraken2 --threads {threads} --db {params.database} {input.read1} --use-names --report {output.report} --output {output.out}")

rule Krona:
    input:
        report = rules.Kraken2.output.report
    output:
        krona_file = protected(krona_out + sample + ".krona"),
        krona_html = protected(krona_out + sample + ".html")
    benchmark:
        os.path.join(benchmark_dir, (str(date) + "_" + sample +".krona.benchmark.txt"))
    log:
        logs_dir + str(date) + ".krona.log"
    threads: config["Kraken2"]["kraken2_cores"]
    shell:
        "kreport2krona.py -r  {input.report} -o {output.krona_file}; ktImportText {output.krona_file} -o {output.krona_html}"

#Extraction of reads classified to a particular taxa id, by default the unclassified reads
rule extract_reads:
    input:
        read1 = read1_selected,
        read2 = read2_selected,
        report = rules.Kraken2.output.report,
        kraken = rules.Kraken2.output.out
    output:
        ext1 = protected(extracted_fa_out + sample + ".1.fastq.gz"),
        ext2 = protected(extracted_fa_out + sample + ".2.fastq.gz")
    params:
        intermediate1 = extracted_fa_out + sample + "1.fastq",
        intermediate2 = extracted_fa_out + sample + "2.fastq"
    benchmark:
        os.path.join(benchmark_dir, (str(date) + "_" + sample +".extracted_reads.benchmark.txt"))
    log:
        logs_dir + str(date) + ".extracted_reads.log"
    threads: 4
    shell:
        "extract_kraken_output_reads -k {input.kraken} -r {input.report} -s1 {input.read1} -s2 {input.read2} -t 0 --include-children -o {params.intermediate1} -o2 {params.intermediate2};"
        "cat {params.intermediate1} | pigz -p {threads} -c > {output.ext1}; cat {params.intermediate2} | pigz -p {threads} -c  > {output.ext2};"


#Optional task in the whole taxonomy rule in order to perform bracken estimation considering the kraken2 assignment

if kmer_dist != None or kmer_dist != "null":
    rule Bracken:
        input:
            report = rules.Kraken2.output.report
        output:
            abundance =  protected(kraken_out + sample + ".bracken_abundance.txt"),
            report =  protected(kraken_out + sample + ".kraken2_bracken_species.report")
        benchmark:
            os.path.join(benchmark_dir, (str(date) + sample +".bracken.benchmark.txt"))
        params:
            kmers=config["Inputs"]["kmer_dist"]
        log:
            logs_dir + str(date) + ".illumina_kraken.log"
        threads: config["Kraken2"]["kraken2_cores"]
        shell:
            "est_abundance.py -i {input.report} -k {params.kmers} -l S -t 10 -o {output.abundance} --out-report {output.report}; "