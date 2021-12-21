#Author: Diego Fuentes
#Contact email: diego.fuentes@bsc.es
#Barcelona
#Date:2021-01-29

######################
# TAXONOMY ABUNDANCE #
######################
#Check whether we use the filtered reads or the unmapped reads when the host genome is provided
if reference_genome == None or reference_genome == "null":
    read1_selected = rules.Concat_reads.output.concat1,
    read2_selected = rules.Concat_reads.output.concat2
else:
    read1_selected = rules.extracted_unmapped_reads.output.fasta,
    read2_selected = rules.extracted_unmapped_reads.output.fasta

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
    conda: os.path.join(workflow.basedir, "WGS_env.yml")
    #Run interpretes the following block as python code, keep python synthax
    run:
        if reference_genome == None or reference_genome == "null":
            #Running kraken2 for paired data
            shell("kraken2 --threads {threads} --db {params.database}  --paired {input.read1} {input.read2} --use-names --report {output.report} --output {output.out}")

        else:
            #Running kraken2 for read extracted from BAM file
            shell("kraken2 --threads {threads} --db {params.database} {input.read1} --use-names --report {output.report} --output {output.out}")

rule Bracken:
    input:
        report = rules.Kraken2.output.report
    output:
        abundance =  protected(kraken_out + sample + ".bracken_abundance.txt")
    benchmark:
        os.path.join(benchmark_dir, (str(date) + sample +".bracken.benchmark.txt"))
    params:
        kmers=config["Inputs"]["kmer_dist"]
    log:
        logs_dir + str(date) + ".illumina_kraken.log"
    threads: config["Kraken2"]["kraken2_cores"]
    conda: os.path.join(workflow.basedir, "WGS_env.yml")
    shell:
        "est_abundance.py -i {input.report} -k {params.kmers} -l S -t 10 -o {output.abundance}; "

rule Krona:
    input:
        report = rules.Kraken2.output.report
    output:
        krona_file = protected(krona_out + sample + ".krona"),
        krona_html = protected(krona_out + sample + ".hmtl")
    benchmark:
        os.path.join(benchmark_dir, (str(date) + "_" + sample +".krona.benchmark.txt"))
    log:
        logs_dir + str(date) + ".krona.log"
    threads: config["Kraken2"]["kraken2_cores"]
    conda: 2
    shell:
        "kreport2krona.py -r  {input.report} -o {output.krona_file}; ktImportText {output.krona_file} -o {output.krona_html}"

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
    conda: os.path.join(workflow.basedir, "WGS_env.yml")
    shell:
        "extract_kraken_reads.py -k {input.kraken} -r {input.report} -s1 {input.read1} -s2 {input.read2}  -t 2 2157 --exclude --include-children  -o {params.intermediate1}  -o2 {params.intermediate2};"
        "cat {params.intermediate1} | pigz -p {threads} -c > {output.ext1}; cat {params.intermediate2} | pigz -p {threads} -c  > {output.ext2};"
        "rm {params.intermediate1}; rm {params.intermediate2};"