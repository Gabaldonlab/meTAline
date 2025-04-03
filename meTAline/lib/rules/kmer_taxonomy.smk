#Author: Diego Fuentes, Dániel Majer and Olfat Khannous Lleiffe
#Contact email: olfat.khannous@bsc.es
#Barcelona
#Date:2025-03-27

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

#Taxonomy assignment with kraken2 (k-mer based approach)
rule Kraken2:
    input:
        read1 = read1_selected,
        read2 = read2_selected
    output:
        report = protected(os.path.join(kraken_out, f"{sample}.kraken2.report")),
        out  =  protected(os.path.join(kraken_out, f"{sample}.kraken2.txt"))
    params:
        database = config["Inputs"]["krakendb"],
        outdir = config["Outputs"]["kraken_out"]

    threads: config["Kraken2"]["kraken2_cores"]
    #Run interpretes the following block as python code, keep python synthax
    run:
        if reference_genome == None or reference_genome == "null":
            #Running kraken2 for paired data
            shell("kraken2 --threads {threads} --db {params.database}  --paired {input.read1} {input.read2} --use-names --report {output.report} --output {output.out}")

        else:
            #Running kraken2 for read extracted from BAM file
            shell("kraken2 --threads {threads} --db {params.database} {input.read1} --use-names --report {output.report} --output {output.out}")

#Krona representation of the kraken2 assignment 
rule Krona:
    input:
        report = rules.Kraken2.output.report
    output:
        krona_file = protected(os.path.join(krona_out, f"{sample}.krona")),
        krona_html = protected(os.path.join(krona_out, f"{sample}.html"))
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
        ext1 = protected(os.path.join(extracted_fa_out, f"{sample}.1.fastq.gz")),
        ext2 = protected(os.path.join(extracted_fa_out, f"{sample}.2.fastq.gz"))
    params:
        intermediate1 = os.path.join(extracted_fa_out, f"{sample}1.fastq"),
        intermediate2 = os.path.join(extracted_fa_out, f"{sample}2.fastq"),
        taxid = config["Inputs"]["taxid"]
    threads: 4
    shell:
        "extract_kraken_output_reads -k {input.kraken} -r {input.report} -s1 {input.read1} -s2 {input.read2} -t {params.taxid} --include-children -o {params.intermediate1} -o2 {params.intermediate2};"
        "cat {params.intermediate1} | pigz -p {threads} -c > {output.ext1}; cat {params.intermediate2} | pigz -p {threads} -c  > {output.ext2};"
