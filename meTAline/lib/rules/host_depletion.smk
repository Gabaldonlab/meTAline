#Author: Diego Fuentes, Dániel Majer and Olfat Khannous Lleiffe
#Contact email: olfat.khannous@bsc.es
#Barcelona
#Date:2025-03-27


###################
# ALIGNMENT RULES #
###################
if reference_genome != None or reference_genome != "null":
    rule hisat2:
        input:
            read1 = rules.Concat_reads.output.concat1,
            read2 = rules.Concat_reads.output.concat2
        output:
            out = protected(os.path.join(alignment_out, f"{sample}.hisat2.bam"))
        params:
            outdir = config["Outputs"]["alignment_out"],
            ref = reference_genome
        benchmark:
            os.path.join(benchmark_dir, (sample +".alignment.benchmark.txt"))
        threads: config["hisat2"]["hisat2_cores"]
        shell:
            #Important, check that the hisat2 index for your reference genome has been generated and provided to the config
            "hisat2 --threads {threads} -x {params.ref} -1 {input.read1} -2 {input.read2} | samtools sort -@ {threads} -O BAM -o {output.out};"

    rule index_bam:
        input:
            BAM = rules.hisat2.output.out
        output:
            BAI = protected(os.path.join(alignment_out, f"{sample}.hisat2.bam.bai"))

        threads:
            config["hisat2"]["hisat2_cores"]

        shell:
            "samtools index {input.BAM} {output.BAI};"

    rule extracted_unmapped_reads:
        input:
            BAM = rules.hisat2.output.out,
            BAI = rules.index_bam.output.BAI
        output:
            fastq = protected(os.path.join(alignment_out, f"{sample}.unmapped.fastq.gz")),
            intermediate=temp(os.path.join(alignment_out, "temp.bam"))
        benchmark:
            os.path.join(benchmark_dir, (sample +".extraction_unmapped_reads.benchmark.txt"))
        threads: config["hisat2"]["hisat2_cores"]
        shell:
            #Generate unmaped bam file, use it and then remove it
            "samtools view -b -f 4 {input.BAM} > {output.intermediate};"
            "samtools fastq {output.intermediate} | pigz -p {threads} -c  > {output.fastq};"
