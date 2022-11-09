#Author: Diego Fuentes
#Contact email: diego.fuentes@bsc.es
#Barcelona
#Date:2022-10-27


###################
# ALIGNMENT RULES #
###################
if reference_genome != None or reference_genome != "null":
    rule hisat2:
        input:
            read1 = rules.Concat_reads.output.concat1,
            read2 = rules.Concat_reads.output.concat2
        output:
            out = protected(alignment_out + sample +".hisat2.bam")
        params:
            outdir = config["Outputs"]["alignment_out"],
            ref = reference_genome

        log:
            logs_dir + str(date) + "_" + sample +".alignment.log"
        benchmark:
            os.path.join(benchmark_dir, (str(date) + "_" + sample +".alignment.benchmark.txt"))
        threads: config["hisat2"]["hisat2_cores"]
        conda: os.path.join(workflow.basedir, "meTAline_env.yml")
        shell:
            #Important, check that the hisat2 index for your reference genome has been generated and provided to the config
            "hisat2 --threads {threads} -x {params.ref} -1 {input.read1} -2 {input.read2} | samtools sort -@ {threads} -O BAM -o {output.out} 2> {log};"

    rule index_bam:
        input:
            BAM = rules.hisat2.output.out
        output:
            BAI = protected(alignment_out + sample +".hisat2.bam.bai")
        log:
            logs_dir + str(date) + "_"+sample+".index_bam.log"  
        threads:
            config["hisat2"]["hisat2_cores"]

        conda: os.path.join(workflow.basedir, "meTAline_env.yml")

        shell:
            "samtools index {input.BAM} {output.BAI} 2> {log}"

    rule extracted_unmapped_reads:
        input:
            BAM = rules.hisat2.output.out,
            BAI = rules.index_bam.output.BAI
        output:
            fasta = protected(alignment_out + sample +".unmapped.fastq.gz"),
            intermediate=temp(alignment_out + "temp.bam")
        log: logs_dir + str(date) + "_"+sample+".bam2fasta.log"  
        threads: config["hisat2"]["hisat2_cores"]
        conda: os.path.join(workflow.basedir, "meTAline_env.yml")
        shell:
            #Generate unmaped bam file, use it and then remove it
            "samtools view -b -f 4 {input.BAM} > {output.intermediate};"
            "samtools fasta {output.intermediate} | pigz -p {threads} -c  > {output.fasta};"
