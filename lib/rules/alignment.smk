#Author: Diego Fuentes
#Contact email: diego.fuentes@bsc.es
#Barcelona
#Date:2021-01-29


###################
# ALIGNMENT RULES #
###################
if reference_genome != None or reference_genome != "null":
    rule bowtie2:
        input:
            read1 = rules.Concat_reads.output.concat1,
            read2 = rules.Concat_reads.output.concat2
        output:
            out = protected(alignment_out + sample +".bowtie2.bam")
        params:
            outdir = config["Outputs"]["alignment_out"],
            D = config["Bowtie2"]["bowtie2_D_param"],
            R = config["Bowtie2"]["bowtie2_R_param"],
            N = config["Bowtie2"]["bowtie2_N_param"],
            L = config["Bowtie2"]["bowtie2_L_param"],
            i = config["Bowtie2"]["bowtie2_i_param"],
            score_min = config["Bowtie2"]["bowtie2_score_min"],
            ref = reference_genome

        log:
            logs_dir + str(date) + "_" + sample +".alignment.log"
        benchmark:
            os.path.join(benchmark_dir, (str(date) + "_" + sample +".alignment.benchmark.txt"))
        threads: config["Bowtie2"]["bowtie2_cores"]
        conda: os.path.join(workflow.basedir, "WGS_env.yml")
        shell:
            #Important, check that the bowtie2 index for your reference genome has been generated and provided to the config
            "bowtie2-align-s --local --threads {threads} -x {params.ref} -1 {input.read1} -2 {input.read2} -D {params.D} -R {params.R} -N {params.N} -L {params.L} -i {params.i} --score-min {params.score_min} | samtools sort -@ {threads} -O BAM -o {output.out} 2> {log};"

    rule index_bam:
        input:
            BAM = rules.bowtie2.output.out
        output:
            BAI = protected(alignment_out + sample +".bowtie2.bam.bai")
        log:
            logs_dir + str(date) + "_"+sample+".index_bam.log"  
        threads:
            config["Bowtie2"]["bowtie2_cores"]

        conda: os.path.join(workflow.basedir, "WGS_env.yml")

        shell:
            "samtools index {input.BAM} {output.BAI} 2> {log}"

    rule extracted_unmapped_reads:
        input:
            BAM = rules.bowtie2.output.out,
            BAI = rules.index_bam.output.BAI
        output:
            fasta = protected(alignment_out + sample +".unmapped.fastq.gz"),
            intermediate=temp(alignment_out + "temp.bam")
        log: logs_dir + str(date) + "_"+sample+".bam2fasta.log"  
        threads: config["Bowtie2"]["bowtie2_cores"]
        conda: os.path.join(workflow.basedir, "WGS_env.yml")
        shell:
            #Generate unmaped bam file, use it and then remove it
            "samtools view -b -f 4 {input.BAM} > {output.intermediate};"
            "samtools fasta {output.intermediate} | pigz -p {threads} -c  > {output.fasta};"
