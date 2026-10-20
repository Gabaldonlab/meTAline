#Author: Laura McCluskey, Diego Fuentes, Dániel Majer and Olfat Khannous Lleiffe
#Contact email: olfat.khannous@bsc.es
#Barcelona
#Date:2026-05-12


###################
# Removal of Host and Phix reads, in case that these are given by the users as reference genomes
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
            BAI = temp(os.path.join(alignment_out, f"{sample}.hisat2.bam.bai"))
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
            #Generate unmapped bam file, use it and then remove it
            "samtools view -b -f 4 {input.BAM} > {output.intermediate};"
            "samtools fastq {output.intermediate} | pigz -p {threads} -c  > {output.fastq};"

    #If the reference genome of the PhiX control is given, add an extra step of alignment to the reference genome to remove the corresponding reads.
    if reference_phix != None or reference_phix != "null":

        rule hisat2_phix:
                input:
                    reads = rules.extracted_unmapped_reads.output.fastq
                output:
                    bam_phix = protected(os.path.join(alignment_out, f"{sample}.phix.hisat2.bam"))
                params:
                    ref_phix = reference_phix
                benchmark:
                    os.path.join(benchmark_dir, (sample + ".phix_alignment.benchmark.txt"))
                threads:
                    config["hisat2"]["hisat2_cores"]
                shell:
                    """
                    hisat2 \
                        --threads {threads} \
                        -x {params.ref_phix} \
                        -U {input.reads} | \
                    samtools sort \
                        -@ {threads} \
                        -O BAM \
                        -o {output.bam_phix}
                    """

        rule index_phix_bam:
                input:
                    BAM_phix = rules.hisat2_phix.output.bam_phix
                output:
                    BAI_phix = protected(os.path.join(alignment_out, f"{sample}.phix.hisat2.bam.bai"))
                threads:
                    config["hisat2"]["hisat2_cores"]
                shell:
                    """
                    samtools index {input.BAM_phix} {output.BAI_phix}
                    """

        rule extract_phix_unmapped_reads:
                input:
                    BAM_phix = rules.hisat2_phix.output.bam_phix,
                    BAI_phix = rules.index_phix_bam.output.BAI_phix
                output:
                    fastq_phix = protected(os.path.join(alignment_out, f"{sample}.clean.fastq.gz")),
                    intermediate_phix = temp(os.path.join(alignment_out, f"{sample}.phix_unmapped.bam"))
                benchmark:
                    os.path.join(benchmark_dir, (sample + ".phix_unmapped_extraction.benchmark.txt"))
                threads:
                    config["hisat2"]["hisat2_cores"]
                shell:
                    """
                    samtools view -b -f 4 {input.BAM_phix} > {output.intermediate_phix}

                    samtools fastq {output.intermediate_phix} | \
                    pigz -p {threads} -c > {output.fastq_phix}
                    """
