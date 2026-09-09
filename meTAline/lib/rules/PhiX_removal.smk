#Author: Laura McCluskey, Diego Fuentes, Dániel Majer and Olfat Khannous Lleiffe
#Contact email: olfat.khannous@bsc.es
#Barcelona
#Date:2026-05-12


###################
# Removal of PhiX reads if host genome was not given (e.g enviromental sample):
################### 1
if reference_phix not in [None, "null"] and reference_genome in [None, "null"]:

    rule remove_phix_reads:
        input:
            read1 = rules.Concat_reads.output.concat1,
            read2 = rules.Concat_reads.output.concat2
        output:
            out_bam_phix = temp(os.path.join(alignment_out, f"{sample}.hisat2_phix.bam")),
            BAI_phix = temp(os.path.join(alignment_out, f"{sample}.hisat2_phix.bam.bai")),
            intermediate_unmapped_phix = temp(os.path.join(alignment_out, "temp_phix.bam")),
            fastq_phix=protected(os.path.join(alignment_out, f"{sample}.unmapped.fastq.gz")),
        params:
            outdir = config["Outputs"]["alignment_out"],
            ref_phix = reference_phix
        threads: 8
        benchmark:
            os.path.join(benchmark_dir, (sample +".alignment_phix.benchmark.txt"))
        shell:
           r"""
           hisat2 --threads {threads} -x {params.ref_phix} -1 {input.read1} -2 {input.read2} | samtools sort -@ {threads} -O BAM -o {output.out_bam_phix}
           samtools index {output.out_bam_phix} {output.BAI_phix}
           samtools view -b -f 4 {output.out_bam_phix} > {output.intermediate_unmapped_phix}
           samtools fastq {output.intermediate_unmapped_phix} | pigz -p {threads} -c  > {output.fastq_phix}
           """
