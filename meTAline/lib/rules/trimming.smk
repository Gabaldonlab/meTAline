#Author: Diego Fuentes, Dániel Majer and Olfat Khannous Lleiffe
#Contact email: olfat.khannous@bsc.es
#Barcelona
#Date:2025-03-27

###################
# TRIMMING AND QC #
###################

#Trimming of raw data
rule Trimmomatic:
    input:
        read1 = os.path.join(config["Inputs"]["reads_directory"], f"{file}{config['Parameters']['prefix_fr']}1.{config['Parameters']['extension']}"),
        read2 = os.path.join(config["Inputs"]["reads_directory"], f"{file}{config['Parameters']['prefix_fr']}2.{config['Parameters']['extension']})"
    output:
        trim1 = protected(os.path.join(trimmomatic_out, "{file}_1_paired.fq.gz")),
        trim2 = protected(os.path.join(trimmomatic_out, "{file}_2_paired.fq.gz")),
        unpaired1 = protected(os.path.join(trimmomatic_out, "{file}_1_unpaired.fq.gz")),
        unpaired2 = protected(os.path.join(trimmomatic_out, "{file}_2_unpaired.fq.gz"))
    params:
        illuminaclip = os.path.join(workflow.basedir, config["Trimmomatic"]["illuminaclip"]),
        leading = config["Trimmomatic"]["leading"],
        trailing = config["Trimmomatic"]["trailing"],
        slidingwindow = config["Trimmomatic"]["slidingwindow"],
        minlen = config["Trimmomatic"]["minlen"]
    threads:
        config["Trimmomatic"]["trimmo_cores"]


    shell:
        "trimmomatic PE -threads {threads} -phred33 {input.read1} {input.read2} {output.trim1} {output.unpaired1} {output.trim2} {output.unpaired2} ILLUMINACLIP:{params.illuminaclip} LEADING:{params.leading} TRAILING:{params.trailing} SLIDINGWINDOW:{params.slidingwindow} MINLEN:{params.minlen} ;"

#Concatenation of reads in case the samples are sequenced in different sequencing lanes
rule Concat_reads:
    input:
        trim1 = lambda wildcards: expand(rules.Trimmomatic.output.trim1, file=files.split(',')),
        trim2 = lambda wildcards: expand(rules.Trimmomatic.output.trim2,  file=files.split(','))
    output:
        concat1 = protected(os.path.join(trimmomatic_out, f"{sample}.1.fastq.gz")),
        concat2 = protected(os.path.join(trimmomatic_out, f"{sample}.2.fastq.gz"))

    threads:
        config["Trimmomatic"]["trimmo_cores"]


    shell:
        "zcat {input.trim1} | pigz -p {threads} -c  > {output.concat1}; zcat {input.trim2} | pigz -p {threads} -c  > {output.concat2}; "

#Quality assessment of the reads, using fastqc
rule fastqc:
    input:
        reads1 = rules.Concat_reads.output.concat1,
        reads2 = rules.Concat_reads.output.concat2
    output:
        DIR = directory(os.path.join(trimmomatic_out, f"{sample}_qc"))
    params:
        mode = "--outdir",
        adapters = os.path.join(workflow.basedir, config["fastqc"]["adapters"])
    threads: 4
    shell:
        #Fastqc called from snakemake requires to have the outdir already created. The (-a) is the list of adapters to be checked.
        "mkdir {output.DIR}; fastqc {input.reads1} {input.reads2} {params.mode}={output.DIR} -t 4 -a {params.adapters}; "
