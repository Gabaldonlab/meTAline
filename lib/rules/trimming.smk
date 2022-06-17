#Author: Diego Fuentes
#Contact email: diego.fuentes@bsc.es
#Barcelona
#Date:2021-01-29

###################
# TRIMMING AND QC #
###################

rule Trimmomatic:
    input:
        read1 = config["Inputs"]["reads_directory"] + "{file}_1."+config["Parameters"]["extension"],
        read2 = config["Inputs"]["reads_directory"] + "{file}_2."+config["Parameters"]["extension"]
    output:
        trim1 = protected(trimmomatic_out + "{file}_1_paired.fq.gz"),
        trim2 = protected(trimmomatic_out + "{file}_2_paired.fq.gz"),
        unpaired1 = temp(trimmomatic_out + "{file}_1_unpaired.fq.gz"),
        unpaired2 = temp(trimmomatic_out + "{file}_2_unpaired.fq.gz")
    params:
        illuminaclip = config["Trimmomatic"]["illuminaclip"],
        leading = config["Trimmomatic"]["leading"],
        trailing = config["Trimmomatic"]["trailing"],
        slidingwindow = config["Trimmomatic"]["slidingwindow"],
        minlen = config["Trimmomatic"]["minlen"]
    benchmark:
        os.path.join(benchmark_dir, (str(date) + "_" + sample +".{file}.trimming.benchmark.txt"))

    log: 
        logs_dir + str(date) + ".{file}.trimmomatic.log"
    threads: 
        config["Trimmomatic"]["trimmo_cores"]

    conda: os.path.join(workflow.basedir, "meTAline_env.yml")

    shell:
        "trimmomatic PE -threads {threads} -phred33 {input.read1} {input.read2} {output.trim1} {output.unpaired1} {output.trim2} {output.unpaired2} ILLUMINACLIP:{params.illuminaclip} LEADING:{params.leading} TRAILING:{params.trailing} SLIDINGWINDOW:{params.slidingwindow} MINLEN:{params.minlen} ;"

rule Concat_reads:
    input:
        trim1 = lambda wildcards: expand(rules.Trimmomatic.output.trim1, file=files.split(',')),
        trim2 = lambda wildcards: expand(rules.Trimmomatic.output.trim2,  file=files.split(','))
    output:
        concat1 = protected(trimmomatic_out + sample + ".1.fastq.gz"),
        concat2 = protected(trimmomatic_out + sample+ ".2.fastq.gz")
    benchmark:
        os.path.join(benchmark_dir, (str(date) + "_" + sample +".concat.benchmark.txt"))

    log:
        logs_dir + str(date) + ".concatenation.log"
    threads: 
        config["Trimmomatic"]["trimmo_cores"]

    conda: os.path.join(workflow.basedir, "meTAline_env.yml")

    shell:
        "zcat {input.trim1} | pigz -p {threads} -c  > {output.concat1}; zcat {input.trim2} | pigz -p {threads} -c  > {output.concat2}; "

rule fastqc:
    input:
        reads1 = rules.Concat_reads.output.concat1,
        reads2 = rules.Concat_reads.output.concat2
    output:
        DIR = directory(trimmomatic_out +sample+"_qc")
    params:
        mode = "--outdir"
    benchmark:
        os.path.join(benchmark_dir, (str(date) + "_" + sample +".fastqc.benchmark.txt"))
    log:
        logs_dir+str(date)+".fastqc.log"
    threads: 4
    conda: os.path.join(workflow.basedir, "meTAline_env.yml")
    shell:
        #Fastqc called from snakemake requires to have the outdir already created
        "mkdir {output.DIR}; fastqc {input.reads1} {input.reads2} {params.mode}={output.DIR} -t 4"
