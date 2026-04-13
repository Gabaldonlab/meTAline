#Author: Diego Fuentes, Dániel Majer and Olfat Khannous Lleiffe
#Contact email: olfat.khannous@bsc.es
#Barcelona
#Date:2026-02-03

###################
# Merging taxonomy assignment of different samples
###################

# extract samples names from the completed samples files and store in a list


####################################################
# Rules that wait until all samples have completed #
####################################################

# Checkpoint to determine if all samples have completed
checkpoint completed_samples_trigger:
    input:
        os.path.join(config["Outputs"]["metaphlan4_out"], sample),
        os.path.join(kraken_out, f"{sample}.kraken2.report"),
        os.path.join(config["Outputs"]["filtered_euk"], f"{sample}.tsv"),
        os.path.join(config["Outputs"]["filtered_bact"], f"{sample}.tsv")
    output:
        directory(os.path.join(workingdir, f"completed_samples/{sample}"))
    run:
        completed_samples_folder = os.path.join(workingdir, f"completed_samples")
        completed_sample_dir = os.path.join(completed_samples_folder, f"{sample}")
        os.makedirs(completed_samples_folder, exist_ok=True)
        os.makedirs(completed_sample_dir, exist_ok=True)
        sample_completed_file=os.path.join(completed_sample_dir, f"{sample}_completed.txt")
        if not os.path.isfile(sample_completed_file):
          open(sample_completed_file, "a").close()

        completed_samples=sum(1 for f in Path(completed_samples_folder).rglob('*completed.txt') if f.is_file())
     
        # Counts the number of overall samples
        joblist_files = glob.glob(f"{workingdir}/joblist*.txt")
        total_samples=0
        for joblist_file in joblist_files:
          with open(joblist_file, "r") as f:
            num_lines = len(f.readlines())
            total_samples=total_samples+num_lines
        
        if total_samples==completed_samples:
          all_samples_completed_file=os.path.join(workingdir, f"completed_samples/{sample}/all_samples_completed.txt")
          if not os.path.isfile(all_samples_completed_file):
            open(all_samples_completed_file, "a").close()

# Merges files depending if the all_samples_completed.txt file is present
def maybe_merge(wildcards):
    ckpt = checkpoints.completed_samples_trigger.get(**wildcards)
    outdir = ckpt.output[0]

    trigger_file = os.path.join(outdir, f"all_samples_completed.txt")

    if os.path.exists(trigger_file):
        return [os.path.join(outdir, f"all_samples_completed.txt")]
    else:
        return []

def get_kraken_output(wildcards):
    # Run the checkpoint first
    ckpt = checkpoints.completed_samples_trigger.get(**wildcards)
    KRAKEN_OUTPUTS = []
    for file_path in glob.glob(os.path.join(workingdir, "configs/*.json")):
      with open(file_path, "r") as f:
       data= json.load(f)
       sample= data["Parameters"]["sample_barcode"]
       kraken_outputs = data["Outputs"]["kraken_out"]
       kraken_file_path=os.path.join(kraken_outputs, f"{sample}.kraken2.report")
       if os.path.exists(kraken_file_path):
         KRAKEN_OUTPUTS.append(kraken_file_path)
    print("Kraken outputs are......")
    print(KRAKEN_OUTPUTS) 
    return KRAKEN_OUTPUTS


def get_metaphlan_output(wildcards):
    # Run the checkpoint first
    ckpt = checkpoints.completed_samples_trigger.get(**wildcards)
    METAPHLAN4_OUTPUTS=[]
    for file_path in glob.glob(os.path.join(workingdir, "configs/*.json")):
      with open(file_path, "r") as f:
       data= json.load(f)
       sample= data["Parameters"]["sample_barcode"] 
       metaphlan4_outputs = data["Outputs"]["metaphlan4_out"]
       metaphlan4_file_path= os.path.join(metaphlan4_outputs, f"{sample}_profiled.txt")
       if os.path.exists(metaphlan4_file_path):
         METAPHLAN4_OUTPUTS.append(metaphlan4_file_path)
    return METAPHLAN4_OUTPUTS
    
def get_sylph_fungi_output(wildcards):
    # Run the checkpoint first
    ckpt = checkpoints.completed_samples_trigger.get(**wildcards)
    SYLPH_FUNGI_OUTPUTS=[]
    for file_path in glob.glob(os.path.join(workingdir, "configs/*.json")):
      with open(file_path, "r") as f:
       data= json.load(f)
       sample= data["Parameters"]["sample_barcode"] 
       sylph_fungi_outputs = data["Outputs"]["filtered_euk"]
       sylph_fungi_file_path= os.path.join(sylph_fungi_outputs, f"{sample}.tsv")
       if os.path.exists(sylph_fungi_file_path):
         SYLPH_FUNGI_OUTPUTS.append(sylph_fungi_file_path)
       print("sylph outputs are......")
       print(SYLPH_FUNGI_OUTPUTS) 
    return SYLPH_FUNGI_OUTPUTS
    
    
def get_sylph_bacteria_output(wildcards):
    # Run the checkpoint first
    ckpt = checkpoints.completed_samples_trigger.get(**wildcards)   
    SYLPH_BACTERIA_OUTPUTS=[] 
    for file_path in glob.glob(os.path.join(workingdir, "configs/*.json")):
      with open(file_path, "r") as f:
       data= json.load(f)
       sample= data["Parameters"]["sample_barcode"] 
       sylph_bacteria_outputs = data["Outputs"]["filtered_bact"]
       sylph_bacteria_file_path= os.path.join(sylph_bacteria_outputs, f"{sample}.tsv")
       if os.path.exists(sylph_bacteria_file_path):
         SYLPH_BACTERIA_OUTPUTS.append(sylph_bacteria_file_path)
       print("sylph outputs are......")
       print(SYLPH_BACTERIA_OUTPUTS) 
    return SYLPH_BACTERIA_OUTPUTS
    

#Kraken2 reports to biom format that will be exported to phyloseq
rule merge_kraken2:
    input:
        merge = maybe_merge,
        report=get_kraken_output
    output:
        protected(os.path.join(merged_out, f"{sample}_kraken_merged_flag.txt"))
    benchmark:
        os.path.join(benchmark_dir, ("conversion_bio.benchmark.txt"))
    threads: config["Kraken2"]["kraken2_cores"]
    params:
        merged_biom_file = protected(os.path.join(merged_out, "merged.Kraken2.biom"))
    run:
        if not input.merge:
          print("merge didn't run")
        else:
          print("input report is....")
          print(input.report)
          shell("""
             echo "Processing report: {input.report}"
             echo "Processing params: {params.merged_biom_file}"
             kraken-biom {input.report} -o {params.merged_biom_file} --fmt json
          """)
        open(output[0], "w").close()


#Metaphlan4 reports to txt table
rule merge_metaphlan4:
    input:
        merge = maybe_merge,
        metaphlan4_report=get_metaphlan_output
    output:
        protected(os.path.join(merged_out, f"{sample}_metaphlan4_merged_flag.txt"))
    benchmark:
        os.path.join(benchmark_dir, ("merging_metaphlan.benchmark.txt"))
    params:
        script = os.path.join(workflow.basedir, "lib/scripts/merge_metaphlan_tables.py"),
        metaphlan4_merged = protected(os.path.join(merged_out, "merged_abundance_table.txt"))
    run:
        if not input.merge:
          print("merge didn't run")
        else:
          shell("""
            set -euo pipefail
            python3 {params.script} {input.metaphlan4_report} > {params.metaphlan4_merged}
            """)
        open(output[0], "w").close()

#Sylph, merge the fungal identification by sylph to tsv table
rule sylph_fungi_merge:
    input:
        merge = maybe_merge,
        fungi_identified = get_sylph_fungi_output
    output:
        protected(os.path.join(merged_out, f"{sample}_sylph_fungi_merged_flag.txt"))
    benchmark:
        os.path.join(benchmark_dir, ("merging_sylph_fungi.benchmark.txt"))
    params:
        script = os.path.join(workflow.basedir, "lib/scripts/merge_sylph_taxprof.py"),
        out_fungi = protected(os.path.join(merged_out, "fungi.tsv"))
    run:
        if not input.merge:
            print("merge didn't run")
        else:
          print("input report is....")
          print(input.fungi_identified)
          shell("""
            python3 {params.script} {input.fungi_identified} --column sequence_abundance -o {params.out_fungi}
            """)
        open(output[0], "w").close()

#Sylph, merge the bacterial identification by sylph to tsv table
rule sylph_bacteria_merge:
   input:
       merge = maybe_merge,
       bacteria_identified = get_sylph_bacteria_output
   output:
        protected(os.path.join(merged_out, f"{sample}_sylph_bacteria_merged_flag.txt"))
   benchmark:
       os.path.join(benchmark_dir, ("merging_sylph_bacteria.benchmark.txt"))
   params:
        script = os.path.join(workflow.basedir, "lib/scripts/merge_sylph_taxprof.py"),
        out_bacteria= protected(os.path.join(merged_out, "bacteria.tsv"))
   run:
       if not input.merge:
          print("merge didn't run")
       else:
          print("input report is....")
          print(input.bacteria_identified)
          shell("""
             python3 {params.script} {input.bacteria_identified} --column sequence_abundance -o {params.out_bacteria}
          """)    
       open(output[0], "w").close()   
