#Author: Diego Fuentes, Dániel Majer and Olfat Khannous Lleiffe
#Contact email: olfat.khannous@bsc.es
#Barcelona
#Date:2026-02-03

###################
# Merging taxonomy assignment of different samples
###################


####################################################
# Rules that wait until all samples have completed #
####################################################
   
          
def checkpoint_code(tool):
    completed_tool_dir = os.path.join(completed_samples_folder, f"{tool}")
    completed_sample_dir = os.path.join(completed_tool_dir, f"{sample}")
    os.makedirs(completed_sample_dir, exist_ok=True)
    sample_completed_file = os.path.join(completed_sample_dir, f"{sample}_completed.txt")
    if not os.path.isfile(sample_completed_file):
      open(sample_completed_file, "a").close()

    # Counts the number of samples that have completed and have the completed txt file
    completed_samples = sum(1 for f in Path(completed_tool_dir).rglob('*completed.txt') if f.is_file())
     
    # Counts the number of overall samples in the joblist files
    joblist_files = glob.glob(f"{workingdir}/joblist*.txt")
    total_samples = 0
    for joblist_file in joblist_files:
      with open(joblist_file, "r") as f:
        num_lines = len(f.readlines())
        total_samples = total_samples+num_lines
        
    #The all_samples_completed.txt file should only be created if all samples have completed
    if total_samples == completed_samples:
      all_samples_completed_file = os.path.join(completed_tool_dir, f"{sample}/all_samples_completed.txt")
      if not os.path.isfile(all_samples_completed_file):
        open(all_samples_completed_file, "a").close()


# Checkpoint to determine if all samples have completed Kraken2
checkpoint completed_samples_trigger_kraken:
    input:
        os.path.join(kraken_out, f"{sample}.kraken2.report"),
    output:
        directory(os.path.join(completed_samples_folder, f"Kraken2/{sample}"))
    run:
      checkpoint_code("Kraken2")
     
# The merged files should only be created if all_samples_completed.txt file has been created
def maybe_merge_kraken(wildcards):
    ckpt = checkpoints.completed_samples_trigger_kraken.get(**wildcards)
    outdir = ckpt.output[0]

    trigger_file = os.path.join(outdir, f"all_samples_completed.txt")

    if os.path.exists(trigger_file):
        return [os.path.join(outdir, f"all_samples_completed.txt")]
    else:
        return []      
        
            
 # Checkpoint to determine if all samples have completed Metaphlan4
checkpoint completed_samples_trigger_metaphlan:
    input:
        os.path.join(config["Outputs"]["metaphlan4_out"], sample),
    output:
        directory(os.path.join(completed_samples_folder, f"Metaphlan4/{sample}"))
    run:
      checkpoint_code("Metaphlan4")

# The merged files should only be created if all_samples_completed.txt file has been created
def maybe_merge_metaphlan(wildcards):
    ckpt = checkpoints.completed_samples_trigger_metaphlan.get(**wildcards)
    outdir = ckpt.output[0]

    trigger_file = os.path.join(outdir, f"all_samples_completed.txt")

    if os.path.exists(trigger_file):
        return [os.path.join(outdir, f"all_samples_completed.txt")]
    else:
        return []
        

# Checkpoint to determine if all samples have completed Sylph
checkpoint completed_samples_trigger_sylph:
    input:
        os.path.join(config["Outputs"]["filtered_euk"], f"{sample}.tsv"),
        os.path.join(config["Outputs"]["filtered_bact"], f"{sample}.tsv")
    output:
        directory(os.path.join(completed_samples_folder, f"Sylph/{sample}"))
    run:
       checkpoint_code("Sylph")
            
# The merged files should only be created if all_samples_completed.txt file has been created
def maybe_merge_sylph(wildcards):
    ckpt = checkpoints.completed_samples_trigger_sylph.get(**wildcards)
    outdir = ckpt.output[0]

    trigger_file = os.path.join(outdir, f"all_samples_completed.txt")

    if os.path.exists(trigger_file):
        return [os.path.join(outdir, f"all_samples_completed.txt")]
    else:
        return []


def get_kraken_output(wildcards):
    # Run the checkpoint first then create an array of all the kraken reports
    ckpt = checkpoints.completed_samples_trigger_kraken.get(**wildcards)
    KRAKEN_OUTPUTS = []
    for file_path in glob.glob(os.path.join(workingdir, "configs/*.json")):
      with open(file_path, "r") as f:
       data = json.load(f)
       sample = data["Parameters"]["sample_barcode"]
       kraken_outputs = data["Outputs"]["kraken_out"]
       kraken_file_path = os.path.join(kraken_outputs, f"{sample}.kraken2.report")
       if os.path.exists(kraken_file_path):
         KRAKEN_OUTPUTS.append(kraken_file_path)
    return KRAKEN_OUTPUTS


def get_metaphlan_output(wildcards):
    # Run the checkpoint first then create an array of all the metaphlan reports
    ckpt = checkpoints.completed_samples_trigger_metaphlan.get(**wildcards)
    METAPHLAN4_OUTPUTS = []
    for file_path in glob.glob(os.path.join(workingdir, "configs/*.json")):
      with open(file_path, "r") as f:
       data = json.load(f)
       sample = data["Parameters"]["sample_barcode"]
       metaphlan4_outputs = data["Outputs"]["metaphlan4_out"]
       metaphlan4_file_path = os.path.join(metaphlan4_outputs, f"{sample}_profiled.txt")
       if os.path.exists(metaphlan4_file_path):
         METAPHLAN4_OUTPUTS.append(metaphlan4_file_path)
    return METAPHLAN4_OUTPUTS
    
def get_sylph_output_fungi(wildcards):
    # Run the checkpoint first then create an array of all the sylph fungi outputs
    ckpt = checkpoints.completed_samples_trigger_sylph.get(**wildcards)
    SYLPH_FUNGI_OUTPUTS = []
    for file_path in glob.glob(os.path.join(workingdir, "configs/*.json")):
      with open(file_path, "r") as f:
       data = json.load(f)
       sample = data["Parameters"]["sample_barcode"] 
       
       sylph_fungi_outputs = data["Outputs"]["filtered_euk"]
       sylph_fungi_file_path = os.path.join(sylph_fungi_outputs, f"{sample}.tsv")
       if os.path.exists(sylph_fungi_file_path):
         SYLPH_FUNGI_OUTPUTS.append(sylph_fungi_file_path)
         
    return SYLPH_FUNGI_OUTPUTS
    
    
def get_sylph_output_bacteria(wildcards):
    # Run the checkpoint first then create an array of all the sylph fungi outputs
    ckpt = checkpoints.completed_samples_trigger_sylph.get(**wildcards)
    SYLPH_BACTERIA_OUTPUTS = [] 
    for file_path in glob.glob(os.path.join(workingdir, "configs/*.json")):
      with open(file_path, "r") as f:
       data = json.load(f)
       sample = data["Parameters"]["sample_barcode"] 
             
       sylph_bacteria_outputs = data["Outputs"]["filtered_bact"]
       sylph_bacteria_file_path = os.path.join(sylph_bacteria_outputs, f"{sample}.tsv")
       if os.path.exists(sylph_bacteria_file_path):
         SYLPH_BACTERIA_OUTPUTS.append(sylph_bacteria_file_path)
         
    return SYLPH_BACTERIA_OUTPUTS
    
   
#Kraken2 reports to biom format that will be exported to phyloseq
rule merge_kraken2:
    input:
        merge = maybe_merge_kraken,
        report = get_kraken_output
    output:
        os.path.join(completed_samples_folder, f"{sample}_kraken_merged_flag.txt")
    benchmark:
        os.path.join(benchmark_dir, ("conversion_bio.benchmark.txt"))
    threads: config["Kraken2"]["kraken2_cores"]
    params:
        merged_biom_file = protected(os.path.join(merged_out, "Kraken2_merged.biom"))
    run:
        if not input.merge:
          print("merge didn't run")
        else:
          shell("kraken-biom {input.report} -o {params.merged_biom_file} --fmt json")
        open(output[0], "w").close()


#Metaphlan4 reports to txt table
rule merge_metaphlan4:
    input:
        merge = maybe_merge_metaphlan,
        metaphlan4_report = get_metaphlan_output
    output:
        os.path.join(completed_samples_folder, f"{sample}_metaphlan4_merged_flag.txt")
    benchmark:
        os.path.join(benchmark_dir, ("merging_metaphlan.benchmark.txt"))
    params:
        script = os.path.join(workflow.basedir, "lib/scripts/merge_metaphlan_tables.py"),
        metaphlan4_merged = protected(os.path.join(merged_out, "metaphlan_merged_abundance_table.txt"))
    run:
        if not input.merge:
          print("merge didn't run")
        else:
          shell("python3 {params.script} {input.metaphlan4_report} > {params.metaphlan4_merged}")
        open(output[0], "w").close()

#Sylph, merge the fungal identification by sylph to tsv table
rule sylph_merge:
    input:
        merge = maybe_merge_sylph,
        sylph_output_bacteria = get_sylph_output_fungi,
        sylph_output_fungi = get_sylph_output_bacteria
        
    output:
        os.path.join(completed_samples_folder, f"{sample}_sylph_merged_flag.txt")
    benchmark:
        os.path.join(benchmark_dir, ("merging_sylph.benchmark.txt"))
    params:
        script = os.path.join(workflow.basedir, "lib/scripts/merge_sylph_taxprof.py"),
        out_fungi = protected(os.path.join(merged_out, "sylph_fungi_merged.tsv")),
        out_bacteria = protected(os.path.join(merged_out, "sylph_bacteria_merged.tsv"))
    run:
        if not input.merge:
            print("merge didn't run")
        else:
          shell("""
            python3 {params.script} {input.sylph_output_fungi} --column sequence_abundance -o {params.out_fungi}
            python3 {params.script} {input.sylph_output_bacteria} --column sequence_abundance -o {params.out_bacteria}
          """)
        open(output[0], "w").close()

