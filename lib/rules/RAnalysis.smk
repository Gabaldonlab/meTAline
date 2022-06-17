#Author: Diego Fuentes
#Contact email: diego.fuentes@bsc.es 
#Barcelona
#Date:2022-05-26

#Rule to convert the bracken report to biom format. 
rule convert_biom:
    input:
        report= rules.Bracken.output.report
    output:
        biom_file = protected(kraken_out + sample + ".Braken.biom")
    benchmark:
        os.path.join(benchmark_dir, (str(date) + "_" + sample +".biom.benchmark.txt"))
    log:
        logs_dir + str(date) + ".biom.log"
    threads: config["Kraken2"]["kraken2_cores"]
    shell:
        "kraken-biom {input.report} -o {output.biom_file} --fmt json" 

#Rule to convert the biom format to a phyloseq object. 
rule Biom_to_Phyloseq:
    input:
        biom= rules.convert_biom.output.biom_file
    output:
        phyloseq_object =  ranalysis_out + sample + ".rds"
    benchmark:
        os.path.join(benchmark_dir, (str(date) + sample +".R_conversion.benchmark.txt"))
    log:
        logs_dir + str(date) + ".R_phyloseq.log"
    shell:
        r"""
        conda activate R_phyloS
        Rscript lib/scripts/Biom_to_Phyloseq.R {input.biom} {output.phyloseq_object}
        """

#Rule to calculate alpha diversity metrics
rule alpha_diversity:
    input:
        phylo_object = rules.Biom_to_Phyloseq.output.phyloseq_object
    output:
        rich_object =  ranalysis_out + sample + "_alpha_div.csv"
    benchmark:
        os.path.join(benchmark_dir, (str(date) + sample +".R_alpha_div.benchmark.txt"))
    log:
        logs_dir + str(date) + ".R_alpha_div.log"
    shell:
        r"""
        conda activate R_phyloS
        Rscript lib/scripts/alpha_diversity.R {input.phylo_object} {output.rich_object}
        """      

#Rule that plots a barplot of relative abundances of the phyla corresponding to the top 25 taxa in the sample  
rule bar_plot_toptaxa:
    input:
        phylo_object = rules.Biom_to_Phyloseq.output.phyloseq_object

    output:
        out_plot =  ranalysis_out + sample + "_Barplot_phyla.jpeg"

    benchmark:
        os.path.join(benchmark_dir, (str(date) + sample +".R_barplot.benchmark.txt"))
    log:
        logs_dir + str(date) + ".R_barplot.log"
    shell:
        r"""
        conda activate R_phyloS
        Rscript lib/scripts/bar_plot_toptaxa.R {input.phylo_object} {output.out_plot}
        """     

