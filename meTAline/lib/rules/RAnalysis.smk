#Author: Diego Fuentes, Dániel Majer and Olfat Khannous Lleiffe
#Contact email: olfat.khannous@bsc.es
#Barcelona
#Date:2025-03-27

#Rule to convert the Kraken2 report to biom format.
rule convert_biom:
    input:
        report= rules.Kraken2.output.report
    output:
        biom_file = protected(os.path.join(kraken_out, f"{sample}.Kraken2.biom"))

    threads: config["Kraken2"]["kraken2_cores"]
    shell:
        "kraken-biom {input.report} -o {output.biom_file} --fmt json"

#Rule to convert the biom format to a phyloseq object.
rule Biom_to_Phyloseq:
    input:
        biom= rules.convert_biom.output.biom_file
    output:
        phyloseq_object =  os.path.join(ranalysis_out, f"{sample}.rds")

    params:
        script = os.path.join(workflow.basedir, "lib/scripts/Biom_to_Phyloseq.R")

    shell:
        r"""
        Rscript {params.script} {input.biom} {output.phyloseq_object}
        """

#Rule to calculate alpha diversity metrics
rule alpha_diversity:
    input:
        phylo_object = rules.Biom_to_Phyloseq.output.phyloseq_object
    output:
        rich_object =  os.path.join(ranalysis_out, f"{sample}_alpha_div.csv")

    params:
        script = os.path.join(workflow.basedir, "lib/scripts/alpha_diversity.R")

    shell:
        r"""
        Rscript {params.script} {input.phylo_object} {output.rich_object}
        """

#Rule that plots a barplot of relative abundances of the phyla corresponding to the top 25 taxa in the sample
rule bar_plot_toptaxa:
    input:
        phylo_object = rules.Biom_to_Phyloseq.output.phyloseq_object

    output:
        out_plot =  os.path.join(ranalysis_out, f"{sample}_Barplot_phyla.jpeg")

    params:
        script = os.path.join(workflow.basedir, "lib/scripts/bar_plot_toptaxa.R")

    shell:
        r"""
        Rscript {params.script} {input.phylo_object} {output.out_plot}
        """

