#Loading packages that are needed for the following functions
library(phyloseq)
library(ggplot2)


#Extraction of command line arguments, given by the snakemake rule 
args <- commandArgs(TRUE)
#First argument corresponds to the input
input <- readRDS(args[1])
#Second argument corresponds to the output
out <- args[2]

## Indicate jpeg to save the resulting plot in the working dir
jpeg(file=out)
#First calculate the relative abundances
data_rel <- transform_sample_counts(input, function(OTU) (OTU/sum(OTU)*100))
#Which are the 25 most represented taxa?
top25_taxa <- names(sort(taxa_sums(data_rel), decreasing=TRUE))[1:25]
#Restrict the phyloseq object to only these taxa
data.top25 <- prune_taxa(top25_taxa, data_rel)
#Barplot at the phylum level (Rank2), considering the top25 taxa:
plot_bar(data.top25, fill="Rank2")  + 
  geom_bar(aes(color=Rank2, fill=Rank2), stat="identity", position="stack")
dev.off()