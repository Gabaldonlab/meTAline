#Loading packages that are needed for the following functions
library(phyloseq)

#Extraction of command line arguments, given by the snakemake rule 
args <- commandArgs(TRUE)
#First argument corresponds to the input
input <- readRDS(args[1])
#Second argument corresponds to the output
out <- args[2]

#Calculate alpha diversity metrics
rich = estimate_richness(input, measures = c("Observed","Shannon", "Simpson"))
#Write the calculated metrics in an excel file
write.csv(rich, out)