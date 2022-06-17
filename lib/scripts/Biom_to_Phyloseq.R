#Loading packages that are needed for the following functions
library(phyloseq)
library(biomformat)

#Extraction of command line arguments, given by the snakemake rule 
args <- commandArgs(TRUE)
#First argument corresponds to the input
file <- args[1]
#Second argument corresponds to the output
out <- args[2]

#Read the report of the taxonomy assignment, in biom format 
biomfilename = read_biom(file)
#Conversion of the biom format into a phyloseq object
data <- import_biom(biomfilename, parseFunction=parse_taxonomy_default)
#Save the phyloseq object
saveRDS(data,out)

