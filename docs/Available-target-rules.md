## Available target rules <a id="target-rules" />

The target rules currently available to use are:

-   rule **all**: From trimming and filtering of the reads, to host removal (if provided the indexed host genome to make an alignment of the reads) to taxonomy assignment based on both a k-mer and gene marker approach. Functional prediction based on Humann. It also includes other functionalities such as representation of the taxonomy assignment by krona, extraction of desired reads, etc.

- rule **trimming**: Trimming of the reads, quality assessment and concatenation (useful when having samples that are sequenced in different sequencing lanes).
- rule **host_depletion**: This rule is to take trimmed reads and align them to an indexed reference genome, for host substraction.
- rule **kmer_taxonomy**: This rule is used if you already have the concatenated/unmapped reads and you want to perform the taxonomic assignment based on k-mers, without running the initial steps of the pipeline.
- rule **RAnalysis**: This rule is used to make a basic R analysis and plotting alpha diversity.
- rule **BioBakery**: This rule is to perform the taxonomy and functional profiling using Biobakery tools based on gene markers: Metaphlan4 and Humann.

---