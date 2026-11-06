## Merge outputs <a id="Merge outputs" />

MeTAline is designed as a single-input, single-output workflow. However, as you may have discovered while navigating our wiki and manuscript, its highly parallelizable architecture makes it particularly well suited for efficiently preprocessing large-scale studies.

Since the output generated for each sample is stored in a separate folder, we provide a set of helper commands below to facilitate further analysis of the dataset as a whole.

There are two ways to use these commands: you can either collect all the individual reports into a single folder and run the commands on that folder, or provide a list of report files directly.

Example of command to find reports and copy them in a folder: 
```sh
find /PROJECT1 -type f -name "*kraken2.report" -exec cp {} /Kraken2_reports \;
```

The commands to be used for merging are (You need to have installed python3, kraken-biom and HUMAnN)

### Merge Sylph outputs

```sh
python3 lib/scripts/merge_sylph_taxprof.py {list of sylph output fungi tsv files} --column sequence_abundance -o sylph_fungi_merged.tsv
python3 lib/scripts/merge_sylph_taxprof.py {list of sylph output bacteria tsv files} --column sequence_abundance -o sylph_bacteria_merged.tsv
```
---

### Merge Kraken outputs
```sh
kraken-biom {list_of_kraken_reports} -o Kraken2_merged.biom --fmt json
```
---

### Merge Metaphlan outputs
```sh
python3 lib/scripts/merge_sylph_taxprof.py {list of metaphalan output profiled.txt files} > metaphlan_merged_abundance_table.txt
```
---

### Merge HUMAnN outputs

```sh
humann_join_tables -i . --file_name pathabundance.tsv --output merged_pathabundance.tsv
```