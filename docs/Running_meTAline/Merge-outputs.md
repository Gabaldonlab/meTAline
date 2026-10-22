## Merge outputs <a id="Merge outputs" />

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
