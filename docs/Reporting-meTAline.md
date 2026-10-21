After the completion of a meTAline run it is possible to automatically generate an html report that includes the workflow topology, runtime statistics as well as the representation of the metrics obtained in the benchmark of each of the rules. In order to obtain it, you need to run, after the completion of the pipeline:

```bash
singularity run --cleanenv ./metaline.sif generate-benchmark-html \
    --output-html-path ./OUT_kraken2_BGUT/Benchmark/plots/report.html \ #This is the path were to store the report and the name to be given
    --config-file config.ERR1018192_BGUT.json #the config file you used for the run
```