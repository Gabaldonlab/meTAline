## Benchmark outputs and reporting meTAline <a id="benchmark-output" />

After the completion of a meTAline run, you might want to know the resources that were used. it is possible to automatically generate an html report that includes the workflow topology, runtime statistics as well as the representation of the metrics obtained in the benchmark of each of the rules. 

In order to obtain it, you need to run, after the completion of the pipeline:

```bash
singularity run --cleanenv ./metaline.sif generate-benchmark-html \
    --output-html-path ./OUT_kraken2_BGUT/Benchmark/plots/report.html \ #This is the path were to store the report and the name to be given
    --config-file config.ERR1018192_BGUT.json #the config file you used for the run
```

For this command to work fully, it is essential that there is a required.txt file in the Benchmark folder. Without this file, the Results section will not be included in the html report.

In addition to the html report, the command above will also create files in the benchmark folder containing information on the following parameters:

`s` --> Running time in seconds

`h: m: s` --> Running time in hours, minutes and seconds

`max_rss` --> "Maximum Resident Set Size". Total amount of physical memory.

`max_vms` --> "Maximum Virtual Memory Size". Total amount of virtual memory.

`max_uss` --> "Unique Set Size". Memory that is unique to a process.

`max_pss` --> "Proportional Set Size". Amount of memory shared with other processes.

`io_in` --> MB read.

`io_out` --> MB written.

`mean_load` --> CPU usage over time / Total running time.

---
