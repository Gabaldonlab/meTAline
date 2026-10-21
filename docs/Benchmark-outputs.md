## Benchmark outputs <a id="benchmark-output" />

At the end of the procedure you might want to know the resources that were used. A folder called "Benchmark" will be created containing for each of some of the rules the following parameters:

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