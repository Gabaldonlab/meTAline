Welcome to the meTAline wiki!

This resource is designed to guide you through the setup, execution, and optimization of meTAline, a Snakemake pipeline tailored for comprehensive metagenomic analysis.
Whether you're a researcher, bioinformatician, or enthusiast delving into metagenomics, this wiki provides detailed instructions and insights to facilitate your journey.

**About meTAline**

meTAline streamlines the preprocessing of short-read metagenomic data. Its workflow encompasses:

`Read Trimming and Filtering`: Enhances data quality by removing low-quality sequences and adapters. 

`Host Read Subtraction`: In case a reference genome is provided,it eliminates host-derived sequences to focus on microbial content.

`Taxonomic Classification`: Employs both k-mer and gene-marker approaches for microbial identification.

`Functional Profiling`: Analyzes the functional potential of microbial communities.

By automating these steps, meTAline ensures a consistent and reproducible analysis pipeline, saving time and reducing potential errors.

---

[**Quick Start**](https://github.com/Gabaldonlab/meTAline/wiki/Quick-Start-&-Try%E2%80%90Out)

To get started with meTAline, refer to our Quick Start & Try-Out guide.
This section provides step-by-step instructions to help you set up and execute the pipeline efficiently for a quick trial.

---

**Deployment Options**

meTAline offers flexibility in deployment to suit various computational environments:

**Docker**: *For containerized execution, learn how to Build and Run meTAline using Docker.*
  - [Build Docker image](https://github.com/Gabaldonlab/meTAline/wiki/Build-meTAline-docker-image)
  - [Run with Docker](https://github.com/Gabaldonlab/meTAline/wiki/Run-meTAline-with-docker)

**Singularity**: *Ideal for HPC clusters, follow the guides to Build and Run meTAline with Singularity.*
  - [Build Singularity image](https://github.com/Gabaldonlab/meTAline/wiki/Build-meTAline-singularity-image)
  - [Run with Singularity](https://github.com/Gabaldonlab/meTAline/wiki/Run-meTAline-with-singularity)

**High-Performance Computing (HPC)**: *For large-scale analyses, see our instructions on Running in HPC environments.*
  - [Running in HPC environment ](https://github.com/Gabaldonlab/meTAline/wiki/Run-meTAline-in-HPC-environment)

---

**Parallelization**

To expedite processing, meTAline supports parallel execution:

**HPC Cluster**: *Discover strategies for Parallelizing in HPC clusters.*
- [In HPC cluster](https://github.com/Gabaldonlab/meTAline/wiki/Parallelize-meTAline-in-HPC-environment)

**Local Machine**: *Learn how to Parallelize on a local machine to leverage multi-core processors.*
- [On local machine](https://github.com/Gabaldonlab/meTAline/wiki/Parallelize-meTAline-on-local-machine)

---

**Outputs and Target Rules**

Understand the expected results and how to customize them:

**Available Target Rules**: *Explore the Available target rules to tailor the pipeline to your specific needs.*
- [Available target rules](https://github.com/Gabaldonlab/meTAline/wiki/Available-target-rules)

---

**Benchmark and reporting**

After the completion of each meTAline run you will obtain a set of metrics per rule that you can include in an html report
- [Benchmark output](https://github.com/Gabaldonlab/meTAline/wiki/Benchmark-outputs)
- [Reporting](https://github.com/Gabaldonlab/meTAline/wiki/Reporting-meTAline)
---

[**Debugging**](https://github.com/Gabaldonlab/meTAline/wiki/Debugging-meTAline)

Want to contribute or encountered an issue? Our Debugging meTAline section offers solutions to common problems and tips for troubleshooting.

[**Citations and Acknowledgments**](https://github.com/Gabaldonlab/meTAline/wiki/Citations---Acknowledgments)

If meTAline contributes to your research, please refer to our Citations / Acknowledgments page for proper attribution and to recognize the efforts of the development team.

