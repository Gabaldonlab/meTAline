### Install Singularity on Debian based Linux. <a id="install-singularity" />
**IMPORTANT: If you are going to run the pipeline in an HPC environment, skip this section!**

*Note: you will always have to build the singularity images on your local machine, because it requires sudo and that you won't have in your remote hpc environment!*

To be able to build the Singularity image of MeTAline, you will need to install Singularity first:

1. Download from [HERE](https://github.com/sylabs/singularity/releases/tag/v4.1.5) the corresponding installation package (_.deb or _.rpm) to your operating system. (E.G.: Ubuntu 24.04 needs the "singularity-ce_4.1.5-noble_amd64.deb" file).
2. Install it from the downloaded package. Example command for the above mentioned version:

```bash
sudo apt install ./singularity-ce_4.1.5-noble_amd64.deb
```

3. Test if the installation is OK.

```bash
singularity --version
```

**NOTE:**
The sanbdox builds after version 4.1.0. has some bugs (like stripped '%runscript' section, not found '/bin/sh' symbolic links, etc..).
These kind of problems can be found in the Singularity-CE and Apptainer versions too!
E.G.: [runscript is being wiped out during build process #2561](https://github.com/apptainer/apptainer/issues/2561)

### Execute the image build and setup in cluster <a id="singularity-image-build" />

1. cd into the root directory of MeTAline (where the Dockerfile and the Makefile is)
2. Run the following command (**REQUIRES "sudo" PERMISSIONS! + This might take several minutes!**):
   **NOTE: This will install a SingularityCE provided tool with pip3 (spython) to translate the Dockerfile to Singularity!**

```bash
make singularity
```

Which translates to:
```bash
./_scripts_singularity/get_singularity_def_file_from_dockerfile.sh
sudo singularity build metaline.sif metaline-singularity.def
```

3. Copy the output **metaline.sif** file to your cluster with the **rsync** command. Example:

```bash
rsync ./metaline.sif <cluster_user>@<cluster_transfer_address>:</path/to/your/destination>
```
