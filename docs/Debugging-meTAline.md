## With Docker

### VsCode

1. Install the Docker extension inside your editor
![image](https://github.com/user-attachments/assets/27b5820b-9b73-4854-b754-6a8843e56200)


2. Spin up an interactive container based on the cgenomics/metaline:latest Docker image
```sh
docker run -dit \
    --name metaline-debug \
    -v ./test_input:/input-directory \
    -v ./test_output:/output-directory \
    cgenomics/metaline:latest
```

3. Attach a VsCode instance to the running container
![image](https://github.com/user-attachments/assets/d1c83d87-241b-43a0-b5f5-f8ccd3332ba6)

4. Install the Python extensions inside the attach VsCode instance
![image](https://github.com/user-attachments/assets/f89063ad-b389-4590-8d35-6d644a9b8583)

5. Then you can set your breakpoint and debug as normally
![image](https://github.com/user-attachments/assets/45038888-ec12-45f1-a12a-39e1033ef9f1)


---

## With Singularity

### Setup Singularity image for debugging <a id="singularity-debugging" />
To debug the **meTAline** pipeline efficiently, use the `--bind` flag in Singularity to mount your modified source code into the container without rebuilding the image.

#### Accessing the meTAline Source Code in the image <a id="enter-into-image" />
The **meTAline** source code inside the Singularity image (`metaline.sif`) is located at `/meTAline`.

You can verify this with:
```bash
[user@host workspace]$ singularity run --cleanenv <path/to/metaline.sif> ls /meTAline
Helper_scripts_MN5  Illumina_MGI_adapters.fa  README.md  adapter_list_new.txt  conda_envs  debug.greasy.job  lib  meTAline.smk
[user@host workspace]$
```

#### Overwriting Source Code with the --bind flag <a id="overwrite-with-bind-flag" />
To test modifications to the pipeline, replace the source code in the image with your local version during execution.
Use the --bind flag to mount your local repository into the container, overwriting the existing directory.

##### Example <a id="example" />
For this example you have the current working directory inside the cloned meTAline repository, which looks like this:

*Note that the meTAline snakemake directory, ./meTAline, is right next to the built ./metaline.sif file!*

```bash
[user@host workspace]$ ls -lisah <path/to/cloned/meTAline/repository>
total 2.7G
75807561554 1.0K drwxrwxr-x  9 user bsc40 4.0K Dec 19 14:27 .
 6980952720 1.0K drwxr-sr-x 25 user bsc40 4.0K Dec 19 10:43 ..
75825726583 1.0K drwxrwxr-x  8 user bsc   4.0K Dec 11 11:22 .git
75825726584 1.0K drwxrwxr-x  3 user bsc   4.0K Dec 11 11:22 .github
75825726590  512 -rw-rw-r--  1 user bsc   3.2K Dec 11 11:22 .gitignore
75825726591  16K -rw-rw-r--  1 user bsc   5.8K Dec 11 11:22 Dockerfile
75825726592  48K -rw-rw-r--  1 user bsc    35K Dec 11 11:22 LICENSE
75825726595  512 -rw-rw-r--  1 user bsc    798 Dec 11 11:22 Makefile
75825726597  32K -rwxrwxr-x  1 user bsc    17K Dec 11 11:22 README.md
75825726599  512 -rw-rw-r--  1 user bsc   1.7K Dec 11 11:22 example_hpc_sbatch.job
75825726600  512 -rw-rw-r--  1 user bsc   1.7K Dec 11 11:22 example_job_def_mn5.job
75825726585 1.0K drwxrwxr-x  3 user bsc   4.0K Dec 11 11:22 external-sources
75825726586 1.0K drwxrwxr-x  5 user bsc   4.0K Dec 19 10:41 meTAline
75825726603  16K -rw-rw-r--  1 user bsc   5.2K Dec 11 11:54 metaline-singularity.def
75810056215 2.7G -rwxr-xr-x  1 user bsc   2.7G Dec 19 11:00 metaline.sif
75825726588 1.0K drwxrwxr-x  2 user bsc   4.0K Dec 11 11:22 scripts
75825726589 1.0K drwxrwxr-x  2 user bsc   4.0K Dec 11 11:22 test_input
[user@host workspace]$
```

Modify the command as follows:

**Original Command:**
```bash
singularity run --cleanenv metaline.sif \
    metaline \
    all \
    -j 16 \
    --configfile my_metaline_config.json
```

**Modified command:**
```bash
singularity run --cleanenv --bind ./meTAline:/meTAline metaline.sif \
    metaline \
    all \
    -j 16 \
    --configfile my_metaline_config.json
```

This allows you to test modifications to the pipeline without needing to rebuild or recompile the Singularity image, significantly speeding up development and debugging iterations.

#### Debugging tools <a id="debugging-tools" />
While IDE debuggers may not always work with Singularity images, most of the pipeline uses Python scripts, allowing you to utilize pdb (Python Debugger) for interactive debugging sessions.

To set a breakpoint in a Python script, use the **breakpoint()** function:

```bash
x = 10
y = 'Hi'
z = 'Hello'
print(y)

breakpoint()

print(z)
```

For further information about pdb refer to the [official Python 3 documentation](https://docs.python.org/3/library/pdb.html).

---
