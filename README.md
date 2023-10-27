# SPARCED

SPARCED is a simple and efficient pipeline for construction, merging, expansion, and simulation of large-scale, single-cell mechanistic models. With minimal set-up, users can configure the model for parallel runs on a Kubernetes cluster (SPARCED-nf), or small-scale experiments on their local machine (SPARCED-jupyter). More information on the model itself can be found [here](https://rdcu.be/cP6tK).


## Dependencies

- [Docker](https://docs.docker.com/get-docker/)
- [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html) *(Optional: SPARCED-nf requirement only)*
- [kubectl](https://kubernetes.io/docs/tasks/tools/install-kubectl/) *(Optional: SPARCED-nf requirement only)* 
- [Anaconda](https://www.anaconda.com/download) *(Optional: Linux only)*

## Instructions
### Setting up your environment

1. Clone this repository from the command-line using `git clone --recursive https://github.com/birtwistlelab/SPARCED.git`
2. Make sure the dependencies listed above are installed
3. (SPARCED-nf only) Ensure you have a Kubernetes `config` file for your chosen cluster located in your `~/.kube` folder
4. (Linux only) In the same directory step 1 was performed in, use the command `conda env create -f ./SPARCED/environment.yml`

### SPARCED-jupyter

#### Setup
Once you have Docker installed, follow the simple steps below.

1. In a terminal window, use the command `docker login` with your account credentials. If you don't have an account yet, head over to hub.docker.com to set one up.
2. Use `docker pull birtwistlelab/sparced-notebook:latest` to download the latest version of the docker image
3. Once the download is complete, use `docker run -p 8888:8888 --name testnb1 -i -t birtwistlelab/sparced-notebook:latest`, and in your browser, go to the last URL produced in your terminal from this command
    - `testnb1` is just a sample name for the container you're creating with this command. If you try to create another container with this name, delete the old one first with `docker rm testnb1`
4. Voila! Once this command finishes, you should see a URL in your terminal that looks similar to `http://127.0.0.1:8888/?token=4a7c71c7a3b0080a4f331b256ae435fbc70`--paste this into your browser. You can now begin stepping through the commands in each of the files in the `jupyter_notebooks` folder to learn more about the model and perform small runs.

To use custom data in the SPARCED-jupyter workflow, start a container with the commands above, and look at the `input_files` directory. Delete the current files and run `docker cp <datafile> <container-name(e.g. testnb1)>:/app/input_files/<datafile>` to replace them.

#### Demo

See the `Demo` folder for an example model compilation and simulation results.

### SPARCED-nf

#### Setting up the model

1. Edit the files in the `input_files` folder as needed. These values will be built into the *creation* of the model. For editing the values present for the model's *simulation*, see the directions accompanying the configuration step.
2. Use `./kube-runner/kube-load.sh <pvc-name> input_files` to load your input data to the PVC of the kube cluster. `kube-load.sh` assumes a `/workspaces` folder as the base of the PVC, and saves this input data at the path `/workspaces/$USER/input_files/`.  __Important__: Before every run where you plan on uploading new data to the PVC, it's important that you run `./kube-runner/kube-login.sh <pvc-name>` and delete the currently existing `/workspaces/$USER/input_files/` folder. If not, you run the risk of your model selecting the wrong files.
3. Edit the values in the `kube-nextflow.config` file in the `configs` folder (for help, see the config README [here](https://github.com/birtwistlelab/SPARCED/blob/master/configs/README.md))

#### Executing the workflow

1. Starting the workflow: 
    - `nextflow kuberun birtwistlelab/SPARCED -v <PVC-name> -c configs/kube-nextflow.config`
2. Retrieving data
    - After the run is finished, save your data from the PVC down to your local machine with `./kube-runner/kube-save.sh <pvc-name> <work-directory>` (`kube-save.sh` will find your `work directory` path as relative to your `workspace/$USER` directory in the PVC. So with the default configurations, it should just be `work`)

#### Final SPARCED-nf notes

- For problems, feel free to consult our [troubleshooting guide](https://github.com/birtwistlelab/SPARCED/blob/master/TROUBLESHOOTING.md) or create issue requests. 
- For brave souls with a lot of local computer power but no Kubernetes cluster, SPARCED-nf is also built to be able to run locally. Simply launch with `nextflow run` instead of `kuberun` and work from the `configs/local-nextflow.config` template.

### Converting SPARCED Docker container to Apptainer

Apptainer has built-in functionality for DockerHub container conversion. 
1. To build the apptainer container, in a terminal, enter the command `apptainer build --sandbox sparced/ docker://jonahhuggins/sparced:latest4`. 
    - Apptainer containers are read-only by default. The `--sandbox` flag enables read & write capabilities. 
2. Once the build is complete, use `apptainer shell --writable --pwd /SPARCED sparced/` to enter the container and work with SPARCED. 
    - the `--writable` flag enables writable execution of commands within the container, which is necessary for Message Passing Interface functionality.
    - The structure of the container binds the user's $HOME into the container by default. The `--pwd` flag allows for the starting directory to be specified by the user. If this flag is not used, one can find the SPARCED directory at root, /SPARCED. 



### SPARCED_Brep

Files to exactly replicate Bouhaddou2018 model. In this format, mRNAs and protein species are tracked in separate variables. The number of species and ratelaws are also the same Bouhaddou2018 model. SPARCED is the new, updated, and cleaned up version.

### Model testing and performance

The Docker image of the model and simulations are tested with multiple machines (see below). The simulation times are:

    - Model file creation time: ~2min
    - Model compilation time: 15-25min (depends on model size, i.e. numbers of species and observables)
    - Model simulation time: ~1min for 24-hour simulation

1. Ubuntu-Desktops:

    - Ubuntu 18.04, Intel Core i7 3930 CPU @ 3.20 GHz, 32 GB DDR3, Nvidia GTX 690 GPU

2. Windows-Desktops:

    - Windows 10 Education, Intel Core i5-3470 CPU @ 3.20 GHz, 8.00 GB RAM, Nvidia GTX 650 GPU, 64-bit operating system
    - Windows 10 Education, Intel Core i7-9700 CPU @ 5.00 GHz, 32.00 GB RAM, Nvidia RTX 2070 GPU, 64-bit operating system

3. Windows-Laptops:

    - Windows 10 Pro, Intel Core i7-8550U CPU @ 2.00 GHz, 16.00 GB RAM, Intel UHD 620 GPU, 64-bit operating system
    - Windows 10 Education, Intel Core i7-10705H CPU @ 2.60 GHz, 16.00 GB RAM, Nvidia GTX 1660ti GPU, 64-bit operating system

## Acknowledgements:

SPARCED is a product of [Birtwistle Lab](http://www.birtwistlelab.com/) and we greatly appreciate the help from multiple colloborators, including [Feltus Lab](https://www.clemson.edu/science/departments/genetics-biochemistry/people/profiles/ffeltus), [Hasenauer Lab](https://www.mathematics-and-life-sciences.uni-bonn.de/en/group-members/jan-hasenauer), and [Robert C. Blake](https://bbs.llnl.gov/RobertBlake.html) from LLNL.

## The Acronym:
The acronym SPARCED is composed of following elements, based on the sub-models in the large-scale mechanistic model.

### S: SBML
### P: Proliferation
### A: Apoptosis
### R: Receptor Signaling
### C: Cell Cycle
### E: Expression
### D: DNA Damage
