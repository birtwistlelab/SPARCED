# SPARCED-nf: A Nextflow Pipeline for SPARCED

SPARCED-nf is a Nextflow pipeline designed to be a more scalable and user-friendly version of the SPARCED model (previously the mechanistic pan-cancer signaling model) by the Birtwistle Lab. With minimal set-up, a user can configure the model for high-intensity runs on a Kubernetes cluster, or small-scale experiments on their local machine. More information on the model itself can be found [here](https://github.com/birtwistlelab/SPARCED).


## Dependencies

- [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html)
- [Docker](https://docs.docker.com/get-docker/)
- (Optional) [kubectl](https://kubernetes.io/docs/tasks/tools/install-kubectl/)

## Walkthrough/Instructions
### Setting up your environment

1. Clone this repository from the command-line using `git clone --recursive https://github.com/ebenz99/SPARCED-nf.git`
2. Make sure the dependencies listed above are installed
3. (Kubernetes only) Ensure you have a Kubernetes `config` file for your chosen cluster located in your `~/.kube` folder

### Setting up the model

1. Edit the files in the `input_data` folder as needed. These values will be built into the *creation* of the model. For editing the values present for the model's *simulation*, see the directions accompanying the configuration step.
2. (Kubernetes only) Use `./kube-runner/kube-load.sh <pvc-name> input_data` to load your input data to the PVC of the kube cluster. `kube-load.sh` assumes a `/workspaces` folder as the base of the PVC, and saves this input data at the path `/workspaces/$USER/input_data/`.  __Important__: If you do this more than once, `kube-login` into the cluster and delete the original copy--otherwise it will not be overwritten.
3. Edit the values in either `kube-nextflow.config` or `local-nextflow.config` section of the `configs` folder (for help, see the config README [here](https://github.com/ebenz99/SPARCED-nf/blob/master/configs/README.md))
4. (Kubernetes only) Use `./kube-runner/kube-load.sh <pvc-name> configs` in the same way you did earlier to move your configuration files to the PVC. __Important__: If you do this more than once, `kube-login` into the cluster and delete the original copy--otherwise it will not be overwritten.

### Running the workflow

1. To start the workflow:
    - If you're running locally: `nextflow run ebenz99/SPARCED-nf -c configs/local-nextflow.config`
    - If you're running with Kubernetes: `nextflow kuberun ebenz99/SPARCED-nf -v <PVC-name> -c configs/kube-nextflow.config`
2. Retrieving data
    - For local runs, all of your data should already be available to you in your new `work` directory.
    - For Kubernetes, after the run is finished, save your data from the PVC down to your local machine with `./kube-runner/kube-save.sh <pvc-name> <work-directory>` (`kube-save.sh` will find your `work directory` path as relative to your `workspace/$USER` directory in the PVC. So with the default configurations, it should just be `work`)

## Debugging

- A known problem with Docker on Mac is running Docker containers that make any use of the `/var/` directory (there's a symlink involved, look it up at your own peril). To fix this and run the pipeline locally on Mac, you have to change your default Docker settings. Go to Docker->Preferences->File Sharing and remove `/private` from your shared directories. Then add `/private/var/folders` and `/var/folders`. This should fix any related issues.
- Verify the path to your input data -> `kube-runner/kube-login.sh <pvc-name>` - find the files and check its path with `pwd`
- `Repository corrupted` or `Newer revision available` -> `kube-runner/kube-login.sh <pvc-name>` - the workflow has been updated since your last run, so manually delete the workflow contents stored in the PVC (which Nextflow stores in the `projects` folder)
 - Using `kube-login.sh` to obtain a shell into the PVC is also helpful for checking workflow output in any of the `.command` files nextflow generates in its `workDir` during runtime


## Acknowledgements:

The SPARCED model and SPARCED-nf pipeline is a product of [Birtwistle Lab](http://www.birtwistlelab.com/) and multiple colloborators, including [Feltus Lab](https://www.clemson.edu/science/departments/genetics-biochemistry/people/profiles/ffeltus), [Hasenauer Lab](https://www.mathematics-and-life-sciences.uni-bonn.de/en/group-members/jan-hasenauer), and [Robert Blake](https://bbs.llnl.gov/RobertBlake.html) from LLNL.



## The Acronym:
The acronym SPARCED is composed of following elements, based on the sub-models in the mechanistic ODE model.

### S: SBML
### P: Proliferation
### A: Apoptosis
###  R: Receptor
###  C: Cell
###  E: Expression
###  D: Death
