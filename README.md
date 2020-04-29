# SPARCED-nf: A Nextflow Pipeline for SPARCED

SPARCED-nf is a Nextflow pipeline designed to be a more scalable and user-friendly version of the SPARCED model (previously the mechanistic pan-cancer signaling model) by the Birtwistle Lab. With minimal set-up, a user can configure the model for high-intensity runs on a Kubernetes cluster, or small-scale experiments on their local machine. More information on the model itself can be found [here](https://github.com/birtwistlelab/SPARCED).


## Dependencies

- [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html)
- [Docker](https://docs.docker.com/get-docker/)
- (Optional) [kubectl](https://kubernetes.io/docs/tasks/tools/install-kubectl/)

## Instructions to run locally
1. Clone this repository
2. Make sure the dependencies listed above are installed
3. Edit the files in the `input_data` folder as needed. These values will be built into the *creation* of the model. For editing the values present for the model's *simulation*, see the directions accompanying the next step.
4. Edit the `local-nextflow.config` file located in the `configs` folder. (for help, see [here](https://github.com/ebenz99/SPARCED-nf/tree/master/configs))
5. Navigate to the base directory of this project, and run the workflow with `nextflow kuberun . -c configs/local-nextflow.config`

## Instructions to run with Kubernetes
1. Clone this repository
2. Make sure the dependencies listed above are installed
3. Ensure you have a Kubernetes `config` file for your chosen cluster located in your `~/.kube` folder
4. Edit the files in the `input_data` folder as needed. These values will be built into the *creation* of the model. For editing the values present for the model's *simulation*, see the directions accompanying the next step.
5. Use `./kube-runner/kube-load.sh <pvc-name> input_data` to load your input data to the PVC of the kube cluster. `kube-load.sh` assumes a `/workspaces` folder as the base of the PVC, and saves this input data at the path `/workspaces/$USER/input_data/`.
6. Edit the `kube-nextflow.config` file located in the `configs` folder. (for help, see [here](https://github.com/ebenz99/SPARCED-nf/tree/master/configs))
8. Run from the command line with `nextflow kuberun ebenz99/SPARCED -c configs/kube-nextflow.config`
9. After the run is finished, save your data from the PVC down to your laptop with `./kube-runner/kube-save.sh <pvc-name> <work-directory>` (`kube-save.sh` will find your `work directory` path as relative to your `workspace/$USER` directory in the PVC. So with the default configurations, it should just be `work`)

## Debugging

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

###  C: Cell Cycle

###  E: Expression

###  D: Death
