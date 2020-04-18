# SPARCED: An SBML Model of Pan-Cancer RTK Signaling and Stochastic Gene Expression

SPARCED is a new version of the mechanistic pan-cancer signaling model by Birtwistle Lab, originally written in MATLAB. The cellular signaling portion of the earlier model is converted into an SBML model. The new scheme still combines stochastic gene expression with deterministic protein signaling. For more details please refer to:

1) The bioxrviv: XXX

2) The original model paper: https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005985

3) Older python version of the original model: https://github.com/birtwistlelab/Mechanistic_Pan-Cancer_Model


## Cloud Use Instructions
1. Clone this repository
3. Install Nextflow
4. Move your input data in tab-separated format into the `input_data` directory (or edit the files already there)
5. Use `./kube-scripts/kube-load.sh <pvc-name> input_data` to load your input data to the PVC
6. Change the `nextflow.config` file to meet your specifications (for more details, check the **Editing your nextflow.config** section of this README)
8. Run from the command line with `nextflow kuberun ebenz99/SPARCED -C nextflow.config`
9. After the run is finished, save your data from the PVC down to your laptop with `./kube-scripts/kube-save.sh <pvc-name> <work-directory>` (this work directory path is relative to your `workspace/$USER` directory in the PVC. So with the default configurations, it should just be `work`)


## Containerized requirements
sudo apt install libatlas-base-dev
sudo apt-get install libhdf5-serial-dev
sudo apt-get install swig
pip3 install requirements.txt


## Acknowledgements:

This work is a product of [Birtwistle Lab](http://www.birtwistlelab.com/) and multiple colloborators, including [Feltus Lab](https://www.clemson.edu/science/departments/genetics-biochemistry/people/profiles/ffeltus), [Hasenauer Lab](https://www.mathematics-and-life-sciences.uni-bonn.de/en/group-members/jan-hasenauer), and [Robert Blake](https://bbs.llnl.gov/RobertBlake.html) from LLNL.



## The Acronym:
The acronym SPARCED is composed of following elements, based on the sub-models in the mechanistic ODE model.

### S: SBML

### P: Proliferation

### A: Apoptosis

###  R: Receptor

###  C: Cell Cycle

###  E: Expression

###  D: Death
