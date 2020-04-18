# SPARCED: An SBML Model of Pan-Cancer RTK Signaling and Stochastic Gene Expression

SPARCED is a new version of the mechanistic pan-cancer signaling model by Birtwistle Lab, originally written in MATLAB. The cellular signaling portion of the earlier model is converted into an SBML model. The new scheme still combines stochastic gene expression with deterministic protein signaling. For more details please refer to:

1) The bioxrviv: XXX

2) The original model paper: https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005985

3) Older python version of the original model: https://github.com/birtwistlelab/Mechanistic_Pan-Cancer_Model


## Use Instructions

1. Clone this repository from Github
2. Install Nextflow
3. Move your input data in tab-separated format into the `input_data` directory
4. Use `./kube-scripts/kube-load.sh <pvc-name> input_data` to load your input data to the PVC
5. Change the `nextflow.config` file to meet your specifications (to check your username, type `echo $USER` in your terminal)




## New installation instructions
sudo apt install libatlas-base-dev
sudo apt-get install libhdf5-serial-dev
sudo apt-get install swig
pip3 install requirements.txt


## Model simulation:

After completing the "Model Creation" part above or by using the provided files, you can use the jupyter-notebook called "RunModel.ipynb" to run the model for the specified conditions.



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
