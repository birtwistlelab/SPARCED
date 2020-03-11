# SPARCED: An SBML Model of Pan-Cancer RTK Signaling and Stochastic Gene Expression 

SPARCED is a new version of the mechanistic pan-cancer signaling model by Birtwistle Lab, originally written in MATLAB. The cellular signaling portion of the earlier model is converted into an SBML model. The new scheme still combines stochastic gene expression with deterministic protein signaling. For more details please refer to:

1) The bioxrviv: XXX

2) The original model paper: https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005985

3) Older python version of the original model: https://github.com/birtwistlelab/Mechanistic_Pan-Cancer_Model



## Installation:

### 1) Install Anaconda
  - Follow the instructions depending on your machine: https://docs.anaconda.com/anaconda/user-guide/getting-started/
  - Install basic packages like python (3.7 or later), numpy, scipy, matplotlib, pandas, jupyter-notebooks...

### 2) Install AMICI
  - Follow the instructions on project home: https://github.com/ICB-DCM/AMICI
  - For WINDOWS machines, we recommend installation via creating a new dedicated conda environment. 

### 3) Install QTAntimony GUI (or Tellurium package): 
  - http://antimony.sourceforge.net/
  - https://tellurium.readthedocs.io/en/latest/installation.html#installation-options



## Model creation (OPTIONAL):

1) Modify the input files for desired model conditions.

2) Use the jupyter-notebook called "SPARCED_ModelCreateWrite.ipynb" to create the model (the Antimony text file) locally.

3) Use QTAntimony GUI to create the SBML xml file from the text file above.

4) (OPTIONAL) Use the second part of the "SPARCED_ModelCreateWrite.ipynb" to add annotations for model elements (i.e. species, compartments, and ratelaws).

5) Use the third part of the "SPARCED_ModelCreateWrite.ipynb" to convert the SBML (.xml) file to an AMICI model.



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