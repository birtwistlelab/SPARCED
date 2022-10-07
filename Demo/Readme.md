# SPARCED model demo

This Demo folder containes pre-created and pre-compiled SPARCED model files. Users can create the model files using "createModel" notebook and can run the model using "runModel" notebook.
The pre-simulated data are also given (GrowthStim_S_demo.txt and GrowthStim_G_demo.txt) and users can compare their simulations to these provided files.

## createModel notebook:

- This notebook builds the 'SPARCED' hybrid model from the given input files. This notebook is already run and the model instance (AMICI required) files are saved in the `SPARCED` subfolder. Users can run the notebook again, with a new model file name to compare to our model provided here.

## createModel_o4a notebook:

- This notebook has functionality analogous to the createModel notebook, for the integrated-SBML version of 'SPARCED' model. Compiled model files are saved in the 'SPARCEDo4a_v1' folder.

## testModel notebook:

- This notebook imports and tests both 'Hybrid' and 'Integrated-SBML' versions of the 'SPARCED' model and draws comparison plots to demonstrate accuracy of outputs. The user needs to make sure both versions of the model have been already built using the above 'createModel' scripts and current version of AMICI available in the environment before this script can be run.

### Simulation conditions:
- Stimulate with ligands EGF = 1nM, HGF = 0.005nM
- 72 hours of simulation
- Deterministic simulation mode (flagD=1)

## runModel notebook:

- The stimulation conditions are reported below, where users can run the notebook and compare to existing results/figures.

### Simulation conditions:
- Stimulate with all ligands
- 12 hours of simulation
- Deterministic simulation mode (flagD=1)


## System properties

- Windows desktop computer (Windows 10 Education, Intel Core i5-3470 CPU @ 3.20GHz, 8.00GB RAM, 64-bit operating system)
- Model creation time: ~2 min
- Model compilation time: ~20min
- Model simulation time <1min
