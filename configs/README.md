# Editing Your Config Files


## Local Config

Make the following edits to the `local-nextflow.config`
1. Change the value of `params.input_dir` to be the absolute path of the input directory on your local machine
2. Edit the `params.speciesVals` value to the string that encompasses any parameter sweeps or changes you want reflected in the run-time environment to the model (Important Note: If you want any of these values built into the *creation* of the model, not just its simulation, you'll need to edit their values in the input_data files directly)
