#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This script is an example of how to build and run the SPARCED model.

# Model creation
from compilation.utils import *
from compilation.create_model import *
# Model simulation
from simulation.ligands import *
from simulation.run import *

# First step is to build the model.
# You may also want to comment this section out and run a model that was
# previously built instead.
# Note that building the model may take several minutes.
launch_model_creation() # Process parsed arguments and launch model creation

# Second step is to set the initial conditions.
# You may want to load a custom extracellular concentration of ligands.
# The basic function is usefull if you don't want to add any ligands or if you
# need to adjust only EGF and/or insulin concentrations.
# TODO: implement reading arguments
ligands = basic_ligands(1.0, 17.21) # WARNING: hard-coded values
# If you need to adjust more species or simply if you're more comfortable
# working with input files, use the load function instead.
# ligands = load_input_ligands(your_file_path)

# Run model
run_model(ligands)
