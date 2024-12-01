#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This script is an example of how to build and run the SPARCED model.

# Model compilation
from compilation.launcher import launch_model_creation
# Model simulation
from simulation.launchers import launch_experiment_simulation

# First step is to build the model.
# You may also want to comment this section out and run a model that was
# previously built instead.
# Note that building the model may take several minutes.
launch_model_creation() # Process parsed arguments and launch model creation

# Finally you can run an experiment, i.e. one or several cell simulations
# within the given initial conditions.
launch_experiment_simulation() # Process parsed arguments and launch experiment
