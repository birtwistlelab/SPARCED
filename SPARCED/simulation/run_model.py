#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys

import importlib
import numpy as np
import pandas as pd

from simulation.modules.RunSPARCED import RunSPARCED
from simulation.utils.initial_conditions import load_species_initial_conditions
from simulation.utils.output import save_simulation_output

# TODO: why to we have to pass again initial concentrations to the function
# when those are already in the SBML file?
# We should just pass the modified ligands and compound concentrations and
# leave the rest as is.

def run_experiment(model_name: str, model_path: str, simulation_name: str,
                   sbml_model: str, is_deterministic: bool, exchange: int,
                   population: int, duration: float, f_species: str, ligands: np.ndarray,
                   verbose: bool, is_SPARCED: bool, compound: str=None,
                   dose: str=None) -> None:
    """
    Run an experiment (one or several cells with the same initial conditions)
    """
    # Load model
    sys.path.insert(0, os.path.abspath(model_path)) # TODO: fix the import on the next line to avoid messing up with paths
    model_module = importlib.import_module(model_name)
    model = model_module.getModel()
    if verbose:
        print("{model_name}: Successfully loaded {module}\n"
               .format(model_name=model_name, module=model_module))
    # Set time points
    model.setTimepoints(np.linspace(0, exchange, 2))
    # Compute initial conditions
    species_initial_conditions = load_species_initial_conditions(f_species, ligands, compound, dose)
    # Run experiment
    cell_number = 0
    while cell_number < int(population):
        run_single_simulation(model, simulation_name, cell_number, sbml_model,
                              is_deterministic, exchange, duration,
                              species_initial_conditions, verbose)
        cell_number += 1 

def run_single_simulation(model, simulation_name: str, simulation_number: int,
                          sbml_model: str, is_deterministic: bool,
                          exchange: int, duration: float,
                          species_initial_conditions: np.ndarray,
                          verbose: bool) -> None:
    """
    Run a single simulation of SPARCED
    """
    # Set initial conditions
    model.setInitialStates(species_initial_conditions)
    if verbose:
        print("{name} n°{number}: Now ready to run\n"
              .format(name=simulation_name, number=simulation_number))
    # Run the simulation
    xoutS_all, xoutG_all, tout_all = RunSPARCED(is_deterministic, float(duration),
                                                species_initial_conditions, [],
                                                sbml_model, model)
    if verbose:
        print("{name} n°{number}: Simulation is over. Now saving results, \
                please do not exit.\n".format(name=simulation_name,
                                              number=simulation_number))
    # Save output
    save_simulation_output(model, simulation_name, simulation_number,
                           xoutS_all, xoutG_all, tout_all)
    if verbose:
        print("{name} n°{number}: Simulation saved.\n"
              .format(name=simulation_name, number=simulation_number))

