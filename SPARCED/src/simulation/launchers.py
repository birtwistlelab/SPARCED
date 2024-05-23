#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys

import numpy as np

from utils.arguments import parse_args
from utils.data_handling import *
from utils.sanitize import *
from simulation.run_model import run_experiment


def launch_experiment_simulation() -> None:
    """Launch experiment simulation function

    Small routine to process the parsed arguments and call run_model()

    Arguments:
        None.

    Returns:
        Nothing.
    """

    args = parse_args()
    # Model
    model_name = sanitize_model_name(args.name)
    model_path = append_subfolder(args.model, model_name, True)
    amici_name = "amici_" + model_name
    amici_path = append_subfolder(model_path, amici_name, True)
    sbml_name = "sbml_" + model_name + ".xml"
    sbml_path = append_subfolder(model_path, sbml_name, True)
    # Load simulation conditions data files
    input_folder = append_subfolder(model_path, args.input_data, True)
    simulation_files = load_simulation_config(input_folder, args.yaml)
    # Load incubation and perturbations
    perturbations = incubation_conditions = incubation_duration = None
    if "perturbations" in simulation_files:
        perturbations = load_perturbations(simulation_files["root"], simulation_files["perturbations"], args.perturbations)
    if "incubation" in simulation_files:
        incubation_conditions = load_perturbations(simulation_files["root"], simulation_files["incubation"], args.perturbations)
        incubation_duration = float(simulation_files["incubation"]["duration"])
    # Population size
    popsize = sanitize_popsize(args.population_size)
    # Runtime booleans
    is_SPARCED = not args.wild  # if it's not wild then it's SPARCED
    verbose = args.verbose
    run_experiment(model_name, args.simulation, args.results, amici_path,
                   sbml_path, simulation_files, args.deterministic, popsize,
                   args.time, args.exchange, args.verbose, is_SPARCED,
                   perturbations, incubation_conditions, incubation_duration)

def load_perturbations(input_folder: str | os.PathLike,
                       config: dict[str, str], file: str=None) -> np.ndarray:
    """Load perturbations

    Load an input perturbations file structured as tab separated, and remove
    the first column (human-readable) to keep only the second (SPARCED symbol)
    and the third (concentration in nM) columns.

    Arguments:
        input_folder
        config
        file

    Returns:
        A numpy array containing ligands symbols and concentrations.
    """
        
    perturbations_root = append_subfolder(input_folder, config["root"])
    # Choose the file to load
    if file != None:
        perturbations_file = append_subfolder(perturbations_root, file)
    else:
        perturbations_file = append_subfolder(perturbations_root, config["default"])
    perturbations = load_input_data_file(perturbations_file)
    # Convert concentrations from strings to floats
    for row in perturbations:
        row[2] = float(row[2])
    return(perturbations[:,1:])

def load_simulation_config(data_folder, f_config):
    # Load the YAML configuration file in a dictionnary
    simulation_files = load_input_data_config(data_folder, f_config)
    # Use only the portion related to simulation data
    simulation_files = simulation_files["simulation"]
    # Reconstruct full path for each input file listed
    simulation_files["root"] = append_subfolder(data_folder, simulation_files["root"], True)
    keys_to_exclude = ['root', 'perturbations', 'incubation']
    for input_file in simulation_files.keys():
        if input_file not in keys_to_exclude:
            simulation_files[input_file] = append_subfolder(simulation_files["root"], simulation_files[input_file], True)
    return(simulation_files) 


if __name__ == '__main__':
    launch_experiment_simulation()

