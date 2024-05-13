#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys

import numpy as np

from utils.arguments import parse_args
from utils.data_handling import *
from utils.sanitize import *
from simulation.run_model import run_experiment


def launch_experiment_simulation(perturbations: np.ndarray) -> None:
    """Launch experiment simulation function

    Small routine to process the parsed arguments and call run_model()

    Arguments:
        perturbations:  The array specifying initial conditions perturbations.

    Returns:
        Nothing.
    """

    args = parse_args()
    # Model
    model_name = sanitize_model_name(args.name)
    model_path = append_subfolder(args.model, model_name, True)
    amici_name = "amici_" + model_name
    amici_path = append_subfolder(args.model, amici_name, True)
    sbml_name = "sbml_" + model_name + ".xml"
    sbml_path = append.subfolder(args.model, sbml_name, True)
    # Input data files
    input_files = load_simulation_input_files(model_path, args.input_data, args.yaml)
    # Population size
    popsize = sanitize_popsize(args.popsize)
    # Runtime booleans
    is_SPARCED = not args.wild  # if it's not wild then it's SPARCED
    verbose = args.verbose
    run_experiment(model_name, args.simulation, args.results, amici_path,
                   sbml_path, input_files, perturbations, args.deterministic,
                   args.popsize, args.time, args.exchange, args.verbose, is_SPARCED)

def load_simulation_input_files(model_path: str | os.PathLike, data_folder: str,
                                yaml_name: stre) -> dict[str, str | os.PathLike]:
    """Load simulation input data files paths dictionnary

    Arguments:
        model_path: The model folder path.
        data_folder: The input data files folder name.
        yaml_name: The YAML configuration file name.

    Returns:
        A dictionnary containing all the input data file paths.
    """

    input_files_configuration = load_input_data_config(model_path, data_folder,
                                                       yaml_name)
    simulation_files = input_files_configuration["simulation"]
    # Load root of simulation input data files
    simulation_root = simulation_files.pop("root", None)
    simulation_data_path = append_subfolder(data_path, simulation_root, True)
    # Reconstruct full path for each input file listed in the dictionnary
    for input_file in simulation_files.keys():
        if input_file != "perturbations":
            simulation_files[input_file] = append_subfolder(simulation_data_path, simulation_files[input_file])
    return(simulation_files)


if __name__ == '__main__':
    launch_experiment_simulation()

