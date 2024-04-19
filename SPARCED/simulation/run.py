#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys

import importlib
import libsbml
import numpy as np

from utils.arguments import parse_args()
from utils.data_handling import load_input_data_file


def initialize_simulation(f_species: str, ligands: np.ndarray,
                          compound: str=None, dose: float=None) -> np.ndarray:
    # Warning: first line is considered as a header and hence skipped
    # Species names should be on first column
    # Species initial conditions, should be on third column
    species = load_input_data_file(f_species)
    species_names = []
    species_initial_conditions = []
    for row in species[1:]:
        species_names = np.append(species_names, str(row[0]))
        species_initial_conditions = np.append(species_initial_conditions, float(row[2]))
    # Any concentration bellow 1e-6 is considered as zero (0)
    species_initial_conditions[np.argwhere(species_initial_conditions <= 1e-6)] = 0.0
    # Adjust ligands concentration
    for l in ligands:
        species_initial_conditons[np.argwhere(species_names == l[0])] = l[1]
    if compound is not None:
        species_initial_conditions[np.argwhere(species_names == compound)] = float(dose)
    return(species_initial_conditions)


def run_model(ligands):
    args = parse_args() 

    # Import model
    model_path = args.outputdir + args.name
    model_module = importlib.import_module(model_path)
    if arg.verbose:
        print(model_module)
    # Set model
    model = model_module.getModel()
    model.setTimepoints(np.linspace(0, args.exchange, 2))
    # -- begin
    # Save list of species
    species_list = list(model.getStateIds())
    # Could save those in an output file but what for???
    # -- end

    # Run simulations
    cell_number = 0
    try:
        assert args.population > 0
    except:
        print("{model_name}-{simulation_name}: Cell population size should be \
               superior to zero (0). Current population size: {size}.\n"
               .format(model_name=args.name, simulation_name=args.simulation,
                        size=args.population))
        sys.exit(0)
    # Don't alter the initial conditions later!!!
    species_initial_conditions = initialize_simulation(args.species, ligands,
                                 args.compound, args.dose)
    while cell_number < int(args.population):
        # INITIALIZATION
        # Set initial conditions
        # Actually you might want to initialize before the while loop to avoid
        # recomputing those initial conditions each time
        model.setInitialStates(species_initial_conditions)
        # SIMULATION
        if args.verbose:
            print("{model_name}-{simulation_name}: Now ready to run the simulation.\n"
                   .format(model_name=args.name, simulation_name=args.simulation))
        # Broken stuff here
        xoutS_all, xoutG_all, tout_all = run_sparced,(args.deterministic, float(args.time), species_initial_conditions, args.name + ".xml", model)
        # SAVE OUTPUT
        save_output(model, args.name, cell_number, xoutS_all, xoutG_all, tout_all)
        if args.verbose:
            print("{model_name}-{simulation_name}: Simulation is finished and saved.\n"
                   .format(model_name = args.name, simulation_name = args.simulation))
        cell_number += 1

