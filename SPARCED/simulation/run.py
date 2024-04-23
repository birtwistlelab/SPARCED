#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys

import importlib
import numpy as np
import pandas as pd

from utils.arguments import parse_args
from utils.data_handling import load_input_data_file
from simulation.modules.RunSPARCED import RunSPARCED

# TODO: why to we have to pass again initial concentrations to the function
# when those are already in the SBML file?
# We should just pass the modified ligands and compound concentrations and
# leave the rest as is.

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
        species_initial_conditions[np.argwhere(species_names == l[0])] = l[1]
    if compound is not None:
        species_initial_conditions[np.argwhere(species_names == compound)] = float(dose)
    return(species_initial_conditions)

def save_output(model, file_prefix, cell_number, xoutS_all, xoutG_all, tout_all):
    #xoutS
    columnsS = [ele for ele in model.getStateIds()]
    condsSDF = pd.DataFrame(data=xoutS_all, columns=columnsS)
    condsSDF.to_csv(file_prefix+'_S_'+str(cell_number)+'.txt',sep="\t")
    condsSDF = None
    #xoutG
    columnsG = [x for n, x in enumerate(columnsS) if 'm_' in x]
    columnsG = columnsG[1:] # Skip header
    resa = [sub.replace('m_', 'ag_') for sub in columnsG]
    resi = [sub.replace('m_', 'ig_') for sub in columnsG]
    columnsG2 = np.concatenate((resa, resi), axis=None)
    condsGDF = pd.DataFrame(data=xoutG_all, columns=columnsG2)
    condsGDF.to_csv(file_prefix+'_G_'+str(cell_number)+'.txt',sep="\t")
    condsGDF = None
    #tout
    np.savetxt(file_prefix+'_T_'+str(cell_number)+'.txt', tout_all, newline="\t", fmt="%s")



def run_model(ligands):
    args = parse_args() 

    # Load model
    model_path = args.outputdir + "amici_" + args.name
    # TODO: fix the import in a clean way
    sys.path.insert(0, os.path.abspath(model_path))
    model_module = importlib.import_module(args.name)
    model = model_module.getModel()
    if args.verbose:
        print("{model_name}: Successfully loaded {module}\n"
              .format(model_name=args.name, module=model_module))
    # Set time points
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
    species_initial_conditions = initialize_simulation(args.inputdir+args.species, ligands,
                                 args.compound, args.dose)
    while cell_number < int(args.population):
        # INITIALIZATION
        # Set initial conditions
        model.setInitialStates(species_initial_conditions)
        # SIMULATION
        if args.verbose:
            print("{model_name}-{simulation_name}: Now ready to run the simulation.\n"
                   .format(model_name=args.name, simulation_name=args.simulation))
        xoutS_all, xoutG_all, tout_all = RunSPARCED(args.deterministic,
                float(args.time), species_initial_conditions, [], args.outputdir + "sbml_" + args.name + ".xml", model)
        print("{model_name}-{simulation_name}: Simulation is over. \
               Now saving results, please do not exit.\n"
               .format(model_name=args.name, simulation_name=args.simulation))
        # SAVE OUTPUT
        save_output(model, args.name, cell_number, xoutS_all, xoutG_all, tout_all)
        if args.verbose:
            print("{model_name}-{simulation_name}: Simulation saved.\n"
                   .format(model_name = args.name, simulation_name = args.simulation))
        cell_number += 1

