#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Sparced Pipeline: Model Simulation"""

import os
import sys

import argparse
import importlib
import libsbml
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from bin.run_sparced import run_sparced

# ---------------------------------------------------------------------------------------------------------------------------------------------------

def save_output(model, file_prefix, cell_number, xoutS_all, xoutG_all, tout_all):
    # xoutS
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
    condsGDF = pd.DataFrame(data=xoutG_all,columns=columnsG2)
    condsGDF.to_csv(file_prefix+'_G_'+str(cell_number)+'.txt',sep="\t") 
    condsGDF = None
    #tout
    np.savetxt(file_prefix+'_T_'+str(cell_number)+'.txt', tout_all, newline="\t", fmt="%s")

# ---------------------------------------------------------------------------------------------------------------------------------------------------

if __name__ == '__main__':
    args = parse_args()

    # Import model
    sbml_model_name = model_output_dir = args.sbml
    sys.path.insert(0, os.path.abspath(model_output_dir))
    model_module = importlib.import_module(sbml_model_name)
    if args.verbose: print(model_module)
    # Set model
    model = model_module.getModel()
    model.setTimepoints(np.linspace(0, args.exchange, 2))
    # Save list of species
    species_all = list(model.getStateIds())
    np.savetxt(args.list, species_all, newline="\t", fmt="%s")

    # Run simulations
    cell_number = 0
    while cell_number < int(args.pop):
        # INITIALIZATION
        # Set initial conditions
        species_sheet = np.array([np.array(line.strip().split("\t")) for line in open(args.species, encoding='latin-1')])
        species_initializations = []
        for row in species_sheet[1:]:
            species_initializations = np.append(species_initializations, float(row[2]))
            species_initializations = np.array(species_initializations)
            species_initializations[np.argwhere(species_initializations <= 1e-6)] = 0.0
        # Input ligand concentrations (in order): EGF, Her, HGF, PDGF, FGF, IGF, INS
        STIMligs = [float(args.egf), 0.0, 0.0, 0.0, 0.0, 0.0, float(args.ins)] # in nM, in extracellular volume
        species_initializations[155:162] = STIMligs
        # Add compound if specified
        if args.compound is not None: species_initializations[species_all.index(args.compound)] = args.dose
        model.setInitialStates(species_initializations)
        # SIMULATION
        if args.verbose: print("SPARCED: Now ready to run a simulation")
        xoutS_all, xoutG_all, tout_all = run_sparced(args.deterministic, float(args.time), species_initializations, sbml_model_name + ".xml", model)
        # SAVE OUTPUT
        save_output(model, args.name, cell_number, xoutS_all, xoutG_all, tout_all)
        if args.verbose: print("SPARCED: Simulation is now over")
        cell_number += 1
