#!/usr/bin/env python3

import libsbml
import importlib
import amici
import amici.plotting
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from datetime import datetime
import scipy.stats
import argparse

from modules.RunSPARCED import RunSPARCED



# SBML model we want to import
sbml_file = 'SPARCED.xml'
# Name of the model that will also be the name of the python module
model_name = sbml_file[0:-4]
# Directory to which the generated model code is written
model_output_dir = model_name



parser = argparse.ArgumentParser(description='Provide arguments to build the SPARCED model')
parser.add_argument('--deterministic', metavar='flagD', type=int, help='0 for deterministic run, 1 for stochastic')
parser.add_argument('--time', metavar='time', type=int, help='experiment run time (in hours)')
parser.add_argument('--Vn', metavar='Vn', help='the volume of the nucleus in liters')
parser.add_argument('--Vc', metavar='Vc', help='the volume of the cytoplasm in liters')
parser.add_argument('--outfile', metavar='outfile', help='the prefix for the name of the output files')
args = parser.parse_args()


if args.time == None or args.deterministic == None or args.Vn == None or args.Vc == None or args.outfile == None:
    print("ERROR: missing arguments. Need to pass --time, --deterministic, --Vn, --Vc, --outfile. Use -h for help.")

flagD = args.deterministic
th = args.time
Vn = float(args.Vn)
Vc = float(args.Vc)
outfile = args.outfile
ts = 30


if flagD == 0:
    flagWr = 1
    nmxlsfile = outfile
    
    sys.path.insert(0, os.path.abspath(model_output_dir))

    species_sheet = np.array([np.array(line.strip().split("\t")) for line in open('Species.txt', encoding='latin-1')])

    species_initializations = []
    for row in species_sheet[1:]:
        species_initializations.append(float(row[2]))
    species_initializations = np.array(species_initializations)

    model_module = importlib.import_module(model_name)
    model = model_module.getModel()
    solver = model.getSolver() # Create solver instance
    solver.setMaxSteps = 1e10
    model.setTimepoints(np.linspace(0,ts)) # np.linspace(0, 30) # set timepoints

    xoutS_all, xoutG_all, tout_all = RunSPARCED(flagD,th,species_initializations,[],sbml_file,model)

    if flagWr==1:
        columnsS=[ele for ele in model.getStateIds()]
        columnsG = [x for n, x in enumerate(columnsS) if 'm_' in x]
        columnsG = columnsG[1:]
        resa = [sub.replace('m_', 'ag_') for sub in columnsG]
        resi = [sub.replace('m_', 'ig_') for sub in columnsG]
        columnsG2 = np.concatenate((resa, resi), axis=None)
        condsSDF = pd.DataFrame(data=xoutS_all,columns=columnsS)
        condsSDF.to_excel(nmxlsfile+'S_0.xlsx')
        condsSDF = None
        condsGDF = pd.DataFrame(data=xoutG_all,columns=columnsG2)
        condsGDF.to_excel(nmxlsfile+'G_0.xlsx')
        condsGDF = None

elif flagD == 1:
    flagWr = 1
    nmxlsfile = outfile

    sys.path.insert(0, os.path.abspath(model_output_dir))
    species_sheet = np.array([np.array(line.strip().split("\t")) for line in open('Species.txt', encoding='latin-1')])

    species_initializations = []
    for row in species_sheet[1:]:
        species_initializations.append(float(row[2]))

    species_initializations = np.array(species_initializations)
    species_initializations[np.argwhere(species_initializations <= 1e-6)] = 0.0

    model_module = importlib.import_module(model_name)
    model = model_module.getModel()
    solver = model.getSolver()          # Create solver instance
    solver.setMaxSteps = 1e10
    model.setTimepoints(np.linspace(0,ts)) # np.linspace(0, 30) # set timepoints

    xoutS_all, xoutG_all, tout_all = RunSPARCED(flagD,th,species_initializations,[],sbml_file,model)

    if flagWr==1:
        columnsS=[ele for ele in model.getStateIds()]
        columnsG = [x for n, x in enumerate(columnsS) if 'm_' in x]
        columnsG = columnsG[1:]
        resa = [sub.replace('m_', 'ag_') for sub in columnsG]
        resi = [sub.replace('m_', 'ig_') for sub in columnsG]
        columnsG2 = np.concatenate((resa, resi), axis=None)
        condsSDF = pd.DataFrame(data=xoutS_all,columns=columnsS)
        condsSDF.to_excel(nmxlsfile+'S_0.xlsx')
        condsSDF = None
        condsGDF = pd.DataFrame(data=xoutG_all,columns=columnsG2)
        condsGDF.to_excel(nmxlsfile+'G_0.xlsx')
        condsGDF = None
