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

from modules.SGEmodule import SGEmodule
from modules.RunPrep import RunPrep
from modules.RunSPARCED import RunSPARCED



# SBML model we want to itimemport
sbml_file = 'SPARCEDv6.xml'
# Name of the model that will also be the name of the python module
model_name = sbml_file[0:-4] # 'BigModel_byparts_v1'
# Directory to which the generated model code is written
model_output_dir = model_name



parser = argparse.ArgumentParser(description='Provide arguments to build the SPARCED model')
parser.add_argument('--deterministic', metavar='flagD', type=int, help='0 for deterministic run, 1 for stochastic')
parser.add_argument('--time', metavar='time', type=int, help='experiment run time (in hours)')
parser.add_argument('--feedTime', metavar='time', type=int, help='time frame for gene to feed into smbl (seconds)')
parser.add_argument('--cells', metavar='cells', type=int, help='number of cells to run with')
parser.add_argument('--Vn', metavar='Vn', help='number of cells to run with')
parser.add_argument('--Vc', metavar='Vc', help='number of cells to run with')
args = parser.parse_args()


if args.time == None or args.feedTime == None or args.cells == None or args.deterministic == None or args.Vn == None or args.Vc == None:
    print("ERRROR: missing arguments. Need to pass --time, --feedTime, --cells, --deterministic. Use -h for help.")

flagD = args.deterministic
numStocCells = args.cells
ts = args.feedTime
th = args.time
Vn = float(args.Vn)
Vc = float(args.Vc)


if flagD == 0:
    STIMligs = [100,100.0,100.0,100.0,100.0,100.0,1721.0] # EGF, Her, HGF, PDGF, FGF, IGF, INS
    flagWr = 1
    nmxlsfile = 'GrowthStim_stoc_'

    sys.path.insert(0, os.path.abspath(model_output_dir))

    species_sheet = np.array([np.array(line.strip().split("\t")) for line in open('Species.txt', encoding='latin-1')])

    species_initializations = []
    for row in species_sheet[1:]:
        species_initializations.append(float(row[2]))
    species_initializations = np.array(species_initializations)

    for nn in range(numStocCells):
        species_initializations[155:162] = STIMligs

        model_module = importlib.import_module(model_name)
        model = model_module.getModel()
        solver = model.getSolver() # Create solver instance
        solver.setMaxSteps = 1e10
        model.setTimepoints(np.linspace(0,ts)) # np.linspace(0, 30) # set timepoints

        xoutS_all, xoutG_all, tout_all = RunSPARCED(flagD,th,species_initializations,[],Vn,Vc,model)

        if flagWr==1:
            columnsS=[ele for ele in model.getStateIds()]
            columnsG = columnsS[773:914]
            resa = [sub.replace('m_', 'ag_') for sub in columnsG]
            resi = [sub.replace('m_', 'ig_') for sub in columnsG]
            columnsG2 = np.concatenate((resa, resi), axis=None)
            condsSDF = pd.DataFrame(data=xoutS_all,columns=columnsS)
            condsSDF.to_excel(nmxlsfile+'S_'+str(nn)+'.xlsx')
            condsSDF = None
            condsGDF = pd.DataFrame(data=xoutG_all,columns=columnsG2)
            condsGDF.to_excel(nmxlsfile+'G_'+str(nn)+'.xlsx')
            condsGDF = None

elif flagD == 1:
    STIMligs = [100,100.0,100.0,100.0,100.0,100.0,1721.0] # EGF, Her, HGF, PDGF, FGF, IGF, INS
    flagWr = 1
    nmxlsfile = 'GrowthStim_det_'

    sys.path.insert(0, os.path.abspath(model_output_dir))
    species_sheet = np.array([np.array(line.strip().split("\t")) for line in open('Species.txt', encoding='latin-1')])

    species_initializations = []
    for row in species_sheet[1:]:
        try:
            species_initializations.append(float(row[2]))
        except:
            print(row)
            exit(1)
    species_initializations = np.array(species_initializations)
    species_initializations[np.argwhere(species_initializations <= 1e-6)] = 0.0

    for nn in range(numStocCells):
        species_initializations[155:162] = STIMligs

        model_module = importlib.import_module(model_name)
        model = model_module.getModel()
        solver = model.getSolver() # Create solver instance
        solver.setMaxSteps = 1e10
        model.setTimepoints(np.linspace(0,ts)) # np.linspace(0, 30) # set timepoints

        xoutS_all, xoutG_all, tout_all = RunSPARCED(flagD,th,species_initializations,[],Vn,Vc,model)

        if flagWr==1:
            columnsS=[ele for ele in model.getStateIds()]
            columnsG = columnsS[773:914]
            resa = [sub.replace('m_', 'ag_') for sub in columnsG]
            resi = [sub.replace('m_', 'ig_') for sub in columnsG]
            columnsG2 = np.concatenate((resa, resi), axis=None)
            condsSDF = pd.DataFrame(data=xoutS_all,columns=columnsS)
            condsSDF.to_excel(nmxlsfile+'S_'+str(nn)+'.xlsx')
            condsSDF = None
            condsGDF = pd.DataFrame(data=xoutG_all,columns=columnsG2)
            condsGDF.to_excel(nmxlsfile+'G_'+str(nn)+'.xlsx')
            condsGDF = None
