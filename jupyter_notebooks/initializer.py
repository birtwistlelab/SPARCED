
import pandas as pd
import numpy as np
import libsbml
import os
import sys
import importlib
import amici
import amici.plotting

import matplotlib.pyplot as plt
import matplotlib as mpl
import time
import yaml
import pypesto
import pypesto.petab
import pypesto.optimize as optimize
import pypesto.visualize as visualize
import petab


mpl.rcParams['figure.dpi'] = 300

#%%

# test module - load model

sbml_file = "SPARCED.xml"
model_name= sbml_file[0:-4]
model_output_dir = model_name
sys.path.insert(0, os.path.abspath(model_output_dir))
model_module = importlib.import_module(model_name)
model = model_module.getModel()



sbml_reader = libsbml.SBMLReader()
sbml_doc = sbml_reader.readSBML(sbml_file)
sbml_model = sbml_doc.getModel()
params_model = []
[params_model.append(p.getId()) for p in sbml_model.getListOfParameters()]


# get kTLids

# fileRatelaws = "Ratelaws.txt"
# Ratelawsf = pd.read_csv(fileRatelaws,header=0,index_col=0,sep='\t')
ratelaw_sheet = np.array([np.array(line.strip().split("\t")) for line in open('Ratelaws.txt')])
ratelaw_data = np.array([line[1:] for line in ratelaw_sheet[1:]])



numRlaws = Ratelawsf.shape[0]






gene_params = pd.read_csv(os.path.join('input_files',"cell_definition_mcf10a_sparced.csv"), index_col="gene")
model_genes = gene_params.index
numberofgenes = len(model_genes)
S_PARCDL = pd.read_csv(os.path.join('input_files',"StoicMat.txt"), header=0, index_col=0, sep='\t')
S_TL = S_PARCDL.iloc[:,2:(2+numberofgenes)]

ObsMat = pd.read_csv(os.path.join('input_files','Observables.txt'), header=0, index_col=0, sep='\t')
NumObs = len(ObsMat.columns)
