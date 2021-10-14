#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 13 21:13:39 2021

@author: arnab
"""
import pandas as pd
import numpy as np
import re
import libsbml
import os
import sys
import importlib
import amici
import amici.plotting


wd = str(os.getcwd()).replace("jupyter_notebooks","")



sbml_file = "SPARCED_au565.xml"
model_name= sbml_file[0:-4]
model_output_dir = model_name
sys.path.insert(0, os.path.join(wd,model_output_dir))
model_module = importlib.import_module(model_name)
model = model_module.getModel()



sbml_reader = libsbml.SBMLReader()
sbml_doc = sbml_reader.readSBML(sbml_file)
sbml_model = sbml_doc.getModel()

#%% create observables
ObsMat = pd.read_csv(os.path.join(wd,'input_files','Observables.txt'), header=0, index_col=0, sep='\t')
species_sheet = np.array([np.array(line.strip().split("\t")) for line in open(os.path.join(wd,'input_files','Species.txt'), encoding='latin-1')])
compartment_sheet = np.array([np.array(line.strip().split("\t")) for line in open(os.path.join(wd,'input_files','Compartments.txt'))])
Vc = float(compartment_sheet[compartment_sheet[:,0]=='Cytoplasm',1])
species_names = np.array([species_sheet[i][0] for i in range(1,len(species_sheet))])
VxPARCDL = np.array([species_sheet[i][1] for i in range(1,len(species_sheet))])
VxPARCDL = [float(compartment_sheet[compartment_sheet[:,0]==VxPARCDL[i],1][0]) for i in range(len(VxPARCDL))]
VxPARCDL = pd.Series(VxPARCDL, index=species_names)



formula_obs = []

for obs in ObsMat.columns:
    sp_obs = ObsMat.index[np.nonzero(np.array(ObsMat.loc[:,obs]>0))[0]]
    sp_obs_id = np.nonzero(np.array(ObsMat.loc[:,obs]>0))[0]
    Vr = VxPARCDL/Vc
    Vf = Vr*ObsMat.loc[:,obs].values
    
    if len(sp_obs) == 1:
        formula_i = sp_obs[0]+'*'+str(Vf[sp_obs_id][0])
    elif len(sp_obs) == 2:
        formula_i = str(sp_obs[0]+'*'+str(Vf[sp_obs_id][0])+'+'+sp_obs[1]+'*'+str(Vf[sp_obs_id][1]))
    elif len(sp_obs) > 2:
        formula_i = ''
        for j in range(len(sp_obs)-1):
            formula_i = formula_i+sp_obs[j]+'*'+str(Vf[sp_obs_id][j])+'+'
        formula_i = formula_i+str(sp_obs[-1])+'*'+str(Vf[sp_obs_id][-1])
    formula_obs.append(formula_i)

observables = {}
obs_names = list(ObsMat.columns)

for i in range(len(obs_names)):
    observables[obs_names[i]] = {}
    observables[obs_names[i]]['formula'] = formula_obs[i]
    
#%% update params from au565 sbml

params_all = pd.read_csv(os.path.join(wd,'input_files','params_all.csv'), header=0, index_col=0, sep=',')
for i in range(np.shape(params_all)[0]):
    params_all.loc[params_all.index[i],'value'] = float(sbml_model.getParameter(str(params_all.index[i])).getValue())
    


#%% compile initialized model

sbml_importer = amici.SbmlImporter(os.path.join(wd,'SPARCED_au565.xml'))


model_path = os.path.join(wd,model_output_dir)

constantParameters = np.array(params_all.index)

sbml_importer.sbml2amici('SPARCED_au565',
                         model_path,
                         verbose=False,
                         observables=observables,
                         constantParameters=constantParameters)
