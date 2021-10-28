#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 21 13:12:04 2021

@author: arnab
"""
import sys
import os
import pandas as pd
sys.path.append(os.getcwd()[0:os.getcwd().rfind('/')]+'/sparced'+'/bin')

wd = str(os.getcwd()).replace("jupyter_notebooks","")

sys.path.append(os.path.join(wd,'bin'))

import libsbml
import importlib
import amici
import amici.plotting
import sys
import numpy as np
import matplotlib.pyplot as plt
import re
from antimony import *
# from modules.paramSweep import paramSweep
# from modules.copyDir import copyDirectory
#%%
#copy input files over to current directory
current_dir = os.getcwd()
input_data_folder = current_dir[0:current_dir.rfind('/')+1]+'/sparced/'+'input_files'
# copyDirectory(input_data_folder, os.getcwd()+"/")

# Input and output file name definitions
fileComps = 'Compartments.txt' # input
fileSpecies = 'Species.txt' # input
fileStoic = 'StoicMat.txt' # input
fileRatelaws = 'Ratelaws.txt' # input
fileAntimony = 'SPARCED.txt'
# Antimony model name and text
fileModel = open(os.path.join(wd,fileAntimony),'w') # file name
fileModel.write("# PanCancer Model by Birtwistle Lab \n") # some explanation
fileModel.write("model SPARCED()\n") # model name

# SBML model we want to import
sbml_file = 'SPARCED.xml'
# Name of the model that will also be the name of the python module
model_name = sbml_file[0:-4] # 'BigModel_byparts_v1'
# Directory to which the generated model code is written
model_output_dir = model_name

# Initializing compartment and volume lists
compartments = []
volumes = []

# Create/write compartments
compartment_sheet = np.array([np.array(line.strip().split("\t")) for line in open(os.path.join(wd,'input_files',fileComps))])

#read in each line minus the header row of compartments file
for row in compartment_sheet[1:]:
    compartments.append(row[0])
    volumes.append(row[1])
    
#Write each compartment to model txt file
fileModel.write("\n  # Compartments and Species:\n")
for idx in range(len(compartments)):
    compName = compartments[idx]
    fileModel.write("  Compartment %s; " % (compName))
fileModel.write("\n")

# Write species and assign compartments
species_sheet = np.array([np.array(line.strip().split("\t")) for line in open(os.path.join(wd,'input_files','Species.txt'), encoding='latin-1')])

#read in each line minus the header row of species file
species_compartments = []
for row in species_sheet[1:]:
    species_compartments.append(row[1])
species_compartments = np.array(species_compartments)

#Write each species to model txt file
fileModel.write("\n")
for idx,val in enumerate(species_sheet[1:]):
    fileModel.write("  Species ")
    fileModel.write("%s in %s" % (val[0], species_compartments[idx]))
    fileModel.write(";\n")
    
    # Write reactions and rate laws
fileModel.write("\n\n  # Reactions:\n")

#reads in file from excel and gets rid of first row and column (they're data labels)
stoic_sheet = np.array([np.array(line.strip().split("\t")) for line in open(os.path.join(wd,'input_files','StoicMat.txt'))])

#gets first column minus blank space at the beginning, adds to stoic data list
stoic_columnnames = stoic_sheet[0]
stoic_rownames = [line[0] for line in stoic_sheet[1:]]
stoic_data = np.array([line[1:] for line in stoic_sheet[1:]])

#creates associated ratelaw data list
ratelaw_sheet = np.array([np.array(line.strip().split("\t")) for line in open(os.path.join(wd,'input_files','Ratelaws.txt'))])
ratelaw_data = np.array([line[1:] for line in ratelaw_sheet[1:]])

# builds the important ratelaw+stoic lines into the txt file 
paramnames = []
paramvals = []
paramrxns = []
paramidxs = []
for rowNum, ratelaw in enumerate(ratelaw_data):
    reactants = []
    products = []
    formula="k"+str(rowNum+1)+"*"

    for i, stoic_rowname in enumerate(stoic_rownames):
        stoic_value = int(stoic_data[i][rowNum])
        if stoic_value < 0:
            for j in range(0,stoic_value*-1):
                reactants.append(stoic_rowname)
                formula=formula+stoic_rowname+"*"
        elif stoic_value > 0:
            for j in range(0,stoic_value):
                products.append(stoic_rowname)

    if "k" not in ratelaw[1]:
        # the mass-action formula
        formula=formula[:-1]
        #the parameter
        paramnames.append("k"+str(rowNum+1))
        paramvals.append(np.double(ratelaw[1]))
        paramrxns.append(ratelaw_sheet[rowNum+1][0])
        paramidxs.append(int(0))
        
    else:
        # specific formula (non-mass-action)
        formula = ratelaw[1]
        j = 1
        params = np.genfromtxt(ratelaw[2:], float) # parameters
        params = params[~np.isnan(params)]
        if len(params) == 1:
            paramnames.append("k"+str(rowNum+1)+"_"+str(j))
            paramvals.append(float(ratelaw[j+1]))
            paramrxns.append(ratelaw_sheet[rowNum+1][0])
            paramidxs.append(int(0))
            pattern = 'k\D*\d*'
            compiled = re.compile(pattern)
            matches = compiled.finditer(formula)
            for ematch in matches:
                formula = formula.replace(ematch.group(),paramnames[-1])
        else:
            # for p in params:
            for q,p in enumerate(params):    
                paramnames.append("k"+str(rowNum+1)+"_"+str(j))
                paramvals.append(float(ratelaw[j+1]))
                paramrxns.append(ratelaw_sheet[rowNum+1][0])
                paramidxs.append(q)
                pattern1 = 'k(\D*)\d*'+'_'+str(j)
                compiled1 = re.compile(pattern1)
                matches1 = compiled1.finditer(formula)
                for ematch in matches1:
                    formula = formula.replace(ematch.group(),paramnames[-1])
                j +=1
    if ratelaw[0] == 'Cytoplasm':
        valcomp = 5.25e-12
    elif ratelaw[0] == 'Extracellular':
        valcomp = 5.00e-5
    elif ratelaw[0] == 'Nucleus':
        valcomp = 1.75e-12
    elif ratelaw[0] == 'Mitochondrion':
        valcomp = 3.675e-13
    #don't include reactions without products or reactants
    if products == [] and reactants == []:
        pass
    else:
        fileModel.write("  %s: %s => %s; (%s)*%.6e;\n" % (stoic_columnnames[rowNum], " + ".join(reactants), " + ".join(products), formula, valcomp))
        
params_all = pd.DataFrame({'value':paramvals,'rxn':paramrxns,'idx':paramidxs},index=paramnames)
params_all.to_csv(os.path.join(wd,'input_files','params_all.csv'),sep=',',header=True, index=True)


# Write compartment ICs
fileModel.write("\n  # Compartment initializations:\n")
for idx in range(len(compartments)):
    fileModel.write("  %s = %.6e;\n" % (compartments[idx], np.double(volumes[idx]))) # '%.3e'  "%.4g"
    fileModel.write("  %s has volume;\n" % (compartments[idx]))
    
# Write species ICs
fileModel.write("\n  # Species initializations:\n")
for idx, val in enumerate(species_sheet[1:]):
    fileModel.write("  %s = %.6e;\n" % (val[0],np.double(val[2]))) # '%.3e'  "%.4g"
    
    
    
# Write parameter ICs
fileModel.write("\n  # Parameter initializations:\n")
count = 0
for param in paramnames:
    fileModel.write("  %s = %.6e;\n" % (param, np.double(paramvals[count]))) # '%.3e'  "%.4g"
    count += 1
    
# Write other declarations
constantVars = ['Cytoplasm','Extracellular','Nucleus','Mitochondrion']

fileModel.write("\n  # Other declarations:\n")
fileModel.write("  const")
for constVar in constantVars[:-1]:
    fileModel.write("  %s," % (constVar))
#last item in row needs semicolon
fileModel.write("  %s;\n" % (constantVars[-1]))

# Write unit definitions
fileModel.write("\n  # Unit definitions:")
fileModel.write("\n  unit time_unit = second;")
fileModel.write("\n  unit volume = litre;")
fileModel.write("\n  unit substance = 1e-9 mole;")
fileModel.write("\n  unit nM = 1e-9 mole / litre;")
fileModel.write("\n")

# End the model file
fileModel.write("\nend")
# Close the file
fileModel.close()

# load model and convert to SBML
if loadFile(os.path.join(wd,fileAntimony)) == 1:
    print("Success loading antimony file")
else:
    print("Failed to load antimony file")
    exit(1)

if writeSBMLFile(os.path.join(wd,sbml_file),"SPARCED") == 1:
    print("Success converting antimony to SBML")
else:
    print("Failure converting antimony to SBML")
    exit(1)
    
#%%
# create observables
ObsMat = pd.read_csv(os.path.join(wd,'input_files',"Observables.txt"), sep='\t',header=0, index_col=0)
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




#%%
# create interaction components
sbml_reader = libsbml.SBMLReader()
sbml_doc = sbml_reader.readSBML(os.path.join(wd,sbml_file))
sbml_model = sbml_doc.getModel()

# Set species annotations
for idx,row in enumerate(species_sheet[1:]):
    Annot=""
    for col in range(4,(len(row))):
        aa=str(row[col].strip())
        if aa=="nan" or aa == "":
            break
        else:
            Annot=Annot+" "+row[col]
    sbml_model.getSpecies(row[0]).setAnnotation(Annot)
    
# Set compartment annotations
for row in compartment_sheet[1:]:
    sbml_model.getCompartment(row[0]).setAnnotation(row[2])
    
# Write with the same name or use the next section instead of below lines
writer = libsbml.SBMLWriter()
writer.writeSBML(sbml_doc, sbml_file)

# prepares to use interaction components to synthesize model
model_name = sbml_file[0:-4]
# model_output_dir = model_name
model_output_dir = os.path.join(wd,model_name)

sbml_reader = libsbml.SBMLReader()
sbml_doc = sbml_reader.readSBML(sbml_file)
sbml_model = sbml_doc.getModel()

# Create an SbmlImporter instance for our SBML model
sbml_importer = amici.SbmlImporter(sbml_file)

#sets important constants for model build
constantParameters = [params.getId() for params in sbml_model.getListOfParameters()]
# observables={'actp53':{'formula': 'p53ac'},'phERK':{'formula': 'ppERK'},'egfLR':{'formula': 'EE1'}}

# uses amici to make mathematical conversions for model (this step can take a while)
sbml_importer.sbml2amici(model_name,
                         model_output_dir,
                         verbose=False,
                         observables=observables,
                         constantParameters=constantParameters)

#%%

# sets environment for model creation
sys.path.insert(0, model_output_dir)
model_module = importlib.import_module(model_name)
model = model_module.getModel()

# Create solver instance
solver = model.getSolver()
solver.setMaxSteps = 1e10

# Run simulation using default model parameters and solver options
rdata = amici.runAmiciSimulation(model, solver)
