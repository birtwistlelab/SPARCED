#!/usr/bin/env python3


import sys
import libsbml
import importlib
import amici
import amici.plotting
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import re
from antimony import *
import argparse

from modules.paramSweep import paramSweep
from modules.copyDir import copyDirectory


parser = argparse.ArgumentParser(description='Provide arguments to build the SPARCED model')
parser.add_argument('--folder', metavar='folder', help='input data folder path')
# parser.add_argument('--paramfile', metavar='paramfile', help='file containing any non-default values')
args = parser.parse_args()


if args.folder == None:
    print("ERROR: No input folder given. Folder must be supplied with --folder option. Use -h for more information.")
    exit(1)

input_data_folder = args.folder

#move input data into working directory
copyDirectory(input_data_folder, os.getcwd()+"/")

# Antimony model name and text
fileModel = open('SPARCED.txt','w') # file name
fileModel.write("# PanCancer Model by Birtwistle Lab \n") # some explanation
fileModel.write("model SPARCED()\n") # model name

# SBML model we want to import
sbml_file = 'SPARCED.xml'
# Name of the model that will also be the name of the python module
model_name = sbml_file[0:-4] # 'BigModel_byparts_v1'
# Directory to which the generated model code is written
model_output_dir = model_name


compartments = []
volumes = []

# Create/write compartments
compartment_sheet = np.array([np.array(line.strip().split("\t")) for line in open('Compartments.txt')])

#read in each line minus the header row
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
species_sheet = np.array([np.array(line.strip().split("\t")) for line in open('Species.txt', encoding='latin-1')])

species_compartments = []
for row in species_sheet[1:]:
    species_compartments.append(row[1])
species_compartments = np.array(species_compartments)

#is this whitespace necessary?
fileModel.write("\n")

for idx,val in enumerate(species_sheet[1:]):
    fileModel.write("  Species ")
    fileModel.write("%s in %s" % (val[0], species_compartments[idx]))
    fileModel.write(";\n")

# Write reactions and rate laws
fileModel.write("\n\n  # Reactions:\n")

#reads in file from excel and gets rid of first row and column (they're data labels)
stoic_sheet = np.array([np.array(line.strip().split("\t")) for line in open('StoicMat.txt')])

#gets first column minus blank space at the beginning
stoic_columnnames = stoic_sheet[0]
stoic_rownames = [line[0] for line in stoic_sheet[1:]]
stoic_data = np.array([line[1:] for line in stoic_sheet[1:]])


ratelaw_sheet = np.array([np.array(line.strip().split("\t")) for line in open('Ratelaws.txt')])
ratelaw_data = np.array([line[1:] for line in ratelaw_sheet[1:]])

paramnames = []
paramvals = []
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
    else:
        # specific formula (non-mass-action)
        formula = ratelaw[1]
        j = 1
        params = np.genfromtxt(ratelaw[2:], float) # parameters
        params = params[~np.isnan(params)]
        if len(params) == 1:
            paramnames.append("k"+str(rowNum+1)+"_"+str(j))
            paramvals.append(float(ratelaw[j+1]))
            pattern = 'k\D*\d*'
            compiled = re.compile(pattern)
            matches = compiled.finditer(formula)
            for ematch in matches:
                formula = formula.replace(ematch.group(),paramnames[-1])
        else:
            for p in params:
                paramnames.append("k"+str(rowNum+1)+"_"+str(j))
                paramvals.append(float(ratelaw[j+1]))
                pattern1 = 'k(\D*)\d*'+'_'+str(j)
                compiled1 = re.compile(pattern1)
                matches1 = compiled1.finditer(formula)
                for ematch in matches1:
                    formula = formula.replace(ematch.group(),paramnames[-1])
                j +=1
    #don't include reactions without products or reactants
    if products == [] and reactants == []:
        pass
    else:
        fileModel.write("  %s: %s => %s; (%s)*%.6e;\n" % (stoic_columnnames[rowNum], " + ".join(reactants), " + ".join(products), formula, ratelaw[0]))

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

if loadFile("SPARCED.txt") == 1:
    print("Success loading antimony file")
else:
    print("Failed to load antimony file")
    exit(1)

if writeSBMLFile("SPARCED.xml","SPARCED") == 1:
    print("Success converting antimony to SBML")
else:
    print("Failure converting antimony to SBML")
    exit(1)


sbml_reader = libsbml.SBMLReader()
sbml_doc = sbml_reader.readSBML(sbml_file)
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

# # Change model name and write with a new name
# sbml_model.setName(model_name+'_Annot')
# sbml_model.setId(model_name+'_Annot')
# sbml_filewAnnot = model_name+'_Annot.xml'
# writer = libsbml.SBMLWriter()
# writer.writeSBML(sbml_doc, sbml_filewAnnot)

# sbml_file = sbml_filewAnnot # Comment-out to use the file name  w/t annotations
model_name = sbml_file[0:-4]
model_output_dir = model_name

sbml_reader = libsbml.SBMLReader()
sbml_doc = sbml_reader.readSBML(sbml_file)
sbml_model = sbml_doc.getModel()

# Create an SbmlImporter instance for our SBML model
sbml_importer = amici.SbmlImporter(sbml_file)

constantParameters = [params.getId() for params in sbml_model.getListOfParameters()]

observables={'actp53':{'formula': 'p53ac'},'phERK':{'formula': 'ppERK'},'egfLR':{'formula': 'EE1'}}

sbml_importer.sbml2amici(model_name,
                         model_output_dir,
                         verbose=False,
                         observables=observables,
                         constantParameters=constantParameters)

sys.path.insert(0, os.path.abspath(model_output_dir))
model_module = importlib.import_module(model_name)
model = model_module.getModel()

# set timepoints for which we want to simulate the model
model.setTimepoints(np.linspace(0, 86400, 2881))
# Create solver instance
solver = model.getSolver()
solver.setMaxSteps = 1e10
# Run simulation using default model parameters and solver options
rdata = amici.runAmiciSimulation(model, solver)
