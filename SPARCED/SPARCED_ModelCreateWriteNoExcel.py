# Antimony model name and text
fileModel = open('SPARCEDv6.txt','w') # file name
fileModel.write("// PanCancer Model by Birtwistle Lab \n") # some explanation
fileModel.write("model SPARCEDv6()\n") # model name

# Input and output file name definitions
fileComps = 'Compartments_v6.xlsx' # input
fileSpecies = 'Species_v6.xlsx' # input
fileStoic = 'StoicMat_v6.xlsx' # input
fileRatelaws = 'Ratelaws_v6.xlsx' # input
fileParams = 'ParamsO_v6.xlsx' # output

# SBML model we want to import
sbml_file = 'SPARCEDv6.xml'
# Name of the model that will also be the name of the python module
model_name = sbml_file[0:-4] # 'BigModel_byparts_v1'
# Directory to which the generated model code is written
model_output_dir = model_name

import libsbml
import importlib
import amici
import amici.plotting
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import csv
import re

compartments = []
volumes = []

# Create/write compartments
compartment_data = np.array([np.array(line.strip().split(",")) for line in open('Compartments_v6.csv')])
compartment_headers = compartment_data[0]

for row in compartment_data[1:]:
    compartments.append(row[0])
    volumes.append(row[1])

#Write each compartment to model
fileModel.write("\n  // Compartments and Species:\n")
for idx in range(len(compartments)):
    compName = compartments[idx]
    print(compName)
    fileModel.write("  Compartment %s; " % (compName))
fileModel.write("\n")


# Write species and assign compartments
species_data = np.array([np.array(line.strip().split(",")) for line in open('Species_v6.csv')])
species_headers = compartment_data[0]
species_compartments = []
for row in species_data[1:]:
    species_compartments.append(row[1])
species_compartments = np.array(species_compartments)

#suspicious minus 1
numSpecies = len(species_compartments) - 1

#is this whitespace necessary?
fileModel.write("\n")

print(numSpecies)

# ##----------------------REFACTORED TO HERE---------------------------
#
# ICf = pd.read_excel(fileSpecies,header=0,index_col=0)
# count = 0
# count2 = 0
#
# species_csv = csv.reader(open("Species_v6.csv"))
# for val in ICf.index:
#     if count2 == 0:
#         fileModel.write("  Species ")
#     if count == numSpecies:
#         fileModel.write("%s in %s;" % (val, SpComps[count]))
#     elif count2 < 4:
#         fileModel.write("%s in %s," % (val, SpComps[count]))
#         count2 = count2 + 1
#     else:
#         fileModel.write("%s in %s;\n" % (val, SpComps[count]))
#         count2 = 0
#     count = count + 1

# for val in ICf.index:
#     if count2 == 0:
#         fileModel.write("  Species ")
#     if count == numSpecies:
#         fileModel.write("%s in %s;" % (val, SpComps[count]))
#     elif count2 < 4:
#         fileModel.write("%s in %s," % (val, SpComps[count]))
#         count2 = count2 + 1
#     else:
#         fileModel.write("%s in %s;\n" % (val, SpComps[count]))
#         count2 = 0
#     count = count + 1
#
# # Write reactions and rate laws
# fileModel.write("\n\n  // Reactions:\n")
#
# Sdf = pd.read_excel(fileStoic,header=0,index_col=0)
# S = Sdf.values
# Ratelawsf = pd.read_excel(fileRatelaws,header=0,index_col=0)
# numRlaws = Ratelawsf.shape[0]
#
# paramnames = []
# paramvals = []
# count = 0
# for rr in range(numRlaws):
#     reactants = []
#     products = []
#     i = 0
#     formula="k"+str(count+1)+"*"
#     line = np.array(Ratelawsf.iloc[rr])
#     for sp in Sdf.index:
#         if S[i,count] < 0:
#             for j in range(0,S[i,count]*-1):
#                 reactants.append(sp)
#                 formula=formula+sp+"*"
#         elif S[i,count] > 0:
#             for j in range(0,S[i,count]):
#                 products.append(sp)
#         i +=1
#     if "k" not in str(line[1]):
#         # the mass-action formula
#         formula=formula[:-1]
#         #the parameter
#         paramnames.append("k"+str(count+1))
#         paramvals.append(np.double(line[1]))
#     else:
#         # specific formula (non-mass-action)
#         formula = line[1]
#         j = 1
#         params = np.array(list(line[2:])) # parameters
#         params = params[~np.isnan(params)]
#         if len(params) == 1:
#             paramnames.append("k"+str(count+1)+"_"+str(j))
#             paramvals.append(float(line[j+1]))
#             pattern = 'k\D*\d*'
#             compiled = re.compile(pattern)
#             matches = compiled.finditer(formula)
#             for ematch in matches:
#                 formula = formula.replace(ematch.group(),paramnames[-1])
#         else:
#             for p in params:
#                 paramnames.append("k"+str(count+1)+"_"+str(j))
#                 paramvals.append(float(line[j+1]))
#                 pattern1 = 'k(\D*)\d*'+'_'+str(j)
#                 compiled1 = re.compile(pattern1)
#                 matches1 = compiled1.finditer(formula)
#                 for ematch in matches1:
#                     formula = formula.replace(ematch.group(),paramnames[-1])
#                 j +=1
#     if line[0] == 'Cytoplasm':
#         valcomp = 5.25e-12
#     elif line[0] == 'Extracellular':
#         valcomp = 5.00e-5
#     elif line[0] == 'Nucleus':
#         valcomp = 1.75e-12
#     elif line[0] == 'Mitochondrion':
#         valcomp = 3.675e-13
#     fileModel.write("  %s: %s => %s; (%s)*%.6e;\n" % (Sdf.columns[count], " + ".join(reactants), " + ".join(products), formula, valcomp))
#     count +=1
# paramsSet =  pd.DataFrame(paramvals, columns = ['param_values'], index=paramnames)
# paramsSet.to_excel(fileParams)
#
#
# # Write compartment ICs
# fileModel.write("\n  // Compartment initializations:\n")
# for inds in range(lenCmps):
#     fileModel.write("  %s = %.6e;\n" % (Comps[inds], np.double(Vols[inds]))) # '%.3e'  "%.4g"
#     fileModel.write("  %s has volume;\n" % (Comps[inds]))
#
# # Write species ICs
# fileModel.write("\n  // Species initializations:\n")
# count = 0
# for sp in ICf.index:
#     vall = np.double(ICf.values[count,1])
#     fileModel.write("  %s = %.6e;\n" % (sp,vall)) # '%.3e'  "%.4g"
#     count += 1
#
# # Write parameter ICs
# fileModel.write("\n  // Parameter initializations:\n")
# count = 0
# for param in paramnames:
#     fileModel.write("  %s = %.6e;\n" % (param, np.double(paramvals[count]))) # '%.3e'  "%.4g"
#     count += 1
#
# # Write other declarations
# constantVars = ['Cytoplasm','Extracellular','Nucleus','Mitochondrion']
# # constantVars = ['Cytoplasm']
# lenConVs = len(constantVars)
#
# fileModel.write("\n  // Other declarations:\n")
# fileModel.write("  const")
# count = 1
# for teconst in constantVars:
#     if count<lenConVs:
#         fileModel.write("  %s," % (teconst))
#         count+=1
#     else:
#         fileModel.write("  %s;\n" % (teconst))
#
# # Write unit definitions
# fileModel.write("\n  // Unit definitions:")
# fileModel.write("\n  unit time_unit = second;")
# fileModel.write("\n  unit volume = litre;")
# fileModel.write("\n  unit substance = 1e-9 mole;")
# fileModel.write("\n  unit nM = 1e-9 mole / litre;")
# fileModel.write("\n")
#
# # End the model file
# fileModel.write("\nend")
# # Close the file
# fileModel.close()
#
# # To add annotations for the species, read-in the xml and Species files, then write a new xml with annotations.
# # This is because Antimony does not support annotations
#
# sbml_reader = libsbml.SBMLReader()
# sbml_doc = sbml_reader.readSBML(sbml_file)
# sbml_model = sbml_doc.getModel()
#
# # Set species annotations
# ICf = pd.read_excel(fileSpecies,header=0,index_col=0)
# for count,val in enumerate(ICf.index):
#     Annot=""
#     for col in range(3,(ICf.shape[1])):
#         aa=str(ICf.values[count,col])
#         if aa=="nan":
#             break
#         else:
#             Annot=Annot+" "+ICf.values[count,col]
#     sbml_model.getSpecies(val).setAnnotation(Annot)
#
# # Set compartment annotations
# CoVols = pd.read_excel(fileComps,header=0,index_col=0)
# for count,val in enumerate(CoVols.index):
#     sbml_model.getCompartment(val).setAnnotation(CoVols.values[count,1])
#
# # Write with the same name or use the next section instead of below lines
# writer = libsbml.SBMLWriter()
# writer.writeSBML(sbml_doc, sbml_file)
#
# # # Change model name and write with a new name
# # sbml_model.setName(model_name+'_Annot')
# # sbml_model.setId(model_name+'_Annot')
# # sbml_filewAnnot = model_name+'_Annot.xml'
# # writer = libsbml.SBMLWriter()
# # writer.writeSBML(sbml_doc, sbml_filewAnnot)
#
# # sbml_file = sbml_filewAnnot # Comment-out to use the file name  w/t annotations
# model_name = sbml_file[0:-4]
# model_output_dir = model_name
#
# sbml_reader = libsbml.SBMLReader()
# sbml_doc = sbml_reader.readSBML(sbml_file)
# sbml_model = sbml_doc.getModel()
#
# # Create an SbmlImporter instance for our SBML model
# sbml_importer = amici.SbmlImporter(sbml_file)
#
# constantParameters = [params.getId() for params in sbml_model.getListOfParameters()]
# constantParameters[-1]
#
# observables={'actp53':{'formula': 'p53ac'},'phERK':{'formula': 'ppERK'},'egfLR':{'formula': 'EE1'}}
#
# sbml_importer.sbml2amici(model_name,
#                          model_output_dir,
#                          verbose=False,
#                          observables=observables,
#                          constantParameters=constantParameters)
#
# sys.path.insert(0, os.path.abspath(model_output_dir))
# model_module = importlib.import_module(model_name)
# model = model_module.getModel()
#
# # set timepoints for which we want to simulate the model
# model.setTimepoints(np.linspace(0, 86400, 2881))
# # Create solver instance
# solver = model.getSolver()
# solver.setMaxSteps = 1e10
# # Run simulation using default model parameters and solver options
# rdata = amici.runAmiciSimulation(model, solver)
#
#
# rdata['x'][:,1]
#
# rdata['y'][:,0]
