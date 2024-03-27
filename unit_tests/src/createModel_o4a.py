#!/usr/bin/env python
# coding: utf-8

# # SPARCED Model Creation (Integrated-SBML)
# 
# This file transforms the input files into Antimony and SBML model files. Here, we do not execute any simulation. 

# **These are the "default" input file names**   
# They are structured, tab-separated text files, located in the "input files" folder by default.

# In[1]:


# Input file name definitions
fnameInput = 'SPARCEDo4a_v1' # model file name = USER input

fileComps = 'Compartments.txt' # input
fileSpecies = 'Species.txt' # input
fileStoic = 'StoicMat.txt' # input
fileRatelaws = 'Ratelaws.txt' # input
fileOmicsdata = 'OmicsData.txt' # input
fileGeneReg = 'GeneReg.txt' # input
fileObservables = 'Observables.txt' # input
fileParamsOut = 'ParamsAll_v1.txt' # output: Lists all parameter names, rxn names, values


# ### Import required packages and scripts

# In[2]:


import sys
import os
import pandas as pd



cd = os.getcwd()
## Get the directory path
wd = os.path.dirname(os.path.abspath(__file__))

# Ensure the SPARCED root and bin directories are in the system path
sparced_root = '/'.join(wd.split(os.path.sep)[:wd.split(os.path.sep).index('SPARCED')+1])
sys.path.append(os.path.join(sparced_root, 'bin'))

input_path = os.path.join(sparced_root +'/input_files')

from modules.copyDir import copyDirectory
copyDirectory(input_path, os.getcwd()+"/")


import libsbml
import importlib
import amici
import numpy as np
import re
import pandas as pd
from antimony import *


# Optional packages to import
import amici.plotting
import matplotlib.pyplot as plt


# ### Define the model file and SBML file name here
# 
# First, the Antimony model file will be created using the input files. Users can add information per the antimony file structure. A default explanation for the SPARCED model is inserted here.
# Then the SBML file will be created using the Antimony file. Its name can also be defined here. 
# For consistency, we use the same name for the Antimony file, SBML file, and the model name.

# In[3]:


fileModel = open(os.path.join(cd,fnameInput+'.txt'),'w') 
fileModel.write("# PanCancer Model by Birtwistle Lab \n") # some explanation
fileModel.write("model "+fnameInput+"()\n") # model name


# SBML file name
sbml_file = fnameInput+'.xml'
# Name of the model that will also be the name of the python module
model_name = sbml_file[0:-4] 
# Directory to which the generated model code is written
model_output_dir = model_name 
# The AMICI package will create this folder while compiling the model, refer to AMICI github page for more details


# ## Input file processing
# 
# Below are steps of input file processing and writing the Antimony file.

# ### First, read-in the "Compartments" and "Volumes"

# In[4]:


#%%
# Initializing compartment and volume lists
compartments = []
volumes = []

# Create/write compartments
compartment_sheet = np.array([np.array(line.strip().split("\t")) for line in open(os.path.join(input_path,fileComps))])

#read in each line minus the header row of compartments file
for row in compartment_sheet[1:]:
    compartments.append(row[0])
    volumes.append(row[1])
    
#

fileModel.write("\n  # Compartments and Species:\n") # Antimony Compartments/Species module title
for idx in range(len(compartments)):
    compName = compartments[idx]
    fileModel.write("  Compartment %s; " % (compName))
fileModel.write("\n")


# ### Second, read-in the "Species" and related information
# 

# In[5]:


model_genes = list(pd.read_csv(os.path.join(input_path,fileOmicsdata),sep='\t',index_col=0,header=0).index)
sp_input = pd.read_csv(os.path.join(input_path,fileSpecies),sep='\t',index_col=None,header=None)
sp_mrna = ['m_'+str(m) for m in model_genes]

sp_input = sp_input.loc[~np.isin(sp_input.iloc[:,0],sp_mrna),:]
species_sheet = [np.array(sp_input.iloc[s,~sp_input.iloc[s,:].isnull().values]).astype(str) for s in range(np.shape(sp_input)[0])]

#%%# Write species and assign compartments

#read in each line minus the header row of species file
species_compartments = [] # Create a list of "home" compartments for each species in the model
for row in species_sheet[1:]:
    species_compartments.append(row[1]) 
species_compartments = np.array(species_compartments)

#Write each species to model txt file
fileModel.write("\n")
for idx,val in enumerate(species_sheet[1:]):
    fileModel.write("  Species ")
    fileModel.write("%s in %s" % (val[0], species_compartments[idx]))
    fileModel.write(";\n")
    


omicsRead = pd.read_csv(os.path.join(input_path,fileOmicsdata),header=0,index_col=0,sep="\t")
gExp_mpc = np.float64(omicsRead.values[:,0])
mExp_mpc = np.float64(omicsRead.values[:,1])
kGin = np.float64(omicsRead.values[:,2])
kGac = np.float64(omicsRead.values[:,3])
kTCleak = np.float64(omicsRead.values[:,4])
kTCmaxs = np.float64(omicsRead.values[:,5])
kTCd = np.float64(omicsRead.values[:,6])
gnames = list(omicsRead.index.values)
ac_genenames = ['ag_' + x for x in gnames]
ia_genenames = ['ig_' + x for x in gnames]
mRNAnames = ['m_' + x for x in gnames]

Vn = np.double(compartment_sheet[3,1])
Vc = np.double(compartment_sheet[1,1])
mpc2nMVc = (1E9/(Vc*6.023E+23))


genesComp = np.array(omicsRead.values[:,12])
# CC mRNA IC exceptions:


cc_mrna = pd.read_csv(os.path.join(input_path,'Initializer.txt'),sep='\t',squeeze=True,usecols=['Step1_mrna','Step1_mrna_mpc'],index_col='Step1_mrna')
cc_mrna = cc_mrna[cc_mrna.index.notnull()]

for c in range(len(cc_mrna)):
    cmrna_sp = 'm_'+str(cc_mrna.index[c])
    mExp_mpc[list(mRNAnames).index(cmrna_sp)] = cc_mrna[c]

##

for idx, val in enumerate(ac_genenames):
    fileModel.write("  Species ")
    fileModel.write("%s in %s" % (val, genesComp[idx]))
    fileModel.write(";\n")
for idx, val in enumerate(ia_genenames):
    fileModel.write("  Species ")
    fileModel.write("%s in %s" % (val, genesComp[idx]))
    fileModel.write(";\n")
for idx, val in enumerate(mRNAnames):
    fileModel.write("  Species ")
    fileModel.write("%s in %s" % (val, genesComp[idx]))
    fileModel.write(";\n")


# ### Read-in the "Reactions" and "StoicMat" files to create reactions of the model

# In[6]:


# Write reactions
fileModel.write("\n\n  # Reactions:\n") # Antimony Reactions module title

#%% Create the reations of the model based on the stoichiometric input and (if available) ratelaw input
stm_input = pd.read_csv(os.path.join(input_path,fileStoic),sep='\t',index_col=None,header=None)
stm_input = stm_input.iloc[~np.isin(stm_input.iloc[:,0],sp_mrna),:]

stoic_sheet = [np.array(stm_input.iloc[s,~stm_input.iloc[s,:].isnull().values]).astype(str) for s in range(np.shape(stm_input)[0])]


#creates associated ratelaw data list
ratelaw_sheet = np.array([np.array(line.strip().split("\t")) for line in open(os.path.join(input_path,fileRatelaws))])
ratelaw_data = np.array([line[1:] for line in ratelaw_sheet[1:]])

#gets first column minus blank space at the beginning, adds to stoic data list
stoic_columnnames = stoic_sheet[0]
stoic_rownames = [line[0] for line in stoic_sheet[1:]]
stoic_data = np.array([line[1:] for line in stoic_sheet[1:]])


# #### Create the reations of the model based on the stoichiometric input and (if available) ratelaw input

# In[7]:


# builds the important ratelaw+stoic lines into the txt file 
paramnames = []
paramvals = []
paramrxns = []
paramidxs = []
for rowNum, ratelaw in enumerate(ratelaw_data):
    reactants = []
    products = []
    formula = "k"+str(rowNum+1)+"*"

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
        formula = formula[:-1]
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

    if products == [] and reactants == []:
        pass
    else:
        fileModel.write("  %s: %s => %s; (%s)*%s;\n" % (stoic_columnnames[rowNum], " + ".join(reactants), " + ".join(products), formula, ratelaw[0]))


rowLastNo = rowNum+1
for rowNum, val in enumerate(ia_genenames):

    formula = "k"+str(rowLastNo+rowNum+1)+"*"+val
    fileModel.write("  %s: %s => %s; (%s);\n" % ('vGac'+str(rowNum+1), val, ac_genenames[rowNum], formula))
    
    paramnames.append("k"+str(rowLastNo+rowNum+1))
    paramvals.append(kGac[rowNum])

    paramrxns.append('vGac'+str(rowNum+1))
    paramidxs.append(0)

rowLastNo = rowLastNo+rowNum+1
for rowNum, val in enumerate(ac_genenames):

    formula = "k"+str(rowLastNo+rowNum+1)+"*"+val
    fileModel.write("  %s: %s => %s; (%s);\n" % ('vGin'+str(rowNum+1), val, ia_genenames[rowNum], formula))
    
    paramnames.append("k"+str(rowLastNo+rowNum+1))
    paramvals.append(kGin[rowNum])

    paramrxns.append('vGin'+str(rowNum+1))
    paramidxs.append(0)
    
#

TARsRead = pd.read_csv(os.path.join(input_path,fileGeneReg),header=0,index_col=0,sep="\t")
TARs0 = (TARsRead.values)
numberofTARs = len(TARsRead.columns)
numberofgenes = int(len(gExp_mpc))

tcnas = np.ones((numberofgenes, numberofTARs))
tck50as = np.zeros((numberofgenes, numberofTARs))
tcnrs = np.ones((numberofgenes, numberofTARs))
tck50rs = np.zeros((numberofgenes, numberofTARs))
for qq in range(numberofgenes):
    for ww in range(numberofTARs):
        pars = str(TARs0[qq,ww]).find(';')
        if pars>0:
            nH = np.float(TARs0[qq,ww][0:pars])
            kH = np.float(TARs0[qq,ww][pars+2::])
            if nH>0:
                tcnas[qq,ww] = nH
                tck50as[qq,ww] = kH
            else:
                tcnrs[qq,ww] = abs(nH)
                tck50rs[qq,ww] = kH
mpc2nmcf_Vn = 1.0E9/(np.float(Vn)*6.023E+23)
mpc2nmcf_Vc = 1.0E9/(np.float(Vc)*6.023E+23)

mrnas2mod = ["m_CCND1","m_CCND2","m_CCND3"]


rowLastNo = rowLastNo+rowNum+1
for rowNum, val in enumerate(mRNAnames):

    tempTAs = []
    tempnAs = []
    tempKAs = []
    tempTRs = []
    tempnRs = []
    tempKRs = []

    for qq in range(numberofTARs):
        if tck50as[rowNum,qq] > 0:
            tempTAs.append(TARsRead.columns[qq])
            tempnAs.append(tcnas[rowNum,qq])
            tempKAs.append(tck50as[rowNum,qq])
        if tck50rs[rowNum,qq] > 0:
            tempTRs.append(TARsRead.columns[qq])
            tempnRs.append(tcnrs[rowNum,qq])
            tempKRs.append(tck50rs[rowNum,qq])


    j = 3
    formulaAs = []
    for pp in range(np.size(tempTAs)):
        formulaAs.append("("+tempTAs[pp]+"/"+"k"+str(rowLastNo+rowNum+1)+"_"+str(j)+")^"+"k"+str(rowLastNo+rowNum+1)+"_"+str(j+1))
        paramnames.append("k"+str(rowLastNo+rowNum+1)+"_"+str(j))
        paramvals.append(tempKAs[pp])
        paramrxns.append('vTC'+str(rowNum+1))
        paramidxs.append(j)
        paramnames.append("k"+str(rowLastNo+rowNum+1)+"_"+str(j+1))
        paramvals.append(tempnAs[pp])
        paramrxns.append('vTC'+str(rowNum+1))
        paramidxs.append(j+1)       
        j = j + 2
    if val not in mrnas2mod: # if not Cyclin D gene
        formulaAs = " + ".join(formulaAs)


    
    formulaRs = []
    for pp in range(np.size(tempTRs)):
        formulaRs.append("("+tempTRs[pp]+"/"+"k"+str(rowLastNo+rowNum+1)+"_"+str(j)+")^"+"k"+str(rowLastNo+rowNum+1)+"_"+str(j+1))
        paramnames.append("k"+str(rowLastNo+rowNum+1)+"_"+str(j))
        paramvals.append(tempKRs[pp])
        paramrxns.append('vTC'+str(rowNum+1))
        paramidxs.append(j)
        paramnames.append("k"+str(rowLastNo+rowNum+1)+"_"+str(j+1))
        paramvals.append(tempnRs[pp])
        paramrxns.append('vTC'+str(rowNum+1))
        paramidxs.append(j+1)  
        j = j + 2
    formulaRs = " + ".join(formulaRs)       

    if val in mrnas2mod: # if Cyclin D gene            
        formulaARs = formulaAs[0]+"/(1.0+"+formulaAs[0]+")"+"*"+formulaAs[1]+"/(1.0+"+formulaAs[1]+")"
    elif formulaAs and formulaRs:
        formulaARs = "("+"("+formulaAs+")"+"/(1.0+"+formulaAs+"+"+formulaRs+"))"
    elif formulaAs:
        formulaARs = "("+"("+formulaAs+")/(1.0+"+formulaAs+"))"
    elif formulaRs:
        formulaARs = "(1.0"+"/(1.0+"+formulaRs+"))"
    else:
        formulaARs = "0.0"

    j = 1


    formula = "("+"k"+str(rowLastNo+rowNum+1)+"_"+str(j)+"*"+ac_genenames[rowNum]+")+"+"("+"k"+str(rowLastNo+rowNum+1)+"_"+str(j+1)+"*"+ac_genenames[rowNum]+"*"+formulaARs+")"

    fileModel.write("  %s: %s => %s; (%s)*%s;\n" % ('vTC'+str(rowNum+1), "", val, formula, genesComp[rowNum])) # for vTC


    paramnames.append("k"+str(rowLastNo+rowNum+1)+"_"+str(j))
    paramvals.append(kTCleak[rowNum])
    paramrxns.append('vTC'+str(rowNum+1))
    paramidxs.append(1)
    paramnames.append("k"+str(rowLastNo+rowNum+1)+"_"+str(j+1))
    paramvals.append(kTCmaxs[rowNum])
    paramrxns.append('vTC'+str(rowNum+1))
    paramidxs.append(2)  

rowLastNo = rowLastNo+rowNum+1
for rowNum, val in enumerate(mRNAnames):


    formula = "("+"k"+str(rowLastNo+rowNum+1)+"*"+val+")"

    fileModel.write("  %s: %s => %s; (%s)*%s;\n" % ('vTCd'+str(rowNum+1), val, "", formula, genesComp[rowNum])) # for vTCd
    paramnames.append("k"+str(rowLastNo+rowNum+1))
    paramvals.append(kTCd[rowNum])
    paramrxns.append('vTCd'+str(rowNum+1))
    paramidxs.append(0)
    
# Export parameters for each reaction, with corresponding order within the ratelaw and its value
params_all = pd.DataFrame({'value':paramvals,'rxn':paramrxns,'idx':paramidxs},index=paramnames)
params_all.to_csv(os.path.join(cd,fileParamsOut),sep='\t',header=True, index=True)



# ### Next, we set the "Initial Conditions" for Compartments, Species, and Parameters

# In[8]:


#% The Compartment initial conditions

# Write compartment ICs
fileModel.write("\n  # Compartment initializations:\n")
for idx in range(len(compartments)):
    fileModel.write("  %s = %.6e;\n" % (compartments[idx], np.double(volumes[idx])))
    fileModel.write("  %s has volume;\n" % (compartments[idx]))
    
# The Species initial conditions
# Write species ICs
fileModel.write("\n  # Species initializations:\n")
for idx, val in enumerate(species_sheet[1:]):
    fileModel.write("  %s = %.6e;\n" % (val[0],np.double(val[2])))
    

xgac_mpc_D = (kGac*gExp_mpc)/(kGin+kGac) #active genes initial condition
xgin_mpc_D = gExp_mpc - xgac_mpc_D

genesIC = np.concatenate((xgac_mpc_D,xgin_mpc_D))

#%%



for idx, val in enumerate(ac_genenames):
    fileModel.write("  %s = %.6e;\n" % (val,np.double(genesIC[idx]*mpc2nmcf_Vc)))
for idx, val in enumerate(ia_genenames):
    fileModel.write("  %s = %.6e;\n" % (val,np.double(genesIC[idx+len(model_genes)]*mpc2nmcf_Vc)))
for idx, val in enumerate(mRNAnames):
    fileModel.write("  %s = %.6e;\n" % (val,np.double(mExp_mpc[idx]*mpc2nmcf_Vc)))
    
# The Parameter (of reactions) initial conditions

# Write parameter ICs
fileModel.write("\n  # Parameter initializations:\n")
count = 0
for param in paramnames:
    fileModel.write("  %s = %.6e;\n" % (param, np.double(paramvals[count])))
    count += 1



# ### Other declarations supported by Antimony

# In[9]:


#%% Other declarations supported by Antimony
# Write other declarations
constantVars = compartments

fileModel.write("\n  # Other declarations:\n")
fileModel.write("  const")
for constVar in constantVars[:-1]:
    fileModel.write("  %s," % (constVar))
#last item in row needs semicolon
fileModel.write("  %s;\n" % (constantVars[-1]))


# ### The Unit Definitions
# 
# We define the time, volume, and substance units

# In[10]:


# Write unit definitions
fileModel.write("\n  # Unit definitions:")
fileModel.write("\n  unit time_unit = second;")
fileModel.write("\n  unit volume = litre;")
fileModel.write("\n  unit substance = 1e-9 mole;")
fileModel.write("\n  unit nM = 1e-9 mole / litre;")
fileModel.write("\n")


# ### End of the Antimony file, close the file

# In[11]:


# End the model file
fileModel.write("\nend")
# Close the file
fileModel.close()


# ### The Antimony file import and conversion to SBML format

# In[12]:


# load model and convert to SBML
if loadFile(os.path.join(cd,fnameInput+".txt")) == 1:
    print("Success loading antimony file")
else:
    print("Failed to load antimony file")
    exit(1)

if writeSBMLFile(os.path.join(cd,fnameInput+".xml"),fnameInput) == 1:
    print("Success converting antimony to SBML")
else:
    print("Failure converting antimony to SBML")


# #### The initial SBML file is imported again, to add annotations for the Compartments and Species
# 
# Currently, the Antimony file format does not support Annotation definitions. So, we re-process the SBML file and add Annotations for the Compartments (GO terms) and Species (HGNC or ENSEMBL identifiers).

# ### Import the SBML file and get the model handle

# In[13]:


# create interaction components
sbml_reader = libsbml.SBMLReader()
sbml_doc = sbml_reader.readSBML(os.path.join(cd,sbml_file))
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
    
#
# Set compartment annotations
for row in compartment_sheet[1:]:
    sbml_model.getCompartment(row[0]).setAnnotation(row[2])
    
# Write with the same name or use the next section instead of below lines
writer = libsbml.SBMLWriter()
writer.writeSBML(sbml_doc, os.path.join(cd,sbml_file))



# ## The Model Compilation

# ### Import the SBML file and compile using the AMICI package

# In[14]:


# prepares to use interaction components to synthesize model
sys.path.insert(0, os.path.abspath(model_output_dir))
model_name = sbml_file[0:-4]
model_output_dir = model_name

sbml_reader = libsbml.SBMLReader()
sbml_doc = sbml_reader.readSBML(os.path.join(cd,sbml_file))
sbml_model = sbml_doc.getModel()

# Create an SbmlImporter instance for our SBML model
sbml_importer = amici.SbmlImporter(os.path.join(cd,sbml_file))


# #### Setting AMICI converter arguments
# 
# We set all the rate parameters as constant for faster compilation. It is also a good idea to define "Observables" for easier comparison to experimental data. An observable in our model corresponds to the totality of a protein sepcies, including all of its post-translationally modified and unmodified states, corrected for compartmental volume differences and reported based on cytoplasmic volume. We defined 102 observables using the "Observables.txt" input file in `AMICI` required format.

# In[15]:


#sets important constants for model build
constantParameters = [params.getId() for params in sbml_model.getListOfParameters()]

# creates observables.

obs_input = pd.read_csv(os.path.join(input_path,fileObservables),sep='\t',header=0,index_col=0)
ObsMat = obs_input.loc[~np.isin(obs_input.index,sp_mrna),:]

# ObsMat = pd.read_csv(os.path.join(wd,'input_files',fileObservables), sep='\t',header=0, index_col=0)
Vc = float(compartment_sheet[compartment_sheet[:,0]=='Cytoplasm',1])

species_names = np.array([species_sheet[i][0] for i in range(1,len(species_sheet))])
Vol_species = np.array([species_sheet[i][1] for i in range(1,len(species_sheet))])
Vol_species = [float(compartment_sheet[compartment_sheet[:,0]==Vol_species[i],1][0]) for i in range(len(Vol_species))]
Vol_species = pd.Series(Vol_species, index=species_names)

formula_obs = []
for obs in ObsMat.columns:
    sp_obs = ObsMat.index[np.nonzero(np.array(ObsMat.loc[:,obs]>0))[0]]
    sp_obs_id = np.nonzero(np.array(ObsMat.loc[:,obs]>0))[0]
    Vr = Vol_species/Vc
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



# #### Compile the model and create required C++ dependencies for simulation. 
# 
# Creates a sub-folder with the same name as the model (i.e. SPARCED)

# In[16]:


# The actual compilation step by AMICI, takes a while to complete for large models
sbml_importer.sbml2amici(model_name,
                         os.path.join(cd,model_output_dir),
                         verbose=False)
                        #  observables=observables,
                        #  constant_parameters=constantParameters)


# ## The model creation is now complete! Enjoy...
