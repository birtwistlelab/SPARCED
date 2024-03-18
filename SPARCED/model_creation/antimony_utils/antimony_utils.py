#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import re

def antimony_init(f_cv, f_s):
    """Load compartments, volumes and species lists from input files

    Args:
        f_cv: the name of the compartments (and volumes) file as a string
        f_s: the name of the species file as a string

    Returns:
        A tuple with the compartments, volumes and species lists.
    """
    # Compartments and Volumes
    comp = []
    vol = []
    sheet = np.array([np.array(line.strip().split("\t")) for line in open(f_cv)])
    for row in sheet[1:]: # Skip header row
        comp.append(row[0])
        vol.append(row[1])
    # Species
    spec = np.array([np.array(line.strip().split("\t")) for line in open(f_s)], dtype="object")
    return((comp, vol, spec))

def antimony_terminal(f):
    """Write constant variables and unit definitions in the given Antimony file
    WARNING: This function contains hard-coded values

    Args:
        f: the Antimony file to write in as an open (w) file

    Returns:
        None.
    """
    # Constant variables
    f.write("# Other declarations:\nconst ")
    constants = ['Cytoplasm', 'Extracellular', 'Nucleus', 'Mitochondrion']
    for c in constants[:-1]:
        f.write("{name}, ".format(name=c))
    f.write("{last_name}\n\n".format(last_name=constants[-1]))
    # Unit definitions
    f.write("# Unit definitions:\n")
    f.write("  unit time_unit = second;\n")
    f.write("  unit volume = litre;\n")
    f.write("  unit substance = 1e-9 mole;\n")
    f.write("  unit nM = 1e-9 mole / litre;\n\n")

def antimony_write_compartments(f, comp):
    """Write compartments' names in the given Antimony file

    Args:
        f: the Antimony file to write in as an open (w) file
        comp: the compartments' names as a list of strings

    Returns:
        None.
    """
    f.write("# Compartments and Species:\n")
    for i in range(len(comp)):
        f.write("Compartment {comp_name}; ".format(comp_name=comp[i]))
    f.write("\n\n")

def antimony_write_init_compartments(f, comp, vol):
    """Write compartments' initial conditions in the given Antimony file

    Args:
        f: the Antimony file to write in as an open (w) file
        comp: the compartments' names as a list of strings
        vol: the compartments' volumes as a list of numericals

    Returns:
        None.
    """
    f.write("# Compartments initialization:\n")
    for i in range(len(comp)):
        f.write("{name} = {volume:.6e};\n{name} has volume;\n".format(name=comp[i], volume=np.double(vol[i])))
    f.write("\n")

def antimony_write_init_reactions(f, p_names, p_vals):
    """Write reactions' parameters' initial conditions in the given Antimony file

    Args:
        f: the Antimony file to write in as an open (w) file
        p_names: the parameters' names as a list of strings
        p_vals: the parameters' values as a list of numericals

    Returns:
        None.
    """
    f.write("# Reactions' parameters initializations:\n ")
    for i, val in enumerate(p_names):
        f.write("{name} = {value:.6e};\n".format(name=val, value=np.double(p_vals[i])))
    f.write("\n")

def antimony_write_init_species(f, spec):
    """Write species' initial concentrations in the given Antimony file

    Args:
        f: the Antimony file to write in as an open (w) file
        spec: the species' names and concentrations as a Numpy array

    Returns:
        None.
    """
    f.write("# Species initialization:\n")
    for i, val in enumerate(spec[1:]): # Skip header
        f.write("{name} = {concentration:.6e};\n".format(name=val[0], concentration=np.double(val[2])))
    f.write("\n")

def antimony_write_reactions(f, f_rl, f_sm, f_outp):
    """Write reactions in the given Antimony file
    WARNING: This function contains a massive copy/paste from some older code

    Args:
        f: the Antimony file to write in as an open (w) file
        f_rl: the name of the ratelaws file as a string
        f_sm: the name of the stoichiometric matrix file as a string
        f_outp: the name of the output parameters file as a string

    Returns:
        A typle with the parameters' names list and the parameters' values list.
    """
    f.write("# Reactions:\n")
    stoic_sheet = np.array([np.array(line.strip().split("\t")) for line in open(f_sm)], dtype="object")
    ratelaw_sheet = np.array([np.array(line.strip().split("\t")) for line in open(f_rl)], dtype="object")
    ratelaw_data = np.array([line[1:] for line in ratelaw_sheet[1:]], dtype="object")
    # ========== COPY/PASTE ==========
    #gets first column minus blank space at the beginning, adds to stoic data list
    stoic_columnnames = stoic_sheet[0]
    stoic_rownames = [line[0] for line in stoic_sheet[1:]]
    stoic_data = np.array([line[1:] for line in stoic_sheet[1:]])
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
            f.write("  %s: %s => %s; (%s)*%.6e;\n" % (stoic_columnnames[rowNum], " + ".join(reactants), " + ".join(products), formula, valcomp))
    
    # Export parameters for each reaction, with corresponding order within the ratelaw and its value
    params_all = pd.DataFrame({'value':paramvals,'rxn':paramrxns,'idx':paramidxs},index=paramnames)
    params_all.to_csv(f_outp,sep='\t',header=True, index=True)
    # ========== END OF COPY/PASTE ==========
    f.write("\n")
    return((paramnames, paramvals))

def antimony_write_species(f, spec):
    """Write species' names and compartments in the given Antimony file

    Args:
        f: the Antimony file to write in as an open (w) file
        spec: the species' names and concentrations as a list of Numpy array

    Returns:
        None.
    """
    for i, val in enumerate(spec[1:]): # Skip header row
        f.write("Species {name} in {compartment};\n".format(name=val[0], compartment=val[1]))
    f.write("\n")
