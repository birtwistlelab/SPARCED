#!/usr/bin/env python
# -*- coding: utf-8 -*-

from typing import IO

import numpy as np
import pandas as pd
import re

from compilation.utils.antimony_scripts.antimony_utils import load_input_data_file


def antimony_write_constant_variables(f_antimony: IO[str], constants: np.ndarray) -> None:
    """Write constant variables in the given Antimony file

    Arguments:
        f_antomony: The open Antimony file.
        constants: The constant variables to declare.

    Returns:
        Nothing.
    """

    f_antimony.write("# Other declarations:\nconst ")
    for const_var in constants[:-1]:
        f_antimony.write("{name}, ".format(name=const_var))
    f_antimony.write("{last_name};\n\n".format(last_name=constants[-1]))

def antimony_write_compartments_names(f_antimony: IO[str], compartments: np.ndarray) -> None:
    """Write compartments names in the given Antimony file

    Warning:
        The first row is considered as a header, and hence it is skipped.

    Note:
        Names should be located on the first column of the array.

    Arguments:
        f_antimony: The open Antimony file.
        compartments: Content of the input compartments file.

    Returns:
        Nothing.
    """

    f_antimony.write("# Compartments:\n")
    for i, value in enumerate(compartments[1:]):
        f_antimony.write("Compartment {name}; ".format(name=value[0]))
    f_antimony.write("\n")

def antimony_write_reactions(f: IO[str], f_ratelaws: str, f_stoichmat: str, f_outp: str):
    """Write reactions in the given Antimony file

    Warning:
        This function contains a massive copy/paste from some older code.
    Todo:
        Clean function and remove hard-coded values.

    Arguments:
        f: The open Antimony file.
        f_ratelaws: The ratelaws file path.
        f_stoichmat: The stoichiometric matrix file path.
        f_outp: The output parameters file path.

    Returns:
        A typle with the parameters' names list and the parameters' values list.
    """

    f.write("# Reactions:\n")

    stoic_sheet = load_input_data_file(f_stoichmat)
    ratelaw_sheet = load_input_data_file(f_ratelaws)
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

def antimony_write_species_names(f_antimony: IO[str], species: np.ndarray) -> None:
    """Write species names and affiliated compartments in the given Antimony file

    Warning:
        The first row is considered as a header, and hence it is skipped.

    Notes:
        Species names should be located on the first column of the array.
        Species compartments should be located on the second column of the array.

    Argurments:
        f_antimony: The open Antimony file.
        species: Content of the input species file.

    Returns:
        Nothing.
    """

    f_antimony.write("# Species:\n")
    for i, value in enumerate(species[1:]):
         f_antimony.write("Species {name} in {compartment};\n"
                          .format(name=value[0], compartment=value[1]))
    f_antimony.write("\n")

def antimony_write_unit_definitions(f_antimony: IO[str]) -> None:
    """Write unit definitions in the given Antimony file

    Warning:
        This function contains hard-coded values

    Arguments:
        f_antimony: The open Antimony file.

    Returns:
        Nothing.
    """

    f_antimony.write("# Unit definitions:\n")
    f_antimony.write("  unit time_unit = second;\n")
    f_antimony.write("  unit volume = litre;\n")
    f_antimony.write("  unit substance = 1e-9 mole;\n")
    f_antimony.write("  unit nM = 1e-9 mole / litre;\n\n")

