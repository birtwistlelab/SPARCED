#!/usr/bin/env python
# -*- coding: utf-8 -*-

import antimony
import numpy as np

from compilation.utils.arguments import parse_args
from compilation.utils.antimony_scripts.antimony_utils import *
from compilation.utils.antimony_scripts.antimony_write import *
from compilation.utils.antimony_scripts.antimony_write_IC import *

def antimony_write_model(f_antimony: str, f_compartments: str,
                         f_stoichmatrix: str, f_output_parameters: str,
                         f_ratelaws: str, f_species: str, is_SPARCED: bool) \
                         -> tuple[np.ndarray, np.ndarray]:

    """Generate an Antimony file

    Provided an empty Antimony file and some SPARCED formatted data on
    compartments, stoichiometric matrix, ratelaws, and species, generate the
    content of the Antimony file.

    Note:
    This process also creates a parameters file as an output.
    Arguments:
        f_antimony: The name of the empty Antimony file to write in.
        f_compartments: The name of the input compartments file.
        f_stoichmatrix: The name of the input schiochiometric matrix file.
        f_output_parameters: The name of the output parameters file.
        f_ratelaws: The name of the input ratelaws file.
        f_species: The name of the input species file.
        is_SPARCED: A flag to raise SPARCED specific behavior or not.
    Returns:
        The contents of the compartments and species files.
    """

    with open(f_antimony,"w") as f:
        # Header
        if (is_SPARCED):
            f.write("# PanCancer Model by Birtwistle Lab\n")
        model_name = extract_antimony_model_name(f_antimony)
        f.write("model {antimony}()\n\n".format(antimony=model_name))
        # Compartments
        compartments = load_input_data_file(f_compartments)
        antimony_write_compartments_names(f, compartments)
        # Species
        species = load_input_data_file(f_species)
        antimony_write_species_names(f, species)
        # Reactions
        param_names, param_vals = antimony_write_reactions(f, f_ratelaws,
                                  f_stoichmatrix, f_output_parameters)
        # Initial conditions
        antimony_write_compartments_IC(f, compartments)
        antimony_write_species_IC(f, species)
        antimony_write_reactions_IC(f, param_names, param_vals)
        # Declare compartments as constant variables
        if (is_SPARCED):
            antimony_write_constant_variables(f, compartments[:,0][1:])
        # Unit definitions
        antimony_write_unit_definitions(f)
        f.write("\nend") 
    return(compartments, species)


if __name__ == '__main__':
    args = parse_args()
    # Add path of input directory to input files names
    antimony_model_name = args.outputdir + "ant_" + args.name + ".txt"
    f_compartments = args.inputdir + args.compartments
    f_stoichmatrix = args.inputdir + args.stoichmatrix
    f_output_params = args.inputdir + args.outputparams
    f_ratelaws = args.inputdir + args.ratelaws
    f_species = args.inputdir + args.species
    # Write model
    antimony_write_model(antimony_model_name, f_compartments, f_stoichmatrix,
                         f_output_params, f_ratelaws, f_species)
