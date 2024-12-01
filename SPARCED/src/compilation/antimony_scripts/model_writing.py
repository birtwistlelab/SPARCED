#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os

import antimony
import numpy as np

from utils.arguments import parse_args
from utils.data_handling import load_input_data_file
from compilation.antimony_scripts.antimony_write import *
from compilation.antimony_scripts.antimony_write_IC import *

def antimony_write_model(model_name: str, f_antimony: str | os.PathLike,
                         input_files: dict[str, str | os.PathLike],
                         f_output_parameters: str | os.PathLike,
                         is_SPARCED: bool) -> tuple[np.ndarray, np.ndarray]:
    """Generate an Antimony file

    Provided an empty Antimony file and paths to the SPARCED formatted input
    data files, generate the content of the Antimony file.

    Note:
        This process also creates a parameters file as an output.

    Arguments:
        f_antimony: The path of the empty Antimony file to write in.
        input_files: The dictionnary containing the input data files paths.
        f_output_parameters: The name of the output parameters file.
        is_SPARCED: A flag to raise SPARCED specific behavior or not.

    Returns:
        The contents of the compartments and species files.
    """

    with f_antimony.open(mode='w') as f:
        # Header
        if (is_SPARCED):
            f.write("# PanCancer Model by Birtwistle Lab\n")
        f.write("model {name}()\n\n".format(name=model_name))
        # Compartments
        compartments = load_input_data_file(str(input_files["compartments"]))
        antimony_write_compartments_names(f, compartments)
        # Species
        species = load_input_data_file(str(input_files["species"]))
        antimony_write_species_names(f, species)
        # Reactions
        param_names, param_vals = antimony_write_reactions(f, str(input_files["ratelaws"]),
                                  str(input_files["stoicmat"]), f_output_parameters)
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

