#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys

import amici
from antimony import *
import libsbml

from compilation.amici_scripts.model_compilation import compile_sbml_to_amici
from compilation.antimony_scripts.model_writing import antimony_write_model
from compilation.sbml_scripts.model_annotation import sbml_annotate_model
from utils.paths_handling import append_subfolder

def create_model(model_name: str, model_path: str | os.PathLike,
                 input_files: dict[str, str | os.PathLike],
                 output_parameters: str | os.PathLike, verbose: bool,
                 is_SPARCED: bool) -> None:
    """
    Generate Antimony, SBML and AMICI models based on given input data files.

    Warning:
        The model name should not contain any dash ('-') characters as those
        are not supported by Antimony.

    Arguments:
        model_name: The name of the model.
        model_path: The path to the model folder.
        input_files: A dictionnary containing paths towards input data file.
        output_parameters: The desired name for the output parameters file.
        verbose: Verbose.
        is_SPARCED: Use SPARCED hard-coded values/behaviors.

    Returns:
        Nothing.
    """
    
    # ------------------------------- ANTIMONY --------------------------------
    # Create and load an Antimony model
    antimony_model_name = "ant_" + model_name + ".txt"
    antimony_file_path = append_subfolder(model_path, antimony_model_name)
    compartments, species = antimony_write_model(model_name, antimony_file_path,
                                                 input_files, output_parameters,
                                                 is_SPARCED)
    try:
        assert not loadFile(str(antimony_file_path)) == -1
    except:
        print("{name}: Failed to load Antimony file".format(name=model_name))
        sys.exit(0)
    else:
         if verbose: print("{name}: Success loading Antimony file"
                         .format(name=model_name))
    # --------------------------------- SBML ----------------------------------
    # Convert the newly created Antimony model into an SBML model
    sbml_file_name = "sbml_" + model_name + ".xml"
    sbml_file_path = append_subfolder(model_path, sbml_file_name)
    try:
        assert not writeSBMLFile(str(sbml_file_path), model_name) == 0
    except:
        print("{name}: Failed to convert Antimony file to SBML"
             .format(name=model_name))
        sys.exit(0)
    else:
        if verbose: print("{name}: Success converting Antimony file to SBML"
                         .format(name=model_name))
    # Annotate the SBML model
    if is_SPARCED:
        sbml_annotate_model(str(sbml_file_path), species, compartments)
    # --------------------------------- AMICI ---------------------------------
    # Compile the SBML model into an AMICI model
    amici_folder_name = "amici_" + model_name
    amici_folder_path = append_subfolder(model_path, amici_folder_name)
    compile_sbml_to_amici(str(sbml_file_path), str(amici_folder_path),
                          input_files["observables"], compartments, species, verbose)
    if verbose: print("{name}: Sucess compiling the model"
                     .format(name=model_name))

