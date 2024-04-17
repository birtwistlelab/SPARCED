#!/usr/bin/env python
# -*- coding: utf-8 -*-

import amici
import libsbml
import numpy as np

from compilation.utils.amici_scripts.amici_utils import *


def compile_sbml_to_amici(sbml_file_name: str, model_output_dir: str,
                          f_observables: str, compartments: np.ndarray,
                          species: np.ndarray) -> None:
    """Compile the given SBML model into an AMICI model

    Warning:
        Passing observables and constant parameters to the importer is
        currently broken, in the original script too.

    Arguments:
        sbml_file_name: The SBML file path.
        model_output_dir: Output directory path for the generated AMICI model.
        f_observables: The observables input file path.
        compartments: Content of the compartments input file.
        species: Content of the species input file.

    Returns:
        Nothing.
    """

    # Import SBML file
    sbml_reader = libsbml.SBMLReader()
    sbml_doc = sbml_reader.readSBML(sbml_file_name)
    sbml_model = sbml_doc.getModel()
    sbml_importer = amici.SbmlImporter(sbml_file_name)
    # Set all the rate parameters as constant (faster compilation)
    constant_parameters = [parameters.getId() \
                           for parameters in sbml_model.getListOfParameters()]
    # Create observables
    # observables = define_observables(f_observables, compartments, species)
    # Compile
    model_name = extract_amici_model_name(model_output_dir)
    sbml_importer.sbml2amici(model_name, model_output_dir, verbose=verbose)
            # observables=observables, constantParameters=constant_parameters)

