#!/usr/bin/env python
# -*- coding: utf-8 -*-

import amici
import libsbml


def compile_sbml_to_amici(sbml_file_name, model_name, output_dir, verbose):
    # Import SBML file
    sbml_reader = libsbml.SBMLReader()
    sbml_doc = sbml_reader.readSBML(sbml_file_name)
    sbml_model = sbml_doc.getModel()
    sbml_importer = amici.SbmlImporter(sbml_file)
    # Set all the parameters as constant
    constant_parameters = [params.getId() for params in sbml_model.getListOfParameters()]
    # Set observables
    # TODO
    # Compile
    model_output_dir = output_dir + "amici_" + model_name
    sbml_importer.sbml2amici(model_name, model_output_dir, verbose=verbose,
                             observables=observables, constantParameters=constant_parameters)

