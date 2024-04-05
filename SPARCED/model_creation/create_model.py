#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys

import amici
from antimony import *
import libsbml

from model_creation.utils import *
from model_creation.utils.antimony_scripts.model_writing import antimony_write_model
from model_creation.utils.sbml_scripts.model_annotation import sbml_annotate_model


def create_model(antimony_model_name, sbml_model_name, f_compartments,
                 f_stoichmatrix, output_dir_path, f_output_parameters,
                 f_ratelaws, f_species, verbose):                                  

    """
    Generate Antimony, SBML and AMICI models based on given input data
    Run from parent folder (use python -m) ;
    OR remove prefix "model_creation." from import statements
    AND update default paths passed as input arguments
    
    :param antimony_model_name: desired name for the generated Antimony model
    :type antimony_model_name: [str]
    :param sbml_model_name: desired name for the generated SBML model
    :type sbml_model_name: [str]
    :param f_compartments: path to the compartments & volumes file
    :type f_compartments: [str]
    :param f_stoichmatrix: path to the stoichiometric matrix file
    :type f_stoichmatrix: [str]
    :param f_output_dir_path: path to the desired output directory
    :type f_output_dir_path: [str]
    :param f_output_parameters: path to the output parameters file
    :type f_output_parameters: [str]
    :param f_ratelaws: path to the ratelaws file
    :param f_ratelaws: [str]
    :param f_species: path to the species file
    :type f_species: [str]
    :param verbose: verbose
    :type verbose: [bool]
    :return: Nothing TODO change this
    :rtype: [void]

    """
    # ------------------------------- ANTIMONY --------------------------------
    # Create and load an Antimony model
    antimony_file_name, compartments, species = \
            antimony_write_model(antimony_model_name, output_dir_path,
                                 f_compartments, f_stoichmatrix,
                                 f_output_parameters, f_ratelaws, f_species)
    try:
        assert not loadFile(antimony_file_name) == -1
    except:
        print("SPARCED: Failed to load Antimony file")
        sys.exit(0)
    else:
        if verbose: print("SPARCED: Success loading Antimony file")
    # --------------------------------- SBML ----------------------------------
    # Convert the newly created Antimony model into an SBML model
    sbml_file_name = output_dir_path + sbml_model_name + ".xml"
    try:
        assert not writeSBMLFile(sbml_file_name, antimony_model_name) == 0
    except:
        print("SPARCED: Failed to convert Antimony file to SBML")
        sys.exit(0)
    else:
        if verbose: print("SPARCED: Success converting Antimony file to SBML")
    # Annotate the SBML model
    sbml_annotate_model(sbml_file_name, species, compartments)
    sys.exit(2)
    # --------------------------------- AMICI ---------------------------------
    # Import
    sys.path.insert(0, os.path.abspath(sbml_model_name))
    sbml_reader = libsbml.SBMLReader()
    sbml_doc = sbml_reader.readSBML(sbml_file_name)
    sbml_model = sbml_doc.getModel()
    sbml_importer = amici.SbmlImporter(sbml_file_name)
    const_params = [params.getId() for params in sbml_model.getListOfParameters()]
    # Compile
    model_output_dir = sbml_model_name
    # Need to add constant parameters and observable's script here
    sbml_importer.sbml2amici(sbml_model_name, model_output_dir, verbose=args.verbose)
    if verbose: print("SPARCED: Sucess compiling the model")

def launch_model_creation():
    """
    Small routine to process parsed arguments and call create_model
    """
    args = parse_args()                                                         
    # Add path of input directory to input files names                          
    f_compartments = args.inputdir + args.compartments                          
    f_stoichmatrix = args.inputdir + args.stoichmatrix                          
    f_output_params = args.inputdir + args.outputparams                         
    f_ratelaws = args.inputdir + args.ratelaws                                  
    f_species = args.inputdir + args.species                                    
    # Create model                                                              
    create_model(args.antimony, args.sbml, f_compartments, f_stoichmatrix,      
                 args.outputdir, f_output_params, f_ratelaws, f_species,        
                 args.verbose)


if __name__ == '__main__':
    launch_model_creation()

