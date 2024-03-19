#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys

import amici
from antimony import *
import argparse
import libsbml

from antimony_utils import *
# from src.copydir import copy_directory
# from src.sbml_utils import *


def create_model(antimony_model_name,sbml_model_name,f_comp,f_stoi,f_outp,f_rate,f_spec,verbose):
    # Create and load Antimony model
    antimony_file_name = antimony_write_model(antimony_model_name,f_comp,f_stoi,f_outp,f_rate,f_spec)
    try:
        assert not loadFile(antimony_file_name) == -1
    except:
        print("SPARCED: Failed to load Antimony file")
        sys.exit(0)
    else:
        if args.verbose: print("SPARCED: Success loading Antimony file")
    # Convert Antimony model into SBML
    sbml_file_name = sbml_model_name + ".xml"
    try:
        assert not writeSBMLFile(sbml_file_name, antimony_model_name) == 0
    except:
        print("SPARCED: Failed to convert Antimony file to SBML")
        sys.exit(0)
    else:
        if args.verbose: print("SPARCED: Success converting Antimony file to SBML")
    # SBML: Annotation
    # Import SBML file
    sbml_reader = libsbml.SBMLReader()
    sbml_doc = sbml_reader.readSBML(sbml_file_name)
    sbml_model = sbml_doc.getModel()
    # Set species annotations
    write_species_annotations(sbml_model, species)
    # Set compartments annotations
    # write_compartments_annotations(sbml_model, compartments)
    # Export the annotated SBML file
    writer = libsbml.SBMLWriter()
    writer.writeSBML(sbml_doc, sbml_file_name)

    # SBML: Compilation
    # Import annotated SBML file
    sys.path.insert(0, os.path.abspath(model_output_dir))
    sbml_reader = libsbml.SBMLReader()
    sbml_doc = sbml_reader.readSBML(sbml_file_name)
    sbml_model = sbml_doc.getModel()
    sbml_importer = amici.SbmlImporter(sbml_file_name)
    const_params = [params.getId() for params in sbml_model.getListOfParameters()]
    # Compile
    sbml_importer.sbml2amici(sbml_model_name, model_output_dir, verbose=args.verbose)
    if args.verbose: print("SPARCED: Sucess compiling the model")


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', '--antimony',     default="SPARCED",          help="name of the antimony model")
    parser.add_argument('-b', '--sbml',         default="SPARCED",          help="name of the SBML model")
    parser.add_argument('-c', '--compartments', default="Compartments.txt", help="name of the compartments' file")
    parser.add_argument('-i', '--inputdir',     default="./../data/",       help="relative path to input files directory")
    parser.add_argument('-m', '--stoichmatrix', default="StoicMat.txt",     help="name of the stoichiometric matrix' file")
    parser.add_argument('-o', '--outputparams', default="ParamsAll.txt",    help="name of the output parameters' file")
    parser.add_argument('-r', '--ratelaws',     default="Ratelaws.txt",     help="name of the rate laws' file")
    parser.add_argument('-s', '--species',      default="Species.txt",      help="name of the species' file")
    parser.add_argument('-v', '--verbose',       action='store_true',       help="display additional details during execution")
    return(parser.parse_args())


if __name__ == '__main__':
    args = parse_args()
    # Update path to files
    f_compartments = args.inputdir + args.compartments
    f_stoichmat = args.inputdir + args.stoichmatrix
    f_output_params = args.inputdir + args.outputparams
    f_ratelaws = args.inputdir + args.ratelaws
    f_species = args.inputdir + args.species
    # Create model
    create_model(args.antimony,args.sbml,f_compartments,f_stoichmat,
                 f_output_params,f_ratelaws,f_species,args.verbose)
