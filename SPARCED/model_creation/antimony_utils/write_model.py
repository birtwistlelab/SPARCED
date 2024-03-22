#!/usr/bin/env python
# -*- coding: utf-8 -*-

import antimony
import argparse

from antimony_utils import *

def antimony_write_model(antimony_model_name,f_compartments,f_stoichmat,
                         f_output_params,f_ratelaws,f_species):
    """
    Generate an Antimony file based on given data

    :param antimony_model_name: name for the generated antimony model
    :type antimony_model_name: [str]
    :param f_compartments: compartments & volumes file
    :type f_compartments: [str]
    :param f_stoichmat: stoichiometric matrix file
    :type f_stoichmat: [str]
    :param f_output_params: output parameters file
    :type f_output_params: [str]
    :param f_ratelaws: ratelaws file
    :type f_ratelaws: [str]
    :param f_species: species file
    :type f_species: [str]
    :return: The file name of the generated Antimony model + compartments and species
    :rtype: ([str], list[str], list[str])
    
    """
    antimony_file = antimony_model_name + ".txt"

    with open(antimony_file,"w") as antimony_model:
        # Write file's header
        antimony_model.write("# PanCancer Model by Birtwistle Lab\n")
        antimony_model.write("model {antimony}()\n\n"
                             .format(antimony=antimony_model_name))
        # Write compartments, species and reactions
        compartments, volumes, species = antimony_init(f_compartments,
                                                       f_species)
        antimony_write_compartments(antimony_model,compartments)
        antimony_write_species(antimony_model,species)
        param_names, param_vals = antimony_write_reactions(antimony_model,
                                                           f_ratelaws,
                                                           f_stoichmat,
                                                           f_output_params)
        # Write initial conditions
        antimony_write_init_compartments(antimony_model,compartments,volumes)
        antimony_write_init_species(antimony_model,species)
        # Write reaction parameters
        antimony_write_init_reactions(antimony_model,param_names,param_vals)
        # Write other declarations and unit definitions
        antimony_terminal(antimony_model)
        antimony_model.write("\nend")
    return(antimony_file, compartments, species)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', '--antimony',     default="SPARCED",
    help="Give a name to the antimony model")
    parser.add_argument('-c', '--compartments', default="Compartments.txt")
    parser.add_argument('-i', '--inputdir',     default="./../../data/")
    parser.add_argument('-m', '--stoichmatrix', default="StoicMat.txt")
    parser.add_argument('-o', '--outputparams', default="ParamsAll.txt")
    parser.add_argument('-r', '--ratelaws',     default="Ratelaws.txt")
    parser.add_argument('-s', '--species',      default="Species.txt")
    return(parser.parse_args())

# --------------------------------------------------------------------------- #

if __name__ == '__main__':
    args = parse_args()
    # Update path to files
    f_compartments = args.inputdir + args.compartments
    f_stoichmat = args.inputdir + args.stoichmatrix
    f_output_params = args.inputdir + args.outputparams
    f_ratelaws = args.inputdir + args.ratelaws
    f_species = args.inputdir + args.species
    # Write model
    antimony_write_model(args.antimony,f_compartments,f_stoichmat,
                         f_output_params,f_ratelaws,f_species)

