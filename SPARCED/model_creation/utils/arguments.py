#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse


def parse_args():
    """
    Retrieve and parse arguments
    :return: a namespace populated with all the attributes
    :rtype: [Namespace]
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', '--antimony',     default="ant_SPARCED",
                        help="desired name for the generated antimony model")
    parser.add_argument('-b', '--sbml',         default="sbml_SPARCED",
                        help="desired name for the generated SBML model")
    parser.add_argument('-c', '--compartments', default="Compartments.txt",
                        help="name of the compartments file")
    parser.add_argument('-i', '--inputdir',     default="./../data/core/",
                        help="relative path to the input directory")
    parser.add_argument('-m', '--stoichmatrix', default="StoicMat.txt",
                        help="name of the stoichiometric matrix file")
    parser.add_argument('-o', '--outputdir',    default="./../models/SPARCED-standard/",
                        help="relative path to the output directory")
    parser.add_argument('-p', '--outputparams', default="out_Parameters.txt",
                        help="desired name for the output parameters file")
    parser.add_argument('-r', '--ratelaws',     default="Ratelaws.txt",
                        help="name of the ratelaws file")
    parser.add_argument('-s', '--species',      default="Species.txt",
                        help="name of the species file")
    parser.add_argument('-v', '--verbose',      action='store_true',
                        help="display additional details during execution")
    return(parser.parse_args())

