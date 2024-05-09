#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse


def parse_args():
    """Retrieve and parse arguments necessary for model creation
    
    Returns:
        A namespace populated with all the attributes.
    """
    
    parser = argparse.ArgumentParser()
    # Capital letter => file & path names + flagD
    parser.add_argument('-D', '--deterministic', action="store_false",
                        help="don't run simulation in deterministic mode")
    parser.add_argument('-I', '--inputdir',      default="./../data/core/",
                        help="relative path to the core input directory (model creation)")
    parser.add_argument('-J', '--inputsim',      default="./../data/simulation/",
                        help="relative path to the simulation input directory")
    parser.add_argument('-O', '--outputdir',     default="./../models/SPARCED-standard/",
                        help="relative path to the output directory")
    parser.add_argument('-P', '--outputparams',  default="out_Parameters.txt",
                        help="desired name for the output parameters file")
    parser.add_argument('-W', '--outputsim',    default="./../results/New-Simulation/",
                        help="relative path to the directory where simulation \
                              results will be saved")
    # Lowcase letter => running options
    parser.add_argument('-c', '--compound',      default=None,
                        help="name of the compound (for example a drug) you wish to alter")
    parser.add_argument('-d', '--dose',          default=0.0,
                         help="concentration of the altered compound (nM)")
    parser.add_argument('-n', '--name',          default="SPARCED",
                        help="desired name for the generated model AND name of \
                              the input model")
    parser.add_argument('-s', '--simulation',    default="GrowthStim",
                        help="desired name for the simulation output files")
    parser.add_argument('-p', '--popsize',    default=1,
                        help="desired cell population size for the simulation")
    parser.add_argument('-t', '--time',          default=1.0,
                        help="desired virtual duration of the simulation (h)")
    parser.add_argument('-v', '--verbose',       action='store_false',
                        help="don't display additional details during execution")
    parser.add_argument('-w', '--wild',          action='store_true',
                        help="UNDER CONSTRUCTION\nrunning wild (without SPARCED \
                              hard-coded values/behaviors")
    parser.add_argument('-x', '--exchange',      default=30,
                         help="timeframe between modules information exchange \
                               during the simulation")
    return(parser.parse_args())

