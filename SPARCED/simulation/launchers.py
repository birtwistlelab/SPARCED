#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys

import numpy as np

from utils.arguments import parse_args
from simulation.run_model import run_experiment


def launch_experiment_simulation(ligands: np.ndarray) -> None:
    """
    Small routine to process parsed arguments and call run_model
    """
    args = parse_args()
    # Process arguments
    is_SPARCED = not args.wild # if it's not wild then it's SPARCED
    f_species = args.inputdir + args.species
    model_path = args.outputdir + "amici_" + args.name
    sbml_model = args.outputdir + "sbml_" + args.name + ".xml"
    try:
        assert args.population > 0
    except:
        print("{model_name}-{simulation_name}: Cell population size should be \
               superior to zero (0). Current population size is: {size}.\n \
               Aborting now.".format(model_name=args.name,
                                     simulation_name=args.simulation,
                                     size = args.population))
        sys.exit(0)
    run_experiment(args.name, model_path, args.simulation, sbml_model,
                   args.deterministic, args.exchange, args.population,
                   args.time, f_species, ligands, args.verbose, is_SPARCED,
                   args.compound, args.dose)


if __name__ == '__main__':
    launch_experiment_simulation()

