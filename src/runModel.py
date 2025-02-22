#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Filename: runModel.py
Created: 2025-02-19
Author(s):
Description:
"""


#<-----------------------------Import Packages------------------------------->
import os
import argparse

from run_sparced import RunSPARCED
from utils.data_handler import Results
from utils.config_manager import ConfigManager
from utils.model_handler import PerturbationHandler
from utils.model_factory import create_model_handler
from singe.singe_engine import SINGEEngine

parser = argparse.ArgumentParser(description='Provide arguments to build the SPARCED model')
parser.add_argument('--config_path', '-p', metavar='config', 
                    help='the path to the configuration file')
args = parser.parse_args()

if args.config_path is None:
    print("ERROR: missing arguments. Need to pass --config_path. Use -h for help.")

class Simulation:
    """
    Represents the simulation entrypoint.
    """
    def __init__(self, config_path):
        """
        Initialize the SimulationEntry class.
        """
        self.config_manager = ConfigManager(config_path)

    def run(self):
        """
        Run the simulation.
        """
        # Get the gene regulation and Omics data paths from the config file
        gene_regulation = self.config_manager.get("simulation.data.gene_regulation")
        omics_data = self.config_manager.get("simulation.data.omics")

        # First, instantiate the ODE model handler
        model_handler = create_model_handler(self.config_manager.get("simulation.model.sbml"),
                                             self.config_manager.get("simulation.model.amici"))

        # Apply perturbations to the model
        perturbants = self.config_manager.get("simulation.perturbations", {})

        ode_model = PerturbationHandler(model_handler)
        ode_model.apply_perturbations(perturbants)

        # Get Model Exceptions from the config file
        model_exceptions = self.config_manager.get("simulation.model.exceptions", None)

        # Get true/false statement on using hybrid setting
        solver_flag = self.config_manager.get("simulation.protocol.hybrid", False)

        # Set the overall simulation time
        duration = self.config_manager.get("simulation.protocol.duration", 0)
        exchange = self.config_manager.get("simulation.exchange", 30)

        # Future: Integrate components of pertrubed model into SINGE model and remove from function. 
        singe_model = SINGEEngine(gene_regulation, 
                                  omics_data, 
                                  solver_flag, 
                                  duration, 
                                  exchange, 
                                  model_exceptions)

        # Instatiate the SPARCED class, run simulation, return results
        results = RunSPARCED(ode_model, singe_model, duration, exchange)

        # Save results to file. Note; save_data method takes !tuple! datatype as input
        Results(os.path.join(self.config_manager.get("simulation.results.directory"),
                self.config_manager.get("simulation.results.filename"))).save_data(results)

if __name__ == "__main__":
    Simulation(args.config_path).run()

# th = args.time
# Vn = float(args.Vn)
# Vc = float(args.Vc)
# outfile = args.outfile
# ts = 30


# if flagD == 0:
#     flagWr = 1
#     nmxlsfile = outfile
    
#     sys.path.insert(0, os.path.abspath(model_output_dir))

#     species_sheet = np.array([np.array(line.strip().split("\t")) for line in open('Species.txt', encoding='latin-1')])

#     species_initializations = []
#     for row in species_sheet[1:]:
#         species_initializations.append(float(row[2]))
#     species_initializations = np.array(species_initializations)

#     model_module = importlib.import_module(model_name)
#     model = model_module.getModel()
#     solver = model.getSolver() # Create solver instance
#     solver.setMaxSteps = 1e10
#     model.setTimepoints(np.linspace(0,ts)) # np.linspace(0, 30) # set timepoints

#     xoutS_all, xoutG_all, tout_all = RunSPARCED(flagD,th,species_initializations,[],sbml_file,model)


# elif flagD == 1:
#     flagWr = 1
#     nmxlsfile = outfile

#     sys.path.insert(0, os.path.abspath(model_output_dir))
#     species_sheet = np.array([np.array(line.strip().split("\t")) for line in open('Species.txt', encoding='latin-1')])

#     species_initializations = []
#     for row in species_sheet[1:]:
#         species_initializations.append(float(row[2]))

#     species_initializations = np.array(species_initializations)
#     species_initializations[np.argwhere(species_initializations <= 1e-6)] = 0.0

#     model_module = importlib.import_module(model_name)
#     model = model_module.getModel()
#     solver = model.getSolver()          # Create solver instance
#     solver.setMaxSteps = 1e10
#     model.setTimepoints(np.linspace(0,ts)) # np.linspace(0, 30) # set timepoints

#     xoutS_all, xoutG_all, tout_all = RunSPARCED(flagD,th,species_initializations,[],sbml_file,model)

