#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Filename: model_handler.py
Created: 2025-02-19
Author(s):
Description:
"""


#<-----------------------------Import Packages------------------------------->
import importlib
import amici
import libsbml
import numpy as np
import tellurium as te

#<------------------------------Parent Class--------------------------------->
# Parent class
class ModelHandler:
    """
    Loads and performs CRUD operations on SBML and AMICI models.
    """
    def __init__(self, model_path):
        self.model_path = model_path
        self.model = None  # Placeholder for the actual model object

    def load_model(self):
        """Method to load a model - must be implemented in child classes."""
        raise NotImplementedError("Subclasses must implement this method.")

    def modify_parameter(self, param_name, value):
        """Modify a model parameter - must be implemented in child classes."""
        raise NotImplementedError("Subclasses must implement this method.")
    
    def modify_species(self, species_name, value):
        """Modify a model species - must be implemented in child classes."""
        raise NotImplementedError("Subclasses must implement this method.")

    def simulate(self):
        """Run model simulation - must be implemented in child classes."""
        raise NotImplementedError("Subclasses must implement this method.")

#<------------------------------Child Classes-------------------------------->
class SBMLModelHandler(ModelHandler):
    """
    Handles loading and simulating SBML models using libSBML.
    """
    def __init__(self, model_path):
        super().__init__(model_path)
        self.load_model()

    def load_model(self):
        """Loads an SBML model using libSBML."""
        reader = libsbml.SBMLReader()
        document = reader.readSBML(self.model_path)
        self.sbml_model = document.getModel()
        if self.sbml_model is None:
            raise ValueError(f"Failed to load SBML model from {self.model_path}")

    def modify_parameter(self, param_name, value):
        """Modify a parameter in the SBML model."""
        param = self.model.getParameter(param_name)
        if param:
            param.setValue(value)
        else:
            raise KeyError(f"Parameter '{param_name}' not found in SBML model.")
       
    def modify_species(self, species_name, value):
        """No need to override this method for SBML models."""

    def module_exchange(self, exchange=30):
        """No need to override this method for SBML models."""

    def simulate(self):
        """Simulation logic would go here (Tellurium, COPASI, etc.)"""
        print("Running SBML simulation... (implement with Tellurium or another solver)")

    def getVolume(self, compartment):
        """
        retrieves the volume of a particular component by name.
        """
        return self.model.getCompartment(compartment).getVolume()

# Child class for AMICI
class AMICIModelHandler(SBMLModelHandler):
    """
    Handles loading and simulating AMICI models.
    """
    def __init__(self, model_path):
        super().__init__(model_path)

        self.load_model()

    def load_model(self):
        """Loads an AMICI model using importlib."""
        self.amici_model = importlib.import_module(self.model_path)
        self.model = self.amici_model.getModel()

    def modify_parameter(self, param_name, value):
        """Modify a parameter in the AMICI model."""
        param_index = self.model.getParameterIndex(param_name)
        if param_index >= 0:
            self.model.setParameterById(param_index, value)
        else:
            raise KeyError(f"Parameter '{param_name}' not found in AMICI model.")
        
    def modify_species(self, species_name, value):
        """Modify a species concentration in the AMICI model."""
        species_index = self.model.getSpeciesIndex(species_name)
        if species_index >= 0:
            self.model.setInitialStatesById(species_index, value)
        else:
            raise KeyError(f"Species '{species_name}' not found in AMICI model.")

    def module_exchange(self, exchange=30):
        """Set the amount of time the AMICI model simulates before exchanging 
        information with the SINGE engine."""
        self.model.setTimepoints(np.linespace(0, exchange)) # np.linspace default points set to 50

    def simulate(self):
        """Runs the AMICI model simulation."""
        solver = self.model.getSolver()
        solver.setMaxSteps(1e10)
        results = amici.runAmiciSimulation(self.model, solver)
        print("AMICI simulation complete.")
        return results

class TelluriumModelHandler(SBMLModelHandler):
    """
    Handles loading and simulating SBML models using Tellurium.
    """
    def __init__(self, model_path):
        super().__init__(model_path)
        self.load_model()

    def load_model(self):
        """Loads an SBML model using Tellurium."""
        self.model = te.loadSBMLModel(self.model_path)
        if self.model is None:
            raise ValueError(f"Failed to load SBML model from {self.model_path}")

    def modify_parameter(self, param_name, value):
        """Modify a parameter in the SBML model."""
        self.model[param_name] = value

    def modify_species(self, species_name, value):
        """Modify a species in the SBML model."""
        self.model[species_name] = value

    def module_exchange(self, exchange=30):
        """Set the amount of time the Tellurium model simulates before exchanging 
        information with the SINGE engine."""
        timestep_results = self.model.simulate(0, exchange, points = 50) # Points set to match AMICI
        return timestep_results

    def simulate(self):
        """Runs the Tellurium model simulation."""
        results = self.model.simulate()
        print("Tellurium simulation complete.")
        return results

class PerturbationHandler:
    """
    Handles applying perturbations to a model instance.
    """
    def __init__(self, model_handler):
        """
        :param model_handler: An instance of a ModelHandler subclass.
        :param perturbations: A dictionary of perturbations to apply to the model.
        """
        self.model_handler = model_handler

    def apply_perturbations(self, perturbations):
        """
        Applies perturbations to the model.
        """

        for param, value in perturbations.get("parameters", {}).items():
            self.model_handler.modify_parameter(param, value)

        for species, value in perturbations.get("species", {}).items():
            self.model_handler.modify_species(species, value)
