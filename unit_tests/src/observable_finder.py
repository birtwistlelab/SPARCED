#!/usr/bin/env python     
# -*- coding: utf-8 -*-

import libsbml
import argparse

# Parse the command line arguments
parser = argparse.ArgumentParser(description='Generate the observable formula for a given observable')
parser.add_argument('--observable','-o', type=str, help='The observable to generate the formula for')
parser.add_argument('--model', '-m',type=str, help='The path to the SBML model file')
args = parser.parse_args()


# Read the SBML model
reader = libsbml.SBMLReader()
document = reader.readSBML(args.model)
model = document.getModel()

def _get_observable_formula(model: libsbml.Model, observable: str):
    """Extract the observable formula as the sum of the product 
        of the species concentration and the compartment volume ratio.
    Input:
        model: libsbml.Model: The SBML model
        observable: str: The observable to generate the formula for

    Output:
        str: The observable formula
    """
    compartments = model.getListOfCompartments()

    cytoplasm_volume = [compartment.getVolume() for compartment in compartments if compartment.getId() == 'Cytoplasm'][0]

    species = model.getListOfSpecies()

    observable_formula = []

    for index, specie in enumerate(species):
    
        specie_name = species.get(index).getId()
    
        if observable in specie_name and 'm_' not in specie_name:

            specie_compartment = species.get(index).getCompartment()
            
            specie_compartment_volume = [compartment.getVolume() for compartment in compartments if compartment.getId() == specie_compartment][0]

            Volume_ratio = specie_compartment_volume / cytoplasm_volume

            observable_formula.append(f'{specie_name} * {Volume_ratio}')

    observable_formula = ' + '.join(observable_formula)
    
    print(observable_formula)

_get_observable_formula(model, args.observable)