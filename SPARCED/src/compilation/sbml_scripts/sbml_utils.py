#/usr/bin/env python
# -*- coding: utf-8 -*-

import libsbml


def write_compartments_annotations(f_sbml, compartments):
    """Set compartments annotations in the given SBML file

    Note:
        The first row is considered as a header, and hence it is skipped.
        Compartment names should be located on the first column of the array.
        Annotations should be located on the third column of the array.

    Arguments:
        f_sbml: The loaded SBML file.
        compartments: Content of the compartments input file.

    Returns:
        Nothing.
    """

    for row in compartments[1:]:
        f_sbml.getCompartment(row[0]).setAnnotation(row[2])

def write_species_annotations(f_sbml, species):
    """Set species annotations in the given SBML file

    Note:
        The first row is considered as a header, and hence it is skipped.
        Annotations should be located between the 4th and the last column of
        the array.

    Arguments:
        f_sbml: The loaded SBML file.
        species: Content of the species input file.

    Returns:
        Nothing.
    """

    for i, row in enumerate(species[1:]):
        annotation = ""
        for col in range(4, len(row)):
            specie_annotation = str(row[col].strip())
            if not (specie_annotation == "nan" or specie_annotation == ""):
                annotation += " " + row[col]
        f_sbml.getSpecies(row[0]).setAnnotation(annotation)

