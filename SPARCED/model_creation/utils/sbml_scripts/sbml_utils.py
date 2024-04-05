#/usr/bin/env python
# -*- coding: utf-8 -*-


def write_compartments_annotations(f, comp):
    """Set compartments annotations using the last column of the compartments input file"""
    for row in comp:
        # For some unexplained reason setAnnotation throws an error
        # Some investigation in the official documentation shows it should work...
        f.getCompartment(row[0]).setAnnotation(row[2])

def write_species_annotations(f, spec):
    """Set species annotations using the last column of the species input file"""
    for i, row in enumerate(spec[1:]): # Skip header
        annot = ""
        for col in range(4, len(row)):
            aa = str(row[col].strip())
            if not (aa == "nan" or aa == ""): annot += " " + row[col]
        f.getSpecies(row[0]).setAnnotation(annot)
