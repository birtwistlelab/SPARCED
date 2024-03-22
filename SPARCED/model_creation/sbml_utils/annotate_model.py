#!/usr/bin/env python
# -*- coding: utf-8 -*-

import libsbml

def sbml_annotate_model(sbml_file_name, species, compartments):
    # Import SBML file
    sbml_reader = libsbml.SBMLReader()
    sbml_doc = sbml_reader.readSBML(sbml_file_name)
    sbml_model = sbml_doc.getModel()
    # Set species annotations
    write_species_annotations(sbml_model, species)
    # Set compartment annotations
    # write_compartments_annotations(sbml_model, compartments)
    # Export the annotated SBML file
    writer = libsbml.SBMLWriter()
    writer.writerSBML(sbml_doc, sbml_file_name)

