#!/usr/bin/env python
# -*- coding: utf-8 -*-

import libsbml
import numpy as np

from compilation.sbml_scripts.sbml_utils import *


def sbml_annotate_model(sbml_file_name: str, species: np.ndarray,
                        compartments: np.ndarray) -> None:
    """Annotate species and compartments of the given SBML model

    Arguments:
        sbml_file_name: The SBML file path.
        species: Content of the species input file.
        compartments: Content of the compartments input file.

    Returns:
        Nothing.
    """

    # Import SBML file
    sbml_reader = libsbml.SBMLReader()
    sbml_doc = sbml_reader.readSBML(sbml_file_name)
    sbml_model = sbml_doc.getModel()
    # Set species annotations
    write_species_annotations(sbml_model, species)
    # Set compartment annotations
    write_compartments_annotations(sbml_model, compartments)
    # Export the annotated SBML file
    writer = libsbml.SBMLWriter()
    writer.writeSBML(sbml_doc, sbml_file_name)

