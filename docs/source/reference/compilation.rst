===============================================================================
Model Creation
===============================================================================

General
-------------------------------------------------------------------------------

Main function for model creation

.. autofunction:: create_model.create_model()

.. autofunction:: create_model.launch_model_creation()

AMICI
-------------------------------------------------------------------------------

Model compilation

.. autofunction:: model_compilation.compile_sbml_to_amici()

Utilities

.. autofunction:: amici_utils.define_observables()

.. autofunction:: amici_utils.extract_amici_model_name()


Antimony
-------------------------------------------------------------------------------

Model writing

.. autofunction:: model_writing.antimony_write_model()

General Utilities

.. autofunction:: antimony_utils.extract_antimony_model_name()

Writing Utilities

.. autofunction:: antimony_write.antimony_write_constant_variables()

.. autofunction:: antimony_write.antimony_write_compartments_names()

.. autofunction:: antimony_write.antimony_write_reactions()

.. autofunction:: antimony_write.antimony_write_species_names()

.. autofunction:: antimony_write.antimony_write_unit_definitions()

Writing Initial Conditions Utilities

.. autofunction:: antimony_write_IC.antimony_write_compartments_IC()

.. autofunction:: antimony_write_IC.antimony_write_reactions_IC()

.. autofunction:: antimony_write_IC.antimony_write_species_IC()


SBML
-------------------------------------------------------------------------------

Model annotation

.. autofunction:: model_annotation.sbml_annotate_model()

Utilities

.. autofunction:: sbml_utils.write_compartments_annotations()

.. autofunction:: sbml_utils.write_species_annotations()

