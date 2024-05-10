"""AMICI-generated module for model SPARCED_standard"""

import amici

# Ensure we are binary-compatible, see #556
if '0.11.12' != amici.__version__:
    raise RuntimeError('Cannot use model SPARCED_standard, generated with AMICI '
                       'version 0.11.12, together with AMICI version'
                       f' {amici.__version__} which is present in your '
                       'PYTHONPATH. Install the AMICI package matching the '
                       'model version or regenerate the model with the AMICI '
                       'currently in your path.')

from SPARCED_standard._SPARCED_standard import *

__version__ = '0.1.0'
