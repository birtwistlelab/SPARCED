"""AMICI-generated module for model SPARCEDo4a_v1"""

import amici

# Ensure we are binary-compatible, see #556
if '0.11.18' != amici.__version__:
    raise RuntimeError('Cannot use model SPARCEDo4a_v1, generated with AMICI '
                       'version 0.11.18, together with AMICI version'
                       f' {amici.__version__} which is present in your '
                       'PYTHONPATH. Install the AMICI package matching the '
                       'model version or regenerate the model with the AMICI '
                       'currently in your path.')

from SPARCEDo4a_v1._SPARCEDo4a_v1 import *

__version__ = '0.1.0'
