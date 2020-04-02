"""AMICI-generated module for model SPARCEDv6"""

import amici

# Ensure we are binary-compatible, see #556
if '0.10.20' != amici.__version__:
    raise RuntimeError('Cannot use model SPARCEDv6, generated with AMICI '
                       'version 0.10.20, together with AMICI version'
                       f' {amici.__version__} which is present in your '
                       'PYTHONPATH. Install the AMICI package matching the '
                       'model version or regenerate the model with the AMICI '
                       'currently in your path.')

from SPARCEDv6._SPARCEDv6 import *

__version__ = '0.1.0'
