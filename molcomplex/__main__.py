#!/usr/bin/env python

###############################################.
#          __main__ file for the code         #
###############################################.

from __future__ import absolute_import

import sys
from molcomplex import molcomplex

# If we are running from a wheel, add the wheel to sys.path
# This allows the usage python pip-*.whl/pip install pip-*.whl

if __package__ != 'molcomplex':
    print('molcomplex is not installed! Use: pip install molcomplex (anywhere, using a terminal) or python setup.py install (from the downloaded /molcomplex/molcomplex folder).')

if __name__ == '__main__':
    molcomplex.main()
    sys.exit()
