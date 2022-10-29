#!/usr/bin/bash

set -e  # stop on error
set -u  # raise error if variable is unset
set -o pipefail  # fail if any prior step failed

# run pycodestyle on new python scripts used in unit and functional testing
pycodestyle funclib.py
pycodestyle segment_skimage.py
pycodestyle test_funclib.py

# run unit testing
python -m unittest test_funclib.py

# run functional testing
./func_tests.sh
