#!/usr/bin/bash

# NOTE: directories are referenced here with the expectation that this script
# will be run from on directory up, using tests/run_tests.sh

set -e  # stop on error
set -u  # raise error if variable is unset
set -o pipefail  # fail if any prior step failed

# run pycodestyle on new python scripts used in unit and functional testing
pycodestyle funclib.py
pycodestyle segment.py
pycodestyle segment_multiple.py
pycodestyle tests/unit_tests.py

# run unit testing
python -m unittest tests/unit_tests.py

# run functional testing
tests/func_tests.sh
