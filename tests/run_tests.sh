#!/usr/bin/bash

# run pycodestyle on new python scripts used in unit and functional testing
pycodestyle funclib.py
pycodestyle segment.py
pycodestyle test_funclib.py

# run unit testing
python -m unittest tests/unit_tests.py

# run functional testing
tests/func_tests.sh
