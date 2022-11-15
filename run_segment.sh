#!/usr/bin/bash

set -e  # stop on error
set -u  # raise error if variable is unset
set -o pipefail  # fail if any prior step failed

python segment.py --input_file 'IBIS.granulation.aligned.25Apr2019.seq56.sav' --skimage_method 'li' --plot_intermed True --out_file 'output.fits' --out_dir 'output/'