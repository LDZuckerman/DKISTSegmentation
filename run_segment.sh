#!/usr/bin/bash

set -e  # stop on error
set -u  # raise error if variable is unset
set -o pipefail  # fail if any prior step failed

# to run on IBIS data
#python segment.py --input_file 'data/IBIS.granulation.aligned.25Apr2019.seq56.sav' --skimage_method 'li' --plot_intermed True --out_file 'output.fits' --out_dir 'output_IBIS/'

# to run on DKIST data
python segment.py --input_file 'data/dkist.cont789nm.scaled.fits' --skimage_method 'li' --plot_intermed True --out_file 'output.fits' --out_dir 'output_DKIST/'