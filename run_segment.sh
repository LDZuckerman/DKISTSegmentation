#!/usr/bin/bash
set -e  # stop on error
set -u  # raise error if variable is unset
set -o pipefail  # fail if any prior step failed

# to run on IBIS data
python segment.py --data_path 'data/IBIS' --skimage_method 'li' --plot_intermed True --out_file 'segmented_data' --vel_comparison_file 'velocity_comparison' --out_dir 'example_outputs/IBIS/'

# to run on DKIST data
python segment.py --data_path 'data/DKIST' --skimage_method 'li' --plot_intermed True --out_file 'segmented_data' --out_dir 'example_outputs/DKIST/'

